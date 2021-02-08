using BioSequences
using CodecZlib
using FASTX

mutable struct BarcodeCorrectionCounter
    candidates::Vector{LongDNASeq}
    likelihoods::Vector{Float64}
    likelihood_sum::Float64
end

UmiCounts = Dict{LongSequence,Int}
SeqOrNothing = Union{LongSequence,Nothing}


@inline function get_barcode_seq(record::FASTQ.Record)::LongSequence
    return @view FASTX.sequence(record)[1:16]
end


@inline function get_barcode_quality(record::FASTQ.Record)::Vector{UInt8}
    return @view FASTQ.quality(record)[1:16]
end


@inline function get_umi_seq(record::FASTQ.Record)::LongSequence
    return @view FASTX.sequence(record)[17:26]
end


@inline function get_barcodes_seqs_as_str(record::FASTQ.Record)::Tuple{String, String}
    sequence = FASTX.sequence(record)
    return (string(@view sequence[1:16]), string(@view sequence[17:26]))
end


@inline function auto_gzopen(filepath)
    return endswith(filepath, ".gz") ? GzipDecompressorStream(open(filepath)) :
           open(filepath)
end


"""
    get_barcode_dist(fastq_fn, seq_whitelist)

Obtains the barcode distribution counting only barcodes present in whitelist
"""
function get_barcode_dist(
    fastq_fn::String,
    barcode_whitelist::Set{S},
)::AbstractDict{S,Float64} where {S<:LongSequence}
    frequencies = Dict{S,Float64}()

    reader = FASTQ.Reader(auto_gzopen(fastq_fn))
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        barcode = get_barcode_seq(record)
        if barcode in barcode_whitelist
            frequencies[barcode] = get(frequencies, barcode, 0) + 1
        end
    end

    # normalize
    nrecords = sum(values(frequencies))
    for (k, v) in frequencies
        frequencies[k] = v / nrecords
    end

    return frequencies
end


"""
    cell_barcode_correction(record, barcode_whitelist, bc_dist[, confidence_threshold])

Performs barcode correction for a single barcode testing each barcode at Hamming-distance 1
"""
function cell_barcode_correction(
    record::FASTQ.Record,
    barcode_whitelist::Set{S},
    bc_dist::AbstractDict{S,Float64};
    confidence_threshold::Float64 = 0.975,
)::Union{Nothing,LongSequence} where {S<:LongSequence}
    dna_bases = [DNA_A, DNA_C, DNA_G, DNA_T]
    sequence = get_barcode_seq(record)
    quality = get_barcode_quality(record)
    counter = BarcodeCorrectionCounter([], [], 0)
    sizehint!(counter.likelihoods, length(sequence) * length(dna_bases) - 1)

    for (i, nucl) in enumerate(sequence)
        # test each possible change
        for dna_base in dna_bases
            if nucl != dna_base
                sequence[i] = dna_base

                p_bc = if sequence in keys(bc_dist)
                    # observed frequency if observed
                    bc_dist[sequence]
                elseif sequence in barcode_whitelist
                    # approximate pseudocount if unobserved but on whitelist
                    1 / length(barcode_whitelist)
                else
                    nothing
                end

                if !isnothing(p_bc)
                    p_edit = 10^(-quality[i] / 10)
                    likelihood = p_bc * p_edit
                    push!(counter.candidates, sequence)
                    push!(counter.likelihoods, likelihood)
                    counter.likelihood_sum += likelihood
                end
            end
        end
        sequence[i] = nucl
    end

    if length(counter.likelihoods) > 0
        max_value, max_idx = findmax(counter.likelihoods)

        if max_value / counter.likelihood_sum > confidence_threshold
            return counter.candidates[max_idx]
        end
    end

    return nothing
end


"""
    load_whitelist_barcodes(filepath)

Loads whitelist barcodes as a set
"""
function load_whitelist_barcodes(filepath::String)::Set{LongSequence}
    f = auto_gzopen(filepath)
    whitelist_barcodes = [LongDNASeq(line) for line in eachline(f)]
    close(f)

    return Set{LongSequence}(whitelist_barcodes)
end


"""
    correct_barcodes_and_count_umis(fastq_fn, barcode_whitelist, barcode_dist)

Performs barcode correction and counts UMI by barcode in the same pass
"""
function correct_barcodes_and_count_umis(
    fastq_fn::String,
    barcode_whitelist::Set{S},
    barcode_dist::AbstractDict{S,Float64},
)::Tuple{
    AbstractDict{S,SeqOrNothing},
    AbstractDict{SeqOrNothing,UmiCounts},
} where {S<:LongSequence}
    barcode_corrections = Dict{S,SeqOrNothing}()
    umi_counts = Dict{SeqOrNothing,UmiCounts}()

    reader = FASTQ.Reader(auto_gzopen(fastq_fn))
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)
        barcode = get_barcode_seq(record)

        # use the dictionary as cache to don't compute the same again
        if !haskey(barcode_corrections, barcode)
            # if barcode is in whitelist, it's assumed that it doesn't need correction
            if !(barcode in barcode_whitelist)
                barcode_corrections[barcode] =
                    cell_barcode_correction(record, barcode_whitelist, barcode_dist)
                # override the variable with the corrected barcode
                barcode = barcode_corrections[barcode]
            else
                # save the barcode when it's in whitelist anyway
                barcode_corrections[barcode] = barcode
            end
        end

        # increase UMI counter
        if !haskey(umi_counts, barcode)
            umi_counts[barcode] = UmiCounts()
        end
        umi = get_umi_seq(record)
        umi_counts[barcode][umi] = get(umi_counts[barcode], umi, 0) + 1
    end

    return barcode_corrections, umi_counts
end


"""
    correct_umis(umi_counts)

Performs UMI corrections testing UMIs at Hamming-distance 1
"""
function correct_umis(umi_counts::UmiCounts)::Dict{LongSequence,LongSequence}
    umi_corrections = Dict{LongSequence,LongSequence}()
    dna_bases = [DNA_A, DNA_C, DNA_G, DNA_T]
    test_umi = dna""

    for (umi, count) in umi_counts
        # make the original as the best by default
        best_dest_umi = umi
        best_dest_count = count
        for (i, nucl) in enumerate(umi)
            for dna_base in dna_bases
                if nucl != dna_base
                    copy!(test_umi, umi)
                    test_umi[i] = dna_base

                    test_count = get(umi_counts, test_umi, 0)
                    # if tested umi has greater count or equal count but is lexicographically larger,
                    # make it the best candidate
                    if test_count > best_dest_count ||
                       (test_count == best_dest_count && test_umi > best_dest_umi)
                        best_dest_umi = copy(test_umi)
                        best_dest_count = test_count
                    end
                end
            end
        end
        if best_dest_umi != umi
            umi_corrections[umi] = best_dest_umi
        end
    end

    return umi_corrections
end
