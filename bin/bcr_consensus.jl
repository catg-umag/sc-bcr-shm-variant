#!/usr/bin/env julia
using ArgParse
using BioAlignments
using BioSequences
using DataStructures
using CSV
using FASTX
using XAM


include("../lib/julia/dna_counter.jl")
include("../lib/julia/others.jl")


function main()
    args = parse_arguments()

    # load data
    reference = get_reference(args["reference"], args["name"])
    info_dict = Dict(
        row.id => Dict(:cell => Symbol(row.cell), :umi => Symbol(row.umi))
        for row in CSV.File(args["index"])
    )

    records = load_records(args["bamfile"], info_dict, args["name"])

    open(args["output"], "w") do f
        write(f, "cell,umi,nrecords,ref_coverage,consensus,filled_consensus,depths\n")

        for (cell, cell_records) in records
            for umi in keys(cell_records)
                nrecords = length(cell_records[umi])
                consensus, depths =
                    make_consensus!(cell_records[umi], reference, fill_with_ref = false)
                coverage = length(filter(x -> x > 0, depths)) / length(depths)

                # fill Ns in consensus with reference
                filled_consensus = copy(consensus)
                for i in eachindex(filled_consensus)
                    if filled_consensus[i] == DNA_N
                        filled_consensus[i] = reference[i]
                    end
                end

                # write result
                write(
                    f,
                    "$(string(cell)),$(string(umi)),$(nrecords),$(coverage)," *
                    "$(string(LongDNASeq(consensus))),$(string(LongDNASeq(filled_consensus))),$(join(depths, ';'))\n",
                )
            end
        end
    end
end

@inline function is_valid(record::BAM.Record)
    return BAM.isfilled(record) &&
           BAM.ismapped(record) &&
           BAM.isnextmapped(record) &&
           BAM.refid(record) == BAM.nextrefid(record)
end


"""
    make_consensus!(records, reference)

Creates a consensus sequence by frequency, filling the gaps with the reference sequence.
Returns the consensus and a vector with the depth for each reference position.
"""
function make_consensus!(
    records::Vector{BAM.Record},
    reference::LongSequence;
    fill_with_ref::Bool = true,
)::Tuple{Vector{DNA},Vector{Int64}}
    sort!(records, by = x -> BAM.position(x))

    pos_counts_start = DefaultOrderedDict{Int,Int}(0)
    for r in records
        pos_counts_start[BAM.position(r)] += 1
    end

    # in reverse order to pop from first to last
    reverse!(records)

    consensus = fill(DNA_N, length(reference))
    depth = fill(0, length(reference))
    current_reads = MutableLinkedList{SeqRecord}()
    dna_counter = DNACounter()
    pos = nothing

    while (!isempty(records))
        if isnothing(pos) || isempty(current_reads)
            pos = first(pos_counts_start).first
        end

        if haskey(pos_counts_start, pos)
            nreads = pop!(pos_counts_start, pos)
            for i = 1:nreads
                push!(current_reads, bam_record_to_seq(pop!(records)))
            end
        end

        reset_dna_counter!(dna_counter)
        idx = 1
        for read in current_reads
            if pos < read.startpos + read.length
                add_to_dna_counter!(
                    dna_counter,
                    read.sequence[pos-read.startpos+1],
                    (10^(read.quality[pos-read.startpos+1] / 10)),
                )

                idx += 1
                depth[pos] += 1
            else
                delete!(current_reads, idx)
            end
        end

        # replace base in consensus
        consensus[pos] = dna_counter.dna_bases[argmax(dna_counter.counters)]
        pos += 1
    end

    if fill_with_ref
        # fill Ns in consensus with reference nucleotide
        for i = 1:length(consensus)
            if consensus[i] == DNA_N
                consensus[i] = reference[i]
            end
        end
    end

    return consensus, depth
end


"""
    read_reference(references_file, name)

Loads and returns a specific sequence from a FASTA file.
"""
function get_reference(references_filename::String, name::String)::LongSequence
    reader = open(FASTA.Reader, references_filename)
    record = FASTA.Record()
    while !eof(reader)
        read!(reader, record)

        if identifier(record) == name
            return sequence(record)
        end
    end
end


"""
    load_records(records_filename, reads_info)

Loads records in a dictionary by Cell and UMI.
"""
function load_records(
    records_filename::String,
    reads_info::Dict,
    reference_name::String,
)::Dict
    reads_dict = Dict()
    for record in BAM.Reader(open(records_filename))
        # only keep reads where both mates are aligned against the same reference
        if is_valid(record) && BAM.refname(record) == reference_name
            cell = Symbol(reads_info[BAM.tempname(record)][:cell])
            umi = Symbol(reads_info[BAM.tempname(record)][:umi])

            cell_reads = get!(reads_dict, cell, Dict())
            if haskey(cell_reads, umi)
                push!(cell_reads[umi], record)
            else
                cell_reads[umi] = [record]
            end
        end
    end

    return reads_dict
end


function parse_arguments()
    s = ArgParseSettings(description = "Make consensus by Cell/UMI")

    @add_arg_table! s begin
        "--output", "-o"
            help = "output file (.csv)"
            default = "consensus.csv"
        "bamfile"
            help = "BAM file containing reads"
            required = true
        "index"
            help = "List containg Cell and UMI correction for each read"
            required = true
        "reference"
            help = "FASTA file containing reference sequence"
            required = true
        "name"
            help = "Reference to use"
            required = true
    end

    return parse_args(s)
end

main()
