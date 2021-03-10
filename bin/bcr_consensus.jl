#!/usr/bin/env julia
using Base.Threads
using ArgParse, BioAlignments, BioSequences, DataStructures, CSV, FASTX, XAM


include("../lib/julia/dna_counter.jl")
include("../lib/julia/others.jl")


function main()
    args = parse_arguments()

    # load data
    references = load_references(args["references"])
    records = load_records(args["bamfile"], Set(keys(references)))
    regions = load_regions(args["reference_regions"])

    open(args["output"], "w") do f
        write(
            f,
            "cell,umi,nreads,ref_vdj_coverage,ref_cdr_coverage,consensus,aligned_consensus,depths\n",
        )

        Threads.@threads for cell in collect(eachindex(records))
            cell_records = records[cell]
            for umi in keys(cell_records)
                nrecords = length(cell_records[umi])
                consensus, depths, refname = make_consensus!(cell_records[umi], references)
                cleaned_consensus = clean_consensus(consensus)
                coverages = get_coverages(depths, regions[refname])

                # write result
                row = [
                    string(cell),
                    string(umi),
                    nrecords,
                    coverages.vdj,
                    coverages.cdr,
                    cleaned_consensus,
                    string(LongDNASeq(consensus)),
                    join(depths, ";"),
                ]
                write(f, join(row, ",") * "\n")
            end
        end
    end
end


"""
    make_consensus!(records, references, references_positions_maps, position_counters, reference_counters)

Creates a consensus sequence by frequency, filling the gaps with the reference sequence.
Returns the consensus and a vector with the depth for each reference position.
"""
function make_consensus!(
    records::Vector{BAM.Record},
    references::Dict{String,LongSequence},
)::Tuple{Vector{DNA},Vector{Int64},String}
    refname = BAM.refname(records[1])
    reference = references[refname]
    reference_length = length(reference)
    position_counters = [DNACounter() for _ = 1:reference_length]

    consensus = fill(DNA_N, reference_length)
    depth = fill(0, reference_length)

    for record in records
        sequence = BAM.sequence(record)
        quality = BAM.quality(record)

        pos = 1
        ref_pos = BAM.position(record)
        for (op, n) in zip(BAM.cigar_rle(record)...)
            if op in (OP_SOFT_CLIP, OP_INSERT)
                pos += n
            elseif op == OP_DELETE
                ref_pos += n
            elseif op == OP_MATCH
                for i = 0:(n-1)
                    add!(
                        position_counters[ref_pos+i],
                        sequence[pos+i],
                        (10^(quality[pos+i] / 10)),
                    )
                    depth[ref_pos+i] += 1
                end
                pos += n
                ref_pos += n
            end

        end
    end

    for pos in eachindex(consensus)
        counter = position_counters[pos]
        if counter.sum > 0
            consensus[pos] = counter.dna_bases[argmax(counter.counters)]
        end
    end

    return consensus, depth, refname
end


"""
    get_coverages(depths, regions)

Gets coverages from V(D)J and CDR regions
"""
function get_coverages(depths::Vector{Int64}, regions::NamedTuple)::NamedTuple
    cdrs = ["cdr1", "cdr2", "cdr3"]
    cdrs = filter(
        x -> all(not_missing_key(regions, "$(x)_$(y)") for y in ("start", "end")),
        cdrs,
    )
    vdj_coverage = cdr_coverage = 0.0
    for (i, d) in enumerate(depths)
        if d > 0
            if all(not_missing_key(regions, "vdj_$(y)") for y in ("start", "end")) &&
               (regions.vdj_start <= i <= regions.vdj_end)
                vdj_coverage += 1

                if any(
                    regions[Symbol("$(cdr)_start")] <= i <= regions[Symbol("$(cdr)_end")]
                    for cdr in cdrs
                )
                    cdr_coverage += 1
                end
            end
        end
    end

    vdj_coverage = vdj_coverage / (regions.vdj_end - regions.vdj_start + 1)
    cdr_coverage =
        cdr_coverage / sum(
            regions[Symbol("$(cdr)_end")] - regions[Symbol("$(cdr)_start")] + 1 for
            cdr in cdrs
        )

    return (vdj = vdj_coverage, cdr = cdr_coverage)
end


"""
    clean_consensus(consensus)

Clean consensus from gaps and undefined bases (N) on both ends
"""
function clean_consensus(consensus::Vector{DNA})::String
    # first, delete gaps
    clean_consensus = string(LongDNASeq([x for x in consensus if x != DNA_Gap]))
    # next, delete N's at both ends
    clean_consensus = replace(clean_consensus, r"^N*|N*$" => "")

    return clean_consensus
end


"""
    load_references(references_file)

Loads and returns a specific sequence from a FASTA file.
"""
function load_references(references_filename::String)::Dict{String,LongSequence}
    references = Dict()
    for record in open(FASTA.Reader, references_filename)
        references[identifier(record)] = FASTX.sequence(record)
    end

    return references
end


"""
    load_records(records_filename, selected_references)

Loads records in a dictionary by Cell and UMI.
"""
function load_records(records_filename::String, selected_references::Set{String})::Dict
    reads_dict = Dict()
    for record in BAM.Reader(open(records_filename))
        if is_valid(record) && BAM.refname(record) in selected_references
            id_splitted = split(BAM.tempname(record), ":")
            cellbc = string(id_splitted[end-1])
            umi = string(id_splitted[end])

            cell_reads = get!(reads_dict, cellbc, Dict())
            if haskey(cell_reads, umi)
                push!(cell_reads[umi], record)
            else
                cell_reads[umi] = [record]
            end
        end
    end

    return reads_dict
end


function load_regions(regions_filename::String)::Dict
    regions = Dict(row.sequence_id => NamedTuple(row) for row in CSV.File(regions_filename))
    return regions
end


@inline function is_valid(record::BAM.Record)
    return BAM.isprimary(record) && BAM.isfilled(record) && BAM.ismapped(record)
end


@inline function not_missing_key(tuple::NamedTuple, key::String)::Bool
    key_sym = Symbol(key)
    return haskey(tuple, key_sym) && !ismissing(tuple[key_sym])
end


function parse_arguments()
    s = ArgParseSettings(description = "Make consensus by Cell/UMI")

    @add_arg_table! s begin
        #! format: off
        "--output", "-o"
            help = "output file (.csv)"
            default = "consensus.csv"
        "bamfile"
            help = "BAM file containing reads"
            required = true
        "references"
            help = "FASTA file containing aligned reference sequences"
            required = true
        "reference_regions"
            help = "CSV file containing VDJ and CDR region positions"
            required = true
        #! format: on
    end

    return parse_args(s)
end

main()
