#!/usr/bin/env julia
using Base.Threads
using ArgParse, BioAlignments, BioSequences, DataStructures, CSV, FASTX, XAM


include("../lib/julia/dna_counter.jl")
include("../lib/julia/others.jl")


function main()
    args = parse_arguments()

    # load data
    references = load_references(args["references"])
    ref_position_maps = Dict(k => get_reference_position_maps(v) for (k, v) in references)
    records = load_records(args["bamfile"], Set(keys(references)))
    regions = load_regions(args["reference_regions"])

    open(args["output"], "w") do f
        write(
            f,
            "cell,umi,nreads,ref_vdj_coverage,ref_cdr_coverage,gapped_consensus,consensus,depths\n",
        )

        Threads.@threads for cell in collect(eachindex(records))
            cell_records = records[cell]
            position_counters = [DNACounter() for _ = 1:length(last(first(references)))]
            reference_counters = Dict(k => 0 for k in keys(references))
            for umi in keys(cell_records)
                # reset counters
                for c in position_counters
                    reset!(c)
                end
                for k in keys(reference_counters)
                    reference_counters[k] = 0
                end

                nrecords = length(cell_records[umi])
                consensus, depths, refname = make_consensus!(
                    cell_records[umi],
                    references,
                    ref_position_maps,
                    position_counters,
                    reference_counters,
                )

                # "ungap" consensus and depth
                ungapped_consensus = clean_consensus(consensus)
                depths = [x for (i, x) in enumerate(depths) if consensus[i] != DNA_Gap]

                coverages = get_coverages(depths, regions[refname])

                # write result
                row = [
                    string(cell),
                    string(umi),
                    nrecords,
                    coverages.vdj,
                    coverages.cdr,
                    string(LongDNASeq(consensus)),
                    ungapped_consensus,
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
    references_positions_maps::Dict{String,Vector{Int64}},
    position_counters::Vector{DNACounter{Float64}},
    reference_counters::Dict{String,Int64},
)::Tuple{Vector{DNA},Vector{Int64},String}
    reference_length = length(position_counters)

    consensus = fill(DNA_N, reference_length)
    depth = fill(0, reference_length)

    for record in records
        sequence = BAM.sequence(record)
        quality = BAM.quality(record)
        refname = BAM.refname(record)
        ref_positions = references_positions_maps[refname]

        pos = 1
        ref_pos = BAM.position(record)
        for (op, n) in zip(BAM.cigar_rle(record)...)
            if op in (OP_SOFT_CLIP, OP_INSERT)
                pos += n
            elseif op == OP_DELETE
                ref_pos += n
            elseif op == OP_MATCH
                for i = 0:(n-1)
                    corr_pos = ref_positions[ref_pos+i]
                    add!(
                        position_counters[corr_pos],
                        sequence[pos+i],
                        (10^(quality[pos+i] / 10)),
                    )
                    depth[corr_pos] += 1
                end
                pos += n
                ref_pos += n
            end

        end
        reference_counters[refname] += 1
    end

    for pos in eachindex(consensus)
        counter = position_counters[pos]
        if counter.sum > 0
            consensus[pos] = counter.dna_bases[argmax(counter.counters)]
        end
    end

    # add reference gaps
    most_used_refname = last(findmax(reference_counters))
    most_used_ref = references[most_used_refname]
    for i in eachindex(most_used_ref)
        if most_used_ref[i] == DNA_Gap
            consensus[i] = DNA_Gap
        end
    end

    return consensus, depth, most_used_refname
end


"""
    get_coverages(depths, regions)

Gets coverages from V(D)J and CDR regions
"""
function get_coverages(depths::Vector{Int64}, regions::NamedTuple)::NamedTuple
    cdrs = ["cdr1", "cdr2", "cdr3"]
    cdrs = filter(
        x -> all(Symbol("$(x)_$(y)") in keys(regions) for y in ("start", "end")),
        cdrs,
    )
    vdj_coverage = cdr_coverage = 0.0
    for (i, d) in enumerate(depths)
        if d > 0
            if all(Symbol("vdj_$(y)") in keys(regions) for y in ("start", "end")) &&
               (regions.vdj_start <= i <= regions.vdj_end)
                vdj_coverage += 1
            end
            if any(
                regions[Symbol("$(cdr)_start")] <= i <= regions[Symbol("$(cdr)_end")] for
                cdr in cdrs
            )
                cdr_coverage += 1
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
    get_reference_position_maps(sequence)

Get "new" positions from old positions using gaps
"""
function get_reference_position_maps(sequence::LongSequence)::Vector{Int64}
    maps = Int64[]
    for (i, x) in enumerate(sequence)
        if x != DNA_Gap
            push!(maps, i)
        end
    end

    return maps
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
