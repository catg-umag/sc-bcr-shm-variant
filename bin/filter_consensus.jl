#!/usr/bin/env julia
using ArgParse, CSV


function main()
    args = parse_arguments()

    region_positions = load_regions_positions(args["region_positions"])
    positions = region_positions[args["reference_name"]]

    passed_consensus = get_filtered_consensus(
        args["input_consensus"],
        positions,
        args["min_vdj_cov"],
        args["min_cdr_cov"],
        args["min_reads"],
    )

    # write fasta
    open(args["output_reads"], "w") do f
        for row in passed_consensus
            write(f, ">$(row.cell)_$(row.umi)\n$(row.filled_consensus)\n")
        end
    end
end


"""
    get_filtered_consensus(consensus_file, positions, min_vdj_cov, min_cdr_cov, min_reads)

Reads consensus from file and checks them against filters,
returning only the ones that passed all filters
"""
function get_filtered_consensus(
    consensus_file::String,
    positions::Dict{String,Int64},
    min_vdj_cov::AbstractFloat,
    min_cdr_cov::AbstractFloat,
    min_reads::Integer,
)::Vector{NamedTuple}
    data = []
    for x in CSV.File(consensus_file)
        # if not enough reads, skip
        (x.nrecords < min_reads) && continue

        # transform depths to numeric
        depths = [parse(Int, x) for x in split(x.depths, ";")]

        # get coverages
        vdj_cov =
            get_coverage(
                depths,
                positions["VDJ-REGION_start"],
                positions["VDJ-REGION_end"],
            ) / (positions["VDJ-REGION_end"] - positions["VDJ-REGION_start"] + 1)

        # no enough VDJ coverage? skip
        (vdj_cov < min_vdj_cov) && continue

        cdr_names = unique(map(
            x -> x.match,
            filter(
                x -> !isnothing(x),
                map(x -> match(r"^CDR.", x), collect(keys(positions))),
            ),
        ))
        cdr_cov =
            sum(
                get_coverage(
                    depths,
                    positions["$(cdr)-IMGT_start"],
                    positions["$(cdr)-IMGT_end"],
                ) for cdr in cdr_names
            ) / sum(
                (positions["$(cdr)-IMGT_end"] - positions["$(cdr)-IMGT_start"] + 1)
                for cdr in cdr_names
            )

        if cdr_cov >= min_cdr_cov
            push!(data, (cell = x.cell, umi = x.umi, filled_consensus = x.filled_consensus))
        end
    end

    return data
end


function load_regions_positions(region_positions_file)::Dict{String,Dict{String,Int64}}
    region_positions = Dict(
        x.Sequence_ID =>
            Dict(String(y) => x[y] for y in eachindex(x) if y != :Sequence_ID)
        for x in CSV.File(region_positions_file)
    )

    return region_positions
end


function get_coverage(depths::Vector{Int}, pos_start::Int, pos_end::Int)::Int
    return sum(x > 0 for x in depths[pos_start:pos_end])
end


function parse_arguments()
    s = ArgParseSettings(
        description = "Filters Cell/UMI consensus (from .csv) by a set of metrics",
    )

    @add_arg_table! s begin
        #! format: off
        "--output", "-o"
            help = "output FASTA file"
            default = "consensus_filtered.fasta"
            dest_name = "output_reads"
        "--min-vdj-coverage", "-v"
            help = "minimum VDJ Region coverage"
            arg_type = Float64
            default = 0.9
            dest_name = "min_vdj_cov"
        "--min-cdr-coverage", "-c"
            help = "minimum CDR Regions coverage"
            arg_type = Float64
            default = 0.95
            dest_name = "min_cdr_cov"
        "--min-reads", "-r"
            help = "minimum reads associated with the consensus"
            arg_type = Int
            default = 2
            dest_name = "min_reads"
        "--reference-name", "-n"
            help = "Reference for the current file, if not defined input basename will be used"
            dest_name = "reference_name"
        "input_consensus"
            help = "Input CSV with consensus table"
            required = true
        "region_positions"
            help = "Positions of each of IG Regions (VDJ and CDRs) in the references"
            required = true
        #! format: on
    end

    return parse_args(s)
end


main()
