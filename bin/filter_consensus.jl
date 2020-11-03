#!/usr/bin/env julia
using ArgParse, CSV, DataFrames, Query


function main()
    args = parse_arguments()

    region_positions = load_regions_positions(args["region_positions"])
    positions = region_positions[args["reference_name"]]

    df_consensus = load_consensus(args["input_consensus"], positions)

    # filter consensus
    df_pass =
        df_consensus |>
        @filter(
            (_.vdj_cov >= args["min_vdj_cov"]) &&
            (_.cdr_cov >= args["min_cdr_cov"]) &&
            (_.nrecords >= args["min_reads"])
        ) |>
        DataFrame

    # write fasta
    open(args["output_reads"], "w") do f
        for row in eachrow(df_pass)
            write(f, ">$(row.cell)_$(row.umi)\n$(row.filled_consensus)\n")
        end
    end
end


function load_consensus(consensus_file, positions)
    data = []
    for x in CSV.File(consensus_file)
        row_dict = Dict(y => x[y] for y in eachindex(x))
        row_dict[:depths] = [parse(Int, x) for x in split(row_dict[:depths], ";")]

        # get coverages
        row_dict[:vdj_cov] =
            get_coverage(
                row_dict[:depths],
                positions["VDJ-REGION_start"],
                positions["VDJ-REGION_end"],
            ) / (positions["VDJ-REGION_end"] - positions["VDJ-REGION_start"] + 1)

        cdr_names = unique(map(
            x -> x.match,
            filter(
                x -> !isnothing(x),
                map(x -> match(r"^CDR.", x), collect(keys(positions))),
            ),
        ))
        row_dict[:cdr_cov] =
            sum(
                get_coverage(
                    row_dict[:depths],
                    positions["$(cdr)-IMGT_start"],
                    positions["$(cdr)-IMGT_end"],
                ) for cdr in cdr_names
            ) / sum(
                (positions["$(cdr)-IMGT_end"] - positions["$(cdr)-IMGT_start"] + 1)
                for cdr in cdr_names
            )

        push!(data, row_dict)
    end

    df = DataFrame([Dict(y => x[y] for y in keys(x) if y != :depths) for x in data])
end


function load_regions_positions(region_positions_file)
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
            arg_type = Number
            default = 0.9
            dest_name = "min_vdj_cov"
        "--min-cdr-coverage", "-c"
            help = "minimum CDR Regions coverage"
            arg_type = Number
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