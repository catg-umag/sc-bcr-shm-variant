#!/usr/bin/env julia
using ArgParse, CSV


function main()
    args = parse_arguments()


    passed_consensus = filter(
        x ->
            x.nreads >= args["min_reads"] &&
                x.ref_vdj_coverage >= args["min_vdj_cov"] &&
                x.ref_cdr_coverage >= args["min_cdr_cov"],
        CSV.File(args["input_consensus"]),
    )

    # write fasta
    open(args["output_reads"], "w") do f
        for row in passed_consensus
            write(f, ">$(row.cell)_$(row.umi)\n$(row.consensus)\n")
        end
    end
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
            default = 0.95
            dest_name = "min_vdj_cov"
        "--min-cdr-coverage", "-c"
            help = "minimum CDR Regions coverage"
            arg_type = Float64
            default = 1.0
            dest_name = "min_cdr_cov"
        "--min-reads", "-r"
            help = "minimum reads associated with the consensus"
            arg_type = Int
            default = 2
            dest_name = "min_reads"
        "input_consensus"
            help = "Input CSV with consensus table"
            required = true
        #! format: on
    end

    return parse_args(s)
end


main()
