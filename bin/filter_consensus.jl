#!/usr/bin/env julia
using ArgMacros, CSV


function main()
    args = parse_arguments()


    passed_consensus = filter(
        x ->
            isa(x.ref_vdj_coverage, Float64) &&
                isa(x.ref_cdr_coverage, Float64) &&
                x.nreads >= args.min_reads &&
                x.ref_vdj_coverage >= args.min_vdj_cov &&
                x.ref_cdr_coverage >= args.min_cdr_cov,
        CSV.File(args.input),
    )

    # write fasta
    open(args.output, "w") do f
        for row in passed_consensus
            write(f, ">$(row.cell)_$(row.umi)\n$(row.consensus)\n")
        end
    end
end


function parse_arguments()
    ratio_validator = x -> 0 <= x <= 1.0
    positive_validator = x -> x > 0

    args = @tuplearguments begin
        @helpusage """
            filter_consensus.jl -i INPUT_CONSENSUS [o OUTPUT]
                [-v MIN_VDJ_COVERAGE] [-c MIN_CDR_COVERAGE] [-r MIN_READS]"""
        @helpdescription "Filters Cell/UMI consensus (from .csv) by a set of metrics"

        @argumentrequired String input "-i" "--input-consensus"
        @arghelp "input file with consensus list (.csv)"
        @argtest input isfile "The input argument must be a valid file"

        @argumentdefault String "consensus_filtered.fasta" output "-o" "--output"
        @arghelp "output FASTA file with consensuses passing the filters"

        @argumentdefault Float64 0.95 min_vdj_cov "-v" "--min-vdj-coverage"
        @arghelp "minimum VDJ region coverage"
        @argtest min_vdj_cov ratio_validator "The min-vdj-coverage must be a valid ratio"

        @argumentdefault Float64 1.0 min_cdr_cov "-c" "--min-cdr-coverage"
        @arghelp "minimum CDR regions coverage"
        @argtest min_cdr_cov ratio_validator "The min-cdr-coverage must be a valid ratio"

        @argumentdefault Int 2 min_reads "-c" "--min-reads"
        @arghelp "minimum reads associated with the consensu"
        @argtest min_cdr_cov positive_validator "The min-reads must be a positive value"

    end

    return args
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
