#!/usr/bin/env julia
using ArgParse
using CSV
using JSON
using FASTX

include("../lib/julia/barcode_correction.jl")

JSON.lower(s::LongSequence) = string(s)


function main()
    args = parse_arguments()

    whitelist_barcodes = load_whitelist_barcodes(args["whitelist"])
    barcode_dist = get_barcode_dist(args["fastq_file"], whitelist_barcodes)

    barcodes_corrected, umi_counts_by_bc = correct_barcodes_and_count_umis(
        args["fastq_file"],
        whitelist_barcodes,
        barcode_dist,
    )

    # correct UMIs for each barcode
    umi_corrections = Dict()
    for (bc, umi_counts) in umi_counts_by_bc
        if length(umi_counts) > 1
            umi_corrector = correct_umis(umi_counts)
            if !isempty(umi_corrector)
                umi_corrections[bc] = umi_corrector
            end
        end
    end

    # write outputs
    open(args["corrections_outfile"], "w") do f
        JSON.print(f, Dict("cells" => barcodes_corrected, "umis" => umi_corrections), 4)
    end
end


function parse_arguments()
    s = ArgParseSettings(
        description = "Perform cell barcode and UMI correction and saves outputs in JSON files",
    )

    @add_arg_table! s begin
        "--corrections-outfile", "-o"
            help = "filename for corrections output (should include .json)"
            default = "barcode_corrections.json"
            dest_name = "corrections_outfile"
        "fastq_file"
            help = "FASTQ file to process (must be R1)"
            required = true
        "whitelist"
            help = "cell barcode whitelist"
            required = true
    end

    return parse_args(s)
end


main()
