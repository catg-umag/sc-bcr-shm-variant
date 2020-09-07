#!/usr/bin/env julia
using ArgParse
using CodecZlib
using CSV
using FASTX
using JSON

include("../lib/julia/barcode_correction.jl")


function main()
    args = parse_arguments()

    offset = 26 # 16 Cell + 10 UMI
    corrections = haskey(args, "corrections") ? JSON.parsefile(args["corrections"]) : nothing
    reads_writer = FASTQ.Writer(GzipCompressorStream(open(args["output_reads"], "w")))
    
    barcodes_writer = open(args["output_barcodes"], "w")
    write(barcodes_writer, "id,cell,umi\n")

    reader = FASTQ.Reader(auto_gzopen(args["input_fastq"]))
    record = FASTQ.Record()
    while !eof(reader)
        read!(reader, record)

        cell = string(get_barcode_seq(record))
        umi = string(get_umi_seq(record))
        if !isnothing(corrections)
            cell = corrections["cells"][cell]
            if haskey(corrections["umis"], cell) && haskey(corrections["umis"][cell], umi)
                umi = corrections["umis"][cell][umi]
            end
        end

        write(
            barcodes_writer, "$(identifier(record)),$(cell),$(umi)\n"
        )

        write(
            reads_writer,
            FASTQ.Record(
                identifier(record),
                description(record),
                sequence(record)[offset+1:end],
                quality(record)[offset+1:end],
            ),
        )


    end

    close(reads_writer)
    close(barcodes_writer)
end

function parse_arguments()
    s = ArgParseSettings(
        description = "Extracts RAW Cell and UMI barcodes in a differente FASTQ file",
    )

    @add_arg_table! s begin
        "--output-reads", "-r"
            help = "output FASTQ file (reads without barcodes)"
            default = "reads_cleaned.fastq.gz"
            dest_name = "output_reads"
        "--output-barcodes", "-b"
            help = "output CSV file with barcodes"
            default = "barcodes.csv"
            dest_name = "output_barcodes"
        "--barcode-corrections", "-c"
            help = "JSON file containing corrections for Cell and UMI barcodes"
            dest_name = "corrections"
        "input_fastq"
            help = "Input FASTQ file to extract barcodes (muste be R1)"
            required = true
    end

    return parse_args(s)
end


main()
