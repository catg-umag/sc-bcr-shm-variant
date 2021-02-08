#!/usr/bin/env julia
using ArgParse, CodecZlib, CSV, FASTX, JSON

include("../lib/julia/barcode_correction.jl")


function main()
    args = parse_arguments()

    offset = 26 # 16 Cell + 10 UMI
    corrections =
        haskey(args, "corrections") ? JSON.parsefile(args["corrections"]) : nothing

    writers = (
        r1 = FASTQ.Writer(GzipCompressorStream(open(args["output_r1"], "w"))),
        r2 = FASTQ.Writer(GzipCompressorStream(open(args["output_r2"], "w"))),
    )
    readers = (
        r1 = FASTQ.Reader(auto_gzopen(args["input_r1"])),
        r2 = FASTQ.Reader(auto_gzopen(args["input_r2"])),
    )
    record = (r1 = FASTQ.Record(), r2 = FASTQ.Record())
    while !eof(readers.r1) && !eof(readers.r2)
        read!(readers.r1, record.r1)
        read!(readers.r2, record.r2)

        cellbc, umi = get_barcodes_seqs_as_str(record.r1)
        if !isnothing(corrections)
            cellbc = corrections["cells"][cellbc]
            if haskey(corrections["umis"], cellbc) &&
               haskey(corrections["umis"][cellbc], umi)
                umi = corrections["umis"][cellbc][umi]
            end
        end

        bc_str = ":$(cellbc):$(umi)"
        write(
            writers.r1,
            FASTQ.Record(
                identifier(record.r1) * bc_str,
                description(record.r1),
                sequence(record.r1)[offset+1:end],
                quality(record.r1)[offset+1:end],
            ),
        )
        write(
            writers.r2,
            FASTQ.Record(
                identifier(record.r2) * bc_str,
                description(record.r2),
                sequence(record.r2),
                quality(record.r2),
            ),
        )
    end

    close(writers.r1)
    close(writers.r2)
end

function parse_arguments()
    s = ArgParseSettings(
        description = "Extracts RAW Cell and UMI barcodes and puts them into the FASTQ header",
    )

    @add_arg_table! s begin
        #! format: off
        "--output-reads-r1", "-o"
            help = "output FASTQ file (reads without barcodes)"
            required = true
            dest_name = "output_r1"
        "--output-reads-r2", "-O"
            help = "output FASTQ file (reads without barcodes)"
            required = true
            dest_name = "output_r2"
        "--input-reads-r1", "-i"
            help = "Input FASTQ file (R1)"
            required = true
            dest_name = "input_r1"
        "--input-reads-r2", "-I"
            help = "Input FASTQ file (R2)"
            required = true
            dest_name = "input_r2"
        "--barcode-corrections", "-c"
            help = "JSON file containing corrections for Cell and UMI barcodes"
            dest_name = "corrections"
        #! format: on
    end

    return parse_args(s)
end


main()
