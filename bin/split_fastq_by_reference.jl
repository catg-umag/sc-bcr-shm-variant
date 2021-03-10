#!/usr/bin/env julia
using ArgParse, CodecZlib, CSV, FASTX

include("../lib/julia/barcode_correction.jl")


function main()
    args = parse_arguments()

    references = Dict(
        x.cell => x.reference for
        x in CSV.File(args["reference_index"]) if x.reference != ""
    )

    mkpath(args["output_dir"])

    readers = (
        r1 = FASTQ.Reader(auto_gzopen(args["input_r1"])),
        r2 = FASTQ.Reader(auto_gzopen(args["input_r2"])),
    )
    writers = Dict()

    record = (r1 = FASTQ.Record(), r2 = FASTQ.Record())
    while !eof(readers.r1) && !eof(readers.r2)
        read!(readers.r1, record.r1)
        read!(readers.r2, record.r2)

        cell = split(FASTQ.identifier(record.r1), ":")[end-1]

        if haskey(references, cell)
            reference = references[cell]
            if !haskey(writers, reference)
                output_base =
                    "$(args["output_dir"])/" *
                    replace(args["input_r1"], r"^.*/|\..*$" => "") *
                    "--$(reference)"
                writers[reference] = (
                    r1 = FASTQ.Writer(open("$(output_base).r1.fastq", "w")),
                    r2 = FASTQ.Writer(open("$(output_base).r2.fastq", "w")),
                )
            end

            write(writers[reference].r1, record.r1)
            write(writers[reference].r2, record.r2)
        end
    end

    for w in values(writers)
        close(w.r1)
        close(w.r2)
    end
end


function parse_arguments()
    s = ArgParseSettings(description = "Splits FASTQ files by best reference")

    @add_arg_table! s begin
        #! format: off
        "input_r1"
            help = "FASTQ file (R1)"
        "input_r2"
            help = "FASTQ file (R2)"
        "reference_index"
            help = "File with each cell and its corresponding reference (.csv)"
        "--output-dir", "-o"
            help = "Output directory"
            default = "."
            dest_name = "output_dir"
        #! format: on
    end

    return parse_args(s)
end


main()