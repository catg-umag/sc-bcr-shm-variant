#!/usr/bin/env julia
using ArgMacros, CodecZlib, CSV, FASTX

include("../lib/julia/barcode_correction.jl")


function main()
    args = parse_arguments()

    references = Dict(
        x.cell => x.reference for
        x in CSV.File(args.index) if x.reference != ""
    )

    mkpath(args.output_dir)

    readers = (
        r1 = FASTQ.Reader(auto_gzopen(args.reads_r1)),
        r2 = FASTQ.Reader(auto_gzopen(args.reads_r2)),
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
                    "$(args.output_dir)/" *
                    replace(args.reads_r1, r"^.*/|\..*$" => "") *
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
    args = @tuplearguments begin
        @helpusage """
        split_fastq_by_reference.jl --r1 FASTQ_R1 --r2 FASTQ_R2
            -i INDEX [-o OUTPUT_DIR]"""
        @helpdescription "Splits FASTQ files by best reference"

        @argumentrequired String reads_r1 "--r1"
        @arghelp "R1 reads"
        @argtest reads_r1 isfile "The R1 argument must be a valid file"

        @argumentrequired String reads_r2 "--r2"
        @arghelp "R2 reads"
        @argtest reads_r2 isfile "The R2 argument must be a valid file"

        @argumentrequired String index "-i" "--reference-index"
        @arghelp "file with each cell and its corresponding reference (.csv)"
        @argtest index isfile "The index argument must be a valid file"

        @argumentdefault String "." output_dir "-o" "--output-dir"
        @arghelp "output directory"
    end

    return args
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
