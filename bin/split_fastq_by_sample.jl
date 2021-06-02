#!/usr/bin/env julia
using ArgMacros, CSV, FASTX, JSON

include("../lib/julia/barcode_correction.jl")
include("../lib/julia/paired_fastq_writer.jl")


function main()
    args = parse_arguments()

    # load data
    subjects_from_cell = load_subject_index(args.index, args.experiment)
    subjects = collect(Set(values(subjects_from_cell)))

    # prepare writers
    mkpath(args.output_dir)
    fastq_writers = Dict(
        subject =>
            PairedFastqWriter("$(args.output_dir)/$(args.experiment)_$(subject)")
        for subject in vcat(subjects, "NA")
    )

    for (r1, r2) in zip(
        FASTQ.Reader(auto_gzopen(args.reads_r1)),
        FASTQ.Reader(auto_gzopen(args.reads_r2)),
    )
        cellbc, umi = split(FASTX.identifier(r1), ":")[end-1: end]

        if !isnothing(cellbc) && haskey(subjects_from_cell, cellbc)
            subject = subjects_from_cell[cellbc]
        else
            subject = "NA"
        end

        write_record(fastq_writers[subject], r1, r2, args.experiment, string(cellbc), string(umi))
    end

    # close writers
    for subject in keys(fastq_writers)
        close_writer(fastq_writers[subject])
    end
end


function load_subject_index(index_filename::String, experiment::String)::Dict{String,String}
    index_dict = Dict()
    for row in CSV.File(index_filename)
        if row.source == experiment && row.cluster != "NA"
            cell = match(r"(?<=:)(.*)(?=-)", row.cell).match
            index_dict[cell] = row.cluster
        end
    end

    return index_dict
end


function parse_arguments()
    args = @tuplearguments begin
        @helpusage """
            split_fastq_by_sample.jl --r1 FASTQ_R1 --r2 FASTQ_R2
                -i INDEX -n EXP_NAME [-o OUTPUT_DIR]"""
        @helpdescription "Splits FASTQ files by subject"

        @argumentrequired String reads_r1 "--r1"
        @arghelp "R1 reads"
        @argtest reads_r1 isfile "The R1 argument must be a valid file"

        @argumentrequired String reads_r2 "--r2"
        @arghelp "R2 reads"
        @argtest reads_r2 isfile "The R2 argument must be a valid file"

        @argumentrequired String index "-i" "--index"
        @arghelp "experiment/cell/cluster index"
        @argtest index isfile "The index argument must be a valid file"

        @argumentrequired String experiment "-n" "--experiment-name"
        @arghelp "experiment name"

        @argumentdefault String "." output_dir "-o" "--output-dir"
        @arghelp "output directory"
    end

    return args
end


main()
