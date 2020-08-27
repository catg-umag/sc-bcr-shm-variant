#!/usr/bin/env julia
using ArgParse
using CSV
using JSON

include("../lib/julia/barcode_correction.jl")
include("../lib/julia/paired_fastq_writer.jl")


function main()
    args = parse_arguments()

    # load data
    corrections = JSON.parsefile(args["corrections"])
    subjects_from_cell = load_subject_index(args["index"], args["experiment"])
    subjects = collect(Set(values(subjects_from_cell)))

    # prepare writers
    mkpath(args["output_dir"])
    fastq_writers = Dict(
        subject =>
            PairedFastqWriter("$(args["output_dir"])/$(args["experiment"])_$(subject)")
        for subject in vcat(subjects, "NA")
    )

    for (r1, r2) in zip(
        FASTQ.Reader(auto_gzopen(args["reads_r1"])),
        FASTQ.Reader(auto_gzopen(args["reads_r2"])),
    )
        cell = corrections["cells"][string(get_barcode_seq(r1))]
        umi = string(get_umi_seq(r1))

        if !isnothing(cell) && haskey(subjects_from_cell, cell)
            subject = subjects_from_cell[cell]
        else
            subject = "NA"
        end

        if haskey(corrections["umis"], cell) && haskey(corrections["umis"][cell], umi)
            umi = corrections["umis"][cell][umi]
        end

        write_record(fastq_writers[subject], r1, r2, args["experiment"], cell, umi)
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
    s = ArgParseSettings(description = "Splits FASTQ files by subject")

    @add_arg_table! s begin
        "--output-dir", "-o"
            help = "output directory"
            default = "."
            dest_name = "output_dir"
        "--r1"
            help = "R1 reads"
            dest_name = "reads_r1"
            required = true
        "--r2"
            help = "R2 reads"
            dest_name = "reads_r2"
            required = true
        "--index", "-i"
            help = "experiment/cell/cluster index"
            required = true
        "--experiment-name", "-n"
            help = "experiment name"
            required = true
            dest_name = "experiment"
        "--corrections", "-c"
            help = "Cell/UMI barcode corrections"
            required = true
    end

    return parse_args(s)
end


main()
