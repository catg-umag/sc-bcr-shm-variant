#!/usr/bin/env julia
using ArgMacros, FASTX, JSON3, StructTypes

include("../lib/julia/barcode_correction.jl")

StructTypes.StructType(::Type{LongDNASeq}) = StructTypes.StringType()


function main()
    args = parse_arguments()

    whitelist_barcodes = load_whitelist_barcodes(args.whitelist)
    barcode_dist = get_barcode_dist(args.input_r1, whitelist_barcodes)

    barcodes_corrected, umi_counts_by_bc =
        correct_barcodes_and_count_umis(args.input_r1, whitelist_barcodes, barcode_dist)

    # get umi corrections
    umi_corrections = get_umi_corrections(umi_counts_by_bc)

    # write corrections (if defined)
    if !isnothing(args.output_corrections)
        open(args.output_corrections, "w") do f
            JSON3.pretty(f, Dict("cells" => barcodes_corrected, "umis" => umi_corrections))
        end
    end

    # write fastq files
    writers = (
        r1 = FASTQ.Writer(auto_gzopen(args.output_r1; mode = "w")),
        r2 = FASTQ.Writer(auto_gzopen(args.output_r2; mode = "w")),
    )
    readers = (
        r1 = FASTQ.Reader(auto_gzopen(args.input_r1)),
        r2 = FASTQ.Reader(auto_gzopen(args.input_r2)),
    )
    record = (r1 = FASTQ.Record(), r2 = FASTQ.Record())
    offset = args.cellbc_size + args.umibc_size

    while !eof(readers.r1) && !eof(readers.r2)
        read!(readers.r1, record.r1)
        read!(readers.r2, record.r2)

        cellbc = get_barcode_seq(record.r1)
        umi = get_umi_seq(record.r1)
        if haskey(barcodes_corrected, cellbc) && isnothing(barcodes_corrected[cellbc])
            cellbc = barcodes_corrected[cellbc]
        end
        if haskey(umi_corrections, cellbc) && haskey(umi_corrections[cellbc], umi)
            umi = umi_corrections[cellbc][umi]
        end

        bc_str = ":$(cellbc):$(umi)"
        write(
            writers.r1,
            FASTQ.Record(
                identifier(record.r1) * bc_str,
                description(record.r1),
                FASTQ.sequence(record.r1)[offset+1:end],
                quality(record.r1)[offset+1:end],
            ),
        )
        write(
            writers.r2,
            FASTQ.Record(
                identifier(record.r2) * bc_str,
                description(record.r2),
                FASTQ.sequence(record.r2),
                quality(record.r2),
            ),
        )
    end

    close(writers.r1)
    close(writers.r2)
end


"""
    get_umi_corrections(umi_counts)
"""
function get_umi_corrections(umi_counts)
    umi_corrections = Dict()
    for (bc, umi_count) in umi_counts
        if length(umi_count) > 1
            umi_corrector = correct_umis(umi_count)
            if !isempty(umi_corrector)
                umi_corrections[bc] = umi_corrector
            end
        end
    end

    return umi_corrections
end


function parse_arguments()
    args = @tuplearguments begin
        @helpusage """
            extract_barcodes.jl -i FASTQ_R1 -I FASTQ_R2 -o FASTQ_R1 -O FASTQ_R2
                -w BARCODE_WHITELIST [-c CORRECTIONS] [-b CELL_BARCODE_SIZE] [-u UMI_BARCODE_SIZE]
        """
        @helpdescription "Adds Cell and UMI barcodes to the read name and deleting them from the sequence"

        @argumentrequired String input_r1 "-i" "--input-r1"
        @arghelp "FASTQ input file (R1, must have the barcodes)"
        @argtest input_r1 isfile "The input_r1 argument must be a valid file"

        @argumentrequired String input_r2 "-I" "--input-r2"
        @arghelp "FASTQ input file (R2)"
        @argtest input_r2 isfile "The input_r2 argument must be a valid file"

        @argumentrequired String output_r1 "-o" "--output-r1"
        @arghelp "FASTQ output file (R1)"

        @argumentrequired String output_r2 "-O" "--output-r2"
        @arghelp "FASTQ output file (R2)"

        @argumentrequired String whitelist "-w" "--barcode-whitelist"
        @arghelp "cell barcode whitelist (one per line)"
        @argtest whitelist isfile "The whitelist argument must be a valid file"

        @argumentoptional String output_corrections "-c" "--output-corrections"
        @arghelp "output corrections in a JSON file (must file the filename desired)"

        @argumentdefault Int 16 cellbc_size "-b" "--cell-barcode-size"
        @arghelp "cell barcode size"

        @argumentdefault Int 10 umibc_size "-u" "--umi-barcode-size"
        @arghelp "UMI barcode size"
    end

    return args
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
