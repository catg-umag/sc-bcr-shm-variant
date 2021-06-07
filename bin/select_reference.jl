#!/usr/bin/env julia
using ArgMacros, XAM

RefCounts = Dict{String,Dict{String,Integer}}


function main()
    args = parse_arguments()

    counts = load_counts(args.bamfile)

    best_references = get_best_references(counts)

    open(args.output, "w") do f
        write(f, "cell,reference\n")
        for row in best_references
            if row[2] != ""
                write(f, "$(row[1]),$(row[2])\n")
            end
        end
    end
end


"""
    load_counts(filename)

Load a sum of MAPQ values by each cell/reference
"""
function load_counts(filename::String)::RefCounts
    cell_ref_scores = Dict{String,Dict}()
    reader = open(BAM.Reader, filename)
    record = BAM.Record()
    while !eof(reader)
        empty!(record)
        read!(reader, record)

        if BAM.isprimary(record) &&
           BAM.isfilled(record) &&
           BAM.ismapped(record) &&
           BAM.mappingquality(record) != 0
            cell = split(BAM.tempname(record), ":")[end-1]
            reference = BAM.refname(record)

            get!(cell_ref_scores, cell, Dict{String,Integer}())
            cell_ref_scores[cell][reference] =
                get(cell_ref_scores[cell], reference, 0) + BAM.mappingquality(record)
        end
    end

    return cell_ref_scores
end


"""
    get_best_references(counts)

For each cell, select the references (heavy and light) with the max count
"""
function get_best_references(counts::RefCounts)
    cell_best_refs = []
    for (cell, counts) in counts
        bests = Dict(:name => "", :count => 0)
        for (refname, count) in counts

            if bests[:count] < count
                bests[:name] = refname
                bests[:count] = count
            end
        end
        push!(cell_best_refs, [cell, bests[:name]])
    end
    
    sort!(cell_best_refs, by = x -> x[1])

    return cell_best_refs
end


function parse_arguments()
    args = @tuplearguments begin
        @helpusage "select_reference.jl -i BAMFILE -o OUTBASE"
        @helpdescription "Selects the best reference for each cell"

        @argumentrequired String bamfile "-i" "--bamfile"
        @arghelp "input file (.bam)"
        @argtest bamfile isfile "The input file must be a valid file"

        @argumentrequired String output "-o" "--output"
        @arghelp "Output file (.csv)"
    end

    return args
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
