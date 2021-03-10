#!/usr/bin/env julia
using ArgParse, XAM

RefCounts = Dict{String,Dict{String,Integer}}


function main()
    args = parse_arguments()

    counts = load_counts(args["bamfile"])

    best_references = get_best_references(counts)

    open(args["output_heavy"], "w") do f
        write(f, "cell,reference\n")
        for row in best_references
            if row[2] != ""
                write(f, "$(row[1]),$(row[2])\n")
            end
        end
    end

    open(args["output_light"], "w") do f
        write(f, "cell,reference\n")
        for row in best_references
            if row[3] != ""
                write(f, "$(row[1]),$(row[3])\n")
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
        bests = Dict(
            :heavy => Dict(:name => "", :count => 0),
            :light => Dict(:name => "", :count => 0),
        )
        for (refname, count) in counts
            chain = split(refname, "_")[3]
            chain_type = chain == "H" ? :heavy : :light

            if bests[chain_type][:count] < count
                bests[chain_type][:name] = refname
                bests[chain_type][:count] = count
            end
        end
        push!(cell_best_refs, [cell, bests[:heavy][:name], bests[:light][:name]])
    end
    
    sort!(cell_best_refs, by = x -> x[1])

    return cell_best_refs
end


function parse_arguments()
    s = ArgParseSettings(description = "Selects the best reference for each cell")

    @add_arg_table! s begin
        #! format: off
        "--bamfile", "-i"
            help = "input file (.bam)"
            required = true
        "--output-heavy", "-H"
            help = "output file for heavy chain (.csv)"
            required = true
            dest_name = "output_heavy"
        "--output-light", "-L"
            help = "output file for light chain (.csv)"
            required = true
            dest_name = "output_light"
        #! format: on
    end

    return parse_args(s)
end


main()