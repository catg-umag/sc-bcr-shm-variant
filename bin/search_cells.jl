#!/usr/bin/env julia
using Base.Threads
using ArgParse, CSV, DataFrames, Edlib, FASTX


function main()
    args = parse_arguments()

    external_sequences = load_sequences_str(args["external_sequences"])
    own_sequences = load_sequences_str(args["own_sequences"])

    found = find_cross_sequences(external_sequences, own_sequences, args["max_dist"])

    open(args["output_summary"], "w") do f
        CSV.write(f, DataFrame(found))
    end
end


function find_cross_sequences(external_seqs, own_seqs, max_dist::Int)
    distances::Matrix{Union{Int8,Missing}} =
        fill(Int8(127), length(external_seqs), length(own_seqs))
    Threads.@threads for i in eachindex(external_seqs)
        @simd for j in eachindex(own_seqs)
            @inbounds distances[i, j] = Edlib.edit_distance(
                external_seqs[i].sequence,
                own_seqs[j].sequence,
                mode = :infix,
                max_distance = 127,
            )
        end
    end

    similar = []
    for i in eachindex(external_seqs)
        cr_dists = @view distances[i, :]

        if any(x -> !ismissing(x), cr_dists)
            cr_dists_no_missing = skipmissing(cr_dists)
            min_dist = minimum(cr_dists_no_missing)
        else
            min_dist = nothing
        end

        if !isnothing(min_dist) && min_dist <= max_dist
            candidates = findall(==(min_dist), cr_dists_no_missing)

            best_candidate = nothing
            # iterate through all the candidates to check if any of them has the same cell
            for c in candidates
                same_cell = check_same_cell(external_seqs[i].id, own_seqs[c].id)
                best_candidate = (
                    cellranger = external_seqs[i].id,
                    ours = own_seqs[c].id,
                    dist = distances[i, c],
                    same_cell = same_cell,
                )
                # if it has the same cell id, then make it automatically the best candidate
                if same_cell
                    break
                end
            end
            push!(similar, best_candidate)
        end
    end

    return similar
end


function load_sequences_str(filename::String)
    open(filename) do f
        return [
            (id = identifier(x), sequence = string(sequence(x))) for x in FASTA.Reader(f)
        ]
    end
end


function check_same_cell(cr_id::String, own_id::String)
    return replace(cr_id, r"-1.*$" => "") == replace(own_id, r"_.*$" => "")
end


function parse_arguments()
    s = ArgParseSettings(description = "Blabla")

    @add_arg_table! s begin
        #! format: off
        "--output-summary", "-o"
            help = "output summary file (CSV)"
            default = "align_summary.csv"
            dest_name = "output_summary"
        "--max-dist", "-d"
            help = "maximum distance against external consensus"
            arg_type = Int
            default = 2
            dest_name = "max_dist"
        "own_sequences"
            help = "FASTA file with sequences (own)"
            required = true
        "external_sequences"
            help = "FASTA file with external sequences"
            required = true
        #! format: on
    end

    return parse_args(s)
end


main()