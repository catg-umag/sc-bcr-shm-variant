#!/usr/bin/env julia
using Base.Threads
using ArgParse, Edlib, FASTX


function main()
    args = parse_arguments()

    external_sequences = load_sequences_str(args["external_sequences"])
    own_sequences = load_sequences_str(args["own_sequences"])

    open(args["output_summary"], "w") do f
        write(f, "cellranger,ours,dist,same_cell\n")
        Threads.@threads for i in eachindex(external_sequences)
            best_candidate = find_best_cantidadate(
                external_sequences[i],
                own_sequences,
                args["max_dist"],
            )

            if !isnothing(best_candidate)
                write(f, join(best_candidate, ",") * "\n")
            else
                # write empty line
                write(f, external_sequences[i].id * ",,,\n")
            end
        end
    end
end


function find_best_cantidadate(external_seq, own_seqs, max_dist)
    distances = get_distances(external_seq, own_seqs)

    if any(x -> !ismissing(x), distances)
        dists_no_missing = skipmissing(distances)
        min_dist = minimum(dists_no_missing)
    else
        min_dist = nothing
    end

    if !isnothing(min_dist) && min_dist <= max_dist
        candidates = findall(==(min_dist), dists_no_missing)

        best_candidate = nothing
        # iterate through all the candidates to check if any of them has the same cell
        for c in candidates
            same_cell = check_same_cell(external_seq.id, own_seqs[c].id)
            best_candidate = (
                cellranger = external_seq.id,
                ours = own_seqs[c].id,
                dist = distances[c],
                same_cell = same_cell,
            )
            # if it has the same cell id, then make it automatically the best candidate
            if same_cell
                break
            end
        end

        return best_candidate
    end
    
    return nothing
end


function get_distances(external_seq, own_seqs)::Vector{Union{Int8,Missing}}
    distances::Vector{Union{Int8,Missing}} = fill(Int8(127), length(own_seqs))
    @simd for i in eachindex(own_seqs)
        @inbounds distances[i] = Edlib.edit_distance(
            external_seq.sequence,
            own_seqs[i].sequence,
            mode = :infix,
            max_distance = 127,
        )
    end

    return distances
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