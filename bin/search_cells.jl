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


"""
    get_distances(external_seq, own_seqs)

Calculates and returns all the distances of `external_seq` against `own_seqs`
"""
function get_distances(external_seq, own_seqs)::Vector{Union{Int8,Missing}}
    distances::Vector{Union{Int8,Missing}} = fill(Int8(127), length(own_seqs))
    @simd for i in eachindex(own_seqs)
        @inbounds distances[i] = get_minimum_distance(
            external_seq.sequence,
            own_seqs[i].sequence
        )
    end

    return distances
end


"""
    get_minimum_distance(s1, s2)

Gets edit distance between `s1` and `s2` without counting edges.
To determine the start/end, a number of consecutive matches are required.
"""
function get_minimum_distance(s1::AbstractString, s2::AbstractString)
    # distance higher than 40 should mean that they are really different,
    # so continuing with the alignment is not worth it
    aln = Edlib.alignment(s1, s2; max_distance=40)
    if !ismissing(aln)
        alignment = aln.alignment
        range = (
            low = find_consecutive_end_matches(alignment),
            high = find_consecutive_end_matches(alignment; direction_fw = false),
        )

        return aln.distance - (range.low[2] + range.high[2])
    end
    return missing
end


"""
    find_consecutive_end_matches(alignment, direction_fw, required_matches)

Finds consecutive matches on some ends alignment.
"""
function find_consecutive_end_matches(
    alignment::Vector{Edlib.Alignment};
    direction_fw::Bool = true,
    required_matches::Integer = 5,
)::Tuple{Integer, Integer}
    range = direction_fw ? (1:length(alignment)) : (length(alignment):-1:1)
    consecutive_matches = 0
    missmatches = 0
    for i in range
        if alignment[i] == Edlib.MATCH
            consecutive_matches += 1
            if consecutive_matches == required_matches
                return (i - (required_matches - 1) * (direction_fw ? 1 : -1), missmatches)
            end
        else
            missmatches += 1
            if consecutive_matches > 0
                consecutive_matches = 0
            end
        end
    end
end


"""
    load_sequences_str(filename)

Loads all sequences from a FASTA file as strings
"""
function load_sequences_str(filename::String)
    open(filename) do f
        return [
            (id = identifier(x), sequence = sequence(String, x)) for x in FASTA.Reader(f)
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