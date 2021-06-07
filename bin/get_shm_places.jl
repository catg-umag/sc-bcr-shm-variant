#!/usr/bin/env julia
using ArgMacros
using BioSequences
using CSV
using FASTX
using Query


function main()
    args = parse_arguments()

    # load cell references
    cell_references = Dict(x.cell => x.reference for x in CSV.File(args.cell_references))

    # load reference regions
    reference_regions = load_reference_regions(args.reference_regions)

    # load reference V-genes aligned positions
    reference_vgene_positions = load_aligned_v_positions(args.aligned_v_regions)

    # load consensuses
    data = CSV.File(args.input) |> @filter(_.ref_vdj_coverage == 1)

    # load consensus productivity
    consensus_producive = Dict(
        row.sequence_id => row.productive for row in CSV.File(args.consensus_productivity)
    )

    # get cells that complies with the min UMI criteria
    selected_cells =
        data |>
        @groupby(_.cell) |>
        @map({cell = key(_), n = length(_)}) |>
        @filter(_.n >= args.min_umis) |>
        @map(_.cell) |>
        collect

    variants = []
    for cell in selected_cells
        cell_data = data |> @filter(_.cell == cell) |> collect
        sequences = cell_data |> @map(_.aligned_consensus) |> collect
        umis = cell_data |> @map(_.umi) |> collect
        depths = cell_data |> @map(parse.(Int64, split(_.depths, ";"))) |> collect
        cell_reference = cell_references[cell]
        cell_regions = reference_regions[cell_reference]

        for i in eachindex(first(sequences))
            # only valid nucleotides
            nucleotides = [x[i] for x in sequences if x[i] in "ACGT"]

            if length(unique(nucleotides)) > 1
                counts = count_nucleotides(nucleotides)
                ratios = sort(collect(values(counts)); rev = true) ./ length(nucleotides)
                cons_nucl = get_consensus_nucleotide(counts)

                if length(nucleotides) >= args.min_umis && ratios[2] >= args.min_ratio
                    region = get_variant_region(i, cell_regions)
                    if region == "V"
                        position_in_v = i - cell_regions["V"][1] + 1
                        position_in_v_aln =
                            reference_vgene_positions[cell_reference][position_in_v]
                    else
                        position_in_v = position_in_v_aln = ""
                    end

                    for (umi, seq, d) in zip(umis, sequences, depths)
                        if seq[i] in "ACGT"
                            productive = get(consensus_producive, "$(cell)_$(umi)", "")
                            variant = (
                                cell,
                                umi,
                                i,
                                cons_nucl,
                                seq[i],
                                get_context(seq, i),
                                d[i],
                                productive,
                                region,
                                position_in_v,
                                position_in_v_aln,
                            )
                            push!(variants, variant)
                        end
                    end
                end
            end
        end
    end

    # write 
    open(args.output, "w") do f
        write(
            f,
            "cell,umi,position,cons_nucl,nucl,context,depth,productive," *
            "region,vgene_position,vgene_position_aligned\n",
        )
        for v in variants
            write(f, join(string.(v), ",") * "\n")
        end
    end
end


"""
    get_variant_region(position, regions)

Uses information contained in regions to get the Ig Region associated with a position
"""
function get_variant_region(position::Int, regions::Dict)::String
    # if regions empty, return empty
    isempty(regions) && return ""

    for (region, range) in regions
        if range[1] <= position <= range[2]
            return region
        end
    end

    # if the position was not covered, the is should be Leader, Constant or N (N1/N2)
    if position < regions["V"][1]
        return "L"
    elseif position > regions["J"][2]
        return "C"
    else
        return "N"
    end
end


"""
    load_reference_regions(filename)

Loads reference regions in a multi-level dict
"""
function load_reference_regions(filename::String)::Dict
    reference_regions = Dict()
    for row in CSV.File(filename)
        reference_regions[row.sequence_id] = Dict()
        for r in ("v", "d", "j")
            rstart = Symbol("$(r)_start")
            rend = Symbol("$(r)_end")
            if !ismissing(row[rstart]) && !ismissing(row[rend])
                reference_regions[row.sequence_id][uppercase(r)] = (row[rstart], row[rend])
            end
        end
    end

    return reference_regions
end


"""
    load_aligned_v_positions(filename)

Loads aligned reference V-genes matching the alined and unaligned position
"""
function load_aligned_v_positions(filename::String)::Dict
    positions = Dict()
    reader = open(FASTA.Reader, filename)
    record = FASTA.Record()
    while !eof(reader)
        read!(reader, record)
        ref_positions = Dict()
        unaligned_pos = 0
        for (i, nucl) in enumerate(FASTA.sequence(record))
            if nucl != DNA_Gap
                unaligned_pos += 1
                ref_positions[unaligned_pos] = i
            end
        end
        positions[FASTA.identifier(record)] = ref_positions
    end

    return positions
end


"""
    get_consensus_nucleotide(counts)

Gets most frequent nucleotide from counts dictionary
"""
function get_consensus_nucleotide(counts::Dict)::Char
    isempty(counts) && return 'N'

    cons_nucl = nothing
    for (n, c) in counts
        if isnothing(cons_nucl) || cons_nucl[2] < c
            cons_nucl = (n, c)
        end
    end

    return cons_nucl[1]
end


"""
    count_nucleotides(bases)

Counts nucleotides in a vector, returning a dict with the counts
"""
function count_nucleotides(bases::Vector{Char})::Dict{Char,Int}
    counts = Dict('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0)
    for b in bases
        counts[b] += 1
    end

    return counts
end


"""
    get_context(sequence, pos[, size])

Get position context from sequence
"""
function get_context(sequence::String, pos::Int, size::Int = 2)
    context = ""
    if pos > 1
        context *= sequence[max(1, pos - size):max(1, pos - 1)]
    end
    context *= "."
    if pos < length(sequence)
        context *=
            sequence[min(length(sequence), pos + 1):min(length(sequence), pos + size)]
    end

    return context
end


function parse_arguments()
    ratio_validator = x -> 0 <= x <= 0.5

    args = @tuplearguments begin
        @helpusage """
            get_shm_places.jl -i INPUT -o OUTPUT -c CELL_REFERENCES -g REFERENCE_REGIONS
                [-u MIN_UMIS] [-r MIN_RATIO]"""
        @helpdescription "Search places of possible SHM occurence"

        @argumentrequired String input "-i" "--input"
        @arghelp "consensus list by cell/UMI (.csv)"
        @argtest input isfile "The input must be a valid file"

        @argumentrequired String output "-o" "--output"
        @arghelp "output file (.csv)"

        @argumentrequired String reference_regions "-g" "--reference-regions"
        @arghelp "CSV file containing the positions of the regions (V, D, J, etc) for each reference"
        @argtest reference_regions isfile "the reference regions must be a valid file"

        @argumentrequired String cell_references "-c" "--cell-references"
        @arghelp "CSV file containing the reference used for each cell"
        @argtest cell_references isfile "the cell references must be a valid file"

        @argumentrequired String aligned_v_regions "-v" "--aligned-v-regions"
        @arghelp "FASTA file with MSA of reference V regions"
        @argtest aligned_v_regions isfile "the aligned V regions must be a valid file"

        @argumentrequired String consensus_productivity "-p" "--consensus-productivity"
        @arghelp "CSV file with the productive status of the consensus"
        @argtest consensus_productivity isfile "the consensus productivity must be a valid file"

        @argumentdefault Int 2 min_umis "-u" "--min-umis"
        @arghelp "minimum amount of UMIs supporting cell"

        @argumentdefault Float64 0.0 min_ratio "-r" "--min-ratio"
        @arghelp "minimum alt allele ratio to consider it a valid case"
        @argtest min_ratio ratio_validator "Value must be between 0 and 0.5"
    end

    return args
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end