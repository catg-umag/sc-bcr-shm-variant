"""
A structure designed to count DNA bases
"""
mutable struct DNACounter
    dna_bases::Vector{DNA}
    counters::Vector{Float64}
    sum::Float64

    DNACounter(;
        dna_bases::Vector{DNA} = [DNA_A, DNA_C, DNA_G, DNA_T],
        default_value::Number = 0.0,
    ) = new(
        dna_bases,
        fill(default_value, length(dna_bases)),
        default_value * length(dna_bases),
    )
end


"""
    add_to_dna_counter!(counter, base, value)

Adds a value to the counter of one nucleotide
"""
@inline function add_to_dna_counter!(counter::DNACounter, base::DNA, value::Number)
    for i = 1:length(counter.dna_bases)
        if counter.dna_bases[i] == base
            @inbounds counter.counters[i] += value
            counter.sum += value
            break
        end
    end
end


"""
    reset_dna_counter!(counter, value)

Sets nucleotide counters and sum to a given value
"""
@inline function reset_dna_counter!(counter::DNACounter, value::Number = 0.0)
    @simd for i = 1:length(counter.dna_bases)
        @inbounds counter.counters[i] = value
    end
    counter.sum = value
end
