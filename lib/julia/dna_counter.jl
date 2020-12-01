using BioSequences

"""
A structure designed to count DNA bases
"""
mutable struct DNACounter{T<:Number}
    dna_bases::Vector{DNA}
    counters::Vector{T}
    sum::T

    DNACounter(;
        dna_bases::Vector{DNA} = [DNA_A, DNA_C, DNA_G, DNA_T],
        default_value::Number = 0.0,
    ) = new{typeof(default_value)}(
        dna_bases,
        fill(default_value, length(dna_bases)),
        default_value * length(dna_bases),
    )
end


"""
    add!(counter, base, value)

Adds a value to the counter of one nucleotide
"""
@inline function add!(
    counter::DNACounter{T},
    base::DNA,
    value::T = convert(T, 1),
) where {T<:Number}
    for i = 1:length(counter.dna_bases)
        if counter.dna_bases[i] == base
            @inbounds counter.counters[i] += value
            counter.sum += value
            break
        end
    end
end


"""
    reset!(counter, value)

Sets nucleotide counters and sum to a given value
"""
@inline function reset!(counter::DNACounter{T}, value::T = convert(T, 0)) where {T<:Number}
    @simd for i = 1:length(counter.dna_bases)
        @inbounds counter.counters[i] = value
    end
    counter.sum = value
end
