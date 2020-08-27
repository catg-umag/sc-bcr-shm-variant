using XAM

SeqRecord = NamedTuple{
    (:contig, :startpos, :length, :sequence, :quality),
    Tuple{String,Int,Int,LongDNASeq,Vector{UInt8}},
}


"""
    bam_record_to_seq(record)

Transforms a BAM record to a tuple containing only essential info.
Also modifies the sequence (and corresponding quality) applying the CIGAR info,
to being able to be compared directly with the reference.
"""
function bam_record_to_seq(record::BAM.Record)::SeqRecord
    sequence = BAM.sequence(record)
    quality = BAM.quality(record)
    ops, ns = BAM.cigar_rle(record)

    # "apply" cigar
    seq_idx = 1
    for (op, n) in zip(ops, ns)
        if op in (OP_SOFT_CLIP, OP_INSERT)
            deleteat!(sequence, seq_idx:seq_idx+n-1)
            deleteat!(quality, seq_idx:seq_idx+n-1)
        elseif op == OP_DELETE
            for _ = 1:n
                insert!(sequence, seq_idx, DNA_Gap)
                insert!(quality, seq_idx, 0x0)
            end
            seq_idx += n
        elseif op == OP_MATCH
            seq_idx += n
        end
    end

    return (
        contig = BAM.refname(record),
        startpos = BAM.position(record),
        length = length(sequence),
        sequence = sequence,
        quality = quality,
    )
end