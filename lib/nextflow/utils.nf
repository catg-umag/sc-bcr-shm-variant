def build_cutadapt_args(primers) {
    cmd = primers.collect {
        primer -> primer.usages.collect {
            format_cutadapt_primer(
                primer.name, primer.seq, it.position, it.read, it.rc, it.anchored
            )
        }.join(' ')
    }.join(' ')

    return cmd
}

def format_cutadapt_primer(name, seq, position, read, rc, anchored) {
    options = [1: [front: 'g', back: 'a'], 2: [front: 'G', back: 'A']]
    
    anchor = anchored ? (position == 'front' ? '^_seq_' : '_seq_$') : '_seq_'
    bc_seq = anchor.replace('_seq_', rc ? reverse_complement(seq) : seq)
    
    return "-${options[read][position]} ${name}${rc ? '_rc' : ''}=${bc_seq}"
}

def reverse_complement(seq) {
    dna_rev = [A: 'T', C: 'G', G: 'C', T: 'A']
    return seq.reverse().split('').collect { dna_rev[it] }.join('')
}
