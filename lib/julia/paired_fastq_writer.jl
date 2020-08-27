using FASTX


mutable struct PairedFastqWriter
    r1::FASTQ.Writer
    r2::FASTQ.Writer
    info::IOStream

    PairedFastqWriter(prefix) = begin
        x = new(
            FASTQ.Writer(open("$(prefix).r1.fastq", "w")),
            FASTQ.Writer(open("$(prefix).r2.fastq", "w")),
            open("$(prefix)_info.csv", "w"),
        )
        write(x.info, "experiment,id,cell,umi\n")
        x
    end
end


@inline function write_record(
    writer::PairedFastqWriter,
    r1::FASTQ.Record,
    r2::FASTQ.Record,
    experiment::String,
    cell::Union{String,Nothing},
    umi::Union{String,Nothing},
    extra_offset::Int = 20,
)
    # delete this from 5' R1 (cell BC, UMI and some extra nucleotides)
    offset = 26 + extra_offset + 1

    write(
        writer.r1,
        FASTQ.Record(
            identifier(r1),
            description(r1),
            sequence(r1)[offset:end],
            quality(r1)[offset:end],
        ),
    )
    write(writer.r2, r2)
    write(writer.info, "$(experiment),$(identifier(r1)),$(cell),$(umi)\n")
end


@inline function close_writer(writer::PairedFastqWriter)
    close(writer.r1)
    close(writer.r2)
    close(writer.info)
end
