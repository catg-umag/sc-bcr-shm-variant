#!/usr/bin/env python3
import csv
import re
from os import makedirs, path
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq

EXP_SHEET = "$FL_PROJ/005_Sample_sheet/Experiment_sheet.txt"
EXT_CONSENSUS = "$FL_PROJ/050_BCR_ext_ref.freeze01/FL_PO.csv"
INT_CONSENSUS = "$FL_PROJ/250_Clonotypes/outs/consensus10X_20200731.fasta"
OUTPUT_DIR = "../input/references/"
IMMCANTATION_IMAGE = "docker://kleinstein/immcantation:3.0.0"


def main():
    makedirs(path.expandvars(OUTPUT_DIR), exist_ok=True)
    process_external_consensus(EXT_CONSENSUS, EXP_SHEET, OUTPUT_DIR)
    process_internal_consensus(INT_CONSENSUS, OUTPUT_DIR)


def load_subjects(filename):
    with open(path.expandvars(filename)) as f:
        d = {
            x["subject"]: x["chip"]
            for x in csv.DictReader(f, delimiter="\t")
            if x["analysis"] == "BCR" and x["chip"] in ("K1B", "K2B")
        }

    return d


def process_external_consensus(filename, subjects_filename, output_dir):
    subjects = load_subjects(subjects_filename)

    with open(path.expandvars(filename)) as f:
        sequences = {x["subject"]: x["sequence"] for x in csv.DictReader(f)}

    with open(
        path.join(path.expandvars(output_dir), "consensus_ext_vdj.fasta"), "w"
    ) as f:
        for seqid, seq in sequences.items():
            subject = seqid.split("_")[0]
            if subject in subjects:
                new_id = f"{subjects[subject]}_" + seqid.replace("_PO", "")
                f.write(f">{new_id}\n{seq}\n")


def process_internal_consensus(filename, output_dir):
    sequences = [x for x in SeqIO.parse(path.expandvars(filename), "fasta")]

    for seq in sequences:
        subject, experiment, chain = re.match(r"^(S\d+)_(K..)_(IG.V)", seq.id).groups()
        seq.id = f"{experiment}_{subject}_{'HC' if 'IGH' in chain else 'LC'}"
        seq.description = seq.id
        seq.seq = Seq(str(seq.seq).replace("?", "N"))

    sequences = sorted(sequences, key=lambda x: get_sortable_id(x.id))
    with open(
        path.join(path.expandvars(output_dir), "consensus_int_vdj.fasta"), "w"
    ) as f:
        SeqIO.write(sequences, f, "fasta-2line")

    with open(
        path.join(path.expandvars(output_dir), "consensus_int_k3b_vdj.fasta"), "w"
    ) as f:
        SeqIO.write([x for x in sequences if "K3B" in x.id], f, "fasta-2line")


def get_sortable_id(seq_id):
    exp, sample, chain = re.match("K(\d)B_S(\d+)_(.C)", seq_id).groups()
    return exp + sample.zfill(6) + chain


def get_consensus_vdj(filename, output_dir):
    """
    Extracts VDJ region from full consensus using IgBlast
    """
    singularity_cmd = ["singularity", "exec", IMMCANTATION_IMAGE]
    exec_cmd = ["changeo-igblast", "-s", path.basename(filename)]
    base_dir = path.dirname(filename)
    basename = path.splitext(path.basename(filename))[0]
    subprocess.run(singularity_cmd + exec_cmd, cwd=base_dir)

    db_filename = path.join(base_dir, basename, f"{basename}_db-pass.tab")
    with open(db_filename) as f:
        sequences = [x for x in csv.DictReader(f, delimiter="\t")]

    with open(
        path.join(path.expandvars(output_dir), "consensus_ext_vdj.fasta"), "w"
    ) as f:
        for seq in sequences:
            f.write(f">{seq['SEQUENCE_ID']}\n{seq['SEQUENCE_VDJ']}\n")


if __name__ == "__main__":
    main()
