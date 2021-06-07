#!/bin/bash
awk -F "\"*,\"*" '
NR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
    }
}
{ if (NR!=1) { printf ">%s\n%s\n", $(f["sequence_id"]), $(f["v_sequence"]) }}
' $1