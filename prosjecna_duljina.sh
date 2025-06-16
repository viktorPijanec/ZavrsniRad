#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Korištenje: $0 <direktorij>"
    exit 1
fi

total_seqs=0
sum_lengths=0

for file in "$1"/*.tfa; do
    if [ ! -f "$file" ]; then
        continue
    fi
    
    read seq_count seq_length <<< $(awk '
    BEGIN {count=0; total=0; current_len=0} 
    /^>/ {
        if(count>0) total += current_len
        count++
        current_len=0
        next
    } 
    {
        current_len += length($0)
    } 
    END {
        if(count>0) total += current_len
        print count, total
    }' "$file")

    total_seqs=$((total_seqs + seq_count))
    sum_lengths=$((sum_lengths + seq_length))
done

if [ $total_seqs -eq 0 ]; then
    echo "Nema sekvenci u datotekama."
    exit 0
fi

average_length=$(echo "scale=2; $sum_lengths / $total_seqs" | bc)
echo "Ukupan broj sekvenci u svim datotekama: $total_seqs"
echo "Prosječna duljina svih sekvenci: $average_length"