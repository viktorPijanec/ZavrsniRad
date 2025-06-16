#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Korištenje: $0 <datoteka.tfa>"
    exit 1
fi

input_file="$1"

if [ ! -f "$input_file" ]; then
    echo "Greška: Datoteka '$input_file' ne postoji ili nije regularna datoteka"
    exit 2
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
}' "$input_file")

if [ $seq_count -eq 0 ]; then
    echo "Nema sekvenci u datoteci."
    exit 0
fi

average_length=$(echo "scale=2; $seq_length / $seq_count" | bc)
echo "Statistika za datoteku: $input_file"
echo "------------------------------"
echo "Ukupan broj sekvenci: $seq_count"
echo "Ukupna duljina svih sekvenci: $seq_length"
echo "Prosječna duljina sekvenci: $average_length"
