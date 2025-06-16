#!/usr/bin/bash

sp_score_uk=0
tc_score_uk=0
sp_score_nj_uk=0
tc_score_nj_uk=0
sp_score_clustal_uk=0
tc_score_clustal_uk=0
brojac_uk=0

for tfa_file in $1/*.tfa; do
    base_name=$(basename "$tfa_file" .tfa)
    echo "$base_name"
    
    build/msa "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "${1}/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_uk=$(echo "$sp_score_uk + $sp_score" | bc)
    tc_score_uk=$(echo "$tc_score_uk + $tc_score" | bc)

    build/msa -n "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "${1}/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_nj_uk=$(echo "$sp_score_nj_uk + $sp_score" | bc)
    tc_score_nj_uk=$(echo "$tc_score_nj_uk + $tc_score" | bc)

    clustalw "$tfa_file" > /dev/null
    sed -i '1a\\//' ${1}/${base_name}.aln
    bb3_release/bali_score_src/bali_score ${1}/${base_name}.msf ${1}/${base_name}.aln > pom.txt
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_clustal_uk=$(echo "$sp_score_clustal_uk + $sp_score" | bc)
    tc_score_clustal_uk=$(echo "$tc_score_clustal_uk + $tc_score" | bc)

    brojac_uk=$((brojac_uk + 1))
done

echo "uk: SP=$sp_score_uk/$brojac_uk, TC=$tc_score_uk/$brojac_uk"
echo "uk_nj: SP=$sp_score_nj_uk/$brojac_uk, TC=$tc_score_nj_uk/$brojac_uk"
echo "uk_clustal: SP=$sp_score_clustal_uk/$brojac_uk, TC=$tc_score_clustal_uk/$brojac_uk"