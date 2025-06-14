#!/usr/bin/bash

sp_score_RV11=0
tc_score_RV11=0
sp_score_nj_RV11=0
tc_score_nj_RV11=0
sp_score_clustal_RV11=0
tc_score_clustal_RV11=0
brojac=0

for tfa_file in bb3_release/RV11/*.tfa; do
    base_name=$(basename "$tfa_file" .tfa)
    echo "$base_name"
    
    # Generate MSA and calculate scores
    build/msa "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV11/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "normalni sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_RV11=$(echo "$sp_score_RV11 + $sp_score" | bc)
    tc_score_RV11=$(echo "$tc_score_RV11 + $tc_score" | bc)

    # Generate MSA and calculate scores
    build/msa -n "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV11/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "nj sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_nj_RV11=$(echo "$sp_score_nj_RV11 + $sp_score" | bc)
    tc_score_nj_RV11=$(echo "$tc_score_nj_RV11 + $tc_score" | bc)

    clustalw "$tfa_file" > /dev/null
    sed -i '1a\\//' bb3_release/RV11/${base_name}.aln
    bb3_release/bali_score_src/bali_score bb3_release/RV11/${base_name}.msf bb3_release/RV11/${base_name}.aln > pom.txt
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "clustalw sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_clustal_RV11=$(echo "$sp_score_clustal_RV11 + $sp_score" | bc)
    tc_score_clustal_RV11=$(echo "$tc_score_clustal_RV11 + $tc_score" | bc)

    brojac=$((brojac + 1))
done
echo "RV11: SP=$sp_score_RV11/$brojac, TC=$tc_score_RV11/$brojac"
echo "RV11_nj: SP=$sp_score_nj_RV11/$brojac, TC=$tc_score_nj_RV11/$brojac"
echo "RV11_clustal: SP=$sp_score_clustal_RV11/$brojac, TC=$tc_score_clustal_RV11/$brojac"
