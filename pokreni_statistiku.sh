#!/usr/bin/bash

sp_score_RV11=0
tc_score_RV11=0
sp_score_nj_RV11=0
tc_score_nj_RV11=0
sp_score_clustal_RV11=0
tc_score_clustal_RV11=0
brojac_RV11=0
sp_score_RV12=0
tc_score_RV12=0
sp_score_nj_RV12=0
tc_score_nj_RV12=0
sp_score_clustal_RV12=0
tc_score_clustal_RV12=0
brojac_RV12=0
sp_score_RV20=0
tc_score_RV20=0
sp_score_nj_RV20=0
tc_score_nj_RV20=0
sp_score_clustal_RV20=0
tc_score_clustal_RV20=0
brojac_RV20=0
sp_score_RV30=0
tc_score_RV30=0
sp_score_nj_RV30=0
tc_score_nj_RV30=0
sp_score_clustal_RV30=0
tc_score_clustal_RV30=0
brojac_RV30=0
sp_score_RV40=0
tc_score_RV40=0
sp_score_nj_RV40=0
tc_score_nj_RV40=0
sp_score_clustal_RV40=0
tc_score_clustal_RV40=0
brojac_RV40=0
sp_score_RV50=0
tc_score_RV50=0
sp_score_nj_RV50=0
tc_score_nj_RV50=0
sp_score_clustal_RV50=0
tc_score_clustal_RV50=0
brojac_RV50=0

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

    brojac_RV11=$((brojac_RV11 + 1))
done

for tfa_file in bb3_release/RV12/*.tfa; do
    base_name=$(basename "$tfa_file" .tfa)
    echo "$base_name"
    
    # Generate MSA and calculate scores
    build/msa "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV12/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "normalni sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_RV12=$(echo "$sp_score_RV12 + $sp_score" | bc)
    tc_score_RV12=$(echo "$tc_score_RV12 + $tc_score" | bc)

    # Generate MSA and calculate scores
    build/msa -n "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV12/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "nj sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_nj_RV12=$(echo "$sp_score_nj_RV12 + $sp_score" | bc)
    tc_score_nj_RV12=$(echo "$tc_score_nj_RV12 + $tc_score" | bc)

    clustalw "$tfa_file" > /dev/null
    sed -i '1a\\//' bb3_release/RV12/${base_name}.aln
    bb3_release/bali_score_src/bali_score bb3_release/RV12/${base_name}.msf bb3_release/RV12/${base_name}.aln > pom.txt
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "clustalw sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_clustal_RV12=$(echo "$sp_score_clustal_RV12 + $sp_score" | bc)
    tc_score_clustal_RV12=$(echo "$tc_score_clustal_RV12 + $tc_score" | bc)

    brojac_RV12=$((brojac_RV12 + 1))
done

for tfa_file in bb3_release/RV20/*.tfa; do
    base_name=$(basename "$tfa_file" .tfa)
    echo "$base_name"
    
    # Generate MSA and calculate scores
    build/msa "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV20/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "normalni sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_RV20=$(echo "$sp_score_RV20 + $sp_score" | bc)
    tc_score_RV20=$(echo "$tc_score_RV20 + $tc_score" | bc)

    # Generate MSA and calculate scores
    build/msa -n "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV20/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "nj sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_nj_RV20=$(echo "$sp_score_nj_RV20 + $sp_score" | bc)
    tc_score_nj_RV20=$(echo "$tc_score_nj_RV20 + $tc_score" | bc)

    clustalw "$tfa_file" > /dev/null
    sed -i '1a\\//' bb3_release/RV20/${base_name}.aln
    bb3_release/bali_score_src/bali_score bb3_release/RV20/${base_name}.msf bb3_release/RV20/${base_name}.aln > pom.txt
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "clustalw sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_clustal_RV20=$(echo "$sp_score_clustal_RV20 + $sp_score" | bc)
    tc_score_clustal_RV20=$(echo "$tc_score_clustal_RV20 + $tc_score" | bc)

    brojac_RV20=$((brojac_RV20 + 1))
done

for tfa_file in bb3_release/RV30/*.tfa; do
    base_name=$(basename "$tfa_file" .tfa)
    echo "$base_name"
    
    # Generate MSA and calculate scores
    build/msa "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV30/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "normalni sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_RV30=$(echo "$sp_score_RV30 + $sp_score" | bc)
    tc_score_RV30=$(echo "$tc_score_RV30 + $tc_score" | bc)

    # Generate MSA and calculate scores
    build/msa -n "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV30/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "nj sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_nj_RV30=$(echo "$sp_score_nj_RV30 + $sp_score" | bc)
    tc_score_nj_RV30=$(echo "$tc_score_nj_RV30 + $tc_score" | bc)

    clustalw "$tfa_file" > /dev/null
    sed -i '1a\\//' bb3_release/RV30/${base_name}.aln
    bb3_release/bali_score_src/bali_score bb3_release/RV30/${base_name}.msf bb3_release/RV30/${base_name}.aln > pom.txt
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "clustalw sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_clustal_RV30=$(echo "$sp_score_clustal_RV30 + $sp_score" | bc)
    tc_score_clustal_RV30=$(echo "$tc_score_clustal_RV30 + $tc_score" | bc)

    brojac_RV30=$((brojac_RV30 + 1))
done

for tfa_file in bb3_release/RV40/*.tfa; do
    base_name=$(basename "$tfa_file" .tfa)
    echo "$base_name"
    
    # Generate MSA and calculate scores
    build/msa "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV40/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "normalni sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_RV40=$(echo "$sp_score_RV40 + $sp_score" | bc)
    tc_score_RV40=$(echo "$tc_score_RV40 + $tc_score" | bc)

    # Generate MSA and calculate scores
    build/msa -n "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV40/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "nj sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_nj_RV40=$(echo "$sp_score_nj_RV40 + $sp_score" | bc)
    tc_score_nj_RV40=$(echo "$tc_score_nj_RV40 + $tc_score" | bc)

    clustalw "$tfa_file" > /dev/null
    sed -i '1a\\//' bb3_release/RV40/${base_name}.aln
    bb3_release/bali_score_src/bali_score bb3_release/RV40/${base_name}.msf bb3_release/RV40/${base_name}.aln > pom.txt
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "clustalw sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_clustal_RV40=$(echo "$sp_score_clustal_RV40 + $sp_score" | bc)
    tc_score_clustal_RV40=$(echo "$tc_score_clustal_RV40 + $tc_score" | bc)

    brojac_RV40=$((brojac_RV40 + 1))
done

for tfa_file in bb3_release/RV50/*.tfa; do
    base_name=$(basename "$tfa_file" .tfa)
    echo "$base_name"
    
    # Generate MSA and calculate scores
    build/msa "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV50/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "normalni sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_RV50=$(echo "$sp_score_RV50 + $sp_score" | bc)
    tc_score_RV50=$(echo "$tc_score_RV50 + $tc_score" | bc)

    # Generate MSA and calculate scores
    build/msa -n "$tfa_file" > /dev/null
    bb3_release/bali_score_src/bali_score "bb3_release/RV50/${base_name}.msf" izgraden_msa.txt > pom.txt
    
    # Extract scores
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "nj sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_nj_RV50=$(echo "$sp_score_nj_RV50 + $sp_score" | bc)
    tc_score_nj_RV50=$(echo "$tc_score_nj_RV50 + $tc_score" | bc)

    clustalw "$tfa_file" > /dev/null
    sed -i '1a\\//' bb3_release/RV50/${base_name}.aln
    bb3_release/bali_score_src/bali_score bb3_release/RV50/${base_name}.msf bb3_release/RV50/${base_name}.aln > pom.txt
    sp_score=$(grep -oP 'SP score=\s*\K\d+\.\d+' pom.txt)
    echo "clustalw sp score $sp_score"
    tc_score=$(grep -oP 'TC score=\s*\K\d+\.\d+' pom.txt)

    sp_score_clustal_RV50=$(echo "$sp_score_clustal_RV50 + $sp_score" | bc)
    tc_score_clustal_RV50=$(echo "$tc_score_clustal_RV50 + $tc_score" | bc)

    brojac_RV50=$((brojac_RV50 + 1))
done

echo "RV11: SP=$sp_score_RV11/$brojac_RV11, TC=$tc_score_RV11/$brojac_RV11"
echo "RV11_nj: SP=$sp_score_nj_RV11/$brojac_RV11, TC=$tc_score_nj_RV11/$brojac_RV11"
echo "RV11_clustal: SP=$sp_score_clustal_RV11/$brojac_RV11, TC=$tc_score_clustal_RV11/$brojac_RV11"

echo "RV12: SP=$sp_score_RV12/$brojac_RV12, TC=$tc_score_RV12/$brojac_RV12"
echo "RV12_nj: SP=$sp_score_nj_RV12/$brojac_RV12, TC=$tc_score_nj_RV12/$brojac_RV12"
echo "RV12_clustal: SP=$sp_score_clustal_RV12/$brojac_RV12, TC=$tc_score_clustal_RV12/$brojac_RV12"

echo "RV20: SP=$sp_score_RV20/$brojac_RV20, TC=$tc_score_RV20/$brojac_RV20"
echo "RV20_nj: SP=$sp_score_nj_RV20/$brojac_RV20, TC=$tc_score_nj_RV20/$brojac_RV20"
echo "RV20_clustal: SP=$sp_score_clustal_RV20/$brojac_RV20, TC=$tc_score_clustal_RV20/$brojac_RV20"

echo "RV30: SP=$sp_score_RV30/$brojac_RV30, TC=$tc_score_RV30/$brojac_RV30"
echo "RV30_nj: SP=$sp_score_nj_RV30/$brojac_RV30, TC=$tc_score_nj_RV30/$brojac_RV30"
echo "RV30_clustal: SP=$sp_score_clustal_RV30/$brojac_RV30, TC=$tc_score_clustal_RV30/$brojac_RV30"

echo "RV40: SP=$sp_score_RV40/$brojac_RV40, TC=$tc_score_RV40/$brojac_RV40"
echo "RV40_nj: SP=$sp_score_nj_RV40/$brojac_RV40, TC=$tc_score_nj_RV40/$brojac_RV40"
echo "RV40_clustal: SP=$sp_score_clustal_RV40/$brojac_RV40, TC=$tc_score_clustal_RV40/$brojac_RV40"

echo "RV50: SP=$sp_score_RV50/$brojac_RV50, TC=$tc_score_RV50/$brojac_RV50"
echo "RV50_nj: SP=$sp_score_nj_RV50/$brojac_RV50, TC=$tc_score_nj_RV50/$brojac_RV50"
echo "RV50_clustal: SP=$sp_score_clustal_RV50/$brojac_RV50, TC=$tc_score_clustal_RV50/$brojac_RV50"
