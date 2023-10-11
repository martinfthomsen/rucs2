#!/bin/bash
# This script will test the following entrypoints: vpcr, pcrs, fucs, fppp   #full, anno, expl, test

# CMD to run: bash entrypoint_tests.sh
# ENTRY point tests:

# Test 1. vpcr
mkdir -p script_tests/test1/inputs && cp testdata/test1_pair_file.tsv script_tests/test1/inputs/pair_file.tsv
docker run --rm -v `pwd`/script_tests/test1:/workdir -v $BLASTDB:/blastdb rucs vpcr --pairs inputs/pair_file.tsv --references CP063056 ARBW00000000 BBXJ00000000 MLAJ00000000 MLHJ00000000 MLHH00000000

# Test 2. vpcr
mkdir -p script_tests/test2/inputs && cp testdata/test1_pair_file.tsv script_tests/test2/inputs/pair_file.tsv
docker run --rm -v `pwd`/script_tests/test2:/workdir -v $BLASTDB:/blastdb rucs vpcr --pairs inputs/pair_file.tsv --references JWIZ01 JADGLC01 CAJUGW ASM228557v1 GCA_002285575.1 CP000672.1

# Test 3. pcrs
mkdir -p script_tests/test3/inputs && mkdir -p script_tests/test3/results && cp testdata/test1_pair_file.tsv script_tests/test3/inputs/pair_file.tsv && cp testdata/test3_template.fna.gz script_tests/test3/inputs/template.fna.gz
docker run --rm -v `pwd`/script_tests/test3:/workdir rucs pcrs --pairs inputs/pair_file.tsv --template inputs/template.fna.gz > script_tests/test3/results/terminal_output.txt

# Test 4. fucs
mkdir -p script_tests/test4/inputs
docker run --rm -v `pwd`/script_tests/test4:/workdir -v $BLASTDB:/blastdb rucs fucs --positives ASM1929502v1 ASM1300761v1 NZ_KU341381.1 --negatives ASM584v2 ASM886v2

# Test 5. fppp
mkdir -p script_tests/test5/inputs && cp testdata/test4_unique_core_sequences.disscafs.fa script_tests/test5/inputs/template.fa  && cp script_tests/test4/inputs/* script_tests/test5/inputs/
docker run --rm -v `pwd`/script_tests/test5:/workdir -v $BLASTDB:/blastdb rucs fppp --template inputs/template.fa --positives ASM1929502v1 ASM1300761v1 --negatives ASM584v2 ASM886v2 -v

# Test 6. fppp with probe
mkdir -p script_tests/test6/inputs && cp testdata/test4_unique_core_sequences.disscafs.fa script_tests/test6/inputs/template.fa  && cp script_tests/test4/inputs/* script_tests/test6/inputs/
docker run --rm -v `pwd`/script_tests/test6:/workdir -v $BLASTDB:/blastdb rucs fppp --template inputs/template.fa --positives ASM1929502v1 ASM1300761v1 --negatives ASM584v2 ASM886v2 -v --pick_probe

# Test 7. fppp with reuse
mkdir -p script_tests/test7/inputs script_tests/test7/work && cp script_tests/test6/inputs/* script_tests/test7/inputs/ && cp script_tests/test6/work/*pairs.pkl script_tests/test7/work/
docker run --rm -v `pwd`/script_tests/test7:/workdir -v $BLASTDB:/blastdb rucs fppp --template inputs/template.fa --positives ASM1929502v1 ASM1300761v1 --negatives ASM584v2 ASM886v2 -v --pick_probe --reuse


# Evaluate tests:
cmp -s testdata/test1_products.tsv script_tests/test1/results/products.tsv && echo -e "\x1B[32mTest 1 - Passed \x1B[0m" || echo -e "\x1B[31mTest 1 - Failed! \x1B[0m"
cmp -s testdata/test2_products.tsv script_tests/test2/results/products.tsv && echo -e "\x1B[32mTest 2 - Passed \x1B[0m" || echo -e "\x1B[31mTest 2 - Failed! \x1B[0m"
cmp -s testdata/test3_terminal_output.txt script_tests/test3/results/terminal_output.txt && echo -e "\x1B[32mTest 3 - Passed \x1B[0m" || echo -e "\x1B[31mTest 3 - Failed! \x1B[0m"
cmp -s testdata/test4_unique_core_sequences.disscafs.fa script_tests/test4/results/unique_core_sequences.disscafs.fa && echo -e "\x1B[32mTest 4 - Passed \x1B[0m" || echo -e "\x1B[31mTest 4 - Failed! \x1B[0m"
cmp -s testdata/test5_results_best.tsv script_tests/test5/results/results_best.tsv && echo -e "\x1B[32mTest 5 - Passed \x1B[0m" || echo -e "\x1B[31mTest 5 - Failed! \x1B[0m"
cmp -s testdata/test6_results_best.tsv script_tests/test6/results/results_best.tsv && echo -e "\x1B[32mTest 6 - Passed \x1B[0m" || echo -e "\x1B[31mTest 6 - Failed! \x1B[0m"
cmp -s testdata/test6_results_best.tsv script_tests/test7/results/results_best.tsv && echo -e "\x1B[32mTest 7 - Passed \x1B[0m" || echo -e "\x1B[31mTest 7 - Failed! \x1B[0m"
