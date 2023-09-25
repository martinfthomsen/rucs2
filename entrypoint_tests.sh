#!/bin/bash
# This script will test the following entrypoints: full, fucs, fppp, vpcr, anno, pcrs, expl, test

# CMD to run: bash entrypoint_tests.sh
# ENTRY point tests:

# Test 1. vpcr
mkdir -p script_tests/test1 && cp testdata/test1_pair_file.tsv  script_tests/test1/
docker run --rm -v `pwd`/script_tests/test1:/workdir -v $BLASTDB:/blastdb rucs vpcr --pairs pair_file.tsv --references CP063056 ARBW00000000 BBXJ00000000 MLAJ00000000 MLHJ00000000 MLHH00000000

# Test 2. vpcr
mkdir -p script_tests/test2 && cp testdata/test1_pair_file.tsv  script_tests/test2/
docker run --rm -v `pwd`/script_tests/test2:/workdir -v $BLASTDB:/blastdb rucs vpcr --pairs pair_file.tsv --references JWIZ01 JADGLC01 CAJUGW ASM228557v1 GCA_002285575.1 CP000672.1

# Evaluate tests:
cmp -s testdata/test1_products.tsv script_tests/test1/results/products.tsv && echo -e "\x1B[32mTest 1 - Passed \x1B[0m" || echo -e "\x1B[31mTest 1 - Failed! \x1B[0m"
cmp -s testdata/test2_products.tsv script_tests/test2/results/products.tsv && echo -e "\x1B[32mTest 2 - Passed \x1B[0m" || echo -e "\x1B[31mTest 2 - Failed! \x1B[0m"
