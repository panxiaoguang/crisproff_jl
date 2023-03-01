#!/bin/bash

julia CRISPRspec_CRISPRoff_pipeline.jl --guides test_data/gRNAs_with_pam.fa --risearch-results-folder test_data/ --CRISPRoff-scores-folder test_data.example.out --specificity-report test_data.example.out/test_CRISPRspec.jl.txt

## compare the two files
diff test_data.example.out/test_CRISPRspec.jl.txt test_data/test_CRISPRspec.txt > /dev/null

if [ $0 -eq 0 ]; then
echo "Test passed"
else
echo "Test failed"
fi

## clean up
rm -rf test_data.example.out/*