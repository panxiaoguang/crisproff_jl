#!/bin/bash

julia CRISPRspec_CRISPRoff_pipeline.jl --guides test_data/grnas_with_pam.fa --risearch-results-folder test_data/ --CRISPRoff-scores-folder test_data.example.out --specificity-report test_data.example.out/test_CRISPRspec.jl.txt

## compare the two files
diff test_data.example.out/test_CRISPRspec.jl.txt test_data/test_CRISPRspec.jl.txt > /dev/null
if [ $? -eq 0 ]; then
    echo "Test passed"
else
    echo "Test failed"
fi

## clean up#
rm -rf test_data.example.out/*
