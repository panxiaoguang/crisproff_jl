### crisproff_jl

This is a Julia implementation for [crisproff](https://github.com/RTH-tools/crisproff). For more details on the algorithm, please see their official repository.

### Usage 

#### Prerequisites

- Julia 
- JLD2
- BioSequences
- GZip
- Fire
- FASTX
- ViennaRNA_jll

You can download Julia from [julialang](https://julialang.org/downloads/) and install these dependencies by simply typing commands like:

```julia
## type  ] to enter Pkg mode 
add JLD2
```

#### method 
In order to process tons of gRNAs, the Julia version comes up here. So the software only supports FASTA sequences as input. In addition, I also think that you already use [search2](https://rth.dk/resources/risearch/) to search the genome to find all off-targets.

```bash
julia CRISPRspec_CRISPRoff_pipeline.jl --guides <guides.fa> --risearch-results-folder <folder> --CRISPRoff-scores-folder <folder> --specificity-report <file>
```

You can find the usage by:

```bash
julia CRISPRspec_CRISPRoff_pipeline.jl --help
```

Thanks for Julia language, the pipeline can be used at multi-theads mode, you can add `-t <threads>` after julia.

```bash
julia -t <Threads> CRISPRspec_CRISPRoff_pipeline.jl --guides <guides.fa> --risearch-results-folder <folder> --CRISPRoff-scores-folder <folder> --specificity-report <file>
```

#### Finally

It should have higher performance than python version!

#### contact

xiaoguang.pan@bio.ku.dk

