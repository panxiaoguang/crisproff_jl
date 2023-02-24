using JLD2
using BioSequences
using GZip
using Fire
using FASTX
using ViennaRNA_jll
##############################################################
#################### Energy computation #########################
##############################################################
struct CRISPRoff_scores
    Chrom::String
    Start::Int64
    End::Int64
    Strand::Char
    Seq::String
    Eng::Float64
end
const RNA_DNA = JLD2.load_object("energy_dics.jld2")
########## READ THE 2mer 3mer 4mer energies ###############

########## gRNA folding #######################

function get_rnafold_eng(seq::String)
    eng = ccall((:vrna_fold, ViennaRNA_jll.libRNA), Cfloat, (Cstring,), seq)
    Base.cconvert(Float64, eng)
end

# Stacking energy is distributed to positions only for interior loops

function calcRNADNAenergy(guideSeq::String, otSeq::String)
    RI_REV_NT_MAP = Dict('-' => ' ', 'a' => 'T', 'A' => 'T', 'c' => 'G', 'C' => 'G', 'g' => 'C', 'G' => 'C',
        't' => 'A', 'T' => 'A', 'u' => 'A', 'U' => 'A', 'n' => 'N', 'N' => 'N')
    RNA_DNA_internal_loop = Dict(4 => 3.2, 5 => 3.555, 6 => 3.725, 7 => 3.975, 8 => 4.16, 9 => 4.33, 10 => 4.495, 11 => 4.6, 12 => 4.7)
    guideSeq = uppercase(guideSeq[1:end-3])
    seq = join([RI_REV_NT_MAP[c] for c in otSeq[1:end-3]])
    spos = -1
    epos = -1
    energy = fill(0.0, length(guideSeq))
    MATCH = Dict('A' => Dict('A' => false, 'C' => false, 'G' => false, 'T' => true),
        'C' => Dict('A' => false, 'C' => false, 'G' => true, 'T' => false),
        'G' => Dict('A' => false, 'C' => true, 'G' => false, 'T' => false),
        'T' => Dict('A' => true, 'C' => false, 'G' => false, 'T' => false))
    for i in 1:length(seq)
        if MATCH[guideSeq[i]][seq[i]]
            if spos == -1
                spos = i
            end
            epos = i
        end
    end
    i = spos
    while i < epos
        j = i + 1
        while !(MATCH[seq[j]][guideSeq[j]])
            j += 1
            if j > epos
                break
            end
        end
        if j > epos
            break
        end
        loop_size = j - i
        eng_con = 0.0
        if loop_size <= 3
            eng_con = RNA_DNA[loop_size][guideSeq[i:j]][seq[i:j]]
            # if there is a stack in the beginning or end AU GU penalty is still needed
            if loop_size == 1
                if (i == spos && (guideSeq[i] == 'T' || seq[i] == 'T')) || (j == epos && (guideSeq[j] == 'T' || seq[j] == 'T'))
                    eng_con += 0.25
                end
            end
        else
            eng_con = RNA_DNA_internal_loop[loop_size] + RNA_DNA[1][guideSeq[i:i+1]][seq[i:i+1]] + RNA_DNA[1][guideSeq[j-1:j]][seq[j-1:j]]
        end
        for k in 1:(loop_size)
            energy[i+k-1] += eng_con / loop_size
        end
        i = j
    end
    energy
end
################ DNA-DNA opening #############

function calcDNAopeningScore(otSeq::String)
    RI_REV_NT_MAP = Dict('-' => ' ', 'a' => 'T', 'A' => 'T', 'c' => 'G', 'C' => 'G', 'g' => 'C', 'G' => 'C',
        't' => 'A', 'T' => 'A', 'u' => 'A', 'U' => 'A', 'n' => 'N', 'N' => 'N')
    RI_DNA_DNA_NN = Dict("AA" => Dict("TT" => -1.00), "TT" => Dict("AA" => -1.00), "AT" => Dict("TA" => -0.88), "TA" => Dict("AT" => -0.58),
        "CA" => Dict("GT" => -1.45), "TG" => Dict("AC" => -1.45), "GT" => Dict("CA" => -1.44), "AC" => Dict("TG" => -1.44),
        "CT" => Dict("GA" => -1.28), "AG" => Dict("TC" => -1.28), "GA" => Dict("CT" => -1.30), "TC" => Dict("AG" => -1.30),
        "CG" => Dict("GC" => -2.17), "GC" => Dict("CG" => -2.24), "GG" => Dict("CC" => -1.84), "CC" => Dict("GG" => -1.84))
    seq = uppercase(otSeq[1:end-3])
    energy = fill(0.0, length(seq))
    for i in 2:length(seq)
        energy[i] = RI_DNA_DNA_NN[seq[i-1]*seq[i]][RI_REV_NT_MAP[seq[i-1]]*RI_REV_NT_MAP[seq[i]]]
    end
    energy
end
############# Get interaction energy ##################
#Employ all the necessary computations on the score vector to get the final free energy
function get_eng(grna_seq::String, off_seq::String; pos_weight=false, pam_corr=false, grna_folding=false, dna_opening=false, dna_pos_wgh=false)
    POS_WGH = [1.80067099242007, 1.95666668400006, 1.90472004401173, 2.13047270152512, 1.37853848098249, 1.46460783730408, 1.0, 1.387220146823, 1.51401000729362, 1.98058344620751, 1.87939168587699, 1.7222593588838, 2.02228445489326, 1.92692086621503, 2.08041972716723, 1.94496755678903, 2.14539112893591, 2.04277109036766, 2.24911493451185, 2.25]
    DNA_POS_WGH = [1.22245576981774, 1.24561578622024, 1.37883177517399, 1.39146340276523, 1.24308180746857, 1.09598194424544, 1.0, 1.11695025382169, 1.11589045394936, 1.22243614188218, 1.21317477033274, 1.07125942316357, 1.25205871414019, 1.21445408158483, 1.20971491326295, 1.21076785001579, 1.2480898972246, 1.40301355270318, 1.41221084925493, 1.4]
    pam_ratios = Dict("GGG" => 1.0, "AGG" => 1.0, "CGG" => 1.0, "TGG" => 1.0, "GAG" => 0.9, "AAG" => 0.9, "CAG" => 0.9, "TAG" => 0.9, "GGA" => 0.8, "AGA" => 0.8, "CGA" => 0.8, "TGA" => 0.8, "OTHERS" => 0.0)
    pam_ratio_count = 3
    scores = calcRNADNAenergy(grna_seq, off_seq)
    if pos_weight
        for i in 1:length(scores)
            if i < 21
                scores[end-(i-1)] = POS_WGH[end-(i-1)] * scores[end-(i-1)]
            end
        end
    end
    off = -sum(scores)
    if grna_folding
        off += get_rnafold_eng(grna_seq[1:20])
    end
    if dna_opening
        dna_scores = calcDNAopeningScore(off_seq)
        if dna_pos_wgh
            for i in 1:length(dna_scores)
                if i < 21
                    dna_scores[end-(i-1)] = DNA_POS_WGH[end-(i-1)] * dna_scores[end-(i-1)]
                end
            end
        end
        off += sum(dna_scores)
    end
    if pam_corr
        if off_seq[end-pam_ratio_count+1:end] in collect(keys(pam_ratios))
            off = off * pam_ratios[off_seq[end-pam_ratio_count+1:end]]
        else
            off = off * pam_ratios["OTHERS"]
        end
    end
    off
end
## Compute POFF for the given grna and its off-targets ##
function compute_CRISPRspec(ontarget::CRISPRoff_scores, offSeqs::Vector{CRISPRoff_scores}; pos_weight=false, pam_corr=false, grna_folding=false, dna_opening=false, dna_pos_wgh=false)
    pf = zero(Float64)
    on = zero(Float64)
    PAR_BETA = 1.0 / (0.001987 * 310.15)
    nums = length(offSeqs) + 1
    scores = Vector{CRISPRoff_scores}(undef, nums)
    i = 1
    for offSeq in offSeqs
        offSeq_eng = get_eng(ontarget.Seq, offSeq.Seq; pos_weight=pos_weight, pam_corr=pam_corr, grna_folding=grna_folding, dna_opening=dna_opening, dna_pos_wgh=dna_pos_wgh)
        pf = pf + exp(PAR_BETA * offSeq_eng)
        scores[i] = CRISPRoff_scores(offSeq.Chrom, offSeq.Start, offSeq.End, offSeq.Strand, offSeq.Seq, offSeq_eng)
        i += 1
    end
    on_eng = get_eng(ontarget.Seq, ontarget.Seq; pos_weight=pos_weight, pam_corr=pam_corr, grna_folding=grna_folding, dna_opening=dna_opening, dna_pos_wgh=dna_pos_wgh)
    on = exp(PAR_BETA * on_eng)
    scores[i] = CRISPRoff_scores(ontarget.Chrom, ontarget.Start, ontarget.End, ontarget.Strand, ontarget.Seq, on_eng)
    pf / (pf + on), scores
end

function comp_seq(seq::AbstractString)
    REV_NT_MAP = Dict('-' => ' ', 'a' => 'T', 'A' => 'T', 'c' => 'G', 'C' => 'G', 'g' => 'C', 'G' => 'C',
        't' => 'A', 'T' => 'A', 'u' => 'A', 'U' => 'A', 'n' => 'N', 'N' => 'N')
    s = ""
    for c in seq
        s = s * REV_NT_MAP[c]
    end
    s
end

function read_risearch_results(guideSeq::String, ris_file::String; noPAM_given::Bool=false, count_mms::Bool=true)
    off_counts = Dict("GG" => fill(0, 7), "AG" => fill(0, 7), "GA" => fill(0, 7))
    offSeqs = CRISPRoff_scores[]
    on_targets = CRISPRoff_scores[]
    inf = GZip.open(ris_file)
    while !eof(inf)
        line = readline(inf)
        cols = split(line, "\t")
        if length(cols) > 11 && length(cols[11]) > 3 && length(cols[12]) > 0
            gid, qs, qe, tc, ts, te, tst, en, ist, iseq, pamseq, preseq = cols[1:12]
            PAM = comp_seq(pamseq[1:3])
            ts = parse(Int64, ts)
            te = parse(Int64, te)
            if PAM[2:3] in ["GG", "AG", "GA"]
                offseq = comp_seq(iseq)
                if tst == "+"
                    ts = ts - 4
                    tst = "-"
                else
                    ts = ts - 1
                    te = te + 3
                    tst = "+"
                end
                # determine and save the number of mismatches for this guide
                mm_count = sum(offseq[i] != guideSeq[i] for i in 1:length(offseq))
                if count_mms
                    if mm_count < 7
                        off_counts[PAM[2:3]][mm_count+1] += 1
                    end
                end
                if mm_count < 7
                    if (noPAM_given && (guideSeq != offseq)) || ((!noPAM_given) && (guideSeq != offseq * PAM))
                        push!(offSeqs, CRISPRoff_scores(tc, ts, te, tst[1], offseq * PAM, 0.0))
                    else
                        push!(on_targets, CRISPRoff_scores(tc, ts, te, tst[1], offseq * PAM, 0.0))
                    end
                end
            end
        else
            @info "WARNING: RIsearch output line is missing columns, or info in columns.\nLINE: $(line)"
        end
    end
    close(inf)
    offSeqs, off_counts, on_targets
end

"""
Usage:

```
julia -t <Threads> CRISPRspec_CRISPRoff_pipeline.jl --guides <guides.fa> --risearch_results_folder <folder> --CRISPRoff_scores_folder <folder> --specificity_report <file>
```
This is a software build by Xiaoguang based on the python version!
"""
@main function compute_CRISPRoff(; guides::String="none", risearch_results_folder::String="none", CRISPRoff_scores_folder::String="none", specificity_report::String="none")
    guideSeqs = collect(open(FASTA.Reader, guides))
    open(specificity_report, "w") do f
        println(f, join(["Guide_ID", "Guide_sequence", "Genomic_position", "CRISPRspec_specificity_score", "MM_counts", "MM_detailed"], "\t"))
        Threads.@threads for guide in guideSeqs
            GuideID = FASTX.identifier(guide)
            GuideSeq = convert(String, FASTA.sequence(guide))
            risearch_file = joinpath(risearch_results_folder, "risearch_$(GuideID).out.gz")
            offSeqs, off_counts, ontargets = read_risearch_results(GuideSeq, risearch_file)
            on_prob, CRISPRoff_scored_offs = compute_CRISPRspec(ontargets[1], offSeqs; pos_weight=true, pam_corr=true, grna_folding=true, dna_opening=true, dna_pos_wgh=false)
            CRISPRspec = on_prob > 0 ? string(-1.0 * log10(on_prob)) : "INF"
            offscore_file = joinpath(CRISPRoff_scores_folder, "$(ontargets[1].Seq).CRISPRoff.tsv")
            open(offscore_file, "w") do IO
                println(IO, join(["chromosome", "start", "end", "off_target_seq", "CRISPRoff_score", "strand"], "\t"))
                for scored in CRISPRoff_scored_offs
                    println(IO, join([scored.Chrom, string(scored.Start), string(scored.End), scored.Seq, string(scored.Eng), string(scored.Strand)], "\t"))
                end
            end
            MM_counts = join(string.(map(x -> sum(x), [map(x -> x[i], values(off_counts)) for i = 1:7])), ",")
            MM_detailed = join([string("N", pam, ":", join(string.(off_counts[pam]), ",")) for pam = ["GG", "AG", "GA"]], ";")
            for ontarget in ontargets
                coord = string(ontarget.Chrom, ":", ontarget.Start, "-", ontarget.End)
                println(f, join([GuideID, ontarget.Seq, coord, CRISPRspec, MM_counts, MM_detailed], "\t"))
            end
        end
    end
    @info "#END: FINISHED RUNNING THE PIPELINE WITH NO ERROR."
end






