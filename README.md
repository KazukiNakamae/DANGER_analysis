# DANGER Analysis
Risk-averse on/off-target assessment for CRISPR editing without reference genome. 

<img src="https://github.com/KazukiNakamae/DANGER_analysis/blob/main/images/DANGERAnalysis.png" alt="DANGERAnalysis_logo" title="DANGERAnalysis_logo" width="440" height="155">

DANGER (**D**eleterious and **AN**ticipatable **G**uides **E**valuated by **R**NA-sequencing) Analysis; a bioinformatics pipeline can elucidate genomic on/off-target sites on mRNA-transcribed regions related to expression changes and then quantify phenotypic risk at the gene ontology (GO) term level using RNA-seq data. 

DANGER Analysis would be helpful for people who want to
- screen off-target sites on mRNA-transcribed regions without reference genome.
- find off-target sites resulting in decreased expression using RNA-seq data.
- annotate off-target genes with GO terms
- get personal transcriptome-based on/off-target profile
- provide transcriptome-aware on/off-target profile without reference genome.
- search a potential on-target on expressed genes without reference genome.
- quantify gene expressions of on/off-target sites.

If you have a question and request, don't hesitate to contact me.

Kazuki Nakamae, Ph.D.
- kazukinakamae[at mark]gmail.com

## References

[1] Nakamae and Bono, DANGER analysis: risk-averse on/off-target assessment for CRISPR editing without a reference genome. Bioinformatics Advances, vbad114. 2023 https://doi.org/10.1093/bioadv/vbad114

[2] Nakamae and Bono, DANGER analysis: Risk-averse on/off-target assessment for CRISPR editing without a reference genome. bioRxiv. 2023 https://doi.org/10.1101/2023.03.11.531115

## Run DANGER Analysis using Docker images (Recommended)

Please prepare the following files in your current directory
- "fltr_lowexpr_dj1_trinity_out_dir.Trinity.fasta": de novo transcriptome assembly (without redundancy)
- "guide_pam.fa": The binding sequence of protospacer & PAM (e.g. GCCGGTTCAGTGCAGCCGTGAGG)
- Each "RSEM.isoforms.results": the expression profiles exported from Trinity:align_and_estimate_abundance.pl

<img src="https://github.com/KazukiNakamae/DANGER_analysis/blob/main/images/example_fileset.png" alt="example_fileset" title="example_fileset" height="400">

The example dataset was deposited in [SourceForge](https://sourceforge.net/projects/danger-analysis-v1/files/example.tar.gz/download).

### Step1: Summarize expression profiles

Run collect_exp_data.py in Docker image.

```
sudo docker run --name example_collectexp --memory 10g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangeranalysis:1.1 python /tmp/collect_exp_data.py \
-o <Output directory> \
-w <"RSEM.isoforms.results" in the expression profiles of WT samples> \
-e <"RSEM.isoforms.results" in the expression profiles of Edited samples> \
-t <Threshold for Edited/WT ratio>;
```

(EXAMPLE)
```bash
sudo docker run --name example_collectexp --memory 10g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangeranalysis:1.1 python /tmp/collect_exp_data.py \
-o exp_collection \
-w rmrna_dj1_ctrl_rep1/RSEM.isoforms.results rmrna_dj1_ctrl_rep2/RSEM.isoforms.results rmrna_dj1_ctrl_rep3/RSEM.isoforms.results \
-e rmrna_dj1_ko_rep1/RSEM.isoforms.results rmrna_dj1_ko_rep2/RSEM.isoforms.results rmrna_dj1_ko_rep3/RSEM.isoforms.results \
-t 2.5; # TPM ratio < 1/2.5 = 0.4
```

The output:

- ctrl_edited_fltrexpr_contig_expected_count_onratio.csv:
Comma-separated text file, including read count values of each sample, Edited/WT ratio, and expression labels

- ctrl_edited_fltrexpr_contig_fpkm_onratio.csv:
Comma-separated text file, including FPKM values of each sample, Edited/WT ratio, and expression labels

- ctrl_edited_fltrexpr_contig_tpm_onratio.csv:
Comma-separated text file, including TPM values of each sample, Edited/WT ratio, and expression labels


#### collect_exp_data.py options

```
-o:     Output directory. The name should be unique.
-w:     "RSEM.isoforms.results" in the expression profiles of WT samples.
-e:     "RSEM.isoforms.results" in the expression profiles of Edited samples
-t:     Threshold for Edited/WT ratio. Edited/WT > (Threshold) is "upregulated." Edited/WT < 1/(Threshold) is "downregulated."
```

#### (Optional): Use expression profiles with DEG

"DANGER Analysis can utilize the DEG profile instead of the aforementioned expression profiles."

DEG profile can be made according to the below example.

[DEG analysis of de novo transcriptome assembly from Trinity](https://github.com/KazukiNakamae/DEG_analysis)

"The DEG profile can be adapted to the input format of DANGER analysis by processing it as follows."
```
cat DEG_profile.csv | tr -d '"' > DEG_profile.csv_clean.csv;
echo 'id,name,<ctrl sample 1>,<ctrl sample 2>,...,<ctrl sample n>,<edited sample 1>,<edited sample 2>,...,<edited sample n>,gene_id,a.value,m.value,p.value,q.value,rank,estimatedDEG,Exp' > DEG_data_header.csv
cp DEG_data_header.csv DEG_profile_fixed.csv
tail -n +2 DEG_profile.csv_clean.csv >> DEG_profile_fixed.csv
### Use p < <alpha>
awk -F, 'BEGIN{OFS=","} {if (NR==1)print $0;else {if($14 < <alpha> && $13 < 0) print $0, "downregulated";else if($14 < <alpha> && $13 > 0) print $0, "upregulated"; else print $0, "unchanged"}}' \
DEG_profile_fixed.csv \
> complete_DEG_profile.csv;
```

### Step2: GO Annotation & D-index Calculation

Run dangeranalysis_v1.sh in Docker image.

```
sudo docker run --name example_danalysis --memory 100g --rm -v `pwd`:/DATA -w /tmp -i kazukinakamae/dangeranalysis:1.1 bash /tmp/dangeranalysis_v1.sh \
<Database for GO annotation> \
<Database Type> \
<Output directory> \
<expression profile> \
<de novo transcriptome assembly (without redundancy)> \
<binding sequence of protospacer & PAM> \
<mismatch number for off-target search> \
<PAM for off-target search>;
```

(EXAMPLE)
```bash
sudo docker run --name example_danalysis --memory 100g --rm -v `pwd`:/DATA -w /tmp -i kazukinakamae/dangeranalysis:1.1 bash /tmp/dangeranalysis_v1.sh \
Dr \
pep \
output \
exp_collection/ctrl_edited_fltrexpr_contig_tpm_onratio.csv \
fltr_lowexpr_dj1_trinity_out_dir.Trinity.fasta \
guide_pam.fa \
11 \
NGG;
```

#### dangeranalysis_v1.sh options

```
Database for GO annotation:    Available database for GO annotation
Users can choose here:
- Hs : Human
- Mm : Mouse
- Dm : Fly
- Ce : C.elegans
- Dr : Zebrafish
- Ag : Mosquito
- At : A.thaliana
- Bt : Cow
- Cf : Dog
- EcK12 : E.coli K-12
- EcSakai E.coli O157:H7 str. Sakai
- Gg : Chicken
- Mmu : Monkey
- Pt : Chimpanzee
- Rn : Rat
- Sc : Yeast
- Ss : Pig
- Xl : Frog
- Mxanthus : M.xanthus

Database Type:     Sequence used for gene annotation
- pep : Protein sequence
- cdna : mRNA sequence

Output directory:     Output directory. The name should be a unique

summary of the TPM:     Comma-separated text file, including TPM values of each sample, Edited/WT ratio, and expression labels

de novo transcriptome assembly (without redundancy):     de novo transcriptome assembly from Trinity

binding sequence of protospacer & PAM:     Fasta file including single sequence of on-target site

mismatch number for off-target search:     mismatch number used in CrisFlash. It can be 1-11 nt

PAM for off-target search:   PAM sequence used in CrisFlash. SpCas9 is generally NGG. If you want to consider minor PAM in SpCas9, you can choose NRR.
```

The result is saved in DANGER_analysis_result.

The table of D-indice will be saved as "DANGER_index_on_***.txt".

<img src="https://github.com/KazukiNakamae/DANGER_analysis/blob/main/images/DANGER_index_on_Molecular_Function.txt.png" alt="Result tabele" title="Result tabele" height="300">

#### Troubleshooting

- If you are a non-root user, the command can fail to download databases. We recommend that the root user runs the command.

### Step3: D-index Evaluation

In DANGER Analysis, it is possible to use permutation data to evaluate the significance of the calculated D-index. The procedure for this is described below.

First, create expression permutation data. If you want to create 100 instances, input as follows

```
mkdir exp_p_data;
sudo docker run --name genarate_exp_p_data --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
bash /tmp/get_multi_exp_permutation_data.sh \
<original summary of TPM> \
<First seed number> \
<Last seed number> \
exp_p_data/exp;
```

(EXAMPLE)
```bash
mkdir example_exp_p_data;
sudo docker run --name genarate_exp_p_data --memory 10g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
bash /tmp/get_multi_exp_permutation_data.sh \
exp_collection/ctrl_edited_fltrexpr_contig_tpm_onratio.csv \
1 \
100 \
example_exp_p_data/exp;
```

The output data will be outputted as a directory, with the seed value added to the prefix(e.g. park7_2_5_exp_p_data/exp1, park7_2_5_exp_p_data/exp2, ...).

Second, create off-target permutation data. If you want to create 100 instances, input as follows

```
mkdir mm_p_data;
sudo docker run --name genarate_mm_p_data --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
bash /tmp/get_multi_mm_permutation_data.sh \
<Original off-target profile> \
<First seed number> \
<Last seed number> \
mm_p_data/mm;
```

(EXAMPLE)
```bash
mkdir example_mm_p_data;
sudo docker run --name genarate_mm_p_data --memory 10g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
bash /tmp/get_multi_mm_permutation_data.sh \
output/Result_offtarget_all.cas-offinder \
1 \
100 \
example_mm_p_data/mm;
```

The output data will be outputted as a directory, with the seed value added to the prefix(e.g. park7_mm8_NGG_data/mm1, park7_mm8_NGG_data/mm2, ...).

Next, calculate the pseudo D-index based on the permutation data.

```bash
mkdir p_result;
for seed in $(seq <First seed number> <Last seed number>)
do
        sudo docker run \
                --name example_dangeranalysis_false_p${seed} --memory 20g --rm -v `pwd`:/DATA -w /tmp -i kazukinakamae/dangertest:1.0 bash /tmp/dangertest_v2.sh \
                <Database for GO annotation> \
                <Database Type> \
                p_result/p${seed} \
                exp_p_data/exp_p${seed}/ctrl_edited_fltrexpr_contig_tpm_onratio.csv \
                mm_p_data/mm_p${seed}/Result_offtarget_all.cas-offinder \
                <Original directory of DANGER analysis> \
                <binding sequence of protospacer & PAM>;
done
```

(EXAMPLE)
```bash
mkdir false_output;
for seed in $(seq 1 100)
do
  sudo docker run \
          --name example_dangeranalysis_false_p${seed} --memory 20g --rm -v `pwd`:/DATA -w /tmp -i kazukinakamae/dangertest:1.0 bash /tmp/dangertest_v2.sh \
          Dr \
          pep \
          false_output/p${seed} \
          example_exp_p_data/exp_p${seed}/ctrl_edited_fltrexpr_contig_tpm_onratio.csv \
          example_mm_p_data/mm_p${seed}/Result_offtarget_all.cas-offinder \
          output \
          guide_pam.fa;
done
```

D-indices greater than the specified confidence interval from the t-distribution calculated from the pseudo D-indices are considered significant

```
sudo docker run \
--name example_dindextest --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
python /tmp/dindex_test.py \
-t <Original directory of DANGER analysis> \
-n <directory of pseudo DANGER analysis> \
-c <Confidence interval> \
-o <Output directory>;
```

(EXAMPLE)
```bash
sudo docker run \
--name example_dindextest --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
python /tmp/dindex_test.py \
-t output \
-n false_output \
-c 0.999999999999999 \
-o example_dindex_test_n100_ci99_9999999999999percent;
```


The table of significant D-indice will be saved as "Significant_DANGER_index_on_***.txt".

(EXAMPLE)

<img src="https://github.com/KazukiNakamae/DANGER_analysis/blob/main/images/Significant_DANGER_index_on_Molecular_Function.txt.png" alt="Result tabele" title="Result tabele" height="100">

The table will be visualized and saved as "Significant_DANGER_index_on_***.tiff."

(EXAMPLE)

<img src="https://github.com/KazukiNakamae/DANGER_analysis/blob/main/images/Significant_DANGER_index_on_Molecular_Function.tiff.png" alt="Result tabele" title="Result tabele" height="300">


### (Optional) Step4: False positive test of D-index Evaluation

In the D-index evaluation method up to Step 3, parameters such as the number of mismatches to consider, the TPM threshold, PAM sequences, confidence intervals, and permutation data can be decided at the user's discretion. However, it is also possible to measure false positives based on these results. The procedure for this is described below.

First, create expression permutation data for the test. If you want to create 10 instances, input as follows
**An important note is that this permutation data must be created in a separate directory from the permutation data created in Step 3, and the seed value must be set to a number that does not duplicate with that in Step 3.**

```bash
mkdir falsepos_exp_p_data;
sudo docker run \
            --name falsepos_park7_2_5_exp_p_data --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
            bash /tmp/get_multi_exp_permutation_data.sh \
            <original expression profile> \
            <First seed number> \
            <Last seed number> \
            falsepos_exp_p_data/<Output prefix>;
```

(EXAMPLE)
```bash
mkdir falsepos_example_exp_p_data;
sudo docker run --name falsepos_example_exp_p_data --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
bash /tmp/get_multi_exp_permutation_data.sh \
exp_collection/ctrl_edited_fltrexpr_contig_tpm_onratio.csv \
1001 \
1010 \
falsepos_example_exp_p_data/exp;
```

Second, create off-target permutation data. If you want to create 100 instances, input as follows
**An important note is that this permutation data must be created in a separate directory from the permutation data created in Step 3, and the seed value must be set to a number that does not duplicate with that in Step 3.**

```bash
mkdir falsepos_mm_p_data;
sudo docker run \
            --name falsepos_park7_mm8_NGG_data --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
            bash /tmp/get_multi_mm_permutation_data.sh \
            <Original off-target profile> \
            <First seed number> \
            <Last seed number> \
            falsepos_mm_p_data/mm;
```

(EXAMPLE)
```bash
mkdir falsepos_example_mm_p_data;
sudo docker run --name falsepos_example_mm_p_data --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
bash /tmp/get_multi_mm_permutation_data.sh \
output/Result_offtarget_all.cas-offinder \
1001 \
1010 \
falsepos_example_mm_p_data/mm;
```

Next, calculate the pseudo D-index based on the above permutation data.
**The paramter shoud be same as Step3**

```bash
mkdir falsepos_p_result;
for seed in $(seq <First seed number> <Last seed number>)
do
        sudo docker run \
                --name falsepos_park7_2_5_random_mm8_NGG_dangeranalysis_p${seed} --memory 20g --rm -v `pwd`:/DATA -w /tmp -i kazukinakamae/dangertest:1.0 bash /tmp/dangertest_v2.sh \
                <Database for GO annotation> \
                <Database Type> \
                falsepos_p_result/p${seed} \
                falsepos_exp_p_data/exp_p${seed}/ctrl_edited_fltrexpr_contig_tpm_onratio.csv \
                falsepos_mm_p_data/mm_p${seed}/Result_offtarget_all.cas-offinder \
                <Original directory of DANGER analysis> \
                <binding sequence of protospacer & PAM>;
done
```

(EXAMPLE)
```bash
mkdir falsepos_outputs;
for seed in $(seq 1001 1010)
do
        sudo docker run \
                --name falsepos_park7_2_5_random_mm8_NGG_dangeranalysis_p${seed} --memory 20g --rm -v `pwd`:/DATA -w /tmp -i kazukinakamae/dangertest:1.0 bash /tmp/dangertest_v2.sh \
                Dr \
                pep \
                falsepos_outputs/p${seed} \
                falsepos_example_exp_p_data/exp_p${seed}/ctrl_edited_fltrexpr_contig_tpm_onratio.csv \
                falsepos_example_mm_p_data/mm_p${seed}/Result_offtarget_all.cas-offinder \
                output \
                guide_pam.fa;
done
```

The pseudo D-indices greater than the specified confidence interval from the t-distribution calculated from the pseudo D-indices of Step3 are considered significant
**The paramter shoud be same as Step3**

```bash
mkdir <directory of D-index test using pseudo DANGER analysis>;
for seed in $(seq <First seed number> <Last seed number>)
do
        sudo docker run \
                --name dangeranalysis_dindextest_${seed} --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
                python /tmp/dindex_test.py \
                -t falsepos_p_result/p${seed} \
                -n <directory of pseudo DANGER analysis> \
                -c <Confidence interval> \
                -o <directory of D-index test using pseudo DANGER analysis>/<user-defined label>_${seed};
done
```

(EXAMPLE)
```bash
mkdir falsepos_Dindex_test_n100;
for seed in $(seq 1001 1010)
do
        sudo docker run \
                --name dangeranalysis_dindextest_${seed} --memory 20g --rm -v `pwd`:/DATA -w /DATA -i kazukinakamae/dangertest:1.0 \
                python /tmp/dindex_test.py \
                -t falsepos_outputs/p${seed} \
                -n false_output \
                -c 0.999999999999999 \
                -o falsepos_Dindex_test_n100/99_9999999999999percent_seed_${seed};
done
```

The proportion of the number of D-index detected in 'Significant_DANGER_index_on_XXX.txt' relative to that in 'All_DANGER_index_on_XXX.txt' is calculated as the false positive rate. 


## Run DANGER Analysis in your local environment using conda (Deprecated)

<details>
<summary>Click here</summary>

### Installation of DANGER Analysis

The DANGER Analysis consists of python and R with various bioinformatics tools. All processes run under Anaconda and Docker environments.

#### 1. Installation of Docker

##### on MacOSX

1. Download and install Docker Desktop: https://docs.docker.com/engine/install/#desktop

2. Enter Docker settings menu to adjust the memory allocation (≥64GB of memory is recommmended)

##### on Linux

1. Download and install Docker Engine: https://docs.docker.com/engine/

2. start Docker daemon service

```bash
sudo systemctl start docker;
```

#### 2. Download of Docker images

Type the following commands in the terminal.

```bash
# Download Docker image for Trinity
sudo docker pull trinityrnaseq/trinityrnaseq:2.12.0

# Download Docker image for BUSCO
sudo docker pull trinityrnaseq/ezlabgva/busco:v5.2.2_cv1
```

#### 3. Create of Anaconda environments

Type the following commands in the terminal.

```bash
# Create cutadapt_env
conda activate cutadapt_env -y;
conda install -c bioconda cutadapt=1.18;
conda deactivate;

# Create bbtools_env
conda activate bbtools_env -y;
conda install -c bioconda bbmap=38.18;
conda deactivate;

# Create transdecoder_env
conda create -n transdecoder_env -y;
conda activate transdecoder_env;
conda install -c bioconda -y TransDecoder=5.5.0;
conda install -c conda-forge pigz=2.6;
conda install -c bioconda blast=2.12.0;
conda install -c bioconda seqkit=2.3.1;
conda deactivate;

# Create matplotlib_venn_env
conda create -n matplotlib_venn_env;
conda activate matplotlib_venn_env;
conda install -c conda-forge matplotlib-venn=0.11.5;
conda deactivate;

# Create topGO
conda create -n topGO -y;
conda activate topGO;
conda install -c conda-forge -c bioconda bioconductor-topgo -y;
conda install -c conda-forge -c bioconda bioconductor-rgraphviz -y;
conda install -c bioconda bioconductor-org.hs.eg.db bioconductor-org.ag.eg.db bioconductor-org.at.tair.db bioconductor-org.bt.eg.db bioconductor-org.ce.eg.db bioconductor-org.cf.eg.db bioconductor-org.dm.eg.db bioconductor-org.dr.eg.db bioconductor-org.eck12.eg.db bioconductor-org.ecsakai.eg.db bioconductor-org.gg.eg.db bioconductor-org.mm.eg.db bioconductor-org.mmu.eg.db bioconductor-org.mxanthus.db bioconductor-org.pt.eg.db bioconductor-org.rn.eg.db bioconductor-org.sc.sgd.db bioconductor-org.ss.eg.db bioconductor-org.xl.eg.db;
conda deactivate;

# Create calcDANGERindex_env
conda create -n calcDANGERindex_env -y;
conda activate calcDANGERindex_env;
conda install -c anaconda pandas=1.5.2;
conda install -c anaconda scipy=1.10.0;
conda install -c conda-forge matplotlib=3.6.3;
conda deactivate;
```

#### 4 Install Crisflash

Download source code of Crisflash from https://github.com/crisflash/crisflash, and then install it according to install instructions.

#### 5 Download scripts of DANGER analysis, SAQE, and the supplemantal script.

Type the following commands in the terminal.

```
git clone https://github.com/KazukiNakamae/DANGER_analysis.git;
git clone https://github.com/bonohu/SAQE.git;
git clone https://github.com/RyoNozu/Sequence_editor.git;
```

#### 6. Prepare de novo transcriptome assembly

We show the examples using park7(dj1) dataset.
Type the following commands in the terminal.

```bash
mkdir working_dir
cd working_dir
mkdir raw_fastq
cd raw_fastq

# Move raw fastq(.gz) files to the current directory
mv XXX.fq.gz ./;
cd ..;

# make log directory
mkdir log;

# make list of sample names
mkdir metadata;
cat << EOF > metadata/sample_name.txt
dj1_Control_1
dj1_Control_2
dj1_Control_3
dj1_KO_1
dj1_KO_2
dj1_KO_3
EOF

### Quality Control & Adapter Trimming
conda activate cutadapt_env;
mkdir trimmed_fq;
mkdir resource;
cat << EOF > resource/illumina_universal.fa
>AGATCGGAAGAG
AGATCGGAAGAG
EOF
# The average length is 150nt. We set 120nt as minimum-length
while read line; do cutadapt -q 30 -a file:resource/illumina_universal.fa -A file:resource/illumina_universal.fa -o trimmed_fq/trimmed_${line}_1.fq.gz -p trimmed_fq/trimmed_${line}_2.fq.gz --minimum-length=120 --pair-filter=any --trim-n raw_fastq/${line}_1.fq.gz raw_fastq/${line}_2.fq.gz &> log/1.cutadapt_${line}.txt; done < metadata/sample_name.txt
conda deactivate;



### Ribosomal RNA (rRNA) Removal
mkdir rmrrna_first_fq;
mkdir rmrrna_first_match_fq;
mkdir rmrrna_second_fq;
mkdir rmrrna_second_match_fq;
# Download rRNA datasets from https://www.arb-silva.de, and then save them to resource directory
# Remove SSU rRNA
while read line; do bbduk.sh K=31 mcf=0.5 in1=trimmed_fq/trimmed_${line}_1.fq.gz in2=trimmed_fq/trimmed_${line}_2.fq.gz \
out1=rmrrna_first_fq/rmssu_${line}_1.fq.gz out2=rmrrna_first_fq/rmssu_${line}_2.fq.gz outm=rmrrna_first_match_fq/rmssu_${line}.fq.gz \
ref=resource/SILVA_119.1_SSURef_Nr99_tax_silva_trunc.fasta.gz overwrite=t -Xmx8g &> log/2.rmrrna_first_${line}.txt; done < metadata/sample_name.txt;
# Remove LSU rRNA
while read line; do bbduk.sh K=31 mcf=0.5 in1=rmrrna_first_fq/rmssu_${line}_1.fq.gz in2=rmrrna_first_fq/rmssu_${line}_2.fq.gz \
out1=rmrrna_second_fq/rmrna_${line}_1.fq.gz out2=rmrrna_second_fq/rmrna_${line}_2.fq.gz outm=rmrrna_second_match_fq/rmlsu_${line}.fq.gz \
ref=resource/SILVA_119_LSURef_tax_silva_trunc.fasta.gz overwrite=t -Xmx8g &> log/3.rmrrna_second_${line}.txt; done < metadata/sample_name.txt;
conda deactivate;



### de novo Transcriptome Assembly
# Merge WT data
# Left
mkdir merge_ctrl_1_fq;
cp rmrrna_second_fq/rmrna_dj1_Control_1_1.fq.gz merge_ctrl_1_fq;
cp rmrrna_second_fq/rmrna_dj1_Control_2_1.fq.gz merge_ctrl_1_fq;
cp rmrrna_second_fq/rmrna_dj1_Control_3_1.fq.gz merge_ctrl_1_fq;
cat merge_ctrl_1_fq/rmrna_dj1_Control_*_1.fq.gz > merge_ctrl_1_fq/MERGED_dj1_Control_1.fq.gz;
# Right
mkdir merge_ctrl_2_fq;
cp rmrrna_second_fq/rmrna_dj1_Control_1_2.fq.gz merge_ctrl_2_fq;
cp rmrrna_second_fq/rmrna_dj1_Control_2_2.fq.gz merge_ctrl_2_fq;
cp rmrrna_second_fq/rmrna_dj1_Control_3_2.fq.gz merge_ctrl_2_fq;
cat merge_ctrl_2_fq/rmrna_dj1_Control_*_2.fq.gz > merge_ctrl_2_fq/MERGED_dj1_Control_2.fq.gz;
# Run de novo Transcriptome Assembly
cp merge_ctrl_1_fq/MERGED_dj1_Control_1.fq.gz ../;
cp merge_ctrl_2_fq/MERGED_dj1_Control_2.fq.gz ../;
cd ..;
sudo docker run -d --memory 128g --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq:2.12.0 Trinity --seqType fq --left `pwd`/MERGED_dj1_Control_1.fq.gz --right `pwd`/MERGED_dj1_Control_2.fq.gz --min_contig_length 100 --CPU 12 --max_memory 128G --output `pwd`/dj1_trinity_out_dir;



### Evaluate assembly using BUSCO
cp dj1_trinity_out_dir/Trinity.fasta ./Trinity_original.fasta;
docker run --name busco_original -d --memory 128g -itv $PWD:/data -w /data --rm ezlabgva/busco:v5.2.2_cv1 busco -m transcriptome -i ./Trinity_original.fasta -o dj1_trinity_transcripts_BUSCO -l actinopterygii_odb10 -c 8;
cat dj1_trinity_transcripts_BUSCO/short_summary.specific.actinopterygii_odb10.dj1_trinity_transcripts_BUSCO.txt
###
# BUSCO version is: 5.2.2 
# The lineage dataset is: actinopterygii_odb10 (Creation date: 2021-02-19, number of genomes: 26, number of BUSCOs: 3640)
# Summarized benchmarking in BUSCO notation for file /data/Trinity_original.fasta
# BUSCO was run in mode: transcriptome

	***** Results: *****

	C:90.9%[S:42.1%,D:48.8%],F:2.5%,M:6.6%,n:3640	   
	3310	Complete BUSCOs (C)			   
	1534	Complete and single-copy BUSCOs (S)	   
	1776	Complete and duplicated BUSCOs (D)	   
	90	Fragmented BUSCOs (F)			   
	240	Missing BUSCOs (M)			   
	3640	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.1
	metaeuk: 4.a0f584d
###



### Removal of redundancy
# Estimate expression of merged WT data
docker run --rm -d -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq:2.12.0 /usr/local/bin/util/align_and_estimate_abundance.pl --transcripts `pwd`/dj1_trinity_out_dir/Trinity.fasta \
--seqType fq \
--left `pwd`/MERGED_dj1_Control_1.fq.gz --right `pwd`/MERGED_dj1_Control_2.fq.gz \
 --est_method RSEM \
 --aln_method bowtie2 \
 --trinity_mode \
 --prep_reference \
 --coordsort_bam \
 --thread_count 20 \
 --output_dir `pwd`/dj1_trinity_out_dir.Trinity_RSEM_outdir;
# removal of low-expression transcripts
docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq:2.12.0 /usr/local/bin/util/filter_low_expr_transcripts.pl --transcripts `pwd`/dj1_trinity_out_dir/Trinity.fasta \
 --highest_iso_only \
 --trinity_mode \
 --matrix `pwd`/dj1_trinity_out_dir.Trinity_RSEM_outdir/RSEM.genes.results \
 > `pwd`/fltr_lowexpr_dj1_trinity_out_dir.Trinity.fasta;
```


#### 7. Expression quantification

Type the following commands in the terminal.

```bash
# Make sample lists
cat << EOF > `pwd`/fastq_sample_fltr_lowexpr.txt
ctrl `pwd`/rmrrna_second_fq/rmrna_dj1_ctrl_rep1 `pwd`/rmrrna_second_fq/rmrna_dj1_Control_1_1.fq.gz `pwd`/rmrrna_second_fq/rmrna_dj1_Control_1_2.fq.gz
ctrl `pwd`/rmrrna_second_fq/rmrna_dj1_ctrl_rep2 `pwd`/rmrrna_second_fq/rmrna_dj1_Control_2_1.fq.gz `pwd`/rmrrna_second_fq/rmrna_dj1_Control_2_2.fq.gz
ctrl `pwd`/rmrrna_second_fq/rmrna_dj1_ctrl_rep3 `pwd`/rmrrna_second_fq/rmrna_dj1_Control_3_1.fq.gz `pwd`/rmrrna_second_fq/rmrna_dj1_Control_3_2.fq.gz
ko `pwd`/rmrrna_second_fq/rmrna_dj1_ko_rep1 `pwd`/rmrrna_second_fq/rmrna_dj1_KO_1_1.fq.gz `pwd`/rmrrna_second_fq/rmrna_dj1_KO_1_2.fq.gz
ko `pwd`/rmrrna_second_fq/rmrna_dj1_ko_rep2 `pwd`/rmrrna_second_fq/rmrna_dj1_KO_2_1.fq.gz `pwd`/rmrrna_second_fq/rmrna_dj1_KO_2_2.fq.gz
ko `pwd`/rmrrna_second_fq/rmrna_dj1_ko_rep3 `pwd`/rmrrna_second_fq/rmrna_dj1_KO_3_1.fq.gz `pwd`/rmrrna_second_fq/rmrna_dj1_KO_3_2.fq.gz
EOF
cat `pwd`/fastq_sample_fltr_lowexpr.txt | tr ' ' '\t' > `pwd`/fastq_sample_fltr_lowexpr_tab.txt;
# Estimate expression of WT/Edited data
sudo docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq:2.12.0 /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts `pwd`/fltr_lowexpr_dj1_trinity_out_dir.Trinity.fasta \
--thread_count 8 \
--prep_reference \
--seqType fq \
--samples_file `pwd`/fastq_sample_fltr_lowexpr_tab.txt \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--coordsort_bam;
``` 

Summarize TPM values into ctrl_ko_fltrexpr_contig_tpm_onratio.csv



#### 8. Detection of on/off-target sites on de novo transcriptome assembly

Type the following commands in the terminal.

``` bash
# Make on-target sequence data
cat << EOT >> guide_pam.fa
>GCCGGTTCAGTGCAGCCGTGAGG
GCCGGTTCAGTGCAGCCGTGAGG
EOT
# on/off-target detection
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' `pwd`/fltr_lowexpr_dj1_trinity_out_dir.Trinity.fasta > Trinity.fa; # Remove \n in sequences
cat Trinity.fa | sed -e 's/ .*//g' > Trinity_simple.fa; # Modify the headers
$HOME/bin/crisflash -g Trinity_simple.fa -s guide_pam.fa -o results_gRNAs.cas-offinder -m 11 -p NRR -t 8 -C;
wc -l results_gRNAs.cas-offinder
# hit site (upto 11mm, NRR PAM) = 4232634 sites

# Get all off-targets 
awk '{ if ($4 != "GCCGGTTCAGTGCAGCCGTGAGG") { print } }' results_gRNAs.cas-offinder > offtarget_all.cas-offinder;
wc -l offtarget_all.cas-offinder;
# off-target site (upto 11mm, NRR PAM) = 4232633 sites

# Get 0-11 off-targets（NRR PAM）
awk -v MM=0 '{ if ($6 == MM) { print } }' offtarget_all.cas-offinder > offtarget_mm0.cas-offinder;
for i in {1..11};do awk -v MM=$i '{ if ($6 == MM) { print } }' offtarget_all.cas-offinder > offtarget_mm"$i".cas-offinder;wc -l offtarget_mm"$i".cas-offinder;done;
###
       0 offtarget_mm1.cas-offinder
       0 offtarget_mm2.cas-offinder
       1 offtarget_mm3.cas-offinder
      25 offtarget_mm4.cas-offinder
     278 offtarget_mm5.cas-offinder
    2179 offtarget_mm6.cas-offinder
   14172 offtarget_mm7.cas-offinder
   73027 offtarget_mm8.cas-offinder
  296222 offtarget_mm9.cas-offinder
 1008070 offtarget_mm10.cas-offinder
 2838659 offtarget_mm11.cas-offinder
###

# Count transcripts have off-target sites
awk '{ print $2 }' offtarget_all.cas-offinder > offtarget_all_list.txt;
awk '!seen[$0]++' offtarget_all_list.txt > offtarget_all_uniq_list.txt;
wc -l offtarget_all_uniq_list.txt;
# 865,452 transcripts have off-target sites
``` 



#### 9. GO enrichment analysis

Type the following commands in the terminal.

```bash
# Predict ORFs
conda activate transdecoder_env;
cd 120221227_dj1_denovo_assembly;
cp ../SAQE/11TransDecoder.sh ./;
cp ./fltr_lowexpr_dj1_trinity_out_dir.Trinity.fasta ./Trinity.fasta;
chmod +x ./11TransDecoder.sh;
./11TransDecoder.sh

# Download Ensembl protein database
mkdir db;
cp /Volumes/denovoseq/SAQE/12GetRefProts.sh ./;
# Prepare modified 12GetRefProts to download zebrafish protein data
# Run script for download 
chmod +x ./12GetRefProts.sh;
./12GetRefProts.sh

# Add dummy word into blank of database
unpigz -c db/Danio_rerio.GRCz11.pep.all.fa.gz > Danio_rerio.GRCz11.pep.all.fa;
cp ../Sequence_editor/00_prepare_faa_4Fanflow.sh ./;
chmod +x 00_prepare_faa_4Fanflow.sh;
./00_prepare_faa_4Fanflow.sh

# Run ggsearch
cp ../SAQE/15ggsearch.sh ./;
chmod +x ./15ggsearch.sh;
./15ggsearch.sh Trinity.fasta.transdecoder.pep Danio_rerio.GRCz11.pep.all.fa2;

# Make annotation table
cp ../SAQE/15parseggsearch.sh ./;
cp ../SAQE/15parseggsearch.pl ./;
# add target=Trinity.fasta.transdecoder.pep into in15parseggsearch.sh
chmod +x 15parseggsearch.sh;
./15parseggsearch.sh Danio_rerio.GRCz11.pep.all.fa2;
perl -nle 'print $1 if(/^\>(\S+)/)' Trinity.fasta.transdecoder.pep > coding-transcript.pid.txt;
cp ../SAQE/15mkannotbl.pl ./;
cat coding-transcript.pid.txt | perl 15mkannotbl.pl Trinity.fasta.transdecoder.pep-Danio_rerio.GRCz11.pep.all.fa2.txt > Trinity.fasta.transdecoder.pep_all.txt;

# Search downrregulated ONratio
awk -F',' '{ if ($NF ~ /downregulated/) { print } }' ctrl_ko_fltrexpr_contig_tpm_onratio.csv > downregulated_ONratio.csv;
awk -F',' '{ print $1 }' downregulated_ONratio.csv > downregulated_ONratio_list.txt;
wc -l downregulated_ONratio_list.txt
# downregulated_ONratio = 439908

# search dTPM (downrregulated ONratio & off-target)
# Make Venn diagram using https://bioinformatics.psb.ugent.be/webtools/Venn/
###
List names	number of elements	number of unique elements
1（off-target）	865452	865452
2（downrregulated ONratio）	439908	439908
Overall number of unique elements	935252
###
all_offtarget_vs_dONrartio.svg
all_offtarget_vs_dONrartio.png
all_offtarget_vs_dONrartio.txt
# downrregulated ONratio & off-target n = 370,108
# Saved it to all_offtarget_vs_dONrartio.id.txt
# Modify all_offtarget_vs_dONrartio.id.txt
cat all_offtarget_vs_dONrartio.txt | sed -n "2,370109p" > all_offtarget_vs_dONrartio.id.txt;

# Rnrichment analysis in group of downrregulated ONratio & off-target
while read line; do grep $line coding-transcript.pid.txt >> all_offtarget_vs_dONrartio.pid.txt;done < all_offtarget_vs_dONrartio.id.txt;
wc -l all_offtarget_vs_dONrartio.pid.txt # ORFありtranscript 2500
cat all_offtarget_vs_dONrartio.pid.txt | perl 15mkannotbl.pl Trinity.fasta.transdecoder.pep-Danio_rerio.GRCz11.pep.all.fa2.txt > all_offtarget_vs_dONrartio_annotation_table.txt;
awk '{ if ($2 ~ /ENSDARP/) { print } }' all_offtarget_vs_dONrartio_annotation_table.txt > all_offtarget_vs_dONrartio_annotation_table.txt_hasENSDARP.txt;
wc -l all_offtarget_vs_dONrartio_annotation_table.txt_hasENSDARP.txt;
# 1246 transcripts have ENSDARP ID

# search downrregulated ONratio for each mismatch number
for i in {0..11};do awk '{ print $2 }' offtarget_mm"$i".cas-offinder > offtarget_mm"$i"_list.txt;awk '!seen[$0]++' offtarget_mm"$i"_list.txt > offtarget_mm"$i"_uniq_list.txt;wc -l offtarget_mm"$i"_uniq_list.txt;done;
###
       0 offtarget_mm0_uniq_list.txt
       0 offtarget_mm1_uniq_list.txt
       0 offtarget_mm2_uniq_list.txt
       1 offtarget_mm3_uniq_list.txt
      25 offtarget_mm4_uniq_list.txt
     278 offtarget_mm5_uniq_list.txt
    2149 offtarget_mm6_uniq_list.txt
   13227 offtarget_mm7_uniq_list.txt
   60910 offtarget_mm8_uniq_list.txt
  198645 offtarget_mm9_uniq_list.txt
  475284 offtarget_mm10_uniq_list.txt
  786219 offtarget_mm11_uniq_list.txt
###

# Search downrregulated ONratio & off-target for each mismatch number
mkdir mm_offtarget_dONratio;
# Draw Venn diagram
conda activate matplotlib_venn_env;
for i in {0..11};do python ../drawVennDiagram.py mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.id.txt mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.tiff offtarget_mm"$i"_uniq_list.txt downregulated_ONratio_list.txt;wc -l mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.id.txt;done;
### ONratio & off-target for each mismatch number
       0 mm_offtarget_dONratio/mm0_offtarget_dONratio.id.txt
       0 mm_offtarget_dONratio/mm1_offtarget_dONratio.id.txt
       0 mm_offtarget_dONratio/mm2_offtarget_dONratio.id.txt
       1 mm_offtarget_dONratio/mm3_offtarget_dONratio.id.txt
       6 mm_offtarget_dONratio/mm4_offtarget_dONratio.id.txt
      62 mm_offtarget_dONratio/mm5_offtarget_dONratio.id.txt
     500 mm_offtarget_dONratio/mm6_offtarget_dONratio.id.txt
    3489 mm_offtarget_dONratio/mm7_offtarget_dONratio.id.txt
   18172 mm_offtarget_dONratio/mm8_offtarget_dONratio.id.txt
   67658 mm_offtarget_dONratio/mm9_offtarget_dONratio.id.txt
  182774 mm_offtarget_dONratio/mm10_offtarget_dONratio.id.txt
  330733 mm_offtarget_dONratio/mm11_offtarget_dONratio.id.txt
###
conda deactivate;

# Search transcript annotated with ENSDARP
for i in {0..11};do while read line; do grep $line coding-transcript.pid.txt >> mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.pid.txt;done < mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.id.txt;wc -l mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.pid.txt;done;
### transcript with ORFs
wc: mm_offtarget_dONratio/mm0_offtarget_dONratio.pid.txt: open: No such file or directory
wc: mm_offtarget_dONratio/mm1_offtarget_dONratio.pid.txt: open: No such file or directory
wc: mm_offtarget_dONratio/mm2_offtarget_dONratio.pid.txt: open: No such file or directory
       0 mm_offtarget_dONratio/mm3_offtarget_dONratio.pid.txt
       0 mm_offtarget_dONratio/mm4_offtarget_dONratio.pid.txt
       1 mm_offtarget_dONratio/mm5_offtarget_dONratio.pid.txt
      19 mm_offtarget_dONratio/mm6_offtarget_dONratio.pid.txt
      85 mm_offtarget_dONratio/mm7_offtarget_dONratio.pid.txt
     412 mm_offtarget_dONratio/mm8_offtarget_dONratio.pid.txt
    1170 mm_offtarget_dONratio/mm9_offtarget_dONratio.pid.txt
    2055 mm_offtarget_dONratio/mm10_offtarget_dONratio.pid.txt
    2469 mm_offtarget_dONratio/mm11_offtarget_dONratio.pid.txt
###
for i in {3..11};do cat mm_offtarget_dONratio/mm"$i"_offtarget_dONratio.pid.txt | perl 15mkannotbl.pl Trinity.fasta.transdecoder.pep-Danio_rerio.GRCz11.pep.all.fa2.txt > mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt;done;
for i in {3..11};do awk '{ if ($2 ~ /ENSDARP/) { print } }' mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt > mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt;wc -l mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt;done;
### transcript with ENSDARP ID
       0 mm_offtarget_dONratio/mm3_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
       0 mm_offtarget_dONratio/mm4_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
       1 mm_offtarget_dONratio/mm5_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
      11 mm_offtarget_dONratio/mm6_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
      50 mm_offtarget_dONratio/mm7_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
     212 mm_offtarget_dONratio/mm8_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
     567 mm_offtarget_dONratio/mm9_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
    1032 mm_offtarget_dONratio/mm10_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
    1228 mm_offtarget_dONratio/mm11_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt
###

### Calculate DENGERindex
# Get gene list for each mismatch number
for i in {3..11};do awk '{ print $3}' mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_annotation_table.txt_hasENSDARP.txt > mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_ENSDARP_genelist.txt;done;
# GO annotation with topGO for each mismatch number
conda activate topGO;
for i in {1..11};do Rscript ../annotateGOv2.R mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_ENSDARP_genelist.txt mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_ENSDARP_geneGO zebrafish;done;
# Enrichment analysis for each mismatch number
for i in {1..11};do Rscript ../makeDANGERenrichmentTablev2.R mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_ENSDARP_geneGO mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_ENSDARP_geneGO_enrichment mm_offtarget_dONratio/mm"$i"_offtarget_dONratio_ENSDARP_genelist.txt $i;done;
conda deactivate;
# DENGERindex Calculation
conda activate calcDANGERindex_env;
python ../calcDANGERindex.py  mm_offtarget_dONratio/mm6_offtarget_dONratio_ENSDARP_geneGO_enrichment mm_offtarget_dONratio/mm7_offtarget_dONratio_ENSDARP_geneGO_enrichment mm_offtarget_dONratio/mm8_offtarget_dONratio_ENSDARP_geneGO_enrichment mm_offtarget_dONratio/mm9_offtarget_dONratio_ENSDARP_geneGO_enrichment mm_offtarget_dONratio/mm10_offtarget_dONratio_ENSDARP_geneGO_enrichment mm_offtarget_dONratio/mm11_offtarget_dONratio_ENSDARP_geneGO_enrichment;
###
Biological Process
Sum of DANGER indexes
224.24425603474646
Upper boundary
0.14470317255959442
Cellular Component
Sum of DANGER indexes
68.85219343072643
Upper boundary
0.16767265997982977
Molecular Function
Sum of DANGER indexes
60.94544580544128
Upper boundary
0.13143181335175058
###
conda deactivate;
```

The result is saved in DANGER_analysis_result.

</details>