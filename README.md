**Analysis of shotgun metagenomics** 

Sixty Duroc pigs (30 immunocastrated males and 30 females) from the same genetic line were housed at the IRTA experimental farm. They were assigned to stress and control groups based on space allowance (stress = 1 m²/pig; control = 1.5 m²/pig). Additionally, pigs in the stress group underwent two mixing events during the growing-finishing period. Data availability in BioProject accession number PRJNA1208715

**1. Computational analysis of metagenomics**

The nf-core/mag 2.5.1 pipeline (https://nf-co.re/mag/2.5.1/) was employed to perform quality control and host decontamination by removing reads that mapped against the porcine Susscrofa11.1 assembly genome. Then, the remaining reads were processed following complementary approaches: 

**1.1 Taxonomic classification with sylph v0.8 against the GTDB-R220 database release 09-RS220, April 24th, 2024**

1.1.1 Sketch multiple samples

sylph sketch -1 /input_files/*R1.fastq.gz -2 /input_files/*R2.fastq.gz

1.1.2 Metagenomic profiling

sylph profile *.paired.sylsp gtdb-r220-c200-dbv1.syldb -t 30 -o results_profile.tsv

conda activate pyt37

python /sylph-utils/sylph_to_taxprof.py -s results_profile.tsv -m sylph-utils/gtdb_r220_metadata.tsv.gz -o profiles

python sylph-utils/merge_sylph_taxprof.py *.sylphmpa --column relative_abundance -o Taxprofiles.tsv



**1.2 Single sample assemblies**

Individual sample assemblies were done with the nf-core/mag 2.5.1 pipeline with Megahit v1.0.2. Binning was performed using MaxBin2 and MetaBAT2. The resulting bins were refined with the bin refinement module of metaWRAP and dereplicated with dRep v3.5.0. The parameters for the refinement module were set at a minimum completion of 70% (-c 70) and a maximum contamination of 10% (-x 10). 

REF="/Database/Sscrofa11.1"

OUTDIR="/MAGs"

INPUT="samples_list.csv"

MAX_CPUS=35


nextflow run nf-core/mag -r 2.5.1 \
  --input "${INPUT}" \
  -profile singularity \
  --outdir "${OUTDIR}" \
  --max_cpus ${MAX_CPUS} \
  --host_fasta "${REF}/Sscrofa11.1.dna.toplevel.fa" \
  --skip_spades \
  --skip_quast \
  --skip_concoct \
  --skip_prokka \
  --skip_prodigal \
  --skip_metaeuk \
  --skip_clipping \
  --skip_binqc \
  --binning_map_mode own


  **2. MAGs Taxonomy and functional annotation**
  
  Taxonomy classification of MAGs was done using GTDB-Tk v2.4.0, and their functional annotation was performed with DRAM against the PFAM-A, KOfam, and dbCAN-V10 databases. 

2.1 Taxonomy classification

conda activate GTDBTKv2.4.0
gtdbtk classify_wf --genome_dir dereplicated_bins/ --out_dir MAGs_gtdbk.dir -x *.fa --cpus 30

2.2 MAGs functional annotation

conda activate DRAM

DRAM.py annotate -i '/MAGs2used/*.fa' -o dram_annotationt --threads 30

**3. Microbial Genes Selection**

Scripts/Lasso_Selected_Microbial_Genes.R

This script applies a LASSO (Least Absolute Shrinkage and Selection Operator) regression model to perform feature selection on microbial annotated genes data. The goal is to identify the most relevant microbial genes associated with the phenotype of interest after correcting for confounding effects (e.g., sex). The script includes data preprocessing, standardization, adjustment for confounding factors on the binary trait via logistic regression, model training with cross-validation, and extraction of selected features based on the optimal regularization parameter.


 
  



