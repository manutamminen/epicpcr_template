#########
# Prepare databases
#########

rule download_tax_db:
  output:
    temp("SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz")
  shell:
    "wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"


rule download_aln_db:
  output:
    "data/external/aln-dbs/SILVA_132_SSURef_NR99_13_12_17_opt.arb",
    "data/external/aln-dbs/SSURef_NR99_128_SILVA_07_09_16_opt.arb"
  shell:
    """
    mkdir -p data/external/aln-dbs &&\
    cd data/external/aln-dbs &&\
    wget https://www.arb-silva.de/fileadmin/arb_web_db/release_132/ARB_files/SILVA_132_SSURef_NR99_13_12_17_opt.arb.gz &&\
    wget https://www.arb-silva.de/fileadmin/arb_web_db/release_128/ARB_files/SSURef_NR99_128_SILVA_07_09_16_opt.arb.gz &&\
    gunzip *.gz
    """


rule prepare_tax_dbs:
  input:
    "data/external/aln-dbs/SILVA_132_SSURef_NR99_13_12_17_opt.arb",
    "data/external/aln-dbs/SSURef_NR99_128_SILVA_07_09_16_opt.arb",
    "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
  output:
    "data/external/tax-dbs/16S_db.fasta",
    "data/external/tax-dbs/18S_db.fasta"
  script:
    "src/data/prepare_taxonomy_dbs.py"



#########
# Prepare sequences
#########

IDs, = glob_wildcards("data/raw/{id}_L001_R1_001.fastq.gz")


rule read_merging:
  input:
    read1="data/raw/{id}_L001_R1_001.fastq.gz",
    read2="data/raw/{id}_L001_R2_001.fastq.gz"
  output:
    merged=temp("data/fastq/{id}_merged.fastq")
  conda:
    "envs/vsearch.yaml"
  shell:
    """
    vsearch --fastq_mergepairs {input.read1} \
        --reverse {input.read2} \
        --fastq_minovlen 50 \
        --fastq_maxdiffs 15 \
        --fastqout {output.merged}
    """

rule quality_filtering:
  input:
    merged="data/fastq/{id}_merged.fastq"
  output:
    filtered=temp("data/fasta/{id}_filtered.fasta")
  conda:
    "envs/vsearch.yaml"
  shell:
    """
    vsearch --fastq_filter {input.merged} \
        --fastq_maxee 1 \
        --fastaout {output.filtered} \
    """

rule combine_fastas:
  input:
    seqs=expand(rules.quality_filtering.output.filtered, id=IDs)
  output:
    "data/combined/16S_seqs.fasta",
    "data/combined/16S_bcs.csv",
    "data/combined/18S_seqs.fasta",
    "data/combined/18S_bcs.csv"
  script:
    "src/data/make_dataset.py"


rule cluster:
  input:
    bact="data/combined/16S_seqs.fasta",
    euk="data/combined/18S_seqs.fasta"
  output:
    bact_clust=protected("data/clustered/16S_seqs_centroids.fasta"),
    bact_uc=protected("data/clustered/16S_seqs.uc"),
    euk_clust=protected("data/clustered/18S_seqs_centroids.fasta"),
    euk_uc=protected("data/clustered/18S_seqs.uc")
  conda:
    "envs/vsearch.yaml"
  shell:
    """
    vsearch --cluster_size {input.bact} \
    --id 0.95 \
    --strand plus \
    --sizeout \
    --uc {output.bact_uc} \
    --centroids {output.bact_clust} &&\
    vsearch --cluster_size {input.euk} \
    --id 0.95 \
    --strand plus \
    --sizeout \
    --uc {output.euk_uc} \
    --centroids {output.euk_clust}
    """


#########
# Prepare taxonomies
#########

rule taxonomy:
  input:
    bact="data/clustered/16S_seqs_centroids.fasta",
    euk="data/clustered/18S_seqs_centroids.fasta",
    bact_db="data/external/tax-dbs/16S_db.fasta",
    euk_db="data/external/tax-dbs/18S_db.fasta"
  output:
    bact_clust="data/taxonomy/16S_seqs_centroids_tax.txt",
    euk_clust="data/taxonomy/18S_seqs_centroids_tax.txt"
  conda:
    "envs/vsearch.yaml"
  shell:
    """
    vsearch \
     --sintax {input.bact} \
     --db {input.bact_db} \
     --tabbedout {output.bact_clust} \
     --strand both &&\
    vsearch \
     --sintax {input.euk} \
     --db {input.euk_db} \
     --tabbedout {output.euk_clust} \
     --strand both
    """


rule collapse_taxonomy:
  input:
    "data/taxonomy/16S_seqs_centroids_tax.txt",
    "data/clustered/16S_seqs_centroids.fasta",
    "data/taxonomy/18S_seqs_centroids_tax.txt",
    "data/clustered/18S_seqs_centroids.fasta"
  output:
    "data/taxonomy/16S_otu_sequences.txt",
    "data/taxonomy/18S_otu_sequences.txt"
  script:
    "src/data/collapse_taxonomy.py"


rule annotate_bcs:
  input:
    "data/clustered/16S_seqs.uc",
    "data/taxonomy/16S_seqs_centroids_tax.txt",
    "data/combined/16S_bcs.csv",
    "data/clustered/18S_seqs.uc",
    "data/taxonomy/18S_seqs_centroids_tax.txt",
    "data/combined/18S_bcs.csv"
  output:
    "data/taxonomy/16S_bc_tax.txt",
    "data/taxonomy/18S_bc_tax.txt"
  script:
    "src/data/annotate_bcs.py"


rule complete_taxonomy:
  input:
    "data/taxonomy/16S_bc_tax.txt",
    "data/taxonomy/18S_bc_tax.txt",
    "data/taxonomy/16S_otu_sequences.txt",
    "data/taxonomy/18S_otu_sequences.txt"


#########
# Prepare alignments
#########

rule add_outgroups:
  input:
    "data/clustered/16S_seqs_centroids.fasta",
    "data/clustered/18S_seqs_centroids.fasta"
  output:
    "data/alignment/16S_seqs_outgroup.fasta",
    "data/alignment/18S_seqs_outgroup.fasta"
  script:
    "src/data/add_outgroups.py"


rule align_outgrouped_sequences:
  input:
    bact_outgrouped = "data/alignment/16S_seqs_outgroup.fasta",
    euk_outgrouped = "data/alignment/18S_seqs_outgroup.fasta",
    bact_db = "data/external/aln-dbs/SSURef_NR99_128_SILVA_07_09_16_opt.arb",
    euk_db = "data/external/aln-dbs/SILVA_132_SSURef_NR99_13_12_17_opt.arb"
  output:
    bact_aligned = "data/alignment/16S_seqs_outgroup_align.fasta",
    euk_aligned = "data/alignment/18S_seqs_outgroup_align.fasta"
  conda:
    "envs/sina.yaml"
  shell:
    """
    sina -i {input.bact_outgrouped} \
    --intype fasta \
    -o {output.bact_aligned} \
    --outtype fasta \
    --db {input.bact_db} &&\
    sina -i {input.euk_outgrouped} \
    --intype fasta \
    -o {output.euk_aligned} \
    --outtype fasta \
    --db {input.euk_db}
    """


rule infer_phylogenies:
  input:
    bact_aligned = "data/alignment/16S_seqs_outgroup_align.fasta",
    euk_aligned = "data/alignment/18S_seqs_outgroup_align.fasta"
  output:
    bact_tree = "data/trees/16S.tre",
    euk_tree = "data/trees/18S.tre"
  conda:
    "envs/fasttree.yaml"
  shell:
    """
    fasttree -nt {input.bact_aligned} > {output.bact_tree} &&\
    fasttree -nt {input.euk_aligned} > {output.euk_tree}
    """
