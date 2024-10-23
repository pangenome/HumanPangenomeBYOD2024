# Pangenome-based haplotype deconvolution

## Preparation

Create a folder for the upcoming analysis:

```shell
cd /cbio/projects/037/$USER
mkdir -p haplotype_deconvolution
cd /cbio/projects/037/$USER/haplotype_deconvolution
```

Quicker! If at the workshop, just use a symlink to get the data (sequencing reads, reference, and pangenome for chr6) in context:

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
ln -s /cbio/projects/037/erikg/haplotype_deconvolution/sequencing_reads .
ln -s /cbio/projects/037/erikg/haplotype_deconvolution/reference .
ln -s /cbio/projects/037/erikg/haplotype_deconvolution/pangenome .
```

(and you can skip the downloads)

ALSO, we need to install gfainject and gafpack.

```shell
module load rust
cargo install --git https://github.com/pangenome/gafpack
cargo install --git https://github.com/chfi/gfainject
export PATH=~/.cargo/bin/:$PATH
```

## Sample reads

We are going to use pre-aligned reads against the human reference genome hg38 in CRAM format. We download 3 samples from 1000G.

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p sequencing_reads
cd sequencing_reads

# Get list of 1000G samples
wget https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index

# Get 3 samples
grep 'HG00438\|HG00673\|HG00735' 1000G_698_related_high_coverage.sequence.index | cut -f 1 > 1000G.selected.txt

# Add also links for CRAM indexes
sed 's/$/.crai/' 1000G.selected.txt >> 1000G.selected.txt

# Download CRAM and index files
wget -i 1000G.selected.txt

cd ..
```

Since we have sequencing reads in CRAM format, a reference genome is needed to fetch reads from the original alignments:

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p reference
cd reference

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

cd ..
```

## Pangenome assemblies

We will use year 1 assemblies from the Human Pangenome Reference Consortium. For simplicity, we will download the `PGGB` pangenome graph of chromosome 6 and get the contigs in FASTA format with `ODGI`:

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p pangenome
cd pangenome

# Download chromosome 6 pangenome graph
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2021_11_16_pggb_wgg.88/chroms/chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa.gz
gunzip chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa.gz

# GFA -> FASTA
module load pggb

odgi paths -i chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.gfa -f -t 16 -P > chr6.pan.fa
samtools faidx chr6.pan.fa

cd ..
```

## Region of interest pangenome

RESTART HERE.

We will genotype samples the C4 region. For that, we need to find the corresponding region on all the assemblies in the pangenome. As easy way is to first, we need to align the pangenome against the reference:

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p wfmash

wfmash \
    reference/GRCh38_full_analysis_set_plus_decoy_hla.fa pangenome/chr6.pan.fa \
    -s 10k -p 95 \
    -t 16 \
    > wfmash/pangenome_to_reference.paf
```

and then project the alignments of the C4 region by using `IMPG`:

```shell
# load impg module
module load impg

# Download BED file containing the C4 region coordinates on hg38's chromosome 6
wget https://raw.githubusercontent.com/pangenome/HumanPangenomeBYOD2024/refs/heads/main/data/hap-deconv.region-of-interest.bed

cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p impg

impg \
    -p wfmash/pangenome_to_reference.paf \
    -b hap-deconv.region-of-interest.bed \
    -x \
    > impg/projected.bedpe
```

Let's extract now the C4 regions in the pangenome:

```shell
# Merge regions
bedtools sort -i impg/projected.bedpe | \
    bedtools merge -d 100000 -i - \
    > impg/merged.bed

samtools faidx \
    -r <(awk -v OFS='\t' '{print($1":"$2+1"-"$3)}' impg/merged.bed) \
    pangenome/chr6.pan.fa \
    > impg/extracted.fasta
samtools faidx impg/extracted.fasta
```

Now we build a C4 region pangenome graph with `PGGB` (n.b. there are 36 haplotypes):

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p pggb

pggb -i impg/extracted.fasta -o pggb -t 16 -n 36

# Rename the final ODGI graph in a more human-friendly way
mv pggb/*smooth.final.og pggb/final.og
```

Let's chop the graph and get a haplotype coverage matrix:

```shell
odgi chop \
    -i pggb/final.og \
    -c 32 \
    -o odgi/chopped.og

odgi paths \
    -i odgi/chopped.og \
    -H | \
    cut -f 1,4- | \
    gzip > odgi/paths_matrix.tsv.gz
```

## Reads-vs-graph alignment

Let's index the C4 region pangenome and align sequencing reads to it with `BWA MEM`:

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p alignment

module load bwa
bwa index impg/extracted.fasta

ls sequencing_reads/*cram | while read CRAM; do
    echo $CRAM

    NAME=$(basename $CRAM .cram)
    
    # Extract reads covering the C4 region and then align them against the pangenome
    samtools view \
        -T reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        -L hap-deconv.region-of-interest.bed \
        -M \
        -b \
        $CRAM | \
        samtools fasta | \
            bwa mem -t 14 impg/extracted.fasta - | \
            samtools view -b -F 4 -@ 2 - \
            > alignment/$NAME.reads_vs_extracted.bam
done
```

Inject the pangenome alignments into the pangenome graph:

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p odgi

odgi view \
    -i odgi/chopped.og \
    -g > odgi/chopped.gfa

ls sequencing_reads/*cram | while read CRAM; do
    echo $CRAM

    NAME=$(basename $CRAM .cram)

    gfainject \
        --gfa odgi/chopped.gfa \
        --bam alignment/$NAME.reads_vs_extracted.bam \
        > alignment/$NAME.injected.gaf
done
```

## Matrixes

Let's get the sample coverage vectors:

```shell
ls sequencing_reads/*cram | while read CRAM; do
    echo $CRAM

    NAME=$(basename $CRAM .cram)

    gafpack \
        -g odgi/chopped.gfa \
        -a alignment/$NAME.injected.gaf \
        --len-scale | \
        gzip > alignment/$NAME.coverage.gafpack.gz
done
```

## Genotyping

**IMPORTANT**: currently compatible with `cosigt` v0.1.0 (commit `92622f6c095a1b43b0e13fef2893c74b9bfee887`).

```shell
cd /cbio/projects/037/$USER/haplotype_deconvolution
mkdir -p cosigt
module load cosigt

# First generate the similarity matrix using odgi
odgi similarity \
    -i odgi/chopped.og \
    > odgi/similarity.tsv

# Download the clustering script from cosigt repository
wget https://raw.githubusercontent.com/davidebolo1993/cosigt/16b18815cf9fdfcbf2afbf588a02740c27941ee3/cosigt_smk/workflow/scripts/cluster.r
chmod +x cluster.r

module load R

# First install required R packages in user's home directory
R --vanilla <<EOF
.libPaths(c("~/R/library", .libPaths()))
dir.create("~/R/library", recursive=TRUE, showWarnings=FALSE)
install.packages(c("reshape2", "NbClust", "rjson", "dendextend", "ggplot2", "data.table"), 
                lib="~/R/library",
                repos="https://cloud.r-project.org")
EOF

# Create a wrapper script that sets the library path before running cluster.r
cat > run_cluster.r <<'EOF'
#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
source("cluster.r")
EOF

# Run the clustering script to generate the JSON
./run_cluster.r odgi/similarity.tsv cosigt/clusters.json

ls sequencing_reads/*cram | while read CRAM; do
    echo $CRAM

    NAME=$(basename $CRAM .cram)

    cosigt \
    -i $NAME \
    -p odgi/paths_matrix.tsv.gz \
    -g alignment/$NAME.coverage.gafpack.gz \
    -o cosigt

    mv cosigt/cosigt_genotype.tsv cosigt/$NAME.cosigt_genotype.tsv
    mv cosigt/sorted_combos.tsv cosigt/$NAME.sorted_combos.tsv
done
```

Check the output:


```shell
grep '^alignment' cosigt/*.cosigt_genotype.tsv | column -t
```
