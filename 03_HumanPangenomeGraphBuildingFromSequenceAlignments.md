# Human Pangenome Bring Your Own Data (BYOD) analysis Workshop 2024

## Human pangenome graph building from sequence alignments

### Learning objectives

- collect and preprocess *de novo* human assemblies
- partition the assembly contigs by chromosome
- built chromosome-specific pangenome graphs

### Getting started

Ask for interactive session (let's ask for a bit more CPUs this round):

    srun --nodes=1 -c32 --mem=250g --time 24:00:00 --job-name "interactive_small" --pty /bin/bash

Make sure you have `pggb`, its tools loaded and `bgzip`:

    module load pggb
    module load htslib

Create a directory to work on for this tutorial:

    cd /cbio/projects/037/$USER
	mkdir human_pangenome_graphs
	cd human_pangenome_graphs

Download 2 human references and 4 diploid human *de novo* assemblies from the Human Pangenome Reference Consortium (HPRC) data:

    DIR_BASE=/cbio/projects/037/$USER
    mkdir -p $DIR_BASE/human_pangenome_graphs/assemblies
    cd $DIR_BASE/human_pangenome_graphs/assemblies

    wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.maternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.maternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.maternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.paternal.f1_assembly_v2_genbank.fa.gz
    wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.maternal.f1_assembly_v2_genbank.fa.gz

Decompress and index the assemblies:

    DIR_BASE=/cbio/projects/037/$USER
    cd $DIR_BASE/human_pangenome_graphs/assemblies
    gunzip *genbank.fa.gz
    ls *genbank.fa | while read f; do echo $f; samtools faidx $f; done

When making pangenome graphs, you should always include at least one reference genome in order to use it as a coordinate system (for projecting variants and/or exploiting its annotations).
We usually use both GRCh38 and CHM13.

### Pangenome Sequence Naming

We follow the [PanSN-spec](https://github.com/pangenome/PanSN-spec) naming to simplify the identification of samples and haplotypes in pangenomes.
The HPRC samples already follow such a convention (`1` is the PATERNAL haplotype, `2` is the MATERNAL haplotype), but we need to apply PanSN naming to the reference genomes.
So, let's add a prefix to their sequence names.
We can do that by using [fastix](https://github.com/ekg/fastix):

    DIR_BASE=/cbio/projects/037/$USER
    cd $DIR_BASE/human_pangenome_graphs/assemblies

    fastix -p 'grch38#1#' <(zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz) | bgzip -@ 16 > grch38_full.fa.gz
    samtools faidx grch38_full.fa.gz
    
    fastix -p 'chm13#1#' <(zcat chm13v2.0.fa.gz) | bgzip -@ 16 > chm13.fa.gz
    samtools faidx chm13.fa.gz

About GRCh38, we remove the unplaced contigs that are (hopefully) properly represented in CHM13:

    DIR_BASE=/cbio/projects/037/$USER
    cd $DIR_BASE/human_pangenome_graphs/assemblies
    samtools faidx grch38_full.fa.gz $(cut -f 1 grch38_full.fa.gz.fai | grep -v _ ) | bgzip -@ 16 > grch38.fa.gz

Cleaning:

    DIR_BASE=/cbio/projects/037/$USER
    cd $DIR_BASE/human_pangenome_graphs/assemblies
    rm chm13v2.0.fa.gz GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz grch38_full.fa.gz.*

Take a look at how sequence names are changed in the FASTA files.

### Sequence partitioning

To reduce analysis complexity, we partition assembly contigs by chromosome and generate chromosome-specific pangenome graphs.
For doing that, we first need to put the two reference genomes together

    DIR_BASE=/cbio/projects/037/$USER
    cd $DIR_BASE/human_pangenome_graphs/assemblies
    zcat chm13.fa.gz grch38.fa.gz | bgzip -@ 16 > chm13+grch38.fa.gz && samtools faidx chm13+grch38.fa.gz

and then map each assembly against the two reference genomes.
**For running this, open a new terminal (do not close the current one), login again and then run:**

    module load pggb
    module load htslib

    DIR_BASE=/cbio/projects/037/$USER
    cd $DIR_BASE/human_pangenome_graphs

    PATH_REFERENCES_FASTA=$DIR_BASE/human_pangenome_graphs/assemblies/chm13+grch38.fa.gz

    mkdir -p $DIR_BASE/human_pangenome_graphs/assemblies/partitioning

    ls $DIR_BASE/human_pangenome_graphs/assemblies/*.f1_assembly_v2_genbank.fa | while read FASTA; do
      NAME=$(basename $FASTA .fa);
      echo $NAME

      PATH_PAF=$DIR_BASE/human_pangenome_graphs/assemblies/partitioning/$NAME.vs.ref.paf
      sbatch -c32 -p Main --wrap "wfmash $PATH_REFERENCES_FASTA $FASTA -m -N -t 32 > $PATH_PAF"
    done

`wfmash` should take 4-5 minutes for each haplotype (each FASTA file).

Why are we using two reference genomes?

<details>
  <summary>Click me for the answer</summary>

A single human genome cannot fully represent the genetic variability of the entire population.
Having multiple reference genomes allows greater genetic variability to be represented, allowing contigs to be better partitioned.
</details>

Recently, a new high-quality human diploid assembly was released at [HG002](https://github.com/marbl/HG002).
You could also try using these 2 new haplotypes (`HG002#MATERNAL` and `HG002#PATERNAL`) to partition assembly contigs and see if you can partition more
(_we haven't seriously tried it yet, so we're very curious how much this new diploid assembly can help with the chromosome partitioning_).

Run `wfmash` without parameters to get information on the meaning of each parameter.

What does `-m` mean?

<details>
  <summary>Click me for the answer</summary>

We ask `wfmash` to compute only the mapping, not the alignment, to save computation time.
To partition chromosomes, we don't need the base-level alignments.
</details>

What does `-N` mean?

<details>
  <summary>Click me for the answer</summary>

With `-N` we ask `wfmash` to generate mappings that completely cover the assembly contigs.
Without `-N`, mappings can only cover parts of the contigs, and different parts of the same contig could map to difference reference chromosomes.
This is frequent for contigs belonging to acrocentric chromosomes or sex chromosomes because of their homologous regions (https://www.nature.com/articles/s41586-023-05976-y, https://doi.org/10.1093/hmg/7.13.1991).
</details>

For each haplotype (for each FASTA file), count how many contigs were partitioned in each reference chromosome.
Which of the two references (GRCh38 and CHM13) do more contigs map to? If there is a clear winner, why?

<details>
  <summary>Click me for the answer</summary>

CHM13 is a complete human genome (GRCh38 is almost complete, but not 100%), so we expect to be able to map more assembly contigs to it.
</details>

What chromosome do most contigs map to in GRCh38 and CHM13? and which chromosome has the least number of contigs mapped to GRCh38 and CHM13? Why?

<details>
  <summary>Click me for the answer</summary>

In GRCh38 the shor arms of the acrocentric chromosomes (chromosome 13, 14, 15, 21 and 22) are missing, so we struggle to align acrocentric contigs to it.
These short arms are available in CHM13 and indeed we are able to map lots of contigs against them.
</details>

It should be noted that also with `wfmash -N`, there can be cases with contigs fully mapping to different cromosomes.
For example:

    DIR_BASE=/cbio/projects/037/$USER
    cd $DIR_BASE/human_pangenome_graphs
    grep 'HG00438#2#JAHBCA010000147.1' $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/*paf

HG00438.maternal.f1_assembly_v2_genbank.vs.ref.paf:HG00438#2#JAHBCA010000147.1  738336  0       738336  +       chm13#1#chr13 113566686       8582789 9321125 245     738336  23      id:f:99.4899    kc:f:0.057824
HG00438.maternal.f1_assembly_v2_genbank.vs.ref.paf:HG00438#2#JAHBCA010000147.1  738336  0       738336  +       chm13#1#chr22 51324926        5299094 6037430 245     738336  23      id:f:99.4899    kc:f:0.057824

For which there is not enough information to determine which is the best chromosome to map against (_acrocentric chromosomes are hard!_).
For these case, we just randomly take one result (_we are working on implementing the random sampling directly in `wfmash`_, to make user life easier).

    DIR_BASE=/cbio/projects/037/$USER
    ls $DIR_BASE/human_pangenome_graphs/assemblies/*.f1_assembly_v2_genbank.fa | while read FASTA; do
      NAME=$(basename $FASTA .fa);
      echo $NAME

      PATH_PAF=$DIR_BASE/human_pangenome_graphs/assemblies/partitioning/$NAME.vs.ref.paf
      cut -f 1,6 $PATH_PAF | sed -e 's/chm13#1#//g' -e 's/grch38#1#//g' | awk '{
          contig = $1;
      
          # If the contig is not already in the data array, add it
          if (!(contig in data)) {
              data[contig] = $2;
          }
      }
      END {
        # Output the result
        for (contig in data) {
          print contig "\t" data[contig];
        }
      }' > $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/$NAME.vs.ref.assignments.tsv
    done

Now we can subset assembly contigs by chromosome:

    DIR_BASE=/cbio/projects/037/$USER
    ( seq 1 22; echo X; echo Y ) | while read i; do
      echo chr$i
      
      grep chr$i --no-filename -w $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/*.assignments.tsv | cut -f 1 > $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr$i.contigs.txt
    done

Then, we create a FASTA file for each chromosome, reference chromosomes included.
To save time and space, let's take only sequences from chromosome 20.
**For running this, return to the first terminal you've opened:**

    DIR_BASE=/cbio/projects/037/$USER
    ( echo 20 ) | while read i; do
      echo chr$i
      samtools faidx $DIR_BASE/human_pangenome_graphs/assemblies/chm13+grch38.fa.gz chm13#1#chr$i grch38#1#chr$i > $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr$i.fa
      
      ls $DIR_BASE/human_pangenome_graphs/assemblies/*.f1_assembly_v2_genbank.fa | while read FASTA; do
        NAME=$(basename $FASTA .fa);
        echo $NAME

        samtools faidx $FASTA $( comm -12 <(cut -f 1 $FASTA.fai | sort) <(sort $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr$i.contigs.txt) ) >> $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr$i.fa
      done

      bgzip -@ 16 $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr$i.fa
      samtools faidx $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr$i.fa.gz
    done

Check that everything went fine:

    DIR_BASE=/cbio/projects/037/$USER
    head $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr20.fa.gz.fai | column -t

### Building chromosome-specific pangenome graphs

Build the pangenome graph for chromosome 20.
**For running this, go to the second terminal you've opened:**

    module load pggb

    DIR_BASE=/cbio/projects/037/$USER
    mkdir -p $DIR_BASE/human_pangenome_graphs/graphs
    cd $DIR_BASE/human_pangenome_graphs/graphs

    sbatch -c32 -p Main --wrap "pggb -i $DIR_BASE/human_pangenome_graphs/assemblies/partitioning/chr20.fa.gz -o $DIR_BASE/human_pangenome_graphs/graphs/pggb.chr20 -p 98 -s 10k -k 79 -V 'chm13,grch38' -D /dev/shm -t 32"

This should take approximately 1 hour and will generate a pangenome graph, several graph visualizations, and 2 variant sets called from the assemblies.

**IMPORTANT**: The `-D` parameter in `pggb` is used to specify the directory used for temporary files.
This directory should be on a high-speed disk (like an SSD or, as here, a [ramdisk](https://en.wikipedia.org/wiki/RAM_drive)) to avoid severe slowdowns during graph construction, sorting and graph layout generation.

<details>
  <summary>Click me for considerations about the `-n` parameter</summary>

Note that we are not specifying the number of haplotypes (`-n` parameter).
Indeed, `pggb` can automatically obtain this information thanks to the fact that the input sequences respect the PanSN naming.
</details>

<details>
  <summary>Click me for considerations about the `-p` and `-s` parameters</summary>

Since human presents a low sequence divergence, we set the mapping identity (`-p` parameter) in `pggb` to `98`.
Additionally, we specify `-s 10k` (equivalent to specifying `-s 10000`) to get a simpler and more linear graph structure, which is easier to work with.
Lower values of `-s` and `-p` lead to more sensitive mappings but to the possibility of having circular graphs due to the sequence similarity of the telomeres.
</details>

<details>
  <summary>Click me for considerations about the `-k` parameter</summary>

The `-k` parameter is used to filter exact matches in the sequence alignments shorter than 79 bps.
Indeed, graph induction with `seqwish` often works better when we filter short matches out of the input alignments.
In practice, these often occur in regions of low alignment quality, which are typical of areas with large indels and structural variations.
Removing short matches can simplify the graph and remove spurious relationships caused by short repeated homologies.
However, this filter might lead to under-alignment, that is resolved in the graph normalization step with `smoothxg`.
By default, `pggb` uses `-k 23`.
With human data, higher values work well.
</details>

<details>
  <summary>Click me for considerations about the `-V` parameter</summary>

The `-V 'chm13,grch38'` parameter specifies to call variants from the assemblies using two different genomes as reference.
This will generate two VCF files: in one the variants will be expressed relative to CHM13, in the other the variants will be expressed relative to GRCh38.
</details>

Use `odgi stats` to obtain the graph length, and the number of nodes, edges, and paths.
Do you think the resulting pangenome graph represents the input sequences well?
Check the length and the number of the input sequences to answer this question.

<details>
  <summary>Click me for the answer</summary>

chr20 is approximately 64-66 Mbp long in reference genomes.
There is no right graph length because it always depends on how much genetic variability we want to represent in the graph.
Lower `-s` and `-p` values will tend to represent more homologies, leading to shorter, but more complex, graphs.
</details>

Take a look at the PNG files in the `pggb.chr20` folder.
Is the layout of the graph roughly linear?
Please note that when pangenome graphs are linear, the `*.draw multiqc.png` files may be difficult to view with the default image viewer.
In these cases, you can view them with alternative software such as [GIMP](https://www.gimp.org/).

Generate another `odgi viz` visualization with

    cd $DIR_BASE/graphs/pggb.chr20
    odgi paths -i chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og -L | cut -f 1,2 -d '#' | uniq > prefixes.txt
    odgi viz -i chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og -o chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og.viz_multiqc.2.png -x 1500 -y 500 -a 10 -M prefixes.txt

What do you think is different between the `chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og.viz_multiqc.png` image
and the newly generated image (`chr20.fa.gz.a8a102b.c2fac19.afc7f52.smooth.final.og.viz_multiqc.2.png`)?

<details>
  <summary>Click me for the answer</summary>

Contigs belonging to the same haplotype are visually collapsed into one colored bar each.
</details>
