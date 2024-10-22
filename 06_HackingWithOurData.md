# Human Pangenome Bring Your Own Data (BYOD) analysis Workshop 2024

## 2024: Understanding H3ABioNet pangenome graphs

### Chromosome 6 pangenome graph with XX haplotypes

You can find `XX` whole-genome assemblyes here: `/cbio/projects/XXXX`:

Let's align the assemblies (our pangenome) against the `CHM13` human reference genome. We use `WFMASH` for the alignment:

```shell
ls /lizardfs/guarracino/pangenomes/washu_pedigree/*.fa.gz | while read FASTA; do
    SAMPLE=$(basename $FASTA .fa.gz)
    sbatch -c 48 -p allnodes --job-name wfmash-$SAMPLE --wrap "wfmash /lizardfs/guarracino/robertsonian_translocation/assemblies/chm13v2.0.fa.gz $FASTA -p 95 -s 10k -t 48 > /lizardfs/guarracino/h3africa/$SAMPLE-vs-chm13v2.aln.paf"
done
```

With `IMPG` we can then project any reference genome onto the assemblies. We will projec the MHC region:

```shell
# Merge all alignments together for convenience
cat /lizardfs/guarracino/h3africa/*-vs-chm13v2.aln.paf > /lizardfs/guarracino/h3africa/pangenome-vs-chm13v2.aln.paf

impg -p /lizardfs/guarracino/h3africa/pangenome-vs-chm13v2.aln.paf -r chr6:29595119-32911317 | bedtools sort | bedtools merge -d 300000 > /lizardfs/guarracino/h3africa/chm13v2+pangenome.MHC.bed
```

Let's extract the MHC from the pangenome (and `CHM13`):

```shell
# Take MHC on the reference and apply PanSN
samtools faidx /lizardfs/guarracino/robertsonian_translocation/assemblies/chm13v2.0.fa.gz chr6:29595119-32911317 | sed 's/>chr/>chm13#1#chr/g' > /lizardfs/guarracino/h3africa/MHC.ref+pan.fa
ls /lizardfs/guarracino/pangenomes/washu_pedigree/*.fa.gz | while read FASTA; do
    samtools faidx $FASTA -r <( grep -Ff <(cut -f 1 $FASTA.fai) chm13v2+pangenome.MHC.bed -w | awk -v OFS='\t' '{print($1":"$2+1"-"$3)}' )
done >> /lizardfs/guarracino/h3africa/MHC.ref+pan.fa
samtools faidx /lizardfs/guarracino/h3africa/MHC.ref+pan.fa
```

Time for [building pangenome graphs](https://doi.org/10.1038/s41592-024-02430-3):

```shell
pggb -i /lizardfs/guarracino/h3africa/MHC.ref+pan.fa -o /lizardfs/guarracino/h3africa/pggb.MHC.ref+pan -D /scratch
```

Let's take a look at 1D and 2D visualization:


![chr6 MHC locus 1D](images/MHC.ref+pan.fa.bf3285f.11fba48.7b43761.smooth.final.og.viz_depth_multiqc.png)

![chr6 MHC locus 2D](images/MHC.ref+pan.fa.bf3285f.11fba48.7b43761.smooth.final.og.lay.draw.png)


Let's try the `IMPG` way with the C4 region:

```shell
impg -p /lizardfs/guarracino/h3africa/pangenome-vs-chm13v2.aln.paf -r chr6:31823879-31909831 | bedtools sort | bedtools merge -d 300000 > /lizardfs/guarracino/h3africa/chm13v2+pangenome.C4.bed

samtools faidx /lizardfs/guarracino/robertsonian_translocation/assemblies/chm13v2.0.fa.gz chr6:31823879-31909831 | sed 's/>chr/>chm13#1#chr/g' > /lizardfs/guarracino/h3africa/C4.ref+pan.fa
ls /lizardfs/guarracino/pangenomes/washu_pedigree/*.fa.gz | while read FASTA; do
    samtools faidx $FASTA -r <( grep -Ff <(cut -f 1 $FASTA.fai) chm13v2+pangenome.C4.bed -w | awk -v OFS='\t' '{print($1":"$2+1"-"$3)}' )
done >> /lizardfs/guarracino/h3africa/C4.ref+pan.fa
samtools faidx /lizardfs/guarracino/h3africa/C4.ref+pan.fa

pggb -i /lizardfs/guarracino/h3africa/C4.ref+pan.fa -o /lizardfs/guarracino/h3africa/pggb.C4.ref+pan -D /scratch
```

![chr6 C4 locus 1D](images/C4.ref+pan.fa.bf3285f.11fba48.7b43761.smooth.final.og.viz_depth_multiqc.png)

![chr6 C4 locus 2D](images/C4.ref+pan.fa.bf3285f.11fba48.7b43761.smooth.final.og.lay.draw.png)

Our [grafiocavallo](https://en.wikipedia.org/wiki/Caciocavallo) is a bit dirty because of short matches that lead to new edges in the graph.

## 2023: Understanding H3ABioNet pangenome graphs

### Chromosome 6 pangenome graph with 14 haplotypes

You can find the `pggb` graphs for all chromosomes here: `/cbio/projects/031/andreaguarracino/6samples`.
Each graph was made with 6 diploid samples (so 12 haplotypes) plus 2 reference genomes (`GRCh38` and `CHM13`).

We extract the MHC locus and color bars by haplotype.
 
![chr6 MHC locus](images/chr6.pan.MHC.png)

The MHC locus is a bit broken:
- the worst sample, `CMI0Q`, covers the locus with multiple contigs;
- the best sample, `FPV0R`, covers the locus with 1 contig for each haplotype.

Let's take a look at the C4 locus.
In 1D, we color by copy number:
- white: 0 copies
- grey: 1 copy
- red: 2 copies
- orange: 3 copies

![chr6 C4 locus 1D](images/chr6.pan.C4.sorted.m.png)

The `EFB0A#1` haplotype has 3 copies of the C4 genes.

The graph layout is the [grafiocavallo](https://en.wikipedia.org/wiki/Caciocavallo) that we expect:

![chr6 C4 locus 2D](images/chr6.pan.C4.sorted.2D.png)

The big loop represents the copy number variation, while the nested loop differentiates the short and long forms of the C4 genes.

If we inject `GRCh38` annotations and untangle the graph, we can see a bit closer the C4 variation:

![chr6 C4 locus untangling](images/chr6.pan.C4.untangling.png)

The locus is inverted in the graph (not a problem).

### Same variant, different representations

The same variant can be expressed differently into a VCF file with respect to the reference sequence to which variation is expressed.

```shell
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CONFLICT,Number=.,Type=String,Description="Sample names for which there are multiple paths in the graph with conflicting alleles">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">

##contig=<ID=chm13#1#chr20,length=66210255>
#CHROM          POS     ID              REF  ALT  QUAL  FILTER  INFO                                                                      FORMAT  AB2563  CMI0Q  CQD0R  EFB0A  FPM0J  FPV0R  grch38
chm13#1#chr20   553892  <338912<338910  CT   C    60    .       AC=9;AF=0.692308;AN=13;AT=>338910>338911>338912,>338910>338912;NS=7;LV=0  GT      0|1     1|1    0|0    1|1    1|1    1|0    1

##contig=<ID=grch38#1#chr20,length=64444167>
#CHROM          POS     ID              REF  ALT  QUAL  FILTER  INFO                                                                      FORMAT  AB2563  CMI0Q  CQD0R  EFB0A  FPM0J  FPV0R  chm13
grch38#1#chr20  510173  <338912<338910  C    CT   60    .       AC=5;AF=0.384615;AN=13;AT=>338910>338912,>338910>338911>338912;NS=7;LV=0  GT      1|0     0|0    1|1    0|0    0|0    0|1    1

#CHROM               POS     ID              REF  ALT  QUAL  FILTER  INFO                                                                      FORMAT  AB2563  CMI0Q  CQD0R  EFB0A  FPM0J  chm13  grch38
FPV0R#1#h1tg000005l  546533  <338912<338910  C    CT   60    .       AC=4;AF=0.333333;AN=12;AT=>338910>338912,>338910>338911>338912;NS=7;LV=0  GT      1|0     0|0    1|1    0|0    0|0    1      0
```

![chr20 bubble Bandage](images/6samples.littleBubble.Bandage.chr20.png)


### Chromosome 6 pangenome graph with 46 haplotypes

You can find the `pggb` graphs for all chromosomes here: `/cbio/projects/031/andreaguarracino/6samples`.
Each graph was made with 6 Baylor diploid samples (12 haplotypes) plus 16 HPRC (32 haplotypes) samples plus 2 reference genomes (`GRCh38` and `CHM13`).
We extract the MHC locus and color bars by haplotype.

![chr6 MHC locus](images/chr6.pan.MHC.46.png)
![chr6 MHC locus](images/chr6.pan.MHC.hack.png)

There are suspicious "indels". Let's try a local `pggb` run:

![chr6 MHC locus](images/chr6.pan.MHC.hack.2.png)
