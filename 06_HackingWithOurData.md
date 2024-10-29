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
https://github.com/AndreaGuarracino/ReferenceGraphPangenomeDataAnalysisHackathon2023/blob/master/06_HackingWithOurData.md
