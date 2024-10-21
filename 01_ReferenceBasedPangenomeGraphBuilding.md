# Reference Graph Pangenome Data Analysis Hackathon 2023

## Reference-based pangenome graph building

### Learning objectives

- construct graphs using `vg construct`
- visualize graphs using `vg view` and `Bandage`
- convert graphs using `vg view`

### Getting started

Ask for interactive session

    srun --nodes=1 --tasks=2 --mem=4g --time 24:00:00 --job-name "interactive_small" --pty /bin/bash

Make sure you have `vg` loaded (it is inside `pggb` tools)

    module load pggb

In this exercise, you will use small toy examples from the `test` directory of `vg`.
So make sure you have checked out `vg` repository:

    cd /cbio/projects/037/$USER
	git clone https://github.com/vgteam/vg.git

Now create a directory to work on for this tutorial:

    cd /cbio/projects/037/$USER
	mkdir ref_based_pangenome_graph_building
	cd ref_based_pangenome_graph_building
	ln -s /cbio/projects/037/$USER/vg/test/tiny

### Constructing and viewing your first graphs

First we will use `vg construct` to build our first graph.
We will construct a reference-based graph from the sequence in file `tiny/tiny.fa`.

To construct a simple graph, run:

	vg construct -r tiny/tiny.fa -m 32 > tiny.ref.vg

This will create a graph that just consists of a linear chain of nodes, each with 32 characters.

The `-m` option tells `vg` to put at most 32 characters into each graph node.

To visualize a graph, you can use `vg view`.
By default, `vg view` will output a graph in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format.

	vg view tiny.ref.vg | column -t

In GFA format, each line is a separate record of some part of the graph.
The lines come in several types, which are indicated by the first character of the line.
What do you think they indicate?

<details>
  <summary>Click me for the answers</summary>

In a GFA file there are these following line types:
- `H`: a header.
- `S`: a "sequence" line, which is the sequence and ID of a node in the graph.
- `L`: a "link" line, which is an edge in the graph.
- `P`: a "path" line, which labels a path of interest in the graph. In this case, the path is the walk that the reference sequence takes through the graph.

Of note, the format does not specify that these lines come in a particular order.
</details>

Try to run `vg construct` with different values of `-m` and observe the different results.

Now let's build a new graph that has some variants built into it.
First, take a look at `tiny/tiny.vcf.gz`, which contains variants in (gzipped) [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format.

    zgrep '^##' -v tiny/tiny.vcf.gz | column -t

then use `vg construct` to build the graph:

	vg construct -r tiny/tiny.fa -v tiny/tiny.vcf.gz -m 32 > tiny.vg

We can write the graph in [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) format:

    vg view tiny.vg > tiny.gfa

A tool for visualizing (not too big) pangenome graphs is [BandageNG](https://github.com/asl/BandageNG).
It supports graphs in GFA format.
To use Bandage, just download it **locally on your computer**.

If you have Linux, run:

    wget -c https://github.com/asl/BandageNG/releases/download/v2022.09/BandageNG-9eb84c2-x86_64.AppImage
    chmod +x BandageNG-9eb84c2-x86_64.AppImage

If you have a Mac, run:

    wget -c https://github.com/asl/BandageNG/releases/download/v2022.09/BandageNG.dmg
    chmod +x BandageNG.dmg

Then download the graph on your computer by executing **locally on your computer**:

    USER="PUT HERE YOUR USER NAME ON ILIFU"
    scp $USER@slurm.ilifu.ac.za:/cbio/projects/037/$USER/ref_based_pangenome_graph_building/tiny.gfa ~/Desktop

and try to visualize it **locally on your computer**.

If you have a Linux, run:

    ./BandageNG-9eb84c2-x86_64.AppImage load ~/Desktop/tiny.gfa

If you have a Mac, run:

    ./BandageNG.dmg load ~/Desktop/tiny.gfa

Click the `Drawn graph` button to visualize the graph.
Click the `Length Node` button under Node labels feature.
Select the longest node (it is in the middle, 19 bp) by clicking on it, then click the `Paths...` button to see the paths that pass through the selected node.
Check that the paths displayed are the same as you expect by checking the GFA. Any problems?

<details>
  <summary>Click me for the answers</summary>

The GFA format does not specify that its lines have to follow a particular order.
However, the downloaded version of Bandage seems to have a bug in reading the paths properly.
Edit the GFA file by moving the `P` at the end of the file.
Now `Bandage` should show correctly the paths traversing the selected node.
Never trust the tools!
</details>

Try to generate a PNG image with `Bandage` (take a look at the `BandageNG-9eb84c2-x86_64.AppImage -h` / `BandageNG.dmg -h` output).
This can be helpful to take at graphs that are too big to be directly loaded with `Bandage`.

Now, let's build a new complex graph.
First, take a look at `tiny/multi.vcf.gz`, which contains multiallelic variants.
Try to generate with that input a graph in GFA format and visualize the output with `Bandage`. What are you able to see?
