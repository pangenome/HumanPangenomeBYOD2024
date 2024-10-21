# Human Pangenome Bring Your Own Data (BYOD) analysis Workshop 2024

## Short read mapping and variant calling

### Simulated short reads from Yeast

We have a `pggb` graph with 3 Yeast whole genomes and we align simulated short reads (250 bp long) from another Yeast whole genome:

```shell
vg autoindex --workflow giraffe -g $GFA -p $NAME.index --tmp-dir /scratch3/users/$USER/ --threads 32
vg giraffe --threads 32 -Z $NAME.index.giraffe.gbz -d $NAME.index.dist -m $NAME.index.min --sample $SAMPLE --read-group $SAMPLE -p -f $READS1 -f $READS2 > $SAMPLE.gam
vg pack -x $NAME.index.giraffe.gbz -g $SAMPLE.gam -o $SAMPLE.pack -Q 5
vg call $NAME.index.giraffe.gbz -k $SAMPLE.pack --ploidy 1 --sample $SAMPLE > $SAMPLE.vcf
```

A few VCF file lines:

```shell
#CHROM        POS   ID            REF  ALT  QUAL     FILTER  INFO                                                        FORMAT                    DBVPG6044
S288C#1#chrI  3991  >29061>29064  G    A    968.529  PASS    AT=>29061>29062>29064,>29061>29063>29064;DP=42              GT:DP:AD:GL:GQ:GP:XD:MAD  1:42:0,42:-97.899314,-1.347271:256:-0.693147:39.141644:42
S288C#1#chrI  4168  >29079>29082  G    A    577.724  PASS    AT=>29079>29080>29082,>29079>29081>29082;DP=25              GT:DP:AD:GL:GQ:GP:XD:MAD  1:25:0,25:-59.874843,-2.403389:256:-0.693147:38.983608:25
S288C#1#chrI  4373  >29093>29096  G    A    554.735  PASS    AT=>29093>29094>29096,>29093>29095>29096;DP=24              GT:DP:AD:GL:GQ:GP:XD:MAD  1:24:0,24:-57.061942,-1.889346:256:-0.693147:34.262753:24
S288C#1#chrI  4867  >29127>29132  C    T    1014.51  PASS    AT=>29127>29128>29130>29132,>29127>29129>29130>29132;DP=44  GT:DP:AD:GL:GQ:GP:XD:MAD  1:44:0,44:-102.468880,-1.319120:256:-1.098612:43.416668:44
S288C#1#chrI  4880  >29132>29135  A    T    140.941  PASS    AT=>29132>29133>29135,>29132>29134>29135;DP=6               GT:DP:AD:GL:GQ:GP:XD:MAD  1:6:0,6:-25.693250,-11.900101:137:-0.693147:43.416668:6
S288C#1#chrI  4908  >29138>29141  T    C    140.941  PASS    AT=>29138>29139>29141,>29138>29140>29141;DP=6               GT:DP:AD:GL:GQ:GP:XD:MAD  1:6:0,6:-25.958225,-12.165076:137:-0.693147:44.123718:6
S288C#1#chrI  6080  >29236>29239  A    G    968.529  PASS    AT=>29236>29237>29239,>29236>29238>29239;DP=42              GT:DP:AD:GL:GQ:GP:XD:MAD  1:42:0,42:-97.874503,-1.322460:256:-0.693147:40.090141:42
S288C#1#chrI  6108  >29239>29242  A    T    715.655  PASS    AT=>29239>29240>29242,>29239>29241>29242;DP=31              GT:DP:AD:GL:GQ:GP:XD:MAD  1:31:0,31:-72.963875,-1.699272:256:-0.693147:40.090141:31
S288C#1#chrI  6140  >29242>29245  C    T    991.517  PASS    AT=>29242>29243>29245,>29242>29244>29245;DP=43              GT:DP:AD:GL:GQ:GP:XD:MAD  1:43:0,43:-100.205969,-1.355068:256:-0.693147:40.090141:43
```

**NOTE:** it worked with `vg 1.40.0`, not with the latest `vg 1.52.0`.
