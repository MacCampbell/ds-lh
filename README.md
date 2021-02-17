# ds-lh
Delta Smelt life history genetic association testing

## Structure of Files

I have used links to point to files that already have been prepared.    

__1__ /genome/ Hypomesus-transpacificus_10X_F_A.pseudohap2.1.fasta -> /group/millermrgrp2/shannon/projects/genome_Hypomesus-transpacificus/data/04-supernova_out/60x/A/Hypomesus-transpacificus_10X_F_A.pseudohap2.1.fasta     

__2__ /bams/ Ht33-95_2017_G12.sort-n.fixmate-m.sort.markdup-r.bam -> /group/millermrgrp2/shannon/projects/DS_history/data/10X_alignments/Ht33-95_2017_G12.sort-n.fixmate-m.sort.markdup-r.bam (and indices)     

__3__ /outputs/ Are meant to have outputs for files in series 1xx, 2xx and so on (outputs/100 and outputs/200 for example)

# Notes as I go

### 101-test-asso.sh
My original subset didn't line up 1:1. For example:
maccamp\@farm:~/ds-lh$ find /group/millermrgrp2/shannon/projects/DS_history/data/ -name Ht20-6_2012_F01.sort-n.fixmate-m.sort.markdup-r.bam

There are several.
Ht20-6_2012_F01.sort-n.fixmate-m.sort.markdup-r.bam
Ht20-3_2012_C01.sort-n.fixmate-m.sort.markdup-r.bam
Ht19-7_2012_G01.sort-n.fixmate-m.sort.markdup-r.bam
I'm going back to redo the bamlist based on the files I have.

### New Reference 20210204 Version
Let's do some comparisons between the original data and the new genome.     
`maccamp@farm:~/genomes/hypomesus-20210204$ srun -p high -t 8:00:00 bwa index Hyp_tra_F_20210204.fa`     

Things to think about:     
__1__ Association    
__2__ Linkage    
__3__ Network Analysis      

Should have 369 samples with sequence data, see "metadata/expanded-meta.csv". Starting 1200 series.    