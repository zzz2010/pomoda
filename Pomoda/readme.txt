#*******************************************************
# Developed by Zhizhuo Zhang (zzz2010@gmail.com), 
# National University of Singapore.
# 
# This software is provided "AS IS".  ZZZ makes no warranties, express
# or implied, including no representation or warranty with respect to
# the performance of the software and derivatives or their safety,
# effectiveness, or commercial viability.  ZZZ does not warrant the
# merchantability or fitness of the software and derivatives for any
# particular purpose, or that they may be exploited without infringing
# the copyrights, patent rights or property rights of others. NIA shall
# not be liable for any claim, demand or action for any loss, harm,
# illness or other damage or injury arising from access to or use of the
# software or associated information, including without limitation any
# direct, indirect, incidental, exemplary, special or consequential
# damages.
# 
# This software program may not be sold, leased, transferred, exported
# or otherwise disclaimed to anyone, in whole or in part, without the
# prior written consent of ZZZ.
#*******************************************************

JPomoda ToolKit is designed to find enriched DNA motifs in ChIP-chip or 
ChIP-seq data and to scan sequences for motif matches.
Disclaimer: these programs may malfunction without generating an error 
message or crash if you use wrong syntax or inappropriate input file format.

1.Installation:
Download the executable jar file, and install JVM 6.0 or above. 

2.Program "JPomoda"
example:
java -jar JPomoda.jar -i ./data/pos.fa  -c ./data/neg.fa   -oops -mask  -prefix ./output/

usage: JPomoda
 -bgmodel <arg>   background model file
 -c <arg>         control fasta file
 -clust <arg>     linkage type of hierachical clustering:[SINGLE,
                  COMPLETE, UPGMA, WPGMA, UPGMC, WPGMC, Ward]
 -FDR <arg>       fasle positive rate
 -i <arg>         input fasta file
 -mask            whether marking the top motif location in order to find
                  co-motif
 -maxlen <arg>    maxmimum length of the motif (default 30)
 -n <arg>         number of motifs in final report (default 5)
 -oops            whether assuming only one occurrence per sequence
 -prefix <arg>    output directory
 -ratio <arg>     sampling ratio (default 0.8)
 -rs <arg>        counting resolution (default 40 bp)
 -seedlen <arg>   kmer seed motif length (default 5)
 -supp <arg>      minimum support ratio, the percentage of peaks contains
                  motif (default 0.05)
                  
                  
                  
                  
3.Program "PWMevaluator"
example:
java -jar PWMevaluator.jar -i ./data/pos.fa  -c ./data/neg.fa  -pwm ./data/test.pwm -convert -roc -match ./data/jaspar.pwm 

usage: PWMevaluator
 -bgmodel <arg>   background model file
 -c <arg>         control fasta file
 -convert         convert input PWM file to the transfac format
 -FDR <arg>       fasle positive rate
 -i <arg>         input fasta file
 -match <arg>     find similar motifs in known PWM library (path to the library, e.g.,
                  jaspar.pwm)
 -maxlen <arg>    maxmimum length of gap (default 8)
 -prefix <arg>    output directory
 -pwm <arg>       input PWM file
 -ratio <arg>     sampling ratio (default 1)
 -roc             compute AUC and draw ROC curve for the given pwm file
 -thresh <arg>    minimum entropy threshold for considering a position as
                  a gap(default 0.5)

4.Program "Pomoda.py"
example:
python Pomoda.py chipseq_peak.txt 200 /home/XXX/genome/hg18  /home/XXX/output

usage:
Pomoda.py peakfile peakwidth genomedir <outputdir>

###before you run it###
please open Pomoda.py and set up the locations for two external programs (by Dr. Chengwei Chang,changcw99@gis.a-star.edu.sg ):
genlogo=""#format("python {rootdir}/PROG/genlogo/genlogo_multi.py",locals())
extractfas=""#format("python {rootdir}/PROG/extractfas.py",locals())


Output:
jpomoda_raw.pwm : a set of redundant enriched PWM motifs (before clustering)
jpomoda_clust.pwm: a set of redundant enriched PWM motifs (after clustering)
Motif_clustX.logo.png: PWM logo
jpomoda_clust.pwm_eval.txt: AUC result, and similar known motifs in JASPAR library
jpomoda_clust.pwm_roc.png: ROC curve figure
jpomoda_clust.pwm_sorted.pwm: PWM motifs sorted by AUC value
