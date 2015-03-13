/ Pomoda: Peak Oriented Motif Discovery Algorithm
**http://compbio.ddns.nus.edu.sg/~chipseq/Pomoda**
**Function: generates position weight matrixes (PWM) for motifs concentrate around the ChIP-Seq peaks**
**Syntax:** POMODA -i inputFasta -o outputDIR [-w weightFile, -ot overlapThreshold -pdt PWM\_Divergence\_Threshold
**-ratio minSupportRatio, -FDR p-valueCutoff, -maxlen maxMotifLength, -seedlen SeedLength,-rs resolution\_bp, -mbr min\_binding\_range** -n numberOfMotifs]
 Comments:
**inputFasta: input data, using fasta format, recommended length per sequence is 10k** outputDIR: the directory for output of one task, different task should have different outputDIR, otherwise, the output position file will be replaced.
**weightFile: the file store the peak intensity, each line corresponding to one sequence in inputFasta, e.g., ":3251" in the 1st line, says the first sequence intensity is 3251, [optional](optional.md)** overlapThreshold: the binding site overlap percentage cutoff for filtering the redundant candidate motif. The higher value, the less motifs will be filtered. [default: 0.02]
**PWM\_Divergence\_Threshold: the dissimilarity cuttof for filtering the redundant candidate motif. The higher value, the more motifs will be filtered. [default:0.18]** minSupportRatio: the minimal proportion of sequences contains the reported moitf [default:0.05].
**FDR: the p-value cuttoff for testing uniformilty of motif distribution. [default:1.e-7]** maxMotifLength: the maximum length of reported motif. [default:26]
**SeedLength: the length of k-mer using in seed selection phase. [default:6]** rs: the resolution level of binding region estimation (unit:bp).  [default:100]
**min\_binding\_range: the minimal size of binding region. [default:100]** numberOfMotifs: the number of output motifs. [default:20]
 Author: Zhizhuo Zhang   10/10/2009
**Computational Biology Lab, School of Computing, National University of Singapore (NUS)** Email: zzz2010@gmail.com