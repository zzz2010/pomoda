#!/usr/bin/python

import os
import sys
#from common import *


genlogo=""#format("python {rootdir}/PROG/genlogo/genlogo_multi.py",locals())
extractfas=""#format("python {rootdir}/PROG/extractfas.py",locals())



peakfile=sys.argv[1]
peakwidth=sys.argv[2]
genomedir=sys.argv[3]
toolpath=sys.path[0]


try:

	outputdir="."
	if len(sys.argv)>4:
		outputdir=sys.argv[4]
	
	ossystem(extractfas+" -genomedir "+genomedir+" -peakfile "+peakfile+" -w "+peakwidth+" > "+outputdir+"/posfa.fa")
	outBG=open(outputdir+"/bg.peak",'w')
	lines=open(peakfile,"r").readlines()
	bias=int(peakwidth)*3
	for line in lines:
		comps=line.strip().split()
		chrom=comps[0]
		leftbg=list()
		rightbg=list()
		leftbg.append(chrom)
		rightbg.append(chrom)
		for el in comps[1:]:
			left1=int(el)-bias
			if(left1>0):
				leftbg.append(str(left1))
			right1=int(el)+bias
			rightbg.append(str(right1))
		if(len(leftbg)==len(rightbg)):
			outBG.write("\t".join(leftbg)+"\n")
		outBG.write("\t".join(rightbg)+"\n")
	ossystem(extractfas+" -genomedir "+genomedir+" -peakfile "+outputdir+"/bg.peak"+" -w "+peakwidth+" > "+outputdir+"/negfa.fa")
	outBG.close()
except Exception,e:
	print("Pomoda.py peakfile peakwidth genomedir <outputdir>")
	print e
	
	
posseqFile=outputdir+"/posfa.fa"
negseqFile=outputdir+"/negfa.fa"
runmode="-mask -oops -clust WPGMA"

print "running program......"
ossystem("time java -d64 -Xmx8g -jar "+toolpath+"/JPomoda.jar -prefix "+outputdir+"/  -i  "+posseqFile+" -c "+negseqFile+" "+runmode+"  > "+outputdir+"/JPomoda.out")

pwmfile=outputdir+"/jpomoda_clust.pwm"

print "converting and evaluating......"
ossystem("time java -Djava.awt.headless=true -jar "+toolpath+"/JEvaluator.jar -i "+posseqFile+" -c "+negseqFile+" -pwm "+pwmfile+" -convert -roc -match "+toolpath+"/jaspar.pwm")


print "creating PWM logo ......"

command=format("{genlogo} -pwmfile {pwmfile} -prefix {outputdir}/",locals())
if ossystem(command):
        raise Exception("genlogo failed")

