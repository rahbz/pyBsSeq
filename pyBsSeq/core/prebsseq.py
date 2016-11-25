"""
  pybsseq
"""
import numpy as np
import argparse
import logging
import sys
import os
import scipy.stats as st
import vcfnp
import subprocess

def setLog(args):
  if args['logDebug']:
    numeric_level = getattr(logging, "DEBUG", None)
  else:
    numeric_level = getattr(logging, "CRITICAL", None)
  logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=numeric_level)

def readVcf(inFile):
  bvcf = vcfnp.variants(inFile, cache=True).view(np.recarray)
  bvcfD = vcfnp.calldata_2d(inFile, cache=True).view(np.recarray)
  return(bvcf, bvcfD)

def BinomTest(per_n, n, p, alternative="greater"):
  tpVal = st.binom_test(per_n * n, n, p, alternative = alternative)
  return tpVal
 
 def getErrorRate(bsCHROM, bsCONTEXT, bsMethPer, bstC, chrs = "ChrC"):
  chrInd = np.where(bsCHROM == chrs)[0]
  contInd = chrInd[np.where((bsCONTEXT[chrInd] == 'CG') | (bsCONTEXT[chrInd] == 'CHG') | (bsCONTEXT[chrInd] == 'CHH'))[0]]
  chrMethPer = bsMethPer[contInd]
  chrDepth = bstC[contInd]
  conv_rate = np.nansum(np.multiply(chrMethPer, chrDepth))/np.nansum(chrDepth)
  return conv_rate

def callMPs(bsMethPer, bstC, error_rate, window=300000):
  bsPval = np.zeros(0,dtype=float)
  npBinomTest = np.vectorize(BinomTest)
  for i in range(0, len(bsMethPer), window):
    pVal = npBinomTest(bsMethPer[i:i+window], bstC[i:i+window], error_rate)
    bsPval = np.append(bsPval, pVal)
  return bsPval

def writeBED(bsCHROM, bsPOS, bsCONTEXT, bstC, bsMethPer, bsPval, bsSTRAND, outBED):
  out = open(outBED, 'w')
  for i in len(bsPOS):
    out.write("%s\t%s\t%s\t%s:%s:%s\t%s\t%s\n" % (bsCHROM[i], bsPOS[i], bsPOS[i]+1, bsCONTEXT[i], bstC[i], bsMethPer[i], bsPval[i], bsSTRAND[i]))
  out.close()

def getMPsfromVCF(args):
  (bVCF, bvcfD) = readVCF(args['inVCF'])
  error_rate = getErrorRate(bvcf.CHROM, bvcf.CX, bvcfD.BT[:,0], bvcfD.CV[:,0], chrs = "ChrC")
  logging.info("Conversion rate: %s", 100 - error_rate * 100)
  ChrsNums = np.array(("Chr1","Chr2","Chr3","Chr4","Chr5"))
  MethInd = np.where((np.in1d(bvcf.CHROM, ChrsNums)) & (np.in1d(bvcf.CX, MethContext)))[0]
  logging.info("Number of positions: %s", len(MethInd))
  bsPval = callMPs(bvcfD.BT[MethInd,0], bvcfD.CV[MethInd,0], error_rate, window=args['window'])
  bsSTRAND = np.core.defchararray.replace(np.core.defchararray.replace(bvcf.REF[MethInd], "C", "+"), "G", "-")
  logging.info("writing MPs info in out bed file")
  writeBED(bvcf.CHROM[MethInd], bvcf.POS[MethInd], bvcf.CX[MethInd], bvcfD.CV[MethInd,0], bvcfD.BT[MethInd,0], bsPval, bsSTRAND, outBED)
  logging.info("finished!")

def getMPsBED(args):
  ## input the split bed file
  










