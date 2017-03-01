"""
  pybsseq
"""
import numpy as np
import argparse
import logging
import sys
import os, os.path
import scipy.stats as st
import vcfnp
import subprocess
from methylpy.call_mc import *
from methylpy.DMRfind import DMRfind
import glob

log = logging.getLogger(__name__)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

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
  log.info("Conversion rate: %s", 100 - error_rate * 100)
  ChrsNums = np.array(("Chr1","Chr2","Chr3","Chr4","Chr5"))
  MethInd = np.where((np.in1d(bvcf.CHROM, ChrsNums)) & (np.in1d(bvcf.CX, MethContext)))[0]
  log.info("Number of positions: %s", len(MethInd))
  bsPval = callMPs(bvcfD.BT[MethInd,0], bvcfD.CV[MethInd,0], error_rate, window=args['window'])
  bsSTRAND = np.core.defchararray.replace(np.core.defchararray.replace(bvcf.REF[MethInd], "C", "+"), "G", "-")
  log.info("writing MPs info in out bed file")
  writeBED(bvcf.CHROM[MethInd], bvcf.POS[MethInd], bvcf.CX[MethInd], bvcfD.CV[MethInd,0], bvcfD.BT[MethInd,0], bsPval, bsSTRAND, outBED)
  log.info("finished!")

def methylpy_callmcs(args):
    files = [args['inFile']]
    libraries = [1]
    sample = args['sample_id']
    f_reference = args['ref_fol'] + '_f'  ## for forward ref
    r_reference = args['ref_fol'] + '_r'
    reference_fasta = args['ref_fasta']
    uMeth = args['unMeth']
    mem = args['memory']
    procs = args['nt']
    run_methylation_pipeline(files,libraries,sample,f_reference,r_reference,reference_fasta,uMeth,sort_mem=mem,num_procs=procs,sig_cutoff=0.01,min_cov=3)

def methylpy_dmrfind(args):
    chromosomes = ["1","2","3","4","5"]
    region_dict = {}
    for chrom in chromosomes:
        region_dict[chrom] = [0, 5000000000]
    mc_type = args['mc_type'].split(",")
    samples = args['sample_id'].split(",")
    if args['sample_cat'] == '0':
        sample_cat=range(0, len(samples))
    else:   ## if samples given some categories
        sample_cat = args['sample_cat'].split(",")
    ### Check which of the samples given have the all c files
    dbacc=np.zeros(0)
    for eacc in glob.glob(args['path_to_allc']+"/*_1.tsv"):
        accID = os.path.basename(eacc).replace("allc_", "").replace("_1.tsv", "")
        dbacc = np.append(dbacc, accID)
    reqAcc=np.array(samples)[np.where(np.in1d(np.array(samples), dbacc))[0]]
    log.info("Number of accessions present in allc: %s", len(reqAcc))
    log.info(reqAcc.tolist())
    DMRfind(mc_type=mc_type, region_dict=region_dict, samples=reqAcc.tolist(), path_to_allc=args['path_to_allc'], num_procs=args['nt'], save_result=args['outDMR'], min_cov=1)
