"""
  pybsseq
"""
import numpy as np
import numpy.ma
import pandas as pd
import argparse
import logging
import sys
import os.path
import scipy.stats as st
import tabix
import h5py

logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

def read_bsdata(inFile):
  BSData = pd.read_table(inFile, header=None, compression = "gzip")
  BSData = BSData.as_matrix(columns=BSData.columns[0:])
#  BS_chrom = np.array(BSData[0], dtype="string")
#  BS_pos = np.array(BSData[1], dtype="int")
#  BS_strand = np.array(BSData[5], dtype="string")
#  BS_depcon =  pd.DataFrame(list(BSData[3].str.split(":")))

  np.core.defchararray.rsplit



  BS_depth = np.array(BS_depcon[1], dtype=int)
  BS_context = np.array(BS_depcon[0], dtype="string")
  BS_perMeth = np.array(BSData[4], dtype="float")
  del BSData
  return (BS_chrom, BS_pos, BS_strand, BS_context, BS_depth, BS_perMeth)

def getErrorMeth(BS_chrom, BS_pos, BS_depth, BS_perMeth, chrName):
  chrInd = np.where(BS_chrom == chrName)[0]
  methPer = np.sum(np.multiply(BS_perMeth[chrInd], BS_depth[chrInd]))/np.sum(BS_depth[chrInd])
  return (methPer, np.mean(BS_depth[chrInd]))

def BinomTest(per_n, n, p, alternative="greater"):
  tpVal = st.binom_test(per_n * n, n, p, alternative = alternative)
  return tpVal

def getMPs_readData(args):
  logging.info("Reading the input file")
  (BS_chrom, BS_pos, BS_strand, BS_context, BS_depth, BS_perMeth) = read_bsdata(args['inFile'])
  logging.info("Calculating the conversion rate")
  ConversionRate = getErrorMeth(BS_chrom, BS_pos, BS_depth, BS_perMeth, "chrC")
  npBinomTest = np.vectorize(BinomTest)
  window = 300000
  BS_pVal = np.zeros(0, dtype=float)
  for i in range(0, len(BS_chrom), window):
    pVal = npBinomTest(BS_perMeth[i:i+window], BS_depth[i:i+window], ConversionRate[0])
    BS_pVal = np.append(BS_pVal, pVal)
    logging.info("Done analyzing %s number of positions", i)
 np.savetxt(args['outFile'], numpy.column_stack((, TotNonMatPos)))


def getMPs(args):
  logging.info("Reading the input file")
  (BS_chrom, BS_pos, BS_strand, BS_context, BS_depth, BS_perMeth) = read_bsdata(args['inFile'])
  BS_no = len(BS_chrom)
  logging.info("Calculating the conversion rate")
  ConversionRate = getErrorMeth(BS_chrom, BS_pos, BS_depth, BS_perMeth, "chrC")
  logging.info("Conversion rate: %s with a average depth in chloroplast at %s", ConversionRate[0], ConversionRate[1])
  h5file = h5py.File(args['outFile'], 'w')
  h5file.create_dataset('chromosome', data=BS_chrom, shape=(BS_no,))
  del BS_chrom
  h5file.create_dataset('position', data=BS_pos, shape=(BS_no,))
  del BS_pos
  h5file.create_dataset('strand', data=BS_strand, shape=(BS_no,))
  del BS_strand
  h5file.create_dataset('context', data=BS_context, shape=(BS_no,))
  del BS_context
  h5file.create_dataset('depth', data=BS_depth, shape=(BS_no,), dtype="int8")
  h5file.create_dataset('perMeth', data=BS_perMeth, shape=(BS_no,))
  logging.info("Reading the file for Methylated Positions")
  h5file.create_dataset('pVal', shape=(BS_no,), dtype='float')
  npBinomTest = np.vectorize(BinomTest)
  window = 300000
  for i in range(0, BS_no, window):
    pVal = npBinomTest(BS_perMeth[i:i+window], BS_depth[i:i+window], ConversionRate[0])
    h5file['pVal'][i:i+window] = np.array(pVal, dtype=float)
    logging.info("Done analyzing %s number of positions", i+window)
  h5file.close()

