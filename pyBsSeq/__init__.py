"""
    pyBsSeq
    ~~~~~~~~~~~~~
    A python toolkit for analysing the bisulfite-sequence data
    :copyright: Rahul Pisupati @ 2016
    :license:   GMI
"""
import os
import argparse
import sys
from pyBsSeq.core import prebsseq
import logging, logging.config

__version__ = '1.0.0'
__updated__ = "1.3.2017"
__date__ = "10.12.2016"

def setLog(logDebug):
  log = logging.getLogger()
  if logDebug:
    numeric_level = getattr(logging, "DEBUG", None)
  else:
    numeric_level = getattr(logging, "ERROR", None)
  log_format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
  lch = logging.StreamHandler()
  lch.setLevel(numeric_level)
  lch.setFormatter(log_format)
  log.setLevel(numeric_level)
  log.addHandler(lch)

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)

def get_options(program_license,program_version_message):
    inOptions = argparse.ArgumentParser(description=program_license)
    inOptions.add_argument('-V', '--version', action='version', version=program_version_message)
    subparsers = inOptions.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')
    callmc_parser = subparsers.add_parser('callmc', help="Call mc using methylpy from fastq file")
    callmc_parser.add_argument("-i", "--input_file", dest="inFile", help="Input fastq file for methylpy")
    callmc_parser.add_argument("-s", "--sample_id", dest="sample_id", help="unique sample ID for allc Files")
    callmc_parser.add_argument("-r", "--ref_fol", dest="ref_fol", help="methylpy reference folder for indices and refid", default="/home/GMI/rahul.pisupati/TAiR10_ARABIDOPSIS/03.methylpy.indices/tair10")
    callmc_parser.add_argument("-f", "--ref_fasta", dest="ref_fasta", help="reference fasta file", default="/home/GMI/rahul.pisupati/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta")
    callmc_parser.add_argument("-n", "--nt", dest="nt", help="number of threads", default=2, type="int")
    callmc_parser.add_argument("-c", "--unMethylatedControl", dest="unMeth", help="unmethylated control", type="string", default="ChrC:")
    callmc_parser.add_argument("-m", "--mem", dest="memory", help="memory for sorting", type="string", default="2G")
    callmc_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    callmc_parser.set_defaults(func=callmcs_onesample)
    dmr_parser = subparsers.add_parser('dmrfind', help="Identify DMR using methylpy")
    dmr_parser.add_argument("-s", "--sample_ids", dest="sample_ids", help="sample ids, comma seperated", type="string")
    dmr_parser.add_argument("-r", "--sample_categories", dest="sample_cat", help="sample categories indicating replicates, comma separated", type="string", default="0")
    dmr_parser.add_argument("-p", "--path", dest="path_to_allc", help="path to allc files", type="string")
    dmr_parser.add_argument("-c", "--context", dest="mc_type", help="methylation context, context separated", type="string")
    dmr_parser.add_argument("-n", "--nt", dest="nt", help="number of threads", default=2, type="int")
    dmr_parser.add_argument("-o", "--outDMR", dest="outDMR", help="output file for DMR", type="string")
    dmr_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
    dmr_parser.set_defaults(func=dmrfind)
    #lowfreq_parser = subparsers.add_parser('callLowFreq', help="Get lowfreq positions from allc files")
    return inOptions

def callMPsfromVCF(args):
  if not args['inFile']:
    die("input file not specified")
  if not args['outFile']:
    die("output file not specified")
  if not os.path.isfile(args['inFile']):
    die("input file does not exist: " + args['inFile'])
  prebsseq.getMPsfromVCF(args)

def callmcs_onesample(args):
    if not args['inFile']:
        die("input file not specified")
    if not os.path.isfile(args['inFile']):
        die("input file does not exist: " + args['inFile'])
    prebsseq.methylpy_callmcs(args)

def dmrfind(args):
    prebsseq.methylpy_dmrfind(args)

def main():
  ''' Command line options '''
  program_version = "v%s" % __version__
  program_build_date = str(__updated__)
  program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
  program_shortdesc = "The main module for pyBsHap"
  program_license = '''%s
  Created by Rahul Pisupati on %s.
  Copyright 2016 Gregor Mendel Institute. All rights reserved.

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
USAGE
''' % (program_shortdesc, str(__date__))

  parser = get_options(program_license,program_version_message)
  args = vars(parser.parse_args())
  setLog(args['logDebug'])
  if 'func' not in args:
    parser.print_help()
    return 0
  try:
    args['func'](args)
    return 0
  except KeyboardInterrupt:
    return 0
  except Exception as e:
    logging.exception(e)
    return 2

if __name__=='__main__':
  sys.exit(main())
