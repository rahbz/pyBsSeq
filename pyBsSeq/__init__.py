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

def die(msg):
  sys.stderr.write('Error: ' + msg + '\n')
  sys.exit(1)


def get_options():
  inOptions = argparse.ArgumentParser()
  subparsers = inOptions.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')

  getmps_parser = subparsers.add_parser('callMPsfromVCF', help="Get the methylated positions in genome")
  getmps_parser.add_argument("-i", "--input_file", dest="inFile", help="VCF file from biscuit after aligning the reads")
  getmps_parser.add_argument("-o", "--output", dest="outFile", help="Output file with the probability scores non corrected")
  getmps_parser.add_argument("-v", "--verbose", action="store_true", dest="logDebug", default=False, help="Show verbose debugging output")
  getmps_parser.set_defaults(func=callMPsfromVCF)
  
  
  return inOptions

def callMPsfromVCF(args):
  if not args['inFile']:
    die("input file not specified")
  if not args['outFile']:
    die("output file not specified")
  if not os.path.isfile(args['inFile']):
    die("input file does not exist: " + args['inFile'])
  prebsseq.getMPsfromVCF(args)

def main():
  parser = get_options()
  args = vars(parser.parse_args())
  if 'func' not in args:
    parser.print_help()
    return 0
  try:
    args['func'](args)
    return 0
  except KeyboardInterrupt:
    return 0
  except Exception as e:
    return 2


if __name__=='__main__':
  sys.exit(main())

