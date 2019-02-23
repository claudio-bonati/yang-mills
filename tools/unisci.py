#!/usr/bin/env python3

import sys
import os

try:
  with open("statanalysis/jackknife.py") as f: 
    pass
except:
  print("Get the package 'statanalysis' from")
  print("https://github.com/claudio-bonati/statanalysis")
  print("and put it in this directory.")
  print()
  print("If you have 'git' installed you can use")
  print("git clone https://github.com/claudio-bonati/statanalysis")
  print()
  sys.exit(1)

sys.path.append(os.getcwd()+"/statanalysis")

import glob
import os.path
import statanalysis.concatenate as conc

#***************************
if __name__=="__main__":

  try:
    indir =sys.argv[1]
  except:
    print("USE: %s dir_name thermalization\n" % sys.argv[0])
    sys.exit(1)

  try:
    therm =int(sys.argv[2])
  except:
    print("USE: %s dir_name thermalization\n" % sys.argv[0])
    sys.exit(1)

  if not os.path.isdir(indir):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)

  if therm<0:
    print("ERROR: 'therm' must be a nonnegative integer\n")
    sys.exit(1)

  conc.concatenate(indir, therm)
