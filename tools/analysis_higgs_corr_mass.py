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

import numpy as np
import os.path
import statanalysis.jackknife as jack

# main
if __name__=="__main__":

  try:
    infile =sys.argv[1]
  except:
    print("USE: %s file_name therm block_size\n" % sys.argv[0])
    sys.exit(1)

  try:
    therm =int(sys.argv[2])
  except:
    print("USE: %s file_name therm block_size\n" % sys.argv[0])
    sys.exit(1)

  try:
    blocksize =int(sys.argv[3])
  except:
    print("USE: %s file_name therm block_size\n" % sys.argv[0])
    sys.exit(1)

  if not os.path.isfile(infile):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)

  f=open(infile, 'r')
  firstline=f.readline()
  f.close()
  L=int(firstline.split()[2])

  # functions to be evaluated
  def id(x):
    return x
  def logratio(x):
    return 1./np.log(x[1]/x[0])

  # data acquisition
  indata=np.loadtxt(infile, skiprows=therm, dtype=np.float)
  data=np.transpose(indata)     #column ordered

  #print(blocksize, end=' ')

  for i in range(L):
    print(i, end=' ')

    ris, err = jack.jackknife_for_secondary(logratio, blocksize, [id, data[9+i+1+L]], [id, data[9+i+L]])
    print(ris, err, end=' ')

    print('')

