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
    print("USE: %s file_name block_size\n" % sys.argv[0])
    sys.exit(1)

  try:
    blocksize =int(sys.argv[2])
  except:
    print("USE: %s file_name block_size\n" % sys.argv[0])
    sys.exit(1)

  if not os.path.isfile(infile):
    print("ERROR: file %s does not exists\n" % infile)
    sys.exit(1)

  # functions to be evaluated
  def id(x):
    return x
  def ratio(x):
    return (x[2]-x[1])/x[0]
  def test(x):
    return x[1]-x[0]

  # data acquisition
  indata=np.loadtxt(infile, skiprows=1, dtype=np.float)
  data=np.transpose(indata)     #column ordered

  ris, err = jack.jackknife_for_secondary(ratio, blocksize, [id, data[0]], [id, data[2]], [id, data[4]])
  print(ris, err)


  print("check on imaginary parts", end = ' ')
  num_columns=len(indata[0])
  for i in range(1, num_columns, 2):
    ris, err = jack.jackknife_for_primary(id, data[i], blocksize)
    print(ris, err, end=' ')
  print('')


