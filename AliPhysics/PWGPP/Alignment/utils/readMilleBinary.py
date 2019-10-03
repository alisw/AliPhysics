#!/usr/bin/python

## \file
# Read millepede binary file and print records
#
# Hardcoded defaults can be replaced by command line arguments for
#    -  Name of binary file
#    -  Number of records to print
#    -  Number of records to skip (optional)
#
# Description of the output from readMilleBinary.py
#    -  Records (tracks) start with \c '===' followed by record number and length 
#       (<0 for binary files containing doubles)
#    -  Measurements: A measurement with global derivatives is called a 'global measurement', 
#       otherwise 'local measurement'. Usually the real measurements from the detectors are 'global'
#       ones and virtual measurements e.g. to describe multiple scattering are 'local'.
#    -  'Global' measurements start with \c '-g-' followed by measurement number, first global label,
#       number of local and global derivatives, measurement value and error. The next lines contain 
#       local and global labels (array('i')) and derivatives (array('f') or array('d')).
#    -  'Local' measurements start with \c '-l-' followed by measurement number, first local label, 
#       number of local and global derivatives, measurement value and error. The next lines contain
#       local labels (array('i')) and derivatives (array('f') or array('d')).
#
# Tested with SL4, SL5, SL6

import array, sys

# ############### read millepede binary file #################
#  
## Binary file type (C or Fortran)
Cfiles = 1  # Cfiles
#Cfiles = 0 # Fortran files
# 
## Integer format 
intfmt = 'i'  # SL5, gcc-4
#intfmt = 'l' # SL4, gcc-3
#
## Binary file name
fname = "milleBinaryISN.dat"
#
## number of records (tracks) to show
mrec = 10
## number of records (track) to skip before 
skiprec = 0
#
# ## C. Kleinwort - DESY ########################

# ## use command line arguments ?
narg = len(sys.argv)
if narg > 1:
  if narg < 3:
    print " usage: readMilleBinary.py <file name> <number of records> [<number of records to skip>]" 
    sys.exit(2)
  else:
    fname = sys.argv[1]
    mrec = int(sys.argv[2])
    if narg > 3:
      skiprec = int(sys.argv[3])

#print " input ", fname, mrec, skiprec
  
f = open(fname, "rb")

nrec = 0
try:
    while (nrec < mrec + skiprec):
# read 1 record    
        if (Cfiles == 0): 
           lenf = array.array(intfmt)
           lenf.fromfile(f, 2)
           
        length = array.array(intfmt)
        length.fromfile(f, 1)
        nr = abs(length[0] / 2)
        nrec += 1

        if length[0] > 0: 
            glder = array.array('f')
        else:
            glder = array.array('d')          
        glder.fromfile(f, nr)

        inder = array.array(intfmt)
        inder.fromfile(f, nr)
        
        if (Cfiles == 0): 
           lenf = array.array(intfmt)
           lenf.fromfile(f, 2)

        if (nrec <= skiprec):  # must be after last fromfile
            continue

        print " === NR ", nrec, length[0] / 2

        i = 0
        nh = 0
        ja = 0
        jb = 0
        jsp = 0
        nsp = 0
        while (i < (nr - 1)):
            i += 1
            while (i < nr) and (inder[i] != 0): i += 1
            ja = i
            i += 1
            while (i < nr) and (inder[i] != 0): i += 1
            jb = i
            i += 1
            while (i < nr) and (inder[i] != 0): i += 1
            i -= 1
# special data ?
            if (ja + 1 == jb) and (glder[jb] < 0.):
               jsp = jb
               nsp = int(-glder[jb])
               i += nsp
               print ' ### spec. ', nsp, inder[jsp + 1:i + 1], glder[jsp + 1:i + 1]
               continue
            nh += 1           
            if (jb < i):
# measurement with global derivatives
               print ' -g- meas. ', nh, inder[jb + 1], jb - ja - 1, i - jb, glder[ja], glder[jb]
            else:
# measurement without global derivatives
               print ' -l- meas. ', nh, inder[ja + 1], jb - ja - 1, i - jb, glder[ja], glder[jb]            
            if (ja + 1 < jb):
               print " local  ", inder[ja + 1:jb]
               print " local  ", glder[ja + 1:jb]
            if (jb + 1 < i + 1):
               print " global ", inder[jb + 1:i + 1]
               print " global ", glder[jb + 1:i + 1]

               
except EOFError:
     pass
#    print "end of file"

f.close()
