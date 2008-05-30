#!/usr/bin/env python

#/**************************************************************************
# * This file is property of and copyright by the ALICE HLT Project        *
# * All rights reserved.                                                   *
# *                                                                        *
# * Primary Authors:                                                       *
# *   Seforo Mohlalisi <seforomohlalisi@yahoo.co.uk>                       *
# *                                                                        *
# * Permission to use, copy, modify and distribute this software and its   *
# * documentation strictly for non-commercial purposes is hereby granted   *
# * without fee, provided that the above copyright notice appears in all   *
# * copies and that both the copyright notice and this permission notice   *
# * appear in the supporting documentation. The authors make no claims     *
# * about the suitability of this software for any purpose. It is          *
# * provided "as is" without express or implied warranty.                  *
# **************************************************************************/

# This script tries to find all dependencies of a library on other libraries.
#
# The first argument to the script should be the directory where your library
# of interest resides in. i.e. for
#  > find_libs directory1 directory2 ... directoryN libXX.so
# The library for which we want to find dependencies, libXX.so, must reside
# in directory1. All other directories are optional, but the last argument must
# be the library name of the library for which we want to find dependencies.

import sys
import os
import commands
import time
libs_names = []
used_symbols = []
unfound_symbols = []
sys.argv = sys.argv[1: ]
check = len(sys.argv) - 1
def libsfind(libname):
    if libname == "sh": #Could be helpful to catch some unwanted behaviour (This line is subject to change)
       sys.exit()
    libs_names.append(libname)
    command = """nm -C """ + libname  +  """ | grep " U " | grep :: | grep -v for | grep -v operator """
    list = []
    for line in os.popen(command).readlines():
        string = line.strip()
        string = (string.lstrip('U ')).strip()
        k = string.find('::')
        string = string[ : k+2]
        if not string in list:
           list.append(string)
           
    whi = 0 # Iterator for the while loop
    for i in range(len(list)):
        if list[i] == "n":
           print "Using n"
        if not list[i] in used_symbols:
           print "Looking for undefined symbol " + list[i] + " found library " + libname
           used_symbols.append(list[i])
           found = 0 #Indicates where the symbol is found or not
           round = 1 #keep track of how many directories we have searched
           while whi < check: #Used basically for changing directories
                 directory = sys.argv[whi]
                 os.chdir(directory)
                 files=os.listdir(".")
                 files=[filename for filename in files if filename[0] != '.']
                 for f in files:
                     com1 = "nm -AC " + f
                     char = "T " + list[i]
                     char = '"' + char + '"'
                     com2 = """ | grep -F """ +  char
                     command = com1 + com2
                     for line in os.popen(command).readlines():
                         string = line.strip()
                         if string[ : 2] != "nm": #Incase nm command returned error
                            k = string.find(':')
                            string = string[ : k] #Take only the library file name
                            if string in libs_names:
                               found = 1
                               print "Found in " + string
                               break
                            found = 1
                            print "Found in " + string
                            libsfind(string)
                            break
                     if found == 1:
                        break
                 if found == 1:
                    break
                 elif round < len(sys.argv):
                      round = round + 1
                 else:
                     unfound_symbols.append(list[i])
                     whi = 0
                     print "We can't find it"
                     break
                 whi = whi + 1
                 if whi == check:
                    whi = 0
                 
os.chdir(sys.argv[0])
starttime = time.time()
libsfind(sys.argv[-1])
stoptime = time.time()
elapsed = stoptime - starttime

print "Found libraries in the order which they are called"
for li in libs_names:
    print li
if len(unfound_symbols) > 0:
   print "Symbols not found are "
   for kk in unfound_symbols:
        print kk
print "Running findlibs took %.3f " % (elapsed)
