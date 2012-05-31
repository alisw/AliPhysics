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
definedSymbols = {}
libraries = {}
foundLibs = []
notFound = []
sys.argv = sys.argv[1: ]
def NM(libname):
        undefinedSymbols = []
        if libname.endswith('.so') == True:
                command = """nm -C """ + libname + """ | awk '/ U /{print ;}/ T /{print ; }' """
                for line in os.popen(command).readlines():
                        if '::' in line and not 'for' in line and not 'operator' in line and not 'std::' in line and not 'gnu_cxx::' in line and not 'non-virtual' in line: #To avoid looking for what we cannot find
                                line = line.strip()
                                if line.startswith('U'):
                                        line.strip()
                                        line = (line.lstrip('U')).strip()
                                        U = line.find('::')
                                        line = line[:U] #Take only the class name in the symbol 
                                        if line not in undefinedSymbols:
                                                undefinedSymbols.append(line)
                                elif 'T ' in line:
                                        line.strip()
                                        k = line.find('T')
                                        line = (line[k+1 : ]).strip()
                                        T = line.find('::')
                                        line = line[:T] #Take only the class name in the symbol
                                        if definedSymbols.has_key(line) == False:
                                                definedSymbols[line] = libname
        libraries[libname] = undefinedSymbols

def findLibs(libname,directory):
        Usymbols = libraries[libname]
        foundLibs.append(libname)
        i = 0
        while i < len(Usymbols):
                found = 0
                symbol = Usymbols[i]	
		files=os.listdir(".")
                files=[filename for filename in files if filename[0] != '.']
                for f in files:
                        if f.endswith('.so'):
                                if libraries.has_key(f) == False:
                                        NM(f)
                                if definedSymbols.has_key(symbol) and definedSymbols[symbol] in foundLibs:
                                        found = 1
                                        break
                                elif definedSymbols.has_key(symbol) and definedSymbols[symbol] not in foundLibs:
                                        print definedSymbols[symbol] #Print found library name
                                        found = 1
                                        findLibs(definedSymbols[symbol], directory)
                                        break
                else:
                        if directory < len(sys.argv) - 1:
                                os.chdir(sys.argv[directory])
                                directory = directory + 1	
				continue # Avoid incrementing the index we are searching
                        else:
				notFound.append("'"+symbol+"'" + " in " + libname)                                
                i = i + 1

starttime = time.time()
os.chdir(sys.argv[0])
print "Libraries " + sys.argv[-1] + " depends on; Printed in the order which they are called"
NM(sys.argv[-1])
findLibs(sys.argv[-1],1)
stoptime = time.time()
elapsed = stoptime - starttime
print "Running findLibDeps.py took %.3f seconds" % (elapsed)

if len(notFound) > 0:
        print "Symbols not found"
        for sim in notFound:
                print sim

