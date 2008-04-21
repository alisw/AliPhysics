#!/usr/bin/python

import sys
import os
import re
import string
import getopt

"""
Given a directory, will look into lib*.pkg files to produce a dependency graph
 of libraries (and DA if there are some)
"""

__author__ = "L. Aphecetche aphecetc_at_in2p3_dot_fr"
__version__ = "$Id$"

#_______________________________________________________________________________
def usage():
  """Describe usage of script
  """
  print "Usage: %s [-h | --help] [-d | --debug] [--noda] directory_to_scan" % sys.argv[0]
  sys.exit(1)
  
#_______________________________________________________________________________
def getSourceFiles(lib):
  """Extract the list of classes from a libXXX.pkg file
  """

  f = open(lib)
  sourcefiles = []
  for line in f:
    l = line.strip()
    if re.search('Ali',l) and re.search('.cxx',l):
      l = re.sub('SRCS',' ',l)
      l = re.sub(':=',' ',l)
      l = re.sub('=',' ',l)
      l = re.sub("\\\\",' ',l)
      l = re.sub(".cxx",' ',l)
      for i in l.split():
        sourcefiles.append(i)
  f.close()
  return sourcefiles
  
#_______________________________________________________________________________
def getIncludeFiles(srcfile):
  """Extract the list of included classes from a class
  """

  includes = []

  try:
    f = open("%s.cxx" % srcfile)
  except:
    print "Could not open file %s.cxx" % srcfile
    return includes
    
  for line in f:
    line = line.strip()
    if re.search("^#",line) and re.search("#include",line) and re.search('Ali',line):
      line = re.sub("#include",' ',line)
      i = line.index(".h")
      line = line[:i]
      line = re.sub("\"",' ',line)
      line = line.strip()
      includes.append(line)
  f.close()
  return includes
  
#_______________________________________________________________________________
def unique(list):
  """Extract a unique list from list
  """
  d = {}
  for l in list:
    d[l] = 1
  return d.keys()
  
#_______________________________________________________________________________
def findLibrary(file,allfiles):
  """Find in which library a given class is defined
  """
  for k,v in allfiles.items():
    for f in v:
      a,f = os.path.split(f)
      if file == f:
        return k
  return "external"
  
#_______________________________________________________________________________
def shorten(libname):
  """From libYYYxxx.pkg to YYYxxx
  """
  s = libname
  if re.search("lib",libname):
    s = re.sub("lib","",s)
    s = re.sub("\.pkg","",s)
  return s
  
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
def main():

  debug = False
  noda = False
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hd",["help", "debug","noda"])
  except getopt.GetoptError:
    print "Error in options"
    usage()

  for o, a in opts:
    if o in ( "-d","--debug" ):
      debug = True
    elif o in ( "-h","--help" ):
      usage()
      sys.exit()
    elif o == "--noda":
      noda = True
    else:
      assert False, "unhandled option"

  dir = args[0]
  
  # find the libraries defined in this directory (looking for libXXXX.pkg files)  
  libraries = []
  
  for file in os.listdir(dir):
    if re.search('^lib',file) and re.search('.pkg$',file):
      libraries.append(file)

  # append fake libraries for DAs
  if not noda:
    libraries.append("libMUONTRKda.pkg")
    libraries.append("libMUONTRGda.pkg")
  
  # srcfiles is a dictionary :
  # srcfiles[libXXX.pkg] -> { list of classes (inferred from list of .cxx files) }
  #
  srcfiles = {}
  
  # allfiles is a dictonary :
  # allfiles[libXXX.pkg] -> { list of all included files of that library }
  allfiles = {}
  
  for lib in libraries:
    if not re.search("da",lib):
      # handle the special case of DAs which are not part of libs, really
      srcfiles[lib] = getSourceFiles(lib)
    else:
      l = lib
      l = re.sub("lib","",l)
      l = re.sub("\.pkg","",l)
      srcfiles[lib] = [ l ]
    files = []
    for src in srcfiles[lib]:
      files.extend(getIncludeFiles(src))
    allfiles[lib] = unique(files)
  
  if debug:
    for lib in libraries:
      print lib
      for f in allfiles[lib]:
        l = findLibrary(f,srcfiles)
        print "   ",f,"(",l,")"
      print
      
  # deps is a dictionary
  # deps[libXXX.pkg] -> { list of libraries libXXX.pkg directly depends upon }
  
  deps = {}
  
  for lib,files in allfiles.items():
    d = []
    for f in files:
      l = findLibrary(f,srcfiles)
      if l != lib:
        d.append(l)
    deps[lib] = unique(d)
  
  if debug:
    for lib in deps:
      print lib, " depends on "
      for l in deps[lib]:
        print "   ",l
      print
        
  ofile = "%s.dot" % os.path.splitext(os.path.basename(sys.argv[0]))[0]
  
  f = open(ofile,"w")
  
  f.write("digraph G {")
  f.write("rankdir=BT;")

  for l,d in deps.items():
    for dl in d:
      if re.search("MUON",dl):
        f.write("node [shape=box,style=filled,fillcolor=yellow];")
      else:
        f.write("node [shape=ellipse,style=filled,fillcolor=lightgray];")
      f.write("%s -> %s;\n" %(shorten(l),shorten(dl)))
      
  f.write("}")
    
  print "You should now do :"
  print "tred %s > %s.bis" % ( ofile, ofile )
  print "dot -Tpng %s.bis > %s.png" % (ofile,ofile)
  
if __name__ == "__main__":
  main()
  