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

notassociatedfiles = {}

#_______________________________________________________________________________
def usage():
  """Describe usage of script
  """
  print "Usage: %s [-h | --help] [-d | --debug] [--da] directory_to_scan" % sys.argv[0]
  sys.exit(1)
  
#_______________________________________________________________________________
def append(list,a):
  """ append a to list, if a not there yet
  """
  if not a in list:
    list.append(a)
    
#_______________________________________________________________________________
def isempty(line):

  return len(line) < 2
  
#_______________________________________________________________________________
def iscomment(line):

  return re.search("^#",line)

#_______________________________________________________________________________
def iscontinuation(line):

  return line.find('\\') > 0

#_______________________________________________________________________________
def compactLines(ilines):
  """ Given an array of lines, remove empty lines, comment lines
  and concatenate lines that should be one (i.e. removing the \ continuation
  marks...
  """
  
  continuation = False

  olines = []
  
  currentline = ""
  
  i = 0
  
  for line in ilines:
    
    i = i + 1
    
    if line.find('\\') > 0:
      line = re.sub("\t|\n|\r|\\\\","",line)
      continuation = True
    else:
      continuation = False
      if isempty(currentline):
        currentline = line
      else:
        l = re.sub("\t","",line)
        currentline = currentline + re.sub(" {1,}"," ",l)
      if not iscomment(currentline) and not isempty(currentline):
        olines.append(re.sub(" {1,}"," ",currentline))
      currentline = ""
          
    if continuation:
      currentline = currentline + line
      
  return olines

#_______________________________________________________________________________
def tokenizeLines(lines):
  """
  Return a dict of keys -> value items.
  keys are the left part of lines supposed to be of the form :
  KEY:=VALUE
  or
  KEY+=VALUE
  or
  KEY=VALUE
  """
  
  tokens = {}
  
  define = ":="
  plus = "+="
  equal = "="
  
  separators = [ define, plus, equal ]
  
  for l in lines:
    sep = False
    for s in separators:
      if s in l:        
        a = l.split(s,1)
        key = a[0].strip()
        value = re.sub("\n|\t|\n","",a[1])
        if len(a) > 2:
          print "Something fishy here !"
        sep = True
        break
    if not sep:
      continue
    if not key in tokens.keys():
      tokens[key] = ""
    if s == plus:
      tokens[key] += value
    else:
      tokens[key] = value
      
  return tokens

#_______________________________________________________________________________
def variableSubstitution(tokenname,alltokens):
  """
  """

  value = alltokens.get(tokenname,"")
      
  for k in alltokens.keys():
    if re.search("\$\(%s\)"%k,value):
      # found something like $(VARIABLE), so expand that variable, 
      # by calling us again
      rep = value.replace("$(%s)"%k,variableSubstitution(k,alltokens))
      return rep
    if re.search("\$\(%s\:"%k,value):
      # found something like $(VARIABLE: 
      # we suppose it's then something like $(VARIABLE:.x=.y)
      # i.e. we replace x by y in VARIABLE's expansion
      t = variableSubstitution(k,alltokens)
      i1 = value.index(":")
      i2 = value.index(")")
      change = value[i1+1:i2].split("=")
      return t.replace(change[0],change[1])

  return value

#_______________________________________________________________________________
def patternSubstitution(value):

  if re.search("\$\(patsubst",value):
    rv = []
    # found the Makefile function $(patsubst %.x, %.y, list)
    i1 = value.index("(")
    i2 = value.rindex(")")
    a = value[i1+1:i2].split(",")
    source = a[0].replace("patsubst ","")
    destination = a[1].replace("%","")
    if source != "%":      
      print "Houston, we have a problem : ",value
      sys.exit(1)
    for l in a[2].split():
      rv.append(destination + l)
    return rv
    
  return value.split()

#_______________________________________________________________________________
def getSourceFiles2(lib,rootsys,alice_root):
  """Extract the list of files from a libXXX.pkg file
    Return a pair of list (sourceFiles,einclude), where einclude
    is the list of directories needed to be included compile the files.
  """
  
  try:
    f = open(lib)
  except:
    print "Cannot open file ", file
    sys.exit(1)

  filelines = f.readlines()
  
  f.close()
  
  lines = compactLines(filelines)
  
  tokens = tokenizeLines(lines)

  sourcesfiles = patternSubstitution(variableSubstitution("SRCS",tokens))
  eincludes = patternSubstitution(variableSubstitution("EINCLUDE",tokens))
  
  pkg = getLibPackage(lib)
  dir = os.path.join(alice_root,pkg)
  
  sourcesfiles = [ os.path.join(dir,x) for x in sourcesfiles ]
  
  return sourcesfiles,eincludes

#_______________________________________________________________________________
def getSourceFiles(lib,rootsys,alice_root):
  """Extract the list of files from a libXXX.pkg file
    Return a pair of list (sourceFiles,einclude), where einclude
    is the list of directories needed to be included compile the files.
  """
  
  # list of possible .pkg variables
  pkgkeys = [ "SRCS","EINCLUDE","HDRS","FSRCS","DHDR","CSRCS","CHDRS","ELIBS","EDEFINE","PACKFFLAGS","PACKCXXFLAGS","PACKCFLAGS","PACKSOFLAGS","EXPORT","EHDRS" ]
  
  keySRCS = pkgkeys[0]
  keyEINCLUDE = pkgkeys[1]
  
  sourcefiles = []
  pkg = getLibPackage(lib)
  einclude = [ "%s/include" % rootsys, "%s/STEER"  % alice_root, "%s/%s" % (alice_root,pkg) ]

  dir = os.path.dirname(lib)
  
  try:
    f = open(lib)
  except:
    print "getSourceFiles : could not open package file %s" % lib
    return sourcefiles, einclude
    
  src = False
  
  for line in f:
    l = line.strip()
    key = False
    for k in pkgkeys:
      if re.search(k,l):
        key = True
    if key:
      if re.search("^%s" % keySRCS,l):
        src = True
      else:
        src = False
      if re.search("^%s" % keyEINCLUDE,l):
        l = re.sub(keyEINCLUDE,' ',l)
        l = re.sub(':',' ',l)
        l = re.sub('=',' ',l)
        l = re.sub('\+',' ',l)
        a = l.split()
        for i in a:
          append(einclude,os.path.join(alice_root,i))
        
    if src:
      if re.search('Ali',l) and ( re.search('.cxx',l) or re.search('.h',l) ):
        l = re.sub(keySRCS,' ',l)
        l = re.sub(':',' ',l)
        l = re.sub('=',' ',l)
        l = re.sub('\+',' ',l)
        l = re.sub("\\\\",' ',l)
        for i in l.split():
          append(sourcefiles,os.path.join(dir,i))

  f.close()
  return sourcefiles,einclude
  
#_______________________________________________________________________________
def getIncludeFiles2(srcfile,alice_root,alice_target,rootsys):
  """Extract the list of included classes from a class, using the dep
   files generated in $ALICE_ROOT/package/tgt_ALICE_TARGET/*.d files
  It is much faster than getIncludeFile, as it reuses the output of 
  previously preprocessing part. Drawback is that it will only work
  on a compiled version of aliroot...
  """

  includes = []

  package = getFilePackage(srcfile,alice_root)
  file = re.sub("%s/%s" % (alice_root,package)," ",srcfile).strip()
  if file[0] == '/':
    file = file[1:]
  depfile = "%s/%s/tgt_%s/%s" % (alice_root,package,alice_target,file)
  depfile = re.sub("\.cxx",".d",depfile)
  
  try:
    f = open(depfile)
  except:
    print "Could not open file %s" % depfile
    print "From",srcfile
    return includes
    
  for line in f:
    line = line.strip()
    i = line.find(":")
    if i > 0:
      line = line[i+1:]
      parts = line.strip().split()
      for p in parts:
        if re.search(rootsys,p):
          p = rootsys
        else:
          if p[0] != '/':
            p = "%s/%s" % (alice_root,p)
          p = re.sub("%s/include" % alice_root,"%s/STEER" % alice_root,p)
        append(includes,p)

  f.close()
  
  return includes

#_______________________________________________________________________________
def getIncludeFiles(srcfile,eincludes,rootsys):
  """Extract the list of included classes from a class, using :
    gcc -MM srcfile -MG 
    and then parses the output...
    This version is quite slow as we're (re-)doing the preprocessing of
    all the files, but the advantage is that it'll work for a fresh checkout
    of aliroot, i.e. even before compilation
  """

  includes = []

#  try:
#    f = open(srcfile)
#  except:
#    print "Could not open file %s" % srcfile
#    return includes
#    
#  f.close()
  

  incdir = ""
  
  for i in eincludes:
    incdir = "%s -I%s" % (incdir,i)
    
  cmd = "gcc %s -MM %s -MG" % (incdir,srcfile)
  
  pre = os.popen(cmd)
  
  for line in pre:
    line = line.strip()
    line = re.sub("\\\\"," ",line)
    i = line.find(":")
    if i > 0:
      line = line[i+1:]
    line = line.strip()
    if len(line) > 0 and line != srcfile:
      if line.find('/') < 0:
        print "Got no path for file",srcfile," line=",line
        print "cmd was",cmd
      if re.search(rootsys,line):
        line = rootsys
      append(includes,line)
  pre.close()
          
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
def libshorten(libname):
  """From libYYYxxx.pkg to YYYxxx
  """
  
  s = os.path.basename(libname)
  if re.search("^lib",s):
    s = re.sub("^lib","",s)
    s = re.sub("\.pkg","",s)
    
  return s

#_______________________________________________________________________________
def fileshorten(file,path):
  """From path/toto/file to toto/file
  """
  
  s = re.sub(path," ",file).strip()
  if s[0] == '/':
    s = s[1:]
    
  return s

#_______________________________________________________________________________
def getFilePackage(file,alice_root):
  """ Get the package in which this file is defined
  """
  
  f = re.sub(alice_root,"/",file)
  while f[0] == '/':
    f = f[1:]
  p = f.split('/')
  return p[0]
  
#_______________________________________________________________________________
def getLibPackage(libname):
  """ Get the package in which this library is defined
  """
  
  p = libname.split('/')
  return p[len(p)-2]

#_______________________________________________________________________________
def tryRecover(f,inc2src,src2lib,alice_root):
  """ This method should try to recover the "father" of file f (most probably
  f is an include file
  The idea would be to find a cxx file that *directly* includes f, and take
  the lib of that cxx file as the source of f...
  Would that work ?
  Is it needed really ?
  """
  
  """
  print "tryRecover:",f
  
  if not f.find('\.h'): 
    return ""
    
  p = getFilePackage(f,alice_root)
  
  cxxfiles = inc2src.get(f,[])
  
  for file in cxxfiles:
    libs = src2lib.get(file,[])
  
    for l in libs:
      pl = getLibPackage(l)
      print f,file,p,pl
  """
  
  return ""
  
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
def main():

  # we cannot work w/o those environement variables, so check them...
  requiredVariables = [ "ROOTSYS", "ALICE_ROOT", "ALICE_TARGET" ]
  
  for r in requiredVariables:
    if not r in os.environ:
      print "%s is not defined. Cannot work." % r
      sys.exit(1)
      
  alice_root = os.environ.get("ALICE_ROOT")
  rootsys = os.environ.get("ROOTSYS")
  alice_target = os.environ.get("ALICE_TARGET")

  debug = 0
  noda = True
  
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hd",["help", "debug","da"])
  except getopt.GetoptError:
    print "Error in options"
    usage()

  for o, a in opts:
    if o in ( "-d","--debug" ):
      debug = debug + 1
    elif o in ( "-h","--help" ):
      usage()
      sys.exit()
    elif o == "--da":
      noda = False
    else:
      assert False, "unhandled option"

  dir = os.path.abspath(args[0])  
  dirs = []
  
  for sd in os.listdir(dir):
    ld = os.path.join(dir,sd)
    if os.path.isdir(ld) and not os.path.islink(ld):
      dirs.append(ld)

  dirs.append(dir)
    
#  requestedPackages = [ "MUON", "STEER", "RAW", "ITS", "TRD", "VZERO", "TPC", "PHOS", "TOF", "ZDC", "EMCAL", "HMPID", "SHUTTLE", "ACORDE", "HLT", "EVE" ];
  requestedPackages = [ "RAW", "STEER", "MUON", "HLT","EVE" ]
  
  # find the libraries defined in this directory (looking for libXXXX.pkg files)  
  libraries = []
  
  for d in dirs:
    for f in os.listdir(d):
      fulllib = os.path.join(d,f)
      p = getLibPackage(fulllib)
      if not p in requestedPackages:
        continue
      if re.search('^lib',f) and re.search('.pkg$',f):
        libraries.append(fulllib)
      if not noda and re.search('da.cxx',f) and not re.search('.svn',f):
        # append fake libraries for DAs
        tmp = re.sub("cxx","pkg",f)
        tmp = "lib%s" % tmp
        libraries.append(os.path.join(d,tmp))
      
  # from list of library files (libXXXyyy.pkg), try to find back the list of
  # "packages" = XXX
  packages = {}
  
  for l in libraries:
    p = getLibPackage(l)
    packages[p] = []
    
  for l in libraries:
    p = getLibPackage(l)
    packages[p].append(l)
    
  # src2inc[file.cxx] -> { all included files of that file }    
  src2inc = {}

  # inc2src[file.h] -> { all files that include that one }
  inc2src = {}
  
  # lib2src[libXXX.pkg] -> { list of files of that library }  
  lib2src = {}

  # src2lib[file.cxx] -> { list of libraries including that file }
  src2lib = {}
  
  # eincludes[libXXX.pkg] -> list of directories to be included to be able to compile the files
  eincludes = {} 
  
  # lib2inc[libXXX.pkg] -> { list of all included files of that library }
  lib2inc = {}
  
  # inc2lib[file.h] -> { list of libraries that include that file }
  inc2lib = {}
  
  for p in packages:

    print "Scanning ",p

    for lib in packages[p]:

      lib2inc[lib] = []
      
      print "  ",libshorten(lib),"..."

      if not re.search("da.pkg",lib):
        # handle the special case of DAs which are not part of libs, really
        lib2src[lib], eincludes[lib] = getSourceFiles2(lib,rootsys,alice_root)        
      else:
        l = lib
        l = re.sub("lib","",l)
        l = re.sub("\.pkg","",l)
        lib2src[lib] = [ "%s.cxx" % l ]
        eincludes[lib] = []
  
      files = []

      for src in lib2src[lib]:
#        inc = getIncludeFiles(src,eincludes[lib],rootsys)
        inc = getIncludeFiles2(src,alice_root,alice_target,rootsys)
        src2inc[src] = inc
        
        if not src in src2lib:
          src2lib[src] = []
          
        append(src2lib[src],lib)
        
        for i in inc:
          if not i in inc2src.keys():
            inc2src[i] = []
          append(inc2src[i],src)
          append(lib2inc[lib],i)
          if not i in inc2lib.keys():
            inc2lib[i] = []
          append(inc2lib[i],lib)
  
  # some debug at this point...
  
  if debug>=1:
    for lib in libraries:
      print lib," is made of "
      for f in lib2src[lib]:
        print "   ",fileshorten(f,alice_root)
      print " and includes "
      for h in lib2inc[lib]:
        print "   ",fileshorten(h,alice_root)
      if len(eincludes[lib]) > 0:
        print " and needs the following directories to be compiled "
        for f in eincludes[lib]:
          print "    ",f
      
  if debug>=2:
    print "src2lib relationship"
    for src,lib in src2lib.items():
      print fileshorten(src,alice_root),"(",
      for l in lib:
        print libshorten(l)," ",
      print ")"
      
  # now fills the ultimate array, lib2lib
  # lib2lib[libXXX.pkg] -> { libYYY.pkg }, list of libraries that lib depends on
  
  lib2lib = {}
  
  for lib in libraries:

    lib2lib[lib] = []

    for hfile in lib2inc[lib]:
      
      l = "external"
      
      # start simple : is f contains ROOTSYS, it's ROOT.
      if re.search(rootsys,hfile):
        l = "ROOT"
      else:
        # not that simple, let's try to find out...
        cxx = re.sub("\.h",".cxx",hfile)
        dl = src2lib.get(cxx,[])
        if len(dl)==1:
          l = dl[0]
        elif len(dl)>1:
          print "Got several libs(",len(dl),"for ",hfile,":"
          print dl
          
      if l =="external":
        notassociatedfiles[hfile] = 1
      else:
        append(lib2lib[lib],l)

  ###################### Debug parts...
  if debug>=1:
    for lib in libraries:
      print libshorten(lib),"depends on"
      for f in lib2lib[lib]:
        print "   ",libshorten(f)
  
  if debug>=2:
      print 
      print "From source files to include files : "
      for cxxfile, hfile in src2inc.items():
        print fileshorten(cxxfile,alice_root)
        for h in hfile:
          print "    ",fileshorten(h,alice_root)

  if debug>=3:
      print
      print "From include files to source files : "
      for i,sources in inc2src.items():
        print fileshorten(i,alice_root), len(sources)
        for s in sources:
          print "   ",fileshorten(s,alice_root),"(",
          for l in src2lib[s]:
            print libshorten(l),
          print ")"
  ###################### 
  
  if len(notassociatedfiles) > 0:
    print "The following files could not be associated with any library..."
    for f in notassociatedfiles.keys():
      print f,
      t=tryRecover(f,inc2src,src2lib,alice_root)
      print t
      
  # output the dot file that will have to be processed by the dot program
  
  ofile = "%s.dot" % os.path.splitext(os.path.basename(sys.argv[0]))[0]
  
  f = open(ofile,"w")
  
  f.write("digraph G {\n")
  f.write("rankdir=BT;\n")

  defaultcolor = "lightblue"
  
  colors = {}
  
  colors["MUON"] = "lightyellow"
  colors["STEER"] = "lightgray"
  
  for l,d in lib2lib.items():
    for dl in d:
      f.write("%s->%s;\n" %(libshorten(l),libshorten(dl)))
      
  for p in packages:
    f.write("subgraph cluster_%s {\n" % p.lower())
    color = colors.get(p,defaultcolor)      
    f.write("style=filled;\n")
    f.write("color=%s;\n" % color)
    f.write('label="%s";\n' % p)
    for lib in packages[p]:
      f.write("%s\n" % libshorten(lib))
    f.write("}\n")
    
  f.write("}\n")
    
  print "You should now do :"
  print "tred %s > %s.bis" % ( ofile, ofile )
  print "dot -Tpng %s.bis -o %s.png" % (ofile,ofile)
  
if __name__ == "__main__":
  main()
  