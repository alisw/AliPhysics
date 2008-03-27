c --- these define the lengths of the path string and
c --- full absolute file names.  If these are too short
c --- they can be increased as desired.  lhaname should be
c --- 40 larger than lhapth to accomodate the larges grid
c --- file name.  
c --- include 'parsetup' must be done first to get a value
c ---  for nmxset.
c ---  do a "make clean" before "make" to force a recompilation
c --- of everything
c
      character*232 lhapath
      character*272 lhaname,lhanames(nmxset)
