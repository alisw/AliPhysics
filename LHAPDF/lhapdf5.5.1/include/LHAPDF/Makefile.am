EXTRA_DIST = FortranWrappers.h.in

if ENABLE_CCWRAP

PKGincludedir = $(includedir)/LHAPDF
PKGinclude_HEADERS = LHAPDF.h LHAPDFfw.h
nodist_PKGinclude_HEADERS = FortranWrappers.h
if ENABLE_OLDCCWRAP
PKGinclude_HEADERS += LHAPDFWrap.h
endif

else 

## Make sure that headers are always bundled
EXTRA_DIST += LHAPDF.h LHAPDFfw.h LHAPDFWrap.h

endif
