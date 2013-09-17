#!/bin/bash

#                                 Save, MCcase, IncludeEWfromTherm, SC,      G,      EW,   GRS, EDBin,  CH,   Mbin,  Ktbin
# without G
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     1)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     2)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     3)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     4)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     5)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     6)'
# with G
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     1)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     2)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     3)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     4)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     5)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     6)'
#
# Include EW from Therminator fits
# without G
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     1)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     2)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     3)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     4)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     5)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kFALSE, kFALSE, kTRUE,     0,  +1,    '$1',     6)'
# with G
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     1)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     2)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     3)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     4)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     5)'
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kTRUE,   kTRUE,  kTRUE, kFALSE, kTRUE,     0,  +1,    '$1',     6)'
#
# 3-particle
# EDbin 1
# GRS
# SC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kTRUE,    0,  +1,    '$1',     10)'
# SC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kTRUE,    0,  -1,    '$1',     10)'
# MC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kTRUE,    0,  +1,    '$1',     10)'
# MC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kTRUE,    0,  -1,    '$1',     10)'
# Omega0
# SC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kFALSE,    0,  +1,    '$1',     10)'
# SC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kFALSE,    0,  -1,    '$1',     10)'
# MC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kFALSE,    0,  +1,    '$1',     10)'
# MC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kFALSE,    0,  -1,    '$1',     10)'
#
# EDbin 2
# GRS
# SC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kTRUE,    1,  +1,    '$1',     10)'
# SC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kTRUE,    1,  -1,    '$1',     10)'
# MC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kTRUE,    1,  +1,    '$1',     10)'
# MC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kTRUE,    1,  -1,    '$1',     10)'
# Omega0
# SC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kFALSE,    1,  +1,    '$1',     10)'
# SC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kTRUE,  kFALSE, kTRUE, kFALSE,    1,  -1,    '$1',     10)'
# MC EW, +
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kFALSE,    1,  +1,    '$1',     10)'
# MC EW, -
root -b -q 'Plot_PDCumulants.C++(kTRUE, kFALSE,        kFALSE,   kFALSE,  kFALSE, kTRUE, kFALSE,    1,  -1,    '$1',     10)'