=========================================================
 Welcome to the tutorials for PWG2/RESONANCES package!!!
=========================================================

1. PRELIMINARY STEP

If you want to be sure that you compiled the most recent PAR libraries,
and clear any eventually resident older version, you should launch a ROOT
job like this:

> root -q -b prepareTutorial.C

which will clear all PARs and re-compile the ones in the working directory.
Of course, make sure that you are in the right working directory and that
the PAR files you have there are up-to-date!!!

2. TUTORIALS

The tutorial is projected to be run in a PROOF session on CAF,
so you have to make sure that you have all necessary rights to use it,
and have obtained an AliEn token, as recently required.

All tutorials are launched using the main front-end macro:

> root -q -b runMyProcess.C

This macro recalls the other macro named "myRsnAnalysis.C" which initializes
all required tasks and the stuff which will produce the output histograms.
To switch from one tutorial to the other, you must do just a simple change:

1) edit macro "myRsnAnalysis.C"
2) set default argument to the name of the tutorial you want to run
3) save the macro
4) launch "runMyProcess.C"

ENJOY!