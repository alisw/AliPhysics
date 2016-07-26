SPD tracklet dNch/deta analysis
===============================

Code in this directory is from Roberto and Valentina's work in Pb-Pb
at 5.02TeV.  The code has been slightly modified so that we can select
the eta binning more easily. 

This code is put here mainly for back-up purposes.  Once there's a
full official task, we can remove this sub-directory. 

Content:
--------

* `AliTrackletTaskMulti.cxx`
* `AliTrackletTaskMulti.h`
   Task code 
   
* `SaveCanvas.C`
* `SaveCanvas.h`
* `CorrectSpectraMultiMCBG.C`
  Modified original script from Ruben.  The script does the final
  calculations. 

* `SimpleCorrect.C`
  A simplified version of the above 
  
* `Extract.C`
  Extract results into GSE's

* `Final.C`
  Run all the final steps 
 
* `FixPaths.C`
  A small piece of code to fix include paths in AliEn 

* `REWEIGHTpid_ka-.root`
* `REWEIGHTpid_ka+.root`
* `REWEIGHTpid_pi-.root`
* `REWEIGHTpid_pi+.root`
* `REWEIGHTpid_pr-.root`
* `REWEIGHTpid_pr+.root`
* `REWEIGHTpid.root`
* `REWEIGHTpt.root`
* `REWEIGHTstr-.root`
* `REWEIGHTstr.root`
* `REWEIGHTstr+.root`
  MC reweighs 

* `run.sh`
  Run the task on AliEn.  Edit for selecting some other runs 

