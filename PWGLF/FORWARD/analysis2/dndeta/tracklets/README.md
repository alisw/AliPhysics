# SPD tracklet dNch/deta analysis

Code in this directory is from Roberto and Valentina's work in Pb-Pb
at 5.02TeV.  The code has been slightly modified so that we can select
the eta binning more easily. 

This code is put here mainly for back-up purposes.  Once there's a
full official task, we can remove this sub-directory. 

## Content:

### Task code: 

* `AliTrackletTaskMulti.cxx`
* `AliTrackletTaskMulti.h`
   Task code 

* `FixPaths.C`
  A small piece of code to fix include paths in AliEn 

* `run.sh`
  Run the task on AliEn.  Edit for selecting some other runs 

### Weights:

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

### Original correct and associated scripts:

* `SaveCanvas.C`
* `SaveCanvas.h`
* `CorrectSpectraMultiMCBG.C`
  Modified original script from Ruben.  The script does the final
  calculations. 

* `Extract.C`
  Extract results into GSE's

* `Final.C`
  Run all the final steps 
 
### Alternative correct and associated scripts:

* `SimpleCorrect.C`
  A simplified version of the above 

* `Correct.C` 
  Simple interface to `SimpleCorrect.C` that inputs and outputs data
  from and to separate directories. 
  
* `DoOne.C` 
  Script to run corrections for left, middle, and right vertex bins. 

* `Compare.C` 
  Compares results using various weights to no weights 
  
* `DrawDeltas.C` 
  (Deprecated) script to draw Delta distributions 
  
* `DrawDeltas2.C` 
  Script to draw a single centrality bin Delta distribution 
  
* `DrawKs.C` 
  Script to draw background scaling factors 
  
<!-- EOF -->

