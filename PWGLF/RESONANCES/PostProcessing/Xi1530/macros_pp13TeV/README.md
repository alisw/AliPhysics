# Postprocessing macro for Xi(1530)0 in pp 13TeV #

## Files ##

 - DrawXi1530.C (core)
 - PlotXi1530.C (default drawing macro)
 - grid.py (grid macro for systematics)
   - extract.jdl
   - extract.sh
 - Xi1530Systematics.C (systematic)
 - YieldMean.C
 - DrawXi1530forXi1820.C

## How to use ##

Before start, all the outputs from train should be in the same folder:
 - AnalysisResults_LHC15fi_16deghijklop_17cefgijklmor_SYS_light_fixed.root (merged output from LF_pp#1201, #1200)
 - AnalysisResults_Xi1530LHC18c6a_RsnMC_Final_part1.root (merged output from LF_pp_MC#1074)
 - AnalysisResults_Xi1530LHC18c6a_RsnMC_Final_part2.root
 - AnalysisResults_Xi1530LHC18c6b_RsnMC_Final.root
 - AnalysisResults_Xi1530LHC18c6b4_RsnMC_Final.root
 - AnalysisResults_Xi1530LHC18c6c_RsnMC_Final.root
 - AnalysisResults_Xi1530LHC16_GenMC_final.root (merged output from LF_pp_MC#1094, #1099)
 
For the RsnMC files, they failed to merge because of the data limit in one bin(true input)

### DrawXi1530.C ###
Core macro to extract the signal, background, efficiency, corrected spectra.

Arguments:
 - $1: integer for the cut systematic variation. eg. 1 = default, 2= "TPCNsigmaXi1530PionLoose" and so on..
 - $2: multiplicity bin start eg. 0
 - $3: multiplicity bin end eg. 0
 - $4: Option parameter eg. "Default" see availble options below
 - $5: sub option parameter eg. 1

Basic usage:
```
root -l -b -q DrawXi1530.C
# will draw with default option: 0-100% bin
```

Advanced usage:
```
root -l -b -q DrawXi1530.C\(1,0,10,\"BinCount\",1\)
# will draw with with "BinCount" option in 0-10% bin

root -l -b -q DrawXi1530.C\(1,0,100,\"FitVarLm\",3\)
# will draw with with "Fit range variation to left -3" in 0-100% bin

# Pararell processing
cat pararell.txt
-l -b -q DrawXi1530.C\(1,0,100,\"FitVarLm\",3\)
-l -b -q DrawXi1530.C\(1,0,100,\"FitVarLm\",2\)
-l -b -q DrawXi1530.C\(1,0,100,\"FitVarLm\",1\)
-l -b -q DrawXi1530.C\(1,0,100,\"FitVarLp\",3\)
-l -b -q DrawXi1530.C\(1,0,100,\"FitVarLp\",2\)
-l -b -q DrawXi1530.C\(1,0,100,\"FitVarLp\",1\)

cat parallel.txt | xargs -P6 -L 1 root
# it will runs every input line in the pararell.txt file in a sametime.
```

#### Available options ####
(n) means it takes suboption value

 - BinCount
 - LikeSignBkg
 - BkgFit
 - MCcheck
 - FitVarLm (n)
 - FitVarLp (n)
 - FitVarRp (n)
 - FitVarRm (n)
 - FitVarBothm (n)
 - FitVarBothp (n)
 - NormVarm (n)
 - NormVarp (n)
 - NormVarLp (n)
 - NormVarLm (n)
 - NormVarRp (n)
 - NormVarRm (n)
 - BinCountLm (n)
 - BinCountLp (n)
 - BinCountRm (n)
 - BinCountRp (n)
 - BinCountBothm (n)
 - BinCountBothp (n)
 - inel (under development)
 

### PlotXi1530.C ###
Ploting macro for the output of DrawXi1530.C

Argument:
filename

Basic usage:
```
root -l -b -q PlotXi1530.C\(\"AnalysisResults_Extracted_1_Multi_0.00-100.00_Default1.root\"\)
# will make a folder and save the deafult outputs.
```

### grid.py ###
Grid macro for the run full systematic variation.
Before use, need to update the path and value in the grid.py file.

Basic usage:
```
python grid.py submit
# will copy the file default files and submit using pararell processing

python grid.py download
# will download the finished files in the grid. if something missing, it will ask for resubmit

python grid.py check
# will check if all the files are in the folder. if something missing, it will ask for resubmit

python grid.py local
# will try to run all remains task in local.
```

### Xi1530Systematics.C ###
Prototype systematic macro
(will be updated to separtated class, "SystematicHelper" later)

It needs all the outputs from grid.py in ./data folder

Basic usage:
```
root -l -b -q Xi1530Systematics.C
```
