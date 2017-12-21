/*! \page READMEjetQA Jet QA

## Running the Jet QA task

The task AliAnalysisTaskPWGJEQA performs basic QA relevant to the jet group. This task is appropriate as a first step for an analysis (as opposed to light QA to be run shortly after data-taking -- for this see the [AnalysisQA train](https://alice.its.cern.ch/jira/browse/PWGPP-281)). See the task for details.

- For an example configuration in PbPb data, see Jets_EMC_PbPb train 1872.
- For an example in a pp Pt-hard production, see Jets_EMC_pp_MC train 1031.

## Plotting Jet QA

There are a handful of scripts used to process the output from the QA task, covering several different use cases:
  - runJetQA.sh
  - plotPWGJEQA.py
  - plotPowerpoint.py

For Pt-hard productions, there are two additional scripts:
  - downloadPtHard.sh
  - scalePtHardHistos.py

The scripts are located at $ALICE_PHYSICS/PWGJE//EMCALJetTasks/macros/JetQA.

Pre-requisites:
  - The script plotPowerpoint.py uses the python package python-pptx. You can install this (via pip) with: pip install python-pptx
  - In order for error bars to be computed correctly, you should use ROOT 6. There is an issue with ROOT 5 versions used by ALICE that calling Sumw2() on a THnSparse after filling does not propagate the errors -- a feature that is extensively used in this QA machinery.

Please read on for simple instructions of the workflow.

### Run-by-run QA

The main script is runJetQA.sh. It downloads the output files from your train, generates plots from these roots files, and constructs a powerpoint presentation. See the script for options that need to be configured.

The script uses two python scripts:
  - plotPWGJEQA.py is a plotting macro that generates plots in an output directory
  - plotPowerpoint.py is a macro that generates a powerpoint presentation using these plots

The run-by-run QA machinery also expects a reference file consisting of all runs merged together. This file should be placed in the location specified in runJetQA.sh.

#### Steps for data or non-Pt-hard productions

1. Download the merged output file from your train (merged over all runs), and place it as the reference file. Run the script plotPWGJEQA.py to generate QA plots for this reference file, and place it in a directory TrainOutput/AllRuns/QAoutput.
2. Configure and execute runJetQA.sh.

#### Steps for Pt-hard productions

Rather than using runJetQA.sh to download data, you should instead use the script downloadPtHard.sh. This script downloads data per Pt-hard bin, and merges and scales appropriately using Pt-hard weights computed from histograms filled in AliAnalysisTaskPWGJEQA, using the script scalePtHardHistos.py.

1. Configure downloadPtHard.sh to download data (but not merge or scale), and execute. This may take awhile since it is a lot of data! You should then copy the data -- we will use one copy to generate the reference file, and one for the run-by-run QA. You may additionally wish to make an additional backup, in case something is misconfigured in later steps, so that you don't have to re-download.
2. Configure downloadPtHard.sh to merge/scale/merge (but not download data) and with doRunByRun=false. This will produce your reference file -- place it appropriately. Run the script plotPWGJEQA.py to generate QA plots for this reference file, and place it in a directory TrainOutput/AllRuns/QAoutput.
3. Configure downloadPtHard.sh to merge/scale/merge (but not download data) and with doRunByRun=true. You should also enable the option useReferenceFile=true, in order to take Pt-hard weights from the reference file. This will produce your run-by-run QA root files.
4. Configure runJetQA.sh to plot QA and generate presentation, but not to download data. Execute the script.

### Non Run-by-run QA

If you do not need run-by-run QA, you can download the the output by hand and directly execute the plotting macro plotPWGJEQA.py. See the script for details.

In case of Pt-hard production, you should use downloadPtHard.sh as described above (but with doRunByRun = false) to download, merge, and scale the output. Then proceed to execute plotPWGJEQA.py.

## Questions / Comments ?

If you encounter problems, or if anything is confusing, or if you would like to make improvements (which are very welcome!) please write to <james.mulligan@yale.edu> or post to [EMCAL-24](https://alice.its.cern.ch/jira/browse/EMCAL-24).

*/
