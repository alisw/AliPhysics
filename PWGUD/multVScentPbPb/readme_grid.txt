# Example of running the job on the grid w/o rebuilding the tracklets (hence
# only MC bg can be used)
# 

mode="full"

# Pb-p
aliroot -b -q MyAnalysisTaskTrackletMultiGRID.C\(\"${mode}\",196433,-1,\"/alice/sim/2014/LHC14i2\"\)
aliroot -b -q MyAnalysisTaskTrackletMultiGRID.C\(\"${mode}\",196433,-1,\"/alice/data/2013/LHC13f\",\"pass2/*ESDs.root\"\)

The merged outputs of the jobs (e.g. trmult.root) from data and MC should be fed to postprocessing macro
aliroot -b -q CorrectSpectraMultiMCBG.C(char* <dataoutput>,char* <mcoutput>,<somePrefix>);

