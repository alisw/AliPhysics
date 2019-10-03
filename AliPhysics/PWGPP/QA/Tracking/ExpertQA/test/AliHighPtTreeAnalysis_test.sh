

#
# Input data for unit test should be provided in standard way within ALICE
# Standard to be discussed in 1-2 weeks
# For the moment inputFile for the unit test hardwired
#
inputFile="/hera/alice/local/filtered/alice/data/2013/LHC13f/000197150/pass2/150/root_archive.zip#FilterEvents_Trees.root"
#inputFile="/hera/alice/local/filtered/alice/data/2013/LHC13f/000196702/pass2/150/root_archive.zip#FilterEvents_Trees.root";

#
# 1. Check code is running and producing list of histograms
#    a.) run the code
#    b.) check the logs
#    c.) in case of problems report and quit
aliroot -b -q $ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/makePlots.C\(\"$inputFile\"\) 2>&1 | tee experQA.log

#
# 2. Check some plots/invariants
#







