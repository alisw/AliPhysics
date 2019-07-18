The TPC dEdx recalibration files are stored in the EOS directory /eos/experiment/alice/analysis-data/OADB/PWGGA which is accessible via lxplus.
Analysis tasks can access these files via, for example: ```AliDataFile::GetFileNameOADB("PWGGA/TPCdEdxRecalibOADB.root")``` which has been implemented into AliConversionPhotonCuts and the corresponding analysis tasks in AliPhysics.
If you want to have the OADB files locally, you can download them from lxplus via:
~~~{.sh}
rsync -av --delete cern_user@lxplus.cern.ch:/eos/experiment/alice/analysis-data/ /path/to/my/local/oadb/
~~~

In order for local tests to work properly, please add the ALICE_DATA global variable to your bashrc or similar.
~~~{.sh}
ALICE_DATA=/path/to/my/local/oadb
~~~

Furthermore, the "export" command should be used in addition to adding the path to the bashrc. This will make the variable available to all processes:
~~~{.sh}
export ALICE_DATA=/path/to/my/local/oadb
~~~

The code for the creation of the OADB files as well as for bookkeeping of input files can be found in the following repository: https://gitlab.cern.ch/alice-pcg/AnalysisSoftware (folder TaskQA/UpdatedEdxRecalib_OADB.C)

In addition, a short history of changes to the files in EOS will be listed here:

- 20180524: Created first file on EOS containing calibrations for LHC13b-f

*/
