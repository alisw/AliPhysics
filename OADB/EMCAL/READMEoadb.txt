/*! \page READMEoadb Accessing the EMCAL OADB from EOS

The EMCal OADB files have been moved to the EOS directory /eos/experiment/alice/analysis-data/OADB/EMCAL which is accessible via lxplus.
Analysis tasks can access these files via, for example: ```AliDataFile::GetFileNameOADB("EMCAL/EMCALTimeL1PhaseCalib.root")``` which has been implemented into the Tender and the Correction Framework as well as most analysis tasks in AliPhysics.
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

In addition, a short history of changes to the files in EOS will be listed here:

- 20180213: Moved all files to EOS.
- 20180220: Update of EMCALBadChannels.root with LHC17o BC maps
- 20180305: Update of EMCALBadChannels.root with LHC17m BC maps
- 20180308: Update of EMCALBadChannels.root with LHC17ghjk BC maps
- 20180308: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with LHC17k,m,h,o,p,q time calibrations
- 20180311: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with LHC17o,p,q time calibrations
- 20180320: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with LHC17j time calibrations
- 20180403: Update of EMCALBadChannels.root with LHC16f BC maps
- 20180406: !! Fix for AliOADBContainer deletion problem uploaded, all files recreated for this !!
- 20180406: Update of EMCALBadChannels.root with new maps for LHC17i,l and additional bad channels for LHC16h,i,j,l,q,t
- 20180423: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with LHC17i time calibrations
- 20180504: Update of EMCALBadChannels.root with additional bad channels for LHC17n
- 20180524: Update of EMCALBadChannels.root with new maps for LHC17r
- 20180530: Update of EMCALTimeCalib.root with updated LHC17m and LHC17k time calibrations
- 20180531: Update of EMCALTimeL1PhaseCalib.root with two missing runs from LHC16kv time calibration
- 20180606: Update of EMCALBadChannels.root with new maps for LHC17c,f
- 20180612: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with LHC17l time calibrations, as well as final calibrations for LHC16i and LHC16j
- 20180615: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with increased run range for LHC17l
- 20180618: Update of EMCALBadChannels.root with runwise maps for LHC18d (except for runs 286129 286350 286341 286313)
- 20180621: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with LHC16o and LHC16p time calibrations
- 20180622: Update of EMCALBadChannels.root with additional cells for LHC16j and LHC16l and a fix for the LHC16h maps which were not properly added in OADB before
- 20180706: Update of EMCALBadChannels.root with additional cells for LHC16p and update of EMCALTimeCalib.root with updated calibrations for LHC17o
- 20180711: Update of EMCALBadChannels.root with new maps for LHC18d
- 20180711: Update of EMCALBadChannels.root with new maps for LHC17d, LHC18e, LHC18f as well as updates to the maps of LHC17m and LHC17o
- 20180713: Update of EMCALTimeCalib.root with fixed run range for period LHC16i
- 20180713: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with LHC18d time calibrations
- 20180718: Update of EMCALBadChannels.root with new maps for LHC18b and LHC18c including special maps for some runs
- 20180730: Update of EMCALBadChannels.root with new maps for LHC16f (nominal and lowB) as well as LHC16g
- 20180713: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with updated LHC16h time calibrations

*/
