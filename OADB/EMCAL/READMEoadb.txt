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

The code for the creation of the OADB files as well as for bookkeeping of input files can be found in the following repository: https://gitlab.cern.ch/alice-EMC/EMCCalib

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
- 20180821: Update of EMCALBadChannels.root with additional cells for LHC17l,m,o
- 20180828: Update of EMCALTemperatureCalibSM.root and EMCALTemperatureCalibParam.root containing new temperature calibration parameters for Run1+Run2 (only run2 active in OADB)
- 20180910: Update of EMCalBadChannels.root with additional cells for LHC16f, LHC16i, LHC17g and LHC18b
- 20180911: Recreated and updated EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with updated LHC16kl and new LHC17r time calibrations
- 20180924: Update of EMCALBadChannels.root with additional cells for LHC16k and LHC16o
- 20181011: Update of EMCALBadChannels.root with additional cells for LHC17h and LHC17i
- 20181106: Update of EMCALRecalib.root with the change: "from run 295275 no calibration is applied to SMs 17/18 since the online voltages were updated"
- 20181013: Update of EMCALBadChannels.root with new maps for LHC18g, LHC18h, LHC18i
- 20181122: Update of EMCALBadChannels.root with new maps for LHC18m
- 20181129: Update of EMCALBadChannels.root with additional cells for LHC17f
- 20181206: Update of EMCALBadChannels.root with new maps for LHC18n and LHC18o
- 20181207: Update of EMCALBadChannels.root with VERY PRELIMINARY maps for LHC18p,q,r
- 20190118: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for LHC17c,d,f,g + fix for LHC17n
- 20190124: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for LHC16g, LHC18f, LHC18g
- 20190129: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for LHC18m
- 20190206: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for LHC18e
- 20190208: Update of EMCALBadChannels.root with new maps for LHC18p and upaded maps for DCal for LHC17h,i,j,k,l,m,o,r
- 20190212: Update of EMCALBadChannels.root with upaded maps for LHC16h,i,j,k,l,q,t
- 20190214: Update of EMCALBadChannels.root with upaded maps for LHC16g,o,p and LHC17c
- 20190215: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for LHC18b/c and LHC18n/o
- 20190301: Update of EMCALBadChannels.root with upaded maps for LHC15o
- 20190308: Update of EMCALBadChannels.root with new maps for LHC18j,k,l and adjusted range of LHC18m(V6) map
- 20190320: Update of EMCALBadChannels.root with new maps for LHC18d that cover larger run ranges and LHC18b,e,f with additional bad channels
- 20190321: Update of EMCALBadChannels.root with additional bad channels for LHC18h,i
- 20190329: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for LHC18j,k,l,p and updates for LHC18m,o
- 20190410: Update of EMCALBadChannels.root with additional bad channels for LHC18b,e,f,h,i,m as well as updated time calibrations with default 600ns shift for uncalibrated cells in LHC18g,n
- 20190411: Update of EMCALBadChannels.root with additional bad channels for LHC18o
- 20190412: Update of EMCALBadChannels.root with new maps for LHC18q and LHC18r
- 20190415: Update of EMCALBadChannels.root with additional bad channels for LHC18c (low B)
- 20190417: Update of EMCALBadChannels.root with additional bad channels for LHC18p
- 20190423: Update of EMCALBadChannels.root with additional bad channels for LHC18d,j,k,l
- 20190429: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for LHC18q,r
- 20190513: Update of EMCALBadChannels.root with remaining run range for LHC18r bad channels
- 20190520: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with calibs for remaining runs of LHC18r
- 20190528: Update of EMCALBadChannels.root with analysis level QA bad channels for LHC18q,r
- 20190612: Update of EMCALBadChannels.root with analysis level QA bad channels for LHC17n (XeXe)
- 20190612: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with updated calibs for LHC17g, LHC18g and LHC18n (lower thresholds)
- 20190617: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with updated calibs for LHC15o (with BC maps and good runs)
- 20190626: Update of EMCALBadChannels.root with analysis level QA bad channels for LHC18q and r (second update)
- 20190716: Update of EMCALBadChannels.root with analysis level QA bad channels for LHC18g (split in three ranges) and LHC18n
- 20190722: Update of EMCALTimeCalib.root and EMCALTimeL1PhaseCalib.root with updated calibs for LHC17f and new calibs for LHC16f_lowB
- 20190722: Update of EMCALBadChannels.root with analysis level QA bad channels for LHC17g and LHC18c
- 20190830: Update of EMCALTemperatureCalibSM.root and EMCALTemperatureCalibParam.root containing new temperature calibration parameters for full Run1 and Run2
- 20191030: Update of EMCALTimeL1PhaseCalib.root with 2016 PAR calibrations
- 20200316: Re-addition of 2016-2018 PAR calibrations in EMCALTimeL1PhaseCalib.root + missing PAR calib for three runs of LHC18q
- 20200414: Update of EMCALTimeL1PhaseCalib.root with additional 2016 and 2018 PAR calibrations

*/
