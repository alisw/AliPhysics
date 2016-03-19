// $Id$

/*! 

\page README_shuttle MUON Shuttle

How to test the Shuttle preprocessor(s) for MUON.

We will get two "logical" MUON preprocessors : one for the tracker and one for the trigger.
Both will manage several subtasks (e.g. the tracker one will handle pedestals,
  and occupancy for instance, while the trigger one will handle masks and trigger lut)
"Physically", only one class will manage both the tracker and the trigger : AliMUONPreprocessor.
Depending on the subsystem and on the task to be performed (based on the run type), this class
 will instanciate the correct set of AliMUONVSubProcessor(s) which does the actual job.
Output of most processors will end up in OCDB (Offine Condition DataBase). A set of helper functions
 to peek at this OCDB are gathered in AiMUONCDB class.
 

\section shuttle_s1 TestMUONPreprocessor.C

This is the master macro used to check the MUON part of the Shuttle.
Depending on what you want to test, you'll have to modify the input files 
(using shuttle->AddInputFile) and/or the run type (using shuttle->AddInputRunParameter())


\section shuttle_s2 AliMUONPreprocessor(const TString& detName)

Depending on how this one is constructed, and depending on the runtype, it will
 perform differents tasks. Note that for the moment the runtypes are "fake", i.e.
 put by hand in the TestMUONPreprocessor.C macro, and might not correspond to
 the final values to be used by the DAQ.

<pre> 
detName   runType                     task to be done           worker class (AliMUONVSubprocessor child)
--------------------------------------------------------------------------------------------------------
MCH       PEDESTAL                    read ASCII ped files      AliMUONPedestalSubprocessor
                                      and put them into OCDB
                        
MCH       GMS                         read GMS alignment files  AliMUONGMSSubprocessor
                                      and put them into OCDB

MCH       PHYSICS                     read DCS HV values and    AliMUONHVSubprocessor
                                      put them into OCDB

MTR       CALIBRATION		      read date files (masks,   AliMUONTriggerSubprocessor
	  			      crates, and LUT)
				      and put them in OCDB           
                                      
MTR       PHYSICS                     read date files as in     AliMUONTriggerDCSSubprocessor
	                              CALIBRATION,
				      read DCS HV and currents
				      and put them in OCDB
				     
</pre>


\section shuttle_s3 Pedestals

Two options here. You can either use a pre-done set of ASCII pedestals files (generated as
 explained below for the 2nd option), located at /afs/cern.ch/user/l/laphecet/public/LDC*.ped, 
 or build you own set.
 
We've written an AliMUONPedestalEventGenerator which creates fake pedestal events. The pedestal values
are taken from the Offline Condition DataBase (OCDB) (which is itself fakely filled
using the static WritePedestals() method of AliMUONCDB class

So first generate a valid pedestal CDB entry by using the AliMUONCDB class. There's one
 little trick : you should first point to the "default" OCDB (local://$ALICE_ROOT/OCDB) in
 order to get the mapping loaded, then only you can play with another OCDB. 
 Or, alternatively, you can put the mapping stuff in the test OCDB, like this :
 
<pre>
root[] AliMpDDLStore::ReadData(); // read mapping from ASCII files
root[] const char* cdbpath="local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB"; // where to put the CDB
root[] AliCDBManager::Instance()->SetDefaultStorage(cdbpath);
root[] AliMpCDB::WriteMpSegmentation();
root[] AliMpCDB::WriteDDLStore();
</pre>

If you've not put the mapping in the test database, then you must start with the default OCDB, load the mapping, and then only switch to the 
 test database :
 
<pre>
root[] AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB"); // only if you've not put the mapping in test OCDB
root[] AliCDBManager::Instance()->SetRun(0); // only if you've not put the mapping in test OCDB
root[] AliMpCDB::LoadDDLStore(); // only if you've not put the mapping in test OCDB
// below are lines to be executed whatever you did with the mapping...
root[] const char* cdbpath="local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB"; // where to put the CDB
root[] AliCDBManager::Instance()->SetDefaultStorage(cdbpath);
root[] Bool_t defaultValues = kFALSE; // to generate random values instead of plain zeros...
root[] Int_t startRun = 80;
root[] Int_t endRun = 80;
root[] AliMUONCDB::WritePedestals(defaultValues, startRun, endRun);
</pre>

Expected output is (don't worry about the warnings, they are harmless) :

<pre>
I-AliMUONCDB::ManuList: Generating ManuList...
I-AliMUONCDB::ManuList: Manu List generated.
I-AliMUONCDB::MakePedestalStore: 16828 Manus and 1064008 channels.
I-AliMUONCDB::WritePedestals: Ngenerated = 1064008
I-AliCDBManager::Init: AliEn classes enabled in Root. AliCDBGrid factory registered.
I-AliCDBManager::SetDefaultStorage: Setting Default storage to: local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestCDB
I-AliCDBLocal::PutEntry: CDB object stored into file ($ALICE_ROOT)/SHUTTLE/TestShuttle/TestCDB/MUON/Calib/Pedestals/Run80_80_v0_s0.root
</pre>

Then use the AliMUONPedestalEventGenerator to produce simulated pedestal events.

Usage (from the Root prompt) :
<pre>
AliMpCDB::LoadDDLStore2(); // load mapping from "default" OCDB=local://$ALICE_ROOT/OCDB
AliCDBManager::Instance()->SetDefaultStorage(cdbpath); // so you will read 
// back pedestals values generated in the previous step
const char* dateFileName = "raw.date"; // base filename for the output
Int_t runNumber = 80; // run number used to fetch the pedestals from the OCDB 
Int_t nevents = 100; // # of events to generate. 100 should be enough
gSystem->Load("libMUONshuttle"); // needed or not depending on whether it's already loaded or not
AliMUONPedestalEventGenerator ped(runNumber,nevents,dateFileName);
ped.Exec("");
</pre>

It *will* take a lot of time (mainly due to the fact that we're writing a 
bunch of ASCII files = DDL files), so please be patient.

The output should be the normal simulated sequence of MUON.Hits.root, MUON.SDigits.root,
 MUON.Digits.root, raw/*.ddl files and raw.date.LDCi where i=0-3 (i.e. one DATE file
per LDC, as will be used in real life), the latter ones being roughly 100 MB each.

// FIXME : instructions below should be replaced with usage of MUONTRKda
//

The raw.date.LDC* files are then processed using the DA online program (which is not built by default, but must be made
 explicitely using make daqDA-MCH from $ALICE_ROOT, and requires some DATE setup..., see \ref README_mchda )
 
<pre>
 MUONTRKda.exe -f raw.date.LCDi -a LDCi.ped (i=0,1,2,3)
 
 (repeat for each LDC)
</pre>

The LDCi.ped files are the input for the pedestal subprocessor,
which is tested using the TestMUONPreprocessor.C macro. 
The output of the pedestal subprocessor (upon success only) is written into the OCDB. 
Difference between the input and the output can be inferred using the diff() function
of MUONCDB.C macro.


\section shuttle_s5 HV

HV DCS values are created in CreateDCSAliasMap() of TestMUONPreprocessor.C
You might want to modify this function to create a given set of error conditions
 in order to test whether the HVSubprocessor is reacting properly to those errors.


\section shuttle_s6 GMS

The GMS alignment data for testing can be generated with
the macro MUONGenerateTestGMS.C:
The matrices of TGeoHMatrix type, with TObject::fUniqueID equal to the geometry
module Id, are put in a TClonesArray and saved in the Root file with a 
key "GMSarray".


\section shuttle_s7 Trigger DCS

HV and currents DCS values are created in CreateDCSAliasMap() of TestMUONPreprocessor.C
As done in Tracker HV, you might want to modify this function to create a given set of 
error conditions in order to test whether the TriggerDCSSubprocessor is reacting 
properly to those errors.

COMMENT: the trigger subprocessor requires trigger .dat files to be present.
In order to test only the DCS values when the remaining trigger files are not present, 
a workaround is put in the TestMUONPreprocessor.C which creates fake date files.
This results in some "ERROR" flags appearing in the LOG file (for masks and LUT), 
which are absolutely expected.
The fake date files are automatically removed at the end of the macro.


This chapter is defined in the READMEshuttle.txt file.

*/
 
