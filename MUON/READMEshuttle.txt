// $Id$

/*! 

\page README_shuttle Shuttle 

How to test the Shuttle preprocessor(s) for MUON.

We will get two "logical" MUON preprocessors : one for the tracker and one for the trigger.
Both will manage several subtasks (e.g. the tracker one will handle pedestals,
 gains and deadchannels, while the trigger one will handle masks and trigger lut)
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
MCH       PEDESTAL_RUN                read ASCII ped files      AliMUONPedestalSubprocessor
                                      and put them into OCDB
                        
MCH       GMS                         read GMS alignment files  AliMUONGMSSubprocessor
                                      and put them into OCDB

MCH       PHYSICS                     read DCS HV values and    AliMUONHVSubprocessor
                                      put them into OCDB
                                      
MCH       ELECTRONICS_CALIBRATION_RUN read ASCII gain files     prototype only = AliMUONGainSubprocessor
                                      and put them into OCDB
                                      
MTR       to be defined               to be defined             to be done
</pre>


\section shuttle_s3 Pedestals

Two options here. You can either use a pre-done set of ASCII pedestals files (generated as
 explained below for the 2nd option), located at /afs/cern.ch/user/l/laphecet/public/LDC*.ped, 
 or build you own set.
 
We've written an AliMUONPedestalEventGenerator which creates fake pedestal events. The pedestal values
are taken from the Offline Condition DataBase (OCDB) (which is itself fakely filled
using the WritePedestals() method of AliMUONCDB class

So first generate a valid pedestal CDB entry by using the AliMUONCDB class

<pre>
root[] const char* cdbpath="local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB"; // where to put the CDB
root[] AliMUONCDB cdb(cdbpath)
root[] Bool_t defaultValues = kFALSE; // to generate random values instead of plain zeros...
root[] Int_t startRun = 80;
root[] Int_t endRun = 80;
root[] cdb.WritePedestals(defaultValues, startRun, endRun);
</pre>

Expected output is (don't worry about the warnings, they are harmless) :

<pre>
I-AliMUONCDB::ManuList: Generating ManuList...
I-AliMUONCDB::ManuList: Manu List generated.
I-AliMUONCDB::MakePedestalStore: 16828 Manus and 1064008 channels.
I-AliMUONCDB::WritePedestals: Ngenerated = 1064008
I-AliCDBManager::Init: AliEn classes enabled in Root. AliCDBGrid factory registered.
I-AliCDBManager::SetDefaultStorage: Setting Default storage to: local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB
I-AliCDBLocal::PutEntry: CDB object stored into file ($ALICE_ROOT)/SHUTTLE/TestShuttle/TestCDB/MUON/Calib/Pedestals/Run80_80_v0_s0.root
</pre>

Then use the AliMUONPedestalEventGenerator to produce simulated pedestal events.

Usage (from the Root prompt) :
<pre>
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

The raw.date.LDC* files are then processed using the makeped online program 
(currently found, pending an agreement on where to put online programs under cvs,
 under /afs/cern.ch/user/a/abaldiss/public/v16; Please contact Alberto to check 
 it's the latest version) which outputs manus-*.ped ASCII files (one per LDC) :
 
<pre>
 makeped -f raw.date.LCDi -a LDCi.ped (i=0,1,2,3)
 
 (repeat for each LDC)
</pre>

The LDCi.ped files are the input for the pedestal subprocessor,
which is tested using the TestMUONPreprocessor.C macro. 
The output of the pedestal subprocessor (upon success only) is written into the OCDB. 
Difference between the input and the output can be inferred using the diff() function
of MUONCDB.C macro.


\section shuttle_s4 Gains

Like pedestals, you have two options here. You can either use a pre-done set of 
ASCII gain files (generated as explained below for the 2nd option), 
located at /afs/cern.ch/user/l/laphecet/public/LDC*.gains,  or build you own set.
 
We've written an AliMUONGainEventGenerator which creates fake gain events. 
The pedestal and gain values are taken from the Offline Condition DataBase (OCDB)
 (which is itself fakely filled using the WritePedestals() and WriteGains() 
 methods of AliMUONCDB class).

So first you need to generate a valid pedestal CDB entry and a valid gain CDB 
entry by using the AliMUONCDB class, from the Root prompt:

<pre>
const char* cdbpath="local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB"; // where to put the CDB
AliMUONCDB cdb(cdbpath)
Bool_t defaultValues = kFALSE; // to generate random values instead of plain zeros...
Int_t gainRun = 80;
Int_t pedRun = 81;
cdb.WritePedestals(defaultValues, pedRun, pedRun);
cdb.WriteGains(defaultValues, gainRun, gainRun);
</pre>

Expected output is (don't worry about the warnings, they are harmless) :

Then use the AliMUONGainEventGenerator to produce simulated gain events : the output 
will be n x 4 date files (n is the number of fake injections, currently 9, and 4
 is the number of LDCs)

Usage (again, from the Root prompt) :

<pre>
const char* cdbpath="local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB"; // where to get the CDB
AliCDBManager::Instance()->SetDefaultStorage(cdbpath); // so you will read 
// back pedestals and gain values generated in the previous step
const char* dateFileName = "raw.date"; // base filename for the output
Int_t gainRunNumber = 80; // run number used to fetch gains from OCDB
Int_t pedRunNumber = 81; // run number used to fetch the pedestals from the OCDB 
// generated ped files will be for r = 81, 83, etc...
Int_t nevents = 100; // # of events to generate. 100 should be enough for testing, but 1000 would be better for prod
gSystem->Load("libMUONshuttle"); // needed or not depending on whether it's already loaded or not
AliMUONGainEventGenerator g(gainRunNumber,pedRunNumber,nevents,dateFileName);
g.Exec("");
</pre>

It *will* take a lot of time (mainly due to the fact that we're writing a 
bunch of ASCII files = DDL files), so please be patient.

The output should be a sequence of directories, RUN81, RUN82, etc..., each 
containing the normal simulated sequence of MUON.Hits.root, MUON.SDigits.root,
 MUON.Digits.root, raw/*.ddl files and raw.date.LDCi where i=0-3 (i.e. one DATE file
per LDC, as will be used in real life), the latter ones being roughly 100 MB each.

<pre>
// FIXME
// Below should follow instructions on how to feed the MUONTRKda with the
// generated files.
</pre>


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

This chapter is defined in the READMEshuttle.txt file.

*/
 
