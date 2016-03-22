/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#if !defined(__CINT__) || defined(__MAKECINT__)

/// \ingroup macros
/// \file TestMUONPreprocessor.C
/// \brief The macro for testing the shuttle preprocessors
///
/// This macro runs the test preprocessor for MUON.
/// It uses AliTestShuttle to simulate a full Shuttle process
///
/// You must make a link of some OCDB entries to have the mapping loaded
/// correctly :
///
/// <pre>
/// cd $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB
/// mkdir -p MUON/Calib/MappingData
/// cd MUON/Calib/MappingData/
/// ln -si $ALICE_ROOT/OCDB/MUON/Calib/MappingData/* .
/// cd $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB
/// mkdir -p MUON/Calib/Config
/// cd MUON/Calib/Config
/// ln -si $ALICE_ROOT/OCDB/MUON/Calib/Config/* .
/// </pre>
///
/// and Align/Baseline if you'd like to test GMS subprocessor :
///
/// <pre>
/// cd $ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB
/// mkdir -p MUON/Align
/// cd MUON/Align
/// ln -si $ALICE_ROOT/OCDB/MUON/Align/Baseline .
/// </pre>
///
/// The input data has to be created first by other processes (or is created
/// here by CreateDCSAliasMap() for tracker and trigger HV).
///
/// To play with it, you'll have to set/modify several lines, to
/// - a) select input files, using shuttle->AddInputFile()
/// - b) select run type, using shuttle->AddInputRunParameter() (the run type
///      dictates which task is really performed by the MUONPreprocessor
///
/// The sourceDirectory is there to "emulate" what the real preprocessor will
/// find on the FXS, and is assumed to have the following structure :
/// <pre>
/// CONFIG/
///    LDC0.config
///    LDC1.config
///    LDC2.config
///    LDC3.config
/// CONFIGPAR/
///    LDC0.config
///    LDC1.config
///    LDC2.config
///    LDC3.config
/// GMS/
///    GMS.root
/// OCCUPANCY/
///    mch.occupancy
/// BPEVO/
///    mchbpevo.root
/// PEDESTALS/
///    LDC0.ped
///    LDC1.ped
///    LDC2.ped
///    LDC3.ped
/// TRIGGER/
///    ExportedFiles.dat (mandatory)
///    MtgGlobalCrate-1.dat
///    MtgLocalLut-1.dat
///    MtgLocalMask-1.dat
///    MtgRegionalCrate-1.dat
/// </pre>
///
///
/// An example set of input files can be found at https://cernbox.cern.ch/public.php?service=files&t=c9363b0a1f92daf4963dd4f8b2a8f72a
///
/// (just get the zip file and unpack it into a directory that will be sourceDirectory)
///
/// IMPORTANT:
/// The trigger files have to be present in order for the algorithm to work correctly.
/// If you want to test the Trigger DCS maps only, but you don't have the .dat trigger files,
/// you have to create dummy files through :
/// <pre>
/// cd sourceDirectory/TRIGGER
/// echo -e "MtgLocalMask-1.dat\nMtgRegionalCrate-1.dat\nMtgGlobalCrate-1.dat\nMtgLocalLut-1.dat" > ExportedFiles.dat
/// touch MtgLocalMask-1.dat MtgRegionalCrate-1.dat MtgGlobalCrate-1.dat MtgLocalLut-1.dat
/// </pre>
///
/// For more information on usage, please see the \ref README_shuttle page.
///
/// \author Laurent Aphecetche, SUBATECH Nantes; \n
///         Diego Stocco, SUBATECH Nantes

#include "TestMUONPreprocessor.h"

#include "AliMUONTrackerPreprocessor.h"
#include "AliMUONTriggerPreprocessor.h"

#include "AliLog.h"

#include "AliMpBusPatch.h"
#include "AliMpExMap.h"
#include "AliMpHelper.h"
#include "AliMpDDLStore.h"
#include "AliMpDCSNamer.h"
#include "AliMpCDB.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "AliShuttleInterface.h"
#include "AliTestShuttle.h"
#include "AliDCSValue.h"

#include "Riostream.h"
#include "TSystem.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TRandom.h"

#endif

//______________________________________________________________________________
void TestMUONPreprocessor(Int_t runNumber=80,
                          const char* runType="PHYSICS",
                          const char* sourceDirectory="$HOME/Downloads/muontestshuttle")
{
  // runType can be :
  //
  // PEDESTAL -> pedestals
  // PHYSICS -> HV
  // GMS

  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime

  AliTestShuttle* shuttle = new AliTestShuttle(runNumber, 0, 1);

  const char* inputCDB = "local://$ALICE_ROOT/../src/SHUTTLE/TestShuttle/TestCDB";

  AliTestShuttle::SetMainCDB(inputCDB);
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/OCDB/SHUTTLE/TestShuttle/TestReference");

  TString rt(runType);
  rt.ToUpper();

  if ( rt.Contains("PHYSICS") || rt.Contains("CALIBRATION") )
  {
    // Create DCS aliases
    TMap* dcsAliasMap = CreateDCSAliasMap(inputCDB, runNumber);

    if ( dcsAliasMap )
    {
      // now give the alias map to the shuttle
      shuttle->SetDCSInput(dcsAliasMap);
    }
  }

  printf("Test Shuttle temp dir: %s\n", AliShuttleInterface::GetShuttleTempDir());
  printf("Test Shuttle log dir: %s\n", AliShuttleInterface::GetShuttleLogDir());
  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

  // The shuttle can process files that originate from DCS, DAQ and HLT.
  // To test it, we provide some local files and locations where these would be found when
  // the online machinery would be there.
  // In real life this functions would be produces by the sub-detectors
  // calibration programs in DCS, DAQ or HLT. These files can then be retrieved using the Shuttle.
  //
  // Files are added with the function AliTestShuttle::AddInputFile. The syntax is:
  // AddInputFile(<system>, <detector>, <id>, <source>, <local-file>)
  // In this example we add 4 files originating from different LDCs but with the same id (PEDESTALS)

  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC0",Form("%s/PEDESTALS/LDC0.ped",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC1",Form("%s/PEDESTALS/LDC1.ped",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC2",Form("%s/PEDESTALS/LDC2.ped",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC3",Form("%s/PEDESTALS/LDC3.ped",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC4",Form("%s/PEDESTALS/LDC4.ped",sourceDirectory));

  if ( rt.Contains("PHYSICS") )
  {
    // simulate a change of configuration during the run
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC0",Form("%s/CONFIGPAR/LDC0.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC1",Form("%s/CONFIGPAR/LDC1.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC2",Form("%s/CONFIGPAR/LDC2.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC3",Form("%s/CONFIGPAR/LDC3.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC4",Form("%s/CONFIGPAR/LDC4.conf",sourceDirectory));
  }
  else
  {
    // configuration as done by the pedestal run
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC0",Form("%s/CONFIG/LDC0.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC1",Form("%s/CONFIG/LDC1.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC2",Form("%s/CONFIG/LDC2.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC3",Form("%s/CONFIG/LDC3.conf",sourceDirectory));
    shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","CONFIG","LDC4",Form("%s/CONFIG/LDC4.conf",sourceDirectory));
  }

  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","OCCUPANCY","MON",Form("%s/OCCUPANCY/mch.occupancy",sourceDirectory));

  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","BPEVO","MON",Form("%s/BPEVO/mchbpevo.root",sourceDirectory));

  // and GMS file
  shuttle->AddInputFile(AliTestShuttle::kDCS,"MCH","GMS","GMS",Form("%s/GMS/GMS.root",sourceDirectory));

  // and then the trigger stuff
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","LOCAL","LDC0",Form("%s/TRIGGER/MtgLocalMask-1.dat",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","REGIONAL","LDC0",Form("%s/TRIGGER/MtgRegionalCrate-1.dat",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","GLOBAL","LDC0",Form("%s/TRIGGER/MtgGlobalCrate-1.dat",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","LUT","LDC0",Form("%s/TRIGGER/MtgLocalLut-1.dat",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","EXPORTED","LDC0",Form("%s/TRIGGER/ExportedFiles.dat",sourceDirectory));
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","TRIGSCAL","LDC0",Form("%s/TRIGGER/MtgTrigScalers.dat",sourceDirectory));

  // The shuttle can read run parameters stored in the DAQ run logbook.
  // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
  // using GetRunParameter function.
  // In real life the parameters will be retrieved automatically from the run logbook;
  shuttle->SetInputRunType(runType);

  shuttle->AddInputRunParameter("totalEvents","20");

  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  new AliMUONTrackerPreprocessor(shuttle);
  new AliMUONTriggerPreprocessor(shuttle);

  shuttle->Print();

  // Test the preprocessor
  shuttle->Process();
}

//______________________________________________________________________________
void GenerateConfig()
{
  /// Generate "fake" configuration files for the tracker. One per LDC.

  Bool_t undefStorage(kFALSE);

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet())
  {
    undefStorage = kTRUE;
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man->SetRun(0);
  }

  // Load mapping
  Bool_t ok = AliMpCDB::LoadDDLStore();

  if (undefStorage)
  {
    man->UnsetDefaultStorage();
  }

  if (!ok)
  {
    AliErrorGeneral("GenerateConfig","Could not load DDLStore from OCDB");
    return;
  }

  ofstream* files[5];
  for ( Int_t i = 0; i < 5; ++i )
  {
    files[i]=0;
  }

  TIter next(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp;

  while ( ( bp = static_cast<AliMpBusPatch*>(next()) ) )
  {
    Int_t ddl = bp->GetDdlId();

    Int_t ldc = ddl/4;

    if (!files[ldc])
    {
      files[ldc] = new ofstream(Form("LDC%d.conf",ldc));
      *(files[ldc]) << "# changed" << endl;
    }

    for ( Int_t imanu = 0; imanu < bp->GetNofManus(); ++imanu )
    {
      *(files[ldc]) << bp->GetId() << " " << bp->GetManuId(imanu) << endl;
    }
  }

  for ( Int_t i = 0; i < 5; ++i )
  {
    if ( files[i] ) files[i]->close();
    delete files[i];
  }
}

//______________________________________________________________________________
TMap* CreateDCSAliasMap(const char* inputCDB, Int_t runNumber)
{
  /// Creates a DCS structure for MUON Tracker HV + LV and Trigger DCS and Currents
  ///
  /// The structure is the following:
  ///   TMap (key --> value)
  ///     <DCSAlias> --> <valueList>
  ///     <DCSAlias> is a string
  ///     <valueList> is a TObjArray of AliDCSValue
  ///     An AliDCSValue consists of timestamp and a value in form of a AliSimpleValue

  Bool_t undefStorage(kFALSE);

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet())
  {
    undefStorage = kTRUE;
    man->SetDefaultStorage(inputCDB);
    man->SetRun(runNumber);
  }

  // Load mapping
  Bool_t ok = AliMpCDB::LoadDDLStore();

  if (undefStorage)
  {
    man->UnsetDefaultStorage();
  }

  if (!ok)
  {
    AliErrorGeneral("CreateDCSAliasMap","Could not load DDLStore from OCDB");
    return 0x0;
  }

  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(kTRUE);

  TRandom random(0);

  const Char_t* detName[2] = { "TRACKER", "TRIGGER" };

  for(Int_t idet=0; idet<2; idet++){

    TString sDetName(detName[idet]);
    sDetName.ToUpper();

    AliMpDCSNamer dcsNamer(detName[idet]);

    TObjArray* aliases = dcsNamer.GenerateAliases();

    for ( Int_t i = 0; i < aliases->GetEntries(); ++i )
    {
      TObjString* alias = static_cast<TObjString*>(aliases->At(i));
      TString& aliasName = alias->String();
      if ( aliasName.Contains("sw") && sDetName.Contains("TRACKER"))
      {
        // HV Switch (St345 only)
        TObjArray* valueSet = new TObjArray;
        valueSet->SetOwner(kTRUE);
        Bool_t bvalue = kTRUE;
        //      Float_t r = random.Uniform();
        //      if ( r < 0.007 ) value = kFALSE;
        //      if ( aliasName.Contains("DE513sw2") ) value = kFALSE;

        for ( UInt_t timeStamp = 0; timeStamp < 60*3; timeStamp += 60 )
        {
          AliDCSValue* dcsValue = new AliDCSValue(bvalue,timeStamp);
          valueSet->Add(dcsValue);
        }
        aliasMap->Add(new TObjString(*alias),valueSet);
      }
      else
      {
        TObjArray* valueSet = new TObjArray;
        valueSet->SetOwner(kTRUE);
        for ( UInt_t timeStamp = 0; timeStamp < 60*15; timeStamp += 120 )
        {
          Float_t value = 0;
          if(sDetName.Contains("TRACKER")){
            if ( aliasName.Contains("vMon"))
            {
              value = random.Gaus(1750,62.5);
              if ( aliasName == "MchHvLvLeft/Chamber00Left/Quad2Sect1.actual.vMon") value = 500;
            }
            else if ( aliasName.Contains("iMon") )
            {
              value = random.Gaus(10,2);
              if ( aliasName == "MchHvLvLeft/Chamber00Left/Quad2Sect1.actual.iMon") value = 50;
            }
            else if ( aliasName.Contains("ann"))
            {
              value = random.Gaus(2.62,0.05);
            }
            else if ( aliasName.Contains("anp"))
            {
              value = random.Gaus(2.72,0.05);
            }
            else if ( aliasName.Contains("dig"))
            {
              value = random.Gaus(3.4,0.05);
            }
            else if ( aliasName.Contains("Crocu"))
            {
              value = random.Gaus(3.35,0.05);
            }
          }
          else if(aliasName.Contains("iMon")){
            value = random.Gaus(2.,0.4);
          }
          else {
            value = random.Gaus(8000.,16.);
          }

          AliDCSValue* dcsValue = new AliDCSValue(value,timeStamp);
          valueSet->Add(dcsValue);
        }
        if ( aliasName == "MchHvLvLeft/Chamber04Left/Slat06.actual.vMon" ) continue;
        if ( aliasName == "MTR_INSIDE_MT22_RPC3_HV.vEff" ) continue;
        if ( aliasName == "MTR_OUTSIDE_MT21_RPC4_HV.actual.iMon" ) continue;
        aliasMap->Add(new TObjString(*alias),valueSet);
      }
    } // loop on aliases

    delete aliases;
  } // loop on detectors (tracker and trigger)

  AliMpCDB::UnloadAll();

  return aliasMap;
}
