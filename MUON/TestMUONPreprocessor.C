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
/// The input data has to be created first by other processes (or is created
/// here by CreateDCSAliasMap() for tracker HV).
///
/// To play with it, you'll have to set/modify several lines, to
/// - a) select input files, using shuttle->AddInputFile()
/// - b) select run type, using shuttle->AddInputRunParameter() (the run type
///      dictates which task is really performed by the MUONPreprocessor
///
/// You must load relevant libraries (besides normal MUON ones) before
/// compiling this macro :
/// <pre>
/// gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle");
/// gSystem->Load("libMUONshuttle.so");
/// </pre>
///
/// For more information on usage, please see the \ref README_shuttle page.
///
/// \author Laurent Aphecetche, SUBATECH Nantes

#include "TestMUONPreprocessor.h"

#include "AliMUONTrackerPreprocessor.h"
#include "AliMUONTriggerPreprocessor.h"

#include "AliLog.h"

#include "AliMpExMap.h"
#include "AliMpHelper.h"
#include "AliMpHVNamer.h"
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

void TestMUONPreprocessor(Int_t runNumber=80, const char* runType="CALIBRATION")
{
  // runType can be :
  //
  // PEDESTAL -> pedestals
  // CALIBRATION -> gains
  // PHYSICS -> HV
  // GMS
  
  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(runNumber, 0, 1);
  
  const char* inputCDB = "local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB";
  //const char* inputCDB = "alien://folder=/alice/testdata/2008/TS08a/OCDB";
  
  AliTestShuttle::SetMainCDB(inputCDB);
  AliTestShuttle::SetMainRefStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestReference");

  TString rt(runType);
  rt.ToUpper();
  
  if ( rt.Contains("PHYSICS") )
  {
    // Create DCS HV aliases
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

//  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC0","$ALICE_ROOT/MUON/data/LDC0.ped");
//  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC1","$ALICE_ROOT/MUON/data/LDC1.ped");
//  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC2","$ALICE_ROOT/MUON/data/LDC2.ped");
//  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC3","$ALICE_ROOT/MUON/data/LDC3.ped");

  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC0","/afs/cern.ch/user/l/laphecet/public/PEDESTALS/LDC0.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC1","/afs/cern.ch/user/l/laphecet/public/PEDESTALS/LDC1.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC2","/afs/cern.ch/user/l/laphecet/public/PEDESTALS/LDC2.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC3","/afs/cern.ch/user/l/laphecet/public/PEDESTALS/LDC3.ped");

  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","GAINS","LDC0","/afs/cern.ch/user/l/laphecet/public/GAINS/LDC0.gain");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","GAINS","LDC1","/afs/cern.ch/user/l/laphecet/public/GAINS/LDC1.gain");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","GAINS","LDC2","/afs/cern.ch/user/l/laphecet/public/GAINS/LDC2.gain");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","GAINS","LDC3","/afs/cern.ch/user/l/laphecet/public/GAINS/LDC3.gain");

  // and GMS file
  shuttle->AddInputFile(AliTestShuttle::kDCS,"MCH","GMS","GMS","/afs/cern.ch/user/l/laphecet/public/GMS/GMS.root");

  // and then the trigger stuff
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","LOCAL","LDC0","/afs/cern.ch/user/l/laphecet/public/TRIGGER/MtgLocalMask-1.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","REGIONAL","LDC0","/afs/cern.ch/user/l/laphecet/public/TRIGGER/MtgRegionalCrate-1.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","GLOBAL","LDC0","/afs/cern.ch/user/l/laphecet/public/TRIGGER/MtgGlobalCrate-1.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","LUT","LDC0","/afs/cern.ch/user/l/laphecet/public/TRIGGER/MtgLocalLut-1.dat");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MTR","EXPORTED","LDC0","/afs/cern.ch/user/l/laphecet/public/TRIGGER/ExportedFiles.dat");
  
  // The shuttle can read run parameters stored in the DAQ run logbook.
  // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
  // using GetRunParameter function.
  // In real life the parameters will be retrieved automatically from the run logbook;
  shuttle->SetInputRunType(runType);
  
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  new AliMUONTrackerPreprocessor(shuttle);
  new AliMUONTriggerPreprocessor(shuttle);
  
  shuttle->Print();
  
  // Test the preprocessor
  shuttle->Process();
}

TMap* CreateDCSAliasMap(const char* inputCDB, Int_t runNumber)
{
  /// Creates a DCS structure for MUON Tracker HV
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
  
  AliMpHVNamer hvNamer;
  
  TObjArray* aliases = hvNamer.GenerateAliases();
  
  for ( Int_t i = 0; i < aliases->GetEntries(); ++i ) 
  {
    TObjString* alias = static_cast<TObjString*>(aliases->At(i));
    TString& aliasName = alias->String();
    if ( aliasName.Contains("sw") ) 
    {
      // HV Switch (St345 only)
      TObjArray* valueSet = new TObjArray;
      valueSet->SetOwner(kTRUE);
      Bool_t value = kTRUE;
//      Float_t r = random.Uniform();
//      if ( r < 0.007 ) value = kFALSE;      
//      if ( aliasName.Contains("DE513sw2") ) value = kFALSE;
      
      for ( UInt_t timeStamp = 0; timeStamp < 60*3; timeStamp += 60 )
      {
        AliDCSValue* dcsValue = new AliDCSValue(value,timeStamp);
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
        Float_t value = random.Gaus(1750,62.5);
        if ( aliasName == "MchHvLvLeft/Chamber00Left/Quad2Sect1.actual.vMon") value = 500;
        AliDCSValue* dcsValue = new AliDCSValue(value,timeStamp);
        valueSet->Add(dcsValue);
      }
      if ( aliasName == "MchHvLvLeft/Chamber04Left/Slat06.actual.vMon" ) continue;
      aliasMap->Add(new TObjString(*alias),valueSet);
    }
  }
  
  delete aliases;
    
  return aliasMap;
}

