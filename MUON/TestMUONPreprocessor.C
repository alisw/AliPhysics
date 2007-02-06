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

/// This macro runs the test preprocessor for MUON.
/// It uses AliTestShuttle to simulate a full Shuttle process
///
/// The input data has to be created first by other processes (or is created
/// here by CreateDCSAliasMap() for tracker HV).
///
/// To play with it, you'll have to set/modify several lines, to
/// a) select input files, using shuttle->AddInputFile()
/// b) select run type, using shuttle->AddInputRunParameter() (the run type
///    dictates which task is really performed by the MUONPreprocessor
///
/// For more information on usage, please see READMEshuttle.
///
// By Laurent Aphecetche, SUBATECH Nantes

#include "TestMUONPreprocessor.h"
#include "AliCDBManager.h"
#include "AliShuttleInterface.h"
#include "AliCDBId.h"
#include "AliTestShuttle.h"
#include "TMap.h"
#include "Riostream.h"
#include "TSystem.h"
#include "AliMpExMap.h"
#include "TMap.h"
#include "TString.h"
#include "TObjArray.h"
#include "AliMpHelper.h"
#include "AliDCSValue.h"
#include "TObjString.h"
#include "TRandom.h"
#include "AliMUONPreprocessor.h"
#include "AliCDBEntry.h"
#include "AliMUONHVNamer.h"
#endif

void TestMUONPreprocessor(Int_t runNumber=1500)
{
  // load library
  gSystem->Load("../SHUTTLE/TestShuttle/libTestShuttle.so");
  gSystem->Load("libMUONshuttle.so");
  
  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/SHUTTLE/TestShuttle/TestCDB");

  // create AliTestShuttle instance
  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(runNumber, 0, 1);

  // Create DCS HV aliases
  TMap* dcsAliasMap = CreateDCSAliasMap();

  // now give the alias map to the shuttle
  shuttle->SetDCSInput(dcsAliasMap);

  // The shuttle can process files that originate from DCS, DAQ and HLT.
  // To test it, we provide some local files and locations where these would be found when
  // the online machinery would be there.
  // In real life this functions would be produces by the sub-detectors
  // calibration programs in DCS, DAQ or HLT. These files can then be retrieved using the Shuttle.
  //
  // Files are added with the function AliTestShuttle::AddInputFile. The syntax is:
  // AddInputFile(<system>, <detector>, <id>, <source>, <local-file>)
  // In this example we add 4 files originating from different LDCs but with the same id (PEDESTALS)

  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC0","$ALICE_ROOT/MUON/data/LDC0.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC1","$ALICE_ROOT/MUON/data/LDC1.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC2","$ALICE_ROOT/MUON/data/LDC2.ped");
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","PEDESTALS","LDC3","$ALICE_ROOT/MUON/data/LDC3.ped");
  
  // and GMS file
  shuttle->AddInputFile(AliTestShuttle::kDAQ,"MCH","GMS","GMS","$ALICE_ROOT/MUON/data/GMS.root");

  // The shuttle can read run parameters stored in the DAQ run logbook.
  // To test it, we must provide the run parameters manually. They will be retrieved in the preprocessor
  // using GetRunParameter function.
  // In real life the parameters will be retrieved automatically from the run logbook;
  shuttle->AddInputRunParameter("RunType", "PEDESTAL_RUN"); 
//  shuttle->AddInputRunParameter("RunType", "PHYSICS"); 
  // PEDESTAL_RUN -> pedestals
  // ELECTRONICS_CALIBRATION_RUN -> gains
  // PHYSICS ? -> HV
  
  // Create the preprocessor that should be tested, it registers itself automatically to the shuttle
  new AliMUONPreprocessor("MCH", shuttle);

  shuttle->Print();
  
  // Test the preprocessor
  shuttle->Process();
}

TMap* CreateDCSAliasMap()
{
  /// Creates a DCS structure for MUON Tracker HV
  ///
  /// The structure is the following:
  ///   TMap (key --> value)
  ///     <DCSAlias> --> <valueList>
  ///     <DCSAlias> is a string
  ///     <valueList> is a TObjArray of AliDCSValue
  ///     An AliDCSValue consists of timestamp and a value in form of a AliSimpleValue
  
  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(kTRUE);
  
  TRandom random(0);
  
  AliMUONHVNamer hvNamer;
  
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
        if ( aliasName == "MchHvLvLeft/Chamber01Left/Quad3Sect2.actual.vMon") value = 500;
//        if ( aliasName == "MchHvLvLeft/Chamber05Left/Slat07.actual.vMon") value = 1300;
//        if ( aliasName == "MchHvLvLeft/Chamber05Left/Slat03.actual.vMon") value = 2500;
        AliDCSValue* dcsValue = new AliDCSValue(value,timeStamp);
        valueSet->Add(dcsValue);
      }
      aliasMap->Add(new TObjString(*alias),valueSet);
    }
  }
  
  delete aliases;
    
  return aliasMap;
}

