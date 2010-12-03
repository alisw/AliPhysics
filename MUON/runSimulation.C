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

/// \ingroup macros
/// \file runSimulation.C
/// \brief Macro for running simulation
///
/// Macro extracted from the MUON test script
///
/// \author Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliSimulation.h"
#include <TRandom.h>
#endif

void runSimulation(int seed, 
                   int nevents, 
                   const char* config,
                   const char* embedwith)
{ 
// Uncoment following lines to run simulation with local residual mis-alignment
// (generated via MUONGenerateGeometryData.C macro)
// AliCDBManager* man = AliCDBManager::Instance();
// man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
// man->SetSpecificStorage("MUON/Align/Data","local://$ALICE_ROOT/OCDB/MUON/ResMisAlignCDB");

  AliSimulation MuonSim(config);
  
  if ( strlen(embedwith) > 0 )
  {
    // setups specific to embedding
    
    gAlice->SetConfigFunction("Config(\"\", \"param\", \"AliMUONDigitStoreV2S\",kTRUE);");
    
    // get the run number from real data
    
    AliRunLoader* runLoader = AliRunLoader::Open(embedwith,"titi");
    if (runLoader == 0x0) 
    {
      AliError(Form("Cannot open file %s",filename));    
      return;
    }
    
    runLoader->LoadHeader();
    
    if ( ! runLoader->GetHeader() ) {
      AliError("Cannot load header.");    
      return;
    }
    else {
      Int_t runNumber = runLoader->GetHeader()->GetRun();
      MuonSim.SetRunNumber(runNumber);
      cout << Form("***** RUN NUMBER SET TO %09d ACCORDING TO %s ",runNumber,embedwith) << endl;
    }  
    runLoader->UnloadHeader(); 
    delete runLoader;
    
    cout << "***** EMBEDDING MODE : USING RAW OCDB" << endl;
    AliCDBManager::Instance()->SetDefaultStorage("raw://");
    AliCDBManager::Instance()->SetSpecificStorage("local://$ALICE_ROOT/OCDB","MUON/Calib/Gains");
    AliCDBManager::Instance()->SetSpecificStorage("local://$ALICE_ROOT/OCDB","MUON/Align/Data");
    
  }
  else
  {
    gAlice->SetConfigFunction("Config(\"\", \"param\", \"AliMUONDigitStoreV2S\",kFALSE);");    
  }
  
  MuonSim.SetSeed(seed);
  MuonSim.SetTriggerConfig("MUON");
  MuonSim.SetWriteRawData("MUON HLT","raw.root",kTRUE);

  MuonSim.SetMakeDigits("MUON");
  MuonSim.SetMakeSDigits("MUON");
  MuonSim.SetMakeDigitsFromHits("");

  MuonSim.SetRunHLT("libAliHLTMUON.so chains=dHLT-sim");

  MuonSim.SetRunQA("MUON:ALL");
  
  if ( strlen(embedwith) > 0 ) 
  {
    MuonSim.MergeWith(embedwith);
  }
  
  MuonSim.Run(nevents);
  //gObjectTable->Print();

}
