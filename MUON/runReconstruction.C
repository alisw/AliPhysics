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
/// \file runReconstruction.C
/// \brief Macro for running reconstruction
///
/// Macro extracted from the MUON test script
///
/// \author Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliMUONReconstructor.h"
#include "AliMUONRecoParam.h"
#include "AliRecoParam.h"
#include "AliCDBManager.h"
#include "AliTracker.h"
#include "AliReconstruction.h"
#include <TRandom.h>
//#include <TObjectTable.h>
#endif

void runReconstruction(int seed, const char* input, const char* recoptions)
{ 
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  gRandom->SetSeed(seed);
  
  AliReconstruction* MuonRec = new AliReconstruction("galice.root");
  MuonRec->SetInput(input);
  MuonRec->SetRunVertexFinder(kFALSE);
  MuonRec->SetRunLocalReconstruction("MUON");
  MuonRec->SetRunTracking("MUON");
  MuonRec->SetFillESD("");
  MuonRec->SetLoadAlignData("MUON");
  MuonRec->SetNumberOfEventsPerFile(1000);
  MuonRec->SetOption("MUON",recoptions);
  MuonRec->SetRunQA("MUON:ALL");
  MuonRec->SetQAWriteExpert(AliQA::kMUON);
  // uncomment the following lines if you want to set custom RecoParam
  // instead of getting them from the OCDB
  //  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  //  muonRecoParam->SaveFullClusterInESD(kTRUE,100.);
  //  MuonRec->SetRecoParam("MUON",muonRecoParam);
  
  MuonRec->Run();
  
  delete MuonRec;
  
  //gObjectTable->Print();
}
