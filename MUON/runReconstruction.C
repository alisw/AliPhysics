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

// Macro extracted from MUON test script
// By Laurent Aphecetche

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliMUONReconstructor.h"
#include "AliMUONRecoParam.h"
#include "AliCDBManager.h"
#include "AliMagFMaps.h"
#include "AliTracker.h"
#include "AliReconstruction.h"
#include <TRandom.h>
//#include <TObjectTable.h>
#endif

void runReconstruction(int run, int seed, const char* input, const char* recoptions)
{ 
  AliCDBManager::Instance()->SetRun(run);
  
  gRandom->SetSeed(seed);
  
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);
  
  AliReconstruction* MuonRec = new AliReconstruction("galice.root");
  MuonRec->SetInput(input);
  MuonRec->SetRunVertexFinder(kFALSE);
  MuonRec->SetRunLocalReconstruction("MUON");
  MuonRec->SetRunTracking("MUON");
  MuonRec->SetFillESD("");
  MuonRec->SetLoadAlignData("MUON");
  MuonRec->SetNumberOfEventsPerFile(1000);
  MuonRec->SetOption("MUON",recoptions);
  //  MuonRec->SetEventRange(319,319);
  MuonRec->SetWriteAOD();
  
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  AliMUONReconstructor::SetRecoParam(muonRecoParam);
  muonRecoParam->Print("FULL");
  
  MuonRec->Run();
  
  delete MuonRec;
  
  //gObjectTable->Print();
}

