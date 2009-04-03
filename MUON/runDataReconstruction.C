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

/* $Id: runReconstruction.C 23207 2007-12-20 09:59:20Z ivana $ */

/// \ingroup macros
/// \file runDataReconstruction.C
/// \brief Macro for running reconstruction
///
/// Macro for running reconstruction on the cosmics run data.
///
/// \author Laurent Aphecetche, Nicole Bastid, Bogdan Vulpescu, ...

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliMUONReconstructor.h"
#include "AliMUONRecoParam.h"
#include "AliRecoParam.h"
#include "AliCDBManager.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliReconstruction.h"
#include <TRandom.h>
#include <TGrid.h>
#include <TSystem.h>
//#include <TObjectTable.h>
#endif

TString caliboption1 = "NOGAIN";
TString caliboption2 = "GAINCONSTANTCAPA";
TString recoptions = "SAVEDIGITS";
Int_t seed = 1234567;

void runDataReconstruction(const char* input = "/Users/laurent/Alice/Data/Raw/09000067495031.10.root",
                           const char* ocdbPath = "alien://folder=/alice/data/2009/OCDB",
                           Int_t calib = 1)
{ 
  TGrid::Connect("alien://");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbPath);

  gRandom->SetSeed(seed);

  TString socdb(ocdbPath);
  if ( socdb.Contains("local://") )
  {
    // no magnetic field
    AliMagF* field = new AliMagF("Maps","Maps",2,0.,0., 10.,AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
    TGeoGlobalMagField::Instance()->Lock();
  }
  
  AliReconstruction *MuonRec = new AliReconstruction();
  
  MuonRec->SetInput(gSystem->ExpandPathName(input));
  MuonRec->SetRunVertexFinder(kFALSE);
  MuonRec->SetRunLocalReconstruction("MUON");
  MuonRec->SetRunTracking("MUON");
  MuonRec->SetFillESD(" ");
  MuonRec->SetLoadAlignData("MUON");
  MuonRec->SetNumberOfEventsPerFile(0);
  MuonRec->SetOption("MUON",recoptions.Data());
  
  // reconstruction parameters
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetCosmicParam();
  
  // digit selection
  muonRecoParam->SetPadGoodnessMask(0x400BE80);
  TString caliboption = caliboption1;
  if ( calib == 2 ) caliboption = caliboption2;
  muonRecoParam->SetCalibrationMode(caliboption.Data());
  
  // chamber resolution (incuding misalignment)
  for (Int_t iCh=0; iCh<10; iCh++) {
    muonRecoParam->SetDefaultNonBendingReso(iCh,0.4);
    muonRecoParam->SetDefaultBendingReso(iCh,0.4);
  }
  muonRecoParam->SetMaxNonBendingDistanceToTrack(10.);
  muonRecoParam->SetMaxBendingDistanceToTrack(10.);
  
  // cut on (non)bending slopes
  //muonRecoParam->SetMaxNonBendingSlope(0.6);
  //muonRecoParam->SetMaxBendingSlope(0.6);
  
  // tracking algorithm
  muonRecoParam->MakeMoreTrackCandidates(kTRUE);
  muonRecoParam->RequestStation(0, kFALSE);
  muonRecoParam->RequestStation(2, kFALSE);
  muonRecoParam->RequestStation(3, kFALSE);
  muonRecoParam->RequestStation(4, kFALSE);
  muonRecoParam->SetSigmaCutForTracking(7.);
  muonRecoParam->ImproveTracks(kTRUE, 7.);
  
  muonRecoParam->Print("FULL");
  
  MuonRec->SetRecoParam("MUON",muonRecoParam);
  
  MuonRec->SetRunQA("MUON:ALL");

  MuonRec->Run();
  
  delete MuonRec;
  
  //gObjectTable->Print();
}

