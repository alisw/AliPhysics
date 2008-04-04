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
#include "AliMagFMaps.h"
#include "AliTracker.h"
#include "AliReconstruction.h"
#include <TRandom.h>
#include <TGrid.h>
//#include <TObjectTable.h>
#endif

// Data file, OCDB on Grid
TString input="alien:///alice/data/2008/LHC08a/000024841/raw/08000024841010.10.root";
//TString input="alien:///alice/data/2008/LHC08a/000021931/raw/08000021931001.50.root";
TString ocdbPath = "alien://folder=/alice/data/2008/LHC08a/OCDB";

// Data file, OCDB locally
//TString input="$ALICE_ROOT/MUON/test_out.100/raw.root";
//TString ocdbPath = "local://$ALICE_ROOT";

TString caliboption1 = "NOGAIN";
TString caliboption2 = "GAINCONSTANTCAPA";
TString recoptions = "SAVEDIGITS";
Int_t seed = 1234567;

void runDataReconstruction(Int_t calib = 1)
{ 
  TGrid::Connect("alien://");
  
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbPath.Data());
  //man->SetSpecificStorage("MUON/Calib/Mapping","local://$ALICE_ROOT");
  //man->SetSpecificStorage("MUON/Calib/DDLStore","local://$ALICE_ROOT");
  //man->SetSpecificStorage("MUON/Calib/Gains","local://$ALICE_ROOT");
  
  gRandom->SetSeed(seed);
  
  // no magnetic field --> factor (4th parameter) = 0
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 0., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);
  AliReconstruction *MuonRec = new AliReconstruction();
  
  
  MuonRec->SetInput(input.Data());
  MuonRec->SetRunVertexFinder(kFALSE);
  MuonRec->SetRunLocalReconstruction("MUON");
  MuonRec->SetRunTracking("MUON");
  MuonRec->SetFillESD(" ");
  MuonRec->SetLoadAlignData("MUON");
  MuonRec->SetNumberOfEventsPerFile(1000);
  MuonRec->SetOption("MUON",recoptions.Data());
  //  MuonRec->SetEventRange(319,319);
  MuonRec->SetWriteAOD();
  //MuonRec.SetEventRange(0,100); 
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetCosmicParam();
  muonRecoParam->BypassSt45(kTRUE);
  muonRecoParam->RequestStation(2,kFALSE);
  TString caliboption = caliboption1;
  if ( calib == 2 ) caliboption = caliboption2;
  muonRecoParam->SetCalibrationMode(caliboption.Data());
  //muonRecoParam->SetClusteringMode("PEAKFIT");
  //muonRecoParam->SetClusteringMode("PEAKCOG");
  muonRecoParam->Print("FULL");
  AliRecoParam::Instance()->RegisterRecoParam(muonRecoParam);
  
  MuonRec->Run();
  
  delete MuonRec;
  
  //gObjectTable->Print();
}

