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
#include <Riostream.h>
//#include <TObjectTable.h>
#endif

void runReconstruction(int seed, const char* input, const char* recoptions, bool rawocdb)
{ 
  AliCDBManager* man = AliCDBManager::Instance();
  
  if ( rawocdb ) 
  {
    cout << "**** WILL USE RAW OCDB" << endl;
    man->SetDefaultStorage("raw://"); //alien://folder=/alice/data/2011/OCDB?cacheFold=/Users/laurent/OCDBcache");
    man->SetSpecificStorage("ITS/Calib/RecoParam","alien://folder=/alice/cern.ch/user/p/ppillot/OCDB_PbPbSim");
  } 
  else
  {
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

    man->SetSpecificStorage("GRP/GRP/Data",
                            Form("local://%s",gSystem->pwd()));

  }
  
  gRandom->SetSeed(seed);
  
  AliReconstruction* MuonRec = new AliReconstruction("galice.root");
  MuonRec->SetInput(gSystem->ExpandPathName(input));
  MuonRec->SetRunReconstruction("MUON ITS");
  MuonRec->SetFillESD("HLT");
  MuonRec->SetOption("HLT", "libAliHLTMUON.so");
  MuonRec->SetNumberOfEventsPerFile(10000);
  MuonRec->SetOption("MUON",recoptions);
  MuonRec->SetRunQA("MUON:ALL");
  MuonRec->SetQAWriteExpert(AliQAv1::kMUON);
  MuonRec->SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  MuonRec->SetWriteESDfriend(kFALSE);
  MuonRec->SetCleanESD(kFALSE);  
  MuonRec->SetStopOnError(kFALSE);
  
  // uncomment the following lines if you want to set custom RecoParam
  // instead of getting them from the OCDB
  //  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  //  muonRecoParam->SaveFullClusterInESD(kTRUE,100.);
  //  MuonRec->SetRecoParam("MUON",muonRecoParam);
  
  MuonRec->Run();
  
  delete MuonRec;
  
  //gObjectTable->Print();
}
