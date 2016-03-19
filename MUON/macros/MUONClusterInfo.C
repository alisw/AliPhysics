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
/// \file MUONClusterInfo.C
/// \brief Macro to fill AliMUONClusterInfo objects
///
/// \author Philippe Pillot, SUBATECH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStopwatch.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TString.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TMath.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliRecoParam.h"
#include "AliCDBManager.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliHeader.h"

// MUON includes
#include "AliMpConstants.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"
#include "AliMUONCDB.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONPadInfo.h"
#include "AliMUONClusterInfo.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#endif

const Int_t printLevel = 1;

TTree* GetESDTree(TFile *esdFile);
UInt_t buildClusterMap(AliMUONTrack &track);

//-----------------------------------------------------------------------
void MUONClusterInfo(Int_t nevents = -1, const char* esdFileName = "AliESDs.root", const char* inFileName = "galice.root",
		     const TString ocdbPath = "local://$ALICE_ROOT/OCDB", const char* outFileName = "clusterInfo.root")
{
  /// 1) if (esdFileName != "")
  /// loop over ESD event and fill AliMUONClusterInfo object with cluster + corresponding track parameters;
  /// 2) if (inFileName != "")
  /// loop over RecPoints and fill AliMUONClusterInfo object with cluster not attached to a track;
  /// 3) write results in a new root file.
  ///
  /// ******************************************* WARNING ******************************************* ///
  /// track parameters at each cluster are recomputed by the interface using Kalman filter + Smoother ///
  /// (It can be changed by resetting the tracker in the interface with a new recoParam object)       ///
  /// and the magnetic field set in the function prepare()                                            ///
  /// ******************************************* WARNING ******************************************* ///
  
  Bool_t useESD = (strcmp(esdFileName,""));
  Bool_t useRecPoints = (strcmp(inFileName,""));
  if (!useESD && !useRecPoints) {
    Error("MUONClusterInfo","you must provide ESD and/or galice root file(s)");
    return;
  }
  
  AliMUONClusterInfo* clusterInfo = new AliMUONClusterInfo();
  AliMUONPadInfo padInfo;
  AliMUONCalibrationData* calibData = 0x0;
  AliMUONESDInterface esdInterface;
  AliMUONVClusterStore* clusterStore = 0x0;
  AliMUONVDigitStore* digitStore = 0x0;
  
  // open the ESD file and tree and connect the ESD event
  TFile* esdFile = 0x0;
  TTree* esdTree = 0x0;
  AliESDEvent* esd = 0x0;
  Int_t runNumber = -1;
  if (useESD) {
    esdFile = TFile::Open(esdFileName);
    esdTree = GetESDTree(esdFile);
    esd = new AliESDEvent();
    esd->ReadFromTree(esdTree);
    if (esdTree->GetEvent(0) <= 0) {
      Error("MUONClusterInfo", "no ESD object found for event 0");
      return;
    }
    runNumber = esd->GetRunNumber();
  }
  
  // get the cluster from RecPoints
  AliRunLoader * rl = 0x0;
  AliLoader* MUONLoader = 0x0;
  if (useRecPoints) {
    rl = AliRunLoader::Open(inFileName,"MUONLoader");
    MUONLoader = rl->GetDetectorLoader("MUON");
    MUONLoader->LoadRecPoints("READ");   
    MUONLoader->LoadDigits("READ");
    rl->LoadHeader();
    if (runNumber < 0) runNumber = rl->GetHeader()->GetRun();
  }
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  AliCDBManager::Instance()->SetRun(runNumber);
  if (!AliMUONCDB::LoadField()) return;
  if (!AliMUONCDB::LoadMapping(kTRUE)) return;
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  // reset tracker for track restoring initial track parameters at cluster
  AliMUONESDInterface::ResetTracker(recoParam);
  
  // prepare access to calibration data
  calibData = new AliMUONCalibrationData(runNumber);
  
  // prepare the output tree
  gROOT->cd();
  TFile* clusterInfoFile = TFile::Open(outFileName, "RECREATE");
  clusterInfoFile->SetCompressionLevel(1);
  
  TTree* clusterInfoTree = new TTree("clusterInfoTree","clusterInfoTree");
  clusterInfoTree->Branch("clusterInfo", &clusterInfo, 32000, 99);
  
  // timer start...
  TStopwatch timer;
  
  // Loop over events
  if (useESD) {
    if (nevents > 0) nevents = TMath::Min(nevents,(Int_t)esdTree->GetEntries());
    else nevents = (Int_t)esdTree->GetEntries();
  } else {
    if (nevents > 0) nevents = TMath::Min(nevents,(Int_t)rl->GetNumberOfEvents());
    else nevents = (Int_t)rl->GetNumberOfEvents();
  }
  for (Int_t iEvent = 0; iEvent < nevents; iEvent++) {
    
    //----------------------------------------------//
    // -------------- process event --------------- //
    //----------------------------------------------//
    // get the ESD of current event
    if (useESD) {
      esdTree->GetEvent(iEvent);
      if (!esd) {
        Error("MUONClusterInfo", "no ESD object found for event %d", iEvent);
        return;
      }
      // load the current esd event
      esdInterface.LoadEvent(*esd);
      // get digit store
      if (!useRecPoints) digitStore = esdInterface.GetDigits();
    }
    
    if (useRecPoints) {
      if (!(rl->GetEvent(iEvent) == 0)) {
        Error("MUONClusterInfo", "unable to load event %d", iEvent);
	return;
      }
      // get the clusters of current event
      TTree* treeR = MUONLoader->TreeR();
      clusterStore = AliMUONVClusterStore::Create(*treeR);
      if ( clusterStore != 0x0 ) {
	clusterStore->Clear();
	clusterStore->Connect(*treeR);
	treeR->GetEvent(0);
      }
      // get the digits of current event
      TTree* treeD = MUONLoader->TreeD();
      if (treeD) {
	digitStore = AliMUONVDigitStore::Create(*treeD);
	if ( digitStore != 0x0 ) {
	  digitStore->Clear();
	  digitStore->Connect(*treeD);
	  treeD->GetEvent(0);
	}
      }
    }
    
    //----------------------------------------------//
    // ------------- fill cluster info ------------ //
    //----------------------------------------------//
    // --- loop over the refitted tracks ---
    if (useESD) {
      
      TIter nextTrack(esdInterface.CreateTrackIterator());
      AliMUONTrack* track;
      while ((track = static_cast<AliMUONTrack*>(nextTrack()))) {
	
	UInt_t muonClusterMap = buildClusterMap(*track);
	
	// loop over clusters
	AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->First());
	while (trackParam) {
	  clusterInfo->Clear("C");
	  
	  // fill cluster info
	  AliMUONVCluster* cluster = trackParam->GetClusterPtr();
	  clusterInfo->SetRunId(esd->GetRunNumber());
	  clusterInfo->SetEventId(iEvent);
	  clusterInfo->SetZ(cluster->GetZ());
	  clusterInfo->SetClusterId(cluster->GetUniqueID());
	  clusterInfo->SetClusterXY(cluster->GetX(), cluster->GetY());
	  clusterInfo->SetClusterXYErr(cluster->GetErrX(), cluster->GetErrY());
	  clusterInfo->SetClusterChi2(cluster->GetChi2());
	  clusterInfo->SetClusterCharge(cluster->GetCharge());
	  
	  // fill track info
	  clusterInfo->SetTrackId(track->GetUniqueID());
	  clusterInfo->SetTrackXY(trackParam->GetNonBendingCoor(), trackParam->GetBendingCoor());
	  clusterInfo->SetTrackThetaXY(TMath::ATan(trackParam->GetNonBendingSlope()), TMath::ATan(trackParam->GetBendingSlope()));
	  clusterInfo->SetTrackP(trackParam->P());
	  const TMatrixD paramCov = trackParam->GetCovariances();
	  clusterInfo->SetTrackXYErr(TMath::Sqrt(paramCov(0,0)), TMath::Sqrt(paramCov(2,2)));
	  clusterInfo->SetTrackChi2(track->GetNormalizedChi2());
	  clusterInfo->SetTrackCharge((Short_t)trackParam->GetCharge());
	  clusterInfo->SetTrackNHits(track->GetNClusters());
	  clusterInfo->SetTrackChamberHitMap(muonClusterMap);
	  
	  // fill pad info if available	  
	  if (digitStore) for (Int_t i=0; i<cluster->GetNDigits(); i++) {
	    AliMUONVDigit* digit = digitStore->FindObject(cluster->GetDigitId(i));
	    if (!digit) continue;
	    
	    // pad location
	    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->
	    GetMpSegmentation(digit->DetElemId(),AliMp::GetCathodType(digit->Cathode()));
	    AliMpPad pad = seg->PadByIndices(digit->PadX(), digit->PadY());
	    
	    // calibration parameters
	    AliMUONVCalibParam* ped = calibData->Pedestals(digit->DetElemId(), digit->ManuId());
	    Int_t manuChannel = digit->ManuChannel();
	    Int_t planeType = 0;
	    if ( digit->ManuId() & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) {
	      planeType = 1;
	    }
	    
	    // fill pad info
	    padInfo.SetPadId(digit->GetUniqueID());
	    padInfo.SetPadPlaneType(planeType);
	    padInfo.SetPadXY(pad.GetPositionX(), pad.GetPositionY());
	    padInfo.SetPadDimXY(pad.GetDimensionX(), pad.GetDimensionY());
	    padInfo.SetPadCharge((Double_t)digit->Charge());
	    padInfo.SetPadADC(digit->ADC());
	    padInfo.SetSaturated(digit->IsSaturated());
	    padInfo.SetCalibrated(digit->IsCalibrated());
	    padInfo.SetPedestal(ped->ValueAsFloatFast(manuChannel,0), ped->ValueAsFloatFast(manuChannel,1));
	    
	    clusterInfo->AddPad(padInfo);
	  }
	  
	  // remove clusters attached to a track
	  if (useRecPoints) {
	    AliMUONVCluster* cl = clusterStore->FindObject(cluster->GetUniqueID());
	    if (cl) clusterStore->Remove(*cl);
	  }
	  
	  // fill cluster info tree
	  clusterInfoTree->Fill();
	  
	  trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->After(trackParam));
	}
	
      }
      
    }
    
    // --- loop over clusters not attached to a track ---
    if (useRecPoints) {
      
      TIter nextCluster(clusterStore->CreateIterator());
      AliMUONVCluster *cluster;
      while ( ( cluster = static_cast<AliMUONVCluster*>(nextCluster()) ) ) {
	
	clusterInfo->Clear("C");
	
	// fill cluster info
	clusterInfo->SetRunId(rl->GetRunNumber());
	clusterInfo->SetEventId(iEvent);
	clusterInfo->SetZ(cluster->GetZ());
	clusterInfo->SetClusterId(cluster->GetUniqueID());
	clusterInfo->SetClusterXY(cluster->GetX(), cluster->GetY());
	clusterInfo->SetClusterXYErr(cluster->GetErrX(), cluster->GetErrY());
	clusterInfo->SetClusterChi2(cluster->GetChi2());
	clusterInfo->SetClusterCharge(cluster->GetCharge());
	
	// fill dummy track info
	clusterInfo->SetTrackId(0);
	clusterInfo->SetTrackXY(0.,0.);
	clusterInfo->SetTrackThetaXY(0.,0.);
	clusterInfo->SetTrackP(0.);
	clusterInfo->SetTrackXYErr(0.,0.);
	clusterInfo->SetTrackChi2(0.);
	clusterInfo->SetTrackCharge(0);
	clusterInfo->SetTrackNHits(0);
	clusterInfo->SetTrackChamberHitMap(0);
	
	// fill pad info if available	  
	if (digitStore) for (Int_t i=0; i<cluster->GetNDigits(); i++) {
	  AliMUONVDigit* digit = digitStore->FindObject(cluster->GetDigitId(i));
	  if (!digit) continue;
	  
	  // pad location
	  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->
	  GetMpSegmentation(digit->DetElemId(),AliMp::GetCathodType(digit->Cathode()));
	  AliMpPad pad = seg->PadByIndices(digit->PadX(), digit->PadY());
	  
	  // calibration parameters
	  AliMUONVCalibParam* ped = calibData->Pedestals(digit->DetElemId(), digit->ManuId());
	  Int_t manuChannel = digit->ManuChannel();
	  Int_t planeType = 0;
	  if ( digit->ManuId() & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) {
	    planeType = 1;
	  }
	  
	  // fill pad info
	  padInfo.SetPadId(digit->GetUniqueID());
	  padInfo.SetPadPlaneType(planeType);
	  padInfo.SetPadXY(pad.GetPositionX(), pad.GetPositionY());
	  padInfo.SetPadDimXY(pad.GetDimensionX(), pad.GetDimensionY());
	  padInfo.SetPadCharge((Double_t)digit->Charge());
	  padInfo.SetPadADC(digit->ADC());
	  padInfo.SetSaturated(digit->IsSaturated());
	  padInfo.SetCalibrated(digit->IsCalibrated());
	  padInfo.SetPedestal(ped->ValueAsFloatFast(manuChannel,0), ped->ValueAsFloatFast(manuChannel,1));
	  
	  clusterInfo->AddPad(padInfo);
	}
	
	// fill cluster info tree
	clusterInfoTree->Fill();
	
      }
      
      delete digitStore;
      delete clusterStore;
      
    }
    
  }
  
  // ...timer stop
  timer.Stop();
  printf("Writing Tree\n");
  // write output tree
  clusterInfoFile->cd();
  clusterInfoTree->Write();
  printf("Deleting Tree\n");
  delete clusterInfoTree;
  printf("Closing File\n");
  clusterInfoFile->Close();
  
  // free memory
  printf("Deleting calibData\n");
  delete calibData;
  printf("Deleting clusterInfo\n");
  delete clusterInfo;
  if (useRecPoints) {
    MUONLoader->UnloadDigits();
    MUONLoader->UnloadRecPoints();
    rl->UnloadHeader();
    delete rl;
  }
  if (useESD) {
    esdFile->Close();
    delete esd;
  }
  cout<<endl<<"time to fill cluster/track info: R:"<<timer.RealTime()<<" C:"<<timer.CpuTime()<<endl<<endl;
}

//-----------------------------------------------------------------------
TTree* GetESDTree(TFile *esdFile)
{
  /// Check that the file is properly open
  /// Return pointer to the ESD Tree
  
  if (!esdFile || !esdFile->IsOpen()) {
    Error("GetESDTree", "opening ESD file failed");
    exit(-1);
  }
  
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("GetESDTree", "no ESD tree found");
    exit(-1);
  }
  
  return tree;
  
}

      
//-----------------------------------------------------------------------
UInt_t buildClusterMap(AliMUONTrack &track)
{
  /// Build the map of clusters in tracking chambers
  
  UInt_t muonClusterMap = 0;
  
  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
  while (trackParam) {
    
    muonClusterMap |= BIT(trackParam->GetClusterPtr()->GetChamberId());
    
    trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
  }
  
  return muonClusterMap;
  
}

