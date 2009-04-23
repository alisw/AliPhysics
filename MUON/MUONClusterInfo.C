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
#include <TClonesArray.h>
#include <TTree.h>
#include <TString.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TMath.h>

// STEER includes
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliRecoParam.h"
#include "AliCDBManager.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

// MUON includes
#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"
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

void Prepare();
TTree* GetESDTree(TFile *esdFile);
UInt_t buildClusterMap(AliMUONTrack &track);

//-----------------------------------------------------------------------
void MUONClusterInfo(Int_t nevents = -1, const char* esdFileName = "AliESDs.root",
		     const char* inFileName = "galice.root", const char* outFileName = "clusterInfo.root")
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
  
  // prepare the refitting during ESD->MUON conversion
  Prepare();
  
  // open the ESD file and tree and connect the ESD event
  TFile* esdFile = 0x0;
  TTree* esdTree = 0x0;
  AliESDEvent* esd = 0x0;
  if (useESD) {
    esdFile = TFile::Open(esdFileName);
    esdTree = GetESDTree(esdFile);
    esd = new AliESDEvent();
    esd->ReadFromTree(esdTree);
  }
  
  // get the cluster from RecPoints
  AliRunLoader * rl = 0x0;
  AliLoader* MUONLoader = 0x0;
  if (useRecPoints) {
    rl = AliRunLoader::Open(inFileName,"MUONLoader");
    MUONLoader = rl->GetDetectorLoader("MUON");
    MUONLoader->LoadRecPoints("READ");   
    MUONLoader->LoadDigits("READ");
  }
  
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
      digitStore = AliMUONVDigitStore::Create(*treeD);
      if ( digitStore != 0x0 ) {
	digitStore->Clear();
	digitStore->Connect(*treeD);
	treeD->GetEvent(0);
      }
    }
    
    // prepare access to calibration data
    if (useESD && !calibData) calibData = new AliMUONCalibrationData(esd->GetESDRun()->GetRunNumber());
    else if (!calibData) calibData = new AliMUONCalibrationData(rl->GetRunNumber());
    
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
	  for (Int_t i=0; i<cluster->GetNDigits(); i++) {
	    AliMUONVDigit* digit = digitStore->FindObject(cluster->GetDigitId(i));
	    if (!digit) continue;
	    
	    // pad location
	    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->
	    GetMpSegmentation(digit->DetElemId(),AliMp::GetCathodType(digit->Cathode()));
	    AliMpPad pad = seg->PadByIndices(digit->PadX(), digit->PadY());
	    
	    // calibration parameters
	    AliMUONVCalibParam* ped = calibData->Pedestals(digit->DetElemId(), digit->ManuId());
	    AliMUONVCalibParam* gain = calibData->Gains(digit->DetElemId(), digit->ManuId());
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
	    padInfo.SetGain(gain->ValueAsFloatFast(manuChannel,0), gain->ValueAsFloatFast(manuChannel,1),
			    gain->ValueAsFloatFast(manuChannel,2), gain->ValueAsFloatFast(manuChannel,3));
	    
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
	for (Int_t i=0; i<cluster->GetNDigits(); i++) {
	  AliMUONVDigit* digit = digitStore->FindObject(cluster->GetDigitId(i));
	  if (!digit) continue;
	  
	  // pad location
	  const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->
	  GetMpSegmentation(digit->DetElemId(),AliMp::GetCathodType(digit->Cathode()));
	  AliMpPad pad = seg->PadByIndices(digit->PadX(), digit->PadY());
	  
	  // calibration parameters
	  AliMUONVCalibParam* ped = calibData->Pedestals(digit->DetElemId(), digit->ManuId());
	  AliMUONVCalibParam* gain = calibData->Gains(digit->DetElemId(), digit->ManuId());
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
	  padInfo.SetGain(gain->ValueAsFloatFast(manuChannel,0), gain->ValueAsFloatFast(manuChannel,1),
			  gain->ValueAsFloatFast(manuChannel,2), gain->ValueAsFloatFast(manuChannel,3));
	  
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
    delete rl;
  }
  if (useESD) {
    esdFile->Close();
    delete esd;
  }
  cout<<endl<<"time to fill cluster/track info: R:"<<timer.RealTime()<<" C:"<<timer.CpuTime()<<endl<<endl;
}

//-----------------------------------------------------------------------
void Prepare()
{
  /// Set the magnetic field, the mapping and the reconstruction parameters
  
  gRandom->SetSeed(0);
  
  // set mag field
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    printf("Loading field map...\n");
    AliMagF* field = new AliMagF("Maps","Maps",2,1.,1., 10.,AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }
  
  // Load mapping
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  if ( ! AliMpCDB::LoadDDLStore() ) {
    Error("Prepare","Could not access mapping from OCDB !");
    exit(-1);
  }
  
  // Reset the reconstruction parameters for track refitting if needed
  // (by default will use Kalman filter + Smoother)
  //  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  //  AliMUONESDInterface::ResetTracker(muonRecoParam);
  
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

