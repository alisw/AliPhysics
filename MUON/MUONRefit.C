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
/// \file MUONRefit.C
/// \brief Macro for refitting ESD tracks from ESD pads
///
/// \author Philippe Pillot, SUBATECH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStopwatch.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Riostream.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TROOT.h>

// STEER includes
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliRecoParam.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliGeomManager.h"

// MUON includes
#include "AliMpCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONRefitter.h"
#include "AliMUONVDigit.h"
#include "AliMUONTrack.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#endif

const Int_t printLevel = 1;

void Prepare();
TTree* GetESDTree(TFile *esdFile);

//-----------------------------------------------------------------------
void MUONRefit(Int_t nevents = -1, const char* esdFileNameIn = "AliESDs.root", const char* esdFileNameOut = "AliESDs_New.root")
{
  /// refit ESD tracks from ESD pads (i.e. re-clusterized the attached ESD clusters);
  /// reset the charge of the digit using their raw charge before refitting;
  /// compare results with original ESD tracks; 
  /// write results in a new ESD file
  
  // prepare the refitting
  gRandom->SetSeed(1);
  Prepare();
  
  // reconstruction parameters for the refitting
  AliMUONRecoParam* recoParam = AliMUONRecoParam::GetLowFluxParam();
  Info("MUONRefit", "\n Reconstruction parameters for refitting:");
  recoParam->Print("FULL");
  
  AliMUONESDInterface esdInterface;
  AliMUONRefitter refitter(recoParam);
  refitter.Connect(&esdInterface);
  
  // open the ESD file and tree
  TFile* esdFile = TFile::Open(esdFileNameIn);
  TTree* esdTree = GetESDTree(esdFile);
  
  // create the ESD output file and tree
  TFile* newESDFile = TFile::Open(esdFileNameOut, "RECREATE");
  newESDFile->SetCompressionLevel(2);
  TTree* newESDTree = esdTree->CloneTree(0);
  
  // connect ESD event to the ESD tree
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);

  // timer start...
  TStopwatch timer;
  
  // Loop over ESD events
  if (nevents > 0) nevents = TMath::Min(nevents,(Int_t)esdTree->GetEntries());
  else nevents = (Int_t)esdTree->GetEntries();
  for (Int_t iEvent = 0; iEvent < nevents; iEvent++) {
    if (printLevel>0) cout<<endl<<"            ****************event #"<<iEvent+1<<"****************"<<endl;
    
    //----------------------------------------------//
    // -------------- process event --------------- //
    //----------------------------------------------//
    // get the ESD of current event
    esdTree->GetEvent(iEvent);
    if (!esd) {
      Error("MUONRefit", "no ESD object found for event %d", iEvent);
      return;
    }
    Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks();
    if (nTracks < 1) continue;
    
    // load the current event
    esdInterface.LoadEvent(*esd);
    
    // loop over digit to modify their charge
    AliMUONVDigit *digit;
    TIter next(esdInterface.CreateDigitIterator());
    while ((digit = static_cast<AliMUONVDigit*>(next()))) {
      digit->SetCharge(digit->ADC());
      digit->Calibrated(kFALSE);
    }
    
    // refit the tracks from digits
    AliMUONVTrackStore* newTrackStore = refitter.ReconstructFromDigits();
    
    //----------------------------------------------//
    // ------ fill new ESD and print results ------ //
    //----------------------------------------------//
    // loop over the list of ESD tracks
    TClonesArray *esdTracks = (TClonesArray*) esd->FindListObject("MuonTracks");
    for (Int_t iTrack = 0; iTrack <  nTracks; iTrack++) {
      
      // get the ESD track
      AliESDMuonTrack* esdTrack = (AliESDMuonTrack*) esdTracks->UncheckedAt(iTrack);
      
      // skip ghost tracks (leave them unchanged in the new ESD file)
      if (!esdTrack->ContainTrackerData()) continue;
      
      // get the corresponding MUON track
      AliMUONTrack* track = esdInterface.FindTrack(esdTrack->GetUniqueID());
      
      // Find the corresponding re-fitted MUON track
      AliMUONTrack* newTrack = (AliMUONTrack*) newTrackStore->FindObject(esdTrack->GetUniqueID());
      
      // replace the content of the current ESD track or remove it if the refitting has failed
      if (newTrack) {
	Double_t vertex[3] = {esdTrack->GetNonBendingCoor(), esdTrack->GetBendingCoor(), esdTrack->GetZ()};
	AliMUONESDInterface::MUONToESD(*newTrack, *esdTrack, vertex, esdInterface.GetDigits());
      } else {
	esdTracks->Remove(esdTrack);
      }
      
      // print initial and re-fitted track parameters at first cluster if any
      if (printLevel>0) {
	cout<<"            ----------------track #"<<iTrack+1<<"----------------"<<endl;
	cout<<"before refit:"<<endl;
	AliMUONTrackParam *param = (AliMUONTrackParam*) track->GetTrackParamAtCluster()->First();
	param->Print("FULL");
	if (printLevel>1) param->GetCovariances().Print();
	if (!newTrack) continue;
	cout<<"after refit:"<<endl;
	param = (AliMUONTrackParam*) newTrack->GetTrackParamAtCluster()->First();
	param->Print("FULL");
	if (printLevel>1) param->GetCovariances().Print();
	cout<<"            ----------------------------------------"<<endl;
      }
      
    }
    
    // free memory
    delete newTrackStore;
    
    // fill new ESD tree with new tracks
    esdTracks->Compress();
    newESDTree->Fill();
    
    if (printLevel>0) cout<<"            ****************************************"<<endl;
  }
  
  // ...timer stop
  timer.Stop();
  
  // write output ESD tree
  newESDFile->cd();
  newESDTree->Write();
  delete newESDTree;
  newESDFile->Close();
  
  // free memory
  esdFile->Close();
  delete esd;
  
  cout<<endl<<"time to refit: R:"<<timer.RealTime()<<" C:"<<timer.CpuTime()<<endl<<endl;
}

//-----------------------------------------------------------------------
void Prepare()
{
  /// Set the geometry, the magnetic field, the mapping and the reconstruction parameters
  
  // Import TGeo geometry (needed by AliMUONTrackExtrap::ExtrapToVertex)
  if (!gGeoManager) {
    AliGeomManager::LoadGeometry("geometry.root");
    if (!gGeoManager) {
      Error("MUONRefit", "getting geometry from file %s failed", "generated/galice.root");
      return;
    }
  }
  
  // set  mag field 
  // waiting for mag field in CDB 
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    printf("Loading field map...\n");
    AliMagF* field = new AliMagF("Maps","Maps",2,1.,1., 10.,AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }
  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();
  
  // Load mapping
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  if ( ! AliMpCDB::LoadDDLStore() ) {
    Error("MUONRefit","Could not access mapping from OCDB !");
    exit(-1);
  }
  
  // Load initial reconstruction parameters from OCDB
  AliMUONRecoParam* recoParam = 0x0;
  AliCDBPath path("MUON","Calib","RecoParam");
  AliCDBEntry *entry=man->Get(path.GetPath());
  if(entry) {
    recoParam = dynamic_cast<AliMUONRecoParam*>(entry->GetObject());
    entry->SetOwner(0);
    AliCDBManager::Instance()->UnloadFromCache(path.GetPath());
  }
  if (!recoParam) {
    printf("Couldn't find RecoParam object in OCDB: create default one");
    recoParam = AliMUONRecoParam::GetLowFluxParam();
  }
  Info("MUONRefit", "\n initial recontruction parameters:");
  recoParam->Print("FULL");
  AliMUONESDInterface::ResetTracker(recoParam);
  
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

