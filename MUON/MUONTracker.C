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

// Macro MUONTracker.C (TO BE COMPILED)
// for testing the C++ reconstruction code
// Output is using aliroot standard output MUON.Tracks.root
// The output is a TClonesArray of AliMUONTracks.
// Allow to make track reconstruction directly from AliTrackReference hits 
// recorded in TrackRefs.root file.

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONEventReconstructor.h"
#endif

void MUONTracker (Int_t FirstEvent = 0, Int_t LastEvent = 9999, Text_t *FileName = "galice.root")
{
  //
  cout << "MUONTracker" << endl;
  cout << "FirstEvent " << FirstEvent << endl;
  cout << "LastEvent " << LastEvent << endl;
  cout << "FileName ``" << FileName << "''" << endl;
  
  // Creating Run Loader and openning file containing Hits, Digits and RecPoints
  AliRunLoader * RunLoader = AliRunLoader::Open(FileName,"MUONLoader","UPDATE");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",FileName);
    return;
  }
  // Loading AliRun master
  if (RunLoader->GetAliRun() == 0x0) RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();
  RunLoader->LoadTrackRefs("READ");
  
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
 
  MUONLoader->LoadTracks("UPDATE");
  
  Int_t nevents;
  nevents = RunLoader->GetNumberOfEvents();
  
  AliMUONEventReconstructor* Reco = new AliMUONEventReconstructor(MUONLoader);
  AliMUONData* muondata = Reco->GetMUONData();

  // Testing if Tracker has already been done
  RunLoader->GetEvent(0);
  if (MUONLoader->TreeT()) {
    if (muondata->IsTrackBranchesInTree()) {
      MUONLoader->UnloadTracks();
      MUONLoader->LoadTracks("RECREATE");
      printf("Recreating tracks files\n");
    }
  }
  
  // The right place for changing AliMUONEventReconstructor parameters
  // with respect to the default ones
  //   Reco->SetMaxSigma2Distance(100.0);
  //   Reco->SetPrintLevel(20);
  Reco->SetPrintLevel(0);
  //   Reco->SetBendingResolution(0.0);
  //   Reco->SetNonBendingResolution(0.0);
  Reco->SetRecTrackRefHits(1); // 1: reconst. from track ref. hits 0: from clusters

  if (Reco->GetRecTrackRefHits() == 0)
    MUONLoader->LoadRecPoints("READ");

  if  (LastEvent > nevents-1) LastEvent=nevents-1;
  
  // Loop over events
  for (Int_t event = FirstEvent; event <= LastEvent; event++) {
    cout << "Event: " << event << endl;
    RunLoader->GetEvent(event);   
    // Test if trigger track has already been done before
    if (MUONLoader->TreeT() == 0x0) {	
      MUONLoader->MakeTracksContainer();
    }      else {
      if (muondata->IsTrackBranchesInTree()){ // Test if track has already been done before
	if (event==FirstEvent) MUONLoader->UnloadTracks();
	MUONLoader->MakeTracksContainer();  // Redoing Tracking
	Info("TrackContainer","Recreating TrackContainer and deleting previous ones");
      }
    }
    
    muondata->MakeBranch("RT");
    muondata->SetTreeAddress("RT");
    Reco->EventReconstruct();
    
    muondata->Fill("RT");
    MUONLoader->WriteTracks("OVERWRITE");  
    muondata->ResetRecTracks();
    if (Reco->GetRecTrackRefHits() == 0)
      muondata->ResetRawClusters();
  } // Event loop
  
  MUONLoader->UnloadRecPoints();
  MUONLoader->UnloadTracks();
  RunLoader->UnloadTrackRefs();
}
