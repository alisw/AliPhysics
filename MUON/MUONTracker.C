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

// Macro MUONTracker.C (TO BE COMPILED)
// for testing the C++ reconstruction code
// Output is using aliroot standard output MUON.Tracks.root
// The output is a TClonesArray of AliMUONTracks.
// Gines MARTINEZ, Subatech, sep 2003

#include <TClonesArray.h>

#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONEventReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackHit.h"
#include "AliMUONTrackParam.h"

void MUONTracker (Text_t *FileName = "galice.root", Int_t FirstEvent = 0, Int_t LastEvent = 9999)
{
  //
  cout << "MUONTracker" << endl;
  cout << "FirstEvent " << FirstEvent << endl;
  cout << "LastEvent " << LastEvent << endl;
  cout << "FileName ``" << FileName << "''" << endl;
  
  // Creating Run Loader and openning file containing Hits, Digits and RecPoints
  AliRunLoader * RunLoader = AliRunLoader::Open(FileName,"Event","UPDATE");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",FileName);
    return;
  }
  // Loading AliRun master
  RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();
  RunLoader->LoadKinematics("READ");
  
  // Loading MUON subsystem
  AliMUON * MUON = (AliMUON *) gAlice->GetDetector("MUON");
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadHits("READ");
  MUONLoader->LoadRecPoints("READ");
  AliMUONData * muondata = MUON->GetMUONData();
  muondata->SetLoader(MUONLoader);
  
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  
  AliMUONEventReconstructor *Reco = new AliMUONEventReconstructor();

  // The right place for changing AliMUONEventReconstructor parameters
  // with respect to the default ones
  //   Reco->SetMaxSigma2Distance(100.0);
  //   Reco->SetPrintLevel(20);
  Reco->SetPrintLevel(1);
  //   Reco->SetBendingResolution(0.0);
  //   Reco->SetNonBendingResolution(0.0);
  cout << "AliMUONEventReconstructor: actual parameters" << endl;
  Reco->Dump();
  //   gObjectTable->Print();

  if  (LastEvent>nevents) LastEvent=nevents;
  // Loop over events
  for (Int_t event = FirstEvent; event < LastEvent; event++) {
    //MUONLoader->LoadHits("READ");
    MUONLoader->LoadRecPoints("READ");
    cout << "Event: " << event << endl;
    RunLoader->GetEvent(event);   
    muondata->SetTreeAddress("RC");
    if (MUONLoader->TreeT() == 0x0) MUONLoader->MakeTree("T");
    muondata->MakeBranch("RT");
    muondata->SetTreeAddress("RT");
    Reco->EventReconstruct();
    // Dump current event
    Reco->EventDump();

    // Duplicating rectrack data in muondata for output
    for(Int_t i=0; i<Reco->GetNRecTracks(); i++) {
      AliMUONTrack * track = (AliMUONTrack*) Reco->GetRecTracksPtr()->At(i);
      muondata->AddRecTrack(*track);
      //printf(">>> TEST TEST Number of hits in the track %d is %d \n",i,track->GetNTrackHits());
    }

    muondata->Fill("RT");
    MUONLoader->WriteTracks("OVERWRITE");
    muondata->ResetRecTracks();
    //MUONLoader->UnloadHits();
    MUONLoader->UnloadRecPoints();
  } // Event loop
}
