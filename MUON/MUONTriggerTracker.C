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

/* */

// Macro MUONTriggerTracker.C (TO BE COMPILED) 
// for testing the C++ trigger reconstruction code
// Output is using aliroot standard output MUON.Tracks.root
// The output is a TClonesArray of AliMUONTriggerTracks.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONEventReconstructor.h"
#endif

void MUONTriggerTracker (Text_t *FileName = "galice.root", Int_t FirstEvent = 0, Int_t LastEvent = 9999)
{
  //
  cout << "MUONTriggerTracker" << endl;
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
  // RunLoader->LoadKinematics("READ");
  
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
 
  Int_t nevents;
  nevents = RunLoader->GetNumberOfEvents();

  MUONLoader->LoadRecPoints("READ");
  MUONLoader->LoadTracks("UPDATE"); 


  AliMUONEventReconstructor *Reco = new AliMUONEventReconstructor(MUONLoader);
  AliMUONData* muondata = Reco->GetMUONData();

  // Testing if Trigger Tracking has already been done
  RunLoader->GetEvent(0);
  if (MUONLoader->TreeT()) {
      if (muondata->IsTriggerTrackBranchesInTree()) {
      MUONLoader->UnloadTracks();
      MUONLoader->LoadTracks("RECREATE");
      printf("Recreating Tracks files\n");
    }
  }

  Reco->SetPrintLevel(0);    
  if  (LastEvent>nevents) LastEvent=nevents;

  // Loop over events
  for (Int_t event = FirstEvent; event < LastEvent; event++) {
      cout << "Event: " << event << endl;
      RunLoader->GetEvent(event);
      if (MUONLoader->TreeT() == 0x0) MUONLoader->MakeTracksContainer();      
      
      muondata->MakeBranch("RL");
      muondata->SetTreeAddress("RL");
      Reco->EventReconstructTrigger();
 
      muondata->Fill("RL");
      MUONLoader->WriteTracks("OVERWRITE");  
      muondata->ResetRecTriggerTracks();
      muondata->ResetTrigger();
  } // Event loop

  MUONLoader->UnloadRecPoints();
  MUONLoader->UnloadTracks();
}
