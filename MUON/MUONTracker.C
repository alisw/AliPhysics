// Macro MUONTracker.C (TO BE COMPILED)
// for testing the C++ reconstruction code
// Output is using aliroot standard output MUON.Tracks.root
// The output is a TClonesArray of AliMUONTracks.
#include <TClonesArray.h>

#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONEventReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackHit.h"
#include "AliMUONTrackParam.h"

void MUONTracker (Int_t FirstEvent = 0, Int_t LastEvent = 0, Text_t *FileName = "galice.root")
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

  // Loop over events
  for (Int_t event = FirstEvent; event < LastEvent; event++) {
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
    }

    MUONLoader->TreeT()->Fill();
    MUONLoader->WriteTracks("OVERWRITE");
    muondata->ResetRecTracks();
  } // Event loop
  MUONLoader->UnloadRecPoints();
}
