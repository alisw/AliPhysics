#include "AliAnalysisTaskNanoPt.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliNanoAODTrack.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskNanoPt)

//--------------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskNanoPt::AliAnalysisTaskNanoPt()
  : AliAnalysisTaskSE(),
    fInputEvent(nullptr),
    fEvent(nullptr),
    fEvtCuts(nullptr),
    fTrack(nullptr),
    fProtonTrack(nullptr),
    fAntiProtonTrack(nullptr),
    fDeuteronTrack(nullptr),
    fAntiDeuteronTrack(nullptr),
    fConfig(nullptr),
    fIsMC(false),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    fGTI(nullptr),
    fProtonRestMass(nullptr),
    fAntiProtonRestMass(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fTrackBufferSize(2500) {}
//-----------------------------------------------------------------------------------------------------------------------

AliAnalysisTaskNanoPt::AliAnalysisTaskNanoPt(
  const char *name, const bool isMC)
  : AliAnalysisTaskSE(name),
    fInputEvent(nullptr),
    fEvent(nullptr),
    fEvtCuts(nullptr),
    fTrack(nullptr),
    fProtonTrack(nullptr),
    fAntiProtonTrack(nullptr),
    fDeuteronTrack(nullptr),
    fAntiDeuteronTrack(nullptr),
    fConfig(nullptr),
    fIsMC(isMC),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    fGTI(nullptr),
    fProtonRestMass(nullptr),
    fAntiProtonRestMass(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fTrackBufferSize(2500) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Dueteron Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiDeuteron Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA
  if (isMC) {
    DefineOutput(8, TList::Class());  //Output for the Proton MC
    DefineOutput(9, TList::Class());  //Output for the AntiProton MC
    DefineOutput(10, TList::Class());  //Output for the Deuteron MC
    DefineOutput(11, TList::Class());  //Output for the AntiDeuteron MC
  }
}



//---------------------------------------------------------------------------------------------------------------------------

AliAnalysisTaskNanoPt::~AliAnalysisTaskNanoPt() {
  delete fEvent;
  delete fTrack;
  delete fProtonTrack;
  delete fAntiProtonTrack;
  delete fDeuteronTrack;
  delete fAntiDeuteronTrack;
  delete fPairCleaner;
  delete fPartColl;
}
//-----------------------------------------------------------------------------------------------------------------

Float_t AliAnalysisTaskNanoPt::GetMass2sq(AliFemtoDreamTrack *track) const {
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (!(beta > 0)) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

//-------------------------------------------------------UserCreateOutPut----------------------------------------------------------------------------------

void AliAnalysisTaskNanoPt::UserCreateOutputObjects() {

  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEvtCuts) {
    AliError("No Event cuts \n");
  } else {
    fEvtCuts->InitQA();
  }

  if (!fProtonTrack) {
    AliError("No Proton cuts \n");
  } else {
    fProtonTrack->Init();
    fProtonRestMass = new TH2F("fProtonRestMass", "Proton", 36, 0.5, 4.05, 180, 0.2, 2.);
    fProtonRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fProtonRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fProtonList = fProtonTrack->GetQAHists();
    fProtonList->Add(fProtonRestMass);
  }

  if (!fAntiProtonTrack) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProtonTrack->Init();
    fAntiProtonRestMass = new TH2F("fAntiProtonRestMass", "AntiProton", 36, 0.5, 4.05, 180, 0.2, 2.);
    fAntiProtonRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiProtonRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiProtonList = fAntiProtonTrack->GetQAHists();
    fAntiProtonList->Add(fAntiProtonRestMass);
  }

  if (!fDeuteronTrack) {
    AliError("No Proton cuts \n");
  } else {
    fDeuteronTrack->Init();
    fDeuteronRestMass = new TH2F("fDeuteronRestMass", "Deuteron", 36, 0.5, 4.05, 180, 0.2, 2.);
    fDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronList = fDeuteronTrack->GetQAHists();
    fDeuteronList->Add(fDeuteronRestMass);
  }

  if (!fAntiDeuteronTrack) {
    AliError("No Proton cuts \n");
  } else {
    fAntiDeuteronTrack->Init();
    fAntiDeuteronRestMass = new TH2F("fAntiDeuteronRestMass", "AntiDeuteron", 36, 0.5, 4.05, 180, 0.2, 2.);
    fAntiDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronList = fAntiDeuteronTrack->GetQAHists();
    fAntiDeuteronList->Add(fAntiDeuteronRestMass);

  }

  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
        fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 0,
        fConfig->GetMinimalBookingME());
  }

  fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates(), false);
  //fMCEvent = new AliMCEvent();
  fTrack = new AliFemtoDreamTrack();
  //fProtonTrack->SetUseMCInfo();
  //fAntiProtonTrack->SetUseMCInfo();
// fDeuteronTrack->SetUseMCInfo();
// fAntiDeuteronTrack->SetUseMCInfo();
   fTrack->SetUseMCInfo(fIsMC);
  if (!fEvtCuts->GetMinimalBooking()) {
    fEvtList = fEvtCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }


  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");

  if (fConfig->GetUseEventMixing()) {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME()) {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fDeuteronList);
  PostData(5, fAntiDeuteronList);
  PostData(6, fResults);
  PostData(7, fResultsQA);

//-----------------------------------------MCTracks------------------------------------------------------------------------------
  if (fProtonTrack->GetIsMonteCarlo()) {
    if (!fProtonTrack->GetMinimalBooking()) {
      fProtonMCList = fProtonTrack->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCProtOnTrackCuts");
      fProtonMCList->SetOwner();
    }
    PostData(8, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    if (!fAntiProtonTrack->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProtonTrack->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiProtonTrackCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(9, fAntiProtonMCList);
  }

  if (fDeuteronTrack->GetIsMonteCarlo()) {
    if (!fDeuteronTrack->GetMinimalBooking()) {
      fDeuteronMCList = fDeuteronTrack->GetMCQAHists();
    } else {
      fDeuteronMCList = new TList();
      fDeuteronMCList->SetName("MCDeuteronTrackCuts");
      fDeuteronMCList->SetOwner();
    }
    PostData(10, fDeuteronMCList);
  }
  if (fAntiDeuteronTrack->GetIsMonteCarlo()) {
    if (!fAntiDeuteronTrack->GetMinimalBooking()) {
      fAntiDeuteronMCList = fAntiDeuteronTrack->GetMCQAHists();
    } else {
      fAntiDeuteronMCList = new TList();
      fAntiDeuteronMCList->SetName("MCAntiDeuteronTrackCuts");
      fAntiDeuteronMCList->SetOwner();
    }
    PostData(11, fAntiDeuteronMCList);
  }

}

//------------------------------------------UserExec()----------------------------------------------------------------------------


void AliAnalysisTaskNanoPt::UserExec(Option_t  *option ) {
  AliVEvent *fInputEvent = InputEvent();


  // PREAMBLE - CHECK EVERYTHING IS THERE
  if (!fInputEvent) {
    AliError("No Input event");
    return;
  }

  if (fIsMC && !fMCEvent) {
    AliError("No MC event");
    return;
  }


  if (!fEvtCuts) {
    AliError("Event Cuts missing");
    return;
  }

  if (!fProtonTrack || !fAntiProtonTrack) {
    AliError("Proton Cuts missing");
    return;
  }

  if (!fDeuteronTrack || !fAntiDeuteronTrack) {
    AliError("Deuteron Cuts missing");
    return;
  }


  // EVENT SELECTION
  fEvent->SetEvent(fInputEvent);
  if (!fEvtCuts->isSelected(fEvent))
    return;

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }

    StoreGlobalTrackReference(track);
  }

  std::vector<AliFemtoDreamBasePart> Proton;
  std::vector<AliFemtoDreamBasePart> AntiProton;
  std::vector<AliFemtoDreamBasePart> Deuteron;
  std::vector<AliFemtoDreamBasePart> AntiDeuteron;

  const int multiplicity = fEvent->GetMultiplicity();
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  //fTrack->SetUseMCInfo(true);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent, multiplicity);

    if (fProtonTrack->isSelected(fTrack)) {
      fProtonRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      Proton.push_back(*fTrack);
    }
    if (fAntiProtonTrack->isSelected(fTrack)) {
      fAntiProtonRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      AntiProton.push_back(*fTrack);
    }
    if (fDeuteronTrack->isSelected(fTrack)) {
      fDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      Deuteron.push_back(*fTrack);
    }
    if (fAntiDeuteronTrack->isSelected(fTrack)) {
      fAntiDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      AntiDeuteron.push_back(*fTrack);
    }
  }
  if (fIsMC) {
    AliAODInputHandler *eventHandler =
      dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()
                                        ->GetInputEventHandler());
    AliMCEvent* fMC = eventHandler->MCEvent();
    for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = static_cast<AliAODMCParticle *>(fMC->GetTrack(iPart));
      if (TMath::Abs(mcPart->Eta()) < 0.8 && mcPart->IsPhysicalPrimary()) {
        if (mcPart->GetPdgCode() == fProtonTrack->GetPDGCode()) {
          fProtonTrack->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiProtonTrack->GetPDGCode()) {
          fAntiProtonTrack->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fDeuteronTrack->GetPDGCode()) {
          fDeuteronTrack->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiDeuteronTrack->GetPDGCode()) {
          fAntiDeuteronTrack->FillGenerated(mcPart->Pt());
        }
      }
    }
  }


  fPairCleaner->CleanTrackAndDecay(&Proton, &Deuteron, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProton, &AntiDeuteron, 1);


  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Proton);
  fPairCleaner->StoreParticle(AntiProton);
  fPairCleaner->StoreParticle(Deuteron);
  fPairCleaner->StoreParticle(AntiDeuteron);


  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());


  // flush the data
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fDeuteronList);
  PostData(5, fAntiDeuteronList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  //-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fProtonTrack->GetIsMonteCarlo()) {
    PostData(8, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    PostData(9, fAntiProtonMCList);
  }

  if (fDeuteronTrack->GetIsMonteCarlo()) {
    PostData(10, fDeuteronMCList);
  }
  if (fAntiDeuteronTrack->GetIsMonteCarlo()) {
    PostData(11, fAntiDeuteronMCList);
  }
}


///------------------------------------------------------------------------
void AliAnalysisTaskNanoPt::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//-------------------------------------------------------------------------

void AliAnalysisTaskNanoPt::StoreGlobalTrackReference(AliVTrack * track) {
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }
  if (fGTI[trackID]) {

    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}


