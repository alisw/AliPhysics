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
    fDeuteronTrackNoTOF(nullptr),
    fAntiDeuteronTrackNoTOF(nullptr),
    fConfig(nullptr),
    fIsMC(false),
    fIsMCTruth(false),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    fDeuteronNoTOFList(nullptr),
    fDeuteronMCNoTOFList(nullptr),
    fAntiDeuteronNoTOFList(nullptr),
    fAntiDeuteronMCNoTOFList(nullptr),
    fGTI(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fProtonDeuteronDump(nullptr),
    fAntiProtonAntiDeuteronDump(nullptr),
    fDumpster(nullptr),
    fUseDumpster(nullptr),
    fTrackBufferSize(2500) {
}
//-----------------------------------------------------------------------------------------------------------------------

AliAnalysisTaskNanoPt::AliAnalysisTaskNanoPt(const char *name, const bool isMC)
  : AliAnalysisTaskSE(name),
    fInputEvent(nullptr),
    fEvent(nullptr),
    fEvtCuts(nullptr),
    fTrack(nullptr),
    fProtonTrack(nullptr),
    fAntiProtonTrack(nullptr),
    fDeuteronTrack(nullptr),
    fAntiDeuteronTrack(nullptr),
    fDeuteronTrackNoTOF(nullptr),
    fAntiDeuteronTrackNoTOF(nullptr),
    fConfig(nullptr),
    fIsMC(isMC),
    fIsMCTruth(false),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fDeuteronList(nullptr),
    fDeuteronMCList(nullptr),
    fAntiDeuteronList(nullptr),
    fAntiDeuteronMCList(nullptr),
    fDeuteronNoTOFList(nullptr),
    fDeuteronMCNoTOFList(nullptr),
    fAntiDeuteronNoTOFList(nullptr),
    fAntiDeuteronMCNoTOFList(nullptr),
    fGTI(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fProtonDeuteronDump(nullptr),
    fAntiProtonAntiDeuteronDump(nullptr),
    fDumpster(nullptr),
    fUseDumpster(nullptr),
    fTrackBufferSize(2500) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Dueteron Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiDeuteron Cuts
  DefineOutput(6, TList::Class());  //Output for the DueteronNoTOF Cuts
  DefineOutput(7, TList::Class());  //Output for the AntiDeuteronNoTOF Cuts
  DefineOutput(8, TList::Class());  //Output for the Results
  DefineOutput(9, TList::Class());  //Output for the Results QA
  DefineOutput(10, TList::Class());  //Output for the Dumpster
  if (fIsMC) {
    DefineOutput(11, TList::Class());  //Output for the Proton MC
    DefineOutput(12, TList::Class());  //Output for the AntiProton MC
    DefineOutput(13, TList::Class());  //Output for the Deuteron MC
    DefineOutput(14, TList::Class());  //Output for the AntiDeuteron MC
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
  delete fDeuteronTrackNoTOF;
  delete fAntiDeuteronTrackNoTOF;
  delete fPairCleaner;
  delete fPartColl;

  if (fProtonDeuteronDump) {
    delete fProtonDeuteronDump;
  }
  if (fAntiProtonAntiDeuteronDump) {
    delete fAntiProtonAntiDeuteronDump;
  }
  if (fDumpster) {
    delete fDumpster;
  }
}
//-----------------------------------------------------------------------------------------------------------------
Float_t AliAnalysisTaskNanoPt::GetMass2sq(AliFemtoDreamTrack *track) const {
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (beta > 0) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}
//-------------------------------------------------------UserCreateOutPut----------------------------------------------------------------------------------
void AliAnalysisTaskNanoPt::UserCreateOutputObjects() {

  fGTI = new AliVTrack*[fTrackBufferSize];

  if (!fEvtCuts) {
    AliError("No Event cuts \n");
  } else {
    fEvtCuts->InitQA();
  }

  if (!fProtonTrack) {
    AliError("No Proton cuts \n");
  } else {
    fProtonTrack->Init();
    fProtonList = fProtonTrack->GetQAHists();
    if (fIsMC) {
      fProtonMCList = fProtonTrack->GetMCQAHists();
    }
  }

  if (!fAntiProtonTrack) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProtonTrack->Init();
    fAntiProtonList = fAntiProtonTrack->GetQAHists();
    if (fIsMC) {
      fAntiProtonMCList = fAntiProtonTrack->GetMCQAHists();
    }
  }

  if (!fDeuteronTrack) {
    AliError("No deuteron cuts \n");
  } else {
    fDeuteronTrack->Init();
    fDeuteronRestMass = new TH2F("fDeuteronRestMass", "Deuteron", 72, 0.5, 8.05,
                                 800, 0.00, 10.0);
    fDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronList = fDeuteronTrack->GetQAHists();
    if (fIsMC) {
      fDeuteronMCList = fDeuteronTrack->GetMCQAHists();
    }
    fDeuteronList->Add(fDeuteronRestMass);
  }

  if (!fAntiDeuteronTrack) {
    AliError("No Proton cuts \n");
  } else {
    fAntiDeuteronTrack->Init();
    fAntiDeuteronRestMass = new TH2F("fAntiDeuteronRestMass", "AntiDeuteron", 72,
                                     0.5, 8.05, 800, 0.00, 10.0);
    fAntiDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronList = fAntiDeuteronTrack->GetQAHists();
    if (fIsMC) {
      fAntiDeuteronMCList = fAntiDeuteronTrack->GetMCQAHists();
    }
    fAntiDeuteronList->Add(fAntiDeuteronRestMass);
  }
//--------------------------------------------------------------------------------------------------------------------
  if (!fDeuteronTrackNoTOF) {
    AliError("No DeuteronNoTOF cuts \n");
  } else {
    fDeuteronTrackNoTOF->Init();
    fDeuteronRestMassNoTOF = new TH2F("fDeuteronRestMassNoTOF", "DeuteronNoTOF",
                                      72, 0.5, 8.05, 800, 0.00, 10.0);
    fDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronNoTOFList = fDeuteronTrackNoTOF->GetQAHists();
    fDeuteronNoTOFList->Add(fDeuteronRestMassNoTOF);
  }

  if (!fAntiDeuteronTrackNoTOF) {
    AliError("No AntiDeuteronNoTOF cuts \n");
  } else {
    fAntiDeuteronTrackNoTOF->Init();
    fAntiDeuteronRestMassNoTOF = new TH2F("fAntiDeuteronRestMassNoTOF",
                                          "AntiDeuteronNoTOF", 72, 0.5, 8.05,
                                          800, 0.00, 10.0);
    fAntiDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronNoTOFList = fAntiDeuteronTrackNoTOF->GetQAHists();
    fAntiDeuteronNoTOFList->Add(fAntiDeuteronRestMassNoTOF);
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
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fDumpster = new TList();
  fDumpster->SetName("Dumpster");
  fDumpster->SetOwner(kTRUE);

  if (fUseDumpster) {
    fProtonDeuteronDump = new AliFemtoDreamDump("pd");
    fDumpster->Add(fProtonDeuteronDump->GetOutput());

    fAntiProtonAntiDeuteronDump = new AliFemtoDreamDump("ApAd");
    fDumpster->Add(fAntiProtonAntiDeuteronDump->GetOutput());
  }

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
  PostData(6, fDeuteronNoTOFList);
  PostData(7, fAntiDeuteronNoTOFList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fDumpster);

  if (fProtonTrack->GetIsMonteCarlo()) {
    PostData(11, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    PostData(12, fAntiProtonMCList);
  }

  if (fDeuteronTrack->GetIsMonteCarlo()) {
    PostData(13, fDeuteronMCList);
  }
  if (fAntiDeuteronTrack->GetIsMonteCarlo()) {
    PostData(14, fAntiDeuteronMCList);
  }
}

//------------------------------------------UserExec()----------------------------------------------------------------------------

void AliAnalysisTaskNanoPt::UserExec(Option_t *option) {
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
  if (!fDeuteronTrackNoTOF || !fAntiDeuteronTrackNoTOF) {
    AliError("DeuteronNoTOF Cuts missing");
    return;
  }

  // EVENT SELECTION
  fEvent->SetEvent(fInputEvent);
  if (!fEvtCuts->isSelected(fEvent))
    return;

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack*>(fInputEvent->GetTrack(iTrack));
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

  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack*>(fInputEvent->GetTrack(iTrack));
    if (!track)
      continue;
    fTrack->SetTrack(track, fInputEvent, multiplicity);

    if (fIsMCTruth && fIsMC) {
      int mcpdg;
      mcpdg = fTrack->GetMCPDGCode();
      if ((mcpdg == 2212) && (fProtonTrack->isSelected(fTrack))) {
        Proton.push_back(*fTrack);
      }
      if ((mcpdg == -2212) && (fAntiProtonTrack->isSelected(fTrack))) {
        AntiProton.push_back(*fTrack);
      }
      if ((mcpdg == 1000010020) && (fDeuteronTrack->isSelected(fTrack))) {
        fDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        Deuteron.push_back(*fTrack);
      }
      if ((mcpdg == -1000010020)
          && (fAntiDeuteronTrack->isSelected(fTrack))) {
        fAntiDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        AntiDeuteron.push_back(*fTrack);
      }
      if ((mcpdg == 1000010020)
          && (fDeuteronTrackNoTOF->isSelected(fTrack))) {
        fDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      }
      if ((mcpdg == -1000010020)
          && (fAntiDeuteronTrackNoTOF->isSelected(fTrack))) {
        fAntiDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
      }
    } else {

      if (fProtonTrack->isSelected(fTrack)) {
        Proton.push_back(*fTrack);
      }
      if (fAntiProtonTrack->isSelected(fTrack)) {
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
      if (fDeuteronTrackNoTOF->isSelected(fTrack)) {
        fDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        // DeuteronNoTOF.push_back(*fTrack);
      }
      if (fAntiDeuteronTrackNoTOF->isSelected(fTrack)) {
        fAntiDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        //AntiDeuteronNoTOF.push_back(*fTrack);
      }
    }
  }
  //loop once over the MC stack to calculate Efficiency/Purity
  if (fIsMC) {
    AliAODInputHandler *eventHandler =
      dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()
                                        ->GetInputEventHandler());
    AliMCEvent* fMC = eventHandler->MCEvent();

    for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
      if (TMath::Abs(mcPart->Eta()) < 0.8 && mcPart->IsPhysicalPrimary()) {
        if (mcPart->GetPdgCode() == fProtonTrack->GetPDGCode()) {
          fProtonTrack->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiProtonTrack->GetPDGCode()) {
          fAntiProtonTrack->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiDeuteronTrack->GetPDGCode()) {
          std::cout << "deuterons " << mcPart->GetPdgCode() << std::endl;
          fAntiDeuteronTrack->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiDeuteronTrack->GetPDGCode()) {
          fAntiDeuteronTrack->FillGenerated(mcPart->Pt());
        }
      }
    }
  }
  fPairCleaner->ResetArray();

  fPairCleaner->CleanTrackAndDecay(&Proton, &Deuteron, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProton, &AntiDeuteron, 1);

  fPairCleaner->StoreParticle(Proton);
  fPairCleaner->StoreParticle(AntiProton);
  fPairCleaner->StoreParticle(Deuteron);
  fPairCleaner->StoreParticle(AntiDeuteron);

  if (fUseDumpster) {
    if (fProtonDeuteronDump) {
      fProtonDeuteronDump->SetEvent(Proton, Deuteron, fEvent, 2212, 1000010020);
    }
    if (fAntiProtonAntiDeuteronDump) {
      fAntiProtonAntiDeuteronDump->SetEvent(AntiProton, AntiDeuteron, fEvent,
                                            -2212, -1000010020);
    }
  }

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
// flush the data
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fDeuteronList);
  PostData(5, fAntiDeuteronList);
  PostData(6, fDeuteronNoTOFList);
  PostData(7, fAntiDeuteronNoTOFList);
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fDumpster);
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fProtonTrack->GetIsMonteCarlo()) {
    PostData(11, fProtonMCList);
  }
  if (fAntiProtonTrack->GetIsMonteCarlo()) {
    PostData(12, fAntiProtonMCList);
  }

  if (fDeuteronTrack->GetIsMonteCarlo()) {
    PostData(13, fDeuteronMCList);
  }
  if (fAntiDeuteronTrack->GetIsMonteCarlo()) {
    PostData(14, fAntiDeuteronMCList);
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
void AliAnalysisTaskNanoPt::StoreGlobalTrackReference(AliVTrack *track) {
// see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);
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
    if (dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

