/*
 * AliAnalysisTaskFemtoDreamDeuteron.cxx
 *
 *  Created on: 21 Mar 2018
 *      Author: bernhardhohlweger
 */

#include "AliAnalysisTaskFemtoDreamDeuteron.h"
#include "AliFemtoDreamBasePart.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliVEvent.h"
ClassImp(AliAnalysisTaskFemtoDreamDeuteron)
AliAnalysisTaskFemtoDreamDeuteron::AliAnalysisTaskFemtoDreamDeuteron()
  : AliAnalysisTaskSE(),
    fIsMC(false),
    fIsMCTruth(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteronDCA(nullptr),
    fTrackCutsDeuteronMass(nullptr),
    fTrackCutsAntiDeuteronDCA(nullptr),
    fTrackCutsAntiDeuteronMass(nullptr),
    fTrackCutsProtonDCA(nullptr),
    fTrackCutsAntiProtonDCA(nullptr),
    fConfig(nullptr),
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
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fProtonDeuteronDump(nullptr),
    fAntiProtonAntiDeuteronDump(nullptr),
    fDumpster(nullptr),
    fUseDumpster(false),
    fTrackBufferSize() {
}

AliAnalysisTaskFemtoDreamDeuteron::AliAnalysisTaskFemtoDreamDeuteron(
  const char *name, bool isMC)
  : AliAnalysisTaskSE(name),
    fIsMC(isMC),
    fIsMCTruth(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteronDCA(nullptr),
    fTrackCutsDeuteronMass(nullptr),
    fTrackCutsAntiDeuteronDCA(nullptr),
    fTrackCutsAntiDeuteronMass(nullptr),
    fTrackCutsProtonDCA(nullptr),
    fTrackCutsAntiProtonDCA(nullptr),
    fConfig(nullptr),
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
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fProtonDeuteronDump(nullptr),
    fAntiProtonAntiDeuteronDump(nullptr),
    fDumpster(nullptr),
    fUseDumpster(false),
    fTrackBufferSize(2000) {
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

AliAnalysisTaskFemtoDreamDeuteron::~AliAnalysisTaskFemtoDreamDeuteron() {
  delete fEvent;
  delete fTrack;
  delete fTrackCutsDeuteronDCA;
  delete fTrackCutsDeuteronMass;
  delete fTrackCutsAntiDeuteronDCA;
  delete fTrackCutsAntiDeuteronMass;
  delete fTrackCutsProtonDCA;
  delete fTrackCutsAntiProtonDCA;
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

Float_t AliAnalysisTaskFemtoDreamDeuteron::GetMass2sq(
  AliFemtoDreamTrack *track) {
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (!(beta < 0)) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

void AliAnalysisTaskFemtoDreamDeuteron::UserCreateOutputObjects() {

  fGTI = new AliAODTrack*[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }

  if (!fTrackCutsProtonDCA) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsProtonDCA->Init();
    fProtonList = fTrackCutsProtonDCA->GetQAHists();
    if (fIsMC) {
      fProtonMCList = fTrackCutsProtonDCA->GetMCQAHists();
    }
  }

  if (!fTrackCutsAntiProtonDCA) {
    AliError("No AntiProton cuts \n");
  } else {
    fTrackCutsAntiProtonDCA->Init();
    fAntiProtonList = fTrackCutsAntiProtonDCA->GetQAHists();
    if (fIsMC) {
      fAntiProtonMCList = fTrackCutsAntiProtonDCA->GetMCQAHists();
    }
  }

  if (!fTrackCutsDeuteronDCA) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsDeuteronDCA->Init();
    fDeuteronRestMass = new TH2F("fDeuteronRestMass", "Deuteron", 72, 0.5, 8.05,
                                 800, 0.00, 10.0);
    fDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronList = fTrackCutsDeuteronDCA->GetQAHists();
    if (fIsMC) {
      fDeuteronMCList = fTrackCutsDeuteronDCA->GetMCQAHists();
    }
    fDeuteronList->Add(fDeuteronRestMass);
  }

  if (!fTrackCutsAntiDeuteronDCA) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsAntiDeuteronDCA->Init();
    fAntiDeuteronRestMass = new TH2F("fAntiDeuteronRestMass", "AntiDeuteron", 72,
                                     0.5, 8.05, 800, 0.00, 10.0);
    fAntiDeuteronRestMass->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMass->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronList = fTrackCutsAntiDeuteronDCA->GetQAHists();
    if (fIsMC) {
      fAntiDeuteronMCList = fTrackCutsAntiDeuteronDCA->GetMCQAHists();
    }
    fAntiDeuteronList->Add(fAntiDeuteronRestMass);
  }
//--------------------------------------------------------------------------------------------------------------------
  if (!fTrackCutsDeuteronMass) {
    AliError("No DeuteronNoTOF cuts \n");
  } else {
    fTrackCutsDeuteronMass->Init();
    fDeuteronRestMassNoTOF = new TH2F("fDeuteronRestMassNoTOF", "DeuteronNoTOF",
                                      72, 0.5, 8.05, 800, 0.00, 10.0);
    fDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fDeuteronNoTOFList = fTrackCutsDeuteronMass->GetQAHists();
    fDeuteronNoTOFList->Add(fDeuteronRestMassNoTOF);
  }

  if (!fTrackCutsAntiDeuteronMass) {
    AliError("No AntiDeuteronNoTOF cuts \n");
  } else {
    fTrackCutsAntiDeuteronMass->Init();
    fAntiDeuteronRestMassNoTOF = new TH2F("fAntiDeuteronRestMassNoTOF",
                                          "AntiDeuteronNoTOF", 72, 0.5, 8.05,
                                          800, 0.00, 10.0);
    fAntiDeuteronRestMassNoTOF->GetXaxis()->SetTitle("pT(GeV)");
    fAntiDeuteronRestMassNoTOF->GetYaxis()->SetTitle("m^2(Gev)^2");
    fAntiDeuteronNoTOFList = fTrackCutsAntiDeuteronMass->GetQAHists();
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

  fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates());
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

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
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

  if (fTrackCutsProtonDCA->GetIsMonteCarlo()) {
    PostData(11, fProtonMCList);
  }
  if (fTrackCutsAntiProtonDCA->GetIsMonteCarlo()) {
    PostData(12, fAntiProtonMCList);
  }

  if (fTrackCutsDeuteronDCA->GetIsMonteCarlo()) {
    PostData(13, fDeuteronMCList);
  }
  if (fTrackCutsAntiDeuteronDCA->GetIsMonteCarlo()) {
    PostData(14, fAntiDeuteronMCList);
  }
}

void AliAnalysisTaskFemtoDreamDeuteron::UserExec(Option_t*) {
  AliAODEvent *Event = static_cast<AliAODEvent*>(fInputEvent);

  if (!Event) {
    AliWarning("No Input Event");
  } else {

    fEvent->SetEvent(Event);

    if (fEventCuts->isSelected(fEvent)) {

      ResetGlobalTrackReference();

      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        StoreGlobalTrackReference(track);
      }

      fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

      static std::vector<AliFemtoDreamBasePart> DCADeuterons;
      DCADeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> DCAAntiDeuterons;
      DCAAntiDeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> DCAProtons;
      DCAProtons.clear();
      static std::vector<AliFemtoDreamBasePart> DCAAntiProtons;
      DCAAntiProtons.clear();


      //Now we loop over all the tracks in the reconstructed event.
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }

        fTrack->SetTrack(track);
        if (fIsMCTruth && fIsMC) {
          int mcpdg;
          mcpdg = fTrack->GetMCPDGCode();
          if ((mcpdg == 2212) && (fTrackCutsProtonDCA->isSelected(fTrack))) {
            DCAProtons.push_back(*fTrack);
          }
          if ((mcpdg == -2212) && (fTrackCutsAntiProtonDCA->isSelected(fTrack))) {
            DCAAntiProtons.push_back(*fTrack);
          }
          if ((mcpdg == 1000010020) && (fTrackCutsDeuteronDCA->isSelected(fTrack))) {
            fDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            DCADeuterons.push_back(*fTrack);
          }
          if ((mcpdg == -1000010020)
              && (fTrackCutsAntiDeuteronDCA->isSelected(fTrack))) {
            fAntiDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            DCAAntiDeuterons.push_back(*fTrack);
          }
          if ((mcpdg == 1000010020)
              && (fTrackCutsDeuteronMass->isSelected(fTrack))) {
            fDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          }
          if ((mcpdg == -1000010020)
              && (fTrackCutsAntiDeuteronMass->isSelected(fTrack))) {
            fAntiDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          }

        } else {
          if (fTrackCutsDeuteronDCA->isSelected(fTrack)) {
            DCADeuterons.push_back(*fTrack);
            fDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          }
          if (fTrackCutsDeuteronMass->isSelected(fTrack)) {
            fDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          }
          if (fTrackCutsAntiDeuteronDCA->isSelected(fTrack)) {
            DCAAntiDeuterons.push_back(*fTrack);
            fAntiDeuteronRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          }
          if (fTrackCutsAntiDeuteronMass->isSelected(fTrack)) {
            fAntiDeuteronRestMassNoTOF->Fill(fTrack->GetPt(),
                                             GetMass2sq(fTrack));
          }
          if (fTrackCutsProtonDCA->isSelected(fTrack)) {
            DCAProtons.push_back(*fTrack);
          }
          if (fTrackCutsAntiProtonDCA->isSelected(fTrack)) {
            DCAAntiProtons.push_back(*fTrack);
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
            if (mcPart->GetPdgCode() == fTrackCutsProtonDCA->GetPDGCode()) {
              fTrackCutsProtonDCA->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsAntiProtonDCA->GetPDGCode()) {
              fTrackCutsAntiProtonDCA->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsDeuteronDCA->GetPDGCode()) {
              std::cout << "deuterons " << mcPart->GetPdgCode() << std::endl;
              fTrackCutsDeuteronDCA->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsAntiDeuteronDCA->GetPDGCode()) {
              fTrackCutsAntiDeuteronDCA->FillGenerated(mcPart->Pt());
            }
          }
        }
      }
      fPairCleaner->CleanTrackAndDecay(&DCAProtons, &DCADeuterons, 0);
      fPairCleaner->CleanTrackAndDecay(&DCAAntiProtons, &DCAAntiDeuterons, 1);
      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(DCAProtons);
      fPairCleaner->StoreParticle(DCAAntiProtons);
      fPairCleaner->StoreParticle(DCADeuterons);
      fPairCleaner->StoreParticle(DCAAntiDeuterons);
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent->GetZVertex(), fEvent->GetRefMult08(),
                          fEvent->GetV0MCentrality());
      void SetEvent(std::vector<AliFemtoDreamBasePart> &vec1,
                    std::vector<AliFemtoDreamBasePart> &vec2,
                    AliFemtoDreamEvent * evt, const int pdg1, const int pdg2);

      if (fUseDumpster) {
        if (fProtonDeuteronDump) {
          fProtonDeuteronDump->SetEvent(DCAProtons, DCADeuterons, fEvent, 2212,
                                        1000010020);
        }

        if (fAntiProtonAntiDeuteronDump) {
          fAntiProtonAntiDeuteronDump->SetEvent(DCAAntiProtons,
                                                DCAAntiDeuterons, fEvent, -2212,
                                                -1000010020);
        }
      }
    }
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
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fTrackCutsProtonDCA->GetIsMonteCarlo()) {
    PostData(11, fProtonMCList);
  }
  if (fTrackCutsAntiProtonDCA->GetIsMonteCarlo()) {
    PostData(12, fAntiProtonMCList);
  }
  if (fTrackCutsDeuteronDCA->GetIsMonteCarlo()) {
    PostData(13, fDeuteronMCList);
  }
  if (fTrackCutsAntiDeuteronDCA->GetIsMonteCarlo()) {
    PostData(14, fAntiDeuteronMCList);
  }
}

void AliAnalysisTaskFemtoDreamDeuteron::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskFemtoDreamDeuteron::StoreGlobalTrackReference(
  AliAODTrack *track) {

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

    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }

    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}
