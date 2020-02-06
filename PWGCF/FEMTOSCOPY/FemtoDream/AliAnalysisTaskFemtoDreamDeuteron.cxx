/*
 * AliAnalysisTaskFemtoDreamDeuteron.cxx
 *
 *  Created on: 21 Mar 2018
 *      Author: bernhardhohlweger
 */

#include "AliAnalysisTaskFemtoDreamDeuteron.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
ClassImp(AliAnalysisTaskFemtoDreamDeuteron)
AliAnalysisTaskFemtoDreamDeuteron::AliAnalysisTaskFemtoDreamDeuteron()
  : AliAnalysisTaskSE(),
    fIsMC(false),
    fOutput(nullptr),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteronDCA(nullptr),
    fTrackCutsDeuteronMass(nullptr),
    fTrackCutsAntiDeuteronDCA(nullptr),
    fTrackCutsAntiDeuteronMass(nullptr),
    fTrackCutsProtonDCA(nullptr),
    fTrackCutsProtonMass(nullptr),
    fTrackCutsAntiProtonDCA(nullptr),
    fTrackCutsAntiProtonMass(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fDRestMass(nullptr),
    fAntiDRestMass(nullptr),
    fPRestMass(nullptr),
    fAntiPRestMass(nullptr),
    fProtonRestMassMC(nullptr),
    fAntiProtonRestMassMC(nullptr),
    fDeuteronRestMassMC(nullptr),
    fAntiDeuteronRestMassMC(nullptr),
    fKaonRestMassMC(nullptr),
    fAntiKaonRestMassMC(nullptr),
    fDProtonRestMassMC(nullptr),
    fDKaonRestMassMC(nullptr),
    fAntiDProtonRestMassMC (nullptr),
    fAntiDKaonRestMassMC(nullptr),
    fPionRestMassMC(nullptr),
    fAntiPionRestMassMC(nullptr),
    fDPionRestMassMC(nullptr),
    fAntiDPionRestMassMC(nullptr),
    fProtonBackgroundMC(nullptr),
    fAntiProtonBackgroundMC(nullptr),
    fDeuteronBackgroundMC(nullptr),
    fAntiDeuteronBackgroundMC(nullptr),

    fTrackBufferSize() {

}

AliAnalysisTaskFemtoDreamDeuteron::AliAnalysisTaskFemtoDreamDeuteron(const char *name, bool isMC)
  : AliAnalysisTaskSE(name),
    fIsMC(isMC),
    fOutput(nullptr),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsDeuteronDCA(nullptr),
    fTrackCutsDeuteronMass(nullptr),
    fTrackCutsAntiDeuteronDCA(nullptr),
    fTrackCutsAntiDeuteronMass(nullptr),
    fTrackCutsProtonDCA(nullptr),
    fTrackCutsProtonMass(nullptr),
    fTrackCutsAntiProtonDCA(nullptr),
    fTrackCutsAntiProtonMass(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fDRestMass(nullptr),
    fAntiDRestMass(nullptr),
    fPRestMass(nullptr),
    fAntiPRestMass(nullptr),
    fProtonRestMassMC(nullptr),
    fAntiProtonRestMassMC(nullptr),
    fDeuteronRestMassMC(nullptr),
    fAntiDeuteronRestMassMC(nullptr),
    fKaonRestMassMC(nullptr),
    fAntiKaonRestMassMC(nullptr),
    fDProtonRestMassMC(nullptr),
    fDKaonRestMassMC(nullptr),
    fAntiDProtonRestMassMC (nullptr),
    fAntiDKaonRestMassMC(nullptr),
    fPionRestMassMC(nullptr),
    fAntiPionRestMassMC(nullptr),
    fDPionRestMassMC(nullptr),
    fAntiDPionRestMassMC(nullptr),
    fProtonBackgroundMC(nullptr),
    fAntiProtonBackgroundMC(nullptr),
    fDeuteronBackgroundMC(nullptr),
    fAntiDeuteronBackgroundMC(nullptr),

    fTrackBufferSize(2000) {
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskFemtoDreamDeuteron::~AliAnalysisTaskFemtoDreamDeuteron() {
  // TODO Auto-generated destructor stub
}

Float_t AliAnalysisTaskFemtoDreamDeuteron::GetMass2sq(AliFemtoDreamTrack *track) {
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (!(beta < 0)) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

void AliAnalysisTaskFemtoDreamDeuteron::InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName) {
  if (!trkCuts) {
    // If the track cuts didn't arrive here, we can go home
    AliFatal("Track Cuts not set!");
  }

  trkCuts->Init();

  trkCuts->SetName(trkCutsName.Data());
  fOutput->Add(trkCuts->GetQAHists());

  if (trkCuts->GetIsMonteCarlo()) {
    trkCuts->SetMCName(MCName.Data());
    fOutput->Add(trkCuts->GetMCQAHists());
  }
}

void AliAnalysisTaskFemtoDreamDeuteron::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output"); // Every output objects needs a name, be careful names can collide!
  fOutput->SetOwner();        // This tells ROOT that this list belongs to the top list / top object


  fDRestMass = new TH2F("fDRestMass", "Deuteron_RestMsq", 72, 0.5, 8.05, 800, 0.00, 10.0);
  fDRestMass->GetXaxis()->SetTitle("pT");
  fDRestMass->GetYaxis()->SetTitle("m^2");
  fAntiDRestMass = new TH2F("fAntiDRestMass", "Antideuteron_RestMsq", 72, 0.5, 8.05, 800, 0.00, 10.0);
  fAntiDRestMass->GetXaxis()->SetTitle("pT");
  fAntiDRestMass->GetYaxis()->SetTitle("m^2");
  fPRestMass = new TH2F("fPRestMass", "Proton_RestMsq", 36, 0.5, 4.05, 400, 0.00, 3);
  fPRestMass->GetXaxis()->SetTitle("pT");
  fPRestMass->GetYaxis()->SetTitle("m^2");
  fAntiPRestMass = new TH2F("fAntiPRestMass", "Antiproton_RestMsq", 36, 0.5, 4.05, 400, 0.00, 3);
  fAntiPRestMass->GetXaxis()->SetTitle("pT");
  fAntiPRestMass->GetYaxis()->SetTitle("m^2");
  fOutput->Add(fDRestMass);
  fOutput->Add(fAntiDRestMass);
  fOutput->Add(fPRestMass);
  fOutput->Add(fAntiPRestMass);

//===================================================================== MC HistogramsfProtonRestMassMCfKaonRestMassMC
  if (fIsMC) {
    if (!fTrackCutsProtonMass) {
      AliError("No Proton cuts \n");
    } else {
      //fProtonTrackNoTOF->Init();
      fProtonRestMassMC = new TH2F("fProtonRestMassMC", "Proton", 36, 0.5, 4.05, 400, 0.0, 3);
      fProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fKaonRestMassMC = new TH2F("fKaonRestMassMC", "kaon", 36, 0.5, 4.05, 400, 0.0, 3);
      fKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fPionRestMassMC = new TH2F("fPionRestMassMC", "pion", 36, 0.5, 4.05, 400, 0.0, 3);
      fPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fProtonBackgroundMC = new TH2F("fProtonBackgroungMC", "background", 36, 0.5, 4.05, 400, 0.0, 3);
      fProtonBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fProtonBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fOutput->Add(fProtonRestMassMC);
      fOutput->Add(fKaonRestMassMC);
      fOutput->Add(fPionRestMassMC);
      fOutput->Add(fProtonBackgroundMC);
    }

    if (!fTrackCutsAntiProtonMass) {
      AliError("No AntiProton cuts \n");
    } else {
      //fAntiProtonTrackNoTOF->Init();
      fAntiProtonRestMassMC = new TH2F("fAntiProtonRestMassMC", "AntiProton", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiKaonRestMassMC = new TH2F("fAntiKaonRestMassMC", "AntiProton", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiPionRestMassMC = new TH2F("fAntiPionRestMassMC", "AntiPion", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiProtonBackgroundMC = new TH2F("fAntiProtonBackgroundMC", "AntiProtonBackgroundMC", 36, 0.5, 4.05, 400, 0.00, 3);
      fAntiProtonBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiProtonBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fOutput->Add(fAntiProtonRestMassMC);
      fOutput->Add(fAntiKaonRestMassMC);
      fOutput->Add(fAntiPionRestMassMC);
      fOutput->Add(fAntiProtonBackgroundMC);
    }

    if (!fTrackCutsDeuteronMass) {
      AliError("No Proton cuts \n");
    } else {
      //fDeuteronTrackNoTOF->Init();
      fDeuteronRestMassMC = new TH2F("fDeuteronRestMassMC", "Deuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDeuteronRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDeuteronRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDProtonRestMassMC = new TH2F("fDProtonRestMassMC", "Proton", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDKaonRestMassMC = new TH2F("fDKaonRestMassMC", "Kaon", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDPionRestMassMC = new TH2F("fDPionRestMassMC", "Pion", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fDPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fDeuteronBackgroundMC = new TH2F("fDeuteronBackgroundMC", "Pion", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fDeuteronBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fDeuteronBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fOutput->Add(fDeuteronRestMassMC);
      fOutput->Add(fDProtonRestMassMC);
      fOutput->Add(fDKaonRestMassMC);
      fOutput->Add(fDPionRestMassMC);
      fOutput->Add(fDeuteronBackgroundMC);
    }

    if (!fTrackCutsAntiDeuteronMass) {
      AliError("No Proton cuts \n");
    } else {
      //fAntiDeuteronTrackNoTOF->Init();
      fAntiDeuteronRestMassMC = new TH2F("fAntiDeuteronRestMassMC", "AntiDeuteron", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDeuteronRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDeuteronRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDProtonRestMassMC = new TH2F("fAntiDProtonRestMassMC", "AntiProton", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDProtonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDProtonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDKaonRestMassMC = new TH2F("fAntiDKaonRestMassMC", "AntiKaon", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDKaonRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDKaonRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDPionRestMassMC = new TH2F("fAntiDPionRestMassMC", "AntiPion", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDPionRestMassMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDPionRestMassMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fAntiDeuteronBackgroundMC = new TH2F("fAntiDeuteronBackgroundMC", "AntiDeuteronBackgroundMC", 72, 0.5, 8.05, 800, 0.00, 10.0);
      fAntiDeuteronBackgroundMC->GetXaxis()->SetTitle("pT(GeV)");
      fAntiDeuteronBackgroundMC->GetYaxis()->SetTitle("m^2(Gev)^2");
      fOutput->Add(fAntiDeuteronRestMassMC);
      fOutput->Add(fAntiDProtonRestMassMC);
      fOutput->Add(fAntiDKaonRestMassMC);
      fOutput->Add(fAntiDPionRestMassMC);
      fOutput->Add(fAntiDeuteronBackgroundMC);
    }

  }


  fEvent = new AliFemtoDreamEvent(false, true, GetCollisionCandidates());
  fOutput->Add(fEvent->GetEvtCutList());
  //Nothing special about the Femto Track, we just initialize it
  fTrack = new AliFemtoDreamTrack();
  //If this is false, the MC information is not included
  fTrack->SetUseMCInfo(fIsMC);

  fGTI = new AliAODTrack*[fTrackBufferSize];

  if (fEventCuts) {
    fEventCuts->InitQA();
    if (fEventCuts->GetHistList()) {
      fOutput->Add(fEventCuts->GetHistList());
    }
  } else {
    AliWarning("Event cuts are missing! \n");
  }

  InitHistograms(fTrackCutsDeuteronDCA, "DCADeuterons", "MCDCADeuterons");
  InitHistograms(fTrackCutsDeuteronMass, "MassDeuterons", "MCMassDeuterons");
  InitHistograms(fTrackCutsAntiDeuteronDCA, "DCAAntiDeuterons", "MCDCAAntiDeuterons");
  InitHistograms(fTrackCutsAntiDeuteronMass, "MassAntiDeuterons", "MCMassAntiDeuterons");
  InitHistograms(fTrackCutsProtonDCA, "DCAProtons", "MCDCAProtons");
  InitHistograms(fTrackCutsProtonMass, "MassProtons", "MCMassProtons");
  InitHistograms(fTrackCutsAntiProtonDCA, "DCAAntiProtons", "MCDCAAntiProtons");
  InitHistograms(fTrackCutsAntiProtonMass, "MassAntiProtons", "MCMassAntiProtons");

  fPairCleaner = new AliFemtoDreamPairCleaner(2, 4, false);
  fOutput->Add(fPairCleaner->GetHistList());
  fPartColl = new AliFemtoDreamPartCollection(fConfig, false);

  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());
  PostData(1, fOutput);
}



void AliAnalysisTaskFemtoDreamDeuteron::UserExec(Option_t *) {
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
      static std::vector<AliFemtoDreamBasePart> MassDeuterons;
      MassDeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> DCAAntiDeuterons;
      DCAAntiDeuterons.clear();
      static std::vector<AliFemtoDreamBasePart> MassAntiDeuterons;
      MassAntiDeuterons.clear();

      static std::vector<AliFemtoDreamBasePart> DCAProtons;
      DCAProtons.clear();
      static std::vector<AliFemtoDreamBasePart> MassProtons;
      MassProtons.clear();
      static std::vector<AliFemtoDreamBasePart> DCAAntiProtons;
      DCAAntiProtons.clear();
      static std::vector<AliFemtoDreamBasePart> MassAntiProtons;
      MassAntiProtons.clear();
      //Now we loop over all the tracks in the reconstructed event.
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        //First we fill all the track information into the track object
        fTrack->SetTrack(track);
        //then we pass it to the track cut objects, which checks if it passes our selection criteria
        //for particle 1 and if it returns true ...
        if (fTrackCutsDeuteronDCA->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          DCADeuterons.push_back(*fTrack);
        }
        if (fTrackCutsDeuteronMass->isSelected(fTrack)) {
          MassDeuterons.push_back(*fTrack);
          fDRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          if (fIsMC) {

            if (fTrack->GetMCPDGCode() == 1000010020) {
              fDeuteronRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            } else if (fTrack->GetMCPDGCode() == 2212) {
              fDProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            } else if (fTrack->GetMCPDGCode() == 321 ) {
              fDKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            } else if (fTrack->GetMCPDGCode() == 211 ) {
              fDPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            } else {
              fDeuteronBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            }
          }
        }

        if (fTrackCutsAntiDeuteronDCA->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          DCAAntiDeuterons.push_back(*fTrack);
        }

        if (fTrackCutsAntiDeuteronMass->isSelected(fTrack)) {
          MassAntiDeuterons.push_back(*fTrack);
          fAntiDRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          if (fIsMC) {

            if (fTrack->GetMCPDGCode() == -1000010020) {

              fAntiDeuteronRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == -2212) {

              fAntiDProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == -321 ) {

              fAntiDKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == -211 ) {

              fAntiDPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else {
              fAntiDeuteronBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            }
          }

        }

        if (fTrackCutsProtonDCA->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          fTrack->SetCPA(gRandom->Uniform());
          DCAProtons.push_back(*fTrack);
        }
        if (fTrackCutsProtonMass->isSelected(fTrack)) {
          MassProtons.push_back(*fTrack);
          fPRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          if (fIsMC) {

            if (fTrack->GetMCPDGCode() == 2212) {

              fProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == 321 ) {

              fKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == 211 ) {

              fPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else {
              fProtonBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
            }
          }
        }

        if (fTrackCutsAntiProtonDCA->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          DCAAntiProtons.push_back(*fTrack);
        }

        if (fTrackCutsAntiProtonMass->isSelected(fTrack)) {
          MassAntiProtons.push_back(*fTrack);
          fAntiPRestMass->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
          if (fIsMC) {
            if (fTrack->GetMCPDGCode() == -2212) {

              fAntiProtonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == -321 ) {

              fAntiKaonRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else if (fTrack->GetMCPDGCode() == -211 ) {

              fAntiPionRestMassMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            } else {
              fAntiProtonBackgroundMC->Fill(fTrack->GetPt(), GetMass2sq(fTrack));

            }
          }
        }
      }


      fPairCleaner->CleanTrackAndDecay(&DCAProtons, &DCADeuterons, 0);
      fPairCleaner->CleanTrackAndDecay(&DCAAntiProtons, &DCAAntiDeuterons, 1);
      fPairCleaner->CleanDecay(&DCAProtons, 0);
      fPairCleaner->CleanDecay(&DCAAntiProtons, 1);
      fPairCleaner->CleanDecay(&DCADeuterons, 2);
      fPairCleaner->CleanDecay(&DCAAntiDeuterons, 3);
      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(DCAProtons);
      fPairCleaner->StoreParticle(DCAAntiProtons);
      fPairCleaner->StoreParticle(DCADeuterons);
      fPairCleaner->StoreParticle(DCAAntiDeuterons);
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                          fEvent->GetRefMult08(), fEvent->GetV0MCentrality());
    }
  }
  PostData(1, fOutput);
}

void AliAnalysisTaskFemtoDreamDeuteron::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}
void AliAnalysisTaskFemtoDreamDeuteron::StoreGlobalTrackReference(AliAODTrack *track) {

  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }

  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
           , trackID, fTrackBufferSize);
    return;
  }

  if (fGTI[trackID]) {

    if ( (!track->GetFilterMap()) && (!track->GetTPCNcls()) ) {
      return;
    }

    if ( fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()  ) {
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
