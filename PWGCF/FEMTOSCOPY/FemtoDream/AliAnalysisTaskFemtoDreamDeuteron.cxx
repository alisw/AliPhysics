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
    fisLightWeight(false),
    fTrackBufferSize(),
    fIsMC(false),
    fdoSideband(false),
    fSigmaUp(0.0),
    fSigmaLow(0.0),
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
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
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
    fResults(nullptr),
    fResultsQA(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr){
}

AliAnalysisTaskFemtoDreamDeuteron::AliAnalysisTaskFemtoDreamDeuteron(
  const char *name, bool isMC)
  : AliAnalysisTaskSE(name),
    fisLightWeight(false),
    fTrackBufferSize(2000),
    fIsMC(isMC),
    fdoSideband(false),
    fSigmaUp(0.0),
    fSigmaLow(0.0),
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
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
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
    fResults(nullptr),
    fResultsQA(nullptr),
    fDeuteronRestMass(nullptr),
    fAntiDeuteronRestMass(nullptr),
    fDeuteronRestMassNoTOF(nullptr),
    fAntiDeuteronRestMassNoTOF(nullptr){
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Dueteron Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiDeuteron Cuts
  DefineOutput(6, TList::Class());  //Output for the DueteronNoTOF Cuts
  DefineOutput(7, TList::Class());  //Output for the AntiDeuteronNoTOF Cuts
  DefineOutput(8, TList::Class());  //Output for the Results
  DefineOutput(9, TList::Class());  //Output for the Results QA
  if (fIsMC) {
    DefineOutput(10, TList::Class());  //Output for the Proton MC
    DefineOutput(11, TList::Class());  //Output for the AntiProton MC
    DefineOutput(12, TList::Class());  //Output for the Deuteron MC
    DefineOutput(13, TList::Class());  //Output for the AntiDeuteron MC
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
}

Float_t AliAnalysisTaskFemtoDreamDeuteron::GetMass2sq(
  AliFemtoDreamTrack *track)const{
  Float_t p = track->GetP();
  Float_t mass2sq = -999;
  Float_t beta = track->GetbetaTOF();
  if (!(beta < 0)) {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

float AliAnalysisTaskFemtoDreamDeuteron::MeanTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const{
  float pTVal = track->GetPt();
  float par0 =  3.55375e+00;
  float par1 = -1.25749e+00;
  float par2 = -3.60444e-01;
  float par3 = -1.00250e-01;
  float par4 = -1.00782e-02;
  return par0 + TMath::Exp(par1 * pTVal + par2 * pTVal * pTVal+ par3 * pTVal * pTVal * pTVal+ par4 * pTVal* pTVal * pTVal * pTVal);
};
float AliAnalysisTaskFemtoDreamDeuteron::SigmaTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const{
  float pTVal = track->GetPt();
  Float_t par0 = 1.19287e-02;
  Float_t par1 = 0.202460e-02;
  Float_t par2 = 1.23058e-02;//par[2];
  Float_t par3 = 30.23644e-04;
  Float_t par4 = 45.80006e-05;
  return 0.088 + 0.1*(par0 * pTVal + par1 * pTVal * pTVal + par2 * pTVal * pTVal* pTVal+ par3 * pTVal * pTVal* pTVal* pTVal+ par4 * pTVal * pTVal* pTVal* pTVal* pTVal);
};

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
    fDeuteronList->Add(fDeuteronRestMass);
    if (fIsMC) {
      fDeuteronMCList = fTrackCutsDeuteronDCA->GetMCQAHists();
    }
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
    fAntiDeuteronList->Add(fAntiDeuteronRestMass);
    if (fIsMC) {
      fAntiDeuteronMCList = fTrackCutsAntiDeuteronDCA->GetMCQAHists();
    }
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

  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

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

  if (fTrackCutsProtonDCA->GetIsMonteCarlo()) {
    PostData(10, fProtonMCList);
  }
  if (fTrackCutsAntiProtonDCA->GetIsMonteCarlo()) {
    PostData(11, fAntiProtonMCList);
  }

  if (fTrackCutsDeuteronDCA->GetIsMonteCarlo()) {
    PostData(12, fDeuteronMCList);
  }
  if (fTrackCutsAntiDeuteronDCA->GetIsMonteCarlo()) {
    PostData(13, fAntiDeuteronMCList);
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
        if (fTrackCutsProtonDCA->isSelected(fTrack)) {
          DCAProtons.push_back(*fTrack);
        }
        if (fTrackCutsAntiProtonDCA->isSelected(fTrack)) {
          DCAAntiProtons.push_back(*fTrack);
        }
        if (fTrackCutsDeuteronMass->isSelected(fTrack)) {
          fDeuteronRestMassNoTOF->Fill(fTrack->GetPt(), GetMass2sq(fTrack));
        }
        if (fTrackCutsAntiDeuteronMass->isSelected(fTrack)) {
          fAntiDeuteronRestMassNoTOF->Fill(fTrack->GetPt(),GetMass2sq(fTrack));
        }

        if (fTrackCutsDeuteronDCA->isSelected(fTrack)){
          float MassSqaured = GetMass2sq(fTrack);
          if(fdoSideband){
            float meanMass = MeanTOFMassSqdDeuteron(fTrack);
            float sigmaMass = SigmaTOFMassSqdDeuteron(fTrack);
            float upMass = meanMass+ (fSigmaUp* sigmaMass);
            float LowMass = meanMass+ (fSigmaLow*sigmaMass);
            if((MassSqaured>= LowMass)&&(MassSqaured<=upMass)) {
              DCADeuterons.push_back(*fTrack);
              fDeuteronRestMass->Fill(fTrack->GetPt(), MassSqaured);
            }
          }else{
            DCADeuterons.push_back(*fTrack);
            fDeuteronRestMass->Fill(fTrack->GetPt(), MassSqaured);
          }
        }
        if (fTrackCutsAntiDeuteronDCA->isSelected(fTrack)){
          float MassSqaured = GetMass2sq(fTrack);
          if(fdoSideband){
            float meanMass = MeanTOFMassSqdDeuteron(fTrack);
            float sigmaMass = SigmaTOFMassSqdDeuteron(fTrack);
            float upMass = meanMass+fSigmaUp*sigmaMass;
            float LowMass = meanMass+fSigmaLow*sigmaMass;
            if((MassSqaured >= LowMass)&&(MassSqaured<=upMass)) {
              DCAAntiDeuterons.push_back(*fTrack);
              fAntiDeuteronRestMass->Fill(fTrack->GetPt(), MassSqaured);
            }
          }else{
            DCAAntiDeuterons.push_back(*fTrack);
            fAntiDeuteronRestMass->Fill(fTrack->GetPt(), MassSqaured);
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
          if (mcPart->IsPhysicalPrimary()) {
            if (mcPart->GetPdgCode() == fTrackCutsProtonDCA->GetPDGCode()) {
              fTrackCutsProtonDCA->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsAntiProtonDCA->GetPDGCode()) {
              fTrackCutsAntiProtonDCA->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsDeuteronDCA->GetPDGCode()) {
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
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fTrackCutsProtonDCA->GetIsMonteCarlo()) {
    PostData(10, fProtonMCList);
  }
  if (fTrackCutsAntiProtonDCA->GetIsMonteCarlo()) {
    PostData(11, fAntiProtonMCList);
  }
  if (fTrackCutsDeuteronDCA->GetIsMonteCarlo()) {
    PostData(12, fDeuteronMCList);
  }
  if (fTrackCutsAntiDeuteronDCA->GetIsMonteCarlo()) {
    PostData(13, fAntiDeuteronMCList);
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
