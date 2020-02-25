/*
 * AliAnalysisTaskNanoBBar.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskNanoBBar.h"
#include "AliNanoAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

ClassImp(AliAnalysisTaskNanoBBar)
AliAnalysisTaskNanoBBar::AliAnalysisTaskNanoBBar()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fIsMC(false),
      fUseDumpster(false),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fProtonMCList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fAntiProtonMCList(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fLambdaMCList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fXiMCList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fAntiXiMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fProtonAntiProtonDump(nullptr),
      fProtonAntiLambdaDump(nullptr),
      fAntiProtonLambdaDump(nullptr),
      fLambdaAntiLambdaDump(nullptr),
      fProtonAntiXiDump(nullptr),
      fAntiProtonXiDump(nullptr),
      fLambdaAntiXiDump(nullptr),
      fAntiLambdaXiDump(nullptr),
      fXiAntiXiDump(nullptr),
      fDumpster(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskNanoBBar::AliAnalysisTaskNanoBBar(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fIsMC(false),
      fUseDumpster(false),
      fQA(nullptr),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fProtonMCList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fAntiProtonMCList(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fLambdaMCList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fAntiLambdaMCList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fXiMCList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fAntiXiMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fProtonAntiProtonDump(nullptr),
      fProtonAntiLambdaDump(nullptr),
      fAntiProtonLambdaDump(nullptr),
      fLambdaAntiLambdaDump(nullptr),
      fProtonAntiXiDump(nullptr),
      fAntiProtonXiDump(nullptr),
      fLambdaAntiXiDump(nullptr),
      fAntiLambdaXiDump(nullptr),
      fXiAntiXiDump(nullptr),
      fDumpster(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class());  //Output for the Event Cuts
  DefineOutput(3, TList::Class());  //Output for the Proton Cuts
  DefineOutput(4, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(5, TList::Class());  //Output for the Lambda Cuts
  DefineOutput(6, TList::Class());  //Output for the AntiLambda Cuts
  DefineOutput(7, TList::Class());  //Output for the Cascade Cuts
  DefineOutput(8, TList::Class());  //Output for the AntiCascade Cuts
  DefineOutput(9, TList::Class());  //Output for the Results
  DefineOutput(10, TList::Class());  //Output for the Results QA
  DefineOutput(11, TList::Class());  //Output for the Results Sample
  DefineOutput(12, TList::Class());  //Output for the Results Sample QA
  DefineOutput(13, TList::Class());  //Output for the Dumpster
  if (isMC) {
    DefineOutput(14, TList::Class());  //Output for the Track MC
    DefineOutput(15, TList::Class());  //Output for the Anti Track MC
    DefineOutput(16, TList::Class());  //Output for the V0 MC
    DefineOutput(17, TList::Class());  //Output for the Anti V0 MC
    DefineOutput(18, TList::Class());  //Output for the Xi MC
    DefineOutput(19, TList::Class());  //Output for the Anti Xi MC
  }
}

AliAnalysisTaskNanoBBar::~AliAnalysisTaskNanoBBar() {
  if (fEvent) {
    delete fEvent;
  }
  if (fEventCuts) {
    delete fEventCuts;
  }
  if (fTrack) {
    delete fTrack;
  }
  if (fProton) {
    delete fProton;
  }
  if (fAntiProton) {
    delete fAntiProton;
  }
  if (fv0) {
    delete fv0;
  }
  if (fLambda) {
    delete fLambda;
  }
  if (fAntiLambda) {
    delete fAntiLambda;
  }
    if (fCascade) {
    delete fCascade;
  }
  if (fXi) {
    delete fXi;
  }
  if (fAntiXi) {
    delete fAntiXi;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if (fSample) {
    delete fSample;
  }
  if (fProtonAntiProtonDump) {
    delete fProtonAntiProtonDump;
  }
  if (fProtonAntiLambdaDump) {
    delete fProtonAntiLambdaDump;
  }
  if (fAntiProtonLambdaDump) {
    delete fAntiProtonLambdaDump;
  }
  if (fLambdaAntiLambdaDump) {
    delete fLambdaAntiLambdaDump;
  }
  if (fProtonAntiXiDump) {
    delete fProtonAntiXiDump;
  }
  if (fAntiProtonXiDump) {
    delete fAntiProtonXiDump;
  }
  if (fLambdaAntiXiDump) {
    delete fLambdaAntiXiDump;
  }
  if (fAntiLambdaXiDump) {
    delete fAntiLambdaXiDump;
  }
  if (fXiAntiXiDump) {
    delete fXiAntiXiDump;
  }
  if (fDumpster){
    delete fDumpster;
  }
}

void AliAnalysisTaskNanoBBar::UserCreateOutputObjects() {
  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }
  if (!fProton) {
    AliError("No Proton cuts \n");
  } else {
    fProton->Init();
  }
  if (!fAntiProton) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProton->Init();
  }
  if (!fLambda) {
    AliError("No Lambda cuts \n");
  } else {
    fLambda->Init();
  }
  if (!fAntiLambda) {
    AliError("No AntiLambda cuts \n");
  } else {
    fAntiLambda->Init();
  }
 if (!fXi) {
    AliError("No Xi cuts \n");
  } else {
    fXi->Init();
  }
  if (!fAntiXi) {
    AliError("No AntiXi cuts \n");
  } else {
    fAntiXi->Init();
  }
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(8, 12,
                                                fConfig->GetMinimalBookingME());
    if (fConfig->GetUsePhiSpinning()) {
      fSample = new AliFemtoDreamControlSample(fConfig);
    }
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  fEvent->SetCalcSpherocity(fEventCuts->GetDoSpherocityCuts());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(
      fProton->GetIsMonteCarlo() || fAntiProton->GetIsMonteCarlo());

  fv0 = new AliFemtoDreamv0();
  fv0->SetUseMCInfo(
      fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());
  //PDG Codes should be set assuming Lambda0 to also work for AntiLambda
  fv0->SetPDGCode(3122);
  fv0->SetPDGDaughterPos(2212);
  fv0->GetPosDaughter()->SetUseMCInfo(
      fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());
  fv0->SetPDGDaughterNeg(211);
  fv0->GetNegDaughter()->SetUseMCInfo(
      fLambda->GetIsMonteCarlo() || fAntiLambda->GetIsMonteCarlo());

  fCascade = new AliFemtoDreamCascade();
  fCascade->SetUseMCInfo(fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  //PDG Codes should be set assuming Xi- to also work for Xi+
  fCascade->SetPDGCode(3312);
  fCascade->SetPDGDaugPos(2212);
  fCascade->GetPosDaug()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugNeg(211);
  fCascade->GetNegDaug()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugBach(211);
  fCascade->GetBach()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->Setv0PDGCode(3122);

  fQA = new TList();
  fQA->SetOwner();
  fQA->SetName("QA");
  fQA->Add(fEvent->GetEvtCutList());
  fQA->Add(fPairCleaner->GetHistList());

  fDumpster = new TList();
  fDumpster->SetName("Dumpster");
  fDumpster->SetOwner(kTRUE);

  if (fUseDumpster) {
  fProtonAntiProtonDump = new AliFemtoDreamDump("pAp");
  fProtonAntiProtonDump->SetkstarThreshold(0.3);
  fDumpster->Add(fProtonAntiProtonDump->GetOutput());

  fProtonAntiLambdaDump = new AliFemtoDreamDump("pAL");
  fProtonAntiLambdaDump->SetkstarThreshold(0.3);
  fDumpster->Add(fProtonAntiLambdaDump->GetOutput());

  fAntiProtonLambdaDump = new AliFemtoDreamDump("ApL");
  fAntiProtonLambdaDump->SetkstarThreshold(0.3);
  fDumpster->Add(fAntiProtonLambdaDump->GetOutput());

  fLambdaAntiLambdaDump = new AliFemtoDreamDump("LAL");
  fLambdaAntiLambdaDump->SetkstarThreshold(0.3);
  fDumpster->Add(fLambdaAntiLambdaDump->GetOutput());

  fProtonAntiXiDump = new AliFemtoDreamDump("pAXi");
  fProtonAntiXiDump->SetkstarThreshold(0.3);
  fDumpster->Add(fProtonAntiXiDump->GetOutput());

  fAntiProtonXiDump = new AliFemtoDreamDump("ApXi");
  fAntiProtonXiDump->SetkstarThreshold(0.3);
  fDumpster->Add(fAntiProtonXiDump->GetOutput());

  fLambdaAntiXiDump = new AliFemtoDreamDump("LAXi");
  fLambdaAntiXiDump->SetkstarThreshold(0.3);
  fDumpster->Add(fLambdaAntiXiDump->GetOutput());

  fAntiLambdaXiDump = new AliFemtoDreamDump("ALXi");
  fAntiLambdaXiDump->SetkstarThreshold(0.3);
  fDumpster->Add(fAntiLambdaXiDump->GetOutput());

  fXiAntiXiDump = new AliFemtoDreamDump("XiAXi");
  fXiAntiXiDump->SetkstarThreshold(0.3);
  fDumpster->Add(fXiAntiXiDump->GetOutput());
  }

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  fProtonList = fProton->GetQAHists();
  fAntiProtonList = fAntiProton->GetQAHists();
  fLambdaList = fLambda->GetQAHists();
  fAntiLambdaList = fAntiLambda->GetQAHists();
  fXiList = fXi->GetQAHists();
  fAntiXiList = fAntiXi->GetQAHists();

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

  fResultsSampleQA = new TList();
  fResultsSampleQA->SetOwner();
  fResultsSampleQA->SetName("ResultsSampleQA");

  if (fConfig->GetUsePhiSpinning()) {
    fResultsSample = fSample->GetHistList();

    if (!fConfig->GetMinimalBookingSample()) {
      fResultsSampleQA->Add(fSample->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResultsSample = new TList();
    fResultsSample->SetOwner();
    fResultsSample->SetName("ResultsSample");
  }

  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fProtonList);
  PostData(4, fAntiProtonList);
  PostData(5, fLambdaList);
  PostData(6, fAntiLambdaList);
  PostData(7, fXiList);
  PostData(8, fAntiXiList);
  PostData(9, fResults);
  PostData(10, fResultsQA);
  PostData(11, fResultsSample);
  PostData(12, fResultsSampleQA);
  PostData(13, fDumpster);

  if (fProton->GetIsMonteCarlo()) {
    if (!fProton->GetMinimalBooking()) {
      fProtonMCList = fProton->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCTrkCuts");
      fProtonMCList->SetOwner();
    }
    PostData(14, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    if (!fAntiProton->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProton->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiTrkCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(15, fAntiProtonMCList);
  }

  if (fLambda->GetIsMonteCarlo()) {
    if (!fLambda->GetMinimalBooking()) {
      fLambdaMCList = fLambda->GetMCQAHists();
    } else {
      fLambdaMCList = new TList();
      fLambdaMCList->SetName("MCv0Cuts");
      fLambdaMCList->SetOwner();
    }
    PostData(16, fLambdaMCList);
  }
  if (fAntiLambda->GetIsMonteCarlo()) {
    if (!fAntiLambda->GetMinimalBooking()) {
      fAntiLambdaMCList = fAntiLambda->GetMCQAHists();
    } else {
      fAntiLambdaMCList = new TList();
      fAntiLambdaMCList->SetName("MCAntiv0Cuts");
      fAntiLambdaMCList->SetOwner();
    }
    PostData(17, fAntiLambdaMCList);
  }
    if (fXi->GetIsMonteCarlo()) {
    if (!fXi->GetMinimalBooking()) {
      fXiMCList = fXi->GetMCQAHists();
    } else {
      fXiMCList = new TList();
      fXiMCList->SetName("MCXiCuts");
      fXiMCList->SetOwner();
    }
    PostData(18, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo()) {
    if (!fAntiXi->GetMinimalBooking()) {
      fAntiXiMCList = fAntiXi->GetMCQAHists();
    } else {
      fAntiXiMCList = new TList();
      fAntiXiMCList->SetName("MCAntiv0Cuts");
      fAntiXiMCList->SetOwner();
    }
    PostData(19, fAntiXiMCList);
  }
}

void AliAnalysisTaskNanoBBar::UserExec(Option_t *option) {
//  AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }

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
  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  const int multiplicity = fEvent->GetMultiplicity();
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent, multiplicity);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
  }

  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;
  AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);
  fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
      iv0 < static_cast<TClonesArray *>(aodEvt->GetV0s())->GetEntriesFast();
      ++iv0) {
    AliAODv0* casc = aodEvt->GetV0(iv0);
    fv0->Setv0(fInputEvent, casc, fEvent->GetMultiplicity());
    if (fLambda->isSelected(fv0)) {
      Lambdas.push_back(*fv0);
    }
    if (fAntiLambda->isSelected(fv0)) {
      AntiLambdas.push_back(*fv0);
    }
  }

  std::vector<AliFemtoDreamBasePart> Xis;
  std::vector<AliFemtoDreamBasePart> AntiXis;
  for (int iCasc = 0;
      iCasc
          < static_cast<TClonesArray *>(aodEvt->GetCascades())->GetEntriesFast();
      ++iCasc) {
    AliAODcascade* casc = aodEvt->GetCascade(iCasc);
    fCascade->SetCascade(fInputEvent, casc);
    if (fXi->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
    }
    if (fAntiXi->isSelected(fCascade)) {
      AntiXis.push_back(*fCascade);
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
        if (mcPart->GetPdgCode() == fProton->GetPDGCode()) {
          fProton->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiProton->GetPDGCode()) {
          fAntiProton->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fLambda->GetPDGv0()) {
          fLambda->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiLambda->GetPDGv0()) {
          fAntiLambda->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fXi->GetPDGv0()) {
          fXi->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiXi->GetPDGv0()) {
          fAntiXi->FillGenerated(mcPart->Pt());
        }
      }
    }
  }

  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Protons, &Lambdas, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiLambdas, 1);
  fPairCleaner->CleanTrackAndDecay(&Protons, &AntiLambdas, 2);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &Lambdas, 3);
  fPairCleaner->CleanTrackAndDecay(&Protons, &Xis, 4);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiXis, 5);
  fPairCleaner->CleanTrackAndDecay(&Protons, &AntiXis, 6);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &Xis, 7);


  fPairCleaner->CleanDecay(&Lambdas, 0);
  fPairCleaner->CleanDecay(&AntiLambdas, 1);
  fPairCleaner->CleanDecay(&Xis, 2);
  fPairCleaner->CleanDecay(&AntiXis, 3);
  fPairCleaner->CleanDecayAndDecay(&Lambdas, &AntiLambdas, 4);
  fPairCleaner->CleanDecayAndDecay(&Xis, &AntiXis, 5);
  fPairCleaner->CleanDecayAndDecay(&Lambdas, &Xis, 6);
  fPairCleaner->CleanDecayAndDecay(&AntiLambdas, &AntiXis, 7);
  fPairCleaner->CleanDecayAndDecay(&Lambdas, &AntiXis, 8);
  fPairCleaner->CleanDecayAndDecay(&AntiLambdas, &Xis, 9);


  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Lambdas);
  fPairCleaner->StoreParticle(AntiLambdas);
  fPairCleaner->StoreParticle(Xis);
  fPairCleaner->StoreParticle(AntiXis);
  if (fPairCleaner->GetCounter() > 0) {
    if (fConfig->GetUseEventMixing()) {
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent->GetZVertex(), fEvent->GetMultiplicity(),
                          fEvent->GetV0MCentrality());
    }
    if (fConfig->GetUsePhiSpinning()) {
      fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
    }
  }
    void SetEvent(std::vector<AliFemtoDreamBasePart> &vec1,
                std::vector<AliFemtoDreamBasePart> &vec2,
                AliFemtoDreamEvent *evt, const int pdg1, const int pdg2);
    if (fUseDumpster) {
    if(fProtonAntiProtonDump) {
    fProtonAntiProtonDump->SetEvent(Protons, AntiProtons, fEvent, 2212, -2212);
  }
    if(fProtonAntiLambdaDump) {
    fProtonAntiLambdaDump->SetEvent(Protons, AntiLambdas, fEvent, 2212, -3122);
  }
    if(fAntiProtonLambdaDump) {
    fAntiProtonLambdaDump->SetEvent(AntiProtons, Lambdas, fEvent, -2212, 3122);
  }
    if(fLambdaAntiLambdaDump) {
    fLambdaAntiLambdaDump->SetEvent(Lambdas, AntiLambdas, fEvent, 3122, -3122);
  }
    if(fProtonAntiXiDump) {
    fProtonAntiXiDump->SetEvent(Protons, AntiXis, fEvent, 2212, -3312);
  }
    if(fAntiProtonXiDump) {
    fAntiProtonXiDump->SetEvent(AntiProtons, Xis, fEvent, -2212, 3312);
  }
    if(fLambdaAntiXiDump) {
    fLambdaAntiXiDump->SetEvent(Lambdas, AntiXis, fEvent, 3122, -3312);
  }
    if(fAntiLambdaXiDump) {
    fAntiLambdaXiDump->SetEvent(AntiLambdas, Xis, fEvent, -3122, 3312);
  }
    if(fXiAntiXiDump) {
    fXiAntiXiDump->SetEvent(Xis, AntiXis, fEvent, 3312, -3312);
  }
    }


  PostData(1, fQA);
  PostData(2, fEvtList);
  PostData(3, fProtonList);
  PostData(4, fAntiProtonList);
  PostData(5, fLambdaList);
  PostData(6, fAntiLambdaList);
  PostData(7, fXiList);
  PostData(8, fAntiXiList);
  PostData(9, fResults);
  PostData(10, fResultsQA);
  PostData(11, fResultsSample);
  PostData(12, fResultsSampleQA);
  PostData(13, fDumpster);

  if (fProton->GetIsMonteCarlo()) {
    PostData(14, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    PostData(15, fAntiProtonMCList);
  }
  if (fLambda->GetIsMonteCarlo()) {
    PostData(16, fLambdaMCList);
  }
  if (fAntiLambda->GetIsMonteCarlo()) {
    PostData(17, fAntiLambdaMCList);
  }
  if (fXi->GetIsMonteCarlo()) {
    PostData(18, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo()) {
    PostData(19, fAntiXiMCList);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoBBar::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoBBar::StoreGlobalTrackReference(AliVTrack *track) {
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
