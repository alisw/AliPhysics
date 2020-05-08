/*
 * AliAnalysisTaskNanoLD.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 *  Modifications for LambdaDeuteron
 *      Author: stheckel
 */

#include "AliAnalysisTaskNanoLD.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskNanoLD)
AliAnalysisTaskNanoLD::AliAnalysisTaskNanoLD()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fDeuteron(nullptr),
      fDeuteronList(nullptr),
      fDeuteronMassSqTOF(nullptr),
      fAntiDeuteron(nullptr),
      fAntiDeuteronList(nullptr),
      fAntiDeuteronMassSqTOF(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskNanoLD::AliAnalysisTaskNanoLD(const char* name)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fDeuteron(nullptr),
      fDeuteronList(nullptr),
      fDeuteronMassSqTOF(nullptr),
      fAntiDeuteron(nullptr),
      fAntiDeuteronList(nullptr),
      fAntiDeuteronMassSqTOF(nullptr),
      fv0(nullptr),
      fLambda(nullptr),
      fLambdaList(nullptr),
      fAntiLambda(nullptr),
      fAntiLambdaList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Deuteron Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiDeuteron Cuts
  DefineOutput(4, TList::Class());  //Output for the Lambda Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiLambda Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA
}

AliAnalysisTaskNanoLD::~AliAnalysisTaskNanoLD() {
  if (fEvent) {
    delete fEvent;
  }
  if (fEventCuts) {
    delete fEventCuts;
  }
  if (fTrack) {
    delete fTrack;
  }
  if (fDeuteron) {
    delete fDeuteron;
  }
  if (fAntiDeuteron) {
    delete fAntiDeuteron;
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
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoLD::UserCreateOutputObjects() {
  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }
  if (!fDeuteron) {
    AliError("No Deuteron cuts \n");
  } else {
    fDeuteron->Init();
  }
  if (!fAntiDeuteron) {
    AliError("No AntiDeuteron cuts \n");
  } else {
    fAntiDeuteron->Init();
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
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2,
                                                fConfig->GetMinimalBookingME());
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), false);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(false);
  
  fv0 = new AliFemtoDreamv0();
  fv0->SetUseMCInfo(false);
  //PDG Codes should be set assuming Lambda0 to also work for AntiLambda
  fv0->SetPDGCode(3122);
  fv0->SetPDGDaughterPos(2212);
  fv0->GetPosDaughter()->SetUseMCInfo(false);
  fv0->SetPDGDaughterNeg(211);
  fv0->GetNegDaughter()->SetUseMCInfo(false);

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  fDeuteronList = fDeuteron->GetQAHists();
  fAntiDeuteronList = fAntiDeuteron->GetQAHists();
  fLambdaList = fLambda->GetQAHists();
  fAntiLambdaList = fAntiLambda->GetQAHists();

  // Mass squared plots
  fDeuteronMassSqTOF = new TH2F("fDeuteronMassSqTOF", "Deuterons", 50, 0. ,5., 400, 0., 8.);
  fDeuteronMassSqTOF->GetXaxis()->SetTitle("p_T (GeV/c)");
  fDeuteronMassSqTOF->GetYaxis()->SetTitle("m^2 (GeV/c^2)^2");
  fDeuteronList->Add(fDeuteronMassSqTOF);

  fAntiDeuteronMassSqTOF = new TH2F("fAntiDeuteronMassSqTOF", "AntiDeuterons", 50, 0. ,5., 400, 0., 8.);
  fAntiDeuteronMassSqTOF->GetXaxis()->SetTitle("p_T (GeV/c)");
  fAntiDeuteronMassSqTOF->GetYaxis()->SetTitle("m^2 (GeV/c^2)^2");
  fAntiDeuteronList->Add(fAntiDeuteronMassSqTOF);

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
  PostData(2, fDeuteronList);
  PostData(3, fAntiDeuteronList);
  PostData(4, fLambdaList);
  PostData(5, fAntiLambdaList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoLD::UserExec(Option_t *option) {
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }

  // DEUTERON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Deuterons;
  std::vector<AliFemtoDreamBasePart> AntiDeuterons;
  const int multiplicity = fEvent->GetMultiplicity();
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent, multiplicity);
    if (fDeuteron->isSelected(fTrack)) {
      Deuterons.push_back(*fTrack);
      fDeuteronMassSqTOF->Fill(fTrack->GetPt(), CalculateMassSqTOF(fTrack));
    }
    if (fAntiDeuteron->isSelected(fTrack)) {
      AntiDeuterons.push_back(*fTrack);
      fAntiDeuteronMassSqTOF->Fill(fTrack->GetPt(), CalculateMassSqTOF(fTrack));
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

  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Deuterons, &Lambdas, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiDeuterons, &AntiLambdas, 1);

  fPairCleaner->CleanDecay(&Lambdas, 0);
  fPairCleaner->CleanDecay(&AntiLambdas, 1);

  fPairCleaner->StoreParticle(Deuterons);
  fPairCleaner->StoreParticle(AntiDeuterons);
  fPairCleaner->StoreParticle(Lambdas);
  fPairCleaner->StoreParticle(AntiLambdas);

  if (fPairCleaner->GetCounter() > 0) {
    if (fConfig->GetUseEventMixing()) {
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent->GetZVertex(), fEvent->GetMultiplicity(),
                          fEvent->GetV0MCentrality());
    }
  }

  PostData(1, fEvtList);
  PostData(2, fDeuteronList);
  PostData(3, fAntiDeuteronList);
  PostData(4, fLambdaList);
  PostData(5, fAntiLambdaList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskNanoLD::CalculateMassSqTOF(AliFemtoDreamTrack *track) {
  // Calculate the mass squared from TOF
  Float_t p = track->GetP();
  Float_t beta = track->GetbetaTOF();
  Float_t massSq = -999;

  if (beta > 0.) {
    massSq = ((1 / (beta * beta)) - 1) * (p * p);
  }

  //printf("p = %f - beta = %f, massSq = %f \n",p,beta,massSq);

  return massSq;
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoLD::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskNanoLD::StoreGlobalTrackReference(AliVTrack *track) {
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
