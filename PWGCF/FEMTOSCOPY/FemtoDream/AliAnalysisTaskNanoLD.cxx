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
      fPairCleanerSettings(1),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fSimpleEventCounter(nullptr),
      fSimpleParticleCounter(nullptr),
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
      fProton(nullptr),
      fProtonList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
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
      fPairCleanerSettings(1),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fSimpleEventCounter(nullptr),
      fSimpleParticleCounter(nullptr),
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
      fProton(nullptr),
      fProtonList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
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
  DefineOutput(8, TList::Class());  //Output for the Proton Cuts
  DefineOutput(9, TList::Class());  //Output for the AntiProton Cuts
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
  if (fProton) {
    delete fProton;
  }
  if (fAntiProton) {
    delete fAntiProton;
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

  // Set number of histograms for clean track and decay
  Int_t nHistCleanTrackDecay = 2;
  Int_t nHistCleanDecayDecay = 2;
  if (fPairCleanerSettings == 0) {
    nHistCleanTrackDecay = 0;
    nHistCleanDecayDecay = 0;
  } else if (fPairCleanerSettings == 1) {
    nHistCleanTrackDecay = 2;
    nHistCleanDecayDecay = 2;
  } else if (fPairCleanerSettings == 2) {
    nHistCleanTrackDecay = 4;
    nHistCleanDecayDecay = 2;
  }

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
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(nHistCleanTrackDecay,
                                                nHistCleanDecayDecay,
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

  // Simple event and particle counters, always active independent of fullBlastQA
  fSimpleEventCounter = new TH1F("fSimpleEventCounter", "Simple event counter", 1, 0., 1.);
  fSimpleEventCounter->GetYaxis()->SetTitle("Number of events");
  fEvtList->Add(fSimpleEventCounter);

  fSimpleParticleCounter = new TH1F("fSimpleParticleCounter", "Simple particle counter", 6, 0., 6.);
  fSimpleParticleCounter->GetYaxis()->SetTitle("Number of particles");
  fSimpleParticleCounter->GetXaxis()->SetBinLabel(1,"Deuterons");
  fSimpleParticleCounter->GetXaxis()->SetBinLabel(2,"AntiDeuterons");
  fSimpleParticleCounter->GetXaxis()->SetBinLabel(3,"Lambdas");
  fSimpleParticleCounter->GetXaxis()->SetBinLabel(4,"Antilambdas");
  fSimpleParticleCounter->GetXaxis()->SetBinLabel(5,"Protons");
  fSimpleParticleCounter->GetXaxis()->SetBinLabel(6,"AntiProtons");
  fSimpleParticleCounter->GetXaxis()->SetTitle(0);
  fEvtList->Add(fSimpleParticleCounter);

  fDeuteronList = fDeuteron->GetQAHists();
  fAntiDeuteronList = fAntiDeuteron->GetQAHists();
  fLambdaList = fLambda->GetQAHists();
  fAntiLambdaList = fAntiLambda->GetQAHists();
  fProtonList = fProton->GetQAHists();
  fAntiProtonList = fAntiProton->GetQAHists();

  // Deuteron and antideuteron TOF mass squared plots
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
  PostData(8, fProtonList);
  PostData(9, fAntiProtonList);
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

  // Fill simple event counter
  fSimpleEventCounter->Fill(0.5);

  // Get and store global track reference
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }

  // Deuteron and Anti-Deuteron selection
  std::vector<AliFemtoDreamBasePart> Deuterons;
  std::vector<AliFemtoDreamBasePart> AntiDeuterons;
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent);
    fTrack->SetInvMass(1.87561); // PDG value, cannot be otained from TDatabasePDG
                                 // in case of deuterons, therefore hard coded here
    if (fDeuteron->isSelected(fTrack)) {
      Deuterons.push_back(*fTrack);
      fSimpleParticleCounter->Fill(0.);
      fDeuteronMassSqTOF->Fill(fTrack->GetPt(), CalculateMassSqTOF(fTrack));
    }
    if (fAntiDeuteron->isSelected(fTrack)) {
      AntiDeuterons.push_back(*fTrack);
      fSimpleParticleCounter->Fill(1.);
      fAntiDeuteronMassSqTOF->Fill(fTrack->GetPt(), CalculateMassSqTOF(fTrack));
    }
  }

  // Lambda and Anti-Lambda selection
  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;
  AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);
  fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
      iv0 < static_cast<TClonesArray *>(aodEvt->GetV0s())->GetEntriesFast();
      ++iv0) {
    AliAODv0* casc = aodEvt->GetV0(iv0);
    fv0->Setv0(fInputEvent, casc);
    if (fLambda->isSelected(fv0)) {
      Lambdas.push_back(*fv0);
      fSimpleParticleCounter->Fill(2.);
    }
    if (fAntiLambda->isSelected(fv0)) {
      AntiLambdas.push_back(*fv0);
      fSimpleParticleCounter->Fill(3.);
    }
  }

  // Proton and Anti-Proton selection (only for pair cleaner)
  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
      fSimpleParticleCounter->Fill(4.);
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
      fSimpleParticleCounter->Fill(5.);
    }
  }

  // Pair cleaner
  fPairCleaner->ResetArray();

  // Clean deuterons and lambda daughters (default = activated)
  if (fPairCleanerSettings > 0) {
    fPairCleaner->CleanTrackAndDecay(&Deuterons, &Lambdas, 0);
    fPairCleaner->CleanTrackAndDecay(&AntiDeuterons, &AntiLambdas, 1);
  }

  // Clean protons and lambda daughters in case this is activated
  if (fPairCleanerSettings == 2) {
    fPairCleaner->CleanTrackAndDecay(&Protons, &Lambdas, 2);
    fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiLambdas, 3);
  }

  // Clean lambdas and lambdas (default = activated)
  if (fPairCleanerSettings > 0) {
    fPairCleaner->CleanDecay(&Lambdas, 0);
    fPairCleaner->CleanDecay(&AntiLambdas, 1);
  }

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
  PostData(8, fProtonList);
  PostData(9, fAntiProtonList);
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
