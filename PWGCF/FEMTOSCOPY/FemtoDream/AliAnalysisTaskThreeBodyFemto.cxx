/*
 * AliAnalysisTaskNanoXioton.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskThreeBodyFemto.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskThreeBodyFemto)
AliAnalysisTaskThreeBodyFemto::AliAnalysisTaskThreeBodyFemto()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
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
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsThreeBody(nullptr),  
      fRunThreeBody(true),
      sameEventDistributionPL(NULL),
      sameEventDistributionPPL(NULL),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskThreeBodyFemto::AliAnalysisTaskThreeBodyFemto(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
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
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsThreeBody(nullptr), 
      fRunThreeBody(true),
      sameEventDistributionPL(NULL),
      sameEventDistributionPPL(NULL),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
        DefineOutput(1, TList::Class());  //Output for the Event Cuts
        DefineOutput(2, TList::Class());  //Output for the Proton Cuts
        DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
        DefineOutput(4, TList::Class());  //Output for the Lambda Cuts
        DefineOutput(5, TList::Class());  //Output for the AntiLambda Cuts
        DefineOutput(6, TList::Class());  //Output for the Results
        DefineOutput(7, TList::Class());  //Output for the Results QA
        DefineOutput(8, TList::Class());  //Output for the Results
        DefineOutput(9, TList::Class());  //Output for the Results QA
        DefineOutput(10, TList::Class());  //Output for the Results Three body
        if (isMC) {
          DefineOutput(11, TList::Class());  //Output for the Track MC
          DefineOutput(12, TList::Class());  //Output for the Anti Track MC
          DefineOutput(13, TList::Class());  //Output for the V0 MC
          DefineOutput(14, TList::Class());  //Output for the Anti V0 MC
        }
      }

AliAnalysisTaskThreeBodyFemto::~AliAnalysisTaskThreeBodyFemto() {
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
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if (fSample) {
    delete fSample;
  }
}

void AliAnalysisTaskThreeBodyFemto::UserCreateOutputObjects() {
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
    AliError("No AntiXi cuts \n");
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
    if (fConfig->GetUsePhiSpinning()) {
      fSample = new AliFemtoDreamControlSample(fConfig);
    }
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), false);
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

  if (fRunThreeBody){
    fResultsThreeBody = new TList();
    fResultsThreeBody->SetOwner();
    fResultsThreeBody->SetName("ResultsTEST");

    sameEventDistributionPL = new TH1F("sameEventDistributionPL","sameEventDistributionPL",1000,0, 1);
    fResultsThreeBody->Add(sameEventDistributionPL);

    sameEventDistributionPPL = new TH1F("sameEventDistributionPPL","sameEventDistributionPPL",8000,0, 8);
    fResultsThreeBody->Add(sameEventDistributionPPL);

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

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fLambdaList);
  PostData(5, fAntiLambdaList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fResultsSample);
  PostData(9, fResultsSampleQA);
  PostData(10, fResultsThreeBody);
  if (fProton->GetIsMonteCarlo()) {
    if (!fProton->GetMinimalBooking()) {
      fProtonMCList = fProton->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCTrkCuts");
      fProtonMCList->SetOwner();
    }
    PostData(11, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    if (!fAntiProton->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProton->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiTrkCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(12, fAntiProtonMCList);
  }

  if (fLambda->GetIsMonteCarlo()) {
    if (!fLambda->GetMinimalBooking()) {
      fLambdaMCList = fLambda->GetMCQAHists();
    } else {
      fLambdaMCList = new TList();
      fLambdaMCList->SetName("MCv0Cuts");
      fLambdaMCList->SetOwner();
    }
    PostData(13, fLambdaMCList);
  }
  if (fAntiLambda->GetIsMonteCarlo()) {
    if (!fAntiLambda->GetMinimalBooking()) {
      fAntiLambdaMCList = fAntiLambda->GetMCQAHists();
    } else {
      fAntiLambdaMCList = new TList();
      fAntiLambdaMCList->SetName("MCAntiv0Cuts");
      fAntiLambdaMCList->SetOwner();
    }
    PostData(14, fAntiLambdaMCList);
  }

    

}

void AliAnalysisTaskThreeBodyFemto::UserExec(Option_t *option) {
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
  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Protons, &Lambdas, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiLambdas, 1);

  fPairCleaner->CleanDecay(&Lambdas, 0);
  fPairCleaner->CleanDecay(&AntiLambdas, 1);

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Lambdas);
  fPairCleaner->StoreParticle(AntiLambdas);

  
  if(fRunThreeBody){
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles();

    // Proton Lambda
    auto ProtonVector = ParticleVector.begin();
    auto LambdaVector = ParticleVector.begin()+2;
    // Loop over particles creating pairs
    for (auto iPart1 = ProtonVector->begin(); iPart1 != ProtonVector->end(); ++iPart1) {
      for (auto iPart2 = LambdaVector->begin(); iPart2 != LambdaVector->end(); ++iPart2) {
        // Set lorentz vectors
        TLorentzVector part1_LorVec, part2_LorVec;
        part1_LorVec.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
          iPart1->GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(2212)->Mass());
        part2_LorVec.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
          iPart2->GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(3122)->Mass());
        // Get momentum
        float RelativeK = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec, part2_LorVec);
        // No need to check pair selection because p lambda
        sameEventDistributionPL->Fill(RelativeK);
      }
    }
  
  
    // Proton Proton Lambda
    auto ProtonVectorPPL1 = ParticleVector.begin();
    auto ProtonVectorPPL2 = ParticleVector.begin();
    auto LambdaVectorPPL = ParticleVector.begin()+2;
    //loop over all lambdas
    for (auto iPart1 = LambdaVectorPPL->begin(); iPart1 != LambdaVectorPPL->end(); ++iPart1) {
      // loop over all protons in the first vector - lambda and proton will be never counted twice
      
      for (auto iPart2 = ProtonVectorPPL1->begin(); iPart2 != ProtonVectorPPL1->end(); ++iPart2) {
        // loop over antiprotons from second vector from antiproton1 +1 to the end

        for (auto iPart3 = iPart2+1; iPart3 != ProtonVectorPPL2->end(); ++iPart3) {

          // From now on: iPart1 is a lambda, iPart2 is a proton, iPart3 is a proton
          // Now we have the three particles, lets create their Lorentz vectors
          auto massProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
          auto massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
          TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
          part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
          iPart1->GetMomentum().Z(), sqrt(iPart1->GetMomentum().X()*iPart1->GetMomentum().X()+iPart1->GetMomentum().Y()*iPart1->GetMomentum().Y() + iPart1->GetMomentum().Z()*iPart1->GetMomentum().Z()+massLambda*massLambda));
          part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
          iPart2->GetMomentum().Z(), sqrt(iPart2->GetMomentum().X()*iPart2->GetMomentum().X()+iPart2->GetMomentum().Y()*iPart2->GetMomentum().Y() + iPart2->GetMomentum().Z()*iPart2->GetMomentum().Z()+massProton*massProton));
          part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
          iPart3->GetMomentum().Z(), sqrt(iPart3->GetMomentum().X()*iPart3->GetMomentum().X()+iPart3->GetMomentum().Y()*iPart3->GetMomentum().Y() + iPart3->GetMomentum().Z()*iPart3->GetMomentum().Z()+massProton*massProton) );
          
          // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
          TLorentzVector q12 = RelativePairMomentum(part1_LorVec,part2_LorVec);
          TLorentzVector q23 = RelativePairMomentum(part2_LorVec,part3_LorVec);
          TLorentzVector q31 = RelativePairMomentum(part3_LorVec,part1_LorVec);
          // The particles in current methodology are put in bins of:
          //                 Q3=sqrt(q12^2+q23^2+q31^2)

          float Q32 = q12*q12+q23*q23+q31*q31;
          // From 3 pion paper, the q must be multiplied by -1 before taking quare root
          float Q3 = sqrt(-Q32); 
          sameEventDistributionPPL->Fill(Q3);
          
        }
      }
    }
  }
  

  if (fPairCleaner->GetCounter() > 0) {
    if (fConfig->GetUseEventMixing()) {
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent);
    }
    if (fConfig->GetUsePhiSpinning()) {
      fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
    }
  }
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fLambdaList);
  PostData(5, fAntiLambdaList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fResultsSample);
  PostData(9, fResultsSampleQA);
  PostData(10, fResultsThreeBody);
  if (fProton->GetIsMonteCarlo()) {
    PostData(11, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    PostData(12, fAntiProtonMCList);
  }
  if (fLambda->GetIsMonteCarlo()) {
    PostData(13, fLambdaMCList);
  }
  if (fAntiLambda->GetIsMonteCarlo()) {
    PostData(14, fAntiLambdaMCList);
  }

}

//____________________________________________________________________________________________________
void AliAnalysisTaskThreeBodyFemto::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskThreeBodyFemto::StoreGlobalTrackReference(AliVTrack *track) {
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


TLorentzVector AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(
    TLorentzVector &PartOne, TLorentzVector &PartTwo) {
  // The q12 components can be calculated as:
  //             q^mu = (p1-p2)^mu /2 - ((p1-p2)*P/(2P^2))*P^mu
  //             where P = p1+p2
  // Reference: https://www.annualreviews.org/doi/pdf/10.1146/annurev.nucl.55.090704.151533
  // In the following code the above written equation will be expressed as:
  //             q = trackDifference/2 -  scaling * trackSum
  // where scaling is a float number:
  //             scaling = trackDifference*trackSum/(2*trackSum^2) = ((p1-p2)*P/(2P^2))
  TLorentzVector trackSum = PartOne + PartTwo;
  TLorentzVector trackDifference = PartOne - PartTwo;
  float scaling = trackDifference*trackSum/(2*trackSum*trackSum);
  TLorentzVector qPartOnePartTwo = trackDifference*0.5 -  scaling * trackSum;

  return qPartOnePartTwo;
}
