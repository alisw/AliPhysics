/*
 * AliAnalysisTaskThreeBodyFemtoPionProton.cxx
 *
 *  Created on: May 13, 2019
 *      Author: Laura Serksnyte 
 */
#include "AliAnalysisTaskThreeBodyFemtoPionProton.h"
#include "AliFemtoDreamHigherPairMath.h"

ClassImp(AliAnalysisTaskThreeBodyFemtoPionProton)
AliAnalysisTaskThreeBodyFemtoPionProton::AliAnalysisTaskThreeBodyFemtoPionProton()
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
      fTrackPion(nullptr),
      fPion(nullptr),
      fPionList(nullptr),
      fPionMCList(nullptr),
      fAntiPion(nullptr),
      fAntiPionList(nullptr),
      fAntiPionMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsThreeBody(nullptr),  
      fSameEvent(nullptr),  
      fMixedEvent(nullptr), 
      fSameEventMult(nullptr), 
      fMixedEventMult(nullptr), 
      fSameEventPhiTheta(nullptr), 
      fMixedEventPhiTheta(nullptr), 
      fRejectedParticles(nullptr),
      fAcceptedParticles(nullptr),
      fOtherHistos(nullptr),    
      fRunThreeBody(true),
      fRejectedParticlesArray(nullptr),
      fAcceptedParticlesArray(nullptr),
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fPartContainer(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fTriggerOn(false),
      fIsMC(false),
      fQ3Limit(0.0),
      fAllEvents(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskThreeBodyFemtoPionProton::AliAnalysisTaskThreeBodyFemtoPionProton(const char* name, bool isMC, bool triggerOn)
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
      fTrackPion(nullptr),
      fPion(nullptr),
      fPionList(nullptr),
      fPionMCList(nullptr),
      fAntiPion(nullptr),
      fAntiPionList(nullptr),
      fAntiPionMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsThreeBody(nullptr),  
      fSameEvent(nullptr),  
      fMixedEvent(nullptr), 
      fSameEventMult(nullptr), 
      fMixedEventMult(nullptr), 
      fSameEventPhiTheta(nullptr), 
      fMixedEventPhiTheta(nullptr), 
      fRejectedParticles(nullptr),
      fAcceptedParticles(nullptr),
      fOtherHistos(nullptr),    
      fRunThreeBody(true),
      fRejectedParticlesArray(nullptr),
      fAcceptedParticlesArray(nullptr),
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fPartContainer(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fTriggerOn(false),
      fIsMC(false),
      fQ3Limit(0.0),
      fAllEvents(nullptr),
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

AliAnalysisTaskThreeBodyFemtoPionProton::~AliAnalysisTaskThreeBodyFemtoPionProton() {
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
  if (fTrackPion) {
    delete fTrackPion;
  }
  if (fPion) {
    delete fPion;
  }
  if (fAntiPion) {
    delete fAntiPion;
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

void AliAnalysisTaskThreeBodyFemtoPionProton::UserCreateOutputObjects() {
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
  if (!fPion) {
    AliError("No Pion cuts \n");
  } else {
    fPion->Init();
  }
  if (!fAntiPion) {
    AliError("No AntiPion cuts \n");
  } else {
    fAntiPion->Init();
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

  float DummyFloat = fEventCuts->GetlowerPtBoundSpherCalc();
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), false, DummyFloat);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  fEvent->SetCalcSpherocity(fEventCuts->GetDoSpherocityCuts());
  fEvent->SetfLowPtSpherCalc(DummyFloat); //quick fix to set the right value of the pT Cut


  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(
      fProton->GetIsMonteCarlo() || fAntiProton->GetIsMonteCarlo());

  fTrackPion = new AliFemtoDreamTrack();
  fTrackPion->SetUseMCInfo(
      fPion->GetIsMonteCarlo() || fAntiPion->GetIsMonteCarlo());

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }
  

  fProtonList = fProton->GetQAHists();
  fAntiProtonList = fAntiProton->GetQAHists();
  fPionList = fPion->GetQAHists();
  fAntiPionList = fAntiPion->GetQAHists();
  

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
    
    
    fOtherHistos = new TList();
    fOtherHistos->SetOwner();
    fOtherHistos->SetName("OtherHistos");

    fAllEvents =  new TH1F("fAllEvents","fAllEvents", 25, 0, 100);

    fOtherHistos->Add(fAllEvents);


    // Rejected
    fRejectedParticles = new TList();
    fRejectedParticles->SetOwner();
    fRejectedParticles->SetName("RejectedParticles");

    fRejectedParticlesArray = new TH1F*[7];
    TString histRejected[7] ={"ProtonPionAntiPion","ProtonPionPion", "AntiProtonPionPion", 
      "ProtonProtonPion", "ProtonProtonAntiPion","CaseProtonPionPion1", "CaseProtonProtonPion2"};
    for (int i = 0; i < 7; ++i) {
      fRejectedParticlesArray[i] = new TH1F(histRejected[i],histRejected[i], 25, 0, 100);
      fRejectedParticles->Add(fRejectedParticlesArray[i]);
     }

    // Accepted
    fAcceptedParticles = new TList();
    fAcceptedParticles->SetOwner();
    fAcceptedParticles->SetName("AcceptedParticles");

    fAcceptedParticlesArray = new TH1F*[7];
    TString histAccepted[7] ={"ProtonPionAntiPion","ProtonPionPion", "AntiProtonPionPion", 
      "ProtonProtonPion", "ProtonProtonAntiPion","CaseProtonPionPion1", "CaseProtonProtonPion2"};
    for (int i = 0; i < 7; ++i) {
      fAcceptedParticlesArray[i] = new TH1F(histAccepted[i],histAccepted[i], 25, 0, 100);
      fAcceptedParticles->Add(fAcceptedParticlesArray[i]);
     }


    // Same event
   // Same event
    fSameEvent = new TList();
    fSameEvent->SetOwner();
    fSameEvent->SetName("SameEvent");

    fSameEventTripletArray = new TH1F*[32];
    TString histTitlesSame[32] ={"ProtonPionAntiPion","AntiProtonPionAntiPion","ProtonPionPion",
      "AntiProtonPionPion", "ProtonAntiPionAntiPion", "AntiProtonAntiPionAntiPion",
      "ProtonProtonPion", "AntiProtonAntiProtonAntiPion",
      "ProtonProtonAntiPion", "AntiProtonAntiProtonPion", // after this are mixed/same
      "ProtonPionSameAntiPionMixed","ProtonMixedPionAntiPionSame","ProtonAntiPionSamePionMixed", 
      "AntiProtonPionSameAntiPionMixed","AntiProtonMixedPionAntiPionSame","AntiProtonAntiPionSamePionMixed",
      "ProtonPionSamePionMixed", "ProtonMixedPionPionSame", 
      "AntiProtonPionSamePionMixed", "AntiProtonMixedPionPionSame", 
      "ProtonAntiPionSameAntiPionMixed", "ProtonMixedAntiPionAntiPionSame",
      "AntiProtonAntiPionSameAntiPionMixed", "AntiProtonMixedAntiPionAntiPionSame", 
      "ProtonProtonSamePionMixed", "ProtonMixedProtonPionSame",
      "AntiProtonAntiProtonSameAntiPionMixed", "AntiProtonMixedAntiProtonAntiPionSame",
      "ProtonProtonSameAntiPionMixed","ProtonMixedProtonAntiPionSame",
      "AntiProtonAntiProtonSamePionMixed", "AntiProtonMixedAntiProtonPionSame"};
    for (int i = 0; i < 32; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }
    // Same event multiplicity dist
    fSameEventMult = new TList();
    fSameEventMult->SetOwner();
    fSameEventMult->SetName("SameEventMult");

    fSameEventTripletMultArray = new TH2F*[32];
    TString histTitlesSameMult[32] ={"ProtonPionAntiPion","AntiProtonPionAntiPion","ProtonPionPion",
      "AntiProtonPionPion", "ProtonAntiPionAntiPion", "AntiProtonAntiPionAntiPion",
      "ProtonProtonPion", "AntiProtonAntiProtonAntiPion",
      "ProtonProtonAntiPion", "AntiProtonAntiProtonPion", // after this are mixed/same
      "ProtonPionSameAntiPionMixed","ProtonMixedPionAntiPionSame","ProtonAntiPionSamePionMixed", 
      "AntiProtonPionSameAntiPionMixed","AntiProtonMixedPionAntiPionSame","AntiProtonAntiPionSamePionMixed",
      "ProtonPionSamePionMixed", "ProtonMixedPionPionSame", 
      "AntiProtonPionSamePionMixed", "AntiProtonMixedPionPionSame", 
      "ProtonAntiPionSameAntiPionMixed", "ProtonMixedAntiPionAntiPionSame",
      "AntiProtonAntiPionSameAntiPionMixed", "AntiProtonMixedAntiPionAntiPionSame", 
      "ProtonProtonSamePionMixed", "ProtonMixedProtonPionSame",
      "AntiProtonAntiProtonSameAntiPionMixed", "AntiProtonMixedAntiProtonAntiPionSame",
      "ProtonProtonSameAntiPionMixed","ProtonMixedProtonAntiPionSame",
      "AntiProtonAntiProtonSamePionMixed", "AntiProtonMixedAntiProtonPionSame"};
    
    for (int i = 0; i < 32; ++i) {
      fSameEventTripletMultArray[i] =  new TH2F(histTitlesSameMult[i],histTitlesSameMult[i], 8000, 0, 8,26,1,27);
      fSameEventMult->Add(fSameEventTripletMultArray[i]);
     }

    // Mixed event
    fMixedEvent = new TList();
    fMixedEvent->SetOwner();
    fMixedEvent->SetName("MixedEvent");

    fMixedEventTripletArray = new TH1F*[10];
    TString histTitlesMixed[10] ={"ProtonPionAntiPion","AntiProtonPionAntiPion","ProtonPionPion",
      "AntiProtonPionPion", "ProtonAntiPionAntiPion", "AntiProtonAntiPionAntiPion",
      "ProtonProtonPion", "AntiProtonAntiProtonAntiPion",
      "ProtonProtonAntiPion", "AntiProtonAntiProtonPion"};
    for (int i = 0; i < 10; ++i) {
      fMixedEventTripletArray[i] = new TH1F(histTitlesMixed[i],histTitlesMixed[i], 8000, 0, 8);
      fMixedEvent->Add(fMixedEventTripletArray[i]);
     }
    // Mixed event multiplicity dist
    fMixedEventMult = new TList();
    fMixedEventMult->SetOwner();
    fMixedEventMult->SetName("MixedEventMult");

    fMixedEventTripletMultArray = new TH2F*[10];
    TString histTitlesMixedMult[10] ={"ProtonPionAntiPion","AntiProtonPionAntiPion","ProtonPionPion",
      "AntiProtonPionPion", "ProtonAntiPionAntiPion", "AntiProtonAntiPionAntiPion",
      "ProtonProtonPion", "AntiProtonAntiProtonAntiPion",
      "ProtonProtonAntiPion", "AntiProtonAntiProtonPion"};
    for (int i = 0; i < 10; ++i) {
      fMixedEventTripletMultArray[i] = new TH2F(histTitlesMixedMult[i],histTitlesMixedMult[i], 8000, 0, 8,26,1,27);
      fMixedEventMult->Add(fMixedEventTripletMultArray[i]);
     }

    // Same event phi theta distribution
    fSameEventPhiTheta = new TList();
    fSameEventPhiTheta->SetOwner();
    fSameEventPhiTheta->SetName("SameEventPhiTheta");

    fSameEventTripletPhiThetaArray = new TH2F*[66];
    TString histTitlesSamePhiEta[33] = {"ProtonPionAntiPion","AntiProtonPionAntiPion","ProtonPionPion",
      "AntiProtonPionPion", "ProtonAntiPionAntiPion", "AntiProtonAntiPionAntiPion",
      "ProtonProtonPion", "AntiProtonAntiProtonAntiPion",
      "ProtonProtonAntiPion", "AntiProtonAntiProtonPion", // after this are mixed/same
      "ProtonPionSameAntiPionMixed","ProtonMixedPionAntiPionSame","ProtonAntiPionSamePionMixed", 
      "AntiProtonPionSameAntiPionMixed","AntiProtonMixedPionAntiPionSame","AntiProtonAntiPionSamePionMixed",
      "ProtonPionSamePionMixed", "ProtonMixedPionPionSame", 
      "AntiProtonPionSamePionMixed", "AntiProtonMixedPionPionSame", 
      "ProtonAntiPionSameAntiPionMixed", "ProtonMixedAntiPionAntiPionSame",
      "AntiProtonAntiPionSameAntiPionMixed", "AntiProtonMixedAntiPionAntiPionSame", 
      "ProtonProtonSamePionMixed", "ProtonMixedProtonPionSame",
      "AntiProtonAntiProtonSameAntiPionMixed", "AntiProtonMixedAntiProtonAntiPionSame",
      "ProtonProtonSameAntiPionMixed","ProtonMixedProtonAntiPionSame",
      "AntiProtonAntiProtonSamePionMixed", "AntiProtonMixedAntiProtonPionSame", "TRASH"};
    for(int i=0;i<33;i++){
      fSameEventTripletPhiThetaArray[i] = new TH2F(histTitlesSamePhiEta[i]+"Before",histTitlesSamePhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventTripletPhiThetaArray[33+i] = new TH2F(histTitlesSamePhiEta[i]+"After",histTitlesSamePhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[i]);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[33+i]);
    }

    // Mixed event phi theta distribution
    fMixedEventPhiTheta = new TList();
    fMixedEventPhiTheta->SetOwner();
    fMixedEventPhiTheta->SetName("MixedEventPhiTheta");

    fMixedEventTripletPhiThetaArray = new TH2F*[20];
    TString histTitlesMixedPhiEta[10] = {"ProtonPionAntiPion","AntiProtonPionAntiPion","ProtonPionPion",
      "AntiProtonPionPion", "ProtonAntiPionAntiPion", "AntiProtonAntiPionAntiPion",
      "ProtonProtonPion", "AntiProtonAntiProtonAntiPion",
      "ProtonProtonAntiPion", "AntiProtonAntiProtonPion"};
    for(int i=0;i<10;i++){
      fMixedEventTripletPhiThetaArray[i] = new TH2F(histTitlesMixedPhiEta[i]+"Before",histTitlesMixedPhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventTripletPhiThetaArray[10+i] = new TH2F(histTitlesMixedPhiEta[i]+"After",histTitlesMixedPhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[i]);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[10+i]);
    }

    fResultsThreeBody = new TList();
    fResultsThreeBody->SetOwner();
    fResultsThreeBody->SetName("ResultsThreeBody");
    fResultsThreeBody->Add(fSameEvent);
    fResultsThreeBody->Add(fMixedEvent);
    fResultsThreeBody->Add(fSameEventMult);
    fResultsThreeBody->Add(fMixedEventMult);
    fResultsThreeBody->Add(fSameEventPhiTheta);
    fResultsThreeBody->Add(fMixedEventPhiTheta);
    fResultsThreeBody->Add(fRejectedParticles);
    fResultsThreeBody->Add(fAcceptedParticles);
    fResultsThreeBody->Add(fOtherHistos);



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
  PostData(4, fPionList);
  PostData(5, fAntiPionList);
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

  if (fPion->GetIsMonteCarlo()) {
    if (!fPion->GetMinimalBooking()) {
      fPionMCList = fPion->GetMCQAHists();
    } else {
      fPionMCList = new TList();
      fPionMCList->SetName("MCTrackPionCuts");
      fPionMCList->SetOwner();
    }
    PostData(13, fPionMCList);
  }
  if (fAntiPion->GetIsMonteCarlo()) {
    if (!fAntiPion->GetMinimalBooking()) {
      fAntiPionMCList = fAntiPion->GetMCQAHists();
    } else {
      fAntiPionMCList = new TList();
      fAntiPionMCList->SetName("MCTrackAntiPionCuts");
      fAntiPionMCList->SetOwner();
    }
    PostData(14, fAntiPionMCList);
  }

 // Mixed event distribution ------------------------------------------------------------------------------
      // Take care of the mixing PartContainer for three particles
    auto ZVtxBinsSize = fConfig->GetNZVtxBins();
    auto MultBinsSize = fConfig->GetNMultBins();

    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    
    for(int iZVtx = 0; iZVtx<ZVtxBinsSize; iZVtx++){
      std::vector<std::vector<AliFemtoDreamPartContainer>> MultContainer;
      for(int iMult = 0; iMult<MultBinsSize; iMult++){
        std::vector<AliFemtoDreamPartContainer> AllUsedParticles;
        for(unsigned int iSpecies = 0; iSpecies<PDGCodes.size(); iSpecies++){
          auto tempPartContainer = new AliFemtoDreamPartContainer(fConfig->GetMixingDepth());
          AllUsedParticles.push_back(*tempPartContainer);
        }
        MultContainer.push_back(AllUsedParticles);
      }
      fPartContainer.push_back(MultContainer);
    }


}

void AliAnalysisTaskThreeBodyFemtoPionProton::UserExec(Option_t *option) {
  //AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);


  // 3-Body Trigger studies
  const int multiplicityForTrigger = fEvent->GetMultiplicity();
  static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
  fAllEvents->Fill(multiplicityForTrigger);

  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }
  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard ");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  std::vector<AliFemtoDreamBasePart> Pions;
  std::vector<AliFemtoDreamBasePart> AntiPions;

  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
    if (fPion->isSelected(fTrack)) {
      Pions.push_back(*fTrack);
    }
    if (fAntiPion->isSelected(fTrack)) {
      AntiPions.push_back(*fTrack);
    }
  }


  fPairCleaner->ResetArray();

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Pions);
  fPairCleaner->StoreParticle(AntiPions);

  
  if(fRunThreeBody){

    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    int bins[2] = { 0, 0 };
    float ZVtx = fEvent->GetZVertex();
    float Mult = fEvent->GetMultiplicity();
    fPartColl->FindBin(ZVtx, Mult, bins);

    // Same event distribution -------------------------------------------------------------------------------
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles();

    if(fTriggerOn){
      MyLovely3BodyTrigger( ParticleVector, fIsMC,PDGCodes, fAcceptedParticlesArray,fRejectedParticlesArray, *fConfig,  Mult, fSameEventTripletPhiThetaArray, 32);
    }


    // Proton Pion AntiPion 
    FillTripletDistribution( ParticleVector, 0, 2, 3, fSameEventTripletArray[0],PDGCodes, bins[1],fSameEventTripletMultArray[0], fSameEventTripletPhiThetaArray,0, *fConfig);
    // AntiProton Pion AntiPion 
    FillTripletDistribution( ParticleVector, 1, 2, 3, fSameEventTripletArray[1],PDGCodes, bins[1],fSameEventTripletMultArray[1], fSameEventTripletPhiThetaArray,1, *fConfig);
    // Proton Pion Pion
    FillTripletDistribution( ParticleVector, 0, 2, 2, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletPhiThetaArray,2, *fConfig);
    // AntiProton Pion Pion
    FillTripletDistribution( ParticleVector, 1, 2, 2, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletPhiThetaArray,3, *fConfig);
    // Proton AntiPion AntiPion
    FillTripletDistribution( ParticleVector, 0, 3, 3, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletPhiThetaArray,4, *fConfig);
    //AntiProton AntiPion AntiPion
    FillTripletDistribution( ParticleVector, 1, 3, 3, fSameEventTripletArray[5],PDGCodes, bins[1],fSameEventTripletMultArray[5], fSameEventTripletPhiThetaArray,5, *fConfig);
    //Proton Proton Pion
    FillTripletDistribution( ParticleVector, 0, 0, 2, fSameEventTripletArray[6],PDGCodes, bins[1],fSameEventTripletMultArray[6], fSameEventTripletPhiThetaArray,6, *fConfig);
    //AntiProton AntiProton AntiPion
    FillTripletDistribution( ParticleVector, 1, 1, 3, fSameEventTripletArray[7],PDGCodes, bins[1],fSameEventTripletMultArray[7], fSameEventTripletPhiThetaArray,7, *fConfig);
    // Proton Proton AntiPion
    FillTripletDistribution( ParticleVector, 0, 0, 3, fSameEventTripletArray[8],PDGCodes, bins[1],fSameEventTripletMultArray[8], fSameEventTripletPhiThetaArray,8, *fConfig);
    // AntiProton AntiProton Pion
    FillTripletDistribution( ParticleVector, 1, 1, 2, fSameEventTripletArray[9],PDGCodes, bins[1],fSameEventTripletMultArray[9], fSameEventTripletPhiThetaArray,9, *fConfig);


     // Mixed event distribution

    // TAKE CARE OF MULT AND ZVtx!!!!!!!!!!1
    if (!(bins[0] == -99 || bins[1] == -99)) {

      auto itZVtx = fPartContainer.begin()+ bins[0];
      auto itMult = itZVtx->begin() + bins[1];

      // Test (pp)L and (Lp)p
      // ProtonPionSame AntiPionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 3, fSameEventTripletArray[10], PDGCodes, bins[1],fSameEventTripletMultArray[10], fSameEventTripletPhiThetaArray, 10, *fConfig);
      // ProtonMixed PionAntiPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 3, 0, fSameEventTripletArray[11], PDGCodes, bins[1],fSameEventTripletMultArray[11], fSameEventTripletPhiThetaArray, 11, *fConfig);
      // ProtonAntiPionSame PionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 3, 2, fSameEventTripletArray[12], PDGCodes, bins[1],fSameEventTripletMultArray[12], fSameEventTripletPhiThetaArray, 12, *fConfig);
      // AntiProtonPionSame AntiPionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 2, 3, fSameEventTripletArray[13], PDGCodes, bins[1],fSameEventTripletMultArray[13], fSameEventTripletPhiThetaArray, 13, *fConfig);
      // AntiProtonMixed PionAntiPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 3, 1, fSameEventTripletArray[14], PDGCodes, bins[1],fSameEventTripletMultArray[14], fSameEventTripletPhiThetaArray, 14, *fConfig);
      // AntiProtonAntiPionSame PionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 2, fSameEventTripletArray[15], PDGCodes, bins[1],fSameEventTripletMultArray[15], fSameEventTripletPhiThetaArray, 15, *fConfig);
      // ProtonPionSame PionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 2, fSameEventTripletArray[16], PDGCodes, bins[1],fSameEventTripletMultArray[16], fSameEventTripletPhiThetaArray, 16, *fConfig);
      // ProtonMixedPionPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 0, fSameEventTripletArray[17], PDGCodes, bins[1],fSameEventTripletMultArray[17], fSameEventTripletPhiThetaArray, 17, *fConfig);
      // AntiProtonPionSame PionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 2, 2, fSameEventTripletArray[18], PDGCodes, bins[1],fSameEventTripletMultArray[18], fSameEventTripletPhiThetaArray, 18, *fConfig);
      // AntiProtonMixedPionPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 1, fSameEventTripletArray[19], PDGCodes, bins[1],fSameEventTripletMultArray[19], fSameEventTripletPhiThetaArray, 19, *fConfig);
      // ProtonAntiPionSame AntiPionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 3, 3, fSameEventTripletArray[20], PDGCodes, bins[1],fSameEventTripletMultArray[20], fSameEventTripletPhiThetaArray, 20, *fConfig);
      // ProtonMixedAntiPionAntiPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 0, fSameEventTripletArray[21], PDGCodes, bins[1],fSameEventTripletMultArray[21], fSameEventTripletPhiThetaArray, 21, *fConfig);
      // AntiProtonAntiPionSameAntiPionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 3, fSameEventTripletArray[22], PDGCodes, bins[1],fSameEventTripletMultArray[22], fSameEventTripletPhiThetaArray, 22, *fConfig);
      // AntiProtonMixedAntiPionAntiPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 1, fSameEventTripletArray[23], PDGCodes, bins[1],fSameEventTripletMultArray[23], fSameEventTripletPhiThetaArray, 23, *fConfig);
      // ProtonProtonSamePionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 2, fSameEventTripletArray[24], PDGCodes, bins[1],fSameEventTripletMultArray[24], fSameEventTripletPhiThetaArray, 24, *fConfig);
      // ProtonMixedProtonPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 0, fSameEventTripletArray[25], PDGCodes, bins[1],fSameEventTripletMultArray[25], fSameEventTripletPhiThetaArray, 25, *fConfig);
      // AntiProtonAntiProtonSame AntiPionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 3, fSameEventTripletArray[26], PDGCodes, bins[1],fSameEventTripletMultArray[26], fSameEventTripletPhiThetaArray, 26, *fConfig);
      // AntiProtonMixed AntiProtonAntiPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 1, fSameEventTripletArray[27], PDGCodes, bins[1],fSameEventTripletMultArray[27], fSameEventTripletPhiThetaArray, 27, *fConfig);
      // ProtonProtonSame AntiPionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 3, fSameEventTripletArray[28], PDGCodes, bins[1],fSameEventTripletMultArray[28], fSameEventTripletPhiThetaArray, 28, *fConfig);
      // ProtonMixed ProtonAntiPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 3, 0, fSameEventTripletArray[29], PDGCodes, bins[1],fSameEventTripletMultArray[29], fSameEventTripletPhiThetaArray, 29, *fConfig);
      // AntiProtonAntiProtonSame PionMixed
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 2, fSameEventTripletArray[30], PDGCodes, bins[1],fSameEventTripletMultArray[30], fSameEventTripletPhiThetaArray, 30, *fConfig);
      // AntiProtonMixedAntiProtonPionSame
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 2, 1, fSameEventTripletArray[31], PDGCodes, bins[1],fSameEventTripletMultArray[31], fSameEventTripletPhiThetaArray, 31, *fConfig);

      // Normal mixing
      FillTripletDistributionMEPP(ParticleVector, *itMult, 0, 2, 3, fMixedEventTripletArray[0], PDGCodes, bins[1],fMixedEventTripletMultArray[0], fMixedEventTripletPhiThetaArray,0, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 1, 2, 3, fMixedEventTripletArray[1], PDGCodes, bins[1],fMixedEventTripletMultArray[1], fMixedEventTripletPhiThetaArray,1, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 0, 2, 2, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletPhiThetaArray,2, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 1, 2, 2, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletPhiThetaArray,3, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 0, 3, 3, fMixedEventTripletArray[4], PDGCodes, bins[1],fMixedEventTripletMultArray[4], fMixedEventTripletPhiThetaArray,4, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 1, 3, 3, fMixedEventTripletArray[5], PDGCodes, bins[1],fMixedEventTripletMultArray[5], fMixedEventTripletPhiThetaArray,5, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 0, 0, 2, fMixedEventTripletArray[6], PDGCodes, bins[1],fMixedEventTripletMultArray[6], fMixedEventTripletPhiThetaArray,6, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 1, 1, 3, fMixedEventTripletArray[7], PDGCodes, bins[1],fMixedEventTripletMultArray[7], fMixedEventTripletPhiThetaArray,7, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 0, 0, 3, fMixedEventTripletArray[8], PDGCodes, bins[1],fMixedEventTripletMultArray[8], fMixedEventTripletPhiThetaArray,8, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 1, 1, 2, fMixedEventTripletArray[9], PDGCodes, bins[1],fMixedEventTripletMultArray[9], fMixedEventTripletPhiThetaArray,9, *fConfig);

      // Update the particle container with current event
      SetMixedEvent(ParticleVector, &(*itMult));
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
  PostData(4, fPionList);
  PostData(5, fAntiPionList);
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
  if (fPion->GetIsMonteCarlo()) {
    PostData(13, fPionMCList);
  }
  if (fAntiPion->GetIsMonteCarlo()) {
    PostData(14, fAntiPionMCList);
  }

}

//____________________________________________________________________________________________________
void AliAnalysisTaskThreeBodyFemtoPionProton::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskThreeBodyFemtoPionProton::StoreGlobalTrackReference(AliVTrack *track) {
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
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap()|| fGTI[trackID]->GetTPCNcls()) {
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


TLorentzVector AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(
    TLorentzVector &PartOne, TLorentzVector &PartTwo) {
  // The q12 components can be calculated as:
  //             q^mu = (p1-p2)^mu /2 - ((p1-p2)*P/(2P^2))*P^mu
  //             where P = p1+p2
  // Reference: https://www.annualreviews.org/doi/pdf/10.1146/annurev.nucl.55.090704.151533
  // In the following code the above written equation will be expressed as:
  //             q = trackDifference/2 -  scaling * trackSum
  // where scaling is a float number:
  //             scaling = trackDifference*trackSum/(2*trackSum^2) = ((p1-p2)*P/(2P^2))
  // PR on 15-06-2020 -> don't use the reduced vector - no division by 2
  TLorentzVector trackSum = PartOne + PartTwo;
  TLorentzVector trackDifference = PartOne - PartTwo;
  float scaling = trackDifference*trackSum/(trackSum*trackSum);
  TLorentzVector qPartOnePartTwo = trackDifference -  scaling * trackSum;

  return qPartOnePartTwo;
}


void AliAnalysisTaskThreeBodyFemtoPionProton::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // This function creates a triplet distribution in Q3 bins (defined lower).
  // It requires the particle vector from PairCleaner() and the three indices of particles of interest. So
  // if you want to get distribution for particles that are saved in particle vector as 1 2 3 element, just 
  // call the function with firstSpecies=1,secondSpecies=2,thirdSpecies=3

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;
  auto Particle3Vector = ParticleVector.begin()+thirdSpecies;

  // Get the PID codes std::vector<int> 
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;
  auto itPDGPar3 = PDGCodes.begin()+thirdSpecies;
  // Get particle masses 
  auto massparticle1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massparticle2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();
  auto massparticle3 = TDatabasePDG::Instance()->GetParticle(*itPDGPar3)->Mass();

  
  unsigned int DoThisPair12 = 11;
  unsigned int DoThisPair23 = 11;
  unsigned int DoThisPair31 = 11;


  // loop over first particle 
  for (auto iPart1 = Particle1Vector->begin(); iPart1 != Particle1Vector->end(); ++iPart1) {
    // if second particle species is different than first - start with the first particle in the vector
    auto iPart2 = Particle2Vector->begin();
    // if second particle  and first are the species, start second loop from the next particle (to not double count)
    if (firstSpecies==secondSpecies) iPart2 = iPart1+1;
    // loop over second particle ...
    for (; iPart2 != Particle2Vector->end(); ++iPart2) {
      auto iPart3 = Particle3Vector->begin();
      if (firstSpecies==thirdSpecies) iPart3 = iPart1+1;
      if (secondSpecies==thirdSpecies) iPart3 = iPart2+1;
      for ( ; iPart3 != Particle3Vector->end(); ++iPart3) {
        bool Pair12 = true;
        bool Pair23 = true;
        bool Pair31 = true;
        Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
        Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
        Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
        if(!Pair12||!Pair23||!Pair31) {continue;}
        // Now we have the three particles, lets create their Lorentz vectors  
        TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
        part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
        iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massparticle1,2)));
        part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
        iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massparticle2,2)));
        part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
        iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massparticle3,2)));
        
        // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
        TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part3_LorVec,part1_LorVec);
        // The particles in current methodology are put in bins of:
        //                 Q3=sqrt(q12^2+q23^2+q31^2)
        float Q32 = q12*q12+q23*q23+q31*q31;
        // From 3 pion paper, the q must be multiplied by -1 before taking quare root
        float Q3 = sqrt(-Q32); // the minus from pion paper
        hist->Fill(Q3); 
        hist2d->Fill(Q3,mult+1); 
        
      }
    }
  }
}



void AliAnalysisTaskThreeBodyFemtoPionProton::SetMixedEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  for(unsigned int iSpecies = 0; iSpecies<ParticleVector.size(); iSpecies++){
    if ((ParticleVector.begin()+iSpecies)->size() > 0) {
      (fPartContainer->begin()+iSpecies)->SetEvent(*(ParticleVector.begin()+iSpecies));
    }
  }
}


void AliAnalysisTaskThreeBodyFemtoPionProton::FillTripletDistributionMEPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // Description of function given in AliAnalysisTaskThreeBodyFemtoPionProton::FillTripletDistribution
  // In this function, only one particle is used from current event, and the other two - from other two events

  // Current behavior with the mixed events:
  //              1) implemented ONLY IF ME1 and ME2 are the same species!!!!!!!!!!! 
  //              2) the first of the two ME particles: takes Nth event from ME, takes every particle in this event
  //              2) the second of the two ME particles: takes (N+1)th event from ME, takes every particle in this event

  
  auto ParticleSE = ParticleVector.begin()+speciesSE;
  auto MixedEvent1Container = fPartContainer.begin()+speciesME1;
  auto MixedEvent2Container = fPartContainer.begin()+speciesME2;

  // Get the PID codes std::vector<int> 
  auto itPDGParSE = PDGCodes.begin()+speciesSE;
  auto itPDGParME1 = PDGCodes.begin()+speciesME1;
  auto itPDGParME2 = PDGCodes.begin()+speciesME2;

  // Get particle masses 
  auto massParticleSE = TDatabasePDG::Instance()->GetParticle(*itPDGParSE)->Mass();
  auto massParticleME1 = TDatabasePDG::Instance()->GetParticle(*itPDGParME1)->Mass();
  auto massParticleME2 = TDatabasePDG::Instance()->GetParticle(*itPDGParME2)->Mass();

  unsigned int DoThisPair12 = 11;
  unsigned int DoThisPair23 = 11;
  unsigned int DoThisPair31 = 11;

  
  // loop over first particle 
  for (auto iPart1 = ParticleSE->begin(); iPart1 != ParticleSE->end(); ++iPart1) {
    // loop over second particle ...
    for (int iDepth1 = 0; iDepth1 < (int) MixedEvent1Container->GetMixingDepth(); ++iDepth1) {    
      std::vector<AliFemtoDreamBasePart> iEvent2 = MixedEvent1Container->GetEvent(iDepth1);
      for ( auto iPart2 = iEvent2.begin(); iPart2 != iEvent2.end(); ++iPart2) {
        int iDepth2 = 0;
        if(speciesME1==speciesME2) iDepth2 = iDepth1+1; 
        for ( ; iDepth2 < (int) MixedEvent2Container->GetMixingDepth(); ++iDepth2) {
          std::vector<AliFemtoDreamBasePart> iEvent3 = MixedEvent2Container->GetEvent(iDepth2);
          for ( auto iPart3 = iEvent3.begin(); iPart3 != iEvent3.end(); ++iPart3) {
            bool Pair12 = true;
            bool Pair23 = true;
            bool Pair31 = true;
            Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config);
            Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,false,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config);
            Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,false,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config);
            if(!Pair12||!Pair23||!Pair31) {continue;}

            // Now we have the three particles, lets create their Lorentz vectors  
            TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
            part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
            iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE,2)));
            part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
            iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleME1,2)));
            part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
            iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME2,2)));
            // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
            TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part3_LorVec,part1_LorVec);
            // The particles in current methodology are put in bins of:
            //                 Q3=sqrt(q12^2+q23^2+q31^2)
            float Q32 = q12*q12+q23*q23+q31*q31;
            // From 3 pion paper, the q must be multiplied by -1 before taking quare root
            float Q3 = sqrt(-Q32); // the minus from pion paper
            hist->Fill(Q3); 
            hist2d->Fill(Q3,mult+1); 
          }
        }  
      }
    }
  }
}



void AliAnalysisTaskThreeBodyFemtoPionProton::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // Description of function given in AliAnalysisTaskThreeBodyFemtoPionProton::FillTripletDistribution
  // In this function, two particles are used from current event, and one - from other event


  auto ParticleSE1 = ParticleVector.begin()+speciesSE1;
  auto ParticleSE2 = ParticleVector.begin()+speciesSE2;
  auto MixedEventContainer = fPartContainer.begin()+speciesME;

  // Get the PID codes std::vector<int> 
  auto itPDGParSE1 = PDGCodes.begin()+speciesSE1;
  auto itPDGParSE2 = PDGCodes.begin()+speciesSE2;
  auto itPDGParME = PDGCodes.begin()+speciesME;

  // Get particle masses 
  auto massParticleSE1 = TDatabasePDG::Instance()->GetParticle(*itPDGParSE1)->Mass();
  auto massParticleSE2 = TDatabasePDG::Instance()->GetParticle(*itPDGParSE2)->Mass();
  auto massParticleME = TDatabasePDG::Instance()->GetParticle(*itPDGParME)->Mass();
  // get info for close pair rejection
  unsigned int DoThisPair12 = 11;
  unsigned int DoThisPair23 = 11;
  unsigned int DoThisPair31 = 11;


  // loop over first particle 
  for (auto iPart1 = ParticleSE1->begin(); iPart1 != ParticleSE1->end(); ++iPart1) {
    auto iPart2 = ParticleSE2->begin();
    if(speciesSE1==speciesSE2) iPart2 = iPart1+1;
    for(;iPart2 != ParticleSE2->end(); ++iPart2)
    {
      for ( int iDepth1 = 0; iDepth1 < (int) MixedEventContainer->GetMixingDepth(); ++iDepth1) {
        std::vector<AliFemtoDreamBasePart> iEvent3 = MixedEventContainer->GetEvent(iDepth1);
        for ( auto iPart3 = iEvent3.begin(); iPart3 != iEvent3.end(); ++iPart3) {
          // Now we have the three particles, lets create their Lorentz vectors  
          bool Pair12 = true;
          bool Pair23 = true;
          bool Pair31 = true;
          Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
          Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,false,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
          Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,false,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
          if(!Pair12||!Pair23||!Pair31) {continue;}


          TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
          part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
          iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE1,2)));
          part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
          iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleSE2,2)));
          part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
          iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME,2)));
          // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
          TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part1_LorVec,part2_LorVec);
          TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part2_LorVec,part3_LorVec);
          TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part3_LorVec,part1_LorVec);
          // The particles in current methodology are put in bins of:
          //                 Q3=sqrt(q12^2+q23^2+q31^2)
          float Q32 = q12*q12+q23*q23+q31*q31;
          // From 3 pion paper, the q must be multiplied by -1 before taking quare root
          float Q3 = sqrt(-Q32); // the minus from pion paper
          hist->Fill(Q3); 
          hist2d->Fill(Q3,mult+1); 
        }
      }  
    }
  }

}

bool AliAnalysisTaskThreeBodyFemtoPionProton::DeltaEtaDeltaPhi(
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist,
                                                   AliFemtoDreamCollConfig Config) {
  // DoThisPair = ij where i is the number of daughters for first particle, j for the second
  
  static const float piHi = TMath::Pi();
  auto fDeltaPhiSqMax = Config.GetDeltaPhiMax() * Config.GetDeltaPhiMax();
  auto fDeltaEtaSqMax = Config.GetDeltaEtaMax() * Config.GetDeltaEtaMax() ;

  bool pass = true;
  // if nDaug == 1 => Single Track, else decay
  unsigned int nDaug1 = (unsigned int) DoThisPair / 10;
  if (nDaug1 > 9) {
    AliWarning("you are doing something wrong \n");
  }
  if (nDaug1 > part1.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 1 (%u) and Radii 1 (%u) do not correspond \n",
            DoThisPair, nDaug1, (unsigned int)part1.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  unsigned int nDaug2 = (unsigned int) DoThisPair % 10;

  if (nDaug2 > part2.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 2 (%u) and Radii 2 (%u) do not correspond \n",
            DoThisPair, nDaug2, (unsigned int)part2.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  std::vector<float> eta1 = part1.GetEta();
  std::vector<float> eta2 = part2.GetEta();

  for (unsigned int iDaug1 = 0; iDaug1 < nDaug1; ++iDaug1) {
    std::vector<float> PhiAtRad1 = part1.GetPhiAtRaidius().at(iDaug1);
    float etaPar1;
    if (nDaug1 == 1) {
      etaPar1 = eta1.at(0);
    } else {
      etaPar1 = eta1.at(iDaug1 + 1);
    }
    for (unsigned int iDaug2 = 0; iDaug2 < nDaug2; ++iDaug2) {
      std::vector<float> phiAtRad2 = part2.GetPhiAtRaidius().at(iDaug2);
      float etaPar2;
      if (nDaug2 == 1) {
        etaPar2 = eta2.at(0);
      } else {
        etaPar2 = eta2.at(iDaug2 + 1);
      }
      float deta = etaPar1 - etaPar2;
      const int size =
          (PhiAtRad1.size() > phiAtRad2.size()) ?
              phiAtRad2.size() : PhiAtRad1.size();
      float dphiAvg = 0;
      for (int iRad = 0; iRad < size; ++iRad) {
        float dphi = PhiAtRad1.at(iRad) - phiAtRad2.at(iRad);
        if (dphi > piHi) {
          dphi += -piHi * 2;
        } else if (dphi < -piHi) {
          dphi += piHi * 2;
        }
        dphi = TVector2::Phi_mpi_pi(dphi);

        dphiAvg += dphi;
      }
      beforeHist->Fill(dphiAvg/ (float) size, deta);
      if (pass) {
        if ((dphiAvg / (float) size) * (dphiAvg / (float) size) / fDeltaPhiSqMax
            + deta * deta / fDeltaEtaSqMax < 1.) {
          pass = false;
        }
        else{
          afterHist->Fill(dphiAvg/ (float) size, deta);
        }
      }
    }
  }
  return pass;
}

void AliAnalysisTaskThreeBodyFemtoPionProton::MyLovely3BodyTrigger(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector ,  bool isMC, std::vector<int> PDGCodes,  TH1F** histAccepted, TH1F** histRejected, AliFemtoDreamCollConfig Config, float mult, TH2F** fEventTripletPhiThetaArray, int phiEtaHistNo){
  
  bool triplet1=false, triplet2=false, triplet3=false, triplet4=false, triplet5=false;
  
  if(CalculatePPLTriggerQ3Min(ParticleVector, 0, 2, 3, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ||  CalculatePPLTriggerQ3Min(ParticleVector, 1, 2, 3, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ){
    triplet1 = true;
    fAcceptedParticlesArray[0]->Fill(mult);
  }
  else{
    fRejectedParticlesArray[0]->Fill(mult);
  }

  if(CalculatePPLTriggerQ3Min(ParticleVector, 0, 2, 2, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ||  CalculatePPLTriggerQ3Min(ParticleVector, 1, 3, 3, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ){
    triplet2 = true;
    fAcceptedParticlesArray[1]->Fill(mult);
  } 
  else{
    fRejectedParticlesArray[1]->Fill(mult);
  } 

  if(CalculatePPLTriggerQ3Min(ParticleVector, 1, 2, 2, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ||  CalculatePPLTriggerQ3Min(ParticleVector, 0, 3, 3, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ){
    triplet3 = true;
    fAcceptedParticlesArray[2]->Fill(mult);
  }
  else{
    fRejectedParticlesArray[2]->Fill(mult);
  }
  
  if(CalculatePPLTriggerQ3Min(ParticleVector, 0, 0, 2, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ||  CalculatePPLTriggerQ3Min(ParticleVector, 1, 1, 3, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ){
    triplet4 = true;
    fAcceptedParticlesArray[3]->Fill(mult);
  }
  else{
    fRejectedParticlesArray[3]->Fill(mult);
  }

  if(CalculatePPLTriggerQ3Min(ParticleVector, 0, 0, 3, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ||  CalculatePPLTriggerQ3Min(ParticleVector, 1, 1, 2, PDGCodes, Config, fEventTripletPhiThetaArray, phiEtaHistNo)<=fQ3Limit ){
    triplet5 = true;
    fAcceptedParticlesArray[4]->Fill(mult);
  }
  else{
    fRejectedParticlesArray[4]->Fill(mult);
  }

  if(triplet1||triplet2||triplet3){
    fAcceptedParticlesArray[5]->Fill(mult);
  }
  else{
    fRejectedParticlesArray[5]->Fill(mult);
  }

  if(triplet4||triplet5){
    fAcceptedParticlesArray[6]->Fill(mult);
  }
  else{
    fRejectedParticlesArray[6]->Fill(mult);
  }
  
}

double AliAnalysisTaskThreeBodyFemtoPionProton::CalculatePPLTriggerQ3Min(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, std::vector<int> PDGCodes, AliFemtoDreamCollConfig Config, TH2F** fEventTripletPhiThetaArray, int phiEtaHistNo){
  // This function creates a triplet distribution in Q3 bins (defined lower).
  // It requires the particle vector from PairCleaner() and the three indices of particles of interest. So
  // if you want to get distribution for particles that are saved in particle vector as 1 2 3 element, just 
  // call the function with firstSpecies=1,secondSpecies=2,thirdSpecies=3

  double minQ3InEvent = 1000.; // min q3 tripelt present in the sample

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;
  auto Particle3Vector = ParticleVector.begin()+thirdSpecies;

  // Get the PID codes std::vector<int> 
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;
  auto itPDGPar3 = PDGCodes.begin()+thirdSpecies;
  // Get particle masses 
  auto massparticle1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massparticle2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();
  auto massparticle3 = TDatabasePDG::Instance()->GetParticle(*itPDGPar3)->Mass();

  unsigned int DoThisPair12 = 11;
  unsigned int DoThisPair23 = 11;
  unsigned int DoThisPair31 = 11;

  // loop over first particle 
  for (auto iPart1 = Particle1Vector->begin(); iPart1 != Particle1Vector->end(); ++iPart1) {
    // if second particle species is different than first - start with the first particle in the vector
    auto iPart2 = Particle2Vector->begin();
    // if second particle  and first are the species, start second loop from the next particle (to not double count)
    if (firstSpecies==secondSpecies) iPart2 = iPart1+1;
    // loop over second particle ...
    for (; iPart2 != Particle2Vector->end(); ++iPart2) {
      auto iPart3 = Particle3Vector->begin();
      if (firstSpecies==thirdSpecies) iPart3 = iPart1+1;
      if (secondSpecies==thirdSpecies) iPart3 = iPart2+1;
      for ( ; iPart3 != Particle3Vector->end(); ++iPart3) {

        bool Pair12 = true;
        bool Pair23 = true;
        bool Pair31 = true;
        Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
        Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,false,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
        Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,false,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[33+phiEtaHistNo],Config);
        if(!Pair12||!Pair23||!Pair31) {continue;}


        // Now we have the three particles, lets create their Lorentz vectors  
        TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
        part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
        iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massparticle1,2)));
        part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
        iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massparticle2,2)));
        part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
        iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massparticle3,2)));
        
        // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
        TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoPionProton::RelativePairMomentum(part3_LorVec,part1_LorVec);
        // The particles in current methodology are put in bins of:
        //                 Q3=sqrt(q12^2+q23^2+q31^2)
        float Q32 = q12*q12+q23*q23+q31*q31;
        // From 3 pion paper, the q must be multiplied by -1 before taking quare root
        float Q3 = sqrt(-Q32); // the minus from pion paper  
        if (Q3<minQ3InEvent)    minQ3InEvent= Q3;   
      }
    }
  }
  return minQ3InEvent;
}
