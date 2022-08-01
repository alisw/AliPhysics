/*
 * AliAnalysisTaskThreeBodyFemtoAOD.cxx
 *
 *  Created on: May 13, 2019
 *      Author: Laura Serksnyte 
 */
#include "AliAnalysisTaskThreeBodyFemtoAOD.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskThreeBodyFemtoAOD)
AliAnalysisTaskThreeBodyFemtoAOD::AliAnalysisTaskThreeBodyFemtoAOD()
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
      fSameEvent(nullptr),    
      fMixedEvent(nullptr),   
      fSameEventMult(nullptr),  
      fMixedEventMult(nullptr),   
      fSameEventPhiTheta(nullptr),  
      fMixedEventPhiTheta(nullptr),       
      fOtherHistos(nullptr),   
      fRandomGen(nullptr),   
      fRunThreeBody(true),
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fTriggerOn(false),
      fTriggerOnSample(false),
      fIsMC(false),
      fRunppp0ppL1(false),
      fQ3Limit(0.0),
      fQ3LimitSample(0.0),
      fQ3LimitSample2(0.0),
      fQ3LimitSampleFraction(0.0),
      fQ3LimitSampleFraction2(0.0),
      fQ3LimitFraction(0.0),
      fRejectedEvents(nullptr),
      fAcceptedEvents(nullptr),
      fAcceptedEventsButNoPPL(nullptr),
      fTripletsPerCollision(nullptr),
      fWhichSample(nullptr),
      fSameEventOnlyLowestQ3(nullptr),
      fSameEventOnlyLowestQ3Anti(nullptr),
      ftotalTripletCount(0),
      fEventCutsTrigger(nullptr),
      fEventCutsTriggerList(nullptr),
      fTrackCutsTrigger(nullptr),
      fTrackCutsTriggerList(nullptr),
      fAntiTrackCutTrigger(nullptr),
      fAntiTrackCutTriggerList(nullptr),
      fv0CutsTrigger(nullptr),
      fv0CutsTriggerList(nullptr),
      fAntiv0CutsTrigger(nullptr),
      fAntiv0CutsTriggerList(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskThreeBodyFemtoAOD::AliAnalysisTaskThreeBodyFemtoAOD(const char* name, bool isMC, bool triggerOn)
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
      fSameEvent(nullptr),    
      fMixedEvent(nullptr),   
      fSameEventMult(nullptr),  
      fMixedEventMult(nullptr),   
      fSameEventPhiTheta(nullptr),  
      fMixedEventPhiTheta(nullptr),       
      fOtherHistos(nullptr),   
      fRandomGen(nullptr),   
      fRunThreeBody(true),  
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fTriggerOn(false),
      fTriggerOnSample(false),
      fIsMC(false),
      fRunppp0ppL1(false),
      fQ3Limit(0.0),
      fQ3LimitSample(0.0),
      fQ3LimitSample2(0.0),
      fQ3LimitSampleFraction(0.0),
      fQ3LimitSampleFraction2(0.0),
      fQ3LimitFraction(0.0),
      fRejectedEvents(nullptr),
      fAcceptedEvents(nullptr),
      fAcceptedEventsButNoPPL(nullptr),
      fTripletsPerCollision(nullptr),
      fWhichSample(nullptr),
      fSameEventOnlyLowestQ3(nullptr),
      fSameEventOnlyLowestQ3Anti(nullptr),
      ftotalTripletCount(0),
      fEventCutsTrigger(nullptr),
      fEventCutsTriggerList(nullptr),
      fTrackCutsTrigger(nullptr),
      fTrackCutsTriggerList(nullptr),
      fAntiTrackCutTrigger(nullptr),
      fAntiTrackCutTriggerList(nullptr),
      fv0CutsTrigger(nullptr),
      fv0CutsTriggerList(nullptr),
      fAntiv0CutsTrigger(nullptr),
      fAntiv0CutsTriggerList(nullptr),
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
        if (triggerOn) {
          DefineOutput(11, TList::Class());  //Output for the event cuts trigger
          DefineOutput(12, TList::Class());  //Output for the p cuts trigger
          DefineOutput(13, TList::Class());  //Output for the pbar cuts trigger
          DefineOutput(14, TList::Class());  //Output for the L cuts trigger
          DefineOutput(15, TList::Class());  //Output for the Lbar cuts trigger
        }
        if (isMC) {
          DefineOutput(16, TList::Class());  //Output for the Track MC
          DefineOutput(17, TList::Class());  //Output for the Anti Track MC
          DefineOutput(18, TList::Class());  //Output for the V0 MC
          DefineOutput(19, TList::Class());  //Output for the Anti V0 MC
        }
      }

AliAnalysisTaskThreeBodyFemtoAOD::~AliAnalysisTaskThreeBodyFemtoAOD() {
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

  if(fEventCutsTrigger){
    delete fEventCutsTrigger;
  }
  if(fTrackCutsTrigger){
    delete fTrackCutsTrigger;
  }
  if(fAntiTrackCutTrigger){
    delete fAntiTrackCutTrigger;
  }
  if(fv0CutsTrigger){
    delete fv0CutsTrigger;
  }
  if(fAntiv0CutsTrigger){
    delete fAntiv0CutsTrigger;
  }
}

void AliAnalysisTaskThreeBodyFemtoAOD::UserCreateOutputObjects() {
  fGTI = new AliAODTrack *[fTrackBufferSize];

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

  if(fTriggerOn){
    if (!fEventCutsTrigger) {
      AliError("No Event cuts \n");
    } else {
      fEventCutsTrigger->InitQA();
    }
    if (!fTrackCutsTrigger) {
      AliError("No Proton cuts \n");
    } else {
      fTrackCutsTrigger->Init();
    }
    if (!fAntiTrackCutTrigger) {
      AliError("No AntiProton cuts \n");
    } else {
      fAntiTrackCutTrigger->Init();
    }
    if (!fv0CutsTrigger) {
      AliError("No Lambda cuts \n");
    } else {
      fv0CutsTrigger->Init();
    }
    if (!fAntiv0CutsTrigger) {
      AliError("No AntiXi cuts \n");
    } else {
      fAntiv0CutsTrigger->Init();
    }
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

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }
  if(fTriggerOn){
    if (!fEventCutsTrigger->GetMinimalBooking()) {
      fEventCutsTriggerList = fEventCutsTrigger->GetHistList();
    } else {
      fEventCutsTriggerList = new TList();
      fEventCutsTriggerList->SetName("EventCutsTrigger");
      fEventCutsTriggerList->SetOwner();
    }
  }

  fProtonList = fProton->GetQAHists();
  fAntiProtonList = fAntiProton->GetQAHists();
  fLambdaList = fLambda->GetQAHists();
  fAntiLambdaList = fAntiLambda->GetQAHists();
  if(fTriggerOn){
    fTrackCutsTriggerList = fTrackCutsTrigger->GetQAHists();
    fAntiTrackCutTriggerList = fAntiTrackCutTrigger->GetQAHists();
    fv0CutsTriggerList = fv0CutsTrigger->GetQAHists();
    fAntiv0CutsTriggerList = fAntiv0CutsTrigger->GetQAHists();
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

  if (fRunThreeBody){
    fOtherHistos = new TList();  
    fOtherHistos->SetOwner(); 
    fOtherHistos->SetName("OtherHistos");

    fRejectedEvents =  new TH1F("fRejectedEvents","fRejectedEvents", 25, 0, 100);
    fAcceptedEvents =  new TH1F("fAcceptedEvents","fAcceptedEvents", 25, 0, 100);
    fAcceptedEventsButNoPPL =  new TH1F("fAcceptedEventsButNoPPL","fAcceptedEventsButNoPPL", 25, 0, 100);
    fTripletsPerCollision =  new TH1F("fTripletsPerCollision","fTripletsPerCollision", 25, 0, 25);
    fWhichSample = new TH1F("fWhichSample","fWhichSample", 4, 0, 4);
    fSameEventOnlyLowestQ3 = new TH1F("fSameEventOnlyLowestQ3","fSameEventOnlyLowestQ3", 8000, 0, 8);
    fSameEventOnlyLowestQ3Anti = new TH1F("fSameEventOnlyLowestQ3Anti","fSameEventOnlyLowestQ3Anti", 8000, 0, 8);

    fOtherHistos->Add(fRejectedEvents);
    fOtherHistos->Add(fAcceptedEvents);
    fOtherHistos->Add(fAcceptedEventsButNoPPL);
    fOtherHistos->Add(fTripletsPerCollision);
    fOtherHistos->Add(fWhichSample);
    fOtherHistos->Add(fSameEventOnlyLowestQ3);
    fOtherHistos->Add(fSameEventOnlyLowestQ3Anti);


    // Same event

    fSameEvent = new TList(); 
    fSameEvent->SetOwner(); 
    fSameEvent->SetName("SameEvent");

    fSameEventTripletArray = new TH1F*[22];
    TString histTitlesSame[22] ={"sameEventDistributionPL","sameEventDistributionPPL","sameEventDistributionAPAPAL",
      "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPLL","sameEventDistributionAPALAL",
      "sameEventDistributionLLL","sameEventDistributionALALAL", "sameEventDistributionppSameLMixed", "sameEventDistributionpLSamepMixed",
      "sameEventDistributionppSamepMixed","sameEventDistributionPP","sameEventDistributionapapSameaLMixed", "sameEventDistributionapaLSameapMixed",
      "sameEventDistributionapapSameapMixed", "sameEventDistributionpLSameLMixed", "sameEventDistributionLLSamepMixed",
      "sameEventDistributionapaLSameaLMixed", "sameEventDistributionaLaLSameapMixed","sameEventDistributionLLSameLMixed", "sameEventDistributionaLaLSameaLMixed"};
    fSameEventTripletArray[0] = new TH1F(histTitlesSame[0],histTitlesSame[0],1000,0, 1);
    fSameEvent->Add(fSameEventTripletArray[0]);
    for (int i = 1; i < 12; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }
    fSameEventTripletArray[12] = new TH1F(histTitlesSame[12],histTitlesSame[12],1000,0, 1);
    fSameEvent->Add(fSameEventTripletArray[12]);
    for (int i = 13; i <22 ; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }

    // Same event multiplicity dist
    fSameEventMult = new TList(); 
    fSameEventMult->SetOwner(); 
    fSameEventMult->SetName("SameEventMult");

    fSameEventTripletMultArray = new TH2F*[21];
    TString histTitlesSameMult[21] ={"sameEventDistributionMultPL","sameEventDistributionMultPPL","sameEventDistributionMultAPAPAL",
      "sameEventDistributionMultPPP", "sameEventDistributionMultAPAPAP", "sameEventDistributionMultPLL","sameEventDistributionMultAPALAL",
      "sameEventDistributionMultLLL","sameEventDistributionMultALALAL",  "sameEventDistributionMultppSameLMixed", "sameEventDistributionMultpLSamepMixed",
      "sameEventDistributionMultppSamepMixed","sameEventDistributionMultapapSameaLMixed", "sameEventDistributionMultapaLSameapMixed",
      "sameEventDistributionMultapapSameapMixed", "sameEventDistributionMultpLSameLMixed", "sameEventDistributionMultLLSamepMixed",
      "sameEventDistributionMultapaLSameaLMixed", "sameEventDistributionMultaLaLSameapMixed","sameEventDistributionMultLLSameLMixed", "sameEventDistributionMultaLaLSameaLMixed"};
    
    fSameEventTripletMultArray[0] = new TH2F(histTitlesSameMult[0],histTitlesSameMult[0],1000,0, 1,26,1,27);
    fSameEventMult->Add(fSameEventTripletMultArray[0]);
    for (int i = 1; i < 21; ++i) {
      fSameEventTripletMultArray[i] =  new TH2F(histTitlesSameMult[i],histTitlesSameMult[i], 8000, 0, 8,26,1,27);
      fSameEventMult->Add(fSameEventTripletMultArray[i]);
     }

    // Mixed event
    fMixedEvent = new TList(); 
    fMixedEvent->SetOwner();  
    fMixedEvent->SetName("MixedEvent");
     
    fMixedEventTripletArray = new TH1F*[14];
    TString histTitlesMixed[14] ={"mixedEventDistributionPL","mixedEventDistributionPPL","mixedEventDistributionAPAPAL",
      "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPLL","mixedEventDistributionAPALAL",
      "mixedEventDistributionLLL","mixedEventDistributionALALAL", "mixedEventDistributionPPLTEST", "mixedEventDistributionAPAPALTEST",
      "mixedEventDistributionPPLTESTppL", "mixedEventDistributionAPAPALTESTppL", "mixedEventDistributionPP"};
    fMixedEventTripletArray[0] = new TH1F(histTitlesMixed[0],histTitlesMixed[0],1000,0, 1);
    fMixedEvent->Add(fMixedEventTripletArray[0]);
    for (int i = 1; i < 13; ++i) {
      fMixedEventTripletArray[i] = new TH1F(histTitlesMixed[i],histTitlesMixed[i], 8000, 0, 8);
      fMixedEvent->Add(fMixedEventTripletArray[i]);
     }
    fMixedEventTripletArray[13] = new TH1F(histTitlesMixed[13],histTitlesMixed[13],1000,0, 1);
    fMixedEvent->Add(fMixedEventTripletArray[13]);

    // Mixed event multiplicity dist
    fMixedEventMult = new TList();  
    fMixedEventMult->SetOwner();  
    fMixedEventMult->SetName("MixedEventMult");

    fMixedEventTripletMultArray = new TH2F*[13];
    TString histTitlesMixedMult[13] ={"mixedEventDistributionPLMult","mixedEventDistributionPPLMult","mixedEventDistributionMultAPAPAL",
      "mixedEventDistributionMultPPP", "mixedEventDistributionMultAPAPAP", "mixedEventDistributionMultPLL","mixedEventDistributionMultAPALAL",
      "mixedEventDistributionMultLLL","mixedEventDistributionMultALALAL", "mixedEventDistributionMultPPLTEST", "mixedEventDistributionMultAPAPALTEST",
      "mixedEventDistributionMultPPLTESTppL", "mixedEventDistributionMultAPAPALTESTppL"};
    fMixedEventTripletMultArray[0] = new TH2F(histTitlesMixedMult[0],histTitlesMixedMult[0],1000,0, 1 ,26,1,27);
    fMixedEventMult->Add(fMixedEventTripletMultArray[0]);
    for (int i = 1; i < 13; ++i) {
      fMixedEventTripletMultArray[i] = new TH2F(histTitlesMixedMult[i],histTitlesMixedMult[i], 8000, 0, 8,26,1,27);
      fMixedEventMult->Add(fMixedEventTripletMultArray[i]);
     }

    // Same event phi theta
    fSameEventPhiTheta = new TList();  
    fSameEventPhiTheta->SetOwner(); 
    fSameEventPhiTheta->SetName("SameEventPhiTheta");

    fSameEventTripletPhiThetaArray = new TH2F*[42];
    TString histTitlesSamePhiEta[20] = {"sameEventPhiEtaPPL","sameEventPhiEtaAPAPAL",
      "sameEventPhiEtaPPP", "sameEventPhiEtaAPAPAP", "sameEventPhiEtaPLL","sameEventPhiEtaAPALAL",
      "sameEventPhiEtaLLL","sameEventPhiEtaALALAL", "sameEventDistributionPhiEtappSameLMixed", "sameEventDistributionPhiEtapLSamepMixed",
      "sameEventDistributionPhiEtappSamepMixed","sameEventDistributionPhiEtaapapSameaLMixed", "sameEventDistributionPhiEtaapaLSameapMixed",
      "sameEventDistributionPhiEtaapapSameapMixed", "sameEventDistributionPhiEtapLSameLMixed", "sameEventDistributionPhiEtaLLSamepMixed",
      "sameEventDistributionPhiEtaapaLSameaLMixed", "sameEventDistributionPhiEtaaLaLSameapMixed","sameEventDistributionPhiEtaLLSameLMixed", "sameEventDistributionPhiEtaaLaLSameaLMixed"};
    for(int i=0;i<20;i++){
      fSameEventTripletPhiThetaArray[i] = new TH2F(histTitlesSamePhiEta[i]+"Before",histTitlesSamePhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventTripletPhiThetaArray[20+i] = new TH2F(histTitlesSamePhiEta[i]+"After",histTitlesSamePhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[i]);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[20+i]);
    }
    fSameEventTripletPhiThetaArray[40] = new TH2F("sameEventPhiEtaPPBefore","sameEventPhiEtaPPBefore", 500, -0.15,0.15,500,-0.15,0.15);
    fSameEventTripletPhiThetaArray[41] = new TH2F("sameEventPhiEtaPPAfter","sameEventPhiEtaPPAfter", 500, -0.15,0.15,500,-0.15,0.15);
    fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[40]);
    fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[41]);

    // Mixed event phi theta

    fMixedEventPhiTheta = new TList();  
    fMixedEventPhiTheta->SetOwner();  
    fMixedEventPhiTheta->SetName("MixedEventPhiTheta");

    fMixedEventTripletPhiThetaArray = new TH2F*[42];
    TString histTitlesMixedPhiEta[20] = {"mixedEventPhiEtaPPL","mixedEventPhiEtaAPAPAL",
      "mixedEventPhiEtaPPP", "mixedEventPhiEtaAPAPAP", "mixedEventPhiEtaPLL","mixedEventPhiEtaAPALAL",
      "mixedEventPhiEtaLLL","mixedEventPhiEtaALALAL", "mixedEventPhiEtaPPLTEST2", "mixedEventPhiEtaAPAPALTEST2",
      "mixedEventPhiEtaPPLTEST3", "mixedEventPhiEtaAPAPALTEST3", "Empty11","Empty12",
      "Empty13","Empty14",
      "Empty1", "Empty2","Empty3", "Empty4"};
    for(int i=0;i<20;i++){
      fMixedEventTripletPhiThetaArray[i] = new TH2F(histTitlesMixedPhiEta[i]+"Before",histTitlesMixedPhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventTripletPhiThetaArray[20+i] = new TH2F(histTitlesMixedPhiEta[i]+"After",histTitlesMixedPhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[i]);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[20+i]);
    }
    fMixedEventTripletPhiThetaArray[40] = new TH2F("mixedEventPhiEtaPPBefore","mixedEventPhiEtaPPBefore", 500, -0.15,0.15,500,-0.15,0.15);
    fMixedEventTripletPhiThetaArray[41] = new TH2F("mixedEventPhiEtaPPAfter","mixedEventPhiEtaPPAfter", 500, -0.15,0.15,500,-0.15,0.15);
    fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[40]);
    fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[41]);

    fResultsThreeBody = new TList();  
    fResultsThreeBody->SetOwner();  
    fResultsThreeBody->SetName("ResultsThreeBody"); 
    fResultsThreeBody->Add(fSameEvent); 
    fResultsThreeBody->Add(fMixedEvent);  
    fResultsThreeBody->Add(fSameEventMult); 
    fResultsThreeBody->Add(fMixedEventMult);  
    fResultsThreeBody->Add(fSameEventPhiTheta); 
    fResultsThreeBody->Add(fMixedEventPhiTheta);  
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
  PostData(4, fLambdaList);
  PostData(5, fAntiLambdaList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fResultsSample);
  PostData(9, fResultsSampleQA);
  PostData(10, fResultsThreeBody);
  if(fTriggerOn){
    PostData(11, fEventCutsTriggerList);
    PostData(12, fTrackCutsTriggerList);
    PostData(13, fAntiTrackCutTriggerList);
    PostData(14, fv0CutsTriggerList);
    PostData(15, fAntiv0CutsTriggerList);
  }

  if (fProton->GetIsMonteCarlo()) {
    if (!fProton->GetMinimalBooking()) {
      fProtonMCList = fProton->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCTrkCuts");
      fProtonMCList->SetOwner();
    }
    PostData(16, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    if (!fAntiProton->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProton->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiTrkCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(17, fAntiProtonMCList);
  }

  if (fLambda->GetIsMonteCarlo()) {
    if (!fLambda->GetMinimalBooking()) {
      fLambdaMCList = fLambda->GetMCQAHists();
    } else {
      fLambdaMCList = new TList();
      fLambdaMCList->SetName("MCv0Cuts");
      fLambdaMCList->SetOwner();
    }
    PostData(18, fLambdaMCList);
  }
  if (fAntiLambda->GetIsMonteCarlo()) {
    if (!fAntiLambda->GetMinimalBooking()) {
      fAntiLambdaMCList = fAntiLambda->GetMCQAHists();
    } else {
      fAntiLambdaMCList = new TList();
      fAntiLambdaMCList->SetName("MCAntiv0Cuts");
      fAntiLambdaMCList->SetOwner();
    }
    PostData(19, fAntiLambdaMCList);
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
      fPartContainerTEST.push_back(MultContainer);
    }

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
      fPartContainerTESTppL.push_back(MultContainer);
    }

    fRandomGen = new TRandom3();
    fRandomGen->SetSeed(69);


}

void AliAnalysisTaskThreeBodyFemtoAOD::UserExec(Option_t *option) {
  AliAODEvent *evt = static_cast<AliAODEvent*>(InputEvent());
  if (!evt) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(evt);


  // 3-Body Trigger studies
  const int multiplicityForTrigger = fEvent->GetMultiplicity();
  static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
  if(fTriggerOn){
    if(!MyLovely3BodyTrigger(evt , fIsMC,PDGCodes)){
      fRejectedEvents->Fill(multiplicityForTrigger);
      return;
    }
    fAcceptedEvents->Fill(multiplicityForTrigger);
  }

  

  if (!fEventCuts->isSelected(fEvent)) {
    fAcceptedEventsButNoPPL->Fill(99);
    return;
  }

  
  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    fTrack->SetTrack(track);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
  }

  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;
  //AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);
  fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
      iv0 < static_cast<TClonesArray *>(evt->GetV0s())->GetEntriesFast();
      ++iv0) {
    AliAODv0* v0 = evt->GetV0(iv0);
    fv0->Setv0(evt, v0);
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
    // for trigger, check contamination
    TAxis *axis = fSameEventTripletArray[1]->GetXaxis(); 
    int bmin = axis->FindBin(0.); //in your case xmin=-1.5
    int bmax = axis->FindBin(fQ3LimitSample2); //in your case xmax=0.8


    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    int bins[2] = { 0, 0 };
    float ZVtx = fEvent->GetZVertex();
    float Mult = fEvent->GetMultiplicity();
    fPartColl->FindBin(ZVtx, Mult, bins);

    // Same event distribution -------------------------------------------------------------------------------
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles();
    // proton lambda, as a test case
    //FillPairDistributionPL(ParticleVector,fSameEventTripletArray[0],bins[1],fSameEventTripletMultArray[0]);
    // pp
    //FillPairDistributionPP(ParticleVector,fSameEventTripletArray[12],  fSameEventTripletPhiThetaArray, *fConfig);

    // proton proton lambda
    int entriesPPLBeforeFill = fSameEventTripletArray[1]->Integral(bmin,bmax);
    if(fRunppp0ppL1==true) FillTripletDistribution( ParticleVector, 0, 2, 0, fSameEventTripletArray[1],PDGCodes, bins[1],fSameEventTripletMultArray[1], fSameEventTripletPhiThetaArray,0, *fConfig);
    if(fRunppp0ppL1==false) FillTripletDistribution( ParticleVector, 0, 0, 0, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletPhiThetaArray,2, *fConfig);
    int entriesPPLAfterFill = fSameEventTripletArray[1]->Integral(bmin,bmax);


    int newEntriesPPL = entriesPPLAfterFill-entriesPPLBeforeFill;

    // antiproton antiproton antilambda
    int entriesAPAPALBeforeFill = fSameEventTripletArray[2]->Integral(bmin,bmax);
    if(fRunppp0ppL1==true) FillTripletDistribution( ParticleVector, 1, 3, 1, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletPhiThetaArray,1, *fConfig);
    if(fRunppp0ppL1==false) FillTripletDistribution( ParticleVector, 1, 1, 1, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletPhiThetaArray,3, *fConfig);
    
    int entriesAPAPALAfterFill = fSameEventTripletArray[2]->Integral(bmin,bmax);


    int newEntriesAPAPAL = entriesAPAPALAfterFill-entriesAPAPALBeforeFill;

    if((newEntriesPPL+newEntriesAPAPAL)==0){
      fAcceptedEventsButNoPPL->Fill(Mult);
    }




    /*// proton proton proton 
    FillTripletDistribution( ParticleVector, 0, 0, 0, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletPhiThetaArray,2, *fConfig);
    // antiproton antiproton antiproton 
    FillTripletDistribution( ParticleVector, 1, 1, 1, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletPhiThetaArray,3, *fConfig);
    // proton lambda lambda 
    FillTripletDistribution( ParticleVector, 0, 2, 2, fSameEventTripletArray[5],PDGCodes, bins[1],fSameEventTripletMultArray[5], fSameEventTripletPhiThetaArray,4, *fConfig);
    // antiproton antilambda antilambda 
    FillTripletDistribution( ParticleVector, 1, 3, 3, fSameEventTripletArray[6],PDGCodes, bins[1],fSameEventTripletMultArray[6], fSameEventTripletPhiThetaArray,5, *fConfig);
    // lambda lambda lambda 
    FillTripletDistribution( ParticleVector, 2, 2, 2, fSameEventTripletArray[7],PDGCodes, bins[1],fSameEventTripletMultArray[7], fSameEventTripletPhiThetaArray,6, *fConfig);
    // antilambda antilambda antilambda 
    FillTripletDistribution( ParticleVector, 3, 3, 3, fSameEventTripletArray[8],PDGCodes, bins[1],fSameEventTripletMultArray[8], fSameEventTripletPhiThetaArray,7, *fConfig);
*/

     // Mixed event distribution

    // TAKE CARE OF MULT AND ZVtx!!!!!!!!!!1
    if (!(bins[0] == -99 || bins[1] == -99)) {

      auto itZVtx = fPartContainer.begin()+ bins[0];
      auto itMult = itZVtx->begin() + bins[1];

      auto itZVtxTEST = fPartContainerTEST.begin()+ bins[0];
      auto itMultTEST = itZVtxTEST->begin() + bins[1];

      auto itZVtxTESTppL = fPartContainerTESTppL.begin()+ bins[0];
      auto itMultTESTppL = itZVtxTESTppL->begin() + bins[1];

      // Test (pp)L and (Lp)p

      /*FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 2, fSameEventTripletArray[9], PDGCodes, bins[1],fSameEventTripletMultArray[9], fSameEventTripletPhiThetaArray, 8, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 0, fSameEventTripletArray[10], PDGCodes, bins[1],fSameEventTripletMultArray[10], fSameEventTripletPhiThetaArray, 9, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 0, fSameEventTripletArray[11], PDGCodes, bins[1],fSameEventTripletMultArray[11], fSameEventTripletPhiThetaArray, 10, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 3, fSameEventTripletArray[13], PDGCodes, bins[1],fSameEventTripletMultArray[12], fSameEventTripletPhiThetaArray, 11, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 1, fSameEventTripletArray[14], PDGCodes, bins[1],fSameEventTripletMultArray[13], fSameEventTripletPhiThetaArray, 12, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 1, fSameEventTripletArray[15], PDGCodes, bins[1],fSameEventTripletMultArray[14], fSameEventTripletPhiThetaArray, 13, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 2, fSameEventTripletArray[16], PDGCodes, bins[1],fSameEventTripletMultArray[15], fSameEventTripletPhiThetaArray, 14, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 0, fSameEventTripletArray[17], PDGCodes, bins[1],fSameEventTripletMultArray[16], fSameEventTripletPhiThetaArray, 15, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 3, fSameEventTripletArray[18], PDGCodes, bins[1],fSameEventTripletMultArray[17], fSameEventTripletPhiThetaArray, 16, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 1, fSameEventTripletArray[19], PDGCodes, bins[1],fSameEventTripletMultArray[18], fSameEventTripletPhiThetaArray, 17, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 2, fSameEventTripletArray[20], PDGCodes, bins[1],fSameEventTripletMultArray[19], fSameEventTripletPhiThetaArray, 18, *fConfig);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 3, fSameEventTripletArray[21], PDGCodes, bins[1],fSameEventTripletMultArray[20], fSameEventTripletPhiThetaArray, 19, *fConfig);

      //Try to reproduce the p-lambda result from FemtoDream
      FillPairDistributionME(ParticleVector, *itMult, 0, 2, fMixedEventTripletArray[0],PDGCodes,bins[1],fMixedEventTripletMultArray[0],  fMixedEventTripletPhiThetaArray, *fConfig);
      FillPairDistributionME(ParticleVector, *itMult, 0, 0, fMixedEventTripletArray[13],PDGCodes,bins[1],fMixedEventTripletMultArray[0],  fMixedEventTripletPhiThetaArray, *fConfig);
    
    
      // Normal mixing
      FillTripletDistributionMEPP(ParticleVector, *itMult, 2, 0, 0, fMixedEventTripletArray[1], PDGCodes, bins[1],fMixedEventTripletMultArray[1], fMixedEventTripletPhiThetaArray,0, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 3, 1, 1, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletPhiThetaArray,1, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 0, 0, 0, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletPhiThetaArray,2, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 1, 1, 1, fMixedEventTripletArray[4], PDGCodes, bins[1],fMixedEventTripletMultArray[4], fMixedEventTripletPhiThetaArray,3, *fConfig);

      FillTripletDistributionMEPP(ParticleVector, *itMult, 0, 2, 2, fMixedEventTripletArray[5], PDGCodes, bins[1],fMixedEventTripletMultArray[5], fMixedEventTripletPhiThetaArray,4, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 1, 3, 3, fMixedEventTripletArray[6], PDGCodes, bins[1],fMixedEventTripletMultArray[6], fMixedEventTripletPhiThetaArray,5, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 2, 2, 2, fMixedEventTripletArray[7], PDGCodes, bins[1],fMixedEventTripletMultArray[7], fMixedEventTripletPhiThetaArray,6, *fConfig);
      FillTripletDistributionMEPP(ParticleVector, *itMult, 3, 3, 3, fMixedEventTripletArray[8], PDGCodes, bins[1],fMixedEventTripletMultArray[8], fMixedEventTripletPhiThetaArray,7, *fConfig);

      // Proton Lambda mixing for both proton and lambda used from same event [lambda_same, proton_mixed, proton_mixed]
      //  and [proton_same, lambda_mixed, proton_mixed]
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTEST, 2, 0, 0, fMixedEventTripletArray[9], PDGCodes, bins[1],fMixedEventTripletMultArray[9], fMixedEventTripletPhiThetaArray,8, *fConfig);
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTEST, 0, 2, 0, fMixedEventTripletArray[9], PDGCodes, bins[1],fMixedEventTripletMultArray[9], fMixedEventTripletPhiThetaArray,8, *fConfig);
      
      // Same for antilambda antiproton
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTEST, 3, 1, 1, fMixedEventTripletArray[10], PDGCodes, bins[1],fMixedEventTripletMultArray[10], fMixedEventTripletPhiThetaArray,9, *fConfig);
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTEST, 1, 3, 1, fMixedEventTripletArray[10], PDGCodes, bins[1],fMixedEventTripletMultArray[10], fMixedEventTripletPhiThetaArray,9, *fConfig);
      
      // Test one more mixing, when events taken for mixing are only the events that have tripplets, mixing bothp same and lambda same
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTESTppL, 2, 0, 0, fMixedEventTripletArray[11], PDGCodes, bins[1],fMixedEventTripletMultArray[11], fMixedEventTripletPhiThetaArray,10, *fConfig);
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTESTppL, 0, 2, 0, fMixedEventTripletArray[11], PDGCodes, bins[1],fMixedEventTripletMultArray[11], fMixedEventTripletPhiThetaArray,10, *fConfig);
      
      // Same for antilambda antiproton
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTESTppL, 3, 1, 1, fMixedEventTripletArray[12], PDGCodes, bins[1],fMixedEventTripletMultArray[12], fMixedEventTripletPhiThetaArray,11, *fConfig);
      FillTripletDistributionMEPPTEST(ParticleVector, *itMultTESTppL, 1, 3, 1, fMixedEventTripletArray[12], PDGCodes, bins[1],fMixedEventTripletMultArray[12], fMixedEventTripletPhiThetaArray,11, *fConfig);
      
      


      // Update the particle container with current event
      SetMixedEvent(ParticleVector, &(*itMult));
      SetMixedEventOnlyPLambdaTEST(ParticleVector, &(*itMultTEST));
      SetMixedEventOnlyPPLambdaTEST( ParticleVector, &(*itMultTESTppL));
      */

    }
  }

  

  /*if (fPairCleaner->GetCounter() > 0) {
    if (fConfig->GetUseEventMixing()) {
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent);
    }
    if (fConfig->GetUsePhiSpinning()) {
      fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
    }
  }
  */
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
void AliAnalysisTaskThreeBodyFemtoAOD::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskThreeBodyFemtoAOD::StoreGlobalTrackReference(AliAODTrack *track) {
  // see AliFemtoDreamAnalysis for details

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
    if ((fGTI[trackID])->GetFilterMap()|| fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}


TLorentzVector AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(
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


void AliAnalysisTaskThreeBodyFemtoAOD::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // This function creates a triplet distribution in Q3 bins (defined lower).
  // It requires the particle vector from PairCleaner() and the three indices of particles of interest. So
  // if you want to get distribution for particles that are saved in particle vector as 1 2 3 element, just 
  // call the function with firstSpecies=1,secondSpecies=2,thirdSpecies=3

  float minQ3Part = CalculatePPLTriggerQ3Min(ParticleVector, firstSpecies, secondSpecies, thirdSpecies, PDGCodes);


  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;
  auto Particle3Vector = ParticleVector.begin()+thirdSpecies;

  if(firstSpecies==0){
    if(minQ3Part>fQ3Limit){
        fSameEventOnlyLowestQ3->Fill(minQ3Part);
    }
  }
  else{
    if(minQ3Part>fQ3Limit){
        fSameEventOnlyLowestQ3Anti->Fill(minQ3Part);
    }

  } 

  // Get the PID codes std::vector<int> 
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;
  auto itPDGPar3 = PDGCodes.begin()+thirdSpecies;
  // Get particle masses 
  auto massparticle1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massparticle2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();
  auto massparticle3 = TDatabasePDG::Instance()->GetParticle(*itPDGPar3)->Mass();
  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;
  unsigned int DaughterPart3 = 0;
  if(abs(*itPDGPar1)==3122) DaughterPart1 = 2;
  if(abs(*itPDGPar1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==3122) DaughterPart2 = 2;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGPar3)==3122) DaughterPart3 = 2;
  if(abs(*itPDGPar3)==2212) DaughterPart3 = 1;
  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;
  unsigned int DoThisPair23 = DaughterPart2*10+DaughterPart3;
  unsigned int DoThisPair31 = DaughterPart3*10+DaughterPart1;

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
        //if(firstSpecies==secondSpecies){
        if(DoThisPair12==11){
          Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
        }
        //if(secondSpecies==thirdSpecies){
        if(DoThisPair23==11){
          Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
        }
        //if(thirdSpecies==firstSpecies){
        if(DoThisPair31==11){
          Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
        }
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
        TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part3_LorVec,part1_LorVec);
        // The particles in current methodology are put in bins of:
        //                 Q3=sqrt(q12^2+q23^2+q31^2)
        float Q32 = q12*q12+q23*q23+q31*q31;
        // From 3 pion paper, the q must be multiplied by -1 before taking quare root
        float Q3 = sqrt(-Q32); // the minus from pion paper
        hist->Fill(Q3); 
        hist2d->Fill(Q3,mult+1); 
        if(firstSpecies==0){
          if(Q3<fQ3Limit){
              fSameEventOnlyLowestQ3->Fill(Q3);
          }
        }
        else{
          if(Q3<fQ3Limit){
              fSameEventOnlyLowestQ3Anti->Fill(Q3);
          }

        }
        
      }
    }
  }
}


void AliAnalysisTaskThreeBodyFemtoAOD::FillPairDistributionPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPL, int mult, TH2F* hist2d){
  // Proton Lambda
  // This function is created to just have a simple check that the looping over vectors in this way works
  // To see if it works, please compare the sameEventDistributionPL distribution with the distribution
  // ontained using SetEvent class from femtoDream for proton lambda 
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
      hist2d->Fill(RelativeK,mult+1); 
    }
  }

}
void AliAnalysisTaskThreeBodyFemtoAOD::FillPairDistributionPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPP, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config){
  // Proton Lambda
  // This function is created to just have a simple check that the looping over vectors in this way works
  // To see if it works, please compare the sameEventDistributionPL distribution with the distribution
  // ontained using SetEvent class from femtoDream for proton lambda 
  auto ProtonVector1 = ParticleVector.begin();
  auto ProtonVector2 = ParticleVector.begin();
  // Loop over particles creating pairs
  for (auto iPart1 = ProtonVector1->begin(); iPart1 != ProtonVector1->end(); ++iPart1) {
    auto iPart2 = iPart1+1;
    for (; iPart2 != ProtonVector2->end(); ++iPart2) {
      // Set lorentz vectors
      if(!DeltaEtaDeltaPhi(*iPart1,*iPart2,true,11, fEventTripletPhiThetaArray[40],fEventTripletPhiThetaArray[41], Config)){continue;}

      TLorentzVector part1_LorVec, part2_LorVec;
      part1_LorVec.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
        iPart1->GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(2212)->Mass());
      part2_LorVec.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
        iPart2->GetMomentum().Z(), TDatabasePDG::Instance()->GetParticle(2212)->Mass());
      // Get momentum
      float RelativeK = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec, part2_LorVec);
      // No need to check pair selection because p lambda
      sameEventDistributionPP->Fill(RelativeK);
    }
  }

}

void AliAnalysisTaskThreeBodyFemtoAOD::SetMixedEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  for(unsigned int iSpecies = 0; iSpecies<ParticleVector.size(); iSpecies++){
    if ((ParticleVector.begin()+iSpecies)->size() > 0) {
      (fPartContainer->begin()+iSpecies)->SetEvent(*(ParticleVector.begin()+iSpecies));
    }
  }
}

void AliAnalysisTaskThreeBodyFemtoAOD::SetMixedEventOnlyPLambdaTEST(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda. 
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 0 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
  if ((ParticleVector.begin()+1)->size() > 0 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
}

void AliAnalysisTaskThreeBodyFemtoAOD::SetMixedEventOnlyPPLambdaTEST(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda. 
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 1 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
  if ((ParticleVector.begin()+1)->size() > 1 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
}


void AliAnalysisTaskThreeBodyFemtoAOD::FillTripletDistributionMEPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // Description of function given in AliAnalysisTaskThreeBodyFemtoAOD::FillTripletDistribution
  // In this function, only one particle is used from current event, and the other two - from other two events

  // Current behavior with the mixed events:
  //              1) implemented ONLY IF ME1 and ME2 are the same species!!!!!!!!!!! 
  //              2) the first of the two ME particles: takes Nth event from ME, takes every particle in this event
  //              2) the second of the two ME particles: takes (N+1)th event from ME, takes every particle in this event

  
  if(speciesME1!=speciesME2) {
    AliError("You chose different species ME1 and ME2 for mixing. This is not yet implemented! \n");
  }
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

  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;
  unsigned int DaughterPart3 = 0;
  if(abs(*itPDGParSE)==3122) DaughterPart1 = 2;
  if(abs(*itPDGParSE)==2212) DaughterPart1 = 1;
  if(abs(*itPDGParME1)==3122) DaughterPart2 = 2;
  if(abs(*itPDGParME1)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParME2)==3122) DaughterPart3 = 2;
  if(abs(*itPDGParME2)==2212) DaughterPart3 = 1;
  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;
  unsigned int DoThisPair23 = DaughterPart2*10+DaughterPart3;
  unsigned int DoThisPair31 = DaughterPart3*10+DaughterPart1;
  
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
            //if(speciesSE==speciesME1){
            if(DoThisPair12==11){
              Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
            }
            //if(speciesME1==speciesME2){
            if(DoThisPair23==11){
              Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,false,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
            }
            //if(speciesME2==speciesSE){
            if(DoThisPair31==11){
              Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,false,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
            }
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
            TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part3_LorVec,part1_LorVec);
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



void AliAnalysisTaskThreeBodyFemtoAOD::FillTripletDistributionMEPPTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // Description of function given in AliAnalysisTaskThreeBodyFemtoAOD::FillTripletDistribution
  // In this function, only one particle is used from current event, and the other two - from other two events

  // Current behavior with the mixed events:
  //              1) implemented ONLY IF ME1 and ME2 are the same species!!!!!!!!!!! 
  //              2) the first of the two ME particles: takes Nth event from ME, takes every particle in this event
  //              2) the second of the two ME particles: takes (N+1)th event from ME, takes every particle in this event

  // THIS USES ALL SAME EVENT PARTICLES TO CREATE MIXED EVENTS! [same1+mixed2+mixed3, same2+mixed1+mixed3,same3+mixed2+mixed1]

  // MUST PASS THE PARTCONTAINER FROM SetMixedEventOnlyPLambdaTEST
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

  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;
  unsigned int DaughterPart3 = 0;
  if(abs(*itPDGParSE)==3122) DaughterPart1 = 2;
  if(abs(*itPDGParSE)==2212) DaughterPart1 = 1;
  if(abs(*itPDGParME1)==3122) DaughterPart2 = 2;
  if(abs(*itPDGParME1)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParME2)==3122) DaughterPart3 = 2;
  if(abs(*itPDGParME2)==2212) DaughterPart3 = 1;
  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;
  unsigned int DoThisPair23 = DaughterPart2*10+DaughterPart3;
  unsigned int DoThisPair31 = DaughterPart3*10+DaughterPart1;

  // loop over first particle 
  for (auto iPart1 = ParticleSE->begin(); iPart1 != ParticleSE->end(); ++iPart1) {
    // loop over second particle ...
    for (int iDepth1 = 0; iDepth1 < (int) MixedEvent1Container->GetMixingDepth(); ++iDepth1) {    
      std::vector<AliFemtoDreamBasePart> iEvent2 = MixedEvent1Container->GetEvent(iDepth1);
      for ( auto iPart2 = iEvent2.begin(); iPart2 != iEvent2.end(); ++iPart2) {
        for ( auto iDepth2 = 0; iDepth2 < (int) MixedEvent2Container->GetMixingDepth(); ++iDepth2) {
          if(iDepth1==iDepth2) continue;
          std::vector<AliFemtoDreamBasePart> iEvent3 = MixedEvent2Container->GetEvent(iDepth2);
          for ( auto iPart3 = iEvent3.begin(); iPart3 != iEvent3.end(); ++iPart3) {
            bool Pair12 = true;
            bool Pair23 = true;
            bool Pair31 = true;
            //if(speciesSE==speciesME1){
            if(DoThisPair12==11){
              Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
            }
            //if(speciesME1==speciesME2){
            if(DoThisPair23==11){
              Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,false,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
            }
            //if(speciesME2==speciesSE){
            if(DoThisPair31==11){
              Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,false,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
            }
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
            TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part3_LorVec,part1_LorVec);
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

void AliAnalysisTaskThreeBodyFemtoAOD::FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config){
  //  Check if reproduces the FemtoDream framework result
  auto ParticleSE = ParticleVector.begin()+speciesSE;
  auto MixedEvent1Container = fPartContainer.begin()+speciesME1;

  // Get the PID codes std::vector<int> 
  auto itPDGParSE = PDGCodes.begin()+speciesSE;
  auto itPDGParME1 = PDGCodes.begin()+speciesME1;

  // Get particle masses 
  auto massParticleSE = TDatabasePDG::Instance()->GetParticle(*itPDGParSE)->Mass();
  auto massParticleME1 = TDatabasePDG::Instance()->GetParticle(*itPDGParME1)->Mass();
  
  // loop over first particle 
  for (auto iPart1 = ParticleSE->begin(); iPart1 != ParticleSE->end(); ++iPart1) {
    // loop over second particle ...
    for (int iDepth1 = 0; iDepth1 < (int) MixedEvent1Container->GetMixingDepth(); ++iDepth1) {    
      std::vector<AliFemtoDreamBasePart> iEvent2 = MixedEvent1Container->GetEvent(iDepth1);
      for ( auto iPart2 = iEvent2.begin(); iPart2 != iEvent2.end(); ++iPart2) {
        if(speciesSE==0 && speciesME1 ==0){
          if(!DeltaEtaDeltaPhi(*iPart1,*iPart2,false,11, fEventTripletPhiThetaArray[40],fEventTripletPhiThetaArray[41], Config)){continue;}
        }
        // Now we have the three particles, lets create their Lorentz vectors  
        TLorentzVector part1_LorVec, part2_LorVec;
        part1_LorVec.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
          iPart1->GetMomentum().Z(),massParticleSE);
        part2_LorVec.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
          iPart2->GetMomentum().Z(),massParticleME1);
        // Get momentum
        float RelativeK = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec, part2_LorVec);
        // No need to check pair selection because p lambda
        hist->Fill(RelativeK);  
        if(speciesSE==0 && speciesME1 ==2) hist2d->Fill(RelativeK,mult+1); 
      }
    }
  }
}


void AliAnalysisTaskThreeBodyFemtoAOD::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // Description of function given in AliAnalysisTaskThreeBodyFemtoAOD::FillTripletDistribution
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
  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;
  unsigned int DaughterPart3 = 0;
  if(abs(*itPDGParSE1)==3122) DaughterPart1 = 2;
  if(abs(*itPDGParSE1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGParSE2)==3122) DaughterPart2 = 2;
  if(abs(*itPDGParSE2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParME)==3122) DaughterPart3 = 2;
  if(abs(*itPDGParME)==2212) DaughterPart3 = 1;
  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;
  unsigned int DoThisPair23 = DaughterPart2*10+DaughterPart3;
  unsigned int DoThisPair31 = DaughterPart3*10+DaughterPart1;

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
          //if(speciesSE==speciesME1){
          if(DoThisPair12==11){
            Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
          }
          //if(speciesME1==speciesME2){
          if(DoThisPair23==11){
            Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,false,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
          }
          //if(speciesME2==speciesSE){
          if(DoThisPair31==11){
            Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,false,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
          }
          if(!Pair12||!Pair23||!Pair31) {continue;}

          TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
          part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
          iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE1,2)));
          part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
          iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleSE2,2)));
          part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
          iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME,2)));
          // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
          TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part1_LorVec,part2_LorVec);
          TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part2_LorVec,part3_LorVec);
          TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part3_LorVec,part1_LorVec);
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

bool AliAnalysisTaskThreeBodyFemtoAOD::DeltaEtaDeltaPhi(
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

bool AliAnalysisTaskThreeBodyFemtoAOD::MyLovely3BodyTrigger(AliAODEvent *evt ,  bool isMC, std::vector<int> PDGCodes){

  if (!fEventCutsTrigger->isSelected(fEvent)) {
    return false;
  }

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return false;
    }
    StoreGlobalTrackReference(track);
  }

  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    fTrack->SetTrack(track);
    if (fTrackCutsTrigger->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fAntiTrackCutTrigger->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
  }

  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;

  fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iv0 = 0;
      iv0 < static_cast<TClonesArray *>(evt->GetV0s())->GetEntriesFast();
      ++iv0) {
    AliAODv0* v0 = evt->GetV0(iv0);
    fv0->Setv0(evt, v0);
    if (fv0CutsTrigger->isSelected(fv0)) {
      Lambdas.push_back(*fv0);
    }
    if (fAntiv0CutsTrigger->isSelected(fv0)) {
      AntiLambdas.push_back(*fv0);
    }
  }

  if(fRunppp0ppL1==true){
    if((Protons.size()<2||Lambdas.size()<1)&&(AntiProtons.size()<2||AntiLambdas.size()<1)){
      return false;
    }
  }
  if(fRunppp0ppL1==false){
    if((Protons.size()<3)&&(AntiProtons.size()<3)){
        return false;
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
  std::vector<std::vector<AliFemtoDreamBasePart>> ParticleVector = fPairCleaner->GetCleanParticles();


  if(fRunppp0ppL1==true){
    if((ParticleVector[0].size()<2||ParticleVector[2].size()<1)&&(ParticleVector[1].size()<2||ParticleVector[3].size()<1)){
      return false;
    }
  }
  if(fRunppp0ppL1==false){
    if((ParticleVector[0].size()<3)&&(ParticleVector[1].size()<3)){
      return false;
    }
  }

  float totalTripletsSumPartAntipart = 0.;
  float minQ3Part = 1e3, minQ3AntiPart = 1e3;

  if(fRunppp0ppL1==true) minQ3Part = CalculatePPLTriggerQ3Min(ParticleVector, 0, 2, 0, PDGCodes);
  if(fRunppp0ppL1==false) minQ3Part = CalculatePPLTriggerQ3Min(ParticleVector, 0, 0, 0, PDGCodes);

  totalTripletsSumPartAntipart = (float)ftotalTripletCount;
  if(fRunppp0ppL1==true) minQ3AntiPart = CalculatePPLTriggerQ3Min(ParticleVector, 1, 3, 1, PDGCodes);
  if(fRunppp0ppL1==false)  minQ3AntiPart = CalculatePPLTriggerQ3Min(ParticleVector, 1, 1, 1, PDGCodes);
  totalTripletsSumPartAntipart +=  (float)ftotalTripletCount-0.01;

  if(minQ3Part<fQ3Limit || minQ3AntiPart <fQ3Limit ){
    float randomNumber = fRandomGen->Uniform(0.,1.);
    if(randomNumber>fQ3LimitFraction){
      return false;
    }
    fTripletsPerCollision->Fill(totalTripletsSumPartAntipart);
    fWhichSample->Fill(0.5);
    return true;
  }


  if(minQ3Part<fQ3LimitSample || minQ3AntiPart <fQ3LimitSample ){
    float randomNumber = fRandomGen->Uniform(0.,1.);
    if(randomNumber>fQ3LimitSampleFraction){
      return false;
    }
    fTripletsPerCollision->Fill(totalTripletsSumPartAntipart);
    fWhichSample->Fill(1.5);
    return true;
  }


  if(minQ3Part<fQ3LimitSample2 || minQ3AntiPart <fQ3LimitSample2 ){
    float randomNumber = fRandomGen->Uniform(0.,1.);
    if(randomNumber>fQ3LimitSampleFraction2){
      return false;
    }
    fTripletsPerCollision->Fill(totalTripletsSumPartAntipart);
    fWhichSample->Fill(2.5);
    return true;
  }


  return false;
}

double AliAnalysisTaskThreeBodyFemtoAOD::CalculatePPLTriggerQ3Min(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, std::vector<int> PDGCodes ){
  // This function creates a triplet distribution in Q3 bins (defined lower).
  // It requires the particle vector from PairCleaner() and the three indices of particles of interest. So
  // if you want to get distribution for particles that are saved in particle vector as 1 2 3 element, just 
  // call the function with firstSpecies=1,secondSpecies=2,thirdSpecies=3


  double minQ3InEvent = 1000.; // min q3 tripelt present in the sample
  ftotalTripletCount = 0;

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

        // Now we have the three particles, lets create their Lorentz vectors  
        TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
        part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
        iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massparticle1,2)));
        part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
        iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massparticle2,2)));
        part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
        iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massparticle3,2)));
        
        // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
        TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoAOD::RelativePairMomentum(part3_LorVec,part1_LorVec);
        // The particles in current methodology are put in bins of:
        //                 Q3=sqrt(q12^2+q23^2+q31^2)
        float Q32 = q12*q12+q23*q23+q31*q31;
        // From 3 pion paper, the q must be multiplied by -1 before taking quare root
        float Q3 = sqrt(-Q32); // the minus from pion paper  
        ftotalTripletCount += 1; 
        if (Q3<minQ3InEvent){
           minQ3InEvent= Q3;   
         }
      }
    }
  }
  
  return minQ3InEvent;
}
