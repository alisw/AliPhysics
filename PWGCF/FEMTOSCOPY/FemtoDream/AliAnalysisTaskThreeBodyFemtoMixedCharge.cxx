/*
 * AliAnalysisTaskThreeBodyFemtoMixedCharge.cxx
 *
 *  Created on: May 13, 2019
 *      Author: Laura Serksnyte 
 */
#include "AliAnalysisTaskThreeBodyFemtoMixedCharge.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliNanoAODTrack.h"
#include "Riostream.h"

ClassImp(AliAnalysisTaskThreeBodyFemtoMixedCharge)
AliAnalysisTaskThreeBodyFemtoMixedCharge::AliAnalysisTaskThreeBodyFemtoMixedCharge()
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
      fRunThreeBody(true),
      fRunPlotPhiTheta(true),
      fRunPlotMult(true),
      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fClosePairRejectionPPPorPPL(false),
      fisMC(false),
      fRun2Body(false),
      fSame2Mixed1Choice(false),
      fQ3LimitForDeltaPhiDeltaEta(1.),
      fMixingChoice(0),  
      fWhichTripletsToRun(11),  
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fPartContainerTESTppp(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskThreeBodyFemtoMixedCharge::AliAnalysisTaskThreeBodyFemtoMixedCharge(const char* name, bool isMC)
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
      fRunThreeBody(true),
      fRunPlotPhiTheta(true),
      fRunPlotMult(true),
      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fClosePairRejectionPPPorPPL(false),
      fisMC(isMC),
      fRun2Body(false),
      fSame2Mixed1Choice(false),
      fQ3LimitForDeltaPhiDeltaEta(1.),
      fMixingChoice(0),
      fWhichTripletsToRun(11),  
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fPartContainerTESTppp(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
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

AliAnalysisTaskThreeBodyFemtoMixedCharge::~AliAnalysisTaskThreeBodyFemtoMixedCharge() {
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

void AliAnalysisTaskThreeBodyFemtoMixedCharge::UserCreateOutputObjects() {
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
    fResultsThreeBody->SetName("ResultsThreeBody");
    // Same event
    fSameEvent = new TList();
    fSameEvent->SetOwner();
    fSameEvent->SetName("SameEvent");

    fSameEventTripletArray = new TH1F*[24];
    TString histTitlesSame[24] ={"sameEventDistributionPL","sameEventDistributionPPL","sameEventDistributionAPAPAL",
      "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPLL","sameEventDistributionAPALAL",
      "sameEventDistributionLLL","sameEventDistributionALALAL", "sameEventDistributionppSameLMixed", "sameEventDistributionpLSamepMixed",
      "sameEventDistributionppSamepMixed","sameEventDistributionPP","sameEventDistributionapapSameaLMixed", "sameEventDistributionapaLSameapMixed",
      "sameEventDistributionapapSameapMixed", "sameEventDistributionpLSameLMixed", "sameEventDistributionLLSamepMixed",
      "sameEventDistributionapaLSameaLMixed", "sameEventDistributionaLaLSameapMixed","sameEventDistributionLLSameLMixed", "sameEventDistributionaLaLSameaLMixed",
      "sameEventDistributionpapSamepMixed", "sameEventDistributionappSameapMixed"};
    fSameEventTripletArray[0] = new TH1F(histTitlesSame[0],histTitlesSame[0],1000,0, 1);
    fSameEvent->Add(fSameEventTripletArray[0]);
    for (int i = 1; i < 12; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }
    fSameEventTripletArray[12] = new TH1F(histTitlesSame[12],histTitlesSame[12],1000,0, 1);
    fSameEvent->Add(fSameEventTripletArray[12]);
    for (int i = 13; i <24 ; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }
    // Same event multiplicity dist
    fSameEventMult = new TList();
    fSameEventMult->SetOwner();
    fSameEventMult->SetName("SameEventMult");

    fSameEventTripletMultArray = new TH2F*[23];
    TString histTitlesSameMult[23] ={"sameEventDistributionMultPL","sameEventDistributionMultPPL","sameEventDistributionMultAPAPAL",
      "sameEventDistributionMultPPP", "sameEventDistributionMultAPAPAP", "sameEventDistributionMultPLL","sameEventDistributionMultAPALAL",
      "sameEventDistributionMultLLL","sameEventDistributionMultALALAL",  "sameEventDistributionMultppSameLMixed", "sameEventDistributionMultpLSamepMixed",
      "sameEventDistributionMultppSamepMixed","sameEventDistributionMultapapSameaLMixed", "sameEventDistributionMultapaLSameapMixed",
      "sameEventDistributionMultapapSameapMixed", "sameEventDistributionMultpLSameLMixed", "sameEventDistributionMultLLSamepMixed",
      "sameEventDistributionMultapaLSameaLMixed", "sameEventDistributionMultaLaLSameapMixed","sameEventDistributionMultLLSameLMixed", "sameEventDistributionMultaLaLSameaLMixed",
      "sameEventDistributionMultpapSamepMixed", "sameEventDistributionMultappSameapMixed"};
    
    fSameEventTripletMultArray[0] = new TH2F(histTitlesSameMult[0],histTitlesSameMult[0],1000,0, 1,26,1,27);
    fSameEventMult->Add(fSameEventTripletMultArray[0]);
    for (int i = 1; i < 23; ++i) {
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

       // end mixed event resolution matrix

    fResultsThreeBody->Add(fSameEvent);
    fResultsThreeBody->Add(fMixedEvent);
    fResultsThreeBody->Add(fSameEventMult);
    fResultsThreeBody->Add(fMixedEventMult);

    // Same event phi theta distribution
    fSameEventPhiTheta = new TList();
    fSameEventPhiTheta->SetOwner();
    fSameEventPhiTheta->SetName("SameEventPhiTheta");

    fSameEventTripletPhiThetaArray = new TH2F*[46];
    TString histTitlesSamePhiEta[23] = {"sameEventPhiEtaPPL","sameEventPhiEtaAPAPAL",
      "sameEventPhiEtaPPP", "sameEventPhiEtaAPAPAP", "sameEventPhiEtaPLL","sameEventPhiEtaAPALAL",
      "sameEventPhiEtaLLL","sameEventPhiEtaALALAL", "sameEventDistributionPhiEtappSameLMixed", "sameEventDistributionPhiEtapLSamepMixed",
      "sameEventDistributionPhiEtappSamepMixed","sameEventDistributionPhiEtaapapSameaLMixed", "sameEventDistributionPhiEtaapaLSameapMixed",
      "sameEventDistributionPhiEtaapapSameapMixed", "sameEventDistributionPhiEtapLSameLMixed", "sameEventDistributionPhiEtaLLSamepMixed",
      "sameEventDistributionPhiEtaapaLSameaLMixed", "sameEventDistributionPhiEtaaLaLSameapMixed","sameEventDistributionPhiEtaLLSameLMixed", "sameEventDistributionPhiEtaaLaLSameaLMixed", "TRASH", "sameEventDistributionPhiEtpapSamepMixed", "sameEventDistributionPhiEtappSameapMixed"};
    for(int i=0;i<23;i++){
      fSameEventTripletPhiThetaArray[i] = new TH2F(histTitlesSamePhiEta[i]+"Before",histTitlesSamePhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventTripletPhiThetaArray[23+i] = new TH2F(histTitlesSamePhiEta[i]+"After",histTitlesSamePhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[i]);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[23+i]);
    }

    // Mixed event phi theta distribution
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



    fResultsThreeBody->Add(fSameEventPhiTheta);
    fResultsThreeBody->Add(fMixedEventPhiTheta);




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
      fPartContainerTESTppp.push_back(MultContainer);
    }


}

void AliAnalysisTaskThreeBodyFemtoMixedCharge::UserExec(Option_t *option) {
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
  }

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
    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    int bins[2] = { 0, 0 };
    float ZVtx = fEvent->GetZVertex();
    float Mult = fEvent->GetMultiplicity();
    fPartColl->FindBin(ZVtx, Mult, bins);
    // Same event distribution -------------------------------------------------------------------------------
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles();
    if(fWhichTripletsToRun ==1||fWhichTripletsToRun ==11){
      // proton proton lambda
      FillTripletDistribution( ParticleVector, 2, 1, 1, fSameEventTripletArray[1],PDGCodes, bins[1],fSameEventTripletMultArray[1], fSameEventTripletPhiThetaArray,0, *fConfig);//, fInvMassSame[0]);
       // antiproton antiproton antilambad
      FillTripletDistribution( ParticleVector, 3, 0, 0, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletPhiThetaArray,1, *fConfig);
    }
    if(fWhichTripletsToRun ==0||fWhichTripletsToRun ==11){
      // proton proton proton 
      FillTripletDistribution( ParticleVector, 1, 0, 0, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletPhiThetaArray,2, *fConfig);
      // antiproton antiproton antiproton 
      FillTripletDistribution( ParticleVector, 0, 1, 1, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletPhiThetaArray,3, *fConfig);
    }

     // Mixed event distribution
    if (!(bins[0] == -99 || bins[1] == -99)) {

      auto itZVtx = fPartContainer.begin()+ bins[0];
      auto itMult = itZVtx->begin() + bins[1];

      // Same2Mixed1
      if(fWhichTripletsToRun ==0||fWhichTripletsToRun ==11){
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 1, fSameEventTripletArray[11], PDGCodes, bins[1],fSameEventTripletMultArray[11], fSameEventTripletPhiThetaArray, 10, *fConfig);//, fQ3VskDistributionsArrayq12[2],fQ3VskDistributionsArrayq23[2]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 0, fSameEventTripletArray[15], PDGCodes, bins[1],fSameEventTripletMultArray[14], fSameEventTripletPhiThetaArray, 13, *fConfig);//, fQ3VskDistributionsArrayq12[5],fQ3VskDistributionsArrayq23[5]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 1, 0, fSameEventTripletArray[22], PDGCodes, bins[1],fSameEventTripletMultArray[21], fSameEventTripletPhiThetaArray, 21, *fConfig);//, fQ3VskDistributionsArrayq12[2],fQ3VskDistributionsArrayq23[2]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 0, 1, fSameEventTripletArray[23], PDGCodes, bins[1],fSameEventTripletMultArray[22], fSameEventTripletPhiThetaArray, 22, *fConfig);//, fQ3VskDistributionsArrayq12[5],fQ3VskDistributionsArrayq23[5]);
      }
      if(fWhichTripletsToRun ==1||fWhichTripletsToRun ==11){
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 2, fSameEventTripletArray[9], PDGCodes, bins[1],fSameEventTripletMultArray[9], fSameEventTripletPhiThetaArray, 8, *fConfig);//, fQ3VskDistributionsArrayq12[0],fQ3VskDistributionsArrayq23[0]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 2, 1, fSameEventTripletArray[10], PDGCodes, bins[1],fSameEventTripletMultArray[10], fSameEventTripletPhiThetaArray, 9, *fConfig);//, fQ3VskDistributionsArrayq12[1],fQ3VskDistributionsArrayq23[1]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 3, fSameEventTripletArray[13], PDGCodes, bins[1],fSameEventTripletMultArray[12], fSameEventTripletPhiThetaArray, 11, *fConfig);//, fQ3VskDistributionsArrayq12[3],fQ3VskDistributionsArrayq23[3]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 3, 0, fSameEventTripletArray[14], PDGCodes, bins[1],fSameEventTripletMultArray[13], fSameEventTripletPhiThetaArray, 12, *fConfig);//, fQ3VskDistributionsArrayq12[4],fQ3VskDistributionsArrayq23[4]);
      }
    
      // Normal mixing

      if(fWhichTripletsToRun ==0||fWhichTripletsToRun ==11){
        FillTripletDistributionME(ParticleVector, *itMult, 1, 0, 0, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletPhiThetaArray,2, *fConfig);
        FillTripletDistributionME(ParticleVector, *itMult, 0, 1, 1, fMixedEventTripletArray[4], PDGCodes, bins[1],fMixedEventTripletMultArray[4], fMixedEventTripletPhiThetaArray,3, *fConfig);
      }
      if(fWhichTripletsToRun ==1||fWhichTripletsToRun ==11){
        FillTripletDistributionME(ParticleVector, *itMult, 3, 0, 0, fMixedEventTripletArray[1], PDGCodes, bins[1],fMixedEventTripletMultArray[1], fMixedEventTripletPhiThetaArray,0, *fConfig);
        FillTripletDistributionME(ParticleVector, *itMult, 2, 1, 1, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletPhiThetaArray,1, *fConfig);
      }



      // Update the particle container with current event
      if(fMixingChoice==0){
        SetMixedEvent(ParticleVector, &(*itMult));
      }
      else if(fMixingChoice==2){
        SetMixedEventOnlyPPLambdaTEST(ParticleVector, &(*itMult));
      }else if(fMixingChoice==3){
        SetMixedEventOnlyPPPTEST(ParticleVector, &(*itMult));
      }
    }
  }

  
  if(fRun2Body){

    if (fPairCleaner->GetCounter() > 0) {
      if (fConfig->GetUseEventMixing()) {
        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                            fEvent);
      }
      if (fConfig->GetUsePhiSpinning()) {
        fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
      }
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
void AliAnalysisTaskThreeBodyFemtoMixedCharge::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskThreeBodyFemtoMixedCharge::StoreGlobalTrackReference(AliVTrack *track) {
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


TLorentzVector AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(
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


void AliAnalysisTaskThreeBodyFemtoMixedCharge::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
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

  bool dontDoIfDiffCharge1, dontDoIfDiffCharge2, dontDoIfDiffCharge3;
  int pairID12=0,  pairID23=0, pairID31=0;
  pairID12 = firstSpecies*10+secondSpecies;
  pairID23 = secondSpecies*10+thirdSpecies;
  pairID31 = thirdSpecies*10+firstSpecies;

  if(pairID12==0||pairID12==11||pairID12==2||pairID12==20||pairID12==13||pairID12==31){
    dontDoIfDiffCharge1 = false;
  } else {
    dontDoIfDiffCharge1 = true;
  }

  if(pairID23==0||pairID23==11||pairID23==2||pairID23==20||pairID23==13||pairID23==31){
    dontDoIfDiffCharge2 = false;
  } else {
    dontDoIfDiffCharge2 = true;
  }

  if(pairID31==0||pairID31==11||pairID31==2||pairID31==20||pairID31==13||pairID31==31){
    dontDoIfDiffCharge3 = false;
  } else {
    dontDoIfDiffCharge3 = true;
  }

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
        TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part3_LorVec,part1_LorVec);
        // The particles in current methodology are put in bins of:
        //                 Q3=sqrt(q12^2+q23^2+q31^2)
        float Q32 = q12*q12+q23*q23+q31*q31;
        // From 3 pion paper, the q must be multiplied by -1 before taking quare root
        float Q3 = sqrt(-Q32); // the minus from pion paper

        bool Pair12 = true;
        bool Pair23 = true;
        bool Pair31 = true;
        
        if(!fturnoffClosePairRejectionCompletely){

          if(fClosePairRejectionForAll){
            if(!dontDoIfDiffCharge1){
              Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
            }
            if(!dontDoIfDiffCharge2){
              Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
            }
            if(!dontDoIfDiffCharge3){
              Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
            }
          }
          if(!fClosePairRejectionForAll){
            if(fClosePairRejectionPPPorPPL){
              if(DoThisPair12==11){ 
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
              if(DoThisPair23==11){ 
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
              if(DoThisPair31==11){ 
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
            }

            if(!fClosePairRejectionPPPorPPL){
              if(DoThisPair12==21||DoThisPair12==12){ 
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
              if(DoThisPair23==21||DoThisPair23==12){ 
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
              if(DoThisPair31==21||DoThisPair31==12){ 
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
            }
          }
        }
        

        if(!Pair12||!Pair23||!Pair31) {continue;}


        hist->Fill(Q3); 
        hist2d->Fill(Q3,mult+1); 

        
      }
    }
  }
}


void AliAnalysisTaskThreeBodyFemtoMixedCharge::SetMixedEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  for(unsigned int iSpecies = 0; iSpecies<ParticleVector.size(); iSpecies++){
    if ((ParticleVector.begin()+iSpecies)->size() > 0) {
      (fPartContainer->begin()+iSpecies)->SetEvent(*(ParticleVector.begin()+iSpecies));
    }
  }
}


void AliAnalysisTaskThreeBodyFemtoMixedCharge::SetMixedEventOnlyPPLambdaTEST(
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

void AliAnalysisTaskThreeBodyFemtoMixedCharge::SetMixedEventOnlyPPPTEST(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda. 
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 2) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
  }
  if ((ParticleVector.begin()+1)->size() > 2) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
  }
}


void AliAnalysisTaskThreeBodyFemtoMixedCharge::FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed){
  // Description of function given in AliAnalysisTaskThreeBodyFemtoMixedCharge::FillTripletDistribution
  // In this function, only one particle is used from current event, and the other two - from other two events

  // Current behavior with the mixed events:
  //              1) implemented ONLY IF ME1 and ME2 are the same species!!!!!!!!!!! 
  //              2) the first of the two ME particles: takes Nth event from ME, takes every particle in this event
  //              2) the second of the two ME particles: takes (N+1)th event from ME, takes every particle in this event

  
  if(speciesME1!=speciesME2) {
    AliError("You chose different species ME1 and ME2 for mixing. This is not yet implemented! \n");
  }

  if(fMixingChoice>0){
    if((speciesSE==0&&speciesME1==0&&speciesME2==0)||(speciesSE==1&&speciesME1==1&&speciesME2==1)){
      if(ParticleVector[speciesSE].size()<3) return;
    }
    if((speciesSE==2&&speciesME1==0&&speciesME2==0)||(speciesSE==3&&speciesME1==1&&speciesME2==1)){
      if(ParticleVector[speciesSE].size()<2||ParticleVector[speciesME1].size()<1) return;
    }
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

  bool dontDoIfDiffCharge1, dontDoIfDiffCharge2, dontDoIfDiffCharge3;
  int pairID12=0,  pairID23=0, pairID31=0;
  pairID12 = speciesSE*10+speciesME1;
  pairID23 = speciesME1*10+speciesME2;
  pairID31 = speciesME2*10+speciesSE;
  if(pairID12==0||pairID12==11||pairID12==2||pairID12==20||pairID12==13||pairID12==31){
    dontDoIfDiffCharge1 = false;
  } else {
    dontDoIfDiffCharge1 = true;
  }

  if(pairID23==0||pairID23==11||pairID23==2||pairID23==20||pairID23==13||pairID23==31){
    dontDoIfDiffCharge2 = false;
  } else {
    dontDoIfDiffCharge2 = true;
  }


  if(pairID31==0||pairID31==11||pairID31==2||pairID31==20||pairID31==13||pairID31==31){
    dontDoIfDiffCharge3 = false;
  } else {
    dontDoIfDiffCharge3 = true;
  }



  
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


            TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
            part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
            iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE,2)));
            part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
            iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleME1,2)));
            part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
            iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME2,2)));
            // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
            TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part3_LorVec,part1_LorVec);
            // The particles in current methodology are put in bins of:
            //                 Q3=sqrt(q12^2+q23^2+q31^2)
            float Q32 = q12*q12+q23*q23+q31*q31;
            // From 3 pion paper, the q must be multiplied by -1 before taking quare root
            float Q3 = sqrt(-Q32); // the minus from pion paper


            bool Pair12 = true;
            bool Pair23 = true;
            bool Pair31 = true;
                    

            if(!fturnoffClosePairRejectionCompletely){

              if(fClosePairRejectionForAll){
                if(!dontDoIfDiffCharge1){
                  Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                 }
                if(!dontDoIfDiffCharge2){
                  Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                }
                if(!dontDoIfDiffCharge3){
                  Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                }
              }
              if(!fClosePairRejectionForAll){
                if(fClosePairRejectionPPPorPPL){
                  if(DoThisPair12==11){ 
                    Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                  }
                  if(DoThisPair23==11){ 
                    Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                  }
                  if(DoThisPair31==11){ 
                    Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                  }
                }
                if(!fClosePairRejectionPPPorPPL){
                  if(DoThisPair12==21||DoThisPair12==12){ 
                    Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                  }
                  if(DoThisPair23==21||DoThisPair23==12){ 
                    Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                  }
                  if(DoThisPair31==21||DoThisPair31==12){ 
                    Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config);
                  }
                }
              }

            }
        
            if(!Pair12||!Pair23||!Pair31) {continue;}
           

            // Now we have the three particles, lets create their Lorentz vectors  

            hist->Fill(Q3); 
            hist2d->Fill(Q3,mult+1); 

          }
        }  
      }
    }
  }
}


void AliAnalysisTaskThreeBodyFemtoMixedCharge::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23){
  // Description of function given in AliAnalysisTaskThreeBodyFemtoMixedCharge::FillTripletDistribution
  // In this function, two particles are used from current event, and one - from other event
  auto ParticleSE1 = ParticleVector.begin()+speciesSE1;
  auto ParticleSE2 = ParticleVector.begin()+speciesSE2;
  auto MixedEventContainer = fPartContainer.begin()+speciesME;

  if(fSame2Mixed1Choice){
    if((speciesSE1==0&&speciesSE2==0&&speciesME==0)||(speciesSE1==1&&speciesSE2==1&&speciesME==1)){
      if(ParticleVector[speciesSE1].size()<3) return;
    }
    if((speciesSE1==0&&speciesSE2==0&&speciesME==2)||(speciesSE1==1&&speciesSE2==1&&speciesME==3)){
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesME].size()<1) return;
    }
    if((speciesSE1==0&&speciesSE2==2&&speciesME==0)||(speciesSE1==1&&speciesSE2==3&&speciesME==1)){
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesSE2].size()<1) return;
    }
  }

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


  bool dontDoIfDiffCharge1, dontDoIfDiffCharge2, dontDoIfDiffCharge3;
  int pairID12=0,  pairID23=0, pairID31=0;
  pairID12 = speciesSE1*10+speciesSE2;
  pairID23 = speciesSE2*10+speciesME;
  pairID31 = speciesME*10+speciesSE1;
  if(pairID12==0||pairID12==11||pairID12==2||pairID12==20||pairID12==13||pairID12==31){
    dontDoIfDiffCharge1 = false;
  } else {
    dontDoIfDiffCharge1 = true;
  }

  if(pairID23==0||pairID23==11||pairID23==2||pairID23==20||pairID23==13||pairID23==31){
    dontDoIfDiffCharge2 = false;
  } else {
    dontDoIfDiffCharge2 = true;
  }


  if(pairID31==0||pairID31==11||pairID31==2||pairID31==20||pairID31==13||pairID31==31){
    dontDoIfDiffCharge3 = false;
  } else {
    dontDoIfDiffCharge3 = true;
  }


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


          TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
          part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(), 
          iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE1,2)));
          part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(), 
          iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleSE2,2)));
          part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(), 
          iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME,2)));
          // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
          TLorentzVector q12 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part1_LorVec,part2_LorVec);
          TLorentzVector q23 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part2_LorVec,part3_LorVec);
          TLorentzVector q31 = AliAnalysisTaskThreeBodyFemtoMixedCharge::RelativePairMomentum(part3_LorVec,part1_LorVec);
          // The particles in current methodology are put in bins of:
          //                 Q3=sqrt(q12^2+q23^2+q31^2)
          float Q32 = q12*q12+q23*q23+q31*q31;
          // From 3 pion paper, the q must be multiplied by -1 before taking quare root
          float Q3 = sqrt(-Q32); // the minus from pion paper


          bool Pair12 = true;
          bool Pair23 = true;
          bool Pair31 = true;
          if(!fturnoffClosePairRejectionCompletely){

            if(fClosePairRejectionForAll){
              if(!dontDoIfDiffCharge1){
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
              if(!dontDoIfDiffCharge2){
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
              if(!dontDoIfDiffCharge3){
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
              }
            }
            if(!fClosePairRejectionForAll){
              if(fClosePairRejectionPPPorPPL){
                if(DoThisPair12==11){ 
                  Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
                }
                if(DoThisPair23==11){ 
                  Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
                }
                if(DoThisPair31==11){ 
                  Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
                }
              }
              if(!fClosePairRejectionPPPorPPL){
                if(DoThisPair12==21||DoThisPair12==12){ 
                  Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
                }
                if(DoThisPair23==21||DoThisPair23==12){ 
                  Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
                }
                if(DoThisPair31==21||DoThisPair31==12){ 
                  Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[23+phiEtaHistNo],Config);
                }
              }
            }

          }

          if(!Pair12||!Pair23||!Pair31) {continue;}


          hist->Fill(Q3); 
          hist2d->Fill(Q3,mult+1); 

        }
      }  
    }
  }
}

bool AliAnalysisTaskThreeBodyFemtoMixedCharge::DeltaEtaDeltaPhi(
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
      if(fRunPlotPhiTheta){
        beforeHist->Fill(dphiAvg/ (float) size, deta);
      }
      if (pass) {
        if ((dphiAvg / (float) size) * (dphiAvg / (float) size) / fDeltaPhiSqMax
            + deta * deta / fDeltaEtaSqMax < 1.) {
          pass = false;
        }
        else{
          if(fRunPlotPhiTheta){
            afterHist->Fill(dphiAvg/ (float) size, deta);
          }
        }
      }
    }
  }
  return pass;
}


