/*
 * AliAnalysisTaskThreeBodyFemto.cxx
 *
 *  Created on: May 13, 2019
 *      Author: Laura Serksnyte 
 */
#include "AliAnalysisTaskThreeBodyFemto.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliNanoAODTrack.h"
#include "Riostream.h"

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
      fSameEvent(nullptr),  
      fMixedEvent(nullptr), 
      fSameEventMult(nullptr), 
      fMixedEventMult(nullptr), 
      fSameEventPhiTheta(nullptr), 
      fMixedEventPhiTheta(nullptr), 
      fSameEventTripletResolution(nullptr), 
      fSameEventTripletResolutionAll(nullptr), 
      fMixedEventTripletResolution(nullptr), 
      fMixedEventTripletResolutionAll(nullptr), 
      fQ3Vskq12(nullptr), 
      fQ3Vskq12Mixed(nullptr), 
      fQ3Vskq23(nullptr), 
      fQ3Vskq23Mixed(nullptr), 
      fOtherHistos(nullptr),   
      fQ3VsPairmTListSE(nullptr),  
      fQ3VsPairmTListME(nullptr),  
      fQ3VsPairmTList2SE1ME(nullptr),  
      fQ3VsTripletmTListSE(nullptr),  
      fQ3VsTripletmTListME(nullptr),  
      fQ3VsTripletmTList2SE1ME(nullptr),  
      fInvMassTripletSame(nullptr),   
      fInvMassTripletMixed(nullptr),  
      fRunThreeBody(true),
      fRunPlotInvMassTriplet(true),
      fRunPlotQ3Vsq(true),
      fRunPlotPhiTheta(true),
      fRunPlotOtherHistos(true),
      fRunPlotMult(true),
      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fClosePairRejectionPPPorPPL(false),
      fisMC(false),
      fRun2Body(false),
      fSame2Mixed1Choice(false),
      fQ3LimitForDeltaPhiDeltaEta(0.4),
      fMixingChoice(0),  
      fWhichTripletsToRun(11),  
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fSameEventTripletArrayResolution(nullptr),
      fSameEventTripletArrayResolutionAll(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fPartContainerTESTppp(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fMixedEventTripletArrayResolution(nullptr),
      fMixedEventTripletArrayResolutionAll(nullptr),
      fQ3VskDistributionsArrayq12(nullptr),
      fQ3VskDistributionsArrayq12Mixed(nullptr),
      fQ3VskDistributionsArrayq23(nullptr),
      fQ3VskDistributionsArrayq23Mixed(nullptr),
      fDoubletVsTrippletPPP(nullptr),
      fInvMassSame(nullptr),
      fInvMassMixed(nullptr),
      fTripletInvMassDet(nullptr), 
      fTripletInvMassPDG(nullptr), 
      fTripletInvMassDetMixed(nullptr), 
      fTripletInvMassPDGMixed(nullptr), 
      fTripletInvMassDetAnti(nullptr), 
      fTripletInvMassPDGAnti(nullptr), 
      fTripletInvMassDetMixedAnti(nullptr), 
      fTripletInvMassPDGMixedAnti(nullptr), 
      fQ3VskmTPairSE(nullptr), 
      fQ3VskmTPairME(nullptr), 
      fQ3VskmTPair2SEME1(nullptr), 
      fQ3VskmTTripletSE(nullptr), 
      fQ3VskmTTripletME(nullptr), 
      fQ3VskmTTriplet2SEME1(nullptr), 
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
      fSameEvent(nullptr),  
      fMixedEvent(nullptr), 
      fSameEventMult(nullptr), 
      fMixedEventMult(nullptr), 
      fSameEventPhiTheta(nullptr), 
      fMixedEventPhiTheta(nullptr), 
      fSameEventTripletResolution(nullptr), 
      fSameEventTripletResolutionAll(nullptr), 
      fMixedEventTripletResolution(nullptr), 
      fMixedEventTripletResolutionAll(nullptr), 
      fQ3Vskq12(nullptr), 
      fQ3Vskq12Mixed(nullptr), 
      fQ3Vskq23(nullptr), 
      fQ3Vskq23Mixed(nullptr), 
      fOtherHistos(nullptr), 
      fQ3VsPairmTListSE(nullptr),  
      fQ3VsPairmTListME(nullptr),  
      fQ3VsPairmTList2SE1ME(nullptr),  
      fQ3VsTripletmTListSE(nullptr),  
      fQ3VsTripletmTListME(nullptr),  
      fQ3VsTripletmTList2SE1ME(nullptr),  
      fInvMassTripletSame(nullptr),   
      fInvMassTripletMixed(nullptr),  
      fRunThreeBody(true),
      fRunPlotInvMassTriplet(true),
      fRunPlotQ3Vsq(true),
      fRunPlotPhiTheta(true),
      fRunPlotOtherHistos(true),
      fRunPlotMult(true),
      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fClosePairRejectionPPPorPPL(false),
      fisMC(isMC),
      fRun2Body(false),
      fSame2Mixed1Choice(false),
      fQ3LimitForDeltaPhiDeltaEta(0.4),
      fMixingChoice(0),
      fWhichTripletsToRun(11),  
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fSameEventTripletArrayResolution(nullptr),
      fSameEventTripletArrayResolutionAll(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fPartContainerTESTppp(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fMixedEventTripletArrayResolution(nullptr),
      fMixedEventTripletArrayResolutionAll(nullptr),
      fQ3VskDistributionsArrayq12(nullptr),
      fQ3VskDistributionsArrayq12Mixed(nullptr),
      fQ3VskDistributionsArrayq23(nullptr),
      fQ3VskDistributionsArrayq23Mixed(nullptr),
      fDoubletVsTrippletPPP(nullptr),
      fInvMassSame(nullptr),
      fInvMassMixed(nullptr),
      fTripletInvMassDet(nullptr), 
      fTripletInvMassPDG(nullptr), 
      fTripletInvMassDetMixed(nullptr), 
      fTripletInvMassPDGMixed(nullptr),  
      fTripletInvMassDetAnti(nullptr), 
      fTripletInvMassPDGAnti(nullptr), 
      fTripletInvMassDetMixedAnti(nullptr), 
      fTripletInvMassPDGMixedAnti(nullptr), 
      fQ3VskmTPairSE(nullptr), 
      fQ3VskmTPairME(nullptr), 
      fQ3VskmTPair2SEME1(nullptr), 
      fQ3VskmTTripletSE(nullptr), 
      fQ3VskmTTripletME(nullptr), 
      fQ3VskmTTriplet2SEME1(nullptr), 
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
    fResultsThreeBody->SetName("ResultsThreeBody");
    if(fRunPlotOtherHistos){
      fOtherHistos = new TList();
      fOtherHistos->SetOwner();
      fOtherHistos->SetName("OtherHistos");
      fDoubletVsTrippletPPP =  new TH1F("fDoubletVsTrippletPPP","fDoubletVsTrippletPPP", 3, 0, 3);

      fTripletInvMassDet = new TH2F("fTripletInvMassDet","fTripletInvMassDet",200,0, 2, 1200, 2.8,4);
      fTripletInvMassPDG = new TH2F("fTripletInvMassPDG","fTripletInvMassPDG",200,0, 2, 1200, 2.8,4);
      fTripletInvMassDetMixed = new TH2F("fTripletInvMassDetMixed","fTripletInvMassDetMixed",200,0, 2, 1200, 2.8,4);
      fTripletInvMassPDGMixed = new TH2F("fTripletInvMassPDGMixed","fTripletInvMassPDGMixed",200,0, 2, 1200, 2.8,4);
      
      fTripletInvMassDetAnti = new TH2F("fTripletInvMassDetAnti","fTripletInvMassDetAnti",200,0, 2, 1200, 2.8,4);
      fTripletInvMassPDGAnti = new TH2F("fTripletInvMassPDGAnti","fTripletInvMassPDGAnti",200,0, 2, 1200, 2.8,4);
      fTripletInvMassDetMixedAnti = new TH2F("fTripletInvMassDetMixedAnti","fTripletInvMassDetMixedAnti",200,0, 2, 1200, 2.8,4);
      fTripletInvMassPDGMixedAnti = new TH2F("fTripletInvMassPDGMixedAnti","fTripletInvMassPDGMixedAnti",200,0, 2, 1200, 2.8,4);


      fOtherHistos->Add(fDoubletVsTrippletPPP);
      fOtherHistos->Add(fTripletInvMassDet);
      fOtherHistos->Add(fTripletInvMassPDG);
      fOtherHistos->Add(fTripletInvMassDetMixed);
      fOtherHistos->Add(fTripletInvMassPDGMixed);
      fOtherHistos->Add(fTripletInvMassDetAnti);
      fOtherHistos->Add(fTripletInvMassPDGAnti);
      fOtherHistos->Add(fTripletInvMassDetMixedAnti);
      fOtherHistos->Add(fTripletInvMassPDGMixedAnti);

      fResultsThreeBody->Add(fOtherHistos);
    }


    // Same event
    fSameEvent = new TList();
    fSameEvent->SetOwner();
    fSameEvent->SetName("SameEvent");

    fSameEventTripletArray = new TH1F*[28];
    TString histTitlesSame[28] ={"sameEventDistributionPL","sameEventDistributionPPL","sameEventDistributionAPAPAL",
      "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPLL","sameEventDistributionAPALAL",
      "sameEventDistributionLLL","sameEventDistributionALALAL", "sameEventDistributionppSameLMixed", "sameEventDistributionpLSamepMixed",
      "sameEventDistributionppSamepMixed","sameEventDistributionPP","sameEventDistributionapapSameaLMixed", "sameEventDistributionapaLSameapMixed",
      "sameEventDistributionapapSameapMixed", "sameEventDistributionpLSameLMixed", "sameEventDistributionLLSamepMixed",
      "sameEventDistributionapaLSameaLMixed", "sameEventDistributionaLaLSameapMixed","sameEventDistributionLLSameLMixed", "sameEventDistributionaLaLSameaLMixed",
      "sameEventDistributionppSameLMixedTEST", "sameEventDistributionpLSamepMixedTEST",
      "sameEventDistributionppSamepMixedTEST","sameEventDistributionapapSameaLMixedTEST", "sameEventDistributionapaLSameapMixedTEST",
      "sameEventDistributionapapSameapMixedTEST"};
    fSameEventTripletArray[0] = new TH1F(histTitlesSame[0],histTitlesSame[0],1000,0, 1);
    fSameEvent->Add(fSameEventTripletArray[0]);
    for (int i = 1; i < 12; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }
    fSameEventTripletArray[12] = new TH1F(histTitlesSame[12],histTitlesSame[12],1000,0, 1);
    fSameEvent->Add(fSameEventTripletArray[12]);
    for (int i = 13; i <28 ; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }
    // Same event multiplicity dist
    fSameEventMult = new TList();
    fSameEventMult->SetOwner();
    fSameEventMult->SetName("SameEventMult");

    fSameEventTripletMultArray = new TH2F*[27];
    TString histTitlesSameMult[27] ={"sameEventDistributionMultPL","sameEventDistributionMultPPL","sameEventDistributionMultAPAPAL",
      "sameEventDistributionMultPPP", "sameEventDistributionMultAPAPAP", "sameEventDistributionMultPLL","sameEventDistributionMultAPALAL",
      "sameEventDistributionMultLLL","sameEventDistributionMultALALAL",  "sameEventDistributionMultppSameLMixed", "sameEventDistributionMultpLSamepMixed",
      "sameEventDistributionMultppSamepMixed","sameEventDistributionMultapapSameaLMixed", "sameEventDistributionMultapaLSameapMixed",
      "sameEventDistributionMultapapSameapMixed", "sameEventDistributionMultpLSameLMixed", "sameEventDistributionMultLLSamepMixed",
      "sameEventDistributionMultapaLSameaLMixed", "sameEventDistributionMultaLaLSameapMixed","sameEventDistributionMultLLSameLMixed", "sameEventDistributionMultaLaLSameaLMixed",
      "sameEventDistributionMultppSameLMixedTEST", "sameEventDistributionMultpLSamepMixedTEST",
      "sameEventDistributionMultppSamepMixedTEST","sameEventDistributionMultapapSameaLMixedTEST", "sameEventDistributionMultapaLSameapMixedTEST",
      "sameEventDistributionMultapapSameapMixedTEST"};
    
    fSameEventTripletMultArray[0] = new TH2F(histTitlesSameMult[0],histTitlesSameMult[0],1000,0, 1,26,1,27);
    fSameEventMult->Add(fSameEventTripletMultArray[0]);
    for (int i = 1; i < 27; ++i) {
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

    // Mixed event resolution Matrix 
    fMixedEventTripletResolution = new TList();
    fMixedEventTripletResolution->SetOwner();
    fMixedEventTripletResolution->SetName("MixedEventResolution");

    fMixedEventTripletResolutionAll = new TList();
    fMixedEventTripletResolutionAll->SetOwner();
    fMixedEventTripletResolutionAll->SetName("MixedEventResolutionAll");

    fMixedEventTripletArrayResolution = new TH2F*[4];
    fMixedEventTripletArrayResolutionAll = new TH2F*[4];
    TString histTitlesMixedResolution[4] ={"mixedEventResolutionPPP","mixedEventResolutionAPAPAP","mixedEventResolutionPPL","mixedEventResolutionAPAPAL"};
    TString histTitlesMixedResolutionAll[4] ={"mixedEventResolutionPPPAll","mixedEventResolutionAPAPAPAll","mixedEventResolutionPPLAll","mixedEventResolutionAPAPALAll"};
    for (int i = 0; i < 4; ++i) {
      fMixedEventTripletArrayResolution[i] = new TH2F(histTitlesMixedResolution[i],histTitlesMixedResolution[i], 1000, 0, 1, 1000, 0, 1);
      fMixedEventTripletArrayResolutionAll[i] = new TH2F(histTitlesMixedResolution[i],histTitlesMixedResolution[i], 1000, 0, 1, 1000, 0, 1);
      fMixedEventTripletResolution->Add(fMixedEventTripletArrayResolution[i]);
      fMixedEventTripletResolutionAll->Add(fMixedEventTripletArrayResolutionAll[i]);
     }

    // Same event resolution Matrix 
    fSameEventTripletResolution = new TList();
    fSameEventTripletResolution->SetOwner();
    fSameEventTripletResolution->SetName("SameEventResolution");

    fSameEventTripletResolutionAll = new TList();
    fSameEventTripletResolutionAll->SetOwner();
    fSameEventTripletResolutionAll->SetName("SameEventResolutionAll");

    fSameEventTripletArrayResolution = new TH2F*[4];
    fSameEventTripletArrayResolutionAll = new TH2F*[4];
    TString histTitlesSameResolution[4] ={"sameEventResolutionPPP","sameEventResolutionAPAPAP","sameEventResolutionPPL","sameEventResolutionAPAPAL"};
    TString histTitlesSameResolutionAll[4] ={"sameEventResolutionPPPAll","sameEventResolutionAPAPAPAll","sameEventResolutionPPLAll","sameEventResolutionAPAPALAll"};
    for (int i = 0; i < 4; ++i) {
      fSameEventTripletArrayResolution[i] = new TH2F(histTitlesSameResolution[i],histTitlesSameResolution[i], 1000, 0, 1, 1000, 0, 1);
      fSameEventTripletArrayResolutionAll[i] = new TH2F(histTitlesSameResolution[i],histTitlesSameResolution[i], 1000, 0, 1, 1000, 0, 1);
      fSameEventTripletResolution->Add(fSameEventTripletArrayResolution[i]);
      fSameEventTripletResolutionAll->Add(fSameEventTripletArrayResolutionAll[i]);
     }

    // end mixed event resolution matrix

    fResultsThreeBody->Add(fSameEvent);
    fResultsThreeBody->Add(fMixedEvent);
    fResultsThreeBody->Add(fSameEventMult);
    fResultsThreeBody->Add(fMixedEventMult);

    fResultsThreeBody->Add(fMixedEventTripletResolution);
    fResultsThreeBody->Add(fMixedEventTripletResolutionAll);

    fResultsThreeBody->Add(fSameEventTripletResolution);
    fResultsThreeBody->Add(fSameEventTripletResolutionAll);


    // Q3 vs pair mT SERKSNYTE
    fQ3VsPairmTListSE = new TList();
    fQ3VsPairmTListSE->SetOwner();
    fQ3VsPairmTListSE->SetName("fQ3VsPairmTListSE");

    fQ3VskmTPairSE =  new TH2F*[8];
    TString histTitlesfQ3VskmTPairSE[8] =  {"Q3vsmTPairsSEPPL","Q3vsmTPairsSEAPAPAL",
      "Q3vsmTPairsSEPPP", "Q3vsmTPairsSEAPAPAP", "Q3vsmTPairsSEPLL","Q3vsmTPairsSEAPALAL",
      "Q3vsmTPairsSELLL","Q3vsmTPairsSEALALAL" };
    for(int i=0;i<8;i++){
      fQ3VskmTPairSE[i] = new TH2F(histTitlesfQ3VskmTPairSE[i],histTitlesfQ3VskmTPairSE[i], 1000, 0., 2.,1000,0.,5.);
      fQ3VsPairmTListSE->Add(fQ3VskmTPairSE[i]);
    }
    fResultsThreeBody->Add(fQ3VsPairmTListSE);


    fQ3VsPairmTListME = new TList();
    fQ3VsPairmTListME->SetOwner();
    fQ3VsPairmTListME->SetName("fQ3VsPairmTListME");

    fQ3VskmTPairME =  new TH2F*[8];
    TString histTitlesfQ3VskmTPairME[8] =  {"Q3vsmTPairsMEPPL","Q3vsmTPairsMEAPAPAL",
      "Q3vsmTPairsMEPPP", "Q3vsmTPairsMEAPAPAP", "Q3vsmTPairsMEPLL","Q3vsmTPairsMEAPALAL",
      "Q3vsmTPairsMELLL","Q3vsmTPairsMEALALAL" };
    for(int i=0;i<8;i++){
      fQ3VskmTPairME[i] = new TH2F(histTitlesfQ3VskmTPairME[i],histTitlesfQ3VskmTPairME[i], 1000, 0., 2.,1000,0.,5.);
      fQ3VsPairmTListME->Add(fQ3VskmTPairME[i]);
    }
    fResultsThreeBody->Add(fQ3VsPairmTListME);


    fQ3VsPairmTList2SE1ME = new TList();
    fQ3VsPairmTList2SE1ME->SetOwner();
    fQ3VsPairmTList2SE1ME->SetName("fQ3VsPairmTList2SE1ME");

    fQ3VskmTPair2SEME1 =  new TH2F*[6];
    TString histTitlesfQ3VskmTPair2SEME1[6] =  {"Q3vsmTPairs2SE1MEPPsameLmixed","Q3vsmTPairs2SE1MEAPAPsameALmixed", 
    "Q3vsmTPairs2SE1MEPLsamePmixed","Q3vsmTPairs2SE1MEAPALsameAPmixed", "Q3vsmTPairs2SE1MEPPsamePmixed","Q3vsmTPairs2SE1MEAPAPsameAPmixed"};
    for(int i=0;i<6;i++){
      fQ3VskmTPair2SEME1[i] = new TH2F(histTitlesfQ3VskmTPair2SEME1[i],histTitlesfQ3VskmTPair2SEME1[i], 1000, 0., 2.,1000,0.,5.);
      fQ3VsPairmTList2SE1ME->Add(fQ3VskmTPair2SEME1[i]);
    }
    fResultsThreeBody->Add(fQ3VsPairmTList2SE1ME);

    // Q3 vs triplet mT SERKSNYTE
    fQ3VsTripletmTListSE = new TList();
    fQ3VsTripletmTListSE->SetOwner();
    fQ3VsTripletmTListSE->SetName("fQ3VsTripletmTListSE");

    fQ3VskmTTripletSE =  new TH2F*[8];
    TString histTitlesfQ3VskmTTripletSE[8] =  {"Q3vsmTTripletsSEPPL","Q3vsmTTripletsSEAPAPAL",
      "Q3vsmTTripletsSEPPP", "Q3vsmTTripletsSEAPAPAP", "Q3vsmTTripletsSEPLL","Q3vsmTTripletsSEAPALAL",
      "Q3vsmTTripletsSELLL","Q3vsmTTripletsSEALALAL" };
    for(int i=0;i<8;i++){
      fQ3VskmTTripletSE[i] = new TH2F(histTitlesfQ3VskmTTripletSE[i],histTitlesfQ3VskmTTripletSE[i], 1000, 0., 2.,1000,0.,5.);
      fQ3VsTripletmTListSE->Add(fQ3VskmTTripletSE[i]);
    }
    fResultsThreeBody->Add(fQ3VsTripletmTListSE);


    fQ3VsTripletmTListME = new TList();
    fQ3VsTripletmTListME->SetOwner();
    fQ3VsTripletmTListME->SetName("fQ3VsTripletmTListME");

    fQ3VskmTTripletME =  new TH2F*[8];
    TString histTitlesfQ3VskmTTripletME[8] =  {"Q3vsmTTripletsMEPPL","Q3vsmTTripletsMEAPAPAL",
      "Q3vsmTTripletsMEPPP", "Q3vsmTTripletsMEAPAPAP", "Q3vsmTTripletsMEPLL","Q3vsmTTripletsMEAPALAL",
      "Q3vsmTTripletsMELLL","Q3vsmTTripletsMEALALAL" };
    for(int i=0;i<8;i++){
      fQ3VskmTTripletME[i] = new TH2F(histTitlesfQ3VskmTTripletME[i],histTitlesfQ3VskmTTripletME[i], 1000, 0., 2.,1000,0.,5.);
      fQ3VsTripletmTListME->Add(fQ3VskmTTripletME[i]);
    }
    fResultsThreeBody->Add(fQ3VsTripletmTListME);


    fQ3VsTripletmTList2SE1ME = new TList();
    fQ3VsTripletmTList2SE1ME->SetOwner();
    fQ3VsTripletmTList2SE1ME->SetName("fQ3VsTripletmTList2SE1ME");

    fQ3VskmTTriplet2SEME1 =  new TH2F*[6];
    TString histTitlesfQ3VskmTTriplet2SEME1[6] =  {"Q3vsmTTriplets2SE1MEPPsameLmixed","Q3vsmTTriplets2SE1MEAPAPsameALmixed", 
    "Q3vsmTTriplets2SE1MEPLsamePmixed","Q3vsmTTriplets2SE1MEAPALsameAPmixed", "Q3vsmTTriplets2SE1MEPPsamePmixed","Q3vsmTTriplets2SE1MEAPAPsameAPmixed"};
    for(int i=0;i<6;i++){
      fQ3VskmTTriplet2SEME1[i] = new TH2F(histTitlesfQ3VskmTTriplet2SEME1[i],histTitlesfQ3VskmTTriplet2SEME1[i], 1000, 0., 2.,1000,0.,5.);
      fQ3VsTripletmTList2SE1ME->Add(fQ3VskmTTriplet2SEME1[i]);
    }
    fResultsThreeBody->Add(fQ3VsTripletmTList2SE1ME);



    // Same event phi theta distribution
    fSameEventPhiTheta = new TList();
    fSameEventPhiTheta->SetOwner();
    fSameEventPhiTheta->SetName("SameEventPhiTheta");

    fSameEventTripletPhiThetaArray = new TH2F*[44];
    TString histTitlesSamePhiEta[21] = {"sameEventPhiEtaPPL","sameEventPhiEtaAPAPAL",
      "sameEventPhiEtaPPP", "sameEventPhiEtaAPAPAP", "sameEventPhiEtaPLL","sameEventPhiEtaAPALAL",
      "sameEventPhiEtaLLL","sameEventPhiEtaALALAL", "sameEventDistributionPhiEtappSameLMixed", "sameEventDistributionPhiEtapLSamepMixed",
      "sameEventDistributionPhiEtappSamepMixed","sameEventDistributionPhiEtaapapSameaLMixed", "sameEventDistributionPhiEtaapaLSameapMixed",
      "sameEventDistributionPhiEtaapapSameapMixed", "sameEventDistributionPhiEtapLSameLMixed", "sameEventDistributionPhiEtaLLSamepMixed",
      "sameEventDistributionPhiEtaapaLSameaLMixed", "sameEventDistributionPhiEtaaLaLSameapMixed","sameEventDistributionPhiEtaLLSameLMixed", "sameEventDistributionPhiEtaaLaLSameaLMixed", "TRASH"};
    for(int i=0;i<21;i++){
      fSameEventTripletPhiThetaArray[i] = new TH2F(histTitlesSamePhiEta[i]+"Before",histTitlesSamePhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventTripletPhiThetaArray[21+i] = new TH2F(histTitlesSamePhiEta[i]+"After",histTitlesSamePhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[i]);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[21+i]);
    }
    fSameEventTripletPhiThetaArray[42] = new TH2F("sameEventPhiEtaPPBefore","sameEventPhiEtaPPBefore", 500, -0.15,0.15,500,-0.15,0.15);
    fSameEventTripletPhiThetaArray[43] = new TH2F("sameEventPhiEtaPPAfter","sameEventPhiEtaPPAfter", 500, -0.15,0.15,500,-0.15,0.15);
    fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[42]);
    fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[43]);

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


    // Q3 vs q12 plot for theory
    if(fRunPlotQ3Vsq){
      fQ3Vskq12 = new TList();
      fQ3Vskq12->SetOwner();
      fQ3Vskq12->SetName("Q3Vskq12");

      fQ3VskDistributionsArrayq12 =  new TH2F*[13];
      TString histTitlesfQ3VskDistributions[13] =  {"Q3vskDistributionppSameLMixed", "Q3vskDistributionpLSamepMixed",
        "Q3vskDistributionppSamepMixed","Q3vskDistributionapapSameaLMixed", "Q3vskDistributionapaLSameapMixed",
        "Q3vskDistributionapapSameapMixed", "Q3vskDistributionpLSameLMixed", "Q3vskDistributionLLSamepMixed",
        "Q3vskDistributionapaLSameaLMixed", "Q3vskDistributionaLaLSameapMixed","Q3vskDistributionLLSameLMixed", "Q3vskDistributionaLaLSameaLMixed", "TRASH"};
      for(int i=0;i<13;i++){
        fQ3VskDistributionsArrayq12[i] = new TH2F(histTitlesfQ3VskDistributions[i],histTitlesfQ3VskDistributions[i], 1000, 0., 1.,1000,0.,1.);
        fQ3Vskq12->Add(fQ3VskDistributionsArrayq12[i]);
      }


      // Q3 vs q12 plot for theory
      fQ3Vskq12Mixed = new TList();
      fQ3Vskq12Mixed->SetOwner();
      fQ3Vskq12Mixed->SetName("Q3Vskq12Mixed");

      fQ3VskDistributionsArrayq12Mixed =  new TH2F*[8];
      TString histTitlesfQ3VskDistributionsMixed[8] =  {"mixedEventDistributionPPL","mixedEventDistributionAPAPAL",
        "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPLL","mixedEventDistributionAPALAL",
        "mixedEventDistributionLLL","mixedEventDistributionALALAL" };
      for(int i=0;i<8;i++){
        fQ3VskDistributionsArrayq12Mixed[i] = new TH2F(histTitlesfQ3VskDistributionsMixed[i],histTitlesfQ3VskDistributionsMixed[i], 1000, 0., 1.,1000,0.,1.);
        fQ3Vskq12Mixed->Add(fQ3VskDistributionsArrayq12Mixed[i]);
      }



      // Q3 vs q23 plot for theory
      fQ3Vskq23 = new TList();
      fQ3Vskq23->SetOwner();
      fQ3Vskq23->SetName("Q3Vskq23");

      fQ3VskDistributionsArrayq23 =  new TH2F*[13];
      for(int i=0;i<13;i++){
        fQ3VskDistributionsArrayq23[i] = new TH2F(histTitlesfQ3VskDistributions[i],histTitlesfQ3VskDistributions[i], 1000, 0., 1.,1000,0.,1.);
        fQ3Vskq23->Add(fQ3VskDistributionsArrayq23[i]);
      }

          // Q3 vs q23 plot for theory
      fQ3Vskq23Mixed = new TList();
      fQ3Vskq23Mixed->SetOwner();
      fQ3Vskq23Mixed->SetName("Q3Vskq23Mixed");

      fQ3VskDistributionsArrayq23Mixed =  new TH2F*[8];
      TString histTitlesfQ3VskDistributionsMixed2[8] =  {"mixedEventDistributionPPL","mixedEventDistributionAPAPAL",
        "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPLL","mixedEventDistributionAPALAL",
        "mixedEventDistributionLLL","mixedEventDistributionALALAL" };
      for(int i=0;i<8;i++){
        fQ3VskDistributionsArrayq23Mixed[i] = new TH2F(histTitlesfQ3VskDistributionsMixed2[i],histTitlesfQ3VskDistributionsMixed2[i], 1000, 0., 1.,1000,0.,1.);
        fQ3Vskq23Mixed->Add(fQ3VskDistributionsArrayq23Mixed[i]);
      }


      fResultsThreeBody->Add(fQ3Vskq12);
      fResultsThreeBody->Add(fQ3Vskq12Mixed);
      fResultsThreeBody->Add(fQ3Vskq23);
      fResultsThreeBody->Add(fQ3Vskq23Mixed);
    }

    if (fRunPlotInvMassTriplet){
      fInvMassTripletSame = new TList();
      fInvMassTripletSame->SetOwner();
      fInvMassTripletSame->SetName("InvMassTripletSame");

      fInvMassSame = new TH2F*[8];
      TString histTitlesInvMassSame[8] ={"sameEventDistributionPPL","sameEventDistributionAPAPAL",
        "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPLL","sameEventDistributionAPALAL",
        "sameEventDistributionLLL","sameEventDistributionALALAL"};
      for(int i=0;i<8;i++){
        fInvMassSame[i] = new TH2F(histTitlesInvMassSame[i],histTitlesInvMassSame[i], 33, 0, 0.99, 400, 1.0, 1.2);
        fInvMassTripletSame->Add(fInvMassSame[i]);
      }


      fInvMassTripletMixed = new TList();
      fInvMassTripletMixed->SetOwner();
      fInvMassTripletMixed->SetName("InvMassTripletMixed");

      fInvMassMixed =new TH2F*[8];
      TString histTitlesInvMassMixed[8] ={"mixedEventDistributionPPL","mixedEventDistributionAPAPAL",
        "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPLL","mixedEventDistributionAPALAL",
        "mixedEventDistributionLLL","mixedEventDistributionALALAL"};
      for(int i=0;i<8;i++){
        fInvMassMixed[i] = new TH2F(histTitlesInvMassMixed[i],histTitlesInvMassMixed[i], 33, 0, 0.99, 400, 1.0, 1.2);
        fInvMassTripletMixed->Add(fInvMassMixed[i]);
      }
    
      fResultsThreeBody->Add(fInvMassTripletSame);
      fResultsThreeBody->Add(fInvMassTripletMixed);
    }

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
    // proton lambda, as a test case
    //FillPairDistributionPL(ParticleVector,fSameEventTripletArray[0],bins[1],fSameEventTripletMultArray[0]);
    // pp
    //FillPairDistributionPP(ParticleVector,fSameEventTripletArray[12],  fSameEventTripletPhiThetaArray, *fConfig);
    if(fWhichTripletsToRun ==1||fWhichTripletsToRun ==11){
      // proton proton lambda
      FillTripletDistribution( ParticleVector, 2, 0, 0, fSameEventTripletArray[1],PDGCodes, bins[1],fSameEventTripletMultArray[1], fSameEventTripletPhiThetaArray,0, *fConfig, fSameEventTripletArrayResolution[2],fSameEventTripletArrayResolutionAll[2] , fisMC, fQ3VskmTPairSE[0], fQ3VskmTTripletSE[0]);//, fInvMassSame[0]);
       // antiproton antiproton antilambad
      FillTripletDistribution( ParticleVector, 3, 1, 1, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletPhiThetaArray,1, *fConfig, fSameEventTripletArrayResolution[3],fSameEventTripletArrayResolutionAll[3] , fisMC, fQ3VskmTPairSE[1], fQ3VskmTTripletSE[1]);//, fInvMassSame[1]);
    }
    if(fWhichTripletsToRun ==0||fWhichTripletsToRun ==11){
      // proton proton proton 
      FillTripletDistribution( ParticleVector, 0, 0, 0, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletPhiThetaArray,2, *fConfig, fSameEventTripletArrayResolution[0],fSameEventTripletArrayResolutionAll[0] , fisMC, fQ3VskmTPairSE[2], fQ3VskmTTripletSE[2]);//, fInvMassSame[2]);
      // antiproton antiproton antiproton 
      FillTripletDistribution( ParticleVector, 1, 1, 1, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletPhiThetaArray,3, *fConfig, fSameEventTripletArrayResolution[1],fSameEventTripletArrayResolutionAll[1] , fisMC, fQ3VskmTPairSE[3], fQ3VskmTTripletSE[3]);//, fInvMassSame[3]);
    }
    /*// proton lambda lambda 
    FillTripletDistribution( ParticleVector, 0, 2, 2, fSameEventTripletArray[5],PDGCodes, bins[1],fSameEventTripletMultArray[5], fSameEventTripletPhiThetaArray,4, *fConfig);//, fInvMassSame[4]);
    // antiproton antilambda antilambda 
    FillTripletDistribution( ParticleVector, 1, 3, 3, fSameEventTripletArray[6],PDGCodes, bins[1],fSameEventTripletMultArray[6], fSameEventTripletPhiThetaArray,5, *fConfig);//, fInvMassSame[5]);
    // lambda lambda lambda 
    FillTripletDistribution( ParticleVector, 2, 2, 2, fSameEventTripletArray[7],PDGCodes, bins[1],fSameEventTripletMultArray[7], fSameEventTripletPhiThetaArray,6, *fConfig);//, fInvMassSame[6]);
    // antilambda antilambda antilambda 
    FillTripletDistribution( ParticleVector, 3, 3, 3, fSameEventTripletArray[8],PDGCodes, bins[1],fSameEventTripletMultArray[8], fSameEventTripletPhiThetaArray,7, *fConfig);//, fInvMassSame[7]);
    */

     // Mixed event distribution

    // TAKE CARE OF MULT AND ZVtx!!!!!!!!!!1
    if (!(bins[0] == -99 || bins[1] == -99)) {

      auto itZVtx = fPartContainer.begin()+ bins[0];
      auto itMult = itZVtx->begin() + bins[1];

      /*auto itZVtxTEST = fPartContainerTEST.begin()+ bins[0];
      auto itMultTEST = itZVtxTEST->begin() + bins[1];

      auto itZVtxTESTppL = fPartContainerTESTppL.begin()+ bins[0];
      auto itMultTESTppL = itZVtxTESTppL->begin() + bins[1];

      auto itZVtxTESTppp = fPartContainerTESTppp.begin()+ bins[0];
      auto itMultTESTppp = itZVtxTESTppp->begin() + bins[1];*/

      // Same2Mixed1
      if(fWhichTripletsToRun ==0||fWhichTripletsToRun ==11){
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 0, fSameEventTripletArray[11], PDGCodes, bins[1],fSameEventTripletMultArray[11], fSameEventTripletPhiThetaArray, 10, *fConfig, fQ3VskmTPair2SEME1[4], fQ3VskmTTriplet2SEME1[4]);//, fQ3VskDistributionsArrayq12[2],fQ3VskDistributionsArrayq23[2]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 1, fSameEventTripletArray[15], PDGCodes, bins[1],fSameEventTripletMultArray[14], fSameEventTripletPhiThetaArray, 13, *fConfig, fQ3VskmTPair2SEME1[5], fQ3VskmTTriplet2SEME1[5]);//, fQ3VskDistributionsArrayq12[5],fQ3VskDistributionsArrayq23[5]);
      }
      if(fWhichTripletsToRun ==1||fWhichTripletsToRun ==11){
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 2, fSameEventTripletArray[9], PDGCodes, bins[1],fSameEventTripletMultArray[9], fSameEventTripletPhiThetaArray, 8, *fConfig, fQ3VskmTPair2SEME1[0], fQ3VskmTTriplet2SEME1[0]);//, fQ3VskDistributionsArrayq12[0],fQ3VskDistributionsArrayq23[0]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 0, fSameEventTripletArray[10], PDGCodes, bins[1],fSameEventTripletMultArray[10], fSameEventTripletPhiThetaArray, 9, *fConfig, fQ3VskmTPair2SEME1[2], fQ3VskmTTriplet2SEME1[2]);//, fQ3VskDistributionsArrayq12[1],fQ3VskDistributionsArrayq23[1]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 3, fSameEventTripletArray[13], PDGCodes, bins[1],fSameEventTripletMultArray[12], fSameEventTripletPhiThetaArray, 11, *fConfig, fQ3VskmTPair2SEME1[1], fQ3VskmTTriplet2SEME1[1]);//, fQ3VskDistributionsArrayq12[3],fQ3VskDistributionsArrayq23[3]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 1, fSameEventTripletArray[14], PDGCodes, bins[1],fSameEventTripletMultArray[13], fSameEventTripletPhiThetaArray, 12, *fConfig, fQ3VskmTPair2SEME1[3], fQ3VskmTTriplet2SEME1[3]);//, fQ3VskDistributionsArrayq12[4],fQ3VskDistributionsArrayq23[4]);
      }
      /*FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 2, fSameEventTripletArray[16], PDGCodes, bins[1],fSameEventTripletMultArray[15], fSameEventTripletPhiThetaArray, 14, *fConfig);//, fQ3VskDistributionsArrayq12[6],fQ3VskDistributionsArrayq23[6]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 0, fSameEventTripletArray[17], PDGCodes, bins[1],fSameEventTripletMultArray[16], fSameEventTripletPhiThetaArray, 15, *fConfig);//, fQ3VskDistributionsArrayq12[7],fQ3VskDistributionsArrayq23[7]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 3, fSameEventTripletArray[18], PDGCodes, bins[1],fSameEventTripletMultArray[17], fSameEventTripletPhiThetaArray, 16, *fConfig);//, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 1, fSameEventTripletArray[19], PDGCodes, bins[1],fSameEventTripletMultArray[18], fSameEventTripletPhiThetaArray, 17, *fConfig);//, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 2, fSameEventTripletArray[20], PDGCodes, bins[1],fSameEventTripletMultArray[19], fSameEventTripletPhiThetaArray, 18, *fConfig);//, fQ3VskDistributionsArrayq12[10],fQ3VskDistributionsArrayq23[10]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 3, fSameEventTripletArray[21], PDGCodes, bins[1],fSameEventTripletMultArray[20], fSameEventTripletPhiThetaArray, 19, *fConfig);//, fQ3VskDistributionsArrayq12[11],fQ3VskDistributionsArrayq23[11]);
      */


      //Test what happens with SE2ME1 if we use only events that have triplets
      /*
      auto protonsvector = ParticleVector[0];
      auto antiprotonsvector = ParticleVector[1];
      auto lambdavector = ParticleVector[2];
      auto antilambdavector = ParticleVector[3];

      int protoncount = protonsvector.size();
      int antiprotonscount = antiprotonsvector.size();
      int lambdacount = lambdavector.size();
      int antilambdacount = antilambdavector.size();
      if(fRunPlotOtherHistos){
        if(protoncount==2){
          fDoubletVsTrippletPPP->Fill(1.5);
        }
        if(protoncount>2){
          fDoubletVsTrippletPPP->Fill(2.5);
        }

        if(antiprotonscount==2){
          fDoubletVsTrippletPPP->Fill(1.5);
        }
        if(antiprotonscount>2){
          fDoubletVsTrippletPPP->Fill(2.5);
        }
      }

      if(protoncount>=3){
        FillTripletDistributionSE2ME1(ParticleVector, *itMultTESTppp, 0, 0, 0, fSameEventTripletArray[24], PDGCodes, bins[1],fSameEventTripletMultArray[23], fSameEventTripletPhiThetaArray, 20, *fConfig, fQ3VskDistributionsArrayq12[12],fQ3VskDistributionsArrayq23[12]);
      }
      if(antiprotonscount>=3){
        FillTripletDistributionSE2ME1(ParticleVector, *itMultTESTppp, 1, 1, 1, fSameEventTripletArray[27], PDGCodes, bins[1],fSameEventTripletMultArray[26], fSameEventTripletPhiThetaArray, 20, *fConfig, fQ3VskDistributionsArrayq12[12],fQ3VskDistributionsArrayq23[12]);
      }
      if(protoncount>=2 && lambdacount>=1){
        FillTripletDistributionSE2ME1(ParticleVector, *itMultTESTppL, 0, 0, 2, fSameEventTripletArray[22], PDGCodes, bins[1],fSameEventTripletMultArray[21], fSameEventTripletPhiThetaArray, 20, *fConfig, fQ3VskDistributionsArrayq12[12],fQ3VskDistributionsArrayq23[12]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMultTESTppL, 0, 2, 0, fSameEventTripletArray[23], PDGCodes, bins[1],fSameEventTripletMultArray[22], fSameEventTripletPhiThetaArray, 20, *fConfig, fQ3VskDistributionsArrayq12[12],fQ3VskDistributionsArrayq23[12]);
      }
      if(antiprotonscount>=2 && antilambdacount>=1){
        FillTripletDistributionSE2ME1(ParticleVector, *itMultTESTppL, 1, 1, 3, fSameEventTripletArray[25], PDGCodes, bins[1],fSameEventTripletMultArray[24], fSameEventTripletPhiThetaArray, 20, *fConfig, fQ3VskDistributionsArrayq12[12],fQ3VskDistributionsArrayq23[12]);
        FillTripletDistributionSE2ME1(ParticleVector, *itMultTESTppL, 1, 3, 1, fSameEventTripletArray[26], PDGCodes, bins[1],fSameEventTripletMultArray[25], fSameEventTripletPhiThetaArray, 20, *fConfig, fQ3VskDistributionsArrayq12[12],fQ3VskDistributionsArrayq23[12]);
      }

      //Try to reproduce the p-lambda result from FemtoDream
      FillPairDistributionME(ParticleVector, *itMult, 0, 2, fMixedEventTripletArray[0],PDGCodes,bins[1],fMixedEventTripletMultArray[0],  fMixedEventTripletPhiThetaArray, *fConfig);
      FillPairDistributionME(ParticleVector, *itMult, 0, 0, fMixedEventTripletArray[13],PDGCodes,bins[1],fMixedEventTripletMultArray[0],  fMixedEventTripletPhiThetaArray, *fConfig);
    */
    
      // Normal mixing

      if(fWhichTripletsToRun ==0||fWhichTripletsToRun ==11){
        FillTripletDistributionME(ParticleVector, *itMult, 0, 0, 0, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletPhiThetaArray,2, *fConfig, fMixedEventTripletArrayResolution[0],fMixedEventTripletArrayResolutionAll[0] , fisMC, fQ3VskmTPairME[2], fQ3VskmTTripletME[2]);//, fInvMassMixed[2], fQ3VskDistributionsArrayq12Mixed[2], fQ3VskDistributionsArrayq23Mixed[2]);
        FillTripletDistributionME(ParticleVector, *itMult, 1, 1, 1, fMixedEventTripletArray[4], PDGCodes, bins[1],fMixedEventTripletMultArray[4], fMixedEventTripletPhiThetaArray,3, *fConfig, fMixedEventTripletArrayResolution[1],fMixedEventTripletArrayResolutionAll[1] , fisMC, fQ3VskmTPairME[3], fQ3VskmTTripletME[3]);//, fInvMassMixed[3], fQ3VskDistributionsArrayq12Mixed[3], fQ3VskDistributionsArrayq23Mixed[3]);
      }
      if(fWhichTripletsToRun ==1||fWhichTripletsToRun ==11){
        FillTripletDistributionME(ParticleVector, *itMult, 2, 0, 0, fMixedEventTripletArray[1], PDGCodes, bins[1],fMixedEventTripletMultArray[1], fMixedEventTripletPhiThetaArray,0, *fConfig, fMixedEventTripletArrayResolution[2],fMixedEventTripletArrayResolutionAll[2] , fisMC, fQ3VskmTPairME[0], fQ3VskmTTripletME[0]);//, fInvMassMixed[0], fQ3VskDistributionsArrayq12Mixed[0], fQ3VskDistributionsArrayq23Mixed[0]);
        FillTripletDistributionME(ParticleVector, *itMult, 3, 1, 1, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletPhiThetaArray,1, *fConfig, fMixedEventTripletArrayResolution[3],fMixedEventTripletArrayResolutionAll[3] , fisMC, fQ3VskmTPairME[1], fQ3VskmTTripletME[1]);//, fInvMassMixed[1], fQ3VskDistributionsArrayq12Mixed[1], fQ3VskDistributionsArrayq23Mixed[1]);
      }
      /*FillTripletDistributionME(ParticleVector, *itMult, 0, 2, 2, fMixedEventTripletArray[5], PDGCodes, bins[1],fMixedEventTripletMultArray[5], fMixedEventTripletPhiThetaArray,4, *fConfig, fMixedEventTripletArrayResolution[4],fMixedEventTripletArrayResolutionAll[4] , isMC);//, fInvMassMixed[4], fQ3VskDistributionsArrayq12Mixed[4], fQ3VskDistributionsArrayq23Mixed[4]);
      FillTripletDistributionME(ParticleVector, *itMult, 1, 3, 3, fMixedEventTripletArray[6], PDGCodes, bins[1],fMixedEventTripletMultArray[6], fMixedEventTripletPhiThetaArray,5, *fConfig, fMixedEventTripletArrayResolution[4],fMixedEventTripletArrayResolutionAll[4] , isMC);//, fInvMassMixed[5], fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5]);
      FillTripletDistributionME(ParticleVector, *itMult, 2, 2, 2, fMixedEventTripletArray[7], PDGCodes, bins[1],fMixedEventTripletMultArray[7], fMixedEventTripletPhiThetaArray,6, *fConfig, fMixedEventTripletArrayResolution[4],fMixedEventTripletArrayResolutionAll[4] , isMC);//, fInvMassMixed[6], fQ3VskDistributionsArrayq12Mixed[6], fQ3VskDistributionsArrayq23Mixed[6]);
      FillTripletDistributionME(ParticleVector, *itMult, 3, 3, 3, fMixedEventTripletArray[8], PDGCodes, bins[1],fMixedEventTripletMultArray[8], fMixedEventTripletPhiThetaArray,7, *fConfig, fMixedEventTripletArrayResolution[4],fMixedEventTripletArrayResolutionAll[4] , isMC);//, fInvMassMixed[7], fQ3VskDistributionsArrayq12Mixed[7], fQ3VskDistributionsArrayq23Mixed[7]);
      */
      // Proton Lambda mixing for both proton and lambda used from same event [lambda_same, proton_mixed, proton_mixed]
      //  and [proton_same, lambda_mixed, proton_mixed]
      /*FillTripletDistributionMETEST(ParticleVector, *itMultTEST, 2, 0, 0, fMixedEventTripletArray[9], PDGCodes, bins[1],fMixedEventTripletMultArray[9], fMixedEventTripletPhiThetaArray,8, *fConfig);
      FillTripletDistributionMETEST(ParticleVector, *itMultTEST, 0, 2, 0, fMixedEventTripletArray[9], PDGCodes, bins[1],fMixedEventTripletMultArray[9], fMixedEventTripletPhiThetaArray,8, *fConfig);
      
      // Same for antilambda antiproton
      FillTripletDistributionMETEST(ParticleVector, *itMultTEST, 3, 1, 1, fMixedEventTripletArray[10], PDGCodes, bins[1],fMixedEventTripletMultArray[10], fMixedEventTripletPhiThetaArray,9, *fConfig);
      FillTripletDistributionMETEST(ParticleVector, *itMultTEST, 1, 3, 1, fMixedEventTripletArray[10], PDGCodes, bins[1],fMixedEventTripletMultArray[10], fMixedEventTripletPhiThetaArray,9, *fConfig);
      
      // Test one more mixing, when events taken for mixing are only the events that have tripplets, mixing bothp same and lambda same
      FillTripletDistributionMETEST(ParticleVector, *itMultTESTppL, 2, 0, 0, fMixedEventTripletArray[11], PDGCodes, bins[1],fMixedEventTripletMultArray[11], fMixedEventTripletPhiThetaArray,10, *fConfig);
      FillTripletDistributionMETEST(ParticleVector, *itMultTESTppL, 0, 2, 0, fMixedEventTripletArray[11], PDGCodes, bins[1],fMixedEventTripletMultArray[11], fMixedEventTripletPhiThetaArray,10, *fConfig);
      
      // Same for antilambda antiproton
      FillTripletDistributionMETEST(ParticleVector, *itMultTESTppL, 3, 1, 1, fMixedEventTripletArray[12], PDGCodes, bins[1],fMixedEventTripletMultArray[12], fMixedEventTripletPhiThetaArray,11, *fConfig);
      FillTripletDistributionMETEST(ParticleVector, *itMultTESTppL, 1, 3, 1, fMixedEventTripletArray[12], PDGCodes, bins[1],fMixedEventTripletMultArray[12], fMixedEventTripletPhiThetaArray,11, *fConfig);
      */
      


      // Update the particle container with current event
      if(fMixingChoice==0){
        SetMixedEvent(ParticleVector, &(*itMult));
      }else if(fMixingChoice==1){
        SetMixedEventOnlyPLambdaTEST(ParticleVector, &(*itMult));
      }else if(fMixingChoice==2){
        SetMixedEventOnlyPPLambdaTEST(ParticleVector, &(*itMult));
      }else if(fMixingChoice==3){
        SetMixedEventOnlyPPPTEST(ParticleVector, &(*itMult));
      }
      //SetMixedEventOnlyPLambdaTEST(ParticleVector, &(*itMultTEST));
      //SetMixedEventOnlyPPLambdaTEST( ParticleVector, &(*itMultTESTppL));
      //SetMixedEventOnlyPPPTEST( ParticleVector, &(*itMultTESTppp));

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
  //PostData(6, fResults);
  //PostData(7, fResultsQA);
  //PostData(8, fResultsSample);
  //PostData(9, fResultsSampleQA);
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
  // PR on 15-06-2020 -> don't use the reduced vector - no division by 2
  TLorentzVector trackSum = PartOne + PartTwo;
  TLorentzVector trackDifference = PartOne - PartTwo;
  float scaling = trackDifference*trackSum/(trackSum*trackSum);
  TLorentzVector qPartOnePartTwo = trackDifference -  scaling * trackSum;

  return qPartOnePartTwo;
}

void AliAnalysisTaskThreeBodyFemto::FillMTDistributionsTriplet(TLorentzVector &PartOne, TLorentzVector &PartTwo, TLorentzVector &PartThree, float m1, float m2, float m3, TH2F* pair, TH2F* triplet, float Q3){
  // first pair
  /*auto ptAverage1 = (pt1+pt2)/2.;
  auto mAverage1 = (m1+m2)/2.;
  auto mTAverage1 = sqrt(pow(ptAverage1,2)+pow(mAverage1,2));
  // second pair
  auto ptAverage2 = (pt2+pt3)/2.;
  auto mAverage2 = (m2+m3)/2.;
  auto mTAverage2 = sqrt(pow(ptAverage2,2)+pow(mAverage2,2));
  // third pair
  auto ptAverage3 = (pt1+pt3)/2.;
  auto mAverage3 = (m1+m3)/2.;
  auto mTAverage3 = sqrt(pow(ptAverage3,2)+pow(mAverage3,2));*/


  TLorentzVector trackSumPair12;
  trackSumPair12 = PartOne + PartTwo;
  auto ptAverage1 = 0.5 * trackSumPair12.Pt();
  auto mAverage1 = (m1+m2)/2.;
  auto mTAverage1 = sqrt(pow(ptAverage1,2)+pow(mAverage1,2));


  
  TLorentzVector trackSumPair23;
  trackSumPair23 = PartThree + PartTwo;
  auto ptAverage2 = 0.5 * trackSumPair23.Pt();
  auto mAverage2 =  (m2+m3)/2.;
  auto mTAverage2 = sqrt(pow(ptAverage2,2)+pow(mAverage2,2));


  TLorentzVector trackSumPair13;
  trackSumPair13 = PartThree + PartOne;
  auto ptAverage3 = 0.5 * trackSumPair13.Pt();
  auto mAverage3 = (m1+m3)/2.;
  auto mTAverage3 = sqrt(pow(ptAverage3,2)+pow(mAverage3,2));



  // average of the average mT over pairs
  auto ultimateAverage = (mTAverage1 + mTAverage2 + mTAverage3)/3.;

  // Fill histos
  pair->Fill(Q3,mTAverage1);
  pair->Fill(Q3,mTAverage2);
  pair->Fill(Q3,mTAverage3);
  triplet->Fill(Q3, ultimateAverage);

}

void AliAnalysisTaskThreeBodyFemto::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Res, TH2F* ResAll , bool isMC, TH2F* pairMT, TH2F* tripletMT){
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
        TLorentzVector q12 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part3_LorVec,part1_LorVec);
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
            Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
            Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
            Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
          }
          if(!fClosePairRejectionForAll){
            if(fClosePairRejectionPPPorPPL){
              if(DoThisPair12==11){ 
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair23==11){ 
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair31==11){ 
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              }
            }

            if(!fClosePairRejectionPPPorPPL){
              if(DoThisPair12==21||DoThisPair12==12){ 
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair23==21||DoThisPair23==12){ 
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair31==21||DoThisPair31==12){ 
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              }
            }


          }

        }
        


        if(!Pair12||!Pair23||!Pair31) {continue;}

        FillMTDistributionsTriplet(part1_LorVec, part2_LorVec, part3_LorVec, massparticle1, massparticle2, massparticle3,  pairMT, tripletMT, Q3);
        hist->Fill(Q3); 
        hist2d->Fill(Q3,mult+1); 

        if(isMC){
          MomentumResolution(ResAll, Res, *iPart1, *itPDGPar1, massparticle1,*iPart2, *itPDGPar2, massparticle2, *iPart3,*itPDGPar3, massparticle3, Q3) ;
        }
        /*if(firstSpecies == 2 || firstSpecies == 3) InvMassSame->Fill(Q3, iPart1->GetInvMass());
        if(secondSpecies == 2 || secondSpecies == 3) InvMassSame->Fill(Q3, iPart2->GetInvMass());
        if(thirdSpecies == 2 || thirdSpecies == 3) InvMassSame->Fill(Q3, iPart3->GetInvMass());*/
        
      }
    }
  }
}

void AliAnalysisTaskThreeBodyFemto::FillTripletDistributionPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassSame, TH2F* InvMassDET,TH2F* InvMassPDG){
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
        TLorentzVector q12 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part3_LorVec,part1_LorVec);
        // The particles in current methodology are put in bins of:
        //                 Q3=sqrt(q12^2+q23^2+q31^2)
        float Q32 = q12*q12+q23*q23+q31*q31;
        // From 3 pion paper, the q must be multiplied by -1 before taking quare root
        float Q3 = sqrt(-Q32); // the minus from pion paper
        hist->Fill(Q3); 
        hist2d->Fill(Q3,mult+1); 

        if(firstSpecies == 2 || firstSpecies == 3) InvMassSame->Fill(Q3, iPart1->GetInvMass());
        if(secondSpecies == 2 || secondSpecies == 3) InvMassSame->Fill(Q3, iPart2->GetInvMass());
        if(thirdSpecies == 2 || thirdSpecies == 3) InvMassSame->Fill(Q3, iPart3->GetInvMass());
        if(fRunPlotOtherHistos){
          FillPairInvMass(*iPart1,*iPart2,*iPart3, InvMassDET, Q3);
          FillPDGPairInvMass(*iPart1, massparticle1, *iPart2,massparticle2, *iPart3, massparticle3, InvMassPDG, Q3);
        }

        
      }
    }
  }
}
void AliAnalysisTaskThreeBodyFemto::FillPairDistributionPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPL, int mult, TH2F* hist2d){
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
void AliAnalysisTaskThreeBodyFemto::FillPairDistributionPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPP, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config){
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
      if(!DeltaEtaDeltaPhi(*iPart1,*iPart2,true,11, fEventTripletPhiThetaArray[42],fEventTripletPhiThetaArray[43], Config)){continue;}

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

void AliAnalysisTaskThreeBodyFemto::SetMixedEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  for(unsigned int iSpecies = 0; iSpecies<ParticleVector.size(); iSpecies++){
    if ((ParticleVector.begin()+iSpecies)->size() > 0) {
      (fPartContainer->begin()+iSpecies)->SetEvent(*(ParticleVector.begin()+iSpecies));
    }
  }
}

void AliAnalysisTaskThreeBodyFemto::SetMixedEventOnlyPLambdaTEST(
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

void AliAnalysisTaskThreeBodyFemto::SetMixedEventOnlyPPLambdaTEST(
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

void AliAnalysisTaskThreeBodyFemto::SetMixedEventOnlyPPPTEST(
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


void AliAnalysisTaskThreeBodyFemto::FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Res, TH2F* ResAll, bool isMC, TH2F* pairMT, TH2F* tripletMT ){//, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed){
  // Description of function given in AliAnalysisTaskThreeBodyFemto::FillTripletDistribution
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
            TLorentzVector q12 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part3_LorVec,part1_LorVec);
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
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
              }
              if(!fClosePairRejectionForAll){
                if(fClosePairRejectionPPPorPPL){
                  if(DoThisPair12==11){ 
                    Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair23==11){ 
                    Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair31==11){ 
                    Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                  }
                }
                if(!fClosePairRejectionPPPorPPL){
                  if(DoThisPair12==21||DoThisPair12==12){ 
                    Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair23==21||DoThisPair23==12){ 
                    Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair31==21||DoThisPair31==12){ 
                    Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[20+phiEtaHistNo],Config, Q3);
                  }
                }
              }

            }
        
            if(!Pair12||!Pair23||!Pair31) {continue;}
           

            // Now we have the three particles, lets create their Lorentz vectors  

            FillMTDistributionsTriplet(part1_LorVec, part2_LorVec, part3_LorVec, massParticleSE, massParticleME1, massParticleME2,  pairMT, tripletMT, Q3);
       

            hist->Fill(Q3); 
            hist2d->Fill(Q3,mult+1); 

            if(isMC){
              MomentumResolution(ResAll, Res, *iPart1, *itPDGParSE, massParticleSE,*iPart2, *itPDGParME1, massParticleME1, *iPart3,*itPDGParME2, massParticleME2, Q3) ;
            }
            /*
            if(fRunPlotQ3Vsq){
              float qSame12= sqrt(-q12*q12);
              float qSame23= sqrt(-q23*q23);
              Q3VskDistribution12Mixed->Fill(Q3, qSame12);
              Q3VskDistribution23Mixed->Fill(Q3, qSame23);
            }
            if (fRunPlotInvMassTriplet){
              if(speciesSE == 2 || speciesSE == 3) InvMassMixed->Fill(Q3, iPart1->GetInvMass());
              if(speciesME1 == 2 || speciesME1 == 3) InvMassMixed->Fill(Q3, iPart2->GetInvMass());
              if(speciesME2 == 2 || speciesME2 == 3) InvMassMixed->Fill(Q3, iPart3->GetInvMass());
            }
            */
          }
        }  
      }
    }
  }
}


void AliAnalysisTaskThreeBodyFemto::FillTripletDistributionMEPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F* InvMassDET,TH2F* InvMassPDG){
  // Description of function given in AliAnalysisTaskThreeBodyFemto::FillTripletDistribution
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
            TLorentzVector q12 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part3_LorVec,part1_LorVec);
            // The particles in current methodology are put in bins of:
            //                 Q3=sqrt(q12^2+q23^2+q31^2)
            float Q32 = q12*q12+q23*q23+q31*q31;
            // From 3 pion paper, the q must be multiplied by -1 before taking quare root
            float Q3 = sqrt(-Q32); // the minus from pion paper
            hist->Fill(Q3); 
            

            hist2d->Fill(Q3,mult+1); 

            if(fRunPlotQ3Vsq){
              float qSame12= sqrt(-q12*q12);
              Q3VskDistribution12Mixed->Fill(Q3, qSame12);
            }
            if (fRunPlotInvMassTriplet){
              if(speciesSE == 2 || speciesSE == 3) InvMassMixed->Fill(Q3, iPart1->GetInvMass());
              if(speciesME1 == 2 || speciesME1 == 3) InvMassMixed->Fill(Q3, iPart2->GetInvMass());
              if(speciesME2 == 2 || speciesME2 == 3) InvMassMixed->Fill(Q3, iPart3->GetInvMass());
            }
            if(fRunPlotOtherHistos){
              FillPairInvMass(*iPart1,*iPart2,*iPart3, InvMassDET, Q3);
              FillPDGPairInvMass(*iPart1, massParticleSE, *iPart2,massParticleME1, *iPart3, massParticleME2, InvMassPDG, Q3);
            }
          }
        }  
      }
    }
  }
}



void AliAnalysisTaskThreeBodyFemto::FillTripletDistributionMETEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // Description of function given in AliAnalysisTaskThreeBodyFemto::FillTripletDistribution
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
            TLorentzVector q12 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part3_LorVec,part1_LorVec);
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

void AliAnalysisTaskThreeBodyFemto::FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray,  AliFemtoDreamCollConfig Config){
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


void AliAnalysisTaskThreeBodyFemto::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* pairMT, TH2F* tripletMT){//, TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23){
  // Description of function given in AliAnalysisTaskThreeBodyFemto::FillTripletDistribution
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
          TLorentzVector q12 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part1_LorVec,part2_LorVec);
          TLorentzVector q23 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part2_LorVec,part3_LorVec);
          TLorentzVector q31 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part3_LorVec,part1_LorVec);
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
              Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
              Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
            }
            if(!fClosePairRejectionForAll){
              if(fClosePairRejectionPPPorPPL){
                if(DoThisPair12==11){ 
                  Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair23==11){ 
                  Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair31==11){ 
                  Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
                }
              }
              if(!fClosePairRejectionPPPorPPL){
                if(DoThisPair12==21||DoThisPair12==12){ 
                  Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair23==21||DoThisPair23==12){ 
                  Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair31==21||DoThisPair31==12){ 
                  Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[21+phiEtaHistNo],Config, Q3);
                }
              }
            }

          }

          if(!Pair12||!Pair23||!Pair31) {continue;}

          FillMTDistributionsTriplet(part1_LorVec, part2_LorVec, part3_LorVec, massParticleSE1, massParticleSE2, massParticleME,  pairMT, tripletMT, Q3);
       

          hist->Fill(Q3); 
          hist2d->Fill(Q3,mult+1); 

          /*if(fRunPlotQ3Vsq){
            float qSame12= sqrt(-q12*q12);
            float qSame23= sqrt(-q23*q23);
            Q3VskDistribution12->Fill(Q3, qSame12);
            Q3VskDistribution23->Fill(Q3, qSame23);
          }*/
        }
      }  
    }
  }
}

bool AliAnalysisTaskThreeBodyFemto::DeltaEtaDeltaPhi(
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



bool AliAnalysisTaskThreeBodyFemto::DeltaEtaDeltaPhi(
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist,
                                                   AliFemtoDreamCollConfig Config, double Q3) {
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


void AliAnalysisTaskThreeBodyFemto::FillPairInvMass( AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2,
                        AliFemtoDreamBasePart &part3, 
                        TH2F* hist, float Q3) {
  TVector3 momPart1 = part1.GetMomentum();
  TVector3 momPart2 = part2.GetMomentum();
  TVector3 momPart3 = part3.GetMomentum();
  TLorentzVector track1, track2, track3;
  track1.SetXYZM(momPart1.Px(), momPart1.Py(), momPart1.Pz(),
                   part1.GetInvMass());
  track2.SetXYZM(momPart2.Px(), momPart2.Py(), momPart2.Pz(),
                   part2.GetInvMass());
  track3.SetXYZM(momPart3.Px(), momPart3.Py(), momPart3.Pz(),
                   part3.GetInvMass());
  TLorentzVector trackSum = track1 + track2 + track3;

  hist->Fill(Q3, trackSum.M());
  
}

void AliAnalysisTaskThreeBodyFemto::FillPDGPairInvMass( AliFemtoDreamBasePart &part1, float massPart1,
                        AliFemtoDreamBasePart &part2, float massPart2,
                        AliFemtoDreamBasePart &part3, float massPart3, 
                        TH2F* hist, float Q3) {

  TVector3 momPart1 = part1.GetMomentum();
  TVector3 momPart2 = part2.GetMomentum();
  TVector3 momPart3 = part3.GetMomentum();
  TLorentzVector track1, track2, track3;
  track1.SetXYZM(momPart1.Px(), momPart1.Py(), momPart1.Pz(),
                   massPart1);
  track2.SetXYZM(momPart2.Px(), momPart2.Py(), momPart2.Pz(),
                   massPart2);
  track3.SetXYZM(momPart3.Px(), momPart3.Py(), momPart3.Pz(),
                   massPart3);
  TLorentzVector trackSum = track1 + track2 + track3;

  hist->Fill(Q3, trackSum.M());
    
}


void AliAnalysisTaskThreeBodyFemto::MomentumResolution( TH2F* histAll, TH2F* hist,
    AliFemtoDreamBasePart& part1, int PDGPart1, float mass1,
    AliFemtoDreamBasePart& part2, int PDGPart2, float mass2,
    AliFemtoDreamBasePart& part3, int PDGPart3, float mass3,
    float Q3Reconstructed) {

  TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;

  float mom1 = sqrt(part1.GetMCMomentum().X()*part1.GetMCMomentum().X()+part1.GetMCMomentum().Y()*part1.GetMCMomentum().Y()+part1.GetMCMomentum().Z()*part1.GetMCMomentum().Z());
  float mom2 = sqrt(part2.GetMCMomentum().X()*part2.GetMCMomentum().X()+part2.GetMCMomentum().Y()*part2.GetMCMomentum().Y()+part2.GetMCMomentum().Z()*part2.GetMCMomentum().Z());
  float mom3 = sqrt(part3.GetMCMomentum().X()*part3.GetMCMomentum().X()+part3.GetMCMomentum().Y()*part3.GetMCMomentum().Y()+part3.GetMCMomentum().Z()*part3.GetMCMomentum().Z());
  
  part1_LorVec.SetPxPyPzE(part1.GetMCMomentum().X(), part1.GetMCMomentum().Y(), 
  part1.GetMCMomentum().Z(), sqrt(pow(mom1,2)+pow(mass1,2)));
  part2_LorVec.SetPxPyPzE(part2.GetMCMomentum().X(), part2.GetMCMomentum().Y(), 
  part2.GetMCMomentum().Z(), sqrt(pow(mom2,2)+pow(mass2,2)));
  part3_LorVec.SetPxPyPzE(part3.GetMCMomentum().X(), part3.GetMCMomentum().Y(), 
  part3.GetMCMomentum().Z(), sqrt(pow(mom3,2)+pow(mass3,2)));

  TLorentzVector q12 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part1_LorVec,part2_LorVec);
  TLorentzVector q23 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part2_LorVec,part3_LorVec);
  TLorentzVector q31 = AliAnalysisTaskThreeBodyFemto::RelativePairMomentum(part3_LorVec,part1_LorVec);
  float Q32 = q12*q12+q23*q23+q31*q31;
  float Q3Real = sqrt(-Q32); 
          

  // check if particles are identified correctly

  histAll->Fill(Q3Real, Q3Reconstructed);
  if ((TMath::Abs(PDGPart1) == TMath::Abs(part1.GetMCPDGCode()))
      && ((TMath::Abs(PDGPart2) == TMath::Abs(part2.GetMCPDGCode())))
      && ((TMath::Abs(PDGPart3) == TMath::Abs(part3.GetMCPDGCode())))) {
      hist->Fill(Q3Real, Q3Reconstructed);
  }
}