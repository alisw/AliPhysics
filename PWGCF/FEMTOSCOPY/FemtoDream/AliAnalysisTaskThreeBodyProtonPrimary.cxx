/*
 * AliAnalysisTaskThreeBodyProtonPrimary.cxx
 *
 *  Created on: May 13, 2019
 *      Authors: Raffaele Del Grande, Marcel Lesch
 *      Based on AliAnalysisTaskThreeBodyFemtoAOD.cxx from Laura Serksnyte
 */
#include "AliAnalysisTaskThreeBodyProtonPrimary.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliNanoAODTrack.h"
#include "Riostream.h"

ClassImp(AliAnalysisTaskThreeBodyProtonPrimary)
AliAnalysisTaskThreeBodyProtonPrimary::AliAnalysisTaskThreeBodyProtonPrimary()
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
      fPrimaryTrack(nullptr), 
      fPrimary(nullptr), 
      fPrimaryList(nullptr), 
      fPrimaryMCList(nullptr),
      fAntiPrimary(nullptr),
      fAntiPrimaryList(nullptr),
      fAntiPrimaryMCList(nullptr),
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
      fQ3Vskq12(nullptr),
      fQ3Vskq12Mixed(nullptr),
      fQ3Vskq23(nullptr),
      fQ3Vskq23Mixed(nullptr),
      fOtherHistos(nullptr),
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
      fQ3LimitForDeltaPhiDeltaEta(0.4),
      fCleanWithLambdas(false),
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fSameEventTripletArray_TwoBody(nullptr),
      fSameEventTripletMultArray_TwoBody(nullptr),
      fSameEventTripletPhiThetaArray_TwoBody(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fPartContainerTESTppp(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fMixedEventTripletArray_TwoBody(nullptr),
      fMixedEventTripletMultArray_TwoBody(nullptr),
      fMixedEventTripletPhiThetaArray_TwoBody(nullptr),
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
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskThreeBodyProtonPrimary::AliAnalysisTaskThreeBodyProtonPrimary(const char* name, bool isMC)
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
      fPrimaryTrack(nullptr),
      fPrimary(nullptr),
      fPrimaryList(nullptr),
      fPrimaryMCList(nullptr),
      fAntiPrimary(nullptr),
      fAntiPrimaryList(nullptr),
      fAntiPrimaryMCList(nullptr),
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
      fQ3Vskq12(nullptr),
      fQ3Vskq12Mixed(nullptr),
      fQ3Vskq23(nullptr),
      fQ3Vskq23Mixed(nullptr),
      fOtherHistos(nullptr),
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
      fQ3LimitForDeltaPhiDeltaEta(0.4),
      fCleanWithLambdas(false),
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPhiThetaArray(nullptr),
      fSameEventTripletArray_TwoBody(nullptr),
      fSameEventTripletMultArray_TwoBody(nullptr),
      fSameEventTripletPhiThetaArray_TwoBody(nullptr),
      fPartContainer(0),
      fPartContainerTEST(0),
      fPartContainerTESTppL(0),
      fPartContainerTESTppp(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPhiThetaArray(nullptr),
      fMixedEventTripletArray_TwoBody(nullptr),
      fMixedEventTripletMultArray_TwoBody(nullptr),
      fMixedEventTripletPhiThetaArray_TwoBody(nullptr),
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
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
        DefineOutput(1, TList::Class());  //Output for the Event Cuts
        DefineOutput(2, TList::Class());  //Output for the Proton Cuts
        DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
        DefineOutput(4, TList::Class());  //Output for the Primary Cuts 
        DefineOutput(5, TList::Class());  //Output for the AntiPrimary Cuts 
        DefineOutput(6, TList::Class());  //Output for the Lambda Cuts 
        DefineOutput(7, TList::Class());  //Output for the AntiLambda Cuts 
        DefineOutput(8, TList::Class());  //Output for the Results
        DefineOutput(9, TList::Class());  //Output for the Results QA
        DefineOutput(10, TList::Class());  //Output for the Results
        DefineOutput(11, TList::Class());  //Output for the Results QA
        DefineOutput(12, TList::Class());  //Output for the Results Three body
        if (isMC) {
          DefineOutput(13, TList::Class());  //Output for the Track MC
          DefineOutput(14, TList::Class());  //Output for the Anti Track MC
          DefineOutput(15, TList::Class());  //Output for the Primary MC
          DefineOutput(16, TList::Class());  //Output for the Anti Primary MC
          DefineOutput(17, TList::Class());  //Output for the V0 MC
          DefineOutput(18, TList::Class());  //Output for the Anti V0 MC
        }
      }


//==================================================================================================================================================

AliAnalysisTaskThreeBodyProtonPrimary::~AliAnalysisTaskThreeBodyProtonPrimary() {
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
  if (fPrimaryTrack) { 
    delete fPrimaryTrack; 
  }
  if (fPrimary) { 
    delete fPrimary;
  }
  if (fAntiPrimary) { 
    delete fAntiPrimary; 
  }
  if (fv0) { 
    delete fv0; 
  }
  if(fCleanWithLambdas){
    if (fLambda) {
      delete fLambda; 
    }
    if (fAntiLambda) { 
      delete fAntiLambda; 
    }
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
} //AliAnalysisTaskThreeBodyProtonPrimary::~AliAnalysisTaskThreeBodyProtonPrimary()

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::UserCreateOutputObjects() {
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
  if (!fPrimary) {   
    AliError("No Primary cuts \n");
  } else {
    fPrimary->Init(); 
  }
  if (!fAntiPrimary) { 
    AliError("No AntiPrimary cuts \n");
  } else {
    fAntiPrimary->Init(); 
  }
  if(fCleanWithLambdas){
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

  fPrimaryTrack = new AliFemtoDreamTrack();
  fPrimaryTrack->SetUseMCInfo(
      fPrimary->GetIsMonteCarlo() || fAntiPrimary->GetIsMonteCarlo()); 

  fv0 = new AliFemtoDreamv0();
  if(fCleanWithLambdas){
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
  fPrimaryList = fPrimary->GetQAHists(); 
  fAntiPrimaryList = fAntiPrimary->GetQAHists(); 
  if(fCleanWithLambdas){
     fLambdaList = fLambda->GetQAHists(); 
     fAntiLambdaList = fAntiLambda->GetQAHists(); 
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
    } //if(fRunPlotOtherHistos)

    // Same event
    fSameEvent = new TList();
    fSameEvent->SetOwner();
    fSameEvent->SetName("SameEvent");

    fSameEventTripletArray = new TH1F*[28];
    TString histTitlesSame[28] = {"sameEventDistributionPPPrim","sameEventDistributionAPAPAPrim", "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPPAPrim","sameEventDistributionAPAPPrim", "sameEventDistributionPPrimPrim", "sameEventDistributionAPAPrimAPrim", "sameEventDistributionPAPrimAPrim", "sameEventDistributionAPPrimPrim","sameEventDistributionPPSamePrimMixed", "sameEventDistributionPPrimSamePMixed", "sameEventDistributionAPAPSameAPrimMixed", "sameEventDistributionAPAPrimSameAPMixed", "sameEventDistributionPPSamePMixed", "sameEventDistributionAPAPSameAPMixed", "sameEventDistributionPPSameAPrimMixed", "sameEventDistributionPAPrimSamePMixed", "sameEventDistributionAPAPSamePrimMixed", "sameEventDistributionAPPrimSameAPMixed", "sameEventDistributionPPrimSamePimMixed", "sameEventDistributionPrimPrimSamePMixed", "sameEventDistributionAPAPrimSameAPMixed", "sameEventDistributionAPrimAPrimSameAPMixed", "sameEventDistributionPAPrimSameAPrimMixed", "sameEventDistributionAPrimAPrimSamePMixed", "sameEventDistributionAPPrimSamePrimMixed", "sameEventDistributionPrimPrimSameAPMixed"}; 

    for (int i = 0; i < 10; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }
    for (int i = 10; i <28 ; ++i) {
      fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray[i]);
     }

    fSameEventTripletArray_TwoBody = new TH1F*[7];
    TString histTitlesSame_TwoBody[7] = {"sameEventDistributionPP","sameEventDistributionAPAP", "sameEventDistributionPPrim", "sameEventDistributionAPAPrim", "sameEventDistributionPAPrim", "sameEventDistributionAPPrim", "sameEventDistributionPrimAPrim"}; 

     for (int i = 0; i < 7; ++i) {
      fSameEventTripletArray_TwoBody[i] =  new TH1F(histTitlesSame_TwoBody[i],histTitlesSame_TwoBody[i], 8000, 0, 8);
      fSameEvent->Add(fSameEventTripletArray_TwoBody[i]);
     } 

    // Same event multiplicity dist
    fSameEventMult = new TList();
    fSameEventMult->SetOwner();
    fSameEventMult->SetName("SameEventMult");


    fSameEventTripletMultArray = new TH2F*[28];
    TString histTitlesSameMult[28] = {"sameEventDistributionMultPPPrim","sameEventDistributionMultAPAPAPrim","sameEventDistributionMultPPP", "sameEventDistributionMultAPAPAP", "sameEventDistributionMultPPAPrim","sameEventDistributionMultAPAPPrim", "sameEventDistributionMultPPrimPrim", "sameEventDistributionMultAPAPrimAPrim", "sameEventDistributionMultPAPrimAPrim", "sameEventDistributionMultAPPrimPrim", "sameEventDistributionMultPPSamePrimMixed", "sameEventDistributionMultPPrimSamePMixed", "sameEventDistributionMultAPAPSameAPrimMixed", "sameEventDistributionMultAPAPrimSameAPMixed", "sameEventDistributionMultPPSamePMixed", "sameEventDistributionMultAPAPSameAPMixed", "sameEventDistributionMultPPSameAPrimMixed", "sameEventDistributionMultPAPrimSamePMixed", "sameEventDistributionMultAPAPSamePrimMixed", "sameEventDistributionMultAPPrimSameAPMixed", "sameEventDistributionMultPPrimSamePimMixed", "sameEventDistributionMultPrimPrimSamePMixed", "sameEventDistributionMultAPAPrimSameAPMixed", "sameEventDistributionMultAPrimAPrimSameAPMixed", "sameEventDistributionMultPAPrimSameAPrimMixed", "sameEventDistributionMultAPrimAPrimSamePMixed", "sameEventDistributionMultAPPrimSamePrimMixed", "sameEventDistributionMultPrimPrimSameAPMixed"}; 

    for (int i = 0; i < 28; ++i) {
      fSameEventTripletMultArray[i] =  new TH2F(histTitlesSameMult[i],histTitlesSameMult[i], 8000, 0, 8,26,1,27); 
      fSameEventMult->Add(fSameEventTripletMultArray[i]);
     }

    fSameEventTripletMultArray_TwoBody = new TH2F*[7];
    TString histTitlesSameMult_TwoBody[7] = {"sameEventDistributionMultPP","sameEventDistributionMultAPAP", "sameEventDistributionMultPPrim", "sameEventDistributionMultAPAPrim", "sameEventDistributionMultPAPrim", "sameEventDistributionMultAPPrim", "sameEventDistributionMultPrimAPrim"}; 

     for (int i = 0; i < 7; ++i) {
      fSameEventTripletMultArray_TwoBody[i] =  new TH2F(histTitlesSameMult_TwoBody[i],histTitlesSameMult_TwoBody[i],8000, 0, 8,26,1,27); 
      fSameEventMult->Add(fSameEventTripletMultArray_TwoBody[i]);
     } 
   
    // Mixed event -------------------------------------------------------------------------------
    fMixedEvent = new TList();
    fMixedEvent->SetOwner();
    fMixedEvent->SetName("MixedEvent");

    fMixedEventTripletArray = new TH1F*[10];
    TString histTitlesMixed[10] = {"mixedEventDistributionPPPrim","mixedEventDistributionAPAPAPrim","mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPPAPrim","mixedEventDistributionAPAPPrim", "mixedEventDistributionPPrimPrim", "mixedEventDistributionAPAPrimAPrim", "mixedEventDistributionPAPrimAPrim", "mixedEventDistributionAPPrimPrim"}; 

    for (int i = 0; i < 10; ++i) {
      fMixedEventTripletArray[i] = new TH1F(histTitlesMixed[i],histTitlesMixed[i], 8000, 0, 8);
      fMixedEvent->Add(fMixedEventTripletArray[i]);
     }

    fMixedEventTripletArray_TwoBody = new TH1F*[7];
    TString histTitlesMixed_TwoBody[7] = {"mixedEventDistributionPP","mixedEventDistributionAPAP", "mixedEventDistributionPPrim", "mixedEventDistributionAPAPrim", "mixedEventDistributionPAPrim", "mixedEventDistributionAPPrim", "mixedEventDistributionPrimAPrim"}; 

     for (int i = 0; i < 7; ++i) {
      fMixedEventTripletArray_TwoBody[i] =  new TH1F(histTitlesMixed_TwoBody[i],histTitlesMixed_TwoBody[i],  8000, 0, 8);
      fMixedEvent->Add(fMixedEventTripletArray_TwoBody[i]);
     } 

    // Mixed event multiplicity dist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fMixedEventMult = new TList();
    fMixedEventMult->SetOwner();
    fMixedEventMult->SetName("MixedEventMult");

    fMixedEventTripletMultArray = new TH2F*[10];
    TString histTitlesMixedMult[10] = {"mixedEventDistributionMultPPPrim","mixedEventDistributionMultAPAPAPrim","mixedEventDistributionMultPPP", "mixedEventDistributionMultAPAPAP", "mixedEventDistributionMultPPAPrim","mixedEventDistributionMultAPAPPrim", "mixedEventDistributionMultPPrimPrim", "mixedEventDistributionMultAPAPrimAPrim", "mixedEventDistributionMultPAPrimAPrim", "mixedEventDistributionMultAPPrimPrim"}; 
    
    for (int i = 0; i < 10; ++i) {
      fMixedEventTripletMultArray[i] = new TH2F(histTitlesMixedMult[i],histTitlesMixedMult[i], 8000, 0, 8,26,1,27);
      fMixedEventMult->Add(fMixedEventTripletMultArray[i]);
     }

     fMixedEventTripletMultArray_TwoBody = new TH2F*[7];
    TString histTitlesMixedMult_TwoBody[7] = {"MixedEventDistributionMultPP","MixedEventDistributionMultAPAP", "MixedEventDistributionMultPPrim", "MixedEventDistributionMultAPAPrim", "MixedEventDistributionMultPAPrim", "MixedEventDistributionMultAPPrim", "MixedEventDistributionMultPrimAPrim"}; 

     for (int i = 0; i < 7; ++i) {
      fMixedEventTripletMultArray_TwoBody[i] =  new TH2F(histTitlesMixedMult_TwoBody[i],histTitlesMixedMult_TwoBody[i], 8000, 0, 8,26,1,27);
      fMixedEventMult->Add(fMixedEventTripletMultArray_TwoBody[i]);
     }

    fResultsThreeBody->Add(fSameEvent);
    fResultsThreeBody->Add(fMixedEvent);
    fResultsThreeBody->Add(fSameEventMult);
    fResultsThreeBody->Add(fMixedEventMult);
    
    // Same event phi theta distribution ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fSameEventPhiTheta = new TList();
    fSameEventPhiTheta->SetOwner();
    fSameEventPhiTheta->SetName("SameEventPhiTheta");

    fSameEventTripletPhiThetaArray = new TH2F*[56];
    TString histTitlesSamePhiEta[28] ={"sameEventPhiEtaMultPPPrim","sameEventPhiEtaMultAPAPAPrim",
      "sameEventPhiEtaMultPPP", "sameEventPhiEtaMultAPAPAP", "sameEventPhiEtaMultPPAPrim","sameEventPhiEtaMultAPAPPrim","sameEventPhiEtaMultPPrimPrim", "sameEventPhiEtaMultAPAPrimAPrim", "sameEventPhiEtaMultPAPrimAPrim", "sameEventPhiEtaMultAPPrimPrim","sameEventPhiEtaMultPPSamePrimMixed", "sameEventPhiEtaMultPPrimSamePMixed", "sameEventPhiEtaMultAPAPSameAPrimMixed", "sameEventPhiEtaMultAPAPrimSameAPMixed", "sameEventPhiEtaMultPPSamePMixed", "sameEventPhiEtaMultAPAPSameAPMixed", "sameEventPhiEtaMultPPSameAPrimMixed", "sameEventPhiEtaMultPAPrimSamePMixed", "sameEventPhiEtaMultAPAPSamePrimMixed", "sameEventPhiEtaMultAPPrimSameAPMixed", "sameEventPhiEtaMultPPrimSamePimMixed", "sameEventPhiEtaMultPrimPrimSamePMixed", "sameEventPhiEtaMultAPAPrimSameAPMixed", "sameEventPhiEtaMultAPrimAPrimSameAPMixed", "sameEventPhiEtaMultPAPrimSameAPrimMixed", "sameEventPhiEtaMultAPrimAPrimSamePMixed", "sameEventPhiEtaMultAPPrimSamePrimMixed", "sameEventPhiEtaMultPrimPrimSameAPMixed"}; 

    for(int i=0;i<28;i++){
      fSameEventTripletPhiThetaArray[i] = new TH2F(histTitlesSamePhiEta[i]+"Before",histTitlesSamePhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventTripletPhiThetaArray[28+i] = new TH2F(histTitlesSamePhiEta[i]+"After",histTitlesSamePhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[i]);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray[28+i]);
    }

    fSameEventTripletPhiThetaArray_TwoBody = new TH2F*[14];
    TString histTitlesSamePhiEta_TwoBody[7] ={"sameEventPhiEtaMultPP", "sameEventPhiEtaMultAPAP", "sameEventPhiEtaMultPPrim", "sameEventPhiEtaMultAPAPrim", "sameEventPhiEtaMultPAPrim", "sameEventPhiEtaMultAPPrim", "sameEventPhiEtaMultPrimAPrim",}; 

    for(int i=0;i<7;i++){
      fSameEventTripletPhiThetaArray_TwoBody[i] = new TH2F(histTitlesSamePhiEta_TwoBody[i]+"Before",histTitlesSamePhiEta_TwoBody[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventTripletPhiThetaArray_TwoBody[7+i] = new TH2F(histTitlesSamePhiEta_TwoBody[i]+"After",histTitlesSamePhiEta_TwoBody[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray_TwoBody[i]);
      fSameEventPhiTheta->Add(fSameEventTripletPhiThetaArray_TwoBody[7+i]);
    } 


    // Mixed event phi theta distribution
    fMixedEventPhiTheta = new TList();
    fMixedEventPhiTheta->SetOwner();
    fMixedEventPhiTheta->SetName("MixedEventPhiTheta");

    fMixedEventTripletPhiThetaArray = new TH2F*[20];
    TString histTitlesMixedPhiEta[10] = {"mixedEventPhiEtaPPPrim","mixedEventPhiEtaAPAPAPrim",
      "mixedEventPhiEtaPPP", "mixedEventPhiEtaAPAPAP", "mixedEventPhiEtaPPAPrim","mixedEventPhiEtaAPAPPrim", "mixedEventPhiEtaPPrimPrim", "mixedEventPhiEtaAPAPrimAPrim", "mixedEventPhiEtaPAPrimAPrim", "mixedEventPhiEtaAPPrimPrim"}; 

    for(int i=0;i<10;i++){
      fMixedEventTripletPhiThetaArray[i] = new TH2F(histTitlesMixedPhiEta[i]+"Before",histTitlesMixedPhiEta[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventTripletPhiThetaArray[10+i] = new TH2F(histTitlesMixedPhiEta[i]+"After",histTitlesMixedPhiEta[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[i]);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray[10+i]);
    }

    fMixedEventTripletPhiThetaArray_TwoBody = new TH2F*[14];
    TString histTitlesMixedPhiEta_TwoBody[7] ={"mixedEventPhiEtaMultPP", "mixedEventPhiEtaMultAPAP", "mixedEventPhiEtaMultPPrim", "mixedEventPhiEtaMultAPAPrim", "mixedEventPhiEtaMultPAPrim", "mixedEventPhiEtaMultAPPrim", "mixedEventPhiEtaMultPrimAPrim",}; 

    for(int i=0;i<7;i++){
      fMixedEventTripletPhiThetaArray_TwoBody[i] = new TH2F(histTitlesMixedPhiEta_TwoBody[i]+"Before",histTitlesMixedPhiEta_TwoBody[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventTripletPhiThetaArray_TwoBody[7+i] = new TH2F(histTitlesMixedPhiEta_TwoBody[i]+"After",histTitlesMixedPhiEta_TwoBody[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray_TwoBody[i]);
      fMixedEventPhiTheta->Add(fMixedEventTripletPhiThetaArray_TwoBody[7+i]);
    } 

    fResultsThreeBody->Add(fSameEventPhiTheta);
    fResultsThreeBody->Add(fMixedEventPhiTheta);

    // Q3 vs q12 plot for theory
    if(fRunPlotQ3Vsq){
      fQ3Vskq12 = new TList();
      fQ3Vskq12->SetOwner();
      fQ3Vskq12->SetName("Q3Vskq12");

      fQ3VskDistributionsArrayq12 =  new TH2F*[13]; //GANESHA CHANGE!!!!
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

      fQ3VskDistributionsArrayq12Mixed =  new TH2F*[6];
      TString histTitlesfQ3VskDistributionsMixed[6] =  {"mixedEventDistributionPPPrim","mixedEventDistributionAPAPAPrim",
        "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPPAPrim","mixedEventDistributionAPAPPrim" }; 
      for(int i=0;i<6;i++){
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

      fQ3VskDistributionsArrayq23Mixed =  new TH2F*[6];
      TString histTitlesfQ3VskDistributionsMixed2[6] =  {"mixedEventDistributionPPPrim","mixedEventDistributionAPAPAPrim",
        "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPPAPrim","mixedEventDistributionAPAPPrim"}; 
      for(int i=0;i<6;i++){
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
  
      fInvMassSame = new TH2F*[10];
      TString histTitlesInvMassSame[10] ={"InvMassSameEventDistributionPPPrim","InvMassSameEventDistributionAPAPAPrim",
        "InvMassSameEventDistributionPPP", "InvMassSameEventDistributionAPAPAP", "InvMassSameEventDistributionPPAPrim","InvMassSameEventDistributionAPAPPrim",  "InvMassSameEventDistributionMultPPrimPrim", "InvMassSameEventDistributionMultAPAPrimAPrim", "InvMassSameEventDistributionMultPAPrimAPrim", "InvMassSameEventDistributionMultAPPrimPrim"}; 
      for(int i=0;i<10;i++){
        fInvMassSame[i] = new TH2F(histTitlesInvMassSame[i],histTitlesInvMassSame[i], 33, 0, 0.99, 400, 1.0, 1.2);
        fInvMassTripletSame->Add(fInvMassSame[i]);
      }

      fInvMassTripletMixed = new TList();
      fInvMassTripletMixed->SetOwner();
      fInvMassTripletMixed->SetName("InvMassTripletMixed");
       
      fInvMassMixed =new TH2F*[10];
      TString histTitlesInvMassMixed[10] = {"InvMassMixedEventDistributionPPPrim","InvMassMixedEventDistributionAPAPAPrim",
        "InvMassMixedEventDistributionPPP", "InvMassMixedEventDistributionAPAPAP", "InvMassMixedEventDistributionPPAPrim","InvMassMixedEventDistributionAPAPPrim",  "InvMassMixedEventDistributionMultPPrimPrim", "InvMassMixedEventDistributionMultAPAPrimAPrim", "InvMassMixedEventDistributionMultPAPrimAPrim", "InvMassMixedEventDistributionMultAPPrimPrim"};

      for(int i=0;i<10;i++){
        fInvMassMixed[i] = new TH2F(histTitlesInvMassMixed[i],histTitlesInvMassMixed[i], 33, 0, 0.99, 400, 1.0, 1.2);
        fInvMassTripletMixed->Add(fInvMassMixed[i]);
      }

      fResultsThreeBody->Add(fInvMassTripletSame);
      fResultsThreeBody->Add(fInvMassTripletMixed);
    }     
  }//if (fRunThreeBody)

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
  PostData(4, fPrimaryList); 
  PostData(5, fAntiPrimaryList); 
  if(fCleanWithLambdas){
    PostData(6, fLambdaList);
    PostData(7, fAntiLambdaList);
  }
  PostData(8, fResults);
  PostData(9, fResultsQA);
  PostData(10, fResultsSample);
  PostData(11, fResultsSampleQA);
  PostData(12, fResultsThreeBody);
  if (fProton->GetIsMonteCarlo()) {
    if (!fProton->GetMinimalBooking()) {
      fProtonMCList = fProton->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCTrkCuts");
      fProtonMCList->SetOwner();
    }
    PostData(13, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    if (!fAntiProton->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProton->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiTrkCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(14, fAntiProtonMCList);
  }

  if (fPrimary->GetIsMonteCarlo()) { 
    if (!fPrimary->GetMinimalBooking()) {
      fPrimaryMCList = fPrimary->GetMCQAHists();
    } else {
      fPrimaryMCList = new TList();
      fPrimaryMCList->SetName("MCPrimaryTrkCuts");
      fPrimaryMCList->SetOwner();
    }
    PostData(15, fPrimaryMCList);
  }
  if (fAntiPrimary->GetIsMonteCarlo()) { 
    if (!fAntiPrimary->GetMinimalBooking()) {
      fAntiPrimaryMCList = fAntiPrimary->GetMCQAHists();
    } else {
      fAntiPrimaryMCList = new TList();
      fAntiPrimaryMCList->SetName("MCAntiPrimaryTrkCuts");
      fAntiPrimaryMCList->SetOwner();
    }
    PostData(16, fAntiPrimaryMCList);
  }

  if(fCleanWithLambdas){
    if (fLambda->GetIsMonteCarlo()) {
      if (!fLambda->GetMinimalBooking()) {
        fLambdaMCList = fLambda->GetMCQAHists();
      } else {
        fLambdaMCList = new TList();
        fLambdaMCList->SetName("MCv0Cuts");
        fLambdaMCList->SetOwner();
      }
      PostData(17, fLambdaMCList);
    }
    if (fAntiLambda->GetIsMonteCarlo()) {
      if (!fAntiLambda->GetMinimalBooking()) {
        fAntiLambdaMCList = fAntiLambda->GetMCQAHists();
      } else {
        fAntiLambdaMCList = new TList();
        fAntiLambdaMCList->SetName("MCAntiv0Cuts");
        fAntiLambdaMCList->SetOwner();
      }
      PostData(18, fAntiLambdaMCList);
    }
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


} //void AliAnalysisTaskThreeBodyProtonPrimary::UserCreateOutputObjects()

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::UserExec(Option_t *option) {
//  AliVEvent *fInputEvent = InputEvent();

//List of actions:
//	a) Particle Selections
//	   a.0) Proton
//	   a.1) General Primary (e.g. K+, K-, Pi+, Pi-)
//	   a.2) General optional V0 (Lambda) Selection
//      b) Optional pair cleaning (might be needed to remove Lambda -> proton pi-)
//      c) Start Three Body Calculus
//	   c.0) Same event distribution
//	   c.1) Mixed event distribution
//		c.1.0) Same 2, Mixed 1
//		c.1.1) Normal mixed
//	   c.2)

  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }
  //a) Particle Selection +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }//for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)

  //a.0) Proton Selection ------------------------------
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
  }//for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)

  //a.1) General Primary Selection ------------------------------
  std::vector<AliFemtoDreamBasePart> Primaries; 
  std::vector<AliFemtoDreamBasePart> AntiPrimaries;
  fPrimaryTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize); 
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fPrimaryTrack->SetTrack(track, fInputEvent);
    if (fPrimary->isSelected(fPrimaryTrack)) {
      Primaries.push_back(*fPrimaryTrack);
    }
    if (fAntiPrimary->isSelected(fPrimaryTrack)) {
      AntiPrimaries.push_back(*fPrimaryTrack);
    }
  }//for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)

  //a.2) General optional V0 (Lambda) Selection ------------------------------
  std::vector<AliFemtoDreamBasePart> Lambdas;
  std::vector<AliFemtoDreamBasePart> AntiLambdas;
  if(fCleanWithLambdas)
  {
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
	  }//for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)
  }

  //b) Optional pair cleaning +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 fPairCleaner->ResetArray();

  if(fCleanWithLambdas)
  {
     fPairCleaner->CleanTrackAndDecay(&Protons, &Lambdas, 0); 
     fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiLambdas, 1); 

     fPairCleaner->CleanDecay(&Lambdas, 0); 
     fPairCleaner->CleanDecay(&AntiLambdas, 1); 
  }

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Primaries); 
  fPairCleaner->StoreParticle(AntiPrimaries); 

  //c) Start Three Body Calculus +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(fRunThreeBody){
    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    int bins[2] = { 0, 0 };
    float ZVtx = fEvent->GetZVertex();
    float Mult = fEvent->GetMultiplicity();
    fPartColl->FindBin(ZVtx, Mult, bins);

    // c.0) Same event distribution -------------------------------- 
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles(); 

    // Proton Proton Primary
    FillTripletDistribution( ParticleVector, 2, 0, 0, fSameEventTripletArray[0],PDGCodes, bins[1],fSameEventTripletMultArray[0], fSameEventTripletPhiThetaArray,0, *fConfig);//, fInvMassSame[0]);
    // Antiproton Antiproton Antiprimary
    FillTripletDistribution( ParticleVector, 3, 1, 1, fSameEventTripletArray[1],PDGCodes, bins[1],fSameEventTripletMultArray[1], fSameEventTripletPhiThetaArray,1, *fConfig);//, fInvMassSame[1]);
    // Proton Proton Proton
    FillTripletDistribution( ParticleVector, 0, 0, 0, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletPhiThetaArray,2, *fConfig);//, fInvMassSame[2]);
    // Antiproton Antiproton Antiproton
    FillTripletDistribution( ParticleVector, 1, 1, 1, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletPhiThetaArray,3, *fConfig);//, fInvMassSame[3]);
    // Proton Proton AntiPrimary
    FillTripletDistribution( ParticleVector, 3, 0, 0, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletPhiThetaArray,4, *fConfig);//, fInvMassSame[4]);
    // Antiproton Antiproton Primary
    FillTripletDistribution( ParticleVector, 2, 1, 1, fSameEventTripletArray[5],PDGCodes, bins[1],fSameEventTripletMultArray[5], fSameEventTripletPhiThetaArray,5, *fConfig);//, fInvMassSame[5]);
    // Proton Primary Primary
    FillTripletDistribution( ParticleVector, 0, 2, 2, fSameEventTripletArray[6],PDGCodes, bins[1],fSameEventTripletMultArray[6], fSameEventTripletPhiThetaArray,6, *fConfig);//, fInvMassSame[5]);
    // Antiproton Antiprimary Antiprimary
    FillTripletDistribution( ParticleVector, 1, 3, 3, fSameEventTripletArray[7],PDGCodes, bins[1],fSameEventTripletMultArray[7], fSameEventTripletPhiThetaArray,7, *fConfig);//, fInvMassSame[5]);
    // Proton Antiprimary Antiprimary
    FillTripletDistribution( ParticleVector, 0, 3, 3, fSameEventTripletArray[8],PDGCodes, bins[1],fSameEventTripletMultArray[8], fSameEventTripletPhiThetaArray,8, *fConfig);//, fInvMassSame[5]);
    // Antiproton Primary Primary
    FillTripletDistribution( ParticleVector, 1, 2, 2, fSameEventTripletArray[9],PDGCodes, bins[1],fSameEventTripletMultArray[9], fSameEventTripletPhiThetaArray,9, *fConfig);//, fInvMassSame[5]);


    //Two Body Analyses...........

    //Proton Proton
    FillPairDistribution( ParticleVector, 0, 0, fSameEventTripletArray_TwoBody[0],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[0], fSameEventTripletPhiThetaArray_TwoBody,0, *fConfig);
    //AntiProton AntiProton
    FillPairDistribution( ParticleVector, 1, 1, fSameEventTripletArray_TwoBody[1],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[1], fSameEventTripletPhiThetaArray_TwoBody,1, *fConfig);
    //Proton Primary
    FillPairDistribution( ParticleVector, 0, 2, fSameEventTripletArray_TwoBody[2],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[2], fSameEventTripletPhiThetaArray_TwoBody,2, *fConfig);
    //Antiproton Antiprimary
    FillPairDistribution( ParticleVector, 1, 3, fSameEventTripletArray_TwoBody[3],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[3], fSameEventTripletPhiThetaArray_TwoBody,3, *fConfig);
    //Proton Antiprimary
    FillPairDistribution( ParticleVector, 0, 3, fSameEventTripletArray_TwoBody[4],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[4], fSameEventTripletPhiThetaArray_TwoBody,4, *fConfig);
    //Antiproton Primary
    FillPairDistribution( ParticleVector, 1, 2, fSameEventTripletArray_TwoBody[5],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[5], fSameEventTripletPhiThetaArray_TwoBody,5, *fConfig);
   //Primary Antiprimary
    FillPairDistribution( ParticleVector, 2, 3, fSameEventTripletArray_TwoBody[6],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[6], fSameEventTripletPhiThetaArray_TwoBody,6, *fConfig);
  
    // c.1) Mixed event distribution --------------------------------

    // TAKE CARE OF MULT AND ZVtx!
    if (!(bins[0] == -99 || bins[1] == -99)) {

      auto itZVtx = fPartContainer.begin()+ bins[0];
      auto itMult = itZVtx->begin() + bins[1];

      //c.1.0) Same 2, Mixed 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      // Proton Proton Primary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 2, fSameEventTripletArray[10], PDGCodes, bins[1],fSameEventTripletMultArray[10], fSameEventTripletPhiThetaArray, 10, *fConfig);//, fQ3VskDistributionsArrayq12[0],fQ3VskDistributionsArrayq23[0]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 0, fSameEventTripletArray[11], PDGCodes, bins[1],fSameEventTripletMultArray[11], fSameEventTripletPhiThetaArray, 11, *fConfig);//, fQ3VskDistributionsArrayq12[1],fQ3VskDistributionsArrayq23[1]);

      // Proton Proton Proton
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 0, fSameEventTripletArray[12], PDGCodes, bins[1],fSameEventTripletMultArray[12], fSameEventTripletPhiThetaArray, 12, *fConfig);//, fQ3VskDistributionsArrayq12[2],fQ3VskDistributionsArrayq23[2]);

      // Antiproton Antiproton Antiprimary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 3, fSameEventTripletArray[13], PDGCodes, bins[1],fSameEventTripletMultArray[13], fSameEventTripletPhiThetaArray, 13, *fConfig);//, fQ3VskDistributionsArrayq12[3],fQ3VskDistributionsArrayq23[3]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 1, fSameEventTripletArray[14], PDGCodes, bins[1],fSameEventTripletMultArray[14], fSameEventTripletPhiThetaArray, 14, *fConfig);//, fQ3VskDistributionsArrayq12[4],fQ3VskDistributionsArrayq23[4]);

      // Antiproton Antiproton Antiproton
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 1, fSameEventTripletArray[15], PDGCodes, bins[1],fSameEventTripletMultArray[15], fSameEventTripletPhiThetaArray, 15, *fConfig);//, fQ3VskDistributionsArrayq12[5],fQ3VskDistributionsArrayq23[5]);

      // Proton Proton AntiPrimary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 0, 3, fSameEventTripletArray[16], PDGCodes, bins[1],fSameEventTripletMultArray[16], fSameEventTripletPhiThetaArray, 16, *fConfig);//, fQ3VskDistributionsArrayq12[6],fQ3VskDistributionsArrayq23[6]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 3, 0, fSameEventTripletArray[17], PDGCodes, bins[1],fSameEventTripletMultArray[17], fSameEventTripletPhiThetaArray, 17, *fConfig);//, fQ3VskDistributionsArrayq12[7],fQ3VskDistributionsArrayq23[7]);

      // Antiproton Antiproton Primary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 1, 2, fSameEventTripletArray[18], PDGCodes, bins[1],fSameEventTripletMultArray[18], fSameEventTripletPhiThetaArray, 18, *fConfig);//, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 2, 1, fSameEventTripletArray[19], PDGCodes, bins[1],fSameEventTripletMultArray[19], fSameEventTripletPhiThetaArray, 19, *fConfig);//, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9]);

      // Proton Primary Primary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 2, 2, fSameEventTripletArray[20], PDGCodes, bins[1],fSameEventTripletMultArray[20], fSameEventTripletPhiThetaArray, 20, *fConfig);//, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 0, fSameEventTripletArray[21], PDGCodes, bins[1],fSameEventTripletMultArray[21], fSameEventTripletPhiThetaArray, 21, *fConfig);//, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9]);

      // Antiproton Antiprimary Antiprimary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 3, 3, fSameEventTripletArray[22], PDGCodes, bins[1],fSameEventTripletMultArray[22], fSameEventTripletPhiThetaArray, 22, *fConfig);//, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 1, fSameEventTripletArray[23], PDGCodes, bins[1],fSameEventTripletMultArray[23], fSameEventTripletPhiThetaArray, 23, *fConfig);//, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9]);

      // Proton Antiprimary Antiprimary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 0, 3, 3, fSameEventTripletArray[24], PDGCodes, bins[1],fSameEventTripletMultArray[24], fSameEventTripletPhiThetaArray, 24, *fConfig);//, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 3, 3, 0, fSameEventTripletArray[25], PDGCodes, bins[1],fSameEventTripletMultArray[25], fSameEventTripletPhiThetaArray, 25, *fConfig);//, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9]);

      // Antiproton Primary Primary
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 1, 2, 2, fSameEventTripletArray[26], PDGCodes, bins[1],fSameEventTripletMultArray[26], fSameEventTripletPhiThetaArray, 26, *fConfig);//, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8]);
      FillTripletDistributionSE2ME1(ParticleVector, *itMult, 2, 2, 1, fSameEventTripletArray[27], PDGCodes, bins[1],fSameEventTripletMultArray[27], fSameEventTripletPhiThetaArray, 27, *fConfig);//, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9]);

      //c.1.1) Normal mixing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // 2 (Anti-)Protons, 1 (Anti-)Primary
      FillTripletDistributionME(ParticleVector, *itMult, 2, 0, 0, fMixedEventTripletArray[0], PDGCodes, bins[1],fMixedEventTripletMultArray[0], fMixedEventTripletPhiThetaArray,0, *fConfig);//, fInvMassMixed[0], fQ3VskDistributionsArrayq12Mixed[0], fQ3VskDistributionsArrayq23Mixed[0]);
      FillTripletDistributionME(ParticleVector, *itMult, 3, 1, 1, fMixedEventTripletArray[1], PDGCodes, bins[1],fMixedEventTripletMultArray[1], fMixedEventTripletPhiThetaArray,1, *fConfig);//, fInvMassMixed[1], fQ3VskDistributionsArrayq12Mixed[1], fQ3VskDistributionsArrayq23Mixed[1]);
      FillTripletDistributionME(ParticleVector, *itMult, 0, 0, 0, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletPhiThetaArray,2, *fConfig);//, fInvMassMixed[2], fQ3VskDistributionsArrayq12Mixed[2], fQ3VskDistributionsArrayq23Mixed[2]);
      FillTripletDistributionME(ParticleVector, *itMult, 1, 1, 1, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletPhiThetaArray,3, *fConfig);//, fInvMassMixed[3], fQ3VskDistributionsArrayq12Mixed[3], fQ3VskDistributionsArrayq23Mixed[3]);
      FillTripletDistributionME(ParticleVector, *itMult, 3, 0, 0, fMixedEventTripletArray[4], PDGCodes, bins[1],fMixedEventTripletMultArray[4], fMixedEventTripletPhiThetaArray,4, *fConfig);//, fInvMassMixed[4], fQ3VskDistributionsArrayq12Mixed[4], fQ3VskDistributionsArrayq23Mixed[4]);
      FillTripletDistributionME(ParticleVector, *itMult, 2, 1, 1, fMixedEventTripletArray[5], PDGCodes, bins[1],fMixedEventTripletMultArray[5], fMixedEventTripletPhiThetaArray,5, *fConfig);//, fInvMassMixed[5], fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5]);

      // 1 (Anti-)Proton, 2 (Anti-)Primary
      FillTripletDistributionME(ParticleVector, *itMult, 0, 2, 2, fMixedEventTripletArray[6], PDGCodes, bins[1],fMixedEventTripletMultArray[6], fMixedEventTripletPhiThetaArray,6, *fConfig);//, fInvMassMixed[5], fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5]);
      FillTripletDistributionME(ParticleVector, *itMult, 1, 3, 3, fMixedEventTripletArray[7], PDGCodes, bins[1],fMixedEventTripletMultArray[7], fMixedEventTripletPhiThetaArray,7, *fConfig);//, fInvMassMixed[5], fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5]);
      FillTripletDistributionME(ParticleVector, *itMult, 0, 3, 3, fMixedEventTripletArray[8], PDGCodes, bins[1],fMixedEventTripletMultArray[8], fMixedEventTripletPhiThetaArray,8, *fConfig);//, fInvMassMixed[5], fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5]);
      FillTripletDistributionME(ParticleVector, *itMult, 1, 2, 2, fMixedEventTripletArray[9], PDGCodes, bins[1],fMixedEventTripletMultArray[9], fMixedEventTripletPhiThetaArray,9, *fConfig);//, fInvMassMixed[5], fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5]);
   
      //Two Body Analyses...........

      //Proton Proton
      FillPairDistributionME( ParticleVector, *itMult, 0, 0, fMixedEventTripletArray_TwoBody[0],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[0], fMixedEventTripletPhiThetaArray_TwoBody,0, *fConfig);
      //AntiProton AntiProton
      FillPairDistributionME( ParticleVector, *itMult, 1, 1, fMixedEventTripletArray_TwoBody[1],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[1], fMixedEventTripletPhiThetaArray_TwoBody,1, *fConfig);
      //Proton Primary
      FillPairDistributionME( ParticleVector, *itMult, 0, 2, fMixedEventTripletArray_TwoBody[2],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[2], fMixedEventTripletPhiThetaArray_TwoBody,2, *fConfig);
      //Antiproton Antiprimary
      FillPairDistributionME( ParticleVector, *itMult, 1, 3, fMixedEventTripletArray_TwoBody[3],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[3], fMixedEventTripletPhiThetaArray_TwoBody,3, *fConfig);
      //Proton Antiprimary
      FillPairDistributionME( ParticleVector, *itMult, 0, 3, fMixedEventTripletArray_TwoBody[4],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[4], fMixedEventTripletPhiThetaArray_TwoBody,4, *fConfig);
      //Antiproton Primary
      FillPairDistributionME( ParticleVector, *itMult, 1, 2, fMixedEventTripletArray_TwoBody[5],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[5], fMixedEventTripletPhiThetaArray_TwoBody,5, *fConfig);
      //Primary Antiprimary
      FillPairDistributionME( ParticleVector, *itMult, 2, 3, fMixedEventTripletArray_TwoBody[6],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[6], fMixedEventTripletPhiThetaArray_TwoBody,6, *fConfig);

      // Update the particle container with current event
      SetMixedEvent(ParticleVector, &(*itMult));
      //SetMixedEventOnlyPLambdaTEST(ParticleVector, &(*itMultTEST));
      //SetMixedEventOnlyPPLambdaTEST( ParticleVector, &(*itMultTESTppL));
      //SetMixedEventOnlyPPPTEST( ParticleVector, &(*itMultTESTppp));

    }//if (!(bins[0] == -99 || bins[1] == -99))

  }//if(fRunThreeBody)



  /*if (fPairCleaner->GetCounter() > 0) {
    if (fConfig->GetUseEventMixing()) {
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent);
    }
    if (fConfig->GetUsePhiSpinning()) {
      fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
    }
  }*/

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPrimaryList); 
  PostData(5, fAntiPrimaryList); 
  if(fCleanWithLambdas){
    PostData(6, fLambdaList);
    PostData(7, fAntiLambdaList);
  }
  //PostData(8, fResults);
  //PostData(9, fResultsQA);
  //PostData(10, fResultsSample);
  //PostData(11, fResultsSampleQA);
  PostData(12, fResultsThreeBody);
  if (fProton->GetIsMonteCarlo()) {
    PostData(13, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    PostData(14, fAntiProtonMCList);
  }
  if (fPrimary->GetIsMonteCarlo()) { 
    PostData(15, fProtonMCList);
  }
  if (fAntiPrimary->GetIsMonteCarlo()) { 
    PostData(16, fAntiProtonMCList);
  }
  if(fCleanWithLambdas){
    if (fLambda->GetIsMonteCarlo()) {
      PostData(17, fLambdaMCList);
    }
    if (fAntiLambda->GetIsMonteCarlo()) {
      PostData(18, fAntiLambdaMCList);
    }
  }

}

//==================================================================================================================================================
void AliAnalysisTaskThreeBodyProtonPrimary::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//==================================================================================================================================================
void AliAnalysisTaskThreeBodyProtonPrimary::StoreGlobalTrackReference(AliVTrack *track) {
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

//==================================================================================================================================================

TLorentzVector AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(
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

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* InvMassSame){
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
        TLorentzVector q12 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);
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
            Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
            Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
            Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
          }
          if(!fClosePairRejectionForAll){
            if(fClosePairRejectionPPPorPPL){
              if(DoThisPair12==11){
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair23==11){
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair31==11){
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              }
            }

            if(!fClosePairRejectionPPPorPPL){
              if(DoThisPair12==21||DoThisPair12==12){
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair23==21||DoThisPair23==12){
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              }
              if(DoThisPair31==21||DoThisPair31==12){
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              }
            }

          }

        }


        if(!Pair12||!Pair23||!Pair31) {continue;}

        hist->Fill(Q3);
        hist2d->Fill(Q3,mult+1);

        /*if(firstSpecies == 2 || firstSpecies == 3) InvMassSame->Fill(Q3, iPart1->GetInvMass());
        if(secondSpecies == 2 || secondSpecies == 3) InvMassSame->Fill(Q3, iPart2->GetInvMass());
        if(thirdSpecies == 2 || thirdSpecies == 3) InvMassSame->Fill(Q3, iPart3->GetInvMass());*/

      }
    }
  }
}

//================================================================================================================================================== //GANESHA CHECK!!
void AliAnalysisTaskThreeBodyProtonPrimary::FillPairDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* InvMassSame){
  // Two Body Same Event 

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;

  // Get the PID codes std::vector<int>
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;

  // Get particle masses
  auto massparticle1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massparticle2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();

  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;

  if(abs(*itPDGPar1)==3122) DaughterPart1 = 2;
  if(abs(*itPDGPar1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==3122) DaughterPart2 = 2;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;

  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;

  // loop over first particle
  for (auto iPart1 = Particle1Vector->begin(); iPart1 != Particle1Vector->end(); ++iPart1) {
    // if second particle species is different than first - start with the first particle in the vector
    auto iPart2 = Particle2Vector->begin();
    // if second particle  and first are the species, start second loop from the next particle (to not double count)
    if (firstSpecies==secondSpecies) iPart2 = iPart1+1;
    // loop over second particle ...
    for (; iPart2 != Particle2Vector->end(); ++iPart2) {

        // Now we have the three particles, lets create their Lorentz vectors 
        TLorentzVector Particle1_LV, Particle2_LV;
        Particle1_LV.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(), massparticle1);
        Particle2_LV.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(), massparticle2);
        // Get momentum
        float RelativeMomentum = AliFemtoDreamHigherPairMath::RelativePairMomentum(Particle1_LV, Particle2_LV);

        bool Pair12 = true;

        if(!fturnoffClosePairRejectionCompletely){

            Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[7+phiEtaHistNo],Config); //GANESHA  Check about DeltaEtaDeltaPhi Function. Might need to change for 2 Body
        }

        if(!Pair12) {continue;}

        hist->Fill(RelativeMomentum);
        hist2d->Fill(RelativeMomentum,mult+1); 
    }
  }
} //void AliAnalysisTaskThreeBodyProtonPrimary::FillPairDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config)

//==================================================================================================================================================

//void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassSame, TH2F* InvMassDET,TH2F* InvMassPDG) 
  
//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  for(unsigned int iSpecies = 0; iSpecies<ParticleVector.size(); iSpecies++){
    if ((ParticleVector.begin()+iSpecies)->size() > 0) {
      (fPartContainer->begin()+iSpecies)->SetEvent(*(ParticleVector.begin()+iSpecies));
    }
  }
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventOnlyPLambdaTEST(
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

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventOnlyPPLambdaTEST(
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

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventOnlyPPPTEST(
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

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed){
  // Description of function given in AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistribution
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


            TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
            part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),
            iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE,2)));
            part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),
            iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleME1,2)));
            part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(),
            iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME2,2)));
            // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
            TLorentzVector q12 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);
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
                Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
              }
              if(!fClosePairRejectionForAll){
                if(fClosePairRejectionPPPorPPL){
                  if(DoThisPair12==11){
                    Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair23==11){
                    Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair31==11){
                    Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                  }
                }
                if(!fClosePairRejectionPPPorPPL){
                  if(DoThisPair12==21||DoThisPair12==12){
                    Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair23==21||DoThisPair23==12){
                    Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                  }
                  if(DoThisPair31==21||DoThisPair31==12){
                    Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config, Q3);
                  }
                }
              }

            }

            if(!Pair12||!Pair23||!Pair31) {continue;}


            // Now we have the three particles, lets create their Lorentz vectors

            hist->Fill(Q3);
            hist2d->Fill(Q3,mult+1);
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

//==================================================================================================================================================

//void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionMEPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F* InvMassDET,TH2F* InvMassPDG)

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionMETEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){
  // Description of function given in AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistribution
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
              Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config);
            }
            //if(speciesME1==speciesME2){
            if(DoThisPair23==11){
              Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,false,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config);
            }
            //if(speciesME2==speciesSE){
            if(DoThisPair31==11){
              Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,false,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[10+phiEtaHistNo],Config);
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
            TLorentzVector q12 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);
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

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){ //GANESHA CHECK

  //  Check if reproduces the FemtoDream framework result
  auto ParticleSE = ParticleVector.begin()+speciesSE;
  auto MixedEvent1Container = fPartContainer.begin()+speciesME1;

  // Get the PID codes std::vector<int>
  auto itPDGParSE = PDGCodes.begin()+speciesSE;
  auto itPDGParME1 = PDGCodes.begin()+speciesME1;

  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;

  if(abs(*itPDGParSE)==3122) DaughterPart1 = 2;
  if(abs(*itPDGParSE)==2212) DaughterPart1 = 1;
  if(abs(*itPDGParME1)==3122) DaughterPart2 = 2;
  if(abs(*itPDGParME1)==2212) DaughterPart2 = 1;

  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;

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

          bool Pair12 = true;

         if(!fturnoffClosePairRejectionCompletely){

            Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[7+phiEtaHistNo],Config); //GANESHA Change number + Check about DeltaEtaDeltaPhi Function. Might need to change for 2 Body
            }

            if(!Pair12) {continue;}  

        }//if(speciesSE==0 && speciesME1 ==0)


        // Now we have the three particles, lets create their Lorentz vectors 
        TLorentzVector part1_LorVec, part2_LorVec;
        part1_LorVec.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(),massParticleSE);
        part2_LorVec.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(),massParticleME1);
        // Get momentum
        float RelativeK = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec, part2_LorVec);
        // No need to check pair selection because p lambda
        hist->Fill(RelativeK);
        hist2d->Fill(RelativeK,mult+1);
      }
    }
  }
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23){
  // Description of function given in AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistribution
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


          TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
          part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),
          iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE1,2)));
          part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),
          iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleSE2,2)));
          part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(),
          iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME,2)));
          // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
          TLorentzVector q12 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
          TLorentzVector q23 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
          TLorentzVector q31 = AliAnalysisTaskThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);
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
              Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
              Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
            }
            if(!fClosePairRejectionForAll){
              if(fClosePairRejectionPPPorPPL){
                if(DoThisPair12==11){
                  Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair23==11){
                  Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair31==11){
                  Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
                }
              }
              if(!fClosePairRejectionPPPorPPL){
                if(DoThisPair12==21||DoThisPair12==12){
                  Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair23==21||DoThisPair23==12){
                  Pair23 = DeltaEtaDeltaPhi(*iPart2,*iPart3,true,  DoThisPair23, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
                }
                if(DoThisPair31==21||DoThisPair31==12){
                  Pair31 = DeltaEtaDeltaPhi(*iPart3,*iPart1,true,  DoThisPair31, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[28+phiEtaHistNo],Config, Q3);
                }
              }
            }

          }

          if(!Pair12||!Pair23||!Pair31) {continue;}


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

//==================================================================================================================================================

bool AliAnalysisTaskThreeBodyProtonPrimary::DeltaEtaDeltaPhi(
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

//==================================================================================================================================================

bool AliAnalysisTaskThreeBodyProtonPrimary::DeltaEtaDeltaPhi(
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
        if(Q3<fQ3LimitForDeltaPhiDeltaEta){
          beforeHist->Fill(dphiAvg/ (float) size, deta);
        }
      }
      if (pass) {
        if ((dphiAvg / (float) size) * (dphiAvg / (float) size) / fDeltaPhiSqMax
            + deta * deta / fDeltaEtaSqMax < 1.) {
          pass = false;
        }
        else{
          if(fRunPlotPhiTheta){
            if(Q3<fQ3LimitForDeltaPhiDeltaEta){
              afterHist->Fill(dphiAvg/ (float) size, deta);
            }
          }
        }
      }
    }
  }
  return pass;
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillPairInvMass( AliFemtoDreamBasePart &part1,
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

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillPDGPairInvMass( AliFemtoDreamBasePart &part1, float massPart1,
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
