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
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"


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
      fSameEventPtPrimaries(nullptr),
      fSameEventPtProtons(nullptr),
      fSameEventPtPrimaries2(nullptr),
      fSameEventPtProtons2(nullptr),
      fMixedEventMult(nullptr),
      fMixedEventPtPrimaries(nullptr),
      fMixedEventPtProtons(nullptr),
      fMixedEventPtPrimaries2(nullptr),
      fMixedEventPtProtons2(nullptr),
      fSameEventPhiTheta_SamePair(nullptr),
      fMixedEventPhiTheta_SamePair(nullptr),
      fSameEventPhiTheta_DifferentPair(nullptr),
      fMixedEventPhiTheta_DifferentPair(nullptr),
      fQ3Vskq12(nullptr),
      fQ3Vskq12Mixed(nullptr),
      fQ3Vskq23(nullptr),
      fQ3Vskq23Mixed(nullptr),
      fOtherHistos(nullptr),
      fInvMassList(nullptr),
      fKinematicsPlots(nullptr),
      fRunThreeBody(true),
      fRunPlotInvMass(false),
      fRunPlotQ3Vsq(true),
      fRunPlotPhiTheta(true),
      fRunPlotOtherHistos(true),
      fRunPlotMult(true),
      fRunPlotPt(false),
      fPlotsMC(false),
      fRunOfficialTwoBody(false), // ADDED BY RAFFA
      fDoKinematicsPlots(false),
      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fClosePairRejectionPPPorPPL(false),
      fQ3LimitForDeltaPhiDeltaEta(1.),
      fDeltaPhiMaxPP(0.017),
      fDeltaEtaMaxPP(0.017),
      fDeltaPhiMaxPPrim(0.04),
      fDeltaEtaMaxPPrim(0.012),
      fDeltaPhiMaxPAPrim(0.04),
      fDeltaEtaMaxPAPrim(0.012),
      fCleanWithLambdas(false),
      fDoOnlyThreeBody(true),
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPtPrimaries(nullptr),
      fSameEventTripletPtProtons2(nullptr),
      fSameEventTripletPtPrimaries2(nullptr),
      fSameEventTripletPtProtons(nullptr),
      fSameEventTripletPtvsQ3Primaries(nullptr),
      fSameEventTripletPtvsQ3Protons(nullptr),
      fSameEventTripletPhiThetaArray_SamePair(nullptr), 
      fSameEventTripletPhiThetaArray_DifferentPair(nullptr),
      fSameEventTripletArray_TwoBody(nullptr),
      fSameEventTripletMultArray_TwoBody(nullptr),
      fSameEventTripletPhiThetaArray_TwoBody(nullptr),
      fPairTranverseMass_TwoBody(nullptr), // ADDED BY RAFFA
      fPairTranverseMassVSkstar_TwoBody(nullptr), // ADDED BY RAFFA
      fPartContainer(0),
      fPartContainerPPP(0),
      fPartContainerPPPrim(0),
      fPartContainerPPAPrim(0),
      fPartContainerPP(0),
      fPartContainerPPrim(0),
      fPartContainerPAPrim(0),
      fVectPartContainers(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPtPrimaries(nullptr),
      fMixedEventTripletPtProtons(nullptr),
      fMixedEventTripletPtPrimaries2(nullptr),
      fMixedEventTripletPtProtons2(nullptr),
      fMixedEventTripletPtvsQ3Primaries(nullptr),
      fMixedEventTripletPtvsQ3Protons(nullptr),
      fMixedEventTripletPhiThetaArray_SamePair(nullptr), 
      fMixedEventTripletPhiThetaArray_DifferentPair(nullptr),
      fMixedEventTripletArray_TwoBody(nullptr),
      fMixedEventTripletMultArray_TwoBody(nullptr),
      fMixedEventTripletPhiThetaArray_TwoBody(nullptr),
      fQ3VskDistributionsArrayq12(nullptr),
      fQ3VskDistributionsArrayq12Mixed(nullptr),
      fQ3VskDistributionsArrayq23(nullptr),
      fQ3VskDistributionsArrayq23Mixed(nullptr),
      fDoubletVsTrippletPPP(nullptr),
      fInvMass(nullptr),
      fKinematics(nullptr),
      fPrimAngles(nullptr),
      fDeta(nullptr),
      fDphi(nullptr),
      fKinematicsME(nullptr),
      fPrimAnglesME(nullptr),
      fDetaME(nullptr),
      fDphiME(nullptr),
      fTripletInvMassDet(nullptr),
      fTripletInvMassPDG(nullptr),
      fTripletInvMassDetMixed(nullptr),
      fTripletInvMassPDGMixed(nullptr),
      fTripletInvMassDetAnti(nullptr),
      fTripletInvMassPDGAnti(nullptr),
      fTripletInvMassDetMixedAnti(nullptr),
      fTripletInvMassPDGMixedAnti(nullptr),
      fpTvsEtaTrueKaons(nullptr),
      fpTvsEtaTrueAntiKaons(nullptr),
      fpTvsEtaTrueProtons(nullptr),
      fpTvsEtaTrueAntiProtons(nullptr),
      fpTvsEtaRecoKaons(nullptr),
      fpTvsEtaRecoAntiKaons(nullptr),
      fpTvsEtaRecoProtons(nullptr),
      fpTvsEtaRecoAntiProtons(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fStandardMixing(true),
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
      fSameEventPtPrimaries(nullptr),
      fSameEventPtProtons(nullptr),
      fSameEventPtPrimaries2(nullptr),
      fSameEventPtProtons2(nullptr),
      fMixedEventMult(nullptr),
      fMixedEventPtPrimaries(nullptr),
      fMixedEventPtProtons(nullptr),
      fMixedEventPtPrimaries2(nullptr),
      fMixedEventPtProtons2(nullptr),
      fSameEventPhiTheta_SamePair(nullptr),
      fMixedEventPhiTheta_SamePair(nullptr),
      fSameEventPhiTheta_DifferentPair(nullptr),
      fMixedEventPhiTheta_DifferentPair(nullptr),
      fQ3Vskq12(nullptr),
      fQ3Vskq12Mixed(nullptr),
      fQ3Vskq23(nullptr),
      fQ3Vskq23Mixed(nullptr),
      fOtherHistos(nullptr),
      fInvMassList(nullptr),
      fKinematicsPlots(nullptr),
      fRunThreeBody(true),
      fRunPlotInvMass(false),
      fRunPlotQ3Vsq(false),
      fRunPlotPhiTheta(true),
      fRunPlotOtherHistos(true),
      fRunPlotMult(true),
      fRunPlotPt(false),
      fPlotsMC(false),
      fRunOfficialTwoBody(false), // ADDED BY RAFFA
      fDoKinematicsPlots(false),
      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fClosePairRejectionPPPorPPL(false),
      fQ3LimitForDeltaPhiDeltaEta(1.),
      fDeltaPhiMaxPP(0.017),
      fDeltaEtaMaxPP(0.017),
      fDeltaPhiMaxPPrim(0.04),
      fDeltaEtaMaxPPrim(0.012),
      fDeltaPhiMaxPAPrim(0.04),
      fDeltaEtaMaxPAPrim(0.012),
      fCleanWithLambdas(false),
      fDoOnlyThreeBody(true),
      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletPtPrimaries(nullptr),
      fSameEventTripletPtProtons(nullptr),
      fSameEventTripletPtPrimaries2(nullptr),
      fSameEventTripletPtProtons2(nullptr),
      fSameEventTripletPtvsQ3Primaries(nullptr),
      fSameEventTripletPtvsQ3Protons(nullptr),
      fSameEventTripletPhiThetaArray_SamePair(nullptr), 
      fSameEventTripletPhiThetaArray_DifferentPair(nullptr),
      fSameEventTripletArray_TwoBody(nullptr),
      fSameEventTripletMultArray_TwoBody(nullptr),
      fSameEventTripletPhiThetaArray_TwoBody(nullptr),
      fPairTranverseMass_TwoBody(nullptr), // ADDED BY RAFFA
      fPairTranverseMassVSkstar_TwoBody(nullptr), // ADDED BY RAFFA
      fPartContainer(0),
      fPartContainerPPP(0),
      fPartContainerPPPrim(0),
      fPartContainerPPAPrim(0),
      fPartContainerPP(0),
      fPartContainerPPrim(0),
      fPartContainerPAPrim(0),
      fVectPartContainers(0),
      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletPtPrimaries(nullptr),
      fMixedEventTripletPtProtons(nullptr),
      fMixedEventTripletPtPrimaries2(nullptr),
      fMixedEventTripletPtProtons2(nullptr),
      fMixedEventTripletPtvsQ3Primaries(nullptr),
      fMixedEventTripletPtvsQ3Protons(nullptr),
      fMixedEventTripletPhiThetaArray_SamePair(nullptr), 
      fMixedEventTripletPhiThetaArray_DifferentPair(nullptr),
      fMixedEventTripletArray_TwoBody(nullptr),
      fMixedEventTripletMultArray_TwoBody(nullptr),
      fMixedEventTripletPhiThetaArray_TwoBody(nullptr),
      fQ3VskDistributionsArrayq12(nullptr),
      fQ3VskDistributionsArrayq12Mixed(nullptr),
      fQ3VskDistributionsArrayq23(nullptr),
      fQ3VskDistributionsArrayq23Mixed(nullptr),
      fDoubletVsTrippletPPP(nullptr),
      fInvMass(nullptr),
      fKinematics(nullptr),
      fPrimAngles(nullptr),
      fDeta(nullptr),
      fDphi(nullptr),
      fKinematicsME(nullptr),
      fPrimAnglesME(nullptr),
      fDetaME(nullptr),
      fDphiME(nullptr),
      fTripletInvMassDet(nullptr),
      fTripletInvMassPDG(nullptr),
      fTripletInvMassDetMixed(nullptr),
      fTripletInvMassPDGMixed(nullptr),
      fTripletInvMassDetAnti(nullptr),
      fTripletInvMassPDGAnti(nullptr),
      fTripletInvMassDetMixedAnti(nullptr),
      fTripletInvMassPDGMixedAnti(nullptr),
      fpTvsEtaTrueKaons(nullptr),
      fpTvsEtaTrueAntiKaons(nullptr),
      fpTvsEtaTrueProtons(nullptr),
      fpTvsEtaTrueAntiProtons(nullptr),
      fpTvsEtaRecoKaons(nullptr),
      fpTvsEtaRecoAntiKaons(nullptr),
      fpTvsEtaRecoProtons(nullptr),
      fpTvsEtaRecoAntiProtons(nullptr),
      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fStandardMixing(true),
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
      fProton->GetIsMonteCarlo() || fAntiProton->GetIsMonteCarlo()||fPrimary->GetIsMonteCarlo() || fAntiPrimary->GetIsMonteCarlo());

/* 
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
*/
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
  /*
  if(fCleanWithLambdas){
     fLambdaList = fLambda->GetQAHists();
     fAntiLambdaList = fAntiLambda->GetQAHists();
  }
*/
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


    //...............................................................................................
    // Same event distributions 1D
    fSameEvent = new TList();
    fSameEvent->SetOwner();
    fSameEvent->SetName("SameEvent");

    fSameEventTripletArray = new TH1F*[28];
    TString histTitlesSame[28] = {"sameEventDistributionPPPrim","sameEventDistributionAPAPAPrim", "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPPAPrim","sameEventDistributionAPAPPrim", "sameEventDistributionPPrimPrim", "sameEventDistributionAPAPrimAPrim", "sameEventDistributionPAPrimAPrim", "sameEventDistributionAPPrimPrim","sameEventDistributionPPSamePrimMixed", "sameEventDistributionPPrimSamePMixed", "sameEventDistributionAPAPSameAPrimMixed", "sameEventDistributionAPAPrimSameAPMixed", "sameEventDistributionPPSamePMixed", "sameEventDistributionAPAPSameAPMixed", "sameEventDistributionPPSameAPrimMixed", "sameEventDistributionPAPrimSamePMixed", "sameEventDistributionAPAPSamePrimMixed", "sameEventDistributionAPPrimSameAPMixed", "sameEventDistributionPPrimSamePrimMixed", "sameEventDistributionPrimPrimSamePMixed", "sameEventDistributionAPAPrimSameAPrimMixed", "sameEventDistributionAPrimAPrimSameAPMixed", "sameEventDistributionPAPrimSameAPrimMixed", "sameEventDistributionAPrimAPrimSamePMixed", "sameEventDistributionAPPrimSamePrimMixed", "sameEventDistributionPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
	for (int i = 0; i < 10; ++i) {
	   fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
	   fSameEvent->Add(fSameEventTripletArray[i]);
	}
	for (int i = 10; i <28 ; ++i) {
	   fSameEventTripletArray[i] =  new TH1F(histTitlesSame[i],histTitlesSame[i], 8000, 0, 8);
	   fSameEvent->Add(fSameEventTripletArray[i]);
	}
    }

    fSameEventTripletArray_TwoBody = new TH1F*[7];
    TString histTitlesSame_TwoBody[7] = {"sameEventDistributionPP","sameEventDistributionAPAP", "sameEventDistributionPPrim", "sameEventDistributionAPAPrim", "sameEventDistributionPAPrim", "sameEventDistributionAPPrim", "sameEventDistributionPrimAPrim"};

      fPairTranverseMass_TwoBody = new TH1F*[7]; // ADDED BY RAFFA
      TString histTitlesPairTranverseMass_TwoBody[7] = {"PairTranverseMassDistributionPP","PairTranverseMassDistributionAPAP", "PairTranverseMassDistributionPPrim", "PairTranverseMassDistributionAPAPrim", "PairTranverseMassDistributionPAPrim", "PairTranverseMassDistributionAPPrim", "PairTranverseMassDistributionPrimAPrim"}; // ADDED BY RAFFA

    if(!fDoOnlyThreeBody){
	for (int i = 0; i < 7; ++i) {
	  fSameEventTripletArray_TwoBody[i] =  new TH1F(histTitlesSame_TwoBody[i],histTitlesSame_TwoBody[i], 8000, 0, 8);
	  fSameEvent->Add(fSameEventTripletArray_TwoBody[i]);

         if(fRunOfficialTwoBody){
           fPairTranverseMass_TwoBody[i] =  new TH1F(histTitlesPairTranverseMass_TwoBody[i],histTitlesPairTranverseMass_TwoBody[i], 8000, 0, 8); // ADDED BY RAFFA
           fSameEvent->Add(fPairTranverseMass_TwoBody[i]); // ADDED BY RAFFA
          }
	}
    }


    //...............................................................................................
    // Same event multiplicity dist
    fSameEventMult = new TList();
    fSameEventMult->SetOwner();
    fSameEventMult->SetName("SameEventMult");


    fSameEventTripletMultArray = new TH2F*[28];
    TString histTitlesSameMult[28] = {"sameEventDistributionMultPPPrim","sameEventDistributionMultAPAPAPrim", "sameEventDistributionMultPPP", "sameEventDistributionMultAPAPAP", "sameEventDistributionMultPPAPrim","sameEventDistributionMultAPAPPrim", "sameEventDistributionMultPPrimPrim", "sameEventDistributionMultAPAPrimAPrim", "sameEventDistributionMultPAPrimAPrim", "sameEventDistributionMultAPPrimPrim","sameEventDistributionMultPPSamePrimMixed", "sameEventDistributionMultPPrimSamePMixed", "sameEventDistributionMultAPAPSameAPrimMixed", "sameEventDistributionMultAPAPrimSameAPMixed", "sameEventDistributionMultPPSamePMixed", "sameEventDistributionMultAPAPSameAPMixed", "sameEventDistributionMultPPSameAPrimMixed", "sameEventDistributionMultPAPrimSamePMixed", "sameEventDistributionMultAPAPSamePrimMixed", "sameEventDistributionMultAPPrimSameAPMixed", "sameEventDistributionMultPPrimSamePrimMixed", "sameEventDistributionMultPrimPrimSamePMixed", "sameEventDistributionMultAPAPrimSameAPrimMixed", "sameEventDistributionMultAPrimAPrimSameAPMixed", "sameEventDistributionMultPAPrimSameAPrimMixed", "sameEventDistributionMultAPrimAPrimSamePMixed", "sameEventDistributionMultAPPrimSamePrimMixed", "sameEventDistributionMultPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 28; ++i) {
	      fSameEventTripletMultArray[i] =  new TH2F(histTitlesSameMult[i],histTitlesSameMult[i], 8000, 0, 8,26,1,27);
              if(fRunPlotMult){fSameEventMult->Add(fSameEventTripletMultArray[i]);}
	     }
    }

    //...............................................................................................
    // Same event Pt Primaries
    fSameEventPtPrimaries = new TList();
    fSameEventPtPrimaries->SetOwner();
    fSameEventPtPrimaries->SetName("SameEventPtPrimaries");


    fSameEventTripletPtPrimaries = new TH1F*[28];
    TString histTitlesSamePtPrimaries[28] = {"sameEventDistributionPtPrimariesPPPrim","sameEventDistributionPtPrimariesAPAPAPrim", "sameEventDistributionPtPrimariesPPP", "sameEventDistributionPtPrimariesAPAPAP", "sameEventDistributionPtPrimariesPPAPrim","sameEventDistributionPtPrimariesAPAPPrim", "sameEventDistributionPtPrimariesPPrimPrim", "sameEventDistributionPtPrimariesAPAPrimAPrim", "sameEventDistributionPtPrimariesPAPrimAPrim", "sameEventDistributionPtPrimariesAPPrimPrim","sameEventDistributionPtPrimariesPPSamePrimMixed", "sameEventDistributionPtPrimariesPPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPSameAPrimMixed", "sameEventDistributionPtPrimariesAPAPrimSameAPMixed", "sameEventDistributionPtPrimariesPPSamePMixed", "sameEventDistributionPtPrimariesAPAPSameAPMixed", "sameEventDistributionPtPrimariesPPSameAPrimMixed", "sameEventDistributionPtPrimariesPAPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPSamePrimMixed", "sameEventDistributionPtPrimariesAPPrimSameAPMixed", "sameEventDistributionPtPrimariesPPrimSamePrimMixed", "sameEventDistributionPtPrimariesPrimPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPrimSameAPrimMixed", "sameEventDistributionPtPrimariesAPrimAPrimSameAPMixed", "sameEventDistributionPtPrimariesPAPrimSameAPrimMixed", "sameEventDistributionPtPrimariesAPrimAPrimSamePMixed", "sameEventDistributionPtPrimariesAPPrimSamePrimMixed", "sameEventDistributionPtPrimariesPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 28; ++i) {
        fSameEventTripletPtPrimaries[i] =  new TH1F(histTitlesSamePtPrimaries[i],histTitlesSamePtPrimaries[i], 4000, 0, 4);
              if(fRunPlotPt){fSameEventPtPrimaries->Add(fSameEventTripletPtPrimaries[i]);}
       }
    }

    fSameEventTripletPtvsQ3Primaries = new TH2F*[10];
    TString histTitlesSamePtvsQ3Primaries[10] = {"sameEventDistributionPtvsQ3PrimariesPPPrim","sameEventDistributionPtvsQ3PrimariesAPAPAPrim", "sameEventDistributionPtvsQ3PrimariesPPP", "sameEventDistributionPtvsQ3PrimariesAPAPAP", "sameEventDistributionPtvsQ3PrimariesPPAPrim","sameEventDistributionPtvsQ3PrimariesAPAPPrim", "sameEventDistributionPtvsQ3PrimariesPPrimPrim", "sameEventDistributionPtvsQ3PrimariesAPAPrimAPrim", "sameEventDistributionPtvsQ3PrimariesPAPrimAPrim", "sameEventDistributionPtvsQ3PrimariesAPPrimPrim"};//,"sameEventDistributionPtvsQ3PrimariesPPSamePrimMixed", "sameEventDistributionPtvsQ3PrimariesPPrimSamePMixed", "sameEventDistributionPtvsQ3PrimariesAPAPSameAPrimMixed", "sameEventDistributionPtvsQ3PrimariesAPAPrimSameAPMixed", "sameEventDistributionPtvsQ3PrimariesPPSamePMixed", "sameEventDistributionPtvsQ3PrimariesAPAPSameAPMixed", "sameEventDistributionPtvsQ3PrimariesPPSameAPrimMixed", "sameEventDistributionPtvsQ3PrimariesPAPrimSamePMixed", "sameEventDistributionPtvsQ3PrimariesAPAPSamePrimMixed", "sameEventDistributionPtvsQ3PrimariesAPPrimSameAPMixed", "sameEventDistributionPtvsQ3PrimariesPPrimSamePrimMixed", "sameEventDistributionPtvsQ3PrimariesPrimPrimSamePMixed", "sameEventDistributionPtvsQ3PrimariesAPAPrimSameAPrimMixed", "sameEventDistributionPtvsQ3PrimariesAPrimAPrimSameAPMixed", "sameEventDistributionPtvsQ3PrimariesPAPrimSameAPrimMixed", "sameEventDistributionPtvsQ3PrimariesAPrimAPrimSamePMixed", "sameEventDistributionPtvsQ3PrimariesAPPrimSamePrimMixed", "sameEventDistributionPtvsQ3PrimariesPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fSameEventTripletPtvsQ3Primaries[i] =  new TH2F(histTitlesSamePtvsQ3Primaries[i],histTitlesSamePtvsQ3Primaries[i],120,0, 1.2, 300, 0, 3);
              if(fRunPlotPt){fSameEventPtPrimaries->Add(fSameEventTripletPtvsQ3Primaries[i]);}
       }
    }

    //...............................................................................................
    // Same event Pt Protons
    fSameEventPtProtons = new TList();
    fSameEventPtProtons->SetOwner();
    fSameEventPtProtons->SetName("SameEventPtProtons");


    fSameEventTripletPtProtons = new TH1F*[28];
    TString histTitlesSamePtProtons[28] = {"sameEventDistributionPtProtonsPPPrim","sameEventDistributionPtProtonsAPAPAPrim", "sameEventDistributionPtProtonsPPP", "sameEventDistributionPtProtonsAPAPAP", "sameEventDistributionPtProtonsPPAPrim","sameEventDistributionPtProtonsAPAPPrim", "sameEventDistributionPtProtonsPPrimPrim", "sameEventDistributionPtProtonsAPAPrimAPrim", "sameEventDistributionPtProtonsPAPrimAPrim", "sameEventDistributionPtProtonsAPPrimPrim","sameEventDistributionPtProtonsPPSamePrimMixed", "sameEventDistributionPtProtonsPPrimSamePMixed", "sameEventDistributionPtProtonsAPAPSameAPrimMixed", "sameEventDistributionPtProtonsAPAPrimSameAPMixed", "sameEventDistributionPtProtonsPPSamePMixed", "sameEventDistributionPtProtonsAPAPSameAPMixed", "sameEventDistributionPtProtonsPPSameAPrimMixed", "sameEventDistributionPtProtonsPAPrimSamePMixed", "sameEventDistributionPtProtonsAPAPSamePrimMixed", "sameEventDistributionPtProtonsAPPrimSameAPMixed", "sameEventDistributionPtProtonsPPrimSamePrimMixed", "sameEventDistributionPtProtonsPrimPrimSamePMixed", "sameEventDistributionPtProtonsAPAPrimSameAPrimMixed", "sameEventDistributionPtProtonsAPrimAPrimSameAPMixed", "sameEventDistributionPtProtonsPAPrimSameAPrimMixed", "sameEventDistributionPtProtonsAPrimAPrimSamePMixed", "sameEventDistributionPtProtonsAPPrimSamePrimMixed", "sameEventDistributionPtProtonsPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 28; ++i) {
        fSameEventTripletPtProtons[i] =  new TH1F(histTitlesSamePtProtons[i],histTitlesSamePtProtons[i], 4000, 0, 4);
              if(fRunPlotPt){fSameEventPtProtons->Add(fSameEventTripletPtProtons[i]);}
       }
    }

    fSameEventTripletPtvsQ3Protons = new TH2F*[10];
    TString histTitlesSamePtvsQ3Protons[10] = {"sameEventDistributionPtvsQ3ProtonsPPPrim","sameEventDistributionPtvsQ3ProtonsAPAPAPrim", "sameEventDistributionPtvsQ3ProtonsPPP", "sameEventDistributionPtvsQ3ProtonsAPAPAP", "sameEventDistributionPtvsQ3ProtonsPPAPrim","sameEventDistributionPtvsQ3ProtonsAPAPPrim", "sameEventDistributionPtvsQ3ProtonsPPrimPrim", "sameEventDistributionPtvsQ3ProtonsAPAPrimAPrim", "sameEventDistributionPtvsQ3ProtonsPAPrimAPrim", "sameEventDistributionPtvsQ3ProtonsAPPrimPrim"};//,"sameEventDistributionPtvsQ3ProtonsPPSamePrimMixed", "sameEventDistributionPtvsQ3ProtonsPPrimSamePMixed", "sameEventDistributionPtvsQ3ProtonsAPAPSameAPrimMixed", "sameEventDistributionPtvsQ3ProtonsAPAPrimSameAPMixed", "sameEventDistributionPtvsQ3ProtonsPPSamePMixed", "sameEventDistributionPtvsQ3ProtonsAPAPSameAPMixed", "sameEventDistributionPtvsQ3ProtonsPPSameAPrimMixed", "sameEventDistributionPtvsQ3ProtonsPAPrimSamePMixed", "sameEventDistributionPtvsQ3ProtonsAPAPSamePrimMixed", "sameEventDistributionPtvsQ3ProtonsAPPrimSameAPMixed", "sameEventDistributionPtvsQ3ProtonsPPrimSamePrimMixed", "sameEventDistributionPtvsQ3ProtonsPrimPrimSamePMixed", "sameEventDistributionPtvsQ3ProtonsAPAPrimSameAPrimMixed", "sameEventDistributionPtvsQ3ProtonsAPrimAPrimSameAPMixed", "sameEventDistributionPtvsQ3ProtonsPAPrimSameAPrimMixed", "sameEventDistributionPtvsQ3ProtonsAPrimAPrimSamePMixed", "sameEventDistributionPtvsQ3ProtonsAPPrimSamePrimMixed", "sameEventDistributionPtvsQ3ProtonsPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fSameEventTripletPtvsQ3Protons[i] =  new TH2F(histTitlesSamePtvsQ3Protons[i],histTitlesSamePtvsQ3Protons[i], 120, 0, 1.2, 400, 0, 4);
             if(fRunPlotPt){fSameEventPtProtons->Add(fSameEventTripletPtvsQ3Protons[i]);}
       }
    }

    //...............................................................................................
    // Same event Pt Primaries L1520
    fSameEventPtPrimaries2 = new TList();
    fSameEventPtPrimaries2->SetOwner();
    fSameEventPtPrimaries2->SetName("SameEventPtPrimariesL1520");


    fSameEventTripletPtPrimaries2 = new TH1F*[28];
    TString histTitlesSamePtPrimaries2[28] = {"sameEventDistributionPtPrimariesPPPrim","sameEventDistributionPtPrimariesAPAPAPrim", "sameEventDistributionPtPrimariesPPP", "sameEventDistributionPtPrimariesAPAPAP", "sameEventDistributionPtPrimariesPPAPrim","sameEventDistributionPtPrimariesAPAPPrim", "sameEventDistributionPtPrimariesPPrimPrim", "sameEventDistributionPtPrimariesAPAPrimAPrim", "sameEventDistributionPtPrimariesPAPrimAPrim", "sameEventDistributionPtPrimariesAPPrimPrim","sameEventDistributionPtPrimariesPPSamePrimMixed", "sameEventDistributionPtPrimariesPPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPSameAPrimMixed", "sameEventDistributionPtPrimariesAPAPrimSameAPMixed", "sameEventDistributionPtPrimariesPPSamePMixed", "sameEventDistributionPtPrimariesAPAPSameAPMixed", "sameEventDistributionPtPrimariesPPSameAPrimMixed", "sameEventDistributionPtPrimariesPAPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPSamePrimMixed", "sameEventDistributionPtPrimariesAPPrimSameAPMixed", "sameEventDistributionPtPrimariesPPrimSamePrimMixed", "sameEventDistributionPtPrimariesPrimPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPrimSameAPrimMixed", "sameEventDistributionPtPrimariesAPrimAPrimSameAPMixed", "sameEventDistributionPtPrimariesPAPrimSameAPrimMixed", "sameEventDistributionPtPrimariesAPrimAPrimSamePMixed", "sameEventDistributionPtPrimariesAPPrimSamePrimMixed", "sameEventDistributionPtPrimariesPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 28; ++i) {
        fSameEventTripletPtPrimaries2[i] =  new TH1F(histTitlesSamePtPrimaries2[i],histTitlesSamePtPrimaries2[i], 4000, 0, 4);
              if(fRunPlotPt){fSameEventPtPrimaries2->Add(fSameEventTripletPtPrimaries2[i]);}
       }
    }

    //...............................................................................................
    // Same event Pt Protons L1520
    fSameEventPtProtons2 = new TList();
    fSameEventPtProtons2->SetOwner();
    fSameEventPtProtons2->SetName("SameEventPtProtonsL1520");


    fSameEventTripletPtProtons2 = new TH1F*[28];
    TString histTitlesSamePtProtons2[28] = {"sameEventDistributionPtProtonsPPPrim","sameEventDistributionPtProtonsAPAPAPrim", "sameEventDistributionPtProtonsPPP", "sameEventDistributionPtProtonsAPAPAP", "sameEventDistributionPtProtonsPPAPrim","sameEventDistributionPtProtonsAPAPPrim", "sameEventDistributionPtProtonsPPrimPrim", "sameEventDistributionPtProtonsAPAPrimAPrim", "sameEventDistributionPtProtonsPAPrimAPrim", "sameEventDistributionPtProtonsAPPrimPrim","sameEventDistributionPtProtonsPPSamePrimMixed", "sameEventDistributionPtProtonsPPrimSamePMixed", "sameEventDistributionPtProtonsAPAPSameAPrimMixed", "sameEventDistributionPtProtonsAPAPrimSameAPMixed", "sameEventDistributionPtProtonsPPSamePMixed", "sameEventDistributionPtProtonsAPAPSameAPMixed", "sameEventDistributionPtProtonsPPSameAPrimMixed", "sameEventDistributionPtProtonsPAPrimSamePMixed", "sameEventDistributionPtProtonsAPAPSamePrimMixed", "sameEventDistributionPtProtonsAPPrimSameAPMixed", "sameEventDistributionPtProtonsPPrimSamePrimMixed", "sameEventDistributionPtProtonsPrimPrimSamePMixed", "sameEventDistributionPtProtonsAPAPrimSameAPrimMixed", "sameEventDistributionPtProtonsAPrimAPrimSameAPMixed", "sameEventDistributionPtProtonsPAPrimSameAPrimMixed", "sameEventDistributionPtProtonsAPrimAPrimSamePMixed", "sameEventDistributionPtProtonsAPPrimSamePrimMixed", "sameEventDistributionPtProtonsPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 28; ++i) {
        fSameEventTripletPtProtons2[i] =  new TH1F(histTitlesSamePtProtons2[i],histTitlesSamePtProtons2[i], 4000, 0, 4);
              if(fRunPlotPt){fSameEventPtProtons2->Add(fSameEventTripletPtProtons2[i]);}
       }
    }




    fSameEventTripletMultArray_TwoBody = new TH2F*[7];
    TString histTitlesSameMult_TwoBody[7] = {"sameEventDistributionMultPP","sameEventDistributionMultAPAP", "sameEventDistributionMultPPrim", "sameEventDistributionMultAPAPrim", "sameEventDistributionMultPAPrim", "sameEventDistributionMultAPPrim", "sameEventDistributionMultPrimAPrim"};

        fPairTranverseMassVSkstar_TwoBody = new TH2F*[7]; 
        TString histTitlesPairTranverseMassVSkstar_TwoBody[7] = {"PairTranverseMassVSkstarDistributionPP","PairTranverseMassVSkstarDistributionAPAP", "PairTranverseMassVSkstarDistributionPPrim", "PairTranverseMassVSkstarDistributionAPAPrim", "PairTranverseMassVSkstarDistributionPAPrim", "PairTranverseMassVSkstarDistributionAPPrim", "PairTranverseMassVSkstarDistributionPrimAPrim"}; 

    if(!fDoOnlyThreeBody){
	for (int i = 0; i < 7; ++i) {
	   fSameEventTripletMultArray_TwoBody[i] =  new TH2F(histTitlesSameMult_TwoBody[i],histTitlesSameMult_TwoBody[i],8000, 0, 8,26,1,27);
           if(fRunPlotMult){fSameEventMult->Add(fSameEventTripletMultArray_TwoBody[i]);}
	   if(fRunOfficialTwoBody){
             fPairTranverseMassVSkstar_TwoBody[i] =  new TH2F(histTitlesPairTranverseMassVSkstar_TwoBody[i],histTitlesPairTranverseMassVSkstar_TwoBody[i], 800, 0, 8, 800, 0, 8); 
             fSameEvent->Add(fPairTranverseMassVSkstar_TwoBody[i]); 
           }
        }
    }

    //...............................................................................................
    // Mixed event 1D
    fMixedEvent = new TList();
    fMixedEvent->SetOwner();
    fMixedEvent->SetName("MixedEvent");

    fMixedEventTripletArray = new TH1F*[10];
    TString histTitlesMixed[10] = {"mixedEventDistributionPPPrim","mixedEventDistributionAPAPAPrim","mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPPAPrim","mixedEventDistributionAPAPPrim", "mixedEventDistributionPPrimPrim", "mixedEventDistributionAPAPrimAPrim", "mixedEventDistributionPAPrimAPrim", "mixedEventDistributionAPPrimPrim"};

    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 10; ++i) {
	      fMixedEventTripletArray[i] = new TH1F(histTitlesMixed[i],histTitlesMixed[i], 8000, 0, 8);
	      fMixedEvent->Add(fMixedEventTripletArray[i]);
	     }
    }

    fMixedEventTripletArray_TwoBody = new TH1F*[7];
    TString histTitlesMixed_TwoBody[7] = {"mixedEventDistributionPP","mixedEventDistributionAPAP", "mixedEventDistributionPPrim", "mixedEventDistributionAPAPrim", "mixedEventDistributionPAPrim", "mixedEventDistributionAPPrim", "mixedEventDistributionPrimAPrim"};

    if(!fDoOnlyThreeBody){
	     for (int i = 0; i < 7; ++i) {
	      fMixedEventTripletArray_TwoBody[i] =  new TH1F(histTitlesMixed_TwoBody[i],histTitlesMixed_TwoBody[i],  8000, 0, 8);
	      fMixedEvent->Add(fMixedEventTripletArray_TwoBody[i]);
	     }
    }


    //...............................................................................................
    // Mixed event multiplicity dist 
    fMixedEventMult = new TList();
    fMixedEventMult->SetOwner();
    fMixedEventMult->SetName("MixedEventMult");

    fMixedEventTripletMultArray = new TH2F*[10];
    TString histTitlesMixedMult[10] = {"mixedEventDistributionMultPPPrim","mixedEventDistributionMultAPAPAPrim","mixedEventDistributionMultPPP", "mixedEventDistributionMultAPAPAP", "mixedEventDistributionMultPPAPrim","mixedEventDistributionMultAPAPPrim", "mixedEventDistributionMultPPrimPrim", "mixedEventDistributionMultAPAPrimAPrim", "mixedEventDistributionMultPAPrimAPrim", "mixedEventDistributionMultAPPrimPrim"};

    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 10; ++i) {
	      fMixedEventTripletMultArray[i] = new TH2F(histTitlesMixedMult[i],histTitlesMixedMult[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fMixedEventMult->Add(fMixedEventTripletMultArray[i]);}
	     }
    }


    fMixedEventPtPrimaries = new TList();
    fMixedEventPtPrimaries->SetOwner();
    fMixedEventPtPrimaries->SetName("MixedEventPtPrimaries");

    fMixedEventTripletPtPrimaries = new TH1F*[10];
    TString histTitlesMixedPtPrimaries[10] = {"mixedEventDistributionPtPrimariesPPPrim","mixedEventDistributionPtPrimariesAPAPAPrim","mixedEventDistributionPtPrimariesPPP", "mixedEventDistributionPtPrimariesAPAPAP", "mixedEventDistributionPtPrimariesPPAPrim","mixedEventDistributionPtPrimariesAPAPPrim", "mixedEventDistributionPtPrimariesPPrimPrim", "mixedEventDistributionPtPrimariesAPAPrimAPrim", "mixedEventDistributionPtPrimariesPAPrimAPrim", "mixedEventDistributionPtPrimariesAPPrimPrim"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fMixedEventTripletPtPrimaries[i] = new TH1F(histTitlesMixedPtPrimaries[i],histTitlesMixedPtPrimaries[i], 4000, 0, 4);
          if(fRunPlotPt){fMixedEventPtPrimaries->Add(fMixedEventTripletPtPrimaries[i]);}
       }
    }


    fMixedEventTripletPtvsQ3Primaries = new TH2F*[10];
    TString histTitlesMixedPtvsQ3Primaries[10] = {"mixedEventDistributionPtvsQ3PrimariesPPPrim","mixedEventDistributionPtvsQ3PrimariesAPAPAPrim","mixedEventDistributionPtvsQ3PrimariesPPP", "mixedEventDistributionPtvsQ3PrimariesAPAPAP", "mixedEventDistributionPtvsQ3PrimariesPPAPrim","mixedEventDistributionPtvsQ3PrimariesAPAPPrim", "mixedEventDistributionPtvsQ3PrimariesPPrimPrim", "mixedEventDistributionPtvsQ3PrimariesAPAPrimAPrim", "mixedEventDistributionPtvsQ3PrimariesPAPrimAPrim", "mixedEventDistributionPtvsQ3PrimariesAPPrimPrim"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fMixedEventTripletPtvsQ3Primaries[i] = new TH2F(histTitlesMixedPtvsQ3Primaries[i],histTitlesMixedPtvsQ3Primaries[i], 120, 0, 1.2, 300, 0, 3);
          if(fRunPlotPt){fMixedEventPtPrimaries->Add(fMixedEventTripletPtvsQ3Primaries[i]);}
       }
    }


    fMixedEventPtProtons = new TList();
    fMixedEventPtProtons->SetOwner();
    fMixedEventPtProtons->SetName("MixedEventPtProtons");


    fMixedEventTripletPtProtons = new TH1F*[10];
    TString histTitlesMixedPtProtons[10] = {"mixedEventDistributionPtProtonsPPPrim","mixedEventDistributionPtProtonsAPAPAPrim","mixedEventDistributionPtProtonsPPP", "mixedEventDistributionPtProtonsAPAPAP", "mixedEventDistributionPtProtonsPPAPrim","mixedEventDistributionPtProtonsAPAPPrim", "mixedEventDistributionPtProtonsPPrimPrim", "mixedEventDistributionPtProtonsAPAPrimAPrim", "mixedEventDistributionPtProtonsPAPrimAPrim", "mixedEventDistributionPtProtonsAPPrimPrim"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fMixedEventTripletPtProtons[i] = new TH1F(histTitlesMixedPtProtons[i],histTitlesMixedPtProtons[i], 4000, 0, 4);
          if(fRunPlotPt){fMixedEventPtProtons->Add(fMixedEventTripletPtProtons[i]);}
       }
    }


    fMixedEventTripletPtvsQ3Protons = new TH2F*[10];
    TString histTitlesMixedPtvsQ3Protons[10] = {"mixedEventDistributionPtvsQ3ProtonsPPPrim","mixedEventDistributionPtvsQ3ProtonsAPAPAPrim","mixedEventDistributionPtvsQ3ProtonsPPP", "mixedEventDistributionPtvsQ3ProtonsAPAPAP", "mixedEventDistributionPtvsQ3ProtonsPPAPrim","mixedEventDistributionPtvsQ3ProtonsAPAPPrim", "mixedEventDistributionPtvsQ3ProtonsPPrimPrim", "mixedEventDistributionPtvsQ3ProtonsAPAPrimAPrim", "mixedEventDistributionPtvsQ3ProtonsPAPrimAPrim", "mixedEventDistributionPtvsQ3ProtonsAPPrimPrim"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fMixedEventTripletPtvsQ3Protons[i] = new TH2F(histTitlesMixedPtvsQ3Protons[i],histTitlesMixedPtvsQ3Protons[i], 120, 0, 1.2, 400, 0, 4);
          if(fRunPlotPt){fMixedEventPtProtons->Add(fMixedEventTripletPtvsQ3Protons[i]);}
       }
    }



    fMixedEventPtPrimaries2 = new TList();
    fMixedEventPtPrimaries2->SetOwner();
    fMixedEventPtPrimaries2->SetName("MixedEventPtPrimariesL1520");

    fMixedEventTripletPtPrimaries2 = new TH1F*[10];
    TString histTitlesMixedPtPrimaries2[10] = {"mixedEventDistributionPtPrimariesPPPrim","mixedEventDistributionPtPrimariesAPAPAPrim","mixedEventDistributionPtPrimariesPPP", "mixedEventDistributionPtPrimariesAPAPAP", "mixedEventDistributionPtPrimariesPPAPrim","mixedEventDistributionPtPrimariesAPAPPrim", "mixedEventDistributionPtPrimariesPPrimPrim", "mixedEventDistributionPtPrimariesAPAPrimAPrim", "mixedEventDistributionPtPrimariesPAPrimAPrim", "mixedEventDistributionPtPrimariesAPPrimPrim"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fMixedEventTripletPtPrimaries2[i] = new TH1F(histTitlesMixedPtPrimaries2[i],histTitlesMixedPtPrimaries2[i], 4000, 0, 4);
          if(fRunPlotPt){fMixedEventPtPrimaries2->Add(fMixedEventTripletPtPrimaries2[i]);}
       }
    }

    fMixedEventPtProtons2 = new TList();
    fMixedEventPtProtons2->SetOwner();
    fMixedEventPtProtons2->SetName("MixedEventPtProtonsL1520");


    fMixedEventTripletPtProtons2 = new TH1F*[10];
    TString histTitlesMixedPtProtons2[10] = {"mixedEventDistributionPtProtonsPPPrim","mixedEventDistributionPtProtonsAPAPAPrim","mixedEventDistributionPtProtonsPPP", "mixedEventDistributionPtProtonsAPAPAP", "mixedEventDistributionPtProtonsPPAPrim","mixedEventDistributionPtProtonsAPAPPrim", "mixedEventDistributionPtProtonsPPrimPrim", "mixedEventDistributionPtProtonsAPAPrimAPrim", "mixedEventDistributionPtProtonsPAPrimAPrim", "mixedEventDistributionPtProtonsAPPrimPrim"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fMixedEventTripletPtProtons2[i] = new TH1F(histTitlesMixedPtProtons2[i],histTitlesMixedPtProtons2[i], 4000, 0, 4);
          if(fRunPlotPt){fMixedEventPtProtons2->Add(fMixedEventTripletPtProtons2[i]);}
       }
    }



     fMixedEventTripletMultArray_TwoBody = new TH2F*[7];
    TString histTitlesMixedMult_TwoBody[7] = {"MixedEventDistributionMultPP","MixedEventDistributionMultAPAP", "MixedEventDistributionMultPPrim", "MixedEventDistributionMultAPAPrim", "MixedEventDistributionMultPAPrim", "MixedEventDistributionMultAPPrim", "MixedEventDistributionMultPrimAPrim"};

    if(!fDoOnlyThreeBody){
	     for (int i = 0; i < 7; ++i) {
	      fMixedEventTripletMultArray_TwoBody[i] =  new TH2F(histTitlesMixedMult_TwoBody[i],histTitlesMixedMult_TwoBody[i], 8000, 0, 8,26,1,27);
             if(fRunPlotMult){ fMixedEventMult->Add(fMixedEventTripletMultArray_TwoBody[i]); }
	     }
    }
    fResultsThreeBody->Add(fSameEvent);
    fResultsThreeBody->Add(fMixedEvent);
    if(fRunPlotMult){
    fResultsThreeBody->Add(fSameEventMult);
    fResultsThreeBody->Add(fMixedEventMult);
    }
    if(fRunPlotPt){
    fResultsThreeBody->Add(fSameEventPtPrimaries);
    fResultsThreeBody->Add(fMixedEventPtPrimaries);
    fResultsThreeBody->Add(fSameEventPtProtons);
    fResultsThreeBody->Add(fMixedEventPtProtons);
    fResultsThreeBody->Add(fSameEventPtPrimaries2);
    fResultsThreeBody->Add(fMixedEventPtPrimaries2);
    fResultsThreeBody->Add(fSameEventPtProtons2);
    fResultsThreeBody->Add(fMixedEventPtProtons2);
    }

    //...............................................................................................
    // Same event phi theta distribution Same Pair
    fSameEventPhiTheta_SamePair = new TList();
    fSameEventPhiTheta_SamePair->SetOwner();
    fSameEventPhiTheta_SamePair->SetName("SameEventPhiTheta");

    fSameEventTripletPhiThetaArray_SamePair = new TH2F*[56];
    TString histTitlesSamePhiEta_SamePair[28] ={"sameEventPhiEtaMultPPPrim","sameEventPhiEtaMultAPAPAPrim",
      "sameEventPhiEtaMultPPP", "sameEventPhiEtaMultAPAPAP", "sameEventPhiEtaMultPPAPrim","sameEventPhiEtaMultAPAPPrim","sameEventPhiEtaMultPPrimPrim", "sameEventPhiEtaMultAPAPrimAPrim", "sameEventPhiEtaMultPAPrimAPrim", "sameEventPhiEtaMultAPPrimPrim","sameEventPhiEtaMultPPSamePrimMixed", "sameEventPhiEtaMultPPrimSamePMixed", "sameEventPhiEtaMultAPAPSameAPrimMixed", "sameEventPhiEtaMultAPAPrimSameAPMixed", "sameEventPhiEtaMultPPSamePMixed", "sameEventPhiEtaMultAPAPSameAPMixed", "sameEventPhiEtaMultPPSameAPrimMixed", "sameEventPhiEtaMultPAPrimSamePMixed", "sameEventPhiEtaMultAPAPSamePrimMixed", "sameEventPhiEtaMultAPPrimSameAPMixed", "sameEventPhiEtaMultPPrimSamePimMixed", "sameEventPhiEtaMultPrimPrimSamePMixed", "sameEventPhiEtaMultAPAPrimSameAPMixed", "sameEventPhiEtaMultAPrimAPrimSameAPMixed", "sameEventPhiEtaMultPAPrimSameAPrimMixed", "sameEventPhiEtaMultAPrimAPrimSamePMixed", "sameEventPhiEtaMultAPPrimSamePrimMixed", "sameEventPhiEtaMultPrimPrimSameAPMixed"};


    if(fDoOnlyThreeBody){
	    for(int i=0;i<28;i++){
	      fSameEventTripletPhiThetaArray_SamePair[i] = new TH2F(histTitlesSamePhiEta_SamePair[i]+"Before",histTitlesSamePhiEta_SamePair[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fSameEventTripletPhiThetaArray_SamePair[28+i] = new TH2F(histTitlesSamePhiEta_SamePair[i]+"After",histTitlesSamePhiEta_SamePair[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);

               if(fRunPlotPhiTheta){
	               fSameEventPhiTheta_SamePair->Add(fSameEventTripletPhiThetaArray_SamePair[i]);
	               fSameEventPhiTheta_SamePair->Add(fSameEventTripletPhiThetaArray_SamePair[28+i]);
               }
	    }
    }


    fSameEventTripletPhiThetaArray_TwoBody = new TH2F*[14];
    TString histTitlesSamePhiEta_TwoBody[7] ={"sameEventPhiEtaMultPP", "sameEventPhiEtaMultAPAP", "sameEventPhiEtaMultPPrim", "sameEventPhiEtaMultAPAPrim", "sameEventPhiEtaMultPAPrim", "sameEventPhiEtaMultAPPrim", "sameEventPhiEtaMultPrimAPrim",};

    if(!fDoOnlyThreeBody){
	    for(int i=0;i<7;i++){
	      fSameEventTripletPhiThetaArray_TwoBody[i] = new TH2F(histTitlesSamePhiEta_TwoBody[i]+"Before",histTitlesSamePhiEta_TwoBody[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fSameEventTripletPhiThetaArray_TwoBody[7+i] = new TH2F(histTitlesSamePhiEta_TwoBody[i]+"After",histTitlesSamePhiEta_TwoBody[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);

              if(fRunPlotPhiTheta){
      	        fSameEventPhiTheta_SamePair->Add(fSameEventTripletPhiThetaArray_TwoBody[i]);
      	        fSameEventPhiTheta_SamePair->Add(fSameEventTripletPhiThetaArray_TwoBody[7+i]);
              }
	    }
    }

    //...............................................................................................
    // Same event phi theta distribution Different Pair
    fSameEventPhiTheta_DifferentPair = new TList();
    fSameEventPhiTheta_DifferentPair->SetOwner();
    fSameEventPhiTheta_DifferentPair->SetName("SameEventPhiTheta_DifferentPair");

    fSameEventTripletPhiThetaArray_DifferentPair = new TH2F*[56];
    TString histTitlesSamePhiEta_DifferentPair[28] ={"sameEventPhiEtaMultPPPrim","sameEventPhiEtaMultAPAPAPrim",
      "sameEventPhiEtaMultPPP", "sameEventPhiEtaMultAPAPAP", "sameEventPhiEtaMultPPAPrim","sameEventPhiEtaMultAPAPPrim","sameEventPhiEtaMultPPrimPrim", "sameEventPhiEtaMultAPAPrimAPrim", "sameEventPhiEtaMultPAPrimAPrim", "sameEventPhiEtaMultAPPrimPrim","sameEventPhiEtaMultPPSamePrimMixed", "sameEventPhiEtaMultPPrimSamePMixed", "sameEventPhiEtaMultAPAPSameAPrimMixed", "sameEventPhiEtaMultAPAPrimSameAPMixed", "sameEventPhiEtaMultPPSamePMixed", "sameEventPhiEtaMultAPAPSameAPMixed", "sameEventPhiEtaMultPPSameAPrimMixed", "sameEventPhiEtaMultPAPrimSamePMixed", "sameEventPhiEtaMultAPAPSamePrimMixed", "sameEventPhiEtaMultAPPrimSameAPMixed", "sameEventPhiEtaMultPPrimSamePimMixed", "sameEventPhiEtaMultPrimPrimSamePMixed", "sameEventPhiEtaMultAPAPrimSameAPMixed", "sameEventPhiEtaMultAPrimAPrimSameAPMixed", "sameEventPhiEtaMultPAPrimSameAPrimMixed", "sameEventPhiEtaMultAPrimAPrimSamePMixed", "sameEventPhiEtaMultAPPrimSamePrimMixed", "sameEventPhiEtaMultPrimPrimSameAPMixed"};


    if(fDoOnlyThreeBody){
      for(int i=0;i<28;i++){
        fSameEventTripletPhiThetaArray_DifferentPair[i] = new TH2F(histTitlesSamePhiEta_DifferentPair[i]+"Before",histTitlesSamePhiEta_DifferentPair[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
        fSameEventTripletPhiThetaArray_DifferentPair[28+i] = new TH2F(histTitlesSamePhiEta_DifferentPair[i]+"After",histTitlesSamePhiEta_DifferentPair[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);

        if(fRunPlotPhiTheta){
          fSameEventPhiTheta_DifferentPair->Add(fSameEventTripletPhiThetaArray_DifferentPair[i]);
          fSameEventPhiTheta_DifferentPair->Add(fSameEventTripletPhiThetaArray_DifferentPair[28+i]);
        }
      }
    }

    //...............................................................................................
    // Mixed event phi theta distribution
    fMixedEventPhiTheta_SamePair = new TList();
    fMixedEventPhiTheta_SamePair->SetOwner();
    fMixedEventPhiTheta_SamePair->SetName("MixedEventPhiTheta");

    fMixedEventTripletPhiThetaArray_SamePair = new TH2F*[20];
    TString histTitlesMixedPhiEta_SamePair[10] = {"mixedEventPhiEtaPPPrim","mixedEventPhiEtaAPAPAPrim",
      "mixedEventPhiEtaPPP", "mixedEventPhiEtaAPAPAP", "mixedEventPhiEtaPPAPrim","mixedEventPhiEtaAPAPPrim", "mixedEventPhiEtaPPrimPrim", "mixedEventPhiEtaAPAPrimAPrim", "mixedEventPhiEtaPAPrimAPrim", "mixedEventPhiEtaAPPrimPrim"};

    if(fDoOnlyThreeBody){
	    for(int i=0;i<10;i++){
	      fMixedEventTripletPhiThetaArray_SamePair[i] = new TH2F(histTitlesMixedPhiEta_SamePair[i]+"Before",histTitlesMixedPhiEta_SamePair[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fMixedEventTripletPhiThetaArray_SamePair[10+i] = new TH2F(histTitlesMixedPhiEta_SamePair[i]+"After",histTitlesMixedPhiEta_SamePair[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);

             if(fRunPlotPhiTheta){
      	        fMixedEventPhiTheta_SamePair->Add(fMixedEventTripletPhiThetaArray_SamePair[i]);
      	        fMixedEventPhiTheta_SamePair->Add(fMixedEventTripletPhiThetaArray_SamePair[10+i]);
             }
	    }
    }

    fMixedEventTripletPhiThetaArray_TwoBody = new TH2F*[14];
    TString histTitlesMixedPhiEta_TwoBody[7] ={"mixedEventPhiEtaMultPP", "mixedEventPhiEtaMultAPAP", "mixedEventPhiEtaMultPPrim", "mixedEventPhiEtaMultAPAPrim", "mixedEventPhiEtaMultPAPrim", "mixedEventPhiEtaMultAPPrim", "mixedEventPhiEtaMultPrimAPrim",};

    if(!fDoOnlyThreeBody){
	    for(int i=0;i<7;i++){
	      fMixedEventTripletPhiThetaArray_TwoBody[i] = new TH2F(histTitlesMixedPhiEta_TwoBody[i]+"Before",histTitlesMixedPhiEta_TwoBody[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fMixedEventTripletPhiThetaArray_TwoBody[7+i] = new TH2F(histTitlesMixedPhiEta_TwoBody[i]+"After",histTitlesMixedPhiEta_TwoBody[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);

              if(fRunPlotPhiTheta){
                fMixedEventPhiTheta_SamePair->Add(fMixedEventTripletPhiThetaArray_TwoBody[i]);
	              fMixedEventPhiTheta_SamePair->Add(fMixedEventTripletPhiThetaArray_TwoBody[7+i]);
              }
	    }
    }

    //...............................................................................................
    // Mixed event phi theta distribution Different Pair
    fMixedEventPhiTheta_DifferentPair = new TList();
    fMixedEventPhiTheta_DifferentPair->SetOwner();
    fMixedEventPhiTheta_DifferentPair->SetName("MixedEventPhiTheta_DifferentPair");

    fMixedEventTripletPhiThetaArray_DifferentPair = new TH2F*[20];
    TString histTitlesMixedPhiEta_DifferentPair[10] = {"mixedEventPhiEtaPPPrim","mixedEventPhiEtaAPAPAPrim",
      "mixedEventPhiEtaPPP", "mixedEventPhiEtaAPAPAP", "mixedEventPhiEtaPPAPrim","mixedEventPhiEtaAPAPPrim", "mixedEventPhiEtaPPrimPrim", "mixedEventPhiEtaAPAPrimAPrim", "mixedEventPhiEtaPAPrimAPrim", "mixedEventPhiEtaAPPrimPrim"};

    if(fDoOnlyThreeBody){
      for(int i=0;i<10;i++){
        fMixedEventTripletPhiThetaArray_DifferentPair[i] = new TH2F(histTitlesMixedPhiEta_DifferentPair[i]+"Before",histTitlesMixedPhiEta_DifferentPair[i]+"Before", 500, -0.15,0.15,500,-0.15,0.15);
        fMixedEventTripletPhiThetaArray_DifferentPair[10+i] = new TH2F(histTitlesMixedPhiEta_DifferentPair[i]+"After",histTitlesMixedPhiEta_DifferentPair[i]+"After", 500, -0.15,0.15,500,-0.15,0.15);

             if(fRunPlotPhiTheta){
                fMixedEventPhiTheta_DifferentPair->Add(fMixedEventTripletPhiThetaArray_DifferentPair[i]);
                fMixedEventPhiTheta_DifferentPair->Add(fMixedEventTripletPhiThetaArray_DifferentPair[10+i]);
             }
      }
    }


    if(fRunPlotPhiTheta){
    fResultsThreeBody->Add(fSameEventPhiTheta_SamePair);
    fResultsThreeBody->Add(fMixedEventPhiTheta_SamePair);
    fResultsThreeBody->Add(fSameEventPhiTheta_DifferentPair);
    fResultsThreeBody->Add(fMixedEventPhiTheta_DifferentPair);
    }

    //...............................................................................................
    // Kinematics plots

    fKinematicsPlots = new TList();
    fKinematicsPlots->SetOwner();
    fKinematicsPlots->SetName("KinematicsPlots");

    fKinematics = new TH2F*[16];
    TString histTitlesfKinematics[16] = {"sameEventKinematicsPPPrim","sameEventKinematicsAPAPAPrim",
        "sameEventKinematicsPPP", "sameEventKinematicsAPAPAP", "sameEventKinematicsPPAPrim","sameEventKinematicsAPAPPrim",
        "sameEventKinematicsPPSamePrimMixed", "sameEventKinematicsPPrimSamePMixed",
        "sameEventKinematicsAPAPSameAPrimMixed", "sameEventKinematicsAPAPrimSameAPMixed",
        "sameEventKinematicsPPSamePMixed", "sameEventKinematicsAPAPSameAPMixed",
        "sameEventKinematicsPPSameAPrimMixed", "sameEventKinematicsPAPrimSamePMixed",
        "sameEventKinematicsAPAPSamePrimMixed", "sameEventKinematicsAPPrimSameAPMixed"};

    fPrimAngles = new TH2F*[16];
    TString histTitlesfPrimAngles[16] = {"sameEventPrimAnglesPPPrim","sameEventPrimAnglesAPAPAPrim",
        "sameEventPrimAnglesPPP", "sameEventPrimAnglesAPAPAP", "sameEventPrimAnglesPPAPrim","sameEventPrimAnglesAPAPPrim",
        "sameEventPrimAnglesPPSamePrimMixed", "sameEventPrimAnglesPPrimSamePMixed",
        "sameEventPrimAnglesAPAPSameAPrimMixed", "sameEventPrimAnglesAPAPrimSameAPMixed",
        "sameEventPrimAnglesPPSamePMixed", "sameEventPrimAnglesAPAPSameAPMixed",
        "sameEventPrimAnglesPPSameAPrimMixed", "sameEventPrimAnglesPAPrimSamePMixed",
        "sameEventPrimAnglesAPAPSamePrimMixed", "sameEventPrimAnglesAPPrimSameAPMixed" };

    fDeta = new TH2F*[16];
    TString histTitlesfDeta[16] = {"sameEventDetaPPPrim","sameEventDetaAPAPAPrim",
        "sameEventDetaPPP", "sameEventDetaAPAPAP", "sameEventDetaPPAPrim","sameEventDetaAPAPPrim",
        "sameEventDetaPPSamePrimMixed", "sameEventDetaPPrimSamePMixed",
        "sameEventDetaAPAPSameAPrimMixed", "sameEventDetaAPAPrimSameAPMixed",
        "sameEventDetaPPSamePMixed", "sameEventDetaAPAPSameAPMixed",
        "sameEventDetaPPSameAPrimMixed", "sameEventDetaPAPrimSamePMixed",
        "sameEventDetaAPAPSamePrimMixed", "sameEventDetaAPPrimSameAPMixed"};

    fDphi = new TH2F*[16];
    TString histTitlesfDphi[16] = {"sameEventDphiPPPrim","sameEventDphiAPAPAPrim",
        "sameEventDphiPPP", "sameEventDphiAPAPAP", "sameEventDphiPPAPrim","sameEventDphiAPAPPrim",
        "sameEventDphiPPSamePrimMixed", "sameEventDphiPPrimSamePMixed",
        "sameEventDphiAPAPSameAPrimMixed", "sameEventDphiAPAPrimSameAPMixed",
        "sameEventDphiPPSamePMixed", "sameEventDphiAPAPSameAPMixed",
        "sameEventDphiPPSameAPrimMixed", "sameEventDphiPAPrimSamePMixed",
        "sameEventDphiAPAPSamePrimMixed", "sameEventDphiAPPrimSameAPMixed"};

    fKinematicsME = new TH2F*[6];
    TString histTitlesfKinematicsME[6] = {"mixedEventKinematicsPPPrim","mixedEventKinematicsAPAPAPrim",
        "mixedEventKinematicsPPP", "mixedEventKinematicsAPAPAP", "mixedEventKinematicsPPAPrim","mixedEventKinematicsAPAPPrim" };

    fPrimAnglesME = new TH2F*[6];
    TString histTitlesfPrimAnglesME[6] = {"mixedEventPrimAnglesPPPrim","mixedEventPrimAnglesAPAPAPrim",
        "mixedEventPrimAnglesPPP", "mixedEventPrimAnglesAPAPAP", "mixedEventPrimAnglesPPAPrim","mixedEventPrimAnglesAPAPPrim" };

    fDetaME = new TH2F*[6];
    TString histTitlesfDetaME[6] = {"mixedEventDetaPPPrim","mixedEventDetaAPAPAPrim",
        "mixedEventDetaPPP", "mixedEventDetaAPAPAP", "mixedEventDetaPPAPrim","mixedEventDetaAPAPPrim" };

    fDphiME = new TH2F*[6];
    TString histTitlesfDphiME[6] = {"mixedEventDphiPPPrim","mixedEventDphiAPAPAPrim",
        "mixedEventDphiPPP", "mixedEventDphiAPAPAP", "mixedEventDphiPPAPrim","mixedEventDphiAPAPPrim" };

    if(fDoKinematicsPlots){

      for(int i=0;i<16;i++){
          fKinematics[i] = new TH2F(histTitlesfKinematics[i],histTitlesfKinematics[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fKinematics[i]);

          fPrimAngles[i] = new TH2F(histTitlesfPrimAngles[i],histTitlesfPrimAngles[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fPrimAngles[i]);

          fDeta[i] = new TH2F(histTitlesfDeta[i],histTitlesfDeta[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fDeta[i]);

          fDphi[i] = new TH2F(histTitlesfDphi[i],histTitlesfDphi[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fDphi[i]);

      }
      for(int i=0;i<6;i++){
          fKinematicsME[i] = new TH2F(histTitlesfKinematicsME[i],histTitlesfKinematicsME[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fKinematicsME[i]);

          fPrimAnglesME[i] = new TH2F(histTitlesfPrimAnglesME[i],histTitlesfPrimAnglesME[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fPrimAnglesME[i]);

          fDetaME[i] = new TH2F(histTitlesfDetaME[i],histTitlesfDetaME[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fDetaME[i]);

          fDphiME[i] = new TH2F(histTitlesfDphiME[i],histTitlesfDphiME[i], 360, 0., 2.*TMath::Pi(),360, 0., 2.*TMath::Pi());
          fKinematicsPlots->Add(fDphiME[i]);

      }


      fpTvsEtaTrueKaons = new TH2F("fpTvsEtaTrueKaons","fpTvsEtaTrueKaons",5000,0, 5, 200, -10,10);
      fpTvsEtaTrueAntiKaons = new TH2F("fpTvsEtaTrueAntiKaons","fpTvsEtaTrueAntiKaons",5000,0, 5, 200, -10,10);
      fpTvsEtaTrueProtons = new TH2F("fpTvsEtaTrueProtons","fpTvsEtaTrueProtons",5000,0, 5, 200, -10,10);
      fpTvsEtaTrueAntiProtons = new TH2F("fpTvsEtaTrueAntiProtons","fpTvsEtaTrueAntiProtons",5000,0, 5, 200, -10,10);
      fpTvsEtaRecoKaons = new TH2F("fpTvsEtaRecoKaons","fpTvsEtaRecoKaons",5000,0, 5, 200, -10,10);
      fpTvsEtaRecoAntiKaons = new TH2F("fpTvsEtaRecoAntiKaons","fpTvsEtaRecoAntiKaons",5000,0, 5, 200, -10,10);
      fpTvsEtaRecoProtons = new TH2F("fpTvsEtaRecoProtons","fpTvsEtaRecoProtons",5000,0, 5, 200, -10,10);
      fpTvsEtaRecoAntiProtons = new TH2F("fpTvsEtaRecoAntiProtons","fpTvsEtaRecoAntiProtons",5000,0, 5, 200, -10,10);
  
     if(fPlotsMC){

      fKinematicsPlots->Add(fpTvsEtaTrueKaons);
      fKinematicsPlots->Add(fpTvsEtaTrueAntiKaons);
      fKinematicsPlots->Add(fpTvsEtaTrueProtons);
      fKinematicsPlots->Add(fpTvsEtaTrueAntiProtons);
      fKinematicsPlots->Add(fpTvsEtaRecoKaons);
      fKinematicsPlots->Add(fpTvsEtaRecoAntiKaons);
      fKinematicsPlots->Add(fpTvsEtaRecoProtons);
      fKinematicsPlots->Add(fpTvsEtaRecoAntiProtons);
      }

      
      fResultsThreeBody->Add(fKinematicsPlots);
    }

    //...............................................................................................
    // Q3 vs q12 plot for theory
      fQ3Vskq12 = new TList();
      fQ3Vskq12->SetOwner();
      fQ3Vskq12->SetName("Q3Vskq12");

      fQ3VskDistributionsArrayq12 =  new TH2F*[10]; 
      TString histTitlesfQ3VskDistributions[10] =  {"Q3vskDistributionPPSamePrimMixed", "Q3vskDistributionPPrimSamePMixed",
        "Q3vskDistributionAPAPSameAPrimMixed","Q3vskDistributionAPAPrimSameaAPMixed", "Q3vskDistributionPPSamePMixed",
        "Q3vskDistributionAPAPSameAPMixed", "Q3vskDistributionPPSameAPrimMixed", "Q3vskDistributionPAPrimSamePMixed",
        "Q3vskDistributionAPAPSameaPrimMixed", "Q3vskDistributionAPPrimSameaAPMixed"};

        //,"Q3vskDistributionLLSameLMixed", "Q3vskDistributionaLaLSameaLMixed", "TRASH"};
    
      for(int i=0;i<10;i++){
        fQ3VskDistributionsArrayq12[i] = new TH2F(histTitlesfQ3VskDistributions[i],histTitlesfQ3VskDistributions[i], 2000, 0., 2.,2000,0.,2.);
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
        fQ3VskDistributionsArrayq12Mixed[i] = new TH2F(histTitlesfQ3VskDistributionsMixed[i],histTitlesfQ3VskDistributionsMixed[i], 2000, 0., 2.,2000,0.,2.);
        fQ3Vskq12Mixed->Add(fQ3VskDistributionsArrayq12Mixed[i]);
      }

      // Q3 vs q23 plot for theory
      fQ3Vskq23 = new TList();
      fQ3Vskq23->SetOwner();
      fQ3Vskq23->SetName("Q3Vskq23");

      fQ3VskDistributionsArrayq23 =  new TH2F*[10];
      for(int i=0;i<10;i++){
        fQ3VskDistributionsArrayq23[i] = new TH2F(histTitlesfQ3VskDistributions[i],histTitlesfQ3VskDistributions[i], 2000, 0., 2.,2000,0.,2.);
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
        fQ3VskDistributionsArrayq23Mixed[i] = new TH2F(histTitlesfQ3VskDistributionsMixed2[i],histTitlesfQ3VskDistributionsMixed2[i], 2000, 0., 2.,2000,0.,2.);
        fQ3Vskq23Mixed->Add(fQ3VskDistributionsArrayq23Mixed[i]);
      }

    if(fRunPlotQ3Vsq){

      fResultsThreeBody->Add(fQ3Vskq12);
      fResultsThreeBody->Add(fQ3Vskq12Mixed);
      fResultsThreeBody->Add(fQ3Vskq23);
      fResultsThreeBody->Add(fQ3Vskq23Mixed);
    }


    //...............................................................................................
    //Invariant Mass
    fInvMass = new TH2F*[16];
    TString histTitlesInvMass[16] = {"InvMassPPPrim","InvMassAPAPAPrim", "InvMassPPP", "InvMassAPAPAP", "InvMassPPAPrim","InvMassAPAPPrim", "InvMassPPSamePrimMixed", "InvMassPPrimSamePMixed", "InvMassAPAPSameAPrimMixed", "InvMassAPAPrimSameAPMixed", "InvMassPPSamePMixed", "InvMassAPAPSameAPMixed", "InvMassPPSameAPrimMixed", "InvMassPAPrimSamePMixed", "InvMassAPAPSamePrimMixed", "InvMassAPPrimSameAPMixed"};


    if (fRunPlotInvMass){
      fInvMassList = new TList();
      fInvMassList->SetOwner();
      fInvMassList->SetName("InvMass");


      for(int i=0;i<16;i++){
        fInvMass[i] = new TH2F(histTitlesInvMass[i],histTitlesInvMass[i], 500, 0., 5., 12, 0., 1.2);
        fInvMassList->Add(fInvMass[i]);
      }

      fResultsThreeBody->Add(fInvMassList);
     
    }
   }//if (fRunThreeBody)


  //...............................................................................................
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
/*
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
*/
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

      if(fStandardMixing){
         fPartContainer.push_back(MultContainer);
      } else {
        if(fDoOnlyThreeBody){
          fPartContainerPPP.push_back(MultContainer);
          fPartContainerPPPrim.push_back(MultContainer);
          fPartContainerPPAPrim.push_back(MultContainer);
        } else {
          fPartContainerPP.push_back(MultContainer);
          fPartContainerPPrim.push_back(MultContainer);
          fPartContainerPAPrim.push_back(MultContainer);
        }
      }
    }



    if(fStandardMixing){

        fVectPartContainers.push_back(fPartContainer);

    }else{

      if(fDoOnlyThreeBody){

        fVectPartContainers.push_back(fPartContainerPPP);
        fVectPartContainers.push_back(fPartContainerPPPrim);
        fVectPartContainers.push_back(fPartContainerPPAPrim);

      } else {
      
        fVectPartContainers.push_back(fPartContainerPP);
        fVectPartContainers.push_back(fPartContainerPPrim);
        fVectPartContainers.push_back(fPartContainerPAPrim);

      }
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
  std::vector<AliFemtoDreamBasePart> Primaries;
  std::vector<AliFemtoDreamBasePart> AntiPrimaries;



  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) { // TO DO: think about double track selection
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
      //if (fPlotsMC) fpTvsEtaRecoProtons->Fill(track->Pt(),track->Eta());
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
      //if (fPlotsMC) fpTvsEtaRecoAntiProtons->Fill(track->Pt(),track->Eta());
    }
    if (fPrimary->isSelected(fTrack)) {
      Primaries.push_back(*fTrack);
      //if (fPlotsMC) fpTvsEtaRecoKaons->Fill(track->Pt(),track->Eta());
    }
    if (fAntiPrimary->isSelected(fTrack)) {
      AntiPrimaries.push_back(*fTrack);
      //if (fPlotsMC) fpTvsEtaRecoAntiKaons->Fill(track->Pt(),track->Eta());
    }
  }



  if (fPlotsMC) {
    AliAODInputHandler *eventHandler =
      dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()
                                        ->GetInputEventHandler());
    AliMCEvent* fMC = eventHandler->MCEvent();

    for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
      if (mcPart->IsPhysicalPrimary()) {
        if (mcPart->GetPdgCode() == fProton->GetPDGCode()) {
          fpTvsEtaTrueProtons ->Fill(mcPart->Pt(),mcPart->Eta());
        } else if (mcPart->GetPdgCode() == fAntiProton->GetPDGCode()) {
          fpTvsEtaTrueAntiProtons ->Fill(mcPart->Pt(),mcPart->Eta());
        } else if (mcPart->GetPdgCode() == fPrimary->GetPDGCode()) {
          fpTvsEtaTrueKaons ->Fill(mcPart->Pt(),mcPart->Eta());
        } else if (mcPart->GetPdgCode() == fAntiPrimary->GetPDGCode()) {
          fpTvsEtaTrueAntiKaons ->Fill(mcPart->Pt(),mcPart->Eta());
        }
      }
    }
  }


  //b) Optional pair cleaning +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 fPairCleaner->ResetArray();
/*
  if(fCleanWithLambdas)
  {
     fPairCleaner->CleanTrackAndDecay(&Protons, &Lambdas, 0);
     fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiLambdas, 1);

     fPairCleaner->CleanDecay(&Lambdas, 0);
     fPairCleaner->CleanDecay(&AntiLambdas, 1);
  }
*/
  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Primaries);
  fPairCleaner->StoreParticle(AntiPrimaries);


  int ContainerIdPPP;
  int ContainerIdPPPrim;
  int ContainerIdPPAPrim;
  int ContainerIdPP;
  int ContainerIdPPrim;
  int ContainerIdPAPrim;


  if(fStandardMixing){

    ContainerIdPPP = 0;
    ContainerIdPPPrim = 0;
    ContainerIdPPAPrim = 0;
    ContainerIdPP = 0;
    ContainerIdPPrim = 0;
    ContainerIdPAPrim = 0;

  }else{

    ContainerIdPPP = 0;
    ContainerIdPPPrim = 1;
    ContainerIdPPAPrim = 2;
    ContainerIdPP = 0;
    ContainerIdPPrim = 1;
    ContainerIdPAPrim = 2;

  }


  //c) Start Three Body Calculus +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(fRunThreeBody){
    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    int bins[2] = { 0, 0 };
    float ZVtx = fEvent->GetZVertex();
    float Mult = fEvent->GetMultiplicity();
    fPartColl->FindBin(ZVtx, Mult, bins);

    // c.0) Same event distribution --------------------------------
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles();

    if(fDoOnlyThreeBody){

       // Proton Proton Primary
       FillTripletDistribution( ParticleVector, 2, 0, 0, fSameEventTripletArray[0],PDGCodes, bins[1],fSameEventTripletMultArray[0], fSameEventTripletPtPrimaries[0], fSameEventTripletPtProtons[0], fSameEventTripletPtPrimaries2[0], fSameEventTripletPtProtons2[0], fSameEventTripletPtvsQ3Primaries[0], fSameEventTripletPtvsQ3Protons[0], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,0, *fConfig, fKinematics[0], fPrimAngles[0], fDeta[0], fDphi[0], fInvMass[0]);
       // Antiproton Antiproton Antiprimary
       FillTripletDistribution( ParticleVector, 3, 1, 1, fSameEventTripletArray[1],PDGCodes, bins[1],fSameEventTripletMultArray[1], fSameEventTripletPtPrimaries[1], fSameEventTripletPtProtons[1], fSameEventTripletPtPrimaries2[1], fSameEventTripletPtProtons2[1], fSameEventTripletPtvsQ3Primaries[1], fSameEventTripletPtvsQ3Protons[1], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,1, *fConfig, fKinematics[1], fPrimAngles[1], fDeta[1], fDphi[1], fInvMass[1]);
       // Proton Proton Proton
       FillTripletDistribution( ParticleVector, 0, 0, 0, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletPtPrimaries[2], fSameEventTripletPtProtons[2], fSameEventTripletPtPrimaries2[2], fSameEventTripletPtProtons2[2], fSameEventTripletPtvsQ3Primaries[2], fSameEventTripletPtvsQ3Protons[2], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,2, *fConfig, fKinematics[2], fPrimAngles[2], fDeta[2], fDphi[2], fInvMass[2]);
       // Antiproton Antiproton Antiproton
       FillTripletDistribution( ParticleVector, 1, 1, 1, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletPtPrimaries[3], fSameEventTripletPtProtons[3], fSameEventTripletPtPrimaries2[3], fSameEventTripletPtProtons2[3], fSameEventTripletPtvsQ3Primaries[3], fSameEventTripletPtvsQ3Protons[3], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,3, *fConfig, fKinematics[3], fPrimAngles[3], fDeta[3], fDphi[3], fInvMass[3]);
       // Proton Proton AntiPrimary
       FillTripletDistribution( ParticleVector, 3, 0, 0, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletPtPrimaries[4], fSameEventTripletPtProtons[4], fSameEventTripletPtPrimaries2[4], fSameEventTripletPtProtons2[4], fSameEventTripletPtvsQ3Primaries[4], fSameEventTripletPtvsQ3Protons[4], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,4, *fConfig, fKinematics[4], fPrimAngles[4], fDeta[4], fDphi[4], fInvMass[4]);
       // Antiproton Antiproton Primary
       FillTripletDistribution( ParticleVector, 2, 1, 1, fSameEventTripletArray[5],PDGCodes, bins[1],fSameEventTripletMultArray[5], fSameEventTripletPtPrimaries[5], fSameEventTripletPtProtons[5], fSameEventTripletPtPrimaries2[5], fSameEventTripletPtProtons2[5], fSameEventTripletPtvsQ3Primaries[5], fSameEventTripletPtvsQ3Protons[5], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,5, *fConfig, fKinematics[5], fPrimAngles[5], fDeta[5], fDphi[5], fInvMass[5]);

    }//if(fDoOnlyThreeBody)
    else {
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

      if(fRunOfficialTwoBody){

        FillPairTransverseMass(ParticleVector, 0, 0, fPairTranverseMass_TwoBody[0], PDGCodes, fPairTranverseMassVSkstar_TwoBody[0]);
        FillPairTransverseMass(ParticleVector, 1, 1, fPairTranverseMass_TwoBody[1], PDGCodes, fPairTranverseMassVSkstar_TwoBody[1]);
        FillPairTransverseMass(ParticleVector, 0, 2, fPairTranverseMass_TwoBody[2], PDGCodes, fPairTranverseMassVSkstar_TwoBody[2]);
        FillPairTransverseMass(ParticleVector, 1, 3, fPairTranverseMass_TwoBody[3], PDGCodes, fPairTranverseMassVSkstar_TwoBody[3]);
        FillPairTransverseMass(ParticleVector, 0, 3, fPairTranverseMass_TwoBody[4], PDGCodes, fPairTranverseMassVSkstar_TwoBody[4]);
        FillPairTransverseMass(ParticleVector, 1, 2, fPairTranverseMass_TwoBody[5], PDGCodes, fPairTranverseMassVSkstar_TwoBody[5]);
        FillPairTransverseMass(ParticleVector, 2, 3, fPairTranverseMass_TwoBody[6], PDGCodes, fPairTranverseMassVSkstar_TwoBody[6]);

      }
  // Two Body Same Event
    }//else

    // c.1) Mixed event distribution --------------------------------

    // TAKE CARE OF MULT AND ZVtx!
    if (!(bins[0] == -99 || bins[1] == -99)) {

      std::vector<std::vector<AliFemtoDreamPartContainer>*> VectItMult;


      if(fStandardMixing){
        auto itZVtx = fVectPartContainers[0].begin()+ bins[0];
        auto itMult = itZVtx->begin() + bins[1];
        VectItMult.push_back(&(*itMult));
      
      } else {
      
        for(int i=0; i<3; i++){
          auto itZVtx = fVectPartContainers[i].begin()+ bins[0];
          auto itMult = itZVtx->begin() + bins[1];
          VectItMult.push_back(&(*itMult));
        }

      }



      //c.1.0) Same 2, Mixed 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(fDoOnlyThreeBody){

         // Proton Proton Primary
        FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 0, 0, 2, fSameEventTripletArray[10], PDGCodes, bins[1],fSameEventTripletMultArray[10], fSameEventTripletPtPrimaries[10], fSameEventTripletPtProtons[10], fSameEventTripletPtPrimaries2[10], fSameEventTripletPtProtons2[10], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 10, *fConfig, fQ3VskDistributionsArrayq12[0],fQ3VskDistributionsArrayq23[0], fKinematics[6], fPrimAngles[6], fDeta[6], fDphi[6], fInvMass[6]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 0, 2, 0, fSameEventTripletArray[11], PDGCodes, bins[1],fSameEventTripletMultArray[11], fSameEventTripletPtPrimaries[11], fSameEventTripletPtProtons[11], fSameEventTripletPtPrimaries2[11], fSameEventTripletPtProtons2[11], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 11, *fConfig, fQ3VskDistributionsArrayq12[1],fQ3VskDistributionsArrayq23[1], fKinematics[7], fPrimAngles[7], fDeta[7], fDphi[7], fInvMass[7]);

         // Antiproton Antiproton Antiprimary
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 1, 1, 3, fSameEventTripletArray[12], PDGCodes, bins[1],fSameEventTripletMultArray[12], fSameEventTripletPtPrimaries[12], fSameEventTripletPtProtons[12], fSameEventTripletPtPrimaries2[12], fSameEventTripletPtProtons2[12], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 12, *fConfig, fQ3VskDistributionsArrayq12[2],fQ3VskDistributionsArrayq23[2], fKinematics[8], fPrimAngles[8], fDeta[8], fDphi[8], fInvMass[8]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 1, 3, 1, fSameEventTripletArray[13], PDGCodes, bins[1],fSameEventTripletMultArray[13], fSameEventTripletPtPrimaries[13], fSameEventTripletPtProtons[13], fSameEventTripletPtPrimaries2[13], fSameEventTripletPtProtons2[13], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 13, *fConfig, fQ3VskDistributionsArrayq12[3],fQ3VskDistributionsArrayq23[3], fKinematics[9], fPrimAngles[9], fDeta[9], fDphi[9], fInvMass[9]);

         // Proton Proton Proton
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPP], 0, 0, 0, fSameEventTripletArray[14], PDGCodes, bins[1],fSameEventTripletMultArray[14], fSameEventTripletPtPrimaries[14], fSameEventTripletPtProtons[14], fSameEventTripletPtPrimaries2[14], fSameEventTripletPtProtons2[14], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 14, *fConfig, fQ3VskDistributionsArrayq12[4],fQ3VskDistributionsArrayq23[4], fKinematics[10], fPrimAngles[10], fDeta[10], fDphi[10], fInvMass[10]);

         // Antiproton Antiproton Antiproton
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPP], 1, 1, 1, fSameEventTripletArray[15], PDGCodes, bins[1],fSameEventTripletMultArray[15], fSameEventTripletPtPrimaries[15], fSameEventTripletPtProtons[15], fSameEventTripletPtPrimaries2[15], fSameEventTripletPtProtons2[15], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 15, *fConfig, fQ3VskDistributionsArrayq12[5],fQ3VskDistributionsArrayq23[5], fKinematics[11], fPrimAngles[11], fDeta[11], fDphi[11], fInvMass[11]);

         // Proton Proton AntiPrimary
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 0, 0, 3, fSameEventTripletArray[16], PDGCodes, bins[1],fSameEventTripletMultArray[16], fSameEventTripletPtPrimaries[16], fSameEventTripletPtProtons[16], fSameEventTripletPtPrimaries2[16], fSameEventTripletPtProtons2[16], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 16, *fConfig, fQ3VskDistributionsArrayq12[6],fQ3VskDistributionsArrayq23[6], fKinematics[12], fPrimAngles[12], fDeta[12], fDphi[12], fInvMass[12]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 0, 3, 0, fSameEventTripletArray[17], PDGCodes, bins[1],fSameEventTripletMultArray[17], fSameEventTripletPtPrimaries[17], fSameEventTripletPtProtons[17], fSameEventTripletPtPrimaries2[17], fSameEventTripletPtProtons2[17], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 17, *fConfig, fQ3VskDistributionsArrayq12[7],fQ3VskDistributionsArrayq23[7], fKinematics[13], fPrimAngles[13], fDeta[13], fDphi[13], fInvMass[13]);

         // Antiproton Antiproton Primary
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 1, 1, 2, fSameEventTripletArray[18], PDGCodes, bins[1],fSameEventTripletMultArray[18], fSameEventTripletPtPrimaries[18], fSameEventTripletPtProtons[18], fSameEventTripletPtPrimaries2[18], fSameEventTripletPtProtons2[18], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 18, *fConfig, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8], fKinematics[14], fPrimAngles[14], fDeta[14], fDphi[14], fInvMass[14]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 1, 2, 1, fSameEventTripletArray[19], PDGCodes, bins[1],fSameEventTripletMultArray[19], fSameEventTripletPtPrimaries[19], fSameEventTripletPtProtons[19], fSameEventTripletPtPrimaries2[19], fSameEventTripletPtProtons2[19], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 19, *fConfig, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9], fKinematics[15], fPrimAngles[15], fDeta[15], fDphi[15], fInvMass[15]);

      }//if(fDoOnlyThreeBody)

      //c.1.1) Normal mixing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(fDoOnlyThreeBody){

           // 2 (Anti-)Protons, 1 (Anti-)Primary
           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPPrim], 2, 0, 0, fMixedEventTripletArray[0], PDGCodes, bins[1],fMixedEventTripletMultArray[0], fMixedEventTripletPtPrimaries[0], fMixedEventTripletPtProtons[0], fMixedEventTripletPtPrimaries2[0], fMixedEventTripletPtProtons2[0], fMixedEventTripletPtvsQ3Primaries[0], fMixedEventTripletPtvsQ3Protons[0], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,0, *fConfig, fQ3VskDistributionsArrayq12Mixed[0], fQ3VskDistributionsArrayq23Mixed[0], fKinematicsME[0], fPrimAnglesME[0], fDetaME[0], fDphiME[0]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPPrim], 3, 1, 1, fMixedEventTripletArray[1], PDGCodes, bins[1],fMixedEventTripletMultArray[1], fMixedEventTripletPtPrimaries[1], fMixedEventTripletPtProtons[1], fMixedEventTripletPtPrimaries2[1], fMixedEventTripletPtProtons2[1], fMixedEventTripletPtvsQ3Primaries[1], fMixedEventTripletPtvsQ3Protons[1], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,1, *fConfig, fQ3VskDistributionsArrayq12Mixed[1], fQ3VskDistributionsArrayq23Mixed[1], fKinematicsME[1], fPrimAnglesME[1], fDetaME[1], fDphiME[1]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPP], 0, 0, 0, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletPtPrimaries[2], fMixedEventTripletPtProtons[2], fMixedEventTripletPtPrimaries2[2], fMixedEventTripletPtProtons2[2], fMixedEventTripletPtvsQ3Primaries[2], fMixedEventTripletPtvsQ3Protons[2], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,2, *fConfig, fQ3VskDistributionsArrayq12Mixed[2], fQ3VskDistributionsArrayq23Mixed[2], fKinematicsME[2], fPrimAnglesME[2], fDetaME[2], fDphiME[2]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPP], 1, 1, 1, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletPtPrimaries[3], fMixedEventTripletPtProtons[3], fMixedEventTripletPtPrimaries2[3], fMixedEventTripletPtProtons2[3], fMixedEventTripletPtvsQ3Primaries[3], fMixedEventTripletPtvsQ3Protons[3], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,3, *fConfig, fQ3VskDistributionsArrayq12Mixed[3], fQ3VskDistributionsArrayq23Mixed[3], fKinematicsME[3], fPrimAnglesME[3], fDetaME[3], fDphiME[3]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPAPrim], 3, 0, 0, fMixedEventTripletArray[4], PDGCodes, bins[1],fMixedEventTripletMultArray[4], fMixedEventTripletPtPrimaries[4], fMixedEventTripletPtProtons[4], fMixedEventTripletPtPrimaries2[4], fMixedEventTripletPtProtons2[4], fMixedEventTripletPtvsQ3Primaries[4], fMixedEventTripletPtvsQ3Protons[4], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,4, *fConfig, fQ3VskDistributionsArrayq12Mixed[4], fQ3VskDistributionsArrayq23Mixed[4], fKinematicsME[4], fPrimAnglesME[4], fDetaME[4], fDphiME[4]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPAPrim], 2, 1, 1, fMixedEventTripletArray[5], PDGCodes, bins[1],fMixedEventTripletMultArray[5], fMixedEventTripletPtPrimaries[5], fMixedEventTripletPtProtons[5], fMixedEventTripletPtPrimaries2[5], fMixedEventTripletPtProtons2[5], fMixedEventTripletPtvsQ3Primaries[5], fMixedEventTripletPtvsQ3Protons[5], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,5, *fConfig, fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5], fKinematicsME[5], fPrimAnglesME[5], fDetaME[5], fDphiME[5]);
  
      } else {
        //Two Body Analyses...........

        //Proton Proton
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPP], 0, 0, fMixedEventTripletArray_TwoBody[0],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[0], fMixedEventTripletPhiThetaArray_TwoBody,0, *fConfig);
        //AntiProton AntiProton
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPP], 1, 1, fMixedEventTripletArray_TwoBody[1],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[1], fMixedEventTripletPhiThetaArray_TwoBody,1, *fConfig);
        //Proton Primary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPPrim], 0, 2, fMixedEventTripletArray_TwoBody[2],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[2], fMixedEventTripletPhiThetaArray_TwoBody,2, *fConfig);
        //Antiproton Antiprimary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPPrim], 1, 3, fMixedEventTripletArray_TwoBody[3],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[3], fMixedEventTripletPhiThetaArray_TwoBody,3, *fConfig);
        //Proton Antiprimary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPAPrim], 0, 3, fMixedEventTripletArray_TwoBody[4],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[4], fMixedEventTripletPhiThetaArray_TwoBody,4, *fConfig);
        //Antiproton Primary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPAPrim], 1, 2, fMixedEventTripletArray_TwoBody[5],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[5], fMixedEventTripletPhiThetaArray_TwoBody,5, *fConfig);
        //Primary Antiprimary
        //FillPairDistributionME( ParticleVector, *itMult, 2, 3, fMixedEventTripletArray_TwoBody[6],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[6], fMixedEventTripletPhiThetaArray_TwoBody,6, *fConfig);

      }


      // Update the particle container with current event
      if(fStandardMixing) {

        SetMixedEvent(ParticleVector, VectItMult[0]);

      }else{

        if(fDoOnlyThreeBody){

          SetMixedEventPPP(ParticleVector, VectItMult[0]);
          SetMixedEventPPPrim(ParticleVector, VectItMult[1]);
          SetMixedEventPPAPrim(ParticleVector, VectItMult[2]);

        } else {
        
          SetMixedEventPP(ParticleVector, VectItMult[0]);
          SetMixedEventPPrim(ParticleVector, VectItMult[1]);
          SetMixedEventPAPrim(ParticleVector, VectItMult[2]);

        }
      }

    }//if (!(bins[0] == -99 || bins[1] == -99))

  }//if(fRunThreeBody)


  //----- OFFICIAL CODE for the 2-body -----------------
  if(fRunOfficialTwoBody){ // ADDED BY RAFFA
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
  //-----------------------------------------------------

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPrimaryList);
  PostData(5, fAntiPrimaryList);
  /*
  if(fCleanWithLambdas){
    PostData(6, fLambdaList);
    PostData(7, fAntiLambdaList);
  }
  */
  if(fRunOfficialTwoBody){ // ADDED BY RAFFA
    PostData(8, fResults);
    PostData(9, fResultsQA);
    PostData(10, fResultsSample);
    PostData(11, fResultsSampleQA);
  }
  PostData(12, fResultsThreeBody);
  if (fProton->GetIsMonteCarlo()) {
    PostData(13, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    PostData(14, fAntiProtonMCList);
  }
  if (fPrimary->GetIsMonteCarlo()) {
    PostData(15, fPrimaryMCList);
  }
  if (fAntiPrimary->GetIsMonteCarlo()) {
    PostData(16, fAntiPrimaryMCList);
  }
  /*
  if(fCleanWithLambdas){
    if (fLambda->GetIsMonteCarlo()) {
      PostData(17, fLambdaMCList);
    }
    if (fAntiLambda->GetIsMonteCarlo()) {
      PostData(18, fAntiLambdaMCList);
    }
  }
  */
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

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F* hPtvsQ3Primaries, TH2F* hPtvsQ3Protons, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* InvMass){
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
  if(abs(*itPDGPar1)==321) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==3122) DaughterPart2 = 2;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGPar2)==321) DaughterPart2 = 1;
  if(abs(*itPDGPar3)==3122) DaughterPart3 = 2;
  if(abs(*itPDGPar3)==2212) DaughterPart3 = 1;
  if(abs(*itPDGPar3)==321) DaughterPart3 = 1;
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
            if(abs(*itPDGPar1)==abs(*itPDGPar2)){
              Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
            }else{
              Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies,*iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
            }
            if(abs(*itPDGPar2)==abs(*itPDGPar3)){
              Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
            }else{
              Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
            }
            if(abs(*itPDGPar3)==abs(*itPDGPar1)){
              Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
            }else{
              Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
            }
          }
          if(!fClosePairRejectionForAll){
           
              if(DoThisPair12==11){
                if(abs(*itPDGPar1)==abs(*itPDGPar2)){
                  Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
                }else{
                  Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
                }
              }
              if(DoThisPair23==11){
                if(abs(*itPDGPar2)==abs(*itPDGPar3)){
                  Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
                }else{
                  Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
                }
              }
              if(DoThisPair31==11){
                if(abs(*itPDGPar3)==abs(*itPDGPar1)){
                  Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
                }else{
                  Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
                }
              }

          }

        }


        if(!Pair12||!Pair23||!Pair31) {continue;}

        hist->Fill(Q3);
        hist2d->Fill(Q3,mult+1);

        if(fRunPlotPt){

          if (Q3<1.) // MODIFIED BY ME
          {
            if(abs(*itPDGPar1)==321) hPtPrimaries -> Fill(part1_LorVec.Pt());
            if(abs(*itPDGPar1)==2212) hPtProtons -> Fill(part1_LorVec.Pt());
            if(abs(*itPDGPar2)==321) hPtPrimaries -> Fill(part2_LorVec.Pt());
            if(abs(*itPDGPar2)==2212) hPtProtons -> Fill(part2_LorVec.Pt());
            if(abs(*itPDGPar3)==321) hPtPrimaries -> Fill(part3_LorVec.Pt());
            if(abs(*itPDGPar3)==2212) hPtProtons -> Fill(part3_LorVec.Pt());
          }
          if (Q3>0.4&&Q3<1.) // REGION OF THE LAMBDA (1520)
          {
            if(abs(*itPDGPar1)==321) hPtPrimaries2 -> Fill(part1_LorVec.Pt());
            if(abs(*itPDGPar1)==2212) hPtProtons2 -> Fill(part1_LorVec.Pt());
            if(abs(*itPDGPar2)==321) hPtPrimaries2 -> Fill(part2_LorVec.Pt());
            if(abs(*itPDGPar2)==2212) hPtProtons2 -> Fill(part2_LorVec.Pt());
            if(abs(*itPDGPar3)==321) hPtPrimaries2 -> Fill(part3_LorVec.Pt());
            if(abs(*itPDGPar3)==2212) hPtProtons2 -> Fill(part3_LorVec.Pt());
          }

           if (Q3<1.2) // MODIFIED BY ME
          {
            if(abs(*itPDGPar1)==321) hPtvsQ3Primaries -> Fill(Q3,part1_LorVec.Pt());
            if(abs(*itPDGPar1)==2212) hPtvsQ3Protons -> Fill(Q3,part1_LorVec.Pt());
            if(abs(*itPDGPar2)==321) hPtvsQ3Primaries -> Fill(Q3,part2_LorVec.Pt());
            if(abs(*itPDGPar2)==2212) hPtvsQ3Protons -> Fill(Q3,part2_LorVec.Pt());
            if(abs(*itPDGPar3)==321) hPtvsQ3Primaries -> Fill(Q3,part3_LorVec.Pt());
            if(abs(*itPDGPar3)==2212) hPtvsQ3Protons -> Fill(Q3,part3_LorVec.Pt());
          }


        }

        if(fDoKinematicsPlots){

          TVector3 momPart1,momPart2,momPart3;

          momPart1.SetX(part1_LorVec.Px());
          momPart1.SetY(part1_LorVec.Py());
          momPart1.SetZ(part1_LorVec.Pz());

          momPart2.SetX(part2_LorVec.Px());
          momPart2.SetY(part2_LorVec.Py());
          momPart2.SetZ(part2_LorVec.Pz());

          momPart3.SetX(part3_LorVec.Px());
          momPart3.SetY(part3_LorVec.Py());
          momPart3.SetZ(part3_LorVec.Pz());

          double theta12 = momPart1.Angle(momPart2);
          double theta23 = momPart2.Angle(momPart3);
          double theta31 = momPart3.Angle(momPart1);

          //cout<<theta12<<"  "<<theta23<<" "<<theta31<<endl;

          double thetaMin=0, thetaMid=0, thetaMax=0,thetaPP=0;
          double deltaetaMin = 0, deltaphiMin = 0, deltaDeta = 0, deltaPhi = 0;


          if(abs(*itPDGPar1)==321){
            if(theta12<theta31){
              thetaMin = theta12;
              thetaMax = theta31;
              auto Eta1 = iPart1->GetEta().at(0);
              auto Eta2 = iPart2->GetEta().at(0);
              deltaetaMin = abs(Eta1-Eta2);
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart2->GetPhi().at(0));
              deltaDeta = abs(iPart1->GetEta().at(0) - iPart3->GetEta().at(0));
              deltaPhi = abs(iPart1->GetPhi().at(0) - iPart3->GetPhi().at(0));
            }else{
              thetaMin = theta31;
              thetaMax = theta12;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart1->GetEta().at(0) - iPart2->GetEta().at(0));
              deltaPhi = abs(iPart1->GetPhi().at(0) - iPart2->GetPhi().at(0));
            }
            thetaPP = theta23;
          }

          if(abs(*itPDGPar2)==321){
            if(theta12<theta23){
              thetaMin = theta12;
              thetaMax = theta23;
              auto Eta1 = iPart1->GetEta().at(0);
              auto Eta2 = iPart2->GetEta().at(0);
              deltaetaMin = abs(Eta1-Eta2);
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart2->GetPhi().at(0));
              deltaDeta = abs(iPart2->GetEta().at(0) - iPart3->GetEta().at(0));
              deltaPhi = abs(iPart2->GetPhi().at(0) - iPart3->GetPhi().at(0));
            }else{
              thetaMin = theta23;
              thetaMax = theta12;
              deltaetaMin = abs(iPart2->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart2->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart2->GetEta().at(0) - iPart1->GetEta().at(0));
              deltaPhi = abs(iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0));
            }
            thetaPP = theta31;
          }

          if(abs(*itPDGPar3)==321){
            if(theta31<theta23){
              thetaMin = theta31;
              thetaMax = theta23;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart3->GetEta().at(0) - iPart2->GetEta().at(0));
              deltaPhi = abs(iPart3->GetPhi().at(0) - iPart2->GetPhi().at(0));
            }else{
              thetaMin = theta23;
              thetaMax = theta31;
              deltaetaMin = abs(iPart2->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart2->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart3->GetEta().at(0) - iPart1->GetEta().at(0));
              deltaPhi = abs(iPart3->GetPhi().at(0) - iPart1->GetPhi().at(0));
            }
            thetaPP = theta12;
          }

          if((abs(*itPDGPar1)==321||abs(*itPDGPar2)==321||abs(*itPDGPar3)==321)&&Q3<1.2){
            hPrimAngles->Fill(thetaMin,thetaPP);
            hDeta->Fill(deltaetaMin,deltaDeta);
            hDphi->Fill(deltaphiMin,deltaPhi);
          }

//          cout<<thetaMin<<"  "<<thetaMax<<" "<<thetaPP<<endl;

          TLorentzVector trackSum = part1_LorVec + part2_LorVec + part3_LorVec;
          double PCM = trackSum.P();
          double ECM = trackSum.E();
          double MassCM = part1_LorVec.M()+part2_LorVec.M()+part3_LorVec.M();//sqrt(pow(ECM,2)-pow(PCM,2));

//          cout<<MassCM<<"  "<<part1_LorVec.M()+part2_LorVec.M()+part3_LorVec.M()<<endl;

          double p1xCM = part1_LorVec.Px()*ECM/MassCM - part1_LorVec.E()*trackSum.Px()/MassCM;
          double p1yCM = part1_LorVec.Py()*ECM/MassCM - part1_LorVec.E()*trackSum.Py()/MassCM;
          double p1zCM = part1_LorVec.Pz()*ECM/MassCM - part1_LorVec.E()*trackSum.Pz()/MassCM;

          double p2xCM = part2_LorVec.Px()*ECM/MassCM - part2_LorVec.E()*trackSum.Px()/MassCM;
          double p2yCM = part2_LorVec.Py()*ECM/MassCM - part2_LorVec.E()*trackSum.Py()/MassCM;
          double p2zCM = part2_LorVec.Pz()*ECM/MassCM - part2_LorVec.E()*trackSum.Pz()/MassCM;

          double p3xCM = part3_LorVec.Px()*ECM/MassCM - part3_LorVec.E()*trackSum.Px()/MassCM;
          double p3yCM = part3_LorVec.Py()*ECM/MassCM - part3_LorVec.E()*trackSum.Py()/MassCM;
          double p3zCM = part3_LorVec.Pz()*ECM/MassCM - part3_LorVec.E()*trackSum.Pz()/MassCM;

          TVector3 momPart1CM,momPart2CM,momPart3CM;

          momPart1CM.SetX(p1xCM);
          momPart1CM.SetY(p1yCM);
          momPart1CM.SetZ(p1zCM);

          momPart2CM.SetX(p2xCM);
          momPart2CM.SetY(p2yCM);
          momPart2CM.SetZ(p2zCM);

          momPart3CM.SetX(p3xCM);
          momPart3CM.SetY(p3yCM);
          momPart3CM.SetZ(p3zCM);

          double theta12CM = momPart1CM.Angle(momPart2CM);
          double theta23CM = momPart2CM.Angle(momPart3CM);
          double theta31CM = momPart3CM.Angle(momPart1CM);

      //   double thetaMin=0,thetaMax=0,thetaPP=0;

          if(abs(*itPDGPar1)==321){
            if(theta12CM<theta31CM){
              thetaMin = theta12CM;
              thetaMax = theta31CM;
            }else{
              thetaMin = theta31CM;
              thetaMax = theta12CM;
            }
            thetaPP = theta23CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGPar2)==321){
            if(theta12CM<theta23CM){
              thetaMin = theta12CM;
              thetaMax = theta23CM;
            }else{
              thetaMin = theta23CM;
              thetaMax = theta12CM;
            }
            thetaPP = theta31CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGPar3)==321){
            if(theta31CM<theta23CM){
              thetaMin = theta31CM;
              thetaMax = theta23CM;
            }else{
              thetaMin = theta23CM;
              thetaMax = theta31CM;
            }
            thetaPP = theta12CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGPar1)==321||abs(*itPDGPar2)==321||abs(*itPDGPar3)==321){
            if(Q3<1.2) hKinematics->Fill(thetaMin+thetaMid,thetaMid-thetaMin);
          }

        }//if(fDoKinematicsPlots)

   
        if(fRunPlotInvMass)
        {
           TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
           TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
           TLorentzVector Sum13 = part1_LorVec + part3_LorVec;

           TLorentzVector Sum123 =  part1_LorVec + part2_LorVec + part3_LorVec;

           InvMass->Fill(Sum12.Mag(),Q3);
           InvMass->Fill(Sum23.Mag(),Q3);
	         InvMass->Fill(Sum13.Mag(),Q3);
           InvMass->Fill(Sum123.Mag(),Q3);

        }//if(fRunPlotInvMass)

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
  if(abs(*itPDGPar1)==321) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==3122) DaughterPart2 = 2;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGPar2)==321) DaughterPart2 = 1;


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

          if(fClosePairRejectionForAll){
             Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[7+phiEtaHistNo],Config); //GANESHA  Check about DeltaEtaDeltaPhi Function. Might need to change for 2 Body
          }
          else{

             if(DoThisPair12==11){
               Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[7+phiEtaHistNo],Config); //GANESHA  Check about DeltaEtaDeltaPhi Function. Might need to change for 2 Body
             }
          } 
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

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventPP(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda.
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 1) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
  }
  if ((ParticleVector.begin()+1)->size() > 1) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventPPrim(
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

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventPAPrim(  
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda.
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 0 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
  if ((ParticleVector.begin()+1)->size() > 0 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventPPP(
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

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventPPPrim(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and primary, 1 and 3 is antiproton antiprimary.
  // Fill the particles only if at least 2 prontons and a primary particle are present in the event
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

void AliAnalysisTaskThreeBodyProtonPrimary::SetMixedEventPPAPrim(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and primary, 1 and 3 is antiproton antiprimary.
  // Fill the particles only if at least 2 prontons and a primary particle are present in the event
  if ((ParticleVector.begin())->size() > 1 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
  if ((ParticleVector.begin()+1)->size() > 1 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F* hPtvsQ3Primaries, TH2F* hPtvsQ3Protons, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi){//, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed){
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
  if(abs(*itPDGParSE)==321) DaughterPart1 = 1;
  if(abs(*itPDGParME1)==3122) DaughterPart2 = 2;
  if(abs(*itPDGParME1)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParME1)==321) DaughterPart2 = 1;
  if(abs(*itPDGParME2)==3122) DaughterPart3 = 2;
  if(abs(*itPDGParME2)==2212) DaughterPart3 = 1;
  if(abs(*itPDGParME2)==321) DaughterPart3 = 1;
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
                if(abs(*itPDGParSE)==abs(*itPDGParME1)){
                  Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[10+phiEtaHistNo],Config, Q3);
                }else{
                  Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[10+phiEtaHistNo],Config, Q3);
                }
                if(abs(*itPDGParME1)==abs(*itPDGParME2)){
                  Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[10+phiEtaHistNo],Config, Q3);
                }else{
                  Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[10+phiEtaHistNo],Config, Q3);
                }
                if(abs(*itPDGParME2)==abs(*itPDGParSE)){
                  Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[10+phiEtaHistNo],Config, Q3);
                }else{
                  Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[10+phiEtaHistNo],Config, Q3);
                }
              }
              if(!fClosePairRejectionForAll){
               
                  if(DoThisPair12==11){
                    if(abs(*itPDGParSE)==abs(*itPDGParME1)){
                      Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[10+phiEtaHistNo],Config, Q3);
                    }else{
                      Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[10+phiEtaHistNo],Config, Q3);
                    }
                  }
                  if(DoThisPair23==11){
                    if(abs(*itPDGParME1)==abs(*itPDGParME2)){
                      Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[10+phiEtaHistNo],Config, Q3);
                    }else{
                      Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[10+phiEtaHistNo],Config, Q3);
                    }
                  }
                  if(DoThisPair31==11){
                    if(abs(*itPDGParME2)==abs(*itPDGParSE)){
                      Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[10+phiEtaHistNo],Config, Q3);
                    }else{
                      Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[10+phiEtaHistNo],Config, Q3);
                    }
                  }  

              }

            }

            if(!Pair12||!Pair23||!Pair31) {continue;}


            // Now we have the three particles, lets create their Lorentz vectors

            hist->Fill(Q3);
            hist2d->Fill(Q3,mult+1);


            if(fRunPlotPt){

              if (Q3<1.) // MODIFIED BY ME
              {
                if(abs(*itPDGParSE)==321) hPtPrimaries -> Fill(part1_LorVec.Pt());
                if(abs(*itPDGParSE)==2212) hPtProtons -> Fill(part1_LorVec.Pt());
                if(abs(*itPDGParME1)==321) hPtPrimaries -> Fill(part2_LorVec.Pt());
                if(abs(*itPDGParME1)==2212) hPtProtons -> Fill(part2_LorVec.Pt());
                if(abs(*itPDGParME2)==321) hPtPrimaries -> Fill(part3_LorVec.Pt());
                if(abs(*itPDGParME2)==2212) hPtProtons -> Fill(part3_LorVec.Pt());
              }
              if (Q3>0.4&&Q3<1.) // REGION OF THE LAMBDA (1520)
              {
                if(abs(*itPDGParSE)==321) hPtPrimaries2 -> Fill(part1_LorVec.Pt());
                if(abs(*itPDGParSE)==2212) hPtProtons2 -> Fill(part1_LorVec.Pt());
                if(abs(*itPDGParME1)==321) hPtPrimaries2 -> Fill(part2_LorVec.Pt());
                if(abs(*itPDGParME1)==2212) hPtProtons2 -> Fill(part2_LorVec.Pt());
                if(abs(*itPDGParME2)==321) hPtPrimaries2 -> Fill(part3_LorVec.Pt());
                if(abs(*itPDGParME2)==2212) hPtProtons2 -> Fill(part3_LorVec.Pt());
              }

              if (Q3<1.2) // MODIFIED BY ME
              {
                if(abs(*itPDGParSE)==321) hPtvsQ3Primaries -> Fill(Q3,part1_LorVec.Pt());
                if(abs(*itPDGParSE)==2212) hPtvsQ3Protons -> Fill(Q3,part1_LorVec.Pt());
                if(abs(*itPDGParME1)==321) hPtvsQ3Primaries -> Fill(Q3,part2_LorVec.Pt());
                if(abs(*itPDGParME1)==2212) hPtvsQ3Protons -> Fill(Q3,part2_LorVec.Pt());
                if(abs(*itPDGParME2)==321) hPtvsQ3Primaries -> Fill(Q3,part3_LorVec.Pt());
                if(abs(*itPDGParME2)==2212) hPtvsQ3Protons -> Fill(Q3,part3_LorVec.Pt());
              }

            }
            
            if(fRunPlotQ3Vsq){
              float qSame12= sqrt(-q12*q12);
              float qSame23= sqrt(-q23*q23);
              Q3VskDistribution12Mixed->Fill(Q3, qSame12);
              Q3VskDistribution23Mixed->Fill(Q3, qSame23);
            }


        if(fDoKinematicsPlots){

          TVector3 momPart1,momPart2,momPart3;

          momPart1.SetX(part1_LorVec.Px());
          momPart1.SetY(part1_LorVec.Py());
          momPart1.SetZ(part1_LorVec.Pz());

          momPart2.SetX(part2_LorVec.Px());
          momPart2.SetY(part2_LorVec.Py());
          momPart2.SetZ(part2_LorVec.Pz());

          momPart3.SetX(part3_LorVec.Px());
          momPart3.SetY(part3_LorVec.Py());
          momPart3.SetZ(part3_LorVec.Pz());

          double theta12 = momPart1.Angle(momPart2);
          double theta23 = momPart2.Angle(momPart3);
          double theta31 = momPart3.Angle(momPart1);

          //cout<<theta12<<"  "<<theta23<<" "<<theta31<<endl;

          double thetaMin=0, thetaMid=0, thetaMax=0,thetaPP=0;
          double deltaetaMin = 0, deltaphiMin = 0, deltaDeta = 0, deltaPhi = 0;


          if(abs(*itPDGParSE)==321){
            if(theta12<theta31){
              thetaMin = theta12;
              thetaMax = theta31;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart2->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart2->GetPhi().at(0));
              deltaDeta = abs(iPart1->GetEta().at(0) - iPart3->GetEta().at(0));
              deltaPhi = abs(iPart1->GetPhi().at(0) - iPart3->GetPhi().at(0));
            }else{
              thetaMin = theta31;
              thetaMax = theta12;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart1->GetEta().at(0) - iPart2->GetEta().at(0));
              deltaPhi = abs(iPart1->GetPhi().at(0) - iPart2->GetPhi().at(0));
            }
            thetaPP = theta23;
          }

          if(abs(*itPDGParME1)==321){
            if(theta12<theta23){
              thetaMin = theta12;
              thetaMax = theta23;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart2->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart2->GetPhi().at(0));
              deltaDeta = abs(iPart2->GetEta().at(0) - iPart3->GetEta().at(0));
              deltaPhi = abs(iPart2->GetPhi().at(0) - iPart3->GetPhi().at(0));
            }else{
              thetaMin = theta23;
              thetaMax = theta12;
              deltaetaMin = abs(iPart2->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart2->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart2->GetEta().at(0) - iPart1->GetEta().at(0));
              deltaPhi = abs(iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0));
            }
            thetaPP = theta31;
          }

          if(abs(*itPDGParME2)==321){
            if(theta31<theta23){
              thetaMin = theta31;
              thetaMax = theta23;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart3->GetEta().at(0) - iPart2->GetEta().at(0));
              deltaPhi = abs(iPart3->GetPhi().at(0) - iPart2->GetPhi().at(0));
            }else{
              thetaMin = theta23;
              thetaMax = theta31;
              deltaetaMin = abs(iPart2->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart2->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart3->GetEta().at(0) - iPart1->GetEta().at(0));
              deltaPhi = abs(iPart3->GetPhi().at(0) - iPart1->GetPhi().at(0));
            }
            thetaPP = theta12;
          }

          if((abs(*itPDGParSE)==321||abs(*itPDGParME1)==321||abs(*itPDGParME2)==321)&&Q3<1.2){
            hPrimAngles->Fill(thetaMin,thetaPP);
            hDeta->Fill(deltaetaMin,deltaDeta);
            hDphi->Fill(deltaphiMin,deltaPhi);
          }


//          cout<<thetaMin<<"  "<<thetaMax<<" "<<thetaPP<<endl;


          TLorentzVector trackSum = part1_LorVec + part2_LorVec + part3_LorVec;
          double PCM = trackSum.P();
          double ECM = trackSum.E();
          double MassCM = part1_LorVec.M()+part2_LorVec.M()+part3_LorVec.M();//sqrt(pow(ECM,2)-pow(PCM,2));


//          cout<<MassCM<<"  "<<part1_LorVec.M()+part2_LorVec.M()+part3_LorVec.M()<<endl;


          double p1xCM = part1_LorVec.Px()*ECM/MassCM - part1_LorVec.E()*trackSum.Px()/MassCM;
          double p1yCM = part1_LorVec.Py()*ECM/MassCM - part1_LorVec.E()*trackSum.Py()/MassCM;
          double p1zCM = part1_LorVec.Pz()*ECM/MassCM - part1_LorVec.E()*trackSum.Pz()/MassCM;

          double p2xCM = part2_LorVec.Px()*ECM/MassCM - part2_LorVec.E()*trackSum.Px()/MassCM;
          double p2yCM = part2_LorVec.Py()*ECM/MassCM - part2_LorVec.E()*trackSum.Py()/MassCM;
          double p2zCM = part2_LorVec.Pz()*ECM/MassCM - part2_LorVec.E()*trackSum.Pz()/MassCM;

          double p3xCM = part3_LorVec.Px()*ECM/MassCM - part3_LorVec.E()*trackSum.Px()/MassCM;
          double p3yCM = part3_LorVec.Py()*ECM/MassCM - part3_LorVec.E()*trackSum.Py()/MassCM;
          double p3zCM = part3_LorVec.Pz()*ECM/MassCM - part3_LorVec.E()*trackSum.Pz()/MassCM;

          TVector3 momPart1CM,momPart2CM,momPart3CM;

          momPart1CM.SetX(p1xCM);
          momPart1CM.SetY(p1yCM);
          momPart1CM.SetZ(p1zCM);

          momPart2CM.SetX(p2xCM);
          momPart2CM.SetY(p2yCM);
          momPart2CM.SetZ(p2zCM);

          momPart3CM.SetX(p3xCM);
          momPart3CM.SetY(p3yCM);
          momPart3CM.SetZ(p3zCM);

          double theta12CM = momPart1CM.Angle(momPart2CM);
          double theta23CM = momPart2CM.Angle(momPart3CM);
          double theta31CM = momPart3CM.Angle(momPart1CM);

      //   double thetaMin=0,thetaMax=0,thetaPP=0;

          if(abs(*itPDGParSE)==321){
            if(theta12CM<theta31CM){
              thetaMin = theta12CM;
              thetaMax = theta31CM;
            }else{
              thetaMin = theta31CM;
              thetaMax = theta12CM;
            }
            thetaPP = theta23CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGParME1)==321){
            if(theta12CM<theta23CM){
              thetaMin = theta12CM;
              thetaMax = theta23CM;
            }else{
              thetaMin = theta23CM;
              thetaMax = theta12CM;
            }
            thetaPP = theta31CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGParME2)==321){
            if(theta31CM<theta23CM){
              thetaMin = theta31CM;
              thetaMax = theta23CM;
            }else{
              thetaMin = theta23CM;
              thetaMax = theta31CM;
            }
            thetaPP = theta12CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGParSE)==321||abs(*itPDGParME1)==321||abs(*itPDGParME2)==321){
            if(Q3<1.2) hKinematics->Fill(thetaMin+thetaMid,thetaMid-thetaMin);
          }

//          cout<<thetaMin<<"  "<<thetaMid<<" "<<thetaMax<<endl;



        }



       
          }
        }
      }
    }
  }
}

//==================================================================================================================================================

//void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionMEPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F* InvMassDET,TH2F* InvMassPDG)


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
  if(abs(*itPDGParSE)==321) DaughterPart1 = 1;
  if(abs(*itPDGParME1)==3122) DaughterPart2 = 2;
  if(abs(*itPDGParME1)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParME1)==321) DaughterPart2 = 1;

  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;

  // Get particle masses
  auto massParticleSE = TDatabasePDG::Instance()->GetParticle(*itPDGParSE)->Mass();
  auto massParticleME1 = TDatabasePDG::Instance()->GetParticle(*itPDGParME1)->Mass();

  // loop over first particle
  for (auto iPart1 = ParticleSE->begin(); iPart1 != ParticleSE->end(); ++iPart1) {
    // loop over second particle ...
    for (int iDepth1 = 0; iDepth1 < (int) MixedEvent1Container->GetMixingDepth(); ++iDepth1) {
      std::vector<AliFemtoDreamBasePart> iEvent2 = MixedEvent1Container->GetEvent(iDepth1); //GANESHA Here check if two species are the same
      for ( auto iPart2 = iEvent2.begin(); iPart2 != iEvent2.end(); ++iPart2) {

        
        bool Pair12 = true;

         if(!fturnoffClosePairRejectionCompletely){

          if(fClosePairRejectionForAll){
             Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[7+phiEtaHistNo],Config); //GANESHA Change number + Check about DeltaEtaDeltaPhi Function. Might need to change for 2 Body
          }
          else{

             if(DoThisPair12==11){
               Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1,*iPart1,*iPart2, *itPDGParSE, *itPDGParME1, false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[7+phiEtaHistNo],Config); //GANESHA Change number + Check about DeltaEtaDeltaPhi Function. Might need to change for 2 Body
             }
          } 
        }

        if(!Pair12) {continue;}

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

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* InvMass){
  // Description of function given in AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistribution
  // In this function, two particles are used from current event, and one - from other event
  auto ParticleSE1 = ParticleVector.begin()+speciesSE1;
  auto ParticleSE2 = ParticleVector.begin()+speciesSE2;
  auto MixedEventContainer = fPartContainer.begin()+speciesME;


  if(!fStandardMixing){
    if((speciesSE1==0&&speciesSE2==0&&speciesME==0)||(speciesSE1==1&&speciesSE2==1&&speciesME==1)){ // (p-p)-p
      if(ParticleVector[speciesSE1].size()<3) return;
    }
    if((speciesSE1==0&&speciesSE2==0&&speciesME==2)||(speciesSE1==1&&speciesSE2==1&&speciesME==3)){ // (p-p)-Prim
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesME].size()<1) return;
    }
    if((speciesSE1==0&&speciesSE2==2&&speciesME==0)||(speciesSE1==1&&speciesSE2==3&&speciesME==1)){ // (p-Prim)-p
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesSE2].size()<1) return;
    }
    if((speciesSE1==0&&speciesSE2==0&&speciesME==3)||(speciesSE1==1&&speciesSE2==1&&speciesME==2)){ // (p-p)-APrim
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesME].size()<1) return;
    }
    if((speciesSE1==0&&speciesSE2==3&&speciesME==0)||(speciesSE1==1&&speciesSE2==2&&speciesME==1)){ // (p-APrim)-p
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
  if(abs(*itPDGParSE1)==321) DaughterPart1 = 1;
  if(abs(*itPDGParSE2)==3122) DaughterPart2 = 2;
  if(abs(*itPDGParSE2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParSE2)==321) DaughterPart2 = 1;
  if(abs(*itPDGParME)==3122) DaughterPart3 = 2;
  if(abs(*itPDGParME)==2212) DaughterPart3 = 1;
  if(abs(*itPDGParME)==321) DaughterPart3 = 1;
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
              if(abs(*itPDGParSE1)==abs(*itPDGParSE2)){
                Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
              }else{
                Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
              }
              if(abs(*itPDGParSE2)==abs(*itPDGParME)){
                Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
              }else{
                Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
              }
              if(abs(*itPDGParSE1)==abs(*itPDGParME)){
                Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
              }else{
                Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
              }
            }
            if(!fClosePairRejectionForAll){
          
                if(DoThisPair12==11){
                  if(abs(*itPDGParSE1)==abs(*itPDGParSE2)){
                    Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
                  }else{
                    Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
                  }
                }
                if(DoThisPair23==11){
                  if(abs(*itPDGParSE2)==abs(*itPDGParME)){
                    Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
                  }else{
                    Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
                  }
                }
                if(DoThisPair31==11){
                  if(abs(*itPDGParSE1)==abs(*itPDGParME)){
                    Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[28+phiEtaHistNo],Config, Q3);
                  }else{
                    Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[28+phiEtaHistNo],Config, Q3);
                  }
                }
              
            }

          }

          if(!Pair12||!Pair23||!Pair31) {continue;}


          hist->Fill(Q3);
          hist2d->Fill(Q3,mult+1);

          if(fRunPlotPt){

            if (Q3<1.) // MODIFIED BY ME
            {
              if(abs(*itPDGParSE1)==321) hPtPrimaries -> Fill(part1_LorVec.Pt());
              if(abs(*itPDGParSE1)==2212) hPtProtons -> Fill(part1_LorVec.Pt());
              if(abs(*itPDGParSE2)==321) hPtPrimaries -> Fill(part2_LorVec.Pt());
              if(abs(*itPDGParSE2)==2212) hPtProtons -> Fill(part2_LorVec.Pt());
              if(abs(*itPDGParME)==321) hPtPrimaries -> Fill(part3_LorVec.Pt());
              if(abs(*itPDGParME)==2212) hPtProtons -> Fill(part3_LorVec.Pt());
            }
            if (Q3>0.4&&Q3<1.) // REGION OF THE LAMBDA (1520)
            {
              if(abs(*itPDGParSE1)==321) hPtPrimaries2 -> Fill(part1_LorVec.Pt());
              if(abs(*itPDGParSE1)==2212) hPtProtons2 -> Fill(part1_LorVec.Pt());
              if(abs(*itPDGParSE2)==321) hPtPrimaries2 -> Fill(part2_LorVec.Pt());
              if(abs(*itPDGParSE2)==2212) hPtProtons2 -> Fill(part2_LorVec.Pt());
              if(abs(*itPDGParME)==321) hPtPrimaries2 -> Fill(part3_LorVec.Pt());
              if(abs(*itPDGParME)==2212) hPtProtons2 -> Fill(part3_LorVec.Pt());
            }

          }
          if(fRunPlotQ3Vsq){
            float qSame12= sqrt(-q12*q12);
            float qSame23= sqrt(-q23*q23);
            Q3VskDistribution12->Fill(Q3, qSame12);
            Q3VskDistribution23->Fill(Q3, qSame23);
          }

        if(fDoKinematicsPlots){

          TVector3 momPart1,momPart2,momPart3;

          momPart1.SetX(part1_LorVec.Px());
          momPart1.SetY(part1_LorVec.Py());
          momPart1.SetZ(part1_LorVec.Pz());

          momPart2.SetX(part2_LorVec.Px());
          momPart2.SetY(part2_LorVec.Py());
          momPart2.SetZ(part2_LorVec.Pz());

          momPart3.SetX(part3_LorVec.Px());
          momPart3.SetY(part3_LorVec.Py());
          momPart3.SetZ(part3_LorVec.Pz());

          double theta12 = momPart1.Angle(momPart2);
          double theta23 = momPart2.Angle(momPart3);
          double theta31 = momPart3.Angle(momPart1);

          //cout<<theta12<<"  "<<theta23<<" "<<theta31<<endl;

          double thetaMin=0, thetaMid=0, thetaMax=0,thetaPP=0;
          double deltaetaMin = 0, deltaphiMin = 0, deltaDeta = 0, deltaPhi = 0;


          if(abs(*itPDGParSE1)==321){
            if(theta12<theta31){
              thetaMin = theta12;
              thetaMax = theta31;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart2->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart2->GetPhi().at(0));
              deltaDeta = abs(iPart1->GetEta().at(0) - iPart3->GetEta().at(0));
              deltaPhi = abs(iPart1->GetPhi().at(0) - iPart3->GetPhi().at(0));
            }else{
              thetaMin = theta31;
              thetaMax = theta12;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart1->GetEta().at(0) - iPart2->GetEta().at(0));
              deltaPhi = abs(iPart1->GetPhi().at(0) - iPart2->GetPhi().at(0));
            }
            thetaPP = theta23;
          }

          if(abs(*itPDGParSE2)==321){
            if(theta12<theta23){
              thetaMin = theta12;
              thetaMax = theta23;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart2->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart2->GetPhi().at(0));
              deltaDeta = abs(iPart2->GetEta().at(0) - iPart3->GetEta().at(0));
              deltaPhi = abs(iPart2->GetPhi().at(0) - iPart3->GetPhi().at(0));
            }else{
              thetaMin = theta23;
              thetaMax = theta12;
              deltaetaMin = abs(iPart2->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart2->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart2->GetEta().at(0) - iPart1->GetEta().at(0));
              deltaPhi = abs(iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0));
            }
            thetaPP = theta31;
          }

          if(abs(*itPDGParME)==321){
            if(theta31<theta23){
              thetaMin = theta31;
              thetaMax = theta23;
              deltaetaMin = abs(iPart1->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart1->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart3->GetEta().at(0) - iPart2->GetEta().at(0));
              deltaPhi = abs(iPart3->GetPhi().at(0) - iPart2->GetPhi().at(0));
            }else{
              thetaMin = theta23;
              thetaMax = theta31;
              deltaetaMin = abs(iPart2->GetEta().at(0)-iPart3->GetEta().at(0));
              deltaphiMin = abs(iPart2->GetPhi().at(0)-iPart3->GetPhi().at(0));
              deltaDeta = abs(iPart3->GetEta().at(0) - iPart1->GetEta().at(0));
              deltaPhi = abs(iPart3->GetPhi().at(0) - iPart1->GetPhi().at(0));
            }
            thetaPP = theta12;
          }

          if((abs(*itPDGParSE1)==321||abs(*itPDGParSE2)==321||abs(*itPDGParME)==321)&&Q3<1.2){
            hPrimAngles->Fill(thetaMin,thetaPP);
            hDeta->Fill(deltaetaMin,deltaDeta);
            hDphi->Fill(deltaphiMin,deltaPhi);
          }



//          cout<<thetaMin<<"  "<<thetaMax<<" "<<thetaPP<<endl;


          TLorentzVector trackSum = part1_LorVec + part2_LorVec + part3_LorVec;
          double PCM = trackSum.P();
          double ECM = trackSum.E();
          double MassCM = part1_LorVec.M()+part2_LorVec.M()+part3_LorVec.M();//sqrt(pow(ECM,2)-pow(PCM,2));


//          cout<<MassCM<<"  "<<part1_LorVec.M()+part2_LorVec.M()+part3_LorVec.M()<<endl;


          double p1xCM = part1_LorVec.Px()*ECM/MassCM - part1_LorVec.E()*trackSum.Px()/MassCM;
          double p1yCM = part1_LorVec.Py()*ECM/MassCM - part1_LorVec.E()*trackSum.Py()/MassCM;
          double p1zCM = part1_LorVec.Pz()*ECM/MassCM - part1_LorVec.E()*trackSum.Pz()/MassCM;

          double p2xCM = part2_LorVec.Px()*ECM/MassCM - part2_LorVec.E()*trackSum.Px()/MassCM;
          double p2yCM = part2_LorVec.Py()*ECM/MassCM - part2_LorVec.E()*trackSum.Py()/MassCM;
          double p2zCM = part2_LorVec.Pz()*ECM/MassCM - part2_LorVec.E()*trackSum.Pz()/MassCM;

          double p3xCM = part3_LorVec.Px()*ECM/MassCM - part3_LorVec.E()*trackSum.Px()/MassCM;
          double p3yCM = part3_LorVec.Py()*ECM/MassCM - part3_LorVec.E()*trackSum.Py()/MassCM;
          double p3zCM = part3_LorVec.Pz()*ECM/MassCM - part3_LorVec.E()*trackSum.Pz()/MassCM;

          TVector3 momPart1CM,momPart2CM,momPart3CM;

          momPart1CM.SetX(p1xCM);
          momPart1CM.SetY(p1yCM);
          momPart1CM.SetZ(p1zCM);

          momPart2CM.SetX(p2xCM);
          momPart2CM.SetY(p2yCM);
          momPart2CM.SetZ(p2zCM);

          momPart3CM.SetX(p3xCM);
          momPart3CM.SetY(p3yCM);
          momPart3CM.SetZ(p3zCM);

          double theta12CM = momPart1CM.Angle(momPart2CM);
          double theta23CM = momPart2CM.Angle(momPart3CM);
          double theta31CM = momPart3CM.Angle(momPart1CM);

      //   double thetaMin=0,thetaMax=0,thetaPP=0;

          if(abs(*itPDGParSE1)==321){
            if(theta12CM<theta31CM){
              thetaMin = theta12CM;
              thetaMax = theta31CM;
            }else{
              thetaMin = theta31CM;
              thetaMax = theta12CM;
            }
            thetaPP = theta23CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGParSE2)==321){
            if(theta12CM<theta23CM){
              thetaMin = theta12CM;
              thetaMax = theta23CM;
            }else{
              thetaMin = theta23CM;
              thetaMax = theta12CM;
            }
            thetaPP = theta31CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGParME)==321){
            if(theta31CM<theta23CM){
              thetaMin = theta31CM;
              thetaMax = theta23CM;
            }else{
              thetaMin = theta23CM;
              thetaMax = theta31CM;
            }
            thetaPP = theta12CM;
            if(thetaPP<thetaMin){
              thetaMid = thetaMin;
              thetaMin = thetaPP;
            }else{
              if(thetaPP>thetaMax){
                thetaMid = thetaMax;
                thetaMax = thetaPP;
              }else{
                thetaMid = thetaPP;
              }
            }
          }

          if(abs(*itPDGParSE1)==321||abs(*itPDGParSE2)==321||abs(*itPDGParME)==321){
            if(Q3<1.2) hKinematics->Fill(thetaMin+thetaMid,thetaMid-thetaMin);
          }

//          cout<<thetaMin<<"  "<<thetaMid<<" "<<thetaMax<<endl;



        } //if(fDoKinematicsPlots)


         if(fRunPlotInvMass)
        {
           TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
           InvMass->Fill(Sum12.Mag(),Q3);
        }//if(fRunPlotInvMass)



        }
      }
    }
  }
}//AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionSE2ME1

//==================================================================================================================================================

bool AliAnalysisTaskThreeBodyProtonPrimary::DeltaEtaDeltaPhi(int species1, int species2,
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   int part1PDGcode,
                                                   int part2PDGcode,
                                                   bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist,
                                                   AliFemtoDreamCollConfig Config) {
  // DoThisPair = ij where i is the number of daughters for first particle, j for the second

  static const float piHi = TMath::Pi();
  //auto fDeltaPhiSqMax = Config.GetDeltaPhiMax() * Config.GetDeltaPhiMax();
  //auto fDeltaEtaSqMax = Config.GetDeltaEtaMax() * Config.GetDeltaEtaMax() ;

  double DeltaPhiSqMaxValue = 0;
  double DeltaEtaSqMaxValue = 0;


  //cout<<part1.GetPDGCode()<<endl;

  
  if ((species1==2&&species2==0)||(species1==0&&species2==2)||
      (species1==3&&species2==1)||(species1==1&&species2==3))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPPrim*fDeltaPhiMaxPPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPPrim*fDeltaEtaMaxPPrim;

  }
  else if ((species1==0&&species2==0)||(species1==1&&species2==1))
  {
    //fDeltaPhiSqMax = 0.017*0.017;
    //fDeltaEtaSqMax = 0.017*0.017;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPP*fDeltaPhiMaxPP;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPP*fDeltaEtaMaxPP;

  }
  else if ((species1==3&&species2==0)||(species1==0&&species2==3)||
      (species1==2&&species2==1)||(species1==1&&species2==2))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPAPrim*fDeltaPhiMaxPAPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPAPrim*fDeltaEtaMaxPAPrim;

  }

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
        if ((dphiAvg / (float) size) * (dphiAvg / (float) size) / DeltaPhiSqMaxValue
            + deta * deta / DeltaEtaSqMaxValue < 1.) {
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

bool AliAnalysisTaskThreeBodyProtonPrimary::DeltaEtaDeltaPhi(int species1, int species2,
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   int part1PDGcode,
                                                   int part2PDGcode,
                                                   bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist,
                                                   AliFemtoDreamCollConfig Config, double Q3) {
  // DoThisPair = ij where i is the number of daughters for first particle, j for the second

  static const float piHi = TMath::Pi();
  //auto fDeltaPhiSqMax = Config.GetDeltaPhiMax() * Config.GetDeltaPhiMax();
  //auto fDeltaEtaSqMax = Config.GetDeltaEtaMax() * Config.GetDeltaEtaMax() ;

  double DeltaPhiSqMaxValue = 0;
  double DeltaEtaSqMaxValue = 0;

  //cout<<part1PDGcode<<endl;

 //bool test = false;

 if ((species1==2&&species2==0)||(species1==0&&species2==2)||
      (species1==3&&species2==1)||(species1==1&&species2==3))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;

    DeltaPhiSqMaxValue = fDeltaPhiMaxPPrim*fDeltaPhiMaxPPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPPrim*fDeltaEtaMaxPPrim;
  //  test = true;
  }
  else if ((species1==0&&species2==0)||(species1==1&&species2==1))
  {
    //fDeltaPhiSqMax = 0.017*0.017;
    //fDeltaEtaSqMax = 0.017*0.017;

    DeltaPhiSqMaxValue = fDeltaPhiMaxPP*fDeltaPhiMaxPP;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPP*fDeltaEtaMaxPP;

  //  test = true;
  }
  else if ((species1==3&&species2==0)||(species1==0&&species2==3)||
      (species1==2&&species2==1)||(species1==1&&species2==2))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPAPrim*fDeltaPhiMaxPAPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPAPrim*fDeltaEtaMaxPAPrim;

  }

  //cout<<species1<<"   "<<species2<<"  "<<test<<endl;

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
        if ((dphiAvg / (float) size) * (dphiAvg / (float) size) / DeltaPhiSqMaxValue
            + deta * deta / DeltaEtaSqMaxValue < 1.) {
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

void AliAnalysisTaskThreeBodyProtonPrimary::FillPairTransverseMass(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector,
                        int firstSpecies, int secondSpecies, TH1F* hist1, std::vector<int> PDGCodes, TH2F* hist2) { // ADDED BY RAFFA

  // Two Body Same Event

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;

  // Get the PID codes std::vector<int>
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;

  // Get particle masses
  auto massPart1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massPart2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();


//  TVector3 momPart1 = part1.GetMomentum();
//  TVector3 momPart2 = part2.GetMomentum();
//  TLorentzVector track1, track2;


  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;

  if(abs(*itPDGPar1)==3122) DaughterPart1 = 2;
  if(abs(*itPDGPar1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGPar1)==321) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==3122) DaughterPart2 = 2;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGPar2)==321) DaughterPart2 = 1;

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
        Particle1_LV.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(), massPart1);
        Particle2_LV.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(), massPart2);
        // Get momentum
        float RelativeMomentum = AliFemtoDreamHigherPairMath::RelativePairMomentum(Particle1_LV, Particle2_LV);


        TLorentzVector trackSum_LV = Particle1_LV + Particle2_LV;

        Double_t pairMT = TMath::Sqrt(pow(0.5*trackSum_LV.Pt(),2.) + pow(0.5*(massPart1+massPart2),2.));

        Double_t Energy1 = sqrt(pow(Particle1_LV.P(),2)+pow(Particle1_LV.M(),2));
        Double_t Energy2 = sqrt(pow(Particle2_LV.P(),2)+pow(Particle2_LV.M(),2));
        Double_t MinvByHand = sqrt(pow(Energy1+Energy2,2) - pow(Particle1_LV.P(),2) - pow(Particle2_LV.P(),2) - 2.*(Particle1_LV.Px()*Particle2_LV.Px() + Particle1_LV.Py()*Particle2_LV.Py() + Particle1_LV.Pz()*Particle2_LV.Pz()));

//        std::cout << "iPart1 " << "  Mass: " << Particle1_LV.M() << "   Mass2: " << massPart1 << "   Momentum: " << Particle1_LV.P() << "  Energy: " << Energy1 << std::endl;
//        std::cout << "iPart2 " << "  Mass: " << Particle2_LV.M() << "   Mass2: " << massPart2 << "   Momentum: " << Particle2_LV.P() << "  Energy: " << Energy2 << std::endl;
//        std::cout << "iPart1-iPart2: " << "   Pt: " << trackSum_LV.Pt() << "  Minv: " << trackSum_LV.M() << "  MinvByHand: " << MinvByHand << std::endl;

        hist2->Fill(RelativeMomentum, pairMT);
        hist1->Fill(pairMT);


//        bool Pair12 = true;

//        if(!fturnoffClosePairRejectionCompletely){

//            Pair12 = DeltaEtaDeltaPhi(*iPart1,*iPart2,true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[7+phiEtaHistNo],Config); //GANESHA  Check about DeltaEtaDeltaPhi Function. Might need to change for 2 Body
//        }

//        if(!Pair12) {continue;}

    }
  }

}
