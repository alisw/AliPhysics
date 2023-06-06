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
      fMixedEventMult(nullptr),
      fSameEventMultDphi(nullptr),
      fMixedEventMultDphi(nullptr),
      fSameEventmT(nullptr),
      fMixedEventmT(nullptr),
      fSameEventPtPrimaries(nullptr),
      fMixedEventPtPrimaries(nullptr),
      fSameEventPtPrimaries2(nullptr),
      fMixedEventPtPrimaries2(nullptr),
      fSameEventPtProtons(nullptr),
      fMixedEventPtProtons(nullptr),
      fSameEventPtProtons2(nullptr),
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
      fP1Histos(nullptr),
      fKHistos(nullptr),
      fInvMassVSmTList(nullptr),

      fRunThreeBody(true),
      fRunPlotInvMass(false),
      fRunPlotQ3Vsq(true),
      fRunPlotPhiTheta(true),
      fRunPlotOtherHistos(true),
      fRunPlotMult(true),
      fRunCorrDeltaPhi(false),
      fRunPlotPt(false),
      fRunOfficialTwoBody(false), // ADDED BY RAFFA
      fDoKinematicsPlots(false),
      fPlotP1(false),
      fCutElectrons(false),
      fRunPlotInvMassVSmT(false),

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
      fQ3cutValue(1.),
      fQ3MinValue(0.),
      fCleanWithLambdas(false),
      fDoOnlyThreeBody(true),
      fStandardMixing(true),
      fPlotsMC(false), //needs to go somewhere else 
      fRunmTPlots(false), 

      fSameEventTripletArray(nullptr),
      fSameEventTripletPtPrimaries(nullptr),
      fSameEventTripletPtProtons(nullptr),
      fSameEventTripletPtPrimaries2(nullptr),
      fSameEventTripletPtProtons2(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletMultArray12(nullptr),
      fSameEventTripletMultArray23(nullptr),
      fSameEventTripletMultArray31(nullptr),
      fSameEventTripletmTArray12(nullptr),
      fSameEventTripletmTArray23(nullptr),
      fSameEventTripletmTArray31(nullptr),
      f2Same1MixedEventTripletmTArray12(nullptr),
      f2Same1MixedEventTripletmTArray23(nullptr),
      f2Same1MixedEventTripletmTArray31(nullptr),
      fSameEventTripletPhiThetaArray_SamePair(nullptr), 
      fSameEventTripletPhiThetaArray_DifferentPair(nullptr),
      fSameEventTripletPtvsQ3Primaries(nullptr),
      fSameEventTripletPtvsQ3Protons(nullptr),

      fSameEventTripletArray_TwoBody(nullptr),
      fSameEventMultArray_TwoBody(nullptr),
      fSameEventDphiArray_TwoBody(nullptr),
      fSameEventTripletMultArray_TwoBody(nullptr),
      fSameEventTripletPhiThetaArray_TwoBody(nullptr),
      fSameEventP1(nullptr),
      fSameEventK(nullptr),

      fPairTranverseMass_TwoBody(nullptr), // ADDED BY RAFFA
      fPairTranverseMassVSkstar_TwoBody(nullptr), // ADDED BY RAFFA
      fKinematics(nullptr),
      fPrimAngles(nullptr),
      fKinematicsME(nullptr),
      fPrimAnglesME(nullptr),
      fDeta(nullptr),
      fDphi(nullptr),
      fDetaSEvsME(nullptr),
      fDphiSEvsME(nullptr),
      fDetaME(nullptr),
      fDphiME(nullptr),

      fPartContainer(0),
      fPartContainerPPP(0),
      fPartContainerPPPrim(0),
      fPartContainerPPAPrim(0),
      fPartContainerPP(0),
      fPartContainerPPrim(0),
      fPartContainerPAPrim(0),
      fVectPartContainers(0),

      fMixedEventTripletArray(nullptr),
      fMixedEventTripletPtPrimaries(nullptr),
      fMixedEventTripletPtProtons(nullptr),
      fMixedEventTripletPtPrimaries2(nullptr),
      fMixedEventTripletPtProtons2(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletMultArray12(nullptr),
      fMixedEventTripletMultArray23(nullptr),
      fMixedEventTripletMultArray31(nullptr),
      fMixedEventTripletmTArray12(nullptr),
      fMixedEventTripletmTArray23(nullptr),
      fMixedEventTripletmTArray31(nullptr),
      fMixedEventTripletPhiThetaArray_SamePair(nullptr), 
      fMixedEventTripletPhiThetaArray_DifferentPair(nullptr),
      fMixedEventTripletPtvsQ3Primaries(nullptr),
      fMixedEventTripletPtvsQ3Protons(nullptr),
      fhistDeltaPhi12(nullptr),
      fhistDeltaPhi23(nullptr),
      fhistDeltaPhi31(nullptr),
      fhistDeltaPhiDeltaEta12(nullptr),
      fhistDeltaPhiDeltaEta23(nullptr),
      fhistDeltaPhiDeltaEta31(nullptr),
      fhistDeltaPhi12ME(nullptr),
      fhistDeltaPhi23ME(nullptr),
      fhistDeltaPhi31ME(nullptr),
      fhistDeltaPhiDeltaEta12ME(nullptr),
      fhistDeltaPhiDeltaEta23ME(nullptr),
      fhistDeltaPhiDeltaEta31ME(nullptr),
      
      fMixedEventTripletArray_TwoBody(nullptr),
      fMixedEventTripletMultArray_TwoBody(nullptr),
      fMixedEventMultArray_TwoBody(nullptr),
      fMixedEventDphiArray_TwoBody(nullptr),
      fMixedEventTripletPhiThetaArray_TwoBody(nullptr),
      fMixedEventP1(nullptr),
      fMixedEventK(nullptr),

      fQ3VskDistributionsArrayq12(nullptr),
      fQ3VskDistributionsArrayq12Mixed(nullptr),
      fQ3VskDistributionsArrayq23(nullptr),
      fQ3VskDistributionsArrayq23Mixed(nullptr),
      fDoubletVsTrippletPPP(nullptr),
      fInvMass12(nullptr),
      fInvMass23(nullptr),
      fInvMass31(nullptr),
      fInvMassVSmT12(nullptr),
      fInvMassVSmT23(nullptr),
      fInvMassVSmT31(nullptr),
      fInvMassVSmT12MixedEvent(nullptr),
      fInvMassVSmT23MixedEvent(nullptr),
      fInvMassVSmT31MixedEvent(nullptr),

      fpTvsEtaTrueKaons(nullptr),
      fpTvsEtaTrueAntiKaons(nullptr),
      fpTvsEtaTrueProtons(nullptr),
      fpTvsEtaTrueAntiProtons(nullptr),
      fpTvsEtaRecoKaons(nullptr),
      fpTvsEtaRecoAntiKaons(nullptr),
      fpTvsEtaRecoProtons(nullptr),
      fpTvsEtaRecoAntiProtons(nullptr),

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
      fSameEventMultDphi(nullptr),
      fMixedEventMultDphi(nullptr),
      fSameEventmT(nullptr),
      fMixedEventmT(nullptr),
      fSameEventPtPrimaries(nullptr),
      fMixedEventPtPrimaries(nullptr),
      fSameEventPtPrimaries2(nullptr),
      fMixedEventPtPrimaries2(nullptr),
      fSameEventPtProtons(nullptr),
      fMixedEventPtProtons(nullptr),
      fSameEventPtProtons2(nullptr),
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
      fP1Histos(nullptr),
      fKHistos(nullptr),
      fInvMassVSmTList(nullptr),

      fRunThreeBody(true),
      fRunPlotInvMass(false),
      fRunPlotQ3Vsq(true),
      fRunPlotPhiTheta(true),
      fRunPlotOtherHistos(true),
      fRunPlotMult(true),
      fRunCorrDeltaPhi(false),
      fRunPlotPt(false),
      fRunOfficialTwoBody(false), // ADDED BY RAFFA
      fDoKinematicsPlots(false),
      fPlotP1(false),
      fCutElectrons(false),
      fRunPlotInvMassVSmT(false),

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
      fQ3cutValue(1.),
      fQ3MinValue(0.),
      fCleanWithLambdas(false),
      fDoOnlyThreeBody(true),
      fStandardMixing(true),
      fPlotsMC(false), //needs to go somewhere else 
      fRunmTPlots(false), 

      fSameEventTripletArray(nullptr),
      fSameEventTripletPtPrimaries(nullptr),
      fSameEventTripletPtProtons(nullptr),
      fSameEventTripletPtPrimaries2(nullptr),
      fSameEventTripletPtProtons2(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletMultArray12(nullptr),
      fSameEventTripletMultArray23(nullptr),
      fSameEventTripletMultArray31(nullptr),
      fSameEventTripletmTArray12(nullptr),
      fSameEventTripletmTArray23(nullptr),
      fSameEventTripletmTArray31(nullptr),
      f2Same1MixedEventTripletmTArray12(nullptr),
      f2Same1MixedEventTripletmTArray23(nullptr),
      f2Same1MixedEventTripletmTArray31(nullptr),
      fSameEventTripletPhiThetaArray_SamePair(nullptr), 
      fSameEventTripletPhiThetaArray_DifferentPair(nullptr),
      fSameEventTripletPtvsQ3Primaries(nullptr),
      fSameEventTripletPtvsQ3Protons(nullptr),

      fSameEventTripletArray_TwoBody(nullptr),
      fSameEventMultArray_TwoBody(nullptr),
      fSameEventDphiArray_TwoBody(nullptr),
      fSameEventTripletMultArray_TwoBody(nullptr),
      fSameEventTripletPhiThetaArray_TwoBody(nullptr),
      fSameEventP1(nullptr),
      fSameEventK(nullptr),

      fPairTranverseMass_TwoBody(nullptr), // ADDED BY RAFFA
      fPairTranverseMassVSkstar_TwoBody(nullptr), // ADDED BY RAFFA
      fKinematics(nullptr),
      fPrimAngles(nullptr),
      fKinematicsME(nullptr),
      fPrimAnglesME(nullptr),
      fDeta(nullptr),
      fDphi(nullptr),
      fDetaSEvsME(nullptr),
      fDphiSEvsME(nullptr),
      fDetaME(nullptr),
      fDphiME(nullptr),

      fPartContainer(0),
      fPartContainerPPP(0),
      fPartContainerPPPrim(0),
      fPartContainerPPAPrim(0),
      fPartContainerPP(0),
      fPartContainerPPrim(0),
      fPartContainerPAPrim(0),
      fVectPartContainers(0),

      fMixedEventTripletArray(nullptr),
      fMixedEventTripletPtPrimaries(nullptr),
      fMixedEventTripletPtProtons(nullptr),
      fMixedEventTripletPtPrimaries2(nullptr),
      fMixedEventTripletPtProtons2(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletMultArray12(nullptr),
      fMixedEventTripletMultArray23(nullptr),
      fMixedEventTripletMultArray31(nullptr),
      fMixedEventTripletmTArray12(nullptr),
      fMixedEventTripletmTArray23(nullptr),
      fMixedEventTripletmTArray31(nullptr),
      fMixedEventTripletPhiThetaArray_SamePair(nullptr), 
      fMixedEventTripletPhiThetaArray_DifferentPair(nullptr),
      fMixedEventTripletPtvsQ3Primaries(nullptr),
      fMixedEventTripletPtvsQ3Protons(nullptr),
      fhistDeltaPhi12(nullptr),
      fhistDeltaPhi23(nullptr),
      fhistDeltaPhi31(nullptr),
      fhistDeltaPhiDeltaEta12(nullptr),
      fhistDeltaPhiDeltaEta23(nullptr),
      fhistDeltaPhiDeltaEta31(nullptr),
      fhistDeltaPhi12ME(nullptr),
      fhistDeltaPhi23ME(nullptr),
      fhistDeltaPhi31ME(nullptr),
      fhistDeltaPhiDeltaEta12ME(nullptr),
      fhistDeltaPhiDeltaEta23ME(nullptr),
      fhistDeltaPhiDeltaEta31ME(nullptr),
      
      fMixedEventTripletArray_TwoBody(nullptr),
      fMixedEventTripletMultArray_TwoBody(nullptr),
      fMixedEventMultArray_TwoBody(nullptr),
      fMixedEventDphiArray_TwoBody(nullptr),
      fMixedEventTripletPhiThetaArray_TwoBody(nullptr),
      fMixedEventP1(nullptr),
      fMixedEventK(nullptr),

      fQ3VskDistributionsArrayq12(nullptr),
      fQ3VskDistributionsArrayq12Mixed(nullptr),
      fQ3VskDistributionsArrayq23(nullptr),
      fQ3VskDistributionsArrayq23Mixed(nullptr),
      fDoubletVsTrippletPPP(nullptr),
      fInvMass12(nullptr),
      fInvMass23(nullptr),
      fInvMass31(nullptr),
      fInvMassVSmT12(nullptr),
      fInvMassVSmT23(nullptr),
      fInvMassVSmT31(nullptr),
      fInvMassVSmT12MixedEvent(nullptr),
      fInvMassVSmT23MixedEvent(nullptr),
      fInvMassVSmT31MixedEvent(nullptr),

      fpTvsEtaTrueKaons(nullptr),
      fpTvsEtaTrueAntiKaons(nullptr),
      fpTvsEtaTrueProtons(nullptr),
      fpTvsEtaTrueAntiProtons(nullptr),
      fpTvsEtaRecoKaons(nullptr),
      fpTvsEtaRecoAntiKaons(nullptr),
      fpTvsEtaRecoProtons(nullptr),
      fpTvsEtaRecoAntiProtons(nullptr),

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

    fP1Histos = new TList();
    fP1Histos->SetOwner();
    fP1Histos->SetName("P1Histos");

    fSameEventP1 = new TH2F*[6];
    TString histTitlesSameP1[6] = {"sameEventDistributionPPPrim","sameEventDistributionAPAPAPrim", "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPPAPrim","sameEventDistributionAPAPPrim"};

    fMixedEventP1 = new TH2F*[6];
    TString histTitlesMixedP1[6] = {"mixedEventDistributionPPPrim","mixedEventDistributionAPAPAPrim", "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPPAPrim","mixedEventDistributionAPAPPrim"};

    if(fPlotP1){
      for (int i = 0; i < 6; ++i) {
        fSameEventP1[i] =  new TH2F(histTitlesSameP1[i],histTitlesSameP1[i], 2000, 0, 2,4000,0,4);
        fMixedEventP1[i] =  new TH2F(histTitlesMixedP1[i],histTitlesMixedP1[i], 2000, 0, 2,4000,0,4);

        fP1Histos->Add(fSameEventP1[i]);
        fP1Histos->Add(fMixedEventP1[i]);

      }

      fResultsThreeBody->Add(fP1Histos);
    }

    fKHistos = new TList();
    fKHistos->SetOwner();
    fKHistos->SetName("KHistos");

    fSameEventK = new TH2F*[6];
    TString histTitlesSameK[6] = {"sameEventDistributionPPPrim","sameEventDistributionAPAPAPrim", "sameEventDistributionPPP", "sameEventDistributionAPAPAP", "sameEventDistributionPPAPrim","sameEventDistributionAPAPPrim"};

    fMixedEventK = new TH2F*[6];
    TString histTitlesMixedK[6] = {"mixedEventDistributionPPPrim","mixedEventDistributionAPAPAPrim", "mixedEventDistributionPPP", "mixedEventDistributionAPAPAP", "mixedEventDistributionPPAPrim","mixedEventDistributionAPAPPrim"};

    if(fPlotP1){
      for (int i = 0; i < 6; ++i) {
        fSameEventK[i] =  new TH2F(histTitlesSameK[i],histTitlesSameK[i], 2000, 0, 2,4000,0,4);
        fMixedEventK[i] =  new TH2F(histTitlesMixedK[i],histTitlesMixedK[i], 2000, 0, 2,4000,0,4);

        fKHistos->Add(fSameEventK[i]);
        fKHistos->Add(fMixedEventK[i]);

      }

      fResultsThreeBody->Add(fKHistos);
    }



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
    fSameEventMultDphi = new TList();
    fSameEventMultDphi->SetOwner();
    fSameEventMultDphi->SetName("SameEventMultDphi");

    fSameEventmT = new TList();
    fSameEventmT->SetOwner();
    fSameEventmT->SetName("fSameEventmT");

    fSameEventTripletMultArray = new TH2F*[28];
    TString histTitlesSameMult[28] = {"sameEventDistributionMultPPPrim","sameEventDistributionMultAPAPAPrim", "sameEventDistributionMultPPP", "sameEventDistributionMultAPAPAP", "sameEventDistributionMultPPAPrim","sameEventDistributionMultAPAPPrim", "sameEventDistributionMultPPrimPrim", "sameEventDistributionMultAPAPrimAPrim", "sameEventDistributionMultPAPrimAPrim", "sameEventDistributionMultAPPrimPrim","sameEventDistributionMultPPSamePrimMixed", "sameEventDistributionMultPPrimSamePMixed", "sameEventDistributionMultAPAPSameAPrimMixed", "sameEventDistributionMultAPAPrimSameAPMixed", "sameEventDistributionMultPPSamePMixed", "sameEventDistributionMultAPAPSameAPMixed", "sameEventDistributionMultPPSameAPrimMixed", "sameEventDistributionMultPAPrimSamePMixed", "sameEventDistributionMultAPAPSamePrimMixed", "sameEventDistributionMultAPPrimSameAPMixed", "sameEventDistributionMultPPrimSamePrimMixed", "sameEventDistributionMultPrimPrimSamePMixed", "sameEventDistributionMultAPAPrimSameAPrimMixed", "sameEventDistributionMultAPrimAPrimSameAPMixed", "sameEventDistributionMultPAPrimSameAPrimMixed", "sameEventDistributionMultAPrimAPrimSamePMixed", "sameEventDistributionMultAPPrimSamePrimMixed", "sameEventDistributionMultPrimPrimSameAPMixed"};

    f2Same1MixedEventTripletmTArray12 = new TH2F*[18];
    TString histTitles2Same1MixedmT12[18] = { "sameEventDistributionmTPPSamePrimMixed12", "sameEventDistributionmTPPrimSamePMixed12", "sameEventDistributionmTAPAPSameAPrimMixed12", "sameEventDistributionmTAPAPrimSameAPMixed12", "sameEventDistributionmTPPSamePMixed12", "sameEventDistributionmTAPAPSameAPMixed12", "sameEventDistributionmTPPSameAPrimMixed12", "sameEventDistributionmTPAPrimSamePMixed12", "sameEventDistributionmTAPAPSamePrimMixed12", "sameEventDistributionmTAPPrimSameAPMixed12", "sameEventDistributionmTPPrimSamePrimMixed12", "sameEventDistributionmTPrimPrimSamePMixed12", "sameEventDistributionmTAPAPrimSameAPrimMixed12", "sameEventDistributionmTAPrimAPrimSameAPMixed12", "sameEventDistributionmTPAPrimSameAPrimMixed12", "sameEventDistributionmTAPrimAPrimSamePMixed12", "sameEventDistributionmTAPPrimSamePrimMixed12", "sameEventDistributionmTPrimPrimSameAPMixed12"};

    f2Same1MixedEventTripletmTArray23 = new TH2F*[18];
    TString histTitles2Same1MixedmT23[18] = { "sameEventDistributionmTPPSamePrimMixed23", "sameEventDistributionmTPPrimSamePMixed23", "sameEventDistributionmTAPAPSameAPrimMixed23", "sameEventDistributionmTAPAPrimSameAPMixed23", "sameEventDistributionmTPPSamePMixed23", "sameEventDistributionmTAPAPSameAPMixed23", "sameEventDistributionmTPPSameAPrimMixed23", "sameEventDistributionmTPAPrimSamePMixed23", "sameEventDistributionmTAPAPSamePrimMixed23", "sameEventDistributionmTAPPrimSameAPMixed23", "sameEventDistributionmTPPrimSamePrimMixed23", "sameEventDistributionmTPrimPrimSamePMixed23", "sameEventDistributionmTAPAPrimSameAPrimMixed23", "sameEventDistributionmTAPrimAPrimSameAPMixed23", "sameEventDistributionmTPAPrimSameAPrimMixed23", "sameEventDistributionmTAPrimAPrimSamePMixed23", "sameEventDistributionmTAPPrimSamePrimMixed23", "sameEventDistributionmTPrimPrimSameAPMixed23"};

    f2Same1MixedEventTripletmTArray31 = new TH2F*[18];
    TString histTitles2Same1MixedmT31[18] = { "sameEventDistributionmTPPSamePrimMixed31", "sameEventDistributionmTPPrimSamePMixed31", "sameEventDistributionmTAPAPSameAPrimMixed31", "sameEventDistributionmTAPAPrimSameAPMixed31", "sameEventDistributionmTPPSamePMixed31", "sameEventDistributionmTAPAPSameAPMixed31", "sameEventDistributionmTPPSameAPrimMixed31", "sameEventDistributionmTPAPrimSamePMixed31", "sameEventDistributionmTAPAPSamePrimMixed31", "sameEventDistributionmTAPPrimSameAPMixed31", "sameEventDistributionmTPPrimSamePrimMixed31", "sameEventDistributionmTPrimPrimSamePMixed31", "sameEventDistributionmTAPAPrimSameAPrimMixed31", "sameEventDistributionmTAPrimAPrimSameAPMixed31", "sameEventDistributionmTPAPrimSameAPrimMixed31", "sameEventDistributionmTAPrimAPrimSamePMixed31", "sameEventDistributionmTAPPrimSamePrimMixed31", "sameEventDistributionmTPrimPrimSameAPMixed31"};


    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 28; ++i) {
	      fSameEventTripletMultArray[i] =  new TH2F(histTitlesSameMult[i],histTitlesSameMult[i], 8000, 0, 8.,26,1,27);
          if(fRunPlotMult){fSameEventMult->Add(fSameEventTripletMultArray[i]);}
	     }

       for (int i = 0; i < 18; ++i) {
        f2Same1MixedEventTripletmTArray12[i] =  new TH2F(histTitles2Same1MixedmT12[i], histTitles2Same1MixedmT12[i], 4000, 0, 8., 100, 0., 5.);
          if(fRunmTPlots && i< 10){fSameEventmT->Add(f2Same1MixedEventTripletmTArray12[i]);}
        f2Same1MixedEventTripletmTArray23[i] =  new TH2F(histTitles2Same1MixedmT23[i], histTitles2Same1MixedmT23[i], 4000, 0, 8., 100, 0., 5.);
          if(fRunmTPlots && i< 10){fSameEventmT->Add(f2Same1MixedEventTripletmTArray23[i]);}
        f2Same1MixedEventTripletmTArray31[i] =  new TH2F(histTitles2Same1MixedmT31[i], histTitles2Same1MixedmT31[i], 4000, 0, 8., 100, 0., 5.);
          if(fRunmTPlots && i< 10){fSameEventmT->Add(f2Same1MixedEventTripletmTArray31[i]);}
	     }
    }

    fSameEventTripletMultArray12 = new TH2F*[10];
    TString histTitlesSameMult12[10] = {"sameEventDistributionMultPPPrim12","sameEventDistributionMultAPAPAPrim12","sameEventDistributionMultPPP12", "sameEventDistributionMultAPAPAP12", "sameEventDistributionMultPPAPrim12","sameEventDistributionMultAPAPPrim12", "sameEventDistributionMultPPrimPrim12", "sameEventDistributionMultAPAPrimAPrim12", "sameEventDistributionMultPAPrimAPrim12", "sameEventDistributionMultAPPrimPrim12"};

    fSameEventTripletMultArray23 = new TH2F*[10];
    TString histTitlesSameMult23[10] = {"sameEventDistributionMultPPPrim23","sameEventDistributionMultAPAPAPrim23","sameEventDistributionMultPPP23", "sameEventDistributionMultAPAPAP23", "sameEventDistributionMultPPAPrim23","sameEventDistributionMultAPAPPrim23", "sameEventDistributionMultPPrimPrim23", "sameEventDistributionMultAPAPrimAPrim23", "sameEventDistributionMultPAPrimAPrim23", "sameEventDistributionMultAPPrimPrim23"};

    fSameEventTripletMultArray31 = new TH2F*[10];
    TString histTitlesSameMult31[10] = {"sameEventDistributionMultPPPrim31","sameEventDistributionMultAPAPAPrim31","sameEventDistributionMultPPP31", "sameEventDistributionMultAPAPAP31", "sameEventDistributionMultPPAPrim31","sameEventDistributionMultAPAPPrim31", "sameEventDistributionMultPPrimPrim31", "sameEventDistributionMultAPAPrimAPrim31", "sameEventDistributionMultPAPrimAPrim31", "sameEventDistributionMultAPPrimPrim31"};

    //To-Do at one point: clean up: last 4 histograms not needed 
    //Q3 vs mT of a pair
    fSameEventTripletmTArray12 = new TH2F*[10];
    TString histTitlesSamemT12[10] = {"sameEventDistributionmTPPPrim12","sameEventDistributionmTAPAPAPrim12","sameEventDistributionmTPPP12", "sameEventDistributionmTAPAPAP12", "sameEventDistributionmTPPAPrim12","sameEventDistributionmTAPAPPrim12", "sameEventDistributionmTPPrimPrim12", "sameEventDistributionmTAPAPrimAPrim12", "sameEventDistributionmTPAPrimAPrim12", "sameEventDistributionmTAPPrimPrim12"};

    fSameEventTripletmTArray23 = new TH2F*[10];
    TString histTitlesSamemT23[10] = {"sameEventDistributionmTPPPrim23","sameEventDistributionmTAPAPAPrim23","sameEventDistributionmTPPP23", "sameEventDistributionmTAPAPAP23", "sameEventDistributionmTPPAPrim23","sameEventDistributionmTAPAPPrim23", "sameEventDistributionmTPPrimPrim23", "sameEventDistributionmTAPAPrimAPrim23", "sameEventDistributionmTPAPrimAPrim23", "sameEventDistributionmTAPPrimPrim23"};

    fSameEventTripletmTArray31 = new TH2F*[10];
    TString histTitlesSamemT31[10] = {"sameEventDistributionmTPPPrim31","sameEventDistributionmTAPAPAPrim31","sameEventDistributionmTPPP31", "sameEventDistributionmTAPAPAP31", "sameEventDistributionmTPPAPrim31","sameEventDistributionmTAPAPPrim31", "sameEventDistributionmTPPrimPrim31", "sameEventDistributionmTAPAPrimAPrim31", "sameEventDistributionmTPAPrimAPrim31", "sameEventDistributionmTAPPrimPrim31"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 10; ++i) {
        fSameEventTripletMultArray12[i] = new TH2F(histTitlesSameMult12[i],histTitlesSameMult12[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fSameEventMult->Add(fSameEventTripletMultArray12[i]);}
        fSameEventTripletMultArray23[i] = new TH2F(histTitlesSameMult23[i],histTitlesSameMult23[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fSameEventMult->Add(fSameEventTripletMultArray23[i]);}
        fSameEventTripletMultArray31[i] = new TH2F(histTitlesSameMult31[i],histTitlesSameMult31[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fSameEventMult->Add(fSameEventTripletMultArray31[i]);}

        fSameEventTripletmTArray12[i] = new TH2F(histTitlesSamemT12[i],histTitlesSamemT12[i], 4000, 0, 8, 100, 0., 5.); 
          if(fRunmTPlots && i<6){fSameEventmT->Add(fSameEventTripletmTArray12[i]);}
        fSameEventTripletmTArray23[i] = new TH2F(histTitlesSamemT23[i],histTitlesSamemT23[i], 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fSameEventmT->Add(fSameEventTripletmTArray23[i]);}
        fSameEventTripletmTArray31[i] = new TH2F(histTitlesSamemT31[i],histTitlesSamemT31[i], 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fSameEventmT->Add(fSameEventTripletmTArray31[i]);}
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
    // Same event Pt Primaries Second Track
    fSameEventPtPrimaries2 = new TList();
    fSameEventPtPrimaries2->SetOwner();
    fSameEventPtPrimaries2->SetName("SameEventPtPrimaries2");


    fSameEventTripletPtPrimaries2 = new TH1F*[28];
    TString histTitlesSamePtPrimaries2[28] = {"sameEventDistributionPtPrimariesPPPrim","sameEventDistributionPtPrimariesAPAPAPrim", "sameEventDistributionPtPrimariesPPP", "sameEventDistributionPtPrimariesAPAPAP", "sameEventDistributionPtPrimariesPPAPrim","sameEventDistributionPtPrimariesAPAPPrim", "sameEventDistributionPtPrimariesPPrimPrim", "sameEventDistributionPtPrimariesAPAPrimAPrim", "sameEventDistributionPtPrimariesPAPrimAPrim", "sameEventDistributionPtPrimariesAPPrimPrim","sameEventDistributionPtPrimariesPPSamePrimMixed", "sameEventDistributionPtPrimariesPPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPSameAPrimMixed", "sameEventDistributionPtPrimariesAPAPrimSameAPMixed", "sameEventDistributionPtPrimariesPPSamePMixed", "sameEventDistributionPtPrimariesAPAPSameAPMixed", "sameEventDistributionPtPrimariesPPSameAPrimMixed", "sameEventDistributionPtPrimariesPAPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPSamePrimMixed", "sameEventDistributionPtPrimariesAPPrimSameAPMixed", "sameEventDistributionPtPrimariesPPrimSamePrimMixed", "sameEventDistributionPtPrimariesPrimPrimSamePMixed", "sameEventDistributionPtPrimariesAPAPrimSameAPrimMixed", "sameEventDistributionPtPrimariesAPrimAPrimSameAPMixed", "sameEventDistributionPtPrimariesPAPrimSameAPrimMixed", "sameEventDistributionPtPrimariesAPrimAPrimSamePMixed", "sameEventDistributionPtPrimariesAPPrimSamePrimMixed", "sameEventDistributionPtPrimariesPrimPrimSameAPMixed"};

    if(fDoOnlyThreeBody){
      for (int i = 0; i < 28; ++i) {
        fSameEventTripletPtPrimaries2[i] =  new TH1F(histTitlesSamePtPrimaries2[i],histTitlesSamePtPrimaries2[i], 4000, 0, 4);
              if(fRunPlotPt){fSameEventPtPrimaries2->Add(fSameEventTripletPtPrimaries2[i]);}
       }
    }

    //...............................................................................................
    // Same event Pt Protons Second Track
    fSameEventPtProtons2 = new TList();
    fSameEventPtProtons2->SetOwner();
    fSameEventPtProtons2->SetName("SameEventPtProtons2");


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
    fSameEventMultArray_TwoBody = new TH2F*[7];
    fSameEventDphiArray_TwoBody = new TH2F*[7];
    TString histTitlesSameDphi_TwoBody[7] = {"sameEventDistributionDphiPP","sameEventDistributionDphiAPAP", "sameEventDistributionDphiPPrim", "sameEventDistributionDphiAPAPrim", "sameEventDistributionDphiPAPrim", "sameEventDistributionDphiAPPrim", "sameEventDistributionDphiPrimAPrim"};

    fPairTranverseMassVSkstar_TwoBody = new TH2F*[7]; 
    TString histTitlesPairTranverseMassVSkstar_TwoBody[7] = {"PairTranverseMassVSkstarDistributionPP","PairTranverseMassVSkstarDistributionAPAP", "PairTranverseMassVSkstarDistributionPPrim", "PairTranverseMassVSkstarDistributionAPAPrim", "PairTranverseMassVSkstarDistributionPAPrim", "PairTranverseMassVSkstarDistributionAPPrim", "PairTranverseMassVSkstarDistributionPrimAPrim"}; 

    if(!fDoOnlyThreeBody){
	   for (int i = 0; i < 7; ++i) {
	      fSameEventTripletMultArray_TwoBody[i] =  new TH2F(histTitlesSameMult_TwoBody[i],histTitlesSameMult_TwoBody[i],8000, 0, 8,26,1,27);
        fSameEventMultArray_TwoBody[i] =  new TH2F(histTitlesSameMult_TwoBody[i],histTitlesSameMult_TwoBody[i],8000, 0, 8,26,1,27);
        fSameEventDphiArray_TwoBody[i] =  new TH2F(histTitlesSameDphi_TwoBody[i],histTitlesSameDphi_TwoBody[i],8000, 0, 8,200,-10,10);
          if(fRunPlotMult) fSameEventMult->Add(fSameEventTripletMultArray_TwoBody[i]);
          if(fRunCorrDeltaPhi){
            fSameEventMultDphi->Add(fSameEventMultArray_TwoBody[i]);
            fSameEventMultDphi->Add(fSameEventDphiArray_TwoBody[i]);
          }
	      if(fRunOfficialTwoBody){
             fPairTranverseMassVSkstar_TwoBody[i] =  new TH2F(histTitlesPairTranverseMassVSkstar_TwoBody[i],histTitlesPairTranverseMassVSkstar_TwoBody[i], 800, 0, 8, 800, 0, 8); 
             fSameEvent->Add(fPairTranverseMassVSkstar_TwoBody[i]); 
        }
      }
    }


    fhistDeltaPhi12 = new TH2F*[6];
    fhistDeltaPhi23 = new TH2F*[6];
    fhistDeltaPhi31 = new TH2F*[6];
    TString histDeltaPhi12[6] = {"sameEventDistributionDphiPPPrim_PPrim12","sameEventDistributionDphiAPAPAPrim_APAPrim12","sameEventDistributionDphiPPP_PP12", "sameEventDistributionDphiAPAPAP_APAP12", "sameEventDistributionDphiPPAPrim_PAPrim12","sameEventDistributionDphiAPAPPrim_APPrim12"};
    TString histDeltaPhi23[6] = {"sameEventDistributionDphiPPPrim_PP23","sameEventDistributionDphiAPAPAPrim_APAP23","sameEventDistributionDphiPPP_PP23", "sameEventDistributionDphiAPAPAP_APAP23", "sameEventDistributionDphiPPAPrim_PP23","sameEventDistributionDphiAPAPPrim_APAP23"};
    TString histDeltaPhi31[6] = {"sameEventDistributionDphiPPPrim_PPrim31","sameEventDistributionDphiAPAPAPrim_APAPrim31","sameEventDistributionDphiPPP_PP31", "sameEventDistributionDphiAPAPAP_APAP31", "sameEventDistributionDphiPPAPrim_PAPrim31","sameEventDistributionDphiAPAPPrim_APPrim31"};
    fhistDeltaPhiDeltaEta12 = new TH2F*[6];
    fhistDeltaPhiDeltaEta23 = new TH2F*[6];
    fhistDeltaPhiDeltaEta31 = new TH2F*[6];
    TString histDeltaPhiDeltaEta12[6] = {"sameEventDistributionDphiDetaPPPrim_PPrim12","sameEventDistributionDphiDetaAPAPAPrim_APAPrim12","sameEventDistributionDphiDetaPPP_PP12", "sameEventDistributionDphiDetaAPAPAP_APAP12", "sameEventDistributionDphiDetaPPAPrim_PAPrim12","sameEventDistributionDphiDetaAPAPPrim_APPrim12"};
    TString histDeltaPhiDeltaEta23[6] = {"sameEventDistributionDphiDetaPPPrim_PP23","sameEventDistributionDphiDetaAPAPAPrim_APAP23","sameEventDistributionDphiDetaPPP_PP23", "sameEventDistributionDphiDetaAPAPAP_APAP23", "sameEventDistributionDphiDetaPPAPrim_PP23","sameEventDistributionDphiDetaAPAPPrim_APAP23"};
    TString histDeltaPhiDeltaEta31[6] = {"sameEventDistributionDphiDetaPPPrim_PPrim31","sameEventDistributionDphiDetaAPAPAPrim_APAPrim31","sameEventDistributionDphiDetaPPP_PP31", "sameEventDistributionDphiDetaAPAPAP_APAP31", "sameEventDistributionDphiDetaPPAPrim_PAPrim31","sameEventDistributionDphiDetaAPAPPrim_APPrim31"};

    if(fDoOnlyThreeBody){
       for (int i = 0; i < 6; ++i) {
        fhistDeltaPhi12[i] =  new TH2F(histDeltaPhi12[i],histDeltaPhi12[i], 8000, 0, 8,200,-10,10);
        fhistDeltaPhi23[i] =  new TH2F(histDeltaPhi23[i],histDeltaPhi23[i], 8000, 0, 8,200,-10,10);
        fhistDeltaPhi31[i] =  new TH2F(histDeltaPhi31[i],histDeltaPhi31[i], 8000, 0, 8,200,-10,10);
        fhistDeltaPhiDeltaEta12[i] =  new TH2F(histDeltaPhiDeltaEta12[i],histDeltaPhiDeltaEta12[i], 200,-10,10, 200,-10,10);
        fhistDeltaPhiDeltaEta23[i] =  new TH2F(histDeltaPhiDeltaEta23[i],histDeltaPhiDeltaEta23[i], 200,-10,10, 200,-10,10);
        fhistDeltaPhiDeltaEta31[i] =  new TH2F(histDeltaPhiDeltaEta31[i],histDeltaPhiDeltaEta31[i], 200,-10,10, 200,-10,10);        
            if(fRunCorrDeltaPhi){ 
              fSameEventMultDphi->Add(fhistDeltaPhi12[i]); 
              fSameEventMultDphi->Add(fhistDeltaPhi23[i]); 
              fSameEventMultDphi->Add(fhistDeltaPhi31[i]); 
              fSameEventMultDphi->Add(fhistDeltaPhiDeltaEta12[i]); 
              fSameEventMultDphi->Add(fhistDeltaPhiDeltaEta23[i]); 
              fSameEventMultDphi->Add(fhistDeltaPhiDeltaEta31[i]); 
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
    fMixedEventMultDphi = new TList();
    fMixedEventMultDphi->SetOwner();
    fMixedEventMultDphi->SetName("MixedEventMultDphi");

    fMixedEventmT = new TList();
    fMixedEventmT->SetOwner();
    fMixedEventmT->SetName("fMixedEventmT");

    fMixedEventTripletMultArray = new TH2F*[10];
    TString histTitlesMixedMult[10] = {"mixedEventDistributionMultPPPrim","mixedEventDistributionMultAPAPAPrim","mixedEventDistributionMultPPP", "mixedEventDistributionMultAPAPAP", "mixedEventDistributionMultPPAPrim","mixedEventDistributionMultAPAPPrim", "mixedEventDistributionMultPPrimPrim", "mixedEventDistributionMultAPAPrimAPrim", "mixedEventDistributionMultPAPrimAPrim", "mixedEventDistributionMultAPPrimPrim"};

    fMixedEventTripletMultArray12 = new TH2F*[10];
    TString histTitlesMixedMult12[10] = {"mixedEventDistributionMultPPPrim12","mixedEventDistributionMultAPAPAPrim12","mixedEventDistributionMultPPP12", "mixedEventDistributionMultAPAPAP12", "mixedEventDistributionMultPPAPrim12","mixedEventDistributionMultAPAPPrim12", "mixedEventDistributionMultPPrimPrim12", "mixedEventDistributionMultAPAPrimAPrim12", "mixedEventDistributionMultPAPrimAPrim12", "mixedEventDistributionMultAPPrimPrim12"};

    fMixedEventTripletMultArray23 = new TH2F*[10];
    TString histTitlesMixedMult23[10] = {"mixedEventDistributionMultPPPrim23","mixedEventDistributionMultAPAPAPrim23","mixedEventDistributionMultPPP23", "mixedEventDistributionMultAPAPAP23", "mixedEventDistributionMultPPAPrim23","mixedEventDistributionMultAPAPPrim23", "mixedEventDistributionMultPPrimPrim23", "mixedEventDistributionMultAPAPrimAPrim23", "mixedEventDistributionMultPAPrimAPrim23", "mixedEventDistributionMultAPPrimPrim23"};

    fMixedEventTripletMultArray31 = new TH2F*[10];
    TString histTitlesMixedMult31[10] = {"mixedEventDistributionMultPPPrim31","mixedEventDistributionMultAPAPAPrim31","mixedEventDistributionMultPPP31", "mixedEventDistributionMultAPAPAP31", "mixedEventDistributionMultPPAPrim31","mixedEventDistributionMultAPAPPrim31", "mixedEventDistributionMultPPrimPrim31", "mixedEventDistributionMultAPAPrimAPrim31", "mixedEventDistributionMultPAPrimAPrim31", "mixedEventDistributionMultAPPrimPrim31"};

    //Mixed Event Q3 vs mT Pair
    fMixedEventTripletmTArray12 = new TH2F*[10];
    TString histTitlesMixedmT12[10] = {"mixedEventDistributionmTPPPrim12","mixedEventDistributionmTAPAPAPrim12","mixedEventDistributionmTPPP12", "mixedEventDistributionmTAPAPAP12", "mixedEventDistributionmTPPAPrim12","mixedEventDistributionmTAPAPPrim12", "mixedEventDistributionmTPPrimPrim12", "mixedEventDistributionmTAPAPrimAPrim12", "mixedEventDistributionmTPAPrimAPrim12", "mixedEventDistributionmTAPPrimPrim12"};
 
    fMixedEventTripletmTArray23 = new TH2F*[10];
    TString histTitlesMixedmT23[10] = {"mixedEventDistributionmTPPPrim23","mixedEventDistributionmTAPAPAPrim23","mixedEventDistributionmTPPP23", "mixedEventDistributionmTAPAPAP23", "mixedEventDistributionmTPPAPrim23","mixedEventDistributionmTAPAPPrim23", "mixedEventDistributionmTPPrimPrim23", "mixedEventDistributionmTAPAPrimAPrim23", "mixedEventDistributionmTPAPrimAPrim23", "mixedEventDistributionmTAPPrimPrim23"};

    fMixedEventTripletmTArray31 = new TH2F*[10];
    TString histTitlesMixedmT31[10] = {"mixedEventDistributionmTPPPrim31","mixedEventDistributionmTAPAPAPrim31","mixedEventDistributionmTPPP31", "mixedEventDistributionmTAPAPAP31", "mixedEventDistributionmTPPAPrim31","mixedEventDistributionmTAPAPPrim31", "mixedEventDistributionmTPPrimPrim31", "mixedEventDistributionmTAPAPrimAPrim31", "mixedEventDistributionmTPAPrimAPrim31", "mixedEventDistributionmTAPPrimPrim31"};


    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 10; ++i) {
	      fMixedEventTripletMultArray[i] = new TH2F(histTitlesMixedMult[i],histTitlesMixedMult[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fMixedEventMult->Add(fMixedEventTripletMultArray[i]);}
        fMixedEventTripletMultArray12[i] = new TH2F(histTitlesMixedMult12[i],histTitlesMixedMult12[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fMixedEventMult->Add(fMixedEventTripletMultArray12[i]);}
        fMixedEventTripletMultArray23[i] = new TH2F(histTitlesMixedMult23[i],histTitlesMixedMult23[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fMixedEventMult->Add(fMixedEventTripletMultArray23[i]);}
        fMixedEventTripletMultArray31[i] = new TH2F(histTitlesMixedMult31[i],histTitlesMixedMult31[i], 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fMixedEventMult->Add(fMixedEventTripletMultArray31[i]);}

        fMixedEventTripletmTArray12[i] = new TH2F(histTitlesMixedmT12[i],histTitlesMixedmT12[i], 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fMixedEventmT->Add(fMixedEventTripletmTArray12[i]);} 
        fMixedEventTripletmTArray23[i] = new TH2F(histTitlesMixedmT23[i],histTitlesMixedmT23[i], 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fMixedEventmT->Add(fMixedEventTripletmTArray23[i]);}
        fMixedEventTripletmTArray31[i] = new TH2F(histTitlesMixedmT31[i],histTitlesMixedmT31[i], 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fMixedEventmT->Add(fMixedEventTripletmTArray31[i]);}
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
    fMixedEventPtPrimaries2->SetName("MixedEventPtPrimaries2");

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
    fMixedEventPtProtons2->SetName("MixedEventPtProtons2");


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
    fMixedEventMultArray_TwoBody = new TH2F*[7];
    fMixedEventDphiArray_TwoBody = new TH2F*[7];
    TString histTitlesMixedDphi_TwoBody[7] = {"MixedEventDistributionDphiPP","MixedEventDistributionDphiAPAP", "MixedEventDistributionDphiPPrim", "MixedEventDistributionDphiAPAPrim", "MixedEventDistributionDphiPAPrim", "MixedEventDistributionDphiAPPrim", "MixedEventDistributionDphiPrimAPrim"};

    if(!fDoOnlyThreeBody){
	     for (int i = 0; i < 7; ++i) {
	      fMixedEventTripletMultArray_TwoBody[i] =  new TH2F(histTitlesMixedMult_TwoBody[i],histTitlesMixedMult_TwoBody[i], 8000, 0, 8,26,1,27);
        fMixedEventMultArray_TwoBody[i] =  new TH2F(histTitlesMixedMult_TwoBody[i],histTitlesMixedMult_TwoBody[i], 8000, 0, 8,26,1,27);
        fMixedEventDphiArray_TwoBody[i] =  new TH2F(histTitlesMixedDphi_TwoBody[i],histTitlesMixedDphi_TwoBody[i], 8000, 0, 8,200,-10,10);
             if(fRunPlotMult){ fMixedEventMult->Add(fMixedEventTripletMultArray_TwoBody[i]); }
             if(fRunCorrDeltaPhi){ 
              fMixedEventMultDphi->Add(fMixedEventMultArray_TwoBody[i]); 
              fMixedEventMultDphi->Add(fMixedEventDphiArray_TwoBody[i]); 
            }
	     }
    }

    fhistDeltaPhi12ME = new TH2F*[6];
    fhistDeltaPhi23ME = new TH2F*[6];
    fhistDeltaPhi31ME = new TH2F*[6];
    TString histDeltaPhi12ME[6] = {"mixedEventDistributionDphiPPPrim_PPrim12","mixedEventDistributionDphiAPAPAPrim_APAPrim12","mixedEventDistributionDphiPPP_PP12", "mixedEventDistributionDphiAPAPAP_APAP12", "mixedEventDistributionDphiPPAPrim_PAPrim12","mixedEventDistributionDphiAPAPPrim_APPrim12"};
    TString histDeltaPhi23ME[6] = {"mixedEventDistributionDphiPPPrim_PP23","mixedEventDistributionDphiAPAPAPrim_APAP23","mixedEventDistributionDphiPPP_PP23", "mixedEventDistributionDphiAPAPAP_APAP23", "mixedEventDistributionDphiPPAPrim_PP23","mixedEventDistributionDphiAPAPPrim_APAP23"};
    TString histDeltaPhi31ME[6] = {"mixedEventDistributionDphiPPPrim_PPrim31","mixedEventDistributionDphiAPAPAPrim_APAPrim31","mixedEventDistributionDphiPPP_PP31", "mixedEventDistributionDphiAPAPAP_APAP31", "mixedEventDistributionDphiPPAPrim_PAPrim31","mixedEventDistributionDphiAPAPPrim_APPrim31"};
    fhistDeltaPhiDeltaEta12ME = new TH2F*[6];
    fhistDeltaPhiDeltaEta23ME = new TH2F*[6];
    fhistDeltaPhiDeltaEta31ME = new TH2F*[6];
    TString histDeltaPhiDeltaEta12ME[6] = {"mixedEventDistributionDphiDetaPPPrim_PPrim12","mixedEventDistributionDphiDetaAPAPAPrim_APAPrim12","mixedEventDistributionDphiDetaPPP_PP12", "mixedEventDistributionDphiDetaAPAPAP_APAP12", "mixedEventDistributionDphiDetaPPAPrim_PAPrim12","mixedEventDistributionDphiDetaAPAPPrim_APPrim12"};
    TString histDeltaPhiDeltaEta23ME[6] = {"mixedEventDistributionDphiDetaPPPrim_PP23","mixedEventDistributionDphiDetaAPAPAPrim_APAP23","mixedEventDistributionDphiDetaPPP_PP23", "mixedEventDistributionDphiDetaAPAPAP_APAP23", "mixedEventDistributionDphiDetaPPAPrim_PP23","mixedEventDistributionDphiDetaAPAPPrim_APAP23"};
    TString histDeltaPhiDeltaEta31ME[6] = {"mixedEventDistributionDphiDetaPPPrim_PPrim31","mixedEventDistributionDphiDetaAPAPAPrim_APAPrim31","mixedEventDistributionDphiDetaPPP_PP31", "mixedEventDistributionDphiDetaAPAPAP_APAP31", "mixedEventDistributionDphiDetaPPAPrim_PAPrim31","mixedEventDistributionDphiDetaAPAPPrim_APPrim31"};

    if(fDoOnlyThreeBody){
       for (int i = 0; i < 6; ++i) {
        fhistDeltaPhi12ME[i] =  new TH2F(histDeltaPhi12ME[i],histDeltaPhi12ME[i], 8000, 0, 8,200,-10,10);
        fhistDeltaPhi23ME[i] =  new TH2F(histDeltaPhi23ME[i],histDeltaPhi23ME[i], 8000, 0, 8,200,-10,10);
        fhistDeltaPhi31ME[i] =  new TH2F(histDeltaPhi31ME[i],histDeltaPhi31ME[i], 8000, 0, 8,200,-10,10);
        fhistDeltaPhiDeltaEta12ME[i] =  new TH2F(histDeltaPhiDeltaEta12ME[i],histDeltaPhiDeltaEta12ME[i], 200,-10,10, 200,-10,10);
        fhistDeltaPhiDeltaEta23ME[i] =  new TH2F(histDeltaPhiDeltaEta23ME[i],histDeltaPhiDeltaEta23ME[i], 200,-10,10, 200,-10,10);
        fhistDeltaPhiDeltaEta31ME[i] =  new TH2F(histDeltaPhiDeltaEta31ME[i],histDeltaPhiDeltaEta31ME[i], 200,-10,10, 200,-10,10);        
            if(fRunCorrDeltaPhi){ 
              fMixedEventMultDphi->Add(fhistDeltaPhi12ME[i]); 
              fMixedEventMultDphi->Add(fhistDeltaPhi23ME[i]); 
              fMixedEventMultDphi->Add(fhistDeltaPhi31ME[i]); 
              fMixedEventMultDphi->Add(fhistDeltaPhiDeltaEta12ME[i]); 
              fMixedEventMultDphi->Add(fhistDeltaPhiDeltaEta23ME[i]); 
              fMixedEventMultDphi->Add(fhistDeltaPhiDeltaEta31ME[i]); 
            }
       }
    }




    if(fRunPlotMult){
    fResultsThreeBody->Add(fSameEventMult);
    fResultsThreeBody->Add(fMixedEventMult);
    }else{
    fResultsThreeBody->Add(fSameEvent);
    fResultsThreeBody->Add(fMixedEvent);
    }
    if(fRunmTPlots){
      fResultsThreeBody->Add(fSameEventmT); //GANEHSA dont forget to add this
      fResultsThreeBody->Add(fMixedEventmT);
    }
    if(fRunCorrDeltaPhi){
    fResultsThreeBody->Add(fSameEventMultDphi);
    fResultsThreeBody->Add(fMixedEventMultDphi);
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


    fDetaSEvsME = new TH2F*[10];
    TString histTitlesfDetaSEvsME[10] = {"SEvsMEDetaPPSamePrimMixed", "SEvsMEDetaPPrimSamePMixed",
        "SEvsMEDetaAPAPSameAPrimMixed", "SEvsMEDetaAPAPrimSameAPMixed",
        "SEvsMEDetaPPSamePMixed", "SEvsMEDetaAPAPSameAPMixed",
        "SEvsMEDetaPPSameAPrimMixed", "SEvsMEDetaPAPrimSamePMixed",
        "SEvsMEDetaAPAPSamePrimMixed", "SEvsMEDetaAPPrimSameAPMixed"};

    fDphiSEvsME = new TH2F*[10];
    TString histTitlesfDphiSEvsME[10] = {"SEvsMEDphiPPSamePrimMixed", "SEvsMEDphiPPrimSamePMixed",
        "SEvsMEDphiAPAPSameAPrimMixed", "SEvsMEDphiAPAPrimSameAPMixed",
        "SEvsMEDphiPPSamePMixed", "SEvsMEDphiAPAPSameAPMixed",
        "SEvsMEDphiPPSameAPrimMixed", "SEvsMEDphiPAPrimSamePMixed",
        "SEvsMEDphiAPAPSamePrimMixed", "SEvsMEDphiAPPrimSameAPMixed"};

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
      for(int i=0;i<10;i++){

          fDetaSEvsME[i] = new TH2F(histTitlesfDetaSEvsME[i],histTitlesfDetaSEvsME[i], 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());
          fKinematicsPlots->Add(fDetaSEvsME[i]);

          fDphiSEvsME[i] = new TH2F(histTitlesfDphiSEvsME[i],histTitlesfDphiSEvsME[i], 360, -TMath::Pi(), TMath::Pi(), 360, -TMath::Pi(), TMath::Pi());
          fKinematicsPlots->Add(fDphiSEvsME[i]);

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

      fQ3VskDistributionsArrayq12 =  new TH2F*[16]; 
      TString histTitlesfQ3VskDistributions[16] =  {"Q3vskDistributionPPPrim","Q3vskDistributionAPAPAPrim",
        "Q3vskDistributionPPP", "Q3vskDistributionAPAPAP", "Q3vskDistributionPPAPrim","Q3vskDistributionAPAPPrim",
        "Q3vskDistributionPPSamePrimMixed", "Q3vskDistributionPPrimSamePMixed",
        "Q3vskDistributionAPAPSameAPrimMixed","Q3vskDistributionAPAPrimSameaAPMixed", "Q3vskDistributionPPSamePMixed",
        "Q3vskDistributionAPAPSameAPMixed", "Q3vskDistributionPPSameAPrimMixed", "Q3vskDistributionPAPrimSamePMixed",
        "Q3vskDistributionAPAPSamePrimMixed", "Q3vskDistributionAPPrimSameAPMixed"};

        //,"Q3vskDistributionLLSameLMixed", "Q3vskDistributionaLaLSameaLMixed", "TRASH"};
    
      for(int i=0;i<16;i++){
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

      fQ3VskDistributionsArrayq23 =  new TH2F*[16];
      for(int i=0;i<16;i++){
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
    fInvMass12 = new TH2F*[16];
    TString histTitlesInvMass12[16] = {"InvMass12_PPPrim","InvMass12_APAPAPrim", "InvMass12_PPP", "InvMass12_APAPAP", "InvMass12_PPAPrim","InvMass12_APAPPrim", "InvMass12_PPSamePrimMixed", "InvMass12_PPrimSamePMixed", "InvMass12_APAPSameAPrimMixed", "InvMass12_APAPrimSameAPMixed", "InvMass12_PPSamePMixed", "InvMass12_APAPSameAPMixed", "InvMass12_PPSameAPrimMixed", "InvMass12_PAPrimSamePMixed", "InvMass12_APAPSamePrimMixed", "InvMass12_APPrimSameAPMixed"};
    fInvMass23 = new TH2F*[16];
    TString histTitlesInvMass23[16] = {"InvMass23_PPPrim","InvMass23_APAPAPrim", "InvMass23_PPP", "InvMass23_APAPAP", "InvMass23_PPAPrim","InvMass23_APAPPrim", "InvMass23_PPSamePrimMixed", "InvMass23_PPrimSamePMixed", "InvMass23_APAPSameAPrimMixed", "InvMass23_APAPrimSameAPMixed", "InvMass23_PPSamePMixed", "InvMass23_APAPSameAPMixed", "InvMass23_PPSameAPrimMixed", "InvMass23_PAPrimSamePMixed", "InvMass23_APAPSamePrimMixed", "InvMass23_APPrimSameAPMixed"};
    fInvMass31 = new TH2F*[16];
    TString histTitlesInvMass31[16] = {"InvMass31_PPPrim","InvMass31_APAPAPrim", "InvMass31_PPP", "InvMass31_APAPAP", "InvMass31_PPAPrim","InvMass31_APAPPrim", "InvMass31_PPSamePrimMixed", "InvMass31_PPrimSamePMixed", "InvMass31_APAPSameAPrimMixed", "InvMass31_APAPrimSameAPMixed", "InvMass31_PPSamePMixed", "InvMass31_APAPSameAPMixed", "InvMass31_PPSameAPrimMixed", "InvMass31_PAPrimSamePMixed", "InvMass31_APAPSamePrimMixed", "InvMass31_APPrimSameAPMixed"};


    if (fRunPlotInvMass){
      fInvMassList = new TList();
      fInvMassList->SetOwner();
      fInvMassList->SetName("InvMass");


      for(int i=0;i<16;i++){
        fInvMass12[i] = new TH2F(histTitlesInvMass12[i],histTitlesInvMass12[i], 500, 0., 5., 12, 0., 1.2);
        fInvMassList->Add(fInvMass12[i]);
        fInvMass23[i] = new TH2F(histTitlesInvMass23[i],histTitlesInvMass23[i], 500, 0., 5., 12, 0., 1.2);
        fInvMassList->Add(fInvMass23[i]);
        fInvMass31[i] = new TH2F(histTitlesInvMass31[i],histTitlesInvMass31[i], 500, 0., 5., 12, 0., 1.2);
        fInvMassList->Add(fInvMass31[i]);
      }

      fResultsThreeBody->Add(fInvMassList);
     
    }

    //Invariant Mass vs mT
    fInvMassVSmT12 = new TH2F*[16];
    TString histTitlesInvMassVSmT12[16] = {"InvMassVSmT12_PPPrim","InvMassVSmT12_APAPAPrim", "InvMassVSmT12_PPP", "InvMassVSmT12_APAPAP", "InvMassVSmT12_PPAPrim","InvMassVSmT12_APAPPrim", "InvMassVSmT12_PPSamePrimMixed", "InvMassVSmT12_PPrimSamePMixed", "InvMassVSmT12_APAPSameAPrimMixed", "InvMassVSmT12_APAPrimSameAPMixed", "InvMassVSmT12_PPSamePMixed", "InvMassVSmT12_APAPSameAPMixed", "InvMassVSmT12_PPSameAPrimMixed", "InvMassVSmT12_PAPrimSamePMixed", "InvMassVSmT12_APAPSamePrimMixed", "InvMassVSmT12_APPrimSameAPMixed"};
    fInvMassVSmT23 = new TH2F*[16];
    TString histTitlesInvMassVSmT23[16] = {"InvMassVSmT23_PPPrim","InvMassVSmT23_APAPAPrim", "InvMassVSmT23_PPP", "InvMassVSmT23_APAPAP", "InvMassVSmT23_PPAPrim","InvMassVSmT23_APAPPrim", "InvMassVSmT23_PPSamePrimMixed", "InvMassVSmT23_PPrimSamePMixed", "InvMassVSmT23_APAPSameAPrimMixed", "InvMassVSmT23_APAPrimSameAPMixed", "InvMassVSmT23_PPSamePMixed", "InvMassVSmT23_APAPSameAPMixed", "InvMassVSmT23_PPSameAPrimMixed", "InvMassVSmT23_PAPrimSamePMixed", "InvMassVSmT23_APAPSamePrimMixed", "InvMassVSmT23_APPrimSameAPMixed"};
    fInvMassVSmT31 = new TH2F*[16];
    TString histTitlesInvMassVSmT31[16] = {"InvMassVSmT31_PPPrim","InvMassVSmT31_APAPAPrim", "InvMassVSmT31_PPP", "InvMassVSmT31_APAPAP", "InvMassVSmT31_PPAPrim","InvMassVSmT31_APAPPrim", "InvMassVSmT31_PPSamePrimMixed", "InvMassVSmT31_PPrimSamePMixed", "InvMassVSmT31_APAPSameAPrimMixed", "InvMassVSmT31_APAPrimSameAPMixed", "InvMassVSmT31_PPSamePMixed", "InvMassVSmT31_APAPSameAPMixed", "InvMassVSmT31_PPSameAPrimMixed", "InvMassVSmT31_PAPrimSamePMixed", "InvMassVSmT31_APAPSamePrimMixed", "InvMassVSmT31_APPrimSameAPMixed"};

    fInvMassVSmT12MixedEvent = new TH2F*[6];
    TString histTitlesInvMassVSmT12MixedEvent[6] = {"InvMassVSmT12MixedEvent_PPPrim","InvMassVSmT12MixedEvent_APAPAPrim", "InvMassVSmT12MixedEvent_PPP", "InvMassVSmT12MixedEvent_APAPAP", "InvMassVSmT12MixedEvent_PPAPrim","InvMassVSmT12MixedEvent_APAPPrim"};
    fInvMassVSmT23MixedEvent = new TH2F*[6];
    TString histTitlesInvMassVSmT23MixedEvent[6] = {"InvMassVSmT23MixedEvent_PPPrim","InvMassVSmT23MixedEvent_APAPAPrim", "InvMassVSmT23MixedEvent_PPP", "InvMassVSmT23MixedEvent_APAPAP", "InvMassVSmT23MixedEvent_PPAPrim","InvMassVSmT23MixedEvent_APAPPrim"};
    fInvMassVSmT31MixedEvent = new TH2F*[6];
    TString histTitlesInvMassVSmT31MixedEvent[6] = {"InvMassVSmT31MixedEvent_PPPrim","InvMassVSmT31MixedEvent_APAPAPrim", "InvMassVSmT31MixedEvent_PPP", "InvMassVSmT31MixedEvent_APAPAP", "InvMassVSmT31MixedEvent_PPAPrim","InvMassVSmT31MixedEvent_APAPPrim"};



    if (fRunPlotInvMassVSmT){
      fInvMassVSmTList = new TList();
      fInvMassVSmTList->SetOwner();
      fInvMassVSmTList->SetName("InvMassVSmT");


      for(int i=0;i<16;i++){
        fInvMassVSmT12[i] = new TH2F(histTitlesInvMassVSmT12[i],histTitlesInvMassVSmT12[i], 500, 0., 5., 100, 0., 5.);
        fInvMassVSmTList->Add(fInvMassVSmT12[i]);
        fInvMassVSmT23[i] = new TH2F(histTitlesInvMassVSmT23[i],histTitlesInvMassVSmT23[i], 500, 0., 5., 100, 0., 5.);
        fInvMassVSmTList->Add(fInvMassVSmT23[i]);
        fInvMassVSmT31[i] = new TH2F(histTitlesInvMassVSmT31[i],histTitlesInvMassVSmT31[i], 500, 0., 5., 100, 0., 5.);
        fInvMassVSmTList->Add(fInvMassVSmT31[i]);
        if(i<6){
          fInvMassVSmT12MixedEvent[i] = new TH2F(histTitlesInvMassVSmT12MixedEvent[i],histTitlesInvMassVSmT12MixedEvent[i], 500, 0., 5., 100, 0., 5.);
          fInvMassVSmTList->Add(fInvMassVSmT12MixedEvent[i]);
          fInvMassVSmT23MixedEvent[i] = new TH2F(histTitlesInvMassVSmT23MixedEvent[i],histTitlesInvMassVSmT23MixedEvent[i], 500, 0., 5., 100, 0., 5.);
          fInvMassVSmTList->Add(fInvMassVSmT23MixedEvent[i]);
          fInvMassVSmT31MixedEvent[i] = new TH2F(histTitlesInvMassVSmT31MixedEvent[i],histTitlesInvMassVSmT31MixedEvent[i], 500, 0., 5., 100, 0., 5.);
          fInvMassVSmTList->Add(fInvMassVSmT31MixedEvent[i]);
        }
      }

      fResultsThreeBody->Add(fInvMassVSmTList);
     
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

  double pTMinP = 0.6;
  double pTMaxP = 0.8;
  double pTMinPrim = 0.3;
  double pTMaxPrim = 0.4;

  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) { // TO DO: think about double track selection
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent);

    if(fCutElectrons){ // this is a strct cut suggested by Ramona to reject electrons

      if (fProton->isSelected(fTrack)&&(fTrack->GetPt()<pTMinP||fTrack->GetPt()>pTMaxP)) {
        Protons.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoProtons->Fill(track->Pt(),track->Eta());
      }
      if (fAntiProton->isSelected(fTrack)&&(fTrack->GetPt()<pTMinP||fTrack->GetPt()>pTMaxP)) {
        AntiProtons.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoAntiProtons->Fill(track->Pt(),track->Eta());
      }
      if (fPrimary->isSelected(fTrack)&&(fTrack->GetPt()<pTMinPrim||fTrack->GetPt()>pTMaxPrim)) {
        Primaries.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoKaons->Fill(track->Pt(),track->Eta());
      }
      if (fAntiPrimary->isSelected(fTrack)&&(fTrack->GetPt()<pTMinPrim||fTrack->GetPt()>pTMaxPrim)) {
        AntiPrimaries.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoAntiKaons->Fill(track->Pt(),track->Eta());
      }

    }else{

      if (fProton->isSelected(fTrack)) {
        Protons.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoProtons->Fill(track->Pt(),track->Eta());
      }
      if (fAntiProton->isSelected(fTrack)) {
        AntiProtons.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoAntiProtons->Fill(track->Pt(),track->Eta());
      }
      if (fPrimary->isSelected(fTrack)) {
        Primaries.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoKaons->Fill(track->Pt(),track->Eta());
      }
      if (fAntiPrimary->isSelected(fTrack)) {
        AntiPrimaries.push_back(*fTrack);
        if (fPlotsMC) fpTvsEtaRecoAntiKaons->Fill(track->Pt(),track->Eta());
      }

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
       FillTripletDistribution( ParticleVector, 2, 0, 0, fSameEventTripletArray[0],PDGCodes, bins[1],fSameEventTripletMultArray[0], fSameEventTripletMultArray12[0] ,fSameEventTripletMultArray23[0], fSameEventTripletMultArray31[0], fSameEventTripletPtPrimaries[0], fSameEventTripletPtProtons[0], fSameEventTripletPtPrimaries2[0], fSameEventTripletPtProtons2[0], fSameEventTripletPtvsQ3Primaries[0], fSameEventTripletPtvsQ3Protons[0], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,0, *fConfig, fQ3VskDistributionsArrayq12[0],fQ3VskDistributionsArrayq23[0], fKinematics[0], fPrimAngles[0], fDeta[0], fDphi[0], fInvMass12[0], fInvMass23[0], fInvMass31[0], fhistDeltaPhi12[0], fhistDeltaPhi23[0], fhistDeltaPhi31[0], fhistDeltaPhiDeltaEta12[0], fhistDeltaPhiDeltaEta23[0], fhistDeltaPhiDeltaEta31[0], fSameEventP1[0], fSameEventK[0], fSameEventTripletmTArray12[0], fSameEventTripletmTArray23[0], fSameEventTripletmTArray31[0], fInvMassVSmT12[0], fInvMassVSmT23[0], fInvMassVSmT31[0]);
       // Antiproton Antiproton Antiprimary
       FillTripletDistribution( ParticleVector, 3, 1, 1, fSameEventTripletArray[1],PDGCodes, bins[1],fSameEventTripletMultArray[1], fSameEventTripletMultArray12[1] ,fSameEventTripletMultArray23[1], fSameEventTripletMultArray31[1], fSameEventTripletPtPrimaries[1], fSameEventTripletPtProtons[1], fSameEventTripletPtPrimaries2[1], fSameEventTripletPtProtons2[1], fSameEventTripletPtvsQ3Primaries[1], fSameEventTripletPtvsQ3Protons[1], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,1, *fConfig, fQ3VskDistributionsArrayq12[1],fQ3VskDistributionsArrayq23[1], fKinematics[1], fPrimAngles[1], fDeta[1], fDphi[1], fInvMass12[1], fInvMass23[1], fInvMass31[1], fhistDeltaPhi12[1], fhistDeltaPhi23[1], fhistDeltaPhi31[1], fhistDeltaPhiDeltaEta12[1], fhistDeltaPhiDeltaEta23[1], fhistDeltaPhiDeltaEta31[1], fSameEventP1[1], fSameEventK[1], fSameEventTripletmTArray12[1], fSameEventTripletmTArray23[1], fSameEventTripletmTArray31[1], fInvMassVSmT12[1], fInvMassVSmT23[1], fInvMassVSmT31[1]);
       // Proton Proton Proton
       FillTripletDistribution( ParticleVector, 0, 0, 0, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletMultArray12[2] ,fSameEventTripletMultArray23[2], fSameEventTripletMultArray31[2], fSameEventTripletPtPrimaries[2], fSameEventTripletPtProtons[2], fSameEventTripletPtPrimaries2[2], fSameEventTripletPtProtons2[2], fSameEventTripletPtvsQ3Primaries[2], fSameEventTripletPtvsQ3Protons[2], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,2, *fConfig, fQ3VskDistributionsArrayq12[2],fQ3VskDistributionsArrayq23[2], fKinematics[2], fPrimAngles[2], fDeta[2], fDphi[2], fInvMass12[2], fInvMass23[2], fInvMass31[2], fhistDeltaPhi12[2], fhistDeltaPhi23[2], fhistDeltaPhi31[2], fhistDeltaPhiDeltaEta12[2], fhistDeltaPhiDeltaEta23[2], fhistDeltaPhiDeltaEta31[2], fSameEventP1[2], fSameEventK[2], fSameEventTripletmTArray12[2], fSameEventTripletmTArray23[2], fSameEventTripletmTArray31[2], fInvMassVSmT12[2], fInvMassVSmT23[2], fInvMassVSmT31[2]);
//       FillTripletDistribution( ParticleVector, 0, 0, 1, fSameEventTripletArray[2],PDGCodes, bins[1],fSameEventTripletMultArray[2], fSameEventTripletMultArray12[2] ,fSameEventTripletMultArray23[2], fSameEventTripletMultArray31[2], fSameEventTripletPtPrimaries[2], fSameEventTripletPtProtons[2], fSameEventTripletPtPrimaries2[2], fSameEventTripletPtProtons2[2], fSameEventTripletPtvsQ3Primaries[2], fSameEventTripletPtvsQ3Protons[2], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,2, *fConfig, fQ3VskDistributionsArrayq12[2],fQ3VskDistributionsArrayq23[2], fKinematics[2], fPrimAngles[2], fDeta[2], fDphi[2], fInvMass12[2], fInvMass23[2], fInvMass31[2], fhistDeltaPhi12[2], fhistDeltaPhi23[2], fhistDeltaPhi31[2], fhistDeltaPhiDeltaEta12[2], fhistDeltaPhiDeltaEta23[2], fhistDeltaPhiDeltaEta31[2], fSameEventP1[2], fSameEventK[2]);
       // Antiproton Antiproton Antiproton
       FillTripletDistribution( ParticleVector, 1, 1, 1, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletMultArray12[3] ,fSameEventTripletMultArray23[3], fSameEventTripletMultArray31[3], fSameEventTripletPtPrimaries[3], fSameEventTripletPtProtons[3], fSameEventTripletPtPrimaries2[3], fSameEventTripletPtProtons2[3], fSameEventTripletPtvsQ3Primaries[3], fSameEventTripletPtvsQ3Protons[3], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,3, *fConfig, fQ3VskDistributionsArrayq12[3],fQ3VskDistributionsArrayq23[3], fKinematics[3], fPrimAngles[3], fDeta[3], fDphi[3], fInvMass12[3], fInvMass23[3], fInvMass31[3], fhistDeltaPhi12[3], fhistDeltaPhi23[3], fhistDeltaPhi31[3], fhistDeltaPhiDeltaEta12[3], fhistDeltaPhiDeltaEta23[3], fhistDeltaPhiDeltaEta31[3], fSameEventP1[3], fSameEventK[3], fSameEventTripletmTArray12[3], fSameEventTripletmTArray23[3], fSameEventTripletmTArray31[3], fInvMassVSmT12[3], fInvMassVSmT23[3], fInvMassVSmT31[3]);
//       FillTripletDistribution( ParticleVector, 1, 1, 0, fSameEventTripletArray[3],PDGCodes, bins[1],fSameEventTripletMultArray[3], fSameEventTripletMultArray12[3] ,fSameEventTripletMultArray23[3], fSameEventTripletMultArray31[3], fSameEventTripletPtPrimaries[3], fSameEventTripletPtProtons[3], fSameEventTripletPtPrimaries2[3], fSameEventTripletPtProtons2[3], fSameEventTripletPtvsQ3Primaries[3], fSameEventTripletPtvsQ3Protons[3], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,3, *fConfig, fQ3VskDistributionsArrayq12[3],fQ3VskDistributionsArrayq23[3], fKinematics[3], fPrimAngles[3], fDeta[3], fDphi[3], fInvMass12[3], fInvMass23[3], fInvMass31[3], fhistDeltaPhi12[3], fhistDeltaPhi23[3], fhistDeltaPhi31[3], fhistDeltaPhiDeltaEta12[3], fhistDeltaPhiDeltaEta23[3], fhistDeltaPhiDeltaEta31[3], fSameEventP1[3], fSameEventK[3]);
       // Proton Proton AntiPrimary
       FillTripletDistribution( ParticleVector, 3, 0, 0, fSameEventTripletArray[4],PDGCodes, bins[1],fSameEventTripletMultArray[4], fSameEventTripletMultArray12[4] ,fSameEventTripletMultArray23[4], fSameEventTripletMultArray31[4], fSameEventTripletPtPrimaries[4], fSameEventTripletPtProtons[4], fSameEventTripletPtPrimaries2[4], fSameEventTripletPtProtons2[4], fSameEventTripletPtvsQ3Primaries[4], fSameEventTripletPtvsQ3Protons[4], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,4, *fConfig, fQ3VskDistributionsArrayq12[4],fQ3VskDistributionsArrayq23[4], fKinematics[4], fPrimAngles[4], fDeta[4], fDphi[4], fInvMass12[4], fInvMass23[4], fInvMass31[4], fhistDeltaPhi12[4], fhistDeltaPhi23[4], fhistDeltaPhi31[4], fhistDeltaPhiDeltaEta12[4], fhistDeltaPhiDeltaEta23[4], fhistDeltaPhiDeltaEta31[4], fSameEventP1[4], fSameEventK[4], fSameEventTripletmTArray12[4], fSameEventTripletmTArray23[4], fSameEventTripletmTArray31[4], fInvMassVSmT12[4], fInvMassVSmT23[4], fInvMassVSmT31[4]);
       // Antiproton Antiproton Primary
       FillTripletDistribution( ParticleVector, 2, 1, 1, fSameEventTripletArray[5],PDGCodes, bins[1],fSameEventTripletMultArray[5], fSameEventTripletMultArray12[5] ,fSameEventTripletMultArray23[5], fSameEventTripletMultArray31[5], fSameEventTripletPtPrimaries[5], fSameEventTripletPtProtons[5], fSameEventTripletPtPrimaries2[5], fSameEventTripletPtProtons2[5], fSameEventTripletPtvsQ3Primaries[5], fSameEventTripletPtvsQ3Protons[5], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,5, *fConfig, fQ3VskDistributionsArrayq12[5],fQ3VskDistributionsArrayq23[5], fKinematics[5], fPrimAngles[5], fDeta[5], fDphi[5], fInvMass12[5], fInvMass23[5], fInvMass31[5], fhistDeltaPhi12[5], fhistDeltaPhi23[5], fhistDeltaPhi31[5], fhistDeltaPhiDeltaEta12[5], fhistDeltaPhiDeltaEta23[5], fhistDeltaPhiDeltaEta31[5], fSameEventP1[5], fSameEventK[5], fSameEventTripletmTArray12[5], fSameEventTripletmTArray23[5], fSameEventTripletmTArray31[5], fInvMassVSmT12[5], fInvMassVSmT23[5], fInvMassVSmT31[5]);

    }//if(fDoOnlyThreeBody)
    else {
      //Two Body Analyses...........
      //Proton Proton
      FillPairDistribution( ParticleVector, 0, 0, fSameEventTripletArray_TwoBody[0],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[0], fSameEventDphiArray_TwoBody[0], fSameEventMultArray_TwoBody[0], fSameEventTripletPhiThetaArray_TwoBody,0, *fConfig);
//      FillPairDistribution( ParticleVector, 0, 1, fSameEventTripletArray_TwoBody[0],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[0], fSameEventDphiArray_TwoBody[0], fSameEventMultArray_TwoBody[0], fSameEventTripletPhiThetaArray_TwoBody,0, *fConfig);
      //AntiProton AntiProton
      FillPairDistribution( ParticleVector, 1, 1, fSameEventTripletArray_TwoBody[1],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[1], fSameEventDphiArray_TwoBody[1], fSameEventMultArray_TwoBody[1], fSameEventTripletPhiThetaArray_TwoBody,1, *fConfig);
      //Proton Primary
      FillPairDistribution( ParticleVector, 0, 2, fSameEventTripletArray_TwoBody[2],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[2], fSameEventDphiArray_TwoBody[2], fSameEventMultArray_TwoBody[2], fSameEventTripletPhiThetaArray_TwoBody,2, *fConfig);
      //Antiproton Antiprimary
      FillPairDistribution( ParticleVector, 1, 3, fSameEventTripletArray_TwoBody[3],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[3], fSameEventDphiArray_TwoBody[3], fSameEventMultArray_TwoBody[3], fSameEventTripletPhiThetaArray_TwoBody,3, *fConfig);
      //Proton Antiprimary
      FillPairDistribution( ParticleVector, 0, 3, fSameEventTripletArray_TwoBody[4],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[4], fSameEventDphiArray_TwoBody[4], fSameEventMultArray_TwoBody[4], fSameEventTripletPhiThetaArray_TwoBody,4, *fConfig);
      //Antiproton Primary
      FillPairDistribution( ParticleVector, 1, 2, fSameEventTripletArray_TwoBody[5],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[5], fSameEventDphiArray_TwoBody[5], fSameEventMultArray_TwoBody[5], fSameEventTripletPhiThetaArray_TwoBody,5, *fConfig);
     //Primary Antiprimary
      FillPairDistribution( ParticleVector, 2, 3, fSameEventTripletArray_TwoBody[6],PDGCodes, bins[1],fSameEventTripletMultArray_TwoBody[6], fSameEventDphiArray_TwoBody[6], fSameEventMultArray_TwoBody[6], fSameEventTripletPhiThetaArray_TwoBody,6, *fConfig);

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
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 0, 0, 2, fSameEventTripletArray[10], PDGCodes, bins[1],fSameEventTripletMultArray[10], fSameEventTripletPtPrimaries[10], fSameEventTripletPtProtons[10], fSameEventTripletPtPrimaries2[10], fSameEventTripletPtProtons2[10], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 10, *fConfig, fQ3VskDistributionsArrayq12[6],fQ3VskDistributionsArrayq23[6], fKinematics[6], fPrimAngles[6], fDeta[6], fDphi[6], fDetaSEvsME[0], fDphiSEvsME[0], fInvMass12[6], fInvMass23[6], fInvMass31[6], f2Same1MixedEventTripletmTArray12[0], f2Same1MixedEventTripletmTArray23[0], f2Same1MixedEventTripletmTArray31[0], fInvMassVSmT12[6], fInvMassVSmT23[6], fInvMassVSmT31[6]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 0, 2, 0, fSameEventTripletArray[11], PDGCodes, bins[1],fSameEventTripletMultArray[11], fSameEventTripletPtPrimaries[11], fSameEventTripletPtProtons[11], fSameEventTripletPtPrimaries2[11], fSameEventTripletPtProtons2[11], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 11, *fConfig, fQ3VskDistributionsArrayq12[7],fQ3VskDistributionsArrayq23[7], fKinematics[7], fPrimAngles[7], fDeta[7], fDphi[7], fDetaSEvsME[1], fDphiSEvsME[1], fInvMass12[7], fInvMass23[7], fInvMass31[7], f2Same1MixedEventTripletmTArray12[1], f2Same1MixedEventTripletmTArray23[1], f2Same1MixedEventTripletmTArray31[1], fInvMassVSmT12[7], fInvMassVSmT23[7], fInvMassVSmT31[7]);

         // Antiproton Antiproton Antiprimary
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 1, 1, 3, fSameEventTripletArray[12], PDGCodes, bins[1],fSameEventTripletMultArray[12], fSameEventTripletPtPrimaries[12], fSameEventTripletPtProtons[12], fSameEventTripletPtPrimaries2[12], fSameEventTripletPtProtons2[12], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 12, *fConfig, fQ3VskDistributionsArrayq12[8],fQ3VskDistributionsArrayq23[8], fKinematics[8], fPrimAngles[8], fDeta[8], fDphi[8], fDetaSEvsME[2], fDphiSEvsME[2], fInvMass12[8], fInvMass23[8], fInvMass31[8], f2Same1MixedEventTripletmTArray12[2], f2Same1MixedEventTripletmTArray23[2], f2Same1MixedEventTripletmTArray31[2], fInvMassVSmT12[8], fInvMassVSmT23[8], fInvMassVSmT31[8]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPPrim], 1, 3, 1, fSameEventTripletArray[13], PDGCodes, bins[1],fSameEventTripletMultArray[13], fSameEventTripletPtPrimaries[13], fSameEventTripletPtProtons[13], fSameEventTripletPtPrimaries2[13], fSameEventTripletPtProtons2[13], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 13, *fConfig, fQ3VskDistributionsArrayq12[9],fQ3VskDistributionsArrayq23[9], fKinematics[9], fPrimAngles[9], fDeta[9], fDphi[9], fDetaSEvsME[3], fDphiSEvsME[3], fInvMass12[9], fInvMass23[9], fInvMass31[9], f2Same1MixedEventTripletmTArray12[3], f2Same1MixedEventTripletmTArray23[3], f2Same1MixedEventTripletmTArray31[3], fInvMassVSmT12[9], fInvMassVSmT23[9], fInvMassVSmT31[9]);

         // Proton Proton Proton
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPP], 0, 0, 0, fSameEventTripletArray[14], PDGCodes, bins[1],fSameEventTripletMultArray[14], fSameEventTripletPtPrimaries[14], fSameEventTripletPtProtons[14], fSameEventTripletPtPrimaries2[14], fSameEventTripletPtProtons2[14], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 14, *fConfig, fQ3VskDistributionsArrayq12[10],fQ3VskDistributionsArrayq23[10], fKinematics[10], fPrimAngles[10], fDeta[10], fDphi[10], fDetaSEvsME[4], fDphiSEvsME[4], fInvMass12[10], fInvMass23[10], fInvMass31[10], f2Same1MixedEventTripletmTArray12[4], f2Same1MixedEventTripletmTArray23[4], f2Same1MixedEventTripletmTArray31[4], fInvMassVSmT12[10], fInvMassVSmT23[10], fInvMassVSmT31[10]);

         // Antiproton Antiproton Antiproton
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPP], 1, 1, 1, fSameEventTripletArray[15], PDGCodes, bins[1],fSameEventTripletMultArray[15], fSameEventTripletPtPrimaries[15], fSameEventTripletPtProtons[15], fSameEventTripletPtPrimaries2[15], fSameEventTripletPtProtons2[15], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 15, *fConfig, fQ3VskDistributionsArrayq12[11],fQ3VskDistributionsArrayq23[11], fKinematics[11], fPrimAngles[11], fDeta[11], fDphi[11], fDetaSEvsME[5], fDphiSEvsME[5], fInvMass12[11], fInvMass23[11], fInvMass31[11], f2Same1MixedEventTripletmTArray12[5], f2Same1MixedEventTripletmTArray23[5], f2Same1MixedEventTripletmTArray31[5], fInvMassVSmT12[11], fInvMassVSmT23[11], fInvMassVSmT31[11]);

         // Proton Proton AntiPrimary
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 0, 0, 3, fSameEventTripletArray[16], PDGCodes, bins[1],fSameEventTripletMultArray[16], fSameEventTripletPtPrimaries[16], fSameEventTripletPtProtons[16], fSameEventTripletPtPrimaries2[16], fSameEventTripletPtProtons2[16], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 16, *fConfig, fQ3VskDistributionsArrayq12[12],fQ3VskDistributionsArrayq23[12], fKinematics[12], fPrimAngles[12], fDeta[12], fDphi[12], fDetaSEvsME[6], fDphiSEvsME[6], fInvMass12[12], fInvMass23[12], fInvMass31[12], f2Same1MixedEventTripletmTArray12[6], f2Same1MixedEventTripletmTArray23[6], f2Same1MixedEventTripletmTArray31[6], fInvMassVSmT12[12], fInvMassVSmT23[12], fInvMassVSmT31[12]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 0, 3, 0, fSameEventTripletArray[17], PDGCodes, bins[1],fSameEventTripletMultArray[17], fSameEventTripletPtPrimaries[17], fSameEventTripletPtProtons[17], fSameEventTripletPtPrimaries2[17], fSameEventTripletPtProtons2[17], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 17, *fConfig, fQ3VskDistributionsArrayq12[13],fQ3VskDistributionsArrayq23[13], fKinematics[13], fPrimAngles[13], fDeta[13], fDphi[13], fDetaSEvsME[7], fDphiSEvsME[7], fInvMass12[13], fInvMass23[13], fInvMass31[13], f2Same1MixedEventTripletmTArray12[7], f2Same1MixedEventTripletmTArray23[7], f2Same1MixedEventTripletmTArray31[7], fInvMassVSmT12[13], fInvMassVSmT23[13], fInvMassVSmT31[13]);

         // Antiproton Antiproton Primary
         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 1, 1, 2, fSameEventTripletArray[18], PDGCodes, bins[1],fSameEventTripletMultArray[18], fSameEventTripletPtPrimaries[18], fSameEventTripletPtProtons[18], fSameEventTripletPtPrimaries2[18], fSameEventTripletPtProtons2[18], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 18, *fConfig, fQ3VskDistributionsArrayq12[14],fQ3VskDistributionsArrayq23[14], fKinematics[14], fPrimAngles[14], fDeta[14], fDphi[14], fDetaSEvsME[8], fDphiSEvsME[8], fInvMass12[14], fInvMass23[14], fInvMass31[14], f2Same1MixedEventTripletmTArray12[8], f2Same1MixedEventTripletmTArray23[8], f2Same1MixedEventTripletmTArray31[8], fInvMassVSmT12[14], fInvMassVSmT23[14], fInvMassVSmT31[14]);

         FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerIdPPAPrim], 1, 2, 1, fSameEventTripletArray[19], PDGCodes, bins[1],fSameEventTripletMultArray[19], fSameEventTripletPtPrimaries[19], fSameEventTripletPtProtons[19], fSameEventTripletPtPrimaries2[19], fSameEventTripletPtProtons2[19], fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, 19, *fConfig, fQ3VskDistributionsArrayq12[15],fQ3VskDistributionsArrayq23[15], fKinematics[15], fPrimAngles[15], fDeta[15], fDphi[15], fDetaSEvsME[9], fDphiSEvsME[9], fInvMass12[15], fInvMass23[15], fInvMass31[15], f2Same1MixedEventTripletmTArray12[9], f2Same1MixedEventTripletmTArray23[9], f2Same1MixedEventTripletmTArray31[9], fInvMassVSmT12[15], fInvMassVSmT23[15], fInvMassVSmT31[15]);

      }//if(fDoOnlyThreeBody)

      //c.1.1) Normal mixing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(fDoOnlyThreeBody){

           // 2 (Anti-)Protons, 1 (Anti-)Primary
           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPPrim], 2, 0, 0, fMixedEventTripletArray[0], PDGCodes, bins[1],fMixedEventTripletMultArray[0], fMixedEventTripletMultArray12[0], fMixedEventTripletMultArray23[0], fMixedEventTripletMultArray31[0], fMixedEventTripletPtPrimaries[0], fMixedEventTripletPtProtons[0], fMixedEventTripletPtPrimaries2[0], fMixedEventTripletPtProtons2[0], fMixedEventTripletPtvsQ3Primaries[0], fMixedEventTripletPtvsQ3Protons[0], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,0, *fConfig, fQ3VskDistributionsArrayq12Mixed[0], fQ3VskDistributionsArrayq23Mixed[0], fKinematicsME[0], fPrimAnglesME[0], fDetaME[0], fDphiME[0], fhistDeltaPhi12ME[0], fhistDeltaPhi23ME[0], fhistDeltaPhi31ME[0], fhistDeltaPhiDeltaEta12ME[0], fhistDeltaPhiDeltaEta23ME[0], fhistDeltaPhiDeltaEta31ME[0], fMixedEventP1[0], fMixedEventK[0], fMixedEventTripletmTArray12[0], fMixedEventTripletmTArray23[0], fMixedEventTripletmTArray31[0], fInvMassVSmT12MixedEvent[0], fInvMassVSmT23MixedEvent[0], fInvMassVSmT31MixedEvent[0]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPPrim], 3, 1, 1, fMixedEventTripletArray[1], PDGCodes, bins[1],fMixedEventTripletMultArray[1], fMixedEventTripletMultArray12[1], fMixedEventTripletMultArray23[1], fMixedEventTripletMultArray31[1], fMixedEventTripletPtPrimaries[1], fMixedEventTripletPtProtons[1], fMixedEventTripletPtPrimaries2[1], fMixedEventTripletPtProtons2[1], fMixedEventTripletPtvsQ3Primaries[1], fMixedEventTripletPtvsQ3Protons[1], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,1, *fConfig, fQ3VskDistributionsArrayq12Mixed[1], fQ3VskDistributionsArrayq23Mixed[1], fKinematicsME[1], fPrimAnglesME[1], fDetaME[1], fDphiME[1], fhistDeltaPhi12ME[1], fhistDeltaPhi23ME[1], fhistDeltaPhi31ME[1], fhistDeltaPhiDeltaEta12ME[1], fhistDeltaPhiDeltaEta23ME[1], fhistDeltaPhiDeltaEta31ME[1], fMixedEventP1[1], fMixedEventK[1], fMixedEventTripletmTArray12[1], fMixedEventTripletmTArray23[1], fMixedEventTripletmTArray31[1], fInvMassVSmT12MixedEvent[1], fInvMassVSmT23MixedEvent[1], fInvMassVSmT31MixedEvent[1]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPP], 0, 0, 0, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletMultArray12[2], fMixedEventTripletMultArray23[2], fMixedEventTripletMultArray31[2], fMixedEventTripletPtPrimaries[2], fMixedEventTripletPtProtons[2], fMixedEventTripletPtPrimaries2[2], fMixedEventTripletPtProtons2[2], fMixedEventTripletPtvsQ3Primaries[2], fMixedEventTripletPtvsQ3Protons[2], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,2, *fConfig, fQ3VskDistributionsArrayq12Mixed[2], fQ3VskDistributionsArrayq23Mixed[2], fKinematicsME[2], fPrimAnglesME[2], fDetaME[2], fDphiME[2], fhistDeltaPhi12ME[2], fhistDeltaPhi23ME[2], fhistDeltaPhi31ME[2], fhistDeltaPhiDeltaEta12ME[2], fhistDeltaPhiDeltaEta23ME[2], fhistDeltaPhiDeltaEta31ME[2], fMixedEventP1[2], fMixedEventK[2], fMixedEventTripletmTArray12[2], fMixedEventTripletmTArray23[2], fMixedEventTripletmTArray31[2], fInvMassVSmT12MixedEvent[2], fInvMassVSmT23MixedEvent[2], fInvMassVSmT31MixedEvent[2]);
//           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPP], 1, 0, 0, fMixedEventTripletArray[2], PDGCodes, bins[1],fMixedEventTripletMultArray[2], fMixedEventTripletMultArray12[2], fMixedEventTripletMultArray23[2], fMixedEventTripletMultArray31[2], fMixedEventTripletPtPrimaries[2], fMixedEventTripletPtProtons[2], fMixedEventTripletPtPrimaries2[2], fMixedEventTripletPtProtons2[2], fMixedEventTripletPtvsQ3Primaries[2], fMixedEventTripletPtvsQ3Protons[2], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,2, *fConfig, fQ3VskDistributionsArrayq12Mixed[2], fQ3VskDistributionsArrayq23Mixed[2], fKinematicsME[2], fPrimAnglesME[2], fDetaME[2], fDphiME[2], fhistDeltaPhi12ME[2], fhistDeltaPhi23ME[2], fhistDeltaPhi31ME[2], fhistDeltaPhiDeltaEta12ME[2], fhistDeltaPhiDeltaEta23ME[2], fhistDeltaPhiDeltaEta31ME[2], fMixedEventP1[2], fMixedEventK[2]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPP], 1, 1, 1, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletMultArray12[3], fMixedEventTripletMultArray23[3], fMixedEventTripletMultArray31[3], fMixedEventTripletPtPrimaries[3], fMixedEventTripletPtProtons[3], fMixedEventTripletPtPrimaries2[3], fMixedEventTripletPtProtons2[3], fMixedEventTripletPtvsQ3Primaries[3], fMixedEventTripletPtvsQ3Protons[3], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,3, *fConfig, fQ3VskDistributionsArrayq12Mixed[3], fQ3VskDistributionsArrayq23Mixed[3], fKinematicsME[3], fPrimAnglesME[3], fDetaME[3], fDphiME[3], fhistDeltaPhi12ME[3], fhistDeltaPhi23ME[3], fhistDeltaPhi31ME[3], fhistDeltaPhiDeltaEta12ME[3], fhistDeltaPhiDeltaEta23ME[3], fhistDeltaPhiDeltaEta31ME[3], fMixedEventP1[3], fMixedEventK[3], fMixedEventTripletmTArray12[3], fMixedEventTripletmTArray23[3], fMixedEventTripletmTArray31[3], fInvMassVSmT12MixedEvent[3], fInvMassVSmT23MixedEvent[3], fInvMassVSmT31MixedEvent[3]);
//           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPP], 0, 1, 1, fMixedEventTripletArray[3], PDGCodes, bins[1],fMixedEventTripletMultArray[3], fMixedEventTripletMultArray12[3], fMixedEventTripletMultArray23[3], fMixedEventTripletMultArray31[3], fMixedEventTripletPtPrimaries[3], fMixedEventTripletPtProtons[3], fMixedEventTripletPtPrimaries2[3], fMixedEventTripletPtProtons2[3], fMixedEventTripletPtvsQ3Primaries[3], fMixedEventTripletPtvsQ3Protons[3], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,3, *fConfig, fQ3VskDistributionsArrayq12Mixed[3], fQ3VskDistributionsArrayq23Mixed[3], fKinematicsME[3], fPrimAnglesME[3], fDetaME[3], fDphiME[3], fhistDeltaPhi12ME[3], fhistDeltaPhi23ME[3], fhistDeltaPhi31ME[3], fhistDeltaPhiDeltaEta12ME[3], fhistDeltaPhiDeltaEta23ME[3], fhistDeltaPhiDeltaEta31ME[3], fMixedEventP1[3], fMixedEventK[3]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPAPrim], 3, 0, 0, fMixedEventTripletArray[4], PDGCodes, bins[1],fMixedEventTripletMultArray[4], fMixedEventTripletMultArray12[4], fMixedEventTripletMultArray23[4], fMixedEventTripletMultArray31[4], fMixedEventTripletPtPrimaries[4], fMixedEventTripletPtProtons[4], fMixedEventTripletPtPrimaries2[4], fMixedEventTripletPtProtons2[4], fMixedEventTripletPtvsQ3Primaries[4], fMixedEventTripletPtvsQ3Protons[4], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,4, *fConfig, fQ3VskDistributionsArrayq12Mixed[4], fQ3VskDistributionsArrayq23Mixed[4], fKinematicsME[4], fPrimAnglesME[4], fDetaME[4], fDphiME[4], fhistDeltaPhi12ME[4], fhistDeltaPhi23ME[4], fhistDeltaPhi31ME[4], fhistDeltaPhiDeltaEta12ME[4], fhistDeltaPhiDeltaEta23ME[4], fhistDeltaPhiDeltaEta31ME[4], fMixedEventP1[4], fMixedEventK[4], fMixedEventTripletmTArray12[4], fMixedEventTripletmTArray23[4], fMixedEventTripletmTArray31[4], fInvMassVSmT12MixedEvent[4], fInvMassVSmT23MixedEvent[4], fInvMassVSmT31MixedEvent[4]);

           FillTripletDistributionME(ParticleVector, *VectItMult[ContainerIdPPAPrim], 2, 1, 1, fMixedEventTripletArray[5], PDGCodes, bins[1],fMixedEventTripletMultArray[5], fMixedEventTripletMultArray12[5], fMixedEventTripletMultArray23[5], fMixedEventTripletMultArray31[5], fMixedEventTripletPtPrimaries[5], fMixedEventTripletPtProtons[5], fMixedEventTripletPtPrimaries2[5], fMixedEventTripletPtProtons2[5], fMixedEventTripletPtvsQ3Primaries[5], fMixedEventTripletPtvsQ3Protons[5], fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair,5, *fConfig, fQ3VskDistributionsArrayq12Mixed[5], fQ3VskDistributionsArrayq23Mixed[5], fKinematicsME[5], fPrimAnglesME[5], fDetaME[5], fDphiME[5], fhistDeltaPhi12ME[5], fhistDeltaPhi23ME[5], fhistDeltaPhi31ME[5], fhistDeltaPhiDeltaEta12ME[5], fhistDeltaPhiDeltaEta23ME[5], fhistDeltaPhiDeltaEta31ME[5], fMixedEventP1[5], fMixedEventK[5], fMixedEventTripletmTArray12[5], fMixedEventTripletmTArray23[5], fMixedEventTripletmTArray31[5], fInvMassVSmT12MixedEvent[5], fInvMassVSmT23MixedEvent[5], fInvMassVSmT31MixedEvent[5]);
  
      } else {
        //Two Body Analyses...........

        //Proton Proton
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPP], 0, 0, fMixedEventTripletArray_TwoBody[0],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[0], fMixedEventDphiArray_TwoBody[0], fMixedEventMultArray_TwoBody[0], fMixedEventTripletPhiThetaArray_TwoBody,0, *fConfig);
//        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPP], 1, 0, fMixedEventTripletArray_TwoBody[0],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[0], fMixedEventDphiArray_TwoBody[0], fMixedEventMultArray_TwoBody[0], fMixedEventTripletPhiThetaArray_TwoBody,0, *fConfig);
        //AntiProton AntiProton
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPP], 1, 1, fMixedEventTripletArray_TwoBody[1],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[1], fMixedEventDphiArray_TwoBody[1], fMixedEventMultArray_TwoBody[1], fMixedEventTripletPhiThetaArray_TwoBody,1, *fConfig);
        //Proton Primary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPPrim], 0, 2, fMixedEventTripletArray_TwoBody[2],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[2], fMixedEventDphiArray_TwoBody[2], fMixedEventMultArray_TwoBody[2], fMixedEventTripletPhiThetaArray_TwoBody,2, *fConfig);
        //Antiproton Antiprimary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPPrim], 1, 3, fMixedEventTripletArray_TwoBody[3],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[3], fMixedEventDphiArray_TwoBody[3], fMixedEventMultArray_TwoBody[3], fMixedEventTripletPhiThetaArray_TwoBody,3, *fConfig);
        //Proton Antiprimary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPAPrim], 0, 3, fMixedEventTripletArray_TwoBody[4],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[4], fMixedEventDphiArray_TwoBody[4], fMixedEventMultArray_TwoBody[4], fMixedEventTripletPhiThetaArray_TwoBody,4, *fConfig);
        //Antiproton Primary
        FillPairDistributionME( ParticleVector, *VectItMult[ContainerIdPAPrim], 1, 2, fMixedEventTripletArray_TwoBody[5],PDGCodes, bins[1],fMixedEventTripletMultArray_TwoBody[5], fMixedEventDphiArray_TwoBody[5], fMixedEventMultArray_TwoBody[5], fMixedEventTripletPhiThetaArray_TwoBody,5, *fConfig);
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

float AliAnalysisTaskThreeBodyProtonPrimary::BoostOneParticle(
    TLorentzVector &PartOne, TLorentzVector &PartTwo, TLorentzVector &PartThree) {
    // The code boost PartThree in the rest frame of PartOne and PartTwo
  TLorentzVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartThreeCMS = PartThree;

  PartThreeCMS.Boost(-betax, -betay, -betaz);

  return PartThreeCMS.P();
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F* hPtvsQ3Primaries, TH2F* hPtvsQ3Protons, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31 , TH2F* histDeltaPhi12, TH2F* histDeltaPhi23, TH2F* histDeltaPhi31, TH2F* histDeltaPhiDeltaEta12, TH2F* histDeltaPhiDeltaEta23, TH2F* histDeltaPhiDeltaEta31, TH2F* histP, TH2F* histK, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31){
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

        float RelativeMomentum12 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec,part2_LorVec);
        float RelativeMomentum23 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part2_LorVec,part3_LorVec);
        float RelativeMomentum31 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part3_LorVec,part1_LorVec);

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

 //       cout<<abs(*itPDGPar1)<<"  "<<abs(*itPDGPar2)<<" "<<abs(*itPDGPar3)<<endl;
 //       cout<<part1_LorVec.Pt()<<"   "<<part2_LorVec.Pt()<<"   "<<part3_LorVec.Pt()<<endl;

        if(fPlotP1){
            float p1boosted = BoostOneParticle(part2_LorVec,part3_LorVec,part1_LorVec);
            histP->Fill(Q3,p1boosted);
            histK->Fill(Q3,RelativeMomentum23);
        }

        double deltaPhi12 = (iPart1->GetPhi().at(0) - iPart2->GetPhi().at(0));
        double deltaPhi23 = (iPart2->GetPhi().at(0) - iPart3->GetPhi().at(0));
        double deltaPhi31 = (iPart3->GetPhi().at(0) - iPart1->GetPhi().at(0));
        double deltaEta12 = (iPart1->GetEta().at(0) - iPart2->GetEta().at(0));
        double deltaEta23 = (iPart2->GetEta().at(0) - iPart3->GetEta().at(0));
        double deltaEta31 = (iPart3->GetEta().at(0) - iPart1->GetEta().at(0));


        if(fQ3MinValue <= Q3 && Q3<fQ3cutValue){
          hist2d12->Fill(RelativeMomentum12,mult+1);
          hist2d23->Fill(RelativeMomentum23,mult+1);
          hist2d31->Fill(RelativeMomentum31,mult+1);
          histDeltaPhi12->Fill(RelativeMomentum12,deltaPhi12);
          histDeltaPhi23->Fill(RelativeMomentum23,deltaPhi23);
          histDeltaPhi31->Fill(RelativeMomentum31,deltaPhi31);
          histDeltaPhiDeltaEta12->Fill(deltaPhi12,deltaEta12);
          histDeltaPhiDeltaEta23->Fill(deltaPhi23,deltaEta23);
          histDeltaPhiDeltaEta31->Fill(deltaPhi31,deltaEta31);
        }

        if(fRunPlotPt){

          if (Q3<1.2) // MODIFIED BY ME
          {
            if(abs(*itPDGPar2)==2212) hPtProtons -> Fill(part2_LorVec.Pt());
          }
          if (Q3<1.2) // Second Proton
          {
            if(abs(*itPDGPar3)==2212) hPtProtons2 -> Fill(part3_LorVec.Pt());
          }


/*
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
*/

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
           TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
           InvMass12->Fill(Sum12.Mag(),Q3);
           InvMass23->Fill(Sum23.Mag(),Q3);
           InvMass31->Fill(Sum31.Mag(),Q3);

        }//if(fRunPlotInvMass)

        if(fRunmTPlots){
          float mT12 = GetmT(part1_LorVec, massparticle1, part2_LorVec, massparticle2);
          float mT23 = GetmT(part2_LorVec, massparticle2, part3_LorVec, massparticle3);
          float mT31 = GetmT(part3_LorVec, massparticle3, part1_LorVec, massparticle1);

          histmTQ312->Fill(Q3,mT12);
          histmTQ323->Fill(Q3,mT23);
          histmTQ331->Fill(Q3,mT31);
        }

        if(fRunPlotInvMassVSmT){
          TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
          TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
          TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
          float mT12 = GetmT(part1_LorVec, massparticle1, part2_LorVec, massparticle2);
          float mT23 = GetmT(part2_LorVec, massparticle2, part3_LorVec, massparticle3);
          float mT31 = GetmT(part3_LorVec, massparticle3, part1_LorVec, massparticle1);
          InvMassVsmT12->Fill(Sum12.Mag(),mT12);
          InvMassVsmT23->Fill(Sum23.Mag(),mT23);
          InvMassVsmT31->Fill(Sum31.Mag(),mT31);

        }

      }
    }
  }
}

//================================================================================================================================================== //GANESHA CHECK!!
void AliAnalysisTaskThreeBodyProtonPrimary::FillPairDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* InvMassSame){
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

        double deltaPhi = (iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0));

        histDeltaPhi->Fill(RelativeMomentum,deltaPhi);
        if(abs(deltaPhi)>0.5&&abs(deltaPhi)<2.5){
          histLowDeltaPhi->Fill(RelativeMomentum,mult+1);
        }



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

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F* hPtvsQ3Primaries, TH2F* hPtvsQ3Protons, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* histDeltaPhi12, TH2F* histDeltaPhi23, TH2F* histDeltaPhi31, TH2F* histDeltaPhiDeltaEta12, TH2F* histDeltaPhiDeltaEta23, TH2F* histDeltaPhiDeltaEta31, TH2F* histP, TH2F* histK, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31){//, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed){
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

            float RelativeMomentum12 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec,part2_LorVec);
            float RelativeMomentum23 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part2_LorVec,part3_LorVec);
            float RelativeMomentum31 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part3_LorVec,part1_LorVec);

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

            if(fPlotP1){
              float p1boosted = BoostOneParticle(part2_LorVec,part3_LorVec,part1_LorVec);
              histP->Fill(Q3,p1boosted);
              histK->Fill(Q3,RelativeMomentum23);
            }

            double deltaPhi12 = (iPart1->GetPhi().at(0) - iPart2->GetPhi().at(0));
            double deltaPhi23 = (iPart2->GetPhi().at(0) - iPart3->GetPhi().at(0));
            double deltaPhi31 = (iPart3->GetPhi().at(0) - iPart1->GetPhi().at(0));
            double deltaEta12 = (iPart1->GetEta().at(0) - iPart2->GetEta().at(0));
            double deltaEta23 = (iPart2->GetEta().at(0) - iPart3->GetEta().at(0));
            double deltaEta31 = (iPart3->GetEta().at(0) - iPart1->GetEta().at(0));


    
            if(fQ3MinValue <= Q3 && Q3<fQ3cutValue){
              hist2d12->Fill(RelativeMomentum12,mult+1);
              hist2d23->Fill(RelativeMomentum23,mult+1);
              hist2d31->Fill(RelativeMomentum31,mult+1);
              histDeltaPhi12->Fill(RelativeMomentum12,deltaPhi12);
              histDeltaPhi23->Fill(RelativeMomentum23,deltaPhi23);
              histDeltaPhi31->Fill(RelativeMomentum31,deltaPhi31);
              histDeltaPhiDeltaEta12->Fill(deltaPhi12,deltaEta12);
              histDeltaPhiDeltaEta23->Fill(deltaPhi23,deltaEta23);
              histDeltaPhiDeltaEta31->Fill(deltaPhi31,deltaEta31);
            }

           if(fRunmTPlots){
             float mT12 = GetmT(part1_LorVec, massParticleSE, part2_LorVec, massParticleME1);
             float mT23 = GetmT(part2_LorVec, massParticleME1, part3_LorVec, massParticleME2);
             float mT31 = GetmT(part3_LorVec, massParticleME2, part1_LorVec, massParticleSE);

             histmTQ312->Fill(Q3,mT12);
             histmTQ323->Fill(Q3,mT23);
             histmTQ331->Fill(Q3,mT31);
           }

            if(fRunPlotPt){


            if (Q3<1.2) // MODIFIED BY ME
            {
              if(abs(*itPDGParME1)==2212) hPtProtons -> Fill(part2_LorVec.Pt());
            }
            if (Q3<1.2) // Second Proton
            {
              if(abs(*itPDGParME2)==2212) hPtProtons2 -> Fill(part3_LorVec.Pt());
            }
/*
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
*/
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


        if(fRunPlotInvMassVSmT){
          TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
          TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
          TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
          float mT12 = GetmT(part1_LorVec, massParticleSE, part2_LorVec, massParticleME1);
          float mT23 = GetmT(part2_LorVec, massParticleME1, part3_LorVec, massParticleME2);
          float mT31 = GetmT(part3_LorVec, massParticleME2, part1_LorVec, massParticleSE);
          InvMassVsmT12->Fill(Sum12.Mag(),mT12);
          InvMassVsmT23->Fill(Sum23.Mag(),mT23);
          InvMassVsmT31->Fill(Sum31.Mag(),mT31);

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

void AliAnalysisTaskThreeBodyProtonPrimary::FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){ //GANESHA CHECK

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


        double deltaPhi = (iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0));

        histDeltaPhi->Fill(RelativeK,deltaPhi);
        if(abs(deltaPhi)>0.5&&abs(deltaPhi)<2.5){
          histLowDeltaPhi->Fill(RelativeK,mult+1);
        }



      }
    }
  }
}

//==================================================================================================================================================

void AliAnalysisTaskThreeBodyProtonPrimary::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH1F* hPtPrimaries, TH1F* hPtProtons, TH1F* hPtPrimaries2, TH1F* hPtProtons2, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* Q3VskDistribution12, TH2F* Q3VskDistribution23, TH2F* hKinematics, TH2F* hPrimAngles, TH2F* hDeta, TH2F* hDphi, TH2F* hDetaSEvsME, TH2F* hDphiSEvsME, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31){
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

            bool fillFirst = false;
            bool fillSecond = false;

            if (Q3<1.2) // MODIFIED BY ME
            {
              if(abs(*itPDGParSE1)==2212){
                hPtProtons -> Fill(part1_LorVec.Pt());
                fillFirst = true;
              }
              if(abs(*itPDGParSE2)==2212&&!fillFirst){
                hPtProtons -> Fill(part2_LorVec.Pt());
                fillFirst = true;
              }else{
                hPtProtons2 -> Fill(part2_LorVec.Pt());
                fillSecond = true;
              }
              if(abs(*itPDGParME)==2212&&!fillSecond) hPtProtons2 -> Fill(part3_LorVec.Pt());
            }
/*
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
*/
          }
          if(fRunPlotQ3Vsq){
            float qSame12= sqrt(-q12*q12);
            float qSame23= sqrt(-q23*q23);
            Q3VskDistribution12->Fill(Q3, qSame12);
            Q3VskDistribution23->Fill(Q3, qSame23);
          }

          if(fRunmTPlots){
            float mT12 = GetmT(part1_LorVec, massParticleSE1, part2_LorVec, massParticleSE2);
            float mT23 = GetmT(part2_LorVec, massParticleSE2, part3_LorVec, massParticleME);
            float mT31 = GetmT(part3_LorVec, massParticleME, part1_LorVec, massParticleSE1);

            histmTQ312->Fill(Q3,mT12);
            histmTQ323->Fill(Q3,mT23);
            histmTQ331->Fill(Q3,mT31);
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
          double deltaEtaSE = 0, deltaEtaME = 0, deltaPhiSE = 0, deltaPhiME = 0;


          if(abs(*itPDGParSE1)==321){
            deltaEtaSE = iPart1->GetEta().at(0) - iPart2->GetEta().at(0);
            deltaPhiSE = iPart1->GetPhi().at(0) - iPart2->GetPhi().at(0);
            deltaEtaME = iPart1->GetEta().at(0) - iPart3->GetEta().at(0);
            deltaPhiME = iPart1->GetPhi().at(0) - iPart3->GetPhi().at(0);
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
            deltaEtaSE = iPart2->GetEta().at(0) - iPart1->GetEta().at(0);
            deltaPhiSE = iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0);
            deltaEtaME = iPart2->GetEta().at(0) - iPart3->GetEta().at(0);
            deltaPhiME = iPart2->GetPhi().at(0) - iPart3->GetPhi().at(0);
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

          if(abs(*itPDGParSE1)==2212&&abs(*itPDGParSE2)==2212){
            deltaEtaSE = iPart1->GetEta().at(0) - iPart2->GetEta().at(0);
            deltaPhiSE = iPart1->GetPhi().at(0) - iPart2->GetPhi().at(0);
            deltaEtaME = iPart1->GetEta().at(0) - iPart3->GetEta().at(0);
            deltaPhiME = iPart1->GetPhi().at(0) - iPart3->GetPhi().at(0);
          }

          if(Q3<1.2){
            hDetaSEvsME->Fill(deltaEtaSE,deltaEtaME);
            hDphiSEvsME->Fill(deltaPhiSE,deltaPhiME);
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
           TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
           TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
           InvMass12->Fill(Sum12.Mag(),Q3);
           InvMass23->Fill(Sum23.Mag(),Q3);
           InvMass31->Fill(Sum31.Mag(),Q3);
        }//if(fRunPlotInvMass)

        if(fRunPlotInvMassVSmT){
          TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
          TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
          TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
          float mT12 = GetmT(part1_LorVec, massParticleSE1, part2_LorVec, massParticleSE2);
          float mT23 = GetmT(part2_LorVec, massParticleSE2, part3_LorVec, massParticleME);
          float mT31 = GetmT(part3_LorVec, massParticleME, part1_LorVec, massParticleSE1);
          InvMassVsmT12->Fill(Sum12.Mag(),mT12);
          InvMassVsmT23->Fill(Sum23.Mag(),mT23);
          InvMassVsmT31->Fill(Sum31.Mag(),mT31);

        }



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

float AliAnalysisTaskThreeBodyProtonPrimary::GetmT(TLorentzVector &Part1, float mass1, TLorentzVector &Part2, float mass2){
  float mT = 0.; 

  TLorentzVector Sum;
  Sum = Part1 + Part2;
  float kT = 0.5 * Sum.Pt();

  float averageMass = 0.5 * (mass1 + mass2);

  mT = TMath::Sqrt(pow(kT, 2.) + pow(averageMass, 2.));
  return mT;
}
