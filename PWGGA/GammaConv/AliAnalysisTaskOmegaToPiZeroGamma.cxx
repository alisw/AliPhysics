/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Daniel Mühlheim, Friederike Bock                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskOmegaToPiZeroGamma.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskOmegaToPiZeroGamma)

//________________________________________________________________________
AliAnalysisTaskOmegaToPiZeroGamma::AliAnalysisTaskOmegaToPiZeroGamma(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fBGHandler(NULL),
  fBGClusHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fPhotonDCAList(NULL),
  fGammaERM02(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fClusterOutputList(NULL),
  fOutputContainer(NULL),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fCutArray(NULL),
  fConversionCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fNeutralPionCutArray(NULL),
  fMesonCutArray(NULL),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaR(NULL),
  fHistoConvGammaEta(NULL),
  fTreeConvGammaPtDcazCat(NULL),
  fPtGamma(0),
  fDCAzPhoton(0),
  fRConvPhoton(0),
  fEtaPhoton(0),
  fCharCatPhoton(0),
  fCharPhotonMCInfo(0),
  fHistoPhotonPairInvMassPt(NULL),
  fHistoPhotonPairMatchedInvMassPt(NULL),
  fHistoPhotonPairYPt(NULL),
  fHistoPhotonPairAlphaPt(NULL),
  fHistoPhotonPairOpenAnglePt(NULL),
  fHistoPhotonPairEtaPhi(NULL),
  fHistoMotherConvPhotonEtaPhi(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherYPt(NULL),
  fHistoMotherAlphaPt(NULL),
  fHistoMotherEtaPhi(NULL),
  fHistoMotherPi0AnglePt(NULL),
  fHistoMotherGammaAnglePt(NULL),
  fHistoPi0GammaAnglePt(NULL),
  fHistoGammaFromMotherPt(NULL),
  fHistoDiff1Diff2SameBackInvMassPt(NULL),
  fHistoDiffSameSameBackInvMassPt(NULL),
  fHistoSameDiffSameBackInvMassPt(NULL),
  fHistoSameSameDiffBackInvMassPt(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllGammaEMCALAccPt(NULL),
  fHistoMCConvGammaPt(NULL),
  fHistoMCConvGammaR(NULL),
  fHistoMCConvGammaEta(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCPi0PtY(NULL),
  fHistoMCPi0PtAlpha(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCGammaFromAllOmegaPt(NULL),
  fHistoMCGammaFromOmegaInAccPt(NULL),
  fHistoMCOmegaInAccInvMassPt(NULL),
  fHistoMCOmegaInvMassPt(NULL),
  fHistoMCAllOmegaYPt(NULL),
  fHistoMCOmegaInAccYPt(NULL),
  fHistoMCAllOmegaInvMassPt(NULL),
  fHistoMCAllOmegaAlphaPt(NULL),
  fHistoMCOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromAllOmegaAlphaPt(NULL),
  fHistoMCPi0FromOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromAllOmegaYPt(NULL),
  fHistoMCPi0FromOmegaInAccYPt(NULL),
  fHistoMCPi0FromAllOmegaEtaPhi(NULL),
  fHistoMCPi0FromOmegaInAccEtaPhi(NULL),
  fHistoMCPi0FromAllOmegaInvMassPt(NULL),
  fHistoMCPi0FromOmegaInAccInvMassPt(NULL),
  fHistoMCAllOmegaEtaPhi(NULL),
  fHistoMCOmegaInAccEtaPhi(NULL),
  fHistoMCAllOmegaPiZeroAnglePt(NULL),
  fHistoMCAllPiZeroGammaAnglePt(NULL),
  fHistoMCAllOmegaGammaAnglePt(NULL),
  fHistoMCInAccOmegaPiZeroAnglePt(NULL),
  fHistoMCInAccPiZeroGammaAnglePt(NULL),
  fHistoMCInAccOmegaGammaAnglePt(NULL),
  fHistoMCAllOmegaPtPi0Pt(NULL),
  fHistoMCInAccOmegaPtPi0Pt(NULL),
  fHistoMCAllOmegaPtGammaPt(NULL),
  fHistoMCInAccOmegaPtGammaPt(NULL),
  fHistoTrueOmegaInvMassPt(NULL),
  fHistoTruePi0FromOmegaInvMassPt(NULL),
  fHistoTruePi0FromOmegaAlphaPt(NULL),
  fHistoTruePi0FromOmegaYPt(NULL),
  fHistoTruePi0FromOmegaEtaPhi(NULL),
  fHistoTruePi0FromOmegaOpenAnglePt(NULL),
  fHistoTrueOmegaPi0AnglePt(NULL),
  fHistoTrueOmegaGammaAnglePt(NULL),
  fHistoTruePi0GammaAnglePt(NULL),
  fHistoTrueOmegaYPt(NULL),
  fHistoTrueOmegaAlphaPt(NULL),
  fHistoTrueOmegaEtaPhi(NULL),
  fHistoTrueOmegaPtPi0Pt(NULL),
  fHistoTrueOmegaPtGammaPt(NULL),
  fHistoTrueGammaFromOmegaPt(NULL),
  fVectorRecTruePi0s(0),
  fVectorDoubleCountTruePi0s(0),
  fHistoMultipleCountTruePi0(NULL),
  fMapMultipleCountTruePi0s(),
  fHistoNEvents(NULL),
  fHistoNEventsMinGamma(NULL),
  fHistoMCOmegaDecayChannels(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNClusterCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fEventPlaneAngle(-100),
  fRandom(0),
  fNGammaCandidates(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMCStackPos(NULL),
  fMCStackNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fnCuts(0),
  fiCut(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fEnableSortForClusMC(kFALSE),
  fReconMethod(0)
{
  
}

//________________________________________________________________________
AliAnalysisTaskOmegaToPiZeroGamma::AliAnalysisTaskOmegaToPiZeroGamma(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fBGHandler(NULL),
  fBGClusHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fPhotonDCAList(NULL),
  fGammaERM02(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fClusterOutputList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fCutArray(NULL),
  fConversionCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fNeutralPionCutArray(NULL),
  fMesonCutArray(NULL),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaR(NULL),
  fHistoConvGammaEta(NULL),
  fTreeConvGammaPtDcazCat(NULL),
  fPtGamma(0),
  fDCAzPhoton(0),
  fRConvPhoton(0),
  fEtaPhoton(0),
  fCharCatPhoton(0),
  fCharPhotonMCInfo(0),
  fHistoPhotonPairInvMassPt(NULL),
  fHistoPhotonPairMatchedInvMassPt(NULL),
  fHistoPhotonPairYPt(NULL),
  fHistoPhotonPairAlphaPt(NULL),
  fHistoPhotonPairOpenAnglePt(NULL),
  fHistoPhotonPairEtaPhi(NULL),
  fHistoMotherConvPhotonEtaPhi(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherYPt(NULL),
  fHistoMotherAlphaPt(NULL),
  fHistoMotherEtaPhi(NULL),
  fHistoMotherPi0AnglePt(NULL),
  fHistoMotherGammaAnglePt(NULL),
  fHistoPi0GammaAnglePt(NULL),
  fHistoGammaFromMotherPt(NULL),
  fHistoDiff1Diff2SameBackInvMassPt(NULL),
  fHistoDiffSameSameBackInvMassPt(NULL),
  fHistoSameDiffSameBackInvMassPt(NULL),
  fHistoSameSameDiffBackInvMassPt(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllGammaEMCALAccPt(NULL),
  fHistoMCConvGammaPt(NULL),
  fHistoMCConvGammaR(NULL),
  fHistoMCConvGammaEta(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0WOWeightInAccPt(NULL),
  fHistoMCPi0PtY(NULL),
  fHistoMCPi0PtAlpha(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCGammaFromAllOmegaPt(NULL),
  fHistoMCGammaFromOmegaInAccPt(NULL),
  fHistoMCPi0FromAllOmegaInvMassPt(NULL),
  fHistoMCAllOmegaYPt(NULL),
  fHistoMCOmegaInAccYPt(NULL),
  fHistoMCAllOmegaInvMassPt(NULL),
  fHistoMCOmegaInvMassPt(NULL),
  fHistoMCOmegaInAccInvMassPt(NULL),
  fHistoMCAllOmegaAlphaPt(NULL),
  fHistoMCOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromAllOmegaAlphaPt(NULL),
  fHistoMCPi0FromOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromOmegaInAccInvMassPt(NULL),
  fHistoMCPi0FromAllOmegaYPt(NULL),
  fHistoMCPi0FromOmegaInAccYPt(NULL),
  fHistoMCPi0FromAllOmegaEtaPhi(NULL),
  fHistoMCPi0FromOmegaInAccEtaPhi(NULL),
  fHistoMCOmegaInAccEtaPhi(NULL),
  fHistoMCAllOmegaEtaPhi(NULL),
  fHistoMCAllOmegaPiZeroAnglePt(NULL),
  fHistoMCAllPiZeroGammaAnglePt(NULL),
  fHistoMCAllOmegaGammaAnglePt(NULL),
  fHistoMCInAccOmegaPiZeroAnglePt(NULL),
  fHistoMCInAccPiZeroGammaAnglePt(NULL),
  fHistoMCInAccOmegaGammaAnglePt(NULL),
  fHistoMCAllOmegaPtPi0Pt(NULL),
  fHistoMCInAccOmegaPtPi0Pt(NULL),
  fHistoMCAllOmegaPtGammaPt(NULL),
  fHistoMCInAccOmegaPtGammaPt(NULL),
  fHistoTrueOmegaInvMassPt(NULL),
  fHistoTrueOmegaYPt(NULL),
  fHistoTrueOmegaAlphaPt(NULL),
  fHistoTruePi0FromOmegaInvMassPt(NULL),
  fHistoTruePi0FromOmegaAlphaPt(NULL),
  fHistoTruePi0FromOmegaYPt(NULL),
  fHistoTruePi0FromOmegaEtaPhi(NULL),
  fHistoTruePi0FromOmegaOpenAnglePt(NULL),
  fHistoTrueOmegaPi0AnglePt(NULL),
  fHistoTrueOmegaGammaAnglePt(NULL),
  fHistoTruePi0GammaAnglePt(NULL),
  fHistoTrueOmegaEtaPhi(NULL),
  fHistoTrueOmegaPtPi0Pt(NULL),
  fHistoTrueOmegaPtGammaPt(NULL),
  fHistoTrueGammaFromOmegaPt(NULL),
  fVectorRecTruePi0s(0),
  fVectorDoubleCountTruePi0s(0),
  fHistoMultipleCountTruePi0(NULL),
  fMapMultipleCountTruePi0s(),
  fHistoNEvents(NULL),
  fHistoNEventsMinGamma(NULL),
  fHistoMCOmegaDecayChannels(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNClusterCandidates(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fEventPlaneAngle(-100),
  fRandom(0),
  fNGammaCandidates(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMCStackPos(NULL),
  fMCStackNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fnCuts(0),
  fiCut(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fSetPlotHistsExtQA(kFALSE),
  fWeightJetJetMC(1),
  fEnableSortForClusMC(kFALSE),
  fReconMethod(0)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskOmegaToPiZeroGamma::~AliAnalysisTaskOmegaToPiZeroGamma()
{
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
  }
  if(fClusterCandidates){
    delete fClusterCandidates;
    fClusterCandidates = 0x0;
  }
  if(fBGHandler){
    delete[] fBGHandler;
    fBGHandler = 0x0;
  }
  if(fBGClusHandler){
    delete[] fBGClusHandler;
    fBGClusHandler = 0x0;
  }
}
//___________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::InitBack(){
  
  const Int_t nDim = 4;
  Int_t nBins[nDim] = {800,300,7,4};
  Double_t xMin[nDim] = {0,0, 0,0};
  Double_t xMax[nDim] = {0.8,30,7,4};

  fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];

  fBGClusHandler = new AliGammaConversionAODBGHandler*[fnCuts];
  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->DoBGCalculation()){
      TString cutstringEvent   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPhoton  = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo   = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson   = ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->GetCutNumber();
      
      Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
      Int_t centMin = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
      Int_t centMax = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));
      
      if(collisionSystem == 1 || collisionSystem == 2 ||
        collisionSystem == 5 || collisionSystem == 8 ||
        collisionSystem == 9){
        centMin = centMin*10;
        centMax = centMax*10;
        if(centMax ==0 && centMax!=centMin) centMax=100;
      }else if(collisionSystem == 3 || collisionSystem == 6){
        centMin = centMin*5;
        centMax = centMax*5;
      }else if(collisionSystem == 4 || collisionSystem == 7){
        centMin = ((centMin*5)+45);
        centMax = ((centMax*5)+45);
      }

      fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                collisionSystem,centMin,centMax,
                                ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->GetNumberOfBGEvents(),
                                ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->UseTrackMultiplicity(),
                                2,8,5);
      fBGClusHandler[iCut] = new AliGammaConversionAODBGHandler(
                                collisionSystem,centMin,centMax,
                                ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->GetNumberOfBGEvents(),
                                ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->UseTrackMultiplicity(),
                                2,8,5);
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::UserCreateOutputObjects(){
  
  if (fIsMC > 1){
    fDoPhotonQA       = 0;
  }
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer  = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer  = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  // Array of current cut's gammas
  fGammaCandidates    = new TList();
  fClusterCandidates  = new TList();
  fClusterCandidates->SetOwner(kTRUE);
  
  fCutFolder          = new TList*[fnCuts];
  fESDList            = new TList*[fnCuts];

  fHistoNEvents               = new TH1F*[fnCuts];
  fHistoNEventsMinGamma       = new TH1F*[fnCuts];
  if (fIsMC > 1){
    fHistoNEventsWOWeight     = new TH1F*[fnCuts];
  }
  if (fIsMC == 2){  
    fProfileJetJetXSection    = new TProfile*[fnCuts];
    fHistoJetJetNTrials       = new TH1F*[fnCuts];
  }
  fHistoNGoodESDTracks        = new TH1F*[fnCuts];
  fHistoVertexZ               = new TH1F*[fnCuts];
  fHistoNClusterCandidates    = new TH1F*[fnCuts];
  fHistoNGammaCandidates      = new TH1F*[fnCuts];
  if(fIsHeavyIon==2)
    fProfileEtaShift          = new TProfile*[fnCuts];
  fHistoNGoodESDTracksVsNGammaCandidates  = new TH2F*[fnCuts];
  fHistoSPDClusterTrackletBackground      = new TH2F*[fnCuts];
  fHistoNV0Tracks                         = new TH1F*[fnCuts];
  fHistoConvGammaPt                       = new TH1F*[fnCuts];

  if (fDoPhotonQA > 0){
    fHistoConvGammaR          = new TH1F*[fnCuts];
    fHistoConvGammaEta        = new TH1F*[fnCuts];
  }

  fHistoPhotonPairInvMassPt             = new TH2F*[fnCuts];
  fHistoMotherInvMassPt                 = new TH2F*[fnCuts];
  fHistoGammaFromMotherPt               = new TH1F*[fnCuts];
  if(fReconMethod<2)
    fHistoPhotonPairMatchedInvMassPt    = new TH2F*[fnCuts];
  if(fReconMethod!=2 && fReconMethod!=5)
    fHistoMotherMatchedInvMassPt        = new TH2F*[fnCuts];

  // BG histograms
  fHistoDiff1Diff2SameBackInvMassPt     = new TH2F*[fnCuts];
  fHistoDiffSameSameBackInvMassPt       = new TH2F*[fnCuts];
  // for other cases SameDiffSame is indistinguishable from DiffSameSame
  if(fReconMethod<2)
    fHistoSameDiffSameBackInvMassPt     = new TH2F*[fnCuts];
  fHistoSameSameDiffBackInvMassPt       = new TH2F*[fnCuts];

  // QA histograms
  if(fDoMesonQA>0){
    fHistoPhotonPairYPt                 = new TH2F*[fnCuts];
    fHistoPhotonPairAlphaPt             = new TH2F*[fnCuts];
    fHistoPhotonPairOpenAnglePt         = new TH2F*[fnCuts];
    fHistoPhotonPairEtaPhi              = new TH2F*[fnCuts];
    if(fReconMethod!=2 && fReconMethod!=5)
      fHistoMotherConvPhotonEtaPhi      = new TH2F*[fnCuts];
    fHistoMotherYPt                     = new TH2F*[fnCuts];
    fHistoMotherAlphaPt                 = new TH2F*[fnCuts];
    fHistoMotherEtaPhi                  = new TH2F*[fnCuts];
    fHistoMotherPi0AnglePt              = new TH2F*[fnCuts];
    fHistoMotherGammaAnglePt            = new TH2F*[fnCuts];
    fHistoPi0GammaAnglePt               = new TH2F*[fnCuts];
  }
  
  if(fReconMethod!=5){
    fClusterOutputList                  = new TList*[fnCuts];
    fHistoClusGammaPt                   = new TH1F*[fnCuts];
    fHistoClusOverlapHeadersGammaPt     = new TH1F*[fnCuts];
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPhoton   = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringNeutralPion    = ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson    = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
    fCutFolder[iCut]          = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s_%s",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringNeutralPion.Data(),cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]            = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s_%s_%s ESD histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringNeutralPion.Data(),cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);
    
    fHistoNEvents[iCut]       = new TH1F("NEvents","NEvents",12,-0.5,11.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){ 
      TString TriggerNames    = "Not Trigger: ";
      TriggerNames            = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    }else {
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    fHistoNEventsMinGamma[iCut]     = new TH1F("NEventsMinGamma","NEventsMinGamma",4,-0.5,3.5);
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(1,"2 EMCal 1 PCM");
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(2,"2 PCM 1 EMCal");
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(3,"3 EMCal");
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(4,"3 PCM");
    fESDList[iCut]->Add(fHistoNEventsMinGamma[iCut]);

    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut]   = new TH1F("NEventsWOWeight","NEventsWOWeight",12,-0.5,11.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){ 
        TString TriggerNames        = "Not Trigger: ";
        TriggerNames                = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      }else {
        fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      }
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);
    }
    
    if (fIsMC == 2){  
      fProfileJetJetXSection[iCut]  = new TProfile("XSection","XSection",1,-0.5,0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials","#sum{NTrials}",1,0,1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }
    
    if(fIsHeavyIon == 1) 
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks","GoodESDTracks",4000,0,4000);
    else if(fIsHeavyIon == 2) 
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks","GoodESDTracks",400,0,400);
    else 
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks","GoodESDTracks",200,0,200);
    fHistoNGoodESDTracks[iCut]->SetXTitle("# TPC tracks");
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    fHistoVertexZ[iCut]             = new TH1F("VertexZ","VertexZ",1000,-50,50);
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);

    if(fIsHeavyIon == 1) 
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates","GammaCandidates",100,0,100);
    else if(fIsHeavyIon == 2) 
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates","GammaCandidates",50,0,50);
    else 
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates","GammaCandidates",50,0,50);
    fHistoNGammaCandidates[iCut]->SetXTitle("# accepted #gamma_{conv}");
    fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);

    fHistoNClusterCandidates[iCut]  = new TH1F("ClusterCandidates","ClusterCandidates",50,0,50);
    fHistoNClusterCandidates[iCut]->SetXTitle("# accepted #gamma_{cluster}");
    fESDList[iCut]->Add(fHistoNClusterCandidates[iCut]);
    
    if(fIsHeavyIon == 1)
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,100,0,100);
    else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,50,0,50);
    else
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,50,0,50);
    fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetXTitle("# TPC tracks");
    fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetYTitle("# accepted #gamma_{conv}");
    fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCandidates[iCut]);

    fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
    fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
    else if(fIsHeavyIon == 2)
      fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
    else
      fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
    fHistoNV0Tracks[iCut]->SetXTitle("VZERO amp [arb. units]");
    fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);

    fHistoConvGammaPt[iCut]         = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",500,0,50);
    fHistoConvGammaPt[iCut]->SetXTitle("p_{T,conv}(GeV/c)");
    fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);

    if(fIsHeavyIon == 2){
      fProfileEtaShift[iCut]          = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNEventsMinGamma[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      fHistoNClusterCandidates[iCut]->Sumw2();
      fHistoNGammaCandidates[iCut]->Sumw2();
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Sumw2();
      fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
      fHistoNV0Tracks[iCut]->Sumw2();
      fHistoConvGammaPt[iCut]->Sumw2();
    }
    
    if (fDoPhotonQA > 0){
      fHistoConvGammaR[iCut]        = new TH1F("ESD_ConvGamma_R","ESD_ConvGamma_R",800,0,200);
      fHistoConvGammaR[iCut]->SetXTitle("R_{conv}(cm)");
      fESDList[iCut]->Add(fHistoConvGammaR[iCut]);
      fHistoConvGammaEta[iCut]      = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",2000,-2,2);
      fHistoConvGammaEta[iCut]->SetXTitle("#eta_{conv}");
      fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
    }

    if(fReconMethod!=5){
      fClusterOutputList[iCut]        = new TList();
      fClusterOutputList[iCut]->SetName(Form("%s_%s_%s_%s_%s Cluster Output",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringNeutralPion.Data(),cutstringMeson.Data()));
      fClusterOutputList[iCut]->SetOwner(1);
      fCutFolder[iCut]->Add(fClusterOutputList[iCut]);

      fHistoClusGammaPt[iCut]         = new TH1F("ClusGamma_Pt","ClusGamma_Pt",200,0,20);
      fHistoClusGammaPt[iCut]->SetXTitle("p_{T,clus}(GeV/c)");
      fClusterOutputList[iCut]->Add(fHistoClusGammaPt[iCut]);
      fHistoClusOverlapHeadersGammaPt[iCut]   = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",200,0,20);
      fHistoClusOverlapHeadersGammaPt[iCut]->SetXTitle("p_{T,clus}(GeV/c), selected header w/ overlap");
      fClusterOutputList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);

      if (fIsMC > 1){
        fHistoClusGammaPt[iCut]->Sumw2();
        fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      }
    }
    
    fHistoPhotonPairInvMassPt[iCut]             = new TH2F("ESD_PhotonPair_InvMass_Pt","ESD_PhotonPair_InvMass_Pt",800,0,0.8,200,0,20);
    fHistoPhotonPairInvMassPt[iCut]->SetXTitle("M_{inv, #pi^{0} cand}(GeV/c^{2})");
    fHistoPhotonPairInvMassPt[iCut]->SetYTitle("p_{T, #pi^{0} cand}(GeV/c)");
    fESDList[iCut]->Add(fHistoPhotonPairInvMassPt[iCut]);

    fHistoMotherInvMassPt[iCut]                 = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0.4,1.2,200,0,20);
    fHistoMotherInvMassPt[iCut]->SetXTitle("M_{inv, #Omega cand}(GeV/c^{2})");
    fHistoMotherInvMassPt[iCut]->SetYTitle("p_{T, #Omega cand}(GeV/c)");
    fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);

    fHistoGammaFromMotherPt[iCut] = new TH1F("ESD_GammaFromMother_Pt","ESD_GammaFromMother_Pt",200,0,20);
    fHistoGammaFromMotherPt[iCut]->SetXTitle("p_{T, #gamma from #Omega cand}(GeV/c)");
    fESDList[iCut]->Add(fHistoGammaFromMotherPt[iCut]);

    if(fReconMethod!=2 && fReconMethod!=5){
      fHistoMotherMatchedInvMassPt[iCut]        = new TH2F("ESD_Mother_Matched_InvMass_Pt","ESD_Mother_Matched_InvMass_Pt",800,0.4,1.2,200,0,20);
      fHistoMotherMatchedInvMassPt[iCut]->SetXTitle("M_{inv, #Omega matched}(GeV/c^{2})");
      fHistoMotherMatchedInvMassPt[iCut]->SetYTitle("p_{T, #Omega matched}(GeV/c)");
      fESDList[iCut]->Add(fHistoMotherMatchedInvMassPt[iCut]);
    }

    fHistoDiff1Diff2SameBackInvMassPt[iCut]     = new TH2F("ESD_Omega_Diff1Diff2Same_InvMass_Pt","ESD_Omega_Diff1Diff2Same_InvMass_Pt",800,0.4,1.2,200,0,20);
    fHistoDiff1Diff2SameBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
    fHistoDiff1Diff2SameBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
    fESDList[iCut]->Add(fHistoDiff1Diff2SameBackInvMassPt[iCut]);

    fHistoDiffSameSameBackInvMassPt[iCut]      = new TH2F("ESD_Omega_DiffSameSame_InvMass_Pt","ESD_Omega_DiffSameSame_InvMass_Pt",800,0.4,1.2,200,0,20);
    fHistoDiffSameSameBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
    fHistoDiffSameSameBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
    fESDList[iCut]->Add(fHistoDiffSameSameBackInvMassPt[iCut]);

    if(fReconMethod<2){
      fHistoSameDiffSameBackInvMassPt[iCut]     = new TH2F("ESD_Omega_SameDiffSame_InvMass_Pt","ESD_Omega_SameDiffSame_InvMass_Pt",800,0.4,1.2,200,0,20);
      fHistoSameDiffSameBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
      fHistoSameDiffSameBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      fESDList[iCut]->Add(fHistoSameDiffSameBackInvMassPt[iCut]);
    }

    fHistoSameSameDiffBackInvMassPt[iCut]     = new TH2F("ESD_Omega_SameSameDiff_InvMass_Pt","ESD_Omega_SameSameDiff_InvMass_Pt",800,0.4,1.2,200,0,20);
    fHistoSameSameDiffBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
    fHistoSameSameDiffBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
    fESDList[iCut]->Add(fHistoSameSameDiffBackInvMassPt[iCut]);

    if(fReconMethod<2){
      fHistoPhotonPairMatchedInvMassPt[iCut]      = new TH2F("ESD_PhotonPair_Matched_InvMass_Pt","ESD_PhotonPair_Matched_InvMass_Pt",800,0,0.8,200,0,20);
      fHistoPhotonPairMatchedInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2}) matched conv e^{+/-}to cluster");
      fHistoPhotonPairMatchedInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      fESDList[iCut]->Add(fHistoPhotonPairMatchedInvMassPt[iCut]);
    }

    if (fIsMC > 1){
      fHistoPhotonPairInvMassPt[iCut]->Sumw2();
      if(fReconMethod<2) fHistoPhotonPairMatchedInvMassPt[iCut]->Sumw2();
    }

    if(fDoMesonQA>0){
      fHistoPhotonPairYPt[iCut]              = new TH2F("ESD_PhotonPair_Y_Pt","ESD_PhotonPair_Y_Pt",300,0.03,30.,150,-1.5,1.5);
      fHistoPhotonPairYPt[iCut]->SetXTitle("p_{T, #pi^{0}cand}(GeV/c)");
      fHistoPhotonPairYPt[iCut]->SetYTitle("y_{#pi^{0}cand}");
      SetLogBinningXTH2(fHistoPhotonPairYPt[iCut]);
      fESDList[iCut]->Add(fHistoPhotonPairYPt[iCut]);

      fHistoPhotonPairAlphaPt[iCut]          = new TH2F("ESD_PhotonPair_Alpha_Pt","ESD_PhotonPair_Alpha_Pt",300,0.03,30.,200,-1,1);
      fHistoPhotonPairAlphaPt[iCut]->SetXTitle("p_{T, #pi^{0}cand}(GeV/c)");
      fHistoPhotonPairAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}cand}");
      SetLogBinningXTH2(fHistoPhotonPairAlphaPt[iCut]);
      fESDList[iCut]->Add(fHistoPhotonPairAlphaPt[iCut]);

      fHistoPhotonPairOpenAnglePt[iCut]      = new TH2F("ESD_PhotonPair_OpenAngle_Pt","ESD_PhotonPair_OpenAngle_Pt",300,0.03,30.,100,0,1);
      fHistoPhotonPairOpenAnglePt[iCut]->SetXTitle("p_{T, #pi^{0}cand}(GeV/c)");
      fHistoPhotonPairOpenAnglePt[iCut]->SetYTitle("#theta_{#pi^{0}cand}");
      SetLogBinningXTH2(fHistoPhotonPairOpenAnglePt[iCut]);
      fESDList[iCut]->Add(fHistoPhotonPairOpenAnglePt[iCut]);

      fHistoPhotonPairEtaPhi[iCut] = new TH2F("ESD_PhotonPair_Eta_Phi","ESD_PhotonPair_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
      fHistoPhotonPairEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}cand}(rad)");
      fHistoPhotonPairEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}cand}");
      fESDList[iCut]->Add(fHistoPhotonPairEtaPhi[iCut]);

      if(fReconMethod!=2 && fReconMethod!=5){
        fHistoMotherConvPhotonEtaPhi[iCut] = new TH2F("ESD_MotherConvPhoton_Eta_Phi","ConvPhoton under #Omega peak",600,0,2*TMath::Pi(),200,-1,1);
        fHistoMotherConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}}(rad)");
        fHistoMotherConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
        fESDList[iCut]->Add(fHistoMotherConvPhotonEtaPhi[iCut]);
      }

      fHistoMotherYPt[iCut]                       = new TH2F("ESD_Mother_Y_Pt","ESD_Mother_Y_Pt",200,0,20,150,-1.5,1.5);
      fHistoMotherYPt[iCut]->SetYTitle("y_{#Omega cand}");
      fHistoMotherYPt[iCut]->SetXTitle("p_{T, #Omega cand}(GeV/c)");
      fESDList[iCut]->Add(fHistoMotherYPt[iCut]);

      fHistoMotherAlphaPt[iCut]                   = new TH2F("ESD_Mother_Alpha_Pt","ESD_Mother_Alpha_Pt",200,0,20,200,-1,1);
      fHistoMotherAlphaPt[iCut]->SetXTitle("p_{T, #Omega cand}(GeV/c)");
      fHistoMotherAlphaPt[iCut]->SetYTitle("#alpha_{#Omega cand}");
      fESDList[iCut]->Add(fHistoMotherAlphaPt[iCut]);

      fHistoMotherEtaPhi[iCut]                  = new TH2F("ESD_Mother_Eta_Phi","ESD_Mother_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
      fHistoMotherEtaPhi[iCut]->SetXTitle("#phi_{#Omega cand}(rad)");
      fHistoMotherEtaPhi[iCut]->SetYTitle("#eta_{#Omega cand}");
      fESDList[iCut]->Add(fHistoMotherEtaPhi[iCut]);

      fHistoMotherPi0AnglePt[iCut] = new TH2F("ESD_MotherPi0_Angle_Pt","ESD_MotherPi0_Angle_Pt",200,0,20,360,0,TMath::Pi());
      fHistoMotherPi0AnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
      fHistoMotherPi0AnglePt[iCut]->SetYTitle("#theta_{#Omega cand, pi^{0}cand}");
      fESDList[iCut]->Add(fHistoMotherPi0AnglePt[iCut]);

      fHistoMotherGammaAnglePt[iCut] = new TH2F("ESD_MotherGamma_Angle_Pt","ESD_MotherGamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
      fHistoMotherGammaAnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
      fHistoMotherGammaAnglePt[iCut]->SetYTitle("#theta_{#Omega cand, #gamma}");
      fESDList[iCut]->Add(fHistoMotherGammaAnglePt[iCut]);

      fHistoPi0GammaAnglePt[iCut] = new TH2F("ESD_Pi0Gamma_Angle_Pt","ESD_Pi0Gamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
      fHistoPi0GammaAnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
      fHistoPi0GammaAnglePt[iCut]->SetYTitle("#theta_{#pi^{0}cand,#gamma}");
      fESDList[iCut]->Add(fHistoPi0GammaAnglePt[iCut]);

      if (fIsMC > 1){
        fHistoPhotonPairYPt[iCut]->Sumw2();
        fHistoPhotonPairAlphaPt[iCut]->Sumw2();
        fHistoPhotonPairOpenAnglePt[iCut]->Sumw2();
        fHistoMotherConvPhotonEtaPhi[iCut]->Sumw2();
        fHistoMotherYPt[iCut]->Sumw2();
        fHistoMotherPi0AnglePt[iCut]->Sumw2();
        fHistoMotherGammaAnglePt[iCut]->Sumw2();
        fHistoPi0GammaAnglePt[iCut]->Sumw2();
        fHistoGammaFromMotherPt[iCut]->Sumw2();
      }
    }
  }

  InitBack(); // Init Background Handler
  
  if(fIsMC>0){
    // MC Histogramms
    fMCList                                         = new TList*[fnCuts];
    // True Histogramms
    fTrueList                                       = new TList*[fnCuts];

    fHistoMCAllGammaPt                              = new TH1F*[fnCuts];
    fHistoMCAllGammaEMCALAccPt                      = new TH1F*[fnCuts];
    fHistoMCConvGammaPt                             = new TH1F*[fnCuts];
    fHistoMCOmegaInAccInvMassPt                     = new TH2F*[fnCuts];
    fHistoMCOmegaInvMassPt                          = new TH2F*[fnCuts];
    fHistoMCPi0FromOmegaInAccInvMassPt              = new TH2F*[fnCuts];
    fHistoMCAllOmegaInvMassPt                       = new TH2F*[fnCuts];
    fHistoMCPi0FromAllOmegaInvMassPt                = new TH2F*[fnCuts];
    fHistoMCPi0Pt                                   = new TH1F*[fnCuts];
    fHistoMCPi0WOWeightPt                           = new TH1F*[fnCuts];
    fHistoMCPi0InAccPt                              = new TH1F*[fnCuts];
    fHistoMCPi0WOWeightInAccPt                      = new TH1F*[fnCuts];
    if (fIsMC > 1){
      fHistoMCPi0WOEvtWeightPt                      = new TH1F*[fnCuts];
    }
    fHistoMCGammaFromAllOmegaPt                     = new TH1F*[fnCuts];
    fHistoMCGammaFromOmegaInAccPt                   = new TH1F*[fnCuts];
    fHistoMCOmegaDecayChannels                      = new TH1F*[fnCuts];

    if(fDoMesonQA>0){
      fHistoMCAllOmegaYPt                             = new TH2F*[fnCuts];
      fHistoMCOmegaInAccYPt                           = new TH2F*[fnCuts];
      fHistoMCAllOmegaAlphaPt                         = new TH2F*[fnCuts];
      fHistoMCOmegaInAccAlphaPt                       = new TH2F*[fnCuts];
      fHistoMCPi0FromAllOmegaAlphaPt                  = new TH2F*[fnCuts];
      fHistoMCPi0FromOmegaInAccAlphaPt                = new TH2F*[fnCuts];
      fHistoMCPi0FromAllOmegaYPt                      = new TH2F*[fnCuts];
      fHistoMCPi0FromOmegaInAccYPt                    = new TH2F*[fnCuts];
      fHistoMCPi0FromAllOmegaEtaPhi                   = new TH2F*[fnCuts];
      fHistoMCPi0FromOmegaInAccEtaPhi                 = new TH2F*[fnCuts];
      fHistoMCAllOmegaEtaPhi                          = new TH2F*[fnCuts];
      fHistoMCOmegaInAccEtaPhi                        = new TH2F*[fnCuts];
      fHistoMCAllOmegaPiZeroAnglePt                   = new TH2F*[fnCuts];
      fHistoMCAllPiZeroGammaAnglePt                   = new TH2F*[fnCuts];
      fHistoMCAllOmegaGammaAnglePt                    = new TH2F*[fnCuts];
      fHistoMCInAccOmegaPiZeroAnglePt                 = new TH2F*[fnCuts];
      fHistoMCInAccPiZeroGammaAnglePt                 = new TH2F*[fnCuts];
      fHistoMCInAccOmegaGammaAnglePt                  = new TH2F*[fnCuts];
      fHistoMCAllOmegaPtPi0Pt                         = new TH2F*[fnCuts];
      fHistoMCInAccOmegaPtPi0Pt                       = new TH2F*[fnCuts];
      fHistoMCAllOmegaPtGammaPt                       = new TH2F*[fnCuts];
      fHistoMCInAccOmegaPtGammaPt                     = new TH2F*[fnCuts];
      fHistoMCPi0PtY                                  = new TH2F*[fnCuts];
      fHistoMCPi0PtAlpha                              = new TH2F*[fnCuts];
      if (fIsMC == 2){
        fHistoMCPi0PtJetPt                            = new TH2F*[fnCuts];
      }
    }

    if (fDoPhotonQA > 0){
      fHistoMCConvGammaR                            = new TH1F*[fnCuts];
      fHistoMCConvGammaEta                          = new TH1F*[fnCuts];
    }

    fHistoTrueOmegaInvMassPt                      = new TH2F*[fnCuts];
    fHistoTruePi0FromOmegaInvMassPt               = new TH2F*[fnCuts];
    fHistoTrueGammaFromOmegaPt                    = new TH1F*[fnCuts];

    if(fDoMesonQA>0){
      fHistoTruePi0FromOmegaYPt                     = new TH2F*[fnCuts];
      fHistoTrueOmegaYPt                            = new TH2F*[fnCuts];
      fHistoTrueOmegaAlphaPt                        = new TH2F*[fnCuts];
      fHistoTruePi0FromOmegaAlphaPt                 = new TH2F*[fnCuts];
      fHistoTrueOmegaPi0AnglePt                     = new TH2F*[fnCuts];
      fHistoTrueOmegaGammaAnglePt                   = new TH2F*[fnCuts];
      fHistoTruePi0GammaAnglePt                     = new TH2F*[fnCuts];
      fHistoTrueOmegaEtaPhi                         = new TH2F*[fnCuts];
      fHistoTruePi0FromOmegaEtaPhi                  = new TH2F*[fnCuts];
      fHistoTruePi0FromOmegaOpenAnglePt             = new TH2F*[fnCuts];
      fHistoTrueOmegaPtPi0Pt                        = new TH2F*[fnCuts];
      fHistoTrueOmegaPtGammaPt                      = new TH2F*[fnCuts];
    }

    fHistoMultipleCountTruePi0                    = new TH1F*[fnCuts];
    
    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPhoton   = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringNeutralPion    = ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson    = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      
      fMCList[iCut]                     = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringNeutralPion.Data(),cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);

      fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",200,0,20);
      fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
      fHistoMCAllGammaEMCALAccPt[iCut]  = new TH1F("MC_AllGammaEMCALAcc_Pt","MC_AllGammaEMCALAcc_Pt",200,0,20);
      fMCList[iCut]->Add(fHistoMCAllGammaEMCALAccPt[iCut]);
      fHistoMCConvGammaPt[iCut]         = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",500,0,50);
      fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);

      if (fIsMC > 1){
        fHistoMCAllGammaPt[iCut]->Sumw2();
        fHistoMCAllGammaEMCALAccPt[iCut]->Sumw2();
        fHistoMCConvGammaPt[iCut]->Sumw2();
      }

      if (fDoPhotonQA > 0){
        fHistoMCConvGammaR[iCut]        = new TH1F("MC_ConvGamma_R","MC_ConvGamma_R",800,0,200);
        fMCList[iCut]->Add(fHistoMCConvGammaR[iCut]);
        fHistoMCConvGammaEta[iCut]      = new TH1F("MC_ConvGamma_Eta","MC_ConvGamma_Eta",2000,-2,2);
        fMCList[iCut]->Add(fHistoMCConvGammaEta[iCut]);
      }
      
      fHistoMCPi0Pt[iCut]               = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",200,0,20);
      fHistoMCPi0Pt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
      fHistoMCPi0WOWeightPt[iCut]       = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",200,0,20);
      fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);

      fHistoMCPi0InAccPt[iCut]          = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",200,0,20);
      fHistoMCPi0InAccPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
      fHistoMCPi0WOWeightInAccPt[iCut]  = new TH1F("MC_Pi0WOWeightInAcc_Pt","MC_Pi0WOWeightInAcc_Pt",200,0,20);
      fMCList[iCut]->Add(fHistoMCPi0WOWeightInAccPt[iCut]);

      fHistoMCGammaFromAllOmegaPt[iCut]    = new TH1F("MC_GammaFromAllOmega_Pt","MC_GammaFromAllOmega_Pt",200,0,20);
      fHistoMCGammaFromAllOmegaPt[iCut]->SetXTitle("p_{T, #gamma}(GeV/c)");
      fHistoMCGammaFromAllOmegaPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCGammaFromAllOmegaPt[iCut]);

      fHistoMCGammaFromOmegaInAccPt[iCut]    = new TH1F("MC_GammaFromOmegaInAcc_Pt","MC_GammaFromOmegaInAcc_Pt",200,0,20);
      fHistoMCGammaFromOmegaInAccPt[iCut]->SetXTitle("p_{T, #gamma}(GeV/c)");
      fHistoMCGammaFromOmegaInAccPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCGammaFromOmegaInAccPt[iCut]);

      fHistoMCOmegaDecayChannels[iCut]     = new TH1F("MC_OmegaDecayChannels","MC_OmegaDecayChannels",3,-0.5,2.5);
      fHistoMCOmegaDecayChannels[iCut]->GetXaxis()->SetBinLabel(1,"PiPlPiMiPiZero");
      fHistoMCOmegaDecayChannels[iCut]->GetXaxis()->SetBinLabel(2,"PiZeroGamma");
      fHistoMCOmegaDecayChannels[iCut]->GetXaxis()->SetBinLabel(3,"Others");
      fMCList[iCut]->Add(fHistoMCOmegaDecayChannels[iCut]);

      fHistoMCAllOmegaInvMassPt[iCut]                = new TH2F("MC_AllOmega_InvMass_Pt","MC_AllOmega_InvMass_Pt",800,0.4,1.2,200,0,20);
      fHistoMCAllOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#Omega}(GeV/c^{2})");
      fHistoMCAllOmegaInvMassPt[iCut]->SetYTitle("#Omega p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCAllOmegaInvMassPt[iCut]);

      fHistoMCPi0FromAllOmegaInvMassPt[iCut]                = new TH2F("MC_Pi0FromAllOmega_InvMass_Pt","MC_Pi0FromAllOmega_InvMass_Pt",800,0,0.8,200,0,20);
      fHistoMCPi0FromAllOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
      fHistoMCPi0FromAllOmegaInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaInvMassPt[iCut]);

      fHistoMCPi0FromOmegaInAccInvMassPt[iCut]              = new TH2F("MC_Pi0FromOmegaInAcc_InvMass_Pt","MC_Pi0FromOmegaInAcc_InvMass_Pt",800,0,0.8,200,0,20);
      fHistoMCPi0FromOmegaInAccInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
      fHistoMCPi0FromOmegaInAccInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccInvMassPt[iCut]);

      fHistoMCOmegaInvMassPt[iCut] = new TH2F("MC_OmegaInvMass_Pt","MC_OmegaInvMass_Pt",800,0.4,1.2,200,0,20);
      fHistoMCOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#Omega}(GeV/c^{2})");
      fHistoMCOmegaInvMassPt[iCut]->SetYTitle("#Omega p_{T}(GeV/c)");
      fHistoMCOmegaInvMassPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCOmegaInvMassPt[iCut]);

      fHistoMCOmegaInAccInvMassPt[iCut] = new TH2F("MC_OmegaInAcc_InvMass_Pt","MC_OmegaInAcc_InvMass_Pt",800,0.4,1.2,200,0,20);
      fHistoMCOmegaInAccInvMassPt[iCut]->SetXTitle("M_{inv,#Omega}(GeV/c^{2})");
      fHistoMCOmegaInAccInvMassPt[iCut]->SetYTitle("#Omega p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCOmegaInAccInvMassPt[iCut]);

      if(fDoMesonQA>0){
        fHistoMCAllOmegaYPt[iCut] = new TH2F("MC_AllOmega_Y_Pt","MC_AllOmega_Y_Pt",200,0,20,150,-1.5,1.5);
        fHistoMCAllOmegaYPt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoMCAllOmegaYPt[iCut]->SetYTitle("Y_{#Omega}");
        fHistoMCAllOmegaYPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCAllOmegaYPt[iCut]);

        fHistoMCOmegaInAccYPt[iCut] = new TH2F("MC_OmegaInAcc_Y_Pt","MC_OmegaInAcc_Y_Pt",200,0,20,150,-1.5,1.5);
        fHistoMCOmegaInAccYPt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoMCOmegaInAccYPt[iCut]->SetYTitle("Y_{#Omega}");
        fHistoMCOmegaInAccYPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCOmegaInAccYPt[iCut]);

        fHistoMCAllOmegaAlphaPt[iCut] = new TH2F("MC_AllOmega_Alpha_Pt","MC_AllOmega_Alpha_Pt",200,0,20,200,-1,1);
        fHistoMCAllOmegaAlphaPt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoMCAllOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#Omega}");
        fMCList[iCut]->Add(fHistoMCAllOmegaAlphaPt[iCut]);

        fHistoMCOmegaInAccAlphaPt[iCut] = new TH2F("MC_OmegaInAcc_Alpha_Pt","MC_OmegaInAcc_Alpha_Pt",200,0,20,200,-1,1);
        fHistoMCOmegaInAccAlphaPt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoMCOmegaInAccAlphaPt[iCut]->SetYTitle("#alpha_{#Omega}");
        fMCList[iCut]->Add(fHistoMCOmegaInAccAlphaPt[iCut]);

        fHistoMCPi0FromAllOmegaAlphaPt[iCut] = new TH2F("MC_Pi0FromAllOmega_Alpha_Pt","MC_Pi0FromAllOmega_Alpha_Pt",200,0,20,200,-1,1);
        fHistoMCPi0FromAllOmegaAlphaPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
        fHistoMCPi0FromAllOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}}");
        fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaAlphaPt[iCut]);

        fHistoMCPi0FromOmegaInAccAlphaPt[iCut] = new TH2F("MC_Pi0FromOmegaInAcc_Alpha_Pt","MC_Pi0FromOmegaInAcc_Alpha_Pt",200,0,20,200,-1,1);
        fHistoMCPi0FromOmegaInAccAlphaPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
        fHistoMCPi0FromOmegaInAccAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}}");
        fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccAlphaPt[iCut]);

        fHistoMCPi0FromAllOmegaYPt[iCut] = new TH2F("MC_Pi0FromAllOmega_Y_Pt","MC_Pi0FromAllOmega_Y_Pt",200,0,20,150,-1.5,1.5);
        fHistoMCPi0FromAllOmegaYPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
        fHistoMCPi0FromAllOmegaYPt[iCut]->SetYTitle("Y_{#pi^{0}}");
        fHistoMCPi0FromAllOmegaYPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaYPt[iCut]);

        fHistoMCPi0FromOmegaInAccYPt[iCut] = new TH2F("MC_Pi0FromOmegaInAcc_Y_Pt","MC_Pi0FromOmegaInAcc_Y_Pt",200,0,20,150,-1.5,1.5);
        fHistoMCPi0FromOmegaInAccYPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
        fHistoMCPi0FromOmegaInAccYPt[iCut]->SetYTitle("Y_{#pi^{0}}");
        fHistoMCPi0FromOmegaInAccYPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccYPt[iCut]);

        fHistoMCPi0FromAllOmegaEtaPhi[iCut] = new TH2F("MC_Pi0FromAllOmega_Eta_Phi","MC_Pi0FromAllOmega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoMCPi0FromAllOmegaEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}}(rad)");
        fHistoMCPi0FromAllOmegaEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}}");
        fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaEtaPhi[iCut]);

        fHistoMCPi0FromOmegaInAccEtaPhi[iCut] = new TH2F("MC_Pi0FromOmegaInAcc_Eta_Phi","MC_Pi0FromOmegaInAcc_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoMCPi0FromOmegaInAccEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}}(rad)");
        fHistoMCPi0FromOmegaInAccEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}}");
        fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccEtaPhi[iCut]);

        fHistoMCAllOmegaEtaPhi[iCut]   = new TH2F("MC_AllOmega_Eta_Phi","MC_AllOmega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoMCAllOmegaEtaPhi[iCut]->SetXTitle("#phi_{#Omega}(rad)");
        fHistoMCAllOmegaEtaPhi[iCut]->SetYTitle("#eta_{#Omega}");
        fMCList[iCut]->Add(fHistoMCAllOmegaEtaPhi[iCut]);

        fHistoMCOmegaInAccEtaPhi[iCut]   = new TH2F("MC_OmegaInAcc_Eta_Phi","MC_OmegaInAcc_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoMCOmegaInAccEtaPhi[iCut]->SetXTitle("#phi_{#Omega}(rad)");
        fHistoMCOmegaInAccEtaPhi[iCut]->SetYTitle("#eta_{#Omega}");
        fMCList[iCut]->Add(fHistoMCOmegaInAccEtaPhi[iCut]);

        fHistoMCAllOmegaPiZeroAnglePt[iCut] = new TH2F("MC_AllOmegaPiZero_Angle_Pt","MC_AllOmegaPiZero_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoMCAllOmegaPiZeroAnglePt[iCut]->SetYTitle("#theta_{#Omega,#pi^{0}}");
        fHistoMCAllOmegaPiZeroAnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fMCList[iCut]->Add(fHistoMCAllOmegaPiZeroAnglePt[iCut]);

        fHistoMCAllPiZeroGammaAnglePt[iCut] = new TH2F("MC_AllPiZeroGamma_Angle_Pt","MC_AllPiZeroGamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoMCAllPiZeroGammaAnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoMCAllPiZeroGammaAnglePt[iCut]->SetYTitle("#theta_{#pi^{0},#gamma}");
        fMCList[iCut]->Add(fHistoMCAllPiZeroGammaAnglePt[iCut]);

        fHistoMCAllOmegaGammaAnglePt[iCut] = new TH2F("MC_AllOmegaGamma_Angle_Pt","MC_AllOmegaGamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoMCAllOmegaGammaAnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoMCAllOmegaGammaAnglePt[iCut]->SetYTitle("#theta_{#Omega,#gamma}");
        fMCList[iCut]->Add(fHistoMCAllOmegaGammaAnglePt[iCut]);

        fHistoMCInAccOmegaPiZeroAnglePt[iCut] = new TH2F("MC_InAccOmegaPiZero_Angle_Pt","MC_InAccOmegaPiZero_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoMCInAccOmegaPiZeroAnglePt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoMCInAccOmegaPiZeroAnglePt[iCut]->SetYTitle("#theta_{#Omega,#pi^{0}}");
        fMCList[iCut]->Add(fHistoMCInAccOmegaPiZeroAnglePt[iCut]);

        fHistoMCInAccPiZeroGammaAnglePt[iCut] = new TH2F("MC_InAccPiZeroGamma_Angle_Pt","MC_InAccPiZeroGamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoMCInAccPiZeroGammaAnglePt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoMCInAccPiZeroGammaAnglePt[iCut]->SetYTitle("#theta_{#pi^{0},#gamma}");
        fMCList[iCut]->Add(fHistoMCInAccPiZeroGammaAnglePt[iCut]);

        fHistoMCInAccOmegaGammaAnglePt[iCut] = new TH2F("MC_InAccOmegaGamma_Angle_Pt","MC_InAccOmegaGamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoMCInAccOmegaGammaAnglePt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoMCInAccOmegaGammaAnglePt[iCut]->SetYTitle("#theta_{#Omega,#gamma}");
        fMCList[iCut]->Add(fHistoMCInAccOmegaGammaAnglePt[iCut]);

        fHistoMCAllOmegaPtPi0Pt[iCut] = new TH2F("MC_All_OmegaPt_Pi0Pt","MC_All_OmegaPt_Pi0Pt",200,0,20,200,0,20);
        fHistoMCAllOmegaPtPi0Pt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoMCAllOmegaPtPi0Pt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
        fMCList[iCut]->Add(fHistoMCAllOmegaPtPi0Pt[iCut]);

        fHistoMCInAccOmegaPtPi0Pt[iCut] = new TH2F("MC_InAcc_OmegaPt_Pi0Pt","MC_InAcc_OmegaPt_Pi0Pt",200,0,20,200,0,20);
        fHistoMCInAccOmegaPtPi0Pt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoMCInAccOmegaPtPi0Pt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
        fMCList[iCut]->Add(fHistoMCInAccOmegaPtPi0Pt[iCut]);

        fHistoMCAllOmegaPtGammaPt[iCut] = new TH2F("MC_All_OmegaPt_GammaPt","MC_All_OmegaPt_GammaPt",200,0,20,200,0,20);
        fHistoMCAllOmegaPtGammaPt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoMCAllOmegaPtGammaPt[iCut]->SetYTitle("#gamma p_{T}(GeV/c)");
        fMCList[iCut]->Add(fHistoMCAllOmegaPtGammaPt[iCut]);

        fHistoMCInAccOmegaPtGammaPt[iCut] = new TH2F("MC_InAcc_OmegaPt_GammaPt","MC_InAcc_OmegaPt_GammaPt",200,0,20,200,0,20);
        fHistoMCInAccOmegaPtGammaPt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoMCInAccOmegaPtGammaPt[iCut]->SetYTitle("#gamma p_{T}(GeV/c)");
        fMCList[iCut]->Add(fHistoMCInAccOmegaPtGammaPt[iCut]);

        fHistoMCPi0PtY[iCut]            = new TH2F("MC_Pi0_Pt_Y","MC_Pi0_Pt_Y",300,0.03,30.,150,-1.5,1.5);
        fHistoMCPi0PtY[iCut]->Sumw2();
        SetLogBinningXTH2(fHistoMCPi0PtY[iCut]);
        fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);

        fHistoMCPi0PtAlpha[iCut]        = new TH2F("MC_Pi0_Pt_Alpha","MC_Pi0_Pt_Alpha",300,0.03,30.,200,-1,1);
        SetLogBinningXTH2(fHistoMCPi0PtAlpha[iCut]);
        fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);

        if (fIsMC >= 2){
          fHistoMCPi0PtAlpha[iCut]->Sumw2();
        }
      }

      if (fIsMC > 1){
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fHistoMCPi0WOWeightInAccPt[iCut]->Sumw2();
        fHistoMCPi0WOEvtWeightPt[iCut]  = new TH1F("MC_Pi0_WOEventWeights_Pt","MC_Pi0_WOEventWeights_Pt",400,0,40);
        fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);

        if (fDoMesonQA > 0 && fIsMC == 2){
          fHistoMCPi0PtJetPt[iCut]      = new TH2F("MC_Pi0_Pt_JetPt","MC_Pi0_Pt_JetPt",300,0.03,30.,200,0,200);
          fHistoMCPi0PtJetPt[iCut]->Sumw2();
          SetLogBinningXTH2(fHistoMCPi0PtJetPt[iCut]);
          fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
        }
      }

      fTrueList[iCut]                           = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringNeutralPion.Data(),cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      fHistoTrueOmegaInvMassPt[iCut]                = new TH2F("True_Omega_InvMass_Pt","True_Omega_InvMass_Pt",800,0.4,1.2,200,0,20);
      fHistoTrueOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#Omega}(GeV/c^{2})");
      fHistoTrueOmegaInvMassPt[iCut]->SetYTitle("#Omega p_{T}(GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueOmegaInvMassPt[iCut]);

      fHistoTruePi0FromOmegaInvMassPt[iCut]          = new TH2F("True_Pi0FromOmega_InvMass_Pt","True_Pi0FromOmega_InvMass_Pt",800,0,0.8,200,0,20);
      fHistoTruePi0FromOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
      fHistoTruePi0FromOmegaInvMassPt[iCut]->SetYTitle("#Omega p_{T}(GeV/c)");
      fTrueList[iCut]->Add(fHistoTruePi0FromOmegaInvMassPt[iCut]);

      fHistoTrueGammaFromOmegaPt[iCut] = new TH1F("True_GammaFromOmega_Pt","True_GammaFromOmega_Pt",200,0,20);
      fHistoTrueGammaFromOmegaPt[iCut]->SetXTitle("p_{T,#gamma}(GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueGammaFromOmegaPt[iCut]);

      if(fDoMesonQA>0){
        fHistoTrueOmegaYPt[iCut] = new TH2F("True_Omega_Y_Pt","True_Omega_Y_Pt",200,0,20,150,-1.5,1.5);
        fHistoTrueOmegaYPt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoTrueOmegaYPt[iCut]->SetYTitle("Y_{#Omega}");
        fHistoTrueOmegaYPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueOmegaYPt[iCut]);

        fHistoTrueOmegaAlphaPt[iCut] = new TH2F("True_Omega_Alpha_Pt","True_Omega_Alpha_Pt",200,0,20,200,-1,1);
        fHistoTrueOmegaAlphaPt[iCut]->SetXTitle("p_{T,#Omega}(GeV/c)");
        fHistoTrueOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#Omega}");
        fHistoTrueOmegaAlphaPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueOmegaAlphaPt[iCut]);

        fHistoTruePi0FromOmegaYPt[iCut] = new TH2F("True_Pi0FromOmega_Y_Pt","True_Pi0FromOmega_Y_Pt",200,0,20,150,-1.5,1.5);
        fHistoTruePi0FromOmegaYPt[iCut]->SetXTitle("p_{T,pi^{0}}(GeV/c)");
        fHistoTruePi0FromOmegaYPt[iCut]->SetYTitle("Y_{#pi^{0}}");
        fHistoTruePi0FromOmegaYPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePi0FromOmegaYPt[iCut]);

        fHistoTruePi0FromOmegaAlphaPt[iCut] = new TH2F("True_Pi0FromOmega_Alpha_Pt","MC_Pi0FromOmega_Alpha_Pt",200,0,20,200,-1,1);
        fHistoTruePi0FromOmegaAlphaPt[iCut]->SetXTitle("p_{T,pi^{0}}(GeV/c)");
        fHistoTruePi0FromOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}}");
        fTrueList[iCut]->Add(fHistoTruePi0FromOmegaAlphaPt[iCut]);

        fHistoTrueOmegaPi0AnglePt[iCut] = new TH2F("True_OmegaPiZero_Angle_Pt","True_OmegaPiZero_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoTrueOmegaPi0AnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoTrueOmegaPi0AnglePt[iCut]->SetYTitle("#theta_{#Omega,#pi^{0}}");
        fTrueList[iCut]->Add(fHistoTrueOmegaPi0AnglePt[iCut]);

        fHistoTruePi0GammaAnglePt[iCut] = new TH2F("True_PiZeroGamma_Angle_Pt","True_PiZeroGamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoTruePi0GammaAnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoTruePi0GammaAnglePt[iCut]->SetYTitle("#theta_{#pi_{0},#gamma}");
        fTrueList[iCut]->Add(fHistoTruePi0GammaAnglePt[iCut]);

        fHistoTrueOmegaGammaAnglePt[iCut] = new TH2F("True_OmegaGamma_Angle_Pt","True_OmegaGamma_Angle_Pt",200,0,20,360,0,TMath::Pi());
        fHistoTrueOmegaGammaAnglePt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoTrueOmegaGammaAnglePt[iCut]->SetYTitle("#theta_{#Omega,#gamma}");
        fTrueList[iCut]->Add(fHistoTrueOmegaGammaAnglePt[iCut]);

        fHistoTrueOmegaEtaPhi[iCut]   = new TH2F("True_Omega_Eta_Phi","True_Omega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoTrueOmegaEtaPhi[iCut]->SetXTitle("#phi_{#Omega}(rad)");
        fHistoTrueOmegaEtaPhi[iCut]->SetYTitle("#eta_{#Omega}");
        fTrueList[iCut]->Add(fHistoTrueOmegaEtaPhi[iCut]);

        fHistoTruePi0FromOmegaEtaPhi[iCut]   = new TH2F("True_Pi0FromOmega_Eta_Phi","True_Pi0FromOmega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoTruePi0FromOmegaEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}}(rad)");
        fHistoTruePi0FromOmegaEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}}");
        fTrueList[iCut]->Add(fHistoTruePi0FromOmegaEtaPhi[iCut]);

        fHistoTruePi0FromOmegaOpenAnglePt[iCut] = new TH2F("True_Pi0FromOmega_OpenAngle_Pt","True_Pi0FromOmega_OpenAngle_Pt",200,0,20,100,0,1);
        fHistoTruePi0FromOmegaOpenAnglePt[iCut]->SetXTitle("p_{T, #pi^{0}}(GeV/c)");
        fHistoTruePi0FromOmegaOpenAnglePt[iCut]->SetYTitle("#theta_{#pi^{0}}");
        fTrueList[iCut]->Add(fHistoTruePi0FromOmegaOpenAnglePt[iCut]);

        fHistoTrueOmegaPtPi0Pt[iCut] = new TH2F("True_OmegaPt_Pi0Pt","True_OmegaPt_Pi0Pt",200,0,20,200,0,20);
        fHistoTrueOmegaPtPi0Pt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoTrueOmegaPtPi0Pt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueOmegaPtPi0Pt[iCut]);

        fHistoTrueOmegaPtGammaPt[iCut] = new TH2F("True_OmegaPt_GammaPt","True_OmegaPt_GammaPt",200,0,20,200,0,20);
        fHistoTrueOmegaPtGammaPt[iCut]->SetXTitle("#Omega p_{T}(GeV/c)");
        fHistoTrueOmegaPtGammaPt[iCut]->SetYTitle("#gamma p_{T}(GeV/c)");
        fTrueList[iCut]->Add(fHistoTrueOmegaPtGammaPt[iCut]);
      }

      fHistoMultipleCountTruePi0[iCut]            = new TH1F("ESD_TrueMultipleCountPi0","ESD_TrueMultipleCountPi0",10,1,11);
      fTrueList[iCut]->Add(fHistoMultipleCountTruePi0[iCut]);

      if (fIsMC > 1){
        fHistoMultipleCountTruePi0[iCut]->Sumw2();
      }
    }
  }

  fVectorDoubleCountTruePi0s.clear();

  fMapMultipleCountTruePi0s.clear();

  fVectorRecTruePi0s.clear();
    
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;}// GetV0Reader
  
  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
    if (fV0Reader->GetV0FindingEfficiencyHistograms())
      fOutputContainer->Add(fV0Reader->GetV0FindingEfficiencyHistograms());


      
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
    if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
    }
    if(fSetPlotHistsExtQA){
      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms()){
        fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms());
      }
    }
    if(!((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))) continue;
    if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
    }
  }
  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskOmegaToPiZeroGamma::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  

    if(fIsHeavyIon == 2){
      if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        continue; // No Eta Shift requested, continue
      }
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
        continue;
      } else{
        printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
            (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      }
    }
  }
  
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  if(fIsMC > 0) fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = 0;

  fInputEvent = InputEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;// incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC > 1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }

  if(fIsMC>0 && fInputEvent->IsA()==AliESDEvent::Class()){
    fMCStack = fMCEvent->Stack();
    if(fMCStack == NULL) fIsMC = 0;
  }
  
  if(fInputEvent->IsA()==AliAODEvent::Class()){
    fInputEvent->InitMagneticField();
  }
  
  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

  // ------------------- BeginEvent ----------------------------
  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;
  
  if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kTRUE);// In case of AODMC relabeling MC
    fV0Reader->RelabelAODs(kTRUE);
  }
  
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    
    fiCut = iCut;
//     cout << ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber() << "_" <<  ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber() << 
//             "_" << ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber() <<endl;
    
    Bool_t isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,isRunningEMCALrelAna);
    
    if(fIsMC==2){
      Float_t xsection      = -1.; 
      Float_t ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials);
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
    }

    if(fIsMC>0){
      fWeightJetJetMC       = 1;
  //     cout << fMCEvent << endl;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC );
      if (fIsMC == 3){
        Double_t weightMult   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }
      if(fIsMC==1) fWeightJetJetMC = 1;
      if (!isMCJet){
        fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(10);
        continue;
      }
    }
    
    Bool_t triggered = kTRUE;
    
    if(eventNotAccepted!= 0){
    // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      if (eventNotAccepted==3 && fIsMC > 0){
        triggered = kFALSE;
      }else {  
        continue;
      }
    }

    if(eventQuality != 0){// Event Not Accepted
      //cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
      continue;
    }

    if (triggered==kTRUE){
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

      fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
      fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
      fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)),fWeightJetJetMC);
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)  fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(), fWeightJetJetMC);
      else fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(), fWeightJetJetMC);
    }
    
    if(fIsMC>0){
      // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fMCEvent);
        }
        else if(fInputEvent->IsA()==AliAODEvent::Class()){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                          ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                          fInputEvent);
        }
      }
    }
    if(fIsMC>0){
      if(fInputEvent->IsA()==AliESDEvent::Class())
        ProcessMCParticles();
      if(fInputEvent->IsA()==AliAODEvent::Class())
        ProcessAODMCParticles();
    }

    if (triggered==kFALSE) continue;
    
    // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
    if(fReconMethod!=5) ProcessClusters();// process calo clusters
    if(fReconMethod!=2) ProcessPhotonCandidates(); // Process this cuts gammas

    fHistoNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries(),fWeightJetJetMC);
    fHistoNClusterCandidates[iCut]->Fill(fClusterCandidates->GetEntries(),fWeightJetJetMC);
    fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fGammaCandidates->GetEntries(),fWeightJetJetMC);
    if(fClusterCandidates->GetEntries()>1 && fGammaCandidates->GetEntries()>0){
      fHistoNEventsMinGamma[iCut]->Fill(0);
    }
    if(fClusterCandidates->GetEntries()>0 && fGammaCandidates->GetEntries()>1){
      fHistoNEventsMinGamma[iCut]->Fill(1);
    }
    if(fClusterCandidates->GetEntries()>2){
      fHistoNEventsMinGamma[iCut]->Fill(2);
    }
    if(fGammaCandidates->GetEntries()>2){
      fHistoNEventsMinGamma[iCut]->Fill(3);
    }

    // Meson Analysis
    if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->UseMCPSmearing() && fIsMC>0){
      fUnsmearedPx = new Double_t[fGammaCandidates->GetEntries()]; // Store unsmeared Momenta
      fUnsmearedPy = new Double_t[fGammaCandidates->GetEntries()];
      fUnsmearedPz = new Double_t[fGammaCandidates->GetEntries()];
      fUnsmearedE =  new Double_t[fGammaCandidates->GetEntries()];

      for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
      fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Px();
      fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Py();
      fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Pz();
      fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->E();
      ((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(gamma)));
      }
    }

    CalculateOmegaCandidates(); // Combine Gammas from conversion and from calo
    CalculateBackground(); // Combinatorial Background
    UpdateEventByEventData(); // Store Event for mixed Events

    if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class()){
      ProcessConversionPhotonsForMissingTagsAOD(); //Count missing tags
    }else if (fIsMC>0 && fInputEvent->IsA()==AliESDEvent::Class()){
      ProcessConversionPhotonsForMissingTags(); //Count missing tags
    }

    if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->UseMCPSmearing() && fIsMC>0){
      for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
      ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
      ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPy(fUnsmearedPy[gamma]);
      ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPz(fUnsmearedPz[gamma]);
      ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetE(fUnsmearedE[gamma]);
      }
      delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
      delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
      delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
      delete[] fUnsmearedE;fUnsmearedE  = 0x0;
    }

    if(fIsMC>0){
      fVectorRecTruePi0s.clear();
      fVectorDoubleCountTruePi0s.clear();
      FillMultipleCountHistoAndClear(fMapMultipleCountTruePi0s,fHistoMultipleCountTruePi0[iCut]);
    }

    fGammaCandidates->Clear(); // delete this cuts good gammas
    fClusterCandidates->Clear(); // delete cluster candidates
  }

  if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }
  
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessClusters(){
  
  Int_t nclus = 0;
  nclus = fInputEvent->GetNumberOfCaloClusters();
  
//   cout << nclus << endl;
  
  if(nclus == 0)  return;

  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);
  
  if(fReconMethod==2){
    // match tracks to clusters
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC);
  }

  // vertex
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  // Loop over EMCal clusters
  for(Int_t i = 0; i < nclus; i++){    
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

    if (!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fIsMC,fWeightJetJetMC)){ delete clus; continue;}

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);
    
    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
    
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}
    
    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton();
    PhotonCandidate->SetCaloClusterRef((Long_t)i);
    // get MC label
    if(fIsMC>0){
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
//       cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k<(Int_t)clus->GetNLabels(); k++){
          if (k<50)PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
//           Int_t pdgCode = fMCStack->Particle(mclabelsCluster[k])->GetPdgCode();
//           cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
        }
      }
    }
    
    fIsFromMBHeader         = kTRUE; 
    fIsOverlappingWithOtherHeader   = kFALSE;
    //TString periodName         = fV0Reader->GetPeriodName();
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCStack, fInputEvent) == 0){
        fIsFromMBHeader = kFALSE;
      }
      if (clus->GetNLabels()>1){
        Int_t* mclabelsCluster = clus->GetLabels();
        for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCStack, fInputEvent) == 0){
            fIsOverlappingWithOtherHeader = kTRUE;
          }
        }
      }
    }
    
    if (fIsOverlappingWithOtherHeader) 
      fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
    
    if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
      if(fIsMC>0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ProcessTrueClusterCandidates(PhotonCandidate);
        }else {
          ProcessTrueClusterCandidatesAOD(PhotonCandidate);
        }
      }
      fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
      fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
    }else{
      delete PhotonCandidate;
    }

    delete clus;
    delete tmpvec;
  } 
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{

  TParticle *Photon = NULL;
  if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");

  if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = fMCStack->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;

  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCStack, fEnableSortForClusMC);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  AliAODMCParticle *Photon = NULL;
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray){
    if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
      else return;
  }else {
    AliInfo("AODMCTrackArray could not be loaded");
    return;
  }

  if(Photon == NULL){
  //  cout << "no photon" << endl;
    return;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fInputEvent, fEnableSortForClusMC);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessPhotonCandidates()
{
  Int_t nV0 = 0;
  TList *GammaCandidatesStepOne = new TList();
  TList *GammaCandidatesStepTwo = new TList();
  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    fIsFromMBHeader = kTRUE;
    if(fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
      if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
      if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
    }
    
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
    !((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
      
      if(fIsFromMBHeader){
        fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
        if (fDoPhotonQA > 0){
          fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
          fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
        }
      }
      if(fIsMC>0){
        if(fInputEvent->IsA()==AliESDEvent::Class())
        ProcessTruePhotonCandidates(PhotonCandidate);
        if(fInputEvent->IsA()==AliAODEvent::Class())
        ProcessTruePhotonCandidatesAOD(PhotonCandidate);
      }

    }else if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
      nV0++;
      GammaCandidatesStepOne->Add(PhotonCandidate);
    }else if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
        ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
      GammaCandidatesStepTwo->Add(PhotonCandidate);
    }
  }
  if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
    for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
      AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
      if(!PhotonCandidate) continue;
      fIsFromMBHeader = kTRUE;
      if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
        fGammaCandidates->Add(PhotonCandidate);
        if(fIsFromMBHeader){
          fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
          if (fDoPhotonQA > 0){
            fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
            fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
          }
        }
      }
      if(fIsMC>0){
        if(fInputEvent->IsA()==AliESDEvent::Class())
          ProcessTruePhotonCandidates(PhotonCandidate);
        if(fInputEvent->IsA()==AliAODEvent::Class())
          ProcessTruePhotonCandidatesAOD(PhotonCandidate);
      }else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
    }
  }
  if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
      if(!PhotonCandidate) continue;
      fIsFromMBHeader = kTRUE;
      if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
      fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
      if(fIsFromMBHeader){
        fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
        if (fDoPhotonQA > 0){
          fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
          fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
        }
      }
      if(fIsMC>0){
        if(fInputEvent->IsA()==AliESDEvent::Class())
          ProcessTruePhotonCandidates(PhotonCandidate);
        if(fInputEvent->IsA()==AliAODEvent::Class())
          ProcessTruePhotonCandidatesAOD(PhotonCandidate);
      }
    }
  }
  
  delete GammaCandidatesStepOne;
  GammaCandidatesStepOne = 0x0;
  delete GammaCandidatesStepTwo;
  GammaCandidatesStepTwo = 0x0;
  
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
  AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};

  if(posDaughter->GetMother() != negDaughter->GetMother()) return;

  else if(posDaughter->GetMother() == -1) return;

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

  AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
  if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

  if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5) return; // check if the daughters come from a conversion

  TruePhotonCandidate->SetIsTrueConvertedPhoton();
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  // Process True Photons
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCStack);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCStack);

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};
  if(posDaughter->GetMother(0) != negDaughter->GetMother(0)) return;
  else if(posDaughter->GetMother(0) == -1) return;

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

  TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCStack);

  if(Photon->GetPdgCode() != 22){
    return; // Mother is no Photon
  }

  if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

  // True Photon
  TruePhotonCandidate->SetIsTrueConvertedPhoton();
  return;
}
//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessAODMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  
  // Loop over all primary MC particle
  for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
    
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
    if (!particle) continue;

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (!isPrimary) continue;
    
    Int_t isMCFromMBHeader = -1;
    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
      if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
    }
    
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kFALSE)){
      fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
      if (fabs(particle->Eta()) < 0.66 ){
        if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
      }
    }
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kTRUE)){
      Double_t rConv = 0;
      for(Int_t daughterIndex=particle->GetDaughter(0);daughterIndex<=particle->GetDaughter(1);daughterIndex++){
        AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(daughterIndex));
        if(!tmpDaughter) continue;
        if(abs(tmpDaughter->GetPdgCode()) == 11){
          rConv = sqrt( (tmpDaughter->Xv()*tmpDaughter->Xv()) + (tmpDaughter->Yv()*tmpDaughter->Yv()) );
        }
      }
      fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
      if (fDoPhotonQA > 0){
        fHistoMCConvGammaR[fiCut]->Fill(rConv);
        fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta());
      }
    }
    // Converted MC Gamma
    if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
      ->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
      AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(0)));
      AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(1)));
      Float_t weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
        if (particle->Pt()>0.005){
          weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
          //                   if(particle->GetPdgCode() == 221){
          //                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
          //                   }
        }
      }
      Double_t mesonY = 10.;
      if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
        mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
      }else{
        mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
      }

      Double_t alpha = -10;
      if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
        alpha = (daughter0->E() - daughter1->E())/(daughter0->E() + daughter1->E());
      }


      if(particle->GetPdgCode() == 111){
        fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
        fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
        if (fDoMesonQA > 0){
          fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC);
          fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);
          if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
        }
      }

      // Check the acceptance for both gammas
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
        // check acceptance of clusters as well, true if one of them points into the Calo acceptance
        if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,AODMCTrackArray) ||
          ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,AODMCTrackArray) ){
          if(particle->GetPdgCode() == 111){
            fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc
            fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc wo weight
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
//   cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;
  
  // Loop over all primary MC particles
  for(UInt_t i = 0; i < fMCStack->GetNtrack(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){ 
      // fill primary histograms
      TParticle* particle = (TParticle *)fMCStack->Particle(i);
      if (!particle) continue;

      // fill histograms for all true omegas
      if(particle->GetPdgCode() == 223){

        Bool_t DecaysToPiPlPiMiPiZero = kFALSE;

        TParticle *posPion = 0x0;
        TParticle *negPion = 0x0;
        TParticle *neutPion = 0x0;

        for(Int_t index= particle->GetFirstDaughter();index <= particle->GetLastDaughter();index++){

          TParticle *temp = (TParticle*)fMCStack->Particle(index);

          switch(temp->GetPdgCode()) {
          case 211:
            posPion      =  temp;
            break;
          case -211:
            negPion      =  temp;
            break;
          case 111:
            neutPion    =  temp;
            break;
          }
        }

        if(posPion && negPion && neutPion) DecaysToPiPlPiMiPiZero = kTRUE;

        if(DecaysToPiPlPiMiPiZero){
          fHistoMCOmegaDecayChannels[fiCut]->Fill(0);
        } else if (particle->GetNDaughters()==2){
          TParticle *gamma2 = 0x0;
          TParticle *pi0 = 0x0;

          for(Int_t index = particle->GetFirstDaughter();index <= particle->GetLastDaughter();index++){
            TParticle *temp = (TParticle*)fMCStack->Particle(index);
            switch(temp->GetPdgCode()){
            case 22:
              gamma2 = temp;
              break;
            case 111:
              pi0   = temp;
              break;
            }
          }

          if(gamma2 && pi0){
            fHistoMCOmegaDecayChannels[fiCut]->Fill(1);
            fHistoMCAllOmegaInvMassPt[fiCut]->Fill(TMath::Sqrt((particle->Energy())*(particle->Energy())-(particle->P())*(particle->P())),particle->Pt(),fWeightJetJetMC);
            fHistoMCGammaFromAllOmegaPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
            fHistoMCPi0FromAllOmegaInvMassPt[fiCut]->Fill(TMath::Sqrt((pi0->Energy())*(pi0->Energy())-(pi0->P())*(pi0->P())),pi0->Pt(),fWeightJetJetMC);
            if(fDoMesonQA>0){

              fHistoMCPi0FromAllOmegaEtaPhi[fiCut]->Fill(pi0->Phi(),pi0->Eta(),fWeightJetJetMC);
              fHistoMCAllOmegaPtPi0Pt[fiCut]->Fill(particle->Pt(),pi0->Pt(),fWeightJetJetMC);
              fHistoMCAllOmegaPtGammaPt[fiCut]->Fill(particle->Pt(),gamma2->Pt(),fWeightJetJetMC);

              Double_t alpha = (pi0->Energy() - gamma2->Energy())/(pi0->Energy() + gamma2->Energy());
              fHistoMCAllOmegaAlphaPt[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);

              Double_t Pi0Y = 10.;
              if(pi0->Energy() - pi0->Pz() == 0 || pi0->Energy() + pi0->Pz() == 0){
                Pi0Y=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                Pi0Y = 0.5*(TMath::Log((pi0->Energy()+pi0->Pz()) / (pi0->Energy()-pi0->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }
              fHistoMCPi0FromAllOmegaYPt[fiCut]->Fill(pi0->Pt(),Pi0Y,fWeightJetJetMC);

              Double_t OmegaY = 10.;
              if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
                OmegaY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                OmegaY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }
              fHistoMCAllOmegaYPt[fiCut]->Fill(particle->Pt(),OmegaY,fWeightJetJetMC);

              fHistoMCAllOmegaEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);

              //get angles and plot
              TVector3 omegavector = TVector3(particle->Px(),particle->Py(),particle->Pz());
              TVector3 pi0vector = TVector3(pi0->Px(),pi0->Py(),pi0->Pz());
              TVector3 gamma2vector = TVector3(gamma2->Px(),gamma2->Py(),gamma2->Pz());
              fHistoMCAllOmegaPiZeroAnglePt[fiCut]->Fill(particle->Pt(),TMath::Pi() - pi0vector.Angle(omegavector));
              fHistoMCAllPiZeroGammaAnglePt[fiCut]->Fill(particle->Pt(),pi0vector.Angle(gamma2vector));
              fHistoMCAllOmegaGammaAnglePt[fiCut]->Fill(particle->Pt(),TMath::Pi() - omegavector.Angle(gamma2vector));

              //check whether pi0 decayed into two gammas
              if (pi0->GetNDaughters()==2){
                TParticle *gamma0 = (TParticle*)fMCStack->Particle(pi0->GetFirstDaughter());
                TParticle *gamma1 = (TParticle*)fMCStack->Particle(pi0->GetLastDaughter());
                if (gamma0->GetPdgCode()==22 && gamma1->GetPdgCode()==22){

                  //plot pi0 alpha
                  Double_t pi0alpha = (gamma0->Energy() - gamma1->Energy())/(gamma0->Energy() + gamma1->Energy());
                  fHistoMCPi0FromAllOmegaAlphaPt[fiCut]->Fill(pi0->Pt(),pi0alpha,fWeightJetJetMC);
                }
              }
            }
          } else{
            fHistoMCOmegaDecayChannels[fiCut]->Fill(2);
          }
        } else{
          fHistoMCOmegaDecayChannels[fiCut]->Fill(2);
        }
      }

      //fill histograms for omegas in acceptance
      Int_t labelNeutPion = -1;
      Int_t labelGamma = -1;

      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
         ->MesonIsSelectedMCPiZeroGamma(particle,fMCStack,labelNeutPion,labelGamma,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        TParticle *neutPion    = fMCStack->Particle(labelNeutPion);

        //fill histograms for acceptance correction
        fHistoMCOmegaInvMassPt[fiCut]->Fill(TMath::Sqrt((particle->Energy())*(particle->Energy())-(particle->P())*(particle->P())),particle->Pt(),fWeightJetJetMC);

        if (neutPion->GetNDaughters()==2){

          TParticle *gamma0 = (TParticle*)fMCStack->Particle(neutPion->GetFirstDaughter());
          TParticle *gamma1 = (TParticle*)fMCStack->Particle(neutPion->GetLastDaughter());
          TParticle *gamma2 = (TParticle*)fMCStack->Particle(labelGamma);

          Bool_t InAcceptance = kFALSE;

          if (fReconMethod%2==1){ //cases 1,3,5 where reconstruction requires gamma2 to be a pcm photon
            if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma2,fMCStack,kFALSE)){ //check that gamma2 is in acceptance
              if(fReconMethod==5){
                //check that both gamma0 and gamma1 are in acceptance
                if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCStack,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCStack,kFALSE)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==3){
                if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCStack) &&
                   ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCStack)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==1){ // both gammas must be in tpc acceptance
                if((((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCStack,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCStack,kFALSE)) &&
                   // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                   (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCStack) ||
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCStack))){
                  InAcceptance = kTRUE;
                }
              }
            }
          } else{ //cases 0,2,4 where reconstruction requires gamma2 to be a calo photon
            if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma2,fMCStack)){
              if(fReconMethod==4){
                //check that both gamma0 and gamma1 are in acceptance
                if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCStack,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCStack,kFALSE)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==2){
                if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCStack) &&
                   ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCStack)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==0){ //either gamma0 is in tpc acc & gamma1 is in emcal acc or vice versa
                if((((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCStack,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCStack,kFALSE)) &&
                   // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                   (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCStack) ||
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCStack))){
                  InAcceptance = kTRUE;
                }
              }
            }
          }
          if(InAcceptance){
            // fill in acceptance histograms
            fHistoMCOmegaInAccInvMassPt[fiCut]->Fill(TMath::Sqrt((particle->Energy())*(particle->Energy())-(particle->P())*(particle->P())),particle->Pt(),fWeightJetJetMC);
            fHistoMCPi0FromOmegaInAccInvMassPt[fiCut]->Fill(TMath::Sqrt((neutPion->Energy())*(neutPion->Energy())-(neutPion->P())*(neutPion->P())),neutPion->Pt(),fWeightJetJetMC);
            fHistoMCGammaFromOmegaInAccPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);

            if(fDoMesonQA>0){

              TVector3 omegavector = TVector3(particle->Px(),particle->Py(),particle->Pz());
              TVector3 pi0vector = TVector3(neutPion->Px(),neutPion->Py(),neutPion->Pz());
              TVector3 gamma2vector = TVector3(gamma2->Px(),gamma2->Py(),gamma2->Pz());

              Double_t OmegaPiZeroAngle = TMath::Pi() - pi0vector.Angle(omegavector);
              Double_t PiZeroGammaAngle = pi0vector.Angle(gamma2vector);
              Double_t OmegaGammaAngle = TMath::Pi() - omegavector.Angle(gamma2vector);

              Double_t OmegaY = 10.;
              if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
                OmegaY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                OmegaY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }
              Double_t Pi0Y = 10.;
              if(neutPion->Energy() - neutPion->Pz() == 0 || neutPion->Energy() + neutPion->Pz() == 0){
                Pi0Y=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                Pi0Y = 0.5*(TMath::Log((neutPion->Energy()+neutPion->Pz()) / (neutPion->Energy()-neutPion->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }

              Double_t OmegaAlpha = (neutPion->Energy() - gamma2->Energy())/(neutPion->Energy() + gamma2->Energy());
              Double_t Pi0Alpha = (gamma0->Energy() - gamma1->Energy())/(gamma0->Energy() + gamma1->Energy());

              fHistoMCOmegaInAccEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
              fHistoMCOmegaInAccYPt[fiCut]->Fill(particle->Pt(),OmegaY,fWeightJetJetMC);
              fHistoMCOmegaInAccAlphaPt[fiCut]->Fill(particle->Pt(),OmegaAlpha,fWeightJetJetMC);
              fHistoMCInAccOmegaPiZeroAnglePt[fiCut]->Fill(particle->Pt(),OmegaPiZeroAngle);
              fHistoMCInAccPiZeroGammaAnglePt[fiCut]->Fill(particle->Pt(),PiZeroGammaAngle);
              fHistoMCInAccOmegaGammaAnglePt[fiCut]->Fill(particle->Pt(),OmegaGammaAngle);
              fHistoMCPi0FromOmegaInAccAlphaPt[fiCut]->Fill(neutPion->Pt(),Pi0Alpha,fWeightJetJetMC);
              fHistoMCPi0FromOmegaInAccYPt[fiCut]->Fill(neutPion->Pt(),Pi0Y,fWeightJetJetMC);
              fHistoMCPi0FromOmegaInAccEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
              fHistoMCInAccOmegaPtPi0Pt[fiCut]->Fill(particle->Pt(),neutPion->Pt(),fWeightJetJetMC);
              fHistoMCInAccOmegaPtGammaPt[fiCut]->Fill(particle->Pt(),gamma2->Pt(),fWeightJetJetMC);
            }
          }
        }
      }


      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kFALSE)){
        fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
        if (fabs(particle->Eta()) < 0.66 ){
          if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        }
      }
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kTRUE)){
        fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        if (fDoPhotonQA > 0){
          fHistoMCConvGammaR[fiCut]->Fill(((TParticle*)fMCStack->Particle(particle->GetFirstDaughter()))->R());
          fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta());
        }
      }// Converted MC Gamma
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
        ->MesonIsSelectedMC(particle,fMCStack,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        TParticle* daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
        TParticle* daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCStack, fInputEvent);
            //                   if(particle->GetPdgCode() == 221){
            //                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
            //                   }
          }
        }
        Double_t mesonY = 10.;
        if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
          mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }else{
          mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }

        Double_t alpha = -10;
        if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
          alpha = (daughter0->Energy() - daughter1->Energy())/(daughter0->Energy() + daughter1->Energy());
        }

        if(particle->GetPdgCode() == 111){
          fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
          fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0){
            fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC); // All MC Pi0
            fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
            if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }
        }
        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t InAcceptance = kFALSE;

        if(fReconMethod/2 == 0){
          if( kDaughter0IsPrim && kDaughter1IsPrim &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCStack,kFALSE) &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCStack,kFALSE)){
            // check acceptance of clusters as well, true if one of them points into the Calo acceptance
            if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCStack) ||
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCStack)){
              InAcceptance = kTRUE;
            }
          }
        } else if(fReconMethod/2 == 1){
          if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCStack) &&
             ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCStack)){
            InAcceptance = kTRUE;
          }
        } else if(fReconMethod/2 == 2){
          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCStack,kFALSE) &&
             ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCStack,kFALSE)){
            InAcceptance = kTRUE;
          }
        }
        if(InAcceptance){
          // fill in acceptance histograms
          if(particle->GetPdgCode() == 111){
            fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Pi0 with gamma in acc
            fHistoMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Pi0 with gamma in acc wo weighting
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::CalculateOmegaCandidates()
{
  switch(fReconMethod){
  //PCM-cal,cal
  case 0:
    if(fClusterCandidates->GetEntries()>1 && fGammaCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;

        for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
          Bool_t matched = kFALSE;
          AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
          matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster,fInputEvent,fWeightJetJetMC);

          AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);

          if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if (matched){
              fHistoPhotonPairMatchedInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            } else{
              fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
              if (pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow() &&
                  pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                if(fDoMesonQA>0){
                  fHistoPhotonPairYPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
                  fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
                  fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand->Phi(),pi0cand->Eta(),fWeightJetJetMC);
                }
                // get third photon from clusters and calculate inv mass of omega
                for(Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
                  if (thirdGammaIndex==secondGammaIndex) continue;
                  Bool_t matchedgamma2wconvgamma = kFALSE;
                  AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
                  if (gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
                  AliVCluster* cluster2 = fInputEvent->GetCaloCluster(gamma2->GetCaloClusterRef());
                  matchedgamma2wconvgamma = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster2, fInputEvent, fWeightJetJetMC);
                  AliAODConversionMother *omegacand = new AliAODConversionMother(pi0cand,gamma2);
                  if (matchedgamma2wconvgamma){
                    fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
                    delete omegacand;
                    omegacand=0x0;
                    continue;
                  }
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                      ProcessTrueMesonCandidates(omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                      ProcessTrueMesonCandidatesAOD(omegacand,gamma0,gamma1,gamma2);
                  }
                  if(omegacand->M()>0.7 && omegacand->M()<0.85){
                    fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                    if(fDoMesonQA>0){
                      fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
                      fHistoMotherYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                      fHistoMotherAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
                      fHistoMotherEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
                      fHistoMotherPi0AnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
                      fHistoMotherGammaAnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
                      fHistoPi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
                    }
                  }
                  delete omegacand;
                  omegacand=0x0;
                }
              }
            }
          }
          delete pi0cand;
          pi0cand=0x0;
        }
      }
    }
    break;
  //PCM-cal,PCM
  case 1:
    if(fGammaCandidates->GetEntries()>1 && fClusterCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;

        for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          Bool_t matched = kFALSE;
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;

          AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
          matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC);

          AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);

          if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if (matched){
              fHistoPhotonPairMatchedInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            } else{
              // fill photon pair histograms
              fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
              if (pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow() &&
                  pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                if(fDoMesonQA>0){
                  fHistoPhotonPairYPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
                  fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
                  fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand->Phi(),pi0cand->Eta(),fWeightJetJetMC);
                }
                // get third photon from conversion photon candidates and calculate inv mass of omega
                for(Int_t thirdGammaIndex=0;thirdGammaIndex<fGammaCandidates->GetEntries();thirdGammaIndex++){
                  if (thirdGammaIndex==firstGammaIndex) continue;
                  Bool_t matchedgamma1wconvgamma2 = kFALSE;
                  AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(thirdGammaIndex));
                  if (gamma2==NULL) continue;
                  matchedgamma1wconvgamma2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2, cluster, fInputEvent, fWeightJetJetMC);
                  AliAODConversionMother *omegacand = new AliAODConversionMother(pi0cand,gamma2);
                  if (matchedgamma1wconvgamma2){
                    fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
                    delete omegacand;
                    omegacand=0x0;
                    continue;
                  }
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                      ProcessTrueMesonCandidates(omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                      ProcessTrueMesonCandidatesAOD(omegacand,gamma0,gamma1,gamma2);
                  }
                  if(omegacand->M()>0.7 && omegacand->M()<0.85){
                    fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                    if(fDoMesonQA>0){
                      fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
                      fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma2->GetPhotonPhi(), gamma2->GetPhotonEta(),fWeightJetJetMC);
                      fHistoMotherYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                      fHistoMotherAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
                      fHistoMotherEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
                      fHistoMotherPi0AnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
                      fHistoMotherGammaAnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
                      fHistoPi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
                    }
                  }
                  delete omegacand;
                  omegacand=0x0;
                }
              }
            }
          }
          delete pi0cand;
          pi0cand=0x0;
        }
      }
    }
    break;
  //cal-cal,cal
  case 2:
    if(fClusterCandidates->GetEntries()>2){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
        if (gamma0==NULL || !(gamma0->GetIsCaloPhoton())) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
          AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0, gamma1);
          if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            if(pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow() ||
               pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){delete pi0cand; pi0cand=0x0; continue;}
            if(fDoMesonQA>0){
              fHistoPhotonPairYPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
              fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
              fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
              fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand->Phi(),pi0cand->Eta(),fWeightJetJetMC);
            }
            for(Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
              if (thirdGammaIndex==secondGammaIndex || thirdGammaIndex==firstGammaIndex) continue;
              AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
              if (gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
              AliAODConversionMother *omegacand = new AliAODConversionMother(pi0cand,gamma2);
              fHistoMotherInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
              if(fIsMC>0){
                if(fInputEvent->IsA()==AliESDEvent::Class())
                  ProcessTrueMesonCandidates(omegacand,gamma0,gamma1,gamma2);
                if(fInputEvent->IsA()==AliAODEvent::Class())
                  ProcessTrueMesonCandidatesAOD(omegacand,gamma0,gamma1,gamma2);
              }
              if(omegacand->M()>0.7 && omegacand->M()<0.85){
                fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                if(fDoMesonQA>0){
                  fHistoMotherYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  fHistoMotherAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
                  fHistoMotherEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
                  fHistoMotherPi0AnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
                  fHistoMotherGammaAnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
                  fHistoPi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
                }
              }
              delete omegacand;
              omegacand=0x0;
            }
            delete pi0cand;
            pi0cand=0x0;
          }
        }
      }
    }
    break;
  //cal-cal,PCM
  case 3:
    if(fClusterCandidates->GetEntries()>1 && fGammaCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
        if (gamma0==NULL || !(gamma0->GetIsCaloPhoton())) continue;
        AliVCluster* cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
          AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0, gamma1);
          if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            if(pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow() ||
               pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){delete pi0cand; pi0cand=0x0; continue;}
            if(fDoMesonQA>0){
              fHistoPhotonPairYPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
              fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
              fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
              fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand->Phi(),pi0cand->Eta(),fWeightJetJetMC);
            }
            AliVCluster* cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            for(Int_t thirdGammaIndex=0;thirdGammaIndex<fGammaCandidates->GetEntries();thirdGammaIndex++){
              AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(thirdGammaIndex));
              Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2,cluster0, fInputEvent, fWeightJetJetMC))
                  || (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2,cluster1, fInputEvent, fWeightJetJetMC));
              AliAODConversionMother *omegacand = new AliAODConversionMother(pi0cand,gamma2);
              if (matched){
                fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
                delete omegacand;
                omegacand=0x0;
                continue;
              }
              fHistoMotherInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
              if(fIsMC>0){
                if(fInputEvent->IsA()==AliESDEvent::Class())
                  ProcessTrueMesonCandidates(omegacand,gamma0,gamma1,gamma2);
                if(fInputEvent->IsA()==AliAODEvent::Class())
                  ProcessTrueMesonCandidatesAOD(omegacand,gamma0,gamma1,gamma2);
              }
              if(omegacand->M()>0.7 && omegacand->M()<0.85){
                fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                if(fDoMesonQA>0){
                  fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma2->GetPhotonPhi(), gamma2->GetPhotonEta(),fWeightJetJetMC);
                  fHistoMotherYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  fHistoMotherAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
                  fHistoMotherEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
                  fHistoMotherPi0AnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
                  fHistoMotherGammaAnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
                  fHistoPi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
                }
              }
              delete omegacand;
              omegacand=0x0;
            }
            delete pi0cand;
            pi0cand=0x0;
          }
        }
      }
    }
    break;
  //PCM-PCM,cal
  case 4:
    if(fGammaCandidates->GetEntries()>1 && fClusterCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if(gamma0==NULL) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
          if(gamma1==NULL) continue;
          AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
          if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            if(pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow() ||
               pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){delete pi0cand; pi0cand=0x0; continue;}
            if(fDoMesonQA>0){
              fHistoPhotonPairYPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
              fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
              fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
              fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand->Phi(),pi0cand->Eta(),fWeightJetJetMC);
            }
            for(Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
              AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
              if(gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
              AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma2->GetCaloClusterRef());
              Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC))
                  || (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma1,cluster, fInputEvent, fWeightJetJetMC));
              AliAODConversionMother *omegacand = new AliAODConversionMother(pi0cand,gamma2);
              if (matched){
                fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
                delete omegacand;
                omegacand=0x0;
                continue;
              }
              fHistoMotherInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
              if(fIsMC>0){
                if(fInputEvent->IsA()==AliESDEvent::Class())
                  ProcessTrueMesonCandidates(omegacand,gamma0,gamma1,gamma2);
                if(fInputEvent->IsA()==AliAODEvent::Class())
                  ProcessTrueMesonCandidatesAOD(omegacand,gamma0,gamma1,gamma2);
              }
              if(omegacand->M()>0.7 && omegacand->M()<0.85){
                fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                if(fDoMesonQA){
                  fHistoMotherYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  fHistoMotherAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
                  fHistoMotherEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
                  fHistoMotherPi0AnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
                  fHistoMotherGammaAnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
                  fHistoPi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
                  fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
                  fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma1->GetPhotonPhi(), gamma1->GetPhotonEta(),fWeightJetJetMC);
                }
              }
              delete omegacand;
              omegacand=0x0;
            }
            delete pi0cand;
            pi0cand=0x0;
          }
        }
      }
    }
    break;
  //PCM-PCM,PCM
  case 5:
    if(fGammaCandidates->GetEntries()>2){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
          if (gamma1==NULL || (secondGammaIndex+1)>=(fGammaCandidates->GetEntries())) continue;
          AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0, gamma1);
          if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
            fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(),fWeightJetJetMC);
            if(pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow() ||
               pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){delete pi0cand; pi0cand=0x0; continue;}
            if(fDoMesonQA>0){
              fHistoPhotonPairYPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
              fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
              fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
              fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand->Phi(),pi0cand->Eta(),fWeightJetJetMC);
            }
            for(Int_t thirdGammaIndex=0;thirdGammaIndex<fGammaCandidates->GetEntries();thirdGammaIndex++){
              if (thirdGammaIndex==firstGammaIndex || thirdGammaIndex==secondGammaIndex) continue;
              AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(thirdGammaIndex));
              if (gamma2==NULL) continue;
              AliAODConversionMother *omegacand = new AliAODConversionMother(pi0cand,gamma2);
              if(fIsMC>0){
                if(fInputEvent->IsA()==AliESDEvent::Class())
                  ProcessTrueMesonCandidates(omegacand,gamma0,gamma1,gamma2);
                if(fInputEvent->IsA()==AliAODEvent::Class())
                  ProcessTrueMesonCandidatesAOD(omegacand,gamma0,gamma1,gamma2);
              }
              fHistoMotherInvMassPt[fiCut]->Fill(omegacand->M(),omegacand->Pt(),fWeightJetJetMC);
              if(omegacand->M()>0.7 && omegacand->M()<0.85){
                fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                if(fDoMesonQA>0){
                  fHistoMotherYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  fHistoMotherAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
                  fHistoMotherEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
                  fHistoMotherPi0AnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
                  fHistoMotherGammaAnglePt[fiCut]->Fill(omegacand->Pt(),omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
                  fHistoPi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
                }
              }
              delete omegacand;
              omegacand=0x0;
            }
            delete pi0cand;
            pi0cand=0x0;
          }
        }
      }
    }
    break;
  }
}
//______________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTrueMesonCandidates(AliAODConversionMother *OmegaCandidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, AliAODConversionPhoton *TrueGammaCandidate2)
{
  Bool_t isTruePi0 = kFALSE;
  Bool_t isTrueOmega = kFALSE;

  switch (fReconMethod){
    // pcm-cal,cal
    case 0:
    {
      // get gamma0MotherLabel
      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCStack);
        TParticle * gammaMC0 = (TParticle*)fMCStack->Particle(gamma0MCLabel);
        gamma0MotherLabel=gammaMC0->GetFirstMother();
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;
      TParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        // get mother or grandmother of gamma1 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC1 = (TParticle*)fMCStack->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=(fMCStack->Particle(gammaMC1->GetMother(0)))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // fill pi0 histograms here if necessary

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;
        TParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          // get mother or grandmother of gamma1 (potentially true omega) depending on whether it is an electron-leading/photon-leading cluster
          gammaMC2 = (TParticle*)fMCStack->Particle(gamma2MCLabel);
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){  // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother(0);
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()) gamma2MotherLabel=fMCStack->Particle(gammaMC2->GetMother(0))->GetMother(0);
            else gamma2MotherLabel=gammaMC2->GetMother(0);
          }
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCStack->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCStack->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){ // also need to consider matching cluster of conversion positron to pcm tracks, not yet implemented
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother *TruePi0 = new AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(fDoMesonQA>0){
            fHistoTrueOmegaYPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTrueOmegaAlphaPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->GetAlpha(),fWeightJetJetMC);
            fHistoTrueOmegaEtaPhi[fiCut]->Fill(OmegaCandidate->Phi(),OmegaCandidate->Eta(),fWeightJetJetMC);
            fHistoTrueOmegaGammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtGammaPt[fiCut]->Fill(OmegaCandidate->Pt(),TrueGammaCandidate2->Pt(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaAlphaPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetAlpha(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaYPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaEtaPhi[fiCut]->Fill(TruePi0->Phi(),TruePi0->Eta(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaOpenAnglePt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetOpeningAngle(),fWeightJetJetMC);
            fHistoTrueOmegaPi0AnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TruePi0->Vect()),fWeightJetJetMC);
            fHistoTruePi0GammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtPi0Pt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Pt(),fWeightJetJetMC);
          }
          delete TruePi0;
          TruePi0=0x0;
        }
      }
    }
      break;
    // pcm-cal,pcm
    case 1:
    {
      // get gamma0MotherLabel
      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCStack);
        TParticle * gammaMC0 = (TParticle*)fMCStack->Particle(gamma0MCLabel);
        gamma0MotherLabel=gammaMC0->GetFirstMother();
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;
      TParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        // get mother or grandmother of gamma1 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC1 = (TParticle*)fMCStack->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=fMCStack->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // fill pi0 histograms here if necessary

        // get mother of gamma2
        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = TrueGammaCandidate2->GetMCParticleLabel(fMCStack);
          TParticle * gammaMC2 = (TParticle*)fMCStack->Particle(gamma2MCLabel);
          gamma2MotherLabel = gammaMC2->GetFirstMother();
        }

        // get pi0MotherLabel
        Int_t pi0MotherLabel = fMCStack->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCStack->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){ // also need to consider matching cluster of conversion positron to pcm tracks, not yet implemented
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother *TruePi0 = new AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(fDoMesonQA>0){
            fHistoTrueOmegaYPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTrueOmegaAlphaPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->GetAlpha(),fWeightJetJetMC);
            fHistoTrueOmegaEtaPhi[fiCut]->Fill(OmegaCandidate->Phi(),OmegaCandidate->Eta(),fWeightJetJetMC);
            fHistoTrueOmegaGammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtGammaPt[fiCut]->Fill(OmegaCandidate->Pt(),TrueGammaCandidate2->Pt(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaAlphaPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetAlpha(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaYPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaEtaPhi[fiCut]->Fill(TruePi0->Phi(),TruePi0->Eta(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaOpenAnglePt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetOpeningAngle(),fWeightJetJetMC);
            fHistoTrueOmegaPi0AnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TruePi0->Vect()),fWeightJetJetMC);
            fHistoTruePi0GammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtPi0Pt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Pt(),fWeightJetJetMC);
          }
          delete TruePi0;
          TruePi0=0x0;
        }
      }
    }
      break;
    // cal-cal,cal
    case 2:
    {
      // get gamma0MotherLabel
      if (!TrueGammaCandidate0->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma0MotherLabel = -1;
      TParticle * gammaMC0 = 0x0;
      if(gamma0MCLabel != -1){
        // get mother or grandmother of gamma0 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC0 = (TParticle*)fMCStack->Particle(gamma0MCLabel);
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother(0);
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion()) gamma0MotherLabel=fMCStack->Particle(gammaMC0->GetMother(0))->GetMother(0);
          else gamma0MotherLabel=gammaMC0->GetMother(0);
        }
      }

      // get gamma1MotherLabel
      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;
      TParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        // get mother or grandmother of gamma1 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC1 = (TParticle*)fMCStack->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=fMCStack->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // fill pi0 histograms here if necessary

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;
        TParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          // get mother or grandmother of gamma2 (potentially true omega) depending on whether it is an electron-leading/photon-leading cluster
          gammaMC2 = (TParticle*)fMCStack->Particle(gamma2MCLabel);
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){  // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother(0);
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()) gamma2MotherLabel=fMCStack->Particle(gammaMC2->GetMother(0))->GetMother(0);
            else gamma2MotherLabel=gammaMC2->GetMother(0);
          }
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCStack->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(isTruePi0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCStack->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother *TruePi0 = new AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(fDoMesonQA>0){
            fHistoTrueOmegaYPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTrueOmegaAlphaPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->GetAlpha(),fWeightJetJetMC);
            fHistoTrueOmegaEtaPhi[fiCut]->Fill(OmegaCandidate->Phi(),OmegaCandidate->Eta(),fWeightJetJetMC);
            fHistoTrueOmegaGammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtGammaPt[fiCut]->Fill(OmegaCandidate->Pt(),TrueGammaCandidate2->Pt(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaAlphaPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetAlpha(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaYPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaEtaPhi[fiCut]->Fill(TruePi0->Phi(),TruePi0->Eta(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaOpenAnglePt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetOpeningAngle(),fWeightJetJetMC);
            fHistoTrueOmegaPi0AnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TruePi0->Vect()),fWeightJetJetMC);
            fHistoTruePi0GammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtPi0Pt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Pt(),fWeightJetJetMC);
          }
          delete TruePi0;
          TruePi0=0x0;
        }
      }
    }
      break;
    // cal-cal,pcm
    case 3:
    {
      // get gamma0MotherLabel
      if (!TrueGammaCandidate0->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma0MotherLabel = -1;
      TParticle * gammaMC0 = 0x0;
      if(gamma0MCLabel != -1){
        // get mother or grandmother of gamma0 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC0 = (TParticle*)fMCStack->Particle(gamma0MCLabel);
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother(0);
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion()) gamma0MotherLabel=fMCStack->Particle(gammaMC0->GetMother(0))->GetMother(0);
          else gamma0MotherLabel=gammaMC0->GetMother(0);
        }
      }

      // get gamma1MotherLabel
      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;
      TParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        // get mother or grandmother of gamma1 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC1 = (TParticle*)fMCStack->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=fMCStack->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // fill pi0 histograms here if necessary

        // get mother of gamma2
        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = TrueGammaCandidate2->GetMCParticleLabel(fMCStack);
          TParticle * gammaMC2 = (TParticle*)fMCStack->Particle(gamma2MCLabel);
          gamma2MotherLabel = gammaMC2->GetFirstMother();
        }

        // get pi0MotherLabel
        Int_t pi0MotherLabel = fMCStack->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(isTruePi0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCStack->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){ // also need to consider matching cluster of conversion positron to pcm tracks, not yet implemented
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother *TruePi0 = new AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(fDoMesonQA>0){
            fHistoTrueOmegaYPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTrueOmegaAlphaPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->GetAlpha(),fWeightJetJetMC);
            fHistoTrueOmegaEtaPhi[fiCut]->Fill(OmegaCandidate->Phi(),OmegaCandidate->Eta(),fWeightJetJetMC);
            fHistoTrueOmegaGammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtGammaPt[fiCut]->Fill(OmegaCandidate->Pt(),TrueGammaCandidate2->Pt(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaAlphaPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetAlpha(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaYPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaEtaPhi[fiCut]->Fill(TruePi0->Phi(),TruePi0->Eta(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaOpenAnglePt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetOpeningAngle(),fWeightJetJetMC);
            fHistoTrueOmegaPi0AnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TruePi0->Vect()),fWeightJetJetMC);
            fHistoTruePi0GammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtPi0Pt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Pt(),fWeightJetJetMC);
          }
          delete TruePi0;
          TruePi0=0x0;
        }
      }
    }
      break;
    // pcm-pcm,cal
    case 4:
    {
      // get gamma0MotherLabel
      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCStack);
        TParticle * gammaMC0 = (TParticle*)fMCStack->Particle(gamma0MCLabel);
        gamma0MotherLabel=gammaMC0->GetFirstMother();
      }

      // get gamma1MotherLabel
      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(fMCStack);
        TParticle * gammaMC1 = (TParticle*)fMCStack->Particle(gamma1MCLabel);
        gamma1MotherLabel=gammaMC1->GetFirstMother();
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma1MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((TParticle*)fMCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // fill pi0 histograms here if necessary

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;
        TParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          // get mother or grandmother of gamma1 (potentially true omega) depending on whether it is an electron-leading/photon-leading cluster
          gammaMC2 = (TParticle*)fMCStack->Particle(gamma2MCLabel);
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){  // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother(0);
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()) gamma2MotherLabel=fMCStack->Particle(gammaMC2->GetMother(0))->GetMother(0);
            else gamma2MotherLabel=gammaMC2->GetMother(0);
          }
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCStack->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(isTruePi0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCStack->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){ // also need to consider matching cluster of conversion positron to pcm tracks, not yet implemented
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother *TruePi0 = new AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(fDoMesonQA>0){
            fHistoTrueOmegaYPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTrueOmegaAlphaPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->GetAlpha(),fWeightJetJetMC);
            fHistoTrueOmegaEtaPhi[fiCut]->Fill(OmegaCandidate->Phi(),OmegaCandidate->Eta(),fWeightJetJetMC);
            fHistoTrueOmegaGammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtGammaPt[fiCut]->Fill(OmegaCandidate->Pt(),TrueGammaCandidate2->Pt(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaAlphaPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetAlpha(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaYPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaEtaPhi[fiCut]->Fill(TruePi0->Phi(),TruePi0->Eta(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaOpenAnglePt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetOpeningAngle(),fWeightJetJetMC);
            fHistoTrueOmegaPi0AnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TruePi0->Vect()),fWeightJetJetMC);
            fHistoTruePi0GammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtPi0Pt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Pt(),fWeightJetJetMC);
          }
          delete TruePi0;
          TruePi0=0x0;
        }
      }
    }
      break;
    // pcm-pcm,pcm
    case 5:
    {
      // get gamma0MotherLabel
      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCStack);
        TParticle * gammaMC0 = (TParticle*)fMCStack->Particle(gamma0MCLabel);
        gamma0MotherLabel=gammaMC0->GetFirstMother();
      }

      // get gamma1MotherLabel
      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(fMCStack);
        TParticle * gammaMC1 = (TParticle*)fMCStack->Particle(gamma1MCLabel);
        gamma1MotherLabel=gammaMC1->GetFirstMother();
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma1MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((TParticle*)fMCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // fill pi0 histograms here if necessary

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = TrueGammaCandidate2->GetMCParticleLabel(fMCStack);
          TParticle * gammaMC2 = (TParticle*)fMCStack->Particle(gamma2MCLabel);
          gamma2MotherLabel=gammaMC2->GetFirstMother();
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCStack->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(isTruePi0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCStack->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){ // also need to consider matching cluster of conversion positron to pcm tracks, not yet implemented
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother *TruePi0 = new AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(fDoMesonQA>0){
            fHistoTrueOmegaYPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTrueOmegaAlphaPt[fiCut]->Fill(OmegaCandidate->Pt(),OmegaCandidate->GetAlpha(),fWeightJetJetMC);
            fHistoTrueOmegaEtaPhi[fiCut]->Fill(OmegaCandidate->Phi(),OmegaCandidate->Eta(),fWeightJetJetMC);
            fHistoTrueOmegaGammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtGammaPt[fiCut]->Fill(OmegaCandidate->Pt(),TrueGammaCandidate2->Pt(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaAlphaPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetAlpha(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaYPt[fiCut]->Fill(TruePi0->Pt(),TruePi0->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaEtaPhi[fiCut]->Fill(TruePi0->Phi(),TruePi0->Eta(),fWeightJetJetMC);
            fHistoTruePi0FromOmegaOpenAnglePt[fiCut]->Fill(TruePi0->Pt(),TruePi0->GetOpeningAngle(),fWeightJetJetMC);
            fHistoTrueOmegaPi0AnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TMath::Pi() - OmegaCandidate->Angle(TruePi0->Vect()),fWeightJetJetMC);
            fHistoTruePi0GammaAnglePt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Angle(TrueGammaCandidate2->Vect()),fWeightJetJetMC);
            fHistoTrueOmegaPtPi0Pt[fiCut]->Fill(OmegaCandidate->Pt(),TruePi0->Pt(),fWeightJetJetMC);
          }
          delete TruePi0;
          TruePi0=0x0;
        }
      }
    }
      break;
  }
}

//______________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *OmegaCandidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, AliAODConversionPhoton *TrueGammaCandidate2)
{
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  Bool_t isTruePi0 = kFALSE;
  Bool_t isTrueOmega = kFALSE;

  switch(fReconMethod){
    // pcm-cal,cal
    case 0:
    {
      AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positiveMC||!negativeMC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positiveMC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0->GetMother();

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;

        AliAODMCParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          gammaMC2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma2MCLabel));
          // get mother label
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){ // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother();
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()){
              AliAODMCParticle * gammaGrandMotherMC2 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC2->GetMother()));
              gamma2MotherLabel=gammaGrandMotherMC2->GetMother();
            }else gamma2MotherLabel=gammaMC2->GetMother();
          }
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)AODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
        }
      }
    }
      break;
    // pcm-cal,pcm
    case 1:
    {
      AliAODMCParticle *positive0MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negative0MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positive0MC||!negative0MC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positive0MC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0->GetMother();

        // get gamma2MotherLabel
        AliAODMCParticle *positive2MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelPositive()));
        AliAODMCParticle *negative2MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelNegative()));

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if(!positive2MC||!negative2MC)
          return;

        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = positive2MC->GetMother();
          AliAODMCParticle * gammaMC2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma2MCLabel));
          gamma2MotherLabel=gammaMC2->GetMother();
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)AODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
        }
      }
    }
      break;
    // cal-cal,cal
    case 2:
    {
      if (!TrueGammaCandidate0->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma0MotherLabel = -1;

      AliAODMCParticle * gammaMC0 = 0x0;
      if(gamma0MCLabel != -1){
        gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
        // get gamma0MotherLabel
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother();
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC0->GetMother()));
            gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
          }else gamma0MotherLabel=gammaMC0->GetMother();
        }
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0->GetMother();

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;

        AliAODMCParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          gammaMC2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma2MCLabel));
          // get gamma2MotherLabel
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){ // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother();
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()){
              AliAODMCParticle * gammaGrandMotherMC2 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC2->GetMother()));
              gamma2MotherLabel=gammaGrandMotherMC2->GetMother();
            }else gamma2MotherLabel=gammaMC2->GetMother();
          }
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)AODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
        }
      }
    }
      break;
    // cal-cal,pcm
    case 3:
    {
      if (!TrueGammaCandidate0->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma0MotherLabel = -1;

      AliAODMCParticle * gammaMC0 = 0x0;
      if(gamma0MCLabel != -1){
        gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
        // get gamma0MotherLabel
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother();
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC0->GetMother()));
            gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
          }else gamma0MotherLabel=gammaMC0->GetMother();
        }
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0->GetMother();

        // get gamma2MotherLabel
        AliAODMCParticle *positive2MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelPositive()));
        AliAODMCParticle *negative2MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelNegative()));

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if(!positive2MC||!negative2MC)
          return;

        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = positive2MC->GetMother();
          AliAODMCParticle * gammaMC2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma2MCLabel));
          gamma2MotherLabel=gammaMC2->GetMother();
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)AODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
        }
      }
    }
      break;
    // pcm-pcm,cal
    case 4:
    {
      AliAODMCParticle *positive0MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negative0MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positive0MC||!negative0MC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positive0MC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      AliAODMCParticle *positive1MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
      AliAODMCParticle *negative1MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));

      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if(!positive1MC||!negative1MC)
        return;

      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = positive1MC->GetMother();
        AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
        gamma1MotherLabel=gammaMC1->GetMother();
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0->GetMother();

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;

        AliAODMCParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          gammaMC2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma2MCLabel));
          // get mother label
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){ // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother();
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()){
              AliAODMCParticle * gammaGrandMotherMC2 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC2->GetMother()));
              gamma2MotherLabel=gammaGrandMotherMC2->GetMother();
            }else gamma2MotherLabel=gammaMC2->GetMother();
          }
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)AODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
        }
      }
    }
      break;
    // pcm-pcm,pcm
    case 5:
    {
      AliAODMCParticle *positive0MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negative0MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positive0MC||!negative0MC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positive0MC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      AliAODMCParticle *positive1MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
      AliAODMCParticle *negative1MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));

      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if(!positive1MC||!negative1MC)
        return;

      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = positive1MC->GetMother();
        AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
        gamma1MotherLabel=gammaMC1->GetMother();
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0->GetMother();

        AliAODMCParticle *positive2MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelPositive()));
        AliAODMCParticle *negative2MC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelNegative()));

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if(!positive2MC||!negative2MC)
          return;

        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = positive2MC->GetMother();
          AliAODMCParticle * gammaMC2 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma2MCLabel));
          gamma2MotherLabel=gammaMC2->GetMother();
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)AODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0->M(),TruePi0->Pt(),fWeightJetJetMC);
        }
      }
    }
      break;
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::CalculateBackground(){

  Int_t zbin = fBGClusHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;

  if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
    mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  }else {
    mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
  }

  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEvent1Vertex = NULL;
  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEvent2Vertex = NULL;

  switch (fReconMethod){
    //pcm-cal,cal
  case 0:
    // we consider first cases where the third cal photon is from the current event
    for(Int_t iCurrent3=0;iCurrent3<fClusterCandidates->GetEntries();iCurrent3++){
      AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent3));
      if(gamma2 == NULL || !(gamma2->GetIsCaloPhoton())) continue;
      AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma2->GetCaloClusterRef());

      // pcm from previous event 1, cal from previous event 2 (Diff1Diff2Same)
      for(Int_t previous1=0;previous1<fBGHandler[fiCut]->GetNBGEvents();previous1++){
        AliGammaConversionAODVector *previous1V0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous1);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous1);
        }
        for(UInt_t igamma0=0;igamma0<previous1V0s->size();igamma0++){
          AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previous1V0s->at(igamma0)));
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
          }
          for(Int_t previous2=previous1+1;previous2<fBGClusHandler[fiCut]->GetNBGEvents();previous2++){
            AliGammaConversionAODVector *previous2clusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous2);
            if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              bgEvent2Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous2);
            }
            for(UInt_t igamma1=0;igamma1<previous2clusters->size();igamma1++){
              AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previous2clusters->at(igamma1)));
              if(!gamma1.GetIsCaloPhoton()) continue;
              if(fMoveParticleAccordingToVertex == kTRUE && bgEvent2Vertex){
                MoveParticleAccordingToVertex(&gamma1,bgEvent2Vertex);
              }
              AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,&gamma1);
              if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                 || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
                 || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                delete BGpi0cand;
                BGpi0cand=0x0;
                continue;
              }
              AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
              fHistoDiff1Diff2SameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
              delete BGpi0cand;
              BGpi0cand=0x0;
              delete BGOmegacand;
              BGOmegacand=0x0;
            }
          }
        }
      }

      // pcm from previous event, cal from current event (DiffSameSame)
      for(Int_t iCurrent2=0;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
        if(iCurrent2==iCurrent3) continue;
        AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
        if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
        for(Int_t previous1=0;previous1<fBGHandler[fiCut]->GetNBGEvents();previous1++){
          AliGammaConversionAODVector *previous1V0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous1);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous1);
          }
          for(UInt_t igamma0=0;igamma0<previous1V0s->size();igamma0++){
            AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previous1V0s->at(igamma0)));
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoDiffSameSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }

      // pcm from current event, cal from previous event (SameDiffSame)
      for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
        AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent1));
        if(gamma0 == NULL) continue;
        Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC));
        if(matched) continue;
        for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
          AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
          }
          for(UInt_t igamma1=0;igamma1<previousclusters->size();igamma1++){
            AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previousclusters->at(igamma1)));
            if(!gamma1.GetIsCaloPhoton()) continue;
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma1,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,&gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoSameDiffSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }

    // gamma0 and gamma1 from the current event while gamma2 is taken from a bg event (SameSameDiff)
    for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
      AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
      if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
        bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
      }
      for(UInt_t igamma2=0;igamma2<previousclusters->size();igamma2++){
        AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousclusters->at(igamma2)));
        if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
          MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
        }
        for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent1));
          if(gamma0 == NULL) continue;
          for(Int_t iCurrent2=0;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
            if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
            AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0, cluster, fInputEvent, fWeightJetJetMC));
            if(matched) continue;
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,&gamma2);
            fHistoSameSameDiffBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }
    break;
  //pcm-cal,pcm
  case 1:
    // we consider first cases where the third pcm photon is from the current event
    for(Int_t iCurrent3=0;iCurrent3<fGammaCandidates->GetEntries();iCurrent3++){
      AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent3));
      if(gamma2 == NULL) continue;

      // pcm from previous event 1, cal from previous event 2 (Diff1Diff2Same)
      for(Int_t previous1=0;previous1<fBGHandler[fiCut]->GetNBGEvents();previous1++){
        AliGammaConversionAODVector *previous1V0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous1);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous1);
        }
        for(UInt_t igamma0=0;igamma0<previous1V0s->size();igamma0++){
          AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previous1V0s->at(igamma0)));
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
          }
          for(Int_t previous2=previous1+1;previous2<fBGClusHandler[fiCut]->GetNBGEvents();previous2++){
            AliGammaConversionAODVector *previous2clusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous2);
            if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              bgEvent2Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous2);
            }
            for(UInt_t igamma1=0;igamma1<previous2clusters->size();igamma1++){
              AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previous2clusters->at(igamma1)));
              if(!gamma1.GetIsCaloPhoton()) continue;
              if(fMoveParticleAccordingToVertex == kTRUE && bgEvent2Vertex){
                MoveParticleAccordingToVertex(&gamma1,bgEvent2Vertex);
              }
              AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,&gamma1);
              if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                 || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
                 || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                delete BGpi0cand;
                BGpi0cand=0x0;
                continue;
              }
              AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
              fHistoDiff1Diff2SameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
              delete BGpi0cand;
              BGpi0cand=0x0;
              delete BGOmegacand;
              BGOmegacand=0x0;
            }
          }
        }
      }

      // pcm from previous event, cal from current event (DiffSameSame)
      for(Int_t iCurrent2=0;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
        AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
        if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
        AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
        Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2,cluster, fInputEvent, fWeightJetJetMC));
        if(matched) continue;
        for(Int_t previous=0;previous<fBGHandler[fiCut]->GetNBGEvents();previous++){
          AliGammaConversionAODVector *previousV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
          }
          for(UInt_t igamma0=0;igamma0<previousV0s->size();igamma0++){
            AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previousV0s->at(igamma0)));
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoDiffSameSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }

      // pcm from current event, cal from previous event (SameDiffSame)
      for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
        if(iCurrent1==iCurrent3) continue;
        AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent1));
        if(gamma0 == NULL) continue;
        for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
          AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
          }
          for(UInt_t igamma1=0;igamma1<previousclusters->size();igamma1++){
            AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previousclusters->at(igamma1)));
            if(!gamma1.GetIsCaloPhoton()) continue;
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma1,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,&gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoSameDiffSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }

    // gamma0 and gamma1 from the current event while gamma2 is taken from a bg event (SameSameDiff)
    for(Int_t previous=0;previous<fBGHandler[fiCut]->GetNBGEvents();previous++){
      AliGammaConversionAODVector *previousV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
      if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
        bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
      }
      for(UInt_t igamma2=0;igamma2<previousV0s->size();igamma2++){
        AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousV0s->at(igamma2)));
        if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
          MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
        }
        for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent1));
          if(gamma0 == NULL) continue;
          for(Int_t iCurrent2=0;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
            if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
            AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0, cluster, fInputEvent, fWeightJetJetMC));
            if(matched) continue;
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,&gamma2);
            fHistoSameSameDiffBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }
    break;
  //cal-cal,cal
  case 2:
    // we consider first cases where the third cal photon is from the current event
    for(Int_t iCurrent3=0;iCurrent3<fClusterCandidates->GetEntries();iCurrent3++){
      AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent3));
      if(gamma2 == NULL || !(gamma2->GetIsCaloPhoton())) continue;

      // pcm from previous event 1, cal from previous event 2 (Diff1Diff2Same)
      for(Int_t previous1=0;previous1<fBGClusHandler[fiCut]->GetNBGEvents();previous1++){
        AliGammaConversionAODVector *previous1clusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous1);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous1);
        }
        for(UInt_t igamma0=0;igamma0<previous1clusters->size();igamma0++){
          AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previous1clusters->at(igamma0)));
          if(!gamma0.GetIsCaloPhoton()) continue;
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
          }
          for(Int_t previous2=previous1+1;previous2<fBGClusHandler[fiCut]->GetNBGEvents();previous2++){
            AliGammaConversionAODVector *previous2clusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous2);
            if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              bgEvent2Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous2);
            }
            for(UInt_t igamma1=0;igamma1<previous2clusters->size();igamma1++){
              AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previous2clusters->at(igamma1)));
              if(!gamma1.GetIsCaloPhoton()) continue;
              if(fMoveParticleAccordingToVertex == kTRUE && bgEvent2Vertex){
                MoveParticleAccordingToVertex(&gamma1,bgEvent2Vertex);
              }
              AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,&gamma1);
              if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                 || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
                 || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                delete BGpi0cand;
                BGpi0cand=0x0;
                continue;
              }
              AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
              fHistoDiff1Diff2SameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
              delete BGpi0cand;
              BGpi0cand=0x0;
              delete BGOmegacand;
              BGOmegacand=0x0;
            }
          }
        }
      }

      // gamma0 from previous event, gamma1 from current event, both cal (DiffSameSame)
      for(Int_t iCurrent2=0;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
        if(iCurrent2==iCurrent3) continue;
        AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
        if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
        for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
          AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
          }
          for(UInt_t igamma0=0;igamma0<previousclusters->size();igamma0++){
            AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previousclusters->at(igamma0)));
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoDiffSameSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }

    // gamma0 and gamma1 from the current event while gamma2 is taken from a bg event (SameSameDiff)
    for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
      AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
      if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
        bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
      }
      for(UInt_t igamma2=0;igamma2<previousclusters->size();igamma2++){
        AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousclusters->at(igamma2)));
        if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
          MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
        }
        for(Int_t iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent1));
          if(gamma0 == NULL || !(gamma0->GetIsCaloPhoton())) continue;
          for(Int_t iCurrent2=iCurrent1+1;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
            if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,&gamma2);
            fHistoSameSameDiffBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }
    break;
  //cal-cal,pcm
  case 3:
    // we consider first cases where the third pcm photon is from the current event
    for(Int_t iCurrent3=0;iCurrent3<fGammaCandidates->GetEntries();iCurrent3++){
      AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent3));
      if(gamma2 == NULL) continue;

      // gamma0 from previous event 1, gamma1 from previous event 2, both cal (Diff1Diff2Same)
      for(Int_t previous1=0;previous1<fBGClusHandler[fiCut]->GetNBGEvents();previous1++){
        AliGammaConversionAODVector *previous1clusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous1);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous1);
        }
        for(UInt_t igamma0=0;igamma0<previous1clusters->size();igamma0++){
          AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previous1clusters->at(igamma0)));
          if(!gamma0.GetIsCaloPhoton()) continue;
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
          }
          for(Int_t previous2=previous1+1;previous2<fBGClusHandler[fiCut]->GetNBGEvents();previous2++){
            AliGammaConversionAODVector *previous2clusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous2);
            if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              bgEvent2Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous2);
            }
            for(UInt_t igamma1=0;igamma1<previous2clusters->size();igamma1++){
              AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previous2clusters->at(igamma1)));
              if(!gamma1.GetIsCaloPhoton()) continue;
              if(fMoveParticleAccordingToVertex == kTRUE && bgEvent2Vertex){
                MoveParticleAccordingToVertex(&gamma1,bgEvent2Vertex);
              }
              AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,&gamma1);
              if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                 || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
                 || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                delete BGpi0cand;
                BGpi0cand=0x0;
                continue;
              }
              AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
              fHistoDiff1Diff2SameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
              delete BGpi0cand;
              BGpi0cand=0x0;
              delete BGOmegacand;
              BGOmegacand=0x0;
            }
          }
        }
      }

      // gamma0 from previous event, gamma1 from current event, both cal (DiffSameSame)
      for(Int_t iCurrent2=0;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
        AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
        if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
        AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
        Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2, cluster, fInputEvent, fWeightJetJetMC));
        if(matched) continue;
        for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
          AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
          }
          for(UInt_t igamma0=0;igamma0<previousclusters->size();igamma0++){
            AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previousclusters->at(igamma0)));
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoDiffSameSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }

    // gamma0 and gamma1 from the current event while gamma2 is taken from a bg event (SameSameDiff)
    for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
      AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
      if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
        bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
      }
      for(UInt_t igamma2=0;igamma2<previousclusters->size();igamma2++){
        AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousclusters->at(igamma2)));
        if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
          MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
        }
        for(Int_t iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent1));
          if(gamma0 == NULL || !(gamma0->GetIsCaloPhoton())) continue;
          for(Int_t iCurrent2=iCurrent1+1;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent2));
            if(gamma1 == NULL || !(gamma1->GetIsCaloPhoton())) continue;
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,&gamma2);
            fHistoSameSameDiffBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }
    break;
  //pcm-pcm,cal
  case 4:
    // we consider first cases where the third cal photon is from the current event
    for(Int_t iCurrent3=0;iCurrent3<fClusterCandidates->GetEntries();iCurrent3++){
      AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent3));
      if(gamma2 == NULL || !(gamma2->GetIsCaloPhoton())) continue;
      AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma2->GetCaloClusterRef());

      // gamma0 from previous event 1, gamma1 from previous event 2, both pcm (Diff1Diff2Same)
      for(Int_t previous1=0;previous1<fBGHandler[fiCut]->GetNBGEvents();previous1++){
        AliGammaConversionAODVector *previous1V0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous1);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous1);
        }
        for(UInt_t igamma0=0;igamma0<previous1V0s->size();igamma0++){
          AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previous1V0s->at(igamma0)));
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
          }
          for(Int_t previous2=previous1+1;previous2<fBGHandler[fiCut]->GetNBGEvents();previous2++){
            AliGammaConversionAODVector *previous2V0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous2);
            if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              bgEvent2Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous2);
            }
            for(UInt_t igamma1=0;igamma1<previous2V0s->size();igamma1++){
              AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previous2V0s->at(igamma1)));
              if(fMoveParticleAccordingToVertex == kTRUE && bgEvent2Vertex){
                MoveParticleAccordingToVertex(&gamma1,bgEvent2Vertex);
              }
              AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,&gamma1);
              if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                 || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
                 || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                delete BGpi0cand;
                BGpi0cand=0x0;
                continue;
              }
              AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
              fHistoDiff1Diff2SameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
              delete BGpi0cand;
              BGpi0cand=0x0;
              delete BGOmegacand;
              BGOmegacand=0x0;
            }
          }
        }
      }

      // gamma0 from previous event, gamma1 from current event, both pcm (DiffSameSame)
      for(Int_t iCurrent2=0;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
        AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent2));
        if(gamma1 == NULL) continue;
        Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma1, cluster, fInputEvent, fWeightJetJetMC));
        if(matched) continue;
        for(Int_t previous=0;previous<fBGHandler[fiCut]->GetNBGEvents();previous++){
          AliGammaConversionAODVector *previousV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
          }
          for(UInt_t igamma0=0;igamma0<previousV0s->size();igamma0++){
            AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previousV0s->at(igamma0)));
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoDiffSameSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }

    // gamma0 and gamma1 (both pcm) from the current event while gamma2 (cal) is taken from a bg event (SameSameDiff)
    for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
      AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
      if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
        bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
      }
      for(UInt_t igamma2=0;igamma2<previousclusters->size();igamma2++){
        AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousclusters->at(igamma2)));
        if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
          MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
        }
        for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent1));
          if(gamma0 == NULL) continue;
          for(Int_t iCurrent2=iCurrent1+1;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent2));
            if(gamma1 == NULL) continue;
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,&gamma2);
            fHistoSameSameDiffBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }
    break;
  //pcm-pcm,pcm
  case 5:
    // we consider first cases where the third pcm photon is from the current event
    for(Int_t iCurrent3=0;iCurrent3<fGammaCandidates->GetEntries();iCurrent3++){
      AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent3));
      if(gamma2 == NULL) continue;

      // gamma0 from previous event 1, gamma1 from previous event 2, both pcm (Diff1Diff2Same)
      for(Int_t previous1=0;previous1<fBGHandler[fiCut]->GetNBGEvents();previous1++){
        AliGammaConversionAODVector *previous1V0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous1);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous1);
        }
        for(UInt_t igamma0=0;igamma0<previous1V0s->size();igamma0++){
          AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previous1V0s->at(igamma0)));
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
          }
          for(Int_t previous2=previous1+1;previous2<fBGHandler[fiCut]->GetNBGEvents();previous2++){
            AliGammaConversionAODVector *previous2V0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous2);
            if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              bgEvent2Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous2);
            }
            for(UInt_t igamma1=0;igamma1<previous2V0s->size();igamma1++){
              AliAODConversionPhoton gamma1 = (AliAODConversionPhoton)(*(previous2V0s->at(igamma1)));
              if(fMoveParticleAccordingToVertex == kTRUE && bgEvent2Vertex){
                MoveParticleAccordingToVertex(&gamma1,bgEvent2Vertex);
              }
              AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,&gamma1);
              if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                 || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
                 || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
                delete BGpi0cand;
                BGpi0cand=0x0;
                continue;
              }
              AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
              fHistoDiff1Diff2SameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
              delete BGpi0cand;
              BGpi0cand=0x0;
              delete BGOmegacand;
              BGOmegacand=0x0;
            }
          }
        }
      }

      // gamma0 from previous event, gamma1 from current event, both pcm (DiffSameSame)
      for(Int_t iCurrent2=0;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
        if(iCurrent2==iCurrent3) continue;
        AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent2));
        if(gamma1 == NULL) continue;
        for(Int_t previous=0;previous<fBGHandler[fiCut]->GetNBGEvents();previous++){
          AliGammaConversionAODVector *previousV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
            bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
          }
          for(UInt_t igamma0=0;igamma0<previousV0s->size();igamma0++){
            AliAODConversionPhoton gamma0 = (AliAODConversionPhoton)(*(previousV0s->at(igamma0)));
            if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
              MoveParticleAccordingToVertex(&gamma0,bgEvent1Vertex);
            }
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(&gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,gamma2);
            fHistoDiffSameSameBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }

    // gamma0 and gamma1 (both pcm) from the current event while gamma2 (pcm) is taken from a bg event (SameSameDiff)
    for(Int_t previous=0;previous<fBGHandler[fiCut]->GetNBGEvents();previous++){
      AliGammaConversionAODVector *previousV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,previous);
      if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
        bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,previous);
      }
      for(UInt_t igamma2=0;igamma2<previousV0s->size();igamma2++){
        AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousV0s->at(igamma2)));
        if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
          MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
        }
        for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton *gamma0 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent1));
          if(gamma0 == NULL) continue;
          for(Int_t iCurrent2=iCurrent1+1;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton *gamma1 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent2));
            if(gamma1 == NULL) continue;
            AliAODConversionMother *BGpi0cand = new AliAODConversionMother(gamma0,gamma1);
            if(!(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
               ->MesonIsSelected(BGpi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
               || BGpi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionLow()
               || BGpi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->GetSelectionHigh()){
              delete BGpi0cand;
              BGpi0cand=0x0;
              continue;
            }
            AliAODConversionMother *BGOmegacand = new AliAODConversionMother(BGpi0cand,&gamma2);
            fHistoSameSameDiffBackInvMassPt[fiCut]->Fill(BGOmegacand->M(),BGOmegacand->Pt(),fWeightJetJetMC);
            delete BGpi0cand;
            BGpi0cand=0x0;
            delete BGOmegacand;
            BGOmegacand=0x0;
          }
        }
      }
    }
    break;
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
  //see header file for documentation
  
  Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
  Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
  Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();
  
  Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
  particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::UpdateEventByEventData(){
  //see header file for documentation
  if(fReconMethod==5){
    if(fGammaCandidates->GetEntries()>0){
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
      }else { // means we use #V0s for multiplicity
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
      }
    }
  } else if(fReconMethod==2){
    if(fClusterCandidates->GetEntries()>0){
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
      }else { // means we use #V0s for multiplicity
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
      }
    }
  } else{
    if(fClusterCandidates->GetEntries()>0 || fGammaCandidates->GetEntries()>0){
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
      }else { // means we use #V0s for multiplicity
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::RelabelAODPhotonCandidates(Bool_t mode){
  
  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  
  if(mode){
    fMCStackPos = new Int_t[fReaderGammas->GetEntries()];
    fMCStackNeg = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayPos = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayNeg = new Int_t[fReaderGammas->GetEntries()];
  }
  
  for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
    if(!PhotonCandidate) continue;
    if(!mode){// Back to ESD Labels
    PhotonCandidate->SetMCLabelPositive(fMCStackPos[iGamma]);
    PhotonCandidate->SetMCLabelNegative(fMCStackNeg[iGamma]);
    PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
    PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
    continue;
    }
    fMCStackPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
    fMCStackNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
    fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
    fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();
    
    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;
    
    for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
      AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if(!AODLabelPos){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
        }
      }
      if(!AODLabelNeg){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelNegative(i);
        AODLabelNeg = kTRUE;
        }
      }
      if(AODLabelNeg && AODLabelPos){
        break;
      }
    }
    if(!AODLabelPos || !AODLabelNeg){
      cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
      if(!AODLabelNeg){
        PhotonCandidate->SetMCLabelNegative(-999999);
        PhotonCandidate->SetLabelNegative(-999999);
      }
      if(!AODLabelPos){
        PhotonCandidate->SetMCLabelPositive(-999999);
        PhotonCandidate->SetLabelPositive(-999999);
      }
    }
  }
  
  
  if(!mode){
    delete[] fMCStackPos;
    delete[] fMCStackNeg;
    delete[] fESDArrayPos;
    delete[] fESDArrayNeg;
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
Int_t AliAnalysisTaskOmegaToPiZeroGamma::GetSourceClassification(Int_t daughter, Int_t pdgCode){
  
  if (daughter == 111) {
    if (abs(pdgCode) == 310) return 1; // k0s
    else if (abs(pdgCode) == 3122) return 2; // Lambda
    else if (abs(pdgCode) == 130) return 3; // K0L
    else if (abs(pdgCode) == 2212) return 4; // proton
    else if (abs(pdgCode) == 2112) return 5; // neutron
    else if (abs(pdgCode) == 211) return 6; // pion
    else if (abs(pdgCode) == 321) return 7; // kaon
    else if (abs(pdgCode) == 113 || abs(pdgCode) == 213 ) return 8; // rho 0,+,-
    else if (abs(pdgCode) == 3222 || abs(pdgCode) == 3212 || abs(pdgCode) == 3112  ) return 9; // Sigma
    else if (abs(pdgCode) == 2224 || abs(pdgCode) == 2214 || abs(pdgCode) == 2114 || abs(pdgCode) == 1114  ) return 10; // Delta
    else if (abs(pdgCode) == 313 || abs(pdgCode) == 323   ) return 11; // K*
    else return 15;
  }
  return 15;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskOmegaToPiZeroGamma::CheckVectorOnly(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else return false;
  }
  return false;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskOmegaToPiZeroGamma::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else{
      vec.push_back(tobechecked);
      return false;
    }
  }
  return false;
}

//_________________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessConversionPhotonsForMissingTags (){

  if (!fMCStack) return;
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
    if (gamma0->IsTrueConvertedPhoton()){
      Int_t gamma0MotherLabel = -1;
      Int_t gamma0MCLabel = gamma0->GetMCParticleLabel(fMCStack);
      if(gamma0MCLabel > -1){
        TParticle * gammaMC0 = (TParticle*)fMCStack->Particle(gamma0MCLabel);
        gamma0MotherLabel = gammaMC0->GetFirstMother();
        if (gamma0MotherLabel>-1){
          if(((TParticle*)fMCStack->Particle(gamma0MotherLabel))->GetPdgCode() == 111){
            if (!CheckVectorForDoubleCount(fVectorRecTruePi0s,gamma0MotherLabel)){
              Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              if (!isPrimary){
                Int_t secMotherLabel = ((TParticle*)fMCStack->Particle(gamma0MotherLabel))->GetMother(0);
                Float_t weightedSec= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCStack, fInputEvent) && fMCStack->Particle(secMotherLabel)->GetPdgCode()==310){
                  weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCStack, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
                }
              }else {
                Float_t weighted= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, fMCStack, fInputEvent)){
                  if (((TParticle*)fMCStack->Particle(gamma0MotherLabel))->Pt()>0.005){
                    weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, fMCStack, fInputEvent);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessConversionPhotonsForMissingTagsAOD (){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    
  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));

    if (gamma0->IsTrueConvertedPhoton()){
      AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0->GetMCLabelPositive()));
      AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0->GetMCLabelNegative()));
      
      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positiveMC||!negativeMC)
        return;
      
      if (gamma0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positiveMC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel = gammaMC0->GetMother();

        if (gamma0MotherLabel>-1){
          if(((AliAODMCParticle*)AODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 111){
            if (!CheckVectorForDoubleCount(fVectorRecTruePi0s,gamma0MotherLabel)){
              Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              if (!isPrimary){
                Int_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel))->GetMother();
                Float_t weightedSec= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
                  weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
                }
              }else {
                Float_t weighted= 1;
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, 0x0, fInputEvent)){
                  if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel))->Pt()>0.005){
                    weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, 0x0, fInputEvent);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
  if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
  else ma[tobechecked] = 2;
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
  map<Int_t, Int_t>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second, fWeightJetJetMC);
  }
  ma.clear();
  return;
}
