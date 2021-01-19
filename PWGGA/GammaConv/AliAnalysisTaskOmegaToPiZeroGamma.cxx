/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Daniel Mühlheim, Ziruo Zhang, Marvin Hemmer                    *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
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
#include "TGeoGlobalMagField.h"
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
#include "AliCaloTrackMatcher.h"
#include <vector>
#include <map>
#include <iostream>

ClassImp(AliAnalysisTaskOmegaToPiZeroGamma)

//________________________________________________________________________
AliAnalysisTaskOmegaToPiZeroGamma::AliAnalysisTaskOmegaToPiZeroGamma(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fBGHandler(NULL),
  fBGClusHandler(NULL),
  fBGPi0Handler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fClusterOutputList(NULL),
  fOutputContainer(NULL),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fPi0Candidates(NULL),
  fPi0CandidatesBackOri(NULL),
  fPi0CandidatesBackAlt(NULL),
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
  fHistoPhotonPairInvMassPt(NULL),
  fHistoPhotonPairMatchedInvMassPt(NULL),
  fHistoPhotonPairYPt(NULL),
  fHistoPhotonPairAlphaPt(NULL),
  fHistoPhotonPairOpenAnglePt(NULL),
  fHistoPhotonPairEtaPhi(NULL),
  fHistoMotherConvPhotonEtaPhi(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherAngleCutRejectedInvMassPt(NULL),
  fHistoMotherYPt(NULL),
  fHistoMotherAlphaPt(NULL),
  fHistoMotherEtaPhi(NULL),
  fHistoMotherPi0AnglePt(NULL),
  fHistoMotherGammaAnglePt(NULL),
  fHistoPi0GammaAnglePt(NULL),
  fHistoMotherRestGammaCosAnglePt(NULL),
  fHistoMotherRestPi0CosAnglePt(NULL),
  fHistoMotherDalitzPlot(NULL),
  fHistoMotherDalitzPlotwoCut(NULL),
  fHistoGammaFromMotherPt(NULL),
  fHistoRecoArmenterosPodolanskiPlot(NULL),
  fHistoRecoArmenterosPodolanskiPlotwoCut(NULL),
  fHistoDiffPi0SameGammaBackInvMassPt(NULL),
  fHistoSamePi0DiffGammaBackInvMassPt(NULL),
  fHistoMotherSwappingBackInvMassPt(NULL),
  fHistoPi0SwappingBackInvMassPt(NULL),
  fHistoBackDalitzPlot(NULL),
  fHistoBackDalitzPlotwoCut(NULL),
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
  fHistoMCOmegaInAccInvMassPt(NULL),
  fHistoMCOmegaInvMassPt(NULL),
  fHistoMCAllOmegaYPt(NULL),
  fHistoMCOmegaInAccYPt(NULL),
  fHistoMCAllOmegaAlphaPt(NULL),
  fHistoMCOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromAllOmegaAlphaPt(NULL),
  fHistoMCPi0FromOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromAllOmegaYPt(NULL),
  fHistoMCPi0FromOmegaInAccYPt(NULL),
  fHistoMCPi0FromOmegaInAccInvMassPt(NULL),
  fHistoMCPi0FromAllOmegaEtaPhi(NULL),
  fHistoMCPi0FromOmegaInAccEtaPhi(NULL),
  fHistoMCAllOmegaEtaPhi(NULL),
  fHistoMCOmegaInAccEtaPhi(NULL),
  fHistoMCAllOmegaPiZeroAnglePt(NULL),
  fHistoMCAllPiZeroGammaAnglePt(NULL),
  fHistoMCAllOmegaGammaAnglePt(NULL),
  fHistoMCInAccOmegaPiZeroAnglePt(NULL),
  fHistoMCInAccPiZeroGammaAnglePt(NULL),
  fHistoMCInAccOmegaGammaAnglePt(NULL),
  fHistoMCAllOmegaInvMassPt(NULL),
  fHistoMCAllOmegaPtPi0Pt(NULL),
  fHistoMCInAccOmegaPtPi0Pt(NULL),
  fHistoMCAllOmegaPtGammaPt(NULL),
  fHistoMCInAccOmegaPtGammaPt(NULL),
  fHistoTrueOmegaInvMassPt(NULL),
  fHistoTrueOmegaYPt(NULL),
  fHistoTrueOmegaAlphaPt(NULL),
  fHistoTruePi0FromOmegaYPt(NULL),
  fHistoTruePi0FromOmegaInvMassPt(NULL),
  fHistoTruePi0FromOmegaAlphaPt(NULL),
  fHistoTruePi0FromOmegaEtaPhi(NULL),
  fHistoTruePi0FromOmegaOpenAnglePt(NULL),
  fHistoTrueOmegaPi0AnglePt(NULL),
  fHistoTrueOmegaGammaAnglePt(NULL),
  fHistoTruePi0GammaAnglePt(NULL),
  fHistoTrueOmegaEtaPhi(NULL),
  fHistoTrueOmegaPtPi0Pt(NULL),
  fHistoTrueOmegaRestGammaCosAnglePt(NULL),
  fHistoTrueOmegaRestPi0CosAnglePt(NULL),
  fHistoTrueDalitzPlot(NULL),
  fHistoTrueOmegaPtGammaPt(NULL),
  fHistoTrueGammaFromOmegaPt(NULL),
  fHistoTrueArmenterosPodolanskiPlot(NULL),
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
  fNGammaCandidates(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMCEventPos(NULL),
  fMCEventNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fnCuts(0),
  fiCut(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fSetPlotHistsExtQA(kFALSE),
  fUseDalitzCut(kFALSE),
  fWeightJetJetMC(1),
  fEnableSortForClusMC(kFALSE),
  fReconMethod(0),
  flowerFactor(0),
  fupperFactor(0),
  fMinPi0Pt(0),
  fmaxfit(NULL),
  fDoPiZeroGammaAngleCut(kFALSE),
  fTrackMatcherRunningMode(0),
  fRandom(0),
  fGenPhaseSpace(),
  fPhotonSelectionMode(0),
  dropOutGammas_CALO{},
  dropOutGammas_CALOBack{},
  dropOutGammas_PCM{},
  fAODMCTrackArray(NULL)
{

}

//________________________________________________________________________
AliAnalysisTaskOmegaToPiZeroGamma::AliAnalysisTaskOmegaToPiZeroGamma(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fBGHandler(NULL),
  fBGClusHandler(NULL),
  fBGPi0Handler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fClusterOutputList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fPi0Candidates(NULL),
  fPi0CandidatesBackOri(NULL),
  fPi0CandidatesBackAlt(NULL),
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
  fHistoPhotonPairInvMassPt(NULL),
  fHistoPhotonPairMatchedInvMassPt(NULL),
  fHistoPhotonPairYPt(NULL),
  fHistoPhotonPairAlphaPt(NULL),
  fHistoPhotonPairOpenAnglePt(NULL),
  fHistoPhotonPairEtaPhi(NULL),
  fHistoMotherConvPhotonEtaPhi(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherAngleCutRejectedInvMassPt(NULL),
  fHistoMotherYPt(NULL),
  fHistoMotherAlphaPt(NULL),
  fHistoMotherEtaPhi(NULL),
  fHistoMotherPi0AnglePt(NULL),
  fHistoMotherGammaAnglePt(NULL),
  fHistoPi0GammaAnglePt(NULL),
  fHistoMotherRestGammaCosAnglePt(NULL),
  fHistoMotherRestPi0CosAnglePt(NULL),
  fHistoMotherDalitzPlot(NULL),
  fHistoMotherDalitzPlotwoCut(NULL),
  fHistoGammaFromMotherPt(NULL),
  fHistoRecoArmenterosPodolanskiPlot(NULL),
  fHistoRecoArmenterosPodolanskiPlotwoCut(NULL),
  fHistoDiffPi0SameGammaBackInvMassPt(NULL),
  fHistoSamePi0DiffGammaBackInvMassPt(NULL),
  fHistoMotherSwappingBackInvMassPt(NULL),
  fHistoPi0SwappingBackInvMassPt(NULL),
  fHistoBackDalitzPlot(NULL),
  fHistoBackDalitzPlotwoCut(NULL),
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
  fHistoMCOmegaInAccInvMassPt(NULL),
  fHistoMCOmegaInvMassPt(NULL),
  fHistoMCAllOmegaYPt(NULL),
  fHistoMCOmegaInAccYPt(NULL),
  fHistoMCAllOmegaAlphaPt(NULL),
  fHistoMCOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromAllOmegaAlphaPt(NULL),
  fHistoMCPi0FromOmegaInAccAlphaPt(NULL),
  fHistoMCPi0FromAllOmegaYPt(NULL),
  fHistoMCPi0FromOmegaInAccYPt(NULL),
  fHistoMCPi0FromOmegaInAccInvMassPt(NULL),
  fHistoMCPi0FromAllOmegaEtaPhi(NULL),
  fHistoMCPi0FromOmegaInAccEtaPhi(NULL),
  fHistoMCAllOmegaEtaPhi(NULL),
  fHistoMCOmegaInAccEtaPhi(NULL),
  fHistoMCAllOmegaPiZeroAnglePt(NULL),
  fHistoMCAllPiZeroGammaAnglePt(NULL),
  fHistoMCAllOmegaGammaAnglePt(NULL),
  fHistoMCInAccOmegaPiZeroAnglePt(NULL),
  fHistoMCInAccPiZeroGammaAnglePt(NULL),
  fHistoMCInAccOmegaGammaAnglePt(NULL),
  fHistoMCAllOmegaInvMassPt(NULL),
  fHistoMCAllOmegaPtPi0Pt(NULL),
  fHistoMCInAccOmegaPtPi0Pt(NULL),
  fHistoMCAllOmegaPtGammaPt(NULL),
  fHistoMCInAccOmegaPtGammaPt(NULL),
  fHistoTrueOmegaInvMassPt(NULL),
  fHistoTrueOmegaYPt(NULL),
  fHistoTrueOmegaAlphaPt(NULL),
  fHistoTruePi0FromOmegaYPt(NULL),
  fHistoTruePi0FromOmegaInvMassPt(NULL),
  fHistoTruePi0FromOmegaAlphaPt(NULL),
  fHistoTruePi0FromOmegaEtaPhi(NULL),
  fHistoTruePi0FromOmegaOpenAnglePt(NULL),
  fHistoTrueOmegaPi0AnglePt(NULL),
  fHistoTrueOmegaGammaAnglePt(NULL),
  fHistoTruePi0GammaAnglePt(NULL),
  fHistoTrueOmegaEtaPhi(NULL),
  fHistoTrueOmegaPtPi0Pt(NULL),
  fHistoTrueOmegaRestGammaCosAnglePt(NULL),
  fHistoTrueOmegaRestPi0CosAnglePt(NULL),
  fHistoTrueDalitzPlot(NULL),
  fHistoTrueOmegaPtGammaPt(NULL),
  fHistoTrueGammaFromOmegaPt(NULL),
  fHistoTrueArmenterosPodolanskiPlot(NULL),
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
  fNGammaCandidates(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMCEventPos(NULL),
  fMCEventNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fnCuts(0),
  fiCut(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fIsHeavyIon(0),
  fDoLightOutput(kFALSE),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsMC(0),
  fSetPlotHistsExtQA(kFALSE),
  fUseDalitzCut(kFALSE),
  fWeightJetJetMC(1),
  fEnableSortForClusMC(kFALSE),
  fReconMethod(0),
  flowerFactor(0),
  fupperFactor(0),
  fMinPi0Pt(0),
  fmaxfit(NULL),
  fDoPiZeroGammaAngleCut(kFALSE),
  fTrackMatcherRunningMode(0),
  fRandom(0),
  fGenPhaseSpace(),
  fPhotonSelectionMode(0),
  dropOutGammas_CALO{},
  dropOutGammas_CALOBack{},
  dropOutGammas_PCM{},
  fAODMCTrackArray(NULL)
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
  if(fPi0Candidates){
    delete fPi0Candidates;
    fPi0Candidates = 0x0;
  }
  if(fPi0CandidatesBackOri){
    delete fPi0CandidatesBackOri;
    fPi0CandidatesBackOri = 0x0;
  }
  if(fPi0CandidatesBackAlt){
    delete fPi0CandidatesBackAlt;
    fPi0CandidatesBackAlt = 0x0;
  }
  if(fBGHandler){
    delete[] fBGHandler;
    fBGHandler = 0x0;
  }
  if(fBGClusHandler){
    delete[] fBGClusHandler;
    fBGClusHandler = 0x0;
  }
  if(fBGPi0Handler){
    delete[] fBGPi0Handler;
    fBGPi0Handler = 0x0;
  }
  if(fDoPiZeroGammaAngleCut){
    delete fmaxfit;
    fmaxfit = 0x0;
  }
}
//___________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::InitBack(){

  fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGClusHandler = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGPi0Handler = new AliGammaConversionAODBGHandler*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConversionMesonCuts*)fNeutralPionCutArray->At(iCut))->DoBGCalculation()){

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
      fBGPi0Handler[iCut] = new AliGammaConversionAODBGHandler(
                                collisionSystem,centMin,centMax,
                                ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                2,8,5);
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::UserCreateOutputObjects(){

  if (fIsMC > 1){
    fDoPhotonQA       = 0;
  }

  // create function for pi0-gamma angle cut
  if(fDoPiZeroGammaAngleCut){
    fmaxfit = new TF1("maxfit", "4.99209 / pow(x + 1.34075, 1.65) + 0.0568024");
  }

  if(fReconMethod > 3) fMinPi0Pt = 0.5; //PCM-PCM
  else if(fReconMethod > 1) fMinPi0Pt = 1.5; //EMCAL-EMCAL
  else fMinPi0Pt = 1.0; //PCM-EMCAL

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
  fPi0Candidates      = new TList();
  fPi0Candidates->SetOwner(kTRUE);
  fPi0CandidatesBackOri = new TList();
  fPi0CandidatesBackOri->SetOwner(kTRUE);
  fPi0CandidatesBackAlt = new TList();
  fPi0CandidatesBackAlt->SetOwner(kTRUE);

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
  if(fDoPiZeroGammaAngleCut)
    fHistoMotherAngleCutRejectedInvMassPt = new TH2F*[fnCuts];

  // BG histograms
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {
    fHistoMotherSwappingBackInvMassPt     = new TH2F*[fnCuts];
    if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg() ) == 2 )
    {
      fHistoPi0SwappingBackInvMassPt        = new TH2F*[fnCuts];
    }
  }
  else {
    fHistoDiffPi0SameGammaBackInvMassPt   = new TH2F*[fnCuts];
    fHistoSamePi0DiffGammaBackInvMassPt   = new TH2F*[fnCuts];
  }


  // QA histograms
  if(!fDoLightOutput){
    if(fDoMesonQA & 0b00000001) {
      fHistoPhotonPairYPt                 = new TH2F*[fnCuts];
      fHistoMotherYPt                     = new TH2F*[fnCuts];
    }
    if(fDoMesonQA & 0b00000010) {
      fHistoPhotonPairAlphaPt             = new TH2F*[fnCuts];
      fHistoMotherAlphaPt                 = new TH2F*[fnCuts];
    }
    if(fDoMesonQA & 0b00000100) {
      fHistoPhotonPairOpenAnglePt         = new TH2F*[fnCuts];
    }
    if(fDoMesonQA & 0b00001000) {
      fHistoPhotonPairEtaPhi              = new TH2F*[fnCuts];
      fHistoMotherEtaPhi                  = new TH2F*[fnCuts];
      if(fReconMethod!=2 && fReconMethod!=5)
      fHistoMotherConvPhotonEtaPhi      = new TH2F*[fnCuts];
    }
    if(fDoMesonQA & 0b00010000) {
      fHistoMotherPi0AnglePt              = new TH2F*[fnCuts];
      fHistoMotherGammaAnglePt            = new TH2F*[fnCuts];
      fHistoPi0GammaAnglePt               = new TH2F*[fnCuts];
    }
    if(fDoMesonQA & 0b00100000) {
      fHistoMotherRestGammaCosAnglePt     = new TH2F*[fnCuts];
      fHistoMotherRestPi0CosAnglePt       = new TH2F*[fnCuts];
    }

    if(fDoMesonQA & 0b01000000) {
      fHistoRecoArmenterosPodolanskiPlot      = new TH2F*[fnCuts];
      fHistoRecoArmenterosPodolanskiPlotwoCut = new TH2F*[fnCuts];
    }
    if(fDoMesonQA & 0b10000000) {
      fHistoMotherDalitzPlot              = new TH2F*[fnCuts];
      fHistoMotherDalitzPlotwoCut         = new TH2F*[fnCuts];
      fHistoBackDalitzPlot                = new TH2F*[fnCuts];
      fHistoBackDalitzPlotwoCut           = new TH2F*[fnCuts];
    }
  }

  if(fReconMethod!=5){
    fClusterOutputList                  = new TList*[fnCuts];
    fHistoClusGammaPt                   = new TH1F*[fnCuts];
    fHistoClusOverlapHeadersGammaPt     = new TH1F*[fnCuts];
  }

  Int_t nBinsPt = 35;
  Float_t maxPt = 50;
  Double_t *arrPtBinning      = new Double_t[100];
  for(Int_t i=0; i<nBinsPt+1;i++){
    if (i < 20) arrPtBinning[i]             = 1.0*i;                            // 00.0 - 20.0 in 1.0 steps
    else if(i<35) arrPtBinning[i]           = 20.0+2.0*(i-20);                  // 20.0 - 50.0 in 2.0 steps
    else  arrPtBinning[i]                   = maxPt;
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

    fHistoNEvents[iCut]       = new TH1F("NEvents","NEvents",14,-0.5,13.5);
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
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    fHistoNEventsMinGamma[iCut]     = new TH1F("NEventsMinGamma","NEventsMinGamma",4,-0.5,3.5);
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(1,"2 EMCal 1 PCM");
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(2,"2 PCM 1 EMCal");
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(3,"3 EMCal");
    fHistoNEventsMinGamma[iCut]->GetXaxis()->SetBinLabel(4,"3 PCM");
    fESDList[iCut]->Add(fHistoNEventsMinGamma[iCut]);

    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut]   = new TH1F("NEventsWOWeight","NEventsWOWeight",14,-0.5,13.5);
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
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
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

      fHistoClusGammaPt[iCut]         = new TH1F("ClusGamma_Pt","ClusGamma_Pt",nBinsPt, arrPtBinning);
      fHistoClusGammaPt[iCut]->SetXTitle("p_{T,clus}(GeV/c)");
      fClusterOutputList[iCut]->Add(fHistoClusGammaPt[iCut]);
      fHistoClusOverlapHeadersGammaPt[iCut]   = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",nBinsPt, arrPtBinning);
      fHistoClusOverlapHeadersGammaPt[iCut]->SetXTitle("p_{T,clus}(GeV/c), selected header w/ overlap");
      fClusterOutputList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);

      if (fIsMC > 1){
        fHistoClusGammaPt[iCut]->Sumw2();
        fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      }
    }

    fHistoPhotonPairInvMassPt[iCut]             = new TH2F("ESD_PhotonPair_InvMass_Pt","ESD_PhotonPair_InvMass_Pt",300,0,0.3,200,0.,20.);
    fHistoPhotonPairInvMassPt[iCut]->SetXTitle("M_{inv, #pi^{0} cand}(GeV/c^{2})");
    fHistoPhotonPairInvMassPt[iCut]->SetYTitle("p_{T, #pi^{0} cand}(GeV/c)");
    fESDList[iCut]->Add(fHistoPhotonPairInvMassPt[iCut]);

    fHistoMotherInvMassPt[iCut]                 = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
    fHistoMotherInvMassPt[iCut]->SetXTitle("M_{inv, #omega cand}(GeV/c^{2})");
    fHistoMotherInvMassPt[iCut]->SetYTitle("p_{T, #omega cand}(GeV/c)");
    fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);

    fHistoGammaFromMotherPt[iCut] = new TH1F("ESD_GammaFromMother_Pt","ESD_GammaFromMother_Pt",nBinsPt, arrPtBinning);
    fHistoGammaFromMotherPt[iCut]->SetXTitle("p_{T, #gamma from #omega cand}(GeV/c)");
    fESDList[iCut]->Add(fHistoGammaFromMotherPt[iCut]);

    if(fReconMethod!=2 && fReconMethod!=5){
      fHistoMotherMatchedInvMassPt[iCut]        = new TH2F("ESD_Mother_Matched_InvMass_Pt","ESD_Mother_Matched_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoMotherMatchedInvMassPt[iCut]->SetXTitle("M_{inv, #omega matched}(GeV/c^{2})");
      fHistoMotherMatchedInvMassPt[iCut]->SetYTitle("p_{T, #omega matched}(GeV/c)");
      fESDList[iCut]->Add(fHistoMotherMatchedInvMassPt[iCut]);
    }

    if(fDoPiZeroGammaAngleCut){
      fHistoMotherAngleCutRejectedInvMassPt[iCut]        = new TH2F("ESD_Mother_AngleCutRejected_InvMass_Pt","ESD_Mother_AngleCutRejected_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoMotherAngleCutRejectedInvMassPt[iCut]->SetXTitle("M_{inv, #omega rejected}(GeV/c^{2})");
      fHistoMotherAngleCutRejectedInvMassPt[iCut]->SetYTitle("p_{T, #omega rejected}(GeV/c)");
      fESDList[iCut]->Add(fHistoMotherAngleCutRejectedInvMassPt[iCut]);
    }

    if(fReconMethod<2){
      fHistoPhotonPairMatchedInvMassPt[iCut]      = new TH2F("ESD_PhotonPair_Matched_InvMass_Pt","ESD_PhotonPair_Matched_InvMass_Pt",300,0,0.3,200,0.,20.);
      fHistoPhotonPairMatchedInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2}) matched conv e^{+/-}to cluster");
      fHistoPhotonPairMatchedInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      fESDList[iCut]->Add(fHistoPhotonPairMatchedInvMassPt[iCut]);
    }

    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {
      fHistoMotherSwappingBackInvMassPt[iCut]     = new TH2F("ESD_Mother_SwappingBack_InvMass_Pt","ESD_Mother_SwappingBack_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoMotherSwappingBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
      fHistoMotherSwappingBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      fESDList[iCut]->Add(fHistoMotherSwappingBackInvMassPt[iCut]);

      if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg() ) == 2 )
      {
        fHistoPi0SwappingBackInvMassPt[iCut]     = new TH2F("ESD_Pi0_SwappingBack_InvMass_Pt","ESD_Pi0_SwappingBack_InvMass_Pt", 300, 0.0, 0.3, 200, 0., 20.);
        fHistoPi0SwappingBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
        fHistoPi0SwappingBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
        fESDList[iCut]->Add(fHistoPi0SwappingBackInvMassPt[iCut]);
      }
    }

    else {
      fHistoDiffPi0SameGammaBackInvMassPt[iCut]     = new TH2F("ESD_Mother_DiffPi0SameGamma_InvMass_Pt","ESD_Mother_DiffPi0SameGamma_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoDiffPi0SameGammaBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
      fHistoDiffPi0SameGammaBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      fESDList[iCut]->Add(fHistoDiffPi0SameGammaBackInvMassPt[iCut]);

      fHistoSamePi0DiffGammaBackInvMassPt[iCut]     = new TH2F("ESD_Mother_SamePi0DiffGamma_InvMass_Pt","ESD_Mother_SamePi0DiffGamma_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoSamePi0DiffGammaBackInvMassPt[iCut]->SetXTitle("M_{inv}(GeV/c^{2})");
      fHistoSamePi0DiffGammaBackInvMassPt[iCut]->SetYTitle("p_{T}(GeV/c)");
      fESDList[iCut]->Add(fHistoSamePi0DiffGammaBackInvMassPt[iCut]);
    }


    if (fIsMC > 1){
      fHistoPhotonPairInvMassPt[iCut]->Sumw2();
      fHistoMotherInvMassPt[iCut]->Sumw2();
      fHistoGammaFromMotherPt[iCut]->Sumw2();
      if(fReconMethod!=2 && fReconMethod!=5) fHistoMotherMatchedInvMassPt[iCut]->Sumw2();
      if(fDoPiZeroGammaAngleCut) fHistoMotherAngleCutRejectedInvMassPt[iCut]->Sumw2();
      if(fReconMethod<2) fHistoPhotonPairMatchedInvMassPt[iCut]->Sumw2();
    }
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {
      fHistoMotherSwappingBackInvMassPt[iCut]->Sumw2();
      if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg() ) == 2 )
      {
        fHistoPi0SwappingBackInvMassPt[iCut]->Sumw2();
      }
    }
    else {
      fHistoDiffPi0SameGammaBackInvMassPt[iCut]->Sumw2();
      fHistoSamePi0DiffGammaBackInvMassPt[iCut]->Sumw2();
    }

    if(!fDoLightOutput){
      if(fDoMesonQA & 0b00000001) {
        fHistoPhotonPairYPt[iCut]              = new TH2F("ESD_PhotonPair_Y_Pt","ESD_PhotonPair_Y_Pt",300,0.03,30.,150,-1.5,1.5);
        fHistoPhotonPairYPt[iCut]->SetXTitle("p_{T, #pi^{0}cand}(GeV/c)");
        fHistoPhotonPairYPt[iCut]->SetYTitle("y_{#pi^{0}cand}");
        SetLogBinningXTH2(fHistoPhotonPairYPt[iCut]);
        fESDList[iCut]->Add(fHistoPhotonPairYPt[iCut]);

        fHistoMotherYPt[iCut]                       = new TH2F("ESD_Mother_Y_Pt","ESD_Mother_Y_Pt",nBinsPt, arrPtBinning, 150,-1.5,1.5);
        fHistoMotherYPt[iCut]->SetYTitle("y_{#omega cand}");
        fHistoMotherYPt[iCut]->SetXTitle("p_{T, #omega cand}(GeV/c)");
        fESDList[iCut]->Add(fHistoMotherYPt[iCut]);
      }
      if(fDoMesonQA & 0b00000010) {
        fHistoPhotonPairAlphaPt[iCut]          = new TH2F("ESD_PhotonPair_Alpha_Pt","ESD_PhotonPair_Alpha_Pt",300,0.03,30.,200,-1,1);
        fHistoPhotonPairAlphaPt[iCut]->SetXTitle("p_{T, #pi^{0}cand}(GeV/c)");
        fHistoPhotonPairAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}cand}");
        SetLogBinningXTH2(fHistoPhotonPairAlphaPt[iCut]);
        fESDList[iCut]->Add(fHistoPhotonPairAlphaPt[iCut]);

        fHistoMotherAlphaPt[iCut]                   = new TH2F("ESD_Mother_Alpha_Pt","ESD_Mother_Alpha_Pt",nBinsPt, arrPtBinning,200,-1,1);
        fHistoMotherAlphaPt[iCut]->SetXTitle("p_{T, #omega cand}(GeV/c)");
        fHistoMotherAlphaPt[iCut]->SetYTitle("#alpha_{#omega cand}");
        fESDList[iCut]->Add(fHistoMotherAlphaPt[iCut]);
      }
      if(fDoMesonQA & 0b00000100) {
        fHistoPhotonPairOpenAnglePt[iCut]      = new TH2F("ESD_PhotonPair_OpenAngle_Pt","ESD_PhotonPair_OpenAngle_Pt",300,0.03,30.,100,0,1);
        fHistoPhotonPairOpenAnglePt[iCut]->SetXTitle("p_{T, #pi^{0}cand}(GeV/c)");
        fHistoPhotonPairOpenAnglePt[iCut]->SetYTitle("#theta_{#pi^{0}cand}");
        SetLogBinningXTH2(fHistoPhotonPairOpenAnglePt[iCut]);
        fESDList[iCut]->Add(fHistoPhotonPairOpenAnglePt[iCut]);
      }
      if(fDoMesonQA & 0b00001000) {
        fHistoPhotonPairEtaPhi[iCut] = new TH2F("ESD_PhotonPair_Eta_Phi","ESD_PhotonPair_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoPhotonPairEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}cand}(rad)");
        fHistoPhotonPairEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}cand}");
        fESDList[iCut]->Add(fHistoPhotonPairEtaPhi[iCut]);

        if(fReconMethod!=2 && fReconMethod!=5){
          fHistoMotherConvPhotonEtaPhi[iCut] = new TH2F("ESD_MotherConvPhoton_Eta_Phi","ConvPhoton under #omega peak",600,0,2*TMath::Pi(),200,-1,1);
          fHistoMotherConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}}(rad)");
          fHistoMotherConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
          fESDList[iCut]->Add(fHistoMotherConvPhotonEtaPhi[iCut]);
        }

        fHistoMotherEtaPhi[iCut]                  = new TH2F("ESD_Mother_Eta_Phi","ESD_Mother_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
        fHistoMotherEtaPhi[iCut]->SetXTitle("#phi_{#omega cand}(rad)");
        fHistoMotherEtaPhi[iCut]->SetYTitle("#eta_{#omega cand}");
        fESDList[iCut]->Add(fHistoMotherEtaPhi[iCut]);
      }
      if(fDoMesonQA & 0b00010000) {
        fHistoMotherPi0AnglePt[iCut] = new TH2F("ESD_MotherPi0_Angle_Pt","ESD_MotherPi0_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
        fHistoMotherPi0AnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
        fHistoMotherPi0AnglePt[iCut]->SetYTitle("#theta_{#omega cand, pi^{0}cand}");
        fESDList[iCut]->Add(fHistoMotherPi0AnglePt[iCut]);

        fHistoMotherGammaAnglePt[iCut] = new TH2F("ESD_MotherGamma_Angle_Pt","ESD_MotherGamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
        fHistoMotherGammaAnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
        fHistoMotherGammaAnglePt[iCut]->SetYTitle("#theta_{#omega cand, #gamma}");
        fESDList[iCut]->Add(fHistoMotherGammaAnglePt[iCut]);

        fHistoPi0GammaAnglePt[iCut] = new TH2F("ESD_Pi0Gamma_Angle_Pt","ESD_Pi0Gamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
        fHistoPi0GammaAnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
        fHistoPi0GammaAnglePt[iCut]->SetYTitle("#theta_{#pi^{0}cand,#gamma}");
        fESDList[iCut]->Add(fHistoPi0GammaAnglePt[iCut]);
      }
      if(fDoMesonQA & 0b00100000) {
        fHistoMotherRestGammaCosAnglePt[iCut] = new TH2F("ESD_MotherRestGamma_CosAngle_Pt","ESD_MotherRestGamma_CosAngle_Pt",nBinsPt, arrPtBinning,200,-1,1);
        fHistoMotherRestGammaCosAnglePt[iCut]->SetXTitle("#it{p}_{T, #omega} (GeV/#it{c})");
        fHistoMotherRestGammaCosAnglePt[iCut]->SetYTitle("cos(#theta^{*}_{#omega#gamma})");
        fESDList[iCut]->Add(fHistoMotherRestGammaCosAnglePt[iCut]);

        fHistoMotherRestPi0CosAnglePt[iCut] = new TH2F("ESD_MotherRestPi0_CosAngle_Pt","ESD_MotherRestPi0_CosAngle_Pt",nBinsPt, arrPtBinning,200,-1,1);
        fHistoMotherRestPi0CosAnglePt[iCut]->SetXTitle("#it{p}_{T, #omega} (GeV/#it{c})");
        fHistoMotherRestPi0CosAnglePt[iCut]->SetYTitle("cos(#theta^{*}_{#omega#pi^{0}})");
        fESDList[iCut]->Add(fHistoMotherRestPi0CosAnglePt[iCut]);
      }
      if(fDoMesonQA & 0b01000000) {
        fHistoRecoArmenterosPodolanskiPlot[iCut] = new TH2F("ESD_Armenteros_Podolanski_Plot", "ESD_Armenteros_Podolanski_Plot", 200, -1.0, 1.0, 100, 0.0, 1.0);
        fHistoRecoArmenterosPodolanskiPlot[iCut]->SetXTitle("#alpha");
        fHistoRecoArmenterosPodolanskiPlot[iCut]->SetYTitle("q_{T} (GeV/#it{c})");
        fESDList[iCut]->Add(fHistoRecoArmenterosPodolanskiPlot[iCut]);

        fHistoRecoArmenterosPodolanskiPlotwoCut[iCut] = new TH2F("ESD_Armenteros_Podolanski_Plot_woCut", "ESD_Armenteros_Podolanski_Plot_woCut", 200, -1.0, 1.0, 100, 0.0, 1.0);
        fHistoRecoArmenterosPodolanskiPlotwoCut[iCut]->SetXTitle("#alpha");
        fHistoRecoArmenterosPodolanskiPlotwoCut[iCut]->SetYTitle("q_{T} (GeV/#it{c})");
        fESDList[iCut]->Add(fHistoRecoArmenterosPodolanskiPlotwoCut[iCut]);
      }
      if(fDoMesonQA & 0b10000000) {
        fHistoMotherDalitzPlot[iCut] = new TH2F("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1","ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1",200,0,0.04,260,0.001,2.6);
        fHistoMotherDalitzPlot[iCut]->SetXTitle("#it{m}^{2}(#gamma_{0}#gamma_{1}) (GeV/#it{c}^{2})^{2}");
        fHistoMotherDalitzPlot[iCut]->SetYTitle("#it{m}^{2}(#gamma_{0/1}#gamma_{2}) (GeV/#it{c}^{2})^{2}");
        SetLogBinningYTH2(fHistoMotherDalitzPlot[iCut]);
        fESDList[iCut]->Add(fHistoMotherDalitzPlot[iCut]);

        fHistoMotherDalitzPlotwoCut[iCut] = new TH2F("ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1_woCut","ESD_Dalitz_Gamma1Gamma2_Gamma0Gamma1_woCut",200,0,0.04,260,0.001,2.6);
        fHistoMotherDalitzPlotwoCut[iCut]->SetXTitle("#it{m}^{2}(#gamma_{0}#gamma_{1}) (GeV/#it{c}^{2})^{2}");
        fHistoMotherDalitzPlotwoCut[iCut]->SetYTitle("#it{m}^{2}(#gamma_{0/1}#gamma_{2}) (GeV/#it{c}^{2})^{2}");
        SetLogBinningYTH2(fHistoMotherDalitzPlotwoCut[iCut]);
        fESDList[iCut]->Add(fHistoMotherDalitzPlotwoCut[iCut]);

        fHistoBackDalitzPlot[iCut] = new TH2F("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1","ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1",200,0,0.04,260,0.001,2.6);
        fHistoBackDalitzPlot[iCut]->SetXTitle("#it{m}^{2}(#gamma_{0}#gamma_{1}) (GeV/#it{c}^{2})^{2}");
        fHistoBackDalitzPlot[iCut]->SetYTitle("#it{m}^{2}(#gamma_{0/1}#gamma_{2}) (GeV/#it{c}^{2})^{2}");
        SetLogBinningYTH2(fHistoBackDalitzPlot[iCut]);
        fESDList[iCut]->Add(fHistoBackDalitzPlot[iCut]);

        fHistoBackDalitzPlotwoCut[iCut] = new TH2F("ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1_woCut","ESD_Dalitz_Back_Gamma1Gamma2_Gamma0Gamma1_woCut",200,0,0.04,260,0.001,2.6);
        fHistoBackDalitzPlotwoCut[iCut]->SetXTitle("#it{m}^{2}(#gamma_{0}#gamma_{1}) (GeV/#it{c}^{2})^{2}");
        fHistoBackDalitzPlotwoCut[iCut]->SetYTitle("#it{m}^{2}(#gamma_{0/1}#gamma_{2}) (GeV/#it{c}^{2})^{2}");
        SetLogBinningYTH2(fHistoBackDalitzPlotwoCut[iCut]);
        fESDList[iCut]->Add(fHistoBackDalitzPlotwoCut[iCut]);
      }


      if (fIsMC > 1){
        if(fDoMesonQA & 0b00000001) {
          fHistoPhotonPairYPt[iCut]->Sumw2();
          fHistoMotherYPt[iCut]->Sumw2();
        }
        if(fDoMesonQA & 0b00000010) {
          fHistoPhotonPairAlphaPt[iCut]->Sumw2();
          fHistoMotherAlphaPt[iCut]->Sumw2();
        }
        if(fDoMesonQA & 0b00000100) {
          fHistoPhotonPairOpenAnglePt[iCut]->Sumw2();
        }
        if(fDoMesonQA & 0b00001000) {
          fHistoPhotonPairEtaPhi[iCut]->Sumw2();
          if(fReconMethod!=2 && fReconMethod!=5) fHistoMotherConvPhotonEtaPhi[iCut]->Sumw2();
          fHistoMotherEtaPhi[iCut]->Sumw2();
        }
        if(fDoMesonQA & 0b00010000) {
          fHistoMotherPi0AnglePt[iCut]->Sumw2();
          fHistoMotherGammaAnglePt[iCut]->Sumw2();
          fHistoPi0GammaAnglePt[iCut]->Sumw2();
        }
        if(fDoMesonQA & 0b00100000) {
          fHistoMotherRestGammaCosAnglePt[iCut]->Sumw2();
          fHistoMotherRestPi0CosAnglePt[iCut]->Sumw2();
        }
        if(fDoMesonQA & 0b01000000) {
          fHistoRecoArmenterosPodolanskiPlot[iCut]->Sumw2();
          fHistoRecoArmenterosPodolanskiPlotwoCut[iCut]->Sumw2();
        }
        if(fDoMesonQA & 0b10000000) {
          fHistoMotherDalitzPlot[iCut]->Sumw2();
          fHistoMotherDalitzPlotwoCut[iCut]->Sumw2();
          fHistoBackDalitzPlot[iCut]->Sumw2();
          fHistoBackDalitzPlotwoCut[iCut]->Sumw2();
        }
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


    if(!fDoLightOutput) {
      if(fDoMesonQA & 0b00000001) {
        fHistoMCAllOmegaYPt                             = new TH2F*[fnCuts];
        fHistoMCOmegaInAccYPt                           = new TH2F*[fnCuts];
        fHistoMCPi0FromOmegaInAccYPt                    = new TH2F*[fnCuts];
        fHistoMCPi0FromAllOmegaYPt                      = new TH2F*[fnCuts];
        fHistoMCPi0PtY                                  = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00000010) {
        fHistoMCAllOmegaAlphaPt                         = new TH2F*[fnCuts];
        fHistoMCOmegaInAccAlphaPt                       = new TH2F*[fnCuts];
        fHistoMCPi0FromAllOmegaAlphaPt                  = new TH2F*[fnCuts];
        fHistoMCPi0FromOmegaInAccAlphaPt                = new TH2F*[fnCuts];
        fHistoMCPi0PtAlpha                              = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00001000) {
        fHistoMCPi0FromAllOmegaEtaPhi                   = new TH2F*[fnCuts];
        fHistoMCPi0FromOmegaInAccEtaPhi                 = new TH2F*[fnCuts];
        fHistoMCAllOmegaEtaPhi                          = new TH2F*[fnCuts];
        fHistoMCOmegaInAccEtaPhi                        = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00010000) {
        fHistoMCAllOmegaPiZeroAnglePt                   = new TH2F*[fnCuts];
        fHistoMCAllPiZeroGammaAnglePt                   = new TH2F*[fnCuts];
        fHistoMCAllOmegaGammaAnglePt                    = new TH2F*[fnCuts];
        fHistoMCInAccOmegaPiZeroAnglePt                 = new TH2F*[fnCuts];
        fHistoMCInAccPiZeroGammaAnglePt                 = new TH2F*[fnCuts];
        fHistoMCInAccOmegaGammaAnglePt                  = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b100000000) {
        fHistoMCAllOmegaPtPi0Pt                         = new TH2F*[fnCuts];
        fHistoMCInAccOmegaPtPi0Pt                       = new TH2F*[fnCuts];
        fHistoMCAllOmegaPtGammaPt                       = new TH2F*[fnCuts];
        fHistoMCInAccOmegaPtGammaPt                     = new TH2F*[fnCuts];
        if (fIsMC == 2){
          fHistoMCPi0PtJetPt                            = new TH2F*[fnCuts];
        }
      }
    }

    if (fDoPhotonQA > 0){
      fHistoMCConvGammaR                            = new TH1F*[fnCuts];
      fHistoMCConvGammaEta                          = new TH1F*[fnCuts];
    }

    fHistoTrueOmegaInvMassPt                      = new TH2F*[fnCuts];
    fHistoTruePi0FromOmegaInvMassPt               = new TH2F*[fnCuts];
    fHistoTrueGammaFromOmegaPt                    = new TH1F*[fnCuts];

    if(!fDoLightOutput){
      if(fDoMesonQA & 0b00000001) {
        fHistoTruePi0FromOmegaYPt                     = new TH2F*[fnCuts];
        fHistoTrueOmegaYPt                            = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00000010) {
        fHistoTrueOmegaAlphaPt                        = new TH2F*[fnCuts];
        fHistoTruePi0FromOmegaAlphaPt                 = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00000100) {
        fHistoTruePi0FromOmegaOpenAnglePt             = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00001000) {
        fHistoTrueOmegaEtaPhi                         = new TH2F*[fnCuts];
        fHistoTruePi0FromOmegaEtaPhi                  = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00010000) {
        fHistoTrueOmegaPi0AnglePt                     = new TH2F*[fnCuts];
        fHistoTrueOmegaGammaAnglePt                   = new TH2F*[fnCuts];
        fHistoTruePi0GammaAnglePt                     = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b00100000) {
        fHistoTrueOmegaRestGammaCosAnglePt            = new TH2F*[fnCuts];
        fHistoTrueOmegaRestPi0CosAnglePt              = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b01000000) {
        fHistoTrueArmenterosPodolanskiPlot            = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b10000000) {
        fHistoTrueDalitzPlot                          = new TH2F*[fnCuts];
      }
      if(fDoMesonQA & 0b100000000) {
        fHistoTrueOmegaPtPi0Pt                        = new TH2F*[fnCuts];
        fHistoTrueOmegaPtGammaPt                      = new TH2F*[fnCuts];
      }
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

      fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
      fHistoMCAllGammaEMCALAccPt[iCut]  = new TH1F("MC_AllGammaEMCALAcc_Pt","MC_AllGammaEMCALAcc_Pt",nBinsPt, arrPtBinning);
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

      fHistoMCPi0Pt[iCut]               = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",nBinsPt, arrPtBinning);
      fHistoMCPi0Pt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
      fHistoMCPi0WOWeightPt[iCut]       = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);

      fHistoMCPi0InAccPt[iCut]          = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",nBinsPt, arrPtBinning);
      fHistoMCPi0InAccPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
      fHistoMCPi0WOWeightInAccPt[iCut]  = new TH1F("MC_Pi0WOWeightInAcc_Pt","MC_Pi0WOWeightInAcc_Pt",nBinsPt, arrPtBinning);
      fMCList[iCut]->Add(fHistoMCPi0WOWeightInAccPt[iCut]);

      fHistoMCGammaFromAllOmegaPt[iCut]    = new TH1F("MC_GammaFromAllOmega_Pt","MC_GammaFromAllOmega_Pt",nBinsPt, arrPtBinning);
      fHistoMCGammaFromAllOmegaPt[iCut]->SetXTitle("p_{T, #gamma}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCGammaFromAllOmegaPt[iCut]);

      fHistoMCGammaFromOmegaInAccPt[iCut]    = new TH1F("MC_GammaFromOmegaInAcc_Pt","MC_GammaFromOmegaInAcc_Pt",nBinsPt, arrPtBinning);
      fHistoMCGammaFromOmegaInAccPt[iCut]->SetXTitle("p_{T, #gamma}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCGammaFromOmegaInAccPt[iCut]);

      fHistoMCOmegaDecayChannels[iCut]     = new TH1F("MC_OmegaDecayChannels","MC_OmegaDecayChannels",3,-0.5,2.5);
      fHistoMCOmegaDecayChannels[iCut]->GetXaxis()->SetBinLabel(1,"PiPlPiMiPiZero");
      fHistoMCOmegaDecayChannels[iCut]->GetXaxis()->SetBinLabel(2,"PiZeroGamma");
      fHistoMCOmegaDecayChannels[iCut]->GetXaxis()->SetBinLabel(3,"Others");
      fMCList[iCut]->Add(fHistoMCOmegaDecayChannels[iCut]);

      fHistoMCAllOmegaInvMassPt[iCut]                = new TH2F("MC_AllOmega_InvMass_Pt","MC_AllOmega_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoMCAllOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#omega}(GeV/c^{2})");
      fHistoMCAllOmegaInvMassPt[iCut]->SetYTitle("#omega p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCAllOmegaInvMassPt[iCut]);

      fHistoMCPi0FromAllOmegaInvMassPt[iCut]                = new TH2F("MC_Pi0FromAllOmega_InvMass_Pt","MC_Pi0FromAllOmega_InvMass_Pt",300,0,0.3,nBinsPt, arrPtBinning);
      fHistoMCPi0FromAllOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
      fHistoMCPi0FromAllOmegaInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaInvMassPt[iCut]);

      fHistoMCPi0FromOmegaInAccInvMassPt[iCut]              = new TH2F("MC_Pi0FromOmegaInAcc_InvMass_Pt","MC_Pi0FromOmegaInAcc_InvMass_Pt",300,0,0.3,nBinsPt, arrPtBinning);
      fHistoMCPi0FromOmegaInAccInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
      fHistoMCPi0FromOmegaInAccInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccInvMassPt[iCut]);

      fHistoMCOmegaInvMassPt[iCut] = new TH2F("MC_OmegaInvMass_Pt","MC_OmegaInvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoMCOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#omega}(GeV/c^{2})");
      fHistoMCOmegaInvMassPt[iCut]->SetYTitle("#omega p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCOmegaInvMassPt[iCut]);

      fHistoMCOmegaInAccInvMassPt[iCut] = new TH2F("MC_OmegaInAcc_InvMass_Pt","MC_OmegaInAcc_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoMCOmegaInAccInvMassPt[iCut]->SetXTitle("M_{inv,#omega}(GeV/c^{2})");
      fHistoMCOmegaInAccInvMassPt[iCut]->SetYTitle("#omega p_{T}(GeV/c)");
      fMCList[iCut]->Add(fHistoMCOmegaInAccInvMassPt[iCut]);

      if(!fDoLightOutput) {
        if(fDoMesonQA & 0b00000001) {
          fHistoMCAllOmegaYPt[iCut] = new TH2F("MC_AllOmega_Y_Pt","MC_AllOmega_Y_Pt",nBinsPt, arrPtBinning,150,-1.5,1.5);
          fHistoMCAllOmegaYPt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoMCAllOmegaYPt[iCut]->SetYTitle("Y_{#omega}");
          fMCList[iCut]->Add(fHistoMCAllOmegaYPt[iCut]);

          fHistoMCOmegaInAccYPt[iCut] = new TH2F("MC_OmegaInAcc_Y_Pt","MC_OmegaInAcc_Y_Pt",nBinsPt, arrPtBinning,150,-1.5,1.5);
          fHistoMCOmegaInAccYPt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoMCOmegaInAccYPt[iCut]->SetYTitle("Y_{#omega}");
          fMCList[iCut]->Add(fHistoMCOmegaInAccYPt[iCut]);

          fHistoMCPi0FromAllOmegaYPt[iCut] = new TH2F("MC_Pi0FromAllOmega_Y_Pt","MC_Pi0FromAllOmega_Y_Pt",200,0,20,150,-1.5,1.5);
          fHistoMCPi0FromAllOmegaYPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
          fHistoMCPi0FromAllOmegaYPt[iCut]->SetYTitle("Y_{#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaYPt[iCut]);

          fHistoMCPi0FromOmegaInAccYPt[iCut] = new TH2F("MC_Pi0FromOmegaInAcc_Y_Pt","MC_Pi0FromOmegaInAcc_Y_Pt",200,0,20,150,-1.5,1.5);
          fHistoMCPi0FromOmegaInAccYPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
          fHistoMCPi0FromOmegaInAccYPt[iCut]->SetYTitle("Y_{#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccYPt[iCut]);

          fHistoMCPi0PtY[iCut]            = new TH2F("MC_Pi0_Pt_Y","MC_Pi0_Pt_Y",300,0.03,30.,150,-1.5,1.5);
          fHistoMCPi0PtY[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
          fHistoMCPi0PtY[iCut]->SetYTitle("Y_{#pi^{0}}");
          SetLogBinningXTH2(fHistoMCPi0PtY[iCut]);
          fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);
        }
        if(fDoMesonQA & 0b00000010) {
          fHistoMCAllOmegaAlphaPt[iCut] = new TH2F("MC_AllOmega_Alpha_Pt","MC_AllOmega_Alpha_Pt",nBinsPt, arrPtBinning,200,-1,1);
          fHistoMCAllOmegaAlphaPt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoMCAllOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#omega}");
          fMCList[iCut]->Add(fHistoMCAllOmegaAlphaPt[iCut]);

          fHistoMCOmegaInAccAlphaPt[iCut] = new TH2F("MC_OmegaInAcc_Alpha_Pt","MC_OmegaInAcc_Alpha_Pt",nBinsPt, arrPtBinning,200,-1,1);
          fHistoMCOmegaInAccAlphaPt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoMCOmegaInAccAlphaPt[iCut]->SetYTitle("#alpha_{#omega}");
          fMCList[iCut]->Add(fHistoMCOmegaInAccAlphaPt[iCut]);

          fHistoMCPi0FromAllOmegaAlphaPt[iCut] = new TH2F("MC_Pi0FromAllOmega_Alpha_Pt","MC_Pi0FromAllOmega_Alpha_Pt",200,0,20,200,-1,1);
          fHistoMCPi0FromAllOmegaAlphaPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
          fHistoMCPi0FromAllOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaAlphaPt[iCut]);

          fHistoMCPi0FromOmegaInAccAlphaPt[iCut] = new TH2F("MC_Pi0FromOmegaInAcc_Alpha_Pt","MC_Pi0FromOmegaInAcc_Alpha_Pt",200,0,20,200,-1,1);
          fHistoMCPi0FromOmegaInAccAlphaPt[iCut]->SetXTitle("p_{T,#pi^{0}}(GeV/c)");
          fHistoMCPi0FromOmegaInAccAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccAlphaPt[iCut]);

          fHistoMCPi0PtAlpha[iCut]        = new TH2F("MC_Pi0_Pt_Alpha","MC_Pi0_Pt_Alpha",300,0.03,30.,200,-1,1);
          SetLogBinningXTH2(fHistoMCPi0PtAlpha[iCut]);
          fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);
        }
        if(fDoMesonQA & 0b00001000) {
          fHistoMCPi0FromAllOmegaEtaPhi[iCut] = new TH2F("MC_Pi0FromAllOmega_Eta_Phi","MC_Pi0FromAllOmega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
          fHistoMCPi0FromAllOmegaEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}}(rad)");
          fHistoMCPi0FromAllOmegaEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCPi0FromAllOmegaEtaPhi[iCut]);

          fHistoMCPi0FromOmegaInAccEtaPhi[iCut] = new TH2F("MC_Pi0FromOmegaInAcc_Eta_Phi","MC_Pi0FromOmegaInAcc_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
          fHistoMCPi0FromOmegaInAccEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}}(rad)");
          fHistoMCPi0FromOmegaInAccEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCPi0FromOmegaInAccEtaPhi[iCut]);

          fHistoMCAllOmegaEtaPhi[iCut]   = new TH2F("MC_AllOmega_Eta_Phi","MC_AllOmega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
          fHistoMCAllOmegaEtaPhi[iCut]->SetXTitle("#phi_{#omega}(rad)");
          fHistoMCAllOmegaEtaPhi[iCut]->SetYTitle("#eta_{#omega}");
          fMCList[iCut]->Add(fHistoMCAllOmegaEtaPhi[iCut]);

          fHistoMCOmegaInAccEtaPhi[iCut]   = new TH2F("MC_OmegaInAcc_Eta_Phi","MC_OmegaInAcc_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
          fHistoMCOmegaInAccEtaPhi[iCut]->SetXTitle("#phi_{#omega}(rad)");
          fHistoMCOmegaInAccEtaPhi[iCut]->SetYTitle("#eta_{#omega}");
          fMCList[iCut]->Add(fHistoMCOmegaInAccEtaPhi[iCut]);
        }
        if(fDoMesonQA & 0b00010000) {
          fHistoMCAllOmegaPiZeroAnglePt[iCut] = new TH2F("MC_AllOmegaPiZero_Angle_Pt","MC_AllOmegaPiZero_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoMCAllOmegaPiZeroAnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoMCAllOmegaPiZeroAnglePt[iCut]->SetYTitle("#theta_{#omega,#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCAllOmegaPiZeroAnglePt[iCut]);

          fHistoMCAllPiZeroGammaAnglePt[iCut] = new TH2F("MC_AllPiZeroGamma_Angle_Pt","MC_AllPiZeroGamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoMCAllPiZeroGammaAnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoMCAllPiZeroGammaAnglePt[iCut]->SetYTitle("#theta_{#pi^{0},#gamma}");
          fMCList[iCut]->Add(fHistoMCAllPiZeroGammaAnglePt[iCut]);

          fHistoMCAllOmegaGammaAnglePt[iCut] = new TH2F("MC_AllOmegaGamma_Angle_Pt","MC_AllOmegaGamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoMCAllOmegaGammaAnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoMCAllOmegaGammaAnglePt[iCut]->SetYTitle("#theta_{#omega,#gamma}");
          fMCList[iCut]->Add(fHistoMCAllOmegaGammaAnglePt[iCut]);

          fHistoMCInAccOmegaPiZeroAnglePt[iCut] = new TH2F("MC_InAccOmegaPiZero_Angle_Pt","MC_InAccOmegaPiZero_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoMCInAccOmegaPiZeroAnglePt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoMCInAccOmegaPiZeroAnglePt[iCut]->SetYTitle("#theta_{#omega,#pi^{0}}");
          fMCList[iCut]->Add(fHistoMCInAccOmegaPiZeroAnglePt[iCut]);

          fHistoMCInAccPiZeroGammaAnglePt[iCut] = new TH2F("MC_InAccPiZeroGamma_Angle_Pt","MC_InAccPiZeroGamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoMCInAccPiZeroGammaAnglePt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoMCInAccPiZeroGammaAnglePt[iCut]->SetYTitle("#theta_{#pi^{0},#gamma}");
          fMCList[iCut]->Add(fHistoMCInAccPiZeroGammaAnglePt[iCut]);

          fHistoMCInAccOmegaGammaAnglePt[iCut] = new TH2F("MC_InAccOmegaGamma_Angle_Pt","MC_InAccOmegaGamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoMCInAccOmegaGammaAnglePt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoMCInAccOmegaGammaAnglePt[iCut]->SetYTitle("#theta_{#omega,#gamma}");
          fMCList[iCut]->Add(fHistoMCInAccOmegaGammaAnglePt[iCut]);
        }
        if(fDoMesonQA & 0b100000000) {
          fHistoMCAllOmegaPtPi0Pt[iCut] = new TH2F("MC_All_OmegaPt_Pi0Pt","MC_All_OmegaPt_Pi0Pt",nBinsPt, arrPtBinning,200,0,20);
          fHistoMCAllOmegaPtPi0Pt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoMCAllOmegaPtPi0Pt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
          fMCList[iCut]->Add(fHistoMCAllOmegaPtPi0Pt[iCut]);

          fHistoMCInAccOmegaPtPi0Pt[iCut] = new TH2F("MC_InAcc_OmegaPt_Pi0Pt","MC_InAcc_OmegaPt_Pi0Pt",nBinsPt, arrPtBinning,200,0,20);
          fHistoMCInAccOmegaPtPi0Pt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoMCInAccOmegaPtPi0Pt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
          fMCList[iCut]->Add(fHistoMCInAccOmegaPtPi0Pt[iCut]);

          fHistoMCAllOmegaPtGammaPt[iCut] = new TH2F("MC_All_OmegaPt_GammaPt","MC_All_OmegaPt_GammaPt",nBinsPt, arrPtBinning,200,0,20);
          fHistoMCAllOmegaPtGammaPt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoMCAllOmegaPtGammaPt[iCut]->SetYTitle("#gamma p_{T}(GeV/c)");
          fMCList[iCut]->Add(fHistoMCAllOmegaPtGammaPt[iCut]);

          fHistoMCInAccOmegaPtGammaPt[iCut] = new TH2F("MC_InAcc_OmegaPt_GammaPt","MC_InAcc_OmegaPt_GammaPt",nBinsPt, arrPtBinning,200,0,20);
          fHistoMCInAccOmegaPtGammaPt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoMCInAccOmegaPtGammaPt[iCut]->SetYTitle("#gamma p_{T}(GeV/c)");
          fMCList[iCut]->Add(fHistoMCInAccOmegaPtGammaPt[iCut]);
        }

        if (fIsMC > 1){
          if(fDoMesonQA & 0b00000001) {
            fHistoMCAllOmegaYPt[iCut]->Sumw2();
            fHistoMCOmegaInAccYPt[iCut]->Sumw2();
            fHistoMCPi0FromAllOmegaYPt[iCut]->Sumw2();
            fHistoMCPi0FromOmegaInAccYPt[iCut]->Sumw2();
            fHistoMCPi0PtY[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00000010) {
            fHistoMCPi0PtAlpha[iCut]->Sumw2();
            fHistoMCAllOmegaAlphaPt[iCut]->Sumw2();
            fHistoMCOmegaInAccAlphaPt[iCut]->Sumw2();
            fHistoMCPi0FromAllOmegaAlphaPt[iCut]->Sumw2();
            fHistoMCPi0FromOmegaInAccAlphaPt[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00001000) {
            fHistoMCPi0FromAllOmegaEtaPhi[iCut]->Sumw2();
            fHistoMCPi0FromOmegaInAccEtaPhi[iCut]->Sumw2();
            fHistoMCAllOmegaEtaPhi[iCut]->Sumw2();
            fHistoMCOmegaInAccEtaPhi[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00010000) {
            fHistoMCAllOmegaPiZeroAnglePt[iCut]->Sumw2();
            fHistoMCAllPiZeroGammaAnglePt[iCut]->Sumw2();
            fHistoMCAllOmegaGammaAnglePt[iCut]->Sumw2();
            fHistoMCInAccOmegaPiZeroAnglePt[iCut]->Sumw2();
            fHistoMCInAccPiZeroGammaAnglePt[iCut]->Sumw2();
            fHistoMCInAccOmegaGammaAnglePt[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b100000000) {
            fHistoMCAllOmegaPtPi0Pt[iCut]->Sumw2();
            fHistoMCInAccOmegaPtPi0Pt[iCut]->Sumw2();
            fHistoMCAllOmegaPtGammaPt[iCut]->Sumw2();
            fHistoMCInAccOmegaPtGammaPt[iCut]->Sumw2();
          }
        }
      }

      if (fIsMC > 1){
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fHistoMCPi0WOWeightInAccPt[iCut]->Sumw2();
        fHistoMCPi0WOEvtWeightPt[iCut]  = new TH1F("MC_Pi0_WOEventWeights_Pt","MC_Pi0_WOEventWeights_Pt",400,0,40);
        fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
        fHistoMCGammaFromAllOmegaPt[iCut]->Sumw2();
        fHistoMCGammaFromOmegaInAccPt[iCut]->Sumw2();
        fHistoMCOmegaInAccInvMassPt[iCut]->Sumw2();
        fHistoMCOmegaInvMassPt[iCut]->Sumw2();
        fHistoMCAllOmegaInvMassPt[iCut]->Sumw2();
        fHistoMCPi0FromAllOmegaInvMassPt[iCut]->Sumw2();
        fHistoMCPi0FromOmegaInAccInvMassPt[iCut]->Sumw2();
        fHistoMCOmegaDecayChannels[iCut]->Sumw2();

        if( (!fDoLightOutput) && (fDoMesonQA & 0b100000000) && (fIsMC == 2) ){
          fHistoMCPi0PtJetPt[iCut]      = new TH2F("MC_Pi0_Pt_JetPt","MC_Pi0_Pt_JetPt",300,0.03,30.,200,0,200);
          fHistoMCPi0PtJetPt[iCut]->SetXTitle("#pi^{0} #it{p}_{T}(GeV/#it{c})");
          fHistoMCPi0PtJetPt[iCut]->SetYTitle("jet #it{p}_{T}(GeV/#it{c})");
          fHistoMCPi0PtJetPt[iCut]->Sumw2();
          SetLogBinningXTH2(fHistoMCPi0PtJetPt[iCut]);
          fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
        }
      }

      fTrueList[iCut]                           = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringNeutralPion.Data(),cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      fHistoTrueOmegaInvMassPt[iCut]                = new TH2F("True_Omega_InvMass_Pt","True_Omega_InvMass_Pt",200,0.0,1.6,nBinsPt, arrPtBinning);
      fHistoTrueOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#omega}(GeV/c^{2})");
      fHistoTrueOmegaInvMassPt[iCut]->SetYTitle("#omega p_{T}(GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueOmegaInvMassPt[iCut]);

      fHistoTruePi0FromOmegaInvMassPt[iCut]          = new TH2F("True_Pi0FromOmega_InvMass_Pt","True_Pi0FromOmega_InvMass_Pt",300,0,0.3,200,0.,20.);
      fHistoTruePi0FromOmegaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
      fHistoTruePi0FromOmegaInvMassPt[iCut]->SetYTitle("p_{T, #pi^{0}}(GeV/c)");
      fTrueList[iCut]->Add(fHistoTruePi0FromOmegaInvMassPt[iCut]);

      fHistoTrueGammaFromOmegaPt[iCut] = new TH1F("True_GammaFromOmega_Pt","True_GammaFromOmega_Pt",nBinsPt, arrPtBinning);
      fHistoTrueGammaFromOmegaPt[iCut]->SetXTitle("p_{T,#gamma}(GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueGammaFromOmegaPt[iCut]);

      if(!fDoLightOutput) {
        if(fDoMesonQA & 0b00000001) {
          fHistoTrueOmegaYPt[iCut] = new TH2F("True_Omega_Y_Pt","True_Omega_Y_Pt",nBinsPt, arrPtBinning,150,-1.5,1.5);
          fHistoTrueOmegaYPt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoTrueOmegaYPt[iCut]->SetYTitle("Y_{#omega}");
          fTrueList[iCut]->Add(fHistoTrueOmegaYPt[iCut]);

          fHistoTruePi0FromOmegaYPt[iCut] = new TH2F("True_Pi0FromOmega_Y_Pt","True_Pi0FromOmega_Y_Pt",nBinsPt, arrPtBinning,150,-1.5,1.5);
          fHistoTruePi0FromOmegaYPt[iCut]->SetXTitle("p_{T,pi^{0}}(GeV/c)");
          fHistoTruePi0FromOmegaYPt[iCut]->SetYTitle("Y_{#pi^{0}}");
          fTrueList[iCut]->Add(fHistoTruePi0FromOmegaYPt[iCut]);
        }
        if(fDoMesonQA & 0b00000010) {
          fHistoTrueOmegaAlphaPt[iCut] = new TH2F("True_Omega_Alpha_Pt","True_Omega_Alpha_Pt",nBinsPt, arrPtBinning,200,-1,1);
          fHistoTrueOmegaAlphaPt[iCut]->SetXTitle("p_{T,#omega}(GeV/c)");
          fHistoTrueOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#omega}");
          fTrueList[iCut]->Add(fHistoTrueOmegaAlphaPt[iCut]);

          fHistoTruePi0FromOmegaAlphaPt[iCut] = new TH2F("True_Pi0FromOmega_Alpha_Pt","MC_Pi0FromOmega_Alpha_Pt",nBinsPt, arrPtBinning,200,-1,1);
          fHistoTruePi0FromOmegaAlphaPt[iCut]->SetXTitle("p_{T,pi^{0}}(GeV/c)");
          fHistoTruePi0FromOmegaAlphaPt[iCut]->SetYTitle("#alpha_{#pi^{0}}");
          fTrueList[iCut]->Add(fHistoTruePi0FromOmegaAlphaPt[iCut]);
        }
        if(fDoMesonQA & 0b00000100) {
          fHistoTruePi0FromOmegaOpenAnglePt[iCut] = new TH2F("True_Pi0FromOmega_OpenAngle_Pt","True_Pi0FromOmega_OpenAngle_Pt",nBinsPt, arrPtBinning,100,0,1);
          fHistoTruePi0FromOmegaOpenAnglePt[iCut]->SetXTitle("p_{T, #pi^{0}}(GeV/c)");
          fHistoTruePi0FromOmegaOpenAnglePt[iCut]->SetYTitle("#theta_{#pi^{0}}");
          fTrueList[iCut]->Add(fHistoTruePi0FromOmegaOpenAnglePt[iCut]);
        }
        if(fDoMesonQA & 0b00001000) {
          fHistoTrueOmegaEtaPhi[iCut]   = new TH2F("True_Omega_Eta_Phi","True_Omega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
          fHistoTrueOmegaEtaPhi[iCut]->SetXTitle("#phi_{#omega}(rad)");
          fHistoTrueOmegaEtaPhi[iCut]->SetYTitle("#eta_{#omega}");
          fTrueList[iCut]->Add(fHistoTrueOmegaEtaPhi[iCut]);

          fHistoTruePi0FromOmegaEtaPhi[iCut]   = new TH2F("True_Pi0FromOmega_Eta_Phi","True_Pi0FromOmega_Eta_Phi",600,0,2*TMath::Pi(),200,-1,1);
          fHistoTruePi0FromOmegaEtaPhi[iCut]->SetXTitle("#phi_{#pi^{0}}(rad)");
          fHistoTruePi0FromOmegaEtaPhi[iCut]->SetYTitle("#eta_{#pi^{0}}");
          fTrueList[iCut]->Add(fHistoTruePi0FromOmegaEtaPhi[iCut]);
        }
        if(fDoMesonQA & 0b00010000) {
          fHistoTrueOmegaPi0AnglePt[iCut] = new TH2F("True_OmegaPiZero_Angle_Pt","True_OmegaPiZero_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoTrueOmegaPi0AnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoTrueOmegaPi0AnglePt[iCut]->SetYTitle("#theta_{#omega,#pi^{0}}");
          fTrueList[iCut]->Add(fHistoTrueOmegaPi0AnglePt[iCut]);

          fHistoTruePi0GammaAnglePt[iCut] = new TH2F("True_PiZeroGamma_Angle_Pt","True_PiZeroGamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoTruePi0GammaAnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoTruePi0GammaAnglePt[iCut]->SetYTitle("#theta_{#pi_{0},#gamma}");
          fTrueList[iCut]->Add(fHistoTruePi0GammaAnglePt[iCut]);

          fHistoTrueOmegaGammaAnglePt[iCut] = new TH2F("True_OmegaGamma_Angle_Pt","True_OmegaGamma_Angle_Pt",nBinsPt, arrPtBinning,360,0,TMath::Pi());
          fHistoTrueOmegaGammaAnglePt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoTrueOmegaGammaAnglePt[iCut]->SetYTitle("#theta_{#omega,#gamma}");
          fTrueList[iCut]->Add(fHistoTrueOmegaGammaAnglePt[iCut]);
        }
        if(fDoMesonQA & 0b00100000) {
          fHistoTrueOmegaRestGammaCosAnglePt[iCut] = new TH2F("True_OmegaRestGamma_CosAngle_Pt","True_OmegaRestGamma_CosAngle_Pt",nBinsPt, arrPtBinning,200,-1,1);
          fHistoTrueOmegaRestGammaCosAnglePt[iCut]->SetXTitle("#it{p}_{T, #omega} (GeV/#it{c})");
          fHistoTrueOmegaRestGammaCosAnglePt[iCut]->SetYTitle("cos(#theta^{*}_{#omega#gamma})");
          fTrueList[iCut]->Add(fHistoTrueOmegaRestGammaCosAnglePt[iCut]);

          fHistoTrueOmegaRestPi0CosAnglePt[iCut] = new TH2F("True_OmegaRestPi0_CosAngle_Pt","True_OmegaRestPi0_CosAngle_Pt",nBinsPt, arrPtBinning,200,-1,1);
          fHistoTrueOmegaRestPi0CosAnglePt[iCut]->SetXTitle("#it{p}_{T, #omega} (GeV/#it{c})");
          fHistoTrueOmegaRestPi0CosAnglePt[iCut]->SetYTitle("cos(#theta^{*}_{#omega#pi^{0}})");
          fTrueList[iCut]->Add(fHistoTrueOmegaRestPi0CosAnglePt[iCut]);
        }
        if(fDoMesonQA & 0b01000000) {
          fHistoTrueArmenterosPodolanskiPlot[iCut] = new TH2F("True_Armenteros_Podolanski_Plot", "True_Armenteros_Podolanski_Plot", 200, -1.0, 1.0, 100, 0.0, 1.0);
          fHistoTrueArmenterosPodolanskiPlot[iCut]->SetXTitle("#alpha");
          fHistoTrueArmenterosPodolanskiPlot[iCut]->SetYTitle("q_{T} (GeV/#it{c})");
          fTrueList[iCut]->Add(fHistoTrueArmenterosPodolanskiPlot[iCut]);
        }
        if(fDoMesonQA & 0b10000000) {
          fHistoTrueDalitzPlot[iCut] = new TH2F("True_Dalitz_Gamma1Gamma2_Gamma0Gamma1","True_Dalitz_Gamma1Gamma2_Gamma0Gamma1",200,0,0.04,260,0.001,2.6);
          fHistoTrueDalitzPlot[iCut]->SetXTitle("#it{m}^{2}(#gamma_{0}#gamma_{1}) (GeV/#it{c}^{2})^{2}");
          fHistoTrueDalitzPlot[iCut]->SetYTitle("#it{m}^{2}(#gamma_{0/1}#gamma_{2}) (GeV/#it{c}^{2})^{2}");
          SetLogBinningYTH2(fHistoTrueDalitzPlot[iCut]);
          fTrueList[iCut]->Add(fHistoTrueDalitzPlot[iCut]);
        }
        if(fDoMesonQA & 0b100000000) {
          fHistoTrueOmegaPtPi0Pt[iCut] = new TH2F("True_OmegaPt_Pi0Pt","True_OmegaPt_Pi0Pt",200,0,20,nBinsPt, arrPtBinning);
          fHistoTrueOmegaPtPi0Pt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoTrueOmegaPtPi0Pt[iCut]->SetYTitle("#pi^{0} p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueOmegaPtPi0Pt[iCut]);

          fHistoTrueOmegaPtGammaPt[iCut] = new TH2F("True_OmegaPt_GammaPt","True_OmegaPt_GammaPt",nBinsPt, arrPtBinning,200,0,20);
          fHistoTrueOmegaPtGammaPt[iCut]->SetXTitle("#omega p_{T}(GeV/c)");
          fHistoTrueOmegaPtGammaPt[iCut]->SetYTitle("#gamma p_{T}(GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueOmegaPtGammaPt[iCut]);
        }

        if (fIsMC > 1){
          if(fDoMesonQA & 0b00000001) {
            fHistoTruePi0FromOmegaYPt[iCut]->Sumw2();
            fHistoTrueOmegaYPt[iCut]->Sumw2();
            fHistoTrueOmegaYPt[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00000010) {
            fHistoTruePi0FromOmegaAlphaPt[iCut]->Sumw2();
            fHistoTrueOmegaAlphaPt[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00000100) {
            fHistoTruePi0FromOmegaOpenAnglePt[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00001000) {
            fHistoTruePi0FromOmegaEtaPhi[iCut]->Sumw2();
            fHistoTrueOmegaEtaPhi[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00010000) {
            fHistoTrueOmegaPi0AnglePt[iCut]->Sumw2();
            fHistoTrueOmegaGammaAnglePt[iCut]->Sumw2();
            fHistoTruePi0GammaAnglePt[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b00100000) {
            fHistoTrueOmegaRestGammaCosAnglePt[iCut]->Sumw2();
            fHistoTrueOmegaRestPi0CosAnglePt[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b01000000) {
            fHistoTrueArmenterosPodolanskiPlot[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b10000000) {
            fHistoTrueDalitzPlot[iCut]->Sumw2();
          }
          if(fDoMesonQA & 0b100000000) {
            fHistoTrueOmegaPtPi0Pt[iCut]->Sumw2();
            fHistoTrueOmegaPtGammaPt[iCut]->Sumw2();
          }
        }
      }

      fHistoMultipleCountTruePi0[iCut]            = new TH1F("ESD_TrueMultipleCountPi0","ESD_TrueMultipleCountPi0",10,1,11);
      fTrueList[iCut]->Add(fHistoMultipleCountTruePi0[iCut]);

      if (fIsMC > 1){
        fHistoMultipleCountTruePi0[iCut]->Sumw2();
        fHistoTrueOmegaInvMassPt[iCut]->Sumw2();
        fHistoTruePi0FromOmegaInvMassPt[iCut]->Sumw2();
        fHistoTrueGammaFromOmegaPt[iCut]->Sumw2();
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

  for(Int_t iMatcherTask = 0; iMatcherTask < 3; iMatcherTask++){
    AliCaloTrackMatcher* temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }


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
  fInputEvent = InputEvent();

  if(fIsMC > 0) fMCEvent = MCEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE || fV0Reader->GetErrorAODRelabeling()) eventQuality = 2;// incomplete event
  // Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete => skip broken event/files
  if(eventQuality == 2 || eventQuality == 3){
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC > 1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }

  if(fInputEvent->IsA()==AliAODEvent::Class()){
    if (!TGeoGlobalMagField::Instance()->GetField()) fInputEvent->InitMagneticField();
    // fInputEvent->InitMagneticField();
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
      Float_t pthard = -1;
  //     cout << fMCEvent << endl;
      Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC, pthard, fInputEvent);
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
      fHistoNEventsMinGamma[iCut]->Fill(0.,fWeightJetJetMC);
    }
    if(fClusterCandidates->GetEntries()>0 && fGammaCandidates->GetEntries()>1){
      fHistoNEventsMinGamma[iCut]->Fill(1,fWeightJetJetMC);
    }
    if(fClusterCandidates->GetEntries()>2){
      fHistoNEventsMinGamma[iCut]->Fill(2,fWeightJetJetMC);
    }
    if(fGammaCandidates->GetEntries()>2){
      fHistoNEventsMinGamma[iCut]->Fill(3,fWeightJetJetMC);
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

    CalculatePi0Candidates();                                                   // Combine Gammas from conversion and from calo to pi0
    CalculateOmegaCandidates();                                                 // Combine Gammas from conversion and from calo with pi0 to omegas if swapping method for background is chosen it will be called inside this function
    if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg() ) == 2 ) {
      CalculatePi0RotationBackground();
    }
    if( !( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg() ) ) {
      CalculateBackground();                                                    // Combinatorial Background
      UpdateEventByEventData();                                                 // Store Event for mixed Events
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
    fPi0Candidates->Clear(); // delete pi0 candidates
    dropOutGammas_CALO.clear();
    dropOutGammas_PCM.clear();
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
  // CHANGED  // TODO
  TClonesArray * arrClustersProcess = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersProcess->GetEntries();
  }
  // nclus = fInputEvent->GetNumberOfCaloClusters();

//   cout << nclus << endl;

  if(nclus == 0)  return;

  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  if(fReconMethod==2){
    // match tracks to clusters
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC, kFALSE); // added kFALSE
  }

  // vertex
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Int_t i = 0; i < nclus; i++){
    AliVCluster* clus = NULL;

    // CHANGED // TODO
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersProcess)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }
    if(!clus) continue;
    // if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    // else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    //
    // if (!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)) {
      delete clus;
      continue;
    }

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate = new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef((Long_t)i);

    // CHANGED // TODO // From AliAnalysisTaskHeavyNeutralMesonToGG
    PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent));

    // get MC label
    if(fIsMC>0){
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
      // cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k<(Int_t)clus->GetNLabels(); k++){
          PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
          // Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
          // cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
        }
      }
    }

    fIsFromMBHeader                 = kTRUE;
    fIsOverlappingWithOtherHeader   = kFALSE;
    //TString periodName         = fV0Reader->GetPeriodName();
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0){
        fIsFromMBHeader = kFALSE;
      }
      if (clus->GetNLabels()>1){
        Int_t* mclabelsCluster = clus->GetLabels();
        for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent) == 0){
            fIsOverlappingWithOtherHeader = kTRUE;
          }
        }
      }
    }

    // CHANGED // TODO from AliAnalysisTaskHeavyNeutralMesonToGG
    // if (fIsOverlappingWithOtherHeader)
    if (fIsFromMBHeader && fIsOverlappingWithOtherHeader && fIsMC > 0 && fHistoClusOverlapHeadersGammaPt[fiCut])
      fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);

    if (fIsFromMBHeader /*&& !fIsOverlappingWithOtherHeader*/){
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
  if (TruePhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return;

  if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;

  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, fEnableSortForClusMC);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
  AliAODMCParticle *Photon = NULL;
  if (!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray){
    if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
      else return;
  }else {
    AliInfo("AODMCTrackArray could not be loaded");
    return;
  }

  if(Photon == NULL){
  //  cout << "no photon" << endl;
    return;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, fEnableSortForClusMC);

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
      Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
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
      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
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
        if(fIsMC>0){
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessTruePhotonCandidates(PhotonCandidate);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessTruePhotonCandidatesAOD(PhotonCandidate);
        }
      } else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
    }
  }
  if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
      if(!PhotonCandidate) continue;
      fIsFromMBHeader = kTRUE;
      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
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
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;
  AliAODMCParticle *posDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
  AliAODMCParticle *negDaughter = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};

  if(posDaughter->GetMother() != negDaughter->GetMother()) return;

  else if(posDaughter->GetMother() == -1) return;

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

  AliAODMCParticle *Photon = (AliAODMCParticle*) fAODMCTrackArray->At(posDaughter->GetMother());
  if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

  if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5) return; // check if the daughters come from a conversion

  TruePhotonCandidate->SetIsTrueConvertedPhoton();
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  // Process True Photons
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};
  if(posDaughter->GetMother(0) != negDaughter->GetMother(0)) return;
  else if(posDaughter->GetMother(0) == -1) return;

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

  TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);

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
//  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
//  Double_t mcProdVtxX   = primVtxMC->GetX();
//  Double_t mcProdVtxY   = primVtxMC->GetY();
//  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if (!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {

    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
    if (!particle) continue;

    // There are nearly no secondary omegas and this check takes a good amount
    // of time. This is way this is commented out for now.
    // Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    // if (!isPrimary) continue;

    // fill histograms for all true omegas
    if(particle->GetPdgCode() == 223){

      Bool_t DecaysToPiPlPiMiPiZero = kFALSE;

      // TParticle *posPion = 0x0;
      // TParticle *negPion = 0x0;
      // TParticle *neutPion = 0x0;
      AliAODMCParticle *posPion = 0x0;
      AliAODMCParticle *negPion = 0x0;
      AliAODMCParticle *neutPion = 0x0;

      for(Int_t index= particle->GetDaughterFirst();index <= particle->GetDaughterLast();index++){
        if(index < 0) continue;
        AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
        // TParticle *temp = (TParticle*)fMCEvent->Particle(index);
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
        fHistoMCOmegaDecayChannels[fiCut]->Fill(0.,fWeightJetJetMC);
      } else if (particle->GetNDaughters()==2){
        // TParticle *gamma2 = 0x0;
        // TParticle *pi0 = 0x0;
        AliAODMCParticle *gamma2 = 0x0;
        AliAODMCParticle *pi0 = 0x0;

        for(Int_t index = particle->GetDaughterFirst();index <= particle->GetDaughterLast();index++){
          if(index < 0) continue;
          AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(index));
          // TParticle *temp = (TParticle*)fMCEvent->Particle(index);
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
          fHistoMCOmegaDecayChannels[fiCut]->Fill(1,fWeightJetJetMC);
          fHistoMCAllOmegaInvMassPt[fiCut]->Fill(((particle->GetCalcMass()) ? particle->GetCalcMass() : particle->M()),particle->Pt(),fWeightJetJetMC);
          fHistoMCGammaFromAllOmegaPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
          fHistoMCPi0FromAllOmegaInvMassPt[fiCut]->Fill( ((pi0->GetCalcMass()) ? pi0->GetCalcMass() : pi0->M()),pi0->Pt(),fWeightJetJetMC);

          if(!fDoLightOutput) {

            if(fDoMesonQA & 0b00000001) {
              Double_t Pi0Y = 10.;
              if(pi0->E() - pi0->Pz() == 0 || pi0->E() + pi0->Pz() == 0){
                Pi0Y=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                Pi0Y = 0.5*(TMath::Log((pi0->E()+pi0->Pz()) / (pi0->E()-pi0->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }
              fHistoMCPi0FromAllOmegaYPt[fiCut]->Fill(pi0->Pt(),Pi0Y,fWeightJetJetMC);

              Double_t OmegaY = 10.;
              if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
                OmegaY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                OmegaY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }
              fHistoMCAllOmegaYPt[fiCut]->Fill(particle->Pt(),OmegaY,fWeightJetJetMC);
            }
            if(fDoMesonQA & 0b00000010) {
              Double_t alpha = (pi0->E() - gamma2->E())/(pi0->E() + gamma2->E());
              fHistoMCAllOmegaAlphaPt[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);
              //check whether pi0 decayed into two gammas
              if (pi0->GetNDaughters()==2 && pi0->GetDaughterFirst()>-1 && pi0->GetDaughterLast()>-1){
                AliAODMCParticle* gamma0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(pi0->GetDaughterFirst()));
                AliAODMCParticle* gamma1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(pi0->GetDaughterLast()));
                // TParticle *gamma0 = (TParticle*)fMCEvent->Particle(pi0->GetFirstDaughter());
                // TParticle *gamma1 = (TParticle*)fMCEvent->Particle(pi0->GetLastDaughter());
                if (gamma0->GetPdgCode()==22 && gamma1->GetPdgCode()==22){

                  //plot pi0 alpha
                  Double_t pi0alpha = (gamma0->E() - gamma1->E())/(gamma0->E() + gamma1->E());
                  fHistoMCPi0FromAllOmegaAlphaPt[fiCut]->Fill(pi0->Pt(),pi0alpha,fWeightJetJetMC);
                }
              }
            }
            if(fDoMesonQA & 0b00001000) {
              fHistoMCPi0FromAllOmegaEtaPhi[fiCut]->Fill(pi0->Phi(),pi0->Eta(),fWeightJetJetMC);
              fHistoMCAllOmegaEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
            }
            if(fDoMesonQA & 0b00010000) {
              //get angles and plot
              TVector3 omegavector = TVector3(particle->Px(),particle->Py(),particle->Pz());
              TVector3 pi0vector = TVector3(pi0->Px(),pi0->Py(),pi0->Pz());
              TVector3 gamma2vector = TVector3(gamma2->Px(),gamma2->Py(),gamma2->Pz());
              fHistoMCAllOmegaPiZeroAnglePt[fiCut]->Fill(particle->Pt(),TMath::Pi() - pi0vector.Angle(omegavector),fWeightJetJetMC);
              fHistoMCAllPiZeroGammaAnglePt[fiCut]->Fill(particle->Pt(),pi0vector.Angle(gamma2vector),fWeightJetJetMC);
              fHistoMCAllOmegaGammaAnglePt[fiCut]->Fill(particle->Pt(),TMath::Pi() - omegavector.Angle(gamma2vector),fWeightJetJetMC);
            }
            if(fDoMesonQA & 0b100000000) {
              fHistoMCAllOmegaPtPi0Pt[fiCut]->Fill(particle->Pt(),pi0->Pt(),fWeightJetJetMC);
              fHistoMCAllOmegaPtGammaPt[fiCut]->Fill(particle->Pt(),gamma2->Pt(),fWeightJetJetMC);
            }
          }
        } else{
          fHistoMCOmegaDecayChannels[fiCut]->Fill(2,fWeightJetJetMC);
        }
      } else{
        fHistoMCOmegaDecayChannels[fiCut]->Fill(2,fWeightJetJetMC);
      }
    }

    //fill histograms for omegas in acceptance
    Int_t labelNeutPion = -1;
    Int_t labelGamma = -1;
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
       ->MesonIsSelectedAODMCPiZeroGamma(particle,fAODMCTrackArray,labelNeutPion,labelGamma,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
     if (labelNeutPion > -1 && labelGamma > -1){
      AliAODMCParticle* neutPion = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelNeutPion));
      // TParticle *neutPion    = fMCEvent->Particle(labelNeutPion);

      //fill histograms for acceptance correction
      fHistoMCOmegaInvMassPt[fiCut]->Fill(TMath::Sqrt((particle->E())*(particle->E())-(particle->P())*(particle->P())),particle->Pt(),fWeightJetJetMC);

      if (neutPion->GetNDaughters()==2 && neutPion->GetDaughterFirst()>-1 && neutPion->GetDaughterLast()>-1){

        AliAODMCParticle* gamma0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(neutPion->GetDaughterFirst()));
        AliAODMCParticle* gamma1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(neutPion->GetDaughterLast()));
        AliAODMCParticle* gamma2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(labelGamma));

        // TParticle *gamma0 = (TParticle*)fMCEvent->Particle(neutPion->GetDaughterFirst());
        // TParticle *gamma1 = (TParticle*)fMCEvent->Particle(neutPion->GetDaughterLast());
        // TParticle *gamma2 = (TParticle*)fMCEvent->Particle(labelGamma);

        Bool_t InAcceptance = kFALSE;

        if (fReconMethod%2==1){ //cases 1,3,5 where reconstruction requires gamma2 to be a pcm photon
          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma2,fAODMCTrackArray,kFALSE)){ //check that gamma2 is in acceptance
            if(fReconMethod==5){
              //check that both gamma0 and gamma1 are in acceptance
              if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma0,fAODMCTrackArray,kFALSE) &&
                 ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma1,fAODMCTrackArray,kFALSE)){
                InAcceptance = kTRUE;
              }
            } else if(fReconMethod==3){
              if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma0,fAODMCTrackArray) &&
                 ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma1,fAODMCTrackArray)){
                InAcceptance = kTRUE;
              }
            } else if(fReconMethod==1){ // both gammas must be in tpc acceptance
              if((((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma0,fAODMCTrackArray,kFALSE) &&
                 ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma1,fAODMCTrackArray,kFALSE)) &&
                 // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                 (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma0,fAODMCTrackArray) ||
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma1,fAODMCTrackArray))){
                InAcceptance = kTRUE;
              }
            }
          }
        } else{ //cases 0,2,4 where reconstruction requires gamma2 to be a calo photon
          if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma2,fAODMCTrackArray)){
            if(fReconMethod==4){
              //check that both gamma0 and gamma1 are in acceptance
              if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma0,fAODMCTrackArray,kFALSE) &&
                 ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma1,fAODMCTrackArray,kFALSE)){
                InAcceptance = kTRUE;
              }
            } else if(fReconMethod==2){
              if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma0,fAODMCTrackArray) &&
                 ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma1,fAODMCTrackArray)){
                InAcceptance = kTRUE;
              }
            } else if(fReconMethod==0){ //either gamma0 is in tpc acc & gamma1 is in emcal acc or vice versa
              if((((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma0,fAODMCTrackArray,kFALSE) &&
                 ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(gamma1,fAODMCTrackArray,kFALSE)) &&
                 // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                 (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma0,fAODMCTrackArray) ||
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma1,fAODMCTrackArray))){
                InAcceptance = kTRUE;
              }
            }
          }
        }
        if(InAcceptance){
          // fill in acceptance histograms
          fHistoMCOmegaInAccInvMassPt[fiCut]->Fill(TMath::Sqrt((particle->E())*(particle->E())-(particle->P())*(particle->P())),particle->Pt(),fWeightJetJetMC);
          fHistoMCPi0FromOmegaInAccInvMassPt[fiCut]->Fill(TMath::Sqrt((neutPion->E())*(neutPion->E())-(neutPion->P())*(neutPion->P())),neutPion->Pt(),fWeightJetJetMC);
          fHistoMCGammaFromOmegaInAccPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);

          if(!fDoLightOutput){
            if(fDoMesonQA & 0b00000001) {
              Double_t OmegaY = 10.;
              if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
                OmegaY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                OmegaY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }
              Double_t Pi0Y = 10.;
              if(neutPion->E() - neutPion->Pz() == 0 || neutPion->E() + neutPion->Pz() == 0){
                Pi0Y=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }else{
                Pi0Y = 0.5*(TMath::Log((neutPion->E()+neutPion->Pz()) / (neutPion->E()-neutPion->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
              }

              fHistoMCOmegaInAccYPt[fiCut]->Fill(particle->Pt(),OmegaY,fWeightJetJetMC);
              fHistoMCPi0FromOmegaInAccYPt[fiCut]->Fill(neutPion->Pt(),Pi0Y,fWeightJetJetMC);
            }
            if(fDoMesonQA & 0b00000010) {
              Double_t OmegaAlpha = (neutPion->E() - gamma2->E())/(neutPion->E() + gamma2->E());
              Double_t Pi0Alpha = (gamma0->E() - gamma1->E())/(gamma0->E() + gamma1->E());
              fHistoMCOmegaInAccAlphaPt[fiCut]->Fill(particle->Pt(),OmegaAlpha,fWeightJetJetMC);
              fHistoMCPi0FromOmegaInAccAlphaPt[fiCut]->Fill(neutPion->Pt(),Pi0Alpha,fWeightJetJetMC);
            }
            if(fDoMesonQA & 0b00001000) {
              fHistoMCOmegaInAccEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
              fHistoMCPi0FromOmegaInAccEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
            }
            if(fDoMesonQA & 0b00010000) {
              TVector3 omegavector = TVector3(particle->Px(),particle->Py(),particle->Pz());
              TVector3 pi0vector = TVector3(neutPion->Px(),neutPion->Py(),neutPion->Pz());
              TVector3 gamma2vector = TVector3(gamma2->Px(),gamma2->Py(),gamma2->Pz());

              Double_t OmegaPiZeroAngle = TMath::Pi() - pi0vector.Angle(omegavector);
              Double_t PiZeroGammaAngle = pi0vector.Angle(gamma2vector);
              Double_t OmegaGammaAngle = TMath::Pi() - omegavector.Angle(gamma2vector);

              fHistoMCInAccOmegaPiZeroAnglePt[fiCut]->Fill(particle->Pt(),OmegaPiZeroAngle,fWeightJetJetMC);
              fHistoMCInAccPiZeroGammaAnglePt[fiCut]->Fill(particle->Pt(),PiZeroGammaAngle,fWeightJetJetMC);
              fHistoMCInAccOmegaGammaAnglePt[fiCut]->Fill(particle->Pt(),OmegaGammaAngle,fWeightJetJetMC);
            }
            if(fDoMesonQA & 0b100000000) {
              fHistoMCInAccOmegaPtPi0Pt[fiCut]->Fill(particle->Pt(),neutPion->Pt(),fWeightJetJetMC);
              fHistoMCInAccOmegaPtGammaPt[fiCut]->Fill(particle->Pt(),gamma2->Pt(),fWeightJetJetMC);
            }
          }
        }
      }
     }
    }

    Int_t isMCFromMBHeader = -1;
    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
      if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
    }

    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kFALSE)){
      fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
      if (TMath::Abs(particle->Eta()) < 0.66 ){
        if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
      }
    }
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,fAODMCTrackArray,kTRUE)){
      Double_t rConv = 0;
      for(Int_t daughterIndex=particle->GetDaughterLabel(0);daughterIndex<=particle->GetDaughterLabel(1);daughterIndex++){
        AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(daughterIndex));
        if(!tmpDaughter) continue;
        if(TMath::Abs(tmpDaughter->GetPdgCode()) == 11){
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
      ->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
      AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
      AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
      Float_t weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
        if (particle->Pt()>0.005){
          weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
          //                   if(particle->GetPdgCode() == 221){
          //                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
          //                   }
        }
      }
      Double_t mesonY = 1.e30;
      Double_t ratio  = 0;
      if (particle->E() != TMath::Abs(particle->Pz())){
        ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
      }
      if( !(ratio <= 0) ){
        mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
      }

      Double_t alpha = -10;
      if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
        alpha = (daughter0->E() - daughter1->E())/(daughter0->E() + daughter1->E());
      }


      if(particle->GetPdgCode() == 111){
        fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
        fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
        if (!fDoLightOutput){
        if (fDoMesonQA & 0b00000001) {
          fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC);
        }
        if (fDoMesonQA & 0b00000010) {
          fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);
        }
          if( (fDoMesonQA & 0b100000000) && (fIsMC == 2) ) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
        }
      }

      // Check the acceptance for both gammas
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,fAODMCTrackArray,kFALSE) &&
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,fAODMCTrackArray,kFALSE)  &&
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
      ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
        // check acceptance of clusters as well, true if one of them points into the Calo acceptance
        if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) ||
          ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray) ){
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
  for(Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      // fill primary histograms
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      // fill histograms for all true omegas
      if(particle->GetPdgCode() == 223){

        Bool_t DecaysToPiPlPiMiPiZero = kFALSE;

        TParticle *posPion = 0x0;
        TParticle *negPion = 0x0;
        TParticle *neutPion = 0x0;

        for(Int_t index= particle->GetFirstDaughter();index <= particle->GetLastDaughter();index++){
          if(index < 0) continue;

          TParticle *temp = (TParticle*)fMCEvent->Particle(index);

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
          fHistoMCOmegaDecayChannels[fiCut]->Fill(0.,fWeightJetJetMC);
        } else if (particle->GetNDaughters()==2){
          TParticle *gamma2 = 0x0;
          TParticle *pi0 = 0x0;

          for(Int_t index = particle->GetFirstDaughter();index <= particle->GetLastDaughter();index++){
            if(index < 0) continue;

            TParticle *temp = (TParticle*)fMCEvent->Particle(index);
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
            fHistoMCOmegaDecayChannels[fiCut]->Fill(1,fWeightJetJetMC);
            fHistoMCAllOmegaInvMassPt[fiCut]->Fill(TMath::Sqrt((particle->Energy())*(particle->Energy())-(particle->P())*(particle->P())),particle->Pt(),fWeightJetJetMC);
            fHistoMCGammaFromAllOmegaPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
            fHistoMCPi0FromAllOmegaInvMassPt[fiCut]->Fill(TMath::Sqrt((pi0->Energy())*(pi0->Energy())-(pi0->P())*(pi0->P())),pi0->Pt(),fWeightJetJetMC);
            if(!fDoLightOutput){
              if(fDoMesonQA & 0b00000001) {
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
              }
              if(fDoMesonQA & 0b00000010) {
                Double_t alpha = (pi0->Energy() - gamma2->Energy())/(pi0->Energy() + gamma2->Energy());
                fHistoMCAllOmegaAlphaPt[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);

                //check whether pi0 decayed into two gammas
                if (pi0->GetNDaughters()==2 && pi0->GetFirstDaughter()>-1 && pi0->GetLastDaughter()>-1){
                  TParticle *gamma0 = (TParticle*)fMCEvent->Particle(pi0->GetFirstDaughter());
                  TParticle *gamma1 = (TParticle*)fMCEvent->Particle(pi0->GetLastDaughter());
                  if (gamma0->GetPdgCode()==22 && gamma1->GetPdgCode()==22){

                    //plot pi0 alpha
                    Double_t pi0alpha = (gamma0->Energy() - gamma1->Energy())/(gamma0->Energy() + gamma1->Energy());
                    fHistoMCPi0FromAllOmegaAlphaPt[fiCut]->Fill(pi0->Pt(),pi0alpha,fWeightJetJetMC);
                  }
                }
              }
              if(fDoMesonQA & 0b00001000) {
                fHistoMCPi0FromAllOmegaEtaPhi[fiCut]->Fill(pi0->Phi(),pi0->Eta(),fWeightJetJetMC);
                fHistoMCAllOmegaEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
              }
              if(fDoMesonQA & 0b00010000) {
                //get angles and plot
                TVector3 omegavector = TVector3(particle->Px(),particle->Py(),particle->Pz());
                TVector3 pi0vector = TVector3(pi0->Px(),pi0->Py(),pi0->Pz());
                TVector3 gamma2vector = TVector3(gamma2->Px(),gamma2->Py(),gamma2->Pz());
                fHistoMCAllOmegaPiZeroAnglePt[fiCut]->Fill(particle->Pt(),TMath::Pi() - pi0vector.Angle(omegavector),fWeightJetJetMC);
                fHistoMCAllPiZeroGammaAnglePt[fiCut]->Fill(particle->Pt(),pi0vector.Angle(gamma2vector),fWeightJetJetMC);
                fHistoMCAllOmegaGammaAnglePt[fiCut]->Fill(particle->Pt(),TMath::Pi() - omegavector.Angle(gamma2vector),fWeightJetJetMC);
              }
              if(fDoMesonQA & 0b100000000) {
                fHistoMCAllOmegaPtPi0Pt[fiCut]->Fill(particle->Pt(),pi0->Pt(),fWeightJetJetMC);
                fHistoMCAllOmegaPtGammaPt[fiCut]->Fill(particle->Pt(),gamma2->Pt(),fWeightJetJetMC);
              }
            }
          } else{
            fHistoMCOmegaDecayChannels[fiCut]->Fill(2,fWeightJetJetMC);
          }
        } else{
          fHistoMCOmegaDecayChannels[fiCut]->Fill(2,fWeightJetJetMC);
        }
      }

      //fill histograms for omegas in acceptance
      Int_t labelNeutPion = -1;
      Int_t labelGamma = -1;

      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
         ->MesonIsSelectedMCPiZeroGamma(particle,fMCEvent,labelNeutPion,labelGamma,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
       if (labelNeutPion > -1 && labelGamma > -1){
        TParticle *neutPion    = fMCEvent->Particle(labelNeutPion);

        //fill histograms for acceptance correction
        fHistoMCOmegaInvMassPt[fiCut]->Fill(TMath::Sqrt((particle->Energy())*(particle->Energy())-(particle->P())*(particle->P())),particle->Pt(),fWeightJetJetMC);

        if (neutPion->GetNDaughters()==2 && neutPion->GetFirstDaughter()>-1 && neutPion->GetLastDaughter()>-1){

          TParticle *gamma0 = (TParticle*)fMCEvent->Particle(neutPion->GetFirstDaughter());
          TParticle *gamma1 = (TParticle*)fMCEvent->Particle(neutPion->GetLastDaughter());
          TParticle *gamma2 = (TParticle*)fMCEvent->Particle(labelGamma);

          Bool_t InAcceptance = kFALSE;

          if (fReconMethod%2==1){ //cases 1,3,5 where reconstruction requires gamma2 to be a pcm photon
            if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma2,fMCEvent,kFALSE)){ //check that gamma2 is in acceptance
              if(fReconMethod==5){
                //check that both gamma0 and gamma1 are in acceptance
                if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCEvent,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==3){
                if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCEvent) &&
                   ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==1){ // both gammas must be in tpc acceptance
                if((((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCEvent,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE)) &&
                   // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                   (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCEvent) ||
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent))){
                  InAcceptance = kTRUE;
                }
              }
            }
          } else{ //cases 0,2,4 where reconstruction requires gamma2 to be a calo photon
            if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma2,fMCEvent)){
              if(fReconMethod==4){
                //check that both gamma0 and gamma1 are in acceptance
                if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCEvent,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==2){
                if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCEvent) &&
                   ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent)){
                  InAcceptance = kTRUE;
                }
              } else if(fReconMethod==0){ //either gamma0 is in tpc acc & gamma1 is in emcal acc or vice versa
                if((((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma0,fMCEvent,kFALSE) &&
                   ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE)) &&
                   // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                   (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma0,fMCEvent) ||
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent))){
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

            if(!fDoLightOutput){
              if(fDoMesonQA & 0b00000001) {
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
                fHistoMCOmegaInAccYPt[fiCut]->Fill(particle->Pt(),OmegaY,fWeightJetJetMC);
                fHistoMCPi0FromOmegaInAccYPt[fiCut]->Fill(neutPion->Pt(),Pi0Y,fWeightJetJetMC);
              }
              if(fDoMesonQA & 0b00000010) {
                Double_t OmegaAlpha = (neutPion->Energy() - gamma2->Energy())/(neutPion->Energy() + gamma2->Energy());
                Double_t Pi0Alpha = (gamma0->Energy() - gamma1->Energy())/(gamma0->Energy() + gamma1->Energy());
                fHistoMCOmegaInAccAlphaPt[fiCut]->Fill(particle->Pt(),OmegaAlpha,fWeightJetJetMC);
                fHistoMCPi0FromOmegaInAccAlphaPt[fiCut]->Fill(neutPion->Pt(),Pi0Alpha,fWeightJetJetMC);
              }
              if(fDoMesonQA & 0b00001000) {
                fHistoMCOmegaInAccEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
                fHistoMCPi0FromOmegaInAccEtaPhi[fiCut]->Fill(particle->Phi(),particle->Eta(),fWeightJetJetMC);
              }
              if(fDoMesonQA & 0b00010000) {
                TVector3 omegavector = TVector3(particle->Px(),particle->Py(),particle->Pz());
                TVector3 pi0vector = TVector3(neutPion->Px(),neutPion->Py(),neutPion->Pz());
                TVector3 gamma2vector = TVector3(gamma2->Px(),gamma2->Py(),gamma2->Pz());

                Double_t OmegaPiZeroAngle = TMath::Pi() - pi0vector.Angle(omegavector);
                Double_t PiZeroGammaAngle = pi0vector.Angle(gamma2vector);
                Double_t OmegaGammaAngle = TMath::Pi() - omegavector.Angle(gamma2vector);

                fHistoMCInAccOmegaPiZeroAnglePt[fiCut]->Fill(particle->Pt(),OmegaPiZeroAngle,fWeightJetJetMC);
                fHistoMCInAccPiZeroGammaAnglePt[fiCut]->Fill(particle->Pt(),PiZeroGammaAngle,fWeightJetJetMC);
                fHistoMCInAccOmegaGammaAnglePt[fiCut]->Fill(particle->Pt(),OmegaGammaAngle,fWeightJetJetMC);
              }
              if(fDoMesonQA & 0b100000000) {
                fHistoMCInAccOmegaPtPi0Pt[fiCut]->Fill(particle->Pt(),neutPion->Pt(),fWeightJetJetMC);
                fHistoMCInAccOmegaPtGammaPt[fiCut]->Fill(particle->Pt(),gamma2->Pt(),fWeightJetJetMC);
              }
            }
          }
        }
       }
      }


      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
        fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
        if (TMath::Abs(particle->Eta()) < 0.66 ){
          if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        }
      }
      if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
        fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
        if (fDoPhotonQA > 0 && particle->GetFirstDaughter()>-1){
          fHistoMCConvGammaR[fiCut]->Fill(((TParticle*)fMCEvent->Particle(particle->GetFirstDaughter()))->R());
          fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta());
        }
      }// Converted MC Gamma
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
        ->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
       if(particle->GetFirstDaughter()>-1 && particle->GetLastDaughter()>-1){
        TParticle* daughter0 = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
        TParticle* daughter1 = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
            //                   if(particle->GetPdgCode() == 221){
            //                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
            //                   }
          }
        }

        Double_t mesonY = 1.e30;
        Double_t ratio  = 0;
        if (particle->Energy() != TMath::Abs(particle->Pz())){
          ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
        }
        if( !(ratio <= 0) ){
          mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }


        Double_t alpha = -10;
        if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
          alpha = (daughter0->Energy() - daughter1->Energy())/(daughter0->Energy() + daughter1->Energy());
        }

        if(particle->GetPdgCode() == 111){
          fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Pi0
          fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (!fDoLightOutput){
            if(fDoMesonQA & 0b00000001) fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC); // All MC Pi0
            if(fDoMesonQA & 0b00000010) fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Pi0
            if( (fDoMesonQA & 0b100000000) && (fIsMC == 2) ) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }
        }
        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t InAcceptance = kFALSE;

        if(fReconMethod/2 == 0){
          if( kDaughter0IsPrim && kDaughter1IsPrim &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)){
            // check acceptance of clusters as well, true if one of them points into the Calo acceptance
            if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) ||
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent)){
              InAcceptance = kTRUE;
            }
          }
        } else if(fReconMethod/2 == 1){
          if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
             ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent)){
            InAcceptance = kTRUE;
          }
        } else if(fReconMethod/2 == 2){
          if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
             ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)){
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
}
//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::CalculatePi0Candidates(){
  switch(fReconMethod){

    // PCM-cal,cal
    case 0:
    {
      if(fClusterCandidates->GetEntries()>1 && fGammaCandidates->GetEntries()>0){
        if (fPhotonSelectionMode >= 1){
          PhotonSelectionCalo();                                                // exclude clusters which would pair with other cluster to a Pi0 to be chosen as photon from the omega
        }
        for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
          AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
          if (gamma0==NULL) continue;
          for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++) {
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
            if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
            Bool_t matched = kFALSE;
            AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster,fInputEvent,fWeightJetJetMC);

            AliAODConversionMother pi0cand = AliAODConversionMother(gamma0,gamma1);

            // Pi0ROtation
            // ...

            if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                && pi0cand.Pt() > fMinPi0Pt) {
              if (matched){
                fHistoPhotonPairMatchedInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
              } else {
                fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
                if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
                  AliAODConversionMother *pi0cand_vec = new AliAODConversionMother(gamma0,gamma1);
                  pi0cand_vec->SetLabel1(firstGammaIndex);                      // safe the Gamma indices in this otherwise empty list
                  pi0cand_vec->SetLabel2(secondGammaIndex);                     // to later get the Gammas back for the MCTrue check
                  fPi0Candidates->Add(pi0cand_vec);
                  if (fPhotonSelectionMode >= 2) {
                    dropOutGammas_CALO.insert(secondGammaIndex);                // exclude clusters which would pair with PCM photon to a Pi0 candidate.
                  }
                  if(!fDoLightOutput){
                    if(fDoMesonQA & 0b00000001){
                      fHistoPhotonPairYPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                    }
                    if(fDoMesonQA & 0b00000010){
                      fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetAlpha(),fWeightJetJetMC);
                    }
                    if(fDoMesonQA & 0b00000100){
                      fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetOpeningAngle(),fWeightJetJetMC);
                    }
                    if(fDoMesonQA & 0b00001000){
                      fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand_vec->Phi(),pi0cand_vec->Eta(),fWeightJetJetMC);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    break;

    // pcm-cal,pcm
    case 1:
    {
      if( (fGammaCandidates->GetEntries() > 1) && (fClusterCandidates->GetEntries() > 0) ) {
        if(fPhotonSelectionMode >= 1) {
          PhotonSelectionPCM();                                                 // exclude PCM photon which would pair with other PCM photon to a Pi0 candidate to be chosen as photon from the omega
        }
        for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
          AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
          if (gamma0==NULL) continue;

          for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
            Bool_t matched = kFALSE;
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
            if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;

            AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC);

            AliAODConversionMother pi0cand = AliAODConversionMother(gamma0,gamma1);

            if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                && pi0cand.Pt() > fMinPi0Pt){
              if (matched){
                fHistoPhotonPairMatchedInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
              } else{
                // fill photon pair histograms
                fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
                if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
                  AliAODConversionMother *pi0cand_vec = new AliAODConversionMother(gamma0,gamma1);
                  pi0cand_vec->SetLabel1(firstGammaIndex);                      // safe the Gamma indices in this otherwise empty list
                  pi0cand_vec->SetLabel2(secondGammaIndex);                     // to later get the Gammas back for the MCTrue check
                  fPi0Candidates->Add(pi0cand_vec);
                  if (fPhotonSelectionMode >=2){
                    dropOutGammas_PCM.insert(firstGammaIndex);                  // exclude PCM photon which would pair with cluster to a Pi0 candidate to be chosen as photon from the omega
                  }
                  if(!fDoLightOutput){
                    if(fDoMesonQA & 0b00000001){
                      fHistoPhotonPairYPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                    }
                    if(fDoMesonQA & 0b00000010){
                      fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetAlpha(),fWeightJetJetMC);
                    }
                    if(fDoMesonQA & 0b00000100){
                      fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetOpeningAngle(),fWeightJetJetMC);
                    }
                    if(fDoMesonQA & 0b00001000){
                      fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand_vec->Phi(),pi0cand_vec->Eta(),fWeightJetJetMC);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    break;

    // cal-cal,cal
    case 2:
    {
      if(fClusterCandidates->GetEntries()>2){
        if ( (fPhotonSelectionMode >= 2) && (fGammaCandidates->GetEntries() > 0) ) {                                        // exclude clusters which would pair with PCM photon to a Pi0 candidate.
          PhotonSelectionMixed();
        }

        for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries()-1;firstGammaIndex++){
          AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
          if (gamma0==NULL || !(gamma0->GetIsCaloPhoton())) continue;
          for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
            if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
            AliAODConversionMother pi0cand = AliAODConversionMother(gamma0, gamma1);
            if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())
                && pi0cand.Pt() > fMinPi0Pt){
              fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
              if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
                AliAODConversionMother *pi0cand_vec = new AliAODConversionMother(gamma0, gamma1);
                pi0cand_vec->SetLabel1(firstGammaIndex);                        // safe the Gamma indices in this otherwise empty list
                pi0cand_vec->SetLabel2(secondGammaIndex);                       // to later get the Gammas back for the MCTrue check
                fPi0Candidates->Add(pi0cand_vec);
                if(fPhotonSelectionMode >= 1){
                  dropOutGammas_CALO.insert(firstGammaIndex);                   // exclude clusters which would pair with PCM photon to a Pi0 candidate
                  dropOutGammas_CALO.insert(secondGammaIndex);                  // from being used as photon coming from an omega
                }
                if(!fDoLightOutput){
                  if(fDoMesonQA & 0b00000001){
                    fHistoPhotonPairYPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000010){
                    fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetAlpha(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000100){
                    fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetOpeningAngle(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00001000){
                    fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand_vec->Phi(),pi0cand_vec->Eta(),fWeightJetJetMC);
                  }
                }
              }
            }
          }
        }
      }
    }
    break;

    // cal-cal,pcm
    case 3:
    {
      if(fClusterCandidates->GetEntries()>1 && fGammaCandidates->GetEntries()>0){ // atleast three photons are needed for one omega!
        if (fPhotonSelectionMode >= 2) {                                        // more strict photon selection, means all possible gamma comibnations are achecked
            PhotonSelectionMixed();                                             // exclude PCM photons which would pair with cluster to a Pi0 candidate from being used as photon coming from an omega
        }
        for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries()-1;firstGammaIndex++){
          AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
          if (gamma0==NULL || !(gamma0->GetIsCaloPhoton())) continue;
          for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
            if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
            AliAODConversionMother pi0cand = AliAODConversionMother(gamma0, gamma1);
            if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())
                && pi0cand.Pt() > fMinPi0Pt){
              fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
              if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
                AliAODConversionMother *pi0cand_vec = new AliAODConversionMother(gamma0, gamma1);
                pi0cand_vec->SetLabel1(firstGammaIndex);                        // safe the Gamma indices in this otherwise empty list
                pi0cand_vec->SetLabel2(secondGammaIndex);                       // to later get the Gammas back for the MCTrue check
                fPi0Candidates->Add(pi0cand_vec);
                if(fPhotonSelectionMode >= 1){
                  PhotonSelectionPCM();                                         // exclude PCM photons which would pair with other PCM photon to a Pi0 candidate from being used as photon coming from an omega
                }
                if(!fDoLightOutput){
                  if(fDoMesonQA & 0b00000001){
                    fHistoPhotonPairYPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000010){
                    fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetAlpha(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000100){
                    fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetOpeningAngle(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00001000){
                    fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand_vec->Phi(),pi0cand_vec->Eta(),fWeightJetJetMC);
                  }
                }
              }
            }
          }
        }
      }
    }
    break;

    // pcm-pcm,cal
    case 4:
    {
      if(fClusterCandidates->GetEntries()>1 && fGammaCandidates->GetEntries()>0){ // atleast three photons are needed for one omega!
        if (fPhotonSelectionMode >= 1){
          PhotonSelectionCalo();                                                // exclude clusters which would pair with PCM photon to a Pi0 candidate from being used as photon coming from an omega
        }
        if (fPhotonSelectionMode >= 2) {                                        // more strict photon selection, means all possible gamma comibnations are achecked
            PhotonSelectionMixed();                                             // exclude PCM photons which would pair with cluster to a Pi0 candidate from being used as photon coming from an omega
        }
      }
      if(fGammaCandidates->GetEntries()>1 && fClusterCandidates->GetEntries()>0){
        for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++)
        {
          AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
          if(gamma0==NULL) continue;
          for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++)
          {
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
            if(gamma1==NULL) continue;
            AliAODConversionMother pi0cand = AliAODConversionMother(gamma0,gamma1);
            if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()) && pi0cand.Pt() > fMinPi0Pt)
            {
              fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
              if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) )
              {
                AliAODConversionMother *pi0cand_vec = new AliAODConversionMother(gamma0,gamma1);
                pi0cand_vec->SetLabel1(firstGammaIndex);                        // safe the Gamma indices in this otherwise empty list
                pi0cand_vec->SetLabel2(secondGammaIndex);                       // to later get the Gammas back for the MCTrue check
                fPi0Candidates->Add(pi0cand_vec);
                if(!fDoLightOutput)
                {
                  if(fDoMesonQA & 0b00000001){
                    fHistoPhotonPairYPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000010){
                    fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetAlpha(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000100){
                    fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetOpeningAngle(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00001000){
                    fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand_vec->Phi(),pi0cand_vec->Eta(),fWeightJetJetMC);
                  }
                }
              }
            }
          }
        }
      }
    }
    break;

    // pcm-pcm,pcm
    case 5:
    {
      if(fGammaCandidates->GetEntries() > 2){
        if ( (fPhotonSelectionMode >= 2) && (fClusterCandidates->GetEntries() > 0) ){
          PhotonSelectionMixed();                                               // exclude PCM photons which would pair with cluster to a Pi0 candidate from being used as photon coming from an omega
        }
        for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
          AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
          if (gamma0==NULL) continue;
          for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
            if (gamma1==NULL || (secondGammaIndex+1)>=(fGammaCandidates->GetEntries())) continue;
            AliAODConversionMother pi0cand = AliAODConversionMother(gamma0, gamma1);
            if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
                 ->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())
                && pi0cand.Pt() > fMinPi0Pt){
              fHistoPhotonPairInvMassPt[fiCut]->Fill(pi0cand.M(),pi0cand.Pt(),fWeightJetJetMC);
              if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
                AliAODConversionMother *pi0cand_vec = new AliAODConversionMother(gamma0, gamma1);
                pi0cand_vec->SetLabel1(firstGammaIndex);                        // safe the Gamma indices in this otherwise empty list
                pi0cand_vec->SetLabel2(secondGammaIndex);                       // to later get the Gammas back for the MCTrue check
                fPi0Candidates->Add(pi0cand_vec);
                if(fPhotonSelectionMode >= 1){
                  dropOutGammas_PCM.insert(firstGammaIndex);                    // exclude PCM photons which would pair with other PCM photon to a Pi0 candidate
                  dropOutGammas_PCM.insert(secondGammaIndex);                   // from being used as photon coming from an omega
                }
                if(!fDoLightOutput){
                  if(fDoMesonQA & 0b00000001){
                    fHistoPhotonPairYPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000010){
                    fHistoPhotonPairAlphaPt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetAlpha(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00000100){
                    fHistoPhotonPairOpenAnglePt[fiCut]->Fill(pi0cand_vec->Pt(),pi0cand_vec->GetOpeningAngle(),fWeightJetJetMC);
                  }
                  if(fDoMesonQA & 0b00001000){
                    fHistoPhotonPairEtaPhi[fiCut]->Fill(pi0cand_vec->Phi(),pi0cand_vec->Eta(),fWeightJetJetMC);
                  }
                }
              }
            }
          }
        }
      }
    }
    break;
  }
}
//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::CalculateOmegaCandidates()
{
  switch(fReconMethod){
  //PCM-cal,cal
  case 0:
    if( (fPi0Candidates->GetEntries() > 0) && (fClusterCandidates->GetEntries() > 0) ) {
      for(Int_t firstPi0Index=0;firstPi0Index<fPi0Candidates->GetEntries();firstPi0Index++){
        AliAODConversionMother *pi0cand_vec=dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(firstPi0Index));
        if (pi0cand_vec==NULL) continue;
        //change energy of pi0 candidate s.t. its mass is the pdg mass
        pi0cand_vec->SetPxPyPzE(pi0cand_vec->Px(),pi0cand_vec->Py(),pi0cand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0cand_vec->P()*pi0cand_vec->P()));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(0)));
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(1)));
        // get third photon from clusters and calculate inv mass of omega
        for(Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALO.find(thirdGammaIndex) != dropOutGammas_CALO.end() ) ) {continue; }
          Bool_t matchedgamma2wconvgamma = kFALSE;
          AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
          if (gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
          AliVCluster* cluster2 = fInputEvent->GetCaloCluster(gamma2->GetCaloClusterRef());
          matchedgamma2wconvgamma = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster2, fInputEvent, fWeightJetJetMC);
          AliAODConversionMother omegacand = AliAODConversionMother(pi0cand_vec,gamma2);
          if (matchedgamma2wconvgamma){
            fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
            continue;
          }
          if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&omegacand, pi0cand_vec, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) ) {
            fHistoRecoArmenterosPodolanskiPlotwoCut[fiCut]->Fill(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec));
            if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec) ) )
              && (omegacand.M() < 1.6) )
            {
              if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&omegacand, gamma0, gamma1, gamma2) )
              {
                Double_t dalitzcand01 = AliAODConversionMother(gamma0,gamma1).M2();
                Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();
                if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                {
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                }
                if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                {
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
                  fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                  if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) == 1 ) {
                    CalculateOmegaRotationBackground(thirdGammaIndex, firstPi0Index);
                  }
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                    ProcessTrueMesonCandidates(&omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                    ProcessTrueMesonCandidatesAOD(&omegacand,gamma0,gamma1,gamma2);
                  }
                  if(!fDoLightOutput){
                    FillQAPlots(&omegacand, pi0cand_vec, gamma0, gamma1, gamma2);
                  }
                }
              }
            }
          } else if(fDoPiZeroGammaAngleCut){
            fHistoMotherAngleCutRejectedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    }
    break;
  //PCM-cal,PCM
  case 1:
    if( (fGammaCandidates->GetEntries() > 0) && (fPi0Candidates->GetEntries() > 0) ) {
      for(Int_t firstPi0Index=0;firstPi0Index<fPi0Candidates->GetEntries();firstPi0Index++){
        AliAODConversionMother *pi0cand_vec=dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(firstPi0Index));
        if (pi0cand_vec==NULL) continue;
        //change energy of pi0 candidate s.t. its mass is the pdg mass
        pi0cand_vec->SetPxPyPzE(pi0cand_vec->Px(),pi0cand_vec->Py(),pi0cand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0cand_vec->P()*pi0cand_vec->P()));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(0)));
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(1)));
        AliVCluster* cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
        // get third photon from conversion photon candidates and calculate inv mass of omega
        for(Int_t thirdGammaIndex=0;thirdGammaIndex<fGammaCandidates->GetEntries();thirdGammaIndex++){
          if ( thirdGammaIndex==pi0cand_vec->GetLabel(0) ) continue;
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_PCM.find(thirdGammaIndex) != dropOutGammas_PCM.end() ) ) {continue; }
          Bool_t matchedgamma1wconvgamma2 = kFALSE;
          AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(thirdGammaIndex));
          if (gamma2==NULL) continue;
          matchedgamma1wconvgamma2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2, cluster1, fInputEvent, fWeightJetJetMC);
          AliAODConversionMother omegacand = AliAODConversionMother(pi0cand_vec,gamma2);
          if (matchedgamma1wconvgamma2){
            fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
            continue;
          }
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&omegacand, pi0cand_vec, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            if( (!fDoLightOutput) && (fDoMesonQA & 0b01000000)) fHistoRecoArmenterosPodolanskiPlotwoCut[fiCut]->Fill(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec));
            if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec) ) )
              && (omegacand.M() < 1.6) )
            {
              if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&omegacand, gamma0, gamma1, gamma2) )
              {
                Double_t dalitzcand01 = AliAODConversionMother(gamma0,gamma1).M2();
                Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();
                if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                {
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                }
                if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                {
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
                  fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                  if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) == 1 ) {
                    CalculateOmegaRotationBackground(thirdGammaIndex, firstPi0Index);
                  }
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                      ProcessTrueMesonCandidates(&omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                      ProcessTrueMesonCandidatesAOD(&omegacand,gamma0,gamma1,gamma2);
                  }
                  if(!fDoLightOutput){
                    FillQAPlots(&omegacand, pi0cand_vec, gamma0, gamma1, gamma2);
                  }
                }
              }
            }
          } else if(fDoPiZeroGammaAngleCut){
              fHistoMotherAngleCutRejectedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    }
    break;
  //cal-cal,cal
  case 2:
    if ( (fClusterCandidates->GetEntries()>2) && (fPi0Candidates->GetEntries() > 0) ){
      for(Int_t firstPi0Index=0;firstPi0Index<fPi0Candidates->GetEntries();firstPi0Index++){
        AliAODConversionMother *pi0cand_vec=dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(firstPi0Index));
        if (pi0cand_vec==NULL) continue;
        //change energy of pi0 candidate s.t. its mass is the pdg mass
        pi0cand_vec->SetPxPyPzE(pi0cand_vec->Px(),pi0cand_vec->Py(),pi0cand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0cand_vec->P()*pi0cand_vec->P()));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(0)));
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(1)));
        for(Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
          if (thirdGammaIndex==pi0cand_vec->GetLabel(0) || thirdGammaIndex==pi0cand_vec->GetLabel(1)) continue;
          if ( (((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALO.find(thirdGammaIndex) != dropOutGammas_CALO.end() ) ) {continue;}
          AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
          if (gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
          AliAODConversionMother omegacand = AliAODConversionMother(pi0cand_vec,gamma2);
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&omegacand, pi0cand_vec, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            if( (!fDoLightOutput) && (fDoMesonQA & 0b01000000)) fHistoRecoArmenterosPodolanskiPlotwoCut[fiCut]->Fill(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec));
            if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec) ) )
              && (omegacand.M() < 1.6) )
            {
              if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&omegacand, gamma0, gamma1, gamma2) )
              {
                Double_t dalitzcand01 = AliAODConversionMother(gamma0,gamma1).M2();
                Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();
                if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                {
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                }
                if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                {
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
                  fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                  if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) == 1 ) {
                    CalculateOmegaRotationBackground(thirdGammaIndex, firstPi0Index);
                  }
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                      ProcessTrueMesonCandidates(&omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                      ProcessTrueMesonCandidatesAOD(&omegacand,gamma0,gamma1,gamma2);
                  }
                  if(!fDoLightOutput){
                    FillQAPlots(&omegacand, pi0cand_vec, gamma0, gamma1, gamma2);
                  }
                }
              }
            }
          } else if(fDoPiZeroGammaAngleCut){
              fHistoMotherAngleCutRejectedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    }
    break;
  //cal-cal,PCM
  case 3:

    if( (fPi0Candidates->GetEntries() > 0) && (fGammaCandidates->GetEntries() > 0) ){
      for(Int_t firstPi0Index=0;firstPi0Index<fPi0Candidates->GetEntries();firstPi0Index++){
        AliAODConversionMother *pi0cand_vec=dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(firstPi0Index));
        if (pi0cand_vec==NULL) continue;
        //change energy of pi0 candidate s.t. its mass is the pdg mass
        pi0cand_vec->SetPxPyPzE(pi0cand_vec->Px(),pi0cand_vec->Py(),pi0cand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0cand_vec->P()*pi0cand_vec->P()));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(0)));
        AliVCluster* cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(1)));
        AliVCluster* cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
        for(Int_t thirdGammaIndex=0;thirdGammaIndex<fGammaCandidates->GetEntries();thirdGammaIndex++){
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_PCM.find(thirdGammaIndex) != dropOutGammas_PCM.end() ) ) continue;
          AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(thirdGammaIndex));
          Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2,cluster0, fInputEvent, fWeightJetJetMC))
              || (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma2,cluster1, fInputEvent, fWeightJetJetMC));
          AliAODConversionMother omegacand = AliAODConversionMother(pi0cand_vec,gamma2);
          if (matched){
            fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
            continue;
          }
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&omegacand, pi0cand_vec, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            if( (!fDoLightOutput) && (fDoMesonQA & 0b01000000)) fHistoRecoArmenterosPodolanskiPlotwoCut[fiCut]->Fill(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec));
            if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec) ) )
              && (omegacand.M() < 1.6) )
            {
              if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&omegacand, gamma0, gamma1, gamma2) )
              {
                Double_t dalitzcand01 = AliAODConversionMother(gamma0,gamma1).M2();
                Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();
                if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                {
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                }
                if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                {
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
                  fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                  if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) == 1 ) {
                    CalculateOmegaRotationBackground(thirdGammaIndex, firstPi0Index);
                  }
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                      ProcessTrueMesonCandidates(&omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                      ProcessTrueMesonCandidatesAOD(&omegacand,gamma0,gamma1,gamma2);
                  }
                  if(!fDoLightOutput){
                    FillQAPlots(&omegacand, pi0cand_vec, gamma0, gamma1, gamma2);
                  }
                }
              }
            }
          } else if(fDoPiZeroGammaAngleCut){
              fHistoMotherAngleCutRejectedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    }
    break;
  //PCM-PCM,cal
  case 4:
    if(fGammaCandidates->GetEntries()>1 && fClusterCandidates->GetEntries()>0){
      for(Int_t firstPi0Index=0;firstPi0Index<fPi0Candidates->GetEntries();firstPi0Index++){
        AliAODConversionMother *pi0cand_vec=dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(firstPi0Index));
        if (pi0cand_vec==NULL) continue;
        //change energy of pi0 candidate s.t. its mass is the pdg mass
        pi0cand_vec->SetPxPyPzE(pi0cand_vec->Px(),pi0cand_vec->Py(),pi0cand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0cand_vec->P()*pi0cand_vec->P()));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(0)));
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(1)));
        for(Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
          AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALO.find(thirdGammaIndex) != dropOutGammas_CALO.end() ) ) {continue; }
          if(gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
          AliVCluster* cluster2 = fInputEvent->GetCaloCluster(gamma2->GetCaloClusterRef());
          Bool_t matched = (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster2, fInputEvent, fWeightJetJetMC))
              || (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma1,cluster2, fInputEvent, fWeightJetJetMC));
          AliAODConversionMother omegacand = AliAODConversionMother(pi0cand_vec,gamma2);
          if (matched){
            fHistoMotherMatchedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
            continue;
          }
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&omegacand, pi0cand_vec, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            if( (!fDoLightOutput) && (fDoMesonQA & 0b01000000)) fHistoRecoArmenterosPodolanskiPlotwoCut[fiCut]->Fill(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec));
            if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec) ) )
              && (omegacand.M() < 1.6) )
            {
              if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&omegacand, gamma0, gamma1, gamma2) )
              {
                Double_t dalitzcand01 = AliAODConversionMother(gamma0,gamma1).M2();
                Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();
                if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                {
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                }
                if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                {
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
                  fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                  if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) == 1 ) {
                    CalculateOmegaRotationBackground(thirdGammaIndex, firstPi0Index);
                  }
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                      ProcessTrueMesonCandidates(&omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                      ProcessTrueMesonCandidatesAOD(&omegacand,gamma0,gamma1,gamma2);
                  }
                  if(!fDoLightOutput){
                    FillQAPlots(&omegacand, pi0cand_vec, gamma0, gamma1, gamma2);
                  }
                }
              }
            }
          } else if(fDoPiZeroGammaAngleCut){
              fHistoMotherAngleCutRejectedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    }
    break;
  //PCM-PCM,PCM
  case 5:
    if(fGammaCandidates->GetEntries()>2){
      for(Int_t firstPi0Index=0;firstPi0Index<fPi0Candidates->GetEntries();firstPi0Index++){
        AliAODConversionMother *pi0cand_vec=dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(firstPi0Index));
        if (pi0cand_vec==NULL) continue;
        //change energy of pi0 candidate s.t. its mass is the pdg mass
        pi0cand_vec->SetPxPyPzE(pi0cand_vec->Px(),pi0cand_vec->Py(),pi0cand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0cand_vec->P()*pi0cand_vec->P()));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(0)));
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(pi0cand_vec->GetLabel(1)));
        for(Int_t thirdGammaIndex=0;thirdGammaIndex<fGammaCandidates->GetEntries();thirdGammaIndex++){
          if (thirdGammaIndex==pi0cand_vec->GetLabel(0) || thirdGammaIndex==pi0cand_vec->GetLabel(1)) continue;
          AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(thirdGammaIndex));
          if (gamma2==NULL) continue;
          AliAODConversionMother omegacand = AliAODConversionMother(pi0cand_vec,gamma2);
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&omegacand, pi0cand_vec, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            if( (!fDoLightOutput) && (fDoMesonQA & 0b01000000)) fHistoRecoArmenterosPodolanskiPlotwoCut[fiCut]->Fill(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec));
            if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&omegacand, pi0cand_vec, gamma2), GetQTPi0(&omegacand, pi0cand_vec) ) )
              && (omegacand.M() < 1.6) )
            {
              if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&omegacand, gamma0, gamma1, gamma2) )
              {
                Double_t dalitzcand01 = AliAODConversionMother(gamma0,gamma1).M2();
                Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();
                if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                {
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                  fHistoMotherDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                }
                if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                {
                  fHistoMotherInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
                  fHistoGammaFromMotherPt[fiCut]->Fill(gamma2->Pt(),fWeightJetJetMC);
                  if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) == 1 ) {
                    CalculateOmegaRotationBackground(thirdGammaIndex, firstPi0Index);
                  }
                  if(fIsMC>0){
                    if(fInputEvent->IsA()==AliESDEvent::Class())
                      ProcessTrueMesonCandidates(&omegacand,gamma0,gamma1,gamma2);
                    if(fInputEvent->IsA()==AliAODEvent::Class())
                      ProcessTrueMesonCandidatesAOD(&omegacand,gamma0,gamma1,gamma2);
                  }
                  if(!fDoLightOutput){
                    FillQAPlots(&omegacand, pi0cand_vec, gamma0, gamma1, gamma2);
                  }
                }
              }
            }
          } else if(fDoPiZeroGammaAngleCut){
              fHistoMotherAngleCutRejectedInvMassPt[fiCut]->Fill(omegacand.M(),omegacand.Pt(),fWeightJetJetMC);
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
  switch (fReconMethod){
    // pcm-cal,cal
    case 0:
    {
      // get gamma0MotherLabel
      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
        if(gamma0MCLabel>-1){
          TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
          gamma0MotherLabel=gammaMC0->GetFirstMother();
        }
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;
      TParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        // get mother or grandmother of gamma1 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=(fMCEvent->Particle(gammaMC1->GetMother(0)))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 111){

        // fill pi0 histograms here if necessary

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;
        TParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          // get mother or grandmother of gamma1 (potentially true omega) depending on whether it is an electron-leading/photon-leading cluster
          gammaMC2 = (TParticle*)fMCEvent->Particle(gamma2MCLabel);
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){  // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother(0);
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion() &&  gammaMC2->GetMother(0)>-1) gamma2MotherLabel=fMCEvent->Particle(gammaMC2->GetMother(0))->GetMother(0);
            else gamma2MotherLabel=gammaMC2->GetMother(0);
          }
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(pi0MotherLabel>=0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCEvent->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          };
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
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
        if(gamma0MCLabel>-1){
          TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
          gamma0MotherLabel=gammaMC0->GetFirstMother();
        }
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;
      TParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        // get mother or grandmother of gamma1 depending on whether it is an electron-leading/photon-leading cluster
        gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){

        // fill pi0 histograms here if necessary

        // get mother of gamma2
        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = TrueGammaCandidate2->GetMCParticleLabel(fMCEvent);
          if(gamma2MCLabel>-1){
            TParticle * gammaMC2 = (TParticle*)fMCEvent->Particle(gamma2MCLabel);
            gamma2MotherLabel = gammaMC2->GetFirstMother();
          }
        }

        // get pi0MotherLabel
        Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(pi0MotherLabel>=0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCEvent->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
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
        gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother(0);
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion() && gammaMC0->GetMother(0)>-1) gamma0MotherLabel=fMCEvent->Particle(gammaMC0->GetMother(0))->GetMother(0);
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
        gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){

        // fill pi0 histograms here if necessary

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;
        TParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          // get mother or grandmother of gamma2 (potentially true omega) depending on whether it is an electron-leading/photon-leading cluster
          gammaMC2 = (TParticle*)fMCEvent->Particle(gamma2MCLabel);
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){  // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother(0);
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion() && gammaMC2->GetMother(0)>-1) gamma2MotherLabel=fMCEvent->Particle(gammaMC2->GetMother(0))->GetMother(0);
            else gamma2MotherLabel=gammaMC2->GetMother(0);
          }
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(pi0MotherLabel>=0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCEvent->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
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
        gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother(0);
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion() && gammaMC0->GetMother(0)>-1) gamma0MotherLabel=fMCEvent->Particle(gammaMC0->GetMother(0))->GetMother(0);
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
        gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){  // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel && ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){

        // fill pi0 histograms here if necessary

        // get mother of gamma2
        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = TrueGammaCandidate2->GetMCParticleLabel(fMCEvent);
          if(gamma2MCLabel>-1){
            TParticle * gammaMC2 = (TParticle*)fMCEvent->Particle(gamma2MCLabel);
            gamma2MotherLabel = gammaMC2->GetFirstMother();
          }
        }

        // get pi0MotherLabel
        Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(pi0MotherLabel>=0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCEvent->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
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
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
        if(gamma0MCLabel>-1){
          TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
          gamma0MotherLabel=gammaMC0->GetFirstMother();
        }
      }

      // get gamma1MotherLabel
      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(fMCEvent);
        if(gamma1MCLabel>-1){
          TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
          gamma1MotherLabel=gammaMC1->GetFirstMother();
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){

        // fill pi0 histograms here if necessary

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;
        TParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          // get mother or grandmother of gamma1 (potentially true omega) depending on whether it is an electron-leading/photon-leading cluster
          gammaMC2 = (TParticle*)fMCEvent->Particle(gamma2MCLabel);
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){  // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother(0);
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){  // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion() && gammaMC2->GetMother(0)>-1) gamma2MotherLabel=fMCEvent->Particle(gammaMC2->GetMother(0))->GetMother(0);
            else gamma2MotherLabel=gammaMC2->GetMother(0);
          }
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(pi0MotherLabel>=0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCEvent->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
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
        gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
        if(gamma0MCLabel>-1){
          TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
          gamma0MotherLabel=gammaMC0->GetFirstMother();
        }
      }

      // get gamma1MotherLabel
      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(fMCEvent);
        if(gamma1MCLabel>-1){
          TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
          gamma1MotherLabel=gammaMC1->GetFirstMother();
        }
      }

      // check if mother of gamma0 and gamma1 is really a pi0
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){

        // fill pi0 histograms here if necessary

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = TrueGammaCandidate2->GetMCParticleLabel(fMCEvent);
          if(gamma2MCLabel>-1){
            TParticle * gammaMC2 = (TParticle*)fMCEvent->Particle(gamma2MCLabel);
            gamma2MotherLabel=gammaMC2->GetFirstMother();
          }
        }

        //get pi0MotherLabel
        Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);
        // check if mother of pi0 and mother of gamma2 are the same particle and that it is an omega
        if(pi0MotherLabel>=0 && gamma2MotherLabel==pi0MotherLabel && ((TParticle*)fMCEvent->Particle(pi0MotherLabel))->GetPdgCode() == 223){
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          // create pi0 candidate and fill histograms
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
        }
      }
    }
      break;
  }
}

//______________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *OmegaCandidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, AliAODConversionPhoton *TrueGammaCandidate2)
{
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;
  Bool_t isTrueOmega = kFALSE;

  switch(fReconMethod){
    // pcm-cal,cal
    case 0:
    {
      AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positiveMC||!negativeMC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positiveMC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      Int_t tmpGammaMotherlabel = gamma0MotherLabel;
      while (tmpGammaMotherlabel > 0) {
        if(((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetPdgCode() != 111)
        {
          tmpGammaMotherlabel = ((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetMother(0);
        }
        else
        {
          gamma0MotherLabel = tmpGammaMotherlabel;
          break;
        }
      }
      tmpGammaMotherlabel = gamma1MotherLabel;
      while (tmpGammaMotherlabel > 0) {
        if(((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetPdgCode() != 111)
        {
          tmpGammaMotherlabel = ((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetMother(0);
        }
        else
        {
          gamma1MotherLabel = tmpGammaMotherlabel;
          break;
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0ForMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0ForMother->GetMother();
        // Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;

        AliAODMCParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          gammaMC2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma2MCLabel));
          // get mother label
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){ // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother();
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()){
              AliAODMCParticle * gammaGrandMotherMC2 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC2->GetMother()));
              gamma2MotherLabel=gammaGrandMotherMC2->GetMother();
            }else gamma2MotherLabel=gammaMC2->GetMother();
          }
        }

        tmpGammaMotherlabel = gamma2MotherLabel;
        while (tmpGammaMotherlabel > 0) {
          if(((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetPdgCode() != 223)
          {
            tmpGammaMotherlabel = ((TParticle*)fMCEvent->Particle(tmpGammaMotherlabel))->GetMother(0);
          }
          else
          {
            gamma2MotherLabel = tmpGammaMotherlabel;
            break;
          }
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)fAODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
        }
      }
    }
      break;
    // pcm-cal,pcm
    case 1:
    {
      AliAODMCParticle *positive0MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negative0MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positive0MC||!negative0MC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positive0MC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0ForMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0ForMother->GetMother();
        // Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);

        // get gamma2MotherLabel
        AliAODMCParticle *positive2MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelPositive()));
        AliAODMCParticle *negative2MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelNegative()));

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if(!positive2MC||!negative2MC)
          return;

        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = positive2MC->GetMother();
          AliAODMCParticle * gammaMC2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma2MCLabel));
          gamma2MotherLabel=gammaMC2->GetMother();
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)fAODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
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
        gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        // get gamma0MotherLabel
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother();
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC0->GetMother()));
            gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
          }else gamma0MotherLabel=gammaMC0->GetMother();
        }
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0ForMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0ForMother->GetMother();
        // Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;
        AliAODMCParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          gammaMC2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma2MCLabel));
          // get gamma2MotherLabel
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){ // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother();
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()){
              AliAODMCParticle * gammaGrandMotherMC2 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC2->GetMother()));
              gamma2MotherLabel=gammaGrandMotherMC2->GetMother();
            }else gamma2MotherLabel=gammaMC2->GetMother();
          }
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)fAODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
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
        gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        // get gamma0MotherLabel
        if (TrueGammaCandidate0->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma0MotherLabel=gammaMC0->GetMother();
        }else if (TrueGammaCandidate0->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate0->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC0->GetMother()));
            gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
          }else gamma0MotherLabel=gammaMC0->GetMother();
        }
      }

      if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
      Int_t gamma1MotherLabel = -1;

      AliAODMCParticle * gammaMC1 = 0x0;
      if(gamma1MCLabel != -1){
        gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
        // get gamma1MotherLabel
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){ // for photons it's the direct mother
          gamma1MotherLabel=gammaMC1->GetMother();
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion()){
            AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC1->GetMother()));
            gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
          }else gamma1MotherLabel=gammaMC1->GetMother();
        }
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0ForMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0ForMother->GetMother();
        // Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);

        // get gamma2MotherLabel
        AliAODMCParticle *positive2MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelPositive()));
        AliAODMCParticle *negative2MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelNegative()));

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if(!positive2MC||!negative2MC)
          return;

        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = positive2MC->GetMother();
          AliAODMCParticle * gammaMC2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma2MCLabel));
          gamma2MotherLabel=gammaMC2->GetMother();
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)fAODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
        }
      }
    }
      break;
    // pcm-pcm,cal
    case 4:
    {
      AliAODMCParticle *positive0MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negative0MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positive0MC||!negative0MC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positive0MC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      AliAODMCParticle *positive1MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
      AliAODMCParticle *negative1MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));

      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if(!positive1MC||!negative1MC)
        return;

      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = positive1MC->GetMother();
        AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
        gamma1MotherLabel=gammaMC1->GetMother();
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0ForMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0ForMother->GetMother();
        // Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);

        if (!TrueGammaCandidate2->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
        Int_t gamma2MCLabel = TrueGammaCandidate2->GetCaloPhotonMCLabel(0); // get most probable MC label
        Int_t gamma2MotherLabel = -1;

        AliAODMCParticle * gammaMC2 = 0x0;
        if(gamma2MCLabel != -1){
          gammaMC2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma2MCLabel));
          // get mother label
          if (TrueGammaCandidate2->IsLargestComponentPhoton()){ // for photons it's the direct mother
            gamma2MotherLabel=gammaMC2->GetMother();
          }else if (TrueGammaCandidate2->IsLargestComponentElectron()){ // for electrons it's either the direct mother or for conversions the grandmother
            if (TrueGammaCandidate2->IsConversion()){
              AliAODMCParticle * gammaGrandMotherMC2 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC2->GetMother()));
              gamma2MotherLabel=gammaGrandMotherMC2->GetMother();
            }else gamma2MotherLabel=gammaMC2->GetMother();
          }
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)fAODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
        }
      }
    }
      break;
    // pcm-pcm,pcm
    case 5:
    {
      AliAODMCParticle *positive0MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
      AliAODMCParticle *negative0MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

      Int_t gamma0MCLabel = -1;
      Int_t gamma0MotherLabel = -1;
      if(!positive0MC||!negative0MC)
        return;

      if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
        gamma0MCLabel = positive0MC->GetMother();
        AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
        gamma0MotherLabel=gammaMC0->GetMother();
      }

      AliAODMCParticle *positive1MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
      AliAODMCParticle *negative1MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));

      Int_t gamma1MCLabel = -1;
      Int_t gamma1MotherLabel = -1;
      if(!positive1MC||!negative1MC)
        return;

      if (TrueGammaCandidate1->IsTrueConvertedPhoton()){
        gamma1MCLabel = positive1MC->GetMother();
        AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
        gamma1MotherLabel=gammaMC1->GetMother();
      }

      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel &&
         ((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){

        // create pi0 and get pi0MotherLabel
        AliAODMCParticle * TruePi0ForMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel));
        Int_t pi0MotherLabel = TruePi0ForMother->GetMother();
        // Int_t pi0MotherLabel = fMCEvent->Particle(gamma0MotherLabel)->GetMother(0);

        AliAODMCParticle *positive2MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelPositive()));
        AliAODMCParticle *negative2MC = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(TrueGammaCandidate2->GetMCLabelNegative()));

        Int_t gamma2MCLabel = -1;
        Int_t gamma2MotherLabel = -1;
        if(!positive2MC||!negative2MC)
          return;

        if (TrueGammaCandidate2->IsTrueConvertedPhoton()){
          gamma2MCLabel = positive2MC->GetMother();
          AliAODMCParticle * gammaMC2 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma2MCLabel));
          gamma2MotherLabel=gammaMC2->GetMother();
        }

        if(pi0MotherLabel==gamma2MotherLabel && ((AliAODMCParticle*)fAODMCTrackArray->At(gamma2MotherLabel))->GetPdgCode() == 223){
          isTrueOmega = kTRUE;
        }
        if(isTrueOmega){
          AliAODConversionMother TruePi0 = AliAODConversionMother(TrueGammaCandidate0, TrueGammaCandidate1);
          fHistoTrueOmegaInvMassPt[fiCut]->Fill(OmegaCandidate->M(),OmegaCandidate->Pt(),fWeightJetJetMC);
          fHistoTruePi0FromOmegaInvMassPt[fiCut]->Fill(TruePi0.M(),TruePi0.Pt(),fWeightJetJetMC);
          fHistoTrueGammaFromOmegaPt[fiCut]->Fill(TrueGammaCandidate2->Pt(),fWeightJetJetMC);
          if(!fDoLightOutput){
            FillQAPlotsMC(OmegaCandidate, &TruePi0, TrueGammaCandidate0, TrueGammaCandidate1, TrueGammaCandidate2);
          }
        }
      }
    }
      break;
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::CalculateBackground(){

  //zbin and mbin from pi0handler
  Int_t Pi0zbin = fBGPi0Handler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t Pi0mbin = 0;

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    Pi0mbin = fBGPi0Handler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  }else {
    Pi0mbin = fBGPi0Handler[fiCut]->GetMultiplicityBinIndex(fPi0Candidates->GetEntries());
  }

  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEvent1Vertex = NULL;

  // background generation using pi0s and gammas
  // pi0 from handler, gamma from current event
  for(Int_t previous=0;previous<fBGPi0Handler[fiCut]->GetNBGEvents();previous++){
    AliGammaConversionMotherAODVector* previousPi0s = fBGPi0Handler[fiCut]->GetBGGoodMesons(Pi0zbin,Pi0mbin,previous);
    if(fMoveParticleAccordingToVertex == kTRUE){
      bgEvent1Vertex = fBGPi0Handler[fiCut]->GetBGEventVertex(Pi0zbin,Pi0mbin,previous);
    }
    for(UInt_t iPi0=0;iPi0<previousPi0s->size();iPi0++){
      AliAODConversionMother BGpi0cand = (AliAODConversionMother)(*(previousPi0s->at(iPi0)));
      if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
        MoveParticleAccordingToVertex(&BGpi0cand, bgEvent1Vertex);
      }
      if(fReconMethod % 2 == 0){ //EMCAL gamma
        for(Int_t iCurrent3=0;iCurrent3<fClusterCandidates->GetEntries();iCurrent3++){
          AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iCurrent3));
          if(gamma2 == NULL || !(gamma2->GetIsCaloPhoton())) continue;
          AliAODConversionMother BGOmegacand = AliAODConversionMother(&BGpi0cand,gamma2);
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&BGOmegacand, &BGpi0cand, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            fHistoDiffPi0SameGammaBackInvMassPt[fiCut]->Fill(BGOmegacand.M(),BGOmegacand.Pt(),fWeightJetJetMC);
          }
        }
      } else{ //PCM gamma
        for(Int_t iCurrent3=0;iCurrent3<fGammaCandidates->GetEntries();iCurrent3++){
          AliAODConversionPhoton *gamma2 = dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent3));
          if(gamma2 == NULL) continue;
          AliAODConversionMother BGOmegacand = AliAODConversionMother(&BGpi0cand,gamma2);
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&BGOmegacand, &BGpi0cand, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            fHistoDiffPi0SameGammaBackInvMassPt[fiCut]->Fill(BGOmegacand.M(),BGOmegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    }
  }

  //pi0 from current event, gamma from handler
  //zbin and mbin from clushandler
  Int_t Cluszbin = fBGClusHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t Clusmbin = 0;

  if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
    Clusmbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  }else {
    Clusmbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
  }

  //zbin and mbin from gammahandler
  Int_t Gammazbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t Gammambin = 0;

  if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
    Gammambin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  }else {
    Gammambin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
  }

  for(Int_t iCurrent=0;iCurrent<fPi0Candidates->GetEntries();iCurrent++){
    AliAODConversionMother* BGpi0cand = dynamic_cast<AliAODConversionMother*>(fPi0Candidates->At(iCurrent));
    if(BGpi0cand == NULL) continue;
    if(fReconMethod % 2 == 0){ //EMCAL gamma
      for(Int_t previous=0;previous<fBGClusHandler[fiCut]->GetNBGEvents();previous++){
        AliGammaConversionAODVector *previousclusters = fBGClusHandler[fiCut]->GetBGGoodV0s(Cluszbin,Clusmbin,previous);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGClusHandler[fiCut]->GetBGEventVertex(Cluszbin,Clusmbin,previous);
        }
        for(UInt_t igamma2=0;igamma2<previousclusters->size();igamma2++){
          AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousclusters->at(igamma2)));
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
          }
          AliAODConversionMother BGOmegacand = AliAODConversionMother(BGpi0cand,&gamma2);
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&BGOmegacand, BGpi0cand, &gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            fHistoSamePi0DiffGammaBackInvMassPt[fiCut]->Fill(BGOmegacand.M(),BGOmegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    } else{ //PCM gamma
      for(Int_t previous=0;previous<fBGHandler[fiCut]->GetNBGEvents();previous++){
        AliGammaConversionAODVector *previousV0s = fBGHandler[fiCut]->GetBGGoodV0s(Gammazbin,Gammambin,previous);
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
          bgEvent1Vertex = fBGHandler[fiCut]->GetBGEventVertex(Gammazbin,Gammambin,previous);
        }
        for(UInt_t igamma2=0;igamma2<previousV0s->size();igamma2++){
          AliAODConversionPhoton gamma2 = (AliAODConversionPhoton)(*(previousV0s->at(igamma2)));
          if(fMoveParticleAccordingToVertex == kTRUE && bgEvent1Vertex){
            MoveParticleAccordingToVertex(&gamma2,bgEvent1Vertex);
          }
          AliAODConversionMother BGOmegacand = AliAODConversionMother(BGpi0cand,&gamma2);
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&BGOmegacand, BGpi0cand, &gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor)){
            fHistoSamePi0DiffGammaBackInvMassPt[fiCut]->Fill(BGOmegacand.M(),BGOmegacand.Pt(),fWeightJetJetMC);
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::CalculateOmegaRotationBackground(Int_t iCurrentGamma, Int_t iCurrentPi0){

  Double_t tempBGCandidateWeight       = fWeightJetJetMC;

  Double_t rotationAngle = TMath::Pi()/2.0; // rotaion angle 90°

  TLorentzVector lvRotationPhoton3;      // photon candidate from omega which gets rotated
  TLorentzVector lvRotationPhoton2;      // photon candidate from pion which gets rotated
  TLorentzVector lvRotationPhoton1;      // photon candidate from pion which gets rotated
  TLorentzVector lvRotationPion;        // pion candidate which gets rotated
  TVector3 lvRotationOmega;             // reconstructed mother particle from photon and pion
  // Needed for TGenPhaseSpace
  TVector3 tvEtaPhiGamma, tvEtaPhiPion, tvEtaPhiGammaDecay, tvEtaPhiPionDecay, tvNormBeforeDecay, tvNormAfterDecay;
  Float_t asymBeforeDecay = 0.;
  Float_t asymAfterDecay = 0.;
  Double_t massPi0Gamma[2] = {0,((AliAODConversionMother*)(fPi0Candidates->At(iCurrentPi0)))->M()};
  Double_t mass2Photons[2] = {0, 0};

  Int_t cellIDRotatedPhoton3 = -1; // cell ID of the theoretical cluster after rotation (photon candidate from omega)
  Int_t cellIDRotatedPhoton2 = -1; // cell ID of the theoretical cluster after rotation (photon candidate from pion)
  Int_t cellIDRotatedPhoton1 = -1; // cell ID of the theoretical cluster after rotation (photon candidate from pion)

  std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
  std::vector<std::array<Double_t, 2>> vSwappingInvMassPTAlphaCut;
  vSwappingInvMassPT.clear();
  vSwappingInvMassPTAlphaCut.clear();
  vSwappingInvMassPT.resize(0);
  vSwappingInvMassPTAlphaCut.resize(0);
  Double_t tempMultWeightSwapping = 1; // weight taking multiplicity of event into account

  // curcial requierment is that the event has at least 4 cluster candidates and one Pi0 candidate
  if( (fClusterCandidates->GetEntries() > 3 ) &&  (fPi0Candidates->GetEntries() > 0) )
  {

    AliAODConversionPhoton* currentEventRotatedPhoton3 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrentGamma));
    if (currentEventRotatedPhoton3==NULL || !(currentEventRotatedPhoton3->GetIsCaloPhoton())) { return;}
    // if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALO.find(iCurrentGamma) != dropOutGammas_CALO.end() ) ) {return; }
    AliAODConversionMother* currentEventGoodV0Temp2 = (AliAODConversionMother*)(fPi0Candidates->At(iCurrentPi0));

    AliAODConversionPhoton* currentEventRotatedPhoton1 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(currentEventGoodV0Temp2->GetLabel(0)));
    AliAODConversionPhoton* currentEventRotatedPhoton2 = dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(currentEventGoodV0Temp2->GetLabel(1)));

    Int_t backClusterIndex[3] = {iCurrentGamma, currentEventGoodV0Temp2->GetLabel(1), currentEventGoodV0Temp2->GetLabel(0)};  // save the indices of the cluster which are rotated

    for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){

      lvRotationPhoton3.SetX(currentEventRotatedPhoton3->Px());
      lvRotationPhoton3.SetY(currentEventRotatedPhoton3->Py());
      lvRotationPhoton3.SetZ(currentEventRotatedPhoton3->Pz());
      lvRotationPhoton3.SetE(currentEventRotatedPhoton3->E());

      lvRotationPhoton2.SetX(currentEventRotatedPhoton2->Px());
      lvRotationPhoton2.SetY(currentEventRotatedPhoton2->Py());
      lvRotationPhoton2.SetZ(currentEventRotatedPhoton2->Pz());
      lvRotationPhoton2.SetE(currentEventRotatedPhoton2->E());

      lvRotationPhoton1.SetX(currentEventRotatedPhoton1->Px());
      lvRotationPhoton1.SetY(currentEventRotatedPhoton1->Py());
      lvRotationPhoton1.SetZ(currentEventRotatedPhoton1->Pz());
      lvRotationPhoton1.SetE(currentEventRotatedPhoton1->E());

      lvRotationPion.SetX(currentEventGoodV0Temp2->Px());
      lvRotationPion.SetY(currentEventGoodV0Temp2->Py());
      lvRotationPion.SetZ(currentEventGoodV0Temp2->Pz());
      lvRotationPion.SetE(currentEventGoodV0Temp2->E());

      lvRotationOmega = (lvRotationPhoton3 + lvRotationPion).Vect();            // create omegacandidate which serves as axis to rotate around

      // rotate both photons around the momentum vector of their hypothetical mother particle
      if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0 || ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1)){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0) {
          rotationAngle = 90. + fRandom.Uniform(-5., 5.); // rotate by 90+-5 degree
          rotationAngle = rotationAngle*TMath::Pi()/180.;
        }
        else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
          Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
          rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
        }
        lvRotationPhoton3.Rotate(rotationAngle, lvRotationOmega);
        lvRotationPhoton2.Rotate(rotationAngle, lvRotationOmega);
        lvRotationPhoton1.Rotate(rotationAngle, lvRotationOmega);
        lvRotationPion.Rotate(rotationAngle, lvRotationOmega);
      } else if (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() >=10){ // generate new decay with TGenPhaseSpace
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
          tvEtaPhiGamma = lvRotationPhoton3.Vect();
          tvEtaPhiPion = lvRotationPion.Vect();
          tvNormBeforeDecay = tvEtaPhiGamma.Cross(tvEtaPhiPion);
          asymBeforeDecay = fabs((lvRotationPhoton3.E()-lvRotationPion.E())/(lvRotationPhoton3.E()+lvRotationPion.E()));
        }

        TLorentzVector lvRotationMother = lvRotationPhoton3 + lvRotationPion;
        fGenPhaseSpace.SetDecay(lvRotationMother, 2, massPi0Gamma);
        fGenPhaseSpace.Generate();
        lvRotationPhoton3 = *fGenPhaseSpace.GetDecay(0);
        lvRotationPion = *fGenPhaseSpace.GetDecay(1);

        fGenPhaseSpace.SetDecay(lvRotationPion, 2,  mass2Photons);
        fGenPhaseSpace.Generate();
        lvRotationPhoton1 = *fGenPhaseSpace.GetDecay(0);
        lvRotationPhoton2 = *fGenPhaseSpace.GetDecay(1);

        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
          tvEtaPhiGammaDecay = lvRotationPhoton3.Vect();
          tvEtaPhiPionDecay = lvRotationPion.Vect();
          tvNormAfterDecay = tvEtaPhiGammaDecay.Cross(tvEtaPhiPionDecay);  // norm vector to decay plane
          asymAfterDecay = fabs((lvRotationPhoton3.E()-lvRotationPion.E())/(lvRotationPhoton3.E()+lvRotationPion.E()));
          // check if decay is nearly the same as original decay: if yes continue with next decay
          if((tvNormAfterDecay.Angle(tvNormBeforeDecay) < 20*TMath::Pi()/180. || tvNormAfterDecay.Angle(tvNormBeforeDecay) > 340*TMath::Pi()/180.) && ( fabs(asymBeforeDecay - asymAfterDecay) < 0.05 )   ) continue;
        }

      }


      cellIDRotatedPhoton3 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton3.Eta(), static_cast<double>((lvRotationPhoton3.Phi()<0) ? lvRotationPhoton3.Phi() + TMath::Pi()*2. : lvRotationPhoton3.Phi()));
      cellIDRotatedPhoton2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()));
      cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));

      std::unique_ptr<AliAODConversionPhoton> currentEventGoodPhotonRotation3 (new AliAODConversionPhoton(&lvRotationPhoton3));
      std::unique_ptr<AliAODConversionPhoton> currentEventGoodPhotonRotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));
      std::unique_ptr<AliAODConversionPhoton> currentEventGoodPhotonRotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));

      if( ( ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton3, lvRotationPhoton3.Phi(), fInputEvent) ) || ( lvRotationPhoton3.E() <= ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut) )->GetMinClusterEnergy() ) )
      {
        backClusterIndex[0] = -1;
      }
      if( ( ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, lvRotationPhoton2.Phi(), fInputEvent) ) || ( lvRotationPhoton2.E() <= ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut) )->GetMinClusterEnergy() ) )
      {
        backClusterIndex[1] = -1;
      }
      if( ( ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent) ) || ( lvRotationPhoton1.E() <= ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut) )->GetMinClusterEnergy() ) )
      {
        backClusterIndex[2] = -1;
      }
      if( (!fDoLightOutput) && (fDoMesonQA & 0b00001000) )
      {
        if(backClusterIndex[0] != -1) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton3.Eta(), static_cast<double>((lvRotationPhoton3.Phi()<0) ? lvRotationPhoton3.Phi() + TMath::Pi()*2. : lvRotationPhoton3.Phi()), 1);
        if(backClusterIndex[1] != -1) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()), 1);
        if(backClusterIndex[2] != -1) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), 1);
      }
      // When no rotated cluster matches with the EMCal/DCal skipp this rotation
      if( (backClusterIndex[0] == -1) && (backClusterIndex[1] == -1) && (backClusterIndex[2] == -1) ) {continue;}
      PhotonSelectionCaloBack(backClusterIndex, currentEventGoodPhotonRotation3.get(), currentEventGoodPhotonRotation2.get(), currentEventGoodPhotonRotation1.get());

      // loop over all possible pi0 which are made of two non altered clusters
      for (Int_t iPi0 = 0; iPi0 < fPi0CandidatesBackOri->GetEntries(); iPi0++)
      {
        AliAODConversionMother* kCurrentPi0Candidate = dynamic_cast<AliAODConversionMother*>(fPi0CandidatesBackOri->At(iPi0));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(0)));
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(1)));

        // check for altered cluster originating from the omega candidate if omega back candidate is possible
        if(backClusterIndex[0] != -1)                                           // chech if the altered cluster even is in acceptance
        {
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALOBack.find(backClusterIndex[0]) == dropOutGammas_CALOBack.end() ) )
          {
            AliAODConversionMother backgroundCandidate = AliAODConversionMother(kCurrentPi0Candidate,currentEventGoodPhotonRotation3.get());
            if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation3.get(), fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) )
            {
              if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation3.get()), GetQTPi0(&backgroundCandidate, kCurrentPi0Candidate) ) )
                && (backgroundCandidate.M() < 1.6) )
              {
                if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&backgroundCandidate, gamma0, gamma1, currentEventGoodPhotonRotation3.get() ) )
                {
                  Double_t dalitzcand01 = kCurrentPi0Candidate->M2();
                  Double_t dalitzcand02 = AliAODConversionMother(gamma0,currentEventGoodPhotonRotation3.get()).M2();
                  Double_t dalitzcand12 = AliAODConversionMother(gamma1,currentEventGoodPhotonRotation3.get()).M2();

                  if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                  {
                    fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                    fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                  }
                  if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                    if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                    {
                      fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                      fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                    }
                    if((!fDoLightOutput) && TMath::Abs(backgroundCandidate.GetAlpha())<0.1)
                    {
                      vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                    }
                  }
                }
              }
            }
          }
        }

        // check for altered cluster 2 originating from the pi0 candidate if omega back candidate is possible
        if(backClusterIndex[1] != -1)                                           // chech if the altered cluster even is in acceptance
        {
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALOBack.find(backClusterIndex[1]) == dropOutGammas_CALOBack.end() ) )
          {
            AliAODConversionMother backgroundCandidate = AliAODConversionMother(kCurrentPi0Candidate,currentEventGoodPhotonRotation2.get());
            if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation2.get(), fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) )
            {
              if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation2.get()), GetQTPi0(&backgroundCandidate, kCurrentPi0Candidate) ) )
                && (backgroundCandidate.M() < 1.6) )
              {
                if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&backgroundCandidate, gamma0, gamma1, currentEventGoodPhotonRotation2.get() ) )
                {
                  Double_t dalitzcand01 = kCurrentPi0Candidate->M2();
                  Double_t dalitzcand02 = AliAODConversionMother(gamma0,currentEventGoodPhotonRotation2.get()).M2();
                  Double_t dalitzcand12 = AliAODConversionMother(gamma1,currentEventGoodPhotonRotation2.get()).M2();

                  if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                  {
                    fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                    fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                  }
                  if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                    if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                    {
                      fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                      fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                    }
                    if((!fDoLightOutput) && TMath::Abs(backgroundCandidate.GetAlpha())<0.1)
                    {
                      vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                    }
                  }
                }
              }
            }
          }
        }

        // check for altered cluster 1 originating from the pi0 candidate if omega back candidate is possible
        if(backClusterIndex[2] != -1)                                           // chech if the altered cluster even is in acceptance
        {
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALOBack.find(backClusterIndex[2]) == dropOutGammas_CALOBack.end() ) )
          {
            AliAODConversionMother backgroundCandidate = AliAODConversionMother(kCurrentPi0Candidate,currentEventGoodPhotonRotation1.get());
            if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation1.get(), fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) )
            {
              if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation1.get()), GetQTPi0(&backgroundCandidate, kCurrentPi0Candidate) ) )
                && (backgroundCandidate.M() < 1.6) )
              {
                if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&backgroundCandidate, gamma0, gamma1, currentEventGoodPhotonRotation1.get() ) )
                {
                  Double_t dalitzcand01 = kCurrentPi0Candidate->M2();
                  Double_t dalitzcand02 = AliAODConversionMother(gamma0,currentEventGoodPhotonRotation1.get()).M2();
                  Double_t dalitzcand12 = AliAODConversionMother(gamma1,currentEventGoodPhotonRotation1.get()).M2();

                  if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                  {
                    fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                    fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                  }
                  if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                    if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                    {
                      fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                      fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                    }
                    if((!fDoLightOutput) && TMath::Abs(backgroundCandidate.GetAlpha())<0.1)
                    {
                      vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                    }
                  }
                }
              }
            }
          }
        }
      }

      // loop over all possible pi0 which are made of one altered clusters and one non altered cluster
      for (Int_t iPi0 = 0; iPi0 < fPi0CandidatesBackAlt->GetEntries(); iPi0++)
      {
        AliAODConversionMother* kCurrentPi0Candidate = dynamic_cast<AliAODConversionMother*>(fPi0CandidatesBackAlt->At(iPi0));
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(0)));
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(1)));
        for(Int_t iOrgGamma = 0;iOrgGamma<fClusterCandidates->GetEntries();iOrgGamma++)
        {
          if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALOBack.find(iOrgGamma) != dropOutGammas_CALOBack.end() ) ) {continue;}
          if( (iOrgGamma == backClusterIndex[0]) || (iOrgGamma == backClusterIndex[1]) || (iOrgGamma == backClusterIndex[2]) ) continue;
          AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iOrgGamma));
          if (gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
          AliAODConversionMother backgroundCandidate = AliAODConversionMother(kCurrentPi0Candidate,gamma2);
          if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&backgroundCandidate, kCurrentPi0Candidate, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) )
          {
            if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&backgroundCandidate, kCurrentPi0Candidate, gamma2), GetQTPi0(&backgroundCandidate, kCurrentPi0Candidate) ) )
              && (backgroundCandidate.M() < 1.6) )
            {
              if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&backgroundCandidate, gamma0, gamma1, gamma2) )
              {
                Double_t dalitzcand01 = kCurrentPi0Candidate->M2();
                Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();

                if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                {
                  fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                  fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                }
                if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                  if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                  {
                    fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand02, tempBGCandidateWeight);
                    fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand12, tempBGCandidateWeight);
                  }
                  if((!fDoLightOutput) && TMath::Abs(backgroundCandidate.GetAlpha())<0.1)
                  {
                    vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate.M(),backgroundCandidate.Pt()});
                  }
                }
              }
            }
          }
        }
      }
      dropOutGammas_CALOBack.clear();
      fPi0CandidatesBackOri->Clear();
      fPi0CandidatesBackAlt->Clear();
    }

    // Fill the histograms
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPT.size() > 0){
      if( (fReconMethod == 0) || (fReconMethod == 2) || (fReconMethod == 4) ) {
        tempMultWeightSwapping = (0.5*(fClusterCandidates->GetEntries()*fClusterCandidates->GetEntries() - fClusterCandidates->GetEntries()))/(vSwappingInvMassPT.size());
      }
      else{
        tempMultWeightSwapping = (0.5*(fGammaCandidates->GetEntries()*fGammaCandidates->GetEntries() - fGammaCandidates->GetEntries()))/(vSwappingInvMassPT.size());
      }
    }
    for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
      fHistoMotherSwappingBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*tempBGCandidateWeight);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::CalculatePi0RotationBackground(){

  Double_t rotationAngle = TMath::Pi()/2.0; //0.78539816339; // rotaion angle 90°

  TLorentzVector lvRotationPhoton1;   // photon candidates which get rotated
  TLorentzVector lvRotationPhoton2;   // photon candidates which get rotated
  TVector3 lvRotationPion;            // reconstructed mother particle from the two photons
  // Needed for TGenPhaseSpace
  TVector3 tvEtaPhigamma1, tvEtaPhigamma2, tvEtaPhigamma1Decay, tvEtaPhigamma2Decay, tvNormBeforeDecay, tvNormAfterDecay;
  Float_t asymBeforeDecay = 0.;
  Float_t asymAfterDecay = 0.;
  Double_t massGamma[2] = {0,0};

  Int_t cellIDRotatedPhoton1 = -1; // cell ID of the cluster after rotation
  Int_t cellIDRotatedPhoton2 = -1; // cell ID of the cluster after rotation

  // curcial requierment is that the event has at least 4 cluster candidates and
  // at least one pi0 candidate
  if(fClusterCandidates->GetEntries() > 3) {

    for(Int_t iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
      AliAODConversionPhoton* currentEventRotatedPhoton1 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent1));

      for(Int_t iCurrent2=iCurrent1+1;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
        AliAODConversionPhoton* currentEventRotatedPhoton2 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent2));

        Int_t backClusterIndex[3] = {-1, iCurrent2, iCurrent1};  // save the indices of the cluster which are rotated

        for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){

          lvRotationPhoton1.SetX(currentEventRotatedPhoton1->Px());
          lvRotationPhoton1.SetY(currentEventRotatedPhoton1->Py());
          lvRotationPhoton1.SetZ(currentEventRotatedPhoton1->Pz());
          lvRotationPhoton1.SetE(currentEventRotatedPhoton1->E());

          lvRotationPhoton2.SetX(currentEventRotatedPhoton2->Px());
          lvRotationPhoton2.SetY(currentEventRotatedPhoton2->Py());
          lvRotationPhoton2.SetZ(currentEventRotatedPhoton2->Pz());
          lvRotationPhoton2.SetE(currentEventRotatedPhoton2->E());

          lvRotationPion = (lvRotationPhoton1 + lvRotationPhoton2).Vect();

          // rotate both photons around the momentum vector of their hypothetical mother particle
          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0 || ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1)){
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0){
              rotationAngle = 90. + fRandom.Uniform(-5., 5.); // rotate by 90+-5 degree
              rotationAngle = rotationAngle*TMath::Pi()/180.;
            }
            else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
               Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
               rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
            }
            lvRotationPhoton1.Rotate(rotationAngle ,lvRotationPion);
            lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);
          } else if (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() >=10){ // generate new decay with TGenPhaseSpace
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
              tvEtaPhigamma1 = lvRotationPhoton1.Vect();
              tvEtaPhigamma2 = lvRotationPhoton2.Vect();
              tvNormBeforeDecay = tvEtaPhigamma1.Cross(tvEtaPhigamma2);
              asymBeforeDecay = fabs((lvRotationPhoton1.E()-lvRotationPhoton2.E())/(lvRotationPhoton1.E()+lvRotationPhoton2.E()));
            }

            TLorentzVector lvRotationMother = lvRotationPhoton1 + lvRotationPhoton2;
            fGenPhaseSpace.SetDecay(lvRotationMother, 2, massGamma);
            fGenPhaseSpace.Generate();
            lvRotationPhoton1 = *fGenPhaseSpace.GetDecay(0);
            lvRotationPhoton2 = *fGenPhaseSpace.GetDecay(1);

            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
              tvEtaPhigamma1Decay = lvRotationPhoton1.Vect();
              tvEtaPhigamma2Decay = lvRotationPhoton2.Vect();
              tvNormAfterDecay = tvEtaPhigamma1Decay.Cross(tvEtaPhigamma2Decay);  // norm vector to decay plane
              asymAfterDecay = fabs((lvRotationPhoton1.E()-lvRotationPhoton2.E())/(lvRotationPhoton1.E()+lvRotationPhoton2.E()));
              // check if decay is nearly the same as original decay: if yes continue with next decay
              if((tvNormAfterDecay.Angle(tvNormBeforeDecay) < 20*TMath::Pi()/180. || tvNormAfterDecay.Angle(tvNormBeforeDecay) > 340*TMath::Pi()/180.) && ( fabs(asymBeforeDecay - asymAfterDecay) < 0.05 )   ) continue;
            }

          }


          cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));
          cellIDRotatedPhoton2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()));

          std::unique_ptr<AliAODConversionPhoton> currentEventGoodPhotonRotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
          std::unique_ptr<AliAODConversionPhoton> currentEventGoodPhotonRotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));

          if( ( ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent)) || ( lvRotationPhoton1.E() <= ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy() ) )
          {
            backClusterIndex[2] = -1;
          }
          if( ( ( (AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, lvRotationPhoton2.Phi(), fInputEvent)) || ( lvRotationPhoton2.E() <= ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy() ) )
          {
            backClusterIndex[1] = -1;
          }
          if( (!fDoLightOutput) && (fDoMesonQA & 0b00001000) )
          {
            if(backClusterIndex[2] != -1) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), 1);
            if(backClusterIndex[1] != -1) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()), 1);
          }
          if( (backClusterIndex[0] == -1) && (backClusterIndex[1] == -1) && (backClusterIndex[2] == -1) ) {continue;}
          PhotonSelectionCaloBack(backClusterIndex, 0, currentEventGoodPhotonRotation2.get(), currentEventGoodPhotonRotation1.get());

          // loop over all possible pi0 which are made of two non altered clusters
          for (Int_t iPi0 = 0; iPi0 < fPi0CandidatesBackOri->GetEntries(); iPi0++)
          {
            AliAODConversionMother* kCurrentPi0Candidate = dynamic_cast<AliAODConversionMother*>(fPi0CandidatesBackOri->At(iPi0));
            AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(0)));
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(1)));
            // check for altered cluster 2 originating from the pi0 candidate if omega back candidate is possible
            if(backClusterIndex[1] != -1)                                       // chech if the altered cluster even is in acceptance
            {
              if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALOBack.find(backClusterIndex[1]) == dropOutGammas_CALOBack.end() ) )
              {
                AliAODConversionMother backgroundCandidate = AliAODConversionMother(kCurrentPi0Candidate,currentEventGoodPhotonRotation2.get());
                if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation2.get(), fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) )
                {
                  if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation2.get()), GetQTPi0(&backgroundCandidate, kCurrentPi0Candidate) ) )
                    && (backgroundCandidate.M() < 1.6) )
                  {
                    if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&backgroundCandidate, gamma0, gamma1, currentEventGoodPhotonRotation2.get() ) )
                    {
                      Double_t dalitzcand01 = kCurrentPi0Candidate->M2();
                      Double_t dalitzcand02 = AliAODConversionMother(gamma0,currentEventGoodPhotonRotation2.get()).M2();
                      Double_t dalitzcand12 = AliAODConversionMother(gamma1,currentEventGoodPhotonRotation2.get()).M2();
                      if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                      {
                        fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                        fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                      }
                      if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                      {
                        fHistoMotherSwappingBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(), backgroundCandidate.Pt(), fWeightJetJetMC);
                        if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                        {
                          fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                          fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                        }
                      }
                    }
                  }
                }
              }
            }

            // check for altered cluster 1 originating from the pi0 candidate if omega back candidate is possible
            if(backClusterIndex[2] != -1)                                       // chech if the altered cluster even is in acceptance
            {
              if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALOBack.find(backClusterIndex[2]) == dropOutGammas_CALOBack.end() ) )
              {
                AliAODConversionMother backgroundCandidate = AliAODConversionMother(kCurrentPi0Candidate,currentEventGoodPhotonRotation1.get());
                if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation1.get(), fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) )
                {
                  if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&backgroundCandidate, kCurrentPi0Candidate, currentEventGoodPhotonRotation1.get()), GetQTPi0(&backgroundCandidate, kCurrentPi0Candidate) ) )
                    && (backgroundCandidate.M() < 1.6) )
                  {
                    if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&backgroundCandidate, gamma0, gamma1, currentEventGoodPhotonRotation1.get() ) )
                    {
                      Double_t dalitzcand01 = kCurrentPi0Candidate->M2();
                      Double_t dalitzcand02 = AliAODConversionMother(gamma0,currentEventGoodPhotonRotation1.get()).M2();
                      Double_t dalitzcand12 = AliAODConversionMother(gamma1,currentEventGoodPhotonRotation1.get()).M2();
                      if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                      {
                        fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                        fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                      }
                      if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                      {
                        fHistoMotherSwappingBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(), backgroundCandidate.Pt(), fWeightJetJetMC);
                        if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                        {
                          fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                          fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          // loop over all possible pi0 which are made of one altered clusters and one non altered cluster
          for (Int_t iPi0 = 0; iPi0 < fPi0CandidatesBackAlt->GetEntries(); iPi0++)
          {
            AliAODConversionMother* kCurrentPi0Candidate = dynamic_cast<AliAODConversionMother*>(fPi0CandidatesBackAlt->At(iPi0));
            AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(0)));
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(kCurrentPi0Candidate->GetLabel(1)));
            for(Int_t iOrgGamma = 0;iOrgGamma<fClusterCandidates->GetEntries();iOrgGamma++)
            {
              if ( ( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseGammaSelection() ) && (dropOutGammas_CALOBack.find(iOrgGamma) == dropOutGammas_CALOBack.end() ) )
              {
                if( (iOrgGamma == backClusterIndex[1]) || (iOrgGamma == backClusterIndex[2]) || (iOrgGamma == kCurrentPi0Candidate->GetLabel(0)) || (iOrgGamma == kCurrentPi0Candidate->GetLabel(1)) ) continue;
                AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(iOrgGamma));
                if (gamma2==NULL || !(gamma2->GetIsCaloPhoton())) continue;
                AliAODConversionMother backgroundCandidate = AliAODConversionMother(kCurrentPi0Candidate,gamma2);
                if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaAngle(&backgroundCandidate, kCurrentPi0Candidate, gamma2, fDoPiZeroGammaAngleCut, fmaxfit, flowerFactor, fupperFactor) )
                {
                  if( ( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->ArmenterosLikeQtCut(GetPodAlpha(&backgroundCandidate, kCurrentPi0Candidate, gamma2), GetQTPi0(&backgroundCandidate, kCurrentPi0Candidate) ) )
                    && (backgroundCandidate.M() < 1.6) )
                  {
                    if( ( (AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedPiZeroGammaOAC(&backgroundCandidate, gamma0, gamma1, gamma2) )
                    {
                      Double_t dalitzcand01 = kCurrentPi0Candidate->M2();
                      Double_t dalitzcand02 = AliAODConversionMother(gamma0,gamma2).M2();
                      Double_t dalitzcand12 = AliAODConversionMother(gamma1,gamma2).M2();
                      if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                      {
                        fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                        fHistoBackDalitzPlotwoCut[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                      }
                      if( (!fUseDalitzCut) || ( (fUseDalitzCut) && (dalitzcand02 >= 0.02) && (dalitzcand12 >= 0.02) ) )
                      {
                        fHistoMotherSwappingBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(), backgroundCandidate.Pt(), fWeightJetJetMC);
                        if( (!fDoLightOutput) && (fDoMesonQA & 0b10000000) )
                        {
                          fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand02, fWeightJetJetMC);
                          fHistoBackDalitzPlot[fiCut]->Fill(dalitzcand01, dalitzcand12, fWeightJetJetMC);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          dropOutGammas_CALOBack.clear();
          fPi0CandidatesBackOri->Clear();
          fPi0CandidatesBackAlt->Clear();
        }
      }
    }
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
void AliAnalysisTaskOmegaToPiZeroGamma::MoveParticleAccordingToVertex(AliAODConversionMother* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
   //see header file for documentation

   Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
   Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
   Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

   Double_t movedPlace[3] = {particle->GetProductionX() - dx,particle->GetProductionY() - dy,particle->GetProductionZ() - dz};
   particle->SetProductionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::UpdateEventByEventData(){
  //see header file for documentation
  if(fPi0Candidates->GetEntries()>0){
    if(fReconMethod % 2 == 1){ //single gamma is PCM
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
      }else { // means we use #V0s for multiplicity
        fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
      }
    } else{ //single gamma is EMCAL
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->UseTrackMultiplicity()){
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
      }else {
        fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),fEventPlaneAngle);
      }
    }
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGPi0Handler[fiCut]->AddMesonEvent(fPi0Candidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
    }else {
      fBGPi0Handler[fiCut]->AddMesonEvent(fPi0Candidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fPi0Candidates->GetEntries(),0);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel

  if(mode){
    fMCEventPos = new Int_t[fReaderGammas->GetEntries()];
    fMCEventNeg = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayPos = new Int_t[fReaderGammas->GetEntries()];
    fESDArrayNeg = new Int_t[fReaderGammas->GetEntries()];
  }

  for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
    if(!PhotonCandidate) continue;
    if(!mode){// Back to ESD Labels
    PhotonCandidate->SetMCLabelPositive(fMCEventPos[iGamma]);
    PhotonCandidate->SetMCLabelNegative(fMCEventNeg[iGamma]);
    PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
    PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
    continue;
    }
    fMCEventPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
    fMCEventNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
    fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
    fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();

    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;

    for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
      AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if(!AODLabelPos){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
        }
      }
      if(!AODLabelNeg){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
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
    delete[] fMCEventPos;
    delete[] fMCEventNeg;
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
void AliAnalysisTaskOmegaToPiZeroGamma::SetLogBinningYTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetYaxis();
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
    if (TMath::Abs(pdgCode) == 310) return 1; // k0s
    else if (TMath::Abs(pdgCode) == 3122) return 2; // Lambda
    else if (TMath::Abs(pdgCode) == 130) return 3; // K0L
    else if (TMath::Abs(pdgCode) == 2212) return 4; // proton
    else if (TMath::Abs(pdgCode) == 2112) return 5; // neutron
    else if (TMath::Abs(pdgCode) == 211) return 6; // pion
    else if (TMath::Abs(pdgCode) == 321) return 7; // kaon
    else if (TMath::Abs(pdgCode) == 113 || TMath::Abs(pdgCode) == 213 ) return 8; // rho 0,+,-
    else if (TMath::Abs(pdgCode) == 3222 || TMath::Abs(pdgCode) == 3212 || TMath::Abs(pdgCode) == 3112  ) return 9; // Sigma
    else if (TMath::Abs(pdgCode) == 2224 || TMath::Abs(pdgCode) == 2214 || TMath::Abs(pdgCode) == 2114 || TMath::Abs(pdgCode) == 1114  ) return 10; // Delta
    else if (TMath::Abs(pdgCode) == 313 || TMath::Abs(pdgCode) == 323   ) return 11; // K*
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

//_________________________________________________________________________________
Double_t AliAnalysisTaskOmegaToPiZeroGamma::GetPodAlpha(AliAODConversionMother* motherPart, AliAODConversionMother* daughter1, AliAODConversionPhoton* daughter2){
  return (cos(motherPart->Angle(daughter1->Vect())) * daughter1->P() - cos(motherPart->Angle(daughter2->Vect())) * daughter2->P()) / (cos(motherPart->Angle(daughter1->Vect())) * daughter1->P() + cos(motherPart->Angle(daughter2->Vect())) * daughter2->P());
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskOmegaToPiZeroGamma::GetPodAlpha(AliAODConversionMother* motherPart, AliAODConversionPhoton* daughter1, AliAODConversionPhoton* daughter2){
  return (cos(motherPart->Angle(daughter1->Vect())) * daughter1->P() - cos(motherPart->Angle(daughter2->Vect())) * daughter2->P()) / (cos(motherPart->Angle(daughter1->Vect())) * daughter1->P() + cos(motherPart->Angle(daughter2->Vect())) * daughter2->P());
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskOmegaToPiZeroGamma::GetQTPi0(AliAODConversionMother* motherPart, AliAODConversionMother* daughter){
  return  sin(motherPart->Angle(daughter->Vect())) * daughter->P();
}

//_________________________________________________________________________________
Double_t AliAnalysisTaskOmegaToPiZeroGamma::GetQTGamma(AliAODConversionMother* motherPart, AliAODConversionPhoton* daughter){
  return  sin(motherPart->Angle(daughter->Vect())) * daughter->P();
}

//_________________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::FillQAPlots(AliAODConversionMother *omegacand, AliAODConversionMother *pi0cand,
  AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1, AliAODConversionPhoton *gamma2){
    if(fDoMesonQA & 0b00000001) {
      fHistoMotherYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00000010) {
      fHistoMotherAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00000100) {
    }
    if(fDoMesonQA & 0b00001000) {
      switch (fReconMethod) {
        case 0: // PCM-Cal-Cal
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
          break;

        case 1: // PCM-Cal-PCM
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma2->GetPhotonPhi(), gamma2->GetPhotonEta(),fWeightJetJetMC);
          break;

        case 3: // Cal-Cal-PCM
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma2->GetPhotonPhi(), gamma2->GetPhotonEta(),fWeightJetJetMC);
          break;

        case 4: // PCM-PCM-Cal
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma1->GetPhotonPhi(), gamma1->GetPhotonEta(),fWeightJetJetMC);
          break;

        case 5: // PCM-PCM-PCM
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma1->GetPhotonPhi(), gamma1->GetPhotonEta(),fWeightJetJetMC);
          fHistoMotherConvPhotonEtaPhi[fiCut]->Fill(gamma2->GetPhotonPhi(), gamma2->GetPhotonEta(),fWeightJetJetMC);
          break;

        default:
          /* do nothing */
          break;
      }
      fHistoMotherEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00010000) {
      fHistoMotherPi0AnglePt[fiCut]->Fill(omegacand->Pt(),TMath::Pi() - omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
      fHistoMotherGammaAnglePt[fiCut]->Fill(omegacand->Pt(),TMath::Pi() - omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
      fHistoPi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00100000) {
      TVector3 Boost = omegacand->BoostVector();
      TLorentzVector Pi0InRestFrame = TLorentzVector(0,0,0,0);
      Pi0InRestFrame.SetPxPyPzE(pi0cand->Px(), pi0cand->Py(), pi0cand->Pz(), pi0cand->E());
      Pi0InRestFrame.Boost(-Boost);
      TLorentzVector GammaInRestFrame = TLorentzVector(0,0,0,0);
      GammaInRestFrame.SetPxPyPzE(gamma2->Px(), gamma2->Py(), gamma2->Pz(), gamma2->E());
      GammaInRestFrame.Boost(-Boost);

      fHistoMotherRestGammaCosAnglePt[fiCut]->Fill(omegacand->Pt(), cos(omegacand->Angle(GammaInRestFrame.Vect()) ),fWeightJetJetMC);
      fHistoMotherRestPi0CosAnglePt[fiCut]->Fill(omegacand->Pt(), cos(omegacand->Angle(Pi0InRestFrame.Vect()) ),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b01000000) {
      fHistoRecoArmenterosPodolanskiPlot[fiCut]->Fill(GetPodAlpha(omegacand, pi0cand, gamma2), GetQTPi0(omegacand, pi0cand));
    }
    if(fDoMesonQA & 0b10000000) {
      TLorentzVector gamma0LV = TLorentzVector(0,0,0,0);
      gamma0LV.SetPxPyPzE(gamma0->Px(), gamma0->Py(), gamma0->Pz(), gamma0->E());
      TLorentzVector gamma1LV = TLorentzVector(0,0,0,0);
      gamma1LV.SetPxPyPzE(gamma1->Px(), gamma1->Py(), gamma1->Pz(), gamma1->E());
      TLorentzVector gamma2LV = TLorentzVector(0,0,0,0);
      gamma2LV.SetPxPyPzE(gamma2->Px(), gamma2->Py(), gamma2->Pz(), gamma2->E());

      fHistoMotherDalitzPlot[fiCut]->Fill( (gamma0LV + gamma1LV).M2(), (gamma1LV + gamma2LV).M2(),fWeightJetJetMC);
      fHistoMotherDalitzPlot[fiCut]->Fill( (gamma0LV + gamma1LV).M2(), (gamma0LV + gamma2LV).M2(),fWeightJetJetMC); // also since 12 and 02 should follow the same structure
    }
  return;
}


//_________________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::FillQAPlotsMC(AliAODConversionMother *omegacand, AliAODConversionMother *pi0cand,
  AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1, AliAODConversionPhoton *gamma2){
    if(fDoMesonQA & 0b00000001) {
      fHistoTrueOmegaYPt[fiCut]->Fill(omegacand->Pt(),omegacand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
      fHistoTruePi0FromOmegaYPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);

    }
    if(fDoMesonQA & 0b00000010) {
      fHistoTrueOmegaAlphaPt[fiCut]->Fill(omegacand->Pt(),omegacand->GetAlpha(),fWeightJetJetMC);
      fHistoTruePi0FromOmegaAlphaPt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha(),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00000100) {
      fHistoTruePi0FromOmegaOpenAnglePt[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00001000) {
      fHistoTrueOmegaEtaPhi[fiCut]->Fill(omegacand->Phi(),omegacand->Eta(),fWeightJetJetMC);
      fHistoTruePi0FromOmegaEtaPhi[fiCut]->Fill(pi0cand->Phi(),pi0cand->Eta(),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00010000) {
      fHistoTrueOmegaGammaAnglePt[fiCut]->Fill(omegacand->Pt(),TMath::Pi() - omegacand->Angle(gamma2->Vect()),fWeightJetJetMC);
      fHistoTrueOmegaPi0AnglePt[fiCut]->Fill(omegacand->Pt(),TMath::Pi() - omegacand->Angle(pi0cand->Vect()),fWeightJetJetMC);
      fHistoTruePi0GammaAnglePt[fiCut]->Fill(omegacand->Pt(),pi0cand->Angle(gamma2->Vect()),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b00100000) {
      TVector3 Boost = omegacand->BoostVector();
      TLorentzVector Pi0InRestFrame = TLorentzVector(0,0,0,0);
      Pi0InRestFrame.SetPxPyPzE(pi0cand->Px(), pi0cand->Py(), pi0cand->Pz(), pi0cand->E());
      Pi0InRestFrame.Boost(-Boost);
      TLorentzVector GammaInRestFrame = TLorentzVector(0,0,0,0);
      GammaInRestFrame.SetPxPyPzE(gamma2->Px(), gamma2->Py(), gamma2->Pz(), gamma2->E());
      GammaInRestFrame.Boost(-Boost);

      fHistoTrueOmegaRestGammaCosAnglePt[fiCut]->Fill(omegacand->Pt(), cos(omegacand->Angle(GammaInRestFrame.Vect()) ),fWeightJetJetMC);
      fHistoTrueOmegaRestPi0CosAnglePt[fiCut]->Fill(omegacand->Pt(), cos(omegacand->Angle(Pi0InRestFrame.Vect()) ),fWeightJetJetMC);
    }
    if(fDoMesonQA & 0b01000000) {
      fHistoTrueArmenterosPodolanskiPlot[fiCut]->Fill(GetPodAlpha(omegacand, pi0cand, gamma2), GetQTPi0(omegacand, pi0cand));
    }
    if(fDoMesonQA & 0b10000000) {
      TLorentzVector gamma0LV = TLorentzVector(0,0,0,0);
      gamma0LV.SetPxPyPzE(gamma0->Px(), gamma0->Py(), gamma0->Pz(), gamma0->E());
      TLorentzVector gamma1LV = TLorentzVector(0,0,0,0);
      gamma1LV.SetPxPyPzE(gamma1->Px(), gamma1->Py(), gamma1->Pz(), gamma1->E());
      TLorentzVector gamma2LV = TLorentzVector(0,0,0,0);
      gamma2LV.SetPxPyPzE(gamma2->Px(), gamma2->Py(), gamma2->Pz(), gamma2->E());

      fHistoTrueDalitzPlot[fiCut]->Fill( (gamma0LV + gamma1LV).M2(), (gamma1LV + gamma2LV).M2(),fWeightJetJetMC);
      fHistoTrueDalitzPlot[fiCut]->Fill( (gamma0LV + gamma1LV).M2(), (gamma0LV + gamma2LV).M2(),fWeightJetJetMC); // also since 12 and 02 should follow the same structure
    }
    if(fDoMesonQA & 0b100000000) {
      fHistoTrueOmegaPtGammaPt[fiCut]->Fill(omegacand->Pt(),gamma2->Pt(),fWeightJetJetMC);
      fHistoTrueOmegaPtPi0Pt[fiCut]->Fill(omegacand->Pt(),pi0cand->Pt(),fWeightJetJetMC);
    }
  return;
}

//______________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::PhotonSelectionCalo() {
  for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
    if (gamma0==NULL || !(gamma0->GetIsCaloPhoton())) continue;
    for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
      AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
      if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
      AliAODConversionMother pi0cand = AliAODConversionMother(gamma0, gamma1);
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
           ->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())
          && pi0cand.Pt() > fMinPi0Pt){
        if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
          dropOutGammas_CALO.insert(firstGammaIndex);
          dropOutGammas_CALO.insert(secondGammaIndex);
        }
      }
    }
  }
  return;
}


//______________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::PhotonSelectionPCM() {
  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
    if (gamma0==NULL) continue;
    for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
      AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
      if (gamma1==NULL || (secondGammaIndex+1)>=(fGammaCandidates->GetEntries())) continue;
      AliAODConversionMother pi0cand = AliAODConversionMother(gamma0, gamma1);
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
           ->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())
          && pi0cand.Pt() > fMinPi0Pt){
        if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
          dropOutGammas_PCM.insert(firstGammaIndex);
          dropOutGammas_PCM.insert(secondGammaIndex);
        }
      }
    }
  }
  return;
}

//______________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::PhotonSelectionMixed() {
  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
    if (gamma0==NULL) continue;
    for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
      Bool_t matched = kFALSE;
      AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
      if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
      AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
      matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0, cluster, fInputEvent, fWeightJetJetMC);
      AliAODConversionMother pi0cand = AliAODConversionMother(gamma0,gamma1);
      if(((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))
           ->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())
          && pi0cand.Pt() > fMinPi0Pt){
        if (!matched){
          if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
            dropOutGammas_PCM.insert(firstGammaIndex);
            dropOutGammas_CALO.insert(secondGammaIndex);
          }
        }
      }
    }
  }
  return;
}

//______________________________________________________________________________
void AliAnalysisTaskOmegaToPiZeroGamma::PhotonSelectionCaloBack(Int_t rotatedGammaIndex[3], AliAODConversionPhoton* rotatedPhoton3, AliAODConversionPhoton* rotatedPhoton2, AliAODConversionPhoton* rotatedPhoton1) {

  // loop over all clusters that are not rotated/from TGPS in the current event in the current background cycle
  for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
    if( (firstGammaIndex == rotatedGammaIndex[0]) || (firstGammaIndex == rotatedGammaIndex[1]) || (firstGammaIndex == rotatedGammaIndex[2]) ) continue;
    AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
    if (gamma1==NULL || !(gamma1->GetIsCaloPhoton())) continue;
    // try to make a pi0 with combining unchnaged cluster with changed photon directly comming from omega candidate
    if(rotatedGammaIndex[0] != -1)                                              // chech if the altered cluster is accepted
    {
      AliAODConversionMother pi0cand3 = AliAODConversionMother(rotatedPhoton3, gamma1);
      if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(&pi0cand3 ,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
          && pi0cand3.Pt() > fMinPi0Pt)
        {
        if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand3, 0) ) {
          if(fPhotonSelectionMode >= 1)
          {
            dropOutGammas_CALOBack.insert(rotatedGammaIndex[0]);
            dropOutGammas_CALOBack.insert(firstGammaIndex);
          }
          AliAODConversionMother *pi0backcand_vec = new AliAODConversionMother(rotatedPhoton3,gamma1);
          pi0backcand_vec->SetLabel1(rotatedGammaIndex[0]);
          pi0backcand_vec->SetLabel2(firstGammaIndex);
          pi0backcand_vec->SetPxPyPzE(pi0backcand_vec->Px(),pi0backcand_vec->Py(),pi0backcand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0backcand_vec->P()*pi0backcand_vec->P()));
          fPi0CandidatesBackAlt->Add(pi0backcand_vec);
        }
      }
    }

    // try to make a pi0 with combining unchnaged cluster with changed photon 2 comming from pi0 from omega candidate
    if(rotatedGammaIndex[1] != -1)                                              // chech if the altered cluster is accepted
    {
      AliAODConversionMother pi0cand2 = AliAODConversionMother(rotatedPhoton2, gamma1);
      if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(&pi0cand2,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
          && pi0cand2.Pt() > fMinPi0Pt)
        {
        if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand2, 0) ) {
          if(fPhotonSelectionMode >= 1)
          {
            dropOutGammas_CALOBack.insert(rotatedGammaIndex[1]);
            dropOutGammas_CALOBack.insert(firstGammaIndex);
          }
          AliAODConversionMother *pi0backcand_vec = new AliAODConversionMother(rotatedPhoton2,gamma1);
          pi0backcand_vec->SetLabel1(rotatedGammaIndex[1]);
          pi0backcand_vec->SetLabel2(firstGammaIndex);
          pi0backcand_vec->SetPxPyPzE(pi0backcand_vec->Px(),pi0backcand_vec->Py(),pi0backcand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0backcand_vec->P()*pi0backcand_vec->P()));
          fPi0CandidatesBackAlt->Add(pi0backcand_vec);
        }
      }
    }

    // try to make a pi0 with combining unchnaged cluster with changed photon 1 comming from pi0 from omega candidate
    if(rotatedGammaIndex[2] != -1)                                              // chech if the altered cluster is accepted
    {
      AliAODConversionMother pi0cand1 = AliAODConversionMother(rotatedPhoton1, gamma1);
      if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(&pi0cand1,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
          && pi0cand1.Pt() > fMinPi0Pt)
        {
        if( ((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand1, 0) ) {
          if(fPhotonSelectionMode >= 1)
          {
            dropOutGammas_CALOBack.insert(rotatedGammaIndex[2]);
            dropOutGammas_CALOBack.insert(firstGammaIndex);
          }
          AliAODConversionMother *pi0backcand_vec = new AliAODConversionMother(rotatedPhoton1,gamma1);
          pi0backcand_vec->SetLabel1(rotatedGammaIndex[2]);
          pi0backcand_vec->SetLabel2(firstGammaIndex);
          pi0backcand_vec->SetPxPyPzE(pi0backcand_vec->Px(),pi0backcand_vec->Py(),pi0backcand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0backcand_vec->P()*pi0backcand_vec->P()));
          fPi0CandidatesBackAlt->Add(pi0backcand_vec);
        }
      }
    }

    // try to make a pi0 with combining unchnaged cluster with another unchanged cluster
    for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
      if( (secondGammaIndex == rotatedGammaIndex[0]) || (secondGammaIndex == rotatedGammaIndex[1]) || (secondGammaIndex == rotatedGammaIndex[2]) ) continue;
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
      if (gamma0==NULL || !(gamma0->GetIsCaloPhoton())) continue;
      AliAODConversionMother pi0cand = AliAODConversionMother(gamma0, gamma1);
      if((((AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelected(&pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
          && pi0cand.Pt() > fMinPi0Pt)
        {
        if( ( (AliConversionMesonCuts*)fNeutralPionCutArray->At(fiCut))->MesonIsSelectedByMassCut(&pi0cand, 0) ) {
          if(fPhotonSelectionMode >= 1)
          {
            dropOutGammas_CALOBack.insert(firstGammaIndex);
            dropOutGammas_CALOBack.insert(secondGammaIndex);
          }
          AliAODConversionMother *pi0backcand_vec = new AliAODConversionMother(gamma0,gamma1);
          pi0backcand_vec->SetLabel1(secondGammaIndex);
          pi0backcand_vec->SetLabel2(firstGammaIndex);
          pi0backcand_vec->SetPxPyPzE(pi0backcand_vec->Px(),pi0backcand_vec->Py(),pi0backcand_vec->Pz(),TMath::Sqrt(0.1349766*0.1349766+pi0backcand_vec->P()*pi0backcand_vec->P()));
          fPi0CandidatesBackOri->Add(pi0backcand_vec);
        }
      }
    }
  }
  return;
}
