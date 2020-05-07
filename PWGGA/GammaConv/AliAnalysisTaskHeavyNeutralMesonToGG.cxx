/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Sandro Wenzel                                 *
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
#include "AliAnalysisTaskHeavyNeutralMesonToGG.h"
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

ClassImp(AliAnalysisTaskHeavyNeutralMesonToGG)

//________________________________________________________________________
AliAnalysisTaskHeavyNeutralMesonToGG::AliAnalysisTaskHeavyNeutralMesonToGG(): AliAnalysisTaskSE(),
  fRandom(0),
  fV0Reader(NULL),
  fBGHandler(NULL),
  fBGHandlerRP(NULL),
  fBGClusHandler(NULL),
  fBGClusHandlerRP(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fEventCuts(NULL),
  fConversionCuts(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCuts(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputContainer(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fClusterCutArray(NULL),
  fMesonCutArray(NULL),
  fReaderGammas(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fFileNameBroken(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fTreeBrokenFiles(NULL),
  fProfileTruePrimaryMesonWeightsInvMassPt(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fHistoMotherMesonPtY(NULL),
  fHistoMotherMesonPtAlpha(NULL),
  fHistoMotherMesonPtOpenAngle(NULL),
  fHistoMotherMesonConvPhotonEtaPhi(NULL),
  fHistoTrueMesonInvMassPt(NULL),
  fHistoTrueMesonMatchedInvMassPt(NULL),
  fHistoTrueMesonCaloPhotonInvMassPt(NULL),
  fHistoTrueMesonCaloConvertedPhotonInvMassPt(NULL),
  fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt(NULL),
  fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt(NULL),
  fHistoTrueMesonCaloElectronInvMassPt(NULL),
  fHistoTrueMesonCaloMergedClusterInvMassPt(NULL),
  fHistoTrueMesonCaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTruePrimaryMesonInvMassPt(NULL),
  fHistoTruePrimaryMesonW0WeightingInvMassPt(NULL),
  fHistoTruePrimaryMesonMCPtResolPt(NULL),
  fHistoTrueMotherMesonConvPhotonEtaPhi(NULL),
  fHistoTrueBckGGInvMassPt(NULL),
  fHistoTrueBckFullMesonContainedInOneClusterInvMassPt(NULL),
  fHistoTrueBckAsymEClustersInvMassPt(NULL),
  fHistoTrueBckContInvMassPt(NULL),
  fHistoTrueMesonPtY(NULL),
  fHistoTrueMesonPtAlpha(NULL),
  fHistoTrueMesonPtOpenAngle(NULL),
  fHistoMCMesonPtY(NULL),
  fHistoMCMesonPtAlpha(NULL),
  fHistoMCMesonPtJetPt(NULL),
  fHistoTrueNLabelsInClusPt(NULL),
  fHistoDoubleCountTrueMesonInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoMCHeaders(NULL),
  fHistoConvGammaPt(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
  fHistoMotherInvMassRejected(NULL),
  fHistoMCMesonPt(NULL),
  fHistoMCMesonWOWeightPt(NULL),
  fHistoMCMesonWOEvtWeightPt(NULL),
  fHistoMCMesonInAccPt(NULL),
  fHistoMCMesonWOWeightInAccPt(NULL),
  fHistoMCMesonWOEvtWeightInAccPt(NULL),
  fHistoTrueConvGammaPt(NULL),
  fHistoTruePrimaryConvGammaPt(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusConvGammaPt(NULL),
  fHistoTrueClusConvGammaFullyPt(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoMultipleCountTrueMeson(NULL),
  fHistoMultipleCountTrueConvGamma(NULL),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoVertexX(NULL),
  fHistoVertexY(NULL),
  fHistoNGammaConvCandidates(NULL),
  fHistoNGammaCaloCandidates(NULL),
  fHistoNV0Tracks(NULL),
  fHistoJetJetNTrials(NULL),
  fMapMultipleCountTrueMesons(),
  fMapMultipleCountTrueConvGammas(),
  fMapMultipleCountTrueClusterGammas(),
  fVectorRecTrueMesons(0),
  fVectorDoubleCountTrueMesons(0),
  fVectorDoubleCountTrueConvGammas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMesonInvMassWindow(NULL),
  fMCEventPos(NULL),
  fMCEventNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fEventPlaneAngle(-100),
  fMesonInvMassMin(0),
  fMesonInvMassMax(0),
  fMesonInvMassNBins(0),
  fWeightJetJetMC(1),
  fNGammaCandidates(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fMesonRecoMode(-1),
  fMesonType(-1),
  fMesonPDG(0),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoClusterQA(0),
  fIsMC(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fDoLightOutput(kFALSE),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fDoConvGammaShowerShapeTree(kFALSE),
  fEnableSortForClusMC(kFALSE),
  fDoPrimaryTrackMatching(kFALSE),
  fDoInvMassShowerShapeTree(kFALSE),
  fAllowOverlapHeaders(kTRUE),
  fEnableClusterCutsForTrigger(kFALSE),
  fGenPhaseSpace()
{

}

//________________________________________________________________________
AliAnalysisTaskHeavyNeutralMesonToGG::AliAnalysisTaskHeavyNeutralMesonToGG(const char *name):
  AliAnalysisTaskSE(name),
  fRandom(0),
  fV0Reader(NULL),
  fBGHandler(NULL),
  fBGHandlerRP(NULL),
  fBGClusHandler(NULL),
  fBGClusHandlerRP(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fEventCuts(NULL),
  fConversionCuts(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCuts(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputContainer(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fClusterCutArray(NULL),
  fMesonCutArray(NULL),
  fReaderGammas(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fFileNameBroken(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fTreeBrokenFiles(NULL),
  fProfileTruePrimaryMesonWeightsInvMassPt(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fHistoMotherMesonPtY(NULL),
  fHistoMotherMesonPtAlpha(NULL),
  fHistoMotherMesonPtOpenAngle(NULL),
  fHistoMotherMesonConvPhotonEtaPhi(NULL),
  fHistoTrueMesonInvMassPt(NULL),
  fHistoTrueMesonMatchedInvMassPt(NULL),
  fHistoTrueMesonCaloPhotonInvMassPt(NULL),
  fHistoTrueMesonCaloConvertedPhotonInvMassPt(NULL),
  fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt(NULL),
  fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt(NULL),
  fHistoTrueMesonCaloElectronInvMassPt(NULL),
  fHistoTrueMesonCaloMergedClusterInvMassPt(NULL),
  fHistoTrueMesonCaloMergedClusterPartConvInvMassPt(NULL),
  fHistoTruePrimaryMesonInvMassPt(NULL),
  fHistoTruePrimaryMesonW0WeightingInvMassPt(NULL),
  fHistoTruePrimaryMesonMCPtResolPt(NULL),
  fHistoTrueMotherMesonConvPhotonEtaPhi(NULL),
  fHistoTrueBckGGInvMassPt(NULL),
  fHistoTrueBckFullMesonContainedInOneClusterInvMassPt(NULL),
  fHistoTrueBckAsymEClustersInvMassPt(NULL),
  fHistoTrueBckContInvMassPt(NULL),
  fHistoTrueMesonPtY(NULL),
  fHistoTrueMesonPtAlpha(NULL),
  fHistoTrueMesonPtOpenAngle(NULL),
  fHistoMCMesonPtY(NULL),
  fHistoMCMesonPtAlpha(NULL),
  fHistoMCMesonPtJetPt(NULL),
  fHistoTrueNLabelsInClusPt(NULL),
  fHistoDoubleCountTrueMesonInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoMCHeaders(NULL),
  fHistoConvGammaPt(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
  fHistoMotherInvMassRejected(NULL),
  fHistoMCMesonPt(NULL),
  fHistoMCMesonWOWeightPt(NULL),
  fHistoMCMesonWOEvtWeightPt(NULL),
  fHistoMCMesonInAccPt(NULL),
  fHistoMCMesonWOWeightInAccPt(NULL),
  fHistoMCMesonWOEvtWeightInAccPt(NULL),
  fHistoTrueConvGammaPt(NULL),
  fHistoTruePrimaryConvGammaPt(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusConvGammaPt(NULL),
  fHistoTrueClusConvGammaFullyPt(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoMultipleCountTrueMeson(NULL),
  fHistoMultipleCountTrueConvGamma(NULL),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoVertexX(NULL),
  fHistoVertexY(NULL),
  fHistoNGammaConvCandidates(NULL),
  fHistoNGammaCaloCandidates(NULL),
  fHistoNV0Tracks(NULL),
  fHistoJetJetNTrials(NULL),
  fMapMultipleCountTrueMesons(),
  fMapMultipleCountTrueConvGammas(),
  fMapMultipleCountTrueClusterGammas(),
  fVectorRecTrueMesons(0),
  fVectorDoubleCountTrueMesons(0),
  fVectorDoubleCountTrueConvGammas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fUnsmearedPx(NULL),
  fUnsmearedPy(NULL),
  fUnsmearedPz(NULL),
  fUnsmearedE(NULL),
  fMesonInvMassWindow(NULL),
  fMCEventPos(NULL),
  fMCEventNeg(NULL),
  fESDArrayPos(NULL),
  fESDArrayNeg(NULL),
  fEventPlaneAngle(-100),
  fMesonInvMassMin(0),
  fMesonInvMassMax(0),
  fMesonInvMassNBins(0),
  fWeightJetJetMC(1),
  fNGammaCandidates(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fMesonRecoMode(-1),
  fMesonType(-1),
  fMesonPDG(0),
  fDoMesonQA(0),
  fDoPhotonQA(0),
  fDoClusterQA(0),
  fIsMC(0),
  fMoveParticleAccordingToVertex(kTRUE),
  fDoLightOutput(kFALSE),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fDoTHnSparse(kTRUE),
  fSetPlotHistsExtQA(kFALSE),
  fDoConvGammaShowerShapeTree(kFALSE),
  fEnableSortForClusMC(kFALSE),
  fDoPrimaryTrackMatching(kFALSE),
  fDoInvMassShowerShapeTree(kFALSE),
  fAllowOverlapHeaders(kTRUE),
  fEnableClusterCutsForTrigger(kFALSE),
  fGenPhaseSpace()
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskHeavyNeutralMesonToGG::~AliAnalysisTaskHeavyNeutralMesonToGG()
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
  if(fBGHandlerRP){
    delete[] fBGHandlerRP;
    fBGHandlerRP = 0x0;
  }
  if(fBGClusHandler){
    delete[] fBGClusHandler;
    fBGClusHandler = 0x0;
  }
  if(fBGClusHandlerRP){
    delete[] fBGClusHandlerRP;
    fBGClusHandlerRP = 0x0;
  }
}
//___________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::InitBack(){

  const Int_t nDim = 4;
  Int_t nBins[nDim] = {800,300,7,4};
  Double_t xMin[nDim] = {0,0, 0,0};
  Double_t xMax[nDim] = {1.2,30,7,4};

  if(fDoTHnSparse){
    fSparseMotherInvMassPtZM      = new THnSparseF*[fnCuts];
    fSparseMotherBackInvMassPtZM  = new THnSparseF*[fnCuts];
  }

  fBGHandler        = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGHandlerRP      = new AliConversionAODBGHandlerRP*[fnCuts];

  fBGClusHandler    = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGClusHandlerRP  = new AliConversionAODBGHandlerRP*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
      TString cutstringEvent        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringConvGamma    = "";
      TString cutstringCaloGamma    = "";
      if (fMesonRecoMode < 2) cutstringConvGamma    = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger ) cutstringCaloGamma    = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      TString fullCutString         = "";
      if (fMesonRecoMode == 0)      fullCutString = Form("%i_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(), cutstringConvGamma.Data(), cutstringMeson.Data());
      else if (fMesonRecoMode == 1) fullCutString = Form("%i_%s_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(), cutstringMeson.Data());
      else if (fMesonRecoMode == 2) fullCutString = Form("%i_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(), cutstringCaloGamma.Data(), cutstringMeson.Data());

      Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
      Int_t centMin = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
      Int_t centMax = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));

      if( collisionSystem == 1 || collisionSystem == 2 ||
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

      if(fDoTHnSparse){
        fBackList[iCut] = new TList();
        fBackList[iCut]->SetName(Form("%s Back histograms",fullCutString.Data()));
        fBackList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fBackList[iCut]);

        fSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m", "Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        fBackList[iCut]->Add(fSparseMotherBackInvMassPtZM[iCut]);

        fMotherList[iCut] = new TList();
        fMotherList[iCut]->SetName(Form("%s Mother histograms",fullCutString.Data()));
        fMotherList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMotherList[iCut]);

        fSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m", "Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
        fMotherList[iCut]->Add(fSparseMotherInvMassPtZM[iCut]);
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
        fBGHandler[iCut] = new AliGammaConversionAODBGHandler(  collisionSystem,centMin,centMax,
                                                                ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                                                ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                                                2,8,5);
        fBGClusHandler[iCut] = new AliGammaConversionAODBGHandler( collisionSystem,centMin,centMax,
                                                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                                                  2,8,5);
        fBGHandlerRP[iCut] = NULL;
      }else if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() != 2){
        fBGHandlerRP[iCut] = new AliConversionAODBGHandlerRP(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsHeavyIon(),
                                                              ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                                              ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
        fBGClusHandlerRP[iCut] = new AliConversionAODBGHandlerRP( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsHeavyIon(),
                                                                  ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
        fBGHandler[iCut] = NULL;
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::UserCreateOutputObjects(){

    fV0Reader = (AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
    if(!fV0Reader){printf("Error: No V0 Reader");return;}// GetV0Reader
    if (fIsMC == 2){
      fDoClusterQA      = 0;
      fDoTHnSparse      = kFALSE;
    } else if (fIsMC == 3){
      fDoTHnSparse      = kFALSE;
    }

    if (fMesonRecoMode == 2)
      fMoveParticleAccordingToVertex  = kFALSE;

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

    if(fDoTHnSparse){
      fBackList         = new TList*[fnCuts];
      fMotherList       = new TList*[fnCuts];
    }

    fHistoNEvents               = new TH1F*[fnCuts];
    if (fIsMC > 1){
      fHistoNEventsWOWeight     = new TH1F*[fnCuts];
    }
    if (fIsMC == 2){
      fProfileJetJetXSection    = new TProfile*[fnCuts];
      fHistoJetJetNTrials       = new TH1F*[fnCuts];
    }
    fHistoNGoodESDTracks        = new TH1F*[fnCuts];
    fHistoVertexZ               = new TH1F*[fnCuts];
    if(!fDoLightOutput){
      fHistoVertexX               = new TH1F*[fnCuts];
      fHistoVertexY               = new TH1F*[fnCuts];
    }
    if (fMesonRecoMode < 2) fHistoNGammaConvCandidates      = new TH1F*[fnCuts];
    if(fIsHeavyIon==2) fProfileEtaShift = new TProfile*[fnCuts];
    if(!fDoLightOutput){
      fHistoSPDClusterTrackletBackground      = new TH2F*[fnCuts];
      fHistoNV0Tracks                         = new TH1F*[fnCuts];
      if (fMesonRecoMode < 2)
        fHistoConvGammaPt                       = new TH1F*[fnCuts];
    }

    if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger) {
      fHistoNGammaCaloCandidates      = new TH1F*[fnCuts];
      if (!fDoLightOutput ){
        fHistoClusGammaPt                   = new TH1F*[fnCuts];
        fHistoClusGammaE                    = new TH1F*[fnCuts];
        if (fIsMC > 0){
          fHistoClusOverlapHeadersGammaPt     = new TH1F*[fnCuts];
          fHistoClusAllHeadersGammaPt         = new TH1F*[fnCuts];
          fHistoClusRejectedHeadersGammaPt    = new TH1F*[fnCuts];
        }
      }
    }

    fHistoMotherInvMassPt             = new TH2F*[fnCuts];
    fHistoMotherInvMassRejected       = new TH1F*[fnCuts];
    fHistoMotherBackInvMassPt         = new TH2F*[fnCuts];
    if(!fDoLightOutput && fMesonRecoMode == 1){
      fHistoMotherMatchedInvMassPt      = new TH2F*[fnCuts];
    }
    if (fDoMesonQA > 0){
      fHistoMotherMesonPtY              = new TH2F*[fnCuts];
      fHistoMotherMesonPtAlpha          = new TH2F*[fnCuts];
      fHistoMotherMesonPtOpenAngle      = new TH2F*[fnCuts];
      fHistoMotherMesonConvPhotonEtaPhi = new TH2F*[fnCuts];
    }

    fMesonInvMassWindow               = new Double_t[2];
    if (fMesonType < 0 || fMesonType > 2){
      if(!fV0Reader){printf("Error: No V0 Reader");return;}// GetV0Reader
    } else if (fMesonType == 0){ // pi0 case 134.9770 ± 0.0005 MeV
      fMesonPDG               = 111;
      fMesonInvMassMin        = 0.;
      fMesonInvMassMax        = 0.400;
      fMesonInvMassNBins      = 400;
      fMesonInvMassWindow[0]  = 0.05;
      fMesonInvMassWindow[1]  = 0.17;
    } else if (fMesonType == 1){ // eta case 547.862 ± 0.017 MeV
      fMesonPDG               = 221;
      fMesonInvMassMin        = 0.300;
      fMesonInvMassMax        = 0.800;
      fMesonInvMassNBins      = 500;
      fMesonInvMassWindow[0]  = 0.45;
      fMesonInvMassWindow[1]  = 0.65;
    } else if (fMesonType == 2){ // eta' case 957.78 ± 0.06 MeV
      fMesonPDG               = 331;
      fMesonInvMassMin        = 0.700;
      fMesonInvMassMax        = 1.200;
      fMesonInvMassNBins      = 500;
      fMesonInvMassWindow[0]  = 0.85;
      fMesonInvMassWindow[1]  = 1.05;
    }
    // set common binning in pT for mesons and photons
    Int_t nBinsPt               = 200;
    Float_t binWidthPt          = 0.1;
    Float_t minPt               = 0;
    Float_t maxPt               = 20;
    Int_t nBinsQAPt             = 170;
    Float_t maxQAPt             = 20;
    Int_t nBinsClusterPt        = 500;
    Float_t minClusterPt        = 0;
    Float_t maxClusterPt        = 50;
    Double_t *arrPtBinning      = new Double_t[1200];
    Double_t *arrQAPtBinning    = new Double_t[1200];
    Double_t *arrClusPtBinning  = new Double_t[1200];
    // Set special pt binning for pp 8TeV
    if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k8TeV ){
      nBinsQAPt                 = 190;
      maxQAPt                   = 40;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<60) arrQAPtBinning[i]              = 0.05*i;
        else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
        else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
        else if(i<190) arrQAPtBinning[i]        = 20.+1.0*(i-170);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsPt                   = 400;
      minPt                     = 0;
      maxPt                     = 40;
      for(Int_t i=0; i<nBinsPt+1;i++){
        arrPtBinning[i]         = ((maxPt-minPt)/nBinsPt)*i;
      }
      nBinsClusterPt            = 800;
      minClusterPt              = 0;
      maxClusterPt              = 80;
      for(Int_t i=0; i<nBinsPt+1;i++){
        arrClusPtBinning[i]     = ((maxClusterPt-minClusterPt)/nBinsClusterPt)*i;
      }
      // Set special pt binning for pPb 5TeV
    } else if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeV ||
                ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeVR2 ||
                ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k5TeV      ){
      binWidthPt                = 0.05;
      nBinsPt                   = 205;
      minPt                     = 0;
      maxPt                     = 60;
      for(Int_t i=0; i<nBinsPt+1;i++){
        if (i < 1) arrPtBinning[i]              = 0.3*i;
        else if(i<55) arrPtBinning[i]           = 0.3+0.05*(i-1);
        else if(i<125) arrPtBinning[i]          = 3.+0.1*(i-55);
        else if(i<165) arrPtBinning[i]          = 10.+0.25*(i-125);
        else if(i<205) arrPtBinning[i]          = 20.+1.0*(i-165);
        else arrPtBinning[i]                    = maxPt;
      }
      nBinsQAPt                 = 210;
      maxQAPt                   = 60;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<60) arrQAPtBinning[i]              = 0.05*i;
        else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
        else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
        else if(i<210) arrQAPtBinning[i]        = 20.+1.0*(i-170);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsClusterPt            = 301;
      minClusterPt              = 0;
      maxClusterPt              = 100;
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        if (i < 1) arrClusPtBinning[i]          = 0.3*i;
        else if(i<55) arrClusPtBinning[i]       = 0.3+0.05*(i-1);
        else if(i<125) arrClusPtBinning[i]      = 3.+0.1*(i-55);
        else if(i<155) arrClusPtBinning[i]      = 10.+0.2*(i-125);
        else if(i<211) arrClusPtBinning[i]      = 16.+0.25*(i-155);
        else if(i<251) arrClusPtBinning[i]      = 30.+0.5*(i-211);
        else if(i<301) arrClusPtBinning[i]      = 50.+1.0*(i-251);
        else arrClusPtBinning[i]                = maxClusterPt;
      }
    // Set special pt binning for pp 13TeV, pPb 8TeV
    } else if ( ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeV ||
                ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeVLowB ||
                ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb8TeV ){
      nBinsPt                   = 285;
      minPt                     = 0;
      maxPt                     = 100;
      binWidthPt                = 0.05;
      for(Int_t i=0; i<nBinsPt+1;i++){
        if (i < 1) arrPtBinning[i]              = 0.3*i;
        else if(i<55) arrPtBinning[i]           = 0.3+0.05*(i-1);
        else if(i<125) arrPtBinning[i]          = 3.+0.1*(i-55);
        else if(i<185) arrPtBinning[i]          = 10.+0.25*(i-125);
        else if(i<235) arrPtBinning[i]          = 25.+0.5*(i-185);
        else if(i<285) arrPtBinning[i]          = 50.+1.0*(i-235);
        else  arrPtBinning[i]                   = maxPt;
      }
      nBinsQAPt                 = 270;
      maxQAPt                   = 100;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<60) arrQAPtBinning[i]              = 0.05*i;
        else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
        else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
        else if(i<210) arrQAPtBinning[i]        = 20.+0.5*(i-170);
        else if(i<270) arrQAPtBinning[i]        = 40.+1.0*(i-210);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsClusterPt            = 301;
      minClusterPt              = 0;
      maxClusterPt              = 100;
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        if (i < 1) arrClusPtBinning[i]          = 0.3*i;
        else if(i<55) arrClusPtBinning[i]       = 0.3+0.05*(i-1);
        else if(i<125) arrClusPtBinning[i]      = 3.+0.1*(i-55);
        else if(i<155) arrClusPtBinning[i]      = 10.+0.2*(i-125);
        else if(i<211) arrClusPtBinning[i]      = 16.+0.25*(i-155);
        else if(i<251) arrClusPtBinning[i]      = 30.+0.5*(i-211);
        else if(i<301) arrClusPtBinning[i]      = 50.+1.0*(i-251);
        else arrClusPtBinning[i]                = maxClusterPt;
      }
    } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kXeXe5440GeV  ){
      nBinsPt                   = 90;
      minPt                     = 0;
      maxPt                     = 20;
      for(Int_t i=0; i<nBinsPt+1;i++){
        if (i < 1) arrPtBinning[i]              = 0.3*i;
        else if(i<58) arrPtBinning[i]           = 0.3+0.1*(i-1);
        else if(i<82) arrPtBinning[i]           = 6.+0.25*(i-58);
        else if(i<90) arrPtBinning[i]           = 12.+1.0*(i-82);
        else arrPtBinning[i]                    = maxPt;
      }
      nBinsQAPt                 = 92;
      maxQAPt                   = 20;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<60) arrQAPtBinning[i]              = 0.1*i;
        else if(i<84) arrQAPtBinning[i]         = 6.+0.25*(i-60);
        else if(i<92) arrQAPtBinning[i]         = 12.+1.0*(i-84);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsClusterPt            = 148;
      minClusterPt              = 0;
      maxClusterPt              = 40;
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        if (i < 1) arrClusPtBinning[i]          = 0.3*i;
        else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
        else if(i<123) arrClusPtBinning[i]      = 10.+0.2*(i-98);
        else if(i<148) arrClusPtBinning[i]      = 15.+1.0*(i-123);
        else arrClusPtBinning[i]                = maxClusterPt;
      }
    } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kPbPb5TeV  ){
      nBinsPt                   = 90;
      minPt                     = 0;
      maxPt                     = 20;
      for(Int_t i=0; i<nBinsPt+1;i++){
        if (i < 1) arrPtBinning[i]              = 0.3*i;
        else if(i<58) arrPtBinning[i]           = 0.3+0.1*(i-1);
        else if(i<82) arrPtBinning[i]           = 6.+0.25*(i-58);
        else if(i<90) arrPtBinning[i]           = 12.+1.0*(i-82);
        else arrPtBinning[i]                    = maxPt;
      }
      nBinsQAPt                 = 92;
      maxQAPt                   = 20;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<60) arrQAPtBinning[i]              = 0.1*i;
        else if(i<84) arrQAPtBinning[i]         = 6.+0.25*(i-60);
        else if(i<92) arrQAPtBinning[i]         = 12.+1.0*(i-84);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsClusterPt            = 148;
      minClusterPt              = 0;
      maxClusterPt              = 40;
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        if (i < 1) arrClusPtBinning[i]          = 0.3*i;
        else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
        else if(i<123) arrClusPtBinning[i]      = 10.+0.2*(i-98);
        else if(i<148) arrClusPtBinning[i]      = 15.+1.0*(i-123);
        else arrClusPtBinning[i]                = maxClusterPt;
      }
              // default binning
    } else {
      for(Int_t i=0; i<nBinsPt+1;i++){
        arrPtBinning[i]         = ((maxPt-minPt)/nBinsPt)*i;
      }
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        arrClusPtBinning[i]     = ((maxClusterPt-minClusterPt)/nBinsClusterPt)*i;
      }
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<60) arrQAPtBinning[i]              = 0.05*i;
        else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
        else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
    }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringConvGamma    = "";
    TString cutstringCaloGamma    = "";
    if (fMesonRecoMode < 2) cutstringConvGamma    = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
    if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger) cutstringCaloGamma    = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    TString fullCutString         = "";
    if (fMesonRecoMode == 0)      fullCutString = Form("%i_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(), cutstringConvGamma.Data(), cutstringMeson.Data());
    else if (fMesonRecoMode == 1) fullCutString = Form("%i_%s_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(), cutstringMeson.Data());
    else if (fMesonRecoMode == 2) fullCutString = Form("%i_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(), cutstringCaloGamma.Data(), cutstringMeson.Data());

    fCutFolder[iCut]          = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s",fullCutString.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]            = new TList();
    fESDList[iCut]->SetName(Form("%s ESD histograms",fullCutString.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

    fHistoNEvents[iCut]       = new TH1F("NEvents", "NEvents", 14, -0.5, 13.5);
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
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut]   = new TH1F("NEventsWOWeight", "NEventsWOWeight", 14, -0.5, 13.5);
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
      fProfileJetJetXSection[iCut]  = new TProfile("XSection", "XSection", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials", "#sum{NTrials}", 1, 0, 1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }

    if(fIsHeavyIon == 1)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks; # TPC tracks", 4000, 0, 4000);
    else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks; # TPC tracks", 400, 0, 400);
    else
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks; # TPC tracks", 200, 0, 200);
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    fHistoVertexZ[iCut]             = new TH1F("VertexZ", "VertexZ", 200, -10, 10);
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);
    if(!fDoLightOutput){
      fHistoVertexX[iCut]             = new TH1F("VertexX", "VertexX", 100, -5, 5);
      fESDList[iCut]->Add(fHistoVertexX[iCut]);
      fHistoVertexY[iCut]             = new TH1F("VertexY", "VertexY", 100, -5, 5);
      fESDList[iCut]->Add(fHistoVertexY[iCut]);
    }

    if (fMesonRecoMode < 2){
      if(fIsHeavyIon == 1)
        fHistoNGammaConvCandidates[iCut]  = new TH1F("GammaConvCandidates", "GammaConvCandidates; # accepted #gamma_{conv}", 100, 0, 100);
      else if(fIsHeavyIon == 2)
        fHistoNGammaConvCandidates[iCut]  = new TH1F("GammaConvCandidates", "GammaConvCandidates; # accepted #gamma_{conv}", 50, 0, 50);
      else
        fHistoNGammaConvCandidates[iCut]  = new TH1F("GammaConvCandidates", "GammaConvCandidates; # accepted #gamma_{conv}", 50, 0, 50);
      fESDList[iCut]->Add(fHistoNGammaConvCandidates[iCut]);
    }
    if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
      if(fIsHeavyIon == 1)
        fHistoNGammaCaloCandidates[iCut]  = new TH1F("GammaCaloCandidates", "GammaCaloCandidates; # accepted #gamma_{conv}", 100, 0, 100);
      else if(fIsHeavyIon == 2)
        fHistoNGammaCaloCandidates[iCut]  = new TH1F("GammaCaloCandidates", "GammaCaloCandidates; # accepted #gamma_{conv}", 50, 0, 50);
      else
        fHistoNGammaCaloCandidates[iCut]  = new TH1F("GammaCaloCandidates", "GammaCaloCandidates; # accepted #gamma_{conv}", 50, 0, 50);
      fESDList[iCut]->Add(fHistoNGammaCaloCandidates[iCut]);
    }

    if(!fDoLightOutput){
      fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", 100, 0, 200, 250, 0, 1000);
      fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity", "V0 Multiplicity; VZERO amp [arb. units]", 30000, 0, 30000);
      else if(fIsHeavyIon == 2)
        fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity", "V0 Multiplicity; VZERO amp [arb. units]", 2500, 0, 2500);
      else
        fHistoNV0Tracks[iCut]         = new TH1F("V0 Multiplicity", "V0 Multiplicity; VZERO amp [arb. units]", 1500, 0, 1500);
      fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);

      if (fMesonRecoMode < 2){
        fHistoConvGammaPt[iCut]         = new TH1F("ESD_ConvGamma_Pt", "ESD_ConvGamma_Pt; p_{T,conv}(GeV/c)", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
      }
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
        fHistoClusGammaPt[iCut]         = new TH1F("ClusGamma_Pt", "ClusGamma_Pt; p_{T,clus} (GeV/c)", nBinsClusterPt, arrClusPtBinning);
        fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
        fHistoClusGammaE[iCut]          = new TH1F("ClusGamma_E", "ClusGamma_E; E_{clus} (GeV)", nBinsClusterPt, arrClusPtBinning);
        fESDList[iCut]->Add(fHistoClusGammaE[iCut]);
        if (fIsMC > 0){
          fHistoClusOverlapHeadersGammaPt[iCut]   = new TH1F("ClusGammaOverlapHeaders_Pt", "ClusGammaOverlapHeaders_Pt; p_{T,clus} (GeV/c), selected header w/ overlap",
                                                             nBinsClusterPt, arrClusPtBinning);
          fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
          fHistoClusAllHeadersGammaPt[iCut]       = new TH1F("ClusGammaAllHeaders_Pt", "ClusGammaAllHeaders_Pt; p_{T,clus} (GeV/c), all headers",
                                                             nBinsClusterPt, arrClusPtBinning);
          fESDList[iCut]->Add(fHistoClusAllHeadersGammaPt[iCut]);
          fHistoClusRejectedHeadersGammaPt[iCut]  = new TH1F("ClusGammaRejectedHeaders_Pt", "ClusGammaRejectedHeaders_Pt; p_{T,clus} (GeV/c), rejected headers",
                                                             nBinsClusterPt, arrClusPtBinning);
          fESDList[iCut]->Add(fHistoClusRejectedHeadersGammaPt[iCut]);
        }
      }
    }

    if(fIsHeavyIon == 2){
      fProfileEtaShift[iCut]          = new TProfile("Eta Shift", "Eta Shift", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      if (fMesonRecoMode < 2) fHistoNGammaConvCandidates[iCut]->Sumw2();
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger)  fHistoNGammaCaloCandidates[iCut]->Sumw2();
      if(!fDoLightOutput){
        fHistoVertexX[iCut]->Sumw2();
        fHistoVertexY[iCut]->Sumw2();
        fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
        fHistoNV0Tracks[iCut]->Sumw2();
        if (fMesonRecoMode < 2)
          fHistoConvGammaPt[iCut]->Sumw2();
        if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
          if (fHistoClusGammaPt[iCut]) fHistoClusGammaPt[iCut]->Sumw2();
          if (fHistoClusGammaE[iCut]) fHistoClusGammaE[iCut]->Sumw2();
          if (fHistoClusOverlapHeadersGammaPt[iCut]) fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
          if (fHistoClusAllHeadersGammaPt[iCut]) fHistoClusAllHeadersGammaPt[iCut]->Sumw2();
          if (fHistoClusRejectedHeadersGammaPt[iCut]) fHistoClusRejectedHeadersGammaPt[iCut]->Sumw2();
        }
      }
    }


    fHistoMotherInvMassPt[iCut]             = new TH2F("ESD_Mother_InvMass_Pt", "ESD_Mother_InvMass_Pt; M_{inv} (GeV/c^{2}); p_{T,pair} (GeV/c)",
                                                       fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
    fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
    fHistoMotherInvMassRejected[iCut]       = new TH1F("ESD_Mother_InvMassRejected", "ESD_Mother_InvMassRejected; M_{inv} (GeV/c^{2})",
                                                       1200, 0, 1.2);
    fESDList[iCut]->Add(fHistoMotherInvMassRejected[iCut]);

    if(!fDoLightOutput){
      if (fMesonRecoMode == 1){
        fHistoMotherMatchedInvMassPt[iCut]      = new TH2F("ESD_MotherMatched_InvMass_Pt", "ESD_MotherMatched_InvMass_Pt; M_{inv} (GeV/c^{2}) matched conv e^{+/-}to cluster; p_{T,pair} (GeV/c)",
                                                         fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
        fESDList[iCut]->Add(fHistoMotherMatchedInvMassPt[iCut]);
      }
    }

    fHistoMotherBackInvMassPt[iCut]         = new TH2F("ESD_Background_InvMass_Pt", "ESD_Background_InvMass_Pt; M_{inv, mxed}(GeV/c^{2}); p_{T,BG pair} (GeV/c)",
                                                       fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
    fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);

    if (fIsMC > 1){
      fHistoMotherInvMassPt[iCut]->Sumw2();
      fHistoMotherBackInvMassPt[iCut]->Sumw2();
      if(!fDoLightOutput && fMesonRecoMode == 1){
        if (fHistoMotherMatchedInvMassPt[iCut]) fHistoMotherMatchedInvMassPt[iCut]->Sumw2();
      }
    }

    if (fDoMesonQA > 0 ){
      fHistoMotherMesonPtY[iCut]              = new TH2F("ESD_MotherMeson_Pt_Y", "ESD_MotherMeson_Pt_Y; p_{T, meson cand} (GeV/c); y_{meson cand}",
                                                         nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
      fESDList[iCut]->Add(fHistoMotherMesonPtY[iCut]);
      fHistoMotherMesonPtAlpha[iCut]          = new TH2F("ESD_MotherMeson_Pt_Alpha", "ESD_MotherMeson_Pt_Alpha; p_{T, meson cand} (GeV/c); #alpha_{meson cand}",
                                                         nBinsQAPt, arrQAPtBinning, 200, -1, 1);
      fESDList[iCut]->Add(fHistoMotherMesonPtAlpha[iCut]);
      fHistoMotherMesonPtOpenAngle[iCut]      = new TH2F("ESD_MotherMeson_Pt_OpenAngle", "ESD_MotherMeson_Pt_OpenAngle; p_{T, meson cand} (GeV/c); #theta_{meson cand}",
                                                         nBinsQAPt, arrQAPtBinning, 100, 0, 1);
      fESDList[iCut]->Add(fHistoMotherMesonPtOpenAngle[iCut]);
      fHistoMotherMesonConvPhotonEtaPhi[iCut] = new TH2F("ESD_MotherMesonConvPhoton_Eta_Phi", "ConvPhoton under meson peak; #phi_{#gamma_{conv}}(rad); #eta_{#gamma_{conv}}",
                                                         600, 0, 2*TMath::Pi(), 200, -1, 1);
      fESDList[iCut]->Add(fHistoMotherMesonConvPhotonEtaPhi[iCut]);
      if (fIsMC > 1){
        fHistoMotherMesonPtY[iCut]->Sumw2();
        fHistoMotherMesonPtAlpha[iCut]->Sumw2();
        fHistoMotherMesonPtOpenAngle[iCut]->Sumw2();
        fHistoMotherMesonConvPhotonEtaPhi[iCut]->Sumw2();
      }
    }
  }

  InitBack(); // Init Background Handler

  if(fIsMC>0){
    // MC Histogramms
    fMCList                                         = new TList*[fnCuts];
    // True Histogramms
    fTrueList                                       = new TList*[fnCuts];

    if(!fDoLightOutput){
      fHistoMCHeaders                                 = new TH1I*[fnCuts];
      if (fMesonRecoMode < 2){
        fHistoTrueConvGammaPt                           = new TH1F*[fnCuts];
        fHistoDoubleCountTrueConvGammaRPt               = new TH2F*[fnCuts];
        fHistoMultipleCountTrueConvGamma                = new TH1F*[fnCuts];
        fHistoTruePrimaryConvGammaPt                    = new TH1F*[fnCuts];
      }

      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
        fHistoTrueClusGammaPt                           = new TH1F*[fnCuts];
        fHistoTrueClusConvGammaPt                       = new TH1F*[fnCuts];
        fHistoTruePrimaryClusGammaPt                    = new TH1F*[fnCuts];
        fHistoTruePrimaryClusConvGammaPt                = new TH1F*[fnCuts];
        fHistoDoubleCountTrueClusterGammaPt             = new TH2F*[fnCuts];
        fHistoMultipleCountTrueClusterGamma             = new TH1F*[fnCuts];
        fHistoTrueNLabelsInClusPt                       = new TH2F*[fnCuts];
        fHistoTrueClusConvGammaFullyPt                = new TH1F*[fnCuts];
      }
    }

    fHistoMCMesonPt                                 = new TH1F*[fnCuts];
    fHistoMCMesonWOWeightPt                         = new TH1F*[fnCuts];
    fHistoMCMesonInAccPt                            = new TH1F*[fnCuts];
    fHistoMCMesonWOWeightInAccPt                    = new TH1F*[fnCuts];
    if (fIsMC > 1){
      fHistoMCMesonWOEvtWeightPt                    = new TH1F*[fnCuts];
      fHistoMCMesonWOEvtWeightInAccPt               = new TH1F*[fnCuts];
    }

    fHistoTrueMesonInvMassPt                        = new TH2F*[fnCuts];
    fHistoDoubleCountTrueMesonInvMassPt             = new TH2F*[fnCuts];
    fHistoMultipleCountTrueMeson                    = new TH1F*[fnCuts];
    fHistoTruePrimaryMesonInvMassPt                 = new TH2F*[fnCuts];
    fHistoTruePrimaryMesonW0WeightingInvMassPt      = new TH2F*[fnCuts];
    fProfileTruePrimaryMesonWeightsInvMassPt        = new TProfile2D*[fnCuts];
    if (fMesonRecoMode == 1){
      fHistoTrueMesonMatchedInvMassPt                 = new TH2F*[fnCuts];
    }

    if (fDoMesonQA > 0){
      fHistoMCMesonPtY                              = new TH2F*[fnCuts];
      fHistoMCMesonPtAlpha                          = new TH2F*[fnCuts];
      if (fIsMC == 2){
        fHistoMCMesonPtJetPt                        = new TH2F*[fnCuts];
      }

      if (fIsMC < 2){
        if (fMesonRecoMode > 0){
          fHistoTrueMesonCaloPhotonInvMassPt                  = new TH2F*[fnCuts];
          fHistoTrueMesonCaloConvertedPhotonInvMassPt         = new TH2F*[fnCuts];
          fHistoTrueMesonCaloElectronInvMassPt                = new TH2F*[fnCuts];
          fHistoTrueMesonCaloMergedClusterInvMassPt           = new TH2F*[fnCuts];
          fHistoTrueMesonCaloMergedClusterPartConvInvMassPt   = new TH2F*[fnCuts];
        }
        if (fMesonRecoMode == 1){
          fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt  = new TH2F*[fnCuts];
          fHistoTrueMotherMesonConvPhotonEtaPhi               = new TH2F*[fnCuts];
        }
        if (fMesonRecoMode == 2){
          fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt   = new TH2F*[fnCuts];
        }
        fHistoTruePrimaryMesonMCPtResolPt                     = new TH2F*[fnCuts];
      }
      if(fDoMesonQA > 1){
        fHistoTrueBckGGInvMassPt                              = new TH2F*[fnCuts];
        fHistoTrueBckFullMesonContainedInOneClusterInvMassPt  = new TH2F*[fnCuts];
        fHistoTrueBckAsymEClustersInvMassPt                   = new TH2F*[fnCuts];
        fHistoTrueBckContInvMassPt                            = new TH2F*[fnCuts];
      }
      fHistoTrueMesonPtY                            = new TH2F*[fnCuts];
      fHistoTrueMesonPtAlpha                        = new TH2F*[fnCuts];
      fHistoTrueMesonPtOpenAngle                    = new TH2F*[fnCuts];
    }

    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringConvGamma    = "";
      TString cutstringCaloGamma    = "";
      if (fMesonRecoMode < 2) cutstringConvGamma    = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger) cutstringCaloGamma    = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      TString fullCutString         = "";
      if (fMesonRecoMode == 0)      fullCutString = Form("%i_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(), cutstringConvGamma.Data(), cutstringMeson.Data());
      else if (fMesonRecoMode == 1) fullCutString = Form("%i_%s_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(), cutstringMeson.Data());
      else if (fMesonRecoMode == 2) fullCutString = Form("%i_%s_%s_%s",fMesonRecoMode, cutstringEvent.Data(), cutstringCaloGamma.Data(), cutstringMeson.Data());

      fMCList[iCut]                     = new TList();
      fMCList[iCut]->SetName(Form("%s MC histograms",fullCutString.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);
      if(!fDoLightOutput){
        fHistoMCHeaders[iCut]             = new TH1I("MC_Headers", "MC_Headers", 20, 0, 20);
        fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
      }

      fHistoMCMesonPt[iCut]               = new TH1F("MC_Meson_Pt", "MC_Meson_Pt; p_{T} (GeV/c)",
                                                      (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
      fHistoMCMesonPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCMesonPt[iCut]);
      fHistoMCMesonWOWeightPt[iCut]       = new TH1F("MC_Meson_WOWeights_Pt", "MC_Meson_WOWeights_Pt; p_{T} (GeV/c)",
                                                      (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
      fMCList[iCut]->Add(fHistoMCMesonWOWeightPt[iCut]);

      fHistoMCMesonInAccPt[iCut]          = new TH1F("MC_MesonInAcc_Pt", "MC_MesonInAcc_Pt; p_{T} (GeV/c)",
                                                      (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
      fHistoMCMesonInAccPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCMesonInAccPt[iCut]);
      fHistoMCMesonWOWeightInAccPt[iCut]  = new TH1F("MC_MesonWOWeightInAcc_Pt", "MC_MesonWOWeightInAcc_Pt; p_{T} (GeV/c)",
                                                      (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
      fMCList[iCut]->Add(fHistoMCMesonWOWeightInAccPt[iCut]);

      if (fIsMC > 1){
        fHistoMCMesonWOWeightPt[iCut]->Sumw2();
        fHistoMCMesonWOWeightInAccPt[iCut]->Sumw2();
        fHistoMCMesonWOEvtWeightPt[iCut]  = new TH1F("MC_Meson_WOEventWeights_Pt", "MC_Meson_WOEventWeights_Pt; p_{T} (GeV/c)",
                                                      (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCMesonWOEvtWeightPt[iCut]);
        fHistoMCMesonWOEvtWeightInAccPt[iCut]  = new TH1F("MC_Meson_WOEventWeightsInAcc_Pt", "MC_Meson_WOEventWeightsInAcc_Pt; p_{T} (GeV/c)",
                                                          (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fMCList[iCut]->Add(fHistoMCMesonWOEvtWeightInAccPt[iCut]);

        if (fDoMesonQA > 0 && fIsMC == 2){
          fHistoMCMesonPtJetPt[iCut]      = new TH2F("MC_Meson_Pt_JetPt", "MC_Meson_Pt_JetPt; p_{T} (GeV/c); p_{T,jet} (GeV/c)",
                                                      nBinsQAPt, arrQAPtBinning, 200, 0, 200);
          fHistoMCMesonPtJetPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCMesonPtJetPt[iCut]);
        }
      }

      // book additional MC QA histograms
      if (fDoMesonQA > 0){
        fHistoMCMesonPtY[iCut]            = new TH2F("MC_Meson_Pt_Y", "MC_Meson_Pt_Y; p_{T} (GeV/c); Y ",
                                                      nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
        fHistoMCMesonPtY[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCMesonPtY[iCut]);
        fHistoMCMesonPtAlpha[iCut]        = new TH2F("MC_Meson_Pt_Alpha", "MC_Meson_Pt_Alpha; p_{T} (GeV/c); #alpha",
                                                      nBinsQAPt, arrQAPtBinning, 200, -1, 1);
        fMCList[iCut]->Add(fHistoMCMesonPtAlpha[iCut]);

        if (fIsMC == 2){
          fHistoMCMesonPtAlpha[iCut]->Sumw2();
        }
      }

      fTrueList[iCut]                           = new TList();
      fTrueList[iCut]->SetName(Form("%s True histograms",fullCutString.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      if(!fDoLightOutput){
        if (fMesonRecoMode < 2){
          fHistoTrueConvGammaPt[iCut]               = new TH1F("ESD_TrueConvGamma_Pt", "ESD_TrueConvGamma_Pt; p_{T,conv} (GeV/c); counts",
                                                            (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);
          fHistoDoubleCountTrueConvGammaRPt[iCut]   = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt", "ESD_TrueDoubleCountConvGamma_R_Pt; R (cm); p_{T,conv} (GeV/c)",
                                                            800, 0, 200, (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);
          fHistoMultipleCountTrueConvGamma[iCut]    = new TH1F("ESD_TrueMultipleCountConvGamma", "ESD_TrueMultipleCountConvGamma", 10, 1, 11);
          fTrueList[iCut]->Add(fHistoMultipleCountTrueConvGamma[iCut]);
          fHistoTruePrimaryConvGammaPt[iCut]        = new TH1F("ESD_TruePrimaryConvGamma_Pt", "ESD_TruePrimaryConvGamma_Pt;p_{T,conv} (GeV/c); counts", (Int_t)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fTrueList[iCut]->Add(fHistoTruePrimaryConvGammaPt[iCut]);
        }
        if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
          fHistoTrueClusGammaPt[iCut]               = new TH1F("TrueClusGamma_Pt", "ESD_TrueClusGamma_Pt; p_{T,clus} (GeV/c); counts", nBinsClusterPt, arrClusPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
          fHistoTruePrimaryClusGammaPt[iCut]        = new TH1F("TruePrimaryClusGamma_Pt", "ESD_TruePrimaryClusGamma_Pt; p_{T,clus} (GeV/c); counts", nBinsClusterPt, arrClusPtBinning);
          fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
          fHistoTrueClusConvGammaPt[iCut]           = new TH1F("TrueClusConvGamma_Pt", "TrueClusConvGamma_Pt; p_{T,clus} (GeV/c); counts", nBinsClusterPt, arrClusPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusConvGammaPt[iCut]);
          fHistoTruePrimaryClusConvGammaPt[iCut]    = new TH1F("TruePrimaryClusConvGamma_Pt", "ESD_TruePrimaryClusConvGamma_Pt; p_{T,clus} (GeV/c); counts", nBinsClusterPt, arrClusPtBinning);
          fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaPt[iCut]);
          fHistoDoubleCountTrueClusterGammaPt[iCut] = new TH2F("TrueDoubleCountClusterGamma_Pt", "TrueDoubleCountClusterGamma_Pt; p_{T,clus} (GeV/c); counts",
                                                              nBinsClusterPt, arrClusPtBinning, 2, 0, 2);
          fTrueList[iCut]->Add(fHistoDoubleCountTrueClusterGammaPt[iCut]);
          fHistoMultipleCountTrueClusterGamma[iCut] = new TH1F("TrueMultipleCountClusterGamma", "TrueMultipleCountClusterGamma", 10, 1, 11);
          fTrueList[iCut]->Add(fHistoMultipleCountTrueClusterGamma[iCut]);
          fHistoTrueClusConvGammaFullyPt[iCut]        = new TH1F("TrueClusConvGammaFullyContained_Pt", "TrueClusConvGammaFullyContained_Pt; p_{T,clus} (GeV/c); counts",
                                                                nBinsClusterPt, arrClusPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusConvGammaFullyPt[iCut]);
          fHistoTrueNLabelsInClusPt[iCut]           = new TH2F("TrueNLabelsInClus_Pt", "TrueNLabelsInClus_Pt; p_{T,clus} (GeV/c); counts",
                                                              100, -0.5, 99.5, nBinsClusterPt, arrClusPtBinning);
          fTrueList[iCut]->Add(fHistoTrueNLabelsInClusPt[iCut]);
        }
      }

      if (fIsMC > 1){
        if(!fDoLightOutput){
          if (fMesonRecoMode < 2){
            fHistoTrueConvGammaPt[iCut]->Sumw2();
            fHistoDoubleCountTrueConvGammaRPt[iCut]->Sumw2();
            fHistoMultipleCountTrueConvGamma[iCut]->Sumw2();
          }
          if (fMesonRecoMode < 2 || fEnableClusterCutsForTrigger){
            fHistoTruePrimaryConvGammaPt[iCut]->Sumw2();
          }
          if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
            fHistoTruePrimaryClusGammaPt[iCut]->Sumw2();
            fHistoTrueNLabelsInClusPt[iCut]->Sumw2();
            fHistoTrueClusConvGammaPt[iCut]->Sumw2();
            fHistoTruePrimaryClusConvGammaPt[iCut]->Sumw2();
            fHistoTrueClusGammaPt[iCut]->Sumw2();
            fHistoDoubleCountTrueClusterGammaPt[iCut]->Sumw2();
            fHistoMultipleCountTrueClusterGamma[iCut]->Sumw2();
            fHistoTrueClusConvGammaFullyPt[iCut]->Sumw2();
          }
        }

      }

      fHistoTrueMesonInvMassPt[iCut]                = new TH2F("ESD_TrueMeson_InvMass_Pt", "ESD_TrueMeson_InvMass_Pt;M_{inv}(GeV/c^{2});p_{T}(GeV/c)",
                                                              fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
      fHistoTrueMesonInvMassPt[iCut]->Sumw2();
      fTrueList[iCut]->Add(fHistoTrueMesonInvMassPt[iCut]);

      fHistoDoubleCountTrueMesonInvMassPt[iCut]     = new TH2F("ESD_TrueDoubleCountMeson_InvMass_Pt", "ESD_TrueDoubleCountMeson_InvMass_Pt;M_{inv}(GeV/c^{2});p_{T}(GeV/c)",
                                                              fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoDoubleCountTrueMesonInvMassPt[iCut]);
      fHistoMultipleCountTrueMeson[iCut]            = new TH1F("ESD_TrueMultipleCountMeson", "ESD_TrueMultipleCountMeson;# number of multiple counts;#", 10, 1, 11);
      fHistoMultipleCountTrueMeson[iCut]->Sumw2();
      fTrueList[iCut]->Add(fHistoMultipleCountTrueMeson[iCut]);

      fHistoTruePrimaryMesonInvMassPt[iCut]         = new TH2F("ESD_TruePrimaryMeson_InvMass_Pt", "ESD_TruePrimaryMeson_InvMass_Pt;M_{inv}(GeV/c^{2});p_{T}(GeV/c)",
                                                              fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
      fHistoTruePrimaryMesonInvMassPt[iCut]->Sumw2();
      fTrueList[iCut]->Add(fHistoTruePrimaryMesonInvMassPt[iCut]);

      fHistoTruePrimaryMesonW0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryMesonW0Weights_InvMass_Pt",
                                                                  "ESD_TruePrimaryMesonW0Weights_InvMass_Pt;M_{inv}(GeV/c^{2});p_{T}(GeV/c)",
                                                                  fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
      fTrueList[iCut]->Add(fHistoTruePrimaryMesonW0WeightingInvMassPt[iCut]);

      fProfileTruePrimaryMesonWeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryMesonWeights_InvMass_Pt",
                                                                        "ESD_TruePrimaryMesonWeights_InvMass_Pt;M_{inv}(GeV/c^{2});p_{T}(GeV/c)",
                                                                        fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
      fProfileTruePrimaryMesonWeightsInvMassPt[iCut]->Sumw2();
      fTrueList[iCut]->Add(fProfileTruePrimaryMesonWeightsInvMassPt[iCut]);

      if (fMesonRecoMode == 1){
        fHistoTrueMesonMatchedInvMassPt[iCut]         = new TH2F("ESD_TrueMeson_Matched_InvMass_Pt", "ESD_TrueMeson_Matched_InvMass_Pt;M_{inv}(GeV/c^{2});p_{T}(GeV/c)",
                                                                 fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
        fHistoTrueMesonMatchedInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueMesonMatchedInvMassPt[iCut]);
      }

      if (fIsMC > 1){
        fHistoTruePrimaryMesonW0WeightingInvMassPt[iCut]->Sumw2();
      }

      if (fDoMesonQA > 0){
        if (fIsMC < 2){
          if (fMesonRecoMode > 0){
            fHistoTrueMesonCaloPhotonInvMassPt[iCut]                  = new TH2F( "ESD_TrueMesonCaloPhoton_InvMass_Pt",
                                                                                  "ESD_TrueMesonCaloPhoton_InvMass_Pt; M_{inv}(GeV/c^{2}) #gamma #gamma;p_{T}(GeV/c) ",
                                                                                  fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueMesonCaloPhotonInvMassPt[iCut]);
            fHistoTrueMesonCaloConvertedPhotonInvMassPt[iCut]         = new TH2F( "ESD_TrueMesonCaloConvertedPhoton_InvMass_Pt",
                                                                                  "ESD_TrueMesonCaloConvertedPhoton_InvMass_Pt;M_{inv}(GeV/c^{2}) #gamma #gamma_{conv};p_{T}(GeV/c)",
                                                                                  fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueMesonCaloConvertedPhotonInvMassPt[iCut]);
            fHistoTrueMesonCaloElectronInvMassPt[iCut]                = new TH2F("ESD_TrueMesonCaloElectron_InvMass_Pt",
                                                                                "ESD_TrueMesonCaloElectron_InvMass_Pt; M_{inv}(GeV/c^{2}) #gamma e^{#pm}; ; p_{T}(GeV/c)",
                                                                                fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueMesonCaloElectronInvMassPt[iCut]);
            fHistoTrueMesonCaloMergedClusterInvMassPt[iCut]           = new TH2F("ESD_TrueMesonCaloMergedCluster_InvMass_Pt",
                                                                                "ESD_TrueMesonCaloMergedCluster_InvMass_Pt; M_{inv}(GeV/c^{2}) #gamma merged cluster; p_{T}(GeV/c)",
                                                                                fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueMesonCaloMergedClusterInvMassPt[iCut]);
            fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[iCut]   = new TH2F( "ESD_TrueMesonCaloMergedClusterPartConv_InvMass_Pt",
                                                                                  "ESD_TrueMesonCaloMergedClusterPartConv_InvMass_Pt; M_{inv}(GeV/c^{2}) #gamma merged cluster, part conv; p_{T}(GeV/c)",
                                                                                  fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[iCut]);
          }
          if (fMesonRecoMode == 1){
            fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt[iCut]  = new TH2F("ESD_TrueMesonCaloConvertedPhotonMatched_InvMass_Pt",
                                                                                 "ESD_TrueMesonCaloConvertedPhotonMatched_InvMass_Pt;M_{inv}(GeV/c^{2}) #gamma #gamma_{conv,matched};p_{T}(GeV/c)",
                                                                                 fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt[iCut]);
            fHistoTrueMotherMesonConvPhotonEtaPhi[iCut]               = new TH2F("ESD_TrueMotherMesonConvPhoton_Eta_Phi", "conv photons for true ; #phi_{#gamma_{conv}}(rad);#eta_{#gamma_{conv}}",
                                                                                  600, 0, 2*TMath::Pi(), 200, -1, 1);
            fTrueList[iCut]->Add(fHistoTrueMotherMesonConvPhotonEtaPhi[iCut]);
          }
          if (fMesonRecoMode == 2){
            fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt[iCut]  = new TH2F("ESD_TrueMesonCaloMixedPhotonConvertedPhoton_InvMass_Pt",
                                                                                "ESD_TrueMesonCaloMixedPhotonConvertedPhoton_InvMass_Pt;M_{inv}(GeV/c^{2}) #gamma #gamma_{conv,matched};p_{T}(GeV/c)",
                                                                                fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
            fTrueList[iCut]->Add(fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt[iCut]);
          }
          fHistoTruePrimaryMesonMCPtResolPt[iCut]                   = new TH2F("ESD_TruePrimaryMeson_MCPt_ResolPt",
                                                                              "ESD_TruePrimaryMeson_ResolPt_MCPt; p_{T,MC}(GeV/c); (p_{T,rec}-p_{T,MC})/p_{T,MC}()",
                                                                              500, 0.03, 25, 1000, -1., 1.);
          fHistoTruePrimaryMesonMCPtResolPt[iCut]->Sumw2();
          SetLogBinningXTH2(fHistoTruePrimaryMesonMCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTruePrimaryMesonMCPtResolPt[iCut]);
        }
        if(fDoMesonQA>1){
          fHistoTrueBckGGInvMassPt[iCut]                          = new TH2F("ESD_TrueBckGG_InvMass_Pt",
                                                                            "ESD_TrueBckGG_InvMass_Pt; M_{inv} (GeV/c^{2}) #gamma #gamma no signal; #pair p_{T}(GeV/c)",
                                                                            fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
          fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut] = new TH2F( "ESD_TrueBckFullMesonContained_InvMass_Pt",
                                                                                "ESD_TrueBckFullMesonContained_InvMass_Pt; M_{inv} (GeV/c^{2}) #gamma #gamma, calo gamma with full pi0; #pair p_{T}(GeV/c)",
                                                                                fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]);
          fHistoTrueBckAsymEClustersInvMassPt[iCut]   = new TH2F("ESD_TrueBckAsymEClus_InvMass_Pt",
                                                                "ESD_TrueBckAsymEClus_InvMass_Pt; M_{inv} (GeV/c^{2}) #gamma #gamma, calo gamma >70% of pi0 energy; #pair p_{T}(GeV/c)",
                                                                fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueBckAsymEClustersInvMassPt[iCut]);
          fHistoTrueBckContInvMassPt[iCut]                        = new TH2F("ESD_TrueBckCont_InvMass_Pt",
                                                                            "ESD_TrueBckCont_InvMass_Pt; M_{inv} (GeV/c^{2}) contamination; #pair p_{T}(GeV/c)",
                                                                            fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);
        }
        fHistoTrueMesonPtY[iCut]                      = new TH2F("ESD_TrueMeson_Pt_Y", "ESD_TrueMeson_Pt_Y; p_{T}(GeV/c); Y",
                                                                nBinsQAPt, arrQAPtBinning, 150, -1.5, 1.5);
        fTrueList[iCut]->Add(fHistoTrueMesonPtY[iCut]);
        fHistoTrueMesonPtAlpha[iCut]                  = new TH2F("ESD_TrueMeson_Pt_Alpha", "ESD_TrueMeson_Pt_Alpha; p_{T}(GeV/c); #alpha",
                                                                nBinsQAPt, arrQAPtBinning, 200, -1, 1);
        fTrueList[iCut]->Add(fHistoTrueMesonPtAlpha[iCut]);
        fHistoTrueMesonPtOpenAngle[iCut]              = new TH2F("ESD_TrueMeson_Pt_OpenAngle", "ESD_TrueMeson_Pt_OpenAngle; p_{T}(GeV/c); #theta",
                                                                nBinsQAPt, arrQAPtBinning, 100, 0, 1);
        fTrueList[iCut]->Add(fHistoTrueMesonPtOpenAngle[iCut]);

        if (fIsMC > 1){
          if(fDoMesonQA>1){
            fHistoTrueBckGGInvMassPt[iCut]->Sumw2();
            fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[iCut]->Sumw2();
            fHistoTrueBckAsymEClustersInvMassPt[iCut]->Sumw2();
            fHistoTrueBckContInvMassPt[iCut]->Sumw2();
          }
          fHistoTrueMesonPtY[iCut]->Sumw2();
          fHistoTrueMesonPtAlpha[iCut]->Sumw2();
          fHistoTrueMesonPtOpenAngle[iCut]->Sumw2();
          if (fHistoTrueMotherMesonConvPhotonEtaPhi[iCut]) fHistoTrueMotherMesonConvPhotonEtaPhi[iCut]->Sumw2();
        }
      }
    }
  }

  fVectorDoubleCountTrueMesons.clear();
  fVectorDoubleCountTrueConvGammas.clear();
  fVectorDoubleCountTrueClusterGammas.clear();

  fMapMultipleCountTrueMesons.clear();
  fMapMultipleCountTrueConvGammas.clear();
  fMapMultipleCountTrueClusterGammas.clear();

  fVectorRecTrueMesons.clear();

  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  if( fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
    if (fV0Reader->GetV0FindingEfficiencyHistograms())
      fOutputContainer->Add(fV0Reader->GetV0FindingEfficiencyHistograms());

  for(Int_t iMatcherTask = 0; iMatcherTask < 3; iMatcherTask++){
    AliCaloTrackMatcher* temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i",iMatcherTask)));
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if (fMesonRecoMode < 2){
      if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
      if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
      }
    }
    if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
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
    }
    if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
    }
  }

  if (fIsMC > 0){
    fTreeBrokenFiles = new TTree("BrokenFiles", "BrokenFiles");
    fTreeBrokenFiles->Branch("fileName",&fFileNameBroken);
    fOutputContainer->Add(fTreeBrokenFiles);
  }

  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskHeavyNeutralMesonToGG::Notify(){

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
void AliAnalysisTaskHeavyNeutralMesonToGG::UserExec(Option_t *){
    //
    // Called for each event
    //
    fInputEvent = InputEvent();

    if(fIsMC > 0) fMCEvent = MCEvent();

    Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
    if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;// incomplete event
    // Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete => abort processing of this event/file
    if(eventQuality == 2 || eventQuality == 3){
      // write out name of broken file for first event
      if (fIsMC > 0){
        if (fInputEvent->IsA()==AliESDEvent::Class()){
          if (((AliESDEvent*)fInputEvent)->GetEventNumberInFile() == 0){
            fFileNameBroken = new TObjString(Form("%s",((TString)fV0Reader->GetCurrentFileName()).Data()));
            if (fTreeBrokenFiles) fTreeBrokenFiles->Fill();
            delete fFileNameBroken;
          }
        }
      }

      for(Int_t iCut = 0; iCut<fnCuts; iCut++){
        fHistoNEvents[iCut]->Fill(eventQuality);
        if (fIsMC > 1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
      }
      return;
    }

//     if(fInputEvent->IsA()==AliAODEvent::Class()){
//       fInputEvent->InitMagneticField();
//     }

    fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

    // ------------------- BeginEvent ----------------------------
    AliEventplane *EventPlane = fInputEvent->GetEventplane();
    if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
    else fEventPlaneAngle=0.0;

    if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled()) && fMesonRecoMode < 2){
      RelabelAODPhotonCandidates(kTRUE);// In case of AODMC relabeling MC
      fV0Reader->RelabelAODs(kTRUE);
    }

    for(Int_t iCut = 0; iCut<fnCuts; iCut++){

      fiCut = iCut;
      Bool_t isRunningEMCALrelAna = kFALSE;
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
        if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
      }

      Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,isRunningEMCALrelAna);

      if(fIsMC==2){
        Float_t xsection      = -1.;
        Float_t ntrials       = -1.;
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials,fInputEvent);
        if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
        fProfileJetJetXSection[iCut]->Fill(0.,xsection);
        fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
      }

      if(fIsMC>0){
        fWeightJetJetMC       = 1;
        //     cout << fMCEvent << endl;
        Float_t pthard = -1;
        Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC, pthard, fInputEvent );
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

      if(fMesonRecoMode < 2){
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetDoElecDeDxPostCalibration()){
            if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->LoadElecDeDxPostCalibration(fInputEvent->GetRunNumber())){
            AliFatal(Form("ERROR: LoadElecDeDxPostCalibration returned kFALSE for %d despite being requested!",fInputEvent->GetRunNumber()));
            }
        }
      }

      if (triggered==kTRUE){
        fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

        fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
        fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
        if(!fDoLightOutput){
          fHistoVertexX[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetX(), fWeightJetJetMC);
          fHistoVertexY[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetY(), fWeightJetJetMC);
          fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)),fWeightJetJetMC);
          if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)  fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(), fWeightJetJetMC);
          else fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(), fWeightJetJetMC);
        }
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

          if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader()){
            for(Int_t i = 0;i<(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
              TString nameBin= fHistoMCHeaders[iCut]->GetXaxis()->GetBinLabel(i+1);
              if (nameBin.CompareTo("")== 0){
                TString nameHeader = ((TObjString*)((TList*)((AliConvEventCuts*)fEventCutArray->At(iCut))
                ->GetAcceptedHeader())->At(i))->GetString();
                fHistoMCHeaders[iCut]->GetXaxis()->SetBinLabel(i+1,nameHeader.Data());
              }
            }
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
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger) // process calo clusters
        ProcessClusters();
      if (fMesonRecoMode < 2) // Process this cuts gammas
        ProcessPhotonCandidates();

      if (fMesonRecoMode < 2) fHistoNGammaConvCandidates[iCut]->Fill(fGammaCandidates->GetEntries(),fWeightJetJetMC);
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger) fHistoNGammaCaloCandidates[iCut]->Fill(fClusterCandidates->GetEntries(),fWeightJetJetMC);

      if (fMesonRecoMode < 2){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC>0){
          fUnsmearedPx = new Double_t[fGammaCandidates->GetEntries()]; // Store unsmeared Momenta
          fUnsmearedPy = new Double_t[fGammaCandidates->GetEntries()];
          fUnsmearedPz = new Double_t[fGammaCandidates->GetEntries()];
          fUnsmearedE =  new Double_t[fGammaCandidates->GetEntries()];

          for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
            fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Px();
            fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Py();
            fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Pz();
            fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->E();
            ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(gamma)));
          }
        }
      }

      // check gamma gamma pairs and veto if necessary
      if (!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetAcceptMassFlag())
        SetPhotonVeto();

      CalculateMesonCandidates(); // Combine Gammas from conversion and from calo

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          CalculateBackground(); // Combinatorial Background
          UpdateEventByEventData(); // Store Event for mixed Events
        } else if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 2){
          CalculateBackgroundSwapp(); // Combinatorial Background
        } else {
          CalculateBackgroundRP(); // Combinatorial Background
          fBGHandlerRP[iCut]->AddEvent(fGammaCandidates,fInputEvent); // Store Event for mixed Events
          fBGClusHandlerRP[iCut]->AddEvent(fClusterCandidates,fInputEvent); // Store Event for mixed Events
        }
      }

      if (fMesonRecoMode < 2){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC>0){
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
      }

      if(fIsMC>0){
        fVectorRecTrueMesons.clear();
        fVectorDoubleCountTrueMesons.clear();
        FillMultipleCountHistoAndClear(fMapMultipleCountTrueMesons,fHistoMultipleCountTrueMeson[iCut]);
        fVectorDoubleCountTrueConvGammas.clear();
        if(!fDoLightOutput && fMesonRecoMode < 2) FillMultipleCountHistoAndClear(fMapMultipleCountTrueConvGammas,fHistoMultipleCountTrueConvGamma[iCut]);
        fVectorDoubleCountTrueClusterGammas.clear();
        if(!fDoLightOutput && fMesonRecoMode > 0) FillMultipleCountHistoAndClear(fMapMultipleCountTrueClusterGammas,fHistoMultipleCountTrueClusterGamma[iCut]);
      }

      fGammaCandidates->Clear(); // delete this cuts good gammas
      fClusterCandidates->Clear(); // delete cluster candidates
    }

    if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled()) && fMesonRecoMode < 2){
      RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
      fV0Reader->RelabelAODs(kFALSE);
    }

    PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessClusters(){
  Int_t nclus = 0;
  TClonesArray * arrClustersProcess = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskHeavyNeutralMesonToGG! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersProcess->GetEntries();
  }

  vector<AliAODConversionPhoton*>         vectorCurrentClusters;
  vector<Int_t>                           vectorRejectCluster;
  vector<Double_t>                        vectorPhotonWeight;
  vector<Double_t>                        vectorClusterM02;

  //   cout << nclus << endl;

  if(nclus == 0)  return;

  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  if(fDoPrimaryTrackMatching) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kFALSE);

  // vertex
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Int_t i = 0; i < nclus; i++){
    AliVCluster* clus = NULL;
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

    if (!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)){
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
    PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent));

    // get MC label
    if(fIsMC>0){
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
      //       cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k<(Int_t)clus->GetNLabels(); k++){
          if (k<50)PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
          //           Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
          //           cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
        }
      }
    }

    fIsFromDesiredHeader            = kTRUE;
    fIsOverlappingWithOtherHeader   = kFALSE;
    //TString periodName         = fV0Reader->GetPeriodName();
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0){
        fIsFromDesiredHeader = kFALSE;
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
    if (fIsMC > 0 && fHistoClusAllHeadersGammaPt[fiCut])
      fHistoClusAllHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
    if (!fIsFromDesiredHeader && fIsMC > 0 && fHistoClusRejectedHeadersGammaPt[fiCut])
      fHistoClusRejectedHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
    if (fIsFromDesiredHeader && fIsOverlappingWithOtherHeader && fIsMC > 0 && fHistoClusOverlapHeadersGammaPt[fiCut])
      fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);

    if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders)){

      vectorCurrentClusters.push_back(PhotonCandidate);
      vectorPhotonWeight.push_back(fWeightJetJetMC);
      vectorClusterM02.push_back(clus->GetM02());
    }else{
      delete PhotonCandidate;
    }

    delete clus;
    delete tmpvec;
  }

  //Bool_t rejected = kFALSE;
  // run conversion recovery in addition
  if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetIsConversionRecovery()){
    /*rejected = */((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckForReconstructedConversionPairs(vectorCurrentClusters,vectorRejectCluster);
  }

  for (Int_t iter = 0; iter < (Int_t)vectorCurrentClusters.size();iter++){
    if (!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckVectorForIndexAndAdd(vectorRejectCluster, iter,kFALSE)){
      fHistoClusGammaPt[fiCut]->Fill(vectorCurrentClusters.at(iter)->Pt(), vectorPhotonWeight.at(iter));
      fHistoClusGammaE[fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
      if(fIsMC> 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ProcessTrueClusterCandidates(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        } else {
          ProcessTrueClusterCandidatesAOD(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        }
      }
      fClusterCandidates->Add(vectorCurrentClusters.at(iter));
    }
  }
  vectorRejectCluster.clear();
  vectorPhotonWeight.clear();
  vectorClusterM02.clear();
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate, Float_t clusM02){

  TParticle *Photon = NULL;
  if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
  if (TruePhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return;
  if(!fDoLightOutput) fHistoTrueNLabelsInClusPt[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(),TruePhotonCandidate->Pt(),fWeightJetJetMC);

  if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
  else return;

  if(Photon == NULL){
    //    cout << "no photon" << endl;
    return;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, fEnableSortForClusMC);

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX           = primVtxMC->GetX();
  Double_t mcProdVtxY           = primVtxMC->GetY();
  Double_t mcProdVtxZ           = primVtxMC->GetZ();
  Bool_t isPrimary              = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

  // to get primary distrinction right put mother of conversion electron as particle to check
  if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
    isPrimary              = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, Photon->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);


  // Fill histograms for inclusive gamma corrections
  // --> a) all clusters with leading real or converted photons
  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()) ){
    fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
    if(!fDoLightOutput){
      // how many of those clusters are from converted photons
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      }
      // --> b) which of these are primary
      if(isPrimary){
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        // --> c) which are from conversions? Just additonal information
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && (Photon->GetMother(0)>-1)){
          fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        }
        // --> d) how do the secondaries look like
      }
    }
  }

  // Some additional QA
  if (fDoClusterQA > 0){
    // how many of the converted photons are fully contained in the cluster
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained())
      fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
  }

  // Check if we are double counting photons
  Int_t motherLab = Photon->GetMother(0);
  if (motherLab > -1){
    if (TruePhotonCandidate->IsLargestComponentPhoton()){
      if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
        fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,fWeightJetJetMC);
        FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
      }
    }
    Int_t grandMotherLab = fMCEvent->Particle(motherLab)->GetMother(0);
    if (grandMotherLab > -1){
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
          fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,fWeightJetJetMC);
          FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,grandMotherLab);
        }
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate, Float_t clusM02){

  AliAODMCParticle *Photon = NULL;
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray){
    if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
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
  if(!fDoLightOutput) fHistoTrueNLabelsInClusPt[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels(),TruePhotonCandidate->Pt(),fWeightJetJetMC);

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

  // to get primary distrinction right put mother of conversion electron as particle to check
  if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
    if (Photon->GetMother()> -1){
      AliAODMCParticle *Mother  = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
      isPrimary                 = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, Mother, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    }
  }

  // True Photon
  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()) ){
    fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
    if(!fDoLightOutput){
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
      }
      if(isPrimary){
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
        if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
        }
      }
    }
  }
  if (fDoClusterQA > 0){
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained())
      fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt(), fWeightJetJetMC);
  }
  Int_t motherLab = Photon->GetMother();
  if (motherLab > -1){
    if (TruePhotonCandidate->IsLargestComponentPhoton()){
      if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,motherLab)){
        fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)0,fWeightJetJetMC);
        FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,motherLab);
      }
    }
    Int_t grandMotherLab = ((AliAODMCParticle*) AODMCTrackArray->At(motherLab))->GetMother();
    if (grandMotherLab > -1){
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueClusterGammas,grandMotherLab)){
          fHistoDoubleCountTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),(Double_t)1,fWeightJetJetMC);
          FillMultipleCountMap(fMapMultipleCountTrueClusterGammas,grandMotherLab);
        }
      }
    }
  }

  return;
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessPhotonCandidates(){

  Int_t nV0 = 0;
  TList *GammaCandidatesStepOne = new TList();
  TList *GammaCandidatesStepTwo = new TList();
  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    fIsFromDesiredHeader = kTRUE;
    if(fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
      if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
    }

    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
        !((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

      if(fIsFromDesiredHeader){
        if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
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
      fIsFromDesiredHeader = kTRUE;
      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
        fGammaCandidates->Add(PhotonCandidate);
        if(fIsFromDesiredHeader){
          if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
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
      fIsFromDesiredHeader = kTRUE;
      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
      fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
      if(fIsFromDesiredHeader){
        if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
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
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
  AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};

  if(posDaughter->GetMother() != negDaughter->GetMother()){
    return;
  }  else if(posDaughter->GetMother() == -1){
    return;
  }

  if(pdgCode[0]!=11 || pdgCode[1]!=11){
    return; //One Particle is not a electron
  }

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()){
    return; // Same Charge
  }

  AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
  if(Photon->GetPdgCode() != 22){
    return; // Mother is no Photon
  }

  if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5){
    return;// check if the daughters come from a conversion
  }
  // STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX

  // True Photon
  if(fIsFromDesiredHeader){
    if(!fDoLightOutput) fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
    if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother())){
      if(!fDoLightOutput) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
      FillMultipleCountMap(fMapMultipleCountTrueConvGammas,posDaughter->GetMother());
    }
  }

  Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if(isPrimary){
    // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
    if(fIsFromDesiredHeader){
      if(!fDoLightOutput){
        fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
      }
    }
    // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
  }
  TruePhotonCandidate->SetIsTrueConvertedPhoton();
  return;
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Photons
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};
  if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){
    return;
  }
  else if(posDaughter->GetMother(0) == -1){
    return;
  }

  if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

  TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);

  if(Photon->GetPdgCode() != 22){
    return; // Mother is no Photon
  }

  if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

  // True Photon
  if(fIsFromDesiredHeader){
    if(!fDoLightOutput) fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
    if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother(0))){
      if(!fDoLightOutput) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt(),fWeightJetJetMC);
      FillMultipleCountMap(fMapMultipleCountTrueConvGammas,posDaughter->GetMother(0));
    }
  }
  Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if(isPrimary){
    // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
    if(fIsFromDesiredHeader){
      if(!fDoLightOutput){
        fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),fWeightJetJetMC);
      }
    }
    // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
  }
  TruePhotonCandidate->SetIsTrueConvertedPhoton();
  return;
}
//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessAODMCParticles(){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {

    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
    if (!particle) continue;

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary) {

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }

      if (fMesonRecoMode < 2){
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      }

      Double_t mesonY = 1.e30;
      Double_t ratio  = 0;
      if (particle->E() != TMath::Abs(particle->Pz())){
        ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
      }
      if( !(ratio <= 0) ){
        mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
      }

      // check neutral mesons
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughterLabel(0)));
        AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughterLabel(1)));
        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
          }
        }
        Double_t alpha = -10;

        if(particle->GetPdgCode() == fMesonPDG){
          alpha = (daughter0->E() - daughter1->E())/(daughter0->E() + daughter1->E());
          if (fHistoMCMesonPt[fiCut])fHistoMCMesonPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Meson
          if (fHistoMCMesonWOWeightPt[fiCut]) fHistoMCMesonWOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          if (fIsMC > 1 && fHistoMCMesonWOEvtWeightPt[fiCut])fHistoMCMesonWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0){
            if (fHistoMCMesonPtY[fiCut])fHistoMCMesonPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC);
            if (fHistoMCMesonPtAlpha[fiCut])fHistoMCMesonPtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC);
            if (fIsMC == 2 && fHistoMCMesonPtJetPt[fiCut])fHistoMCMesonPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }

          // Check the acceptance for both gammas
          if (fMesonRecoMode == 0){
            if( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
            // check acceptance of clusters as well, true if one of them points into the Calo acceptance
              if (fHistoMCMesonInAccPt[fiCut])fHistoMCMesonInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Meson with gamma in acc
              if (fHistoMCMesonWOWeightInAccPt[fiCut])fHistoMCMesonWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Meson with gamma in acc wo weight
              if(fIsMC > 1 && fHistoMCMesonWOEvtWeightInAccPt[fiCut])fHistoMCMesonWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Meson with gamma in acc wo any weight
            }
          } else if (fMesonRecoMode == 1){
            if( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
                ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
              // check acceptance of clusters as well, true if one of them points into the Calo acceptance
              if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,AODMCTrackArray) ||
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,AODMCTrackArray) ){
                if (fHistoMCMesonInAccPt[fiCut])fHistoMCMesonInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Meson with gamma in acc
                if (fHistoMCMesonWOWeightInAccPt[fiCut])fHistoMCMesonWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Meson with gamma in acc wo weight
                if(fIsMC > 1 && fHistoMCMesonWOEvtWeightInAccPt[fiCut])fHistoMCMesonWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Meson with gamma in acc wo any weight
              }
            }
          } else if (fMesonRecoMode == 2){
            if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,AODMCTrackArray) &&
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,AODMCTrackArray) ){
              if (fHistoMCMesonInAccPt[fiCut])fHistoMCMesonInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Meson with gamma in acc
              if (fHistoMCMesonWOWeightInAccPt[fiCut])fHistoMCMesonWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Meson with gamma in acc wo weight
              if(fIsMC > 1 && fHistoMCMesonWOEvtWeightInAccPt[fiCut])fHistoMCMesonWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Meson with gamma in acc wo any weight
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessMCParticles(){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  //   cout << mcProdVtxX <<"\t" << mcProdVtxY << "\t" << mcProdVtxZ << endl;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      // fill primary histograms
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }

      if (fMesonRecoMode < 2){
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      }

      // Fill histograms for other particles
      Double_t mesonY = 1.e30;
      Double_t ratio  = 0;
      if (particle->Energy() != TMath::Abs(particle->Pz())){
        ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
      }
      if( !(ratio <= 0) ){
        mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
      }

      // check neutral mesons
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        TParticle* daughter0 = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
        TParticle* daughter1 = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
          }
        }

        Double_t alpha = -10;
        if(particle->GetPdgCode() == fMesonPDG){
          alpha = (daughter0->Energy() - daughter1->Energy())/(daughter0->Energy() + daughter1->Energy());
          if (fHistoMCMesonPt[fiCut])fHistoMCMesonPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // All MC Meson
          if (fHistoMCMesonWOWeightPt[fiCut]) fHistoMCMesonWOWeightPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC);
          if (fIsMC > 1 && fHistoMCMesonWOEvtWeightPt[fiCut])fHistoMCMesonWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0){
            if (fHistoMCMesonPtY[fiCut])fHistoMCMesonPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted*fWeightJetJetMC); // All MC Meson
            if (fHistoMCMesonPtAlpha[fiCut])fHistoMCMesonPtAlpha[fiCut]->Fill(particle->Pt(),alpha,fWeightJetJetMC); // All MC Meson
            if (fIsMC == 2 && fHistoMCMesonPtJetPt[fiCut])fHistoMCMesonPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
          }

          // Check the acceptance for both gammas & whether they are counted as primaries as well
          Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          if( kDaughter0IsPrim && kDaughter1IsPrim ){
            if (fMesonRecoMode == 0){
              if( ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)  &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
                // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                if (fHistoMCMesonInAccPt[fiCut])fHistoMCMesonInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Meson with gamma in acc
                if (fHistoMCMesonWOWeightInAccPt[fiCut])fHistoMCMesonWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Meson with gamma in acc wo weighting
                if(fIsMC > 1 && fHistoMCMesonWOEvtWeightInAccPt[fiCut])fHistoMCMesonWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Meson with gamma in acc wo any weight
              }
            } else if (fMesonRecoMode == 1){
              if (((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCEvent,kFALSE) &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCEvent,kFALSE)  &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
                // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) ||
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ){
                  if (fHistoMCMesonInAccPt[fiCut])fHistoMCMesonInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Meson with gamma in acc
                  if (fHistoMCMesonWOWeightInAccPt[fiCut])fHistoMCMesonWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Meson with gamma in acc wo weighting
                  if(fIsMC > 1 && fHistoMCMesonWOEvtWeightInAccPt[fiCut])fHistoMCMesonWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Meson with gamma in acc wo any weight
                }
              }
            } else if (fMesonRecoMode == 2){
              if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ){
                if (fHistoMCMesonInAccPt[fiCut])fHistoMCMesonInAccPt[fiCut]->Fill(particle->Pt(),weighted*fWeightJetJetMC); // MC Meson with gamma in acc
                if (fHistoMCMesonWOWeightInAccPt[fiCut])fHistoMCMesonWOWeightInAccPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // MC Meson with gamma in acc wo weighting
                if(fIsMC > 1 && fHistoMCMesonWOEvtWeightInAccPt[fiCut])fHistoMCMesonWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Meson with gamma in acc wo any weight
              }
            }
          }
        }
      }
    }
  }
}


//________________________________________________________________________
// function to reject photons in specific invariant mass window
//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::SetPhotonVeto(){
  TClonesArray * arrClustersMesonCand = NULL;
  if(fCorrTaskSetting.CompareTo("") && fMesonRecoMode > 0)
    arrClustersMesonCand = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));

  if (fMesonRecoMode == 0){
    if(fGammaCandidates->GetEntries()>1){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
          //Check for same Electron ID
          if (gamma1==NULL) continue;
          if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
            gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
            gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
            gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;

          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);
          mesonCand->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if (!(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(mesonCand, 0))){
              gamma0->SetUseForMesonPair(kFALSE);
              gamma1->SetUseForMesonPair(kFALSE);
              fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
            }
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  } else if (fMesonRecoMode == 1){
    if(fGammaCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;

        for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          Bool_t matched = kFALSE;
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL) continue;

          if (gamma1->GetIsCaloPhoton() > 0){
            AliVCluster* cluster = NULL;
            if(fInputEvent->IsA()==AliESDEvent::Class()){
              if(arrClustersMesonCand)
                cluster = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
              else
                cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            } else if(fInputEvent->IsA()==AliAODEvent::Class()){
              if(arrClustersMesonCand)
                cluster = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
              else
                cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            }

            matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC);
            if(arrClustersMesonCand) delete cluster;
          }

          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if (!matched){
              if (!(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(mesonCand, 0))){
                gamma0->SetUseForMesonPair(kFALSE);
                gamma1->SetUseForMesonPair(kFALSE);
                fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
              }
            }
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  } else if (fMesonRecoMode == 2){
    if(fClusterCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL) continue;

          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand, kTRUE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), gamma0->GetLeadingCellID(), gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton()))){
            if (!(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(mesonCand, 0))){
              gamma0->SetUseForMesonPair(kFALSE);
              gamma1->SetUseForMesonPair(kFALSE);
              fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
            }
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::CalculateMesonCandidates(){
  TClonesArray * arrClustersMesonCand = NULL;
  if(fCorrTaskSetting.CompareTo("") && fMesonRecoMode > 0)
    arrClustersMesonCand = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));

  if (fMesonRecoMode == 0){
    if(fGammaCandidates->GetEntries()>1){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
          //Check for same Electron ID
          if (gamma1==NULL) continue;
          if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
            gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
            gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
            gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;


          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);
          mesonCand->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if (!gamma0->GetUseForMesonPair() || !gamma1->GetUseForMesonPair()){
              fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
              delete mesonCand;
              mesonCand=0x0;
              continue;
            }
            if (fHistoMotherInvMassPt[fiCut]) fHistoMotherInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);

            if (fDoMesonQA > 0){
              if (mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1]){
                  if (fHistoMotherMesonPtY[fiCut]) fHistoMotherMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  if (fHistoMotherMesonPtAlpha[fiCut]) fHistoMotherMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(),fWeightJetJetMC);
                  if (fHistoMotherMesonPtOpenAngle[fiCut]) fHistoMotherMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(),fWeightJetJetMC);
                }
            }
            if(fDoTHnSparse && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
              Int_t zbin = 0;
              Int_t mbin = 0;
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
                zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                }else {
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                }
              }else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() != 2){
                zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                  mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                }else {
                  mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                }
              }

              Double_t sparesFill[4] = {mesonCand->M(),mesonCand->Pt(),(Double_t)zbin,(Double_t)mbin};
              if (fSparseMotherInvMassPtZM[fiCut]) fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
            }


            if( fIsMC > 0 ){
              if(fInputEvent->IsA()==AliESDEvent::Class())
                ProcessTrueMesonCandidatesConv(mesonCand, gamma0, gamma1);
              if(fInputEvent->IsA()==AliAODEvent::Class())
                ProcessTrueMesonCandidatesConvAOD(mesonCand, gamma0, gamma1);
            }
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  } else if (fMesonRecoMode == 1){
    if(fGammaCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;

        for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          Bool_t matched = kFALSE;
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL) continue;

          if (gamma1->GetIsCaloPhoton() > 0){
            AliVCluster* cluster = NULL;
            if(fInputEvent->IsA()==AliESDEvent::Class()){
              if(arrClustersMesonCand)
                cluster = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
              else
                cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            } else if(fInputEvent->IsA()==AliAODEvent::Class()){
              if(arrClustersMesonCand)
                cluster = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
              else
                cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
            }

            matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC);
            if(arrClustersMesonCand) delete cluster;
          }

          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if (matched){
              if(!fDoLightOutput && fMesonRecoMode == 1){
                if(fHistoMotherMatchedInvMassPt[fiCut]) fHistoMotherMatchedInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
              }
            }else {
              if (!gamma0->GetUseForMesonPair() || !gamma1->GetUseForMesonPair()){
                fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
                delete mesonCand;
                mesonCand=0x0;
                continue;
              }
              if (fHistoMotherInvMassPt[fiCut]) fHistoMotherInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
            }

            // fill new histograms
            if (!matched){
              if (fDoMesonQA > 0){
                if (mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1]){
                  if (fHistoMotherMesonPtY[fiCut]) fHistoMotherMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  if (fHistoMotherMesonPtAlpha[fiCut]) fHistoMotherMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(),fWeightJetJetMC);
                  if (fHistoMotherMesonPtOpenAngle[fiCut]) fHistoMotherMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(),fWeightJetJetMC);
                  if (fHistoMotherMesonConvPhotonEtaPhi[fiCut]) fHistoMotherMesonConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
                }
              }
              if(fDoTHnSparse && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
                Int_t zbin = 0;
                Int_t mbin = 0;

                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
                  zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                    mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                  }else {
                    mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                  }
                }else if (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() != 2){
                  zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                    mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                  }else {
                    mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                  }
                }
                Double_t sparesFill[4] = {mesonCand->M(),mesonCand->Pt(),(Double_t)zbin,(Double_t)mbin};
                if (fSparseMotherInvMassPtZM[fiCut]) fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
              }
            }

            if(fIsMC>0){
              if(fInputEvent->IsA()==AliESDEvent::Class())
                ProcessTrueMesonCandidatesConvCalo(mesonCand, gamma0, gamma1, matched);
              if(fInputEvent->IsA()==AliAODEvent::Class())
                ProcessTrueMesonCandidatesConvCaloAOD(mesonCand, gamma0, gamma1, matched);
            }
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  } else if (fMesonRecoMode == 2){
    if(fClusterCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL) continue;

          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand, kTRUE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), gamma0->GetLeadingCellID(), gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton()))){
            if (!gamma0->GetUseForMesonPair() || !gamma1->GetUseForMesonPair()){
              fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
              delete mesonCand;
              mesonCand=0x0;
              continue;
            }
            if (fHistoMotherInvMassPt[fiCut]) fHistoMotherInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
            // fill new histograms
            if (fDoMesonQA > 0 && fDoMesonQA < 3){
              if ( mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1]){
                fHistoMotherMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
                fHistoMotherMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(), fWeightJetJetMC);
                fHistoMotherMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(), fWeightJetJetMC);
              }
            }
            if(fDoTHnSparse && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
              Int_t zbin = 0;
              Int_t mbin = 0;

              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
                zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                } else {
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
                }
              }
              Double_t sparesFill[4] = {mesonCand->M(),mesonCand->Pt(),(Double_t)zbin,(Double_t)mbin};
              if (fSparseMotherInvMassPtZM[fiCut]) fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
            }

            if(fIsMC> 0){
              if(fInputEvent->IsA()==AliESDEvent::Class())
                ProcessTrueMesonCandidatesCalo(mesonCand,gamma0,gamma1);
              if(fInputEvent->IsA()==AliAODEvent::Class())
                ProcessTrueMesonCandidatesCaloAOD(mesonCand,gamma0,gamma1);
            }
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueMesonCandidatesConvCalo(  AliAODConversionMother *mesonCand,
                                                                                AliAODConversionPhoton *TrueGammaCandidate0,
                                                                                AliAODConversionPhoton *TrueGammaCandidate1,
                                                                                Bool_t matched )
{
  // obtain MC vertex
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTrueMeson = kFALSE;
    Int_t gamma0MCLabel = -1;
    Int_t gamma0MotherLabel = -1;
    if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
      gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
      if(gamma0MCLabel>-1){
        TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
        gamma0MotherLabel=gammaMC0->GetFirstMother();
      }
    }
    if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

    Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
    Int_t gamma1MotherLabel = -1;
    // check if

    TParticle * gammaMC1 = 0x0;
    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
      if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
        // get mother of interest (pi0 or eta)
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
          if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }
    }

    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fMesonPDG){
        isTrueMeson=kTRUE;
      }
    }

    if(isTrueMeson ){// True Pion or Eta
      if (!matched){
        if (fHistoTrueMesonInvMassPt[fiCut]) fHistoTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
      }else{
        if (fHistoTrueMesonMatchedInvMassPt[fiCut]) fHistoTrueMesonMatchedInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
      }
      if (fDoMesonQA > 0  && fIsMC < 2){
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && !matched){
          if (fHistoTrueMesonCaloPhotonInvMassPt[fiCut]) fHistoTrueMesonCaloPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        }
        if (TrueGammaCandidate1->IsLargestComponentElectron() && !matched){
          if (fHistoTrueMesonCaloElectronInvMassPt[fiCut]) fHistoTrueMesonCaloElectronInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        }
        if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() ){
          if ( !matched){
            fHistoTrueMesonCaloConvertedPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
          }

          if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTrueMeson){
            if (fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt[fiCut])fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
          }
        }
        if ((TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged()) && !matched ){
          if (fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        }
        if (TrueGammaCandidate1->IsMergedPartConv() && !matched){
          if (fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        }
      }
      if (!matched){
        if (fDoMesonQA > 0){
          if (isTrueMeson){
            if (mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1] ){
              if (fHistoTrueMesonPtY[fiCut]) fHistoTrueMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
              if (fHistoTrueMesonPtAlpha[fiCut]) fHistoTrueMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(),fWeightJetJetMC);
              if (fHistoTrueMesonPtOpenAngle[fiCut]) fHistoTrueMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(),fWeightJetJetMC);
              if (fHistoTrueMotherMesonConvPhotonEtaPhi[fiCut])fHistoTrueMotherMesonConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta(),fWeightJetJetMC);
            }
          }
        }
        Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if (isPrimary) {
          Float_t weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
              if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
                weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
                  //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
              }
          }
          if (isTrueMeson){
              if (fHistoTruePrimaryMesonInvMassPt[fiCut]) fHistoTruePrimaryMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
              if (fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]) fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
              if (fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]) fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
              if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMesons,gamma0MotherLabel)){
                if (fHistoDoubleCountTrueMesonInvMassPt[fiCut]) fHistoDoubleCountTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
                FillMultipleCountMap(fMapMultipleCountTrueMesons,gamma0MotherLabel);
              }
          }

          if (fDoMesonQA > 0 && fIsMC < 2){
            if(isTrueMeson){ // Only primary pi0 for resolution
              if (fHistoTruePrimaryMesonMCPtResolPt[fiCut]) fHistoTruePrimaryMesonMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(mesonCand->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted*fWeightJetJetMC);
            }
          }
        }
      }
    }else if(!isTrueMeson ){ // Background
      if (fDoMesonQA > 1){
        if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Meson or Eta
          if (fHistoTrueBckGGInvMassPt[fiCut])fHistoTrueBckGGInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
          if(
            ( ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fMesonPDG
              && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
          ){
            if (fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]) fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
          }else if( TrueGammaCandidate1->E()/mesonCand->E() > 0.7 ){
            if (fHistoTrueBckAsymEClustersInvMassPt[fiCut]) fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
          }
        }else { // No photon or without mother
          if (fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        }
      }
    }
    if (isTrueMeson && !matched){
      fVectorRecTrueMesons.push_back(gamma0MotherLabel);
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueMesonCandidatesConvCaloAOD(   AliAODConversionMother *mesonCand,
                                                                                    AliAODConversionPhoton *TrueGammaCandidate0,
                                                                                    AliAODConversionPhoton *TrueGammaCandidate1,
                                                                                    Bool_t matched )
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  Bool_t isTrueMeson = kFALSE;

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

  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");
  Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); // get most probable MC label
  Int_t gamma1MotherLabel = -1;
  // check if

  AliAODMCParticle * gammaMC1 = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother();
      }else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
          gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
        }else gamma1MotherLabel=gammaMC1->GetMother();
      }
    }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fMesonPDG){
      isTrueMeson=kTRUE;
    }
  }

  if(isTrueMeson ){// True Pion or Eta
    if (!matched){
      if (fHistoTrueMesonInvMassPt[fiCut]) fHistoTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);

    }else{
      if (fHistoTrueMesonMatchedInvMassPt[fiCut]) fHistoTrueMesonMatchedInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
    }
    if (fDoMesonQA > 0 && fIsMC < 2){
      if (TrueGammaCandidate1->IsLargestComponentPhoton() && !matched){
        if (fHistoTrueMesonCaloPhotonInvMassPt[fiCut]) fHistoTrueMesonCaloPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      if (TrueGammaCandidate1->IsLargestComponentElectron() && !matched) {
        if (fHistoTrueMesonCaloElectronInvMassPt[fiCut]) fHistoTrueMesonCaloElectronInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()){
        if (!matched)fHistoTrueMesonCaloConvertedPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTrueMeson)
            if (fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt[fiCut])fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      if ((TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged()) && !matched ){
        if (fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      if (TrueGammaCandidate1->IsMergedPartConv() && !matched) {
        if (fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
    }

    if ( !matched){
      if (fDoMesonQA > 0){
        if (mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1] ){
          if (fHistoTrueMesonPtY[fiCut]) fHistoTrueMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
          if (fHistoTrueMesonPtAlpha[fiCut]) fHistoTrueMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(),fWeightJetJetMC);
          if (fHistoTrueMesonPtOpenAngle[fiCut]) fHistoTrueMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(),fWeightJetJetMC);
          if (fHistoTrueMotherMesonConvPhotonEtaPhi[fiCut]) fHistoTrueMotherMesonConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta(),fWeightJetJetMC);
        }
      }

      Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      if(isPrimary){
        // Only primary pi0 for efficiency calculation
        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
          if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
              //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
          }
        }
        if (isTrueMeson){
          if (fHistoTruePrimaryMesonInvMassPt[fiCut]) fHistoTruePrimaryMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
          if (fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]) fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
          if (fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]) fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMesons,gamma0MotherLabel)){
            if (fHistoDoubleCountTrueMesonInvMassPt[fiCut]) fHistoDoubleCountTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
            FillMultipleCountMap(fMapMultipleCountTrueMesons,gamma0MotherLabel);
          }
        }
        if (fDoMesonQA > 0 && fIsMC < 2){
          if(isTrueMeson){ // Only primary pi0 for resolution
            if (fHistoTruePrimaryMesonMCPtResolPt[fiCut]) fHistoTruePrimaryMesonMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                                            (mesonCand->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted*fWeightJetJetMC);

          }
        }
      }
    }
  }else if(!isTrueMeson ) { // Background
    if (fDoMesonQA > 1){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Meson or Eta
        if (fHistoTrueBckGGInvMassPt[fiCut]) fHistoTrueBckGGInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());

        if( (((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fMesonPDG
            && ((TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv())))
        ){
          if (fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]) fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        }else if( TrueGammaCandidate1->E()/mesonCand->E() > 0.7 ){
          if (fHistoTrueBckAsymEClustersInvMassPt[fiCut]) fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        }
      }else { // No photon or without mother
        if (fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
    }
  }
  if (isTrueMeson && !matched){
    fVectorRecTrueMesons.push_back(gamma0MotherLabel);
  }
}


//______________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueMesonCandidatesCalo(  AliAODConversionMother *mesonCand,
                                                                            AliAODConversionPhoton *TrueGammaCandidate0,
                                                                            AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Bool_t isTrueMeson            = kFALSE;
  Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma0MotherLabel = -1;

  TParticle * gammaMC0 = 0x0;
  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother(0);
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion() && (gammaMC0->GetMother(0) > -1)){
          gamma0MotherLabel=fMCEvent->Particle(gammaMC0->GetMother(0))->GetMother(0);
        } else {
          gamma0MotherLabel=gammaMC0->GetMother(0);
        }
      }
    }
  }
  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

  Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma1MotherLabel = -1;
  // check if

  TParticle * gammaMC1 = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma1MotherLabel      = gammaMC1->GetMother(0);
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          if(gammaMC1->GetMother(0) > -1) gamma1MotherLabel    = fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
        } else {
          gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }
    }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fMesonPDG){
      isTrueMeson=kTRUE;
    }
  }

  if(isTrueMeson ){// True Meson
    if (fHistoTrueMesonInvMassPt[fiCut]) fHistoTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),  fWeightJetJetMC);
    if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2){
      // both gammas are real gammas
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentPhoton()) {
        if (fHistoTrueMesonCaloPhotonInvMassPt[fiCut]) fHistoTrueMesonCaloPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // both particles are electrons
      if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
        if (fHistoTrueMesonCaloElectronInvMassPt[fiCut]) fHistoTrueMesonCaloElectronInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // both particles are converted electrons
      if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        fHistoTrueMesonCaloConvertedPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // 1 gamma is converted the other one is normals
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
        (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
      ) {
        fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }

      // at least one of the photon is merged
      if (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv() || TrueGammaCandidate0->IsDalitzMerged() || TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() ){
        if (fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // at least one of the photon is merged and part conv
      if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {
        if (fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
    }

    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if ( mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1] ){
        if (fHistoTrueMesonPtY[fiCut]) fHistoTrueMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
        if (fHistoTrueMesonPtAlpha[fiCut]) fHistoTrueMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(), fWeightJetJetMC);
        if (fHistoTrueMesonPtOpenAngle[fiCut]) fHistoTrueMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(), fWeightJetJetMC);
      }

    }
    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(isPrimary){ // Only primary pi0 for efficiency calculation
      // filling primary histograms
      Float_t weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
        if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
          weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
          //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
        }
      }
      if (isTrueMeson){
        if (fHistoTruePrimaryMesonInvMassPt[fiCut]) fHistoTruePrimaryMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted* fWeightJetJetMC);
        if (fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]) fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        if (fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]) fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted* fWeightJetJetMC);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMesons,gamma0MotherLabel) && fHistoDoubleCountTrueMesonInvMassPt[fiCut]) fHistoDoubleCountTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), weighted*fWeightJetJetMC);
      }
      if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC<2){
        if(isTrueMeson){ // Only primary pi0 for resolution
          if (fHistoTruePrimaryMesonMCPtResolPt[fiCut]) fHistoTruePrimaryMesonMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(mesonCand->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
        }
      }
    }
  } else if(!isTrueMeson ){ // Background
    if (fDoMesonQA > 1 && fDoMesonQA < 3){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
        if (fHistoTrueBckGGInvMassPt[fiCut]) fHistoTrueBckGGInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);

        if( (((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == fMesonPDG
          && (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv()))
          ||
          (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fMesonPDG
          && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
        ){
          if (fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]) fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        }else if( (TrueGammaCandidate0->E()/mesonCand->E() > 0.7) || (TrueGammaCandidate1->E()/mesonCand->E() > 0.7) ){
          if (fHistoTrueBckAsymEClustersInvMassPt[fiCut]) fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        }
      } else { // No photon or without mother
        if (fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueMesonCandidatesCaloAOD( AliAODConversionMother *mesonCand,
                                                                              AliAODConversionPhoton *TrueGammaCandidate0,
                                                                              AliAODConversionPhoton *TrueGammaCandidate1)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;

  Bool_t isTrueMeson              = kFALSE;
  Int_t gamma0MCLabel             = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma0MotherLabel         = -1;

  // check if
  AliAODMCParticle * gammaMC0 = 0x0;
  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 0
    gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother();
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion()){
          AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC0->GetMother()));
          gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
        } else gamma0MotherLabel=gammaMC0->GetMother();
      }
    }
  }

  Int_t gamma1MCLabel         = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);   // get most probable MC label
  Int_t gamma1MotherLabel     = -1;

  // check if
  AliAODMCParticle *gammaMC1  = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother();
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
          gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
        } else gamma1MotherLabel=gammaMC1->GetMother();
      }
    }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fMesonPDG){
      isTrueMeson=kTRUE;
    }
  }

  if(isTrueMeson ){// True Meson
    if (fHistoTrueMesonInvMassPt[fiCut]) fHistoTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
    if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC < 2){
      // both gammas are real gammas
      if (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentPhoton()) {
        if (fHistoTrueMesonCaloPhotonInvMassPt[fiCut]) fHistoTrueMesonCaloPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // both particles are electrons
      if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
        if (fHistoTrueMesonCaloElectronInvMassPt[fiCut]) fHistoTrueMesonCaloElectronInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // both particles are converted electrons
      if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
        fHistoTrueMesonCaloConvertedPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // 1 gamma is converted the other one is normals
      if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
        (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
      ) {
        fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }

      // at least one of the photon is merged
      if (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv() || TrueGammaCandidate0->IsDalitzMerged() || TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() ){
        if (fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
      // at least one of the photon is merged and part conv
      if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {
        if (fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]) fHistoTrueMesonCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
    }

    if (fDoMesonQA > 0 && fDoMesonQA < 3){
      if ( mesonCand->M() > fMesonInvMassWindow[0]  && mesonCand->M() < fMesonInvMassWindow[1] ){
        if (fHistoTrueMesonPtY[fiCut]) fHistoTrueMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
        if (fHistoTrueMesonPtAlpha[fiCut]) fHistoTrueMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(), fWeightJetJetMC);
        if (fHistoTrueMesonPtOpenAngle[fiCut]) fHistoTrueMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(), fWeightJetJetMC);
      }
    }

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(isPrimary){  // Only primary pi0 for efficiency calculation
      Float_t weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
        if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
          weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
        }
      }

      if (fHistoTruePrimaryMesonInvMassPt[fiCut]) fHistoTruePrimaryMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted* fWeightJetJetMC);
      if (fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]) fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
      if (fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]) fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted* fWeightJetJetMC);
      if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMesons,gamma0MotherLabel) && fHistoDoubleCountTrueMesonInvMassPt[fiCut]) fHistoDoubleCountTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), weighted*fWeightJetJetMC);
      if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC<2){
        if (fHistoTruePrimaryMesonMCPtResolPt[fiCut]) fHistoTruePrimaryMesonMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                                      (mesonCand->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
      }
    }
  } else if(!isTrueMeson ) { // Background
    if (fDoMesonQA > 1 && fDoMesonQA < 3){
      if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
        if (fHistoTrueBckGGInvMassPt[fiCut]) fHistoTrueBckGGInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);

        if( (((AliAODMCParticle*)AODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == fMesonPDG
          && (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv()))
          ||
          (((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fMesonPDG
          && (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv()))
        ){
          if (fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]) fHistoTrueBckFullMesonContainedInOneClusterInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        }else if( (TrueGammaCandidate0->E()/mesonCand->E() > 0.7) || (TrueGammaCandidate1->E()/mesonCand->E() > 0.7) ){
          if (fHistoTrueBckAsymEClustersInvMassPt[fiCut]) fHistoTrueBckAsymEClustersInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
        }
      } else { // No photon or without mother
        if (fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
      }
    }
  }
}


//______________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueMesonCandidatesConv(  AliAODConversionMother *mesonCand,
                                                                            AliAODConversionPhoton *TrueGammaCandidate0,
                                                                            AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTrueMeson = kFALSE;
    Bool_t isTrueMesonDalitz = kFALSE;
    Bool_t gamma0DalitzCand = kFALSE;
    Bool_t gamma1DalitzCand = kFALSE;
    Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
    Int_t gamma0MotherLabel = -1;
    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(fMCEvent);
      TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(fMCEvent);
      TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetFirstMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==fMesonPDG){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-fMesonPDG;
        }
      }
    }
    if(TrueGammaCandidate1->GetV0Index()<fInputEvent->GetNumberOfV0s()){
      Int_t gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(fMCEvent);
      Int_t gamma1MotherLabel = -1;
      if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
        // Daughters Gamma 1
        TParticle * negativeMC = (TParticle*)TrueGammaCandidate1->GetNegativeMCDaughter(fMCEvent);
        TParticle * positiveMC = (TParticle*)TrueGammaCandidate1->GetPositiveMCDaughter(fMCEvent);
        TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
        if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
          if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
            if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
              gamma1MotherLabel=gammaMC1->GetFirstMother();
            }
          }
          if(gammaMC1->GetPdgCode() ==fMesonPDG ){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-fMesonPDG;
          }
        }
      }
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
        if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == fMesonPDG){
          isTrueMeson=kTRUE;
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMesons,gamma0MotherLabel)){
            if (fHistoDoubleCountTrueMesonInvMassPt[fiCut]) fHistoDoubleCountTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
            FillMultipleCountMap(fMapMultipleCountTrueMesons,gamma0MotherLabel);
          }
        }
      }

      //Identify Dalitz candidate
      if (gamma1DalitzCand || gamma0DalitzCand){
        if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
          if (gamma0MotherLabel == -fMesonPDG) isTrueMesonDalitz = kTRUE;
        }
        if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
          if (gamma1MotherLabel == -fMesonPDG) isTrueMesonDalitz = kTRUE;
        }
      }


      if(isTrueMeson ){// True Meosn
        if (fHistoTrueMesonInvMassPt[fiCut])fHistoTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
        if (fDoMesonQA > 0){
          if ( mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1]){
            if(fIsMC < 2){
              if (fHistoTrueMesonPtY[fiCut]) fHistoTrueMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
              if (fHistoTrueMesonPtOpenAngle[fiCut]) fHistoTrueMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle());
              if (fHistoTrueMesonPtAlpha[fiCut]) fHistoTrueMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(),fWeightJetJetMC);
            }
          }
        }
        Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if(isPrimary){ // Only primary pi0 for efficiency calculation
          Float_t weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
            if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
              //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
            }
          }
          if (fHistoTruePrimaryMesonInvMassPt[fiCut]) fHistoTruePrimaryMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
          if (fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]) fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
          if (fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]) fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);

          if (fDoMesonQA > 0 && fIsMC < 2){
            if (fHistoTruePrimaryMesonMCPtResolPt[fiCut]) fHistoTruePrimaryMesonMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(mesonCand->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted);
          }
        }
      } else if(!isTrueMeson ){ // Background
        if (fDoMesonQA > 1 && fIsMC < 2){
          if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
            if (fHistoTrueBckGGInvMassPt[fiCut])fHistoTrueBckGGInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
          } else { // No photon or without mother
            if (fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
          }
        }
        if (!isTrueMesonDalitz && (gamma0DalitzCand || gamma1DalitzCand)){
          if (fDoMesonQA > 1 && fIsMC < 2 && fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::ProcessTrueMesonCandidatesConvAOD( AliAODConversionMother *mesonCand,
                                                                              AliAODConversionPhoton *TrueGammaCandidate0,
                                                                              AliAODConversionPhoton *TrueGammaCandidate1)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  Bool_t isTrueMeson = kFALSE;
  Bool_t isTrueMesonDalitz = kFALSE;
  Bool_t gamma0DalitzCand = kFALSE;
  Bool_t gamma1DalitzCand = kFALSE;

  if (AODMCTrackArray!=NULL && TrueGammaCandidate0 != NULL){
    AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
    AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

    Int_t gamma0MCLabel = -1;
    Int_t gamma0MotherLabel = -1;
    if(!positiveMC||!negativeMC)
      return;

    if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gamma0MCLabel = positiveMC->GetMother();
    }

    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==fMesonPDG){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-fMesonPDG;
        }
      }
    }
    positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
    negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));

    Int_t gamma1MCLabel = -1;
    Int_t gamma1MotherLabel = -1;
    if(!positiveMC||!negativeMC)
      return;

    if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gamma1MCLabel = positiveMC->GetMother();
    }
    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...
          if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma1MotherLabel=gammaMC1->GetMother();
          }
        }
        if(gammaMC1->GetPdgCode() ==fMesonPDG ){ // Dalitz candidate
          gamma1DalitzCand = kTRUE;
          gamma1MotherLabel=-fMesonPDG;
        }
      }
    }
    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == fMesonPDG){
        isTrueMeson=kTRUE;
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMesons,gamma0MotherLabel)){
          if (fHistoDoubleCountTrueMesonInvMassPt[fiCut]) fHistoDoubleCountTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
          FillMultipleCountMap(fMapMultipleCountTrueMesons,gamma0MotherLabel);
        }
      }
    }

    //Identify Dalitz candidate
    if (gamma1DalitzCand || gamma0DalitzCand){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
        if (gamma0MotherLabel == -fMesonPDG) isTrueMesonDalitz = kTRUE;
      }
      if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
        if (gamma1MotherLabel == -fMesonPDG) isTrueMesonDalitz = kTRUE;
      }
    }

    if(isTrueMeson ){// True Meson
      if (fHistoTrueMesonInvMassPt[fiCut])fHistoTrueMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
      if (fDoMesonQA > 0){
        if ( mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1]){
          if(fIsMC < 2){
            if (fHistoTrueMesonPtY[fiCut]) fHistoTrueMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
            if (fHistoTrueMesonPtOpenAngle[fiCut]) fHistoTrueMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle());
          }
          if (fHistoTrueMesonPtAlpha[fiCut]) fHistoTrueMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(),fWeightJetJetMC);
        }

      }
      Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

      if(isPrimary){ // Only primary pi0 for efficiency calculation
        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
          if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
            //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
          }
        }
        if (fHistoTruePrimaryMesonInvMassPt[fiCut]) fHistoTruePrimaryMesonInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);
        if (fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]) fHistoTruePrimaryMesonW0WeightingInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
        if (fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]) fProfileTruePrimaryMesonWeightsInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),weighted*fWeightJetJetMC);

        if (fDoMesonQA > 0 && fIsMC < 2){
          if (fHistoTruePrimaryMesonMCPtResolPt[fiCut])fHistoTruePrimaryMesonMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                                        (mesonCand->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted);
        }
      }
    } else if(!isTrueMeson ) { // Background
      if (fDoMesonQA > 1 && fIsMC < 2){
        if(!(gamma0MotherLabel>-1 && gamma1MotherLabel>-1)){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
          if (fHistoTrueBckGGInvMassPt[fiCut])fHistoTrueBckGGInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        } else { // No photon or without mother
          if (fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
        }
      }
      if (!isTrueMesonDalitz && (gamma0DalitzCand || gamma1DalitzCand)){
        if (fDoMesonQA > 1 && fIsMC < 2)if (fHistoTrueBckContInvMassPt[fiCut]) fHistoTrueBckContInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt());
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::CalculateBackground(){

  // set current BG handler
  AliGammaConversionAODBGHandler* currBGHandler = NULL;
  TList* currPhotonList                         = NULL;
  if (fMesonRecoMode == 0){
    currBGHandler         = fBGHandler[fiCut];
    currPhotonList        = fGammaCandidates;
  } else if (fMesonRecoMode == 1){
    currBGHandler         = fBGClusHandler[fiCut];
    currPhotonList        = fGammaCandidates;
  } else if (fMesonRecoMode == 2){
    currBGHandler         = fBGClusHandler[fiCut];
    currPhotonList        = fClusterCandidates;
  }

  Int_t zbin = currBGHandler->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    mbin = currBGHandler->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  } else {
    mbin = currBGHandler->GetMultiplicityBinIndex(currPhotonList->GetEntries());
  }

  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    for(Int_t nEventsInBG=0;nEventsInBG<currBGHandler->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = currBGHandler->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if (fMesonRecoMode < 2){
        if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0 ){
          bgEventVertex = currBGHandler->GetBGEventVertex(zbin,mbin,nEventsInBG);
        }
      }
      for(Int_t iCurrent=0;iCurrent<currPhotonList->GetEntries();iCurrent++){
        AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(currPhotonList->At(iCurrent));
        if (!currentEventGoodV0.GetUseForMesonPair() )continue;
        for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
          AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
          if(fMoveParticleAccordingToVertex == kTRUE){
            if (bgEventVertex){
              MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
            }
          }
          if (fMesonRecoMode < 2){
            if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
              if (bgEventVertex){
                RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
              }
            }
          }
          AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
          backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate, kFALSE,
            ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()), currentEventGoodV0.GetLeadingCellID(), previousGoodV0.GetLeadingCellID(), currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton())){
            if (fHistoMotherBackInvMassPt[fiCut]) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
            if(fDoTHnSparse){
              Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
              if (fSparseMotherBackInvMassPtZM[fiCut]) fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
            }
          }
        }
      }
    }
  }else {
    for(Int_t nEventsInBG=0;nEventsInBG <currBGHandler->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = currBGHandler->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        if (fMesonRecoMode < 2){
          if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0 ){
            bgEventVertex = currBGHandler->GetBGEventVertex(zbin,mbin,nEventsInBG);
          }
        }
        for(Int_t iCurrent=0;iCurrent<currPhotonList->GetEntries();iCurrent++){
          AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(currPhotonList->At(iCurrent));
          if (!currentEventGoodV0.GetUseForMesonPair() )continue;
          for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
            AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
            if(fMoveParticleAccordingToVertex == kTRUE){
              if (bgEventVertex){
                MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
              }
            }
            if (fMesonRecoMode < 2){
              if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
                if (bgEventVertex){
                  RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
                }
              }
            }
            AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
            backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,
                ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()), currentEventGoodV0.GetLeadingCellID(), previousGoodV0.GetLeadingCellID(), currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton())){
              if (fHistoMotherBackInvMassPt[fiCut]) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt(),fWeightJetJetMC);
              if(fDoTHnSparse){
                Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                if (fSparseMotherBackInvMassPtZM[fiCut]) fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
              }
            }
            delete backgroundCandidate;
            backgroundCandidate = 0x0;
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::CalculateBackgroundSwapp(){

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {

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

    std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
    std::vector<std::array<Double_t, 2>> vSwappingInvMassPTAlphaCut;
    vSwappingInvMassPT.clear();
    vSwappingInvMassPTAlphaCut.clear();
    vSwappingInvMassPT.resize(0);
    vSwappingInvMassPTAlphaCut.resize(0);
    Double_t tempMultWeightSwapping = 1; // weight taking multiplicity of event into account

    if (fMesonRecoMode == 0){ // PCM - PCM
      // curcial requierment is that the event has at least 3 cluster candidates
      if(fGammaCandidates->GetEntries() > 2 ){
        for(Int_t iCurrent1=0;iCurrent1<fGammaCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton* currentEventGoodV0Temp1 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent1));
          if(!currentEventGoodV0Temp1) continue;
          for(Int_t iCurrent2=iCurrent1+1; iCurrent2 < fGammaCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton* currentEventGoodV0Temp2 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent2));
            if(!currentEventGoodV0Temp2) continue;
            for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){

              lvRotationPhoton1.SetX(currentEventGoodV0Temp1->Px());
              lvRotationPhoton1.SetY(currentEventGoodV0Temp1->Py());
              lvRotationPhoton1.SetZ(currentEventGoodV0Temp1->Pz());
              lvRotationPhoton1.SetE(currentEventGoodV0Temp1->E());

              lvRotationPhoton2.SetX(currentEventGoodV0Temp2->Px());
              lvRotationPhoton2.SetY(currentEventGoodV0Temp2->Py());
              lvRotationPhoton2.SetZ(currentEventGoodV0Temp2->Pz());
              lvRotationPhoton2.SetE(currentEventGoodV0Temp2->E());
              lvRotationPion = (lvRotationPhoton1 + lvRotationPhoton2).Vect();

              // rotate both photons around the momentum vector of their hypothetical mother particle
              if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0 || ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1)){
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0) rotationAngle = TMath::Pi()/2.0; // rotate by 90 degree
                else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
                   Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
                   rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
                }
                lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
                lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);
              }
              std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
              std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));
              for(auto const& kCurrentGammaCandidates  : *fGammaCandidates){
                if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentGammaCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentGammaCandidates)) continue;

                std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));
                std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));
                if( fabs(currentEventGoodV0Temp1->Eta()) <= ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetEtaCut())
                {
                  if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                  }
                }
                if( fabs(currentEventGoodV0Temp2->Eta()) <= ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetEtaCut())
                {
                  if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                  }
                }
              }
            }
          }
        }
        // Fill the histograms
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPT.size() > 0){
          tempMultWeightSwapping = (0.5*(fGammaCandidates->GetEntries()*fGammaCandidates->GetEntries() - fGammaCandidates->GetEntries()))/(vSwappingInvMassPT.size());
        }
        for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
          fHistoMotherBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
        }
      }

    } else if (fMesonRecoMode == 1){ // PCM - CALO
      for(Int_t iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
        AliAODConversionPhoton* currentEventGoodV0Temp1 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent1));

        for(Int_t iCurrent2=0;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
          AliAODConversionPhoton *currentEventGoodV0Temp2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(iCurrent2));
          if (currentEventGoodV0Temp2==NULL) continue;

          for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){

            lvRotationPhoton1.SetX(currentEventGoodV0Temp1->Px());
            lvRotationPhoton1.SetY(currentEventGoodV0Temp1->Py());
            lvRotationPhoton1.SetZ(currentEventGoodV0Temp1->Pz());
            lvRotationPhoton1.SetE(currentEventGoodV0Temp1->E());

            lvRotationPhoton2.SetX(currentEventGoodV0Temp2->Px());
            lvRotationPhoton2.SetY(currentEventGoodV0Temp2->Py());
            lvRotationPhoton2.SetZ(currentEventGoodV0Temp2->Pz());
            lvRotationPhoton2.SetE(currentEventGoodV0Temp2->E());

            lvRotationPion = (lvRotationPhoton1 + lvRotationPhoton2).Vect();

            // rotate both photons around the momentum vector of their hypothetical mother particle
            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0 || ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1)){
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0){ // rotate by 90 degree
                rotationAngle = TMath::Pi()/2.0;
              } else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
                 Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
                 rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
              }
              lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
              lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);
            }

            // Fill Eta Phi Map for Calo Photon
            cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), lvRotationPhoton1.Phi());
            if(!fDoLightOutput){
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent))){
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), 1);
              }
            }

            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));

            for(auto const& kCurrentClusterCandidates  : *fClusterCandidates){
              if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentClusterCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentClusterCandidates)) continue;

              std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));

              if( fabs(currentEventGoodV0Temp2->Eta()) <= ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetEtaCut() )
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                }
              }
            }
            for(auto const& kCurrentGammaCandidates  : *fGammaCandidates){
              if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentGammaCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentGammaCandidates)) continue;

              std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));
              // std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));

              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent)))
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                }
              }
            }
          }
        }
      }
      // Fill the histograms
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPT.size() > 0){
        tempMultWeightSwapping = (fGammaCandidates->GetEntries()*fClusterCandidates->GetEntries())/(vSwappingInvMassPT.size());
      }
      for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
        fHistoMotherBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
      }

    } else if (fMesonRecoMode == 2){ // CALO - CALO
      if(fClusterCandidates->GetEntries() > 2 ){
        for(Int_t iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
          AliAODConversionPhoton* currentEventGoodV0Temp1 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent1));

          for(Int_t iCurrent2=iCurrent1+1;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
            AliAODConversionPhoton* currentEventGoodV0Temp2 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent2));

            for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){
              lvRotationPhoton1.SetX(currentEventGoodV0Temp1->Px());
              lvRotationPhoton1.SetY(currentEventGoodV0Temp1->Py());
              lvRotationPhoton1.SetZ(currentEventGoodV0Temp1->Pz());
              lvRotationPhoton1.SetE(currentEventGoodV0Temp1->E());

              lvRotationPhoton2.SetX(currentEventGoodV0Temp2->Px());
              lvRotationPhoton2.SetY(currentEventGoodV0Temp2->Py());
              lvRotationPhoton2.SetZ(currentEventGoodV0Temp2->Pz());
              lvRotationPhoton2.SetE(currentEventGoodV0Temp2->E());

              lvRotationPion = (lvRotationPhoton1 + lvRotationPhoton2).Vect();

              // rotate both photons around the momentum vector of their hypothetical mother particle
              if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0 || ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1)){
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0) rotationAngle = TMath::Pi()/2.0; // rotate by 90 degree
                else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
                   Double_t temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
                   rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
                }
                lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
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


              cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));
              cellIDRotatedPhoton2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()));

              if(!fDoLightOutput){
                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent))){
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), 1);
                }
                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, lvRotationPhoton2.Phi(), fInputEvent))){
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()), 1);
                }
              }

              std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
              std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));

              for(auto const& kCurrentClusterCandidates  : *fClusterCandidates){
                if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentClusterCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentClusterCandidates)){ continue;}

                std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));
                std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));

                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, lvRotationPhoton1.Phi(), fInputEvent)) && lvRotationPhoton1.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
                {
                  if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton1, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID()))
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                    // if((!fDoLightOutput || fDoPi0Only) && TMath::Abs(backgroundCandidate1->GetAlpha())<0.1){
                    //   vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                    // }
                  }
                }
                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, lvRotationPhoton2.Phi(), fInputEvent)) && lvRotationPhoton2.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
                {
                  if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton2, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID()))
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                    // if((!fDoLightOutput || fDoPi0Only) && TMath::Abs(backgroundCandidate2->GetAlpha())<0.1){
                    //   vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                    // }
                  }
                }
              }
            }
          }
        }
        // Fill the histograms
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPT.size() > 0){
          tempMultWeightSwapping = (0.5*(fClusterCandidates->GetEntries()*fClusterCandidates->GetEntries() - fClusterCandidates->GetEntries()))/(vSwappingInvMassPT.size());
        }
        for(Int_t i = 0; i < (Int_t)vSwappingInvMassPT.size(); i++){
          fHistoMotherBackInvMassPt[fiCut]->Fill(vSwappingInvMassPT.at(i)[0], vSwappingInvMassPT.at(i)[1], tempMultWeightSwapping*fWeightJetJetMC);
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::CalculateBackgroundRP(){

  // set current BG handler
  AliConversionAODBGHandlerRP* currBGHandler    = NULL;
  TList* currPhotonList                         = NULL;
  if (fMesonRecoMode == 0){
    currBGHandler         = fBGHandlerRP[fiCut];
    currPhotonList        = fGammaCandidates;
  } else if (fMesonRecoMode == 1){
    currBGHandler         = fBGClusHandlerRP[fiCut];
    currPhotonList        = fGammaCandidates;
  } else if (fMesonRecoMode == 2){
    currBGHandler         = fBGClusHandlerRP[fiCut];
    currPhotonList        = fClusterCandidates;
  }

  Int_t zbin = currBGHandler->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    mbin = currBGHandler->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  } else {
    mbin = currBGHandler->GetMultiplicityBinIndex(currPhotonList->GetEntries());
  }

  //Rotation Method
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseRotationMethod()){
    // Correct for the number of rotations
    // BG is for rotation the same, except for factor NRotations
    Double_t weight=1./Double_t(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents());
    for(Int_t firstGammaIndex=0;firstGammaIndex<currPhotonList->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(currPhotonList->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      if (!gamma0->GetUseForMesonPair() ) continue;
      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<currPhotonList->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(currPhotonList->At(secondGammaIndex));
        if (gamma1 == NULL) continue;
        if (!gamma1->GetUseForMesonPair() ) continue;
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(gamma1,fInputEvent))continue;
        for(Int_t nRandom=0;nRandom<((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents();nRandom++){
          RotateParticle(gamma1);
          AliAODConversionMother backgroundCandidate(gamma0,gamma1);
          backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(&backgroundCandidate, kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), gamma0->GetLeadingCellID(), gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton() )){

            if (fHistoMotherBackInvMassPt[fiCut]) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
            if(fDoTHnSparse){
              Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
              if (fSparseMotherBackInvMassPtZM[fiCut]) fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
            }
          }
        }
      }
    }

  }else {
    // Do Event Mixing
    for(Int_t nEventsInBG=0;nEventsInBG <currBGHandler->GetNBGEvents(currPhotonList,fInputEvent);nEventsInBG++){
      AliGammaConversionPhotonVector *previousEventGammas = currBGHandler->GetBGGoodGammas(currPhotonList,fInputEvent,nEventsInBG);
      if(previousEventGammas){
        // test weighted background
        Double_t weight=1.0;
        // Correct for the number of eventmixing:
        // N gammas -> (N-1) + (N-2) +(N-3) ...+ (N-(N-1))  using sum formula sum(i)=N*(N-1)/2  -> N*(N-1)/2
        // real combinations (since you cannot combine a photon with its own)
        // but BG leads to N_{a}*N_{b}combinations
        weight*=0.5*(Double_t(currPhotonList->GetEntries()-1))/Double_t(previousEventGammas->size());

        for(Int_t iCurrent=0;iCurrent<currPhotonList->GetEntries();iCurrent++){
          AliAODConversionPhoton *gamma0 = (AliAODConversionPhoton*)(currPhotonList->At(iCurrent));
          if (!gamma0->GetUseForMesonPair() ) continue;
          for(UInt_t iPrevious=0;iPrevious<previousEventGammas->size();iPrevious++){
            AliAODConversionPhoton *gamma1 = (AliAODConversionPhoton*)(previousEventGammas->at(iPrevious));
            AliAODConversionMother backgroundCandidate(gamma0,gamma1);
            backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
            if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                ->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), gamma0->GetLeadingCellID(), gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton())){
              if (fHistoMotherBackInvMassPt[fiCut]) fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt(),fWeightJetJetMC);
              if(fDoTHnSparse){
                  Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
                  if (fSparseMotherBackInvMassPtZM[fiCut]) fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
              }
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::RotateParticle(AliAODConversionPhoton *gamma){
  Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
  Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
  Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){
  previousEventEP=previousEventEP+TMath::Pi();
  thisEventEP=thisEventEP+TMath::Pi();
  Double_t rotationValue= thisEventEP-previousEventEP;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
  //see header file for documentation

  Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
  Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
  Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

  Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
  particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::UpdateEventByEventData(){
  //see header file for documentation
  if(fGammaCandidates->GetEntries() >1 && fMesonRecoMode == 0 ){
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                  fV0Reader->GetNumberOfPrimaryTracks(), fEventPlaneAngle);
    }else { // means we use #V0s for multiplicity
      fBGHandler[fiCut]->AddEvent(fGammaCandidates, fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                  fGammaCandidates->GetEntries(), fEventPlaneAngle);
    }
  } else if((fGammaCandidates->GetEntries() >0 || fClusterCandidates->GetEntries() > 0 ) && fMesonRecoMode == 1 ){
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                  fV0Reader->GetNumberOfPrimaryTracks(), fEventPlaneAngle);
      fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                      fV0Reader->GetNumberOfPrimaryTracks(), fEventPlaneAngle);
    }else { // means we use #V0s for multiplicity
      fBGHandler[fiCut]->AddEvent(fGammaCandidates, fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                  fGammaCandidates->GetEntries(), fEventPlaneAngle);
      fBGClusHandler[fiCut]->AddEvent(fClusterCandidates, fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                      fGammaCandidates->GetEntries(), fEventPlaneAngle);
    }
  } else if(fClusterCandidates->GetEntries() >1 && fMesonRecoMode == 2 ){
    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                      fV0Reader->GetNumberOfPrimaryTracks(), fEventPlaneAngle);
    }else { // means we use #V0s for multiplicity
      fBGClusHandler[fiCut]->AddEvent(fClusterCandidates, fInputEvent->GetPrimaryVertex()->GetX(), fInputEvent->GetPrimaryVertex()->GetY(), fInputEvent->GetPrimaryVertex()->GetZ(),
                                      fClusterCandidates->GetEntries(), fEventPlaneAngle);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::RelabelAODPhotonCandidates(Bool_t mode){

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
void AliAnalysisTaskHeavyNeutralMesonToGG::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskHeavyNeutralMesonToGG::Terminate(const Option_t *)
{
  //fOutputContainer->Print(); // Will crash on GRID
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskHeavyNeutralMesonToGG::CheckVectorOnly(vector<Int_t> &vec, Int_t tobechecked){
  if(tobechecked > -1){
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else return false;
  }
  return false;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskHeavyNeutralMesonToGG::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked){
  if(tobechecked > -1){
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
void AliAnalysisTaskHeavyNeutralMesonToGG::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
  if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
  else ma[tobechecked] = 2;
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskHeavyNeutralMesonToGG::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
  map<Int_t, Int_t>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second, fWeightJetJetMC);
  }
  ma.clear();
  return;
}
