/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock, Joshua König                                  *
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
#include "AliAnalysisTaskConvCaloCalibration.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliGAKFVertex.h"
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

ClassImp(AliAnalysisTaskConvCaloCalibration)

//________________________________________________________________________
AliAnalysisTaskConvCaloCalibration::AliAnalysisTaskConvCaloCalibration(): AliAnalysisTaskSE(),
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
  fGeomEMCAL(NULL),
  fElecSelector(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fOutputContainer(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fClusterCutArray(NULL),
  fMesonCutArray(NULL),
  fReaderGammas(NULL),
  fSelectorElectronIndex(),
  fSelectorPositronIndex(),
  fV0Electrons(),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fFileNameBroken(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fTreeBrokenFiles(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fHistoMotherMesonPtY(NULL),
  fHistoMotherMesonPtAlpha(NULL),
  fHistoMotherMesonPtOpenAngle(NULL),
  fHistoMotherMesonConvPhotonEtaPhi(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoMotherInvMassECalibSM(NULL),
  fHistoMotherBackInvMassECalibSM(NULL),
  fHistoMotherInvMassECalib(NULL),
  fHistoMotherBackInvMassECalib(NULL),
  fHistoEVsNCellsInPiMass(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  // fHistoMCHeaders(NULL),
  fHistoConvGammaPt(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusGammaPtSM(NULL),
  fHistoClusGammaESM(NULL),
  fHistoClusGammaERx(NULL),
  fHistoClusGammaERxSM(NULL),
  fHistoClusGammaERxNCellCrit(NULL),
  fHistoClusGammaERxNCellCritSM(NULL),
  fHistoClusTrackdEtaSM(NULL),
  fHistoClusTrackdPhiSM(NULL),
  fHistoClusHighPtTrackdEtaSM(NULL),
  fHistoClusHighPtTrackdPhiSM(NULL),
  fHistoEVsM02(NULL),
  fHistoEVsM02NCell4(NULL),
  fHistoMotherInvMassRejected(NULL),
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
  fnModules(20),
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
  fTrackMatcherRunningMode(0),
  fUseEletronMatchingCalibration(kFALSE),
  fGenPhaseSpace()
{

}

//________________________________________________________________________
AliAnalysisTaskConvCaloCalibration::AliAnalysisTaskConvCaloCalibration(const char *name):
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
  fGeomEMCAL(NULL),
  fElecSelector(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fOutputContainer(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fCutArray(NULL),
  fClusterCutArray(NULL),
  fMesonCutArray(NULL),
  fReaderGammas(NULL),
  fSelectorElectronIndex(),
  fSelectorPositronIndex(),
  fV0Electrons(),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fFileNameBroken(NULL),
  fSparseMotherInvMassPtZM(NULL),
  fSparseMotherBackInvMassPtZM(NULL),
  fTreeBrokenFiles(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherMatchedInvMassPt(NULL),
  fHistoMotherBackInvMassPt(NULL),
  fHistoMotherMesonPtY(NULL),
  fHistoMotherMesonPtAlpha(NULL),
  fHistoMotherMesonPtOpenAngle(NULL),
  fHistoMotherMesonConvPhotonEtaPhi(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoMotherInvMassECalibSM(NULL),
  fHistoMotherBackInvMassECalibSM(NULL),
  fHistoMotherInvMassECalib(NULL),
  fHistoMotherBackInvMassECalib(NULL),
  fHistoEVsNCellsInPiMass(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  // fHistoMCHeaders(NULL),
  fHistoConvGammaPt(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusGammaPtSM(NULL),
  fHistoClusGammaESM(NULL),
  fHistoClusGammaERx(NULL),
  fHistoClusGammaERxSM(NULL),
  fHistoClusGammaERxNCellCrit(NULL),
  fHistoClusGammaERxNCellCritSM(NULL),
  fHistoClusTrackdEtaSM(NULL),
  fHistoClusTrackdPhiSM(NULL),
  fHistoClusHighPtTrackdEtaSM(NULL),
  fHistoClusHighPtTrackdPhiSM(NULL),
  fHistoEVsM02(NULL),
  fHistoEVsM02NCell4(NULL),
  fHistoMotherInvMassRejected(NULL),
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
  fnModules(20),
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
  fTrackMatcherRunningMode(0),
  fUseEletronMatchingCalibration(kFALSE),
  fGenPhaseSpace()
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskConvCaloCalibration::~AliAnalysisTaskConvCaloCalibration()
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
void AliAnalysisTaskConvCaloCalibration::InitBack(){

  const Int_t nDim = 4;
  Int_t nBins[nDim] = {800,300,7,4};
  Double_t xMin[nDim] = {0,0, 0,0};
  Double_t xMax[nDim] = {0.8,30,7,4};

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
      }else{
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
void AliAnalysisTaskConvCaloCalibration::UserCreateOutputObjects(){

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

      fHistoNGammaCaloCandidates          = new TH1F*[fnCuts];
      fHistoClusGammaPt                   = new TH1F*[fnCuts];
      fHistoClusGammaE                    = new TH1F*[fnCuts];
      fHistoClusGammaPtSM                 = new TH1F**[fnCuts];
      fHistoClusGammaESM                  = new TH1F**[fnCuts];
      for(Int_t icuts = 0; icuts < fnCuts; icuts++){
        fHistoClusGammaPtSM[icuts]        = new TH1F*[fnModules];
        fHistoClusGammaESM[icuts]         = new TH1F*[fnModules];
      }
      if(fDoMesonQA > 0){
        fHistoClusGammaERx             = new TH1F*[fnCuts];
        fHistoClusGammaERxNCellCrit    = new TH1F*[fnCuts];
        fHistoClusGammaERxSM           = new TH1F**[fnCuts];
        fHistoClusGammaERxNCellCritSM  = new TH1F**[fnCuts];
        fHistoClusGammaERxNCellCritSM  = new TH1F**[fnCuts];
        fHistoClusHighPtTrackdEtaSM    = new TH1F**[fnCuts];
        fHistoClusHighPtTrackdPhiSM    = new TH1F**[fnCuts];
        fHistoClusTrackdEtaSM          = new TH1F**[fnCuts];
        fHistoClusTrackdPhiSM          = new TH1F**[fnCuts];
        fHistoEVsM02                   = new TH2F*[fnCuts];
        fHistoEVsM02NCell4             = new TH2F*[fnCuts];
        for(Int_t icuts = 0; icuts < fnCuts; icuts++){
          fHistoClusGammaERxSM[icuts]           = new TH1F*[fnModules];
          fHistoClusGammaERxNCellCritSM[icuts]  = new TH1F*[fnModules];
          fHistoClusHighPtTrackdEtaSM[icuts]    = new TH1F*[fnModules];
          fHistoClusHighPtTrackdPhiSM[icuts]    = new TH1F*[fnModules];
          fHistoClusTrackdEtaSM[icuts]          = new TH1F*[fnModules];
          fHistoClusTrackdPhiSM[icuts]          = new TH1F*[fnModules];
        }
      }
    }

    fHistoMotherInvMassPt             = new TH2F*[fnCuts];
    fHistoMotherInvMassRejected       = new TH1F*[fnCuts];
    fHistoMotherBackInvMassPt         = new TH2F*[fnCuts];
    fHistoMotherInvMassECalibSM         = new TH2F**[fnCuts];
    fHistoMotherBackInvMassECalibSM     = new TH2F**[fnCuts];
    for(Int_t icuts = 0; icuts < fnCuts; icuts++){
      fHistoMotherInvMassECalibSM[icuts] = new TH2F*[fnModules];
      fHistoMotherBackInvMassECalibSM[icuts] = new TH2F*[fnModules];
    }
    fHistoMotherInvMassECalib           = new TH2F*[fnCuts];
    fHistoMotherBackInvMassECalib       = new TH2F*[fnCuts];
    if(!fDoLightOutput && fMesonRecoMode == 1){
      fHistoMotherMatchedInvMassPt      = new TH2F*[fnCuts];
    }

    if(!fDoLightOutput && fMesonRecoMode > 0){
      fHistoEVsNCellsInPiMass     = new TH2F*[fnCuts];
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
      fMesonInvMassMax        = 0.300;
      fMesonInvMassNBins      = 300;
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
    Float_t maxClusterPt        = 50;
    Double_t *arrPtBinning      = new Double_t[1200];
    Double_t *arrQAPtBinning    = new Double_t[1200];
    Double_t *arrClusPtBinning  = new Double_t[1200];
    if(fMesonRecoMode < 2 && fDoMesonQA == 0){
      nBinsPt                   = 235;
      minPt                     = 0;
      maxPt                     = 50;
      binWidthPt                = 0.05;
      for(Int_t i=0; i<nBinsPt+1;i++){
        if (i < 1) arrPtBinning[i]              = 0.3*i;
        else if(i<55) arrPtBinning[i]           = 0.3+0.05*(i-1);
        else if(i<125) arrPtBinning[i]          = 3.+0.1*(i-55);
        else if(i<185) arrPtBinning[i]          = 10.+0.25*(i-125);
        else if(i<235) arrPtBinning[i]          = 25.+0.5*(i-185);
        else  arrPtBinning[i]                   = maxPt;
      }
      nBinsQAPt                 = 221;
      maxQAPt                   = 50;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<60) arrQAPtBinning[i]              = 0.05*i;
        else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
        else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
        else if(i<210) arrQAPtBinning[i]        = 20.+0.5*(i-170);
        else if(i<221) arrQAPtBinning[i]        = 40.+1.0*(i-210);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsClusterPt            = 351;
      maxClusterPt              = 50;
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        if (i < 1) arrClusPtBinning[i]          = 0.3*i;
        else if(i<55) arrClusPtBinning[i]       = 0.3+0.05*(i-1);
        else if(i<125) arrClusPtBinning[i]      = 3.+0.1*(i-55);
        else if(i<155) arrClusPtBinning[i]      = 10.+0.2*(i-125);
        else if(i<211) arrClusPtBinning[i]      = 16.+0.25*(i-155);
        else if(i<251) arrClusPtBinning[i]      = 30.+0.5*(i-211);
        else arrClusPtBinning[i]                = maxClusterPt;
    }

  } else if(fMesonRecoMode == 2 && fDoMesonQA == 0){
        nBinsPt                   = 165;
        minPt                     = 0;
        maxPt                     = 20;
        binWidthPt                = 0.05;
        for(Int_t i=0; i<nBinsPt+1;i++){
          if (i < 1) arrPtBinning[i]              = 0.3*i;
          else if(i<55) arrPtBinning[i]           = 0.3+0.05*(i-1);
          else if(i<125) arrPtBinning[i]          = 3.+0.1*(i-55);
          else if(i<165) arrPtBinning[i]          = 10.+0.25*(i-125);
          else  arrPtBinning[i]                   = maxPt;
        }
        nBinsQAPt                 = 210;
        maxQAPt                   = 20;
        for(Int_t i=0; i<nBinsQAPt+1;i++){
          if(i<60) arrQAPtBinning[i]              = 0.05*i;
          else if(i<130) arrQAPtBinning[i]        = 3.+0.1*(i-60);
          else if(i<170) arrQAPtBinning[i]        = 10.+0.25*(i-130);
          else arrQAPtBinning[i]                  = maxQAPt;
        }
        nBinsClusterPt            = 171;
        maxClusterPt              = 20;
        for(Int_t i=0; i<nBinsClusterPt+1;i++){
          if (i < 1) arrClusPtBinning[i]          = 0.3*i;
          else if(i<55) arrClusPtBinning[i]       = 0.3+0.05*(i-1);
          else if(i<125) arrClusPtBinning[i]      = 3.+0.1*(i-55);
          else if(i<155) arrClusPtBinning[i]      = 10.+0.2*(i-125);
          else if(i<171) arrClusPtBinning[i]      = 16.+0.25*(i-155);
          else arrClusPtBinning[i]                = maxClusterPt;
        }
              // default binning
    } else {
      nBinsPt                   = 310;
      minPt                     = 0;
      maxPt                     = 100;
      binWidthPt                = 0.1;
      for(Int_t i=0; i<nBinsPt+1;i++){
        if (i < 1) arrPtBinning[i]              = 0.3*i;
        else if(i<198) arrPtBinning[i]          = 0.3+0.1*(i-1);
        else if(i<238) arrPtBinning[i]          = 20.+0.25*(i-198);
        else if(i<278) arrPtBinning[i]          = 30.+0.5*(i-238);
        else if(i<298) arrPtBinning[i]          = 50.+1.0*(i-278);
        else if(i<310) arrPtBinning[i]          = 70.+2.5*(i-298);
        else  arrPtBinning[i]                   = maxPt;
      }
      nBinsQAPt                 = 240;
      maxQAPt                   = 100;
      for(Int_t i=0; i<nBinsQAPt+1;i++){
        if(i<100) arrQAPtBinning[i]             = 0.1*i;
        else if(i<140) arrQAPtBinning[i]        = 10.+0.25*(i-100);
        else if(i<180) arrQAPtBinning[i]        = 20.+0.5*(i-140);
        else if(i<240) arrQAPtBinning[i]        = 40.+1.0*(i-180);
        else arrQAPtBinning[i]                  = maxQAPt;
      }
      nBinsClusterPt            = 310;
      maxClusterPt              = 100;
      for(Int_t i=0; i<nBinsClusterPt+1;i++){
        if (i < 1) arrClusPtBinning[i]          = 0.3*i;
        else if(i<198) arrClusPtBinning[i]      = 0.3+0.1*(i-1);
        else if(i<238) arrClusPtBinning[i]      = 20.+0.25*(i-198);
        else if(i<278) arrClusPtBinning[i]      = 30.+0.5*(i-238);
        else if(i<298) arrClusPtBinning[i]      = 50.+1.0*(i-278);
        else if(i<310) arrClusPtBinning[i]      = 70.+2.5*(i-298);
        else  arrClusPtBinning[i]               = maxClusterPt;
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
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
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
    if (fMesonRecoMode > 0){ //if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
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
        }
      }
    }

    for(Int_t iModules = 0; iModules < fnModules; iModules++ ){
      fHistoMotherInvMassECalibSM[iCut][iModules]         = new TH2F(Form("ESD_Mother_InvMass_E_Calib_SM%i",iModules), Form("ESD_Mother_InvMass_E_Calib_SM%i",iModules), fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
      fHistoMotherInvMassECalibSM[iCut][iModules]->SetXTitle("M_{inv} (GeV/c^{2})");
      fHistoMotherInvMassECalibSM[iCut][iModules]->SetYTitle("E_{cluster}(GeV)");
      fESDList[iCut]->Add(fHistoMotherInvMassECalibSM[iCut][iModules]);
      fHistoMotherBackInvMassECalibSM[iCut][iModules]     = new TH2F(Form("ESD_Back_InvMass_E_Calib_SM%i",iModules), Form("ESD_Back_InvMass_E_Calib_SM%i",iModules), fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
      fHistoMotherBackInvMassECalibSM[iCut][iModules]->SetXTitle("M_{inv} (GeV/c^{2})");
      fHistoMotherBackInvMassECalibSM[iCut][iModules]->SetYTitle("E_{cluster}(GeV)");
      fESDList[iCut]->Add(fHistoMotherBackInvMassECalibSM[iCut][iModules]);

      fHistoClusGammaPtSM[iCut][iModules]              = new TH1F(Form("ClusGamma_Pt_SM%i",iModules), Form("ClusGamma_Pt_SM%i; p_{T,clus} (GeV/c)",iModules), nBinsClusterPt, arrClusPtBinning);
      fESDList[iCut]->Add(fHistoClusGammaPtSM[iCut][iModules]);
      fHistoClusGammaESM[iCut][iModules]               = new TH1F(Form("ClusGamma_E_SM%i",iModules), Form("ClusGamma_E_SM%i; p_{T,clus} (GeV/c)",iModules), nBinsClusterPt, arrClusPtBinning);
      fESDList[iCut]->Add(fHistoClusGammaESM[iCut][iModules]);
    }

    if(fDoMesonQA > 0){
      fHistoClusGammaERx[iCut]                     = new TH1F("ClusGamma_E_Rx", "ClusGamma_E_Rx; E_{clus} (GeV/c)", nBinsClusterPt, arrClusPtBinning);
      fESDList[iCut]->Add(fHistoClusGammaERx[iCut]);
      fHistoClusGammaERxNCellCrit[iCut]            = new TH1F("ClusGamma_E_Rx_NCellCrit", "ClusGamma_E_Rx_NCellCrit; E_{clus} (GeV/c)", nBinsClusterPt, arrClusPtBinning);
      fESDList[iCut]->Add(fHistoClusGammaERx[iCut]);
      for(Int_t iModules = 0; iModules < fnModules; iModules++ ){
        fHistoClusGammaERxSM[iCut][iModules]          = new TH1F(Form("ClusGamma_E_Rx_SM%i", iModules), Form("ClusGamma_E_Rx_SM%i; E_{clus} (GeV/c)", iModules), nBinsClusterPt, arrClusPtBinning);
        fESDList[iCut]->Add(fHistoClusGammaERxSM[iCut][iModules]);
        fHistoClusGammaERxNCellCritSM[iCut][iModules] = new TH1F(Form("ClusGamma_E_Rx_NCellCrit_SM%i", iModules), Form("ClusGamma_E_Rx_NCellCrit_SM%i; E_{clus} (GeV/c)", iModules), nBinsClusterPt, arrClusPtBinning);
        fESDList[iCut]->Add(fHistoClusGammaERxNCellCritSM[iCut][iModules]);


        fHistoClusTrackdEtaSM[iCut][iModules] = new TH1F(Form("Clus_Track_dEta_SM%i", iModules), Form("Clus_Track_dEta_SM%i; d#eta", iModules), 200, -0.01, 0.01);
        fESDList[iCut]->Add(fHistoClusTrackdEtaSM[iCut][iModules]);
        fHistoClusTrackdPhiSM[iCut][iModules] = new TH1F(Form("Clus_Track_dPhi_SM%i", iModules), Form("Clus_Track_dPhi_SM%i; d#eta", iModules), 200, -0.01, 0.01);
        fESDList[iCut]->Add(fHistoClusTrackdPhiSM[iCut][iModules]);
        fHistoClusHighPtTrackdEtaSM[iCut][iModules] = new TH1F(Form("Clus_HighPtTrack_dEta_SM%i", iModules), Form("Clus_HighPtTrack_dEta_SM%i; d#eta", iModules), 200, -0.01, 0.01);
        fESDList[iCut]->Add(fHistoClusHighPtTrackdEtaSM[iCut][iModules]);
        fHistoClusHighPtTrackdPhiSM[iCut][iModules] = new TH1F(Form("Clus_HighPtTrack_dPhi_SM%i", iModules), Form("Clus_HighPtTrack_dPhi_SM%i; d#eta", iModules), 200, -0.01, 0.01);
        fESDList[iCut]->Add(fHistoClusHighPtTrackdPhiSM[iCut][iModules]);
      }
      fHistoEVsM02[iCut]                       = new TH2F("ESD_ClusterE_M02", "ESD_ClusterE_M02", nBinsClusterPt, arrClusPtBinning, 200, 0, 2);
      fHistoEVsM02[iCut]->SetXTitle("E_{cluster}(GeV)");
      fHistoEVsM02[iCut]->SetYTitle("M_{02}");
      fESDList[iCut]->Add(fHistoEVsM02[iCut]);
      fHistoEVsM02NCell4[iCut]                       = new TH2F("ESD_ClusterE_M02_NCell4", "ESD_ClusterE_M02_NCell4", nBinsClusterPt, arrClusPtBinning, 200, 0, 2);
      fHistoEVsM02NCell4[iCut]->SetXTitle("E_{cluster}(GeV)");
      fHistoEVsM02NCell4[iCut]->SetYTitle("M_{02}");
      fESDList[iCut]->Add(fHistoEVsM02NCell4[iCut]);

    }

    fHistoMotherInvMassECalib[iCut]                       = new TH2F("ESD_Mother_InvMass_E_Calib", "ESD_Mother_InvMass_E_Calib", fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
    fHistoMotherInvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
    fHistoMotherInvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
    fESDList[iCut]->Add(fHistoMotherInvMassECalib[iCut]);
    fHistoMotherBackInvMassECalib[iCut]                   = new TH2F("ESD_Back_InvMass_E_Calib", "ESD_Back_InvMass_E_Calib", fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
    fHistoMotherBackInvMassECalib[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
    fHistoMotherBackInvMassECalib[iCut]->SetYTitle("E_{cluster}(GeV)");
    fESDList[iCut]->Add(fHistoMotherBackInvMassECalib[iCut]);


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

      if(fMesonRecoMode > 0){
        fHistoEVsNCellsInPiMass[iCut]       = new TH2F("ESD_EVsNcell_Pi0Tagged", "ESD_EVsNcell_Pi0Tagged; #it{E}_{clus} (GeV); #it{N}_{cells}",
                                                         nBinsPt, arrPtBinning, 25, 0, 25);
       fESDList[iCut]->Add(fHistoEVsNCellsInPiMass[iCut]);
      }
    }

    fHistoMotherBackInvMassPt[iCut]         = new TH2F("ESD_Background_InvMass_Pt", "ESD_Background_InvMass_Pt; M_{inv, mxed}(GeV/c^{2}); p_{T,BG pair} (GeV/c)",
                                                       fMesonInvMassNBins, fMesonInvMassMin, fMesonInvMassMax, nBinsPt, arrPtBinning);
    fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);

    if (fIsMC > 1){
      fHistoMotherInvMassPt[iCut]->Sumw2();
      fHistoMotherBackInvMassPt[iCut]->Sumw2();
      if (fMesonRecoMode == 1){
        if (fHistoMotherMatchedInvMassPt[iCut]) fHistoMotherMatchedInvMassPt[iCut]->Sumw2();
      }
      if (fMesonRecoMode > 0){
        if (fHistoMotherInvMassECalib[iCut]) fHistoMotherInvMassECalib[iCut]->Sumw2();
        if (fHistoMotherBackInvMassECalib[iCut]) fHistoMotherBackInvMassECalib[iCut]->Sumw2();
        if (fHistoEVsNCellsInPiMass[iCut]) fHistoEVsNCellsInPiMass[iCut]->Sumw2();
        if (fHistoClusGammaERxNCellCrit[iCut])  fHistoClusGammaERxNCellCrit[iCut]->Sumw2();
        if (fHistoEVsM02[iCut])  fHistoEVsM02[iCut]->Sumw2();
        if (fHistoEVsM02NCell4[iCut])  fHistoEVsM02NCell4[iCut]->Sumw2();
        for(Int_t iModules = 0; iModules < fnModules; iModules++ ){
          if(fHistoClusGammaERxSM[iCut][iModules])fHistoClusGammaERxSM[iCut][iModules]->Sumw2();
          if(fHistoClusGammaERxNCellCritSM[iCut][iModules])fHistoClusGammaERxNCellCritSM[iCut][iModules]->Sumw2();
          if (fHistoMotherInvMassECalibSM[iCut][iModules]) fHistoMotherInvMassECalibSM[iCut][iModules]->Sumw2();
          if (fHistoMotherBackInvMassECalibSM[iCut][iModules]) fHistoMotherBackInvMassECalibSM[iCut][iModules]->Sumw2();
          if(fHistoClusTrackdEtaSM[iCut][iModules]) fHistoClusTrackdEtaSM[iCut][iModules]->Sumw2();
          if(fHistoClusTrackdPhiSM[iCut][iModules]) fHistoClusTrackdPhiSM[iCut][iModules]->Sumw2();
          if(fHistoClusHighPtTrackdEtaSM[iCut][iModules]) fHistoClusHighPtTrackdEtaSM[iCut][iModules]->Sumw2();
          if(fHistoClusHighPtTrackdPhiSM[iCut][iModules]) fHistoClusHighPtTrackdPhiSM[iCut][iModules]->Sumw2();
        }
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

  for(Int_t iMatcherTask = 0; iMatcherTask < 5; iMatcherTask++){
    AliCaloTrackMatcher* temp = 0x0;
    if(!fCorrTaskSetting.CompareTo("")){
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
    } else {
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i_%s",iMatcherTask,fTrackMatcherRunningMode,fCorrTaskSetting.Data())));
    }
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }

  // for electron selection
  if(fUseEletronMatchingCalibration == 1){
    fElecSelector=(AliDalitzElectronSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ElectronSelector");
    if(!fElecSelector){AliFatal("Error: No ElectronSelector");}

    if( fElecSelector ){
      if ( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() ){
        fOutputContainer->Add( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() );
      }
    }
  }


    //********************************************************************************************************//
    //*****************************  NOT NEEDED FOR CALIBRATION???    ****************************************//
    //********************************************************************************************************//
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


  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskConvCaloCalibration::Notify(){

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
void AliAnalysisTaskConvCaloCalibration::UserExec(Option_t *){
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

      for(Int_t iCut = 0; iCut<fnCuts; iCut++){
        fHistoNEvents[iCut]->Fill(eventQuality);
      }
      return;
    }
//     if(fInputEvent->IsA()==AliAODEvent::Class()){
//       fInputEvent->InitMagneticField();
//     }
    fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

    if(fUseEletronMatchingCalibration == 1){
      fSelectorElectronIndex = fElecSelector->GetReconstructedElectronsIndex(); // Electrons from default Cut
      fSelectorPositronIndex = fElecSelector->GetReconstructedPositronsIndex(); // Positrons from default Cut
    } else if(fUseEletronMatchingCalibration == 2){
      GetV0Electrons();
    }

    // ------------------- BeginEvent ----------------------------
    AliEventplane *EventPlane = fInputEvent->GetEventplane();
    if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
    else fEventPlaneAngle=0.0;

    for(Int_t iCut = 0; iCut<fnCuts; iCut++){

      fiCut = iCut;
      Bool_t isRunningEMCALrelAna = kFALSE;
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger){
        if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
      }

      if (fIsMC > 0){
        fWeightJetJetMC       = 1;
        Float_t pthard = -1;
        Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC , pthard, fInputEvent);
        if (fIsMC == 3){
          Double_t weightMult   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
          fWeightJetJetMC       = fWeightJetJetMC*weightMult;
        }

        if (!isMCJet){
          fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
          if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(10);
          continue;
        }
      }

      Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,isRunningEMCALrelAna);

      if(eventNotAccepted!= 0){
        fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
        // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
        continue;
      }
      if(eventQuality != 0){// Event Not Accepted
        //cout << "event rejected due to: " <<eventQuality << endl;
        fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
        continue;
      }

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

      if (fMesonRecoMode < 2){
        if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetDoElecDeDxPostCalibration()){
          if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->LoadElecDeDxPostCalibration(fInputEvent->GetRunNumber())){
            AliFatal(Form("ERROR: LoadElecDeDxPostCalibration returned kFALSE for %d despite being requested!",fInputEvent->GetRunNumber()));
          }
        }
      }


      // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger) // process calo clusters
        ProcessClusters();
      if (fMesonRecoMode < 2) // Process this cuts gammas
        ProcessPhotonCandidates();
      if (fMesonRecoMode < 2) fHistoNGammaConvCandidates[iCut]->Fill(fGammaCandidates->GetEntries(),fWeightJetJetMC);
      if (fMesonRecoMode > 0 || fEnableClusterCutsForTrigger) fHistoNGammaCaloCandidates[iCut]->Fill(fClusterCandidates->GetEntries(),fWeightJetJetMC);


      // check gamma gamma pairs and veto if necessary
      if (!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetAcceptMassFlag())
        SetPhotonVeto();

      CalculateMesonCandidates(); // Combine Gammas from conversion and from calo
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          CalculateBackground(); // Combinatorial Background
          UpdateEventByEventData(); // Store Event for mixed Events
        } else if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 2){
          CalculateBackgroundSwapp(); // Combinatorial Background
        } else{
          CalculateBackgroundRP(); // Combinatorial Background
          fBGHandlerRP[iCut]->AddEvent(fGammaCandidates,fInputEvent); // Store Event for mixed Events
          fBGClusHandlerRP[iCut]->AddEvent(fClusterCandidates,fInputEvent); // Store Event for mixed Events
        }
      }
      fGammaCandidates->Clear(); // delete this cuts good gammas
      fClusterCandidates->Clear(); // delete cluster candidates
    }

    // if(fIsMC>0 && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled()) && fMesonRecoMode < 2){
    //   RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    //   fV0Reader->RelabelAODs(kFALSE);
    // }
    fSelectorElectronIndex.clear();
    fSelectorPositronIndex.clear();
    fV0Electrons.clear();
    PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::ProcessClusters(){
  Int_t nclus = 0;
  TClonesArray * arrClustersProcess = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskConvCaloCalibration! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersProcess->GetEntries();
  }

  // check if esd or aod
  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(fInputEvent);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return;
    }
  }

    // cout << nclus << endl;
  vector<AliAODConversionPhoton*>         vectorCurrentClusters;
  vector<Int_t>                           vectorRejectCluster;
  vector<Int_t>                           vectorClusterSM;
  vector<Double_t>                        vectorPhotonWeight;

  if(nclus == 0)  return;

  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  if(fDoPrimaryTrackMatching && fUseEletronMatchingCalibration == 0) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kFALSE);
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
    Bool_t IsClusSelected = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i);
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedBeforeTrackMatch()){
      // get track matching histos for eta and phi for every SM
      AliCaloTrackMatcher* trackMatcher = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloTrackMatcherInstance();
      if(!trackMatcher) AliFatal("AliAnalysisTaskConvCaloCalibration::ProcessClusters()  Track matcher instance could not be retrieved...\n");
      // loop over all tracks
      for (Int_t itr=0;itr<fInputEvent->GetNumberOfTracks();itr++){
        AliVTrack *inTrack = 0x0;
        if(esdev){
          inTrack = esdev->GetTrack(itr);
          if(!inTrack) continue;
          AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);
          // if(!isEMCalOnly){ //match only primaries for hybrid reconstruction schemes
          //    if(!EsdTrackCuts->AcceptTrack(esdt)) continue;
          // }
          const AliExternalTrackParam *in = esdt->GetInnerParam();
          if (!in){AliDebug(2, "Could not get InnerParam of Track, continue");continue;}
        } else if(aodev) {
          inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
          if(!inTrack) continue;
          AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);
          // if(!isEMCalOnly){ //match only primaries for hybrid reconstruction schemes
          if(!aodt->IsHybridGlobalConstrainedGlobal()) continue;
          if(TMath::Abs(aodt->Eta())>0.8) continue;
          if(aodt->Pt()<0.15) continue;
          // }
          if(trackMatcher->GetRunningMode()==7){
            // if you considered negative ID hybrid tracks, one needs to make sure
            // that one only asks fCaloTrackMatcher->GetTrackClusterMatchingResidual()
            // for hybrid tracks as well. Otherwise you will find multiple matches
            // since this loop contains dublicates when not requiring Is
            if(!aodt->IsHybridGlobalConstrainedGlobal()) continue;
            if(TMath::Abs(aodt->Eta())>0.8) continue;
            if(aodt->Pt()<0.15) continue;
          }
        } // end condition AOD/ESD

        Float_t dEta, dPhi;
        if(!trackMatcher->GetTrackClusterMatchingResidual(inTrack->GetID(),clus->GetID(),dEta,dPhi)){
          continue;
        }
        if(!fGeomEMCAL)fGeomEMCAL = AliEMCALGeometry::GetInstance();
        Int_t SupMod = fGeomEMCAL->GetSuperModuleNumber(clus->GetCellAbsId(0));
        fHistoClusTrackdEtaSM[fiCut][SupMod]->Fill(dEta, fWeightJetJetMC);
        fHistoClusTrackdPhiSM[fiCut][SupMod]->Fill(dPhi, fWeightJetJetMC);
        if(inTrack->Pt()>5){
          fHistoClusHighPtTrackdEtaSM[fiCut][SupMod]->Fill(dEta, fWeightJetJetMC);
          fHistoClusHighPtTrackdPhiSM[fiCut][SupMod]->Fill(dPhi, fWeightJetJetMC);
        }
      } // end loop over tracks
    }

    if(!IsClusSelected){
      delete clus;
      continue;
    }

    if(fUseEletronMatchingCalibration == 1){
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchElectronTracksToClusters(fInputEvent, fMCEvent, clus, fIsMC, fSelectorElectronIndex, fWeightJetJetMC);
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchElectronTracksToClusters(fInputEvent, fMCEvent, clus, fIsMC, fSelectorPositronIndex, fWeightJetJetMC);
    } else if(fUseEletronMatchingCalibration == 2){
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchElectronTracksToClusters(fInputEvent, fMCEvent, clus, fIsMC, fV0Electrons, fWeightJetJetMC);
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


    fIsFromDesiredHeader            = kTRUE;
    fIsOverlappingWithOtherHeader   = kFALSE;
    if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders)){

      vectorCurrentClusters.push_back(PhotonCandidate);
      vectorPhotonWeight.push_back(fWeightJetJetMC);

      if(PhotonCandidate->GetIsCaloPhoton() == 1 || PhotonCandidate->GetIsCaloPhoton() == 3 ){
        fGeomEMCAL = AliEMCALGeometry::GetInstance();
        vectorClusterSM.push_back(fGeomEMCAL->GetSuperModuleNumber(PhotonCandidate->GetLeadingCellID()));
      }else {
        vectorClusterSM.push_back(0);
      }

      // Calculate R for cross talk studies
      if(fDoMesonQA > 0){
        if(clus->GetM02() > 0.1 && clus->GetM02() < 0.3){
          Int_t iSM = fGeomEMCAL->GetSuperModuleNumber(PhotonCandidate->GetLeadingCellID());
          if(iSM >= 0 && iSM < 20){
            fHistoClusGammaERxSM[fiCut][iSM]->Fill(PhotonCandidate->E(), fWeightJetJetMC);
          }
          fHistoClusGammaERx[fiCut]->Fill(PhotonCandidate->E(), fWeightJetJetMC);
          if(clus->GetNCells() > 4){
            fHistoClusGammaERxNCellCrit[fiCut]->Fill(PhotonCandidate->E(), fWeightJetJetMC);
            if(iSM >= 0 && iSM < 20){
              fHistoClusGammaERxNCellCritSM[fiCut][iSM]->Fill(PhotonCandidate->E(), fWeightJetJetMC);
            }
          }
        }
        fHistoEVsM02[fiCut]->Fill(PhotonCandidate->E(), clus->GetM02(), fWeightJetJetMC);
        if(clus->GetNCells() > 4){
          fHistoEVsM02NCell4[fiCut]->Fill(PhotonCandidate->E(), clus->GetM02(), fWeightJetJetMC);
        }
      }

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

      if(vectorClusterSM.at(iter) >= 0 && vectorClusterSM.at(iter) < 20){
        fHistoClusGammaPtSM[fiCut][vectorClusterSM.at(iter)]->Fill(vectorCurrentClusters.at(iter)->Pt(),vectorPhotonWeight.at(iter));
        fHistoClusGammaESM[fiCut][vectorClusterSM.at(iter)]->Fill(vectorCurrentClusters.at(iter)->Pt(),vectorPhotonWeight.at(iter));
      }


      fClusterCandidates->Add(vectorCurrentClusters.at(iter));
    }
  }
  vectorRejectCluster.clear();
  vectorClusterSM.clear();
  vectorPhotonWeight.clear();
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::ProcessPhotonCandidates(){

  Int_t nV0 = 0;
  TList *GammaCandidatesStepOne = new TList();
  TList *GammaCandidatesStepTwo = new TList();
  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    fIsFromDesiredHeader = kTRUE;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
        !((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

      if(fIsFromDesiredHeader){
        if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
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
      // if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      //   Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      //   Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
      //   if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromDesiredHeader = kFALSE;
      // }
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
      if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
        fGammaCandidates->Add(PhotonCandidate);
        if(fIsFromDesiredHeader){
          if(!fDoLightOutput) fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt(),fWeightJetJetMC);
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
    }
  }

  delete GammaCandidatesStepOne;
  GammaCandidatesStepOne = 0x0;
  delete GammaCandidatesStepTwo;
  GammaCandidatesStepTwo = 0x0;

}

//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::GetV0Electrons(){

  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;

    //apply cuts to maximize electron purity
    if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent))continue;

    for (Int_t iElec = 0;iElec < 2;iElec++){
      Int_t tracklabel = PhotonCandidate->GetLabel(iElec);
      AliVTrack *inTrack = 0x0;
      if( fInputEvent->IsA()==AliESDEvent::Class() ) {
        if( tracklabel > fInputEvent->GetNumberOfTracks() ) continue;
        inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(tracklabel));
      } else {
        if( ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled() ){
          inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(tracklabel));
        } else {
          for( Int_t ii=0;ii<fInputEvent->GetNumberOfTracks();ii++ ) {
            inTrack = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(ii));
            if(inTrack){
              if(inTrack->GetID() == tracklabel) {
                break;
              }
            }
          }
        }
      }
      fV0Electrons.push_back(tracklabel);
    }

  }

}

//________________________________________________________________________
// function to reject photons in specific invariant mass window
//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::SetPhotonVeto(){
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
void AliAnalysisTaskConvCaloCalibration::CalculateMesonCandidates(){
  TClonesArray * arrClustersMesonCand = NULL;
  if(fCorrTaskSetting.CompareTo("") && fMesonRecoMode > 0)
    arrClustersMesonCand = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));

  // AliVCaloCells* cells = 0x0;
  // cells = fInputEvent->GetEMCALCells();
  Int_t nMod = -1;
  Short_t cellNumber = 0;
  AliVCluster* cluster = NULL;
  AliVCluster* Cluster0 = NULL;
  AliVCluster* Cluster1 = NULL;
  fGeomEMCAL = AliEMCALGeometry::GetInstance();
  if (fMesonRecoMode == 1){ //PCM-CALO
    if(fGammaCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;

        for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          Bool_t matched = kFALSE;
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL) continue;

          if (gamma1->GetIsCaloPhoton() > 0){

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
          }

          matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent, fWeightJetJetMC);
          // if(arrClustersMesonCand) delete cluster;

          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
            if (matched){
              if(fHistoMotherMatchedInvMassPt[fiCut]) fHistoMotherMatchedInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
            }


            if (!matched){

              if (!gamma0->GetUseForMesonPair() || !gamma1->GetUseForMesonPair()){
                fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
                delete mesonCand;
                mesonCand=0x0;
                continue;
              }
              if (fHistoMotherInvMassPt[fiCut]){
                fHistoMotherInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(),fWeightJetJetMC);
                cellNumber  = cluster->GetCellsAbsId()[0];
                nMod = fGeomEMCAL->GetSuperModuleNumber(cellNumber);
                fHistoMotherInvMassECalibSM[fiCut][nMod]->Fill(mesonCand->M(),gamma1->E(),fWeightJetJetMC);
                fHistoMotherInvMassECalib[fiCut]->Fill(mesonCand->M(),gamma1->E(),fWeightJetJetMC);
               }


              if (fDoMesonQA > 0){
                if (mesonCand->M() > fMesonInvMassWindow[0] && mesonCand->M() < fMesonInvMassWindow[1]){
                  if (fHistoMotherMesonPtY[fiCut]) fHistoMotherMesonPtY[fiCut]->Fill(mesonCand->Pt(),mesonCand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),fWeightJetJetMC);
                  if (fHistoMotherMesonPtAlpha[fiCut]) fHistoMotherMesonPtAlpha[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetAlpha(),fWeightJetJetMC);
                  if (fHistoMotherMesonPtOpenAngle[fiCut]) fHistoMotherMesonPtOpenAngle[fiCut]->Fill(mesonCand->Pt(),mesonCand->GetOpeningAngle(),fWeightJetJetMC);
                  if (fHistoMotherMesonConvPhotonEtaPhi[fiCut]) fHistoMotherMesonConvPhotonEtaPhi[fiCut]->Fill(gamma0->GetPhotonPhi(), gamma0->GetPhotonEta(),fWeightJetJetMC);
                }
              }

              if(!fDoLightOutput && (mesonCand->M() > 0.09 && mesonCand->M() < 0.17)) {
                fHistoEVsNCellsInPiMass[fiCut]->Fill(gamma1->E(), cluster->GetNCells(), fWeightJetJetMC);
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
                }else{
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
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  } else if (fMesonRecoMode == 2){   //CALO-CALO
    if(fClusterCandidates->GetEntries()>0){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
        AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
        if (gamma0==NULL) continue;
        for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
          AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
          if (gamma1==NULL) continue;

          AliAODConversionMother *mesonCand = new AliAODConversionMother(gamma0,gamma1);
          mesonCand->SetLabels(firstGammaIndex,secondGammaIndex);

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesonCand, kTRUE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), gamma0->GetLeadingCellID(), gamma1->GetLeadingCellID(),gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton() ))){
            if (!gamma0->GetUseForMesonPair() || !gamma1->GetUseForMesonPair()){
              fHistoMotherInvMassRejected[fiCut]->Fill(mesonCand->M());
              delete mesonCand;
              mesonCand=0x0;
              continue;
            }
            if (fHistoMotherInvMassPt[fiCut]) fHistoMotherInvMassPt[fiCut]->Fill(mesonCand->M(),mesonCand->Pt(), fWeightJetJetMC);
            // fill new histograms



            if(fInputEvent->IsA()==AliESDEvent::Class()){
              if(arrClustersMesonCand){
                Cluster0 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(gamma0->GetCaloClusterRef()));
                Cluster1 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
              }else{
                Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              }
            } else if(fInputEvent->IsA()==AliAODEvent::Class()){
              if(arrClustersMesonCand){
                Cluster0 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(gamma0->GetCaloClusterRef()));
                Cluster1 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMesonCand->At(gamma1->GetCaloClusterRef()));
              }else{
                Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              }
            }

            Double_t tempmesonCandWeight = 1;   // for JetJet MC
            //Get SM for Cluster0 & Cluster1
            Int_t iCellClus0 = Cluster0->GetCellsAbsId()[0];
            Int_t iCellClus1 = Cluster1->GetCellsAbsId()[0];
            if(TMath::Abs(mesonCand->GetAlpha())<0.1 ){
              fHistoMotherInvMassECalib[fiCut]->Fill(mesonCand->M(),mesonCand->E(),tempmesonCandWeight);
              if(fGeomEMCAL->GetSuperModuleNumber(iCellClus0) == fGeomEMCAL->GetSuperModuleNumber(iCellClus1)){
                fHistoMotherInvMassECalibSM[fiCut][fGeomEMCAL->GetSuperModuleNumber(iCellClus0)]->Fill(mesonCand->M(),mesonCand->E(),tempmesonCandWeight);
              }
            }

            if(!fDoLightOutput && (mesonCand->M() > 0.09 && mesonCand->M() < 0.17)) {
              fHistoEVsNCellsInPiMass[fiCut]->Fill(gamma0->E(), Cluster0->GetNCells(), fWeightJetJetMC);
              fHistoEVsNCellsInPiMass[fiCut]->Fill(gamma1->E(), Cluster1->GetNCells(), fWeightJetJetMC);
            }

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
          }
          delete mesonCand;
          mesonCand=0x0;
        }
      }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::CalculateBackground(){

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


  Int_t nMod = -1;
  fGeomEMCAL = AliEMCALGeometry::GetInstance();

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
            // Double_t tempBGCandidateWeight = 1;

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

              if( currentEventGoodV0.GetIsCaloPhoton() > 0 && previousGoodV0.GetIsCaloPhoton() > 0 ){
                if( TMath::Abs(backgroundCandidate->GetAlpha())<0.1){
                  fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->E(),fWeightJetJetMC);
                  if(fGeomEMCAL->GetSuperModuleNumber(currentEventGoodV0.GetLeadingCellID()) == fGeomEMCAL->GetSuperModuleNumber(previousGoodV0.GetLeadingCellID())){
                    nMod = fGeomEMCAL->GetSuperModuleNumber(currentEventGoodV0.GetLeadingCellID());
                    fHistoMotherBackInvMassECalibSM[fiCut][nMod]->Fill(backgroundCandidate->M(),backgroundCandidate->E(),fWeightJetJetMC);
                  }
                }
               } else{
                  if( currentEventGoodV0.GetIsCaloPhoton() > 0 ){
                    nMod = fGeomEMCAL->GetSuperModuleNumber(currentEventGoodV0.GetLeadingCellID());
                    fHistoMotherBackInvMassECalibSM[fiCut][nMod]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(),fWeightJetJetMC);
                    fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),currentEventGoodV0.E(),fWeightJetJetMC);
                  }
                  else if(previousGoodV0.GetIsCaloPhoton() > 0){
                    nMod = fGeomEMCAL->GetSuperModuleNumber(previousGoodV0.GetLeadingCellID());
                    fHistoMotherBackInvMassECalibSM[fiCut][nMod]->Fill(backgroundCandidate->M(),previousGoodV0.E(),fWeightJetJetMC);
                    fHistoMotherBackInvMassECalib[fiCut]->Fill(backgroundCandidate->M(),previousGoodV0.E(),fWeightJetJetMC);
                  }
               }



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
void AliAnalysisTaskConvCaloCalibration::CalculateBackgroundSwapp(){

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {
    // Calo mode
    if (fMesonRecoMode == 2){

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

      // curcial requierment is that the event has at least 3 cluster candidates
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


              cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));
              cellIDRotatedPhoton2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()));

              if(!fDoLightOutput){
                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg()))){
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), 1);
                }
                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg()))){
                  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()), 1);
                }
              }

              std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
              std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));

              for(auto const& kCurrentClusterCandidates  : *fClusterCandidates){
                if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentClusterCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentClusterCandidates)){ continue;}

                std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));
                std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));

                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg())) && lvRotationPhoton1.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
                {
                  if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton1, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID()))
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                    vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                  }
                }
                if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg())) && lvRotationPhoton2.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
                {
                  if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton2, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID()))
                  {
                    vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                    vSwappingInvMassPTAlphaCut.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
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
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoWeightingInSwappBg() && vSwappingInvMassPTAlphaCut.size() > 0){
          tempMultWeightSwapping = (0.5*(fClusterCandidates->GetEntries()*fClusterCandidates->GetEntries() - fClusterCandidates->GetEntries()))/(vSwappingInvMassPTAlphaCut.size());
        }
        for(Int_t i = 0; i < (Int_t)vSwappingInvMassPTAlphaCut.size(); i++){
          fHistoMotherBackInvMassECalib[fiCut]->Fill(vSwappingInvMassPTAlphaCut.at(i)[0], vSwappingInvMassPTAlphaCut.at(i)[1],tempMultWeightSwapping*fWeightJetJetMC);
        }
      }
      // PCM-Calo case
    } else if(fMesonRecoMode == 1){
      Double_t rotationAngle = TMath::Pi()/2.0;

      TLorentzVector lvRotationPhoton1;   // photon candidates which get rotated
      TLorentzVector lvRotationPhoton2;   // photon candidates which get rotated
      TVector3 lvRotationPion;            // reconstructed mother particle from the two photons

      Int_t cellIDRotatedPhoton = -1; // cell ID of the cluster after rotation

      std::vector<std::array<Double_t, 2>> vSwappingInvMassPT;
      std::vector<std::array<Double_t, 2>> vSwappingInvMassPTECalib;
      std::vector<std::array<Double_t, 2>> vSwappingInvMassPTAlphaCut;
      vSwappingInvMassPT.clear();
      vSwappingInvMassPTAlphaCut.clear();
      vSwappingInvMassPTECalib.clear();
      vSwappingInvMassPT.resize(0);
      vSwappingInvMassPTAlphaCut.resize(0);
      vSwappingInvMassPTECalib.resize(0);

      Double_t tempMultWeightSwapping = 1; // weight taking multiplicity of event into account

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
            cellIDRotatedPhoton = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), lvRotationPhoton1.Phi());
            if(!fDoLightOutput){
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg()))){
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
                  vSwappingInvMassPTECalib.push_back({backgroundCandidate2->M(),currentEventGoodV0Temp1->E()});
                }
              }
            }
            for(auto const& kCurrentGammaCandidates  : *fGammaCandidates){
              if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentGammaCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentGammaCandidates)) continue;

              std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));
              // std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentGammaCandidates)));

              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg())))
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                  vSwappingInvMassPTECalib.push_back({backgroundCandidate1->M(),currentEventGoodV0Temp1->E()});
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
        fHistoMotherBackInvMassECalib[fiCut]->Fill(vSwappingInvMassPTECalib.at(i)[0], vSwappingInvMassPTECalib.at(i)[1],tempMultWeightSwapping*fWeightJetJetMC);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::CalculateBackgroundRP(){

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
          if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(&backgroundCandidate, kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), gamma0->GetLeadingCellID(), gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton())){
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
                ->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),gamma0->GetLeadingCellID(), gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton())){
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
void AliAnalysisTaskConvCaloCalibration::RotateParticle(AliAODConversionPhoton *gamma){
  Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
  Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
  Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){
  previousEventEP=previousEventEP+TMath::Pi();
  thisEventEP=thisEventEP+TMath::Pi();
  Double_t rotationValue= thisEventEP-previousEventEP;
  gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
  //see header file for documentation

  Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
  Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
  Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

  Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
  particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::UpdateEventByEventData(){
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
void AliAnalysisTaskConvCaloCalibration::RelabelAODPhotonCandidates(Bool_t mode){

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
void AliAnalysisTaskConvCaloCalibration::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskConvCaloCalibration::Terminate(const Option_t *)
{
  //fOutputContainer->Print(); // Will crash on GRID
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskConvCaloCalibration::CheckVectorOnly(vector<Int_t> &vec, Int_t tobechecked){
  if(tobechecked > -1){
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else return false;
  }
  return false;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskConvCaloCalibration::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked){
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
void AliAnalysisTaskConvCaloCalibration::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
  if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
  else ma[tobechecked] = 2;
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskConvCaloCalibration::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
  map<Int_t, Int_t>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second, fWeightJetJetMC);
  }
  ma.clear();
  return;
}
