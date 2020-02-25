/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock 									              *
 * Version 1                                                              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
// WARNING!: This task is obsolete. Please use AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson instead!
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------


#include <vector>
#include "TParticle.h"
#include "TPDGCode.h"
#include "TMCProcess.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TH2F.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDpidCuts.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliKFParticle.h"
#include "AliMCEventHandler.h"
#include "AliKFVertex.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero.h"
#include "AliCaloTrackMatcher.h"
#include <vector>

ClassImp( AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero():
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fPionSelector(NULL),
  fBGHandlerPiPl(NULL),
  fBGHandlerPiMi(NULL),
  fESDEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fTrueTreeList(NULL),
  fMCList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fSelectorNegPionIndex(0),
  fSelectorPosPionIndex(0),
  fGoodConvGammas(NULL),
  fClusterCandidates(NULL),
  fNeutralPionCandidates(NULL),
  fNeutralPionSidebandCandidates(NULL),
  fPosPionCandidates(NULL),
  fNegPionCandidates(NULL),
  fGoodVirtualParticles(NULL),
  fEventCutArray(NULL),
  fGammaCutArray(NULL),
  fClusterCutArray(NULL),
  fPionCutArray(NULL),
  fNeutralPionMesonCutArray(NULL),
  fMesonCutArray(NULL),
  fEventCuts(NULL),
  fConversionCuts(NULL),
  fClusterCuts(NULL),
  fTreePiPiSameMother(NULL),
  fTreePiPiPiSameMother(NULL),
  fTreeEventInfoOmega(NULL),
  fTreeEventInfoEta(NULL),
  fCasePiPi(-1),
  fSamePiPiMotherID(-1),
  fSamePiPiMotherInvMass(-1),
  fSamePiPiMotherPt(-1),
  fSamePiPiPiMotherID(-1),
  fSamePiPiPiMotherInvMass(-1),
  fSamePiPiPiMotherPt(-1),
  fV0MultiplicityOmegaEvent(-1),
  fTrackMultiplicityOmegaEvent(-1),
  fZVertexOmegaEvent(-1),
  fPtOmega(-1),
  fV0MultiplicityEtaEvent(-1),
  fTrackMultiplicityEtaEvent(-1),
  fZVertexEtaEvent(-1),
  fPtEta(-1),
  fPDGMassPi0(-1),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaEta(NULL),
  fHistoClusterGammaPt(NULL),
  fHistoClusterGammaEta(NULL),
  fHistoNegPionPt(NULL),
  fHistoPosPionPt(NULL),
  fHistoNegPionPhi(NULL),
  fHistoPosPionPhi(NULL),
  fHistoNegPionEta(NULL),
  fHistoPosPionEta(NULL),
  fHistoNegPionClsTPC(NULL),
  fHistoPosPionClsTPC(NULL),
  fHistoPionDCAxy(NULL),
  fHistoPionDCAz(NULL),
  fHistoPionTPCdEdxNSigma(NULL),
  fHistoPionTPCdEdx(NULL),
  fHistoPionPionInvMassPt(NULL),
  fHistoGammaGammaInvMassPt(NULL),
  fHistoGammaGammaInvMassPtBeforeCuts(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherInvMassPtRejectedKinematic(NULL),
  fHistoBackInvMassPtGroup1(NULL),
  fHistoBackInvMassPtGroup2(NULL),
  fHistoBackInvMassPtGroup3(NULL),
  fHistoBackInvMassPtGroup4(NULL),
  fHistoMotherLikeSignBackInvMassPt(NULL),
  fHistoAngleOmegaPiPlPiMi(NULL),
  fHistoAngleOmegaPiZero(NULL),
  fHistoAngleOmegaPiPl(NULL),
  fHistoAngleOmegaPiMi(NULL),
  fHistoAnglePiPlPiMi(NULL),
  fHistoAnglePiZeroPiMi(NULL),
  fHistoAnglePiPlPiZero(NULL),
  fHistoAngleSum(NULL),
  fHistoTrueAngleSum(NULL),
  fHistoMotherInvMassSubPi0(NULL),
  fHistoBackInvMassPtGroup1SubPi0(NULL),
  fHistoBackInvMassPtGroup2SubPi0(NULL),
  fHistoBackInvMassPtGroup3SubPi0(NULL),
  fHistoBackInvMassPtGroup4SubPi0(NULL),
  fHistoMotherLikeSignBackInvMassSubPi0Pt(NULL),
  fHistoMotherInvMassFixedPzPi0(NULL),
  fHistoBackInvMassPtGroup1FixedPzPi0(NULL),
  fHistoBackInvMassPtGroup2FixedPzPi0(NULL),
  fHistoBackInvMassPtGroup3FixedPzPi0(NULL),
  fHistoBackInvMassPtGroup4FixedPzPi0(NULL),
  fHistoMotherLikeSignBackInvMassFixedPzPi0Pt(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCConvGammaPt(NULL),
  fHistoMCAllPosPionsPt(NULL),
  fHistoMCAllNegPionsPt(NULL),
  fHistoMCGammaFromNeutralMesonPt(NULL),
  fHistoMCPosPionsFromNeutralMesonPt(NULL),
  fHistoMCNegPionsFromNeutralMesonPt(NULL),
  fHistoMCEtaPiPlPiMiPiZeroPt(NULL),
  fHistoMCEtaPiPlPiMiPiZeroInAccPt(NULL),
  fHistoMCOmegaPiPlPiMiPiZeroPt(NULL),
  fHistoMCOmegaPiPlPiMiPiZeroInAccPt(NULL),
  fHistoTrueMotherPiPlPiMiPiZeroInvMassPt(NULL),
  fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt(NULL),
  fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt(NULL),
  fHistoTrueMotherGammaGammaInvMassPt(NULL),
  fHistoTrueMotherGammaGammaFromEtaInvMassPt(NULL),
  fHistoTrueMotherGammaGammaFromOmegaInvMassPt(NULL),
  fHistoTrueConvGammaPt(NULL),
  fHistoTrueConvGammaFromNeutralMesonPt(NULL),
  fHistoTrueClusterGammaPt(NULL),
  fHistoTrueClusterGammaFromNeutralMesonPt(NULL),
  fHistoTruePosPionPt(NULL),
  fHistoTruePosPionFromNeutralMesonPt(NULL),
  fHistoTrueNegPionPt(NULL),
  fHistoTrueNegPionFromNeutralMesonPt(NULL),
  fHistoTruePionPionInvMassPt(NULL),
  fHistoTruePionPionFromSameMotherInvMassPt(NULL),
  fHistoTruePionPionFromEtaInvMassPt(NULL),
  fHistoTruePionPionFromOmegaInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt(NULL),
  fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt(NULL),
  fHistoTruePiPlPiMiPiZeroContaminationInvMassPt(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueOmegaInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueOmegas(0),
  fVectorDoubleCountTrueConvGammas(0),
  fHistoNEvents(NULL),
  fHistoNGoodESDTracks(NULL),
  fProfileEtaShift(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fNumberOfESDTracks(0),
  fMoveParticleAccordingToVertex(kFALSE),
  fIsHeavyIon(0),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsMC(kFALSE),
  fDoLightOutput(kFALSE),
  fNeutralPionMode(0),
  fTolerance(-1),
  fTrackMatcherRunningMode(0)
{

}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero( const char* name ):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fPionSelector(NULL),
  fBGHandlerPiPl(NULL),
  fBGHandlerPiMi(NULL),
  fESDEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fTrueTreeList(NULL),
  fMCList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fSelectorNegPionIndex(0),
  fSelectorPosPionIndex(0),
  fGoodConvGammas(NULL),
  fClusterCandidates(NULL),
  fNeutralPionCandidates(NULL),
  fNeutralPionSidebandCandidates(NULL),
  fPosPionCandidates(NULL),
  fNegPionCandidates(NULL),
  fGoodVirtualParticles(NULL),
  fEventCutArray(NULL),
  fGammaCutArray(NULL),
  fClusterCutArray(NULL),
  fPionCutArray(NULL),
  fNeutralPionMesonCutArray(NULL),
  fMesonCutArray(NULL),
  fEventCuts(NULL),
  fConversionCuts(NULL),
  fClusterCuts(NULL),
  fTreePiPiSameMother(NULL),
  fTreePiPiPiSameMother(NULL),
  fTreeEventInfoOmega(NULL),
  fTreeEventInfoEta(NULL),
  fCasePiPi(-1),
  fSamePiPiMotherID(-1),
  fSamePiPiMotherInvMass(-1),
  fSamePiPiMotherPt(-1),
  fSamePiPiPiMotherID(-1),
  fSamePiPiPiMotherInvMass(-1),
  fSamePiPiPiMotherPt(-1),
  fV0MultiplicityOmegaEvent(-1),
  fTrackMultiplicityOmegaEvent(-1),
  fZVertexOmegaEvent(-1),
  fPtOmega(-1),
  fV0MultiplicityEtaEvent(-1),
  fTrackMultiplicityEtaEvent(-1),
  fZVertexEtaEvent(-1),
  fPtEta(-1),
  fPDGMassPi0(-1),
  fHistoConvGammaPt(NULL),
  fHistoConvGammaEta(NULL),
  fHistoClusterGammaPt(NULL),
  fHistoClusterGammaEta(NULL),
  fHistoNegPionPt(NULL),
  fHistoPosPionPt(NULL),
  fHistoNegPionPhi(NULL),
  fHistoPosPionPhi(NULL),
  fHistoNegPionEta(NULL),
  fHistoPosPionEta(NULL),
  fHistoNegPionClsTPC(NULL),
  fHistoPosPionClsTPC(NULL),
  fHistoPionDCAxy(NULL),
  fHistoPionDCAz(NULL),
  fHistoPionTPCdEdxNSigma(NULL),
  fHistoPionTPCdEdx(NULL),
  fHistoPionPionInvMassPt(NULL),
  fHistoGammaGammaInvMassPt(NULL),
  fHistoGammaGammaInvMassPtBeforeCuts(NULL),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherInvMassPtRejectedKinematic(NULL),
  fHistoBackInvMassPtGroup1(NULL),
  fHistoBackInvMassPtGroup2(NULL),
  fHistoBackInvMassPtGroup3(NULL),
  fHistoBackInvMassPtGroup4(NULL),
  fHistoMotherLikeSignBackInvMassPt(NULL),
  fHistoAngleOmegaPiPlPiMi(NULL),
  fHistoAngleOmegaPiZero(NULL),
  fHistoAngleOmegaPiPl(NULL),
  fHistoAngleOmegaPiMi(NULL),
  fHistoAnglePiPlPiMi(NULL),
  fHistoAnglePiZeroPiMi(NULL),
  fHistoAnglePiPlPiZero(NULL),
  fHistoAngleSum(NULL),
  fHistoTrueAngleSum(NULL),
  fHistoMotherInvMassSubPi0(NULL),
  fHistoBackInvMassPtGroup1SubPi0(NULL),
  fHistoBackInvMassPtGroup2SubPi0(NULL),
  fHistoBackInvMassPtGroup3SubPi0(NULL),
  fHistoBackInvMassPtGroup4SubPi0(NULL),
  fHistoMotherLikeSignBackInvMassSubPi0Pt(NULL),
  fHistoMotherInvMassFixedPzPi0(NULL),
  fHistoBackInvMassPtGroup1FixedPzPi0(NULL),
  fHistoBackInvMassPtGroup2FixedPzPi0(NULL),
  fHistoBackInvMassPtGroup3FixedPzPi0(NULL),
  fHistoBackInvMassPtGroup4FixedPzPi0(NULL),
  fHistoMotherLikeSignBackInvMassFixedPzPi0Pt(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCConvGammaPt(NULL),
  fHistoMCAllPosPionsPt(NULL),
  fHistoMCAllNegPionsPt(NULL),
  fHistoMCGammaFromNeutralMesonPt(NULL),
  fHistoMCPosPionsFromNeutralMesonPt(NULL),
  fHistoMCNegPionsFromNeutralMesonPt(NULL),
  fHistoMCEtaPiPlPiMiPiZeroPt(NULL),
  fHistoMCEtaPiPlPiMiPiZeroInAccPt(NULL),
  fHistoMCOmegaPiPlPiMiPiZeroPt(NULL),
  fHistoMCOmegaPiPlPiMiPiZeroInAccPt(NULL),
  fHistoTrueMotherPiPlPiMiPiZeroInvMassPt(NULL),
  fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt(NULL),
  fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt(NULL),
  fHistoTrueMotherGammaGammaInvMassPt(NULL),
  fHistoTrueMotherGammaGammaFromEtaInvMassPt(NULL),
  fHistoTrueMotherGammaGammaFromOmegaInvMassPt(NULL),
  fHistoTrueConvGammaPt(NULL),
  fHistoTrueConvGammaFromNeutralMesonPt(NULL),
  fHistoTrueClusterGammaPt(NULL),
  fHistoTrueClusterGammaFromNeutralMesonPt(NULL),
  fHistoTruePosPionPt(NULL),
  fHistoTruePosPionFromNeutralMesonPt(NULL),
  fHistoTrueNegPionPt(NULL),
  fHistoTrueNegPionFromNeutralMesonPt(NULL),
  fHistoTruePionPionInvMassPt(NULL),
  fHistoTruePionPionFromSameMotherInvMassPt(NULL),
  fHistoTruePionPionFromEtaInvMassPt(NULL),
  fHistoTruePionPionFromOmegaInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt(NULL),
  fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt(NULL),
  fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt(NULL),
  fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt(NULL),
  fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt(NULL),
  fHistoTruePiPlPiMiPiZeroContaminationInvMassPt(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueOmegaInvMassPt(NULL),
  fHistoDoubleCountTrueConvGammaRPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueOmegas(0),
  fVectorDoubleCountTrueConvGammas(0),
  fHistoNEvents(NULL),
  fHistoNGoodESDTracks(NULL),
  fProfileEtaShift(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fNumberOfESDTracks(0),
  fMoveParticleAccordingToVertex(kFALSE),
  fIsHeavyIon(0),
  fDoMesonAnalysis(kTRUE),
  fDoMesonQA(0),
  fIsFromMBHeader(kTRUE),
  fIsMC(kFALSE),
  fDoLightOutput(kFALSE),
  fNeutralPionMode(0),
  fTolerance(-1),
  fTrackMatcherRunningMode(0)
{
  DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::~AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero()
{
  //
  // virtual destructor
  //
  cout<<"Destructor"<<endl;
  if(fGoodConvGammas){
    delete fGoodConvGammas;
    fGoodConvGammas = 0x0;
  }
  if(fClusterCandidates){
    delete fClusterCandidates;
    fClusterCandidates = 0x0;
  }

  if(fNeutralPionCandidates){
    delete fNeutralPionCandidates;
    fNeutralPionCandidates = 0x0;
  }

  if(fNeutralPionSidebandCandidates){
    delete fNeutralPionSidebandCandidates;
    fNeutralPionSidebandCandidates = 0x0;
  }

  if(fPosPionCandidates){
    delete fPosPionCandidates;
    fPosPionCandidates = 0x0;
  }

  if(fNegPionCandidates){
    delete fNegPionCandidates;
    fNegPionCandidates = 0x0;
  }

  if(fGoodVirtualParticles){
    delete fGoodVirtualParticles;
    fGoodVirtualParticles = 0x0;
  }

  if(fBGHandlerPiPl){
    delete[] fBGHandlerPiPl;
    fBGHandlerPiPl = 0x0;
  }

  if(fBGHandlerPiMi){
    delete[] fBGHandlerPiMi;
    fBGHandlerPiMi = 0x0;
  }
}
//___________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::InitBack(){

  fBGHandlerPiPl = new AliGammaConversionAODBGHandler*[fnCuts];
  fBGHandlerPiMi = new AliGammaConversionAODBGHandler*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){

    TString cutstringEvent		= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPion		= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
    TString cutstringConvGamma = "";
    if (fNeutralPionMode < 2)  cutstringConvGamma = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
    TString cutstringCaloGamma = "";
    if (fNeutralPionMode > 0)  cutstringCaloGamma = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringNeutralPion= ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson		= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    TString fullCutString = "";
    if (fNeutralPionMode == 0) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNeutralPionMode == 1) fullCutString = Form("%i_%s_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNeutralPionMode == 2) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());

    Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
    Int_t centMin = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
    Int_t centMax = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));

    if(collisionSystem == 1 || collisionSystem == 2 ||
      collisionSystem == 5 || collisionSystem == 8 ||
      collisionSystem == 9){
      centMin = centMin*10;
      centMax = centMax*10;
    }
    else if(collisionSystem == 3 || collisionSystem == 6){
      centMin = centMin*5;
      centMax = centMax*5;
    }
    else if(collisionSystem == 4 || collisionSystem == 7){
      centMin = ((centMin*5)+45);
      centMax = ((centMax*5)+45);
    }

    fBGHandlerPiPl[iCut] = new AliGammaConversionAODBGHandler(	collisionSystem,centMin,centMax,
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents(),
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                4,8,5);

    fBGHandlerPiMi[iCut] = new AliGammaConversionAODBGHandler(	collisionSystem,centMin,centMax,
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents(),
                                ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                4,8,5);
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::UserCreateOutputObjects()
{
  //
  // Create ouput objects
  //

  // Set pT and mass ranges
  Double_t HistoMassRange[2]            = {0.4,1.0};
  Double_t HistoMassRangeSub[2]         = {0.4,1.0};
  Double_t HistoPtRange[2]              = {0.,25.};
  Double_t HistoNMassBins               = 600;
  Double_t HistoNMassBinsSub            = 600;
  Double_t HistoNPtBins                 = 250;
  fPDGMassPi0                           = 0.1349766; // hard coded PDG value to keep results reproducable later

  // Create the output container
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer            = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer            = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fGoodConvGammas               = new TList();
  fClusterCandidates            = new TList();
  fClusterCandidates->SetOwner(kTRUE);

  fNeutralPionCandidates        = new TList();
  fNeutralPionCandidates->SetOwner(kTRUE);

  fNeutralPionSidebandCandidates        = new TList();
  fNeutralPionSidebandCandidates->SetOwner(kTRUE);


  fPosPionCandidates            = new TList();
  fPosPionCandidates->SetOwner(kTRUE);
  fNegPionCandidates            = new TList();
  fNegPionCandidates->SetOwner(kTRUE);
  fGoodVirtualParticles         = new TList();
  fGoodVirtualParticles->SetOwner(kTRUE);

  fCutFolder                    = new TList*[fnCuts];
  fESDList                      = new TList*[fnCuts];
  fHistoNEvents                 = new TH1I*[fnCuts];
  fHistoNGoodESDTracks          = new TH1I*[fnCuts];
  if(!fDoLightOutput){
      fProfileEtaShift              = new TProfile*[fnCuts];
      fHistoSPDClusterTrackletBackground = new TH2F*[fnCuts];

      if (fNeutralPionMode < 2){
          fHistoConvGammaPt           = new TH1F*[fnCuts];
          fHistoConvGammaEta          = new TH1F*[fnCuts];
      }
      if (fNeutralPionMode > 0){
          fHistoClusterGammaPt        = new TH1F*[fnCuts];
          fHistoClusterGammaEta       = new TH1F*[fnCuts];
      }
      fHistoNegPionPt               = new TH1F*[fnCuts];
      fHistoPosPionPt               = new TH1F*[fnCuts];
      fHistoNegPionPhi              = new TH1F*[fnCuts];
      fHistoPosPionPhi              = new TH1F*[fnCuts];
      fHistoPionPionInvMassPt       = new TH2F*[fnCuts];

      if( fDoMesonQA>0 ) {
          fHistoNegPionEta            = new TH1F*[fnCuts];
          fHistoPosPionEta            = new TH1F*[fnCuts];
          fHistoNegPionClsTPC         = new TH2F*[fnCuts];
          fHistoPosPionClsTPC         = new TH2F*[fnCuts];
          fHistoPionDCAxy             = new TH2F*[fnCuts];
          fHistoPionDCAz              = new TH2F*[fnCuts];
          fHistoPionTPCdEdxNSigma     = new TH2F*[fnCuts];
          fHistoPionTPCdEdx           = new TH2F*[fnCuts];
      }


      fHistoAngleOmegaPiPlPiMi    = new TH2F*[fnCuts];
      fHistoAngleOmegaPiZero      = new TH2F*[fnCuts];
      fHistoAngleOmegaPiPl        = new TH2F*[fnCuts];
      fHistoAngleOmegaPiMi        = new TH2F*[fnCuts];
      fHistoAnglePiZeroPiMi       = new TH2F*[fnCuts];
      fHistoAnglePiPlPiMi         = new TH2F*[fnCuts];
      fHistoAnglePiPlPiZero       = new TH2F*[fnCuts];
      fHistoAngleSum              = new TH2F*[fnCuts];
  }

  fHistoGammaGammaInvMassPt               = new TH2F*[fnCuts];
  fHistoGammaGammaInvMassPtBeforeCuts      = new TH2F*[fnCuts];
  fHistoMotherInvMassPt                   = new TH2F*[fnCuts];
  fHistoMotherInvMassPtRejectedKinematic  = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup1 = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup2 = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup3  = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup4  = new TH2F*[fnCuts];

  fHistoMotherLikeSignBackInvMassPt       = new TH2F*[fnCuts];

  fHistoMotherInvMassSubPi0                               = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup1SubPi0 = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup2SubPi0 = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup3SubPi0  = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup4SubPi0  = new TH2F*[fnCuts];

  fHistoMotherLikeSignBackInvMassSubPi0Pt       = new TH2F*[fnCuts];

  fHistoMotherInvMassFixedPzPi0                     = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup1FixedPzPi0 = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup2FixedPzPi0 = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup3FixedPzPi0  = new TH2F*[fnCuts];
  fHistoBackInvMassPtGroup4FixedPzPi0  = new TH2F*[fnCuts];

  fHistoMotherLikeSignBackInvMassFixedPzPi0Pt       = new TH2F*[fnCuts];

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPion         = ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
    TString cutstringConvGamma    = "";
    if (fNeutralPionMode < 2)
      cutstringConvGamma          = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
    TString cutstringCaloGamma    = "";
    if (fNeutralPionMode > 0)
      cutstringCaloGamma          = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringNeutralPion  = ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    TString fullCutString         = "";
    if (fNeutralPionMode == 0)
      fullCutString               = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNeutralPionMode == 1)
      fullCutString               = Form("%i_%s_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    else if (fNeutralPionMode == 2)
      fullCutString               = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
    TString nameCutFolder         = Form("Cut Number %s", fullCutString.Data());
    TString nameESDList           = Form("%s ESD histograms", fullCutString.Data());

    fCutFolder[iCut]              = new TList();
    fCutFolder[iCut]->SetName(nameCutFolder.Data());
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);

    fESDList[iCut]                = new TList();
    fESDList[iCut]->SetName(nameESDList.Data());
    fESDList[iCut]->SetOwner(kTRUE);

    fHistoNEvents[iCut]           = new TH1I("NEvents","NEvents",14,-0.5,13.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fHistoNEvents[iCut]->GetYaxis()->SetTitle("N_{events}");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if(fIsHeavyIon>0)
      fHistoNGoodESDTracks[iCut]  = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
    else
      fHistoNGoodESDTracks[iCut]  = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
    fHistoNGoodESDTracks[iCut]->GetXaxis()->SetTitle("N_{good ESD tracks}");
    fHistoNGoodESDTracks[iCut]->GetYaxis()->SetTitle("N_{events}");
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);
    if(!fDoLightOutput){
        fProfileEtaShift[iCut]        = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
        fESDList[iCut]->Add(fProfileEtaShift[iCut]);
        fHistoSPDClusterTrackletBackground[iCut] = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
        fHistoSPDClusterTrackletBackground[iCut]->GetXaxis()->SetTitle("N_{SPD tracklets}");
        fHistoSPDClusterTrackletBackground[iCut]->GetYaxis()->SetTitle("N_{SPD clusters}");
        fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);
        if (fNeutralPionMode < 2){
            fHistoConvGammaPt[iCut]     = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
            fHistoConvGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            fHistoConvGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
            fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
            fHistoConvGammaEta[iCut]    = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",600,-1.5,1.5);
            fHistoConvGammaEta[iCut]->GetXaxis()->SetTitle("#eta");
            fHistoConvGammaEta[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
            fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
        }
        if (fNeutralPionMode > 0){
            fHistoClusterGammaPt[iCut]  = new TH1F("ESD_ClusterGamma_Pt","ESD_ClusterGamma_Pt",250,0,25);
            fHistoClusterGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
            fHistoClusterGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
            fESDList[iCut]->Add(fHistoClusterGammaPt[iCut]);
            fHistoClusterGammaEta[iCut] = new TH1F("ESD_ClusterGamma_Eta","ESD_ClusterGamma_Eta",600,-1.5,1.5);
            fHistoClusterGammaEta[iCut]->GetXaxis()->SetTitle("#eta");
            fHistoClusterGammaEta[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
            fESDList[iCut]->Add(fHistoClusterGammaEta[iCut]);
        }
        fHistoNegPionPt[iCut]         = new TH1F("ESD_PrimaryNegPions_Pt","ESD_PrimaryNegPions_Pt",1000,0,25);
        fHistoNegPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoNegPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
        fESDList[iCut]->Add(fHistoNegPionPt[iCut]);
        fHistoPosPionPt[iCut]         = new TH1F("ESD_PrimaryPosPions_Pt","ESD_PrimaryPosPions_Pt",1000,0,25);
        fHistoPosPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoPosPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
        fESDList[iCut]->Add(fHistoPosPionPt[iCut]);
        fHistoNegPionPhi[iCut]        = new TH1F("ESD_PrimaryNegPions_Phi","ESD_PrimaryNegPions_Phi",360,0,2*TMath::Pi());
        fHistoNegPionPhi[iCut]->GetXaxis()->SetTitle("#phi");
        fHistoNegPionPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
        fESDList[iCut]->Add(fHistoNegPionPhi[iCut]);
        fHistoPosPionPhi[iCut]        = new TH1F("ESD_PrimaryPosPions_Phi","ESD_PrimaryPosPions_Phi",360,0,2*TMath::Pi());
        fHistoPosPionPhi[iCut]->GetXaxis()->SetTitle("#phi");
        fHistoPosPionPhi[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
        fESDList[iCut]->Add(fHistoPosPionPhi[iCut]);
        fHistoPionPionInvMassPt[iCut] = new TH2F("ESD_PiPlusPiNeg_InvMassPt","ESD_PiPlusPiNeg_InvMassPt",2000,0.,2.,200,0.,20.);
        fHistoPionPionInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
        fHistoPionPionInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoPionPionInvMassPt[iCut]);

        if ( fDoMesonQA>0 ) {
            fHistoNegPionEta[iCut]        = new TH1F("ESD_PrimaryNegPions_Eta","ESD_PrimaryNegPions_Eta",600,-1.5,1.5);
            fHistoNegPionEta[iCut]->GetXaxis()->SetTitle("#eta");
            fHistoNegPionEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
            fESDList[iCut]->Add(fHistoNegPionEta[iCut]);
            fHistoPosPionEta[iCut]        = new TH1F("ESD_PrimaryPosPions_Eta","ESD_PrimaryPosPions_Eta",600,-1.5,1.5);
            fHistoPosPionEta[iCut]->GetXaxis()->SetTitle("#eta");
            fHistoPosPionEta[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
            fESDList[iCut]->Add(fHistoPosPionEta[iCut]);
            fHistoNegPionClsTPC[iCut]     = new TH2F("ESD_PrimaryNegPions_ClsTPC","ESD_PrimaryNegPions_ClsTPC",100,0,1,400,0.,10.);
            fHistoNegPionClsTPC[iCut]->GetXaxis()->SetTitle("N_{findable cls. TPC #pi^{-}}");
            fHistoNegPionClsTPC[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
            fESDList[iCut]->Add(fHistoNegPionClsTPC[iCut]);
            fHistoPosPionClsTPC[iCut]     = new TH2F("ESD_PrimaryPosPions_ClsTPC","ESD_PrimaryPosPions_ClsTPC",100,0,1,400,0.,10.);
            fHistoPosPionClsTPC[iCut]->GetXaxis()->SetTitle("N_{findable cls. TPC #pi^{+}}");
            fHistoPosPionClsTPC[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
            fESDList[iCut]->Add(fHistoPosPionClsTPC[iCut]);
            fHistoPionDCAxy[iCut]         = new TH2F("ESD_PrimaryPions_DCAxy","ESD_PrimaryPions_DCAxy",800,-4.0,4.0,400,0.,10.);
            fHistoPionDCAxy[iCut]->GetXaxis()->SetTitle("DCA_{xy}");
            fHistoPionDCAxy[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
            fESDList[iCut]->Add(fHistoPionDCAxy[iCut]);
            fHistoPionDCAz[iCut]          = new TH2F("ESD_PrimaryPions_DCAz","ESD_PrimaryPions_DCAz",800,-4.0,4.0,400,0.,10.);
            fHistoPionDCAz[iCut]->GetXaxis()->SetTitle("DCA_{z}");
            fHistoPionDCAz[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
            fESDList[iCut]->Add(fHistoPionDCAz[iCut]);
            fHistoPionTPCdEdxNSigma[iCut] = new TH2F("ESD_PrimaryPions_TPCdEdx","ESD_PrimaryPions_TPCdEdx",150,0.05,20,400,-10,10);
            fHistoPionTPCdEdxNSigma[iCut]->GetXaxis()->SetTitle("p (GeV/c)");
            fHistoPionTPCdEdxNSigma[iCut]->GetYaxis()->SetTitle("#sigma_{PID,TPC}");
            fESDList[iCut]->Add(fHistoPionTPCdEdxNSigma[iCut]);
            fHistoPionTPCdEdx[iCut]       = new TH2F("ESD_PrimaryPions_TPCdEdxSignal","ESD_PrimaryPions_TPCdEdxSignal" ,150,0.05,20.0,800,0.0,200);
            fHistoPionTPCdEdx[iCut]->GetXaxis()->SetTitle("p (GeV/c)");
            fHistoPionTPCdEdx[iCut]->GetYaxis()->SetTitle("dE/dx signal (au)");
            fESDList[iCut]->Add(fHistoPionTPCdEdx[iCut]);
        }
    fHistoGammaGammaInvMassPt[iCut]               = new TH2F("ESD_GammaGamma_InvMass_Pt","ESD_GammaGamma_InvMass_Pt",600,0.,0.6,250,0,25);
    fHistoGammaGammaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
    fHistoGammaGammaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fESDList[iCut]->Add(fHistoGammaGammaInvMassPt[iCut]);
    fHistoGammaGammaInvMassPtBeforeCuts[iCut]               = new TH2F("ESD_GammaGamma_InvMass_Pt_Before_Cuts","ESD_GammaGamma_InvMass_Pt_Before_Cuts",600,0.,0.6,250,0,25);
    fHistoGammaGammaInvMassPtBeforeCuts[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
    fHistoGammaGammaInvMassPtBeforeCuts[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fESDList[iCut]->Add(fHistoGammaGammaInvMassPtBeforeCuts[iCut]);
  }
    fHistoMotherInvMassPt[iCut]                   = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoMotherInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
    fHistoMotherInvMassPtRejectedKinematic[iCut]  = new TH2F("ESD_Mother_InvMass_Pt_KinematicRejected","ESD_Mother_InvMass_Pt_KinematicRejected",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassPtRejectedKinematic[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoMotherInvMassPtRejectedKinematic[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fESDList[iCut]->Add(fHistoMotherInvMassPtRejectedKinematic[iCut]);

    fHistoBackInvMassPtGroup1[iCut] = new TH2F("ESD_Background_1_InvMass_Pt","ESD_Background_1_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup2[iCut] = new TH2F("ESD_Background_2_InvMass_Pt","ESD_Background_2_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup3[iCut]  = new TH2F("ESD_Background_3_InvMass_Pt","ESD_Background_3_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup4[iCut]  = new TH2F("ESD_Background_4_InvMass_Pt","ESD_Background_4_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup1[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup1[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup2[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup2[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup3[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup3[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup4[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup4[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    if(!(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing())){
       fESDList[iCut]->Add(fHistoBackInvMassPtGroup1[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup2[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup3[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup4[iCut]);
    }

    fHistoMotherLikeSignBackInvMassPt[iCut]  = new TH2F("ESD_Background_LikeSign_InvMass_Pt","ESD_Background_LikeSign_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherLikeSignBackInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{#pm} #pi^{#pm} #pi^{0}} (GeV/c^{2})");
    fHistoMotherLikeSignBackInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()){
        fESDList[iCut]->Add(fHistoMotherLikeSignBackInvMassPt[iCut]);
    }
    fHistoMotherInvMassSubPi0[iCut]                       = new TH2F("ESD_InvMass_Mother_Sub_InvMass(NeutralPion)_Pt","ESD_Mother_InvMass_Sub_InvMass(NeutralPion)_Pt",HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassSubPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} - (M_{#pi^{0}}-M_{#pi^{0},PDG}) (GeV/c^{2})");
    fHistoMotherInvMassSubPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fESDList[iCut]->Add(fHistoMotherInvMassSubPi0[iCut]);

    fHistoBackInvMassPtGroup1SubPi0[iCut]   = new TH2F("ESD_Background_1_InvMass_Sub_InvMass(NeutralPion)_Pt","ESD_Background_1_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                                     HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup2SubPi0[iCut]   = new TH2F("ESD_Background_2_InvMass_Sub_InvMass(NeutralPion)_Pt","ESD_Background_2_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                                     HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup3SubPi0[iCut]    = new TH2F("ESD_Background_3_InvMass_Sub_InvMass(NeutralPion)_Pt","ESD_Background_3_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                                     HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup4SubPi0[iCut]    = new TH2F("ESD_Background_4_InvMass_Sub_InvMass(NeutralPion)_Pt","ESD_Background_4_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                                     HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);

    fHistoBackInvMassPtGroup1SubPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} - (M_{#pi^{0}}-M_{#pi^{0},PDG}) (GeV/c^{2})");
    fHistoBackInvMassPtGroup1SubPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup2SubPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} - (M_{#pi^{0}}-M_{#pi^{0},PDG}) (GeV/c^{2})");
    fHistoBackInvMassPtGroup2SubPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup3SubPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} - (M_{#pi^{0}}-M_{#pi^{0},PDG}) (GeV/c^{2})");
    fHistoBackInvMassPtGroup3SubPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup4SubPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} - (M_{#pi^{0}}-M_{#pi^{0},PDG}) (GeV/c^{2})");
    fHistoBackInvMassPtGroup4SubPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    if(!(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing())){
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup4SubPi0[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup1SubPi0[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup2SubPi0[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup3SubPi0[iCut]);
    }
    fHistoMotherLikeSignBackInvMassSubPi0Pt[iCut]    = new TH2F("ESD_Background_LikeSign_InvMass_Sub_InvMass(NeutralPion)_Pt","ESD_Background_LikeSign_InvMass_Sub_InvMass(NeutralPion)_Pt",
                                                                HistoNMassBinsSub,HistoMassRangeSub[0],HistoMassRangeSub[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherLikeSignBackInvMassSubPi0Pt[iCut]->GetXaxis()->SetTitle("M_{#pi^{#pm} #pi^{#pm} #pi^{0}} - M_{#pi^{0}} (GeV/c^{2})");
    fHistoMotherLikeSignBackInvMassSubPi0Pt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()){
        fESDList[iCut]->Add(fHistoMotherLikeSignBackInvMassSubPi0Pt[iCut]);
    }
    fHistoMotherInvMassFixedPzPi0[iCut]                     = new TH2F("ESD_InvMass_Mother_FixedPz(NeutralPion)_Pt","ESD_Mother_InvMass_FixedPz(NeutralPion)_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherInvMassFixedPzPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoMotherInvMassFixedPzPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fESDList[iCut]->Add(fHistoMotherInvMassFixedPzPi0[iCut]);

    fHistoBackInvMassPtGroup1FixedPzPi0[iCut] = new TH2F("ESD_Background_1_InvMass_FixedPz(NeutralPion)_Pt","ESD_Background_1_InvMass_FixedPz(NeutralPion)_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup2FixedPzPi0[iCut] = new TH2F("ESD_Background_2_InvMass_FixedPz(NeutralPion)_Pt","ESD_Background_2_InvMass_FixedPz(NeutralPion)_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup3FixedPzPi0[iCut]  = new TH2F("ESD_Background_3_InvMass_FixedPz(NeutralPion)_Pt","ESD_Background_3_InvMass_FixedPz(NeutralPion)_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup4FixedPzPi0[iCut]  = new TH2F("ESD_Background_4_InvMass_FixedPz(NeutralPion)_Pt","ESD_Background_4_InvMass_FixedPz(NeutralPion)_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoBackInvMassPtGroup1FixedPzPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup1FixedPzPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup2FixedPzPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup2FixedPzPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup3FixedPzPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup3FixedPzPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    fHistoBackInvMassPtGroup4FixedPzPi0[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
    fHistoBackInvMassPtGroup4FixedPzPi0[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    if(!(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing())){
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup1FixedPzPi0[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup2FixedPzPi0[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup3FixedPzPi0[iCut]);
        fESDList[iCut]->Add(fHistoBackInvMassPtGroup4FixedPzPi0[iCut]);
    }
    fHistoMotherLikeSignBackInvMassFixedPzPi0Pt[iCut]  = new TH2F("ESD_Background_LikeSign_InvMass_FixedPz(NeutralPion)_Pt","ESD_Background_LikeSign_InvMass_FixedPz(NeutralPion)_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
    fHistoMotherLikeSignBackInvMassFixedPzPi0Pt[iCut]->GetXaxis()->SetTitle("M_{#pi^{#pm} #pi^{#pm} #pi^{0}} (GeV/c^{2})");
    fHistoMotherLikeSignBackInvMassFixedPzPi0Pt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
    if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseLikeSignMixing()){
        fESDList[iCut]->Add(fHistoMotherLikeSignBackInvMassFixedPzPi0Pt[iCut]);
    }
    if(!fDoLightOutput){
        fHistoAngleOmegaPiPlPiMi[iCut]      = new TH2F("ESD_Mother_AngleOmegaNegPionsPosPions_Pt","ESD_Mother_AngleOmegaNegPionsPosPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
        fHistoAngleOmegaPiPlPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAngleOmegaPiPlPiMi[iCut]->GetYaxis()->SetTitle("#angle (meson,#pi^{+}#pi^{-})");
        fESDList[iCut]->Add(fHistoAngleOmegaPiPlPiMi[iCut]);
        fHistoAngleOmegaPiMi[iCut]          = new TH2F("ESD_Mother_AngleOmegaNegPions_Pt","ESD_Mother_AngleOmegaNegPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
        fHistoAngleOmegaPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAngleOmegaPiMi[iCut]->GetYaxis()->SetTitle("#angle (meson,#pi^{-})");
        fESDList[iCut]->Add(fHistoAngleOmegaPiMi[iCut]);
        fHistoAngleOmegaPiPl[iCut]          = new TH2F("ESD_Mother_AngleOmegaPosPions_Pt","ESD_Mother_AngleOmegaPosPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
        fHistoAngleOmegaPiPl[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAngleOmegaPiPl[iCut]->GetYaxis()->SetTitle("#angle (meson,#pi^{+})");
        fESDList[iCut]->Add(fHistoAngleOmegaPiPl[iCut]);
        fHistoAngleOmegaPiZero[iCut]        = new TH2F("ESD_Mother_AngleOmegaNeutralPion_Pt","ESD_Mother_AngleOmegaNeutralPion_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
        fHistoAngleOmegaPiZero[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAngleOmegaPiZero[iCut]->GetYaxis()->SetTitle("#angle (meson,#pi^{0})");
        fESDList[iCut]->Add(fHistoAngleOmegaPiZero[iCut]);
        fHistoAnglePiPlPiZero[iCut]         = new TH2F("ESD_Mother_AnglePosPionsNeutralPion_Pt","ESD_Mother_AnglePosPionsNeutralPion_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
        fHistoAnglePiPlPiZero[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAnglePiPlPiZero[iCut]->GetYaxis()->SetTitle("#angle (#pi^{+},#pi^{0})");
        fESDList[iCut]->Add(fHistoAnglePiPlPiZero[iCut]);
        fHistoAnglePiPlPiMi[iCut]           = new TH2F("ESD_Mother_AnglePosPionsNegPions_Pt","ESD_Mother_AnglePosPionsNegPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
        fHistoAnglePiPlPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAnglePiPlPiMi[iCut]->GetYaxis()->SetTitle("#angle (#pi^{+},#pi^{-})");
        fESDList[iCut]->Add(fHistoAnglePiPlPiMi[iCut]);
        fHistoAnglePiZeroPiMi[iCut]         = new TH2F("ESD_Mother_AngleNeutralPionNegPions_Pt","ESD_Mother_AngleNeutralPionNegPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],360,0,TMath::Pi());
        fHistoAnglePiZeroPiMi[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAnglePiZeroPiMi[iCut]->GetYaxis()->SetTitle("#angle (#pi^{0},#pi^{-})");
        fESDList[iCut]->Add(fHistoAnglePiZeroPiMi[iCut]);
        fHistoAngleSum[iCut]                = new TH2F("ESD_Mother_AngleSum_Pt","ESD_Mother_AngleSum_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],720,0,2*TMath::Pi());
        fHistoAngleSum[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        fHistoAngleSum[iCut]->GetYaxis()->SetTitle("#sum #angle");
        fESDList[iCut]->Add(fHistoAngleSum[iCut]);
    }
    if ( fDoMesonQA>0 && (!fDoLightOutput) ) {
      TAxis *AxisAfter        = fHistoPionTPCdEdxNSigma[iCut]->GetXaxis();
      Int_t bins              = AxisAfter->GetNbins();
      Double_t from           = AxisAfter->GetXmin();
      Double_t to             = AxisAfter->GetXmax();
      Double_t *newBins       = new Double_t[bins+1];
      newBins[0]              = from;
      Double_t factor         = TMath::Power(to/from, 1./bins);
      for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];

      AxisAfter->Set(bins, newBins);
      AxisAfter               = fHistoPionTPCdEdx[iCut]->GetXaxis();
      AxisAfter->Set(bins, newBins);
      delete [] newBins;
    }

    fCutFolder[iCut]->Add(fESDList[iCut]);

  }

  if( fIsMC ){
    // MC Histogramms
    fMCList                                 = new TList*[fnCuts];
    // True Histogramms
    fTrueList                               = new TList*[fnCuts];
    if(!fDoLightOutput){
        if (fNeutralPionMode < 2){
            fHistoTrueConvGammaPt                 = new TH1F*[fnCuts];
            fHistoDoubleCountTrueConvGammaRPt     = new TH2F*[fnCuts];
            fHistoTrueConvGammaFromNeutralMesonPt = new TH1F*[fnCuts];
        }
        if (fNeutralPionMode > 0){
            fHistoTrueClusterGammaPt                  = new TH1F*[fnCuts];
            fHistoTrueClusterGammaFromNeutralMesonPt  = new TH1F*[fnCuts];
        }
        fHistoTruePosPionPt                     = new TH1F*[fnCuts];
        fHistoTrueNegPionPt                     = new TH1F*[fnCuts];
        fHistoTruePosPionFromNeutralMesonPt     = new TH1F*[fnCuts];
        fHistoTrueNegPionFromNeutralMesonPt     = new TH1F*[fnCuts];


        fHistoMCAllGammaPt                      = new TH1F*[fnCuts];
        if (fNeutralPionMode < 2){
            fHistoMCConvGammaPt                   = new TH1F*[fnCuts];
        }
        fHistoMCAllPosPionsPt                   = new TH1F*[fnCuts];
        fHistoMCAllNegPionsPt                   = new TH1F*[fnCuts];
        fHistoMCGammaFromNeutralMesonPt         = new TH1F*[fnCuts];
        fHistoMCPosPionsFromNeutralMesonPt      = new TH1F*[fnCuts];
        fHistoMCNegPionsFromNeutralMesonPt      = new TH1F*[fnCuts];
    }
    fHistoMCEtaPiPlPiMiPiZeroPt                     = new TH1F*[fnCuts];
    fHistoMCEtaPiPlPiMiPiZeroInAccPt                = new TH1F*[fnCuts];
    fHistoMCOmegaPiPlPiMiPiZeroPt                   = new TH1F*[fnCuts];
    fHistoMCOmegaPiPlPiMiPiZeroInAccPt              = new TH1F*[fnCuts];
    if(!fDoLightOutput){
        fHistoDoubleCountTruePi0InvMassPt               = new TH2F*[fnCuts];
        fHistoDoubleCountTrueEtaInvMassPt               = new TH2F*[fnCuts];
        fHistoDoubleCountTrueOmegaInvMassPt             = new TH2F*[fnCuts];
    }
    fHistoTrueMotherPiPlPiMiPiZeroInvMassPt         = new TH2F*[fnCuts];
    fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt    = new TH2F*[fnCuts];
    fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt      = new TH2F*[fnCuts];
    if(!fDoLightOutput){
        fHistoTrueMotherGammaGammaInvMassPt             = new TH2F*[fnCuts];
        fHistoTrueMotherGammaGammaFromEtaInvMassPt      = new TH2F*[fnCuts];
        fHistoTrueMotherGammaGammaFromOmegaInvMassPt    = new TH2F*[fnCuts];
        fHistoTrueAngleSum                              = new TH2F*[fnCuts];
    }
    if(!fDoLightOutput){
        if (fDoMesonQA>0){
            fHistoTruePionPionInvMassPt                               = new TH2F*[fnCuts];
            fHistoTruePionPionFromSameMotherInvMassPt                 = new TH2F*[fnCuts];
            fHistoTruePionPionFromEtaInvMassPt                        = new TH2F*[fnCuts];
            fHistoTruePionPionFromOmegaInvMassPt                      = new TH2F*[fnCuts];

            fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt              = new TH2F*[fnCuts];
            fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt            = new TH2F*[fnCuts];
            fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt              = new TH2F*[fnCuts];
            fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt         = new TH2F*[fnCuts];
            fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt              = new TH2F*[fnCuts];
            fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt              = new TH2F*[fnCuts];
            fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt            = new TH2F*[fnCuts];
            fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt          = new TH2F*[fnCuts];
            fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt            = new TH2F*[fnCuts];
            fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt            = new TH2F*[fnCuts];
            fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt            = new TH2F*[fnCuts];
            fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt          = new TH2F*[fnCuts];
            fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt            = new TH2F*[fnCuts];
            fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt            = new TH2F*[fnCuts];
            fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt       = new TH2F*[fnCuts];
            fHistoTruePiPlPiMiPiZeroContaminationInvMassPt            = new TH2F*[fnCuts];
            if (fDoMesonQA>1){
                fTrueTreeList                                           = new TList*[fnCuts];
                fTreePiPiSameMother                                     = new TTree*[fnCuts];
                fTreePiPiPiSameMother                                   = new TTree*[fnCuts];
                fTreeEventInfoOmega                                     = new TTree*[fnCuts];
                fTreeEventInfoEta                                       = new TTree*[fnCuts];
            }
        }
    }

    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent            = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringPion             = ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
      TString cutstringConvGamma        = "";
      if (fNeutralPionMode < 2)
        cutstringConvGamma              = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
      TString cutstringCaloGamma        = "";
      if (fNeutralPionMode > 0)
        cutstringCaloGamma              = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringNeutralPion      = ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson            = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      TString fullCutString             = "";
      if (fNeutralPionMode == 0)
        fullCutString                   = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),
                                               cutstringMeson.Data());
      else if (fNeutralPionMode == 1)
        fullCutString                   = Form("%i_%s_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(), cutstringNeutralPion.Data(),
                                               cutstringPion.Data(), cutstringMeson.Data());
      else if (fNeutralPionMode == 2)
        fullCutString                   = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(), cutstringCaloGamma.Data(), cutstringNeutralPion.Data(), cutstringPion.Data(),
                                               cutstringMeson.Data());
      TString nameMCList                = Form("%s MC histograms", fullCutString.Data());
      TString nameTrueRecList           = Form("%s True histograms", fullCutString.Data());
      TString nameTrueRecTTreeList      = Form("%s True TTrees", fullCutString.Data());

      fMCList[iCut]                     = new TList();
      fMCList[iCut]->SetName(nameMCList.Data());
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);

      if(!fDoLightOutput){
          fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoMCAllGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCAllGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma}");
          fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
          if (fNeutralPionMode < 2){
              fHistoMCConvGammaPt[iCut]       = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
              fHistoMCConvGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
              fHistoMCConvGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
              fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);
          }

          fHistoMCAllPosPionsPt[iCut]               = new TH1F("MC_AllPosPions_Pt","MC_AllPosPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoMCAllPosPionsPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCAllPosPionsPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fMCList[iCut]->Add(fHistoMCAllPosPionsPt[iCut]);
          fHistoMCAllNegPionsPt[iCut]               = new TH1F("MC_AllNegPions_Pt","MC_AllNegPions_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoMCAllNegPionsPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCAllNegPionsPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fMCList[iCut]->Add(fHistoMCAllNegPionsPt[iCut]);
          fHistoMCGammaFromNeutralMesonPt[iCut]     = new TH1F("MC_GammaFromNeutralMeson_Pt","MC_GammaFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoMCGammaFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCGammaFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#gamma}");
          fMCList[iCut]->Add(fHistoMCGammaFromNeutralMesonPt[iCut]);
          fHistoMCPosPionsFromNeutralMesonPt[iCut]  = new TH1F("MC_PosPionsFromNeutralMeson_Pt","MC_PosPionsFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoMCPosPionsFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCPosPionsFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fMCList[iCut]->Add(fHistoMCPosPionsFromNeutralMesonPt[iCut]);
          fHistoMCNegPionsFromNeutralMesonPt[iCut]  = new TH1F("MC_NegPionsFromNeutralMeson_Pt","MC_NegPionsFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoMCNegPionsFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoMCNegPionsFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fMCList[iCut]->Add(fHistoMCNegPionsFromNeutralMesonPt[iCut]);
      }
      fHistoMCEtaPiPlPiMiPiZeroPt[iCut]         = new TH1F("MC_Eta_Pt","MC_Eta_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoMCEtaPiPlPiMiPiZeroPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoMCEtaPiPlPiMiPiZeroPt[iCut]->GetYaxis()->SetTitle("N_{#eta}");
      fHistoMCEtaPiPlPiMiPiZeroPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCEtaPiPlPiMiPiZeroPt[iCut]);

      fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut]    = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut]->GetYaxis()->SetTitle("A #times N_{#eta}");
      fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut]);

      fHistoMCOmegaPiPlPiMiPiZeroPt[iCut]       = new TH1F("MC_Omega_Pt","MC_Omega_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoMCOmegaPiPlPiMiPiZeroPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoMCOmegaPiPlPiMiPiZeroPt[iCut]->GetYaxis()->SetTitle("N_{#omega}");
      fHistoMCOmegaPiPlPiMiPiZeroPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCOmegaPiPlPiMiPiZeroPt[iCut]);

      fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut]  = new TH1F("MC_OmegaInAcc_Pt","MC_OmegaInAcc_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut]->GetYaxis()->SetTitle("A #times N_{#omega}");
      fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut]);

      fTrueList[iCut]                           = new TList();
      fTrueList[iCut]->SetName(nameTrueRecList.Data());
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      if(!fDoLightOutput){
          if (fNeutralPionMode < 2){
              fHistoTrueConvGammaPt[iCut]                 = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
              fHistoTrueConvGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
              fHistoTrueConvGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
              fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);
              fHistoDoubleCountTrueConvGammaRPt[iCut]     = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt","ESD_TrueDoubleCountConvGamma_R_Pt",800,0,200,300,0,30);
              fHistoDoubleCountTrueConvGammaRPt[iCut]->GetXaxis()->SetTitle("R_{conv} (cm)");
              fHistoDoubleCountTrueConvGammaRPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);
              fHistoTrueConvGammaFromNeutralMesonPt[iCut] = new TH1F("ESD_TrueConvGammaFromNeutralMeson_Pt","ESD_TrueConvGammaFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
              fHistoTrueConvGammaFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
              fHistoTrueConvGammaFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,conv}");
              fTrueList[iCut]->Add(fHistoTrueConvGammaFromNeutralMesonPt[iCut]);
          }
          if (fNeutralPionMode > 0){
              fHistoTrueClusterGammaPt[iCut]                  = new TH1F("ESD_TrueClusterGamma_Pt","ESD_TrueClusterGamma_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
              fHistoTrueClusterGammaPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
              fHistoTrueClusterGammaPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
              fTrueList[iCut]->Add(fHistoTrueClusterGammaPt[iCut]);
              fHistoTrueClusterGammaFromNeutralMesonPt[iCut]  = new TH1F("ESD_TrueClusterGammaFromNeutralMeson_Pt","ESD_TrueClusterGammaFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
              fHistoTrueClusterGammaFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
              fHistoTrueClusterGammaFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#gamma,cluster}");
              fTrueList[iCut]->Add(fHistoTrueClusterGammaFromNeutralMesonPt[iCut]);
          }
          fHistoTruePosPionPt[iCut]                       = new TH1F("ESD_TruePosPion_Pt","ESD_TruePosPion_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePosPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTruePosPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fTrueList[iCut]->Add(fHistoTruePosPionPt[iCut]);
          fHistoTrueNegPionPt[iCut]                       = new TH1F("ESD_TrueNegPion_Pt","ESD_TrueNegPion_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueNegPionPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTrueNegPionPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fTrueList[iCut]->Add(fHistoTrueNegPionPt[iCut]);

          fHistoTrueNegPionFromNeutralMesonPt[iCut]       = new TH1F("ESD_TrueNegPionFromNeutralMeson_Pt","ESD_TrueNegPionFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueNegPionFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTrueNegPionFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{-}}");
          fTrueList[iCut]->Add(fHistoTrueNegPionFromNeutralMesonPt[iCut]);
          fHistoTruePosPionFromNeutralMesonPt[iCut]       = new TH1F("ESD_TruePosPionFromNeutralMeson_Pt","ESD_TruePosPionFromNeutralMeson_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTruePosPionFromNeutralMesonPt[iCut]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTruePosPionFromNeutralMesonPt[iCut]->GetYaxis()->SetTitle("N_{#pi^{+}}");
          fTrueList[iCut]->Add(fHistoTruePosPionFromNeutralMesonPt[iCut]);

          fHistoDoubleCountTruePi0InvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt","ESD_TrueDoubleCountPi0_InvMass_Pt",800,0,0.8,300,0,30);
          fHistoDoubleCountTruePi0InvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
          fHistoDoubleCountTruePi0InvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
          fHistoDoubleCountTrueEtaInvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt","ESD_TrueDoubleCountEta_InvMass_Pt",800,0,0.8,300,0,30);
          fHistoDoubleCountTrueEtaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#eta} (GeV/c^{2})");
          fHistoDoubleCountTrueEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);
          fHistoDoubleCountTrueOmegaInvMassPt[iCut]       = new TH2F("ESD_TrueDoubleCountOmega_InvMass_Pt","ESD_TrueDoubleCountOmega_InvMass_Pt",800,0,0.8,300,0,30);
          fHistoDoubleCountTrueOmegaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#omega} (GeV/c^{2})");
          fHistoDoubleCountTrueOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoDoubleCountTrueOmegaInvMassPt[iCut]);
      }
      fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut]       = new TH2F("ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt","ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
      fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut]->Sumw2();
      fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
      fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut]);

          fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt[iCut]  = new TH2F("ESD_TrueMotherOmegaPiPlPiMiPiZero_InvMass_Pt","ESD_TrueMotherOmegaPiPlPiMiPiZero_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt[iCut]    = new TH2F("ESD_TrueMotherEtaPiPlPiMiPiZero_InvMass_Pt","ESD_TrueMotherEtaPiPlPiMiPiZero_InvMass_Pt",HistoNMassBins,HistoMassRange[0],HistoMassRange[1],HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt[iCut]->Sumw2();
          fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt[iCut]->Sumw2();
          fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
          fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+} #pi^{-} #pi^{0}} (GeV/c^{2})");
          fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt[iCut]);
          fTrueList[iCut]->Add(fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt[iCut]);

      if(!fDoLightOutput){
          fHistoTrueMotherGammaGammaInvMassPt[iCut]           = new TH2F("ESD_TrueMotherGG_InvMass_Pt","ESD_TrueMotherGG_InvMass_Pt",450,0.,0.45,HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueMotherGammaGammaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
          fHistoTrueMotherGammaGammaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaInvMassPt[iCut]);
          fHistoTrueMotherGammaGammaFromEtaInvMassPt[iCut]    = new TH2F("ESD_TrueMotherGGFromEta_InvMass_Pt","ESD_TrueMotherGGFromEta_InvMass_Pt",450,0.,0.45,HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueMotherGammaGammaFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
          fHistoTrueMotherGammaGammaFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaFromEtaInvMassPt[iCut]);
          fHistoTrueMotherGammaGammaFromOmegaInvMassPt[iCut]  = new TH2F("ESD_TrueMotherGGFromOmega_InvMass_Pt","ESD_TrueMotherGGFromOmega_InvMass_Pt",450,0.,0.45,HistoNPtBins,HistoPtRange[0],HistoPtRange[1]);
          fHistoTrueMotherGammaGammaFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV/c^{2})");
          fHistoTrueMotherGammaGammaFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaFromOmegaInvMassPt[iCut]);
          fHistoTrueAngleSum[iCut]                            = new TH2F("ESD_TrueMother_AngleSum_Pt","ESD_TrueMother_AngleSum_Pt",HistoNPtBins,HistoPtRange[0],HistoPtRange[1],720,0,2*TMath::Pi());
          fHistoTrueAngleSum[iCut]->GetXaxis()->SetTitle("#sum #angle");
          fHistoTrueAngleSum[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
          fTrueList[iCut]->Add(fHistoTrueAngleSum[iCut]);

          if (fDoMesonQA>0){
              fHistoTruePionPionInvMassPt[iCut]                 = new TH2F("ESD_TruePiPlusPiNeg_InvMassPt","ESD_TruePiPlusPiNeg_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePionPionInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePionPionInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePionPionInvMassPt[iCut]);
              fHistoTruePionPionFromSameMotherInvMassPt[iCut]   = new TH2F("ESD_TruePiPlusPiNegFromSameMother_InvMassPt","ESD_TruePiPlusPiNegFromSameMother_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePionPionFromSameMotherInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePionPionFromSameMotherInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePionPionFromSameMotherInvMassPt[iCut]);
              fHistoTruePionPionFromEtaInvMassPt[iCut]          = new TH2F("ESD_TruePiPlusPiNegFromEta_InvMassPt","ESD_TruePiPlusPiNegFromEta_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePionPionFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePionPionFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePionPionFromEtaInvMassPt[iCut]);
              fHistoTruePionPionFromOmegaInvMassPt[iCut]        = new TH2F("ESD_TruePiPlusPiNegFromOmega_InvMassPt","ESD_TruePiPlusPiNegFromOmega_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePionPionFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePionPionFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePionPionFromOmegaInvMassPt[iCut]);

              fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromEta_InvMassPt","ESD_TruePiPlPiMiSameMotherFromEta_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[iCut]);
              fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiSameMotherFromOmega_InvMassPt","ESD_TruePiPlPiMiSameMotherFromOmega_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[iCut]);
              fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromRho_InvMassPt","ESD_TruePiPlPiMiSameMotherFromRho_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[iCut]);
              fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut] = new TH2F("ESD_TruePiPlPiMiSameMotherFromEtaPrime_InvMassPt","ESD_TruePiPlPiMiSameMotherFromEtaPrime_InvMassPt",
                                                                                 2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[iCut]);
              fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromK0s_InvMassPt","ESD_TruePiPlPiMiSameMotherFromK0s_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[iCut]);
              fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]    = new TH2F("ESD_TruePiPlPiMiSameMotherFromK0l_InvMassPt","ESD_TruePiPlPiMiSameMotherFromK0l_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} (GeV/c^{2})");
              fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[iCut]);

              fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]  = new TH2F("ESD_TruePiMiPiZeroSameMotherFromEta_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromEta_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{-}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[iCut]);
              fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut] = new TH2F("ESD_TruePiMiPiZeroSameMotherFromOmega_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromOmega_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{-}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[iCut]);
              fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]  = new TH2F("ESD_TruePiMiPiZeroSameMotherFromRho_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromRho_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{-}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]);
              fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]  = new TH2F("ESD_TruePiMiPiZeroSameMotherFromK0l_InvMassPt","ESD_TruePiMiPiZeroSameMotherFromK0l_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{-}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[iCut]);

              fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiZeroSameMotherFromEta_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromEta_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[iCut]);
              fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut] = new TH2F("ESD_TruePiPlPiZeroSameMotherFromOmega_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromOmega_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[iCut]);
              fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiZeroSameMotherFromRho_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromRho_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[iCut]);
              fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiZeroSameMotherFromK0l_InvMassPt","ESD_TruePiPlPiZeroSameMotherFromK0l_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[iCut]);
              fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiPiZeroPureCombinatorical_InvMassPt","ESD_TruePiPlPiMiPiZeroPureCombinatorical_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt[iCut]);
              fHistoTruePiPlPiMiPiZeroContaminationInvMassPt[iCut]  = new TH2F("ESD_TruePiPlPiMiPiZeroContamination_InvMassPt","ESD_TruePiPlPiMiPiZeroContamination_InvMassPt",2000,0.,2.,200,0.,20.);
              fHistoTruePiPlPiMiPiZeroContaminationInvMassPt[iCut]->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}#pi^{0}} (GeV/c^{2})");
              fHistoTruePiPlPiMiPiZeroContaminationInvMassPt[iCut]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
              fTrueList[iCut]->Add(fHistoTruePiPlPiMiPiZeroContaminationInvMassPt[iCut]);
              if(fDoMesonQA>1){
                  fTrueTreeList[iCut]                               = new TList();
                  fTrueTreeList[iCut]->SetName(nameTrueRecTTreeList.Data());
                  fTrueTreeList[iCut]->SetOwner(kTRUE);
                  fCutFolder[iCut]->Add(fTrueTreeList[iCut]);

                  fTreePiPiSameMother[iCut]                         = new TTree("TreePiPiSameMother","TreePiPiSameMother");
                  fTreePiPiSameMother[iCut]->Branch("fCasePiPi", &fCasePiPi, "fCasePiPi/S");
                  fTreePiPiSameMother[iCut]->Branch("fSamePiPiMotherID", &fSamePiPiMotherID, "fSamePiPiMotherID/F");
                  fTreePiPiSameMother[iCut]->Branch("fSamePiPiMotherInvMass", &fSamePiPiMotherInvMass, "fSamePiPiMotherInvMass/F");
                  fTreePiPiSameMother[iCut]->Branch("fSamePiPiMotherPt", &fSamePiPiMotherPt, "fSamePiPiMotherPt/F");
                  fTrueTreeList[iCut]->Add(fTreePiPiSameMother[iCut]);

                  fTreePiPiPiSameMother[iCut]                         = new TTree("TreePiPiPiSameMother","TreePiPiPiSameMother");
                  fTreePiPiPiSameMother[iCut]->Branch("fSamePiPiPiMotherID", &fSamePiPiPiMotherID, "fSamePiPiPiMotherID/F");
                  fTreePiPiPiSameMother[iCut]->Branch("fSamePiPiPiMotherInvMass", &fSamePiPiPiMotherInvMass, "fSamePiPiPiMotherInvMass/F");
                  fTreePiPiPiSameMother[iCut]->Branch("fSamePiPiPiMotherPt", &fSamePiPiPiMotherPt, "fSamePiPiPiMotherPt/F");
                  fTrueTreeList[iCut]->Add(fTreePiPiPiSameMother[iCut]);

                  fTreeEventInfoOmega[iCut]                         = new TTree("TreeEventInfoOmega","TreeEventInfoOmega");
                  fTreeEventInfoOmega[iCut]->Branch("fV0MultiplicityOmegaEvent", &fV0MultiplicityOmegaEvent, "fV0MultiplicityOmegaEvent/F");
                  fTreeEventInfoOmega[iCut]->Branch("fTrackMultiplicityOmegaEvent", &fTrackMultiplicityOmegaEvent, "fTrackMultiplicityOmegaEvent/F");
                  fTreeEventInfoOmega[iCut]->Branch("fZVertexOmegaEvent", &fZVertexOmegaEvent, "fZVertexOmegaEvent/F");
                  fTreeEventInfoOmega[iCut]->Branch("fPtOmega", &fPtOmega, "fPtOmega/F");
                  fTrueTreeList[iCut]->Add(fTreeEventInfoOmega[iCut]);

                  fTreeEventInfoEta[iCut]                         = new TTree("TreeEventInfoEta","TreeEventInfoEta");
                  fTreeEventInfoEta[iCut]->Branch("fV0MultiplicityEtaEvent", &fV0MultiplicityEtaEvent, "fV0MultiplicityEtaEvent/F");
                  fTreeEventInfoEta[iCut]->Branch("fTrackMultiplicityEtaEvent", &fTrackMultiplicityEtaEvent, "fTrackMultiplicityEtaEvent/F");
                  fTreeEventInfoEta[iCut]->Branch("fZVertexEtaEvent", &fZVertexEtaEvent, "fZVertexEtaEvent/F");
                  fTreeEventInfoEta[iCut]->Branch("fPtEta", &fPtEta, "fPtEta/F");
                  fTrueTreeList[iCut]->Add(fTreeEventInfoEta[iCut]);
              }
          }
      }
    }
  }

  fVectorDoubleCountTruePi0s.clear();
  fVectorDoubleCountTrueEtas.clear();
  fVectorDoubleCountTrueOmegas.clear();
  fVectorDoubleCountTrueConvGammas.clear();

  InitBack(); // Init Background Handler

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  if(fV0Reader){
    if((AliConvEventCuts*)fV0Reader->GetEventCuts()){
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms()){
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
      }
    }

    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts()){
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms()){
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
      }
    }

  }

  for(Int_t iMatcherTask = 0; iMatcherTask < 3; iMatcherTask++){
    AliCaloTrackMatcher* temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }

  fPionSelector=(AliPrimaryPionSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("PionSelector");
  if(!fPionSelector){printf("Error: No PionSelector");return;} // GetV0Reader

  if( fPionSelector && (!fDoLightOutput)){
    if ( ((AliPrimaryPionCuts*)fPionSelector->GetPrimaryPionCuts())->GetCutHistograms() ){
      fOutputContainer->Add( ((AliPrimaryPionCuts*)fPionSelector->GetPrimaryPionCuts())->GetCutHistograms() );
    }
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if( fEventCutArray) {
      if( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
      }
    }

    if( fPionCutArray){
      if( ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutHistograms() );
      }
    }
    if (fNeutralPionMode < 2){
      if( fGammaCutArray ) {
        if( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms() ) {
          fCutFolder[iCut]->Add( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms()  );
        }
      }
    }
    if (fNeutralPionMode > 0){
      if( fClusterCutArray ) {
        if( ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms() ) {
          fCutFolder[iCut]->Add( ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()  );
        }
      }
    }
    if( fNeutralPionMesonCutArray ) {
      if( ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutHistograms());
      }
    }
    if( fMesonCutArray ) {
      if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
      }
    }
  }

  PostData(1, fOutputContainer);

}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::UserExec(Option_t *){

  //
  // Execute analysis for current event
  //

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(InputEvent()->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
    }
    return;
  }

  fPionSelector=(AliPrimaryPionSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("PionSelector");
  if(!fPionSelector){printf("Error: No PionSelector");return;} // GetV0Reader

  if(fIsMC) fMCEvent     =  MCEvent();
  fESDEvent        = (AliESDEvent*)InputEvent();
  fReaderGammas    = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
  fSelectorNegPionIndex = fPionSelector->GetReconstructedNegPionIndex(); // Electrons from default Cut
  fSelectorPosPionIndex = fPionSelector->GetReconstructedPosPionIndex(); // Positrons from default Cut

  fNumberOfESDTracks = fV0Reader->GetNumberOfPrimaryTracks();
  //AddTaskContainers(); //Add conatiner

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    fiCut = iCut;

    Bool_t isRunningEMCALrelAna = kFALSE;
    if (fNeutralPionMode > 0){
      if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
    }

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

    if(eventNotAccepted){
      // 			cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      fHistoNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      continue;
    }

    if(eventQuality != 0){// Event Not Accepted
      // 			cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality);
      continue;
    }

    fHistoNEvents[iCut]->Fill(eventQuality);
    fHistoNGoodESDTracks[iCut]->Fill(fNumberOfESDTracks);
    if(!fDoLightOutput){
        fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)));
    }
    if(fMCEvent){ // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                                            ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                                            fMCEvent);
      }
      ProcessMCParticles();
    }

    if (fNeutralPionMode < 2){
      ProcessConversionPhotonCandidates(); // Process this cuts conversion gammas
    }
    if (fNeutralPionMode > 0){
      ProcessCaloPhotonCandidates(); // Process this cuts calo gammas
    }

    if (fNeutralPionMode == 0 ){
      ProcessNeutralPionCandidatesPureConversions(); // Process neutral pion candidates purely from conversions
    }
    if (fNeutralPionMode == 1){
      ProcessNeutralPionCandidatesMixedConvCalo(); // Process neutral pion candidates mixed conv and calo
    }
    if (fNeutralPionMode == 2){
      ProcessNeutralPionCandidatesPureCalo(); // Process neutral pion candidates purely from calo
    }

    ProcessPionCandidates(); // Process this cuts gammas

    CalculateMesonCandidates();
    CalculateBackground();
    UpdateEventByEventData();

    fVectorDoubleCountTruePi0s.clear();
    fVectorDoubleCountTrueEtas.clear();
    fVectorDoubleCountTrueOmegas.clear();
    fVectorDoubleCountTrueConvGammas.clear();

    fGoodConvGammas->Clear();
    fClusterCandidates->Clear();
    fNeutralPionCandidates->Clear();
    if(fNeutralPionSidebandCandidates) fNeutralPionSidebandCandidates->Clear();
    fPosPionCandidates->Clear();
    fNegPionCandidates->Clear();
    fGoodVirtualParticles->Clear(); // delete this cuts good gammas
  }

  fSelectorNegPionIndex.clear();
  fSelectorPosPionIndex.clear();

  PostData( 1, fOutputContainer );
}
//________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::Notify(){
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }

    if( !((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift() ){
      if(!fDoLightOutput){
        fProfileEtaShift[iCut]->Fill(0.,0.);
      }
      continue; // No Eta Shift requested, continue
    }
    if( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() );
      if(!fDoLightOutput){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      }
      continue;
    } else {
      printf(" Eta t PiPlusPiMinus Gamma Task %s :: Eta Shift Manually Set to %f \n\n",
      (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() );
      if(!fDoLightOutput){
        fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      }
    }
  }
  return kTRUE;
}


void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::Terminate(const Option_t *){
///Grid
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessCaloPhotonCandidates()
{

  Int_t nclus = 0;
  nclus = fInputEvent->GetNumberOfCaloClusters();

// 	cout << nclus << endl;

  if(nclus == 0)	return;

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Long_t i = 0; i < nclus; i++){

    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

    if (!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,1.,i)){ delete clus; continue;}
    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef(i);
    // get MC label
    if(fIsMC){
      Int_t* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
// 			cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
          if (k< 50)PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
// 					Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
// 					cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
        }
      }
    }

    fIsFromMBHeader = kTRUE;
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0) fIsFromMBHeader = kFALSE;

    if (fIsFromMBHeader && (!fDoLightOutput)){
      fHistoClusterGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
      fHistoClusterGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
    }
    fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

    if(fIsMC){
// 			if(fInputEvent->IsA()==AliESDEvent::Class()){
        ProcessTrueCaloPhotonCandidates(PhotonCandidate);
// 			} else {
// 				ProcessTrueClusterCandidatesAOD(PhotonCandidate);
// 			}
    }

    delete clus;
    delete tmpvec;
  }

}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueCaloPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  TParticle *Photon = NULL;
  if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TruePhotonCandidate->GetCaloPhotonMCLabel(0)<0) return;
// 	fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;

  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }

// 	Int_t pdgCodeParticle = Photon->GetPdgCode();
  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, kFALSE);

  // True Photon
  if(fIsFromMBHeader && (!fDoLightOutput)){
    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(isPrimary){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
        if (GammaIsNeutralMesonPiPlPiMiPiZeroDaughter(TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
          fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt());
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
          fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
          if (GammaIsNeutralMesonPiPlPiMiPiZeroDaughter(TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
          fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt());
        }
      }
    }
  }
  return;
}



//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessConversionPhotonCandidates(){
  Int_t nV0 = 0;
  TList *GoodGammasStepOne = new TList();
  TList *GoodGammasStepTwo = new TList();
  // Loop over Photon Candidates allocated by ReaderV1

  for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;

    fIsFromMBHeader = kTRUE;

    if( fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0 ){
      Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
      if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
      if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
    }

    if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

    if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
      !((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas

      fGoodConvGammas->Add(PhotonCandidate);

      if(fIsFromMBHeader && (!fDoLightOutput)){
        fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
        fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
      }

      if(fMCEvent){
        ProcessTrueConversionPhotonCandidates(PhotonCandidate);
      }
    } else if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
      ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
      nV0++;
      GoodGammasStepOne->Add(PhotonCandidate);
    } else if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
        ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
      GoodGammasStepTwo->Add(PhotonCandidate);
    }
  }


  if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){
    for(Int_t i = 0;i<GoodGammasStepOne->GetEntries();i++){
      AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GoodGammasStepOne->At(i);
      if(!PhotonCandidate) continue;
      fIsFromMBHeader = kTRUE;
      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
        Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }
      if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
      if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
        fGoodConvGammas->Add(PhotonCandidate);
        if(fIsFromMBHeader && (!fDoLightOutput)){
          fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
          fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
        }
        if(fMCEvent){
          ProcessTrueConversionPhotonCandidates(PhotonCandidate);
        }
      }
      else GoodGammasStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
    }
  }
  if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GoodGammasStepTwo->GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GoodGammasStepTwo->At(i);
      if(!PhotonCandidate) continue;

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
        Int_t isNegFromMBHeader
        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
        if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }

      if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
      fGoodConvGammas->Add(PhotonCandidate); // Add gamma to current cut TList

      if(fIsFromMBHeader && (!fDoLightOutput)){
        fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt()); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
        fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
      }

      if(fMCEvent){
        ProcessTrueConversionPhotonCandidates(PhotonCandidate);
      }
    }
  }

  delete GoodGammasStepOne;
  GoodGammasStepOne = 0x0;
  delete GoodGammasStepTwo;
  GoodGammasStepTwo = 0x0;
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueConversionPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
  // Process True Photons
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();


  if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
  if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){  // Not Same Mother == Combinatorial Bck
    return;
  }

  else if (posDaughter->GetMother(0) == -1){
    return;
  }

  if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11) return; //One Particle is not electron
  if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge
  if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

  TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);
  if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

  // True Photon

  if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother(0)) && (!fDoLightOutput)) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt());

  Int_t labelGamma = TruePhotonCandidate->GetMCParticleLabel(fMCEvent);
  Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelGamma, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
  if( gammaIsPrimary ){
    if( fIsFromMBHeader && (!fDoLightOutput) ){
      fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
      if (GammaIsNeutralMesonPiPlPiMiPiZeroDaughter(labelGamma)){
        fHistoTrueConvGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt());
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessNeutralPionCandidatesPureConversions(){
  // Conversion Gammas
  if(fGoodConvGammas->GetEntries()>1){
    for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodConvGammas->GetEntries()-1;firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGoodConvGammas->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(secondGammaIndex));
        //Check for same Electron ID
        if (gamma1==NULL) continue;
        if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
        gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
        gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
        gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;

        AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
        pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

        pi0cand->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

        if(!fDoLightOutput){
            fHistoGammaGammaInvMassPtBeforeCuts[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
        }
        if((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if(fIsMC){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueNeutralPionCandidatesPureConversions(pi0cand,gamma0,gamma1);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueNeutralPionCandidatesPureConversionsAOD(pi0cand,gamma0,gamma1);
          }
          if (((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 0)){
            fNeutralPionCandidates->Add(pi0cand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) &&
                    (((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 1))){
            fNeutralPionSidebandCandidates->Add(pi0cand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides()) &&
                    ((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 2)) ||
                     ((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 3))))){
            fNeutralPionSidebandCandidates->Add(pi0cand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
            }
          } else{
            delete pi0cand;
            pi0cand=0x0;
          }
        }else{
          delete pi0cand;
          pi0cand=0x0;
        }
      }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessNeutralPionCandidatesPureCalo(){

  // Conversion Gammas
  if(fClusterCandidates->GetEntries()>0){

    // vertex
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;

      for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
        if (firstGammaIndex == secondGammaIndex) continue;
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
        if (gamma1==NULL) continue;

        AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
        pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

        if(!fDoLightOutput){
            fHistoGammaGammaInvMassPtBeforeCuts[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
        }
        if((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if(fIsMC){
            ProcessTrueNeutralPionCandidatesPureCalo(pi0cand,gamma0,gamma1);
          }

          if (((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 0)){
            fNeutralPionCandidates->Add(pi0cand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) &&
                    (((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 1))){
            fNeutralPionSidebandCandidates->Add(pi0cand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
            }
          } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides()) &&
                    ((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 2)) ||
                      ((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 3))))){
            fNeutralPionSidebandCandidates->Add(pi0cand);

            if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
            }
          }else {
            delete pi0cand;
            pi0cand=0x0;
          }
        } else{
          delete pi0cand;
          pi0cand=0x0;
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesPureCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons

  Bool_t isTruePi0 = kFALSE;
  Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); 	// get most probable MC label
  Int_t gamma0MotherLabel = -1;
  Int_t motherRealLabel = -1;

  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){		// largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){														// for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother(0);
        motherRealLabel=gammaMC0->GetFirstMother();
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
                if (TrueGammaCandidate0->IsConversion() && gammaMC0->GetMother(0)>-1){
          gamma0MotherLabel=fMCEvent->Particle(gammaMC0->GetMother(0))->GetMother(0);
          motherRealLabel=fMCEvent->Particle(gammaMC0->GetMother(0))->GetMother(0);
        } else {
          gamma0MotherLabel=gammaMC0->GetMother(0);
          motherRealLabel=gammaMC0->GetMother(0);
        }
      }
    }
  }

  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

  Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
  Int_t gamma1MotherLabel = -1;
  // check if
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother(0);
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
                if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
        else gamma1MotherLabel=gammaMC1->GetMother(0);
      }
    }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
      isTruePi0=kTRUE;
      if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
    }
  }

  if(isTruePi0){// True Pion
    Pi0Candidate->SetTrueMesonValue(1);
    Pi0Candidate->SetMCLabel(motherRealLabel);
    if(!fDoLightOutput){
        fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
            fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
        if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
            fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
    }
  }
}



//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesPureConversions(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTruePi0 = kFALSE;
    Bool_t isTruePi0Dalitz = kFALSE;
    Bool_t gamma0DalitzCand = kFALSE;
    Bool_t gamma1DalitzCand = kFALSE;
    Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
    Int_t gamma0MotherLabel = -1;
    Int_t motherRealLabel = -1;
    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(fMCEvent);
      TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(fMCEvent);
      TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetFirstMother();
            motherRealLabel=gammaMC0->GetFirstMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
          motherRealLabel=gamma0MCLabel;
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
          if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-111;
          }
        }
      }
      if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
        if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
          isTruePi0=kTRUE;
          if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
        }
      }

      //Identify Dalitz candidate
      if (gamma1DalitzCand || gamma0DalitzCand){
        if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
          if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
        }
        if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
          if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
        }
      }


      if(isTruePi0 || isTruePi0Dalitz){// True Pion
        Pi0Candidate->SetTrueMesonValue(1);
        Pi0Candidate->SetMCLabel(motherRealLabel);
        if(!fDoLightOutput){
            fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
                fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
            if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
                fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesPureConversionsAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{

  // Process True Mesons
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  Bool_t isTruePi0 = kFALSE;
  Bool_t isTruePi0Dalitz = kFALSE;
  Bool_t gamma0DalitzCand = kFALSE;
  Bool_t gamma1DalitzCand = kFALSE;
  Int_t motherRealLabel = -1;

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
            motherRealLabel=gammaMC0->GetMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
          motherRealLabel=gamma0MCLabel;
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
        if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-111;
        }
      }
    }
    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) &&(!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }

    //Identify Dalitz candidate
    if (gamma1DalitzCand || gamma0DalitzCand){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
        if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
      }
      if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
        if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
      }
    }

    if(isTruePi0 || isTruePi0Dalitz){// True Pion
      Pi0Candidate->SetTrueMesonValue(1);
      Pi0Candidate->SetMCLabel(motherRealLabel);
      if(!fDoLightOutput){
          fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
              fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
          if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
              fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessNeutralPionCandidatesMixedConvCalo(){

  // Conversion Gammas
  if(fGoodConvGammas->GetEntries()>0){
    // vertex
    Double_t vertex[3] = {0};
    InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

    for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodConvGammas->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(firstGammaIndex));
      if (gamma0==NULL) continue;

      for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
        Bool_t matched = kFALSE;
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
        if (gamma1==NULL) continue;

        if (gamma1->GetIsCaloPhoton() > 0){
          AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
          matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent );
        }

        AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
        pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

        if(!fDoLightOutput){
          fHistoGammaGammaInvMassPtBeforeCuts[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
        }

        if((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
          if (!matched){
            if(fIsMC){
              ProcessTrueNeutralPionCandidatesMixedConvCalo(pi0cand,gamma0,gamma1);
            }
            if (((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 0)){
              fNeutralPionCandidates->Add(pi0cand);

              if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
              }
            } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) &&
                      (((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 1))){
              fNeutralPionSidebandCandidates->Add(pi0cand);

              if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
              }
            } else if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides()) &&
                      ((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 2)) ||
                      ((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelectedByMassCut(pi0cand, 3))))){
              fNeutralPionSidebandCandidates->Add(pi0cand);

              if(!fDoLightOutput){
                fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
              }
            } else{
              delete pi0cand;
              pi0cand=0x0;
            }
          }else{
            delete pi0cand;
            pi0cand=0x0;
          }
        }else{
          delete pi0cand;
          pi0cand=0x0;
        }
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesMixedConvCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
  // Process True Mesons
  if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
    Bool_t isTruePi0 = kFALSE;
    Bool_t isTruePi0Dalitz = kFALSE;
    Bool_t gamma0DalitzCand = kFALSE;

    Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(fMCEvent);
    Int_t gamma0MotherLabel = -1;
    Int_t motherRealLabel = -1;
    if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(fMCEvent);
      TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(fMCEvent);
      TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
      if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
        if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
          if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
            gamma0MotherLabel=gammaMC0->GetFirstMother();
            motherRealLabel=gammaMC0->GetFirstMother();
          }
        }
        if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
          gamma0DalitzCand = kTRUE;
          gamma0MotherLabel=-111;
          motherRealLabel=gamma0MCLabel;
        }

      }
    }

    if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

    Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
    Int_t gamma1MotherLabel = -1;
    // check if

    if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
      if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
        // get mother of interest (pi0 or eta)
        if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother
          gamma1MotherLabel=gammaMC1->GetMother(0);
        } else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
                    if (TrueGammaCandidate1->IsConversion() && gammaMC1->GetMother(0)>-1) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
          else gamma1MotherLabel=gammaMC1->GetMother(0);
        }
      }
    }

    if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
        isTruePi0=kTRUE;
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
    }

    if (gamma0DalitzCand ){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
        if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
      }
    }

    if(isTruePi0 || isTruePi0Dalitz ){
      Pi0Candidate->SetTrueMesonValue(1);
      Pi0Candidate->SetMCLabel(motherRealLabel);
      if(!fDoLightOutput){
          fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
              fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
          if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) {
              fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
          }
      }
    }
  }
}



//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessPionCandidates(){

  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }

  vector<Int_t> lGoodNegPionIndexPrev(0);
  vector<Int_t> lGoodPosPionIndexPrev(0);

    for(UInt_t i = 0; i < fSelectorNegPionIndex.size(); i++){
    AliESDtrack* negPionCandidate = fESDEvent->GetTrack(fSelectorNegPionIndex[i]);
    if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelected(negPionCandidate) ) continue;
    lGoodNegPionIndexPrev.push_back(   fSelectorNegPionIndex[i] );

    TLorentzVector *negPionforHandler = new TLorentzVector();
    negPionforHandler->SetPxPyPzE(negPionCandidate->Px(), negPionCandidate->Py(), negPionCandidate->Pz(), negPionCandidate->E());

    AliAODConversionPhoton *negPionHandler = new AliAODConversionPhoton(negPionforHandler);
    delete negPionforHandler;

    fNegPionCandidates->Add(negPionHandler);
    if(!fDoLightOutput){
        fHistoNegPionPt[fiCut]->Fill(negPionCandidate->Pt());
        fHistoNegPionPhi[fiCut]->Fill(negPionCandidate->Phi());
    }

    if( fMCEvent ) {
      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();

      Int_t labelNegPion = TMath::Abs( negPionCandidate->GetLabel() );
      Bool_t negPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            if( labelNegPion>-1 && labelNegPion < fMCEvent->GetNumberOfTracks() ){
        TParticle* negPion = fMCEvent->Particle(labelNegPion);
        if( negPion->GetPdgCode() ==  -211 ){
            if(!fDoLightOutput){

                if( negPionIsPrimary ){
                    fHistoTrueNegPionPt[fiCut]->Fill(negPionCandidate->Pt());    //primary negPion
                }
                if( IsEtaPiPlPiMiPiZeroDaughter(labelNegPion) || IsOmegaPiPlPiMiPiZeroDaughter(labelNegPion) ) {
                    if( negPionIsPrimary ) {
                        fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt());
                    }
                }
            }
        }
      }
    }
  }

  for(UInt_t i = 0; i < fSelectorPosPionIndex.size(); i++){
    AliESDtrack* posPionCandidate = fESDEvent->GetTrack( fSelectorPosPionIndex[i] );
    if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelected(posPionCandidate) ) continue;
    lGoodPosPionIndexPrev.push_back(   fSelectorPosPionIndex[i]  );

    TLorentzVector *posPionforHandler = new TLorentzVector();
    posPionforHandler->SetPxPyPzE(posPionCandidate->Px(), posPionCandidate->Py(), posPionCandidate->Pz(), posPionCandidate->E());

    AliAODConversionPhoton *posPionHandler = new AliAODConversionPhoton(posPionforHandler);
    delete posPionforHandler;

    fPosPionCandidates->Add(posPionHandler);
    if(!fDoLightOutput){
        fHistoPosPionPt[fiCut]->Fill( posPionCandidate->Pt() );
        fHistoPosPionPhi[fiCut]->Fill( posPionCandidate->Phi() );
    }
    if( fMCEvent ) {
      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();

      Int_t labelPosPion = TMath::Abs( posPionCandidate->GetLabel() );
      Bool_t posPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            if( labelPosPion>-1 && labelPosPion < fMCEvent->GetNumberOfTracks() ) {
        TParticle* posPion = fMCEvent->Particle(labelPosPion);
        if( posPion->GetPdgCode() ==  211 ){
            if(!fDoLightOutput){
                if( posPionIsPrimary ){
                    fHistoTruePosPionPt[fiCut]->Fill(posPionCandidate->Pt());
                }
                if( IsEtaPiPlPiMiPiZeroDaughter(labelPosPion) || IsOmegaPiPlPiMiPiZeroDaughter(labelPosPion) ) {
                    if(posPionIsPrimary){
                        fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt());
                    }
                }
            }
        }
      }
    }
  }


  for(UInt_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){
    AliESDtrack *negPionCandidate = fESDEvent->GetTrack(lGoodNegPionIndexPrev[i]);
    AliKFParticle negPionCandidateKF( *negPionCandidate->GetConstrainedParam(), 211 );

    for(UInt_t j = 0; j < lGoodPosPionIndexPrev.size(); j++){
      AliESDtrack *posPionCandidate = fESDEvent->GetTrack(lGoodPosPionIndexPrev[j]);
      AliKFParticle posPionCandidateKF( *posPionCandidate->GetConstrainedParam(), 211 );

      AliKFConversionPhoton* virtualPhoton = NULL;
      virtualPhoton = new AliKFConversionPhoton(negPionCandidateKF,posPionCandidateKF);
      AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
// 			primaryVertexImproved+=*virtualPhoton;
      virtualPhoton->SetProductionVertex(primaryVertexImproved);
      virtualPhoton->SetTrackLabels( lGoodPosPionIndexPrev[j], lGoodNegPionIndexPrev[i]);

      Int_t labeln=0;
      Int_t labelp=0;
      Int_t motherlabelp = 0;
      Int_t motherlabeln = 0;
      TParticle *fNegativeMCParticle =NULL;
      TParticle *fPositiveMCParticle =NULL;
      if( fMCEvent ) {
        labeln=TMath::Abs(negPionCandidate->GetLabel());
        labelp=TMath::Abs(posPionCandidate->GetLabel());
                if(labeln>-1) fNegativeMCParticle = fMCEvent->Particle(labeln);
                if(labelp>-1) fPositiveMCParticle = fMCEvent->Particle(labelp);
        // check whether MC particles exist, else abort
        if (fNegativeMCParticle == NULL || fPositiveMCParticle == NULL) return;

        motherlabeln = fNegativeMCParticle->GetMother(0);
        motherlabelp = fPositiveMCParticle->GetMother(0);
        virtualPhoton->SetMCLabelPositive(labelp);
        virtualPhoton->SetMCLabelNegative(labeln);

      }

      AliAODConversionPhoton *vParticle = new AliAODConversionPhoton(virtualPhoton); //To apply mass 2 pion mass cut
      if(!fDoLightOutput){
        if (fMCEvent &&(fDoMesonQA>0)){
          if (fPositiveMCParticle && fNegativeMCParticle ) {
            if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
              if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
                if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
                  fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                  if (motherlabeln == motherlabelp){
                    fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                    if( IsEtaPiPlPiMiPiZeroDaughter(labeln) ) { //|| IsOmegaPiPlPiMiPiZeroDaughter(labeln)
                      fHistoTruePionPionFromEtaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                    }
                    if( IsOmegaPiPlPiMiPiZeroDaughter(labeln) ) { //||
                      fHistoTruePionPionFromOmegaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                    }
                  }
                }
              }
            } else {
              if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
                fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                if (motherlabeln == motherlabelp){
                  fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                  if( IsEtaPiPlPiMiPiZeroDaughter(labeln) ) { //|| IsOmegaPiPlPiMiPiZeroDaughter(labeln)
                    fHistoTruePionPionFromEtaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                  }
                  if( IsOmegaPiPlPiMiPiZeroDaughter(labeln) ) { //||
                    fHistoTruePionPionFromOmegaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
                  }
                }
              }
            }
          }
        }
      }

      if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
        if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
                    fGoodVirtualParticles->Add( vParticle );
                    if(!fDoLightOutput){
                      fHistoPionPionInvMassPt[fiCut]->Fill( vParticle->GetMass(),vParticle->Pt());
                    }
        }else{
          delete vParticle;
          vParticle=0x0;
        }
      } else {
                fGoodVirtualParticles->Add( vParticle );
                if(!fDoLightOutput){
                  fHistoPionPionInvMassPt[fiCut]->Fill( vParticle->GetMass(),vParticle->Pt());
                }
      }

      Double_t clsToFPos = -1.0;
      Double_t clsToFNeg = -1.0;

      Float_t dcaToVertexXYPos = -1.0;
      Float_t dcaToVertexZPos  = -1.0;
      Float_t dcaToVertexXYNeg = -1.0;
      Float_t dcaToVertexZNeg  = -1.0;

      if ( fDoMesonQA>0 ) {
        clsToFPos = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(posPionCandidate);
        clsToFNeg = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(negPionCandidate);

        Float_t bPos[2];
        Float_t bCovPos[3];
        posPionCandidate->GetImpactParameters(bPos,bCovPos);
        if (bCovPos[0]<=0 || bCovPos[2]<=0) {
          AliDebug(1, "Estimated b resolution lower or equal zero!");
          bCovPos[0]=0; bCovPos[2]=0;
        }

        Float_t bNeg[2];
        Float_t bCovNeg[3];
        posPionCandidate->GetImpactParameters(bNeg,bCovNeg);
        if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
          AliDebug(1, "Estimated b resolution lower or equal zero!");
          bCovNeg[0]=0; bCovNeg[2]=0;
        }

        dcaToVertexXYPos = bPos[0];
        dcaToVertexZPos  = bPos[1];
        dcaToVertexXYNeg = bNeg[0];
        dcaToVertexZNeg  = bNeg[1];

        if(!fDoLightOutput){
          fHistoNegPionEta[fiCut]->Fill( negPionCandidate->Eta() );
          fHistoPosPionEta[fiCut]->Fill( posPionCandidate->Eta() );

          fHistoNegPionClsTPC[fiCut]->Fill(clsToFNeg,negPionCandidate->Pt());
          fHistoPosPionClsTPC[fiCut]->Fill(clsToFPos,posPionCandidate->Pt());

          fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionCandidate->Pt() );
          fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionCandidate->Pt() );
          fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionCandidate->Pt() );
          fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionCandidate->Pt() );

          fHistoPionTPCdEdxNSigma[fiCut]->Fill( posPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(posPionCandidate, AliPID::kPion) );
          fHistoPionTPCdEdxNSigma[fiCut]->Fill( negPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(negPionCandidate, AliPID::kPion) );

          fHistoPionTPCdEdx[fiCut]->Fill( posPionCandidate->P(), TMath::Abs(posPionCandidate->GetTPCsignal()));
          fHistoPionTPCdEdx[fiCut]->Fill( negPionCandidate->P(), TMath::Abs(negPionCandidate->GetTPCsignal()));
        }
      }

      delete virtualPhoton;
      virtualPhoton=NULL;
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessMCParticles(){

  // Loop over all primary MC particle
  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  for(Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){

      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader
          = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }

      if(!fDoLightOutput){
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
          // find MC photons
          if (fNeutralPionMode < 2){
            if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
              fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==111){
                  if (fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1){
                    if ( fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 221 ||
                         fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 223 ){
                      if ( fMCEvent->Particle(particle->GetMother(0))->GetNDaughters()==3 )
                        fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All photons from eta or omega via pi0
                    }
                  }
                }
              }
            }
          } else if (fNeutralPionMode == 2){
            if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
              fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==111){
                  if (fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1){
                    if ( fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 221 ||
                         fMCEvent->Particle((fMCEvent->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 223 ){
                      if ( fMCEvent->Particle(particle->GetMother(0))->GetNDaughters()==3 )
                        fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All photons from eta or omega via pi0
                    }
                  }
                }
              }
            }
          }
          if (fNeutralPionMode < 2){
            if (((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
              fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt());
            } // Converted MC Gamma
          }
          if(((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(i,fMCEvent)){
            if( particle->GetPdgCode() == 211){
              fHistoMCAllPosPionsPt[fiCut]->Fill(particle->Pt()); // All pos pions
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==221 || fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==223)
                  fHistoMCPosPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All pos from eta or omega
              }
            }
            if( particle->GetPdgCode() == -211){
              fHistoMCAllNegPionsPt[fiCut]->Fill(particle->Pt()); // All neg pions
              if(particle->GetMother(0) >-1){
                if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==221 || fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==223 )
                  fHistoMCNegPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All pos from eta or omega
              }
            }
          }
        }
      }

        // \eta -> pi+ pi- \gamma
        Int_t labelNeutPion = -1;
        Int_t labelNegPion = -1;
        Int_t labelPosPion = -1;

        if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCPiPlPiMiPiZero(particle,fMCEvent,labelNegPion,labelPosPion,labelNeutPion,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
          Float_t weighted= 1;
          if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) {
            if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
              if (particle->Pt()>0.005){
                weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
              }
            }
          }
          if(particle->GetPdgCode() == 221)fHistoMCEtaPiPlPiMiPiZeroPt[fiCut]->Fill(particle->Pt(), weighted); 						// All MC Eta in respective decay channel
          if(particle->GetPdgCode() == 223)fHistoMCOmegaPiPlPiMiPiZeroPt[fiCut]->Fill(particle->Pt(), weighted); 						// All MC Omega in respective decay channel

          if(labelNeutPion>-1){
            TParticle *neutPion    = fMCEvent->Particle(labelNeutPion);
            if(neutPion->GetDaughter(0)>-1 && neutPion->GetDaughter(1)>-1){
              TParticle *gamma1 = fMCEvent->Particle(neutPion->GetDaughter(0));
              TParticle *gamma2 = fMCEvent->Particle(neutPion->GetDaughter(1));
              Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, neutPion->GetDaughter(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, neutPion->GetDaughter(1), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              Bool_t kNegPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
              Bool_t kPosPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

              if (fNeutralPionMode == 0){
                if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                    ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE) &&					// test first daugther of pi0
                    ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma2,fMCEvent,kFALSE) &&					// test second daughter of pi0
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCEvent) &&								// test negative pion
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCEvent) 								// test positive pion
                ) {
                  if(particle->GetPdgCode() == 221) fHistoMCEtaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Eta pi+ pi- pi0 with gamma's and e+e- in acc
                  if(particle->GetPdgCode() == 223) fHistoMCOmegaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Omega pi+ pi- pi0 with gamma's and e+e- in acc
                }
              } else if (fNeutralPionMode == 1){ // mixed mode
                  // check if within PCM acceptance first
                  if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                          ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCEvent,kFALSE) &&					// test first daugther of pi0
                          ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma2,fMCEvent,kFALSE) &&					// test second daughter of pi0
                          ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCEvent) &&								// test negative pion
                          ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCEvent) 								// test positive pion
                          ) {
                      // check acceptance of clusters as well, true if one of them points into the Calo acceptance
                      if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent) ||
                              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma2,fMCEvent) ){
                          if(particle->GetPdgCode() == 221) fHistoMCEtaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Eta pi+ pi- pi0 with gamma's and e+e- in acc
                          if(particle->GetPdgCode() == 223) fHistoMCOmegaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Omega pi+ pi- pi0 with gamma's and e+e- in acc
                      }
                  }
              } else if (fNeutralPionMode == 2){
                if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCEvent) &&					// test first daugther of pi0
                    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma2,fMCEvent) &&					// test second daughter of pi0
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCEvent) &&								// test negative pion
                    ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCEvent) 								// test positive pion
                ) {
                  if(particle->GetPdgCode() == 221) fHistoMCEtaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Eta pi+ pi- pi0 with gamma's and e+e- in acc
                  if(particle->GetPdgCode() == 223) fHistoMCOmegaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Omega pi+ pi- pi0 with gamma's and e+e- in acc
                }
              }
            }
          }
        }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::CalculateMesonCandidates(){

  // Conversion Gammas
  if( fNeutralPionCandidates->GetEntries() > 0 && fGoodVirtualParticles->GetEntries() > 0 ){
    for(Int_t mesonIndex=0; mesonIndex<fNeutralPionCandidates->GetEntries(); mesonIndex++){
      AliAODConversionMother *neutralPion=dynamic_cast<AliAODConversionMother*>(fNeutralPionCandidates->At(mesonIndex));
      if (neutralPion==NULL) continue;

      for(Int_t virtualParticleIndex=0;virtualParticleIndex<fGoodVirtualParticles->GetEntries();virtualParticleIndex++){

                AliAODConversionPhoton *vParticle=dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualParticles->At(virtualParticleIndex));
        if (vParticle==NULL) continue;
        //Check for same Electron ID

        AliAODConversionMother *mesoncand = new AliAODConversionMother(neutralPion,vParticle);
        mesoncand->SetLabels(mesonIndex,virtualParticleIndex);
        if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesoncand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){

          AliESDtrack *negPionCandidatetmp = (AliESDtrack*) fESDEvent->GetTrack(vParticle->GetTrackLabel(1));
          if(negPionCandidatetmp==NULL){ delete mesoncand; continue;}
          AliAODConversionMother *NegPiontmp = new AliAODConversionMother();
          NegPiontmp->SetPxPyPzE(negPionCandidatetmp->Px(), negPionCandidatetmp->Py(), negPionCandidatetmp->Pz(), negPionCandidatetmp->E());

          AliESDtrack *posPionCandidatetmp = (AliESDtrack*) fESDEvent->GetTrack(vParticle->GetTrackLabel(0));
          if(posPionCandidatetmp==NULL){ delete NegPiontmp; delete mesoncand; continue;}
          AliAODConversionMother *PosPiontmp = new AliAODConversionMother();
          PosPiontmp->SetPxPyPzE(posPionCandidatetmp->Px(), posPionCandidatetmp->Py(), posPionCandidatetmp->Pz(), posPionCandidatetmp->E());

          if(KinematicCut(NegPiontmp, PosPiontmp, neutralPion, mesoncand)){
              if(!fDoLightOutput){
                  fHistoAngleOmegaPiZero[fiCut]->Fill(mesoncand->Pt(),neutralPion->Angle(mesoncand->Vect()));
                  fHistoAngleOmegaPiPl[fiCut]->Fill(mesoncand->Pt(),PosPiontmp->Angle(mesoncand->Vect()));
                  fHistoAngleOmegaPiMi[fiCut]->Fill(mesoncand->Pt(),NegPiontmp->Angle(mesoncand->Vect()));
                  fHistoAnglePiZeroPiMi[fiCut]->Fill(mesoncand->Pt(),NegPiontmp->Angle(neutralPion->Vect()));
                  fHistoAnglePiPlPiMi[fiCut]->Fill(mesoncand->Pt(),NegPiontmp->Angle(PosPiontmp->Vect()));
                  fHistoAnglePiPlPiZero[fiCut]->Fill(mesoncand->Pt(),PosPiontmp->Angle(neutralPion->Vect()));
                  fHistoAngleOmegaPiPlPiMi[fiCut]->Fill(mesoncand->Pt(),vParticle->Angle(mesoncand->Vect()));
                  fHistoAngleSum[fiCut]->Fill(mesoncand->Pt(),((PosPiontmp->Angle(mesoncand->Vect()))+(NegPiontmp->Angle(PosPiontmp->Vect()))+(PosPiontmp->Angle(neutralPion->Vect()))));

                  //Double_t sparesFill[4] = {mesoncand->M(),mesoncand->Pt(),(Double_t)zbin,(Double_t)mbin};
                  //fTHnSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
              }

              // Subtract mass of used pi0 candidate and then add PDG mass to get to right range again
              fHistoMotherInvMassSubPi0[fiCut]->Fill(mesoncand->M()-(neutralPion->M()-fPDGMassPi0),mesoncand->Pt());

              // Fix Pz of pi0 candidate to match pi0 PDG mass
              AliAODConversionMother *Pi0tmp = new AliAODConversionMother();
              Pi0tmp->SetPxPyPzE(neutralPion->Px(), neutralPion->Py(), neutralPion->Pz(), neutralPion->Energy());
              FixPzToMatchPDGInvMassPi0(Pi0tmp);
              AliAODConversionMother *mesontmp = new AliAODConversionMother(Pi0tmp,vParticle);
              fHistoMotherInvMassFixedPzPi0[fiCut]->Fill(mesontmp->M(),mesontmp->Pt());
              delete Pi0tmp;
              delete mesontmp;
              fHistoMotherInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt());
            if(fMCEvent){
              ProcessTrueMesonCandidates(mesoncand,neutralPion,vParticle);
            }
          }else{
            if(!fDoLightOutput){
              fHistoMotherInvMassPtRejectedKinematic[fiCut]->Fill(mesoncand->M(),mesoncand->Pt());
            }
          }
          if(!fDoLightOutput){
            delete NegPiontmp;
            delete PosPiontmp;
          }
        }
        delete mesoncand;
        mesoncand=0x0;
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::CalculateBackground(){

  /* Event mixing histo explanation
  *
  * fHistoBackInvMassPtGroup1 => pi+ and pi- from same event
  * fHistoBackInvMassPtGroup2 => pi+ and pi0 from same event
  * fHistoBackInvMassPtGroup3 => pi- and pi0 from same event
  * fHistoBackInvMassPtGroup4 => no pions from same event
  */

  // Get multiplicity and zbin from fBGHandler
  Int_t zbin= fBGHandlerPiMi[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
  Int_t mbin = 0;

  // Multiplicity can be determined either by number of cluster candidates or track mulitiplicity
  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    mbin = fBGHandlerPiMi[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
  } else {
    if (fNeutralPionMode < 2) mbin = fBGHandlerPiMi[fiCut]->GetMultiplicityBinIndex(fGoodConvGammas->GetEntries());
    else mbin = fBGHandlerPiMi[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
  }

  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertexPl = NULL;
  AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertexMi = NULL;

  // Get N of Pi0 according to chosen mix mode
  Int_t NPi0Candidates = 0;
  if( (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) || (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides())){
      NPi0Candidates = fNeutralPionSidebandCandidates->GetEntries();
   }else{
      NPi0Candidates = fNeutralPionCandidates->GetEntries();
  }
  // Begin loop over all Pi0 candidates
 for(Int_t iCurrentPi0=0; iCurrentPi0<NPi0Candidates; iCurrentPi0++){
     AliAODConversionMother* EventPiZeroGoodMeson;
     if( (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixing()) ||  (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseSidebandMixingBothSides())){
         EventPiZeroGoodMeson = (AliAODConversionMother*)(fNeutralPionSidebandCandidates->At(iCurrentPi0));
     }else{
         EventPiZeroGoodMeson = (AliAODConversionMother*)(fNeutralPionCandidates->At(iCurrentPi0));
     }

if(!(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseLikeSignMixing())){
    // Begin loop over BG events for Pi+
    for(Int_t nEventsInBGPl=0;nEventsInBGPl<fBGHandlerPiPl[fiCut]->GetNBGEvents();nEventsInBGPl++){

      // Store all Pi+ of current event in right binning in vector
      AliGammaConversionMotherAODVector *EventPiPlMeson = fBGHandlerPiPl[fiCut]->GetBGGoodMesons(zbin,mbin,nEventsInBGPl);

      // Begin loop over BG events for Pi-
      for(Int_t nEventsInBGMi=0;nEventsInBGMi<fBGHandlerPiMi[fiCut]->GetNBGEvents();nEventsInBGMi++){
        AliGammaConversionMotherAODVector *EventPiMiMeson = fBGHandlerPiMi[fiCut]->GetBGGoodMesons(zbin,mbin,nEventsInBGMi);

        // If one of the events isn't found skip to next one
        if((EventPiMiMeson && EventPiPlMeson) == kFALSE) continue;

        // Determine Background event vertex
        if(fMoveParticleAccordingToVertex == kTRUE){
          bgEventVertexPl = fBGHandlerPiPl[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBGPl);
          bgEventVertexMi = fBGHandlerPiMi[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBGMi);
        }
        // Loop over all Pi+
        for(UInt_t iCurrentPiPl = 0; iCurrentPiPl<EventPiPlMeson->size();iCurrentPiPl++){
            AliAODConversionMother EventPiPlGoodMeson= (AliAODConversionMother)(*(EventPiPlMeson->at(iCurrentPiPl)));

            // Move Vertex
            if(fMoveParticleAccordingToVertex == kTRUE){
              MoveParticleAccordingToVertex(&EventPiPlGoodMeson, bgEventVertexPl);
            }

            // Combine Pi+ and Pi0
            AliAODConversionMother *PiPlPiZeroBackgroundCandidate = new AliAODConversionMother(EventPiZeroGoodMeson, &EventPiPlGoodMeson);

            for(UInt_t iCurrentPiMi = 0; iCurrentPiMi<EventPiMiMeson->size();iCurrentPiMi++){
              AliAODConversionMother EventPiMiGoodMeson = (AliAODConversionMother)(*(EventPiMiMeson->at(iCurrentPiMi)));

              // Move Vertex
              if(fMoveParticleAccordingToVertex == kTRUE){
                MoveParticleAccordingToVertex(&EventPiMiGoodMeson, bgEventVertexMi);
              }

              // Mass cut (pi+pi-)
              if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
                AliAODConversionMother *backPiPlPiMiCandidate = new AliAODConversionMother(&EventPiPlGoodMeson,&EventPiMiGoodMeson);
                if (backPiPlPiMiCandidate->M() >= ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
                  delete backPiPlPiMiCandidate;
                  backPiPlPiMiCandidate = 0x0;
                  continue;
                }
                delete backPiPlPiMiCandidate;
                backPiPlPiMiCandidate = 0x0;
              }

              // Create (final) Candidate
              AliAODConversionMother *PiPlPiMiPiZeroBackgroundCandidate = new AliAODConversionMother(PiPlPiZeroBackgroundCandidate,&EventPiMiGoodMeson);

              // Check if candidate survives meson cut
              if(  ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(PiPlPiMiPiZeroBackgroundCandidate,kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){

                // Check if candidate survives kinematic cut
                if(KinematicCut(&EventPiMiGoodMeson, &EventPiPlGoodMeson,EventPiZeroGoodMeson,PiPlPiMiPiZeroBackgroundCandidate)){
                  // Create temporary mesons to be able to fix pz
                  AliAODConversionMother *Pi0tmp = new AliAODConversionMother();
                  Pi0tmp->SetPxPyPzE(EventPiZeroGoodMeson->Px(), EventPiZeroGoodMeson->Py(), EventPiZeroGoodMeson->Pz(), EventPiZeroGoodMeson->Energy());
                  FixPzToMatchPDGInvMassPi0(Pi0tmp);
                  AliAODConversionMother *PiMiPiZerotmp = new AliAODConversionMother(&EventPiMiGoodMeson,Pi0tmp);
                  AliAODConversionMother *PiPlPiMiPiZerotmp = new AliAODConversionMother(&EventPiPlGoodMeson,PiMiPiZerotmp);

                  if (nEventsInBGMi != nEventsInBGPl){
                    // Pi+ and Pi- don't come from the same event (but different than pi0 event)
                    // Fill histograms
                    fHistoBackInvMassPtGroup4[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M(),PiPlPiMiPiZeroBackgroundCandidate->Pt());
                    fHistoBackInvMassPtGroup4SubPi0[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M()-(EventPiZeroGoodMeson->M()-fPDGMassPi0),PiPlPiMiPiZeroBackgroundCandidate->Pt());
                    fHistoBackInvMassPtGroup4FixedPzPi0[fiCut]->Fill(PiPlPiMiPiZerotmp->M(),PiPlPiMiPiZerotmp->Pt());

                  } else if(nEventsInBGMi==nEventsInBGPl){
                    // Pi+ and Pi- come from the same event (but different than pi0 event)
                    fHistoBackInvMassPtGroup1[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M(),PiPlPiMiPiZeroBackgroundCandidate->Pt());
                    fHistoBackInvMassPtGroup1SubPi0[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M()-(EventPiZeroGoodMeson->M()-fPDGMassPi0),PiPlPiMiPiZeroBackgroundCandidate->Pt());
                    fHistoBackInvMassPtGroup1FixedPzPi0[fiCut]->Fill(PiPlPiMiPiZerotmp->M(),PiPlPiMiPiZerotmp->Pt());
                  }

                  delete Pi0tmp;
                  delete PiMiPiZerotmp;
                  delete PiPlPiMiPiZerotmp;

                  delete PiPlPiMiPiZeroBackgroundCandidate;
                  PiPlPiMiPiZeroBackgroundCandidate = 0x0;
                }
              }
              if(PiPlPiMiPiZeroBackgroundCandidate!=0x0){
                  delete PiPlPiMiPiZeroBackgroundCandidate;
                  PiPlPiMiPiZeroBackgroundCandidate = 0x0;
              }
            } // end pi- loop
            if(PiPlPiZeroBackgroundCandidate!=0x0){
                delete PiPlPiZeroBackgroundCandidate;
                PiPlPiZeroBackgroundCandidate = 0x0;
            }
        } // end pi+ loop
      } // end loop over all pi- event
    } // end loop over pi+ events

    // Loop over all pi+ events(from Handler)
    for(Int_t nEventsInBGPl=0;nEventsInBGPl<fBGHandlerPiPl[fiCut]->GetNBGEvents();nEventsInBGPl++){
      // Store all Pi+ of current event in right binning in vector
      AliGammaConversionMotherAODVector *EventPiPlMeson = fBGHandlerPiPl[fiCut]->GetBGGoodMesons(zbin,mbin,nEventsInBGPl);

      // Determine Vertex
      if(fMoveParticleAccordingToVertex == kTRUE){
        bgEventVertexPl = fBGHandlerPiPl[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBGPl);
      }
      // Begin loop over all pi+ in ecent
      for(UInt_t iCurrentPiPl = 0; iCurrentPiPl<EventPiPlMeson->size();iCurrentPiPl++){
        AliAODConversionMother EventPiPlGoodMeson= (AliAODConversionMother)(*(EventPiPlMeson->at(iCurrentPiPl)));

        // Move vertex
        if(fMoveParticleAccordingToVertex == kTRUE){
          MoveParticleAccordingToVertex(&EventPiPlGoodMeson, bgEventVertexPl);
        }
        // Combine Pi+ and Pi0
        AliAODConversionMother *PiPlPiZeroBackgroundCandidate = new AliAODConversionMother(EventPiZeroGoodMeson, &EventPiPlGoodMeson);
        // Loop over all pi- (from current event)
        for(Int_t iCurrentPiMi=0; iCurrentPiMi<fNegPionCandidates->GetEntries(); iCurrentPiMi++){
          AliAODConversionMother EventPiNegGoodMeson = *(AliAODConversionMother*)(fNegPionCandidates->At(iCurrentPiMi));

          // Mass cut on pi+pi-
          if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
            AliAODConversionMother *backPiPlPiMiCandidate = new AliAODConversionMother(&EventPiPlGoodMeson,&EventPiNegGoodMeson);
            if (backPiPlPiMiCandidate->M() >= ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
              delete backPiPlPiMiCandidate;
              backPiPlPiMiCandidate = 0x0;
              continue;
            }
            delete backPiPlPiMiCandidate;
            backPiPlPiMiCandidate = 0x0;
          }

          // Create (final) Candidate
          AliAODConversionMother *PiPlPiMiPiZeroBackgroundCandidate = new AliAODConversionMother(PiPlPiZeroBackgroundCandidate,&EventPiNegGoodMeson);

          // Check if candidate survives meson cut
          if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(PiPlPiMiPiZeroBackgroundCandidate,kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){

            // Check if candidate survives kinematic cut
            if(KinematicCut(&EventPiNegGoodMeson, &EventPiPlGoodMeson, EventPiZeroGoodMeson,PiPlPiMiPiZeroBackgroundCandidate)){

              // Create temporary mesons to be able to fix pz
              AliAODConversionMother *Pi0tmp = new AliAODConversionMother();
              Pi0tmp->SetPxPyPzE(EventPiZeroGoodMeson->Px(), EventPiZeroGoodMeson->Py(), EventPiZeroGoodMeson->Pz(), EventPiZeroGoodMeson->Energy());
              FixPzToMatchPDGInvMassPi0(Pi0tmp);
              AliAODConversionMother *PiMiPiZerotmp = new AliAODConversionMother(&EventPiNegGoodMeson,Pi0tmp);
              AliAODConversionMother *PiPlPiMiPiZerotmp = new AliAODConversionMother(&EventPiPlGoodMeson,PiMiPiZerotmp);

              // Fill histograms (pi- and pi0 from same event)
              fHistoBackInvMassPtGroup3[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M(),PiPlPiMiPiZeroBackgroundCandidate->Pt());
              fHistoBackInvMassPtGroup3SubPi0[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M()-(EventPiZeroGoodMeson->M()-fPDGMassPi0),PiPlPiMiPiZeroBackgroundCandidate->Pt());
              fHistoBackInvMassPtGroup3FixedPzPi0[fiCut]->Fill(PiPlPiMiPiZerotmp->M(),PiPlPiMiPiZerotmp->Pt());

              delete Pi0tmp;
              delete PiMiPiZerotmp;
              delete PiPlPiMiPiZerotmp;

              delete PiPlPiMiPiZeroBackgroundCandidate;
              PiPlPiMiPiZeroBackgroundCandidate = 0x0;
            }
          }
          if(PiPlPiMiPiZeroBackgroundCandidate!=0x0){
              delete PiPlPiMiPiZeroBackgroundCandidate;
              PiPlPiMiPiZeroBackgroundCandidate = 0x0;
          }
        } // End loop pi- (from current event)
        if(PiPlPiZeroBackgroundCandidate!=0x0){
            delete PiPlPiZeroBackgroundCandidate;
            PiPlPiZeroBackgroundCandidate = 0x0;
        }
      } // End loop pi+
    } // end loop over pi+ events

    // Loop over all pi- events(from Handler)
    for(Int_t nEventsInBGMi=0;nEventsInBGMi<fBGHandlerPiPl[fiCut]->GetNBGEvents();nEventsInBGMi++){
      // Store all Pi- of current event in right binning in vector
      AliGammaConversionMotherAODVector *EventPiMiMeson = fBGHandlerPiMi[fiCut]->GetBGGoodMesons(zbin,mbin,nEventsInBGMi);

      // Determine vertex
      if(fMoveParticleAccordingToVertex == kTRUE){
        bgEventVertexMi = fBGHandlerPiMi[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBGMi);
      }

      // Begin loop over all pi- in event
      for(UInt_t iCurrentPiMi = 0; iCurrentPiMi<EventPiMiMeson->size();iCurrentPiMi++){
        AliAODConversionMother EventPiMiGoodMeson= (AliAODConversionMother)(*(EventPiMiMeson->at(iCurrentPiMi)));

        // move vertex
        if(fMoveParticleAccordingToVertex == kTRUE){
          MoveParticleAccordingToVertex(&EventPiMiGoodMeson, bgEventVertexMi);
        }


        // Combine Pi- and Pi0
        AliAODConversionMother *PiMiPiZeroBackgroundCandidate = new AliAODConversionMother(EventPiZeroGoodMeson, &EventPiMiGoodMeson);

        // Loop over all pi+ (from current event)
        for(Int_t iCurrentPiPl=0; iCurrentPiPl<fPosPionCandidates->GetEntries(); iCurrentPiPl++){
          AliAODConversionMother EventPiPlGoodMeson = *(AliAODConversionMother*)(fPosPionCandidates->At(iCurrentPiPl));

          // Mass cut on pi+pi-
          if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
            AliAODConversionMother *backPiPlPiMiCandidate = new AliAODConversionMother(&EventPiPlGoodMeson,&EventPiMiGoodMeson);
            if (backPiPlPiMiCandidate->M() >= ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
              delete backPiPlPiMiCandidate;
              backPiPlPiMiCandidate = 0x0;
              continue;
            }
            delete backPiPlPiMiCandidate;
            backPiPlPiMiCandidate = 0x0;
          }

          // Create (final) Candidate
          AliAODConversionMother *PiPlPiMiPiZeroBackgroundCandidate = new AliAODConversionMother(PiMiPiZeroBackgroundCandidate,&EventPiPlGoodMeson);

          // Check if candidate survives meson cut
          if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(PiMiPiZeroBackgroundCandidate,kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){

            // Check if candidate survives kinematic cut
            if(KinematicCut(&EventPiMiGoodMeson, &EventPiPlGoodMeson, EventPiZeroGoodMeson,PiPlPiMiPiZeroBackgroundCandidate)){

              // Create temporary mesons to be able to fix pz
              AliAODConversionMother *Pi0tmp = new AliAODConversionMother();
              Pi0tmp->SetPxPyPzE(EventPiZeroGoodMeson->Px(), EventPiZeroGoodMeson->Py(), EventPiZeroGoodMeson->Pz(), EventPiZeroGoodMeson->Energy());
              FixPzToMatchPDGInvMassPi0(Pi0tmp);
              AliAODConversionMother *PiMiPiZerotmp = new AliAODConversionMother(&EventPiMiGoodMeson,Pi0tmp);
              AliAODConversionMother *PiPlPiMiPiZerotmp = new AliAODConversionMother(&EventPiPlGoodMeson,PiMiPiZerotmp);

              // Fill histograms (pi+ and pi0 from same event)
              fHistoBackInvMassPtGroup2[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M(),PiPlPiMiPiZeroBackgroundCandidate->Pt());
              fHistoBackInvMassPtGroup2SubPi0[fiCut]->Fill(PiPlPiMiPiZeroBackgroundCandidate->M()-(EventPiZeroGoodMeson->M()-fPDGMassPi0),PiPlPiMiPiZeroBackgroundCandidate->Pt());
              fHistoBackInvMassPtGroup2FixedPzPi0[fiCut]->Fill(PiPlPiMiPiZerotmp->M(),PiPlPiMiPiZerotmp->Pt());

              delete Pi0tmp;
              delete PiMiPiZerotmp;
              delete PiPlPiMiPiZerotmp;

              delete PiPlPiMiPiZeroBackgroundCandidate;
              PiPlPiMiPiZeroBackgroundCandidate = 0x0;
            }
          }
          if(PiPlPiMiPiZeroBackgroundCandidate!=0x0){
              delete PiPlPiMiPiZeroBackgroundCandidate;
              PiPlPiMiPiZeroBackgroundCandidate = 0x0;
          }
        } // End loop pi+ (from current event)
        if(PiMiPiZeroBackgroundCandidate!=0x0){
            delete PiMiPiZeroBackgroundCandidate;
            PiMiPiZeroBackgroundCandidate = 0x0;
        }
      } // End loop pi-
    } // end loop over pi+ events
 /*
  * LikeSign Mixing
  */
 } else if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseLikeSignMixing()){
     // Loops for Pi0Pi+Pi+ LikeSign mixing
     for(Int_t iCurrentPiPl=0; iCurrentPiPl<fPosPionCandidates->GetEntries(); iCurrentPiPl++){

         AliAODConversionMother EventPiPlGoodMeson = *(AliAODConversionMother*)(fPosPionCandidates->At(iCurrentPiPl));

          for(Int_t iCurrentPiPl2=0; iCurrentPiPl2<fPosPionCandidates->GetEntries(); iCurrentPiPl2++){

              if(iCurrentPiPl!=iCurrentPiPl2){ // dont mix same particle
                  AliAODConversionMother EventPiPlGoodMeson2 = *(AliAODConversionMother*)(fPosPionCandidates->At(iCurrentPiPl2));

                  // Combine Pi+ and Pi0
                  AliAODConversionMother *PiPlPiZeroBackgroundCandidate = new AliAODConversionMother(&EventPiPlGoodMeson, EventPiZeroGoodMeson);

                  // Mass cut on pi+pi+
                  if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
                    AliAODConversionMother *backPiPlPiPlCandidate = new AliAODConversionMother(&EventPiPlGoodMeson,&EventPiPlGoodMeson2);
                    if (backPiPlPiPlCandidate->M() >= ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
                      delete backPiPlPiPlCandidate;
                      backPiPlPiPlCandidate = 0x0;
                      continue;
                    }
                    delete backPiPlPiPlCandidate;
                    backPiPlPiPlCandidate = 0x0;
                  }

                  // Create (final) Candidate
                  AliAODConversionMother *PiPlPiPlPiZeroBackgroundCandidate = new AliAODConversionMother(PiPlPiZeroBackgroundCandidate, &EventPiPlGoodMeson2);

                  // Check if candidate survives meson cut
                  if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(PiPlPiZeroBackgroundCandidate,kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){

                    // Check if candidate survives kinematic cut
                    if(KinematicCut(&EventPiPlGoodMeson, &EventPiPlGoodMeson2, EventPiZeroGoodMeson,PiPlPiPlPiZeroBackgroundCandidate)){

                      // Create temporary mesons to be able to fix pz
                      AliAODConversionMother *Pi0tmp = new AliAODConversionMother();
                      Pi0tmp->SetPxPyPzE(EventPiZeroGoodMeson->Px(), EventPiZeroGoodMeson->Py(), EventPiZeroGoodMeson->Pz(), EventPiZeroGoodMeson->Energy());
                      FixPzToMatchPDGInvMassPi0(Pi0tmp);
                      AliAODConversionMother *PiPlPiZerotmp = new AliAODConversionMother(&EventPiPlGoodMeson,Pi0tmp);
                      AliAODConversionMother *PiPlPiPlPiZerotmp = new AliAODConversionMother(&EventPiPlGoodMeson2,PiPlPiZerotmp);

                      // Fill histograms (likesign)
                      fHistoMotherLikeSignBackInvMassPt[fiCut]->Fill(PiPlPiPlPiZeroBackgroundCandidate->M(),PiPlPiPlPiZeroBackgroundCandidate->Pt());
                      fHistoMotherLikeSignBackInvMassSubPi0Pt[fiCut]->Fill(PiPlPiPlPiZeroBackgroundCandidate->M()-(EventPiZeroGoodMeson->M()-fPDGMassPi0),PiPlPiPlPiZeroBackgroundCandidate->Pt());
                      fHistoMotherLikeSignBackInvMassFixedPzPi0Pt[fiCut]->Fill(PiPlPiPlPiZerotmp->M(),PiPlPiPlPiZerotmp->Pt());

                      delete Pi0tmp;
                      delete PiPlPiZerotmp;
                      delete PiPlPiPlPiZerotmp;

                      delete PiPlPiPlPiZeroBackgroundCandidate;
                      PiPlPiPlPiZeroBackgroundCandidate = 0x0;
                    }
                  }


              }
          } // end of iCurrentPiPl2
     }// end of iCurrenPiPl

     // Loops for Pi0Pi-Pi- LikeSign mixing
     for(Int_t iCurrentPiMi=0; iCurrentPiMi<fNegPionCandidates->GetEntries(); iCurrentPiMi++){

         AliAODConversionMother EventPiMiGoodMeson = *(AliAODConversionMother*)(fNegPionCandidates->At(iCurrentPiMi));

          for(Int_t iCurrentPiMi2=0; iCurrentPiMi2<fNegPionCandidates->GetEntries(); iCurrentPiMi2++){

              if(iCurrentPiMi!=iCurrentPiMi2){ // dont mix same particle
                  AliAODConversionMother EventPiMiGoodMeson2 = *(AliAODConversionMother*)(fNegPionCandidates->At(iCurrentPiMi2));

                  // Combine Pi- and Pi0
                  AliAODConversionMother *PiMiPiZeroBackgroundCandidate = new AliAODConversionMother(&EventPiMiGoodMeson, EventPiZeroGoodMeson);

                  // Mass cut on pi-pi-
                  if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
                    AliAODConversionMother *backPiMiPiMiCandidate = new AliAODConversionMother(&EventPiMiGoodMeson,&EventPiMiGoodMeson2);
                    if (backPiMiPiMiCandidate->M() >= ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
                      delete backPiMiPiMiCandidate;
                      backPiMiPiMiCandidate = 0x0;
                      continue;
                    }
                    delete backPiMiPiMiCandidate;
                    backPiMiPiMiCandidate = 0x0;
                  }

                  // Create (final) Candidate
                  AliAODConversionMother *PiMiPiMiPiZeroBackgroundCandidate = new AliAODConversionMother(PiMiPiZeroBackgroundCandidate, &EventPiMiGoodMeson2);

                  // Check if candidate survives meson cut
                  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(PiMiPiZeroBackgroundCandidate,kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){

                    // Check if candidate survives kinematic cut
                    if(KinematicCut(&EventPiMiGoodMeson, &EventPiMiGoodMeson2, EventPiZeroGoodMeson,PiMiPiMiPiZeroBackgroundCandidate)){

                      // Create temporary mesons to be able to fix pz
                      AliAODConversionMother *Pi0tmp = new AliAODConversionMother();
                      Pi0tmp->SetPxPyPzE(EventPiZeroGoodMeson->Px(), EventPiZeroGoodMeson->Py(), EventPiZeroGoodMeson->Pz(), EventPiZeroGoodMeson->Energy());
                      FixPzToMatchPDGInvMassPi0(Pi0tmp);
                      AliAODConversionMother *PiMiPiZerotmp = new AliAODConversionMother(&EventPiMiGoodMeson,Pi0tmp);
                      AliAODConversionMother *PiMiPiMiPiZerotmp = new AliAODConversionMother(&EventPiMiGoodMeson2,PiMiPiZerotmp);

                      // Fill histograms (likesign)
                      fHistoMotherLikeSignBackInvMassPt[fiCut]->Fill(PiMiPiMiPiZeroBackgroundCandidate->M(),PiMiPiMiPiZeroBackgroundCandidate->Pt());
                      fHistoMotherLikeSignBackInvMassSubPi0Pt[fiCut]->Fill(PiMiPiMiPiZeroBackgroundCandidate->M()-(EventPiZeroGoodMeson->M()-fPDGMassPi0),PiMiPiMiPiZeroBackgroundCandidate->Pt());
                      fHistoMotherLikeSignBackInvMassFixedPzPi0Pt[fiCut]->Fill(PiMiPiMiPiZerotmp->M(),PiMiPiMiPiZerotmp->Pt());

                      delete Pi0tmp;
                      delete PiMiPiZerotmp;
                      delete PiMiPiMiPiZerotmp;

                      delete PiMiPiMiPiZeroBackgroundCandidate;
                      PiMiPiMiPiZeroBackgroundCandidate = 0x0;
                    }
                  }


              }
          } // end of iCurrentPiMi2
     }// end of iCurrenPiMi
   } // end of LikeSign if
 } //end loop pi0 candidates
}

//______________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::KinematicCut(AliAODConversionMother *negpion, AliAODConversionMother *pospion, AliAODConversionMother *neutpion, AliAODConversionMother *omega){

  if(fTolerance == -1) return kTRUE;
  if((omega->Pt())<=5.){
    if( (omega->Angle(pospion->Vect()))    < ((2.78715*(TMath::Exp(-0.589934*(omega->Pt()))+0.0519574))*fTolerance) &&
        (omega->Angle(negpion->Vect()))    < ((5.94216*(TMath::Exp(-0.444428*(omega->Pt()))-0.0574076))*fTolerance) &&
        (omega->Angle(neutpion->Vect()))   < ((2.79529*(TMath::Exp(-0.565999*(omega->Pt()))+0.0413576))*fTolerance) &&
        (pospion->Angle(negpion->Vect()))  < ((3.14446*(TMath::Exp(-0.666433*(omega->Pt()))+0.0964309))*fTolerance) &&
        (pospion->Angle(neutpion->Vect())) < ((3.08241*(TMath::Exp(-0.650657*(omega->Pt()))+0.0997539))*fTolerance) &&
        (negpion->Angle(neutpion->Vect())) < ((3.18536*(TMath::Exp(-0.752847*(omega->Pt()))+0.1262780))*fTolerance)
      ){
        return kTRUE;
    }
  }else{
    if( (omega->Angle(pospion->Vect()))    < ((0.459270*(TMath::Exp(-0.126007*(omega->Pt()))+0.100475))*fTolerance) &&
        (omega->Angle(negpion->Vect()))    < ((0.521250*(TMath::Exp(-0.152532*(omega->Pt()))+0.114617))*fTolerance) &&
        (omega->Angle(neutpion->Vect()))   < ((0.409766*(TMath::Exp(-0.108566*(omega->Pt()))+0.103594))*fTolerance) &&
        (pospion->Angle(negpion->Vect()))  < ((0.709206*(TMath::Exp(-0.149072*(omega->Pt()))+0.111345))*fTolerance) &&
        (pospion->Angle(neutpion->Vect())) < ((0.662184*(TMath::Exp(-0.123397*(omega->Pt()))+0.104675))*fTolerance) &&
        (negpion->Angle(neutpion->Vect())) < ((0.730228*(TMath::Exp(-0.120859*(omega->Pt()))+0.105522))*fTolerance)
      ){
        return kTRUE;
    }
  }
  return kFALSE;
}


//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueMesonCandidates(AliAODConversionMother *mesoncand, AliAODConversionMother *TrueNeutralPionCandidate, AliAODConversionPhoton *TrueVirtualParticleCandidate)
{

  // Process True Mesons

  Bool_t isSameMotherPiPlPiMiPiZero   = kFALSE;   // pi+ pi- and pi0 have the same mother
  Bool_t isSameMotherPiPlPiMi         = kFALSE;   // pi+ and pi- have the same mother
  Bool_t isSameMotherPiPlPiZero       = kFALSE;   // pi+ and pi0 have the same mother
  Bool_t isSameMotherPiMiPiZero       = kFALSE;   // pi- and pi0 have the same mother
  Bool_t isNoSameMother               = kFALSE;   // none of the pions have the same mother
  Bool_t isNoPiPiPi                   = kFALSE;   // the decay is not a 3 pion decay


  Int_t virtualParticleMCLabel = TrueVirtualParticleCandidate->GetMCParticleLabel(fMCEvent);
  Int_t virtualParticleMotherLabel = -1;
  Int_t trueMesonFlag  = TrueNeutralPionCandidate->GetTrueMesonValue();
  Int_t pi0MCLabel     = TrueNeutralPionCandidate->GetMCLabel();

  Float_t weighted= 1;

  if ( !(trueMesonFlag == 1 && pi0MCLabel != -1)){
      if((fDoMesonQA>0 ) && (!fDoLightOutput)){
          fHistoTruePiPlPiMiPiZeroContaminationInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      }
      return;
  }
  Int_t pi0MotherLabel =  fMCEvent->Particle(pi0MCLabel)->GetMother(0);

  TParticle * negativeMC = (TParticle*)TrueVirtualParticleCandidate->GetNegativeMCDaughter(fMCEvent);
  TParticle * positiveMC = (TParticle*)TrueVirtualParticleCandidate->GetPositiveMCDaughter(fMCEvent);

  Int_t posMotherLabelMC = positiveMC->GetMother(0);
  Int_t negMotherLabelMC = negativeMC->GetMother(0);

  // Check case present
  if((TMath::Abs(negativeMC->GetPdgCode())==211) && (TMath::Abs(positiveMC->GetPdgCode())==211) && (fMCEvent->Particle(pi0MCLabel)->GetPdgCode()==111)){
    // three pion decay
    if(virtualParticleMCLabel!=-1){
      // pi+ pi- have same mother
      virtualParticleMotherLabel  = virtualParticleMCLabel;
      if(virtualParticleMotherLabel==pi0MotherLabel){
        // all pions from same mother
        if(fMCEvent->Particle(pi0MotherLabel)->GetStatusCode()!=21) isSameMotherPiPlPiMiPiZero  = kTRUE;
      } else{
        // only pi+ pi- from same mother
        if(fMCEvent->Particle(virtualParticleMotherLabel)->GetStatusCode()!=21) isSameMotherPiPlPiMi = kTRUE;
      }
    } else{
      if(pi0MotherLabel==negMotherLabelMC && negMotherLabelMC != -1){
        // pi0 and pi- same mother
        if(fMCEvent->Particle(negMotherLabelMC)->GetStatusCode()!=21) isSameMotherPiMiPiZero      = kTRUE;
      } else if(pi0MotherLabel==posMotherLabelMC && posMotherLabelMC != -1){
        // pi0 and pi+ same mother
        if(fMCEvent->Particle(posMotherLabelMC)->GetStatusCode()!=21)  isSameMotherPiPlPiZero      = kTRUE;
      } else{
        // all pions different mother
        isNoSameMother              = kTRUE;
      }
    }
  } else{
    // not a three pion decay
    isNoPiPiPi = kTRUE;
  }

  // Do things for each case
  if(isSameMotherPiPlPiMiPiZero){
    if(fMCEvent->Particle(pi0MotherLabel)->GetPdgCode()                        == 221){
      // eta was found
      fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      fHistoTrueMotherEtaPiPlPiMiPiZeroInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      AliAODConversionMother *PosPiontmp = new AliAODConversionMother();
      PosPiontmp->SetPxPyPzE(positiveMC->Px(), positiveMC->Py(), positiveMC->Pz(), positiveMC->Energy());
      AliAODConversionMother *NegPiontmp = new AliAODConversionMother();
      NegPiontmp->SetPxPyPzE(negativeMC->Px(), negativeMC->Py(), negativeMC->Pz(), negativeMC->Energy());
      if(!fDoLightOutput) fHistoTrueAngleSum[fiCut]->Fill(mesoncand->Pt(),((PosPiontmp->Angle(mesoncand->Vect()))+(NegPiontmp->Angle(PosPiontmp->Vect()))+(PosPiontmp->Angle(TrueNeutralPionCandidate->Vect()))));

      delete PosPiontmp; PosPiontmp = 0x0;
      delete NegPiontmp; NegPiontmp = 0x0;

      // Fill tree to get info about event that the eta was found in
      if(fDoMesonQA>1 && (!fDoLightOutput)){
         fV0MultiplicityEtaEvent = fMCEvent->GetNumberOfV0s();
         fTrackMultiplicityEtaEvent = fMCEvent->GetNumberOfTracks();
         fZVertexEtaEvent = fMCEvent->GetPrimaryVertex()->GetZ();
         fPtEta = mesoncand->Pt();

         fTreeEventInfoEta[fiCut]->Fill();
      }
      if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,pi0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt());
    } else if(fMCEvent->Particle(pi0MotherLabel)->GetPdgCode()                 == 223){
      // omega was found
      fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
      fHistoTrueMotherOmegaPiPlPiMiPiZeroInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);

      AliAODConversionMother *PosPiontmp = new AliAODConversionMother();
      PosPiontmp->SetPxPyPzE(positiveMC->Px(), positiveMC->Py(), positiveMC->Pz(), positiveMC->Energy());
      AliAODConversionMother *NegPiontmp = new AliAODConversionMother();
      NegPiontmp->SetPxPyPzE(negativeMC->Px(), negativeMC->Py(), negativeMC->Pz(), negativeMC->Energy());
      if(!fDoLightOutput) fHistoTrueAngleSum[fiCut]->Fill(mesoncand->Pt(),((PosPiontmp->Angle(mesoncand->Vect()))+(NegPiontmp->Angle(PosPiontmp->Vect()))+(PosPiontmp->Angle(TrueNeutralPionCandidate->Vect()))));

      // Fill tree to get info about event that the omega was found in
      if(fDoMesonQA>1 && (!fDoLightOutput)){
         fV0MultiplicityOmegaEvent = fMCEvent->GetNumberOfV0s();
         fTrackMultiplicityOmegaEvent = fMCEvent->GetNumberOfTracks();
         fZVertexOmegaEvent = fMCEvent->GetPrimaryVertex()->GetZ();
         fPtOmega = mesoncand->Pt();

         fTreeEventInfoOmega[fiCut]->Fill();
      }

      if (CheckVectorForDoubleCount(fVectorDoubleCountTrueOmegas,pi0MotherLabel) && (!fDoLightOutput)) fHistoDoubleCountTrueOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt());
    } else{
      if(fDoMesonQA>1 && (!fDoLightOutput)){
        // Write "unknown" mother to TTree
        fSamePiPiPiMotherID       = fMCEvent->Particle(posMotherLabelMC)->GetPdgCode();
        fSamePiPiPiMotherInvMass  = mesoncand->M();
        fSamePiPiPiMotherPt       = mesoncand->Pt();

        fTreePiPiPiSameMother[fiCut]->Fill();
      }
    }
  } else if(isSameMotherPiPlPiMi &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
    if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()                     == 221){
      // pi+pi- come from eta
      fHistoTruePiPlPiMiSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 223){
      // pi+pi- come from omega
      fHistoTruePiPlPiMiSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 113){
      // pi+pi- come from rho0
      fHistoTruePiPlPiMiSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 331){
      // pi+pi- come from eta prime
      fHistoTruePiPlPiMiSameMotherFromEtaPrimeInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 310){
      // pi+pi- come from K0 short
      fHistoTruePiPlPiMiSameMotherFromK0sInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 130){
      // pi+pi- come from K0 short
      fHistoTruePiPlPiMiSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else{
      // pi+pi- come from something else
      if(fDoMesonQA>1 && (!fDoLightOutput)){
        fCasePiPi = 0;
        // Write "unknown" mother to TTree
        fSamePiPiMotherID = fMCEvent->Particle(posMotherLabelMC)->GetPdgCode();
        fSamePiPiMotherInvMass = mesoncand->M();
        fSamePiPiMotherPt = mesoncand->Pt();

        fTreePiPiSameMother[fiCut]->Fill();
      }
    }
  } else if(isSameMotherPiMiPiZero  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
    if(fMCEvent->Particle(pi0MotherLabel)->GetPdgCode()                       == 221){
      // pi0pi- come from eta
      fHistoTruePiMiPiZeroSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(pi0MotherLabel)->GetPdgCode()                == 223){
      // pi0pi- come from omega
      fHistoTruePiMiPiZeroSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(pi0MotherLabel)->GetPdgCode()                ==-213){
      // pi0pi- come from rho-
      fHistoTruePiMiPiZeroSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(pi0MotherLabel)->GetPdgCode()                == 130){
      // pi0pi- come from rho-
      fHistoTruePiMiPiZeroSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else{
      // pi0pi- come from something else
      if(fDoMesonQA>1){
        fCasePiPi = 1;
        // Write "unknown" mother to TTree
        fSamePiPiMotherID = fMCEvent->Particle(pi0MotherLabel)->GetPdgCode();
        fSamePiPiMotherInvMass = mesoncand->M();
        fSamePiPiMotherPt = mesoncand->Pt();

        fTreePiPiSameMother[fiCut]->Fill();
      }
    }
  } else if(isSameMotherPiPlPiZero  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
    if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()                     == 221){
      // pi+pi0 come from eta
      fHistoTruePiPlPiZeroSameMotherFromEtaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 223){
      // pi+pi0 come from omega
      fHistoTruePiPlPiZeroSameMotherFromOmegaInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 213) {
      // pi+pi0 come from rho+
      fHistoTruePiPlPiZeroSameMotherFromRhoInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else if(fMCEvent->Particle(posMotherLabelMC)->GetPdgCode()              == 130) {
      // pi+pi0 come from rho+
      fHistoTruePiPlPiZeroSameMotherFromK0lInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    } else{
      // pi+pi0 come from something else
      if(fDoMesonQA>1){
        fCasePiPi = 2;
        // Write "unknown" mother to TTree
        fSamePiPiMotherID = fMCEvent->Particle(pi0MotherLabel)->GetPdgCode();
        fSamePiPiMotherInvMass = mesoncand->M();
        fSamePiPiMotherPt = mesoncand->Pt();

        fTreePiPiSameMother[fiCut]->Fill();
      }
    }
  } else if(isNoSameMother  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
    // no same mother purecombinatorical
    fHistoTruePiPlPiMiPiZeroPureCombinatoricalInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
  } else if(isNoPiPiPi  &&  (fDoMesonQA>0 ) && (!fDoLightOutput)){
    // no pi pi pi decay contamination
    fHistoTruePiPlPiMiPiZeroContaminationInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
    // investigate here what was missmatched (?)

  }
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::UpdateEventByEventData(){
  //see header file for documentation

  Int_t method = 1;
  if( method == 1 ) {
    if(fPosPionCandidates->GetEntries() >0 && fNegPionCandidates->GetEntries() >0){
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
        fBGHandlerPiPl[fiCut]->AddMesonEvent(fPosPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
        fBGHandlerPiMi[fiCut]->AddMesonEvent(fNegPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
      } else { // means we use #V0s for multiplicity
        if (fNeutralPionMode < 2){
          fBGHandlerPiPl[fiCut]->AddMesonEvent(fPosPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodConvGammas->GetEntries(),0);
          fBGHandlerPiMi[fiCut]->AddMesonEvent(fNegPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodConvGammas->GetEntries(),0);
        }else {
          fBGHandlerPiPl[fiCut]->AddMesonEvent(fPosPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),0);
          fBGHandlerPiMi[fiCut]->AddMesonEvent(fNegPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),0);
        }
      }
    }
  }
//	else if ( method == 2 ){
//		if(fGoodVirtualParticles->GetEntries() > 0 ){
//			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
//				fBGHandler[fiCut]->AddEvent(fGoodVirtualParticles,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
//			} else{ // means we use #V0s for multiplicity
//				fBGHandler[fiCut]->AddEvent(fGoodVirtualParticles,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodVirtualParticles->GetEntries(),0);
//			}
//		}
//	}
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::MoveParticleAccordingToVertex(AliAODConversionMother* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
  //see header file for documentation

  Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
  Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
  Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

  Double_t movedPlace[3] = {particle->GetProductionX() - dx,particle->GetProductionY() - dy,particle->GetProductionZ() - dz};
  particle->SetProductionPoint(movedPlace);
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::FixPzToMatchPDGInvMassPi0(AliAODConversionMother* particle){

Double_t px = particle->Px();
Double_t py = particle->Py();
Int_t signPz = particle->Pz()<0?-1:1;
Double_t energy = particle->Energy();
Double_t pz = signPz*TMath::Sqrt(TMath::Abs(pow(fPDGMassPi0,2)-pow(energy,2)+pow(px,2)+pow(py,2)));
particle->SetPxPyPzE(px,py,pz,energy);

return;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::IsEtaPiPlPiMiPiZeroDaughter( Int_t label ) const {
//
// Returns true if the particle comes from eta -> pi+ pi- gamma
//
    if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
// 	cout << "found eta? " << endl;
  if( mother->GetPdgCode() != 221 ) return kFALSE;
// 		else cout << "YES" << endl;
  if( IsPiPlPiMiPiZeroDecay( mother ) ) return kTRUE;
  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::IsOmegaPiPlPiMiPiZeroDaughter( Int_t label ) const {
//
// Returns true if the particle comes from eta -> pi+ pi- gamma
//
    if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
// 	cout << "found omega? " << endl;
  if( mother->GetPdgCode() != 223 ) return kFALSE;
// 		else cout << "YES" << endl;
  if( IsPiPlPiMiPiZeroDecay( mother ) ) return kTRUE;
  return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::IsPiPlPiMiPiZeroDecay(TParticle *fMCMother) const
{
// 	cout << fMCMother->GetNDaughters() << endl;
  if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
// 	cout << fMCMother->GetPdgCode() << endl;
  if( !(fMCMother->GetPdgCode() == 221 || fMCMother->GetPdgCode() == 223)  ) return kFALSE;
// 	cout << "made it til here" << endl;

  TParticle *posPion = 0x0;
  TParticle *negPion = 0x0;
  TParticle *neutPion    = 0x0;

  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
    if(index<0) continue;
    TParticle* temp = (TParticle*)fMCEvent->Particle( index );

    switch( temp->GetPdgCode() ) {
    case 211:
      posPion =  temp;
      break;
    case -211:
      negPion =  temp;
      break;
    case 111:
      neutPion = temp;
      break;
    }
  }
  if( posPion && negPion && neutPion) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::GammaIsNeutralMesonPiPlPiMiPiZeroDaughter( Int_t label ) const {
//
// Returns true if the particle comes from eta -> pi+ pi- gamma
//
    if(label<0) return kFALSE;
  Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
  if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;

  TParticle* mother = fMCEvent->Particle( motherLabel );
// 	cout << "found omega? " << endl;
  if( mother->GetPdgCode() != 111 ) return kFALSE;
// 		else cout << "YES" << endl;
  Int_t grandMotherLabel = mother->GetMother(0);
  if( grandMotherLabel < 0 || grandMotherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;
  TParticle* grandmother = fMCEvent->Particle( grandMotherLabel );

  if( IsPiPlPiMiPiZeroDecay( grandmother ) ) return kTRUE;
  return kFALSE;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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


