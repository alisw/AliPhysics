/****************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.   *
 *                                                                          *
 * Authors: Friederike Bock                                                 *
 * Version 1.0                                                              *
 *                                                                          *
 * Permission to use, copy, modify and distribute this software and its     *
 * documentation strictly for non-commercial purposes is hereby granted     *
 * without fee, provided that the above copyright notice appears in all     *
 * copies and that both the copyright notice and this permission notice     *
 * appear in the supporting documentation. The authors make no claims       *
 * about the suitability of this software for any purpose. It is            *
 * provided "as is" without express or implied warranty.                    *
 ***************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling photon selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliConversionPhotonCuts.h"

#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "AliStack.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliV0ReaderV1.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliTRDTriggerAnalysis.h"

class iostream;

using namespace std;

ClassImp(AliConversionPhotonCuts)


const char* AliConversionPhotonCuts::fgkCutNames[AliConversionPhotonCuts::kNCuts] = {
  "V0FinderType",           // 0
  "EtaCut",                 // 1
  "MinRCut",                // 2
  "EtaForPhiCut",           // 3
  "MinPhiCut",              // 4
  "MaxPhiCut",              // 5
  "SinglePtCut",            // 6
  "ClsTPCCut",              // 7
  "ededxSigmaCut",          // 8
  "pidedxSigmaCut",         // 9
  "piMomdedxSigmaCut",      // 10
  "piMaxMomdedxSigmaCut",   // 11
  "LowPRejectionSigmaCut",  // 12
  "TOFelectronPID",         // 13
  "ITSelectronPID",         // 14 -- new ITS PID
  "TRDelectronPID",         // 15 -- new TRD PID
  "QtMaxCut",               // 16
  "Chi2GammaCut",           // 17
  "PsiPair",                // 18
  "DoPhotonAsymmetryCut",   // 19
  "CosinePointingAngle",    // 20
  "SharedElectronCuts",     // 21
  "RejectToCloseV0s",       // 22
  "DcaRPrimVtx",            // 23
  "DcaZPrimVtx",            // 24
  "EvetPlane"               // 25
};


//________________________________________________________________________
AliConversionPhotonCuts::AliConversionPhotonCuts(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fDoLightOutput(kFALSE),
  fV0ReaderName("V0ReaderV1"),
  fMaxR(200),
  fMinR(0),
  fEtaCut(0.9),
  fEtaCutMin(-0.1),
  fEtaForPhiCutMin(-10.),
  fEtaForPhiCutMax(10.),
  fMinPhiCut(0.),
  fMaxPhiCut(100.),
  fDoShrinkTPCAcceptance(kFALSE),
  fPtCut(0.02),
  fSinglePtCut(0),
  fMaxZ(1000),
  fMinClsTPC(0.),
  fMinClsTPCToF(0.),
  fLineCutZRSlope(0.),
  fLineCutZValue(0),
  fLineCutZRSlopeMin(0.),
  fLineCutZValueMin(0),
  fChi2CutConversion(1000),
  fPIDProbabilityCutNegativeParticle(0),
  fPIDProbabilityCutPositiveParticle(0),
  fDodEdxSigmaCut(kTRUE),
  fDoTOFsigmaCut(kFALSE),
  fPIDTRDEfficiency(1),
  fDoTRDPID(kFALSE),
  fPIDnSigmaAboveElectronLine(100),
  fPIDnSigmaBelowElectronLine(-100),
  fTofPIDnSigmaAboveElectronLine(100),
  fTofPIDnSigmaBelowElectronLine(-100),
  fPIDnSigmaAbovePionLine(0),
  fPIDnSigmaAbovePionLineHighPt(-100),
  fPIDMinPnSigmaAbovePionLine(0),
  fPIDMaxPnSigmaAbovePionLine(0),
  fDoKaonRejectionLowP(kFALSE),
  fDoProtonRejectionLowP(kFALSE),
  fDoPionRejectionLowP(kFALSE),
  fPIDnSigmaAtLowPAroundKaonLine(0),
  fPIDnSigmaAtLowPAroundProtonLine(0),
  fPIDnSigmaAtLowPAroundPionLine(0),
  fPIDMinPKaonRejectionLowP(1.5),
  fPIDMinPProtonRejectionLowP(2),
  fPIDMinPPionRejectionLowP(0),	
  fDoQtGammaSelection(kTRUE),
  fDo2DQt(kFALSE),
  fQtMax(100),
  fNSigmaMass(0.),
  fUseEtaMinCut(kFALSE),
  fUseOnFlyV0Finder(kTRUE),
  fUseOnFlyV0FinderSameSign(0),
  fDoPhotonAsymmetryCut(kTRUE),
  fDoPhotonPDependentAsymCut(kFALSE),
  fFAsymmetryCut(0),
  fMinPPhotonAsymmetryCut(100.),
  fMinPhotonAsymmetry(0.),
  fUseCorrectedTPCClsInfo(kFALSE),
  fUseTOFpid(kFALSE),
  fOpeningAngle(0.005),
  fPsiPairCut(10000),
  fDo2DPsiPairChi2(kFALSE),
  fIncludeRejectedPsiPair(kFALSE),
  fCosPAngleCut(10000),
  fDoToCloseV0sCut(kFALSE),
  fminV0Dist(200.),
  fDoSharedElecCut(kFALSE),
  fDoPhotonQualitySelectionCut(kFALSE),
  fPhotonQualityCut(0),
  fRandom(0),
  fElectronArraySize(500),
  fElectronLabelArray(NULL),
  fDCAZPrimVtxCut(1000),
  fDCARPrimVtxCut(1000),
  fInPlaneOutOfPlane(0),
  fConversionPointXArray(0.0),
  fConversionPointYArray(0.0),
  fConversionPointZArray(0.0),
  fCutString(NULL),
  fIsHeavyIon(0),
  fUseITSpid(kFALSE),
  fITSPIDnSigmaAboveElectronLine(100),
  fITSPIDnSigmaBelowElectronLine(-100),
  fMaxPtPIDITS(1.5),
  fTRDPIDAboveCut(100),
  fTRDPIDBelowCut(-100),
  fDoDoubleCountingCut(kFALSE),
  fMinRDC(0.),
  fDeltaR(0.),
  fOpenAngle(0.),
  fSwitchToKappa(kFALSE),
  fKappaMinCut(-1),
  fKappaMaxCut(1000),
  fHistoEtaDistV0s(NULL),
  fHistoEtaDistV0sAfterdEdxCuts(NULL),
  fHistodEdxCuts(NULL),
  fHistoTPCdEdxbefore(NULL),
  fHistoTPCdEdxafter(NULL),
  fHistoTPCdEdxSigbefore(NULL),
  fHistoTPCdEdxSigafter(NULL),
  fHistoKappaafter(NULL),
  fHistoTOFbefore(NULL),
  fHistoTOFSigbefore(NULL),
  fHistoTOFSigafter(NULL),
  fHistoITSSigbefore(NULL),
  fHistoITSSigafter(NULL),
  fHistoPsiPairDeltaPhiafter(NULL),
  fHistoTrackCuts(NULL),
  fHistoPhotonCuts(NULL),
  fHistoInvMassbefore(NULL),
  fHistoArmenterosbefore(NULL),
  fHistoInvMassafter(NULL),
  fHistoArmenterosafter(NULL),
  fHistoAsymmetryafter(NULL),
  fHistoAcceptanceCuts(NULL),
  fHistoCutIndex(NULL),
  fHistoEventPlanePhi(NULL),
  fPreSelCut(kFALSE),
  fProcessAODCheck(kFALSE),
  fProfileContainingMaterialBudgetWeights(NULL),
  fMaterialBudgetWeightsInitialized(kFALSE)
{
  InitPIDResponse();
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());

  fElectronLabelArray = new Int_t[fElectronArraySize];
}

//________________________________________________________________________
AliConversionPhotonCuts::AliConversionPhotonCuts(const AliConversionPhotonCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fDoLightOutput(ref.fDoLightOutput),
  fV0ReaderName("V0ReaderV1"),
  fMaxR(ref.fMaxR),
  fMinR(ref.fMinR),
  fEtaCut(ref.fEtaCut),
  fEtaCutMin(ref.fEtaCutMin),
  fEtaForPhiCutMin(ref.fEtaForPhiCutMin),
  fEtaForPhiCutMax(ref.fEtaForPhiCutMax),
  fMinPhiCut(ref.fMinPhiCut),
  fMaxPhiCut(ref.fMaxPhiCut),
  fDoShrinkTPCAcceptance(ref.fDoShrinkTPCAcceptance),
  fPtCut(ref.fPtCut),
  fSinglePtCut(ref.fSinglePtCut),
  fMaxZ(ref.fMaxZ),
  fMinClsTPC(ref.fMinClsTPC),
  fMinClsTPCToF(ref.fMinClsTPCToF),
  fLineCutZRSlope(ref.fLineCutZRSlope),
  fLineCutZValue(ref.fLineCutZValue),
  fLineCutZRSlopeMin(ref.fLineCutZRSlopeMin),
  fLineCutZValueMin(ref.fLineCutZValueMin),
  fChi2CutConversion(ref.fChi2CutConversion),
  fPIDProbabilityCutNegativeParticle(ref.fPIDProbabilityCutNegativeParticle),
  fPIDProbabilityCutPositiveParticle(ref.fPIDProbabilityCutPositiveParticle),
  fDodEdxSigmaCut(ref. fDodEdxSigmaCut),
  fDoTOFsigmaCut(ref.fDoTOFsigmaCut),
  fPIDTRDEfficiency(ref.fPIDTRDEfficiency),
  fDoTRDPID(ref.fDoTRDPID),
  fPIDnSigmaAboveElectronLine(ref.fPIDnSigmaAboveElectronLine),
  fPIDnSigmaBelowElectronLine(ref.fPIDnSigmaBelowElectronLine),
  fTofPIDnSigmaAboveElectronLine(ref.fTofPIDnSigmaAboveElectronLine),
  fTofPIDnSigmaBelowElectronLine(ref.fTofPIDnSigmaBelowElectronLine),
  fPIDnSigmaAbovePionLine(ref.fPIDnSigmaAbovePionLine),
  fPIDnSigmaAbovePionLineHighPt(ref.fPIDnSigmaAbovePionLineHighPt),
  fPIDMinPnSigmaAbovePionLine(ref.fPIDMinPnSigmaAbovePionLine),
  fPIDMaxPnSigmaAbovePionLine(ref.fPIDMaxPnSigmaAbovePionLine),
  fDoKaonRejectionLowP(ref.fDoKaonRejectionLowP),
  fDoProtonRejectionLowP(ref.fDoProtonRejectionLowP),
  fDoPionRejectionLowP(ref.fDoPionRejectionLowP),
  fPIDnSigmaAtLowPAroundKaonLine(ref.fPIDnSigmaAtLowPAroundKaonLine),
  fPIDnSigmaAtLowPAroundProtonLine(ref.fPIDnSigmaAtLowPAroundProtonLine),
  fPIDnSigmaAtLowPAroundPionLine(ref.fPIDnSigmaAtLowPAroundPionLine),
  fPIDMinPKaonRejectionLowP(ref.fPIDMinPKaonRejectionLowP),
  fPIDMinPProtonRejectionLowP(ref.fPIDMinPProtonRejectionLowP),
  fPIDMinPPionRejectionLowP(ref.fPIDMinPPionRejectionLowP),
  fDoQtGammaSelection(ref.fDoQtGammaSelection),
  fDo2DQt(ref.fDo2DQt),
  fQtMax(ref.fQtMax),
  fNSigmaMass(ref.fNSigmaMass),
  fUseEtaMinCut(ref.fUseEtaMinCut),
  fUseOnFlyV0Finder(ref.fUseOnFlyV0Finder),
  fUseOnFlyV0FinderSameSign(ref.fUseOnFlyV0FinderSameSign),
  fDoPhotonAsymmetryCut(ref.fDoPhotonAsymmetryCut),
  fDoPhotonPDependentAsymCut(ref.fDoPhotonPDependentAsymCut),
  fFAsymmetryCut(ref.fFAsymmetryCut),
  fMinPPhotonAsymmetryCut(ref.fMinPPhotonAsymmetryCut),
  fMinPhotonAsymmetry(ref.fMinPhotonAsymmetry),
  fUseCorrectedTPCClsInfo(ref.fUseCorrectedTPCClsInfo),
  fUseTOFpid(ref.fUseTOFpid),
  fOpeningAngle(ref.fOpeningAngle),
  fPsiPairCut(ref.fPsiPairCut),
  fDo2DPsiPairChi2(ref.fDo2DPsiPairChi2),
  fIncludeRejectedPsiPair(ref.fIncludeRejectedPsiPair),
  fCosPAngleCut(ref.fCosPAngleCut),
  fDoToCloseV0sCut(ref.fDoToCloseV0sCut),
  fminV0Dist(ref.fminV0Dist),
  fDoSharedElecCut(ref.fDoSharedElecCut),
  fDoPhotonQualitySelectionCut(ref.fDoPhotonQualitySelectionCut),
  fPhotonQualityCut(ref.fPhotonQualityCut),
  fRandom(ref.fRandom),
  fElectronArraySize(ref.fElectronArraySize),
  fElectronLabelArray(NULL),
  fDCAZPrimVtxCut(ref.fDCAZPrimVtxCut),
  fDCARPrimVtxCut(ref.fDCAZPrimVtxCut),
  fInPlaneOutOfPlane(ref.fInPlaneOutOfPlane),
  fConversionPointXArray(ref.fConversionPointXArray),
  fConversionPointYArray(ref.fConversionPointYArray),
  fConversionPointZArray(ref.fConversionPointZArray),
  fCutString(NULL),
  fIsHeavyIon(ref.fIsHeavyIon),
  fUseITSpid(ref.fUseITSpid),
  fITSPIDnSigmaAboveElectronLine(ref.fITSPIDnSigmaAboveElectronLine),
  fITSPIDnSigmaBelowElectronLine(ref.fITSPIDnSigmaBelowElectronLine),
  fMaxPtPIDITS(ref.fMaxPtPIDITS),  
  fTRDPIDAboveCut(ref.fTRDPIDAboveCut),
  fTRDPIDBelowCut(ref.fTRDPIDBelowCut),
  fDoDoubleCountingCut(ref.fDoDoubleCountingCut),
  fMinRDC(ref.fMinRDC),
  fDeltaR(ref.fDeltaR),
  fOpenAngle(ref.fOpenAngle),
  fSwitchToKappa(ref.fSwitchToKappa),
  fKappaMinCut(ref.fKappaMinCut),
  fKappaMaxCut(ref.fKappaMaxCut),
  fHistoEtaDistV0s(NULL),
  fHistoEtaDistV0sAfterdEdxCuts(NULL),
  fHistodEdxCuts(NULL),
  fHistoTPCdEdxbefore(NULL),
  fHistoTPCdEdxafter(NULL),
  fHistoTPCdEdxSigbefore(NULL),
  fHistoTPCdEdxSigafter(NULL),
  fHistoKappaafter(NULL),
  fHistoTOFbefore(NULL),
  fHistoTOFSigbefore(NULL),
  fHistoTOFSigafter(NULL),
  fHistoITSSigbefore(NULL),
  fHistoITSSigafter(NULL),
  fHistoPsiPairDeltaPhiafter(NULL),
  fHistoTrackCuts(NULL),
  fHistoPhotonCuts(NULL),
  fHistoInvMassbefore(NULL),
  fHistoArmenterosbefore(NULL),
  fHistoInvMassafter(NULL),
  fHistoArmenterosafter(NULL),
  fHistoAsymmetryafter(NULL),
  fHistoAcceptanceCuts(NULL),
  fHistoCutIndex(NULL),
  fHistoEventPlanePhi(NULL),
  fPreSelCut(ref.fPreSelCut),
  fProcessAODCheck(ref.fProcessAODCheck),
  fProfileContainingMaterialBudgetWeights(ref.fProfileContainingMaterialBudgetWeights),
  fMaterialBudgetWeightsInitialized(ref.fMaterialBudgetWeightsInitialized)
{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());
  fElectronLabelArray = new Int_t[fElectronArraySize];
  // dont copy histograms (if you like histograms, call InitCutHistograms())
}


//________________________________________________________________________
AliConversionPhotonCuts::~AliConversionPhotonCuts() {
  // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  // if(fHistograms)
  //    delete fHistograms;
  // fHistograms = NULL;
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }
  if(fElectronLabelArray){
    delete fElectronLabelArray;
    fElectronLabelArray = NULL;
  }
  
  if(fFAsymmetryCut != NULL){
    delete fFAsymmetryCut;
    fFAsymmetryCut = NULL;
  }
  if(fProfileContainingMaterialBudgetWeights){
      delete fProfileContainingMaterialBudgetWeights;
      fProfileContainingMaterialBudgetWeights = 0x0;
  }  
}

//________________________________________________________________________
void AliConversionPhotonCuts::InitCutHistograms(TString name, Bool_t preCut){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);

  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }
  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("ConvCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  // IsPhotonSelected
  fHistoCutIndex=new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",12,-0.5,11.5);
  fHistoCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kOnFly+1,"onfly");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kNoV0+1,"miss. V0 in AOD");
  if (!fSwitchToKappa)fHistoCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"PID");
  else fHistoCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"Kappa+[TOF,ITS,TRD] PID");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kTrackCuts+1,"Track cuts");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kConvPointFail+1,"ConvPoint fail");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kPhotonCuts+1,"PhotonCuts");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kEventPlane+1,"EventPlane");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
  fHistograms->Add(fHistoCutIndex);

  // Track Cuts
  fHistoTrackCuts=new TH1F(Form("TrackCuts %s",GetCutNumber().Data()),"TrackCuts",9,-0.5,8.5);
  fHistoTrackCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(2,"likesign");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(3,"ntpccl");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(4,"acceptance");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(5,"singlept");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(6,"TPCrefit");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(7,"kink");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(8,"out");
  fHistograms->Add(fHistoTrackCuts);

  // Photon Cuts
  fHistoPhotonCuts=new TH2F(Form("PhotonCuts %s",GetCutNumber().Data()),"PhotonCuts vs p_{T,#gamma}",15,-0.5,14.5,250,0,50);
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(2,"qtcut");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(3,"chi2");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(4,"acceptance");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(5,"asymmetry");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(6,"pidprob");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(7,"cortpcclinfo");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(8,"PsiPair");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(9,"CosPAngle");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(10,"DCA R");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(11,"DCA Z");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(12,"Photon Quality");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(13,"out");
  fHistograms->Add(fHistoPhotonCuts);
  
  if (fProfileContainingMaterialBudgetWeights){
      fProfileContainingMaterialBudgetWeights->SetName("InputMaterialBudgetWeightsPerGamma");
      fHistograms->Add(fProfileContainingMaterialBudgetWeights);
  }

  if(!fDoLightOutput){
    if(preCut){
      fHistoInvMassbefore=new TH1F(Form("InvMass_before %s",GetCutNumber().Data()),"InvMass_before",1000,0,0.3);
      fHistograms->Add(fHistoInvMassbefore);
      fHistoArmenterosbefore=new TH2F(Form("Armenteros_before %s",GetCutNumber().Data()),"Armenteros_before",200,-1,1,1000,0,1.);
      fHistograms->Add(fHistoArmenterosbefore);
      fHistoEtaDistV0s = new TH1F(Form("Eta_before %s",GetCutNumber().Data()),"Eta_before",2000,-2,2);
      fHistograms->Add(fHistoEtaDistV0s);

    }
    fHistoInvMassafter=new TH1F(Form("InvMass_after %s",GetCutNumber().Data()),"InvMass_after",1000,0,0.3);
    fHistograms->Add(fHistoInvMassafter);
    fHistoArmenterosafter=new TH2F(Form("Armenteros_after %s",GetCutNumber().Data()),"Armenteros_after",200,-1,1,250,0,0.25);
    fHistograms->Add(fHistoArmenterosafter);
    if(fDoPhotonAsymmetryCut){
      fHistoAsymmetryafter=new TH2F(Form("Asymmetry_after %s",GetCutNumber().Data()),"Asymmetry_after",150,0.03,20.,200,0,1.);
      fHistograms->Add(fHistoAsymmetryafter);
    }
  }

  fHistoAcceptanceCuts=new TH2F(Form("PhotonAcceptanceCuts %s",GetCutNumber().Data()),"PhotonAcceptanceCuts vs p_{T,#gamma}",11,-0.5,10.5,250,0,50);
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(2,"maxR");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(3,"minR");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(4,"line");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(5,"maxZ");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(6,"eta");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(7,"phisector");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(8,"minpt");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(9,"out");
  fHistograms->Add(fHistoAcceptanceCuts);

  // dEdx Cuts 
  fHistodEdxCuts=new TH2F(Form("dEdxCuts %s",GetCutNumber().Data()),"dEdxCuts vs p_{T,e}",11,-0.5,10.5,250,0,50);
  fHistodEdxCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(2,"TPCelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(3,"TPCpion");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(4,"TPCpionhighp");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(5,"TPCkaonlowprej");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(6,"TPCprotonlowprej");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(7,"TPCpionlowprej");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(8,"TOFelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(9,"ITSelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(10,"TRDelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(11,"out");
  fHistograms->Add(fHistodEdxCuts);
  
  if(!fDoLightOutput){
    TAxis *AxisBeforedEdx = NULL;
    TAxis *AxisBeforedEdxSig = NULL;
    TAxis *AxisBeforeTOF = NULL;
    TAxis *AxisBeforeTOFSig = NULL;
    TAxis *AxisBeforeITSSig = NULL;
    if(preCut){
      fHistoTPCdEdxbefore=new TH2F(Form("Gamma_dEdx_before %s",GetCutNumber().Data()),"dEdx Gamma before" ,150,0.03,20,800,0,200);
      fHistograms->Add(fHistoTPCdEdxbefore);
      AxisBeforedEdx = fHistoTPCdEdxbefore->GetXaxis();
      fHistoTPCdEdxSigbefore=new TH2F(Form("Gamma_dEdxSig_before %s",GetCutNumber().Data()),"dEdx Sigma Gamma before" ,150,0.03,20,400,-10,10);
      fHistograms->Add(fHistoTPCdEdxSigbefore);
      AxisBeforedEdxSig = fHistoTPCdEdxSigbefore->GetXaxis();

      fHistoTOFbefore=new TH2F(Form("Gamma_TOF_before %s",GetCutNumber().Data()),"TOF Gamma before" ,150,0.03,20,11000,-1000,10000);
      fHistograms->Add(fHistoTOFbefore);
      AxisBeforeTOF = fHistoTOFbefore->GetXaxis();
      fHistoTOFSigbefore=new TH2F(Form("Gamma_TOFSig_before %s",GetCutNumber().Data()),"TOF Sigma Gamma before" ,150,0.03,20,400,-6,10);
      fHistograms->Add(fHistoTOFSigbefore);
      AxisBeforeTOFSig = fHistoTOFSigbefore->GetXaxis();

      fHistoITSSigbefore=new TH2F(Form("Gamma_ITSSig_before %s",GetCutNumber().Data()),"ITS Sigma Gamma before" ,150,0.03,20,400,-10,10);
      fHistograms->Add(fHistoITSSigbefore);
      AxisBeforeITSSig = fHistoITSSigbefore->GetXaxis();
    }

    fHistoTPCdEdxSigafter=new TH2F(Form("Gamma_dEdxSig_after %s",GetCutNumber().Data()),"dEdx Sigma Gamma after" ,150,0.03,20,400, -10,10);
    fHistograms->Add(fHistoTPCdEdxSigafter);

    fHistoTPCdEdxafter=new TH2F(Form("Gamma_dEdx_after %s",GetCutNumber().Data()),"dEdx Gamma after" ,150,0.03,20,800,0,200);
    fHistograms->Add(fHistoTPCdEdxafter);

    fHistoKappaafter=new TH2F(Form("Gamma_Kappa_after %s",GetCutNumber().Data()),"Kappa Gamma after" ,150,0.03,20,200,-20,20);
    fHistograms->Add(fHistoKappaafter);

    fHistoTOFSigafter=new TH2F(Form("Gamma_TOFSig_after %s",GetCutNumber().Data()),"TOF Sigma Gamma after" ,150,0.03,20,400,-6,10);
    fHistograms->Add(fHistoTOFSigafter);

    fHistoITSSigafter=new TH2F(Form("Gamma_ITSSig_after %s",GetCutNumber().Data()),"ITS Sigma Gamma after" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoITSSigafter);

    fHistoEtaDistV0sAfterdEdxCuts = new TH1F(Form("Eta_afterdEdx %s",GetCutNumber().Data()),"Eta_afterdEdx",2000,-2,2);
    fHistograms->Add(fHistoEtaDistV0sAfterdEdxCuts);

    fHistoPsiPairDeltaPhiafter=new TH2F(Form("Gamma_PsiPairDeltaPhi_after %s",GetCutNumber().Data()),"Psi Pair vs Delta Phi Gamma after" ,200,-2,2,200,-2,2);
    fHistograms->Add(fHistoPsiPairDeltaPhiafter);

    TAxis *AxisAfter = fHistoTPCdEdxSigafter->GetXaxis();
    Int_t bins = AxisAfter->GetNbins();
    Double_t from = AxisAfter->GetXmin();
    Double_t to = AxisAfter->GetXmax();
    Double_t *newBins = new Double_t[bins+1];
    newBins[0] = from;
    Double_t factor = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoTOFSigafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoTPCdEdxafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoKappaafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoITSSigafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    if(fDoPhotonAsymmetryCut){
      AxisAfter = fHistoAsymmetryafter->GetXaxis();
      AxisAfter->Set(bins, newBins);
    }
    if(preCut){
      AxisBeforedEdx->Set(bins, newBins);
      AxisBeforeTOF->Set(bins, newBins);
      AxisBeforedEdxSig->Set(bins, newBins);
      AxisBeforeTOFSig->Set(bins, newBins);
      AxisBeforeITSSig->Set(bins, newBins);
    }
    delete [] newBins;

    // Event Cuts and Info
    if(!preCut){
      fHistoEventPlanePhi=new TH1F(Form("EventPlaneMinusPhotonAngle %s",GetCutNumber().Data()),"EventPlaneMinusPhotonAngle",360,-TMath::Pi(),TMath::Pi());
      fHistograms->Add(fHistoEventPlanePhi);
    }
  }

  TH1::AddDirectory(kTRUE);
}

//________________________________________________________________________
Bool_t AliConversionPhotonCuts::InitPIDResponse(){
  // Set Pointer to AliPIDResponse

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if(fPIDResponse)return kTRUE;

  }


  return kFALSE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonIsSelectedMC(TParticle *particle,AliStack *fMCStack,Bool_t checkForConvertedGamma){
  // MonteCarlo Photon Selection

  if(!fMCStack)return kFALSE;

  if (particle->GetPdgCode() == 22){


    if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) )
      return kFALSE;
    if(fEtaCutMin>-0.1){
      if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) )
        return kFALSE;
    }

    if(particle->GetMother(0) >-1 && fMCStack->Particle(particle->GetMother(0))->GetPdgCode() == 22){
      return kFALSE; // no photon as mothers!
    }
    
    // removed, decision on primary and secondary taken in main task
// 		if(particle->GetMother(0) >= fMCStack->GetNprimary()){
// 			return kFALSE; // the gamma has a mother, and it is not a primary particle
// 		}

    if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

    // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
    TParticle* ePos = NULL;
    TParticle* eNeg = NULL;

    if(particle->GetNDaughters() >= 2){
      for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
        if(daughterIndex<0) continue;
        TParticle *tmpDaughter = fMCStack->Particle(daughterIndex);
        if(tmpDaughter->GetUniqueID() == 5){
        if(tmpDaughter->GetPdgCode() == 11){
          eNeg = tmpDaughter;
        } else if(tmpDaughter->GetPdgCode() == -11){
          ePos = tmpDaughter;
        }
        }
      }
    }

    if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
      return kFALSE;
    }

    if(ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut){
      return kFALSE; // no reconstruction below the Pt cut
    }

    if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ||
      eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) )
      return kFALSE;

    if(fEtaCutMin > -0.1){
      if( (ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin)) ||
        (eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin)) )
        return kFALSE;
    }

    if(ePos->R()>fMaxR){
      return kFALSE; // cuts on distance from collision point
    }

    if(TMath::Abs(ePos->Vz()) > fMaxZ){
      return kFALSE;  // outside material
    }
    if(TMath::Abs(eNeg->Vz()) > fMaxZ){
      return kFALSE;  // outside material
    }

    if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE;  // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    if( eNeg->R() <= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE; // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   eNeg->R() >= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    return kTRUE;
    //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
  }
  return kFALSE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray,Bool_t checkForConvertedGamma){
  // MonteCarlo Photon Selection

  if(!aodmcArray)return kFALSE;

  if (particle->GetPdgCode() == 22){
    if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) )
      return kFALSE;
    if(fEtaCutMin>-0.1){
      if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) )
        return kFALSE;
    }

    if(particle->GetMother() > -1 && (static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 22){
        return kFALSE; // no photon as mothers!
    }
      // removed, decision on primary and secondary taken in main task
//			Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
//			if(!isPrimary){
//				return kFALSE; // the gamma has a mother, and it is not a primary particle
//			}

    if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

    // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
    AliAODMCParticle* ePos = NULL;
    AliAODMCParticle* eNeg = NULL;

    if(particle->GetNDaughters() >= 2){
      for(Int_t daughterIndex=particle->GetDaughter(0);daughterIndex<=particle->GetDaughter(1);daughterIndex++){
        AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(aodmcArray->At(daughterIndex));
        if(!tmpDaughter) continue;
        if(((tmpDaughter->GetMCProcessCode())) == 5){    // STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX
        if(tmpDaughter->GetPdgCode() == 11){
          eNeg = tmpDaughter;
        } else if(tmpDaughter->GetPdgCode() == -11){
          ePos = tmpDaughter;
        }
        }
      }
    }

    if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
      return kFALSE;
    }

    if(ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut){
      return kFALSE; // no reconstruction below the Pt cut
    }

    if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ||
      eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) )
      return kFALSE;

    if(fEtaCutMin > -0.1){
      if( (ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin)) ||
        (eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin)) )
        return kFALSE;
    }

    Double_t rPos = sqrt( (ePos->Xv()*ePos->Xv()) + (ePos->Yv()*ePos->Yv()) );
    Double_t rNeg = sqrt( (eNeg->Xv()*eNeg->Xv()) + (eNeg->Yv()*eNeg->Yv()) );

    if(rPos>fMaxR){
      return kFALSE; // cuts on distance from collision point
    }
    if(TMath::Abs(ePos->Zv()) > fMaxZ){
      return kFALSE;  // outside material
    }
    if(TMath::Abs(eNeg->Zv()) > fMaxZ){
      return kFALSE;  // outside material
    }

    if( rPos <= ((TMath::Abs(ePos->Zv()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE;  // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   rPos >= ((TMath::Abs(ePos->Zv()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    if( rNeg <= ((TMath::Abs(eNeg->Zv()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE; // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   rNeg >= ((TMath::Abs(eNeg->Zv()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    return kTRUE;
    //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
  }
  return kFALSE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event){   // Specific Photon Cuts

  Int_t cutIndex = 0;
  if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt());
  cutIndex++;

  // Fill Histos before Cuts
  if(fHistoInvMassbefore)fHistoInvMassbefore->Fill(photon->GetMass());
  if(fHistoArmenterosbefore)fHistoArmenterosbefore->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());

  // Gamma selection based on QT from Armenteros
  if(fDoQtGammaSelection == kTRUE){
    if(!ArmenterosQtCut(photon)){
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //1
      return kFALSE;
    }
  }
  cutIndex++; //2

  // Chi Cut
  if(photon->GetChi2perNDF() > fChi2CutConversion || photon->GetChi2perNDF() <=0){
    {
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //2
      return kFALSE;
    }
  }
  cutIndex++;//3

  // Reconstruction Acceptance Cuts
  if(!AcceptanceCuts(photon)){
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //3
    return kFALSE;
  }

  cutIndex++; //4
  // Asymmetry Cut
  if(fDoPhotonAsymmetryCut == kTRUE){
    if(!AsymmetryCut(photon,event)){
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //4
      return kFALSE;
    }
  }

  //Check the pid probability
  cutIndex++; //5
  if(!PIDProbabilityCut(photon, event)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //5
    return kFALSE;
  }

  cutIndex++; //6
  if(!CorrectedTPCClusterCut(photon, event)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //6
    return kFALSE;
  }

  Double_t magField = event->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }
  
  AliVTrack * electronCandidate = GetTrack(event,photon->GetTrackLabelNegative());
  AliVTrack * positronCandidate = GetTrack(event,photon->GetTrackLabelPositive());
  Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

  cutIndex++; //7
  if(!PsiPairCut(photon)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //7
    return kFALSE;
  }

  cutIndex++; //8
  if(!CosinePAngleCut(photon, event)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //8
    return kFALSE;
  }

  AliAODConversionPhoton* photonAOD = dynamic_cast<AliAODConversionPhoton*>(photon);
  if (photonAOD){
    photonAOD->CalculateDistanceOfClossetApproachToPrimVtx(event->GetPrimaryVertex());

    cutIndex++; //9
    if(photonAOD->GetDCArToPrimVtx() > fDCARPrimVtxCut) { //DCA R cut of photon to primary vertex
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //9
      return kFALSE;
    }

    cutIndex++; //10
    if(TMath::Abs(photonAOD->GetDCAzToPrimVtx()) > fDCAZPrimVtxCut) { //DCA Z cut of photon to primary vertex
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //10
      return kFALSE;
    }
  } else {
    cutIndex++; //9
    cutIndex++; //10
  }
  cutIndex++; //11

  if (photonAOD){
    UChar_t photonQuality = 0;
      AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*>(event);
      if(aodEvent) {
        photonQuality = DeterminePhotonQualityAOD(photonAOD, event);
      } else {
        photonQuality = photonAOD->GetPhotonQuality();
      }	
      if (fDoPhotonQualitySelectionCut && photonQuality != fPhotonQualityCut){
        if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //11
        return kFALSE;
      }	
  } 
  cutIndex++; //12
  if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //11

  // Histos after Cuts
  if(fHistoInvMassafter)fHistoInvMassafter->Fill(photon->GetMass());
  if(fHistoArmenterosafter)fHistoArmenterosafter->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());
  if(fHistoPsiPairDeltaPhiafter)fHistoPsiPairDeltaPhiafter->Fill(deltaPhi,photon->GetPsiPair());
  if(fHistoKappaafter)fHistoKappaafter->Fill(photon->GetPhotonPt(), GetKappaTPC(photon, event));
  if(fHistoAsymmetryafter){
    if(photon->GetPhotonP()!=0 && electronCandidate->P()!=0)fHistoAsymmetryafter->Fill(photon->GetPhotonP(),electronCandidate->P()/photon->GetPhotonP());
  }
  return kTRUE;

}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event){   //Cut on corrected TPC Cluster Info

  AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

  if(!negTrack||!posTrack)return kFALSE;

  Double_t negclsToF=0;

  if (!fUseCorrectedTPCClsInfo ){
    if(negTrack->GetTPCNclsF()!=0){
      negclsToF = (Double_t)negTrack->GetNcls(1)/(Double_t)negTrack->GetTPCNclsF();}// Ncluster/Nfindablecluster
  }
  else {
    negclsToF = negTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
  }

  Double_t posclsToF = 0.;
  if (!fUseCorrectedTPCClsInfo ){
    if(posTrack->GetTPCNclsF()!=0){
      posclsToF = (Double_t)posTrack->GetNcls(1)/(Double_t)posTrack->GetTPCNclsF();
    }
  }else{
    posclsToF = posTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
  }

  if( negclsToF < fMinClsTPCToF || posclsToF < fMinClsTPCToF ){
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonIsSelected(AliConversionPhotonBase *photon, AliVEvent * event){
  //Selection of Reconstructed Photons

  FillPhotonCutIndex(kPhotonIn);

  if(event->IsA()==AliESDEvent::Class()) {
    if(!SelectV0Finder( ( ((AliESDEvent*)event)->GetV0(photon->GetV0Index()))->GetOnFlyStatus() ) ){
      FillPhotonCutIndex(kOnFly);
      return kFALSE;
    }
  }

  // Get Tracks
  AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

  if(!negTrack || !posTrack) {
    FillPhotonCutIndex(kNoTracks);
    return kFALSE;
  }

  // check if V0 from AliAODGammaConversion.root is actually contained in AOD by checking if V0 exists with same tracks
  if(event->IsA()==AliAODEvent::Class() && fPreSelCut && ( fIsHeavyIon != 1 || (fIsHeavyIon == 1 && fProcessAODCheck) )) {
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);

    Bool_t bFound = kFALSE;
    Int_t v0PosID = posTrack->GetID();
    Int_t v0NegID = negTrack->GetID();
    AliAODv0* v0 = NULL;
    for(Int_t iV=0; iV<aodEvent->GetNumberOfV0s(); iV++){
      v0 = aodEvent->GetV0(iV);
      if(!v0) continue;
      if( (v0PosID == v0->GetPosID() && v0NegID == v0->GetNegID()) || (v0PosID == v0->GetNegID() && v0NegID == v0->GetPosID()) ){
        bFound = kTRUE;
        break;
      }
    }
    if(!bFound){
      FillPhotonCutIndex(kNoV0);
      return kFALSE;
    }
  }

  photon->DeterminePhotonQuality(negTrack,posTrack);
  // Track Cuts
  if(!TracksAreSelected(negTrack, posTrack)){
    FillPhotonCutIndex(kTrackCuts);
    return kFALSE;
  }
  if (fHistoEtaDistV0s)fHistoEtaDistV0s->Fill(photon->GetPhotonEta());
  // dEdx Cuts
  
  if(!KappaCuts(photon, event) || !dEdxCuts(negTrack) || !dEdxCuts(posTrack)) {
    FillPhotonCutIndex(kdEdxCuts);
    return kFALSE;
  }
    
  if (fHistoEtaDistV0sAfterdEdxCuts)fHistoEtaDistV0sAfterdEdxCuts->Fill(photon->GetPhotonEta());
  // Photon Cuts
  if(!PhotonCuts(photon,event)){
    FillPhotonCutIndex(kPhotonCuts);
    return kFALSE;
  }
  
  // Photon passed cuts
  FillPhotonCutIndex(kPhotonOut);
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::ArmenterosQtCut(AliConversionPhotonBase *photon){   // Armenteros Qt Cut
  if(fDo2DQt){
    if ( !(TMath::Power(photon->GetArmenterosAlpha()/0.95,2)+TMath::Power(photon->GetArmenterosQt()/fQtMax,2) < 1) ){
      return kFALSE;
    }
  } else {
    if(photon->GetArmenterosQt()>fQtMax){
      return kFALSE;
    }
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::AcceptanceCuts(AliConversionPhotonBase *photon) {
  // Exclude certain areas for photon reconstruction

  Int_t cutIndex=0;
  if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
  cutIndex++;

  if(photon->GetConversionRadius()>fMaxR){ // cuts on distance from collision point
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(photon->GetConversionRadius()<fMinR){ // cuts on distance from collision point
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(photon->GetConversionRadius() <= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlope)-fLineCutZValue)){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  else if (fUseEtaMinCut &&  photon->GetConversionRadius() >= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(TMath::Abs(photon->GetConversionZ()) > fMaxZ ){ // cuts out regions where we do not reconstruct
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;


  if( photon->GetPhotonEta() > (fEtaCut)    || photon->GetPhotonEta() < (-fEtaCut) ){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  if(fEtaCutMin>-0.1){
    if( photon->GetPhotonEta() < (fEtaCutMin) && photon->GetPhotonEta() > (-fEtaCutMin) ){
      if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
      return kFALSE;
    }
  }
  cutIndex++;
  
  if (fDoShrinkTPCAcceptance){
    if(photon->GetPhotonEta() > fEtaForPhiCutMin && photon->GetPhotonEta() < fEtaForPhiCutMax ){
      if (fMinPhiCut < fMaxPhiCut){
        if( photon->GetPhotonPhi() > fMinPhiCut && photon->GetPhotonPhi() < fMaxPhiCut ) {
          if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
          return kFALSE;
        }
      } else {
        Double_t photonPhi = photon->GetPhotonPhi();
        if (photon->GetPhotonPhi() < TMath::Pi()) photonPhi = photon->GetPhotonPhi() + 2*TMath::Pi();
        if( photonPhi > fMinPhiCut && photonPhi < fMaxPhiCut+2*TMath::Pi() ) {
          if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
          return kFALSE;
        }	
      }	
    }
  }	
  cutIndex++;
  

  
  if(photon->GetPhotonPt()<fPtCut){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());

  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex) {
  // Track Cuts which require AOD/ESD specific implementation

  if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  AliAODVertex * NegVtxType=negTrack->GetProdVertex();
  AliAODVertex * PosVtxType=posTrack->GetProdVertex();
  if( (NegVtxType->GetType())==AliAODVertex::kKink || (PosVtxType->GetType())==AliAODVertex::kKink) {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  return kTRUE;

}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex) {
  // Track Cuts which require AOD/ESD specific implementation

  if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  if(negTrack->GetKinkIndex(0) > 0  || posTrack->GetKinkIndex(0) > 0 ) {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  return kTRUE;
}



///________________________________________________________________________
Bool_t AliConversionPhotonCuts::TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack) {
  // Track Selection for Photon Reconstruction

  Int_t cutIndex=0;
  if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
  cutIndex++;

  // avoid like sign
  if(fUseOnFlyV0FinderSameSign==0){
    if(negTrack->Charge() == posTrack->Charge()) {
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
      return kFALSE;
    }
  }else if(fUseOnFlyV0FinderSameSign==1){
    if(negTrack->Charge() != posTrack->Charge()) {
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
      return kFALSE;
    }
  }
  cutIndex++;

  // Number of TPC Clusters


  if( negTrack->GetNcls(1) < fMinClsTPC || posTrack->GetNcls(1) < fMinClsTPC ) {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  // Acceptance
  if( posTrack->Eta() > (fEtaCut) || posTrack->Eta() < (-fEtaCut) ||
    negTrack->Eta() > (fEtaCut) || negTrack->Eta() < (-fEtaCut) ){
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  if(fEtaCutMin>-0.1){
    if( (posTrack->Eta() < (fEtaCutMin) && posTrack->Eta() > (-fEtaCutMin)) ||
      (negTrack->Eta() < (fEtaCutMin) && negTrack->Eta() > (-fEtaCutMin)) ){
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
      return kFALSE;
    }
  }
  cutIndex++;

  // Single Pt Cut
  if( negTrack->Pt()< fSinglePtCut || posTrack->Pt()< fSinglePtCut){
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  // AOD ESD specific cuts
  Bool_t passCuts = kTRUE;

  if(negTrack->IsA()==AliAODTrack::Class()) {
    passCuts = passCuts * SpecificTrackCuts(static_cast<AliAODTrack*>(negTrack), static_cast<AliAODTrack*>(posTrack),cutIndex);
  } else {
    passCuts = passCuts * SpecificTrackCuts(static_cast<AliESDtrack*>(negTrack), static_cast<AliESDtrack*>(posTrack),cutIndex);
  }

  if(!passCuts){
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);

  return kTRUE;

}
///________________________________________________________________________
Float_t AliConversionPhotonCuts::GetKappaTPC(AliConversionPhotonBase *gamma, AliVEvent * event){
  
  if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
  if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error
  
  AliVTrack * negTrack = GetTrack(event, gamma->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, gamma->GetTrackLabelPositive());
  
  Float_t KappaPlus, KappaMinus, Kappa;
  KappaMinus = fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kElectron);
  KappaPlus  = fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kElectron);
  Kappa = ( TMath::Abs(KappaMinus) + TMath::Abs(KappaPlus) ) / 2.0 + 2.0*(KappaMinus+KappaPlus);
  
  return Kappa;
  
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::dEdxCuts(AliVTrack *fCurrentTrack){
  // Electron Identification Cuts for Photon reconstruction
  if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
  if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error

  Int_t cutIndex=0;
  if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
  if(fHistoTPCdEdxSigbefore)fHistoTPCdEdxSigbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));
  if(fHistoTPCdEdxbefore)fHistoTPCdEdxbefore->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
  cutIndex++;
  if(fDodEdxSigmaCut == kTRUE && !fSwitchToKappa){
    // TPC Electron Line
    if( fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaBelowElectronLine ||
      fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaAboveElectronLine){

      if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
      return kFALSE;
    }
    cutIndex++;

    // TPC Pion Line
    if( fCurrentTrack->P()>fPIDMinPnSigmaAbovePionLine && fCurrentTrack->P()<fPIDMaxPnSigmaAbovePionLine ){
      if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
        fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine&&
        fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){

        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
    cutIndex++;

    // High Pt Pion rej
    if( fCurrentTrack->P()>fPIDMaxPnSigmaAbovePionLine ){
      if(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)>fPIDnSigmaBelowElectronLine &&
        fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron)<fPIDnSigmaAboveElectronLine &&
        fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){

        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
    cutIndex++;
  }
  else{cutIndex+=3;}

  if(fDoKaonRejectionLowP == kTRUE && !fSwitchToKappa){
    if(fCurrentTrack->P()<fPIDMinPKaonRejectionLowP ){
      if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){

        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
  }
  cutIndex++;
  if(fDoProtonRejectionLowP == kTRUE && !fSwitchToKappa){
    if( fCurrentTrack->P()<fPIDMinPProtonRejectionLowP ){
      if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){

        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
  }
  cutIndex++;

  if(fDoPionRejectionLowP == kTRUE && !fSwitchToKappa){
    if( fCurrentTrack->P()<fPIDMinPPionRejectionLowP ){
      if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){

        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
  }
  cutIndex++;


  // cout<<"Start"<<endl;
  // AliPIDResponse::EDetPidStatus status=fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,fCurrentTrack);

  // if( ( (status & AliVTrack::kTOFout) == AliVTrack::kTOFout ) && ( (status & AliVTrack::kTIME) == AliVTrack::kTIME ))
  //    {cout<<"TOF DA"<<endl;}
  // if(status == AliPIDResponse::kDetPidOk){
  //    Float_t probMis = fPIDResponse->GetTOFMismatchProbability(fCurrentTrack);
  //    cout<<"--> "<<probMis<<endl;
  //    if(probMis > 0.01){

  //    }
  // }

  if((fCurrentTrack->GetStatus() & AliESDtrack::kTOFpid) && !(fCurrentTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
    if(fHistoTOFbefore){
      Double_t t0 = fPIDResponse->GetTOFResponse().GetStartTime(fCurrentTrack->P());
      Double_t  times[AliPID::kSPECIESC];
      fCurrentTrack->GetIntegratedTimes(times,AliPID::kSPECIESC);
      Double_t TOFsignal = fCurrentTrack->GetTOFsignal();
      Double_t dT = TOFsignal - t0 - times[0];
      fHistoTOFbefore->Fill(fCurrentTrack->P(),dT);
    }
    if(fHistoTOFSigbefore) fHistoTOFSigbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
    if(fUseTOFpid){
      if(fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)>fTofPIDnSigmaAboveElectronLine ||
        fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)<fTofPIDnSigmaBelowElectronLine ){
        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
    if(fHistoTOFSigafter)fHistoTOFSigafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
  }
  cutIndex++;
  
  if((fCurrentTrack->GetStatus() & AliESDtrack::kITSpid)){
    if(fHistoITSSigbefore) fHistoITSSigbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron));
    if(fUseITSpid){
      if(fCurrentTrack->Pt()<=fMaxPtPIDITS){
        if(fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron)>fITSPIDnSigmaAboveElectronLine || fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron)<fITSPIDnSigmaBelowElectronLine ){
          if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
          return kFALSE;
        }
      }  
    }
    if(fHistoITSSigafter)fHistoITSSigafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron));
  }
  
  cutIndex++;
  
  // Apply TRD PID
  if(fDoTRDPID){
    if(!fPIDResponse->IdentifiedAsElectronTRD(fCurrentTrack,fPIDTRDEfficiency)){
      if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
      return kFALSE;
    }
  }
  cutIndex++;

  if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
  if(fHistoTPCdEdxSigafter)fHistoTPCdEdxSigafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTPC(fCurrentTrack, AliPID::kElectron));
  if(fHistoTPCdEdxafter)fHistoTPCdEdxafter->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
  
  return kTRUE;
}

Bool_t AliConversionPhotonCuts::KappaCuts(AliConversionPhotonBase * photon,AliVEvent *event) {
  // abort if Kappa selection not enabled
  if (!fSwitchToKappa) return kTRUE;
  
  Float_t kappa = GetKappaTPC(photon, event);
  if (kappa < fKappaMinCut) return kFALSE;
  if (kappa > fKappaMaxCut) return kFALSE;
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::AsymmetryCut(AliConversionPhotonBase * photon,AliVEvent *event) {
  // Cut on Energy Assymetry

  for(Int_t ii=0;ii<2;ii++){

    AliVTrack *track=GetTrack(event,photon->GetTrackLabel(ii));

    if(fDoPhotonPDependentAsymCut){
      Double_t trackNegAsy=0;
      if (photon->GetPhotonP()!=0.){
          trackNegAsy= track->P()/photon->GetPhotonP();
      }

      if( trackNegAsy > fFAsymmetryCut->Eval(photon->GetPhotonP()) || trackNegAsy < 1.-fFAsymmetryCut->Eval(photon->GetPhotonP()) ){
        return kFALSE;
      }
    
    } else {
      if( track->P() > fMinPPhotonAsymmetryCut ){
        Double_t trackNegAsy=0;
        if (photon->GetPhotonP()!=0.){
          trackNegAsy= track->P()/photon->GetPhotonP();
        }

        if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
          return kFALSE;
        }
      }
    }

  }
  return kTRUE;
}

///________________________________________________________________________
AliVTrack *AliConversionPhotonCuts::GetTrack(AliVEvent * event, Int_t label){
  //Returns pointer to the track with given ESD label
  //(Important for AOD implementation, since Track array in AOD data is different
  //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);
  if(esdEvent) {
    if(label > event->GetNumberOfTracks() ) return NULL;
    AliESDtrack * track = esdEvent->GetTrack(label);
    return track;

  } else {
    if(label == -999999) return NULL; // if AOD relabelling goes wrong, immediately return NULL
    AliVTrack * track = 0x0;
    if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
      if(event->GetTrack(label)) track = dynamic_cast<AliVTrack*>(event->GetTrack(label));
      return track;
    }
    else{
      for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
        if(event->GetTrack(ii)) track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
        if(track){
          if(track->GetID() == label) {
            return track;
          }
        }
      }
    }
  }
  //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
  return NULL;
}

///________________________________________________________________________
AliESDtrack *AliConversionPhotonCuts::GetESDTrack(AliESDEvent * event, Int_t label){
  //Returns pointer to the track with given ESD label
  //(Important for AOD implementation, since Track array in AOD data is different
  //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  if(event) {
    if(label > event->GetNumberOfTracks() ) return NULL;
    AliESDtrack * track = event->GetTrack(label);
    return track;
  }
  //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
  return NULL;
}



///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event){
  // Cut on Electron Probability for Photon Reconstruction

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);

  if(esdEvent){

    Bool_t iResult=kFALSE;

    Double_t *posProbArray = new Double_t[AliPID::kSPECIES];
    Double_t *negProbArray = new Double_t[AliPID::kSPECIES];

    AliESDtrack* negTrack   = esdEvent->GetTrack(photon->GetTrackLabelNegative());
    AliESDtrack* posTrack   = esdEvent->GetTrack(photon->GetTrackLabelPositive());

    if(negProbArray && posProbArray){

      negTrack->GetTPCpid(negProbArray);
      posTrack->GetTPCpid(posProbArray);

      if(negProbArray[AliPID::kElectron]>=fPIDProbabilityCutNegativeParticle && posProbArray[AliPID::kElectron]>=fPIDProbabilityCutPositiveParticle){
        iResult=kTRUE;
      }
    }

    delete [] posProbArray;
    delete [] negProbArray;
    return iResult;

  } else {
    ///Not possible for AODs
    return kTRUE;
  }
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg){
  // MC Acceptance Cuts
  //(Certain areas were excluded for photon reconstruction)

  if(particle->R()>fMaxR){
    return kFALSE;}

  if(ePos->R()>fMaxR){
    return kFALSE;
  }

  if(ePos->R()<fMinR){
    return kFALSE;
  }

  if( ePos->R() <= ((TMath::Abs(ePos->Vz())*fLineCutZRSlope)-fLineCutZValue)){
    return kFALSE;
  }
  else if (fUseEtaMinCut &&  ePos->R() >= ((TMath::Abs(ePos->Vz())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
    return kFALSE;
  }

  if(TMath::Abs(eNeg->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
    return kFALSE;
  }

  if(eNeg->Vz()!=ePos->Vz()||eNeg->R()!=ePos->R()){
    return kFALSE;
  }

  if(TMath::Abs(ePos->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
    return kFALSE;
  }


  if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) ){
    return kFALSE;
  }
  if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ){
    return kFALSE;
  }
  if( eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) ){
    return kFALSE;
  }
  if(fEtaCutMin>-0.1){
    if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) ){
      return kFALSE;
    }
    if( ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin) ){
      return kFALSE;
    }
    if( eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin) ){
      return kFALSE;
    }
  }

  if( ePos->Pt()< fSinglePtCut ||  eNeg->Pt()< fSinglePtCut){
    return kFALSE;
  }

  if(particle->Pt()<fPtCut){
    return kFALSE;
  }

  return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::UpdateCutString() {
  ///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
    fCutString->SetString(GetCutNumber());
  } else {
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  // Initialize Cuts from a given Cut string
  AliInfo(Form("Set Photoncut Number: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsDigit()){
    AliError("Cut selection contains characters");
    return kFALSE;
  }

  const char *cutSelection = analysisCutSelection.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = cutSelection[i] - '0'
  for(Int_t ii=0;ii<kNCuts;ii++){
    ASSIGNARRAY(ii);
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
    if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  PrintCutsWithValues();

  return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  switch (cutID) {

    case kv0FinderType:
      if( SetV0Finder(value)) {
        fCuts[kv0FinderType] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case ketaCut:
      if( SetEtaCut(value)) {
        fCuts[ketaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kRCut:
      if( SetRCut(value)) {
        fCuts[kRCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
  
    case kEtaForPhiSector:
      if( SetEtaForPhiCut(value)) {
        fCuts[kEtaForPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    case kMinPhiSector:
      if( SetMinPhiSectorCut(value)) {
        fCuts[kMinPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    case kMaxPhiSector:
      if( SetMaxPhiSectorCut(value)) {
        fCuts[kMaxPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case ksinglePtCut:
      if( SetSinglePtCut(value)) {
        fCuts[ksinglePtCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kclsTPCCut:
      if( SetTPCClusterCut(value)) {
        fCuts[kclsTPCCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
      
    case kededxSigmaCut:
      if (!fSwitchToKappa){
        if( SetTPCdEdxCutElectronLine(value)) {
          fCuts[kededxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        if( SetKappaTPCCut(value)) {
          fCuts[kededxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      }  
    case kpidedxSigmaCut:
      if (!fSwitchToKappa){
        if( SetTPCdEdxCutPionLine(value)) {
          fCuts[kpidedxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kpidedxSigmaCut] = 0;
        return kTRUE;
      }  
    case kpiMomdedxSigmaCut:
      if (!fSwitchToKappa){
        if( SetMinMomPiondEdxCut(value)) {
          fCuts[kpiMomdedxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kpiMomdedxSigmaCut] = 0;
        return kTRUE;
      } 
    case kpiMaxMomdedxSigmaCut:
      if (!fSwitchToKappa){
        if( SetMaxMomPiondEdxCut(value)) {
          fCuts[kpiMaxMomdedxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kpiMaxMomdedxSigmaCut] = 0;
        return kTRUE;
      } 
    case kLowPRejectionSigmaCut:
      if (!fSwitchToKappa){
        if( SetLowPRejectionCuts(value)) {
          fCuts[kLowPRejectionSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kLowPRejectionSigmaCut] = 0;
        return kTRUE;
      } 
    case kTOFelectronPID:
      if( SetTOFElectronPIDCut(value)) {
        fCuts[kTOFelectronPID] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
      
    case kQtMaxCut:
      if( SetQtMaxCut(value)) {
        fCuts[kQtMaxCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    
    case kchi2GammaCut:
      if( SetChi2GammaCut(value)) {
        fCuts[kchi2GammaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kPsiPair:
      if( SetPsiPairCut(value)) {
        fCuts[kPsiPair] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
      
    case kdoPhotonAsymmetryCut:
      if( SetPhotonAsymmetryCut(value)) {
        fCuts[kdoPhotonAsymmetryCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kCosPAngle:
      if( SetCosPAngleCut(value)) {
        fCuts[kCosPAngle] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kElecShare:
      if( SetSharedElectronCut(value)) {
        fCuts[kElecShare] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kToCloseV0s:
      if( SetToCloseV0sCut(value)) {
        fCuts[kToCloseV0s] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDcaRPrimVtx:
      if( SetDCARPhotonPrimVtxCut(value)) {
        fCuts[kDcaRPrimVtx] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDcaZPrimVtx:
      if( SetDCAZPhotonPrimVtxCut(value)) {
        fCuts[kDcaZPrimVtx] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kInPlaneOutOfPlane:
      if( SetInPlaneOutOfPlane(value)) {
        fCuts[kInPlaneOutOfPlane] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
      
    case kITSelectronPID:
    if( SetITSElectronPIDCut(value)) {
      fCuts[kITSelectronPID] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

    case kTRDelectronPID:
    if( SetTRDElectronPIDCut(value)) {
      fCuts[kTRDelectronPID] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
    
    case kNCuts:
      AliError("Cut id out of range");
      return kFALSE;
  }

  AliError("Cut id %d not recognized");
  return kFALSE;
}
///________________________________________________________________________
void AliConversionPhotonCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

void AliConversionPhotonCuts::PrintCutsWithValues() {
  // Print out current Cut Selection with value
  printf("\nConversion cutnumber \n");
  for(Int_t ic = 0; ic < kNCuts; ic++) {
    printf("%d",fCuts[ic]);
  }
  printf("\n\n");	
  printf("Electron cuts & Secondary Track Cuts - only track from secondaries enter analysis: \n");
  printf("\t no like sign pairs from V0s \n");
  if (!fUseCorrectedTPCClsInfo) printf("\t # TPC clusters > %3.2f \n", fMinClsTPC);
  if (fEtaCutMin > -0.1) printf("\t %3.2f < eta_{e} < %3.2f\n", fEtaCutMin, fEtaCut );
    else printf("\t eta_{e} < %3.2f\n", fEtaCut );
  printf("\t reject: %3.2f < phi < %3.2f with %3.2f < eta < %3.2f  \n", fMinPhiCut, fMaxPhiCut, fEtaForPhiCutMin, fEtaForPhiCutMax);	
  printf("\t p_{T,e} > %3.2f\n", fSinglePtCut );
  printf("\t TPC refit \n");
  printf("\t no kinks \n");
  if (!fSwitchToKappa){
    printf("\t accept: %3.2f < n sigma_{e,TPC} < %3.2f\n", fPIDnSigmaBelowElectronLine, fPIDnSigmaAboveElectronLine );
    printf("\t reject: %3.2f < p_{e,T} < %3.2f, n sigma_{pi,TPC} < %3.2f\n", fPIDMinPnSigmaAbovePionLine, fPIDMaxPnSigmaAbovePionLine, fPIDnSigmaAbovePionLine );
    printf("\t reject: p_{e,T} > %3.2f, n sigma_{pi,TPC} < %3.2f\n", fPIDMaxPnSigmaAbovePionLine, fPIDnSigmaAbovePionLineHighPt );
    if (fDoPionRejectionLowP) printf("\t reject: p_{e,T} < %3.2f, -%3.2f < n sigma_{pi,TPC} < %3.2f\n", fPIDMinPPionRejectionLowP, fPIDnSigmaAtLowPAroundPionLine, fPIDnSigmaAtLowPAroundPionLine );
    if (fDoKaonRejectionLowP) printf("\t reject: -%3.2f < n sigma_{K,TPC} < %3.2f\n", fPIDnSigmaAtLowPAroundKaonLine, fPIDnSigmaAtLowPAroundKaonLine );
    if (fDoProtonRejectionLowP) printf("\t reject: -%3.2f < n sigma_{p,TPC} < %3.2f\n", fPIDnSigmaAtLowPAroundProtonLine, fPIDnSigmaAtLowPAroundProtonLine );
  } else {
    printf("\t accept: %3.2f <= Kappa_{TPC} < %3.2f\n", fKappaMinCut, fKappaMaxCut );
  }  
  if (fUseTOFpid) printf("\t accept: %3.2f < n sigma_{e,TOF} < %3.2f\n", fTofPIDnSigmaBelowElectronLine, fTofPIDnSigmaAboveElectronLine);
  if (fUseITSpid) printf("\t accept: %3.2f < n sigma_{e,ITS} < %3.2f\n -- up to pT %3.2f", fITSPIDnSigmaBelowElectronLine, fITSPIDnSigmaAboveElectronLine, fMaxPtPIDITS);
  
  printf("Photon cuts: \n");
  if (fUseOnFlyV0Finder) printf("\t using Onfly V0 finder \n");
  else printf("\t using Offline V0 finder \n");
  if (fDo2DQt){
    printf("\t 2 dimensional q_{T} cut applied with maximum of %3.2f \n", fQtMax );
  } else {
    printf("\t 1 dimensional q_{T} cut applied with maximum of %3.2f \n", fQtMax );
  }
  if (fDo2DPsiPairChi2){
    printf("\t 2 dimensional triangle chi^{2} and psi_{pair} cut applied with maximum of chi^{2} = %3.2f and |psi_{pair}| = %3.2f \n", fChi2CutConversion, fPsiPairCut ); 
  } else {
    printf("\t chi^{2} max cut chi^{2} < %3.2f \n", fChi2CutConversion ); 
    printf("\t psi_{pair} max cut |psi_{pair}| < %3.2f \n", fPsiPairCut ); 
  }	      
  printf("\t %3.2f < R_{conv} < %3.2f\n", fMinR, fMaxR );
  printf("\t Z_{conv} < %3.2f\n", fMaxZ );
  if (fEtaCutMin > -0.1) printf("\t %3.2f < eta_{conv} < %3.2f\n", fEtaCutMin, fEtaCut );
    else printf("\t eta_{conv} < %3.2f\n", fEtaCut );
  if (fDoPhotonAsymmetryCut) printf("\t for p_{T,track} > %3.2f,  A_{gamma} < %3.2f \n", fMinPPhotonAsymmetryCut, fMinPhotonAsymmetry  );
  if (fDoPhotonPDependentAsymCut && fDoPhotonAsymmetryCut) printf("\t p-dependent asymmetry cut \n");
  if (fUseCorrectedTPCClsInfo) printf("\t #cluster TPC/ #findable clusters TPC (corrected for radius) > %3.2f\n", fMinClsTPCToF );
  printf("\t p_{T,gamma} > %3.2f\n", fPtCut );	 
  printf("\t cos(Theta_{point}) > %3.2f \n", fCosPAngleCut );
  printf("\t dca_{R} < %3.2f \n", fDCARPrimVtxCut );
  printf("\t dca_{Z} < %3.2f \n", fDCAZPrimVtxCut );
  if (fDoPhotonQualitySelectionCut) printf("\t selection based on photon quality with quality %d \n", fPhotonQualityCut );
  if (fDoDoubleCountingCut) printf("\t Reject doubly counted photons with R > %3.2f, DeltaR < %3.2f, OpenAngle < %3.2f  \n", fMinRDC, fDeltaR,fOpenAngle );

}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetV0Finder(Int_t v0FinderType){   // Set Cut
  switch (v0FinderType){
  case 0:  // on fly V0 finder
    cout << "have chosen onfly V0" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=0;
    break;
  case 1:  // offline V0 finder
    cout << "have chosen offline V0" << endl;
    fUseOnFlyV0Finder=kFALSE;
    fUseOnFlyV0FinderSameSign=0;
    break;
  case 2:  // on fly V0 finder with same signs
    cout << "have chosen onfly V0 same sign pairing" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=1;
    break;
  case 3:  // on fly V0 finder with unlike signs and same signs
    cout << "have chosen onfly V0 unlike sign and same signs pairing" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=2;
    break;
  default:
    AliError(Form(" v0FinderType not defined %d",v0FinderType));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetEtaCut(Int_t etaCut){   // Set Cut

  //Set Standard LineCutZValues
  fLineCutZValueMin = -2;
  fLineCutZValue = 7.;

  switch(etaCut){
  case 0: // 0.9
    fEtaCut     = 0.9;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 1:  // 0.6  // changed from 1.2 to 0.6 on 2013.06.10
    fEtaCut     = 0.6;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 2:  // 1.4
    fEtaCut     = 1.4;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 3: // 0.65
    fEtaCut     = 0.65;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 4: // 0.75
    fEtaCut     = 0.75;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 5: // 0.5
    fEtaCut     = 0.5;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 6: // 5.
    fEtaCut     = 5.;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 7:
    if (fIsHeavyIon==1){
      fEtaCut     = 0.7;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
    } else {   
      fEtaCut     = 0.3;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
    }
  // case 8: // 0.1 - 0.8
  //    fEtaCut     = 0.9;
  //    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
  //    fEtaCutMin     = 0.1;
  //    fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
  //    break;
  case 8: // 0.4
    fEtaCut     = 0.4;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 9: // 10
    fEtaCut     = 10;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  default:
    AliError(Form(" EtaCut not defined %d",etaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetRCut(Int_t RCut){
  // Set Cut
  switch(RCut){
  case 0:
    fMinR=0;
    fMaxR = 180.;
    break;
  case 1:
    fMinR=2.8;
    fMaxR = 180.;
    break;
  case 2:
    fMinR=5.;
    fMaxR = 180.;
    break;
  case 3:
    fMaxR = 70.;
    fMinR = 10.;
    break;
  case 4:
    fMaxR = 70.;
    fMinR = 5.;
    break;
  case 5:
    fMaxR = 180.;
    fMinR = 10.;
    break;
  case 6:
    fMaxR = 180.;
    fMinR = 20.;
    break;
  case 7:
    fMaxR = 180.;
    fMinR = 35.; //old 26.
    break;
  case 8:
    fMaxR = 180.;
    fMinR = 12.5;
    break;
  case 9:
    fMaxR = 180.;
    fMinR = 7.5;
    break;

  default:
    AliError("RCut not defined");
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetEtaForPhiCut(Int_t etaPhiCut) {

  switch(etaPhiCut) {
  case 0: //no specific eta range selected, full eta range
    fEtaForPhiCutMin = -fEtaCut;
    fEtaForPhiCutMax = fEtaCut;
    break;
  case 1:  //eta < 0 only
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fEtaForPhiCutMin = -fEtaCut;
    fEtaForPhiCutMax = 0.;
    break;
  case 2://eta > 0 only
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fEtaForPhiCutMin = 0.;
    fEtaForPhiCutMax = fEtaCut;
    break;
  default:
    AliError(Form("EtaForPhiCut not defined %d",etaPhiCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
// This are exclusion cuts, everything between fMinPhiCut & fMaxPhiCut will be excluded
Bool_t AliConversionPhotonCuts::SetMinPhiSectorCut(Int_t minPhiCut) {

  switch(minPhiCut) {
  case 0:
    fDoShrinkTPCAcceptance = kFALSE;
    fMinPhiCut = 0;
    break;
  case 1:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMinPhiCut = 1.7; //OROC C08 large cut
    break;
  case 2:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMinPhiCut = 4.4; //EMCal
    break;
  case 3:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMinPhiCut = 1.0; //PHOS
    break;
  case 4:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMinPhiCut = 3.4; //EMCal tight
    break;
  case 5:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMinPhiCut = 2.0; //OROC C08 medium cut 
    break;
  case 6:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMinPhiCut = 2.2; //OROC C08 small cut
    break;
  case 7:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMinPhiCut = 2.4; //OROC C08 tightest cut
    break;
    
  default:
    AliError(Form("MinPhiCut not defined %d",minPhiCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
// This are exclusion cuts, everything between fMinPhiCut & fMaxPhiCut will be excluded
Bool_t AliConversionPhotonCuts::SetMaxPhiSectorCut(Int_t maxPhiCut) {

  switch(maxPhiCut) {
  case 0:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kFALSE;
    fMaxPhiCut = 2*TMath::Pi()+0.00001;
    break;
  case 1:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMaxPhiCut = 4.3; //OROC C08 large cut
    break;
  case 2:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMaxPhiCut = 5.8; //EMCal
    break;
  case 3:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMaxPhiCut = 3.0; //PHOS
    break;
  case 4:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMaxPhiCut = 1.; //EMCal
    break;
  case 5:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMaxPhiCut = 4.0; //OROC C08 medium cut 
    break;
  case 6:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMaxPhiCut = 3.8; //OROC C08 small cut
    break;
  case 7:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = kTRUE;
    fMaxPhiCut = 3.6; //OROC C08 tighest cut
    break;
    
  default:
    AliError(Form("MaxPhiCut not defined %d",maxPhiCut));
    return kFALSE;
  }

  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetSinglePtCut(Int_t singlePtCut){   // Set Cut
  switch(singlePtCut){
  case 0: // 0.050 GeV
    fSinglePtCut = 0.050;
    break;
  case 1:  // 0.100 GeV
    fSinglePtCut = 0.100;
    break;
  case 2:  // 0.150 GeV
    fSinglePtCut = 0.150;
    break;
  case 3:  // 0.200 GeV
    fSinglePtCut = 0.200;
    break;
  case 4:  // 0.075 GeV
    fSinglePtCut = 0.075;
    break;
  case 5:  // 0.125 GeV
    fSinglePtCut = 0.125;
    break;
  case 6:  // 0.04 GeV
    fSinglePtCut = 0.040;
    break;
  case 7:  // 0.0 GeV
    fSinglePtCut = 0.0;
    break;
  case 8:  // 0.02 GeV ; equivalent to .05 for the low B field runs 
    fSinglePtCut = 0.02;
    break;
  default:
    AliError(Form("singlePtCut not defined %d",singlePtCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTPCClusterCut(Int_t clsTPCCut){   // Set Cut
  switch(clsTPCCut){
  case 0: // 0
    fMinClsTPC= 0.;
    break;
  case 1:  // 60
    fMinClsTPC= 60.;
    break;
  case 2:  // 80
    fMinClsTPC= 80.;
    break;
  case 3:  // 100
    fMinClsTPC= 100.;
    break;
  case 4:  // 95% of findable clusters
    fMinClsTPCToF= 0.95;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 5:  // 0% of findable clusters
    fMinClsTPCToF= 0.0;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 6:  // 70% of findable clusters
    fMinClsTPCToF= 0.7;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 7:  // 0% of findable clusters
    fMinClsTPCToF= 0.35;
    fUseCorrectedTPCClsInfo=0;
    break;
  case 8:
    fMinClsTPCToF= 0.35;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 9:
    fMinClsTPCToF= 0.6;
    fUseCorrectedTPCClsInfo=1;
    break;
  default:
    AliError(Form("Warning: clsTPCCut not defined %d",clsTPCCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut){   // Set Cut
  switch(ededxSigmaCut){
  case 0: // -10,10
    fPIDnSigmaBelowElectronLine=-10;
    fPIDnSigmaAboveElectronLine=10;
    break;
  case 1: // -5,5
    fPIDnSigmaBelowElectronLine=-5;
    fPIDnSigmaAboveElectronLine=5;
    break;
  case 2: // -3,5
    fPIDnSigmaBelowElectronLine=-3;
    fPIDnSigmaAboveElectronLine=5;
    break;
  case 3: // -4,5
    fPIDnSigmaBelowElectronLine=-4;
    fPIDnSigmaAboveElectronLine=5;
    break;
  case 4: // -6,7
    fPIDnSigmaBelowElectronLine=-6;
    fPIDnSigmaAboveElectronLine=7;
    break;
  case 5: // -4,4
    fPIDnSigmaBelowElectronLine=-4;
    fPIDnSigmaAboveElectronLine=4;
    break;
  case 6: // -2.5,4
    fPIDnSigmaBelowElectronLine=-2.5;
    fPIDnSigmaAboveElectronLine=4;
    break;
  case 7: // -2,3.5
    fPIDnSigmaBelowElectronLine=-2;
    fPIDnSigmaAboveElectronLine=3.5;
    break;
  case 8: // -2.5,3.
    fPIDnSigmaBelowElectronLine=-2.5;
    fPIDnSigmaAboveElectronLine=3;
    break;
  default:
    AliError("TPCdEdxCutElectronLine not defined");
    return kFALSE;

  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut){   // Set Cut

  switch(pidedxSigmaCut){
  case 0:  // -10
    fPIDnSigmaAbovePionLine=-10;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 1:   // 0
    fPIDnSigmaAbovePionLine=0;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 2:  // 1
    fPIDnSigmaAbovePionLine=1;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 3:  // 1
    fPIDnSigmaAbovePionLine=2.5;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 4:  // 3.0sigma, 1.0 sigma at high pt
    fPIDnSigmaAbovePionLine=3.0;
    fPIDnSigmaAbovePionLineHighPt=1.;
    break;
  case 5:  // 1
    fPIDnSigmaAbovePionLine=2.;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 6:  // 1
    fPIDnSigmaAbovePionLine=2.;
    fPIDnSigmaAbovePionLineHighPt=0.5;
    break;
  case 7:  // 1
    fPIDnSigmaAbovePionLine=3.5;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 8:  // 1
    fPIDnSigmaAbovePionLine=2.;
    fPIDnSigmaAbovePionLineHighPt=1.;
    break;
  case 9:
    fPIDnSigmaAbovePionLine=1; // We need a bit less tight cut on dE/dx
    fPIDnSigmaAbovePionLineHighPt=0.5;
    break;
  default:
    AliError(Form("Warning: pidedxSigmaCut not defined %d",pidedxSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetMinMomPiondEdxCut(Int_t piMomdedxSigmaCut){   // Set Cut
  switch(piMomdedxSigmaCut){
  case 0:  // 0.5 GeV
    fPIDMinPnSigmaAbovePionLine=0.5;
    break;
  case 1:  // 1. GeV
    fPIDMinPnSigmaAbovePionLine=1.;
    break;
  case 2:  // 1.5 GeV
    fPIDMinPnSigmaAbovePionLine=1.5;
    break;
  case 3:  // 20.0 GeV
    fPIDMinPnSigmaAbovePionLine=20.;
    break;
  case 4:  // 50.0 GeV
    fPIDMinPnSigmaAbovePionLine=50.;
    break;
  case 5:  // 0.3 GeV
    fPIDMinPnSigmaAbovePionLine=0.3;
    break;
  case 6:  // 0.25 GeV
    fPIDMinPnSigmaAbovePionLine=0.25;
    break;
  case 7:  // 0.4 GeV
    fPIDMinPnSigmaAbovePionLine=0.4;
    break;
  case 8:  // 0.2 GeV
    fPIDMinPnSigmaAbovePionLine=0.2;
    break;
  default:
    AliError(Form("piMomdedxSigmaCut not defined %d",piMomdedxSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut){   // Set Cut
  switch(piMaxMomdedxSigmaCut){
  case 0:  // 100. GeV
    fPIDMaxPnSigmaAbovePionLine=100.;
    break;
  case 1:  // 5. GeV
    fPIDMaxPnSigmaAbovePionLine=5.;
    break;
  case 2:  // 4. GeV
    fPIDMaxPnSigmaAbovePionLine=4.;
    break;
  case 3:  // 3.5 GeV
    fPIDMaxPnSigmaAbovePionLine=3.5;
    break;
  case 4:  // 3. GeV
    fPIDMaxPnSigmaAbovePionLine=3.;
    break;
  case 5:  // 7. GeV
    fPIDMaxPnSigmaAbovePionLine=7.;
    break;
  case 6:  // 2. GeV
    fPIDMaxPnSigmaAbovePionLine=2.;
    break;
  default:
    AliError(Form("piMaxMomdedxSigmaCut not defined %d",piMaxMomdedxSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut){   // Set Cut
  switch(LowPRejectionSigmaCut){
  case 0:  //
    fPIDnSigmaAtLowPAroundKaonLine=0;
    fPIDnSigmaAtLowPAroundProtonLine=0;
    fPIDnSigmaAtLowPAroundPionLine=0;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kFALSE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 1:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.5;
    fPIDnSigmaAtLowPAroundProtonLine=0.5;
    fPIDnSigmaAtLowPAroundPionLine=0.5;
    fDoKaonRejectionLowP = kTRUE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 2:  //
    fPIDnSigmaAtLowPAroundKaonLine=1;
    fPIDnSigmaAtLowPAroundProtonLine=1;
    fPIDnSigmaAtLowPAroundPionLine=1;
    fDoKaonRejectionLowP = kTRUE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 3:  //
    fPIDnSigmaAtLowPAroundKaonLine=2.;
    fPIDnSigmaAtLowPAroundProtonLine=2.;
    fPIDnSigmaAtLowPAroundPionLine=2.;
    fDoKaonRejectionLowP = kTRUE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 4:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=1;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 5:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=1.5;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 6:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=2.;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 7:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=0.5;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 8:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.5;
    fPIDnSigmaAtLowPAroundPionLine=0.;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kFALSE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  default:
    AliError(Form("LowPRejectionSigmaCut not defined %d",LowPRejectionSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetKappaTPCCut(Int_t kappaCut){   // Set Cut
  switch(kappaCut){
  case 0: // completely open
    fKappaMaxCut=200;
    fKappaMinCut=-200;
    break;
  case 1: // mainly pi pi 
    fKappaMaxCut=-13;
    fKappaMinCut=-20;
    break;
  case 2: // mainly pi e
    fKappaMaxCut=-6;
    fKappaMinCut=-11;
    break;
  case 3: // signal
    fKappaMaxCut=5;
    fKappaMinCut=-3;
    break;
  case 4: // remaining
    fKappaMaxCut=20;
    fKappaMinCut=11;
    break;
  case 5: // -5-10 full signal peak(including background)
    fKappaMaxCut=10;
    fKappaMinCut=-5;
    break;
  case 6: //
    fKappaMaxCut=10;
    fKappaMinCut=-3;
    break;
  case 7: //
    fKappaMaxCut=10;
    fKappaMinCut=0;
    break;
  default:
    AliError("KappaTPCCut not defined");
    return kFALSE;

  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTOFElectronPIDCut(Int_t TOFelectronPID){
  // Set Cut
  switch(TOFelectronPID){
  case 0: // no cut
    fUseTOFpid = kFALSE;
    fTofPIDnSigmaBelowElectronLine=-100;
    fTofPIDnSigmaAboveElectronLine=100;
    break;
  case 1: // -7,7
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-7;
    fTofPIDnSigmaAboveElectronLine=7;
    break;
  case 2: // -5,5
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-5;
    fTofPIDnSigmaAboveElectronLine=5;
    break;
  case 3: // -3,5
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-3;
    fTofPIDnSigmaAboveElectronLine=5;
    break;
  case 4: // -2,3
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-2;
    fTofPIDnSigmaAboveElectronLine=3;
    break;
  case 5: // -3,3
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-3;
    fTofPIDnSigmaAboveElectronLine=3;
    break;
  default:
    AliError(Form("TOFElectronCut not defined %d",TOFelectronPID));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetITSElectronPIDCut(Int_t ITSelectronPID){
  // Set Cut
  switch(ITSelectronPID){
  case 0: // no cut
    fUseITSpid = kFALSE;
    fITSPIDnSigmaBelowElectronLine=-100;
    fITSPIDnSigmaAboveElectronLine=100;
    fMaxPtPIDITS = 1.5;
    break;
  case 1: // -3,3
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=3;
    fMaxPtPIDITS = 1.5;
    break;
  case 2: // -2,2
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-2;
    fITSPIDnSigmaAboveElectronLine=2;
    fMaxPtPIDITS = 1.5;
    break;
  case 3: // -1,1
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-1;
    fITSPIDnSigmaAboveElectronLine=1;
    fMaxPtPIDITS = 1.5;
    break;
  case 4: // -3,5
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=5;
    fMaxPtPIDITS = 1.5;
    break;
  case 5: // -5,5
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-5;
    fITSPIDnSigmaAboveElectronLine=5;
    fMaxPtPIDITS = 1.5;
    break;
  case 6: // -3,3
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=3;
    fMaxPtPIDITS = 2;
    break;
  case 7: // -2,2
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-2;
    fITSPIDnSigmaAboveElectronLine=2;
    fMaxPtPIDITS = 2;
    break;
  case 8: // -1,1
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-1;
    fITSPIDnSigmaAboveElectronLine=1;
    fMaxPtPIDITS = 2;
    break;
  case 9: // -3,5
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=5;
    fMaxPtPIDITS = 2;
    break;
  default:
    AliError(Form("ITSelectronPID not defined %d",ITSelectronPID));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTRDElectronPIDCut(Int_t TRDelectronPID){
  // Set Cut
  switch(TRDelectronPID){
  case 0: // no cut
    fDoTRDPID = kFALSE;
    fTRDPIDBelowCut=-100;
    fTRDPIDAboveCut=100;
    break;
  default:
    AliError(Form("TRDelectronPID not defined %d",TRDelectronPID));
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetQtMaxCut(Int_t QtMaxCut){   // Set Cut
  switch(QtMaxCut){
  case 0: //
    fQtMax=1.;
    fDoQtGammaSelection=kFALSE;
    fDo2DQt=kFALSE;
    break;
  case 1:
    fQtMax=0.1;
    fDo2DQt=kFALSE;
    break;
  case 2:
    fQtMax=0.06;
    fDo2DQt=kTRUE;
    break;
  case 3:
    fQtMax=0.05;
    fDo2DQt=kFALSE;
    break;
  case 4:
    fQtMax=0.03;
    fDo2DQt=kFALSE;
    break;
  case 5:
    fQtMax=0.02;
    fDo2DQt=kFALSE;
    break;
  case 6:
    fQtMax=0.02;
    fDo2DQt=kTRUE;
    break;
  case 7:
    fQtMax=0.15;
    fDo2DQt=kFALSE;
    break;
  case 8:
    fQtMax=0.05;
    fDo2DQt=kTRUE;
    break;   
  case 9:
    fQtMax=0.03;
    fDo2DQt=kTRUE;
    break;      
  default:
    AliError(Form("Warning: QtMaxCut not defined %d",QtMaxCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetChi2GammaCut(Int_t chi2GammaCut){   // Set Cut

  switch(chi2GammaCut){
  case 0: // 100
    fChi2CutConversion = 100.;
    break;
  case 1:  // 50
    fChi2CutConversion = 50.;
    break;
  case 2:  // 30
    fChi2CutConversion = 30.;
    break;
  case 3:
    fChi2CutConversion = 200.;
    break;
  case 4:
    if (fIsHeavyIon==1){
      fChi2CutConversion = 7.;
    } else {
      fChi2CutConversion = 500.;
    }
    break;
  case 5:
    fChi2CutConversion = 100000.;
    break;
  case 6:
    fChi2CutConversion = 5.;
    break;
  case 7:
    fChi2CutConversion = 10.;
    break;
  case 8:
    fChi2CutConversion = 20.;
    break;
  case 9:
    fChi2CutConversion = 15.;
    break;
  default:
    AliError(Form("Warning: Chi2GammaCut not defined %d",chi2GammaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetPsiPairCut(Int_t psiCut) {
  switch(psiCut) {
  case 0:
    fPsiPairCut = 10000; //
    break;
  case 1:
    fPsiPairCut = 0.1; //
    break;
  case 2:
    fPsiPairCut = 0.05; // Standard
    break;
  case 3:
    fPsiPairCut = 0.035; //
    break;
  case 4:
    fPsiPairCut = 0.2; //
    break;   
  case 5:
    fPsiPairCut = 0.1; //
    fDo2DPsiPairChi2 = kTRUE;
    break;
  case 6:
    fPsiPairCut = 0.05; //
    fDo2DPsiPairChi2 = kTRUE;
    break;
  case 7:
    if (fIsHeavyIon==1){
      fPsiPairCut = 0.07; //
    } else {
      fPsiPairCut = 0.035; //
    }
    fDo2DPsiPairChi2 = kTRUE;
    break;
  case 8:
    fPsiPairCut = 0.2; //
    fDo2DPsiPairChi2 = kTRUE; //
    break;
  case 9:
    //   if (fIsHeavyIon==1){ //AM 2016-05-13
      fPsiPairCut = 0.1; //
      fDo2DPsiPairChi2 = kTRUE;
      fIncludeRejectedPsiPair = kTRUE;
      break;
    // } else {
    //   fPsiPairCut = 0.5; //
    //   break;
    // }
  default:
    AliError(Form("PsiPairCut not defined %d",psiCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut){
  // Set Cut
  switch(doPhotonAsymmetryCut){
  case 0:
    fDoPhotonAsymmetryCut=0;
    fMinPPhotonAsymmetryCut=100.;
    fMinPhotonAsymmetry=0.;
    break;
  case 1:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=3.5;
    fMinPhotonAsymmetry=0.04;
    break;
  case 2:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=3.5;
    fMinPhotonAsymmetry=0.06;
    break;
  case 3:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=0.0;
    fMinPhotonAsymmetry=0.05;
    break; 
  case 4:
    fDoPhotonAsymmetryCut=1;
    fDoPhotonPDependentAsymCut=1;
    fFAsymmetryCut = new TF1("fFAsymmetryCut","[0] + [1]*tanh(2*TMath::Power(x,[2]))",0.,100.);
    fFAsymmetryCut->SetParameter(0,0.3);
    fFAsymmetryCut->SetParameter(1,0.66);
    fFAsymmetryCut->SetParameter(2,0.7);
    fMinPPhotonAsymmetryCut=0.0;
    fMinPhotonAsymmetry=0.;
    break; 
  case 5:
    fDoPhotonAsymmetryCut=1;
    fDoPhotonPDependentAsymCut=1;
    fFAsymmetryCut = new TF1("fFAsymmetryCut","[0] + [1]*tanh(2*TMath::Power(x,[2]))",0.,100.);
    fFAsymmetryCut->SetParameter(0,0.14);
    fFAsymmetryCut->SetParameter(1,0.66);
    fFAsymmetryCut->SetParameter(2,0.5);
    fMinPPhotonAsymmetryCut=0.0;
    fMinPhotonAsymmetry=0.;
    break;
  case 6:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=6.;
    fMinPhotonAsymmetry=0.05;
    break; 
  case 7:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=8.;
    fMinPhotonAsymmetry=0.05;
    break;       
  default:
    AliError(Form("PhotonAsymmetryCut not defined %d",doPhotonAsymmetryCut));
    return kFALSE;
  }
  fCuts[kdoPhotonAsymmetryCut]=doPhotonAsymmetryCut;
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetCosPAngleCut(Int_t cosCut) {

  switch(cosCut){
  case 0:
    fCosPAngleCut = -1; 
    break;
  case 1:
    fCosPAngleCut = 0; 
    break;
  case 2:
    fCosPAngleCut = 0.5; 
    break;
  case 3:
    fCosPAngleCut = 0.75; 
    break;
  case 4:
    fCosPAngleCut = 0.85; 
    break;
  case 5:
    fCosPAngleCut = 0.88; 
    break;
  case 6:
    fCosPAngleCut = 0.9;
    break;
  case 7:
    fCosPAngleCut = 0.95;
    break;
  default:
    AliError(Form("Cosine Pointing Angle cut not defined %d",cosCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetSharedElectronCut(Int_t sharedElec) {

    switch(sharedElec){
    case 0:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kFALSE;
      fPhotonQualityCut = 0;
      break;
    case 1:
      fDoSharedElecCut = kTRUE;
      fDoPhotonQualitySelectionCut = kFALSE;
      fPhotonQualityCut = 0;
      break;
    case 2:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kTRUE;
      fPhotonQualityCut = 1;
      break;
    case 3:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kTRUE;	  
      fPhotonQualityCut = 2;
      break;
    case 4:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kTRUE;	  
      fPhotonQualityCut = 3;
      break;
    default:
      AliError(Form("Shared Electron Cut not defined %d",sharedElec));	
      return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetToCloseV0sCut(Int_t toClose) {

  switch(toClose){
  case 0:
    fDoToCloseV0sCut = kFALSE;
    fminV0Dist = 250;
    break;
  case 1:
    fDoToCloseV0sCut = kTRUE;
    fminV0Dist = 1;
    break;
  case 2:
    fDoToCloseV0sCut = kTRUE;
    fminV0Dist = 2;
    break;
  case 3:
    fDoToCloseV0sCut = kTRUE;
    fminV0Dist = 3;
    break;
  case 4:
    fDoToCloseV0sCut = kTRUE;
    fDoDoubleCountingCut = kTRUE;
    fMinRDC=0.;
    fDeltaR=6.;
    fOpenAngle=0.02;
    break;
  case 5:
    fDoToCloseV0sCut = kTRUE;
    fDoDoubleCountingCut = kTRUE;
    fMinRDC=0.;
    fDeltaR=6.;
    fOpenAngle=0.03;
    break;
  case 6:
    fDoToCloseV0sCut = kTRUE;
    fDoDoubleCountingCut = kTRUE;
    fMinRDC=0.;
    fDeltaR=6.;
    fOpenAngle=0.04;
    break;

  default:
    AliError(Form("Shared Electron Cut not defined %d",toClose));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTRDElectronCut(Int_t TRDElectronCut){   // Set Cut
  switch(TRDElectronCut){
  case 0:
    fDoTRDPID=kFALSE;
    break;
  case 1:
    fDoTRDPID=kTRUE;
    fPIDTRDEfficiency=0.1;
    break;
  case 8:
    fDoTRDPID=kTRUE;
    fPIDTRDEfficiency=0.8;
    break;
  case 9:
    fDoTRDPID=kTRUE;
    fPIDTRDEfficiency=0.9;
    break;
  default:
    AliError(Form("TRDElectronCut not defined %d",TRDElectronCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetDCAZPhotonPrimVtxCut(Int_t DCAZPhotonPrimVtx){
  // Set Cut
  switch(DCAZPhotonPrimVtx){
  case 0:  //
    fDCAZPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCAZPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCAZPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCAZPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCAZPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCAZPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCAZPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCAZPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCAZPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCAZPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCAZPhotonPrimVtx not defined "<<DCAZPhotonPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetDCARPhotonPrimVtxCut(Int_t DCARPhotonPrimVtx){
  // Set Cut
  switch(DCARPhotonPrimVtx){
  case 0:  //
    fDCARPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCARPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCARPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCARPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCARPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCARPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCARPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCARPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCARPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCARPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCARPhotonPrimVtx not defined "<<DCARPhotonPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetInPlaneOutOfPlane(Int_t inOutPlane){
  // Set Cut
  switch(inOutPlane){
  case 0:  //
    fInPlaneOutOfPlane = 0; // No Event Plane
    break;
  case 1:  //
    fInPlaneOutOfPlane = 1; // In-Plane
    break;
  case 2:  //
    fInPlaneOutOfPlane = 2; // Out-Of-Plane
    break;
  default:
    cout<<"Warning: In-Plane or Out-Of-Plane not defined "<<inOutPlane<<endl;
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
Int_t AliConversionPhotonCuts::GetFirstTPCRow(Double_t radius){
  // Get first TPC row
  Int_t firstTPCRow = 0;
  Double_t radiusI = 84.8;
  Double_t radiusO = 134.6;
  Double_t radiusOB = 198.;
  Double_t rSizeI = 0.75;
  Double_t rSizeO = 1.;
  Double_t rSizeOB = 1.5;
  Int_t nClsI = 63;
  Int_t nClsIO = 127;

  if(radius <= radiusI){
    return firstTPCRow;
  }
  if(radius>radiusI && radius<=radiusO){
    firstTPCRow = (Int_t)((radius-radiusI)/rSizeI);
  }
  if(radius>radiusO && radius<=radiusOB){
    firstTPCRow = (Int_t)(nClsI+(radius-radiusO)/rSizeO);
  }

  if(radius>radiusOB){
    firstTPCRow =(Int_t)(nClsIO+(radius-radiusOB)/rSizeOB);
  }

  return firstTPCRow;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::CosinePAngleCut(const AliConversionPhotonBase * photon, AliVEvent * event) const {
  ///Check if passes cosine of pointing angle cut
  if(GetCosineOfPointingAngle(photon, event) < fCosPAngleCut){
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Double_t AliConversionPhotonCuts::GetCosineOfPointingAngle( const AliConversionPhotonBase * photon, AliVEvent * event) const{
  // calculates the pointing angle of the recalculated V0

  Double_t momV0[3] = {0,0,0};
  if(event->IsA()==AliESDEvent::Class()){
    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(event);
    if(!esdEvent) return -999;
    AliESDv0 *v0 = esdEvent->GetV0(photon->GetV0Index());
    if(!v0) return -999;
    v0->GetPxPyPz(momV0[0],momV0[1],momV0[2]);
  }
  if(event->IsA()==AliAODEvent::Class()){
    momV0[0] = photon->GetPx();
    momV0[1] = photon->GetPy();
    momV0[2] = photon->GetPz();
  }

  //Double_t momV0[3] = { photon->GetPx(), photon->GetPy(), photon->GetPz() }; //momentum of the V0
  Double_t PosV0[3] = { photon->GetConversionX() - event->GetPrimaryVertex()->GetX(),
              photon->GetConversionY() - event->GetPrimaryVertex()->GetY(),
              photon->GetConversionZ() - event->GetPrimaryVertex()->GetZ() }; //Recalculated V0 Position vector

  Double_t momV02 = momV0[0]*momV0[0] + momV0[1]*momV0[1] + momV0[2]*momV0[2];
  Double_t PosV02 = PosV0[0]*PosV0[0] + PosV0[1]*PosV0[1] + PosV0[2]*PosV0[2];


  Double_t cosinePointingAngle = -999;
  if(momV02*PosV02 > 0.0)
    cosinePointingAngle = (PosV0[0]*momV0[0] +  PosV0[1]*momV0[1] + PosV0[2]*momV0[2] ) / TMath::Sqrt(momV02 * PosV02);

  return cosinePointingAngle;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PsiPairCut(const AliConversionPhotonBase * photon) const {

  if (fDo2DPsiPairChi2){
    if(fIncludeRejectedPsiPair){
      if (TMath::Abs(photon->GetPsiPair()) < -fPsiPairCut/fChi2CutConversion*photon->GetChi2perNDF() + fPsiPairCut || (photon->GetPsiPair()) == 4){
        return kTRUE;
      } else {
        return kFALSE;
      }

    } else {
      if (TMath::Abs(photon->GetPsiPair()) < -fPsiPairCut/fChi2CutConversion*photon->GetChi2perNDF() + fPsiPairCut ){
        return kTRUE;
      } else {
        return kFALSE;
      }
    }
  } else {
    if(fIncludeRejectedPsiPair){
      if(TMath::Abs(photon->GetPsiPair()) > fPsiPairCut || (photon->GetPsiPair()) != 4){
        return kFALSE;
      } else {
        return kTRUE;
      }
    } else {
      if(TMath::Abs(photon->GetPsiPair()) > fPsiPairCut){
        return kFALSE;
      } else {
        return kTRUE;
      }
    }
  } 
}

///________________________________________________________________________
TString AliConversionPhotonCuts::GetCutNumber(){
  // returns TString with current cut number
  TString a(kNCuts);
  for(Int_t ii=0;ii<kNCuts;ii++){
    a.Append(Form("%d",fCuts[ii]));
  }
  return a;
}

///________________________________________________________________________
void AliConversionPhotonCuts::FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0){

  Int_t posLabel = photon->GetTrackLabelPositive();
  Int_t negLabel = photon->GetTrackLabelNegative();

  fElectronLabelArray[nV0*2] = posLabel;
  fElectronLabelArray[(nV0*2)+1] = negLabel;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s){

  Int_t posLabel = photon->GetTrackLabelPositive();
  Int_t negLabel = photon->GetTrackLabelNegative();

  for(Int_t i = 0; i<nV0s*2;i++){
    if(i==nV0*2)     continue;
    if(i==(nV0*2)+1) continue;
    if(fElectronLabelArray[i] == posLabel){
      return kFALSE;}
    if(fElectronLabelArray[i] == negLabel){
      return kFALSE;}
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0){

  if (fDoDoubleCountingCut && photon->GetConversionRadius() < fMinRDC) return kTRUE;

  Double_t posX = photon->GetConversionX();
  Double_t posY = photon->GetConversionY();
  Double_t posZ = photon->GetConversionZ();

  for(Int_t i = 0;i<photons->GetEntries();i++){
    if(nV0 == i) continue;
    AliAODConversionPhoton *photonComp = (AliAODConversionPhoton*) photons->At(i);
    Double_t posCompX = photonComp->GetConversionX();
    Double_t posCompY = photonComp->GetConversionY();
    Double_t posCompZ = photonComp->GetConversionZ();

    if (!fDoDoubleCountingCut){
      Double_t dist = pow((posX - posCompX),2)+pow((posY - posCompY),2)+pow((posZ - posCompZ),2);
  
      if(dist < fminV0Dist*fminV0Dist){
        if(photon->GetChi2perNDF() > photonComp->GetChi2perNDF()) return kFALSE;
      }
    }else{         
      TVector3 v1(photon->Px(),photon->Py(),photon->Pz());
      TVector3 v2(photonComp->Px(),photonComp->Py(),photonComp->Pz());
      Double_t OpeningAngle=v1.Angle(v2);
      if( OpeningAngle < fOpenAngle && TMath::Abs(photon->GetConversionRadius()-photonComp->GetConversionRadius()) < fDeltaR){
        if(photon->GetChi2perNDF() > photonComp->GetChi2perNDF()) return kFALSE;
      }      
    }

  }
  return kTRUE;
}


///________________________________________________________________________
AliConversionPhotonCuts* AliConversionPhotonCuts::GetStandardCuts2010PbPb(){
  //Create and return standard 2010 PbPb cuts
  AliConversionPhotonCuts *cuts=new AliConversionPhotonCuts("StandardCuts2010PbPb","StandardCuts2010PbPb");
  if(!cuts->InitializeCutsFromCutString("04209297002322000000")){
    cout<<"Warning: Initialization of Standardcuts2010PbPb failed"<<endl;}
  return cuts;
}

///________________________________________________________________________
AliConversionPhotonCuts* AliConversionPhotonCuts::GetStandardCuts2010pp(){
  //Create and return standard 2010 PbPb cuts
  AliConversionPhotonCuts *cuts=new AliConversionPhotonCuts("StandardCuts2010pp","StandardCuts2010pp");
  if(!cuts->InitializeCutsFromCutString("00209366300380000000")){
    cout<<"Warning: Initialization of Standardcuts2010pp failed"<<endl;}
  return cuts;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::InPlaneOutOfPlaneCut(Double_t photonPhi, Double_t eventPlaneAngle, Bool_t fill){
  
  //GetPhotonPhi() 0-2 Pi  //eventPlaneAngle -1pi-1pi
  eventPlaneAngle=eventPlaneAngle+TMath::Pi();
  Double_t gammaToEPAngle = eventPlaneAngle-photonPhi;
  if(gammaToEPAngle < 0) gammaToEPAngle=gammaToEPAngle+2*TMath::Pi();
  gammaToEPAngle = gammaToEPAngle-TMath::Pi(); // angle from -pi +pi

  if(!fInPlaneOutOfPlane){
    if(fill&&fHistoEventPlanePhi)fHistoEventPlanePhi->Fill(gammaToEPAngle);
    return kTRUE;
  }
  else if(fInPlaneOutOfPlane == 1){
    if(TMath::Abs(gammaToEPAngle)<=0.25*TMath::Pi() || TMath::Abs(gammaToEPAngle)>=0.75*TMath::Pi()){
      if(fill&&fHistoEventPlanePhi)fHistoEventPlanePhi->Fill(gammaToEPAngle);
      return kTRUE;
    }
    else return kFALSE;
  }
  else if(fInPlaneOutOfPlane == 2){
    if(TMath::Abs(gammaToEPAngle)>0.25*TMath::Pi() && TMath::Abs(gammaToEPAngle)<0.75*TMath::Pi()){
      if(fill&&fHistoEventPlanePhi)fHistoEventPlanePhi->Fill(gammaToEPAngle);
      return kTRUE;
    }
    else return kFALSE;
  }
  return kFALSE;
}

///________________________________________________________________________
UChar_t AliConversionPhotonCuts::DeterminePhotonQualityAOD(AliAODConversionPhoton* photon, AliVEvent* eventDummy){

  AliAODTrack * negTrack = static_cast<AliAODTrack*>(GetTrack(eventDummy, photon->GetTrackLabelNegative()));
  AliAODTrack * posTrack = static_cast<AliAODTrack*>(GetTrack(eventDummy, photon->GetTrackLabelPositive()));

  if(!negTrack || !posTrack) {
      return 0;
  }
  if(negTrack->Charge() == posTrack->Charge()){
      return 0;
  }   
  Int_t nClusterITSneg = negTrack->GetITSNcls();
  Int_t nClusterITSpos = posTrack->GetITSNcls();
  //    cout << nClusterITSneg << "\t" << nClusterITSpos <<endl;
  
  if (nClusterITSneg > 1 && nClusterITSpos > 1){
    return 3;
  } else if (nClusterITSneg > 1 || nClusterITSpos > 1){
    return 2;
  } else {
    return 1;
  }
  return 0;
}

///__________________________________________________________________________________________
Bool_t AliConversionPhotonCuts::InitializeMaterialBudgetWeights(Int_t flag, TString filename){
    
    TString nameProfile;
    if      (flag==1){    
                nameProfile = "profileContainingMaterialBudgetWeights_fewRadialBins";}
    else if (flag==2){
                nameProfile = "profileContainingMaterialBudgetWeights_manyRadialBins";}
    else {
        AliError(Form("%d not a valid flag for InitMaterialBudgetWeightingOfPi0Candidates()",flag));
        return kFALSE;
    }
    TFile* file = TFile::Open(filename.Data());
    if (!file) {
        AliError(Form("File %s for materialbudgetweights not found",filename.Data()));
        return kFALSE;
    }
    fProfileContainingMaterialBudgetWeights = (TProfile*)file->Get(nameProfile.Data());
    if (!fProfileContainingMaterialBudgetWeights){
        AliError(Form("Histogram %s not found in file",nameProfile.Data()));
        return kFALSE;
    }
    fProfileContainingMaterialBudgetWeights->SetDirectory(0);
    file->Close();
    delete file;
    
    fMaterialBudgetWeightsInitialized = kTRUE;
    AliInfo(Form("MaterialBudgetWeightingOfPi0Candidates initialized with flag %d. This means %d radial bins will be used for the weighting. File used: %s.",flag, fProfileContainingMaterialBudgetWeights->GetNbinsX(), filename.Data()));
    return kTRUE;
}

///___________________________________________________________________________________________________
Float_t AliConversionPhotonCuts::GetMaterialBudgetCorrectingWeightForTrueGamma(AliAODConversionPhoton* gamma){
    
    Float_t weight = 1.0;
    Float_t gammaConversionRadius = gamma->GetConversionRadius();
    Int_t bin = fProfileContainingMaterialBudgetWeights->FindBin(gammaConversionRadius);
    if (bin > 0 && bin <= fProfileContainingMaterialBudgetWeights->GetNbinsX()){
        weight = fProfileContainingMaterialBudgetWeights->GetBinContent(bin);
    }
    return weight;
} 
