/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                           *
* Authors: Svein Lindal, Daniel Lohner            *
* Version 1.0                  *
*                    *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims    *
* about the suitability of this software for any purpose. It is    *
* provided "as is" without express or implied warranty.      *
**************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliConversionMesonCuts.h"

#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "AliMCEvent.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TPDGCode.h"
#include "TDatabasePDG.h"
#include "AliAODMCParticle.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"

class iostream;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliConversionMesonCuts)
/// \endcond


const char* AliConversionMesonCuts::fgkCutNames[AliConversionMesonCuts::kNCuts] = {
  "MesonKind", //0
  "BackgroundScheme", //1
  "NumberOfBGEvents", //2
  "DegreesForRotationMethod", //3
  "RapidityMesonCut", //4
  "PtCut", //5
  "AlphaMesonCut", //6
  "SelectionWindow", //7
  "SharedElectronCuts", //8
  "RejectToCloseV0s", //9
  "UseMCPSmearing", //10
  "DcaGammaGamma", //11
  "DcaRPrimVtx", //12
  "DcaZPrimVtx", //13
  "MinOpanMesonCut", //14
  "MaxOpanMesonCut" //15
};


//________________________________________________________________________
AliConversionMesonCuts::AliConversionMesonCuts(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fRandom(0),
  fCaloPhotonCuts(NULL),
  fHistograms(NULL),
  fCutString(NULL),
  fCutStringRead(""),
  fHistoMesonCuts(NULL),
  fHistoMesonBGCuts(NULL),
  fHistoDCAGGMesonBefore(NULL),
  fHistoDCAZMesonPrimVtxBefore(NULL),
  fHistoDCARMesonPrimVtxBefore(NULL),
  fHistoDCAGGMesonAfter(NULL),
  fHistoDCAZMesonPrimVtxAfter(NULL),
  fHistoDCARMesonPrimVtxAfter(NULL),
  fHistoInvMassBefore(NULL),
  fHistoInvMassAfter(NULL),
  fBrem(NULL),
  fFAlphaCut(NULL),
  fFMinOpanCut(NULL),
  fFMaxOpanCut(NULL),
  fMaxR(180),
  fMinPt(0.),
  fSelectionLow(0.0),
  fSelectionHigh(4),
  fSelectionNSigmaLow(999),
  fSelectionNSigmaHigh(999),
  fAlphaMinCutMeson(0),
  fAlphaCutMeson(1),
  fRapidityCutMesonMin(-1),
  fRapidityCutMesonMax(1),
  fMinV0Dist(200.),
  fAPlikeSigma(0),
  fMesonQualityMin(0),
  fPBremSmearing(0),
  fPSigSmearing(0),
  fPSigSmearingCte(0),
  fDCAGammaGammaCut(1000),
  fDCAZMesonPrimVtxCut(1000),
  fDCARMesonPrimVtxCut(1000),
  fMinOpanCutMeson(0),
  fMaxOpanCutMeson(TMath::Pi()),
  fSidebandMixingLow(0.180),
  fSidebandMixingHigh(0.300),
  fSidebandMixingLeftLow(0.05),
  fSidebandMixingLeftHigh(0.100),
  fSidebandMixingRightLow(0.180),
  fSidebandMixingRightHigh(0.300),
  fOpeningAngle(0.005),
  fMode(0),
  fMesonKind(0),
  fIsMergedClusterCut(0),
  fUsePtDepSelectionWindow(kFALSE),
  fUseGammaSelection(kFALSE),
  fSelectionWindowCut(-1),
  fNDegreeRotationPMForBG(0),
  fNumberOfBGEvents(0),
  fElectronLabelArraySize(500),
  fElectronLabelArray(NULL),
  fBackgroundHandler(0),
  fMassParamFunction(0),
  fDoLightOutput(0),
  fDoMinPtCut(kFALSE),
  fEnableMassCut(kFALSE),
  fAcceptMesonMass(kTRUE),
  fUseRotationMethodInBG(kFALSE),
  fUsePtmaxMethodForBG(kFALSE),
  fDoBG(kTRUE),
  fDoBGProbability(kFALSE),
  fDoConvCaloMixing(kFALSE),
  fDoSectorMixing(kFALSE),
  fDoSectorJetMixing(kFALSE),
  fDoJetMixing(kFALSE),
  fDoJetRotateMixing(kFALSE),
  fDoJetPtMixing(kFALSE),
  fDoSphericityMixing(kFALSE),
  fUseTrackMultiplicityForBG(kFALSE),
  fDoGammaSwappForBg(0),
  fDoWeightingInSwappBg(kFALSE),
  fGammaSwappMethodBg(0),
  fNumberOfSwappsForBg(1),
  fEnableMinOpeningAngleCut(kTRUE),
  fEnableOneCellDistCut(kFALSE),
  fAllowCombOnlyInSameRecMethod(kFALSE),
  fDoToCloseV0sCut(kFALSE),
  fDoSharedElecCut(kFALSE),
  fDoMesonQualitySelection(kFALSE),
  fUseMCPSmearing(0),
  fAlphaPtDepCut(kFALSE),
  fDCAGammaGammaCutOn(kFALSE),
  fDCAZMesonPrimVtxCutOn(kFALSE),
  fDCARMesonPrimVtxCutOn(kFALSE),
  fMinOpanPtDepCut(kFALSE),
  fMaxOpanPtDepCut(kFALSE),
  fBackgroundUseSideband(kFALSE),
  fBackgroundUseSidebandBothSides(kFALSE),
  fBackgroundUseLikeSign(kFALSE),
  fBackgroundMode(4),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fDoOutOfJet(0),
  fDoIsolatedAnalysis(kFALSE),
  fDoHighPtHadronAnalysis(kFALSE),
  fEnableOmegaAPlikeCut(kFALSE),
  fDoGammaMinEnergyCut(kFALSE),
  fNDaughterEnergyCut(0),
  fSingleDaughterMinE(0.)
{
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());
  fElectronLabelArray = new Int_t[fElectronLabelArraySize];
  if (fBrem == NULL){
    fBrem = new TF1("fBrem","pow(-log(x),[0]/log(2.0)-1.0)/TMath::Gamma([0]/log(2.0))",0.00001,0.999999999);
    // tests done with 1.0e-14
    fBrem->SetParameter(0,fPBremSmearing);
    fBrem->SetNpx(100000);
  }
}

//________________________________________________________________________
AliConversionMesonCuts::AliConversionMesonCuts(const AliConversionMesonCuts &ref) :
  AliAnalysisCuts(ref),
  fRandom(ref.fRandom),
  fCaloPhotonCuts(ref.fCaloPhotonCuts),
  fHistograms(NULL),
  fCutString(NULL),
  fCutStringRead(""),
  fHistoMesonCuts(NULL),
  fHistoMesonBGCuts(NULL),
  fHistoDCAGGMesonBefore(NULL),
  fHistoDCAZMesonPrimVtxBefore(NULL),
  fHistoDCARMesonPrimVtxBefore(NULL),
  fHistoDCAGGMesonAfter(NULL),
  fHistoDCAZMesonPrimVtxAfter(NULL),
  fHistoDCARMesonPrimVtxAfter(NULL),
  fHistoInvMassBefore(NULL),
  fHistoInvMassAfter(NULL),
  fBrem(NULL),
  fFAlphaCut(NULL),
  fFMinOpanCut(NULL),
  fFMaxOpanCut(NULL),
  fMaxR(ref.fMaxR),
  fMinPt(ref.fMinPt),
  fSelectionLow(ref.fSelectionLow),
  fSelectionHigh(ref.fSelectionHigh),
  fSelectionNSigmaLow(ref.fSelectionNSigmaLow),
  fSelectionNSigmaHigh(ref.fSelectionNSigmaHigh),
  fAlphaMinCutMeson(ref.fAlphaMinCutMeson),
  fAlphaCutMeson(ref.fAlphaCutMeson),
  fRapidityCutMesonMin(ref.fRapidityCutMesonMin),
  fRapidityCutMesonMax(ref.fRapidityCutMesonMax),
  fMinV0Dist(ref.fMinV0Dist),
  fAPlikeSigma(ref.fAPlikeSigma),
  fMesonQualityMin(ref.fMesonQualityMin),
  fPBremSmearing(ref.fPBremSmearing),
  fPSigSmearing(ref.fPSigSmearing),
  fPSigSmearingCte(ref.fPSigSmearingCte),
  fDCAGammaGammaCut(ref.fDCAGammaGammaCut),
  fDCAZMesonPrimVtxCut(ref.fDCAZMesonPrimVtxCut),
  fDCARMesonPrimVtxCut(ref.fDCARMesonPrimVtxCut),
  fMinOpanCutMeson(ref.fMinOpanCutMeson),
  fMaxOpanCutMeson(ref.fMaxOpanCutMeson),
  fSidebandMixingLow(ref.fSidebandMixingLow),
  fSidebandMixingHigh(ref.fSidebandMixingHigh),
  fSidebandMixingLeftLow(ref.fSidebandMixingLeftLow),
  fSidebandMixingLeftHigh(ref.fSidebandMixingLeftHigh),
  fSidebandMixingRightLow(ref.fSidebandMixingRightLow),
  fSidebandMixingRightHigh(ref.fSidebandMixingRightHigh),
  fOpeningAngle(ref.fOpeningAngle),
  fMode(ref.fMode),
  fMesonKind(ref.fMesonKind),
  fIsMergedClusterCut(ref.fIsMergedClusterCut),
  fUsePtDepSelectionWindow(ref.fUsePtDepSelectionWindow),
  fUseGammaSelection(ref.fUseGammaSelection),
  fSelectionWindowCut(ref.fSelectionWindowCut),
  fNDegreeRotationPMForBG(ref.fNDegreeRotationPMForBG),
  fNumberOfBGEvents(ref. fNumberOfBGEvents),
  fElectronLabelArraySize(ref.fElectronLabelArraySize),
  fElectronLabelArray(NULL),
  fBackgroundHandler(ref.fBackgroundHandler),
  fMassParamFunction(ref.fMassParamFunction),
  fDoLightOutput(ref.fDoLightOutput),
  fDoMinPtCut(ref.fDoMinPtCut),
  fEnableMassCut(ref.fEnableMassCut),
  fAcceptMesonMass(ref.fAcceptMesonMass),
  fUseRotationMethodInBG(ref.fUseRotationMethodInBG),
  fUsePtmaxMethodForBG(ref.fUsePtmaxMethodForBG),
  fDoBG(ref.fDoBG),
  fDoBGProbability(ref.fDoBGProbability),
  fDoConvCaloMixing(ref.fDoConvCaloMixing),
  fDoSectorMixing(ref.fDoSectorMixing),
  fDoSectorJetMixing(ref.fDoSectorJetMixing),
  fDoJetMixing(ref.fDoJetMixing),
  fDoJetRotateMixing(ref.fDoJetRotateMixing),
  fDoJetPtMixing(ref.fDoJetPtMixing),
  fDoSphericityMixing(ref.fDoSphericityMixing),
  fUseTrackMultiplicityForBG(ref.fUseTrackMultiplicityForBG),
  fDoGammaSwappForBg(ref.fDoGammaSwappForBg),
  fDoWeightingInSwappBg(ref.fDoWeightingInSwappBg),
  fGammaSwappMethodBg(ref.fGammaSwappMethodBg),
  fNumberOfSwappsForBg(ref.fNumberOfSwappsForBg),
  fEnableMinOpeningAngleCut(ref.fEnableMinOpeningAngleCut),
  fEnableOneCellDistCut(ref.fEnableOneCellDistCut),
  fAllowCombOnlyInSameRecMethod(ref.fAllowCombOnlyInSameRecMethod),
  fDoToCloseV0sCut(ref.fDoToCloseV0sCut),
  fDoSharedElecCut(ref.fDoSharedElecCut),
  fDoMesonQualitySelection(ref.fDoMesonQualitySelection),
  fUseMCPSmearing(ref.fUseMCPSmearing),
  fAlphaPtDepCut(ref.fAlphaPtDepCut),
  fDCAGammaGammaCutOn(ref.fDCAGammaGammaCutOn),
  fDCAZMesonPrimVtxCutOn(ref.fDCAZMesonPrimVtxCutOn),
  fDCARMesonPrimVtxCutOn(ref.fDCARMesonPrimVtxCutOn),
  fMinOpanPtDepCut(ref.fMinOpanPtDepCut),
  fMaxOpanPtDepCut(ref.fMaxOpanPtDepCut),
  fBackgroundUseSideband(ref.fBackgroundUseSideband),
  fBackgroundUseSidebandBothSides(ref.fBackgroundUseSidebandBothSides),
  fBackgroundUseLikeSign(ref.fBackgroundUseLikeSign),
  fBackgroundMode(4),
  fDoJetAnalysis(ref.fDoJetAnalysis),
  fDoJetQA(ref.fDoJetQA),
  fDoOutOfJet(ref.fDoOutOfJet),
  fDoIsolatedAnalysis(ref.fDoIsolatedAnalysis),
  fDoHighPtHadronAnalysis(ref.fDoHighPtHadronAnalysis),
  fEnableOmegaAPlikeCut(ref.fEnableOmegaAPlikeCut),
  fDoGammaMinEnergyCut(kFALSE),
  fNDaughterEnergyCut(0),
  fSingleDaughterMinE(0.)

{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());
  fElectronLabelArray = new Int_t[fElectronLabelArraySize];
  if (fBrem == NULL)fBrem = (TF1*)ref.fBrem->Clone("fBrem");
  // Histograms are not copied, if you need them, call InitCutHistograms
}


//________________________________________________________________________
AliConversionMesonCuts::~AliConversionMesonCuts() {
  // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  // if(fHistograms)
  //   delete fHistograms;
  // fHistograms = NULL;
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }
  if(fElectronLabelArray){
    delete fElectronLabelArray;
    fElectronLabelArray = NULL;
  }

  if(fFAlphaCut != NULL){
    delete fFAlphaCut;
    fFAlphaCut = NULL;
  }
  if(fBrem != NULL){
    delete fBrem;
    fBrem = NULL;
  }
  if(fFMinOpanCut != NULL){
    delete fFMinOpanCut;
    fFMinOpanCut = NULL;
  }
  if(fFMaxOpanCut != NULL){
    delete fFMaxOpanCut;
    fFMaxOpanCut = NULL;
  }
}

//________________________________________________________________________
void AliConversionMesonCuts::InitCutHistograms(TString name, Bool_t additionalHists){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);

  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }

  if(fDoLightOutput==2) {
      AliInfo("Minimal output chosen");
      return;
  }

  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("ConvMesonCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  // Meson Cuts
  if (fIsMergedClusterCut == 1){
    fHistoMesonCuts=new TH2F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts vs Pt",9,-0.5,8.5, 500, 0, 100);
    fHistoMesonCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(4,"mass cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(5,"opening angle");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(6,"alpha max");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(7,"alpha min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(8,"pT min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(9,"out");
    fHistograms->Add(fHistoMesonCuts);
  } else if (fIsMergedClusterCut == 2){
    fHistoMesonCuts=new TH2F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts vs Pt",9,-0.5,8.5, 250, 0, 50);
    fHistoMesonCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(4,"1 cell distance");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(5,"opening angle");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(6,"alpha max");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(7,"alpha min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(8,"pT min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(9,"out");
    fHistograms->Add(fHistoMesonCuts);

    fHistoMesonBGCuts=new TH2F(Form("MesonBGCuts %s",GetCutNumber().Data()),"MesonBGCuts vs Pt",9,-0.5,8.5, 250, 0, 50);
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(4,"1 cell distance");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(5,"opening angle");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(6,"alpha max");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(7,"alpha min");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(8,"pT min");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(9,"out");
    fHistograms->Add(fHistoMesonBGCuts);
  } else {
    fHistoMesonCuts=new TH2F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts vs Pt",12,-0.5,11.5, 250, 0, 50);
    fHistoMesonCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(4,"opening angle");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(5,"alpha max");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(6,"alpha min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(7,"dca gamma gamma");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(8,"dca R prim Vtx");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(9,"dca Z prim Vtx");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(10,"pT min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(11,"Meson Quality min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(12,"out");
    fHistograms->Add(fHistoMesonCuts);

    fHistoMesonBGCuts=new TH2F(Form("MesonBGCuts %s",GetCutNumber().Data()),"MesonBGCuts vs Pt",12,-0.5,11.5, 250, 0, 50);
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(4,"opening angle");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(5,"alpha max");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(6,"alpha min");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(7,"dca gamma gamma");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(8,"dca R prim Vtx");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(9,"dca Z prim Vtx");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(10,"pT min");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(11,"Meson Quality min");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(12,"out");
    fHistograms->Add(fHistoMesonBGCuts);
  }

  if(!fDoLightOutput){
    if (fIsMergedClusterCut == 1){
      fHistoInvMassBefore=new TH1F(Form("InvMassMeson Before %s",GetCutNumber().Data()),"InvMassMeson Before",1000,0,1);
      fHistograms->Add(fHistoInvMassBefore);
      fHistoInvMassAfter=new TH1F(Form("InvMassMeson After %s",GetCutNumber().Data()),"InvMassMeson After",1000,0,1);
      fHistograms->Add(fHistoInvMassAfter);
    }

    if (additionalHists && fIsMergedClusterCut== 0){
      fHistoDCAGGMesonBefore=new TH1F(Form("DCAGammaGammaMeson Before %s",GetCutNumber().Data()),"DCAGammaGammaMeson Before",200,0,10);
      fHistograms->Add(fHistoDCAGGMesonBefore);

      fHistoDCARMesonPrimVtxBefore=new TH1F(Form("DCARMesonPrimVtx Before %s",GetCutNumber().Data()),"DCARMesonPrimVtx Before",200,0,10);
      fHistograms->Add(fHistoDCARMesonPrimVtxBefore);

      fHistoDCAZMesonPrimVtxBefore=new TH1F(Form("DCAZMesonPrimVtx Before %s",GetCutNumber().Data()),"DCAZMesonPrimVtx Before",401,-10,10);
      fHistograms->Add(fHistoDCAZMesonPrimVtxBefore);
    }

    if (fIsMergedClusterCut == 0){
      fHistoDCAGGMesonAfter=new TH1F(Form("DCAGammaGammaMeson After %s",GetCutNumber().Data()),"DCAGammaGammaMeson After",200,0,10);
      fHistograms->Add(fHistoDCAGGMesonAfter);

      fHistoDCAZMesonPrimVtxAfter=new TH2F(Form("InvMassDCAZMesonPrimVtx After %s",GetCutNumber().Data()),"InvMassDCAZMesonPrimVtx After",800,0,0.8,401,-10,10);
      fHistograms->Add(fHistoDCAZMesonPrimVtxAfter);

      fHistoDCARMesonPrimVtxAfter=new TH1F(Form("DCARMesonPrimVtx After %s",GetCutNumber().Data()),"DCARMesonPrimVtx After",200,0,10);
      fHistograms->Add(fHistoDCARMesonPrimVtxAfter);
    }
  }

  TH1::AddDirectory(kTRUE);
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMC(TParticle *fMCMother,AliMCEvent *mcEvent, Double_t fRapidityShift){
  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if(!mcEvent)return kFALSE;

  if(fMCMother->GetPdgCode()==111 || fMCMother->GetPdgCode()==221 || fMCMother->GetPdgCode()==331 ){
    if(fMCMother->R()>fMaxR)  return kFALSE; // cuts on distance from collision point

    Double_t rapidity = 10.;
    if(fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0){
      rapidity=8.-fRapidityShift;
    } else{
      rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
    }

    // Rapidity Cut
    if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

    // min Pt Cut
    if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

    // Select only -> 2y decay channel
    if(fMCMother->GetNDaughters()!=2)return kFALSE;

    for(Int_t i=0;i<2;i++){
      if(fMCMother->GetDaughter(i) < 0) return kFALSE;
      TParticle *MDaughter=mcEvent->Particle(fMCMother->GetDaughter(i));
      // Is Daughter a Photon?
      if(MDaughter->GetPdgCode()!=22)return kFALSE;
      // Is Photon in Acceptance?
      //   if(bMCDaughtersInAcceptance){
      //  if(!PhotonIsSelectedMC(MDaughter,mcEvent)){return kFALSE;}
      //   }
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMC(AliAODMCParticle *MCMother,TClonesArray *AODMCArray, Double_t fRapidityShift){
  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if(!AODMCArray)return kFALSE;

  if(MCMother->GetPdgCode()==111 || MCMother->GetPdgCode()==221 || MCMother->GetPdgCode()==331 ){
    Double_t rMeson = sqrt( (MCMother->Xv()*MCMother->Xv()) + (MCMother->Yv()*MCMother->Yv()) ) ;
    if(rMeson>fMaxR)  return kFALSE; // cuts on distance from collision point

    Double_t rapidity = 10.;
    if(MCMother->E() - MCMother->Pz() == 0 || MCMother->E() + MCMother->Pz() == 0){
      rapidity=8.-fRapidityShift;
    } else{
      rapidity = 0.5*(TMath::Log((MCMother->E()+MCMother->Pz()) / (MCMother->E()-MCMother->Pz())))-fRapidityShift;
    }

    // Rapidity Cut
    if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

    // min Pt Cut
    if(fDoMinPtCut && (MCMother->Pt() < fMinPt)) return kFALSE;

    // Select only -> 2y decay channel
    if(MCMother->GetNDaughters()!=2)return kFALSE;

    for(Int_t i=0;i<2;i++){
      AliAODMCParticle *MDaughter=static_cast<AliAODMCParticle*>(AODMCArray->At(MCMother->GetDaughterLabel(i)));
      // Is Daughter a Photon?
      if(MDaughter->GetPdgCode()!=22)return kFALSE;
      // Is Photon in Acceptance?
      //   if(bMCDaughtersInAcceptance){
      //  if(!PhotonIsSelectedMC(MDaughter,mcEvent)){return kFALSE;}
      //   }
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCAODESD(AliDalitzAODESDMC *fMCMother,AliDalitzEventMC *mcEvent, Double_t fRapidityShift) const{
// Returns true for all pions within acceptance cuts for decay into 2 photons
// If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts
    if(!mcEvent)return kFALSE;

    if(fMCMother->GetPdgCodeG()==111 || fMCMother->GetPdgCodeG()==221 || fMCMother->GetPdgCodeG()==331 ){
        if(fMCMother->GetRatioVxyG()>fMaxR)  return kFALSE; // cuts on distance from collision point
        Double_t rapidity = 10.;
        if(fMCMother->EnergyG() - fMCMother->PzG() == 0 || fMCMother->EnergyG() + fMCMother->PzG() == 0){
            rapidity=8.-fRapidityShift;
        } else{
            rapidity = 0.5*(TMath::Log((fMCMother->EnergyG()+fMCMother->PzG()) / (fMCMother->EnergyG()-fMCMother->PzG())))-fRapidityShift;
        }

        // Rapidity Cut
        if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

        // min Pt Cut
        if(fDoMinPtCut && (fMCMother->PtG() < fMinPt)) return kFALSE;

        // Select only -> 2y decay channel
        if(fMCMother->GetNDaughtersG()!=2)return kFALSE;
        if((fMCMother->GetLastDaughterG()<0) || (fMCMother->GetFirstDaughterG()<0)) return kFALSE;
            std::unique_ptr<AliDalitzAODESDMC> MDaughter0=std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(fMCMother->GetFirstDaughterG()));
            std::unique_ptr<AliDalitzAODESDMC> MDaughter1=std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(fMCMother->GetLastDaughterG()));
        if ((MDaughter0->GetPdgCodeG()!=22) || (MDaughter1->GetPdgCodeG()!=22)) return kFALSE;
            return kTRUE;
    }
    return kFALSE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCDalitz(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !mcEvent )return kFALSE;

  if(  fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 && fMCMother->GetPdgCode() != 331) return kFALSE;

  if(  fMCMother->R()>fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> Dalitz decay channel
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  TParticle *positron = 0x0;
  TParticle *electron = 0x0;
  TParticle    *gamma = 0x0;

  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
    if(index < 0) continue;
    TParticle* temp = (TParticle*)mcEvent->Particle( index );

    switch( temp->GetPdgCode() ) {
    case ::kPositron:
      positron      =  temp;
      labelpositron = index;
      break;
    case ::kElectron:
      electron      =  temp;
      labelelectron = index;
      break;
    case ::kGamma:
      gamma         =  temp;
      labelgamma    = index;
      break;
    }
  }

  if( positron && electron && gamma) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMCDalitz(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !AODMCArray )return kFALSE;

  if(  fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 && fMCMother->GetPdgCode() != 331 ) return kFALSE;

  Double_t rMeson = sqrt( (fMCMother->Xv()*fMCMother->Xv()) + (fMCMother->Yv()*fMCMother->Yv()) ) ;
  if(rMeson>fMaxR)  return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->E() - fMCMother->Pz() == 0 || fMCMother->E() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->E()+fMCMother->Pz()) / (fMCMother->E()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> Dalitz decay channel
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  AliAODMCParticle *positron = 0x0;
  AliAODMCParticle *electron = 0x0;
  AliAODMCParticle    *gamma = 0x0;

  for(Int_t index= fMCMother->GetDaughterFirst();index<= fMCMother->GetDaughterLast();index++){
    if(index < 0) continue;
    AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(AODMCArray->At(index));
    if (!temp) continue;

    switch( temp->GetPdgCode() ) {
    case ::kPositron:
      positron      =  temp;
      labelpositron = index;
      break;
    case ::kElectron:
      electron      =  temp;
      labelelectron = index;
      break;
    case ::kGamma:
      gamma         =  temp;
      labelgamma    = index;
      break;
    }
  }

  if( positron && electron && gamma) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
 Bool_t AliConversionMesonCuts::MesonIsSelectedMCDalitzAODESD(AliDalitzAODESDMC* fMCMother,AliDalitzEventMC *mcEvent, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma, Double_t fRapidityShift) const{
// Returns true for all pions within acceptance cuts for decay into 2 photons
// If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts
    if( !mcEvent )return kFALSE;

    if(  fMCMother->GetPdgCodeG() != 111 && fMCMother->GetPdgCodeG() != 221 && fMCMother->GetPdgCodeG() != 331) return kFALSE;

    if(  fMCMother->GetRatioVxyG()>fMaxR ) return kFALSE; // cuts on distance from collision point

    Double_t rapidity = 10.;

    if( fMCMother->EnergyG() - fMCMother->PzG() == 0 || fMCMother->EnergyG() + fMCMother->PzG() == 0 ){
        rapidity=8.-fRapidityShift;
    }
    else{
        rapidity = 0.5*(TMath::Log((fMCMother->EnergyG()+fMCMother->PzG()) / (fMCMother->EnergyG()-fMCMother->PzG())))-fRapidityShift;
    }

    // Rapidity Cut
    if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

    // min Pt Cut
    if(fDoMinPtCut && (fMCMother->PtG() < fMinPt)) return kFALSE;

    // Select only -> Dalitz decay channel
    if( fMCMother->GetNDaughtersG() != 3 )return kFALSE;

    Bool_t flagpositron=kFALSE;
    Bool_t flagelectron=kFALSE;
    Bool_t flaggamma=kFALSE;
    for(Int_t index= fMCMother->GetFirstDaughterG();index<= fMCMother->GetLastDaughterG();index++){
        if(index < 0) continue;
        std::unique_ptr<AliDalitzAODESDMC> temp = std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(index));
        switch( temp->GetPdgCodeG() ) {
        case ::kPositron:{
            flagpositron=kTRUE;
            labelpositron = index; }
            break;
        case ::kElectron:{
            flagelectron=kTRUE;
            labelelectron = index;}
            break;
        case ::kGamma:{
            flaggamma=kTRUE;
            labelgamma    = index;}
            break;
        }
    }
    if( flagelectron && flagpositron && flaggamma) return kTRUE;
    return kFALSE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCEtaPiPlPiMiGamma(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelGamma, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !mcEvent )return kFALSE;

  if( fMCMother->GetPdgCode() != 221 ) return kFALSE;

  if( fMCMother->R()>fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> Dalitz decay channel
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  TParticle *posPion = 0x0;
  TParticle *negPion = 0x0;
  TParticle    *gamma = 0x0;

  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
    if(index < 0) continue;
    TParticle* temp = (TParticle*)mcEvent->Particle( index );

    switch( temp->GetPdgCode() ) {
    case 211:
      posPion      =  temp;
      labelPosPion = index;
      break;
    case -211:
      negPion      =  temp;
      labelNegPion = index;
      break;
    case ::kGamma:
      gamma         =  temp;
      labelGamma    = index;
      break;
    }
  }

  if( posPion && negPion && gamma) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMCEtaPiPlPiMiGamma(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelGamma, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !AODMCArray )return kFALSE;

  if( fMCMother->GetPdgCode() != 221 ) return kFALSE;

  Double_t rMeson = sqrt( (fMCMother->Xv()*fMCMother->Xv()) + (fMCMother->Yv()*fMCMother->Yv()) ) ;
  if(rMeson>fMaxR)  return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->E() - fMCMother->Pz() == 0 || fMCMother->E() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->E()+fMCMother->Pz()) / (fMCMother->E()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> Dalitz decay channel
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  AliAODMCParticle *posPion = 0x0;
  AliAODMCParticle *negPion = 0x0;
  AliAODMCParticle    *gamma = 0x0;

  for(Int_t index= fMCMother->GetDaughterFirst();index<= fMCMother->GetDaughterLast();index++){
    if(index < 0) continue;
    AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(AODMCArray->At(index));

    switch( temp->GetPdgCode() ) {
    case 211:
      posPion      =  temp;
      labelPosPion = index;
      break;
    case -211:
      negPion      =  temp;
      labelNegPion = index;
      break;
    case ::kGamma:
      gamma         =  temp;
      labelGamma    = index;
      break;
    }
  }

  if( posPion && negPion && gamma) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCPiPlPiMiEta(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelEtaMeson, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !mcEvent )return kFALSE;

  if( !(fMCMother->GetPdgCode() == 331 ) ) return kFALSE;

  if( fMCMother->R()>fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> pi+ pi- pi0
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  TParticle *posPion = 0x0;
  TParticle *negPion = 0x0;
  TParticle *etaMeson = 0x0;

//   cout << "\n"<< fMCMother->GetPdgCode() << "\n" << endl;
  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
    if(index < 0) continue;
    TParticle* temp = (TParticle*)mcEvent->Particle( index );
//     cout << temp->GetPdgCode() << endl;
    switch( temp->GetPdgCode() ) {
    case 211:
      posPion      =  temp;
      labelPosPion = index;
      break;
    case -211:
      negPion      =  temp;
      labelNegPion = index;
      break;
    case 221:
      etaMeson         =  temp;
      labelEtaMeson    = index;
      break;
    }
  }

  if( posPion && negPion && etaMeson ) return kTRUE;
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMCPiPlPiMiEta(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelEtaMeson, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !AODMCArray )return kFALSE;

  if( !(fMCMother->GetPdgCode() == 331 ) ) return kFALSE;

  Double_t rMeson = sqrt( (fMCMother->Xv()*fMCMother->Xv()) + (fMCMother->Yv()*fMCMother->Yv()) ) ;
  if( rMeson >fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->E() - fMCMother->Pz() == 0 || fMCMother->E() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->E()+fMCMother->Pz()) / (fMCMother->E()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> pi+ pi- pi0
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  AliAODMCParticle *posPion = 0x0;
  AliAODMCParticle *negPion = 0x0;
  AliAODMCParticle *etaMeson = 0x0;

//   cout << "\n"<< fMCMother->GetPdgCode() << "\n" << endl;
  for(Int_t index= fMCMother->GetDaughterFirst();index<= fMCMother->GetDaughterLast();index++){
    if(index < 0) continue;
    AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(AODMCArray->At(index));
//     cout << temp->GetPdgCode() << endl;
    switch( temp->GetPdgCode() ) {
    case 211:
      posPion      =  temp;
      labelPosPion = index;
      break;
    case -211:
      negPion      =  temp;
      labelNegPion = index;
      break;
    case 221:
      etaMeson         =  temp;
      labelEtaMeson    = index;
      break;
    }
  }

  if( posPion && negPion && etaMeson ) return kTRUE;
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCPiPlPiMiPiZero(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !mcEvent )return kFALSE;

  if( !(fMCMother->GetPdgCode() == 221 || fMCMother->GetPdgCode() == 223 || fMCMother->GetPdgCode() == 421) ) return kFALSE;

  if( fMCMother->R()>fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> pi+ pi- pi0
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  TParticle *posPion = 0x0;
  TParticle *negPion = 0x0;
  TParticle *neutPion = 0x0;

//   cout << "\n"<< fMCMother->GetPdgCode() << "\n" << endl;
  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
    if(index < 0) continue;
    TParticle* temp = (TParticle*)mcEvent->Particle( index );
//     cout << temp->GetPdgCode() << endl;
    switch( temp->GetPdgCode() ) {
    case 211:
      posPion      =  temp;
      labelPosPion = index;
      break;
    case -211:
      negPion      =  temp;
      labelNegPion = index;
      break;
    case 111:
      neutPion         =  temp;
      labelNeutPion    = index;
      break;
    }
  }

  if( posPion && negPion && neutPion ) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMCPiPlPiMiPiZero(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !AODMCArray )return kFALSE;

  if( !(fMCMother->GetPdgCode() == 221 || fMCMother->GetPdgCode() == 223 || fMCMother->GetPdgCode() == 421) ) return kFALSE;

  Double_t rMeson = sqrt( (fMCMother->Xv()*fMCMother->Xv()) + (fMCMother->Yv()*fMCMother->Yv()) ) ;
  if( rMeson >fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->E() - fMCMother->Pz() == 0 || fMCMother->E() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->E()+fMCMother->Pz()) / (fMCMother->E()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  // Select only -> pi+ pi- pi0
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  AliAODMCParticle *posPion  = 0x0;
  AliAODMCParticle *negPion  = 0x0;
  AliAODMCParticle *neutPion = 0x0;

//   cout << "\n"<< fMCMother->GetPdgCode() << "\n" << endl;
  for(Int_t index= fMCMother->GetDaughterFirst();index<= fMCMother->GetDaughterLast();index++){
    if(index < 0) continue;
    AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(AODMCArray->At(index));
//     cout << temp->GetPdgCode() << endl;
    switch( temp->GetPdgCode() ) {
    case 211:
      posPion      =  temp;
      labelPosPion = index;
      break;
    case -211:
      negPion      =  temp;
      labelNegPion = index;
      break;
    case 111:
      neutPion         =  temp;
      labelNeutPion    = index;
      break;
    }
  }

  if( posPion && negPion && neutPion ) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCPiZeroGamma(TParticle *fMCMother, AliMCEvent *mcEvent, Int_t &labelNeutPion, Int_t &labelGamma, Double_t fRapidityShift){
  // returns true for omegas decaying into pi0 + gamma within the rapidity window

  if(!mcEvent) return kFALSE;

  if(fMCMother->GetPdgCode()!=223) return kFALSE; // we only want omegas

  Double_t rapidity = 10.;

  if(fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  if(fMCMother->GetNDaughters()!=2) return kFALSE;

  TParticle *gamma = 0x0;
  TParticle *pi0 = 0x0;

  for(Int_t index = fMCMother->GetFirstDaughter();index <= fMCMother->GetLastDaughter();index++){
    if(index < 0) continue;
    TParticle *temp = (TParticle*)mcEvent->Particle(index);
    switch(temp->GetPdgCode()){
    case 22:
      gamma = temp;
      labelGamma = index;
      break;
    case 111:
      pi0   = temp;
      labelNeutPion = index;
      break;
    }
  }

  if(gamma && pi0) return kTRUE;
  return kFALSE;

}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMCPiZeroGamma(AliAODMCParticle *fMCMother, TClonesArray *AODMCArray, Int_t &labelNeutPion, Int_t &labelGamma, Double_t fRapidityShift){
  // returns true for omegas decaying into pi0 + gamma within the rapidity window

  if(!AODMCArray) return kFALSE;

  if(fMCMother->GetPdgCode()!=223) return kFALSE; // we only want omegas

  Double_t rapidity = 10.;

  if(fMCMother->E() - fMCMother->Pz() == 0 || fMCMother->E() + fMCMother->Pz() == 0){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->E()+fMCMother->Pz()) / (fMCMother->E()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

  // min Pt Cut
  if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

  if(fMCMother->GetNDaughters()!=2) return kFALSE;

  AliAODMCParticle *gamma = 0x0;
  AliAODMCParticle *pi0 = 0x0;

  for(Int_t index = fMCMother->GetDaughterFirst();index <= fMCMother->GetDaughterLast();index++){
    if(index < 0) continue;
    AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(AODMCArray->At(index));
    switch(temp->GetPdgCode()){
    case 22:
      gamma = temp;
      labelGamma = index;
      break;
    case 111:
      pi0   = temp;
      labelNeutPion = index;
      break;
    }
  }

  if(gamma && pi0) return kTRUE;
  return kFALSE;

}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedPiZeroGammaAngle(AliAODConversionMother *omega, AliAODConversionMother *pi0, AliAODConversionPhoton *gamma, Bool_t DoPiZeroAngleCut, TF1 *maxfit, Double_t lowerFactor, Double_t upperFactor){

  if(!DoPiZeroAngleCut) return kTRUE;

  Double_t PiZeroGammaAngle = pi0->Angle(gamma->Vect());
  Double_t omegaPt = omega->Pt();

  if(PiZeroGammaAngle > lowerFactor * maxfit->Eval(omegaPt) && PiZeroGammaAngle < upperFactor * maxfit->Eval(omegaPt)) return kTRUE;
  return kFALSE;

}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedPiZeroGammaOAC(AliAODConversionMother *omega, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1, AliAODConversionPhoton *gamma2){

  TVector3 vecPi0Gamma0  = TVector3(gamma0->Px(), gamma0->Py(), gamma0->Pz());
  TVector3 vecPi0Gamma1  = TVector3(gamma1->Px(), gamma1->Py(), gamma1->Pz());
  TVector3 vecOmegaGamma = TVector3(gamma2->Px(), gamma2->Py(), gamma2->Pz());

  // Opening Angle Cut
  if( ( fEnableMinOpeningAngleCut ) && (vecPi0Gamma0.Angle(vecPi0Gamma1) < fOpeningAngle) &&
    (vecPi0Gamma0.Angle(vecOmegaGamma) < fOpeningAngle) && (vecPi0Gamma1.Angle(vecOmegaGamma) < fOpeningAngle) )
  {
    return kFALSE;
  }

  // Min Opening Angle
  if (fMinOpanPtDepCut == kTRUE) fMinOpanCutMeson = fFMinOpanCut->Eval(omega->Pt());

  if( (vecPi0Gamma0.Angle(vecPi0Gamma1) < fMinOpanCutMeson) &&
    (vecPi0Gamma0.Angle(vecOmegaGamma) < fMinOpanCutMeson) && (vecPi0Gamma1.Angle(vecOmegaGamma) < fMinOpanCutMeson) )
  {
    return kFALSE;
  }

  // Max Opening Angle
  if (fMaxOpanPtDepCut == kTRUE) fMaxOpanCutMeson = fFMaxOpanCut->Eval(omega->Pt());

  if( (vecPi0Gamma0.Angle(vecPi0Gamma1) > fMaxOpanCutMeson) &&
    (vecPi0Gamma0.Angle(vecOmegaGamma) > fMaxOpanCutMeson) && (vecPi0Gamma1.Angle(vecOmegaGamma) > fMaxOpanCutMeson) )
  {
    return kFALSE;
  }

  return kTRUE;

}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCChiC(TParticle *fMCMother,AliMCEvent *mcEvent,Int_t & labelelectronChiC, Int_t & labelpositronChiC, Int_t & labelgammaChiC, Double_t fRapidityShift){
  // Returns true for all ChiC within acceptance cuts for decay into JPsi + gamma -> e+ + e- + gamma
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if(!mcEvent)return kFALSE;
  // if(fMCMother->GetPdgCode()==20443 ){
  //    return kFALSE;
  // }
  if(fMCMother->GetPdgCode()==10441 || fMCMother->GetPdgCode()==10443 || fMCMother->GetPdgCode()==445 ){
    if(fMCMother->R()>fMaxR)  return kFALSE; // cuts on distance from collision point

    Double_t rapidity = 10.;
    if(fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0){
      rapidity=8.-fRapidityShift;
    }
    else{
      rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
    }

    // Rapidity Cut
    if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

    // min Pt Cut
    if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

    // Select only -> ChiC radiative (JPsi+gamma) decay channel
    if(fMCMother->GetNDaughters()!=2)return kFALSE;

    TParticle *jpsi   = 0x0;
    TParticle *gamma   = 0x0;
    TParticle *positron = 0x0;
    TParticle *electron = 0x0;

    //Int_t labeljpsiChiC = -1;

    for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
      if(index < 0) continue;
      TParticle* temp = (TParticle*)mcEvent->Particle( index );

      switch( temp->GetPdgCode() ) {
      case 443:
        jpsi =  temp;
        //labeljpsiChiC = index;
        break;
      case 22:
        gamma    =  temp;
        labelgammaChiC = index;
        break;
      }
    }

    if ( !jpsi || ! gamma) return kFALSE;
    if(jpsi->GetNDaughters()!=2)return kFALSE;


    for(Int_t index= jpsi->GetFirstDaughter();index<= jpsi->GetLastDaughter();index++){
      if(index < 0) continue;
      TParticle* temp = (TParticle*)mcEvent->Particle( index );
      switch( temp->GetPdgCode() ) {
      case -11:
        electron =  temp;
        labelelectronChiC = index;
        break;
      case 11:
        positron =  temp;
        labelpositronChiC = index;
        break;
      }
    }
    if( !electron || !positron) return kFALSE;
    if( positron && electron && gamma) return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedAODMCChiC(AliAODMCParticle *fMCMother,TClonesArray *AODMCArray,Int_t & labelelectronChiC, Int_t & labelpositronChiC, Int_t & labelgammaChiC, Double_t fRapidityShift){
  // Returns true for all ChiC within acceptance cuts for decay into JPsi + gamma -> e+ + e- + gamma
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if(!AODMCArray)return kFALSE;
  // if(fMCMother->GetPdgCode()==20443 ){
  //    return kFALSE;
  // }
  if(fMCMother->GetPdgCode()==10441 || fMCMother->GetPdgCode()==10443 || fMCMother->GetPdgCode()==445 ){
    Double_t rMeson = sqrt( (fMCMother->Xv()*fMCMother->Xv()) + (fMCMother->Yv()*fMCMother->Yv()) ) ;
    if( rMeson >fMaxR ) return kFALSE; // cuts on distance from collision point

    Double_t rapidity = 10.;
    if(fMCMother->E() - fMCMother->Pz() == 0 || fMCMother->E() + fMCMother->Pz() == 0){
      rapidity=8.-fRapidityShift;
    }
    else{
      rapidity = 0.5*(TMath::Log((fMCMother->E()+fMCMother->Pz()) / (fMCMother->E()-fMCMother->Pz())))-fRapidityShift;
    }

    // Rapidity Cut
    if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

    // min Pt Cut
    if(fDoMinPtCut && (fMCMother->Pt() < fMinPt)) return kFALSE;

    // Select only -> ChiC radiative (JPsi+gamma) decay channel
    if(fMCMother->GetNDaughters()!=2)return kFALSE;

    AliAODMCParticle *jpsi   = 0x0;
    AliAODMCParticle *gamma   = 0x0;
    AliAODMCParticle *positron = 0x0;
    AliAODMCParticle *electron = 0x0;

    //Int_t labeljpsiChiC = -1;

    for(Int_t index= fMCMother->GetDaughterFirst();index<= fMCMother->GetDaughterLast();index++){
      if(index < 0) continue;
      AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(AODMCArray->At(index));

      switch( temp->GetPdgCode() ) {
      case 443:
        jpsi =  temp;
        //labeljpsiChiC = index;
        break;
      case 22:
        gamma    =  temp;
        labelgammaChiC = index;
        break;
      }
    }

    if ( !jpsi || ! gamma) return kFALSE;
    if(jpsi->GetNDaughters()!=2)return kFALSE;


    for(Int_t index= jpsi->GetDaughterFirst();index<= jpsi->GetDaughterLast();index++){
      if(index < 0) continue;
      AliAODMCParticle* temp = static_cast<AliAODMCParticle*>(AODMCArray->At(index));
      switch( temp->GetPdgCode() ) {
      case -11:
        electron =  temp;
        labelelectronChiC = index;
        break;
      case 11:
        positron =  temp;
        labelpositronChiC = index;
        break;
      }
    }
    if( !electron || !positron) return kFALSE;
    if( positron && electron && gamma) return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelectedMCChiCAODESD(AliDalitzAODESDMC* fMCMother,AliDalitzEventMC *mcEvent,Int_t & labelelectronChiC, Int_t & labelpositronChiC, Int_t & labelgammaChiC, Double_t fRapidityShift) const{
// Returns true for all ChiC within acceptance cuts for decay into JPsi + gamma -> e+ + e- + gamma
// If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts
    if(!mcEvent)return kFALSE;
    // if(fMCMother->GetPdgCode()==20443 ){
    //    return kFALSE;
    // }
    if(fMCMother->GetPdgCodeG()==10441 || fMCMother->GetPdgCodeG()==10443 || fMCMother->GetPdgCodeG()==445 ){
        if(fMCMother->GetRatioVxyG()>fMaxR)  return kFALSE; // cuts on distance from collision point

        Double_t rapidity = 10.;
        if(fMCMother->EnergyG() - fMCMother->PzG() == 0 || fMCMother->EnergyG() + fMCMother->PzG() == 0){
            rapidity=8.-fRapidityShift;
        }
        else{
            rapidity = 0.5*(TMath::Log((fMCMother->EnergyG()+fMCMother->PzG()) / (fMCMother->EnergyG()-fMCMother->PzG())))-fRapidityShift;
        }

        // Rapidity Cut
        if( rapidity<fRapidityCutMesonMin || rapidity>fRapidityCutMesonMax)return kFALSE;

        // min Pt Cut
        if(fDoMinPtCut && (fMCMother->PtG() < fMinPt)) return kFALSE;

        // Select only -> ChiC radiative (JPsi+gamma) decay channel
        if(fMCMother->GetNDaughtersG()!=2)return kFALSE;

            AliDalitzAODESDMC* jpsi=0x0;
            AliDalitzAODESDMC* gamma=0x0;
            AliDalitzAODESDMC* positron=0x0;
            AliDalitzAODESDMC* electron=0x0;

        //Int_t labeljpsiChiC = -1;

        for(Int_t index= fMCMother->GetFirstDaughterG();index<= fMCMother->GetLastDaughterG();index++){
            if(index < 0) continue;
            std::unique_ptr<AliDalitzAODESDMC> temp = std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(index));
            //TParticle* temp = (TParticle*)mcEvent->Particle( index );

            switch( temp->GetPdgCodeG() ) {
            case 443:
                jpsi =  temp.get();
                //labeljpsiChiC = index;
                break;
            case 22:
                gamma    =  temp.get();
            labelgammaChiC = index;
            break;
            }
        }
    if ( !jpsi || ! gamma ) return kFALSE;
    if(jpsi->GetNDaughtersG()!=2)return kFALSE;

    for(Int_t index= jpsi->GetFirstDaughterG();index<= jpsi->GetLastDaughterG();index++){
        if(index < 0) continue;
        std::unique_ptr<AliDalitzAODESDMC> temp = std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(index));
        //TParticle* temp = (TParticle*)mcEvent->Particle( index );
        switch( temp->GetPdgCodeG() ) {
        case -11:
            electron = temp.get();
            labelelectronChiC = index;
        break;
        case 11:
            positron =  temp.get();
            labelpositronChiC = index;
        break;
        }
    }
    if( !electron || !positron) return kFALSE;
    if( positron && electron && gamma) return kTRUE;
    }
    return kFALSE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal, Double_t fRapidityShift, Int_t leadingCellID1, Int_t leadingCellID2, Char_t recoMeth1, Char_t recoMeth2)
{

  // Selection of reconstructed Meson candidates
  // Use flag IsSignal in order to fill Fill different
  // histograms for Signal and Background
  TH2 *hist=0x0;

  if(IsSignal){hist=fHistoMesonCuts;}
  else{hist=fHistoMesonBGCuts;}

  Int_t cutIndex=0;

  if(hist)hist->Fill(cutIndex, pi0->Pt());
  cutIndex++;

  // Undefined Rapidity -> Floating Point exception (also catch E==pz case)
  if(pi0->E()==pi0->Pz() || (pi0->E()+pi0->Pz())/(pi0->E()-pi0->Pz())<=0){
    if(hist)hist->Fill(cutIndex, pi0->Pt());
    cutIndex++;
    if (!IsSignal)cout << "undefined rapidity" << endl;
    return kFALSE;
  }
  else{
    // PseudoRapidity Cut --> But we cut on Rapidity !!!
    cutIndex++;
    if( (pi0->Rapidity()-fRapidityShift)<fRapidityCutMesonMin || (pi0->Rapidity()-fRapidityShift)>fRapidityCutMesonMax){

      if(hist)hist->Fill(cutIndex, pi0->Pt());
      return kFALSE;
    }
  }
  cutIndex++;

  if (fHistoInvMassBefore) fHistoInvMassBefore->Fill(pi0->M());
  // Mass cut
  if (fIsMergedClusterCut == 1 ){
    if (fEnableMassCut){
      Double_t massMin = FunctionMinMassCut(pi0->E());
      Double_t massMax = FunctionMaxMassCut(pi0->E());
  //     cout << "Min mass: " << massMin << "\t max Mass: " << massMax << "\t mass current: " <<  pi0->M()<< "\t E current: " << pi0->E() << endl;
      if (pi0->M() > massMax || pi0->M() < massMin ){
        if(hist)hist->Fill(cutIndex, pi0->Pt());
        return kFALSE;
      }
    }
    cutIndex++;
  }else if(fIsMergedClusterCut == 2){
    if(fEnableOneCellDistCut && ((leadingCellID1 == leadingCellID2) || fCaloPhotonCuts->AreNeighbours(leadingCellID1,leadingCellID2)) ){
      if(hist)hist->Fill(cutIndex, pi0->Pt());
      return kFALSE;
    }
    cutIndex++;
  }

  // Opening Angle Cut
  //fOpeningAngle=2*TMath::ATan(0.134/pi0->P());// physical minimum opening angle

  if (fAllowCombOnlyInSameRecMethod) {
    if (recoMeth1 != recoMeth2){
      if(hist)hist->Fill(cutIndex, pi0->Pt());
      return kFALSE;
    }
  }

  if( fEnableMinOpeningAngleCut && pi0->GetOpeningAngle() < fOpeningAngle){
    if(hist)hist->Fill(cutIndex, pi0->Pt());
    return kFALSE;
  }

  // Min Opening Angle
  if (fMinOpanPtDepCut == kTRUE) fMinOpanCutMeson = fFMinOpanCut->Eval(pi0->Pt());

  if (pi0->GetOpeningAngle() < fMinOpanCutMeson){
    if(hist)hist->Fill(cutIndex, pi0->Pt());
    return kFALSE;
  }

  // Max Opening Angle
  if (fMaxOpanPtDepCut == kTRUE) fMaxOpanCutMeson = fFMaxOpanCut->Eval(pi0->Pt());

  if( pi0->GetOpeningAngle() > fMaxOpanCutMeson){
    if(hist)hist->Fill(cutIndex, pi0->Pt());
    return kFALSE;
  }

  cutIndex++;

  // Alpha Max Cut
  if (fIsMergedClusterCut == 1 && fAlphaPtDepCut) fAlphaCutMeson = fFAlphaCut->Eval(pi0->E());
  else if (fAlphaPtDepCut == kTRUE) fAlphaCutMeson = fFAlphaCut->Eval(pi0->Pt());

  if(TMath::Abs(pi0->GetAlpha())>fAlphaCutMeson){
    if(hist)hist->Fill(cutIndex, pi0->Pt());
    return kFALSE;
  }
  cutIndex++;

  // Alpha Min Cut
  if(TMath::Abs(pi0->GetAlpha())<fAlphaMinCutMeson){
    if(hist)hist->Fill(cutIndex, pi0->Pt());
    return kFALSE;
  }
  cutIndex++;

  if (fHistoInvMassAfter) fHistoInvMassAfter->Fill(pi0->M());

  if (fIsMergedClusterCut == 0){
    if (fHistoDCAGGMesonBefore)fHistoDCAGGMesonBefore->Fill(pi0->GetDCABetweenPhotons());
    if (fHistoDCARMesonPrimVtxBefore)fHistoDCARMesonPrimVtxBefore->Fill(pi0->GetDCARMotherPrimVtx());

    if (fDCAGammaGammaCutOn){
      if (pi0->GetDCABetweenPhotons() > fDCAGammaGammaCut){
        if(hist)hist->Fill(cutIndex, pi0->Pt());
        return kFALSE;
      }
    }
    cutIndex++;

    if (fDCARMesonPrimVtxCutOn){
      if (pi0->GetDCARMotherPrimVtx() > fDCARMesonPrimVtxCut){
        if(hist)hist->Fill(cutIndex, pi0->Pt());
        return kFALSE;
      }
    }
    cutIndex++;

    if (fHistoDCAZMesonPrimVtxBefore)fHistoDCAZMesonPrimVtxBefore->Fill(pi0->GetDCAZMotherPrimVtx());

    if (fDCAZMesonPrimVtxCutOn){
      if (TMath::Abs(pi0->GetDCAZMotherPrimVtx()) > fDCAZMesonPrimVtxCut){
        if(hist)hist->Fill(cutIndex, pi0->Pt());
        return kFALSE;
      }
    }
    cutIndex++;

    if (fHistoDCAGGMesonAfter)fHistoDCAGGMesonAfter->Fill(pi0->GetDCABetweenPhotons());
    if (fHistoDCARMesonPrimVtxAfter)fHistoDCARMesonPrimVtxAfter->Fill(pi0->GetDCARMotherPrimVtx());
    if (fHistoDCAZMesonPrimVtxAfter)fHistoDCAZMesonPrimVtxAfter->Fill(pi0->M(),pi0->GetDCAZMotherPrimVtx());
  }

  //PtCut
  if(fDoMinPtCut){
      if(pi0->Pt()< fMinPt){
          if(hist)hist->Fill(cutIndex, pi0->Pt());
          return kFALSE;
      }
  }
  cutIndex++;

  //Meson Quality Selection
  if(fDoMesonQualitySelection){
    if( pi0->GetMesonQuality() < fMesonQualityMin){
         if(hist)hist->Fill(cutIndex, pi0->Pt());
          return kFALSE;
      }
  }
  cutIndex++;


  if(hist)hist->Fill(cutIndex, pi0->Pt());
  return kTRUE;
}



//________________________________________________________________________
//________________________________________________________________________
Bool_t AliConversionMesonCuts::UpdateCutString() {
  ///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
    fCutString->SetString(GetCutNumber());
  } else {
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());

  // Initialize Cuts from a given Cut string
  // AliInfo(Form("Set Meson Cutnumber: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsAlnum()){
    AliError("Cut selection is not alphanumeric");
    return kFALSE;
  }

  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
  const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
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
//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  //cout << "Updating cut  " << fgkCutNames[cutID] << " (" << cutID << ") to " << value << endl;
  switch (cutID) {
  case kMesonKind:
    if( SetMesonKind(value)) {
      fCuts[kMesonKind] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kSelectionCut:
    if (fIsMergedClusterCut == 1){
      if( SetSelectionWindowMergedCut(value)) {
        fCuts[kSelectionCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    } else if(fUsePtDepSelectionWindow){
      if( SetSelectionWindowCutPtDep(value)) {
        fCuts[kSelectionCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    } else {
      if( SetSelectionWindowCut(value)) {
        fCuts[kSelectionCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    }
  case kalphaMesonCut:
    if (fIsMergedClusterCut == 1){
      if( SetAlphaMesonMergedCut(value)) {
        fCuts[kalphaMesonCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    } else {
      if( SetAlphaMesonCut(value)) {
        fCuts[kalphaMesonCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    }
  case kPtCut:
    if( SetMinPtCut(value)) {
      fCuts[kPtCut] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kRapidityMesonCut:
    if( SetRapidityMesonCut(value)) {
      fCuts[kRapidityMesonCut] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kBackgroundScheme:
    if( SetBackgroundScheme(value)) {
      fCuts[kBackgroundScheme] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kDegreesForRotationMethod:
    if( SetNDegreesForRotationMethod(value)) {
      fCuts[kDegreesForRotationMethod] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kNumberOfBGEvents:
    if( SetNumberOfBGEvents(value)) {
      fCuts[kNumberOfBGEvents] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kuseMCPSmearing:
    if( SetMCPSmearing(value)) {
      fCuts[kuseMCPSmearing] = value;
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
  case kDcaGammaGamma:
    if( SetDCAGammaGammaCut(value)) {
      fCuts[kDcaGammaGamma] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kDcaZPrimVtx:
    if( SetDCAZMesonPrimVtxCut(value)) {
      fCuts[kDcaZPrimVtx] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kDcaRPrimVtx:
    if( SetDCARMesonPrimVtxCut(value)) {
      fCuts[kDcaRPrimVtx] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kMinOpanMesonCut:
    if( SetMinOpanMesonCut(value)) {
      fCuts[kMinOpanMesonCut] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kMaxOpanMesonCut:
    if( SetMaxOpanMesonCut(value)) {
      fCuts[kMaxOpanMesonCut] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

  case kNCuts:
    cout << "Error:: Cut id out of range"<< endl;
    return kFALSE;
  }

  cout << "Error:: Cut id " << cutID << " not recognized "<< endl;
  return kFALSE;

}


//________________________________________________________________________
void AliConversionMesonCuts::PrintCuts() {
  // Print out current Cut Selection
  for(Int_t ic = 0; ic < kNCuts; ic++) {
    printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
  }
}

//________________________________________________________________________
void AliConversionMesonCuts::PrintCutsWithValues() {
  // Print out current Cut Selection with values
  printf("\nMeson cutnumber \n");
  for(Int_t ic = 0; ic < kNCuts; ic++) {
    printf("%d",fCuts[ic]);
  }
  printf("\n\n");

  printf("Meson cuts \n");
  printf("\t %3.2f < y < %3.2f \n", fRapidityCutMesonMin, fRapidityCutMesonMax);
  if (fEnableOneCellDistCut)  printf("\t Only valid for GammaCalo: one cell distance cut enabled");
  if (fEnableMinOpeningAngleCut) printf("\t theta_{open} > %3.4f\n", fOpeningAngle);
  if (!fAlphaPtDepCut) printf("\t %3.2f < alpha < %3.2f\n", fAlphaMinCutMeson, fAlphaCutMeson);
  else printf("\t alpha pT-dep cut active\n");
  if (!fIsMergedClusterCut){
    if (fDCAGammaGammaCutOn)printf("\t dca_{gamma,gamma} < %3.2f\n", fDCAGammaGammaCut);
    if (fDCARMesonPrimVtxCutOn)printf("\t dca_{R, prim Vtx} < %3.2f\n", fDCARMesonPrimVtxCut);
    if (fDCAZMesonPrimVtxCutOn)printf("\t dca_{Z, prim Vtx} < %3.2f\n\n", fDCAZMesonPrimVtxCut);
  }
  if (fIsMergedClusterCut == 1 && fEnableMassCut){
    printf("\t Meson selection energy dependent\n\n");
  } else {
      if (fAcceptMesonMass)
        printf("\t Meson selection window for further analysis %3.3f > M_{gamma,gamma} > %3.3f\n\n", fSelectionLow, fSelectionHigh);
      else
        printf("\t Meson rejection window has been applied %3.3f > M_{gamma,gamma} > %3.3f\n\n", fSelectionLow, fSelectionHigh);
  }
  if(fDoMinPtCut) printf("\t pT_{min} > %3.4f\n", fMinPt);
  if (!fMinOpanPtDepCut) printf("\t theta_{open} > %3.4f\n", fMinOpanCutMeson);
  else printf("\t Min theta_{open} pT-dep cut active\n");
  if (!fMaxOpanPtDepCut) printf("\t %3.4f < theta_{open}\n", fMaxOpanCutMeson);
  else printf("\t Max theta_{open} pT-dep cut active\n");
  printf("\t Running mode for cutselection (0 std, 2 PCM-Calo): %d\n", fMode);

  printf("Meson BG settings \n");
  if (!fDoBG){
    printf("\t No BG estimation \n");
  } else {
    if (!fUseRotationMethodInBG  & !fUseTrackMultiplicityForBG & !fBackgroundHandler) printf("\t BG scheme: mixing V0 mult \n");
    if (!fUseRotationMethodInBG  & fUseTrackMultiplicityForBG & !fBackgroundHandler) printf("\t BG scheme: mixing track mult \n");
    if (fUseRotationMethodInBG )printf("\t BG scheme: rotation \n");
    if (fUsePtmaxMethodForBG )printf("\t BG scheme: Ptmax \n");
    if (fDoBGProbability) printf("\t -> use BG probability \n");
    if (fBackgroundHandler) printf("\t -> use new BG handler \n");
    printf("\t depth of pool: %d\n", fNumberOfBGEvents);
    if (fUseRotationMethodInBG )printf("\t degree's for BG rotation: %d\n", fNDegreeRotationPMForBG);
    if (!fUseRotationMethodInBG  & !fUseTrackMultiplicityForBG & fBackgroundHandler) printf("\t BG scheme: event plane angle with V0 mult \n");
    if (fDoGammaSwappForBg && fGammaSwappMethodBg == 0) printf("\t BG scheme: new in event rotation with 90 degree rotation \n");
    if (fDoGammaSwappForBg && fGammaSwappMethodBg == 1) printf("\t BG scheme: new in event rotation with random rotation: %d rotations \n", fNumberOfSwappsForBg);
    if (fDoGammaSwappForBg && fGammaSwappMethodBg == 10) printf("\t BG scheme: new in event rotation with TGPS: %d rotations \n", fNumberOfSwappsForBg);
  }
}


//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetMesonKind(Int_t mesonKind){
  // Set Cut
  switch(mesonKind){
  case 0:
    fMesonKind = 0;
    break;
  case 1:
    fMesonKind = 1;
  case 2:
    fMesonKind = 0;
    fDoJetAnalysis = kTRUE;
    break;
  case 3:
    fMesonKind = 0;
    fDoJetAnalysis = kTRUE;
    fDoJetQA = kTRUE;
    break;
  case 4:
    fMesonKind = 0;
    fDoIsolatedAnalysis = kTRUE;
    break;
  case 5:
    fMesonKind = 0;
    fDoHighPtHadronAnalysis = kTRUE;
    break;
  case 6: // out of jet
    fMesonKind = 0;
    fDoJetAnalysis = kTRUE;
    fDoOutOfJet = 1;
    break;
  case 7: // out of jet, only on opposite side of Jet
    fMesonKind = 0;
    fDoJetAnalysis = kTRUE;
    fDoOutOfJet = 2;
    break;
  case 8: // out of jet, in "donut shape" [R, R + 0.2] around Jet
    fMesonKind = 0;
    fDoJetAnalysis = kTRUE;
    fDoOutOfJet = 3;
    break;
  default:
    cout<<"Warning: Meson kind not defined"<<mesonKind<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetMinPtCut(Int_t PtCut){
  // Set Cut on min pT of meson
  switch(PtCut){
  case 0: // no cut on pT
    fMinPt = 0.;
    fDoMinPtCut = kFALSE;
    break;
  case 1:
    fMinPt = 0.4;
    fDoMinPtCut = kTRUE;
    break;
  case 2:
    fMinPt = 0.7;
    fDoMinPtCut = kTRUE;
    break;
  case 3:
    fMinPt = 0.9;
    fDoMinPtCut = kTRUE;
    break;
  case 4:
    fMinPt = 1.0;
    fDoMinPtCut = kTRUE;
    break;
  case 5:
    fMinPt = 1.2;
    fDoMinPtCut = kTRUE;
    break;
  case 6:
    fMinPt = 1.5;
    fDoMinPtCut = kTRUE;
    break;
  case 7:
    fMinPt = 0.5;
    fDoMinPtCut = kTRUE;
    break;
  case 8: // for triggered omega
    fMinPt = 5.0;
    fDoMinPtCut = kTRUE;
    break;
  case 9: // for triggered omega
    fMinPt = 3.0;
    fDoMinPtCut = kTRUE;
    break;
  case 10: // a for triggered omega
    fMinPt = 4.0;
    fDoMinPtCut = kTRUE;
    break;
  case 11: // b for triggered omega
    fMinPt = 6.0;
    fDoMinPtCut = kTRUE;
    break;
  case 12: // c for triggered omega
    fMinPt = 8.0;
    fDoMinPtCut = kTRUE;
    break;
  case 13: // d for triggered omega
    fMinPt = 10.0;
    fDoMinPtCut = kTRUE;
    break;

  // Instead of applying pt cut, apply a min energy cut on daughters
  // (needs to be treated in analysis task)
  case 14: // e
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 5.;
    break;
  case 15: // f
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 7.5;
    break;
  case 16: // g
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 6.;
    break;
  case 17: // h
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 10.;
    break;
  case 18: // i
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 12.;
    break;
  case 19: //j
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 2.;
    break;
  case 20: //k
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 4.;
    break;
  case 21: //l
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 8.0;
    break;
  case 22: //m
    fDoGammaMinEnergyCut = kTRUE;
    fNDaughterEnergyCut  = 1;
    fSingleDaughterMinE  = 9.0;
      break;
  case 23: //n
    cout<<"Warning: pT cut not defined"<<PtCut<<endl;
    return kFALSE;
    break;
  case 24: //o
    cout<<"Warning: pT cut not defined"<<PtCut<<endl;
    return kFALSE;
    break;
  case 25: //p
    cout<<"Warning: pT cut not defined"<<PtCut<<endl;
    return kFALSE;
    break;
    //continue normal pT cuts on meson
  case 26: //q
      fMinPt = 12.0;
      fDoMinPtCut = kTRUE;
      break;
  case 27: //r
      fMinPt = 14.0;
      fDoMinPtCut = kTRUE;
      break;
  default:
    cout<<"Warning: pT cut not defined"<<PtCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetSelectionWindowCut(Int_t selectionCut){
  // Set Cut
  switch(selectionCut){
    case 0:
      fSelectionLow       = 0.0;
      fSelectionHigh      = 4.;
      fAcceptMesonMass    = kTRUE;
      break;
    case 1:
      fSelectionLow       = 0.1;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kTRUE;
      break;
    case 2:
      fSelectionLow       = 0.11;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kTRUE;
      break;
    case 3:
      fSelectionLow       = 0.12;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kTRUE;
      break;
    case 4:
      fSelectionLow       = 0.1;
      fSelectionHigh      = 0.15;
      fAcceptMesonMass    = kTRUE;
      break;
    case 5:
      fSelectionLow       = 0.11;
      fSelectionHigh      = 0.15;
      fAcceptMesonMass    = kTRUE;
      break;
    case 6:
      fSelectionLow       = 0.12;
      fSelectionHigh      = 0.15;
      fAcceptMesonMass    = kTRUE;
      break;
    case 7:
      fSelectionLow       = 0.1;
      fSelectionHigh      = 0.155;
      fAcceptMesonMass    = kTRUE;
      break;
    case 8:
      fSelectionLow       = 0.125;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kTRUE;
      break;
    case 9:
      fSelectionLow       = 0.11;
      fSelectionHigh      = 0.155;
      fAcceptMesonMass    = kTRUE;
      break;
    case 10: //a
      fSelectionLow       = 0.08;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kTRUE;
      break;
    case 11: //b
      fSelectionLow       = 0.08;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kFALSE;
      break;
    case 12: //c
      fSelectionLow       = 0.1;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kFALSE;
      break;
    case 13: //d
      fSelectionLow       = 0.11;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kFALSE;
      break;
    case 14: //e
      fSelectionLow       = 0.12;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kFALSE;
      break;
    case 15: //f
      fSelectionLow       = 0.1;
      fSelectionHigh      = 0.15;
      fAcceptMesonMass    = kFALSE;
      break;
    case 16: //g
      fSelectionLow       = 0.11;
      fSelectionHigh      = 0.15;
      fAcceptMesonMass    = kFALSE;
      break;
    case 17: //h
      fSelectionLow       = 0.12;
      fSelectionHigh      = 0.15;
      fAcceptMesonMass    = kFALSE;
      break;
    case 18: //i
      fSelectionLow       = 0.1;
      fSelectionHigh      = 0.155;
      fAcceptMesonMass    = kFALSE;
      break;
    case 19: //j
      fSelectionLow       = 0.125;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kFALSE;
      break;
    case 20: //k
      fSelectionLow       = 0.11;
      fSelectionHigh      = 0.155;
      fAcceptMesonMass    = kFALSE;
      break;
    case 21: //l
      fSelectionLow       = 0.500;
      fSelectionHigh      = 0.600;
      fAcceptMesonMass    = kTRUE;
      break;
    case 22: //m
      fSelectionLow       = 0.400;
      fSelectionHigh      = 0.700;
      fAcceptMesonMass    = kTRUE;
      break;
    case 23: //n
      fSelectionLow       = 0.120;
      fSelectionHigh      = 0.160;
      fAcceptMesonMass    = kTRUE;
      break;
    // Pt dependent around pi0 mass
    case 24: //o // EMC-EMC
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 0;
      break;
    case 25: //p // PCM-EMC
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 1;
      break;
    case 26: //q // PHOS-PHOS
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 2;
      break;
    case 27: //r // PCM-PHOS
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 3;
      break;
    case 28: //s // PCM-PCM
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 4;
      break;
    case 29: // t // EMC-EMC  gamma selection
      fSelectionLow       = 0.1;
      fSelectionHigh      = 0.145;
      fAcceptMesonMass    = kTRUE;
      fUseGammaSelection  = kTRUE;

      break;
    default:
      cout<<"Warning: SelectionCut not defined "<<selectionCut<<endl;
      return kFALSE;
  }
  return kTRUE;
}

Bool_t AliConversionMesonCuts::SetSelectionWindowMergedCut(Int_t selectionCut){
  // Set Cut
  fSelectionWindowCut = selectionCut;
  switch(fSelectionWindowCut){
    case 0:
      fEnableMassCut = kFALSE;
      break;
    case 1:   //NLM 1
      fEnableMassCut = kTRUE;
      break;
    case 2:   //NLM 2
      fEnableMassCut = kTRUE;
      break;
    case 3:   //NLM 1
      fEnableMassCut = kTRUE;
      break;
    case 4:   //NLM 2
      fEnableMassCut = kTRUE;
      break;
    case 5:   //NLM 1
      fEnableMassCut = kTRUE;
      break;
    case 6:   //NLM 2
      fEnableMassCut = kTRUE;
      break;
    case 7:   //min mass cut around 0
      fEnableMassCut = kTRUE;
      break;
    case 8:   //min mass cut around 0
      fEnableMassCut = kTRUE;
      break;
    case 9:   //min mass cut around 0
      fEnableMassCut = kTRUE;
      break;
    default:
      cout<<"Warning: SelectionCut merged not defined "<<selectionCut<<endl;
      return kFALSE;
  }

  return kTRUE;

}

Bool_t AliConversionMesonCuts::SetSelectionWindowCutPtDep(Int_t selectionCut){
  // Set Cut
  fSelectionWindowCut = selectionCut;
  switch(fSelectionWindowCut){
    case 0: // EMC-EMC - 99 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 99.;
      fSelectionNSigmaHigh = 99.;
      fMassParamFunction   = 0;
      break;
    case 1: // EMC-EMC - 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 0;
      break;
    case 2: // PCM-EMC 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 1;
      break;
    case 3: // PHOS-PHOS 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 2;
      break;
    case 4: // PCM-PHOS 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 3;
    case 5: //PCM-PCM 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 4;
      break;
    // Variations for 5 TeV
    case 6: // EMC-EMC - 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 5;
      break;
    case 7: // PCM-EMC 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 6;
      break;
    case 8: // PHOS-PHOS 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 7;
      break;
    case 9: // PCM-PHOS 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 8;
    case 10: // a PCM-PCM 2 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 9;
      break;
    case 11: // b EMC-EMC - 1 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 5;
      break;
    case 12: // c EMC-EMC - 3 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 5;
      break;
    case 13: // d EMC-EMC - 4 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 4.;
      fSelectionNSigmaHigh = 4.;
      fMassParamFunction   = 5;
      break;
    case 14: // e PCM-EMC 1 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 6;
      break;
    case 15: // f PCM-EMC 3 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 6;
      break;
    case 16: // g PCM-EMC 4 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 4.;
      fSelectionNSigmaHigh = 4.;
      fMassParamFunction   = 6;
      break;
    case 17: // h PHOS-PHOS 1 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 7;
      break;
    case 18: // i PHOS-PHOS 3 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 7;
      break;
    case 19: // j PHOS-PHOS 4 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 4.;
      fSelectionNSigmaHigh = 4.;
      fMassParamFunction   = 7;
      break;
    case 20: // k PCM-PHOS 1 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 8;
      break;
    case 21: // l PCM-PHOS 3 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 8;
      break;
    case 22: // m PCM-PHOS 4 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 4.;
      fSelectionNSigmaHigh = 4.;
      fMassParamFunction   = 8;
      break;
    case 23: // n// PCM-PCM 1 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 9;
      break;
    case 24: // o // PCM-PCM 3 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 9;
      break;
    case 25: // p // PCM-PCM 4 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 4.;
      fSelectionNSigmaHigh = 4.;
      fMassParamFunction   = 9;
      break;
    case 26: // q // PHOS-PHOS 1 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 2;
      break;
    case 27: // r // PCM-PHOS 1 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 3;
      break;
    case 28: // s //PHOS-PHOS 3 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 2;
      break;
    case 29: // t // PCM-PHOS 3 sigma
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 3;
      break;
      // Varation for 13 TeV EDC
    case 30: // u // EMC-EMC - 1 sigma - gamma selection
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 1.;
      fSelectionNSigmaHigh = 1.;
      fMassParamFunction   = 0;
      break;
    case 31: // v // EMC-EMC - 3 sigma - gamma selection
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 3.;
      fSelectionNSigmaHigh = 3.;
      fMassParamFunction   = 0;
      break;
    case 32: // w // EMC-EMC - 4 sigma - gamma selection
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 4.;
      fSelectionNSigmaHigh = 4.;
      fMassParamFunction   = 0;
      break;
    case 33: // x // EMC-EMC - 2 sigma - gamma selection
      fAcceptMesonMass     = kFALSE;
      fUsePtDepSelectionWindow = kTRUE;
      fUseGammaSelection   = kTRUE;
      fSelectionNSigmaLow  = 2.;
      fSelectionNSigmaHigh = 2.;
      fMassParamFunction   = 0;
      break;
    default:
      cout<<"Warning: SelectionCut merged not defined "<<selectionCut<<endl;
      return kFALSE;
  }

  return kTRUE;

}

Float_t AliConversionMesonCuts::FunctionMaxMassCut(Float_t e){

  Float_t switchMass    = 0;
  Float_t aMassLow      = 0;
  Float_t bMassLow      = 0;
  Float_t aMassHigh     = 0;
  Float_t bMassHigh     = 0;
  Float_t aMass         = 0;
  Float_t bMass         = 0;
  Float_t switchSigma   = 0.;
  Float_t nSigma        = 0;
  Float_t aSigmaLow     = 0.;
  Float_t bSigmaLow     = 0;
  Float_t aSigmaHigh    = 0.;
  Float_t bSigmaHigh    = 0;
  Float_t mass          = 0;
  Float_t sigma         = 0;

  switch(fSelectionWindowCut){
    case 0:
      fEnableMassCut = kFALSE;
      break;
    case 1:   //NLM 1
      aMass         = 0.044;
      bMass         = 0.0049;
      switchSigma   = 19.;
      nSigma        = 3;
      aSigmaLow     = 0.012;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0012;
      bSigmaHigh    = 6e-4;

      mass          = aMass + bMass*e;
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
      return mass + nSigma*sigma;
      break;
    case 2:  //NLM 2
      switchMass    = 21;
      aMassLow      = 0.115;
      bMassLow      = 9.6e-4;
      aMassHigh     = 0.1;
      bMassHigh     = 0.0017;
      switchSigma   = 10.;
      nSigma        = 3;
      aSigmaLow     = 0.009;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0023;
      bSigmaHigh    = 6.7e-4;

      mass          = 0;
      if (e < switchMass){
        mass        = aMassLow + bMassLow*e;
      } else {
        mass        = aMassHigh + bMassHigh*e;
      }
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
      return mass + nSigma*sigma;
      break;
    case 3:   //NLM 1
      aMass         = 0.044;
      bMass         = 0.0049;
      switchSigma   = 19.;
      nSigma        = 2;
      aSigmaLow     = 0.012;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0012;
      bSigmaHigh    = 6e-4;

      mass          = aMass + bMass*e;
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
      return mass + nSigma*sigma;
      break;
    case 4:  //NLM 2
      switchMass    = 21;
      aMassLow      = 0.115;
      bMassLow      = 9.6e-4;
      aMassHigh     = 0.1;
      bMassHigh     = 0.0017;
      switchSigma   = 10.;
      nSigma        = 2;
      aSigmaLow     = 0.009;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0023;
      bSigmaHigh    = 6.7e-4;

      mass          = 0;
      if (e < switchMass){
        mass        = aMassLow + bMassLow*e;
      } else {
        mass        = aMassHigh + bMassHigh*e;
      }
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
      return mass + nSigma*sigma;
      break;
    case 5:   //NLM 1
      aMass         = 0.044;
      bMass         = 0.0049;
      switchSigma   = 19.;
      nSigma        = 4;
      aSigmaLow     = 0.012;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0012;
      bSigmaHigh    = 6e-4;

      mass          = aMass + bMass*e;
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
      return mass + nSigma*sigma;
      break;
    case 6:  //NLM 2
      switchMass    = 21;
      aMassLow      = 0.115;
      bMassLow      = 9.6e-4;
      aMassHigh     = 0.1;
      bMassHigh     = 0.0017;
      switchSigma   = 10.;
      nSigma        = 4;
      aSigmaLow     = 0.009;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0023;
      bSigmaHigh    = 6.7e-4;

      mass          = 0;
      if (e < switchMass){
        mass        = aMassLow + bMassLow*e;
      } else {
        mass        = aMassHigh + bMassHigh*e;
      }
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
      return mass + nSigma*sigma;
      break;
    case 7: // maximum mass
      return 10000.;
      break;
    case 8: // maximum mass
      return 10000.;
      break;
    case 9: // maximum mass
      return 10000.;
      break;
    default:
      cout<<"Warning: SelectionCut merged not defined "<<fSelectionWindowCut<<endl;
      return -1;
  }
  return -1;

}

Float_t AliConversionMesonCuts::FunctionMinMassCut(Float_t e){

  Float_t switchMass      = 0;
  Float_t aMassLow        = 0;
  Float_t bMassLow        = 0;
  Float_t aMassHigh       = 0;
  Float_t bMassHigh       = 0;
  Float_t aMass           = 0;
  Float_t bMass           = 0;
  Float_t switchSigma     = 0.;
  Float_t nSigma          = 0;
  Float_t aSigmaLow       = 0.;
  Float_t bSigmaLow       = 0;
  Float_t aSigmaHigh      = 0.;
  Float_t bSigmaHigh      = 0;
  Float_t mass            = 0;
  Float_t sigma           = 0;

  switch(fSelectionWindowCut){
    case 0:
      fEnableMassCut      = kFALSE;
      break;
    case 1:   //NLM 1
      aMass         = 0.044;
      bMass         = 0.0049;
      switchSigma   = 19.;
      nSigma        = 3;
      aSigmaLow     = 0.012;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0012;
      bSigmaHigh    = 6e-4;

      mass          = aMass + bMass*e;
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
  //     cout << "E: "<< e << "\t mass: " << mass << "\t sigma: " << sigma<<  endl;
      return mass - nSigma*sigma;
      break;
    case 2:  //NLM 2
      switchMass    = 21;
      aMassLow      = 0.115;
      bMassLow      = 9.6e-4;
      aMassHigh     = 0.1;
      bMassHigh     = 0.0017;
      switchSigma   = 10.;
      nSigma        = 3;
      aSigmaLow     = 0.009;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0023;
      bSigmaHigh    = 6.7e-4;

      mass          = 0;
      if (e < switchMass){
        mass        = aMassLow + bMassLow*e;
      } else {
        mass        = aMassHigh + bMassHigh*e;
      }
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
  //     cout << "E: "<< e << "\t mass: " << mass << "\t sigma: "<< sigma << endl;

      return mass - nSigma*sigma;
      break;
    case 3:   //NLM 1
      aMass         = 0.044;
      bMass         = 0.0049;
      switchSigma   = 19.;
      nSigma        = 2;
      aSigmaLow     = 0.012;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0012;
      bSigmaHigh    = 6e-4;

      mass          = aMass + bMass*e;
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
  //     cout << "E: "<< e << "\t mass: " << mass << "\t sigma: " << sigma<<  endl;
      return mass - nSigma*sigma;
      break;
    case 4:  //NLM 2
      switchMass    = 21;
      aMassLow      = 0.115;
      bMassLow      = 9.6e-4;
      aMassHigh     = 0.1;
      bMassHigh     = 0.0017;
      switchSigma   = 10.;
      nSigma        = 2;
      aSigmaLow     = 0.009;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0023;
      bSigmaHigh    = 6.7e-4;

      mass          = 0;
      if (e < switchMass){
        mass        = aMassLow + bMassLow*e;
      } else {
        mass        = aMassHigh + bMassHigh*e;
      }
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
  //     cout << "E: "<< e << "\t mass: " << mass << "\t sigma: "<< sigma << endl;

      return mass - nSigma*sigma;
      break;
    case 5:   //NLM 1
      aMass         = 0.044;
      bMass         = 0.0049;
      switchSigma   = 19.;
      nSigma        = 4;
      aSigmaLow     = 0.012;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0012;
      bSigmaHigh    = 6e-4;

      mass          = aMass + bMass*e;
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
  //     cout << "E: "<< e << "\t mass: " << mass << "\t sigma: " << sigma<<  endl;
      return mass - nSigma*sigma;
      break;
    case 6:  //NLM 2
      switchMass    = 21;
      aMassLow      = 0.115;
      bMassLow      = 9.6e-4;
      aMassHigh     = 0.1;
      bMassHigh     = 0.0017;
      switchSigma   = 10.;
      nSigma        = 4;
      aSigmaLow     = 0.009;
      bSigmaLow     = 0;
      aSigmaHigh    = 0.0023;
      bSigmaHigh    = 6.7e-4;

      mass          = 0;
      if (e < switchMass){
        mass        = aMassLow + bMassLow*e;
      } else {
        mass        = aMassHigh + bMassHigh*e;
      }
      sigma         = 0;
      if (e < switchSigma){
        sigma       = aSigmaLow + bSigmaLow*e;
      } else {
        sigma       = aSigmaHigh + bSigmaHigh*e;
      }
  //     cout << "E: "<< e << "\t mass: " << mass << "\t sigma: "<< sigma << endl;
      return mass - nSigma*sigma;
      break;

    case 7: // just exclude band at 0
      return 0.005+0.004*e;
      break;
    case 8: // just exclude band at 0 looser
      return 0.004+0.004*e;
      break;
    case 9: // just exclude band at 0 tighter
      return 0.006+0.004*e;
      break;

    default:
      cout<<"Warning: SelectionCut merged not defined "<<fSelectionWindowCut<<endl;
      return -1;
  }
  return -1;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetAlphaMesonCut(Int_t alphaMesonCut)
{ // Set Cut
  switch(alphaMesonCut){
  case 0:  // 0- 0.7
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.7;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 1:  // Updated 15 May 2015
    if (fIsMergedClusterCut == 0){
      if( fFAlphaCut ) delete fFAlphaCut;
      fFAlphaCut        = new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
      fFAlphaCut->SetParameter(0,0.65);
      fFAlphaCut->SetParameter(1,1.8);
      fAlphaMinCutMeson =  0.0;
      fAlphaCutMeson    = -1.0;
      fAlphaPtDepCut    = kTRUE;
    } else {
      fAlphaPtDepCut    = kFALSE;
      fAlphaMinCutMeson = 0.5;
      fAlphaCutMeson    = 1;
    }
    break;
  case 2:  // Updated 31 October 2013 before 0.5-1
    if (fIsMergedClusterCut == 0){
      if( fFAlphaCut ) delete fFAlphaCut;
      fFAlphaCut        = new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
      fFAlphaCut->SetParameter(0,0.8);
      fFAlphaCut->SetParameter(1,1.2);
      fAlphaMinCutMeson =  0.0;
      fAlphaCutMeson    = -1.0;
      fAlphaPtDepCut    = kTRUE;
    } else {
      fAlphaPtDepCut    = kFALSE;
      fAlphaMinCutMeson = 0.6;
      fAlphaCutMeson    = 1;
    }
    break;
  case 3:  // 0.0-1
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 1.;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 4:  // 0-0.65
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.65;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 5:  // 0-0.75
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.75;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 6:  // 0-0.8
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.8;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 7:  // 0.0-0.85
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.85;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 8:  // 0.0-0.6
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.6;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 9: // Updated 11 November 2013 before 0.0 - 0.3
    if (fIsMergedClusterCut == 0){
      if( fFAlphaCut ) delete fFAlphaCut;
      fFAlphaCut        = new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
      fFAlphaCut->SetParameter(0,0.65);
      fFAlphaCut->SetParameter(1,1.2);
      fAlphaMinCutMeson =  0.0;
      fAlphaCutMeson    = -1.0;
      fAlphaPtDepCut    = kTRUE;
    } else {
      fAlphaPtDepCut    = kFALSE;
      fAlphaMinCutMeson = 0.4;
      fAlphaCutMeson    = 1;
    }
    break;
  case 10:  //a 0-0.2
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.2;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 11:  //b 0.2-0.6
    fAlphaMinCutMeson   = 0.2;
    fAlphaCutMeson      = 0.6;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 12:  //c 0.6-1.0
    fAlphaMinCutMeson   = 0.6;
    fAlphaCutMeson      = 1.0;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 13:  //d 0-0.1
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.1;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 14:  //e 0-0.3
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.3;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 15:  //f 0-0.4
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.4;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 16:  //g 0-0.5
    fAlphaMinCutMeson   = 0.0;
    fAlphaCutMeson      = 0.5;
    fAlphaPtDepCut      = kFALSE;
    break;
  case 17:  // h (for lowB)
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.8);
    fFAlphaCut->SetParameter(1,2.9);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;
  case 18:  // i (for lowB)
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]*tanh([1]*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.75);
    fFAlphaCut->SetParameter(1,3.5);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;


  default:
    cout<<"Warning: AlphaMesonCut not defined "<<alphaMesonCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetAlphaMesonMergedCut(Int_t alphaMesonCut)
{ // Set Cut
  switch(alphaMesonCut){
  case 0:  // 0- 1
    fAlphaMinCutMeson = 0.0;
    fAlphaCutMeson    = 1;
    fAlphaPtDepCut    = kFALSE;
    break;
  case 1:  // cut for NLM 1
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]+[1]*x+[2]/(x*x*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.96);
    fFAlphaCut->SetParameter(1,0);
    fFAlphaCut->SetParameter(2,-879);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;
  case 2:  // cut for NLM 2
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]+[1]*x+[2]/(x*x*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.95);
    fFAlphaCut->SetParameter(1,0.0015);
    fFAlphaCut->SetParameter(2,-233);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;
  case 3:  // cut for NLM 1 larger
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]+[1]*x+[2]/(x*x*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.975);
    fFAlphaCut->SetParameter(1,0);
    fFAlphaCut->SetParameter(2,-800);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;
  case 4:  // cut for NLM 2 larger
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]+[1]*x+[2]/(x*x*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.97);
    fFAlphaCut->SetParameter(1,0.0015);
    fFAlphaCut->SetParameter(2,-200);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;
  case 5:  // cut for NLM 1 smaller
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]+[1]*x+[2]/(x*x*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.94);
    fFAlphaCut->SetParameter(1,0);
    fFAlphaCut->SetParameter(2,-970);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;
  case 6:  // cut for NLM 2 smaller
    if( fFAlphaCut ) delete fFAlphaCut;
    fFAlphaCut        = new TF1("fFAlphaCut","[0]+[1]*x+[2]/(x*x*x)",0.,100.);
    fFAlphaCut->SetParameter(0,0.935);
    fFAlphaCut->SetParameter(1,0.0015);
    fFAlphaCut->SetParameter(2,-273);
    fAlphaMinCutMeson =  0.0;
    fAlphaCutMeson    = -1.0;
    fAlphaPtDepCut    = kTRUE;
    break;

  default:
    cout<<"Warning: AlphaMesonCut for merged clusters not defined "<<alphaMesonCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetRapidityMesonCut(Int_t RapidityMesonCut){
  // Set Cut
  switch(RapidityMesonCut){
  case 0:  // changed from 0.9 to 1.35
    fRapidityCutMesonMax   = 1.35;
    fRapidityCutMesonMin   = -1.35;
    break;
  case 1:  //
    fRapidityCutMesonMax   = 0.8;
    fRapidityCutMesonMin   = -0.8;
    break;
  case 2:  //
    fRapidityCutMesonMax   = 0.7;
    fRapidityCutMesonMin   = -0.7;
    break;
  case 3:  //
    fRapidityCutMesonMax   = 0.6;
    fRapidityCutMesonMin   = -0.6;
    break;
  case 4:  //
    fRapidityCutMesonMax   = 0.5;
    fRapidityCutMesonMin   = -0.5;
    break;
  case 5:  //
    fRapidityCutMesonMax   = 0.85;
    fRapidityCutMesonMin   = -0.85;
    break;
  case 6:  //
    fRapidityCutMesonMax   = 0.75;
    fRapidityCutMesonMin   = -0.75;
    break;
  case 7:  //
    fRapidityCutMesonMax   = 0.3;
    fRapidityCutMesonMin   = -0.3;
    break;
  case 8:  //changed, before 0.35
    fRapidityCutMesonMax   = 0.25;
    fRapidityCutMesonMin   = -0.25;
    break;
  case 9:  //
    fRapidityCutMesonMax   = 0.4;
    fRapidityCutMesonMin   = -0.4;
    break;
  case 10:  // a
    fRapidityCutMesonMax   = 0.8;
    fRapidityCutMesonMin   = 0.0;
    break;
  case 11:  // b
    fRapidityCutMesonMax   = 0.0;
    fRapidityCutMesonMin   = -0.8;
    break;
  case 12:  // c
    fRapidityCutMesonMax   = 0.6;
    fRapidityCutMesonMin   = 0.0;
    break;
  case 13:  // d
    fRapidityCutMesonMax   = 0.0;
    fRapidityCutMesonMin   = -0.6;
    break;
  case 14:  // e
    fRapidityCutMesonMax   = 1.0;
    fRapidityCutMesonMin   = 0.0;
    break;
  case 15:  // f
    fRapidityCutMesonMax   = 0.0;
    fRapidityCutMesonMin   = -1.0;
    break;
  default:
    cout<<"Warning: RapidityMesonCut not defined "<<RapidityMesonCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetBackgroundScheme(Int_t BackgroundScheme){
  // Set Cut
  switch(BackgroundScheme){
  case 0: //Rotation
    fUseRotationMethodInBG      = kTRUE;
    fDoBGProbability            = kFALSE;
    break;
  case 1: // mixed event with V0 multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;

    // only used for heavy meson task
    // no pions from same event
    fBackgroundMode             = 4;
    break;
  case 2: // mixed event with track multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kTRUE;
    fDoBGProbability            = kFALSE;
    break;
  case 3: //Rotation
    fUseRotationMethodInBG      = kTRUE;
    fDoBGProbability            = kTRUE;
    break;
  case 4: //No BG calculation
    cout << "no BG calculation should be done" << endl;
    fUseRotationMethodInBG      = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoBG                       = kFALSE;
    break;
  case 5: //Rotation
    fUseRotationMethodInBG      = kTRUE;
    fDoBGProbability            = kFALSE;
    fBackgroundHandler          = 1;
    break;
  case 6: // mixed event with V0 multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundHandler          = 1;
    break;
  case 7: // mixed event with track multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kTRUE;
    fDoBGProbability            = kFALSE;
    fBackgroundHandler          = 1;
    break;
  case 8: //Rotation
    fUseRotationMethodInBG      = kTRUE;
    fDoBGProbability            = kTRUE;
    fBackgroundHandler          = 1;
    break;
  case 9: // mixed event with Ptmax method
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fUsePtmaxMethodForBG        = kTRUE;
    break;
  case 10: // a mixed event with likesign mixing
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundUseLikeSign      = kTRUE;
    fBackgroundUseSideband      = kFALSE;
    break;
  case 11: // b same event pi0 sideband candidates (right side of pi0 peak)
    fBackgroundMode             = 6;
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundUseLikeSign      = kFALSE;
    fBackgroundUseSideband      = kTRUE;
    fSidebandMixingLow          = 0.180;
    fSidebandMixingHigh         = 0.220;
    break;
  case 12: // c same event with pi0 sideband candidates (left side of pi0 peak)
    fBackgroundMode             = 6;
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundUseLikeSign      = kFALSE;
    fBackgroundUseSideband      = kTRUE;
    fSidebandMixingLow          = 0.01;
    fSidebandMixingHigh         = 0.05;
    break;
  case 13: // d same event with pi0 sideband candidates (both sides of pi0 peak)
    fBackgroundMode                  = 6;
    fUseRotationMethodInBG           = kFALSE;
    fUseTrackMultiplicityForBG       = kFALSE;
    fDoBGProbability                 = kFALSE;
    fBackgroundUseLikeSign           = kFALSE;
    fBackgroundUseSideband           = kFALSE;
    fBackgroundUseSidebandBothSides  = kTRUE;
    fSidebandMixingLeftLow           = 0.01;
    fSidebandMixingLeftHigh          = 0.05;
    fSidebandMixingRightLow          = 0.180;
    fSidebandMixingRightHigh         = 0.220;
    break;
  case 14: //e same event with pi0 sideband candidates (right side of pi0 peak)
    fBackgroundMode             = 6;
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundUseLikeSign      = kFALSE;
    fBackgroundUseSideband      = kTRUE;
    fSidebandMixingLow          = 0.600;
    fSidebandMixingHigh         = 0.650;
    break;
  case 15: //f same event with pi0 sideband candidates (left side of pi0 peak)
    fBackgroundMode             = 6;
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundUseLikeSign      = kFALSE;
    fBackgroundUseSideband      = kTRUE;
    fSidebandMixingLow          = 0.42;
    fSidebandMixingHigh         = 0.47;
    break;
  case 16: //g same event with pi0 sideband candidates (both sides of pi0 peak)
    fBackgroundMode                  = 6;
    fUseRotationMethodInBG           = kFALSE;
    fUseTrackMultiplicityForBG       = kFALSE;
    fDoBGProbability                 = kFALSE;
    fBackgroundUseLikeSign           = kFALSE;
    fBackgroundUseSideband           = kFALSE;
    fBackgroundUseSidebandBothSides  = kTRUE;
    fSidebandMixingLeftLow           = 0.42;
    fSidebandMixingLeftHigh          = 0.47;
    fSidebandMixingRightLow          = 0.600;
    fSidebandMixingRightHigh         = 0.650;
    break;
  case 17: //h mixed event with V0 multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoConvCaloMixing           = kTRUE;
    break;
  case 18: //i mixed event with V0 multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoSectorMixing             = kTRUE;
    break;
  case 19: //j mixed event with V0 multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoSectorJetMixing          = kTRUE;
    break;
  case 20: //k mixed by jet distance
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoJetMixing                = kTRUE;
    break;
  case 21: //l mixed by jet distance and rotation
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoJetMixing                = kTRUE;
    fDoJetRotateMixing          = kTRUE;
    break;
  case 22: //m mixed by jet distance and jet pt
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoJetMixing                = kTRUE;
    fDoJetPtMixing              = kTRUE;
    break;
  case 23: //n mixed by jet distance, rotation and jet pt
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fDoJetMixing                = kTRUE;
    fDoJetRotateMixing          = kTRUE;
    fDoJetPtMixing              = kTRUE;
    break;
  case 24: // o mixed event for three pions (pi- and pi+ same event)
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundMode             = 1;
    break;
  case 25: // p mixed event for three pions (pi+ and pi0 same event)
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundMode             = 2;
    break;
  case 26: // q mixed event for three pions (pi- and pi0 same event)
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fDoBGProbability            = kFALSE;
    fBackgroundMode             = 3;
    break;
  case 27: // r cluster swapping method with 90 degree rotation angle
    fDoGammaSwappForBg          = 1;
    fDoWeightingInSwappBg       = kFALSE;
    fGammaSwappMethodBg         = 0;
    fNumberOfSwappsForBg        = 1;
    fBackgroundHandler          = 2;
    fBackgroundMode             = 7;
    break;
  case 28: // s cluster swapping method with 90 degree rotation angle with weighting
    fDoGammaSwappForBg          = 1;
    fDoWeightingInSwappBg       = kTRUE;
    fGammaSwappMethodBg         = 0;
    fNumberOfSwappsForBg        = 1;
    fBackgroundHandler          = 2;
    break;
  case 29: // t cluster swapping method with random (between 60 & 120 + 240 & 300) rotation angle
    fDoGammaSwappForBg          = 1;
    fDoWeightingInSwappBg       = kTRUE;
    fGammaSwappMethodBg         = 1;
    fNumberOfSwappsForBg        = 1;
    fBackgroundHandler          = 2;
  case 30: // u cluster swapping method with 4 random (between 60 & 120 + 240 & 300) rotation angle
    fDoGammaSwappForBg          = 1;
    fDoWeightingInSwappBg       = kTRUE;
    fGammaSwappMethodBg         = 1;
    fNumberOfSwappsForBg        = 4;
    fBackgroundHandler          = 2;
    break;
  case 31: // v cluster swapping method with 20 random with TGenPhaseSpace no evt weighting
    fDoGammaSwappForBg          = 1;
    fDoWeightingInSwappBg       = kFALSE;
    fGammaSwappMethodBg         = 10;
    fNumberOfSwappsForBg        = 20;
    fBackgroundHandler          = 2;
    break;
  case 32: // w cluster swapping method with 20 random with TGenPhaseSpace with event weighting
    fDoGammaSwappForBg          = 1;
    fDoWeightingInSwappBg       = kTRUE;
    fGammaSwappMethodBg         = 10;
    fNumberOfSwappsForBg        = 20;
    fBackgroundHandler          = 2;
    break;
  case 33: // x cluster swapping method with 20 random with TGenPhaseSpace with event weighting & forbid decays that are similar to original decay
    fDoGammaSwappForBg          = 1;
    fDoWeightingInSwappBg       = kTRUE;
    fGammaSwappMethodBg         = 11;
    fNumberOfSwappsForBg        = 20;
    fBackgroundHandler          = 2;
    break;
  case 34: // y cluster swapping method with 90 degree rotation angle (around Pi0 for the omega analyses)
    fDoGammaSwappForBg          = 2;
    fDoWeightingInSwappBg       = kFALSE;
    fGammaSwappMethodBg         = 0;
    fNumberOfSwappsForBg        = 1;
    fBackgroundHandler          = 2;
    break;
  case 35: // z cluster swapping method with 20 random with TGenPhaseSpace with event weighting (around Pi0 for the omega analyses)
    fDoGammaSwappForBg          = 2;
    fDoWeightingInSwappBg       = kTRUE;
    fGammaSwappMethodBg         = 10;
    fNumberOfSwappsForBg        = 20;
    fBackgroundHandler          = 2;
    break;
  default:
    cout<<"Warning: BackgroundScheme not defined "<<BackgroundScheme<<endl;
    return kFALSE;
  }
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetNDegreesForRotationMethod(Int_t DegreesForRotationMethod){
  // Set Cut
  switch(DegreesForRotationMethod){
  case 0:
    fNDegreeRotationPMForBG = 5;
    break;
  case 1:
    fNDegreeRotationPMForBG = 10;
    break;
  case 2:
    fNDegreeRotationPMForBG = 15;
    break;
  case 3:
    fNDegreeRotationPMForBG = 20;
    break;
  default:
    cout<<"Warning: DegreesForRotationMethod not defined "<<DegreesForRotationMethod<<endl;
    return kFALSE;
  }
  fCuts[kDegreesForRotationMethod]=DegreesForRotationMethod;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetNumberOfBGEvents(Int_t NumberOfBGEvents){
  // Set Cut
  switch(NumberOfBGEvents){
  case 0:
    fNumberOfBGEvents = 5;
    break;
  case 1:
    fNumberOfBGEvents = 10;
    break;
  case 2:
    fNumberOfBGEvents = 15;
    break;
  case 3:
    fNumberOfBGEvents = 20;
    break;
  case 4:
    fNumberOfBGEvents = 2;
    break;
  case 5:
    fNumberOfBGEvents = 50;
    break;
  case 6:
    fNumberOfBGEvents = 80;
    break;
  case 7:
    fNumberOfBGEvents = 100;
    break;
  default:
    cout<<"Warning: NumberOfBGEvents not defined "<<NumberOfBGEvents<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetSharedElectronCut(Int_t sharedElec) {

  switch(sharedElec){
  case 0:
    fDoSharedElecCut = kFALSE;
    fDoMesonQualitySelection= kFALSE;
    fMesonQualityMin=0;
    break;
  case 1:
    fDoSharedElecCut = kTRUE;
    fDoMesonQualitySelection= kFALSE;
    fMesonQualityMin=0;
    break;
  case 2:
    fDoSharedElecCut = kFALSE;
    fDoMesonQualitySelection= kTRUE;
    fMesonQualityMin=2;
    break;

  default:
    cout<<"Warning: Shared Electron Cut not defined "<<sharedElec<<endl;
    return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetToCloseV0sCut(Int_t toClose) {

  if(!fEnableOmegaAPlikeCut){                                                   // this cut is overloaded for the omega to pizero gamma analysis
    switch(toClose){
    case 0:
      fDoToCloseV0sCut  = kFALSE;
      fMinV0Dist        = 250;
      break;
    case 1:
      fDoToCloseV0sCut  = kTRUE;
      fMinV0Dist        = 1;
      break;
    case 2:
      fDoToCloseV0sCut  = kTRUE;
      fMinV0Dist        = 2;
      break;
    case 3:
      fDoToCloseV0sCut  = kTRUE;
      fMinV0Dist        = 3;
      break;
    default:
      cout<<"Warning: To Close To V0 Cut not defined "<<toClose<<endl;
      return kFALSE;
    }
    return kTRUE;
  }
  else{
    switch(toClose){
    case 0:
      fAPlikeSigma        = 0;
      break;
    case 1:
      fAPlikeSigma        = 1;
      break;
    case 2:
      fAPlikeSigma        = 2;
      break;
    case 3:
      fAPlikeSigma        = 3;
      break;
    default:
      cout<<"Warning: Overloaded To Close To V0 Cut not defined "<<toClose<<endl;
      return kFALSE;
    }
    return kTRUE;
  }
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetMCPSmearing(Int_t useMCPSmearing)
{// Set Cut
  if(fMode == 2){ //PCM-EMCal running
    switch(useMCPSmearing){
    case 0:
      fUseMCPSmearing   = 0;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.;
      fPSigSmearingCte  = 0.;
      break;
    case 1:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1;
      fPSigSmearing     = 0.010;
      fPSigSmearingCte  = 0.010;
      break;
    case 2:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1;
      fPSigSmearing     = 0.015;
      fPSigSmearingCte  = 0.010;
      break;
    case 3:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.020;
      fPSigSmearingCte  = 0.010;
      break;
    case 4:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.020;
      fPSigSmearingCte  = 0.020;
      break;
    case 5:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.030;
      fPSigSmearingCte  = 0.020;
      break;
    case 6:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.030;
      fPSigSmearingCte  = 0.030;
      break;
    case 7:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1;
      fPSigSmearing     = 0.030;
      fPSigSmearingCte  = 0.050;
      break;
    case 8:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.030;
      fPSigSmearingCte  = 0.060;
      break;
    case 9:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.075;
      fPSigSmearingCte  = 0.050;
      break;
    case 10:     //a
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.0275;
      fPSigSmearingCte  = 0.025;
      break;
    case 11:     //b
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.025;
      fPSigSmearingCte  = 0.030;
      break;
    case 12:     //c
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.0275;
      fPSigSmearingCte  = 0.020;
      break;
    case 13:     //d
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.0275;
      fPSigSmearingCte  = 0.035;
      break;
    case 14:     //e
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.0275;
      fPSigSmearingCte  = 0.04;
      break;
    case 15:     //f
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.0275;
      fPSigSmearingCte  = 0.015;
      break;
    case 16:     //g
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.025;
      fPSigSmearingCte  = 0.020;
      break;
    case 17:     //h
      fUseMCPSmearing   = 2;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.025;
      fPSigSmearingCte  = 0.020;
      break;

    default:
      AliError("Warning: UseMCPSmearing not defined");
      return kFALSE;
    }
  }else{
    switch(useMCPSmearing){
    case 0:
      fUseMCPSmearing   = 0;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.;
      fPSigSmearingCte  = 0.;
      break;
    case 1:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.0e-14;
      fPSigSmearing     = 0.;
      fPSigSmearingCte  = 0.;
      break;
    case 2:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.0e-15;
      fPSigSmearing     = 0.0;
      fPSigSmearingCte  = 0.;
      break;
    case 3:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.003;
      fPSigSmearingCte  = 0.002;
      break;
    case 4:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.003;
      fPSigSmearingCte  = 0.007;
      break;
    case 5:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.003;
      fPSigSmearingCte  = 0.016;
      break;
    case 6:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.007;
      fPSigSmearingCte  = 0.016;
      break;
    case 7:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.0e-16;
      fPSigSmearing     = 0.0;
      fPSigSmearingCte  = 0.;
      break;
    case 8:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.007;
      fPSigSmearingCte  = 0.014;
      break;
    case 9:
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.007;
      fPSigSmearingCte  = 0.011;
      break;
   case 10:     //a
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.00;
      fPSigSmearingCte  = 0.02;
      break;
    case 11:     //b
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.00;
      fPSigSmearingCte  = 0.025;
      break;
    case 12:     //c
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.00;
      fPSigSmearingCte  = 0.030;
      break;
    case 13:     //d
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.00;
      fPSigSmearingCte  = 0.01;
      break;
    case 14:     //e
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.00;
      fPSigSmearingCte  = 0.008;
      break;
    case 15:     //f
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.002;
      fPSigSmearingCte  = 0.008;
      break;
    case 16:     //g
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.001;
      fPSigSmearingCte  = 0.008;
      break;
    case 17:     //h
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.00;
      fPSigSmearingCte  = 0.011;
      break;
    case 18:     //i
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.00;
      fPSigSmearingCte  = 0.012;
      break;
    case 19:     //j
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.002;
      fPSigSmearingCte  = 0.01;
      break;
    case 20:     //k
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.002;
      fPSigSmearingCte  = 0.011;
      break;
    case 21:     //l
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.002;
      fPSigSmearingCte  = 0.012;
      break;
    case 22:     //m
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.004;
      fPSigSmearingCte  = 0.011;
      break;
    case 23:     //n
      fUseMCPSmearing   = 1;
      fPBremSmearing    = 1.;
      fPSigSmearing     = 0.009;
      fPSigSmearingCte  = 0.011;
    case 24:     //o
      fUseMCPSmearing   = 2;
      fPSigSmearing     = 0.000222327;
      fPSigSmearingCte  = -0.000123638;

      break;
    default:
      AliError("Warning: UseMCPSmearing not defined");
      return kFALSE;
    }
  }
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetDCAGammaGammaCut(Int_t DCAGammaGamma){
  // Set Cut
  switch(DCAGammaGamma){
  case 0:  //
    fDCAGammaGammaCutOn = kFALSE;
    fDCAGammaGammaCut   = 1000;
    break;
  case 1:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 10;
    break;
  case 2:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 5;
    break;
  case 3:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 4;
    break;
  case 4:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 3;
    break;
  case 5:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 2.5;
    break;
  case 6:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 2;
    break;
  case 7:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 1.5;
    break;
  case 8:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 1;
    break;
  case 9:  //
    fDCAGammaGammaCutOn = kTRUE;
    fDCAGammaGammaCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCAGammaGamma not defined "<<DCAGammaGamma<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetDCAZMesonPrimVtxCut(Int_t DCAZMesonPrimVtx){
  // Set Cut
  switch(DCAZMesonPrimVtx){
  case 0:  //
    fDCAZMesonPrimVtxCutOn = kFALSE;
    fDCAZMesonPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCAZMesonPrimVtxCutOn = kTRUE;
    fDCAZMesonPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCAZMesonPrimVtx not defined "<<DCAZMesonPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetDCARMesonPrimVtxCut(Int_t DCARMesonPrimVtx){
  // Set Cut
  switch(DCARMesonPrimVtx){
  case 0:  //
    fDCARMesonPrimVtxCutOn = kFALSE;
    fDCARMesonPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCARMesonPrimVtxCutOn = kTRUE;
    fDCARMesonPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCARMesonPrimVtx not defined "<<DCARMesonPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetMinOpanMesonCut(Int_t minOpanMesonCut){
  // Set Cut

    switch(minOpanMesonCut){
    case 0:      //
      fMinOpanCutMeson  = 0;
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 1:      //
      fMinOpanCutMeson  = 0.005;
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 2:
      if( fFMinOpanCut ) delete fFMinOpanCut;
      fFMinOpanCut      = new TF1("fFMinOpanCut","[0]*exp(-[1]*x)+[2]",0.,100.);
      fFMinOpanCut->SetParameter(0,1.5);
      fFMinOpanCut->SetParameter(1,1.35);
      fFMinOpanCut->SetParameter(2,0.02);
      fMinOpanCutMeson  = 0;
      fMinOpanPtDepCut  = kTRUE;
      break;
    case 3:      //
      fMinOpanCutMeson  = 0.01;
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 4:      //
      fMinOpanCutMeson  = 0.0152; // minimum 0.75 EMCal cell diagonals
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 5:      //
      fMinOpanCutMeson  = 0.0202; // minimum 1 EMCal cell diagonal
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 6:      //
      fMinOpanCutMeson  = 0.017; // new standard cut for EMCal analyses as of 17.05.2017
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 7:      //
      fMinOpanCutMeson  = 0.016;
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 8:      //
      fMinOpanCutMeson  = 0.018;
      fMinOpanPtDepCut  = kFALSE;
      break;
    case 9:      //
      fMinOpanCutMeson  = 0.019;
      fMinOpanPtDepCut  = kFALSE;
      break;

    //cuts with one cell dist for GammaCalo only
    case 10:      //a
      fMinOpanCutMeson  = 0.;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      break;
    case 11:      //b
      fMinOpanCutMeson  = 0.0152;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      break;
    case 12:      //c
      fMinOpanCutMeson  = 0.016;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      break;
    case 13:      //d
      fMinOpanCutMeson  = 0.017;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      break;
    case 14:      //e
      fMinOpanCutMeson  = 0.018;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      break;
    case 15:      //f
      fMinOpanCutMeson  = 0.019;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      break;
    case 16:      //g
      fMinOpanCutMeson  = 0.0202;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      break;
    case 17:      //h - Min opening angle cut for eta reconstruction - EMCEMC
      if( fFMinOpanCut ) delete fFMinOpanCut;
      fFMinOpanCut      = new TF1("fFMinOpanCut","(exp([0]*(x+1))+[1]*(x-1)+[2]-0.05)*(x<[3])+0.017*(x>[3])",0.,100.);
      fFMinOpanCut->SetParameter(0,-0.530209);
      fFMinOpanCut->SetParameter(1,-0.00536687);
      fFMinOpanCut->SetParameter(2,0.168845);
      fFMinOpanCut->SetParameter(3,20);
      fMinOpanCutMeson  = 0;
      fMinOpanPtDepCut  = kTRUE;
      break;
    case 18:      //i
      fMinOpanCutMeson  = 0.017;
      fEnableMinOpeningAngleCut = kFALSE;
      fMinOpanPtDepCut  = kFALSE;
      fEnableOneCellDistCut = kTRUE;
      fAllowCombOnlyInSameRecMethod = kTRUE;
      break;
    // opening angle cut variations for EMCal related analyses up to 17. May 2017
//    case 5:      //
//      fMinOpanCutMeson  = 0.0202; // minimum 1 EMCal cell diagonal
//      fMinOpanPtDepCut  = kFALSE;
//      break;
//    case 6:      //
//      fMinOpanCutMeson  = 0.0404; // minimum 2 EMCal cell diagonals
//      fMinOpanPtDepCut  = kFALSE;
//      break;
//    case 7:      //
//      fMinOpanCutMeson  = 0.0303; // minimum 1.5 EMCal cell diagonal
//      fMinOpanPtDepCut  = kFALSE;
//      break;
//    case 8:      //
//      fMinOpanCutMeson  = 0.02525; // minimum 1.25 EMCal cell diagonal
//      fMinOpanPtDepCut  = kFALSE;
//      break;
//    case 9:      //
//      fMinOpanCutMeson  = 0.03535; // minimum 1.75 EMCal cell diagonal
//      fMinOpanPtDepCut  = kFALSE;
//      break;
    default:
      cout<<"Warning:minOpanMesonCut  not defined "<<minOpanMesonCut<<endl;
      return kFALSE;
    }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetMaxOpanMesonCut(Int_t maxOpanMesonCut){
  // Set Cut
  switch(maxOpanMesonCut){
  case 0:      //
    fMaxOpanCutMeson  = TMath::Pi();
    fMaxOpanPtDepCut  = kFALSE;
    break;
  case 1:
    if( fFMaxOpanCut ) delete fFMaxOpanCut;
    fFMaxOpanCut      = new TF1("fFMaxOpanCut","[0]*exp(-[1]*x)+[2]",0.,100.);
    fFMaxOpanCut->SetParameter(0,2.5);
    fFMaxOpanCut->SetParameter(1,0.85);
    fFMaxOpanCut->SetParameter(2,0.35);
    fMaxOpanPtDepCut  = kTRUE;
    fMaxOpanCutMeson  = TMath::Pi();
    break;
  case 2:
    if( fFMaxOpanCut ) delete fFMaxOpanCut;
    fFMaxOpanCut      = new TF1("fFMaxOpanCut","[0]*exp(-[1]*x)+[2]",0.,100.);
    fFMaxOpanCut->SetParameter(0,2.3);
    fFMaxOpanCut->SetParameter(1,0.85);
    fFMaxOpanCut->SetParameter(2,0.35);
    fMaxOpanPtDepCut  = kTRUE;
    fMaxOpanCutMeson  = TMath::Pi();
    break;
  case 3:      // - Max opening angle cut for pi0 & eta reconstruction - EMCEMC
      if( fFMaxOpanCut ) delete fFMaxOpanCut;
      fFMaxOpanCut      = new TF1("fFMaxOpanCut","exp([0]*(x-0.5))+[1]*(x-0.5)+[2]+0.15",0.,100.);
      fFMaxOpanCut->SetParameter(0,-0.530209);
      fFMaxOpanCut->SetParameter(1,-0.00536687);
      fFMaxOpanCut->SetParameter(2,0.168845);
      fMaxOpanPtDepCut  = kTRUE;
      fMaxOpanCutMeson  = TMath::Pi();
      break;
  default:
    cout<<"Warning: maxOpanMesonCut not defined "<< maxOpanMesonCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
TString AliConversionMesonCuts::GetCutNumber(){
  // returns TString with current cut number
  return fCutStringRead;
}

//________________________________________________________________________
void AliConversionMesonCuts::FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0){

  Int_t posLabel = photon->GetTrackLabelPositive();
  Int_t negLabel = photon->GetTrackLabelNegative();

  fElectronLabelArray[nV0*2] = posLabel;
  fElectronLabelArray[(nV0*2)+1] = negLabel;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s){

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
//________________________________________________________________________
Bool_t AliConversionMesonCuts::RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0){
  Double_t posX = photon->GetConversionX();
  Double_t posY = photon->GetConversionY();
  Double_t posZ = photon->GetConversionZ();

  for(Int_t i = 0;i<photons->GetEntries();i++){
    if(nV0 == i) continue;
    AliAODConversionPhoton *photonComp = (AliAODConversionPhoton*) photons->At(i);
    Double_t posCompX = photonComp->GetConversionX();
    Double_t posCompY = photonComp->GetConversionY();
    Double_t posCompZ = photonComp->GetConversionZ();

    Double_t dist = pow((posX - posCompX),2)+pow((posY - posCompY),2)+pow((posZ - posCompZ),2);

    if(dist < fMinV0Dist*fMinV0Dist){
      if(photon->GetChi2perNDF() < photonComp->GetChi2perNDF()) return kTRUE;
      else {
        return kFALSE;}
    }

  }
  return kTRUE;
}

//________________________________________________________________________
void AliConversionMesonCuts::SmearParticle(AliAODConversionPhoton* photon)
{

  if (photon==NULL) return;
  Double_t facPBrem = 1.;
  Double_t facPSig = 0.;

  Double_t phi=0.;
  Double_t theta=0.;
  Double_t P=0.;


  P=photon->P();
  phi=photon->Phi();
  if( photon->P()!=0){
    theta=acos( photon->Pz()/ photon->P());
  }

  if (fUseMCPSmearing == 1){
    if( fPSigSmearing != 0. || fPSigSmearingCte!=0. ){
      facPSig = TMath::Sqrt(fPSigSmearingCte*fPSigSmearingCte+fPSigSmearing*fPSigSmearing*P*P)*fRandom.Gaus(0.,1.);
    }

    if( fPBremSmearing != 1.){
      if(fBrem!=NULL){
        facPBrem = fBrem->GetRandom();
      }
    }

    photon->SetPx(facPBrem* (1+facPSig)* P*sin(theta)*cos(phi)) ;
    photon->SetPy(facPBrem* (1+facPSig)* P*sin(theta)*sin(phi)) ;
    photon->SetPz(facPBrem* (1+facPSig)* P*cos(theta)) ;
    photon->SetE(photon->P());
  } else if (fUseMCPSmearing == 2) {
    TF1* fScaling = new TF1("fScaling", "pol2", 0, 100);
    fScaling->SetParameters(1.47693, -0.954702, 0.47791);    
    facPSig = fRandom.Gaus(P,fScaling->GetX(fPSigSmearingCte+fPSigSmearing*P));
    
    photon->SetPx(facPSig*sin(theta)*cos(phi)) ;
    photon->SetPy(facPSig*sin(theta)*sin(phi)) ;
    photon->SetPz(facPSig*cos(theta));

    photon->SetE(photon->P());
  }
  

}
//________________________________________________________________________
void AliConversionMesonCuts::SmearVirtualPhoton(AliAODConversionPhoton* photon)
{

  if (photon==NULL) return;
  Double_t facPBrem = 1.;
  Double_t facPSig = 0.;

  Double_t phi=0.;
  Double_t theta=0.;
  Double_t P=0.;


  P=photon->P();
  phi=photon->Phi();
  if( photon->P()!=0){
    theta=acos( photon->Pz()/ photon->P());
  }

  if( fPSigSmearing != 0. || fPSigSmearingCte!=0. ){
    facPSig = TMath::Sqrt(fPSigSmearingCte*fPSigSmearingCte+fPSigSmearing*fPSigSmearing*P*P)*fRandom.Gaus(0.,1.);
  }

  if( fPBremSmearing != 1.){
    if(fBrem!=NULL){
      facPBrem = fBrem->GetRandom();
    }
  }

  photon->SetPx(facPBrem* (1+facPSig)* P*sin(theta)*cos(phi)) ;
  photon->SetPy(facPBrem* (1+facPSig)* P*sin(theta)*sin(phi)) ;
  photon->SetPz(facPBrem* (1+facPSig)* P*cos(theta)) ;

}
//________________________________________________________________________
TLorentzVector AliConversionMesonCuts::SmearElectron(TLorentzVector particle)
{

  //if (particle==0) return;
  Double_t facPBrem = 1.;
  Double_t facPSig = 0.;

  Double_t phi=0.;
  Double_t theta=0.;
  Double_t P=0.;

  P=particle.P();
  phi=particle.Phi();
  if (phi < 0.) phi += 2. * TMath::Pi();

  if( particle.P()!=0){
    theta=acos( particle.Pz()/ particle.P());
  }


  Double_t fPSigSmearingHalf    =  fPSigSmearing  / 2.0;  //The parameter was set for gammas with 2 particles and here we have just one electron
  Double_t sqrtfPSigSmearingCteHalf =  fPSigSmearingCte / 2.0 ;  //The parameter was set for gammas with 2 particles and here we have just one electron



  if( fPSigSmearingHalf != 0. || sqrtfPSigSmearingCteHalf!=0. ){
    facPSig = TMath::Sqrt(sqrtfPSigSmearingCteHalf*sqrtfPSigSmearingCteHalf+fPSigSmearingHalf*fPSigSmearingHalf*P*P)*fRandom.Gaus(0.,1.);
  }

  if( fPBremSmearing != 1.){
    if(fBrem!=NULL){
      facPBrem = fBrem->GetRandom();
    }
  }

  TLorentzVector SmearedParticle;

  SmearedParticle.SetXYZM( facPBrem* (1+facPSig)* P*sin(theta)*cos(phi) , facPBrem* (1+facPSig)* P*sin(theta)*sin(phi)  ,
        facPBrem* (1+facPSig)* P*cos(theta) , TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass()) ;

  //particle.SetPx(facPBrem* (1+facPSig)* P*sin(theta)*cos(phi)) ;
  //particle.SetPy(facPBrem* (1+facPSig)* P*sin(theta)*sin(phi)) ;
  //particle.SetPz(facPBrem* (1+facPSig)* P*cos(theta)) ;

  return SmearedParticle;

}

//________________________________________________________________________
// function to determine whether meson was selected by mass range
Bool_t AliConversionMesonCuts::MesonIsSelectedByMassCut(AliAODConversionMother *meson, Int_t nominalRange = 0){

  if (fAcceptMesonMass){
      if (nominalRange == 0){
        if (meson->M() > fSelectionLow && meson->M() < fSelectionHigh)
          return kTRUE;
        else
          return kFALSE;
      } else if (nominalRange == 1){
        if (meson->M() > fSidebandMixingLow && meson->M() < fSidebandMixingHigh)
          return kTRUE;
        else
          return kFALSE;
      } else if (nominalRange == 2){
        if (meson->M() > fSidebandMixingLeftLow && meson->M() < fSidebandMixingLeftHigh)
          return kTRUE;
        else
          return kFALSE;
      } else if (nominalRange == 3){
        if (meson->M() > fSidebandMixingRightLow && meson->M() < fSidebandMixingRightHigh)
          return kTRUE;
        else
          return kFALSE;
      }
  } else if (fUsePtDepSelectionWindow){
      // Determine correct mass parametrisation depending on what method is used
      Float_t pt   = meson->Pt();
      Float_t mass = 0;
      Float_t sigma = 999;
      switch(fMassParamFunction){
        case 0: // EMC-EMC
          mass = 0.125306 + 0.001210 * pt;
          sigma =   0.0136138 + ( (-0.00104914) * pt ) + (7.61163e-05 * pt * pt);
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 1: // PCM-EMC
          mass = 0.129756 + 0.000660514 * pt;
          sigma =   0.00990291 + ( ( -0.00114665) * pt ) + (0.000128015 * pt * pt);
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 2: // PHOS-PHOS
          mass = 0.132298 + ( 9.84713e-05 * pt );
          sigma =   4.15716e-03 + ( 3.35025e-03 / pt );
          if (pt>5.){
              sigma+=1.59496e-04*(pt-5);
          }
          if (sigma < 0.004 ) {sigma =0.004;}
          else if (sigma > 0.02) {sigma =0.02;}
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 3: // PCM-PHOS
          mass = 0.13344 - ( (-1.26101e-05) * pt );
          sigma =   3.56197e-03 + ( 7.31591e-04 / pt );
          if (pt>5.){
              sigma+=4.77381e-04*(pt-5);
          }
          if (sigma < 0.0025 ) {sigma =0.0025;}
          else if (sigma > 0.02) {sigma =0.02;}
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 4: // PCM-PCM
          mass = 0.134613 + (-0.000154418 * pt);
          sigma =   0.00223215 + ( (0.000349362) * pt ) + (-1.13689e-05 * pt * pt);
          if (sigma < 0.001 ) {sigma =0.001;}
          else if (sigma > 0.007) {sigma =0.007;}
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 5: // EMC-EMC (optimized for 7 TeV with 31 NonLin)
          mass = 0.122498 + (0.004978 * pt) + (-0.000850 * pow(pt,2.)) + (0.000061 * pow(pt,3)) + (-13e-07 * pow(pt,4));
          sigma =  0.006118 + 0.000373 * pt+ 0.010567 /pt;
          if(mass < 0.123) mass = 0.123;
          if(sigma > 0.017) sigma = 0.017;
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 6: // PCM-EMC (optimized for 7 TeV with 31 NonLin)
          mass = 0.131116 + 0.000420 * pt;
          sigma =0.006128 + 0.003981 / pt;
          if(sigma>0.012) sigma = 0.012;
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 7: // PHOS-PHOS (optimized for 5 TeV with correct nonlin Run2)
          mass = 0.134392 + ( 3.26311e-05  * pt );
          sigma =   4.90919e-03 + (5.27357e-03 * exp(-8.65337e-01 * pt));
          if (sigma < 0.002 ) {sigma =0.002;}
          else if (sigma > 0.02) {sigma =0.02;}
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 8: // PCM-PHOS (optimized for 5 TeV with correct nonlin Run2)
          mass =  1.34073e-01  + (7.79463e-03 * exp(-2.77158 * pt));
          sigma =  3.75986e-03 + (6.41965e-03  * exp(-2.38570 * pt));
          if (mass>0.138) mass = 0.138;
          if (sigma < 0.002 ) {sigma =0.002;}
          else if (sigma > 0.02) {sigma =0.02;}
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        case 9: // PCM-PCM (optimized for 7 TeV)
          mass = 0.133587 + (0.002594 * exp(-1.189987 * pt)) + -0.000035 * pt;
          sigma =   0.00223215 + ( (0.000349362) * pt ) + (-1.13689e-05 * pt * pt);
          if (mass>0.137) mass = 0.137;
          if (sigma < 0.001 ) {sigma =0.001;}
          else if (sigma > 0.007) {sigma =0.007;}
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
        default:
          mass = 0.125306 + 0.001210 * pt;
          sigma =   0.0136138 + ( (-0.00104914) * pt ) + (7.61163e-05 * pt * pt);
          fSelectionLow = mass - (fSelectionNSigmaLow * sigma);
          fSelectionHigh = mass + (fSelectionNSigmaHigh * sigma);
          break;
      }

      if (nominalRange == 0){
        if (meson->M() > fSelectionLow && meson->M() < fSelectionHigh)
          return kTRUE;
        else
          return kFALSE;
      } else if (nominalRange == 1){
        if (meson->M() > fSidebandMixingLow && meson->M() < fSidebandMixingHigh)
          return kTRUE;
        else
          return kFALSE;
      } else if (nominalRange == 2){
        if (meson->M() > fSidebandMixingLeftLow && meson->M() < fSidebandMixingLeftHigh)
          return kTRUE;
        else
          return kFALSE;
      } else if (nominalRange == 3){
        if (meson->M() > fSidebandMixingRightLow && meson->M() < fSidebandMixingRightHigh)
          return kTRUE;
        else
          return kFALSE;
      }
  } else {
    if (!(meson->M() > fSelectionLow && meson->M() < fSelectionHigh))
      return kTRUE;
    else
      return kFALSE;
  }
  return kTRUE;
}

///_____________________________________________________________________________
// Armenteros Qt Cut
Bool_t AliConversionMesonCuts::ArmenterosLikeQtCut(Double_t alpha, Double_t qT){
  if(fAPlikeSigma){
    if( (alpha < 1.0) && (alpha > -0.97) ){
      if(fAPlikeSigma == 1){
        if( (qT < sqrt( (1.-pow(alpha-0.0298, 2.)/pow(1., 2.))*pow(0.43, 2) ) ) &&
        (qT > 0.313 - 0.084*alpha - 0.327 * pow(alpha, 2.) + 0.056 * pow(alpha, 3.) + 0.022 *pow(alpha, 4.) ) )
          return kTRUE;
        else
          return kFALSE;
      }
      else if(fAPlikeSigma == 2){
        if( (qT < sqrt( (1.-pow(alpha-0.0298, 2.)/pow(1., 2.))*pow(0.475, 2) ) ) &&
        (qT > 0.260 - 0.042*alpha - 0.255 * pow(alpha, 2.) + 0.045 * pow(alpha, 3.) - 0.023 *pow(alpha, 4.) ) )
          return kTRUE;
        else
          return kFALSE;
      }
      else if(fAPlikeSigma == 3){
        if( (qT < sqrt( (1.-pow(alpha-0.0298, 2.)/pow(1., 2.))*pow(0.52, 2) ) ) &&
        (qT > 0.210 - 0.052*alpha - 0.316 * pow(alpha, 2.) + 0.063 * pow(alpha, 3.) + 0.049 *pow(alpha, 4.) ) )
          return kTRUE;
        else
          return kFALSE;
      }
      else{
        return kTRUE;
      }
    }
    else{
      return kFALSE;
    }
  }
  else{
    return kTRUE;
  }
}
