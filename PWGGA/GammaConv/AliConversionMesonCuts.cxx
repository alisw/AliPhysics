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
  "RCut", //5
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
  fHistograms(NULL),
  fDoLightOutput(kFALSE),
  fMode(0),
  fCaloPhotonCuts(NULL),
  fMesonKind(0),
  fIsMergedClusterCut(0),
  fMaxR(200),
  fEnableMassCut(kFALSE),
  fSelectionLow(0.08),
  fSelectionHigh(0.145),
  fSelectionWindowCut(-1),
  fAlphaMinCutMeson(0),
  fAlphaCutMeson(1),
  fRapidityCutMeson(1),
  fUseRotationMethodInBG(kFALSE),
  fUsePtmaxMethodForBG(kFALSE),
  fDoBG(kTRUE),
  fdoBGProbability(kFALSE),
  fUseTrackMultiplicityForBG(kFALSE),
  fnDegreeRotationPMForBG(0),
  fNumberOfBGEvents(0),
  fOpeningAngle(0.005),
  fEnableMinOpeningAngleCut(kTRUE),
  fEnableOneCellDistCut(kFALSE),
  fDoToCloseV0sCut(kFALSE),
  fminV0Dist(200.),
  fDoSharedElecCut(kFALSE),
  fUseMCPSmearing(kFALSE),
  fPBremSmearing(0),
  fPSigSmearing(0),
  fPSigSmearingCte(0),
  fBrem(NULL),
  fRandom(0),
  fFAlphaCut(0),
  fAlphaPtDepCut(kFALSE),
  fElectronLabelArraySize(500),
  fElectronLabelArray(NULL),
  fDCAGammaGammaCut(1000),
  fDCAZMesonPrimVtxCut(1000),
  fDCARMesonPrimVtxCut(1000),
  fDCAGammaGammaCutOn(kFALSE),
  fDCAZMesonPrimVtxCutOn(kFALSE),
  fDCARMesonPrimVtxCutOn(kFALSE),
  fMinOpanCutMeson(0),
  fFMinOpanCut(0),
  fMinOpanPtDepCut(kFALSE),
  fMaxOpanCutMeson(TMath::Pi()),
  fFMaxOpanCut(0),
  fMaxOpanPtDepCut(kFALSE),
  fCutString(NULL),
  fCutStringRead(""),
  fBackgroundHandler(0),
  fHistoMesonCuts(NULL),
  fHistoMesonBGCuts(NULL),
  fHistoDCAGGMesonBefore(NULL),
  fHistoDCAZMesonPrimVtxBefore(NULL),
  fHistoDCARMesonPrimVtxBefore(NULL),
  fHistoDCAGGMesonAfter(NULL),
  fHistoDCAZMesonPrimVtxAfter(NULL),
  fHistoDCARMesonPrimVtxAfter(NULL),
  fHistoInvMassBefore(NULL),
  fHistoInvMassAfter(NULL)
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
  fHistograms(NULL),
  fDoLightOutput(ref.fDoLightOutput),
  fMode(ref.fMode),
  fCaloPhotonCuts(ref.fCaloPhotonCuts),
  fMesonKind(ref.fMesonKind),
  fIsMergedClusterCut(ref.fIsMergedClusterCut),
  fMaxR(ref.fMaxR),
  fEnableMassCut(ref.fEnableMassCut),
  fSelectionLow(ref.fSelectionLow),
  fSelectionHigh(ref.fSelectionHigh),
  fSelectionWindowCut(-1),
  fAlphaMinCutMeson(ref.fAlphaMinCutMeson),
  fAlphaCutMeson(ref.fAlphaCutMeson),
  fRapidityCutMeson(ref.fRapidityCutMeson),
  fUseRotationMethodInBG(ref.fUseRotationMethodInBG),
  fUsePtmaxMethodForBG(ref.fUsePtmaxMethodForBG),
  fDoBG(ref.fDoBG),
  fdoBGProbability(ref.fdoBGProbability),
  fUseTrackMultiplicityForBG(ref.fUseTrackMultiplicityForBG),
  fnDegreeRotationPMForBG(ref.fnDegreeRotationPMForBG),
  fNumberOfBGEvents(ref. fNumberOfBGEvents),
  fOpeningAngle(ref.fOpeningAngle),
  fEnableMinOpeningAngleCut(ref.fEnableMinOpeningAngleCut),
  fEnableOneCellDistCut(ref.fEnableOneCellDistCut),
  fDoToCloseV0sCut(ref.fDoToCloseV0sCut),
  fminV0Dist(ref.fminV0Dist),
  fDoSharedElecCut(ref.fDoSharedElecCut),
  fUseMCPSmearing(ref.fUseMCPSmearing),
  fPBremSmearing(ref.fPBremSmearing),
  fPSigSmearing(ref.fPSigSmearing),
  fPSigSmearingCte(ref.fPSigSmearingCte),
  fBrem(NULL),
  fRandom(ref.fRandom),
  fFAlphaCut(NULL),
  fAlphaPtDepCut(ref.fAlphaPtDepCut),
  fElectronLabelArraySize(ref.fElectronLabelArraySize),
  fElectronLabelArray(NULL),
  fDCAGammaGammaCut(ref.fDCAGammaGammaCut),
  fDCAZMesonPrimVtxCut(ref.fDCAZMesonPrimVtxCut),
  fDCARMesonPrimVtxCut(ref.fDCARMesonPrimVtxCut),
  fDCAGammaGammaCutOn(ref.fDCAGammaGammaCutOn),
  fDCAZMesonPrimVtxCutOn(ref.fDCAZMesonPrimVtxCutOn),
  fDCARMesonPrimVtxCutOn(ref.fDCARMesonPrimVtxCutOn),
  fBackgroundHandler(ref.fBackgroundHandler),
  fFMinOpanCut(0),
  fMinOpanPtDepCut(kFALSE),
  fMaxOpanCutMeson(TMath::Pi()),
  fFMaxOpanCut(0),
  fMaxOpanPtDepCut(kFALSE),
  fMinOpanCutMeson(0),
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
  fHistoInvMassAfter(NULL)
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
  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("ConvMesonCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  // Meson Cuts
  if (fIsMergedClusterCut == 1){
    fHistoMesonCuts=new TH2F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts vs Pt",8,-0.5,7.5, 500, 0, 100);
    fHistoMesonCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(4,"mass cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(5,"opening angle");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(6,"alpha max");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(7,"alpha min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(8,"out");
    fHistograms->Add(fHistoMesonCuts);
  } else if (fIsMergedClusterCut == 2){
    fHistoMesonCuts=new TH2F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts vs Pt",8,-0.5,7.5, 250, 0, 50);
    fHistoMesonCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(4,"1 cell distance");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(5,"opening angle");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(6,"alpha max");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(7,"alpha min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(8,"out");
    fHistograms->Add(fHistoMesonCuts);

    fHistoMesonBGCuts=new TH2F(Form("MesonBGCuts %s",GetCutNumber().Data()),"MesonBGCuts vs Pt",8,-0.5,7.5, 250, 0, 50);
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(4,"1 cell distance");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(5,"opening angle");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(6,"alpha max");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(7,"alpha min");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(8,"out");
    fHistograms->Add(fHistoMesonBGCuts);
  } else {
    fHistoMesonCuts=new TH2F(Form("MesonCuts %s",GetCutNumber().Data()),"MesonCuts vs Pt",10,-0.5,9.5, 250, 0, 50);
    fHistoMesonCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(4,"opening angle");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(5,"alpha max");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(6,"alpha min");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(7,"dca gamma gamma");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(8,"dca R prim Vtx");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(9,"dca Z prim Vtx");
    fHistoMesonCuts->GetXaxis()->SetBinLabel(10,"out");
    fHistograms->Add(fHistoMesonCuts);

    fHistoMesonBGCuts=new TH2F(Form("MesonBGCuts %s",GetCutNumber().Data()),"MesonBGCuts vs Pt",10,-0.5,9.5, 250, 0, 50);
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(1,"in");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(2,"undef rapidity");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(3,"rapidity cut");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(4,"opening angle");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(5,"alpha max");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(6,"alpha min");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(7,"dca gamma gamma");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(8,"dca R prim Vtx");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(9,"dca Z prim Vtx");
    fHistoMesonBGCuts->GetXaxis()->SetBinLabel(10,"out");
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

  if(fMCMother->GetPdgCode()==111 || fMCMother->GetPdgCode()==221){
    if(fMCMother->R()>fMaxR)  return kFALSE; // cuts on distance from collision point

    Double_t rapidity = 10.;
    if(fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0){
      rapidity=8.-fRapidityShift;
    } else{
      rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
    }

    // Rapidity Cut
    if(TMath::Abs(rapidity)>fRapidityCutMeson)return kFALSE;

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

  if(MCMother->GetPdgCode()==111 || MCMother->GetPdgCode()==221){
    Double_t rMeson = sqrt( (MCMother->Xv()*MCMother->Xv()) + (MCMother->Yv()*MCMother->Yv()) ) ;
    if(rMeson>fMaxR)  return kFALSE; // cuts on distance from collision point

    Double_t rapidity = 10.;
    if(MCMother->E() - MCMother->Pz() == 0 || MCMother->E() + MCMother->Pz() == 0){
      rapidity=8.-fRapidityShift;
    } else{
      rapidity = 0.5*(TMath::Log((MCMother->E()+MCMother->Pz()) / (MCMother->E()-MCMother->Pz())))-fRapidityShift;
    }

    // Rapidity Cut
    if(TMath::Abs(rapidity)>fRapidityCutMeson)return kFALSE;

    // Select only -> 2y decay channel
    if(MCMother->GetNDaughters()!=2)return kFALSE;

    for(Int_t i=0;i<2;i++){
      AliAODMCParticle *MDaughter=static_cast<AliAODMCParticle*>(AODMCArray->At(MCMother->GetDaughter(i)));
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
Bool_t AliConversionMesonCuts::MesonIsSelectedMCDalitz(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelelectron, Int_t &labelpositron, Int_t &labelgamma, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !mcEvent )return kFALSE;

  if(  fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 ) return kFALSE;

  if(  fMCMother->R()>fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( TMath::Abs(rapidity) > fRapidityCutMeson )return kFALSE;

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

  if(  fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 ) return kFALSE;

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
  if( TMath::Abs(rapidity) > fRapidityCutMeson )return kFALSE;

  // Select only -> Dalitz decay channel
  if( fMCMother->GetNDaughters() != 3 )return kFALSE;

  AliAODMCParticle *positron = 0x0;
  AliAODMCParticle *electron = 0x0;
  AliAODMCParticle    *gamma = 0x0;

  for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
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
  if( TMath::Abs(rapidity) > fRapidityCutMeson )return kFALSE;

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
Bool_t AliConversionMesonCuts::MesonIsSelectedMCPiPlPiMiPiZero(TParticle *fMCMother,AliMCEvent *mcEvent, Int_t &labelNegPion, Int_t &labelPosPion, Int_t &labelNeutPion, Double_t fRapidityShift){

  // Returns true for all pions within acceptance cuts for decay into 2 photons
  // If bMCDaughtersInAcceptance is selected, it requires in addition that both daughter photons are within acceptance cuts

  if( !mcEvent )return kFALSE;

  if( !(fMCMother->GetPdgCode() == 221 || fMCMother->GetPdgCode() == 223) ) return kFALSE;
  
  if( fMCMother->R()>fMaxR ) return kFALSE; // cuts on distance from collision point

  Double_t rapidity = 10.;

  if( fMCMother->Energy() - fMCMother->Pz() == 0 || fMCMother->Energy() + fMCMother->Pz() == 0 ){
    rapidity=8.-fRapidityShift;
  }
  else{
    rapidity = 0.5*(TMath::Log((fMCMother->Energy()+fMCMother->Pz()) / (fMCMother->Energy()-fMCMother->Pz())))-fRapidityShift;
  }

  // Rapidity Cut
  if( TMath::Abs(rapidity) > fRapidityCutMeson )return kFALSE;

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
  if(TMath::Abs(rapidity) > fRapidityCutMeson)return kFALSE;

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
Bool_t AliConversionMesonCuts::MesonIsSelectedPiZeroGammaAngle(AliAODConversionMother *omega, AliAODConversionMother *pi0, AliAODConversionPhoton *gamma, Bool_t DoPiZeroAngleCut, TF1 *maxfit, Double_t lowerFactor, Double_t upperFactor){

  if(!DoPiZeroAngleCut) return kTRUE;

  Double_t PiZeroGammaAngle = pi0->Angle(gamma->Vect());
  Double_t omegaPt = omega->Pt();

  if(PiZeroGammaAngle > lowerFactor * maxfit->Eval(omegaPt) && PiZeroGammaAngle < upperFactor * maxfit->Eval(omegaPt)) return kTRUE;
  return kFALSE;

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
    if(TMath::Abs(rapidity)>fRapidityCutMeson)return kFALSE;

    // Select only -> ChiC radiative (JPsi+gamma) decay channel
    if(fMCMother->GetNDaughters()!=2)return kFALSE;

    TParticle *jpsi   = 0x0;
    TParticle *gamma   = 0x0;
    TParticle *positron = 0x0;
    TParticle *electron = 0x0;

    Int_t labeljpsiChiC = -1;

    for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){
      if(index < 0) continue;
      TParticle* temp = (TParticle*)mcEvent->Particle( index );

      switch( temp->GetPdgCode() ) {
      case 443:
        jpsi =  temp;
        labeljpsiChiC = index;
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
Bool_t AliConversionMesonCuts::MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal, Double_t fRapidityShift, Int_t leadingCellID1, Int_t leadingCellID2)
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

  // Undefined Rapidity -> Floating Point exception
  if((pi0->E()+pi0->Pz())/(pi0->E()-pi0->Pz())<=0){
    if(hist)hist->Fill(cutIndex, pi0->Pt());
    cutIndex++;
    if (!IsSignal)cout << "undefined rapidity" << endl;
    return kFALSE;
  }
  else{
    // PseudoRapidity Cut --> But we cut on Rapidity !!!
    cutIndex++;
    if(TMath::Abs(pi0->Rapidity()-fRapidityShift)>fRapidityCutMeson){
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
  AliInfo(Form("Set Meson Cutnumber: %s",analysisCutSelection.Data()));
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
  case kRCut:
    if( SetRCut(value)) {
      fCuts[kRCut] = value;
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
  printf("\t |y| < %3.2f \n", fRapidityCutMeson);
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
    printf("\t Meson selection window for further analysis %3.3f > M_{gamma,gamma} > %3.3f\n\n", fSelectionLow, fSelectionHigh);
  }
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
    if (fdoBGProbability) printf("\t -> use BG probability \n");
    if (fBackgroundHandler) printf("\t -> use new BG handler \n");
    printf("\t depth of pool: %d\n", fNumberOfBGEvents);
    if (fUseRotationMethodInBG )printf("\t degree's for BG rotation: %d\n", fnDegreeRotationPMForBG);
    if (!fUseRotationMethodInBG  & !fUseTrackMultiplicityForBG & fBackgroundHandler) printf("\t BG scheme: event plane angle with V0 mult \n");
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
    break;
  default:
    cout<<"Warning: Meson kind not defined"<<mesonKind<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetRCut(Int_t rCut){
  // Set Cut
  switch(rCut){
  case 0:
    fMaxR = 180.;
    break;
  case 1:
    fMaxR = 180.;
    break;
  case 2:
    fMaxR = 180.;
    break;
  case 3:
    fMaxR = 70.;
    break;
  case 4:
    fMaxR = 70.;
    break;
  case 5:
    fMaxR = 180.;
    break;
    // High purity cuts for PbPb
  case 9:
    fMaxR = 180.;
    break;
  default:
    cout<<"Warning: rCut not defined"<<rCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetSelectionWindowCut(Int_t selectionCut){
  // Set Cut
  switch(selectionCut){
  case 0:  
    fSelectionLow   = 0.08;
    fSelectionHigh  = 0.145;
    break;
  case 1:  
    fSelectionLow   = 0.1;
    fSelectionHigh  = 0.145;
    break;
  case 2:  
    fSelectionLow   = 0.11;
    fSelectionHigh  = 0.145;
    break;
  case 3:
    fSelectionLow   = 0.12;
    fSelectionHigh  = 0.145;
    break;
  case 4:
    fSelectionLow   = 0.1;
    fSelectionHigh  = 0.15;
    break;
  case 5:
    fSelectionLow   = 0.11;
    fSelectionHigh  = 0.15;
    break;
  case 6:
    fSelectionLow   = 0.12;
    fSelectionHigh  = 0.15;
    break;
  case 7:
    fSelectionLow   = 0.1;
    fSelectionHigh  = 0.155;
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
    fRapidityCutMeson   = 1.35;
    break;
  case 1:  //
    fRapidityCutMeson   = 0.8;
    break;
  case 2:  //
    fRapidityCutMeson   = 0.7;
    break;
  case 3:  //
    fRapidityCutMeson   = 0.6;
    break;
  case 4:  //
    fRapidityCutMeson   = 0.5;
    break;
  case 5:  //
    fRapidityCutMeson   = 0.85;
    break;
  case 6:  //
    fRapidityCutMeson   = 0.75;
    break;
  case 7:  //
    fRapidityCutMeson   = 0.3;
    break;
  case 8:  //changed, before 0.35
    fRapidityCutMeson   = 0.25;
    break;
  case 9:  //
    fRapidityCutMeson   = 0.4;
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
    fdoBGProbability            = kFALSE;
    break;
  case 1: // mixed event with V0 multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fdoBGProbability            = kFALSE;
    break;
  case 2: // mixed event with track multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kTRUE;
    fdoBGProbability            = kFALSE;
    break;
  case 3: //Rotation
    fUseRotationMethodInBG      = kTRUE;
    fdoBGProbability            = kTRUE;
    break;
  case 4: //No BG calculation
    cout << "no BG calculation should be done" << endl;
    fUseRotationMethodInBG      = kFALSE;
    fdoBGProbability            = kFALSE;
    fDoBG                       = kFALSE;
    break;
  case 5: //Rotation
    fUseRotationMethodInBG      = kTRUE;
    fdoBGProbability            = kFALSE;
    fBackgroundHandler          = 1;
    break;
  case 6: // mixed event with V0 multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fdoBGProbability            = kFALSE;
    fBackgroundHandler          = 1;
    break;
  case 7: // mixed event with track multiplicity
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kTRUE;
    fdoBGProbability            = kFALSE;
    fBackgroundHandler          = 1;
    break;
  case 8: //Rotation
    fUseRotationMethodInBG      = kTRUE;
    fdoBGProbability            = kTRUE;
    fBackgroundHandler          = 1;
    break;
  case 9: // mixed event with Ptmax method
    fUseRotationMethodInBG      = kFALSE;
    fUseTrackMultiplicityForBG  = kFALSE;
    fdoBGProbability            = kFALSE;
    fUsePtmaxMethodForBG        = kTRUE;
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
    fnDegreeRotationPMForBG = 5;
    break;
  case 1:
    fnDegreeRotationPMForBG = 10;
    break;
  case 2:
    fnDegreeRotationPMForBG = 15;
    break;
  case 3:
    fnDegreeRotationPMForBG = 20;
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
    break;
  case 1:
    fDoSharedElecCut = kTRUE;
    break;
  default:
    cout<<"Warning: Shared Electron Cut not defined "<<sharedElec<<endl;
    return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliConversionMesonCuts::SetToCloseV0sCut(Int_t toClose) {

  switch(toClose){
  case 0:
    fDoToCloseV0sCut  = kFALSE;
    fminV0Dist        = 250;
    break;
  case 1:
    fDoToCloseV0sCut  = kTRUE;
    fminV0Dist        = 1;
    break;
  case 2:
    fDoToCloseV0sCut  = kTRUE;
    fminV0Dist        = 2;
    break;
  case 3:
    fDoToCloseV0sCut  = kTRUE;
    fminV0Dist        = 3;
    break;
  default:
    cout<<"Warning: Shared Electron Cut not defined "<<toClose<<endl;
    return kFALSE;
  }
  return kTRUE;
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

    if(dist < fminV0Dist*fminV0Dist){
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
