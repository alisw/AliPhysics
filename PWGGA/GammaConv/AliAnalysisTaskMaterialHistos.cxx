/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*																								*
* Authors: Friederike Bock															*
* Version 1.0																				*
*																								*
* Permission to use, copy, modify and distribute this software and its	 *
* documentation strictly for non-commercial purposes is hereby granted	 *
* without fee, provided that the above copyright notice appears in all	 *
* copies and that both the copyright notice and this permission notice	 *
* appear in the supporting documentation. The authors make no claims	 *
* about the suitability of this software for any purpose. It is			*
* provided "as is" without express or implied warranty.						*
**************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskMaterialHistos.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "TFile.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskMaterialHistos)

AliAnalysisTaskMaterialHistos::AliAnalysisTaskMaterialHistos() : AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fConversionGammas(NULL),
  fGammaCandidates(NULL),
  fConversionCutArray(NULL),
  fEventCutArray(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputList(NULL),
  fAllMCGammaList(NULL),
  fAllMCConvGammaList(NULL),
  fPrimVtxZ(0.),
  fNContrVtx(0),
  fNESDtracksEta09(0),
  fNESDtracksEta0914(0),
  fNESDtracksEta14(0),
  fGammaMCPt(0.),
  fGammaMCTheta(0.),
  fGammaMCConvPt(0.),
  fGammaMCConvTheta(0.),
  fGammaPt(0.),
  fGammaTheta(0.),
  fGammaChi2NDF(0.),
  fKind(0),
  fIsHeavyIon(0),
  fIsMC(0),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fnCuts(0),
  fiCut(0),
  hNEvents(NULL),
  hNGoodESDTracksEta09(NULL),
  hNGoodESDTracksEta14(NULL),
  hNGoodESDTracksEta09_14(NULL),
  hESDConversionRPhi(NULL),
  hESDConversionRZ(NULL),
  hESDConversionRPt(NULL),
  hESDConversionREta(NULL),
  hESDConversionDCA(NULL),
  hESDConversionPsiPair(NULL),
  hESDConversionChi2(NULL),
  hESDConversionMass(NULL),
  hESDConversionRRejSmall(NULL),
  hESDConversionRRejLarge(NULL),
  hElectronRdEdx(NULL),
  hElectronRNSigmadEdx(NULL),
  hPositronRdEdx(NULL),
  hPositronRNSigmadEdx(NULL),
  hMCConversionRPhi(NULL),
  hMCConversionRPt(NULL),
  hMCConversionREta(NULL),
  hMCConversionRRejSmall(NULL),
  hMCConversionRRejLarge(NULL),
  hMCAllGammaPt(NULL),
  hMCTrueConversionRPhi(NULL),
  hMCTrueConversionRZ(NULL),
  hMCTrueConversionRPt(NULL),
  hMCTrueConversionREta(NULL),
  hMCTrueConversionDCA(NULL),
  hMCTrueConversionPsiPair(NULL),
  hMCTrueConversionChi2(NULL),
  hMCTrueConversionMass(NULL),
  hMCTrueConversionRRejSmall(NULL),
  hMCTrueConversionRRejLarge(NULL),
  hMCTruePi0DalConversionRPt(NULL),
  hMCTruePi0DalConversionEta(NULL),
  hMCTrueEtaDalConversionRPt(NULL),
  hMCTrueEtaDalConversionEta(NULL),
  hMCTrueCombinatorialConversionRPt(NULL),
  hMCTrueCombinatorialConversionEta(NULL)
{

}


//________________________________________________________________________
AliAnalysisTaskMaterialHistos::AliAnalysisTaskMaterialHistos(const char *name) : AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fConversionGammas(NULL),
  fGammaCandidates(NULL),
  fConversionCutArray(NULL),
  fEventCutArray(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputList(NULL),
  fAllMCGammaList(NULL),
  fAllMCConvGammaList(NULL),
  fPrimVtxZ(0.),
  fNContrVtx(0),
  fNESDtracksEta09(0),
  fNESDtracksEta0914(0),
  fNESDtracksEta14(0),
  fGammaMCPt(0.),
  fGammaMCTheta(0.),
  fGammaMCConvPt(0.),
  fGammaMCConvTheta(0.),
  fGammaPt(0.),
  fGammaTheta(0.),
  fGammaChi2NDF(0.),
  fKind(0),
  fIsHeavyIon(0),
  fIsMC(0),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fnCuts(0),
  fiCut(0),
  hNEvents(NULL),
  hNGoodESDTracksEta09(NULL),
  hNGoodESDTracksEta14(NULL),
  hNGoodESDTracksEta09_14(NULL),
  hESDConversionRPhi(NULL),
  hESDConversionRZ(NULL),
  hESDConversionRPt(NULL),
  hESDConversionREta(NULL),
  hESDConversionDCA(NULL),
  hESDConversionPsiPair(NULL),
  hESDConversionChi2(NULL),
  hESDConversionMass(NULL),
  hESDConversionRRejSmall(NULL),
  hESDConversionRRejLarge(NULL),
  hElectronRdEdx(NULL),
  hElectronRNSigmadEdx(NULL),
  hPositronRdEdx(NULL),
  hPositronRNSigmadEdx(NULL),
  hMCConversionRPhi(NULL),
  hMCConversionRPt(NULL),
  hMCConversionREta(NULL),
  hMCConversionRRejSmall(NULL),
  hMCConversionRRejLarge(NULL),
  hMCAllGammaPt(NULL),
  hMCTrueConversionRPhi(NULL),
  hMCTrueConversionRZ(NULL),
  hMCTrueConversionRPt(NULL),
  hMCTrueConversionREta(NULL),
  hMCTrueConversionDCA(NULL),
  hMCTrueConversionPsiPair(NULL),
  hMCTrueConversionChi2(NULL),
  hMCTrueConversionMass(NULL),
  hMCTrueConversionRRejSmall(NULL),
  hMCTrueConversionRRejLarge(NULL),
  hMCTruePi0DalConversionRPt(NULL),
  hMCTruePi0DalConversionEta(NULL),
  hMCTrueEtaDalConversionRPt(NULL),
  hMCTrueEtaDalConversionEta(NULL),
  hMCTrueCombinatorialConversionRPt(NULL),
  hMCTrueCombinatorialConversionEta(NULL)
{
  // Default constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMaterialHistos::~AliAnalysisTaskMaterialHistos()
{
  // default deconstructor
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskMaterialHistos::UserCreateOutputObjects()
{
  // Create User Output Objects

  if(fOutputList != NULL){
    delete fOutputList;
    fOutputList = NULL;
  }
  if(fOutputList == NULL){
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }

  // Array of current cut's gammas

  fGammaCandidates          = new TList();
  fCutFolder                = new TList*[fnCuts];
  fESDList                  = new TList*[fnCuts];
  fMCList                   = new TList*[fnCuts];
  fTrueList                 = new TList*[fnCuts];

  hNEvents                  = new TH1F*[fnCuts];
  hNGoodESDTracksEta09      = new TH1F*[fnCuts];
  hNGoodESDTracksEta14      = new TH1F*[fnCuts];
  hNGoodESDTracksEta09_14   = new TH1F*[fnCuts];
  hESDConversionRPhi        = new TH2F*[fnCuts];
  hESDConversionRZ          = new TH2F*[fnCuts];
  hESDConversionRPt         = new TH2F*[fnCuts];
  hESDConversionREta        = new TH2F*[fnCuts];
  hESDConversionDCA         = new TH1F*[fnCuts];
  hESDConversionPsiPair     = new TH1F*[fnCuts];
  hESDConversionChi2        = new TH1F*[fnCuts];
  hESDConversionMass        = new TH1F*[fnCuts];
  hESDConversionRRejLarge   = new TH1F*[fnCuts];
  hESDConversionRRejSmall   = new TH1F*[fnCuts];

  hElectronRdEdx            = new TH2F*[fnCuts];
  hElectronRNSigmadEdx      = new TH2F*[fnCuts];
  hPositronRdEdx            = new TH2F*[fnCuts];
  hPositronRNSigmadEdx      = new TH2F*[fnCuts];

  hMCConversionRPhi         = new TH2F*[fnCuts];
  hMCConversionRPt          = new TH2F*[fnCuts];
  hMCConversionREta         = new TH2F*[fnCuts];
  hMCConversionRRejLarge    = new TH1F*[fnCuts];
  hMCConversionRRejSmall    = new TH1F*[fnCuts];
  hMCAllGammaPt             = new TH1F*[fnCuts];

  hMCTrueConversionRPhi     = new TH2F*[fnCuts];
  hMCTrueConversionRZ       = new TH2F*[fnCuts];
  hMCTrueConversionRPt      = new TH2F*[fnCuts];
  hMCTrueConversionREta     = new TH2F*[fnCuts];
  hMCTrueConversionDCA      = new TH1F*[fnCuts];
  hMCTrueConversionPsiPair  = new TH1F*[fnCuts];
  hMCTrueConversionChi2     = new TH1F*[fnCuts];
  hMCTrueConversionMass     = new TH1F*[fnCuts];
  hMCTrueConversionRRejLarge = new TH1F*[fnCuts];
  hMCTrueConversionRRejSmall = new TH1F*[fnCuts];

  hMCTruePi0DalConversionRPt = new TH2F*[fnCuts];
  hMCTruePi0DalConversionEta = new TH1F*[fnCuts];
  hMCTrueEtaDalConversionRPt = new TH2F*[fnCuts];
  hMCTrueEtaDalConversionEta = new TH1F*[fnCuts];
  hMCTrueCombinatorialConversionRPt = new TH2F*[fnCuts];
  hMCTrueCombinatorialConversionEta = new TH1F*[fnCuts];


  for(Int_t iCut = 0; iCut<fnCuts;iCut++){


    TString cutstringEvent      = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPhoton     = ((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))->GetCutNumber();
    fCutFolder[iCut]            = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s",cutstringEvent.Data() ,cutstringPhoton.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputList->Add(fCutFolder[iCut]);

    fESDList[iCut]              = new TList();
    fESDList[iCut]->SetName(Form("%s_%s ESD histograms",cutstringEvent.Data() ,cutstringPhoton.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

    Int_t nBinsR=400;
//     Int_t nBinsX=2000;
//     Int_t nBinsY=2000;
    Int_t nBinsZ=750;
    Int_t nBinsPhi=750;
    Int_t nBinsEta=2000;
    Int_t nBinsPt=400;

    hNEvents[iCut]              = new TH1F("NEvents","NEvents",14,-0.5,13.5);
    hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames      = "Not Trigger: ";
        TriggerNames              = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        hNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
        hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(hNEvents[iCut]);


    hNGoodESDTracksEta09[iCut]      = new TH1F("GoodESDTracksEta09","GoodESDTracksEta09",4000,0,4000);
    fESDList[iCut]->Add(hNGoodESDTracksEta09[iCut]);
    hNGoodESDTracksEta14[iCut]      = new TH1F("GoodESDTracksEta14","GoodESDTracksEta14",4000,0,4000);
    fESDList[iCut]->Add(hNGoodESDTracksEta14[iCut]);
    hNGoodESDTracksEta09_14[iCut]   = new TH1F("GoodESDTracksEta09_14","GoodESDTracksEta09_14",4000,0,4000);
    fESDList[iCut]->Add(hNGoodESDTracksEta09_14[iCut]);

    hESDConversionRPhi[iCut]        = new TH2F("ESD_Conversion_RPhi","ESD_Conversion_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRPhi[iCut]);
    hESDConversionREta[iCut]        = new TH2F("ESD_Conversion_REta","ESD_Conversion_REta",nBinsEta,-2.,2.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionREta[iCut]);
    hESDConversionRPt[iCut]         = new TH2F("ESD_Conversion_RPt","ESD_Conversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRPt[iCut]);
    hESDConversionRZ[iCut]          = new TH2F("ESD_Conversion_RZ","ESD_Conversion_RZ",nBinsZ,-180.,180.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRZ[iCut]);

    hElectronRdEdx[iCut]            = new TH2F("Electron_RdEdx","Electron_RdEdx",200,0.,200.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hElectronRdEdx[iCut]);
    hElectronRNSigmadEdx[iCut]      = new TH2F("Electron_RNSigmadEdx","Electron_RNSigmadEdx",200,-10.,10.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hElectronRNSigmadEdx[iCut]);
    hPositronRdEdx[iCut]            = new TH2F("Positron_RdEdx","Positron_RdEdx",200,0.,200.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hPositronRdEdx[iCut]);
    hPositronRNSigmadEdx[iCut]      = new TH2F("Positron_RNSigmadEdx","Positron_RNSigmadEdx",200,-10.,10.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hPositronRNSigmadEdx[iCut]);

    hESDConversionDCA[iCut]         = new TH1F("ESD_Conversion_DCA","ESD_Conversion_DCA",400,0.,5.);
    fESDList[iCut]->Add(hESDConversionDCA[iCut]);
    hESDConversionPsiPair[iCut]     = new TH1F("ESD_Conversion_PsiPair","ESD_Conversion_PsiPair",400,0.,5.);
    fESDList[iCut]->Add(hESDConversionPsiPair[iCut]);
    hESDConversionChi2[iCut]        = new TH1F("ESD_Conversion_Chi2","ESD_Conversion_Chi2",400,0.,50.);
    fESDList[iCut]->Add(hESDConversionChi2[iCut]);
    hESDConversionMass[iCut]        = new TH1F("ESD_Conversion_Mass","ESD_Conversion_Mass",400,0.,1.);
    fESDList[iCut]->Add(hESDConversionMass[iCut]);

    hESDConversionRRejLarge[iCut]   = new TH1F("ESD_Conversion_RLarge","ESD_Conversion_RLarge",nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRRejLarge[iCut]);
    hESDConversionRRejSmall[iCut]   = new TH1F("ESD_Conversion_RSmall","ESD_Conversion_RSmall",nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRRejSmall[iCut]);

    if (fIsMC>0) {

        fMCList[iCut]               = new TList();
        fMCList[iCut]->SetName(Form("%s_%s MC histograms",cutstringEvent.Data() ,cutstringPhoton.Data()));
        fMCList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMCList[iCut]);

        fTrueList[iCut]             = new TList();
        fTrueList[iCut]->SetName(Form("%s_%s True histograms",cutstringEvent.Data() ,cutstringPhoton.Data()));
        fTrueList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fTrueList[iCut]);

        hMCAllGammaPt[iCut]         = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",nBinsPt,0.,20.);
        fMCList[iCut]->Add(hMCAllGammaPt[iCut]);

        hMCConversionRPhi[iCut]     = new TH2F("MC_Conversion_RPhi","MC_Conversion_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionRPhi[iCut]);
        hMCConversionREta[iCut]     = new TH2F("MC_Conversion_REta","MC_Conversion_REta",nBinsEta,-2.,2.,nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionREta[iCut]);
        hMCConversionRPt[iCut]      = new TH2F("MC_Conversion_RPt","MC_Conversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionRPt[iCut]);
        hMCConversionRRejLarge[iCut] = new TH1F("MC_Conversion_RLarge","MC_Conversion_RLarge",nBinsR,0.,200.);
        fESDList[iCut]->Add(hMCConversionRRejLarge[iCut]);
        hMCConversionRRejSmall[iCut] = new TH1F("MC_Conversion_RSmall","MC_Conversion_RSmall",nBinsR,0.,200.);
        fESDList[iCut]->Add(hMCConversionRRejSmall[iCut]);

        hMCTrueConversionRPhi[iCut] = new TH2F("ESD_TrueConversion_RPhi","ESD_TrueConversion_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionRPhi[iCut]);
        hMCTrueConversionREta[iCut] = new TH2F("ESD_TrueConversion_REta","ESD_TrueConversion_REta",nBinsEta,-2.,2.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionREta[iCut]);
        hMCTrueConversionRPt[iCut]  = new TH2F("ESD_TrueConversion_RPt","ESD_TrueConversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionRPt[iCut]);
        hMCTrueConversionRZ[iCut]   = new TH2F("ESD_TrueConversion_RZ","ESD_TrueConversion_RZ",nBinsZ,-180.,180.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionRZ[iCut]);

        hMCTrueConversionDCA[iCut]  = new TH1F("ESD_TrueConversion_DCA","ESD_TrueConversion_DCA",400,0.,5.);
        fTrueList[iCut]->Add(hMCTrueConversionDCA[iCut]);
        hMCTrueConversionPsiPair[iCut] = new TH1F("ESD_TrueConversion_PsiPair","ESD_TrueConversion_PsiPair",400,0.,5.);
        fTrueList[iCut]->Add(hMCTrueConversionPsiPair[iCut]);
        hMCTrueConversionChi2[iCut] = new TH1F("ESD_TrueConversion_Chi2","ESD_TrueConversion_Chi2",400,0.,50.);
        fTrueList[iCut]->Add(hMCTrueConversionChi2[iCut]);
        hMCTrueConversionMass[iCut] = new TH1F("ESD_TrueConversion_Mass","ESD_TrueConversion_Mass",400,0.,1.);
        fTrueList[iCut]->Add(hMCTrueConversionMass[iCut]);

        hMCTrueConversionRRejLarge[iCut] = new TH1F("ESD_TrueConversion_RLarge","ESD_TrueConversion_RLarge",nBinsR,0.,200.);
        fESDList[iCut]->Add(hMCTrueConversionRRejLarge[iCut]);
        hMCTrueConversionRRejSmall[iCut] = new TH1F("ESD_TrueConversion_RSmall","ESD_TrueConversion_RSmall",nBinsR,0.,200.);
        fESDList[iCut]->Add(hMCTrueConversionRRejSmall[iCut]);

        hMCTruePi0DalConversionRPt[iCut]   = new TH2F("ESD_TruePi0DalConversion_RPt","ESD_TruePi0DalConversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTruePi0DalConversionRPt[iCut]);
        hMCTruePi0DalConversionEta[iCut] = new TH1F("ESD_TruePi0DalConversion_Eta","ESD_TruePi0DalConversion_Eta",nBinsEta,-2.,2.);
        fTrueList[iCut]->Add(hMCTruePi0DalConversionEta[iCut]);

        hMCTrueEtaDalConversionRPt[iCut]   = new TH2F("ESD_TrueEtaDalConversion_RPt","ESD_TrueEtaDalConversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueEtaDalConversionRPt[iCut]);
        hMCTrueEtaDalConversionEta[iCut] = new TH1F("ESD_TrueEtaDalConversion_Eta","ESD_TrueEtaDalConversion_Eta",nBinsEta,-2.,2.);
        fTrueList[iCut]->Add(hMCTrueEtaDalConversionEta[iCut]);

        hMCTrueCombinatorialConversionRPt[iCut]     = new TH2F("ESD_TrueCombinatorialConversion_RPt","ESD_TrueCombinatorialConversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueCombinatorialConversionRPt[iCut]);
        hMCTrueCombinatorialConversionEta[iCut]   = new TH1F("ESD_TrueCombinatorialConversion_Eta","ESD_TrueCombinatorialConversion_Eta",nBinsEta,-2.,2.);
        fTrueList[iCut]->Add(hMCTrueCombinatorialConversionEta[iCut]);


    }
  }


  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
    if (fV0Reader->GetV0FindingEfficiencyHistograms())
      fOutputList->Add(fV0Reader->GetV0FindingEfficiencyHistograms());

  if(fV0Reader){
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputList->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputList->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());

  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))) continue;
    if(((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))->GetCutHistograms());
    }
  }

  PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskMaterialHistos::UserExec(Option_t *){

  fInputEvent = InputEvent();
  if (fInputEvent==NULL) return;

  if(fIsMC>0) fMCEvent = MCEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  // Event Not Accepted due to MC event missing or because it is incomplere or  wrong trigger for V0ReaderV1 => skip broken event/file
  if(eventQuality == 2 || eventQuality == 3){
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      hNEvents[iCut]->Fill(eventQuality);
    }
    return;
  }

  fConversionGammas=fV0Reader->GetReconstructedGammas();// Gammas from default Cut

  // ------------------- BeginEvent ----------------------------

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    fiCut = iCut;
    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
    if(eventNotAccepted){
      // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      hNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      continue;
    }

    if(eventQuality != 0){// Event Not Accepted
      // cout << "event rejected due to: " <<eventQuality << endl;
      hNEvents[iCut]->Fill(eventQuality);
      continue;
    }

    hNEvents[iCut]->Fill(eventQuality); // Should be 0 here

    if(fIsMC > 0){
      // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
            ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                            ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                            fMCEvent);
        }
      }
      ProcessMCPhotons();
    }

    fNESDtracksEta09 = CountTracks09(); // Estimate Event Multiplicity
    fNESDtracksEta0914 = CountTracks0914(); // Estimate Event Multiplicity
    fNESDtracksEta14 = fNESDtracksEta09 + fNESDtracksEta0914;

    hNGoodESDTracksEta09[iCut]->Fill(fNESDtracksEta09);
    hNGoodESDTracksEta14[iCut]->Fill(fNESDtracksEta14);
    hNGoodESDTracksEta09_14[iCut]->Fill(fNESDtracksEta0914);


    if(fInputEvent){
      if(fInputEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
        fNContrVtx = fInputEvent->GetPrimaryVertexTracks()->GetNContributors();
      } else {
        fNContrVtx = 0;
      }
    }
    ProcessPhotons();
    fGammaCandidates->Clear(); // delete this cuts good gammas
  }

  //cout<<" done with the event"<<endl;

  PostData(1, fOutputList);
}

///________________________________________________________________________
void AliAnalysisTaskMaterialHistos::FillMCHistograms(Int_t eventPos){
  TParticle* candidate = (TParticle *)fMCEvent->Particle(eventPos);

  if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(candidate,fMCEvent,kFALSE)){

    fGammaMCPt = candidate->Pt();
    fGammaMCTheta = candidate->Theta();

    hMCAllGammaPt[fiCut]->Fill(candidate->Pt());

  }

  if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(candidate,fMCEvent,kTRUE)){

    fGammaMCConvPt = candidate->Pt();
    fGammaMCConvTheta = candidate->Theta();

    TParticle* daughter1 = (TParticle *)fMCEvent->Particle(candidate->GetFirstDaughter());
//     TParticle* daughter2 = (TParticle *)fMCEvent->Particle(candidate->GetLastDaughter());

    hMCConversionRPhi[fiCut]->Fill(candidate->Phi(),daughter1->R());
    hMCConversionREta[fiCut]->Fill(candidate->Eta(),daughter1->R());
    hMCConversionRPt[fiCut]->Fill(candidate->Pt(),daughter1->R());

    if(daughter1->R() < 75. || daughter1->R() > 85.) hMCConversionRRejSmall[fiCut]->Fill(daughter1->R());
    if(daughter1->R() < 70. || daughter1->R() > 90.) hMCConversionRRejLarge[fiCut]->Fill(daughter1->R());

  } // Converted MC Gamma
}

///________________________________________________________________________
void AliAnalysisTaskMaterialHistos::ProcessMCPhotons(){
  // Loop over all primary MC particle
  for(Int_t i = 0; i < fMCEvent->GetNumberOfPrimaries(); i++) {
    TParticle* particle = (TParticle *)fMCEvent->Particle(i);
    if (!particle) continue;


    if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
      Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
      Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
      if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
    }


    if (particle->GetPdgCode() == 111 && particle->GetFirstDaughter() >= fMCEvent->GetNumberOfPrimaries()){
      //cout << "Undecayed pi0 found with mother: " << particle->GetMother(0) << endl;
      for (Int_t j = 0; j < 2 ; j++){
        FillMCHistograms(particle->GetDaughter(j));
      }
    } else {
        FillMCHistograms(i);
    }

  }
}

///________________________________________________________________________
void AliAnalysisTaskMaterialHistos::ProcessPhotons(){

  // Fill Histograms for QA and MC
  TList *GammaCandidatesStepTwo = new TList();

  for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
    AliAODConversionPhoton *gamma= (AliAODConversionPhoton*)fConversionGammas->At(firstGammaIndex);

    if (gamma == NULL) continue;

    if(!((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelected(gamma,fInputEvent))continue;

    if( ! ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->UseToCloseV0sCut()){
      fGammaCandidates->Add(gamma); // if no second loop is required add to events good gammas
    }else if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->UseToCloseV0sCut()) { // shared electron is disabled, step one not needed -> step two
      GammaCandidatesStepTwo->Add(gamma);
    }
  }

  if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
      if(!PhotonCandidate) continue;
      if(!((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
      fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
    }
  }

  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
    AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
    if (gamma==NULL) continue;

    fGammaPt        = gamma->GetPhotonPt();
    fGammaTheta     = gamma->GetPhotonTheta();
    fGammaChi2NDF   = gamma->GetChi2perNDF();

    AliPIDResponse* pidResonse = ((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();

    AliVTrack * negTrack = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetTrack(fInputEvent, gamma->GetTrackLabelNegative());
    AliVTrack * posTrack = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetTrack(fInputEvent, gamma->GetTrackLabelPositive());

    hESDConversionRPhi[fiCut]->Fill(gamma->GetPhotonPhi(),gamma->GetConversionRadius());
    hESDConversionRZ[fiCut]->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());
    hESDConversionREta[fiCut]->Fill(gamma->GetPhotonEta(),gamma->GetConversionRadius());
    hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius());

    if(negTrack->GetTPCsignal()){
        hElectronRdEdx[fiCut]->Fill(negTrack->GetTPCsignal(),gamma->GetConversionRadius());
        hElectronRNSigmadEdx[fiCut]->Fill(pidResonse->NumberOfSigmasTPC(negTrack, AliPID::kElectron),gamma->GetConversionRadius());
    }
    if(posTrack->GetTPCsignal()){
        hPositronRdEdx[fiCut]->Fill(posTrack->GetTPCsignal(),gamma->GetConversionRadius());
        hPositronRNSigmadEdx[fiCut]->Fill(pidResonse->NumberOfSigmasTPC(posTrack, AliPID::kElectron),gamma->GetConversionRadius());
    }

    if(gamma->GetConversionRadius() < 75. || gamma->GetConversionRadius() > 85.) hESDConversionRRejSmall[fiCut]->Fill(gamma->GetConversionRadius());
    if(gamma->GetConversionRadius() < 70. || gamma->GetConversionRadius() > 90.) hESDConversionRRejLarge[fiCut]->Fill(gamma->GetConversionRadius());

    if(fInputEvent->IsA()==AliESDEvent::Class()){
      AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
      if(esdEvent){
        AliESDv0 *v0 = esdEvent->GetV0(gamma->GetV0Index());
        hESDConversionDCA[fiCut]->Fill(v0->GetDcaV0Daughters());
      }
    }
    hESDConversionPsiPair[fiCut]->Fill(gamma->GetPsiPair());
    hESDConversionChi2[fiCut]->Fill(gamma->GetChi2perNDF());
    hESDConversionMass[fiCut]->Fill(gamma->GetInvMassPair());

    fKind = 9;
    Int_t pdgCodePos = 0.;
    Int_t pdgCodeNeg = 0.;

    if(fIsMC>0){

      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();

      TParticle *posDaughter = gamma->GetPositiveMCDaughter(fMCEvent);
      TParticle *negDaughter = gamma->GetNegativeMCDaughter(fMCEvent);
      //cout << "generate Daughters: "<<posDaughter << "\t" << negDaughter << endl;

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
            Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCEvent, fInputEvent);
            Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCEvent, fInputEvent);
            if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
      }

      if(posDaughter == NULL || negDaughter == NULL){

        fKind = 9; // garbage

      } else if(posDaughter->GetMother(0) != negDaughter->GetMother(0) || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)){

        fKind = 1; //Not Same Mother == Combinatorial Bck
        pdgCodePos = posDaughter->GetPdgCode();
        pdgCodeNeg = negDaughter->GetPdgCode();

        if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11)
            fKind = 10; //Electron Combinatorial
        if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11 && (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1))
            fKind = 15; //direct Electron Combinatorial
        if(TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==211)
            fKind = 11; //Pion Combinatorial
        if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==2212) ||
            (TMath::Abs(pdgCodePos)==2212 && TMath::Abs(pdgCodeNeg)==211))
            fKind = 12; //Pion, Proton Combinatorics
        if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==11) ||
            (TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==211))
            fKind = 13; //Pion, Electron Combinatorics
        if (TMath::Abs(pdgCodePos)==321 || TMath::Abs(pdgCodeNeg)==321)
            fKind = 14; //Kaon combinatorics

      } else {
        //cout << "same mother" << endl;
	pdgCodePos = posDaughter->GetPdgCode();
        pdgCodeNeg = negDaughter->GetPdgCode();
 	Int_t pdgCode;
        pdgCode = gamma->GetMCParticle(fMCEvent)->GetPdgCode();
        if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11)
            fKind = 2; // combinatorics from hadronic decays
        else if ( !(pdgCodeNeg==pdgCodePos)){
            Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
            if(pdgCode == 111)      fKind = 3; // pi0 Dalitz
            else if (pdgCode == 221)fKind = 4; // eta Dalitz
            else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
                if(pdgCode == 22 && gammaIsPrimary){
                    fKind = 0; // primary photons
                } else if (pdgCode == 22){
                    fKind = 5; //secondary photons
                }
            } else 	fKind = 9; //garbage
        } else fKind = 9; //garbage
      }

      if(fKind==0 || fKind==5){
        hMCTrueConversionRPhi[fiCut]->Fill(gamma->GetPhotonPhi(),gamma->GetConversionRadius());
        hMCTrueConversionRZ[fiCut]->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());
        hMCTrueConversionREta[fiCut]->Fill(gamma->GetPhotonEta(),gamma->GetConversionRadius());
        hMCTrueConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius());

        if(gamma->GetConversionRadius() < 75. || gamma->GetConversionRadius() > 85.) hMCTrueConversionRRejSmall[fiCut]->Fill(gamma->GetConversionRadius());
        if(gamma->GetConversionRadius() < 70. || gamma->GetConversionRadius() > 90.) hMCTrueConversionRRejLarge[fiCut]->Fill(gamma->GetConversionRadius());

        hMCTrueConversionPsiPair[fiCut]->Fill(gamma->GetPsiPair());
        hMCTrueConversionChi2[fiCut]->Fill(gamma->GetChi2perNDF());
        hMCTrueConversionMass[fiCut]->Fill(gamma->GetInvMassPair());
        if(fInputEvent->IsA()==AliESDEvent::Class()){
            AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
            if(esdEvent){
                AliESDv0 *v0 = esdEvent->GetV0(gamma->GetV0Index());
                hMCTrueConversionDCA[fiCut]->Fill(v0->GetDcaV0Daughters());
            }
        }

      } else if(fKind==3){
        hMCTruePi0DalConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius());
        hMCTruePi0DalConversionEta[fiCut]->Fill(gamma->GetPhotonEta());
      } else if(fKind==4){
        hMCTrueEtaDalConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius());
        hMCTrueEtaDalConversionEta[fiCut]->Fill(gamma->GetPhotonEta());
      } else {
        hMCTrueCombinatorialConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius());
        hMCTrueCombinatorialConversionEta[fiCut]->Fill(gamma->GetPhotonEta());
      }
    }
  }

  delete GammaCandidatesStepTwo;
  GammaCandidatesStepTwo = 0x0;


}

//________________________________________________________________________
Int_t AliAnalysisTaskMaterialHistos::CountTracks09(){

  Int_t fNumberOfESDTracks = 0;
  if(fInputEvent->IsA()==AliESDEvent::Class()){
  // Using standard function for setting Cuts

//     Bool_t selectPrimaries = kTRUE;
    static AliESDtrackCuts *EsdTrackCuts = 0x0;
    static int prevRun = -1;
    // Using standard function for setting Cuts
    Int_t runNumber = fInputEvent->GetRunNumber();
    if (prevRun!=runNumber) {
      delete EsdTrackCuts;
      EsdTrackCuts = 0;
      prevRun = runNumber;
    }
//     AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
    if (!EsdTrackCuts) {
        // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
        if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

        } else if (runNumber>=209122){ // else if run2 data use 2015 PbPb cuts
            //EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
            // hard coded track cuts for the moment, because AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb() gives spams warnings
            EsdTrackCuts = new AliESDtrackCuts();
            // TPC; clusterCut = 1, cutAcceptanceEdges = kTRUE, removeDistortedRegions = kFALSE
            EsdTrackCuts->AliESDtrackCuts::SetMinNCrossedRowsTPC(70);
            EsdTrackCuts->AliESDtrackCuts::SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
            EsdTrackCuts->SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0);  // only dead zone and not clusters per length
            //EsdTrackCuts->AliESDtrackCuts::SetCutOutDistortedRegionsTPC(kTRUE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterTPC(4);
            EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
            // ITS; selPrimaries = 1
            EsdTrackCuts->AliESDtrackCuts::SetRequireITSRefit(kTRUE);
            EsdTrackCuts->AliESDtrackCuts::SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                        AliESDtrackCuts::kAny);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
            EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);

        } else { // else use 2011 version of track cuts
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
        }
        EsdTrackCuts->SetMaxDCAToVertexZ(2);
        EsdTrackCuts->SetEtaRange(-0.9, 0.9);
        EsdTrackCuts->SetPtRange(0.15);
    }

    for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ){
        if (fMCEvent){
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
                        Int_t isFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()), fMCEvent, fInputEvent);
            if( (isFromMBHeader < 1) ) continue;
          }
        }
        fNumberOfESDTracks++;
      }
    }
    delete EsdTrackCuts;
    EsdTrackCuts=0x0;

  }

  return fNumberOfESDTracks;

}

//________________________________________________________________________
Int_t AliAnalysisTaskMaterialHistos::CountTracks0914(){

  Int_t fNumberOfESDTracks = 0;
  if(fInputEvent->IsA()==AliESDEvent::Class()){
    // Using standard function for setting Cuts

//     Bool_t selectPrimaries = kTRUE;
    static AliESDtrackCuts *EsdTrackCuts = 0x0;
    static int prevRun = -1;
    // Using standard function for setting Cuts
    Int_t runNumber = fInputEvent->GetRunNumber();
    if (prevRun!=runNumber) {
      delete EsdTrackCuts;
      EsdTrackCuts = 0;
      prevRun = runNumber;
    }
//     AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
    if (!EsdTrackCuts) {
        // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
        if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

        } else if (runNumber>=209122){ // else if run2 data use 2015 PbPb cuts
            //EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
            // hard coded track cuts for the moment, because AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb() gives spams warnings
            EsdTrackCuts = new AliESDtrackCuts();
            // TPC; clusterCut = 1, cutAcceptanceEdges = kTRUE, removeDistortedRegions = kFALSE
            EsdTrackCuts->AliESDtrackCuts::SetMinNCrossedRowsTPC(70);
            EsdTrackCuts->AliESDtrackCuts::SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
            EsdTrackCuts->SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0);  // only dead zone and not clusters per length
            //EsdTrackCuts->AliESDtrackCuts::SetCutOutDistortedRegionsTPC(kTRUE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterTPC(4);
            EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
            // ITS; selPrimaries = 1
            EsdTrackCuts->AliESDtrackCuts::SetRequireITSRefit(kTRUE);
            EsdTrackCuts->AliESDtrackCuts::SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                        AliESDtrackCuts::kAny);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
            EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);

        } else { // else use 2011 version of track cuts
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
        }
        EsdTrackCuts->SetMaxDCAToVertexZ(2);
        EsdTrackCuts->SetPtRange(0.15);
    }

    EsdTrackCuts->SetEtaRange(0.9, 1.4);
    for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ){
        if (fMCEvent){
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
                        Int_t isFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()), fMCEvent, fInputEvent);
            if( (isFromMBHeader < 1) ) continue;
          }
        }
        fNumberOfESDTracks++;
      }
    }

    EsdTrackCuts->SetEtaRange(-1.4, -0.9);
    for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack =(AliESDtrack*) fInputEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ){
        if (fMCEvent){
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
                        Int_t isFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()), fMCEvent, fInputEvent);
            if( (isFromMBHeader < 1) ) continue;
          }
        }
        fNumberOfESDTracks++;
      }
    }
    delete EsdTrackCuts;
    EsdTrackCuts=0x0;

  }

  return fNumberOfESDTracks;
}


//________________________________________________________________________
void AliAnalysisTaskMaterialHistos::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskMaterialHistos::Terminate(Option_t *)
{
//    if (fStreamMaterial){
//       fStreamMaterial->GetFile()->Write();
//    }
}
