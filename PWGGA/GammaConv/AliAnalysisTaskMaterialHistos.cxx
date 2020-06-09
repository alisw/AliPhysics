/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Authors: Ana Marin                                                      *
* Version 1.0                                                             *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
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
#include "TMath.h"

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
  fDeDxMapList(NULL),
  fOutputList(NULL),
  fAllMCGammaList(NULL),
  fAllMCConvGammaList(NULL),
  fPrimVtxZ(0.),
  fNContrVtx(0),
  fNESDtracksEta08(0),
  fNESDtracksEta08pt200(0),
  fNESDtracksEta08pt300(0),
  fNESDtracksEta08pt400(0),
  fNESDtracksEta0814(0),
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
  fDoDeDxMaps(0),
  fDoMultWeights(0),
  fWeightMultMC(1),
  hNEvents(NULL),
  hBCNumber(NULL),
  hBCNumberSelected(NULL),
  hNGoodESDTracksEta08(NULL),
  hNGoodESDTracksWeightedEta08(NULL),
  hNGoodESDTracksEta08pt200(NULL),
  hNGoodESDTracksWeightedEta08pt200(NULL),
  hNGoodESDTracksEta08pt300(NULL),
  hNGoodESDTracksWeightedEta08pt300(NULL),
  hNGoodESDTracksEta08pt400(NULL),
  hNGoodESDTracksWeightedEta08pt400(NULL),
  hNGoodESDTracksEta14(NULL),
  hNGoodESDTracksEta08_14(NULL),
  fHistoNV0Tracks(NULL),
  fHistoNV0TracksWeighted(NULL),
  hESDConversionRPhi(NULL),
  hESDConversionRPhiFromConv(NULL),
  hESDConversionRZ(NULL),
  hESDConversionRPt(NULL),
  hESDConversionWOWeightRPt(NULL),
  hESDConversionREta(NULL),
  hESDConversionDCA(NULL),
  hESDConversionPsiPair(NULL),
  hESDConversionChi2(NULL),
  hESDConversionMass(NULL),
  hESDConversionRRejSmall(NULL),
  hESDConversionRRejLarge(NULL),
  hESDConversionAsymP(NULL),
  hElectronRdEdx(NULL),
  hElectronRNSigmadEdx(NULL),
  hPositronRdEdx(NULL),
  hPositronRNSigmadEdx(NULL),
  hMCConversionRPhi(NULL),
  hMCConversionRPhiFromConv(NULL),
  hMCConversionRPt(NULL),
  hMCConversionWOWeightRPt(NULL),
  hMCConversionREta(NULL),
  hMCConversionRRejSmall(NULL),
  hMCConversionRRejLarge(NULL),
  hMCAllGammaPt(NULL),
  hMCAllGammaWOWeightPt(NULL),
  hMCAllSecondaryGammaPt(NULL),
  hMCSecondaryConvGammaPtR(NULL),
  hMCTrueConversionRPhi(NULL),
  hMCTrueConversionRPhiFromConv(NULL),
  hMCTrueConversionRZ(NULL),
  hMCTrueConversionRPt(NULL),
  hMCTrueConversionWOWeightRPt(NULL),
  hMCTrueConversionRPtMCRPt(NULL),
  hMCTrueConversionWOWeightRPtMCRPt(NULL),
  hMCTrueConversionREta(NULL),
  hMCTrueConversionDCA(NULL),
  hMCTrueConversionPsiPair(NULL),
  hMCTrueConversionChi2(NULL),
  hMCTrueConversionMass(NULL),
  hMCTrueConversionAsymP(NULL),
  hMCTrueConversionRRejSmall(NULL),
  hMCTrueConversionRRejLarge(NULL),
  hMCTruePrimConversionRPt(NULL),
  hMCTruePrimConversionWOWeightRPt(NULL),
  hMCTrueSecConversionRPt(NULL),
  hMCTrueSecondaryConvGammaRPt(NULL),
  hMCTrueSecondaryConvGammaMCRPt(NULL),
  hMCTruePi0DalConversionRPt(NULL),
  hMCTruePi0DalConversionEta(NULL),
  hMCTrueEtaDalConversionRPt(NULL),
  hMCTrueEtaDalConversionEta(NULL),
  hMCTrueCombinatorialConversionRPt(NULL),
  hMCTrueCombinatorialConversionEta(NULL),
  hPositrondEdxMapsR0(NULL),
  hElectrondEdxMapsR0(NULL),
  hPositrondEdxMapsR1(NULL),
  hElectrondEdxMapsR1(NULL),
  hPositrondEdxMapsR2(NULL),
  hElectrondEdxMapsR2(NULL),
  hPositrondEdxMapsR3(NULL),
  hElectrondEdxMapsR3(NULL),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  fDoSelectBCNumber(kFALSE),
  fBCNumber(0),
  fRunNumber(0)
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
  fDeDxMapList(NULL),
  fOutputList(NULL),
  fAllMCGammaList(NULL),
  fAllMCConvGammaList(NULL),
  fPrimVtxZ(0.),
  fNContrVtx(0),
  fNESDtracksEta08(0),
  fNESDtracksEta08pt200(0),
  fNESDtracksEta08pt300(0),
  fNESDtracksEta08pt400(0),
  fNESDtracksEta0814(0),
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
  fDoDeDxMaps(0),
  fDoMultWeights(0),
  fWeightMultMC(1),
  hNEvents(NULL),
  hBCNumber(NULL),
  hBCNumberSelected(NULL),
  hNGoodESDTracksEta08(NULL),
  hNGoodESDTracksWeightedEta08(NULL),
  hNGoodESDTracksEta08pt200(NULL),
  hNGoodESDTracksWeightedEta08pt200(NULL),
  hNGoodESDTracksEta08pt300(NULL),
  hNGoodESDTracksWeightedEta08pt300(NULL),
  hNGoodESDTracksEta08pt400(NULL),
  hNGoodESDTracksWeightedEta08pt400(NULL),
  hNGoodESDTracksEta14(NULL),
  hNGoodESDTracksEta08_14(NULL),
  fHistoNV0Tracks(NULL),
  fHistoNV0TracksWeighted(NULL),
  hESDConversionRPhi(NULL),
  hESDConversionRPhiFromConv(NULL),
  hESDConversionRZ(NULL),
  hESDConversionRPt(NULL),
  hESDConversionWOWeightRPt(NULL),
  hESDConversionREta(NULL),
  hESDConversionDCA(NULL),
  hESDConversionPsiPair(NULL),
  hESDConversionChi2(NULL),
  hESDConversionMass(NULL),
  hESDConversionRRejSmall(NULL),
  hESDConversionRRejLarge(NULL),
  hESDConversionAsymP(NULL),
  hElectronRdEdx(NULL),
  hElectronRNSigmadEdx(NULL),
  hPositronRdEdx(NULL),
  hPositronRNSigmadEdx(NULL),
  hMCConversionRPhi(NULL),
  hMCConversionRPhiFromConv(NULL),
  hMCConversionRPt(NULL),
  hMCConversionWOWeightRPt(NULL),
  hMCConversionREta(NULL),
  hMCConversionRRejSmall(NULL),
  hMCConversionRRejLarge(NULL),
  hMCAllGammaPt(NULL),
  hMCAllGammaWOWeightPt(NULL),
  hMCAllSecondaryGammaPt(NULL),
  hMCSecondaryConvGammaPtR(NULL),
  hMCTrueConversionRPhi(NULL),
  hMCTrueConversionRPhiFromConv(NULL),
  hMCTrueConversionRZ(NULL),
  hMCTrueConversionRPt(NULL),
  hMCTrueConversionWOWeightRPt(NULL),
  hMCTrueConversionRPtMCRPt(NULL),
  hMCTrueConversionWOWeightRPtMCRPt(NULL),
  hMCTrueConversionREta(NULL),
  hMCTrueConversionDCA(NULL),
  hMCTrueConversionPsiPair(NULL),
  hMCTrueConversionChi2(NULL),
  hMCTrueConversionMass(NULL),
  hMCTrueConversionAsymP(NULL),
  hMCTrueConversionRRejSmall(NULL),
  hMCTrueConversionRRejLarge(NULL),
  hMCTruePrimConversionRPt(NULL),
  hMCTruePrimConversionWOWeightRPt(NULL),
  hMCTrueSecConversionRPt(NULL),
  hMCTrueSecondaryConvGammaRPt(NULL),
  hMCTrueSecondaryConvGammaMCRPt(NULL),
  hMCTruePi0DalConversionRPt(NULL),
  hMCTruePi0DalConversionEta(NULL),
  hMCTrueEtaDalConversionRPt(NULL),
  hMCTrueEtaDalConversionEta(NULL),
  hMCTrueCombinatorialConversionRPt(NULL),
  hMCTrueCombinatorialConversionEta(NULL),
  hPositrondEdxMapsR0(NULL),
  hElectrondEdxMapsR0(NULL),
  hPositrondEdxMapsR1(NULL),
  hElectrondEdxMapsR1(NULL),
  hPositrondEdxMapsR2(NULL),
  hElectrondEdxMapsR2(NULL),
  hPositrondEdxMapsR3(NULL),
  hElectrondEdxMapsR3(NULL),
  fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE),
  fDoSelectBCNumber(kFALSE),
  fBCNumber(0),
  fRunNumber(0)
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

  Int_t nTracks                 = 200;
  if(fIsHeavyIon == 1){
    nTracks                     = 4000;
  } else if(fIsHeavyIon == 2){
    nTracks                     = 400;
  }
  // Array of current cut's gammas

  fGammaCandidates          = new TList();
  fCutFolder                = new TList*[fnCuts];
  fESDList                  = new TList*[fnCuts];
  fMCList                   = new TList*[fnCuts];
  fTrueList                 = new TList*[fnCuts];
  if (fDoDeDxMaps>0) {
    fDeDxMapList            = new TList*[fnCuts];
  }
  hNEvents                  = new TH1F*[fnCuts];
  hBCNumber                 = new TH1F*[fnCuts];
  if(fDoSelectBCNumber){
  hBCNumberSelected         = new TH1F*[fnCuts];
  }
  hNGoodESDTracksEta08      = new TH1F*[fnCuts];
  hNGoodESDTracksEta08pt200 = new TH1F*[fnCuts];
  hNGoodESDTracksEta08pt300 = new TH1F*[fnCuts];
  hNGoodESDTracksEta08pt400 = new TH1F*[fnCuts];

  if (fDoMultWeights>0 && fIsMC>0 ) {
    hNGoodESDTracksWeightedEta08  = new TH1F*[fnCuts];
    hNGoodESDTracksWeightedEta08pt200  = new TH1F*[fnCuts];
    hNGoodESDTracksWeightedEta08pt300  = new TH1F*[fnCuts];
    hNGoodESDTracksWeightedEta08pt400  = new TH1F*[fnCuts];
  }
  hNGoodESDTracksEta14      = new TH1F*[fnCuts];
  hNGoodESDTracksEta08_14   = new TH1F*[fnCuts];
  fHistoNV0Tracks           = new TH1F*[fnCuts];
  fHistoNV0TracksWeighted   = new TH1F*[fnCuts];
  hESDConversionRPhi        = new TH2F*[fnCuts];
  hESDConversionRPhiFromConv= new TH2F*[fnCuts];
  hESDConversionRZ          = new TH2F*[fnCuts];
  hESDConversionRPt         = new TH2F*[fnCuts];
  hESDConversionWOWeightRPt = new TH2F*[fnCuts];
  hESDConversionREta        = new TH2F*[fnCuts];
  hESDConversionDCA         = new TH1F*[fnCuts];
  hESDConversionPsiPair     = new TH1F*[fnCuts];
  hESDConversionChi2        = new TH1F*[fnCuts];
  hESDConversionMass        = new TH1F*[fnCuts];
  hESDConversionRRejLarge   = new TH1F*[fnCuts];
  hESDConversionRRejSmall   = new TH1F*[fnCuts];
  hESDConversionAsymP       = new TH2F*[fnCuts];

  hElectronRdEdx            = new TH2F*[fnCuts];
  hElectronRNSigmadEdx      = new TH2F*[fnCuts];
  hPositronRdEdx            = new TH2F*[fnCuts];
  hPositronRNSigmadEdx      = new TH2F*[fnCuts];

  if (fDoDeDxMaps>0) {
    hElectrondEdxMapsR0  =   new TH3F*[fnCuts];
    hPositrondEdxMapsR0  =   new TH3F*[fnCuts];

    hElectrondEdxMapsR1  =   new TH3F*[fnCuts];
    hPositrondEdxMapsR1  =   new TH3F*[fnCuts];

    hElectrondEdxMapsR2  =   new TH3F*[fnCuts];
    hPositrondEdxMapsR2  =   new TH3F*[fnCuts];

    hElectrondEdxMapsR3  =   new TH3F*[fnCuts];
    hPositrondEdxMapsR3  =   new TH3F*[fnCuts];
  }

  hMCConversionRPhi         = new TH2F*[fnCuts];
  hMCConversionRPhiFromConv = new TH2F*[fnCuts];
  hMCConversionRPt          = new TH2F*[fnCuts];
  hMCConversionWOWeightRPt  = new TH2F*[fnCuts];
  hMCConversionREta         = new TH2F*[fnCuts];
  hMCConversionRRejLarge    = new TH1F*[fnCuts];
  hMCConversionRRejSmall    = new TH1F*[fnCuts];
  hMCAllGammaPt             = new TH1F*[fnCuts];
  hMCAllGammaWOWeightPt     = new TH1F*[fnCuts];
  hMCAllSecondaryGammaPt    = new TH2F*[fnCuts];
  hMCSecondaryConvGammaPtR  = new TH3F*[fnCuts];

  hMCTrueConversionRPhi        = new TH2F*[fnCuts];
  hMCTrueConversionRPhiFromConv= new TH2F*[fnCuts];
  hMCTrueConversionRZ          = new TH2F*[fnCuts];
  hMCTrueConversionRPt         = new TH2F*[fnCuts];
  hMCTrueConversionWOWeightRPt = new TH2F*[fnCuts];
  hMCTrueConversionRPtMCRPt         = new TH2F*[fnCuts];
  hMCTrueConversionWOWeightRPtMCRPt = new TH2F*[fnCuts];
  hMCTrueConversionREta     = new TH2F*[fnCuts];
  hMCTrueConversionDCA      = new TH1F*[fnCuts];
  hMCTrueConversionPsiPair  = new TH1F*[fnCuts];
  hMCTrueConversionChi2     = new TH1F*[fnCuts];
  hMCTrueConversionMass     = new TH1F*[fnCuts];
  hMCTrueConversionAsymP    = new TH2F*[fnCuts];
  hMCTrueConversionRRejLarge = new TH1F*[fnCuts];
  hMCTrueConversionRRejSmall = new TH1F*[fnCuts];

  hMCTruePrimConversionRPt         = new TH2F*[fnCuts];
  hMCTruePrimConversionWOWeightRPt = new TH2F*[fnCuts];
  hMCTrueSecConversionRPt    = new TH2F*[fnCuts];

  hMCTrueSecondaryConvGammaRPt    = new TH3F*[fnCuts];
  hMCTrueSecondaryConvGammaMCRPt  = new TH3F*[fnCuts];


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
    hNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(hNEvents[iCut]);

    hBCNumber[iCut]                      = new TH1F("BCNumber","BCNumber",3564,-0.5,3563.5);
    fESDList[iCut]->Add(hBCNumber[iCut]);
    if(fDoSelectBCNumber){
    hBCNumberSelected[iCut]              = new TH1F("BCNumberSel","BCNumberSel",3564,-0.5,3563.5);
    fESDList[iCut]->Add(hBCNumberSelected[iCut]);
    }
    hNGoodESDTracksEta08[iCut]      = new TH1F("GoodESDTracksEta08","GoodESDTracksEta08",nTracks, -0.5, nTracks-0.5);
    fESDList[iCut]->Add(hNGoodESDTracksEta08[iCut]);

    hNGoodESDTracksEta08pt200[iCut]      = new TH1F("GoodESDTracksEta08pt200","GoodESDTracksEta08pt200",nTracks, -0.5, nTracks-0.5);
    fESDList[iCut]->Add(hNGoodESDTracksEta08pt200[iCut]);

    hNGoodESDTracksEta08pt300[iCut]      = new TH1F("GoodESDTracksEta08pt300","GoodESDTracksEta08pt300",nTracks, -0.5, nTracks-0.5);
    fESDList[iCut]->Add(hNGoodESDTracksEta08pt300[iCut]);

    hNGoodESDTracksEta08pt400[iCut]      = new TH1F("GoodESDTracksEta08pt400","GoodESDTracksEta08pt400",nTracks, -0.5, nTracks-0.5);
    fESDList[iCut]->Add(hNGoodESDTracksEta08pt400[iCut]);




    if(fDoMultWeights && fIsMC>0) {
      hNGoodESDTracksWeightedEta08[iCut]      = new TH1F("GoodESDTracksWeightedEta08","GoodESDTracksWeigthedEta08",nTracks, -0.5, nTracks-0.5);
      hNGoodESDTracksWeightedEta08[iCut]->Sumw2();
      fESDList[iCut]->Add(hNGoodESDTracksWeightedEta08[iCut]);

      hNGoodESDTracksWeightedEta08pt200[iCut]      = new TH1F("GoodESDTracksWeightedEta08pt200","GoodESDTracksWeigthedEta08pt200",nTracks, -0.5, nTracks-0.5);
      hNGoodESDTracksWeightedEta08pt200[iCut]->Sumw2();
      fESDList[iCut]->Add(hNGoodESDTracksWeightedEta08pt200[iCut]);


      hNGoodESDTracksWeightedEta08pt300[iCut]      = new TH1F("GoodESDTracksWeightedEta08pt300","GoodESDTracksWeigthedEta08pt300",nTracks, -0.5, nTracks-0.5);
      hNGoodESDTracksWeightedEta08pt300[iCut]->Sumw2();
      fESDList[iCut]->Add(hNGoodESDTracksWeightedEta08pt300[iCut]);


      hNGoodESDTracksWeightedEta08pt400[iCut]      = new TH1F("GoodESDTracksWeightedEta08pt400","GoodESDTracksWeigthedEta08pt300",nTracks, -0.5, nTracks-0.5);
      hNGoodESDTracksWeightedEta08pt400[iCut]->Sumw2();
      fESDList[iCut]->Add(hNGoodESDTracksWeightedEta08pt400[iCut]);


    }

    hNGoodESDTracksEta14[iCut]      = new TH1F("GoodESDTracksEta14","GoodESDTracksEta14",nTracks, -0.5, nTracks-0.5);
    fESDList[iCut]->Add(hNGoodESDTracksEta14[iCut]);
    hNGoodESDTracksEta08_14[iCut]   = new TH1F("GoodESDTracksEta08_14","GoodESDTracksEta08_14",nTracks, -0.5, nTracks-0.5);
    fESDList[iCut]->Add(hNGoodESDTracksEta08_14[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNV0Tracks[iCut]            = new TH1F("V0 Multiplicity", "V0 Multiplicity", 20000, 0, 40000);
    else if(fIsHeavyIon == 2)
      fHistoNV0Tracks[iCut]            = new TH1F("V0 Multiplicity", "V0 Multiplicity", 2500, 0, 2500);
    else
      fHistoNV0Tracks[iCut]            = new TH1F("V0 Multiplicity", "V0 Multiplicity", 1500, 0, 1500);

    fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);

    if(fDoMultWeights && fIsMC>0) {
      fHistoNV0TracksWeighted[iCut]            = new TH1F("V0 Multiplicity Weighted", "V0 Multiplicity Weighted", 1500, 0, 1500);
      fHistoNV0TracksWeighted[iCut]->Sumw2();
      fESDList[iCut]->Add(fHistoNV0TracksWeighted[iCut]);
    }

    hESDConversionRPhi[iCut]        = new TH2F("ESD_Conversion_RPhi","ESD_Conversion_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRPhi[iCut]);

    hESDConversionRPhiFromConv[iCut]        = new TH2F("ESD_Conversion_RPhi_FromConv","ESD_Conversion_RPhi_FromConv",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRPhiFromConv[iCut]);


    hESDConversionREta[iCut]        = new TH2F("ESD_Conversion_REta","ESD_Conversion_REta",nBinsEta,-2.,2.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionREta[iCut]);

    hESDConversionRPt[iCut]         = new TH2F("ESD_Conversion_RPt","ESD_Conversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionRPt[iCut]);
    hESDConversionRPt[iCut]->Sumw2();

    hESDConversionWOWeightRPt[iCut]         = new TH2F("ESD_ConversionWOWeight_RPt","ESD_ConversionWOWeight_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
    fESDList[iCut]->Add(hESDConversionWOWeightRPt[iCut]);

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

    hESDConversionAsymP[iCut]       = new TH2F("ESD_ConversionMapping_AsymP","ESD_ConversionMapping_AsymP",nBinsPt,0.01,20.,500,0.,1.);
    fESDList[iCut]->Add(hESDConversionAsymP[iCut]);
    TAxis *AxisAfter = hESDConversionAsymP[iCut]->GetXaxis();
    Int_t bins = AxisAfter->GetNbins();
    Double_t from = AxisAfter->GetXmin();
    Double_t to = AxisAfter->GetXmax();
    Double_t *newBins = new Double_t[bins+1];
    newBins[0] = from;
    Double_t factor = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
    AxisAfter->Set(bins, newBins);


    if ((fDoMultWeights>0 && fIsMC>0) || ( fDoMaterialBudgetWeightingOfGammasForTrueMesons>0 && fIsMC>0) ) {
      hESDConversionRPt[iCut] ->Sumw2();
      hESDConversionDCA[iCut] ->Sumw2();
      hESDConversionChi2[iCut] ->Sumw2();
      hESDConversionPsiPair[iCut] ->Sumw2();
      hESDConversionMass[iCut] ->Sumw2();
      hESDConversionAsymP[iCut] ->Sumw2();
    }

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
        hMCAllGammaWOWeightPt[iCut] = new TH1F("MC_AllGammaWOWeight_Pt","MC_AllGammaWOWeight_Pt",nBinsPt,0.,20.);
        fMCList[iCut]->Add(hMCAllGammaWOWeightPt[iCut]);

	      hMCAllSecondaryGammaPt[iCut]    = new TH2F("MC_AllSecondaryGamma_Pt", "MC_AllSecondaryGamma_Pt", nBinsPt, 0., 20., 4, -0.5, 3.5);
	      hMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
      	hMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
	      hMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
	      hMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"rest");
	      fMCList[iCut]->Add(hMCAllSecondaryGammaPt[iCut]);

      	hMCSecondaryConvGammaPtR[iCut]  = new TH3F("MC_SecondaryConvGamma_PtR", "MC_SecondaryConvGamma_PtR", nBinsPt, 0., 20., nBinsR,0.,200.,4, -0.5, 3.5);
        hMCSecondaryConvGammaPtR[iCut]->GetZaxis()->SetBinLabel(1,"K0s");
        hMCSecondaryConvGammaPtR[iCut]->GetZaxis()->SetBinLabel(2,"K0l");
        hMCSecondaryConvGammaPtR[iCut]->GetZaxis()->SetBinLabel(3,"Lambda");
        hMCSecondaryConvGammaPtR[iCut]->GetZaxis()->SetBinLabel(4,"rest");
        fMCList[iCut]->Add(hMCSecondaryConvGammaPtR[iCut]);


        hMCConversionRPhi[iCut]     = new TH2F("MC_Conversion_RPhi","MC_Conversion_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionRPhi[iCut]);
        hMCConversionRPhiFromConv[iCut]     = new TH2F("MC_Conversion_RPhi_FromConv","MC_Conversion_RPhi_FromConv",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionRPhiFromConv[iCut]);

        hMCConversionREta[iCut]     = new TH2F("MC_Conversion_REta","MC_Conversion_REta",nBinsEta,-2.,2.,nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionREta[iCut]);
        hMCConversionRPt[iCut]      = new TH2F("MC_Conversion_RPt","MC_Conversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionRPt[iCut]);
        hMCConversionWOWeightRPt[iCut] = new TH2F("MC_ConversionWOWeight_RPt","MC_ConversionWOWeight_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fMCList[iCut]->Add(hMCConversionWOWeightRPt[iCut]);

        hMCConversionRRejLarge[iCut] = new TH1F("MC_Conversion_RLarge","MC_Conversion_RLarge",nBinsR,0.,200.);
        fESDList[iCut]->Add(hMCConversionRRejLarge[iCut]);
        hMCConversionRRejSmall[iCut] = new TH1F("MC_Conversion_RSmall","MC_Conversion_RSmall",nBinsR,0.,200.);
        fESDList[iCut]->Add(hMCConversionRRejSmall[iCut]);

        if (fDoMultWeights>0 && fIsMC>0 ) {
          hMCConversionRPt[iCut] ->Sumw2();
        }

        hMCTrueConversionRPhi[iCut] = new TH2F("ESD_TrueConversion_RPhi","ESD_TrueConversion_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionRPhi[iCut]);
        hMCTrueConversionRPhiFromConv[iCut] = new TH2F("ESD_TrueConversion_RPhi_FromConv","ESD_TrueConversion_RPhi_FromConv",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionRPhiFromConv[iCut]);

        hMCTrueConversionREta[iCut] = new TH2F("ESD_TrueConversion_REta","ESD_TrueConversion_REta",nBinsEta,-2.,2.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionREta[iCut]);
        hMCTrueConversionRPt[iCut]  = new TH2F("ESD_TrueConversion_RPt","ESD_TrueConversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionRPt[iCut]);
        hMCTrueConversionWOWeightRPt[iCut]  = new TH2F("ESD_TrueConversionWOWeight_RPt","ESD_TrueConversionWOWeight_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionWOWeightRPt[iCut]);
      	hMCTrueConversionRPtMCRPt[iCut]  = new TH2F("ESD_TrueConversion_RPtMCRPt","ESD_TrueConversion_RPtMCRPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionRPtMCRPt[iCut]);
       	hMCTrueConversionWOWeightRPtMCRPt[iCut]  = new TH2F("ESD_TrueConversionWOWeight_RPtMCRPt","ESD_TrueConversionWOWeight_RPtMCRPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueConversionWOWeightRPtMCRPt[iCut]);
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

      	hMCTruePrimConversionRPt[iCut]  = new TH2F("ESD_TruePrimConversion_RPt","ESD_TruePrimConversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTruePrimConversionRPt[iCut]);
      	hMCTruePrimConversionWOWeightRPt[iCut]  = new TH2F("ESD_TruePrimConversionWOWeight_RPt","ESD_TruePrimConversionWOWeight_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTruePrimConversionWOWeightRPt[iCut]);


        hMCTrueSecConversionRPt[iCut]  = new TH2F("ESD_TrueSecConversion_RPt","ESD_TrueSecConversion_RPt",nBinsPt,0.,20.,nBinsR,0.,200.);
        fTrueList[iCut]->Add(hMCTrueSecConversionRPt[iCut]);

        hMCTrueSecondaryConvGammaRPt[iCut]  = new TH3F("ESD_TrueSecondaryConvGamma_Pt", "ESD_TrueSecondaryConvGamma_Pt", nBinsPt,0.,20.,nBinsR,0.,200. , 4, -0.5, 3.5);
        hMCTrueSecondaryConvGammaRPt[iCut]->GetZaxis()->SetBinLabel(1,"K0s");
        hMCTrueSecondaryConvGammaRPt[iCut]->GetZaxis()->SetBinLabel(2,"K0l");
        hMCTrueSecondaryConvGammaRPt[iCut]->GetZaxis()->SetBinLabel(3,"Lambda");
        hMCTrueSecondaryConvGammaRPt[iCut]->GetZaxis()->SetBinLabel(4,"rest");
        fTrueList[iCut]->Add(hMCTrueSecondaryConvGammaRPt[iCut]);
        hMCTrueSecondaryConvGammaMCRPt[iCut]  = new TH3F("ESD_TrueSecondaryConvGamma_MCPt", "ESD_TrueSecondaryConvGamma_MCPt", nBinsPt, 0.,20.,nBinsR,0.,200., 4, -0.5, 3.5);
        hMCTrueSecondaryConvGammaMCRPt[iCut]->GetZaxis()->SetBinLabel(1,"K0s");
        hMCTrueSecondaryConvGammaMCRPt[iCut]->GetZaxis()->SetBinLabel(2,"K0l");
        hMCTrueSecondaryConvGammaMCRPt[iCut]->GetZaxis()->SetBinLabel(3,"Lambda");
        hMCTrueSecondaryConvGammaMCRPt[iCut]->GetZaxis()->SetBinLabel(4,"rest");
        fTrueList[iCut]->Add(hMCTrueSecondaryConvGammaMCRPt[iCut]);




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

	hMCTrueConversionAsymP[iCut]               = new TH2F("ESD_TrueConversionMapping_AsymP","ESD_TrueConversionMapping_AsymP",nBinsPt,0.01,20.,500,0.,1.);
	fTrueList[iCut]->Add(hMCTrueConversionAsymP[iCut]);

	AxisAfter = hMCTrueConversionAsymP[iCut]->GetXaxis();
	AxisAfter->Set(bins, newBins);

	if ((fDoMultWeights>0 && fIsMC>0) || ( fDoMaterialBudgetWeightingOfGammasForTrueMesons>0 && fIsMC>0) ) {
	  hMCTrueConversionRPt[iCut] ->Sumw2();
	  hMCTrueConversionRPtMCRPt[iCut] -> Sumw2();
	  hMCTrueConversionDCA[iCut] ->Sumw2();
	  hMCTrueConversionChi2[iCut] ->Sumw2();
	  hMCTrueConversionAsymP[iCut] ->Sumw2();
	  hMCTruePrimConversionRPt[iCut] ->Sumw2();
	  hMCTrueSecConversionRPt[iCut] ->Sumw2();
	  hMCTrueSecondaryConvGammaRPt[iCut]->Sumw2();
	  hMCTrueSecondaryConvGammaMCRPt[iCut]->Sumw2();
	}


    }

    Int_t nPBins =12;
    Int_t nEtaBins =20;
    Int_t nSigmaDeDxBins=100;
    Double_t *arrPBinning = new Double_t[13];
    for( Int_t i=0;i<nPBins+1;i++){
      if(i==0){
        arrPBinning[i]= 0.05;
      }else if(i>0 && i<11){
        arrPBinning[i]= 0.1*i;
      }else if(i==11){
        arrPBinning[i]= 2.0;
      }else if(i==12){
        arrPBinning[i]= 10.0;
      }
      //cout<< "pbins::"<< i << " " <<  arrPBinning[i]<< endl;
    }
    Double_t *arrEtaBinning      = new Double_t[21];
    for( Int_t i=0;i<nEtaBins+1;i++){
      arrEtaBinning[i]= -1.+0.1*i;
      //cout<< "Etabins::"<< i << " " <<  arrEtaBinning[i]<< endl;
    }
    Double_t *arrSigmaDeDxBinning      = new Double_t[101];
    for( Int_t i=0;i<nSigmaDeDxBins+1;i++){
      arrSigmaDeDxBinning[i]= -5.+0.1*i;
      //cout<< "dedx::"<< i << " " <<  arrSigmaDeDxBinning[i]<< endl;
    }

    if (fDoDeDxMaps>0) {
      fDeDxMapList[iCut] = new TList();
      fDeDxMapList[iCut] ->SetName(Form("%s_%s  dEdx Maps",cutstringEvent.Data() ,cutstringPhoton.Data()));
      fDeDxMapList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fDeDxMapList[iCut]);

      hElectrondEdxMapsR0[iCut]= new TH3F("R0 electron sigma dEdx P Eta","R0 electron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      hPositrondEdxMapsR0[iCut]= new TH3F("R0 positron sigma dEdx P Eta","R0 positron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      fDeDxMapList[iCut]->Add( hPositrondEdxMapsR0[iCut]);
      fDeDxMapList[iCut]->Add( hElectrondEdxMapsR0[iCut]);

      hElectrondEdxMapsR1[iCut]= new TH3F("R1 electron sigma dEdx P Eta","R1 electron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      hPositrondEdxMapsR1[iCut]= new TH3F("R1 positron sigma dEdx P Eta","R1 positron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      fDeDxMapList[iCut]->Add( hPositrondEdxMapsR1[iCut]);
      fDeDxMapList[iCut]->Add( hElectrondEdxMapsR1[iCut]);

      hElectrondEdxMapsR2[iCut]= new TH3F("R2 electron sigma dEdx P Eta","R2 electron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      hPositrondEdxMapsR2[iCut]= new TH3F("R2 positron sigma dEdx P Eta","R2 positron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      fDeDxMapList[iCut]->Add( hPositrondEdxMapsR2[iCut]);
      fDeDxMapList[iCut]->Add( hElectrondEdxMapsR2[iCut]);

      hElectrondEdxMapsR3[iCut]= new TH3F("R3 electron sigma dEdx P Eta", "R3 electron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      hPositrondEdxMapsR3[iCut]= new TH3F("R3 positron sigma dEdx P Eta","R3 positron sigma dEdx P Eta", nSigmaDeDxBins, arrSigmaDeDxBinning, nEtaBins,arrEtaBinning, nPBins, arrPBinning);
      fDeDxMapList[iCut]->Add( hPositrondEdxMapsR3[iCut]);
      fDeDxMapList[iCut]->Add( hElectrondEdxMapsR3[iCut]);


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

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMaterialHistos::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }
  }
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskMaterialHistos::UserExec(Option_t *){

  fInputEvent = InputEvent();
  if (fInputEvent==NULL) return;


  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    hBCNumber[iCut]->Fill(fBCNumber);
  }

  if(fDoSelectBCNumber) DoSelectBCNumbers();

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
    if(fDoSelectBCNumber) hBCNumberSelected[iCut]->Fill(fBCNumber);
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


    fNESDtracksEta08 = CountTracks08(); // Estimate Event Multiplicity
    fNESDtracksEta08pt200 = CountTracks08pt200(); // Estimate Event Multiplicity
    fNESDtracksEta08pt300 = CountTracks08pt300(); // Estimate Event Multiplicity
    fNESDtracksEta08pt400 = CountTracks08pt400(); // Estimate Event Multiplicity
    fNESDtracksEta0814 = CountTracks0814(); // Estimate Event Multiplicity
    fNESDtracksEta14 = fNESDtracksEta08 + fNESDtracksEta0814;

    // cout<< " pt dependent::  " <<  fNESDtracksEta08 << "  " << fNESDtracksEta08pt200 << "  "<< fNESDtracksEta08pt300 << "  " << fNESDtracksEta08pt400<< endl;

    if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetDoElecDeDxPostCalibration()){
      if(!((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->LoadElecDeDxPostCalibration(fInputEvent->GetRunNumber())){
        AliFatal(Form("ERROR: LoadElecDeDxPostCalibration returned kFALSE for %d despite being requested!",fInputEvent->GetRunNumber()));
      }
    }

    hNGoodESDTracksEta08[iCut]->Fill(fNESDtracksEta08);
    hNGoodESDTracksEta08pt200[iCut]->Fill(fNESDtracksEta08pt200);
    hNGoodESDTracksEta08pt300[iCut]->Fill(fNESDtracksEta08pt300);
    hNGoodESDTracksEta08pt400[iCut]->Fill(fNESDtracksEta08pt400);
    hNGoodESDTracksEta14[iCut]->Fill(fNESDtracksEta14);
    hNGoodESDTracksEta08_14[iCut]->Fill(fNESDtracksEta0814);

    if(fDoMultWeights && fIsMC > 0) {
      fWeightMultMC = 1.;

      if ( fDoMultWeights==1) {
        fWeightMultMC = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fNESDtracksEta08);
      }else if ( fDoMultWeights==2) {
        fWeightMultMC = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());
      }else {
        fWeightMultMC = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());
      }

      hNGoodESDTracksWeightedEta08[iCut]->Fill(fNESDtracksEta08, fWeightMultMC);
      hNGoodESDTracksWeightedEta08pt200[iCut]->Fill(fNESDtracksEta08pt200, fWeightMultMC);
      hNGoodESDTracksWeightedEta08pt300[iCut]->Fill(fNESDtracksEta08pt300, fWeightMultMC);
      hNGoodESDTracksWeightedEta08pt400[iCut]->Fill(fNESDtracksEta08pt400, fWeightMultMC);
 
     fHistoNV0TracksWeighted[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(),fWeightMultMC);
    }

    fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());


    // Calculation of Multiplicity weight moved before ProcessMCPhotons
    // fWeightMultMC shuld also be inserted to input pT distributions .

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
void AliAnalysisTaskMaterialHistos::ProcessMCPhotons(){

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();


  // Loop over all primary MC particle

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

    //  for(Int_t i = 0; i < fMCEvent->GetNumberOfPrimaries(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      // fill primary histogram
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
      }

      Float_t weighted= 1;
      if (particle->Pt()>0.005){
        weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForGamma(i, fMCEvent, fInputEvent);
        //	cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
      }


      if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
        hMCAllGammaPt[fiCut]->Fill(particle->Pt(),weighted*fWeightMultMC); // All MC Gamma
        hMCAllGammaWOWeightPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
      }
      if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
        TParticle* daughter1 = (TParticle *)fMCEvent->Particle(particle->GetFirstDaughter());

        Double_t phiFromConv = TMath::ATan2(daughter1->Vy(),daughter1->Vx());
        if (phiFromConv<0) phiFromConv+=TMath::TwoPi();

        hMCConversionRPhi[fiCut]->Fill(particle->Phi(),daughter1->R());
        hMCConversionRPhiFromConv[fiCut]->Fill(phiFromConv,daughter1->R());
        hMCConversionREta[fiCut]->Fill(particle->Eta(),daughter1->R());
        hMCConversionWOWeightRPt[fiCut]->Fill(particle->Pt(),daughter1->R());
        hMCConversionRPt[fiCut]->Fill(particle->Pt(),daughter1->R(),weighted*fWeightMultMC);

        if(daughter1->R() < 75. || daughter1->R() > 85.) hMCConversionRRejSmall[fiCut]->Fill(daughter1->R());
        if(daughter1->R() < 70. || daughter1->R() > 90.) hMCConversionRRejLarge[fiCut]->Fill(daughter1->R());
      } // Converted MC Gamma
    } else {
      // fill secondary histograms
      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
      }
      if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
        if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
          if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
            hMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,fWeightMultMC);
          } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
            hMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,fWeightMultMC);
          } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
            hMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,fWeightMultMC);
          } else {
            //  if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 &&
            // fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
            hMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightMultMC);
          }
        } else {
          hMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,fWeightMultMC);
        }
        if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
          TParticle* tmpDaughter1 = (TParticle *)fMCEvent->Particle(particle->GetFirstDaughter());
          if (particle->GetMother(0) > -1 && fMCEvent->Particle(particle->GetMother(0))->GetMother(0) > -1) {
            if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
              hMCSecondaryConvGammaPtR[fiCut]->Fill(particle->Pt(),tmpDaughter1->R(),0.,fWeightMultMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
              hMCSecondaryConvGammaPtR[fiCut]->Fill(particle->Pt(),tmpDaughter1->R(),1.,fWeightMultMC);
            } else if (fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
              hMCSecondaryConvGammaPtR[fiCut]->Fill(particle->Pt(),tmpDaughter1->R(),2.,fWeightMultMC);
            } else {
            //              if ( !(TMath::Abs(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()) == 11 &&
            //   fMCEvent->Particle(fMCEvent->Particle(particle->GetMother(0))->GetMother(0))->GetPdgCode() == 22) )
              hMCSecondaryConvGammaPtR[fiCut]->Fill(particle->Pt(),tmpDaughter1->R(),3.,fWeightMultMC);
            }
          } else {
            hMCSecondaryConvGammaPtR[fiCut]->Fill(particle->Pt(),tmpDaughter1->R(),3.,fWeightMultMC);
          }
        }
      }
    }
  }

}

///________________________________________________________________________
void AliAnalysisTaskMaterialHistos::ProcessPhotons(){

  Double_t magField = fInputEvent->GetMagneticField();

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

    AliPIDResponse* pidResponse = ((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();

    AliVTrack * negTrack = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetTrack(fInputEvent, gamma->GetTrackLabelNegative());
    AliVTrack * posTrack = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetTrack(fInputEvent, gamma->GetTrackLabelPositive());

    Short_t Charge    = 1;
    Double_t electronNSigmaTPC = pidResponse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
    Double_t electronNSigmaTPCCor=0.;
    Double_t positronNSigmaTPC = pidResponse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
    Double_t positronNSigmaTPCCor=0.;
    Double_t P=0.;
    Double_t Eta=0.;

    if( ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetDoElecDeDxPostCalibration() ){
      Charge = negTrack->Charge();
      P = negTrack->P();
      Eta = negTrack->Eta();
      electronNSigmaTPCCor = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetCorrectedElectronTPCResponse(Charge,electronNSigmaTPC,P,Eta,negTrack->GetTPCNcls(),gamma->GetConversionRadius());

      Charge = posTrack->Charge();
      P = posTrack->P();
      Eta = posTrack->Eta();
      positronNSigmaTPCCor = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetCorrectedElectronTPCResponse(Charge,positronNSigmaTPC,P,Eta,posTrack->GetTPCNcls(),gamma->GetConversionRadius());
    }

    hESDConversionWOWeightRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC);
    if(fIsMC==0) hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC);
    //In case of MC, this histogram is filled with pT weights for primary photons. Weights not applied for secondaries and also not for contaminations

    Double_t phiFromConv = TMath::ATan2(gamma->GetConversionY(),gamma->GetConversionX());
    if (phiFromConv<0) phiFromConv+=TMath::TwoPi();

    hESDConversionRPhi[fiCut]->Fill(gamma->GetPhotonPhi(),gamma->GetConversionRadius());
    hESDConversionRPhiFromConv[fiCut]->Fill(phiFromConv,gamma->GetConversionRadius());
    hESDConversionRZ[fiCut]->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());
    hESDConversionREta[fiCut]->Fill(gamma->GetPhotonEta(),gamma->GetConversionRadius());


    if( ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetDoElecDeDxPostCalibration() ){
      if(negTrack->GetTPCsignal()){
        hElectronRdEdx[fiCut]->Fill(negTrack->GetTPCsignal(),gamma->GetConversionRadius());
        hElectronRNSigmadEdx[fiCut]->Fill( electronNSigmaTPCCor, gamma->GetConversionRadius());
      }
      if(posTrack->GetTPCsignal()){
        hPositronRdEdx[fiCut]->Fill(posTrack->GetTPCsignal(),gamma->GetConversionRadius());
        hPositronRNSigmadEdx[fiCut]->Fill( positronNSigmaTPCCor, gamma->GetConversionRadius());
      }
    }else{
     if(negTrack->GetTPCsignal()){
        hElectronRdEdx[fiCut]->Fill(negTrack->GetTPCsignal(),gamma->GetConversionRadius());
        hElectronRNSigmadEdx[fiCut]->Fill( electronNSigmaTPC, gamma->GetConversionRadius());
      }
      if(posTrack->GetTPCsignal()){
        hPositronRdEdx[fiCut]->Fill(posTrack->GetTPCsignal(),gamma->GetConversionRadius());
        hPositronRNSigmadEdx[fiCut]->Fill( positronNSigmaTPC, gamma->GetConversionRadius());
      }

    }

    if(gamma->GetConversionRadius() < 75. || gamma->GetConversionRadius() > 85.) hESDConversionRRejSmall[fiCut]->Fill(gamma->GetConversionRadius());
    if(gamma->GetConversionRadius() < 70. || gamma->GetConversionRadius() > 90.) hESDConversionRRejLarge[fiCut]->Fill(gamma->GetConversionRadius());

    if(fInputEvent->IsA()==AliESDEvent::Class()){
      AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
      if(esdEvent){
        AliESDv0 *v0 = esdEvent->GetV0(gamma->GetV0Index());
        hESDConversionDCA[fiCut]->Fill(v0->GetDcaV0Daughters(),fWeightMultMC);
      }
    }
    hESDConversionPsiPair[fiCut]->Fill(gamma->GetPsiPair(),fWeightMultMC);
    hESDConversionChi2[fiCut]->Fill(gamma->GetChi2perNDF(),fWeightMultMC);

    hESDConversionMass[fiCut]->Fill(gamma->GetInvMassPair(),fWeightMultMC);


    if(gamma->GetPhotonP()!=0 && negTrack->P()!=0) {
      if(gamma->GetConversionRadius() > 5. ){
        hESDConversionAsymP[fiCut]->Fill(gamma->GetPhotonP(),negTrack->P()/gamma->GetPhotonP(),fWeightMultMC);
      }
    }

    if(fDoDeDxMaps > 0 ) {
      if( ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetDoElecDeDxPostCalibration()){
        if(gamma->GetConversionRadius() < 33.5){
          hElectrondEdxMapsR0[fiCut]->Fill(electronNSigmaTPCCor, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR0[fiCut]->Fill(positronNSigmaTPCCor, gamma->GetPhotonEta(), posTrack->P());
        }else if (gamma->GetConversionRadius() > 33.5  && gamma->GetConversionRadius() < 72.){
          hElectrondEdxMapsR1[fiCut]->Fill(electronNSigmaTPCCor, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR1[fiCut]->Fill(positronNSigmaTPCCor, gamma->GetPhotonEta(), posTrack->P());
        }else if (gamma->GetConversionRadius() > 72.  && gamma->GetConversionRadius() < 145.){
          hElectrondEdxMapsR2[fiCut]->Fill(electronNSigmaTPCCor, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR2[fiCut]->Fill(positronNSigmaTPCCor, gamma->GetPhotonEta(), posTrack->P());
        }else if (gamma->GetConversionRadius() > 145.  && gamma->GetConversionRadius() < 180.){
          hElectrondEdxMapsR3[fiCut]->Fill(electronNSigmaTPCCor, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR3[fiCut]->Fill(positronNSigmaTPCCor, gamma->GetPhotonEta(), posTrack->P());
        }
      }else{
        if(gamma->GetConversionRadius() < 33.5){
          hElectrondEdxMapsR0[fiCut]->Fill(electronNSigmaTPC, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR0[fiCut]->Fill(positronNSigmaTPC, gamma->GetPhotonEta(), posTrack->P());
        }else if (gamma->GetConversionRadius() > 33.5  && gamma->GetConversionRadius() < 72.){
          hElectrondEdxMapsR1[fiCut]->Fill(electronNSigmaTPC, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR1[fiCut]->Fill(positronNSigmaTPC, gamma->GetPhotonEta(), posTrack->P());
        }else if (gamma->GetConversionRadius() > 72.  && gamma->GetConversionRadius() < 145.){
          hElectrondEdxMapsR2[fiCut]->Fill(electronNSigmaTPC, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR2[fiCut]->Fill(positronNSigmaTPC, gamma->GetPhotonEta(), posTrack->P());
        }else if (gamma->GetConversionRadius() > 145.  && gamma->GetConversionRadius() < 180.){
          hElectrondEdxMapsR3[fiCut]->Fill(electronNSigmaTPC, gamma->GetPhotonEta(), negTrack->P());
          hPositrondEdxMapsR3[fiCut]->Fill(positronNSigmaTPC, gamma->GetPhotonEta(), posTrack->P());
        }
      }
    }


    fKind = 9;
    Int_t pdgCodePos = 0.;
    Int_t pdgCodeNeg = 0.;

    if(fIsMC>0){

      const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX   = primVtxMC->GetX();
      Double_t mcProdVtxY   = primVtxMC->GetY();
      Double_t mcProdVtxZ   = primVtxMC->GetZ();

      TParticle *posDaughter = gamma->GetPositiveMCDaughter(fMCEvent);
      TParticle *negDaughter = gamma->GetNegativeMCDaughter(fMCEvent);
      TParticle *Photon = gamma->GetMCParticle(fMCEvent);
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
        hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC);
      } else {
        //cout << "same mother" << endl;
        pdgCodePos = posDaughter->GetPdgCode();
        pdgCodeNeg = negDaughter->GetPdgCode();
        Int_t pdgCode;
        pdgCode = gamma->GetMCParticle(fMCEvent)->GetPdgCode();
        if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11){
          fKind = 2; // combinatorics from hadronic decays
          hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC);
        }else if ( !(pdgCodeNeg==pdgCodePos)){
          Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          if(pdgCode == 111) {
            fKind = 3; // pi0 Dalitz
            hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC);
          }else if (pdgCode == 221){
            fKind = 4; // eta Dalitz
            hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC);
          }else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
            if(pdgCode == 22 && gammaIsPrimary){
              fKind = 0; // primary photons
            } else if (pdgCode == 22){
              fKind = 5; //secondary photons

              //----------Splitting of secondaries. Part taken from AliAnalysisTaskGammaConvV1-----------------
              if( Photon->GetMother(0) > -1 && fMCEvent->Particle(Photon->GetMother(0))->GetMother(0) > -1){
                if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
                  hMCTrueSecondaryConvGammaRPt[fiCut]->Fill(gamma->Pt(),gamma->GetConversionRadius(),0.,fWeightMultMC);
                  hMCTrueSecondaryConvGammaMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R(),0.,fWeightMultMC);
                  // hMCTrueSecondaryConvGammaFromXFromK0sMCPtESDPtR[fiCut]->Fill(Photon->Pt(),gamma->Pt());
                } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 130) {
                  hMCTrueSecondaryConvGammaRPt[fiCut]->Fill(gamma->Pt(),gamma->GetConversionRadius(),1.,fWeightMultMC);
                  hMCTrueSecondaryConvGammaMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R(),1.,fWeightMultMC);
                  // hMCTrueSecondaryConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),gamma->Pt());
                } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 3122) {
                  hMCTrueSecondaryConvGammaRPt[fiCut]->Fill(gamma->Pt(),gamma->GetConversionRadius(),2.,fWeightMultMC);
                  hMCTrueSecondaryConvGammaMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R(),2.,fWeightMultMC);
                  // hMCTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),gamma->Pt());
                } else if (fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 221) {
                  hMCTrueSecondaryConvGammaRPt[fiCut]->Fill(Photon->Pt(),gamma->GetConversionRadius(),3.,fWeightMultMC);
                  hMCTrueSecondaryConvGammaMCRPt[fiCut]->Fill(gamma->Pt(),negDaughter->R(),3.,fWeightMultMC);
                } else {
                  // if ( !(TMath::Abs(fMCEvent->Particle(Photon->GetMother(0))->GetPdgCode()) == 11 && fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 22) ) {
                  hMCTrueSecondaryConvGammaRPt[fiCut]->Fill(gamma->Pt(),gamma->GetConversionRadius(),3.,fWeightMultMC);
                  hMCTrueSecondaryConvGammaMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R(),3.,fWeightMultMC);
                  //         }
                }
              } else {
                hMCTrueSecondaryConvGammaRPt[fiCut]->Fill(gamma->Pt(),gamma->GetConversionRadius(),3.,fWeightMultMC);
                hMCTrueSecondaryConvGammaMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R(),3.,fWeightMultMC);
              }
              // End spliting of secondaries
            }
          } else   fKind = 9; //garbage
        } else fKind = 9; //garbage
      }
      Float_t weighted = 1.;
      weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForGamma(posDaughter->GetMother(0), fMCEvent, fInputEvent);

      Double_t phiFromConv = TMath::ATan2(gamma->GetConversionY(),gamma->GetConversionX());
      if (phiFromConv<0) phiFromConv+=TMath::TwoPi();

      Float_t weightMatBudget = 1.;
      if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetMaterialBudgetWeightsInitialized()) {
	weightMatBudget = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(gamma,magField);
      }



      if(fKind==0 || fKind==5){
        hMCTrueConversionRPhi[fiCut]->Fill(gamma->GetPhotonPhi(),gamma->GetConversionRadius());
        hMCTrueConversionRPhiFromConv[fiCut]->Fill(phiFromConv,gamma->GetConversionRadius());
        hMCTrueConversionRZ[fiCut]->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());
        hMCTrueConversionREta[fiCut]->Fill(gamma->GetPhotonEta(),gamma->GetConversionRadius());

        if(fKind==0) hMCTrueConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),weighted*fWeightMultMC*weightMatBudget);
        if(fKind==5) hMCTrueConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC*weightMatBudget);

        hMCTrueConversionWOWeightRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius());


	if(fKind==0) hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),weighted*fWeightMultMC*weightMatBudget);
	if(fKind==5) hESDConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC*weightMatBudget);

	if(fKind==0) hMCTruePrimConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),weighted*fWeightMultMC*weightMatBudget);
	if(fKind==0) hMCTruePrimConversionWOWeightRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius());
	if(fKind==5) hMCTrueSecConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),fWeightMultMC*weightMatBudget);

	if(fKind==0)hMCTrueConversionRPtMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R(),weighted*fWeightMultMC*weightMatBudget);
	if(fKind==5)hMCTrueConversionRPtMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R(),fWeightMultMC*weightMatBudget);


        hMCTrueConversionWOWeightRPtMCRPt[fiCut]->Fill(Photon->Pt(),negDaughter->R());

        if(gamma->GetConversionRadius() < 75. || gamma->GetConversionRadius() > 85.) hMCTrueConversionRRejSmall[fiCut]->Fill(gamma->GetConversionRadius());
        if(gamma->GetConversionRadius() < 70. || gamma->GetConversionRadius() > 90.) hMCTrueConversionRRejLarge[fiCut]->Fill(gamma->GetConversionRadius());

        hMCTrueConversionPsiPair[fiCut]->Fill(gamma->GetPsiPair(),fWeightMultMC);
        hMCTrueConversionChi2[fiCut]->Fill(gamma->GetChi2perNDF(),weighted*fWeightMultMC);

        hMCTrueConversionMass[fiCut]->Fill(gamma->GetInvMassPair(),fWeightMultMC);
	      if(gamma->GetPhotonP()!=0 && negTrack->P()!=0) {
	       if(gamma->GetConversionRadius() > 5.){
	         hMCTrueConversionAsymP[fiCut]->Fill(gamma->GetPhotonP(),negTrack->P()/gamma->GetPhotonP(),fWeightMultMC);
      	  }

        }



        if(fInputEvent->IsA()==AliESDEvent::Class()){
            AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
            if(esdEvent){
                AliESDv0 *v0 = esdEvent->GetV0(gamma->GetV0Index());
                hMCTrueConversionDCA[fiCut]->Fill(v0->GetDcaV0Daughters());
            }
        }

      } else if(fKind==3){
        hMCTruePi0DalConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),weighted*fWeightMultMC);
        hMCTruePi0DalConversionEta[fiCut]->Fill(gamma->GetPhotonEta());
      } else if(fKind==4){
        hMCTrueEtaDalConversionRPt[fiCut]->Fill(gamma->GetPhotonPt(),gamma->GetConversionRadius(),weighted*fWeightMultMC);
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
Int_t AliAnalysisTaskMaterialHistos::CountTracks08(){

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
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/TMath::Power(pt,1.1)");
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
            EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);

        } else { // else use 2011 version of track cuts
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
        }
        EsdTrackCuts->SetMaxDCAToVertexZ(2);
        EsdTrackCuts->SetEtaRange(-0.8, 0.8);
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
Int_t AliAnalysisTaskMaterialHistos::CountTracks08pt200(){

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
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/TMath::Power(pt,1.1)");
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
            EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);

        } else { // else use 2011 version of track cuts
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
        }
        EsdTrackCuts->SetMaxDCAToVertexZ(2);
        EsdTrackCuts->SetEtaRange(-0.8, 0.8);
        EsdTrackCuts->SetPtRange(0.2);
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
Int_t AliAnalysisTaskMaterialHistos::CountTracks08pt300(){

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
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/TMath::Power(pt,1.1)");
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
            EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);

        } else { // else use 2011 version of track cuts
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
        }
        EsdTrackCuts->SetMaxDCAToVertexZ(2);
        EsdTrackCuts->SetEtaRange(-0.8, 0.8);
        EsdTrackCuts->SetPtRange(0.3);
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
Int_t AliAnalysisTaskMaterialHistos::CountTracks08pt400(){

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
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/TMath::Power(pt,1.1)");
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
            EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);

        } else { // else use 2011 version of track cuts
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
        }
        EsdTrackCuts->SetMaxDCAToVertexZ(2);
        EsdTrackCuts->SetEtaRange(-0.8, 0.8);
        EsdTrackCuts->SetPtRange(0.4);
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
Int_t AliAnalysisTaskMaterialHistos::CountTracks0814(){

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

    EsdTrackCuts->SetEtaRange(0.8, 1.4);
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

    EsdTrackCuts->SetEtaRange(-1.4, -0.8);
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
void AliAnalysisTaskMaterialHistos::DoSelectBCNumbers(){
  if(fIsMC == 0 ){
    fRunNumber  = fInputEvent->GetRunNumber();
    fBCNumber = fInputEvent->GetBunchCrossNumber();
    if(fRunNumber >= 265332 && fRunNumber <= 265344){//Fill 5506
      if(1390 > fBCNumber || fBCNumber > 1540){
	return;
      }
    }else if(fRunNumber >= 265377 && fRunNumber <= 265388){//Fill 5507
      if(3010 > fBCNumber || fBCNumber > 3140){
	return;
      }
    }
  }
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
