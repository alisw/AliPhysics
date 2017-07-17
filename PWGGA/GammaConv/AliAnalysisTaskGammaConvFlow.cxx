/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: Andrea Dubla, Redmer Alexander Bertens, Friederike Bock	      *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its	  *
* documentation strictly for non-commercial purposes is hereby granted	  *
* without fee, provided that the above copyright notice appears in all	  *
* copies and that both the copyright notice and this permission notice	  *
* appear in the supporting documentation. The authors make no claims	  *
* about the suitability of this software for any purpose. It is	      *
* provided "as is" without express or implied warranty.       		      *
**************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------

// Class used to do analysis on conversion pairs
//---------------------------------------------
///////////////////////////////////////////////
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
#include "AliAnalysisTaskGammaConvFlow.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"

#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"

class AliFlowTrackCuts;

using namespace std;

ClassImp(AliAnalysisTaskGammaConvFlow)

//________________________________________________________________________
AliAnalysisTaskGammaConvFlow::AliAnalysisTaskGammaConvFlow(): AliAnalysisTaskSE(),
fV0Reader(NULL),
fV0ReaderName("V0ReaderV1"),
fBGHandler(NULL),
fBGHandlerRP(NULL),
fInputEvent(NULL),

fCutFolder(NULL),
fESDList(NULL),
fBackList(NULL),
fMotherList(NULL),
fPhotonDCAList(NULL),
fMesonDCAList(NULL),
//fTrueList(NULL),
//fMCList(NULL),
fHeaderNameList(NULL),
fOutputContainer(0),
fReaderGammas(NULL),
fGammaCandidates(NULL),
fEventCutArray(NULL),
fEventCuts(NULL),
fCutArray(NULL),
fConversionCuts(NULL),
fMesonCutArray(NULL),
fMesonCuts(NULL),
hESDConvGammaPt(NULL),
hInvMassPair(NULL),
hLTMPt(NULL),
hLTMPt_MC(NULL),
hPt_TruePt(NULL),
hdPhidRcandidates(NULL),
hdPhidRcandidates_MCsigsig(NULL),
hdPhidRcandidates_MCbkgsig(NULL),
hdPhidRcandidates_MCbkgbkg(NULL),
hKappaTPC(NULL),
hKappaTPC_after(NULL),
hKappaTPC_Temp0(NULL),
hKappaTPC_Temp1(NULL),
hKappaTPC_Temp2(NULL),
hKappaTPC_Temp3(NULL),
hKappaTPC_Temp4(NULL),
hKappaTPC_Temp5(NULL),
hKappaTPC_Temp6(NULL),
hKappaTPC_Temp7(NULL),
hKappaTPC_Temp8(NULL),
hKappaTPC_Temp9(NULL),
hKappaTPC_Temp10(NULL),
hESDConvGammaR(NULL),
hESDConvGammaEta(NULL),
//tESDConvGammaPtDcazCat(NULL),
fPtGamma(0),
fDCAzPhoton(0),
fRConvPhoton(0),
fEtaPhoton(0),
iCatPhoton(0),
iPhotonMCInfo(0),
hESDMotherInvMassPt(NULL),
//sESDMotherInvMassPtZM(NULL),
hESDMotherBackInvMassPt(NULL),
//sESDMotherBackInvMassPtZM(NULL),
hESDMotherInvMassEalpha(NULL),
hESDMotherPi0PtY(NULL),
hESDMotherEtaPtY(NULL),
hESDMotherPi0PtAlpha(NULL),
hESDMotherEtaPtAlpha(NULL),
hESDMotherPi0PtOpenAngle(NULL),
hESDMotherEtaPtOpenAngle(NULL),

hNEvents(NULL),
hNGoodESDTracks(NULL),
hNGammaCandidates(NULL),
hNGoodESDTracksVsNGammaCanditates(NULL),
hNV0Tracks(NULL),
hEtaShift(NULL),
tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
fInvMass(0),
fPt(0),
fDCAzGammaMin(0),
fDCAzGammaMax(0),
iFlag(0),
iMesonMCInfo(0),
fEventPlaneAngle(-100),
fRandom(0),
fnGammaCandidates(0),
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
fDoMesonAnalysis(kTRUE),
fDoMesonQA(0),
fDoPhotonQA(0),
fIsFromMBHeader(kTRUE),
fhistoEPVZ(NULL),
fMinMass(-1),
fMaxMass(10),
fMinKappa(-1),
fMaxKappa(100),
fFilterVariable(1),
fMinFilter(0),
fMaxFilter(0.2),
fApplydPhidRCut(0),
fPerformExtraStudies(0),
fDebug(0),
fCutsRP(0),
fNullCuts(0), 
fFlowEvent(NULL),
fIsMC(0),
fMCEvent(NULL)
{
  // DefineOutput(1, TList::Class());
  // DefineOutput(2, AliFlowEventSimple::Class());
}

//________________________________________________________________________
AliAnalysisTaskGammaConvFlow::AliAnalysisTaskGammaConvFlow(const char *name):
AliAnalysisTaskSE(name),
fV0Reader(NULL),
fV0ReaderName("V0ReaderV1"),
fBGHandler(NULL),
fBGHandlerRP(NULL),
fInputEvent(NULL),
fCutFolder(NULL),
fESDList(NULL),
fBackList(NULL),
fMotherList(NULL),
fPhotonDCAList(NULL),
fMesonDCAList(NULL),
//fTrueList(NULL),
//fMCList(NULL),
fHeaderNameList(NULL),
fOutputContainer(0),
fReaderGammas(NULL),
fGammaCandidates(NULL),
fEventCutArray(NULL),
fEventCuts(NULL),
fCutArray(NULL),
fConversionCuts(NULL),
fMesonCutArray(NULL),
fMesonCuts(NULL),
hESDConvGammaPt(NULL),
hInvMassPair(NULL),
hLTMPt(NULL),
hLTMPt_MC(NULL),
hPt_TruePt(NULL),
hdPhidRcandidates(NULL),
hdPhidRcandidates_MCsigsig(NULL),
hdPhidRcandidates_MCbkgsig(NULL),
hdPhidRcandidates_MCbkgbkg(NULL),
hKappaTPC(NULL),
hKappaTPC_after(NULL),
hKappaTPC_Temp0(NULL),
hKappaTPC_Temp1(NULL),
hKappaTPC_Temp2(NULL),
hKappaTPC_Temp3(NULL),
hKappaTPC_Temp4(NULL),
hKappaTPC_Temp5(NULL),
hKappaTPC_Temp6(NULL),
hKappaTPC_Temp7(NULL),
hKappaTPC_Temp8(NULL),
hKappaTPC_Temp9(NULL),
hKappaTPC_Temp10(NULL),
hESDConvGammaR(NULL),
hESDConvGammaEta(NULL),
//tESDConvGammaPtDcazCat(NULL),
fPtGamma(0),
fDCAzPhoton(0),
fRConvPhoton(0),
fEtaPhoton(0),
iCatPhoton(0),
iPhotonMCInfo(0),
hESDMotherInvMassPt(NULL),
//sESDMotherInvMassPtZM(NULL),
hESDMotherBackInvMassPt(NULL),
//sESDMotherBackInvMassPtZM(NULL),
hESDMotherInvMassEalpha(NULL),
hESDMotherPi0PtY(NULL),
hESDMotherEtaPtY(NULL),
hESDMotherPi0PtAlpha(NULL),
hESDMotherEtaPtAlpha(NULL),
hESDMotherPi0PtOpenAngle(NULL),
hESDMotherEtaPtOpenAngle(NULL),
hNEvents(NULL),
hNGoodESDTracks(NULL),
hNGammaCandidates(NULL),
hNGoodESDTracksVsNGammaCanditates(NULL),
hNV0Tracks(NULL),
hEtaShift(NULL),
tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
fInvMass(0),
fPt(0),
fDCAzGammaMin(0),
fDCAzGammaMax(0),
iFlag(0),
iMesonMCInfo(0),
fEventPlaneAngle(-100),
fRandom(0),
fnGammaCandidates(0),
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
fDoMesonAnalysis(kTRUE),
fDoMesonQA(0),
fDoPhotonQA(0),
fIsFromMBHeader(kTRUE),
fhistoEPVZ(NULL),
fMinMass(-1),
fMaxMass(10),
fMinKappa(-1),
fMaxKappa(100),
fFilterVariable(1),
fMinFilter(0),
fMaxFilter(0.2),
fApplydPhidRCut(0),
fPerformExtraStudies(0),
fDebug(0),
fCutsRP(0), 
fNullCuts(0), 
fFlowEvent(NULL),
fIsMC(0),
fMCEvent(NULL)

{
  // Define output slots here
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliFlowEventSimple::Class());
  DefineOutput(3, AliFlowEventSimple::Class());
  DefineOutput(4, AliFlowEventSimple::Class());
  DefineOutput(5, AliFlowEventSimple::Class());
  DefineOutput(6, AliFlowEventSimple::Class());

}


//________________________________________________________________________
AliAnalysisTaskGammaConvFlow::AliAnalysisTaskGammaConvFlow(const char *name, Int_t nCuts):
AliAnalysisTaskSE(name),
fV0Reader(NULL),
fV0ReaderName("V0ReaderV1"),
fBGHandler(NULL),
fBGHandlerRP(NULL),
fInputEvent(NULL),
fCutFolder(NULL),
fESDList(NULL),
fBackList(NULL),
fMotherList(NULL),
fPhotonDCAList(NULL),
fMesonDCAList(NULL),
//fTrueList(NULL),
//fMCList(NULL),
fHeaderNameList(NULL),
fOutputContainer(0),
fReaderGammas(NULL),
fGammaCandidates(NULL),
fEventCutArray(NULL),
fEventCuts(NULL),
fCutArray(NULL),
fConversionCuts(NULL),
fMesonCutArray(NULL),
fMesonCuts(NULL),
hESDConvGammaPt(NULL),
hInvMassPair(NULL),
hLTMPt(NULL),
hLTMPt_MC(NULL),
hPt_TruePt(NULL),
hdPhidRcandidates(NULL),
hdPhidRcandidates_MCsigsig(NULL),
hdPhidRcandidates_MCbkgsig(NULL),
hdPhidRcandidates_MCbkgbkg(NULL),
hKappaTPC(NULL),
hKappaTPC_after(NULL),
hKappaTPC_Temp0(NULL),
hKappaTPC_Temp1(NULL),
hKappaTPC_Temp2(NULL),
hKappaTPC_Temp3(NULL),
hKappaTPC_Temp4(NULL),
hKappaTPC_Temp5(NULL),
hKappaTPC_Temp6(NULL),
hKappaTPC_Temp7(NULL),
hKappaTPC_Temp8(NULL),
hKappaTPC_Temp9(NULL),
hKappaTPC_Temp10(NULL),
hESDConvGammaR(NULL),
hESDConvGammaEta(NULL),
//tESDConvGammaPtDcazCat(NULL),
fPtGamma(0),
fDCAzPhoton(0),
fRConvPhoton(0),
fEtaPhoton(0),
iCatPhoton(0),
iPhotonMCInfo(0),
hESDMotherInvMassPt(NULL),
//sESDMotherInvMassPtZM(NULL),
hESDMotherBackInvMassPt(NULL),
//sESDMotherBackInvMassPtZM(NULL),
hESDMotherInvMassEalpha(NULL),
hESDMotherPi0PtY(NULL),
hESDMotherEtaPtY(NULL),
hESDMotherPi0PtAlpha(NULL),
hESDMotherEtaPtAlpha(NULL),
hESDMotherPi0PtOpenAngle(NULL),
hESDMotherEtaPtOpenAngle(NULL),
hNEvents(NULL),
hNGoodESDTracks(NULL),
hNGammaCandidates(NULL),
hNGoodESDTracksVsNGammaCanditates(NULL),
hNV0Tracks(NULL),
hEtaShift(NULL),
tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
fInvMass(0),
fPt(0),
fDCAzGammaMin(0),
fDCAzGammaMax(0),
iFlag(0),
iMesonMCInfo(0),
fEventPlaneAngle(-100),
fRandom(0),
fnGammaCandidates(0),
fUnsmearedPx(NULL),
fUnsmearedPy(NULL),
fUnsmearedPz(NULL),
fUnsmearedE(NULL),
fMCEventPos(NULL),
fMCEventNeg(NULL),
fESDArrayPos(NULL),
fESDArrayNeg(NULL),
fnCuts(nCuts),
fiCut(0),
fMoveParticleAccordingToVertex(kTRUE),
fIsHeavyIon(0),
fDoMesonAnalysis(kTRUE),
fDoMesonQA(0),
fDoPhotonQA(0),
fIsFromMBHeader(kTRUE),
fhistoEPVZ(NULL),
fMinMass(-1),
fMaxMass(10),
fMinKappa(-1),
fMaxKappa(100),
fFilterVariable(1),
fMinFilter(0),
fMaxFilter(0.2),
fApplydPhidRCut(0),
fPerformExtraStudies(0),
fDebug(0),
fCutsRP(0), 
fNullCuts(0), 
fFlowEvent(NULL),
fIsMC(0),
fMCEvent(NULL)

{
  // Define output slots here
  DefineOutput(1, TList::Class());
  for (Int_t i = 0; i<nCuts; i++){
    DefineOutput(2+i, AliFlowEventSimple::Class());
  }

}
//___________________________________________________________
AliAnalysisTaskGammaConvFlow::~AliAnalysisTaskGammaConvFlow()
{
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
  }
  if(fBGHandler){
    delete[] fBGHandler;
    fBGHandler = 0x0;
  }
  if(fBGHandlerRP){
    delete[] fBGHandlerRP;
    fBGHandlerRP = 0x0;
  }
  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (fFlowEvent[iCut]) delete fFlowEvent[iCut];
  }

}

//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::UserCreateOutputObjects(){
  
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  //========================= again flow setting==========================
  //Create histograms
  //----------hfe initialising begin---------
  fNullCuts = new AliFlowTrackCuts("null_cuts");
  
  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(10000);
  cc->SetMultMin(0);
  cc->SetMultMax(10000);
  
  cc->SetNbinsPt(200);
  cc->SetPtMin(0);
  cc->SetPtMax(20);
  
  cc->SetNbinsPhi(180);
  cc->SetPhiMin(0.0);
  cc->SetPhiMax(TMath::TwoPi());
  
  cc->SetNbinsEta(40);
  cc->SetEtaMin(-2.0);
  cc->SetEtaMax(+2.0);
  
  cc->SetNbinsQ(500);
  cc->SetQMin(0.0);
  cc->SetQMax(3.0);
  
  cc->SetMassMin(fMinFilter);
  cc->SetMassMax(fMaxFilter);
  cc->SetNbinsMass(100);

  // Array of current cut's gammas
  fGammaCandidates = new TList();
  
  fCutFolder = new TList*[fnCuts];
  fESDList = new TList*[fnCuts];
  fBackList = new TList*[fnCuts];
  fMotherList = new TList*[fnCuts];
  hNEvents = new TH1I*[fnCuts];
  hNGoodESDTracks = new TH1I*[fnCuts];
  hNGammaCandidates = new TH1I*[fnCuts];
  hNGoodESDTracksVsNGammaCanditates = new TH2F*[fnCuts];
  hNV0Tracks = new TH1I*[fnCuts];
  hEtaShift = new TProfile*[fnCuts];
  hESDConvGammaPt = new TH1F*[fnCuts];
  hInvMassPair = new TH2F*[fnCuts];
  hLTMPt = new TH2F*[fnCuts];
  hLTMPt_MC = new TH2F*[fnCuts];
  hPt_TruePt = new TH2F*[fnCuts];
  hdPhidRcandidates = new TH2F*[fnCuts];
  hdPhidRcandidates_MCsigsig = new TH2F*[fnCuts];
  hdPhidRcandidates_MCbkgsig = new TH2F*[fnCuts];
  hdPhidRcandidates_MCbkgbkg = new TH2F*[fnCuts];
  hKappaTPC = new TH2F*[fnCuts];
  hKappaTPC_after = new TH2F*[fnCuts];
  hKappaTPC_Temp0 = new TH2F*[fnCuts];
  hKappaTPC_Temp1 = new TH2F*[fnCuts];
  hKappaTPC_Temp2 = new TH2F*[fnCuts];
  hKappaTPC_Temp3 = new TH2F*[fnCuts];
  hKappaTPC_Temp4 = new TH2F*[fnCuts];
  hKappaTPC_Temp5 = new TH2F*[fnCuts];
  hKappaTPC_Temp6 = new TH2F*[fnCuts];
  hKappaTPC_Temp7 = new TH2F*[fnCuts];
  hKappaTPC_Temp8 = new TH2F*[fnCuts];
  hKappaTPC_Temp9 = new TH2F*[fnCuts];
  hKappaTPC_Temp10 = new TH2F*[fnCuts];
  
  if (fDoPhotonQA == 2){
    fPhotonDCAList = new TList*[fnCuts];
//      tESDConvGammaPtDcazCat = new TTree*[fnCuts];
  }
  if (fDoPhotonQA > 0){
    hESDConvGammaR = new TH1F*[fnCuts];
    hESDConvGammaEta = new TH1F*[fnCuts];
  }
  
  if(fDoMesonAnalysis){
    hESDMotherInvMassPt = new TH2F*[fnCuts];
    hESDMotherBackInvMassPt = new TH2F*[fnCuts];
    hESDMotherInvMassEalpha = new TH2F*[fnCuts];
    if (fDoMesonQA == 2){
      fMesonDCAList = new TList*[fnCuts];
      tESDMesonsInvMassPtDcazMinDcazMaxFlag = new TTree*[fnCuts];
    }
    if (fDoMesonQA > 0){
      hESDMotherPi0PtY =  new TH2F*[fnCuts];
      hESDMotherEtaPtY =  new TH2F*[fnCuts];
      hESDMotherPi0PtAlpha =  new TH2F*[fnCuts];
      hESDMotherEtaPtAlpha =  new TH2F*[fnCuts];
      hESDMotherPi0PtOpenAngle =  new TH2F*[fnCuts];
      hESDMotherEtaPtOpenAngle =  new TH2F*[fnCuts];
    }
  }
  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    
    TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPhoton = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson = "NoMesonCut";
    if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
    fCutFolder[iCut] = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s",cutstringEvent.Data() ,cutstringPhoton.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut] = new TList();
    fESDList[iCut]->SetName(Form("%s_%s ESD histograms",cutstringEvent.Data() ,cutstringPhoton.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);
    
    hNEvents[iCut] = new TH1I("NEvents","NEvents",13,-0.5,12.5);
    hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
      TString TriggerNames = "Not Trigger: ";
      TriggerNames = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
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
    hNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fESDList[iCut]->Add(hNEvents[iCut]);
    
    if(fIsHeavyIon == 1) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",4000,0,4000);
    else if(fIsHeavyIon == 2) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",400,0,400);
    else hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
    fESDList[iCut]->Add(hNGoodESDTracks[iCut]);
    if(fIsHeavyIon == 1) hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",100,0,100);
    else if(fIsHeavyIon == 2) hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
    else hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
    fESDList[iCut]->Add(hNGammaCandidates[iCut]);
    if(fIsHeavyIon == 1) hNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,100,0,100);
    else if(fIsHeavyIon == 2) hNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,50,0,50);
    else hNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,50,0,50);
    fESDList[iCut]->Add(hNGoodESDTracksVsNGammaCanditates[iCut]);
    
    
    if(fIsHeavyIon == 1) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
    else if(fIsHeavyIon == 2) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
    else hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
    fESDList[iCut]->Add(hNV0Tracks[iCut]);
    hEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
    fESDList[iCut]->Add(hEtaShift[iCut]);
    hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
    fESDList[iCut]->Add(hESDConvGammaPt[iCut]);
    
    hInvMassPair[iCut]= new TH2F("InvMassPair_Pt","Gamma invariant mass vs Pt",200,0,0.2,250,0,25);	
    fESDList[iCut]->Add(hInvMassPair[iCut]);
    
    hLTMPt[iCut]= new TH2F("LTM_Pt","LTM vs Pt",200,0,200,250,0,25); 
    fESDList[iCut]->Add(hLTMPt[iCut]);
    hLTMPt_MC[iCut]= new TH2F("LTM_Pt_MCgen","LTM vs Pt (MC)",200,0,200,250,0,25); 
    fESDList[iCut]->Add(hLTMPt_MC[iCut]);
    
    hPt_TruePt[iCut]= new TH2F("hPt_TruePt","Pt vs TruePt (MC)",250,0,25,250,0,25);
    fESDList[iCut]->Add(hPt_TruePt[iCut]);
    hdPhidRcandidates[iCut]= new TH2F("hdPhidRcandidates","dPhi vs dRconvVtx",400,0,4,300,0,300); 
    fESDList[iCut]->Add(hdPhidRcandidates[iCut]);
    hdPhidRcandidates_MCsigsig[iCut]= new TH2F("hdPhidRcandidates_MCsigsig","dPhi vs dRconvVtx",400,0,4,300,0,300); 
    fESDList[iCut]->Add(hdPhidRcandidates_MCsigsig[iCut]);
    hdPhidRcandidates_MCbkgsig[iCut]= new TH2F("hdPhidRcandidates_MCbkgsig","dPhi vs dRconvVtx",400,0,4,300,0,300); 
    fESDList[iCut]->Add(hdPhidRcandidates_MCbkgsig[iCut]);
    hdPhidRcandidates_MCbkgbkg[iCut]= new TH2F("hdPhidRcandidates_MCbkgbkg","dPhi vs dRconvVtx",400,0,4,300,0,300); 
    fESDList[iCut]->Add(hdPhidRcandidates_MCbkgbkg[iCut]);
    
    
    
    hKappaTPC[iCut]= new TH2F("KappaTPC_Pt","Gamma KappaTPC vs Pt",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC[iCut]);
    
    hKappaTPC_after[iCut]= new TH2F("KappaTPC_Pt_after","Gamma KappaTPC vs Pt after cuts",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_after[iCut]);
    
    hKappaTPC_Temp0[iCut]= new TH2F("hKappaTPC_Temp0_ee","Gamma KappaTPC vs Pt Template 0",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp0[iCut]);
    
    hKappaTPC_Temp1[iCut]= new TH2F("hKappaTPC_Temp1_pipi","Gamma KappaTPC vs Pt Template 1",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp1[iCut]);
    
    hKappaTPC_Temp2[iCut]= new TH2F("hKappaTPC_Temp2_pie","Gamma KappaTPC vs Pt Template 2",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp2[iCut]);
    
    hKappaTPC_Temp3[iCut]= new TH2F("hKappaTPC_Temp3_piK","Gamma KappaTPC vs Pt Template 3",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp3[iCut]);
    
    hKappaTPC_Temp4[iCut]= new TH2F("hKappaTPC_Temp4_pip","Gamma KappaTPC vs Pt Template 4",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp4[iCut]);
    
    hKappaTPC_Temp5[iCut]= new TH2F("hKappaTPC_Temp5_eK","Gamma KappaTPC vs Pt Template 5",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp5[iCut]);
    
    hKappaTPC_Temp6[iCut]= new TH2F("hKappaTPC_Temp6_ep","Gamma KappaTPC vs Pt Template 6",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp6[iCut]);
    
    hKappaTPC_Temp7[iCut]= new TH2F("hKappaTPC_Temp7_KK","Gamma KappaTPC vs Pt Template 7",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp7[iCut]);
    
    hKappaTPC_Temp8[iCut]= new TH2F("hKappaTPC_Temp8_had","Gamma KappaTPC vs Pt Template 8",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp8[iCut]);
    
    hKappaTPC_Temp9[iCut]= new TH2F("hKappaTPC_Temp9_rem4","Gamma KappaTPC vs Pt Template 9",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp9[iCut]);
    
    hKappaTPC_Temp10[iCut]= new TH2F("hKappaTPC_Temp10_rem10","Gamma KappaTPC vs Pt Template 10",200,-20,20,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_Temp10[iCut]);    
    
    //2d histogram filling the cut and value - control check for selections
    
    if (fDoPhotonQA == 2){
      fPhotonDCAList[iCut] = new TList();
      fPhotonDCAList[iCut]->SetName(Form("%s_%s_%s Photon DCA tree",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringMeson.Data()));
      fPhotonDCAList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fPhotonDCAList[iCut]);
      
//            tESDConvGammaPtDcazCat[iCut] = new TTree("ESD_ConvGamma_Pt_Dcaz_R_Eta","ESD_ConvGamma_Pt_Dcaz_R_Eta_Cat");
//            tESDConvGammaPtDcazCat[iCut]->Branch("Pt",&fPtGamma,"fPtGamma/F");
//            tESDConvGammaPtDcazCat[iCut]->Branch("DcaZPhoton",&fDCAzPhoton,"fDCAzPhoton/F");
      //          tESDConvGammaPtDcazCat[iCut]->Branch("R",&fRConvPhoton,"fRConvPhoton/F");
      //          tESDConvGammaPtDcazCat[iCut]->Branch("Eta",&fEtaPhoton,"fEtaPhoton/F");
      
    //    tESDConvGammaPtDcazCat[iCut]->Branch("cat",&iCatPhoton,"iCatPhoton/b");

  //        fPhotonDCAList[iCut]->Add(tESDConvGammaPtDcazCat[iCut]);
    }
    
    if (fDoPhotonQA > 0){
      hESDConvGammaR[iCut] = new TH1F("ESD_ConvGamma_R","ESD_ConvGamma_R",800,0,200);
      fESDList[iCut]->Add(hESDConvGammaR[iCut]);
      hESDConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",2000,-2,2);
      fESDList[iCut]->Add(hESDConvGammaEta[iCut]);
    }
    
    if(fDoMesonAnalysis){
      hESDMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,250,0,25);
      fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);
      hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,250,0,25);
      fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);
      hESDMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",800,0,0.8,250,0,25);
      fESDList[iCut]->Add(hESDMotherInvMassEalpha[iCut]);
      if (fDoMesonQA == 2){
        fMesonDCAList[iCut] = new TList();
        fMesonDCAList[iCut]->SetName(Form("%s_%s_%s Meson DCA tree",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
        fMesonDCAList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMesonDCAList[iCut]);
        
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut] = new TTree("ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag","ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("InvMass",&fInvMass,"fInvMass/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("Pt",&fPt,"fPt/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMin",&fDCAzGammaMin,"fDCAzGammaMin/F");
        tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMax",&fDCAzGammaMax,"fDCAzGammaMax/F");
        
        fMesonDCAList[iCut]->Add(tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]);
        
      }
      if (fDoMesonQA > 0 ){
        hESDMotherPi0PtY[iCut] = new TH2F("ESD_MotherPi0_Pt_Y","ESD_MotherPi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
        SetLogBinningXTH2(hESDMotherPi0PtY[iCut]);
        fESDList[iCut]->Add(hESDMotherPi0PtY[iCut]);
        hESDMotherEtaPtY[iCut] = new TH2F("ESD_MotherEta_Pt_Y","ESD_MotherEta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
        SetLogBinningXTH2(hESDMotherEtaPtY[iCut]);
        fESDList[iCut]->Add(hESDMotherEtaPtY[iCut]);
        hESDMotherPi0PtAlpha[iCut] = new TH2F("ESD_MotherPi0_Pt_Alpha","ESD_MotherPi0_Pt_Alpha",150,0.03,15.,100,0,1);
        SetLogBinningXTH2(hESDMotherPi0PtAlpha[iCut]);
        fESDList[iCut]->Add(hESDMotherPi0PtAlpha[iCut]);
        hESDMotherEtaPtAlpha[iCut] = new TH2F("ESD_MotherEta_Pt_Alpha","ESD_MotherEta_Pt_Alpha",150,0.03,15.,100,0,1);
        SetLogBinningXTH2(hESDMotherEtaPtAlpha[iCut]);
        fESDList[iCut]->Add(hESDMotherEtaPtAlpha[iCut]);
        hESDMotherPi0PtOpenAngle[iCut] = new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",150,0.03,15.,100,0,TMath::Pi());
        SetLogBinningXTH2(hESDMotherPi0PtOpenAngle[iCut]);
        fESDList[iCut]->Add(hESDMotherPi0PtOpenAngle[iCut]);
        hESDMotherEtaPtOpenAngle[iCut] = new TH2F("ESD_MotherEta_Pt_OpenAngle","ESD_MotherEta_Pt_OpenAngle",150,0.03,15.,100,0,TMath::Pi());
        SetLogBinningXTH2(hESDMotherEtaPtOpenAngle[iCut]);
        fESDList[iCut]->Add(hESDMotherEtaPtOpenAngle[iCut]);
      }
      
      
    }
    
    
  }
//     if(fDoMesonAnalysis){
//         InitBack(); // Init Background Handler
//     }
  

  
    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
  
  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
    if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
    }
    if(fDoMesonAnalysis){
      if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
      }
    }
  }
  
  fhistoEPVZ = new TH1D("EPVZ", "EPVZ", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputContainer->Add(fhistoEPVZ);


  
  PostData(1, fOutputContainer);
  
  fFlowEvent = new AliFlowEvent*[fnCuts];
  
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      fFlowEvent[iCut] = new AliFlowEvent(10000);
      if (fIsHeavyIon == 1){
        PostData(2+iCut, fFlowEvent[iCut]);
      }   
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvFlow::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
        
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
      hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      continue; // No Eta Shift requested, continue
    }
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
      hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      continue;
    }
    else{
      printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
        (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
      hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
    }
  }
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      hNEvents[iCut]->Fill(eventQuality);
    }
    return;
  }
  
  fInputEvent = InputEvent();
  if(fIsMC) fMCEvent = MCEvent();

  fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
  
  // ------------------- BeginEvent ----------------------------
  
  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;
  
  SetNullCuts(fInputEvent);
  
  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
//     cout << "event cut array =  " << fEventCutArray->At(iCut) << endl;
    PrepareFlowEvent(fInputEvent->GetNumberOfTracks(),fFlowEvent[iCut]);    //Calculate event plane Qvector and EP resolution for inclusive
// 		Int_t iCut = 0;
    fiCut = iCut;
    Int_t eventNotAccepted =
    ((AliConvEventCuts*)fEventCutArray->At(iCut))
    ->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
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
    hNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks());
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)	hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A());
    else hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());

    ProcessPhotonCandidates(); // Process this cuts gammas
    ProcessPhotonCandidatesforV2(); // Process this cuts gammas and do v2
    
    if(fPerformExtraStudies){
      ProcessPhotonCandidatesforLTM(); // calculate the Local Track Multiplicity in a jet cone around the candidates
      GetdPhidRtoCandidate(); // calculate the distances to other conversions in the event for each selected candidate
    }

    hNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries());
    hNGoodESDTracksVsNGammaCanditates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fGammaCandidates->GetEntries());
    if (fIsHeavyIon == 1){
        PostData(2+iCut, fFlowEvent[iCut]);
    }    
    fGammaCandidates->Clear(); // delete this cuts good gammas
  }
  
  fhistoEPVZ->Fill(fEventPlaneAngle);

  PostData(1, fOutputContainer);
  

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::ProcessPhotonCandidates()
{
    Int_t nV0 = 0;
    TList *GammaCandidatesStepOne = new TList();
    TList *GammaCandidatesStepTwo = new TList();
    Int_t PhotonTemplateID;
    // Loop over Photon Candidates allocated by ReaderV1
    for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
        AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
        if(!PhotonCandidate) continue;
        fIsFromMBHeader = kTRUE;
        
        hKappaTPC[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
        if (PhotonCandidate->GetInvMassPair() < fMinMass || PhotonCandidate->GetInvMassPair() > fMaxMass) continue;
        if (((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent) < fMinKappa || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent) > fMaxKappa) continue;
        if(fApplydPhidRCut==1){ if(GetdPhidRtoCandidate(PhotonCandidate,i)==1) continue; }
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
          !((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
            fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
//             cout << "IS IT MC? " << fIsMC << endl;
            if(fIsMC){
              PhotonTemplateID = GetTemplateID(PhotonCandidate);
              if(PhotonTemplateID == 0) hKappaTPC_Temp0[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 1) hKappaTPC_Temp1[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 2) hKappaTPC_Temp2[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 3) hKappaTPC_Temp3[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 4) hKappaTPC_Temp4[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 5) hKappaTPC_Temp5[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 6) hKappaTPC_Temp6[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 7) hKappaTPC_Temp7[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID == 8) hKappaTPC_Temp8[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID != 0 && PhotonTemplateID != 1 && PhotonTemplateID != 2) hKappaTPC_Temp9[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
              if(PhotonTemplateID != 0 && PhotonTemplateID != 1 && PhotonTemplateID != 2 &&
                PhotonTemplateID != 3 && PhotonTemplateID != 4 && PhotonTemplateID != 5 &&
                PhotonTemplateID != 6 && PhotonTemplateID != 7 && PhotonTemplateID != 8) hKappaTPC_Temp10[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());

              TParticle *TRUEPhoton = PhotonCandidate->GetMCParticle(fMCEvent);
              if(TRUEPhoton) hPt_TruePt[fiCut]->Fill(PhotonCandidate->Pt(),TRUEPhoton->Pt());
            }
            
            if(fIsFromMBHeader){
                hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
                hInvMassPair[fiCut]->Fill(PhotonCandidate->GetInvMassPair(),PhotonCandidate->Pt());
                hKappaTPC_after[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
                if (fDoPhotonQA > 0){
                    hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
                    hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
                }
            }

            
            if (fIsFromMBHeader && fDoPhotonQA == 2){
                if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                  // tESDConvGammaPtDcazCat[fiCut]->Fill();
                } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                  //  tESDConvGammaPtDcazCat[fiCut]->Fill();
                }
            }
        } else if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
            nV0++;
            GammaCandidatesStepOne->Add(PhotonCandidate);
        } else if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
            GammaCandidatesStepTwo->Add(PhotonCandidate);
        }
    }
    
    
    
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
        for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
            AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
            if(!PhotonCandidate) continue;
            fIsFromMBHeader = kTRUE;

            if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
            if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
                fGammaCandidates->Add(PhotonCandidate);
                if(fIsFromMBHeader){
                    hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
                    hInvMassPair[fiCut]->Fill(PhotonCandidate->GetInvMassPair(),PhotonCandidate->Pt());
                    hKappaTPC_after[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
                    if (fDoPhotonQA > 0){
                        hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
                        hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
                    }
                }
            }

                
                GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
            
            if (fIsFromMBHeader && fDoPhotonQA == 2){
                if (fIsHeavyIon ==1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                //    tESDConvGammaPtDcazCat[fiCut]->Fill();
                } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                //    tESDConvGammaPtDcazCat[fiCut]->Fill();
                }
            }
        }
    }
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
        for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
            AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
            if(!PhotonCandidate) continue;
            fIsFromMBHeader = kTRUE;

            if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
            fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
            if(fIsFromMBHeader){
                hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
                hInvMassPair[fiCut]->Fill(PhotonCandidate->GetInvMassPair(),PhotonCandidate->Pt());
                hKappaTPC_after[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
                if (fDoPhotonQA > 0){
                    hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
                    hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
                }
            }

            if (fIsFromMBHeader){
                if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
              //     tESDConvGammaPtDcazCat[fiCut]->Fill();
                } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                  //  tESDConvGammaPtDcazCat[fiCut]->Fill();
                }
            }
        }
    }
    
    delete GammaCandidatesStepOne;
    GammaCandidatesStepOne = 0x0;
    delete GammaCandidatesStepTwo;
    GammaCandidatesStepTwo = 0x0;
    
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::UpdateEventByEventData(){
    //see header file for documentation
    if(fGammaCandidates->GetEntries() >0 ){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
            fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
        }
        else{ // means we use #V0s for multiplicity
            fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
        }
    }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskGammaConvFlow::Terminate(const Option_t *)
{
    
    //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::ProcessPhotonCandidatesforV2()
{

  // Loop over Photon Candidates allocated by ReaderV1	
//   cout << "number of gamma's: " << fGammaCandidates->GetEntries() << endl;
  for(Int_t i = 0; i < fGammaCandidates->GetEntries(); i++){
    
    AliAODConversionPhoton *gammaForv2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(i));
    if (gammaForv2 == NULL) return;
    AliFlowTrack *sTrack = new AliFlowTrack();
    sTrack->SetForRPSelection(kFALSE);
    sTrack->SetForPOISelection(kTRUE);
    
    if(fFilterVariable==1){//using mass for POI selections
      sTrack->SetMass(gammaForv2->GetInvMassPair());
    }else if(fFilterVariable==2){//using Kappa for POI selections
      sTrack->SetMass(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(gammaForv2,fInputEvent));
    }else if(fFilterVariable==3 && fIsMC){//MC ElectronElectron + mass
      if(!MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(gammaForv2->GetInvMassPair());
    }else if(fFilterVariable==4 && fIsMC){//MC ElectronElectron + kappa
      if(!MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(gammaForv2,fInputEvent));
    }else if(fFilterVariable==5 && fIsMC){//MC Not ElectronElectron + mass
      if(MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(gammaForv2->GetInvMassPair());
    }else if(fFilterVariable==6 && fIsMC){//MC Not ElectronElectron + kappa
      if(MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(gammaForv2,fInputEvent));
    }else{//no additional POI selection
      sTrack->SetMass(123456);
    }
    sTrack->SetPt(gammaForv2->Pt());
    sTrack->SetPhi(gammaForv2->GetPhotonPhi());
    sTrack->SetEta(gammaForv2->GetPhotonEta());

  /*     for(int iRPs=0; iRPs!=fFlowEvent->NumberOfTracks(); ++iRPs)
    {
      //   cout << " no of rps " << iRPs << endl;
      AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack( iRPs ));
      if (!iRP) continue;
      if (!iRP->InRPSelection()) continue;
      if( sTrack->GetID() == iRP->GetID())
      {
        if(fDebug) printf(" was in RP set");
        //  cout << sTrack->GetID() <<"   ==  " << iRP->GetID() << " was in RP set====REMOVED" <<endl;
        iRP->SetForRPSelection(kFALSE);
        // fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
      }
    } //end of for loop on RPs*/
    fFlowEvent[fiCut]->InsertTrack(((AliFlowTrack*) sTrack));
    fFlowEvent[fiCut]->SetNumberOfPOIs(fFlowEvent[fiCut]->GetNumberOfPOIs()+1);
//     cout << "cutnumber: " << fiCut << " nPoi " << fFlowEvent[fiCut]->GetNumberOfPOIs() << " ntracks " << fFlowEvent[fiCut]->NumberOfTracks() << endl;
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::ProcessPhotonCandidatesforLTM()
{
  
  Float_t gamma_Eta = 0;
  Float_t gamma_Phi = 0;
  Float_t gamma_Pt = 0;
  
  Float_t LTMpart_Eta = 0;
  Float_t LTMpart_Phi = 0;
  Float_t LTMpart_Pt = 0;
  
  Float_t dPhi = 0;
  
  Float_t nCloseByTracks = 0;
  
  if(fIsMC){
    // Loop over Photon Candidates allocated by MCstack
    for(Int_t i = 0; i < fMCEvent->GetNumberOfPrimaries(); i++){
      TParticle* gammaForLTM = (TParticle *)fMCEvent->Particle(i);
      if(!gammaForLTM) return;
      if(!MCConversionPhotonCheck(gammaForLTM)) continue;
      gamma_Eta = gammaForLTM->Eta(); gamma_Phi = gammaForLTM->Phi(); gamma_Pt = gammaForLTM->Pt();
      if( gamma_Eta > 0.9 || gamma_Eta < -0.9 ) continue;
      if(gamma_Eta==0 || gamma_Phi==0 || gamma_Pt==0)continue;
      nCloseByTracks = 0;
      for(Int_t j = 0; j < fInputEvent->GetNumberOfTracks(); j++){
        AliVParticle *LTMpart = fInputEvent->GetTrack(j);
        if (LTMpart == NULL) return;
        LTMpart_Eta = LTMpart->Eta(); LTMpart_Phi = LTMpart->Phi(); LTMpart_Pt = LTMpart->Pt();
        if(LTMpart_Eta==0 || LTMpart_Phi==0 || LTMpart_Pt==0)continue;
        dPhi = TMath::Abs(LTMpart_Phi-gamma_Phi);
        if(dPhi > TMath::Pi()) dPhi = TMath::Abs(dPhi-2.0*TMath::Pi());
        if(TMath::Sqrt(pow((LTMpart_Eta-gamma_Eta),2)+pow(dPhi,2))<0.2) nCloseByTracks+=1;
      }
//       cout << "nCloseByTracks MCgen= " << nCloseByTracks << " with pt= " << gamma_Pt <<  endl;
      hLTMPt_MC[fiCut]->Fill(nCloseByTracks,gamma_Pt);
    }
  }
  // Loop over Photon Candidates allocated by ReaderV1
  for(Int_t i = 0; i < fGammaCandidates->GetEntries(); i++){
    AliAODConversionPhoton *gammaForLTM=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(i));
    if(!gammaForLTM) return;
    gamma_Eta = gammaForLTM->GetPhotonEta(); gamma_Phi = gammaForLTM->GetPhotonPhi(); gamma_Pt = gammaForLTM->GetPhotonPt();
    if(gamma_Eta==0 || gamma_Phi==0 || gamma_Pt==0)continue;
    nCloseByTracks = 0;
    for(Int_t j = 0; j < fInputEvent->GetNumberOfTracks(); j++){
      AliVParticle *LTMpart = fInputEvent->GetTrack(j);
      if (LTMpart == NULL) return;
      LTMpart_Eta = LTMpart->Eta(); LTMpart_Phi = LTMpart->Phi(); LTMpart_Pt = LTMpart->Pt();
      if(LTMpart_Eta==0 || LTMpart_Phi==0 || LTMpart_Pt==0)continue;
      dPhi = TMath::Abs(LTMpart_Phi-gamma_Phi);
      if(dPhi > TMath::Pi()) dPhi = TMath::Abs(dPhi-2.0*TMath::Pi());
      if(TMath::Sqrt(pow((LTMpart_Eta-gamma_Eta),2)+pow(dPhi,2))<0.2) nCloseByTracks+=1;
    }
//     cout << "nCloseByTracks= " << nCloseByTracks << " with pt= " << gamma_Pt <<  endl;
    hLTMPt[fiCut]->Fill(nCloseByTracks,gamma_Pt);
  }
}

//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskGammaConvFlow::SetNullCuts(T* event)
{
  // Set null cuts
  if (fDebug) cout << " fCutsRP " << fCutsRP << endl;
  fCutsRP->SetEvent(event, MCEvent());
  fNullCuts->SetParamType(AliFlowTrackCuts::kGlobal);
  fNullCuts->SetPtRange(+1, -1); // select nothing QUICK
  fNullCuts->SetEtaRange(+1, -1); // select nothing VZERO
  fNullCuts->SetEvent(event, MCEvent());
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const
{
  //Prepare flow events
    FlowEv->ClearFast();
    FlowEv->Fill(fCutsRP, fNullCuts);
    FlowEv->SetReferenceMultiplicity(iMulti);
    FlowEv->DefineDeadZone(0, 0, 0, 0);
    //  FlowEv->TagSubeventsInEta(-0.7, 0, 0, 0.7);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvFlow::MCGammaSignal( AliAODConversionPhoton *MCPhoton ){
  
  
  TParticle *posDaughter = MCPhoton->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = MCPhoton->GetNegativeMCDaughter(fMCEvent);
  if(posDaughter==NULL || negDaughter==NULL) return kFALSE;
  Int_t pdgCodePos = 0; 
  Int_t pdgCodeNeg = 0;
  Bool_t IsItGammaSignal = kFALSE;
  
  if( (posDaughter->GetMother(0) == negDaughter->GetMother(0)) ) {
    pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
    pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
    if(pdgCodePos==11 && pdgCodeNeg==11) IsItGammaSignal = kTRUE;
  }else{
    IsItGammaSignal = kFALSE;
  }
  
  return IsItGammaSignal;
    
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvFlow::MCElectronElectron( AliAODConversionPhoton *MCPhoton ){
  
  
  TParticle *posDaughter = MCPhoton->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = MCPhoton->GetNegativeMCDaughter(fMCEvent);
  if(posDaughter==NULL || negDaughter==NULL) return kFALSE;
  Int_t pdgCodePos = 0; 
  Int_t pdgCodeNeg = 0;
  Bool_t IsItElectronElectron = kFALSE;
  
  if( (posDaughter->GetMother(0) != negDaughter->GetMother(0))  || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1) ) {
    pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
    pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
    if(pdgCodePos==11 && pdgCodeNeg==11) IsItElectronElectron = kTRUE;
  }else{
    IsItElectronElectron = kFALSE;
  }
  
  return IsItElectronElectron;
    
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvFlow::MCConversionPhotonCheck( TParticle *MCPhoton ){
  
  if(MCPhoton->GetPdgCode()!=22) return kFALSE;
  if(MCPhoton->GetNDaughters() != 2) return kFALSE;
  TParticle *posDaughter = (TParticle*)fMCEvent->Particle(MCPhoton->GetFirstDaughter());
  TParticle *negDaughter = (TParticle*)fMCEvent->Particle(MCPhoton->GetLastDaughter());
  if(posDaughter==NULL || negDaughter==NULL) return kFALSE;
  if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() != 5) return kFALSE;
  Int_t pdgCodePos = 0; 
  Int_t pdgCodeNeg = 0;
  Bool_t IsItPhoton = kFALSE;
  
  if(posDaughter->GetMother(0) == negDaughter->GetMother(0)) {
    pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
    pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
    if(pdgCodePos==11 && pdgCodeNeg==11) IsItPhoton = kTRUE;
  }else{
    IsItPhoton = kFALSE;
  }
  
  return IsItPhoton;
    
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaConvFlow::GetTemplateID( AliAODConversionPhoton *MCPhoton ){
  
  
  TParticle *posDaughter = MCPhoton->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = MCPhoton->GetNegativeMCDaughter(fMCEvent);
  if(posDaughter==NULL || negDaughter==NULL) return kFALSE;
  
  Int_t pdgCodePos = 0; 
  Int_t pdgCodeNeg = 0;
  Int_t TemplateID = -1;
  
  pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
  pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
  
  if(pdgCodePos==11 && pdgCodeNeg==11)                                              TemplateID = 0; //signal -> e+ e-
  if(pdgCodePos==211 && pdgCodeNeg==211)                                            TemplateID = 1; //pi pi 211 211
  if((pdgCodePos==211 && pdgCodeNeg==11) || (pdgCodePos==11 && pdgCodeNeg==211))    TemplateID = 2; //pi e 211 11
  if((pdgCodePos==211 && pdgCodeNeg==321) || (pdgCodePos==321 && pdgCodeNeg==211))  TemplateID = 3; //pi K 211 321
  if((pdgCodePos==211 && pdgCodeNeg==2212) || (pdgCodePos==2212 && pdgCodeNeg==211))TemplateID = 4; //pi p 211 2212
  if((pdgCodePos==11 && pdgCodeNeg==321) || (pdgCodePos==321 && pdgCodeNeg==11))    TemplateID = 5; //e K 11 321
  if((pdgCodePos==11 && pdgCodeNeg==2212) || (pdgCodePos==2212 && pdgCodeNeg==11))  TemplateID = 6; //e p 11 2212
  if(pdgCodePos==321 && pdgCodeNeg==321)                                            TemplateID = 7; //K K 321 321
  if(pdgCodePos!=11 && pdgCodeNeg!=11 && TemplateID == -1)                          TemplateID = 8; //hadronic not 11 11
  
//   cout << "TEMPLATE ID IS: " << TemplateID << endl;
  return TemplateID;
    
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::GetdPhidRtoCandidate(){
  
  Float_t gamma1_Eta = 0; Float_t gamma1_Phi = 0; Float_t gamma1_Pt = 0;
  Float_t gamma2_Eta = 0; Float_t gamma2_Phi = 0; Float_t gamma2_Pt = 0;
  
  Float_t dPhi = 0; Float_t dRconvVtx = 0;
  
  Float_t gamma1_Vtx_x = 0; Float_t gamma1_Vtx_y = 0; Float_t gamma1_Vtx_z = 0;
  Float_t gamma2_Vtx_x = 0; Float_t gamma2_Vtx_y = 0; Float_t gamma2_Vtx_z = 0;
  
  for(Int_t i = 0; i < fGammaCandidates->GetEntries(); i++){
    AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(i));
    if(!gamma1) return;
    gamma1_Eta = gamma1->GetPhotonEta(); gamma1_Phi = gamma1->GetPhotonPhi(); gamma1_Pt = gamma1->GetPhotonPt();
    if(gamma1_Eta==0 || gamma1_Phi==0 || gamma1_Pt==0)continue;
    gamma1_Vtx_x = gamma1->GetConversionX(); gamma1_Vtx_y = gamma1->GetConversionY(); gamma1_Vtx_z = gamma1->GetConversionZ();
    for(Int_t j = i+1; j < fGammaCandidates->GetEntries(); j++){
      AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(j));
      if(!gamma2) continue;
      gamma2_Eta = gamma2->GetPhotonEta(); gamma2_Phi = gamma2->GetPhotonPhi(); gamma2_Pt = gamma2->GetPhotonPt();
      if(gamma2_Eta==0 || gamma2_Phi==0 || gamma2_Pt==0)continue;
      gamma2_Vtx_x = gamma2->GetConversionX(); gamma2_Vtx_y = gamma2->GetConversionY(); gamma2_Vtx_z = gamma2->GetConversionZ();
      dPhi = TMath::Abs(gamma2_Phi-gamma1_Phi);
      if(dPhi > TMath::Pi()) dPhi = TMath::Abs(dPhi-2.0*TMath::Pi());
      dRconvVtx = TMath::Sqrt(pow(gamma2_Vtx_x-gamma1_Vtx_x,2)+pow(gamma2_Vtx_y-gamma1_Vtx_y,2)+pow(gamma2_Vtx_z-gamma1_Vtx_z,2));
      hdPhidRcandidates[fiCut]->Fill(dPhi,dRconvVtx);
      if(fIsMC){
        // gamma1 = signal, gamma2 = signal  => sigsig
        if( MCGammaSignal(gamma1) && MCGammaSignal(gamma2) ) hdPhidRcandidates_MCsigsig[fiCut]->Fill(dPhi,dRconvVtx);
        // gamma1 = bkg, gamma2 = signal     => bkgsig
        if( (MCElectronElectron(gamma1) && MCGammaSignal(gamma2)) || (MCElectronElectron(gamma2) && MCGammaSignal(gamma1)) ) hdPhidRcandidates_MCbkgsig[fiCut]->Fill(dPhi,dRconvVtx);
        // gamma1 = bkg, gamma2 = bkg        => bkgbkg
        if( MCElectronElectron(gamma1) && MCElectronElectron(gamma2) ) hdPhidRcandidates_MCbkgbkg[fiCut]->Fill(dPhi,dRconvVtx);
      }
    }
  }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskGammaConvFlow::GetdPhidRtoCandidate( AliAODConversionPhoton* gamma, Int_t PhotonID ){
  
  Float_t gamma1_Eta = 0; Float_t gamma1_Phi = 0; Float_t gamma1_Pt = 0;
  Float_t gamma2_Eta = 0; Float_t gamma2_Phi = 0; Float_t gamma2_Pt = 0;
  
  Float_t dPhi = 0; Float_t dRconvVtx = 0;
  
  Float_t gamma1_Vtx_x = 0; Float_t gamma1_Vtx_y = 0; Float_t gamma1_Vtx_z = 0;
  Float_t gamma2_Vtx_x = 0; Float_t gamma2_Vtx_y = 0; Float_t gamma2_Vtx_z = 0;
  
  if(!gamma) return -1;
  gamma1_Eta = gamma->GetPhotonEta(); gamma1_Phi = gamma->GetPhotonPhi(); gamma1_Pt = gamma->GetPhotonPt();
  if(gamma1_Eta==0 || gamma1_Phi==0 || gamma1_Pt==0) return -1;
  gamma1_Vtx_x = gamma->GetConversionX(); gamma1_Vtx_y = gamma->GetConversionY(); gamma1_Vtx_z = gamma->GetConversionZ();
  for(Int_t i = 0; i < fReaderGammas->GetEntries(); i++){
    if(i==PhotonID) continue;
    AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fReaderGammas->At(i));
    if(!gamma2) continue;
    gamma2_Eta = gamma2->GetPhotonEta(); gamma2_Phi = gamma2->GetPhotonPhi(); gamma2_Pt = gamma2->GetPhotonPt();
    if(gamma2_Eta==0 || gamma2_Phi==0 || gamma2_Pt==0)continue;
    gamma2_Vtx_x = gamma2->GetConversionX(); gamma2_Vtx_y = gamma2->GetConversionY(); gamma2_Vtx_z = gamma2->GetConversionZ();
    dPhi = TMath::Abs(gamma2_Phi-gamma1_Phi);
    if(dPhi > TMath::Pi()) dPhi = TMath::Abs(dPhi-2.0*TMath::Pi());
    dRconvVtx = TMath::Sqrt(pow(gamma2_Vtx_x-gamma1_Vtx_x,2)+pow(gamma2_Vtx_y-gamma1_Vtx_y,2)+pow(gamma2_Vtx_z-gamma1_Vtx_z,2));
    if(dPhi < 0.04 && dRconvVtx < 7.5) return 1;
  }
  return 0;
}
