/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in thce supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**********************************
 * analysis task for CRC with ZDC *
 *                                *
 * author: Jacopo Margutti        *
 *         (margutti@nikhef.nl)   *
 **********************************/

#include "Riostream.h" //needed as include
#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TTimeStamp.h"
#include "TStopwatch.h"
#include "TProfile.h"
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TParticle.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskCRCZDC.h"
#include "AliMultSelection.h"

// ALICE Correction Framework
#include "AliCFManager.h"

// Interface to Event generators to get Reaction Plane Angle
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliGenEposEventHeader.h"

// Interface to Load short life particles
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

// Interface to make the Flow Event Simple used in the flow analysis methods
#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowCommonConstants.h"

ClassImp(AliAnalysisTaskCRCZDC)

//________________________________________________________________________
AliAnalysisTaskCRCZDC::AliAnalysisTaskCRCZDC():
AliAnalysisTaskSE(""),
fAnalysisType("AUTOMATIC"),
fRPType(""),
fCFManager1(NULL),
fCFManager2(NULL),
fCutsEvent(NULL),
fCutsRP(NULL),
fCutsPOI(NULL),
fCutContainer(new TList()),
fQAList(NULL),
fMinMult(0),
fMaxMult(10000000),
fMinA(-1.0),
fMaxA(-0.01),
fMinB(0.01),
fMaxB(1.0),
fGenHeader(NULL),
fPythiaGenHeader(NULL),
fHijingGenHeader(NULL),
fFlowTrack(NULL),
fAnalysisUtil(NULL),
fQAon(kFALSE),
fLoadCandidates(kFALSE),
fNbinsMult(10000),
fNbinsPt(100),
fNbinsPhi(100),
fNbinsEta(200),
fNbinsQ(500),
fNbinsMass(1),
fMultMin(0.),
fMultMax(10000.),
fPtMin(0.),
fPtMax(10.),
fPhiMin(0.),
fPhiMax(TMath::TwoPi()),
fEtaMin(-5.),
fEtaMax(5.),
fQMin(0.),
fQMax(3.),
fMassMin(-1.),
fMassMax(0.),
fHistWeightvsPhiMin(0.),
fHistWeightvsPhiMax(3.),
fExcludedEtaMin(0.),
fExcludedEtaMax(0.),
fExcludedPhiMin(0.),
fExcludedPhiMax(0.),
fAfterburnerOn(kFALSE),
fNonFlowNumberOfTrackClones(0),
fV1(0.),
fV2(0.),
fV3(0.),
fV4(0.),
fV5(0.),
fDifferentialV2(0),
fFlowEvent(NULL),
fShuffleTracks(kFALSE),
fMyTRandom3(NULL),
fAnalysisInput(kAOD),
fIsMCInput(kFALSE),
fUseMCCen(kTRUE),
fRejectPileUp(kTRUE),
fCentrLowLim(0.),
fCentrUpLim(100.),
fCentrEstimator("V0M"),
fOutput(0x0),
fhZNCvsZNA(0x0),
fhZPCvsZPA(0x0),
fhZDCCvsZDCCA(0x0),
fhZNCvsZPC(0x0),
fhZNAvsZPA(0x0),
fhZNvsZP(0x0),
fhZNvsVZERO(0x0),
fhZDCvsVZERO(0x0),
fhZDCvsTracklets(0x0),
fhZDCvsNclu1(0x0),
fhDebunch(0x0),
fhZNCcentroid(0x0),
fhZNAcentroid(0x0),
fhAsymm(0x0),
fhZNAvsAsymm(0x0),
fhZNCvsAsymm(0x0),
fhZNCvscentrality(0x0),
fhZNAvscentrality(0x0),
fhZPCvscentrality(0x0),
fhZPAvscentrality(0x0),
fhZNCpmcvscentr(0x0),
fhZNApmcvscentr(0x0),
fhZPCpmcvscentr(0x0),
fhZPApmcvscentr(0x0),
fhZNCpmcLR(0x0),
fhZNApmcLR(0x0),
fhZPCpmcLR(0x0),
fhZPApmcLR(0x0),
fCRCnRun(0),
fZDCGainAlpha(0.395),
fDataSet("2010"),
fStack(0x0),
fCutTPC(kFALSE),
fCenDis(0x0),
fMultSelection(0x0),
fCentrality(0x0),
fTowerEqList(NULL),
fCachedRunNum(0)
{
 for(int i=0; i<5; i++){
  fhZNCPM[i] = 0x0;
  fhZNAPM[i] = 0x0;
  fhZPCPM[i] = 0x0;
  fhZPAPM[i] = 0x0;
  fhZNCPMlg[i] = 0x0;
  fhZNAPMlg[i] = 0x0;
  fhZPCPMlg[i] = 0x0;
  fhZPAPMlg[i] = 0x0;
 }
 for(int i=0; i<2; i++) fhZEM[i] = 0x0;
 for(int i=0; i<6; i++){
  fhTDCraw[i] = 0x0;
  fhTDC[i] = 0x0;
 }
 for(int i=0; i<4; i++){
  fhZNCPMQiPMC[i] = 0x0;
  fhZNAPMQiPMC[i] = 0x0;
  fhZPCPMQiPMC[i] = 0x0;
  fhZPAPMQiPMC[i] = 0x0;
  fhPMCvsPMQ[i] = 0x0;
 }
 for(Int_t r=0; r<fCRCMaxnRun; r++) {
  fRunList[r] = 0;
 }
  for(Int_t k=0; k<8; k++) {
    fTowerGainEq[k] =  NULL;
  }
 this->InitializeRunArrays();
  fMyTRandom3 = new TRandom3(1);
  gRandom->SetSeed(fMyTRandom3->Integer(65539));
  for(Int_t c=0; c<10; c++) {
    fPtSpecGen[c] = NULL;
    fPtSpecFB32[c] = NULL;
    fPtSpecFB96[c] = NULL;
    fPtSpecFB128[c] = NULL;
    fPtSpecFB768[c] = NULL;
  }
}

//________________________________________________________________________
AliAnalysisTaskCRCZDC::AliAnalysisTaskCRCZDC(const char *name, TString RPtype, Bool_t on, TString DataSet, UInt_t iseed, Bool_t bCandidates):
AliAnalysisTaskSE(name),
fAnalysisType("AUTOMATIC"),
fRPType(RPtype),
fCFManager1(NULL),
fCFManager2(NULL),
fCutsEvent(NULL),
fCutsRP(NULL),
fCutsPOI(NULL),
fCutContainer(new TList()),
fQAList(NULL),
fMinMult(0),
fMaxMult(10000000),
fMinA(-1.0),
fMaxA(-0.01),
fMinB(0.01),
fMaxB(1.0),
fQAon(on),
fLoadCandidates(bCandidates),
fNbinsMult(10000),
fNbinsPt(100),
fNbinsPhi(100),
fNbinsEta(200),
fNbinsQ(500),
fNbinsMass(1),
fMultMin(0.),
fMultMax(10000.),
fPtMin(0.),
fPtMax(10.),
fPhiMin(0.),
fPhiMax(TMath::TwoPi()),
fEtaMin(-5.),
fEtaMax(5.),
fQMin(0.),
fQMax(3.),
fMassMin(-1.),
fMassMax(0.),
fHistWeightvsPhiMin(0.),
fHistWeightvsPhiMax(3.),
fExcludedEtaMin(0.),
fExcludedEtaMax(0.),
fExcludedPhiMin(0.),
fExcludedPhiMax(0.),
fAfterburnerOn(kFALSE),
fNonFlowNumberOfTrackClones(0),
fV1(0.),
fV2(0.),
fV3(0.),
fV4(0.),
fV5(0.),
fDifferentialV2(0),
fFlowEvent(NULL),
fShuffleTracks(kFALSE),
fMyTRandom3(NULL),
fAnalysisInput(kAOD),
fIsMCInput(kFALSE),
fUseMCCen(kTRUE),
fRejectPileUp(kTRUE),
fCentrLowLim(0.),
fCentrUpLim(100.),
fCentrEstimator("V0M"),
fOutput(0x0),
fhZNCvsZNA(0x0),
fhZPCvsZPA(0x0),
fhZDCCvsZDCCA(0x0),
fhZNCvsZPC(0x0),
fhZNAvsZPA(0x0),
fhZNvsZP(0x0),
fhZNvsVZERO(0x0),
fhZDCvsVZERO(0x0),
fhZDCvsTracklets(0x0),
fhZDCvsNclu1(0x0),
fhDebunch(0x0),
fhZNCcentroid(0x0),
fhZNAcentroid(0x0),
fhAsymm(0x0),
fhZNAvsAsymm(0x0),
fhZNCvsAsymm(0x0),
fhZNCvscentrality(0x0),
fhZNAvscentrality(0x0),
fhZPCvscentrality(0x0),
fhZPAvscentrality(0x0),
fhZNCpmcvscentr(0x0),
fhZNApmcvscentr(0x0),
fhZPCpmcvscentr(0x0),
fhZPApmcvscentr(0x0),
fhZNCpmcLR(0x0),
fhZNApmcLR(0x0),
fhZPCpmcLR(0x0),
fhZPApmcLR(0x0),
fDataSet(DataSet),
fCRCnRun(0),
fZDCGainAlpha(0.395),
fGenHeader(NULL),
fPythiaGenHeader(NULL),
fHijingGenHeader(NULL),
fFlowTrack(NULL),
fAnalysisUtil(NULL),
fStack(0x0),
fCutTPC(kFALSE),
fCenDis(0x0),
fMultSelection(0x0),
fCentrality(0x0),
fTowerEqList(NULL),
fCachedRunNum(0)
{
 
 for(int i=0; i<5; i++){
  fhZNCPM[i] = 0x0;
  fhZNAPM[i] = 0x0;
  fhZPCPM[i] = 0x0;
  fhZPAPM[i] = 0x0;
  fhZNCPMlg[i] = 0x0;
  fhZNAPMlg[i] = 0x0;
  fhZPCPMlg[i] = 0x0;
  fhZPAPMlg[i] = 0x0;
 }
 for(int i=0; i<2; i++) fhZEM[i] = 0x0;
 for(int i=0; i<6; i++){
  fhTDCraw[i] = 0x0;
  fhTDC[i] = 0x0;
 }
 for(int i=0; i<4; i++){
  fhZNCPMQiPMC[i] = 0x0;
  fhZNAPMQiPMC[i] = 0x0;
  fhZPCPMQiPMC[i] = 0x0;
  fhZPAPMQiPMC[i] = 0x0;
  fhPMCvsPMQ[i] = 0x0;
 }
 for(Int_t r=0; r<fCRCMaxnRun; r++) {
  fRunList[r] = 0;
 }
  for(Int_t k=0; k<8; k++) {
    fTowerGainEq[k] =  NULL;
  }
 this->InitializeRunArrays();
 fMyTRandom3 = new TRandom3(iseed);
 gRandom->SetSeed(fMyTRandom3->Integer(65539));
 
 DefineInput(0, TChain::Class());
 // Define output slots here
 // Define here the flow event output
 DefineOutput(1, AliFlowEventSimple::Class());
 DefineOutput(2, TList::Class());
 
  for(Int_t c=0; c<10; c++) {
    fPtSpecGen[c] = NULL;
    fPtSpecFB32[c] = NULL;
    fPtSpecFB96[c] = NULL;
    fPtSpecFB128[c] = NULL;
    fPtSpecFB768[c] = NULL;
  }
 
}

//________________________________________________________________________
AliAnalysisTaskCRCZDC::~AliAnalysisTaskCRCZDC()
{
 // Destructor
 if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
  delete fOutput; fOutput=0;
 }
 delete fMyTRandom3;
 delete fFlowEvent;
 delete fFlowTrack;
 delete fCutsEvent;
 if (fAnalysisUtil) delete fAnalysisUtil;
 if (fQAList) delete fQAList;
 if (fCutContainer) fCutContainer->Delete(); delete fCutContainer;
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::InitializeRunArrays()
{
 for(Int_t r=0;r<fCRCMaxnRun;r++) {
  fCRCQVecListRun[r] = NULL;
  for(Int_t i=0;i<fCRCnTow;i++) {
   fhnTowerGain[r][i] = NULL;
  }
 }
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::UserCreateOutputObjects()
{
 // Create the output containers
 
 if (!(fAnalysisType == "AOD" || fAnalysisType == "MCkine" || fAnalysisType == "MCAOD" || fAnalysisType == "AUTOMATIC" || fAnalysisType == "MCESD"))
 {
  AliError("WRONG ANALYSIS TYPE! only MCESD, MCkine, MCAOD, AOD and AUTOMATIC are allowed.");
  exit(1);
 }
 
 //set the common constants
 AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
 cc->SetNbinsMult(fNbinsMult);
 cc->SetNbinsPt(fNbinsPt);
 cc->SetNbinsPhi(fNbinsPhi);
 cc->SetNbinsEta(fNbinsEta);
 cc->SetNbinsQ(fNbinsQ);
 cc->SetNbinsMass(fNbinsMass);
 cc->SetMultMin(fMultMin);
 cc->SetMultMax(fMultMax);
 cc->SetPtMin(fPtMin);
 cc->SetPtMax(fPtMax);
 cc->SetPhiMin(fPhiMin);
 cc->SetPhiMax(fPhiMax);
 cc->SetEtaMin(fEtaMin);
 cc->SetEtaMax(fEtaMax);
 cc->SetQMin(fQMin);
 cc->SetQMax(fQMax);
 cc->SetMassMin(fMassMin);
 cc->SetMassMax(fMassMax);
 cc->SetHistWeightvsPhiMax(fHistWeightvsPhiMax);
 cc->SetHistWeightvsPhiMin(fHistWeightvsPhiMin);
 
 fFlowEvent = new AliFlowEvent(20000);
 fFlowTrack = new AliFlowTrack();
 
 //printf("  AliAnalysisTaskCRCZDC::UserCreateOutputObjects()\n\n");
 fOutput = new TList();
 fOutput->SetOwner(kTRUE);
 //fOutput->SetName("output");
 
 if (fQAon) {
  fQAList = new TList();
  fQAList->SetOwner(kTRUE);
  fQAList->SetName("AliFlowEventCuts QA");
  if (fCutsEvent->GetQA()) fQAList->Add(fCutsEvent->GetQA()); //0
  if (fCutsRP->GetQA()) fQAList->Add(fCutsRP->GetQA());       //1
  if (fCutsPOI->GetQA())fQAList->Add(fCutsPOI->GetQA());      //2
  fOutput->Add(fQAList);
 }
  
  fCenDis = new TH1F("fCenDis", "fCenDis", 100, 0., 100.);
  fOutput->Add(fCenDis);
  
  for(Int_t k=0; k<8; k++) {
    fTowerGainEq[k] =  new TH1D(Form("fTowerGainEq[%d]",k),Form("fTowerGainEq[%d]",k),20,0.,100.);
    fOutput->Add(fTowerGainEq[k]);
  }
 
 Float_t xmin[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.33,2.66,3.,3.5,4.,4.5,5.,6.,7.,8.};
 for(Int_t c=0; c<10; c++) {
   fPtSpecGen[c] = new TH1F(Form("fPtSpecGen[%d]",c), Form("fPtSpecGen[%d]",c), 24, xmin);
   fOutput->Add(fPtSpecGen[c]);
   fPtSpecFB32[c] = new TH1F(Form("fPtSpecFB32[%d]",c), Form("fPtSpecFB32[%d]",c), 24, xmin);
   fOutput->Add(fPtSpecFB32[c]);
   fPtSpecFB96[c] = new TH1F(Form("fPtSpecFB96[%d]",c), Form("fPtSpecFB96[%d]",c), 24, xmin);
   fOutput->Add(fPtSpecFB96[c]);
   fPtSpecFB128[c] = new TH1F(Form("fPtSpecFB128[%d]",c), Form("fPtSpecFB128[%d]",c), 24, xmin);
   fOutput->Add(fPtSpecFB128[c]);
   fPtSpecFB768[c] = new TH1F(Form("fPtSpecFB768[%d]",c), Form("fPtSpecFB768[%d]",c), 24, xmin);
   fOutput->Add(fPtSpecFB768[c]);
 }
 
 fAnalysisUtil = new AliAnalysisUtils();
 
 for(int i=0; i<5; i++){
  char hname[20];
  sprintf(hname,"hZNCPM%d",i);
  fhZNCPM[i] = new TH1F(hname, hname, 200, -50., 140000);
  fOutput->Add(fhZNCPM[i]);
  //
  sprintf(hname,"hZNAPM%d",i);
  fhZNAPM[i] = new TH1F(hname, hname, 200, -50., 140000);
  fOutput->Add(fhZNAPM[i]);
  //
  sprintf(hname,"hZPCPM%d",i);
  fhZPCPM[i] = new TH1F(hname, hname, 200, -50., 50000);
  fOutput->Add(fhZPCPM[i]);
  //
  sprintf(hname,"hZPAPM%d",i);
  fhZPAPM[i] = new TH1F(hname, hname, 200, -50., 50000);
  fOutput->Add(fhZPAPM[i]);
  //
  sprintf(hname,"hZNCPMlg%d",i);
  fhZNCPMlg[i] = new TH1F(hname, hname, 200, -50., 140000);
  fOutput->Add(fhZNCPMlg[i]);
  //
  sprintf(hname,"hZNAPMlg%d",i);
  fhZNAPMlg[i] = new TH1F(hname, hname, 200, -50., 140000);
  fOutput->Add(fhZNAPMlg[i]);
  //
  sprintf(hname,"hZPCPMlg%d",i);
  fhZPCPMlg[i] = new TH1F(hname, hname, 200, -50., 50000);
  fOutput->Add(fhZPCPMlg[i]);
  //
  sprintf(hname,"hZPAPMlg%d",i);
  fhZPAPMlg[i] = new TH1F(hname, hname, 200, -50., 50000);
  fOutput->Add(fhZPAPMlg[i]);
  //
  if(i<4){
   //
   char hnamenc[20];
   sprintf(hnamenc, "hZNCPMQ%dPMC",i+1);
   fhZNCPMQiPMC[i] = new TH1F(hnamenc, hnamenc, 100, 0., 1.);
   fOutput->Add(fhZNCPMQiPMC[i]);
   //
   char hnamena[20];
   sprintf(hnamena, "hZNAPMQ%dPMC",i+1);
   fhZNAPMQiPMC[i] = new TH1F(hnamena, hnamena, 100, 0., 1.);
   fOutput->Add(fhZNAPMQiPMC[i]);
   //
   char hnamepc[20];
   sprintf(hnamepc, "hZPCPMQ%dPMC",i+1);
   fhZPCPMQiPMC[i] = new TH1F(hnamepc, hnamepc, 100, 0., 1.);
   fOutput->Add(fhZPCPMQiPMC[i]);
   //
   char hnamepa[20];
   sprintf(hnamepa, "hZPAPMQ%dPMC",i+1);
   fhZPAPMQiPMC[i] = new TH1F(hnamepa, hnamepa, 100, 0., 1.);
   fOutput->Add(fhZPAPMQiPMC[i]);
  }
 }
 for(int i=0; i<6; i++){
  if(i==0){
   fhZEM[i] = new TH1F("hZEM1","hZEM1",200,-10.,1190.);
   fhTDCraw[i] = new TH1F("hTDCZEM1raw", "hTDCZEM1raw", 200, -200., 0.);
   fhTDC[i] = new TH1F("hTDCZEM1", "hTDCZEM1", 200, -150., 50.);
   fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZNC","hPMCvsPMQZNC",200,-10.,140000,200,-10.,140000);
   //
   fOutput->Add(fhZEM[i]);
   fOutput->Add(fhPMCvsPMQ[i]);
  }
  else if(i==1){
   fhZEM[i] = new TH1F("hZEM2","hZEM2",200,-10.,1190.);
   fhTDCraw[i] = new TH1F("hTDCZEM2raw", "hTDCZEM2raw", 200, -200., 0.);
   fhTDC[i] = new TH1F("hTDCZEM2", "hTDCZEM2", 200, -150., 50.);
   fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZPC","hPMCvsPMQZPC",200,-10.,50000,200,-10.,50000);
   //
   fOutput->Add(fhZEM[i]);
   fOutput->Add(fhPMCvsPMQ[i]);
  }
  else if(i==2){
   fhTDCraw[i] = new TH1F("hTDCZNCraw", "hTDCZNCraw", 200, -200., 0.);
   fhTDC[i] = new TH1F("hTDCZNC", "hTDCZNC", 200, -150., 50.);
   fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZNA","hPMCvsPMQZNA",200,-10.,140000,200,-10.,140000);
   //
   fOutput->Add(fhPMCvsPMQ[i]);
  }
  else if(i==3){
   fhTDCraw[i] = new TH1F("hTDCZPCraw", "hTDCZPCraw", 200, -200., 0.);
   fhTDC[i] = new TH1F("hTDCZPC", "hTDCZPC", 200, -150., 50.);
   fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZPA","hPMCvsPMQZPA",200,-10.,50000,200,-10.,50000);
   //
   fOutput->Add(fhPMCvsPMQ[i]);
  }
  else if(i==4){
   fhTDCraw[i] = new TH1F("hTDCZNAraw", "hTDCZNAraw", 200, -200., 0.);
   fhTDC[i] = new TH1F("hTDCZNA", "hTDCZNA", 200, -150., 50.);
  }
  else if(i==5){
   fhTDCraw[i] = new TH1F("hTDCZPAraw", "hTDCZPAraw", 200, -200., 0.);
   fhTDC[i] = new TH1F("hTDCZPA", "hTDCZPA", 200, -150., 50.);
  }
  //
  fOutput->Add(fhTDC[i]);
 }

 fhZNCvsZNA = new TH2F("hZNCvsZNA","hZNCvsZNA",200,-50.,140000,200,-50.,140000);
 fOutput->Add(fhZNCvsZNA);
 fhZPCvsZPA = new TH2F("hZPCvsZPA","hZPCvsZPA",200,-50.,50000,200,-50.,50000);
 fOutput->Add(fhZPCvsZPA);
 fhZDCCvsZDCCA = new TH2F("hZDCCvsZDCCA","hZDCCvsZDCCA",200,0.,180000.,200,0.,200000.);
 fOutput->Add(fhZDCCvsZDCCA);
 fhZNCvsZPC = new TH2F("hZNCvsZPC","hZNCvsZPC",200,-50.,50000,200,-50.,140000);
 fOutput->Add(fhZNCvsZPC);
 fhZNAvsZPA = new TH2F("hZNAvsZPA","hZNAvsZPA",200,-50.,50000,200,-50.,140000);
 fOutput->Add(fhZNAvsZPA);
 fhZNvsZP = new TH2F("hZNvsZP","hZNvsZP",200,-50.,80000,200,-50.,200000);
 fOutput->Add(fhZNvsZP);
 fhZNvsVZERO = new TH2F("hZNvsVZERO","hZNvsVZERO",250,0.,25000.,200,0.,200000.);
 fOutput->Add(fhZNvsVZERO);
 fhZDCvsVZERO = new TH2F("hZDCvsVZERO","hZDCvsVZERO",250,0.,25000.,250,0.,250000.);
 fOutput->Add(fhZDCvsVZERO);
 fhZDCvsTracklets = new TH2F("hZDCvsTracklets","hZDCvsTracklets",200,0.,4000.,250,0.,250000.);
 fOutput->Add(fhZDCvsTracklets);
 fhZDCvsNclu1 = new TH2F("hZDCvsNclu1", "hZDCvsNclu1", 200, 0.,8000.,200,0.,250000.);
 fOutput->Add(fhZDCvsNclu1);
 fhDebunch = new TH2F("hDebunch","hDebunch",240,-100.,-40.,240,-30.,30.);
 fOutput->Add(fhDebunch);
 fhZNCcentroid = new TH2F("hZNCcentroid","hZNCcentroid",100,-3.5,3.5,100,-3.5,3.5);
 fOutput->Add(fhZNCcentroid);
 fhZNAcentroid = new TH2F("hZNAcentroid","hZNAcentroid",100,-3.5,3.5,100,-3.5,3.5);
 fOutput->Add(fhZNAcentroid);
 
 fhAsymm = new TH1F("hAsymm" , "Asimmetry ",200,-1.,1.);
 fOutput->Add(fhAsymm);
 fhZNAvsAsymm = new TH2F("hZNAvsAsymm","ZNA vs. asymm.",200,-1.,1.,200,0.,80.);
 fOutput->Add(fhZNAvsAsymm);
 fhZNCvsAsymm = new TH2F("hZNCvsAsymm","ZNC vs. asymm.",200,-1.,1.,200,0.,80.);
 fOutput->Add(fhZNCvsAsymm);
 
 fhZNCvscentrality = new TH2F("hZNCvscentrality","ZNC vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZNCvscentrality);
 fhZNAvscentrality = new TH2F("hZNAvscentrality","ZNA vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZNAvscentrality);
 fhZPCvscentrality = new TH2F("hZPCvscentrality","ZPC vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZPCvscentrality);
 fhZPAvscentrality = new TH2F("hZPAvscentrality","ZPA vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZPAvscentrality);
 
 fhZNCpmcvscentr = new TH2F("hZNCpmcvscentr","ZNC PMC vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZNCpmcvscentr);
 fhZNApmcvscentr = new TH2F("hZNApmcvscentr","ZNA PMC vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZNApmcvscentr);
 fhZPCpmcvscentr = new TH2F("hZPCpmcvscentr","ZPC PMC vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZPCpmcvscentr);
 fhZPApmcvscentr = new TH2F("hZPApmcvscentr","ZPA PMC vs. centrality",100,0.,100.,100,0.,100.);
 fOutput->Add(fhZPApmcvscentr);
 
 fhZNCpmcLR = new TH1F("hZNCpmcLR","ZNC PMC lr", 100, 0., 10.);
 fOutput->Add(fhZNCpmcLR);
 fhZNApmcLR = new TH1F("hZNApmcLR","ZNA PMC lr", 100, 0., 10.);
 fOutput->Add(fhZNApmcLR);
 fhZPCpmcLR = new TH1F("hZPCpmcLR","ZPC PMC lr", 100, 0., 10.);
 fOutput->Add(fhZPCpmcLR);
 fhZPApmcLR = new TH1F("hZPApmcLR","ZPA PMC lr", 100, 0., 10.);
 fOutput->Add(fhZPApmcLR);
 
 //********************************************************************
 
 Int_t dRun10h[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};
 
 Int_t dRun11h[] = {167902, 167903, 167915, 167920, 167985, 167987, 167988, 168066, 168068, 168069, 168076, 168104, 168105, 168107, 168108, 168115, 168212, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168461, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168984, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169143, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 169965, 170027, 170036,170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593};
  
 Int_t dRun15h[] = {244917, 244918, 244975, 244980, 244982, 244983, 245061, 245064, 245066, 245068};
 
 if(fDataSet.EqualTo("2010")) {fCRCnRun=92;}
 if(fDataSet.EqualTo("2011")) {fCRCnRun=119;}
 if(fDataSet.EqualTo("2015")) {fCRCnRun=10;}
 if(fDataSet.EqualTo("MCkine")) {fCRCnRun=1;}
 
 Int_t d=0;
 for(Int_t r=0; r<fCRCnRun; r++) {
  if(fDataSet.EqualTo("2010"))   fRunList[d] = dRun10h[r];
  if(fDataSet.EqualTo("2011"))   fRunList[d] = dRun11h[r];
  if(fDataSet.EqualTo("2015"))   fRunList[d] = dRun15h[r];
  if(fDataSet.EqualTo("MCkine")) fRunList[d] = 1;
  d++;
 }

 for(Int_t r=0;r<fCRCnRun;r++) {
  fCRCQVecListRun[r] = new TList();
  fCRCQVecListRun[r]->SetName(Form("Run %d",fRunList[r]));
  fCRCQVecListRun[r]->SetOwner(kTRUE);
  fOutput->Add(fCRCQVecListRun[r]);
  for(Int_t i=0;i<fCRCnTow;i++) {
   fhnTowerGain[r][i] = new TH2D(Form("fhnTowerGain[%d][%d]",fRunList[r],i),
                                 Form("fhnTowerGain[%d][%d]",fRunList[r],i),20,0.,100.,100,0.,100);
   fhnTowerGain[r][i]->Sumw2();
   fCRCQVecListRun[r]->Add(fhnTowerGain[r][i]);
  }
 }
 
 PostData(2, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::UserExec(Option_t */*option*/)
{
 // Execute analysis for current event:
 AliMCEvent* McEvent = MCEvent();
 AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
 AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
// AliMultiplicity* myTracklets = NULL;
// AliESDPmdTrack* pmdtracks = NULL;
// int availableINslot=1;
 
 if (!(fCutsRP&&fCutsPOI&&fCutsEvent)) {
  AliError("cuts not set");
  return;
 }
 
 //DEFAULT - automatically takes care of everything
 if (fAnalysisType == "AUTOMATIC") {
  
  //check event cuts
  if (InputEvent()) {
   if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
   if(fRejectPileUp && fAnalysisUtil->IsPileUpEvent(InputEvent())) return;
  }
  
  //first attach all possible information to the cuts
  fCutsRP->SetEvent( InputEvent(), MCEvent() );  //attach event
  fCutsPOI->SetEvent( InputEvent(), MCEvent() );
  
  //then make the event
  fFlowEvent->Fill( fCutsRP, fCutsPOI );
  
  fFlowEvent->SetReferenceMultiplicity(fCutsEvent->GetReferenceMultiplicity(InputEvent(),McEvent));
  fFlowEvent->SetCentrality(fCutsEvent->GetCentrality(InputEvent(),McEvent));
   
  fFlowEvent->SetCentralityCL1(((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1"));
  fFlowEvent->SetCentralityTRK(((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK"));
  fFlowEvent->SetNITSCL1(((AliVAODHeader*)aod->GetHeader())->GetNumberOfITSClusters(1));
   
   Double_t vtxpos[3]={0.,0.,0.};
   vtxpos[0] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetX();
   vtxpos[1] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetY();
   vtxpos[2] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetZ();
   fFlowEvent->SetVertexPosition(vtxpos);
   
  if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
  
 }
 
  if (fAnalysisType == "MCAOD") {
    
    //check event cuts
    if (InputEvent()) {
      if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
      if(fRejectPileUp && fAnalysisUtil->IsPileUpEvent(InputEvent())) return;
    }
    
    fFlowEvent->ClearFast();
    
    if(!McEvent) {
      AliError("ERROR: Could not retrieve MCEvent");
      return;
    }
    fStack = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
    if(!fStack){
      AliError("ERROR: Could not retrieve MCStack");
      return;
    }
    
    // get centrality (from AliMultSelection or AliCentrality)
    Float_t centr = 300;
    if(fDataSet == "2015") {
      fMultSelection = (AliMultSelection*)aod->FindListObject("MultSelection");
      if(!fMultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
      } else {
        centr = fMultSelection->GetMultiplicityPercentile("V0M");
      }
    } else {
      centr = (((AliVAODHeader*)aod->GetHeader())->GetCentralityP())->GetCentralityPercentile("V0M");
    }
    fCenDis->Fill(centr);
    // centrality bin
    if (centr<fCentrLowLim || centr>=fCentrUpLim ) return;
    Int_t CenBin = -1;
    if (centr>0. && centr<5.) CenBin=0;
    if (centr>5. && centr<10.) CenBin=1;
    if (centr>10. && centr<20.) CenBin=2;
    if (centr>20. && centr<30.) CenBin=3;
    if (centr>30. && centr<40.) CenBin=4;
    if (centr>40. && centr<50.) CenBin=5;
    if (centr>50. && centr<60.) CenBin=6;
    if (centr>60. && centr<70.) CenBin=7;
    if (centr>70. && centr<80.) CenBin=8;
    if (centr>80. && centr<90.) CenBin=9;
    if(CenBin==-1) return;
    
    // reconstructed
    for(Int_t jTracks = 0; jTracks<aod->GetNumberOfTracks(); jTracks++){
      
      AliAODTrack* track = (AliAODTrack*)aod->GetTrack(jTracks);
      if(!track) continue;
      
      // select primaries
      Int_t lp = TMath::Abs(track->GetLabel());
      if (!((AliAODMCParticle*)fStack->At(lp))->IsPhysicalPrimary()) continue;
      
      // general kinematic & quality cuts
      if (track->Pt() < .2 || track->Pt() > 8. || TMath::Abs(track->Eta()) > .8 || track->GetTPCNcls() < 70)  continue;
      
      // test filter bits
      if (track->TestFilterBit(32)) fPtSpecFB32[CenBin]->Fill(track->Pt());
      if (track->TestFilterBit(96)) fPtSpecFB96[CenBin]->Fill(track->Pt());
      if (track->TestFilterBit(128)) fPtSpecFB128[CenBin]->Fill(track->Pt());
      if (track->TestFilterBit(768)) fPtSpecFB768[CenBin]->Fill(track->Pt());
    }
    
    // generated (physical primaries)xw
    
    for(Int_t jTracks = 0; jTracks<fStack->GetEntriesFast(); jTracks++) {
      AliAODMCParticle *MCpart = (AliAODMCParticle*)fStack->At(jTracks);
      if (!MCpart) {
        printf("ERROR: Could not receive MC track %d\n", jTracks);
        continue;
      }
      // kinematic cuts
      if ( MCpart->Pt() < 0.2 || MCpart->Pt() > 8. || TMath::Abs(MCpart->Eta()) > .8 ) continue;
      // select charged primaries
      if ( MCpart->Charge() == 0. || !MCpart->IsPhysicalPrimary()) continue;
      
      fPtSpecGen[CenBin]->Fill(MCpart->Pt());
    }
    
//    fGenHeader = McEvent->GenEventHeader();
//    if(fGenHeader) fPythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(fGenHeader);
    //  printf("#reconstructed : %d (rejected from cuts %d), #MC primaries : %d (rejected from cuts %d) \n",AODPOIs,AODbads,MCPrims,MCSecos);
    fFlowEvent->SetReferenceMultiplicity(aod->GetNumberOfTracks());
    fFlowEvent->SetCentrality(centr);
//    if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
    fFlowEvent->SetRun(aod->GetRunNumber());
    //  printf("Run : %d, RefMult : %d, Cent : %f \n",fFlowEvent->GetRun(),fFlowEvent->GetReferenceMultiplicity(),fFlowEvent->GetCentrality());
  }
  
 if(fAnalysisType ==  "MCESD") {
  
  fFlowEvent->ClearFast();
  
  if(!esd) {
   AliError("ERROR: Could not retrieve ESDEvent");
   return;
  }
  if(!McEvent) {
   AliError("ERROR: Could not retrieve MCEvent");
   return;
  }
  AliStack* fStack = fMCEvent->Stack();
  if(!fStack) {
   AliError("ERROR: Could not retrieve MCStack");
   return;
  }
  
  AliESDVertex *vertex = (AliESDVertex*) esd->GetPrimaryVertex();
  if (!vertex) return;
  if (TMath::Abs(vertex->GetZ()) > 10. ) return;
  if (vertex->GetNContributors() < 1 ) return;
  AliCentrality *centrality = esd->GetCentrality();
  if (!centrality) return;
  Double_t centr = centrality->GetCentralityPercentile("V0M");
  if (centr<fCentrLowLim || centr>=fCentrUpLim ) return;
  Int_t CenBin = -1;
  if (centr>0. && centr<5.) CenBin=0;
  if (centr>5. && centr<10.) CenBin=1;
  if (centr>10. && centr<20.) CenBin=2;
  if (centr>20. && centr<30.) CenBin=3;
  if (centr>30. && centr<40.) CenBin=4;
  if (centr>40. && centr<50.) CenBin=5;
  if (centr>50. && centr<60.) CenBin=6;
  if (centr>60. && centr<70.) CenBin=7;
  if (centr>70. && centr<80.) CenBin=8;
  if (centr>80. && centr<90.) CenBin=9;
  if(CenBin==-1) return;
  
  //Generated
  Int_t MCPrims = 0;
  for ( Int_t i=0 ; i<fStack->GetNtrack() ; i++ )  {
   
   //Primaries Selection
   TParticle *particle = (TParticle*)fStack->Particle(i);
   if (!particle) continue;
   if (!fStack->IsPhysicalPrimary(i)) continue;
   if ( particle->GetPDG()->Charge() == 0.) continue;
   
   //Kinematic Cuts
   if ( particle->Pt()<0.2 || particle->Pt()>10. ) continue;
   if ( TMath::Abs(particle->Eta())>0.8 ) continue;
   
   fFlowTrack->SetPhi(particle->Phi());
   fFlowTrack->SetEta(particle->Eta());
   fFlowTrack->SetPt(particle->Pt());
   fFlowTrack->SetSource(AliFlowTrack::kFromMC);
   fFlowTrack->SetForRPSelection(kTRUE);
   fFlowEvent->IncrementNumberOfPOIs(0);
   fFlowTrack->SetForPOISelection(kFALSE);
   fFlowEvent->InsertTrack(fFlowTrack);
   MCPrims++;
   
   fPtSpecGen[CenBin]->Fill(particle->Pt());
   
  }
  
  //Reconstructed
  Int_t ESDPrims = 0;
  for (Int_t i=0 ; i<esd->GetNumberOfTracks() ; i++)  {
   
   //Get reconstructed track
   AliVTrack *vtrack = static_cast<AliVTrack*>(esd->GetTrack(i));
   AliESDtrack *track = dynamic_cast<AliESDtrack*>(vtrack);
   if (!track) continue;
   
   //Primaries selection
   Int_t lp = TMath::Abs(track->GetLabel());
   if (!fStack->IsPhysicalPrimary(lp)) continue;
   TParticle *particle = (TParticle*)fStack->Particle(lp);
   if (!particle) continue;
   if (particle->GetPDG()->Charge() == 0.) continue;
   
//   if(!fCutsPOI->PassesESDcuts(track)) continue;
   
   Bool_t pass = kTRUE;
   
   if(fCutTPC) {
//    printf("******* cutting TPC ******** \n");
    UShort_t ntpccls = track->GetTPCNcls();
    Double_t tpcchi2 = track->GetTPCchi2();
    if (tpcchi2<0.2 || tpcchi2 >=4.) {
//     printf("TPCchi2 : %e %e ",tpcchi2,track->GetTPCchi2Iter1());
     pass=kFALSE;
    }
    if (ntpccls < 70) {
//     printf("#TPCcluster : %u %u %u %u ",ntpccls,track->GetTPCNclsF(),track->GetTPCNclsFIter1(),track->GetTPCNclsIter1());
     pass=kFALSE;
    }
   }
   
   Float_t dcaxy=0.0;
   Float_t dcaz=0.0;
   track->GetImpactParameters(dcaxy,dcaz);
   if (dcaxy > 0.3 || dcaz > 0.3) {
//    printf("DCA : %e %e ",dcaxy,dcaz);
    pass=kFALSE;
   }
   if(!pass) continue;
   
   //Kinematic Cuts
   if ( track->Pt()<0.2 || track->Pt()>10. ) continue;
   if ( TMath::Abs(track->Eta())>0.8 ) continue;
   
   fFlowTrack->SetPhi(track->Phi());
   fFlowTrack->SetEta(track->Eta());
   fFlowTrack->SetPt(track->Pt());
   fFlowTrack->SetSource(AliFlowTrack::kFromESD);
   fFlowTrack->SetForRPSelection(kFALSE);
   fFlowTrack->SetForPOISelection(kTRUE);
   fFlowEvent->IncrementNumberOfPOIs(1);
   fFlowEvent->InsertTrack(fFlowTrack);
   ESDPrims++;
   
  }
  
//  printf("#reconstructed : %d , #MC primaries : %d \n",ESDPrims,MCPrims);
  fFlowEvent->SetReferenceMultiplicity(esd->GetNumberOfTracks());
  fFlowEvent->SetCentrality(centr);
  if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
  fFlowEvent->SetRun(esd->GetRunNumber());
//  printf("Run : %d, RefMult : %d, Cent : %f \n",fFlowEvent->GetRun(),fFlowEvent->GetReferenceMultiplicity(),fFlowEvent->GetCentrality());
  
 } // end of if(fAnalysisType ==  "MCESD")
 
 if(fAnalysisType ==  "MCkine") {
  
  fFlowEvent->ClearFast();
  
  AliInputEventHandler* McHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!McHandler) {
   AliError("ERROR: Could not retrieve MCtruthEventHandler");
   return;
  }
  McEvent = McHandler->MCEvent();
  if(!McEvent) {
   AliError("ERROR: Could not retrieve MC event");
   return;
  }
  
  Int_t nTracks = McEvent->GetNumberOfTracks();
//  Int_t nPrimTr = McEvent->GetNumberOfPrimaries();
  
  //loop over tracks
  for (Int_t itrkN=0; itrkN<nTracks; itrkN++) {
   //get input particle
   AliMCParticle* pParticle = dynamic_cast<AliMCParticle*>(McEvent->GetTrack(itrkN));
   if (!pParticle) continue;
   
   //check if track passes the cuts
   if (McEvent->IsPhysicalPrimary(itrkN) && pParticle->Charge()!=0) {
    fFlowTrack->Set(pParticle);
    fFlowTrack->SetSource(AliFlowTrack::kFromMC);
    fFlowTrack->SetForRPSelection(kTRUE);
    fFlowEvent->IncrementNumberOfPOIs(0);
    fFlowTrack->SetForPOISelection(kTRUE);
    fFlowEvent->IncrementNumberOfPOIs(1);
    fFlowEvent->InsertTrack(fFlowTrack);
   }
  }// for all tracks
  
  // if monte carlo event get reaction plane from monte carlo (depends on generator)
  if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
  // set reference multiplicity
  fFlowEvent->SetReferenceMultiplicity(McEvent->GetNumberOfTracks());
  // tag subevents
  fFlowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);
  // set centrality from impact parameter
  Double_t ImpPar=0., CenPer=0.;
  fGenHeader = McEvent->GenEventHeader();
  if(fGenHeader){
   fPythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(fGenHeader);
   if(fPythiaGenHeader) ImpPar = fPythiaGenHeader->GetImpactParameter();
   fHijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(fGenHeader);
   if(fHijingGenHeader) ImpPar = fHijingGenHeader->ImpactParameter();
   if(ImpPar) CenPer = 0.3859796743103508*pow(ImpPar,2.);
   if(CenPer>0. && CenPer<100.) fFlowEvent->SetCentrality(CenPer);
   else return;
   fFlowEvent->SetRun(1);
  }
  
 } // end of if(fAnalysisType ==  "MCkine")
 
 if (!fFlowEvent) return; //shuts up coverity
 
 //check final event cuts
 Int_t mult = fFlowEvent->NumberOfTracks();
 //  AliInfo(Form("FlowEvent has %i tracks",mult));
 if (mult<fMinMult || mult>fMaxMult) {
  AliWarning("FlowEvent cut on multiplicity"); return;
 }
 
 //define dead zone
 fFlowEvent->DefineDeadZone(fExcludedEtaMin, fExcludedEtaMax, fExcludedPhiMin, fExcludedPhiMax );
 
 //////////////////////////////////////////////////////////////////////////////
 ///////////////////////////AFTERBURNER
 if (fAfterburnerOn)
 {
  //if reaction plane not set from elsewhere randomize it before adding flow
  if (!fFlowEvent->IsSetMCReactionPlaneAngle())
   fFlowEvent->SetMCReactionPlaneAngle(gRandom->Uniform(0.0,TMath::TwoPi()));
  
  if(fDifferentialV2)
   fFlowEvent->AddV2(fDifferentialV2);
  else
   fFlowEvent->AddFlow(fV1,fV2,fV3,fV4,fV5);     //add flow
  fFlowEvent->CloneTracks(fNonFlowNumberOfTrackClones); //add nonflow by cloning tracks
 }
 //////////////////////////////////////////////////////////////////////////////
 
 //tag subEvents
 fFlowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);
 
 //do we want to serve shullfed tracks to everybody?
 fFlowEvent->SetShuffleTracks(fShuffleTracks);
 
 // associate the mother particles to their daughters in the flow event (if any)
 fFlowEvent->FindDaughters();
 
 //fListHistos->Print();
 //fOutputFile->WriteObject(fFlowEvent,"myFlowEventSimple");
 
//  printf("event : ntr %d, cen %f **********************************************************************************************************\n"
 
 //********************************************************************************************************************************
 
  if(fAnalysisType == "AOD" || fAnalysisType == "AUTOMATIC") {
    
    // PHYSICS SELECTION
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *hdr = (AliInputEventHandler*)am->GetInputEventHandler();
    
    if(hdr->IsEventSelected() && AliVEvent::kAny) {
      
      AliCentrality* centrality = aod->GetCentrality();
      Float_t centrperc = centrality->GetCentralityPercentile(fCentrEstimator.Data());
      
      AliAODTracklets *trackl = aod->GetTracklets();
      Int_t nTracklets = trackl->GetNumberOfTracklets();
      
      AliAODVZERO *vzeroAOD = aod->GetVZEROData();
      Float_t multV0A = vzeroAOD->GetMTotV0A();
      Float_t multV0C = vzeroAOD->GetMTotV0C();
      
      AliAODZDC *aodZDC = aod->GetZDCData();
      
      Float_t energyZNC  = (Float_t) (aodZDC->GetZNCEnergy());
      Float_t energyZPC  = (Float_t) (aodZDC->GetZPCEnergy());
      Float_t energyZNA  = (Float_t) (aodZDC->GetZNAEnergy());
      Float_t energyZPA  = (Float_t) (aodZDC->GetZPAEnergy());
      Float_t energyZEM1 = (Float_t) (aodZDC->GetZEM1Energy());
      Float_t energyZEM2 = (Float_t) (aodZDC->GetZEM2Energy());
      
      const Double_t * towZNC = aodZDC->GetZNCTowerEnergy();
      const Double_t * towZPC = aodZDC->GetZPCTowerEnergy();
      const Double_t * towZNA = aodZDC->GetZNATowerEnergy();
      const Double_t * towZPA = aodZDC->GetZPATowerEnergy();
      //
      const Double_t * towZNClg = aodZDC->GetZNCTowerEnergyLR();
      const Double_t * towZNAlg = aodZDC->GetZNATowerEnergyLR();
      //
      Double_t towZPClg[5], towZPAlg[5]={0.};
      for(Int_t it=0; it<5; it++){
        towZPClg[it] = 8*towZPC[it];
        towZPAlg[it] = 8*towZNA[it];
      }
      
      // Get centroid from ZDCs *******************************************************
      
      Double_t xyZNC[2]={999.,999.}, xyZNA[2]={999.,999.};
      Float_t zncEnergy=0., znaEnergy=0.;
      for(Int_t i=0; i<5; i++){
        zncEnergy += towZNC[i];
        znaEnergy += towZNA[i];
      }
      fFlowEvent->SetZNCEnergy(towZNC[0]);
      fFlowEvent->SetZNAEnergy(towZNA[0]);
      
      Int_t RunNum=fFlowEvent->GetRun();
      if(fTowerEqList) {
        if(RunNum!=fCachedRunNum) {
          for(Int_t i=0; i<8; i++) {
            fTowerGainEq[i] = (TH1D*)(fTowerEqList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fhnTowerGainEqFactor[%d][%d]",RunNum,i)));
          }
        }
      }
      Double_t AvTowerGain[8] = {1., 1., 1., 1., 1., 1., 1., 1.};
      
      if (fUseMCCen) {
        if(aod->GetRunNumber() < 209122) aodZDC->GetZNCentroidInPbPb(1380., xyZNC, xyZNA);
        else                             aodZDC->GetZNCentroidInPbPb(2510., xyZNC, xyZNA);
      } else {
        // set tower gain equalization, if available
        if(fTowerEqList) {
          for(Int_t i=0; i<8; i++) {
            if(fTowerGainEq[i]) AvTowerGain[i] = fTowerGainEq[i]->GetBinContent(fTowerGainEq[i]->FindBin(centrperc));
          }
        }
        const Float_t x[4] = {-1.75, 1.75, -1.75, 1.75};
        const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};
        Float_t numXZNC=0., numYZNC=0., denZNC=0., wZNC;
        Float_t numXZNA=0., numYZNA=0., denZNA=0., wZNA;
        for(Int_t i=0; i<4; i++) {
          if(towZNC[i+1]>0.) {
            wZNC = TMath::Power(towZNC[i+1], fZDCGainAlpha)*AvTowerGain[i];
            numXZNC += x[i]*wZNC;
            numYZNC += y[i]*wZNC;
            denZNC += wZNC;
          }
          if(towZNA[i+1]>0.) {
            wZNA = TMath::Power(towZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
            numXZNA += x[i]*wZNA;
            numYZNA += y[i]*wZNA;
            denZNA += wZNA;
          }
        }
        if(denZNC!=0) {
          xyZNC[0] = numXZNC/denZNC;
          xyZNC[1] = numYZNC/denZNC;
        }
        else{
          xyZNC[0] = xyZNC[1] = 999.;
          zncEnergy = 0.;
        }
        if(denZNA!=0) {
          xyZNA[0] = numXZNA/denZNA;
          xyZNA[1] = numYZNA/denZNA;
        }
        else{
          xyZNA[0] = xyZNA[1] = 999.;
          znaEnergy = 0.;
        }
      }
      
      Float_t MulA=0., MulC=0.;
      for(Int_t i=0; i<4; i++) {
        if(towZNC[i+1]>0.) {
          MulC += TMath::Power(towZNC[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
        }
        if(towZNA[i+1]>0.) {
          MulA += TMath::Power(towZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
        }
      }
      
      fhZNCcentroid->Fill(xyZNC[0], xyZNC[1]);
      fhZNAcentroid->Fill(xyZNA[0], xyZNA[1]);
      fFlowEvent->SetZDC2Qsub(xyZNC,MulC,xyZNA,MulA);
    
      Int_t RunBin=-1, bin=0;
      for(Int_t c=0;c<fCRCnRun;c++) {
        if(fRunList[c]==RunNum) RunBin=bin;
        else bin++;
      }
      if(fDataSet.EqualTo("MCkine")) RunBin=0;
      
      if(RunBin!=-1) {
        for(Int_t i=0; i<4; i++) {
          if(towZNC[i+1]>0.) {
            fhnTowerGain[RunBin][i]->Fill(centrperc,TMath::Power(towZNC[i+1], fZDCGainAlpha)*AvTowerGain[i]);
          }
          if(towZNA[i+1]>0.) {
            fhnTowerGain[RunBin][i+4]->Fill(centrperc,TMath::Power(towZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4]);
          }
        }
      }
      
      // ******************************************************************************
      
      Float_t tdcSum = aodZDC->GetZDCTimeSum();
      Float_t tdcDiff = aodZDC->GetZDCTimeDiff();
      fhDebunch->Fill(tdcDiff, tdcSum);
      
      for(int i=0; i<5; i++){
        fhZNCPM[i]->Fill(towZNC[i]);
        fhZNCPMlg[i]->Fill(towZNClg[i]);
        if((i<4) && (towZNC[0]>0.)) fhZNCPMQiPMC[i]->Fill(towZNC[i+1]/towZNC[0]);
        fhZNCpmcLR->Fill(towZNClg[0]/1000.);
      }
      fhPMCvsPMQ[0]->Fill(towZNC[1]+towZNC[2]+towZNC[3]+towZNC[4], towZNC[0]);
      for(int i=0; i<5; i++){
        fhZPCPM[i]->Fill(towZPC[i]);
        fhZPCPMlg[i]->Fill(towZPClg[i]);
        if(((i<4) && towZPC[0]>0.)) fhZPCPMQiPMC[i]->Fill(towZPC[i+1]/towZPC[0]);
        fhZPCpmcLR->Fill(towZPClg[0]/1000.);
      }
      fhPMCvsPMQ[1]->Fill(towZPC[1]+towZPC[2]+towZPC[3]+towZPC[4], towZPC[0]);
      for(int i=0; i<5; i++){
        fhZNAPM[i]->Fill(towZNA[i]);
        fhZNAPMlg[i]->Fill(towZNAlg[i]);
        if(((i<4) && towZNA[0]>0.)) fhZNAPMQiPMC[i]->Fill(towZNA[i+1]/towZNA[0]);
        fhZNApmcLR->Fill(towZNAlg[0]/1000.);
      }
      fhPMCvsPMQ[2]->Fill(towZNA[1]+towZNA[2]+towZNA[3]+towZNA[4], towZNA[0]);
      for(int i=0; i<5; i++){
        fhZPAPM[i]->Fill(towZPA[i]);
        fhZPAPMlg[i]->Fill(towZPAlg[i]);
        if(((i<4) && towZPA[0]>0.)) fhZPAPMQiPMC[i]->Fill(towZPA[i+1]/towZPA[0]);
        fhZPApmcLR->Fill(towZPAlg[0]/1000.);
      }
      fhPMCvsPMQ[3]->Fill(towZPA[1]+towZPA[2]+towZPA[3]+towZPA[4], towZPA[0]);
      fhZEM[0]->Fill(energyZEM1);
      fhZEM[1]->Fill(energyZEM2);
      
      fhZNCvsZNA->Fill(energyZNA, energyZNC);
      fhZPCvsZPA->Fill(energyZPA, energyZPC);
      fhZDCCvsZDCCA->Fill(energyZNA+energyZPA, energyZNC+energyZPC);
      fhZNCvsZPC->Fill(energyZPC, energyZNC);
      fhZNAvsZPA->Fill(energyZPA, energyZNA);
      fhZNvsZP->Fill(energyZPA+energyZPC, energyZNA+energyZNC);
      fhZNvsVZERO->Fill(multV0A+multV0C, energyZNC+energyZNA);
      fhZDCvsVZERO->Fill(multV0A+multV0C, energyZNA+energyZPA+energyZNC+energyZPC);
      fhZDCvsTracklets->Fill((Float_t) (nTracklets), energyZNA+energyZPA+energyZNC+energyZPC);
      
      Double_t asymmetry = -999.;
      if((energyZNC+energyZNA)>0.) asymmetry = (energyZNC-energyZNA)/(energyZNC+energyZNA);
      fhAsymm->Fill(asymmetry);
      fhZNAvsAsymm->Fill(asymmetry, energyZNA/1000.);
      fhZNCvsAsymm->Fill(asymmetry, energyZNC/1000.);
      
      fhZNCvscentrality->Fill(centrperc, energyZNC/1000.);
      fhZNAvscentrality->Fill(centrperc, energyZNA/1000.);
      fhZPCvscentrality->Fill(centrperc, energyZPC/1000.);
      fhZPAvscentrality->Fill(centrperc, energyZPA/1000.);
      
      fhZNCpmcvscentr->Fill(centrperc, towZNC[0]/1380.);
      fhZNApmcvscentr->Fill(centrperc, towZNA[0]/1380.);
      fhZPCpmcvscentr->Fill(centrperc, towZPC[0]/1380.);
      fhZPApmcvscentr->Fill(centrperc, towZPA[0]/1380.);
      
    } // PHYSICS SELECTION
    
  }

  // p) cache run number
  fCachedRunNum = fFlowEvent->GetRun();
  
 PostData(1, fFlowEvent);
 
 PostData(2, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::SetCutsRP(AliFlowTrackCuts* cutsRP) {
 fCutContainer->AddAt(cutsRP,0); fCutsRP=cutsRP; cutsRP->SetPOItype(0);
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::SetCutsPOI(AliFlowTrackCuts* cutsPOI) {
 fCutContainer->AddAt(cutsPOI,1); fCutsPOI=cutsPOI; cutsPOI->SetPOItype(1);
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::Terminate(Option_t */*option*/)
{
 // Terminate analysis
 //
 /*  if(fDebug > 1) printf(" **** AliAnalysisTaskCRCZDC::Terminate() \n");
  
  //fOutput = dynamic_cast<TList*> (GetOutputData(1));
  //if(!fOutput) printf("ERROR: fOutput not available\n");
  */
}
