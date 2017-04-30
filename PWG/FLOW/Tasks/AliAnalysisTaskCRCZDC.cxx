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
#include <TH3D.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TParticle.h>
#include "TProfile2D.h"
#include "TProfile3D.h"

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
#include "AliLumiTools.h"

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
fRejectPileUpTight(kFALSE),
fResetNegativeZDC(kFALSE),
fCentrLowLim(0.),
fCentrUpLim(100.),
fCentrEstimator(kV0M),
fOutput(0x0),
fhZNCvsZNA(0x0),
fhZDCCvsZDCCA(0x0),
fhZNCvsZPC(0x0),
fhZNAvsZPA(0x0),
fhZNvsZP(0x0),
fhZNvsVZERO(0x0),
fhZDCvsVZERO(0x0),
fhZDCvsTracklets(0x0),
fhZDCvsNclu1(0x0),
fhDebunch(0x0),
fhAsymm(0x0),
fhZNAvsAsymm(0x0),
fhZNCvsAsymm(0x0),
fhZNCvscentrality(0x0),
fhZNAvscentrality(0x0),
fCRCnRun(0),
fZDCGainAlpha(0.395),
fDataSet(kAny),
fStack(0x0),
fCutTPC(kFALSE),
fCenDis(0x0),
fVZEROMult(0x0),
fMultSelection(0x0),
fPileUpCount(0x0),
fPileUpMultSelCount(0x0),
fMultTOFLowCut(0x0),
fMultTOFHighCut(0x0),
fUseTowerEq(kFALSE),
fTowerEqList(NULL),
fUseBadTowerCalib(kFALSE),
fBadTowerCalibList(NULL),
fVZEROGainEqList(NULL),
fVZEROQVecRecList(NULL),
fUseZDCSpectraCorr(kFALSE),
fZDCSpectraCorrList(NULL),
fSpectraMCList(NULL),
fBadTowerStuffList(NULL),
fVZEROStuffList(NULL),
fCachedRunNum(0),
fhZNSpectra(0x0),
fhZNSpectraCor(0x0),
fhZNSpectraPow(0x0),
fhZNBCCorr(0x0)
{
  for(int i=0; i<5; i++){
    fhZNCPM[i] = 0x0;
    fhZNAPM[i] = 0x0;
  }
  for(int i=0; i<4; i++){
    fhZNCPMQiPMC[i] = 0x0;
    fhZNAPMQiPMC[i] = 0x0;
  }
  for(Int_t r=0; r<fCRCMaxnRun; r++) {
    fRunList[r] = 0;
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for(Int_t c=0; c<100; c++) {
    fBadTowerCalibHist[c] = NULL;
  }
  fVZEROGainEqHist = NULL;
  for (Int_t k=0; k<fkVZEROnHar; k++) {
//    fVZEROQVectorRecQx[k] = NULL;
//    fVZEROQVectorRecQy[k] = NULL;
    fVZEROQVectorRecQxStored[k] = NULL;
    fVZEROQVectorRecQyStored[k] = NULL;
    for (Int_t t=0; t<fkVZEROnQAplots; t++) {
      fVZEROQVectorRecFinal[k][t] = NULL;
    }
  }
  for(Int_t i=0; i<8; i++) {
    SpecCorMu1[i] = NULL;
    SpecCorMu2[i] = NULL;
    SpecCorSi[i] = NULL;
    SpecCorAv[i] = NULL;
  }
  this->InitializeRunArrays();
  fMyTRandom3 = new TRandom3(1);
  gRandom->SetSeed(fMyTRandom3->Integer(65539));
  for(Int_t j=0; j<2; j++) {
    for(Int_t c=0; c<10; c++) {
      fPtSpecGen[j][c] = NULL;
      fPtSpecFB32[j][c] = NULL;
      fPtSpecFB96[j][c] = NULL;
      fPtSpecFB128[j][c] = NULL;
      fPtSpecFB768[j][c] = NULL;
    }
  }
  for (Int_t c=0; c<2; c++) {
    fhZNCenDis[c] = NULL;
  }
  fMinRingVZC=1;
  fMaxRingVZC=4;
  fMinRingVZA=5;
  fMaxRingVZA=8;
}

//________________________________________________________________________
AliAnalysisTaskCRCZDC::AliAnalysisTaskCRCZDC(const char *name, TString RPtype, Bool_t on, UInt_t iseed, Bool_t bCandidates):
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
fRejectPileUpTight(kFALSE),
fResetNegativeZDC(kFALSE),
fCentrLowLim(0.),
fCentrUpLim(100.),
fCentrEstimator(kV0M),
fOutput(0x0),
fhZNCvsZNA(0x0),
fhZDCCvsZDCCA(0x0),
fhZNCvsZPC(0x0),
fhZNAvsZPA(0x0),
fhZNvsZP(0x0),
fhZNvsVZERO(0x0),
fhZDCvsVZERO(0x0),
fhZDCvsTracklets(0x0),
fhZDCvsNclu1(0x0),
fhDebunch(0x0),
fhAsymm(0x0),
fhZNAvsAsymm(0x0),
fhZNCvsAsymm(0x0),
fhZNCvscentrality(0x0),
fhZNAvscentrality(0x0),
fDataSet(kAny),
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
fVZEROMult(0x0),
fMultSelection(0x0),
fPileUpCount(0x0),
fPileUpMultSelCount(0x0),
fMultTOFLowCut(0x0),
fMultTOFHighCut(0x0),
fUseTowerEq(kFALSE),
fTowerEqList(NULL),
fUseBadTowerCalib(kFALSE),
fBadTowerCalibList(NULL),
fVZEROGainEqList(NULL),
fVZEROQVecRecList(NULL),
fUseZDCSpectraCorr(kFALSE),
fZDCSpectraCorrList(NULL),
fSpectraMCList(NULL),
fBadTowerStuffList(NULL),
fVZEROStuffList(NULL),
fCachedRunNum(0),
fhZNSpectra(0x0),
fhZNSpectraCor(0x0),
fhZNSpectraPow(0x0),
fhZNBCCorr(0x0)
{
  
  for(int i=0; i<5; i++){
    fhZNCPM[i] = 0x0;
    fhZNAPM[i] = 0x0;
  }
  for(int i=0; i<4; i++){
    fhZNCPMQiPMC[i] = 0x0;
    fhZNAPMQiPMC[i] = 0x0;
  }
  for(Int_t r=0; r<fCRCMaxnRun; r++) {
    fRunList[r] = 0;
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for(Int_t c=0; c<100; c++) {
    fBadTowerCalibHist[c] = NULL;
  }
  fVZEROGainEqHist = NULL;
  for (Int_t k=0; k<fkVZEROnHar; k++) {
    //    fVZEROQVectorRecQx[k] = NULL;
    //    fVZEROQVectorRecQy[k] = NULL;
    fVZEROQVectorRecQxStored[k] = NULL;
    fVZEROQVectorRecQyStored[k] = NULL;
    for (Int_t t=0; t<fkVZEROnQAplots; t++) {
      fVZEROQVectorRecFinal[k][t] = NULL;
    }
  }
  for(Int_t i=0; i<8; i++) {
    SpecCorMu1[i] = NULL;
    SpecCorMu2[i] = NULL;
    SpecCorSi[i] = NULL;
    SpecCorAv[i] = NULL;
  }
  this->InitializeRunArrays();
  fMyTRandom3 = new TRandom3(iseed);
  gRandom->SetSeed(fMyTRandom3->Integer(65539));
  
  DefineInput(0, TChain::Class());
  // Define output slots here
  // Define here the flow event output
  DefineOutput(1, AliFlowEventSimple::Class());
  DefineOutput(2, TList::Class());
  
  for(Int_t j=0; j<2; j++) {
    for(Int_t c=0; c<10; c++) {
      fPtSpecGen[j][c] = NULL;
      fPtSpecFB32[j][c] = NULL;
      fPtSpecFB96[j][c] = NULL;
      fPtSpecFB128[j][c] = NULL;
      fPtSpecFB768[j][c] = NULL;
    }
  }
  for (Int_t c=0; c<2; c++) {
    fhZNCenDis[c] = NULL;
  }
  fMinRingVZC=1;
  fMaxRingVZC=4;
  fMinRingVZA=5;
  fMaxRingVZA=8;
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
  if (fTowerEqList)  delete fTowerEqList;
  if (fBadTowerCalibList) delete fBadTowerCalibList;
  if (fVZEROGainEqList) delete fVZEROGainEqList;
  if (fVZEROQVecRecList) delete fVZEROQVecRecList;
  if (fZDCSpectraCorrList) delete fZDCSpectraCorrList;
  if (fAnalysisUtil) delete fAnalysisUtil;
  if (fQAList) delete fQAList;
  if (fCutContainer) fCutContainer->Delete(); delete fCutContainer;
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::InitializeRunArrays()
{
  for(Int_t r=0;r<fCRCMaxnRun;r++) {
    fCRCQVecListRun[r] = NULL;
    for(Int_t k=0;k<fCRCnTow;k++) {
      fZNCTower[r][k] = NULL;
      fZNATower[r][k] = NULL;
    }
//    fhZNSpectraRbR[r] = NULL;
  }
  //   for(Int_t i=0;i<fnCen;i++) {
  //     fPtPhiEtaRbRFB128[r][i] = NULL;
  //     fPtPhiEtaRbRFB768[r][i] = NULL;
  //   }
  // }
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
  
  fVZEROStuffList = new TList();
  fVZEROStuffList->SetOwner(kTRUE);
  fVZEROStuffList->SetName("VZERO stuff");
  fOutput->Add(fVZEROStuffList);
  
  fBadTowerStuffList = new TList();
  fBadTowerStuffList->SetOwner(kTRUE);
  fBadTowerStuffList->SetName("BadTowerCalib");
  fOutput->Add(fBadTowerStuffList);
  
  fCenDis = new TH1F("fCenDis", "fCenDis", 100, 0., 100.);
  fOutput->Add(fCenDis);
  fPileUpCount = new TH1F("fPileUpCount", "fPileUpCount", 9, 0., 9.);
  fPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultiplicityComb08");
  fPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>7.5");
  fPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(8,"multESDTPCDif");
  fPileUpCount->GetXaxis()->SetBinLabel(9,"extraPileUpMultSel");
  fOutput->Add(fPileUpCount);
  fPileUpMultSelCount = new TH1F("fPileUpMultSelCount", "fPileUpMultSelCount", 8, 0., 8.);
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(1,"IsNotPileup");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(2,"IsNotPileupMV");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(3,"IsNotPileupInMultBins");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(4,"InconsistentVertices");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(6,"AsymmetricInVZERO");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(7,"IncompleteDAQ");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(8,"GoodVertex2016");
  fOutput->Add(fPileUpMultSelCount);
  
  fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
  fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
  fOutput->Add(fMultTOFLowCut);
  fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
  fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
  fOutput->Add(fMultTOFHighCut);
  
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] = new TH1D();
      fOutput->Add(fTowerGainEq[c][i]);
    }
  }
  
  for (Int_t c=0; c<2; c++) {
    fhZNCenDis[c] = new TH3D(Form("fhZNCenDis[%d]",c), Form("fhZNCenDis[%d]",c), 100, 0., 100., 100, -2., 2. , 100., -2., 2.);
    fOutput->Add(fhZNCenDis[c]);
  }
  
  if(fBadTowerCalibList) {
    for(Int_t c=0; c<100; c++) {
      fBadTowerCalibHist[c] = (TH2D*)fBadTowerCalibList->FindObject(Form("TH2Resp[%d]",c));
      fBadTowerStuffList->Add(fBadTowerCalibHist[c]);
    }
  }
  if(fZDCSpectraCorrList) {
    for(Int_t i=0; i<8; i++) {
      SpecCorMu1[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorMu1[%d]",i));
      fOutput->Add(SpecCorMu1[i]);
      SpecCorMu2[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorMu2[%d]",i));
      fOutput->Add(SpecCorMu2[i]);
      SpecCorAv[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorAv[%d]",i));
      fOutput->Add(SpecCorAv[i]);
      SpecCorSi[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorSi[%d]",i));
      fOutput->Add(SpecCorSi[i]);
    }
  }
  
  fhZNSpectra = new TH3D("fhZNSpectra","fhZNSpectra",100,0.,100.,8,0.,8.,1000,0.,1.E5);
  fOutput->Add(fhZNSpectra);
  fhZNSpectraCor = new TH3D("fhZNSpectraCor","fhZNSpectraCor",100,0.,100.,8,0.,8.,1000,0.,1.E5);
  fOutput->Add(fhZNSpectraCor);
  fhZNSpectraPow = new TH3D("fhZNSpectraPow","fhZNSpectraPow",100,0.,100.,8,0.,8.,1000,0.,TMath::Power(1.E5,fZDCGainAlpha));
  fOutput->Add(fhZNSpectraPow);
  fhZNBCCorr = new TH3D("fhZNBCCorr","fhZNBCCorr",100,0.,100.,500,0.,1.E5,500,0.,1.E5);
  fOutput->Add(fhZNBCCorr);
  
  if(fAnalysisType == "MCAOD") {
    
    fSpectraMCList = new TList();
    fSpectraMCList->SetOwner(kTRUE);
    fSpectraMCList->SetName("Spectra");
    fOutput->Add(fSpectraMCList);
    
    Double_t xmin[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.33,2.66,3.,3.5,4.,5.,6.,9.,20.};
    for(Int_t j=0; j<2; j++) {
      for(Int_t c=0; c<10; c++) {
        fPtSpecGen[j][c] = new TH1F(Form("fPtSpecGen[%d][%d]",j,c), Form("fPtSpecGen[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecGen[j][c]);
        fPtSpecFB32[j][c] = new TH1F(Form("fPtSpecFB32[%d][%d]",j,c), Form("fPtSpecFB32[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB32[j][c]);
        fPtSpecFB96[j][c] = new TH1F(Form("fPtSpecFB96[%d][%d]",j,c), Form("fPtSpecFB96[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB96[j][c]);
        fPtSpecFB128[j][c] = new TH1F(Form("fPtSpecFB128[%d][%d]",j,c), Form("fPtSpecFB128[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB128[j][c]);
        fPtSpecFB768[j][c] = new TH1F(Form("fPtSpecFB768[%d][%d]",j,c), Form("fPtSpecFB768[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB768[j][c]);
      }
    }
  }
  
  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);
  
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
    }
  }
  
  fhZNCvsZNA = new TH2F("hZNCvsZNA","hZNCvsZNA",200,-50.,140000,200,-50.,140000);
  fOutput->Add(fhZNCvsZNA);
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
  
  //********************************************************************
  
  Int_t dRun10h[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};
  
  Int_t dRun11h[] = {167902, 167903, 167915, 167920, 167985, 167987, 167988, 168066, 168068, 168069, 168076, 168104, 168105, 168107, 168108, 168115, 168212, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168461, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168984, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169143, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 169965, 170027, 170036,170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593};
  
  // 12 low IR: 244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392
  // 80 high IR ("CentralBarrelTracking" good runs): 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246676, 246675, 246540, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683
  
  Int_t dRun15o[] = {244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};
  
  Int_t dRun15ov6[] = {244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};
  
  if(fDataSet==k2010) {fCRCnRun=92;}
  if(fDataSet==k2011) {fCRCnRun=119;}
  if(fDataSet==k2015) {fCRCnRun=90;}
  if(fDataSet==k2015v6) {fCRCnRun=91;}
  if(fDataSet==kAny) {fCRCnRun=1;}
  
  Int_t d=0;
  for(Int_t r=0; r<fCRCnRun; r++) {
    if(fDataSet==k2010)   fRunList[d] = dRun10h[r];
    if(fDataSet==k2011)   fRunList[d] = dRun11h[r];
    if(fDataSet==k2015)   fRunList[d] = dRun15o[r];
    if(fDataSet==k2015v6) fRunList[d] = dRun15ov6[r];
    if(fDataSet==kAny) fRunList[d] = 1;
    d++;
  }
  
  fVZEROMult = new TProfile2D("fVZEROMult","fVZEROMult",fCRCnRun,0.,1.*fCRCnRun,64,0.,64.);
  for (Int_t i=0; i<fCRCnRun; i++) {
    fVZEROMult->GetXaxis()->SetBinLabel(i+1,Form("%d",fRunList[i]));
  }
  fVZEROStuffList->Add(fVZEROMult);
  
  if(fVZEROGainEqList) {
    fVZEROGainEqHist = (TH2D*)fVZEROGainEqList->FindObject("VZEROEqGain");
    fVZEROStuffList->Add(fVZEROGainEqHist);
  }
  if(fVZEROQVecRecList) {
    for (Int_t k=0; k<fkVZEROnHar; k++) {
      fVZEROQVectorRecQxStored[k] = (TProfile3D*)fVZEROQVecRecList->FindObject(Form("fVZEROQVectorRecQx[%d]",k));
      fVZEROQVectorRecQxStored[k]->SetTitle(Form("fVZEROQVectorRecQxStored[%d]",k));
      fVZEROQVectorRecQxStored[k]->SetName(Form("fVZEROQVectorRecQxStored[%d]",k));
      fVZEROStuffList->Add(fVZEROQVectorRecQxStored[k]);
      fVZEROQVectorRecQyStored[k] = (TProfile3D*)fVZEROQVecRecList->FindObject(Form("fVZEROQVectorRecQy[%d]",k));
      fVZEROQVectorRecQyStored[k]->SetTitle(Form("fVZEROQVectorRecQyStored[%d]",k));
      fVZEROQVectorRecQyStored[k]->SetName(Form("fVZEROQVectorRecQyStored[%d]",k));
      fVZEROStuffList->Add(fVZEROQVectorRecQyStored[k]);
      for (Int_t t=0; t<fkVZEROnQAplots; t++) {
        fVZEROQVectorRecFinal[k][t] = new TProfile2D(Form("fVZEROQVectorRecFinal[%d][%d]",k,t),Form("fVZEROQVectorRecFinal[%d][%d]",k,t),fCRCnRun,0.,1.*fCRCnRun,100,0.,100.,"s");
        fVZEROQVectorRecFinal[k][t]->Sumw2();
        fVZEROStuffList->Add(fVZEROQVectorRecFinal[k][t]);
      }
    }
  }
  
//  for (Int_t k=0; k<fkVZEROnHar; k++) {
//    fVZEROQVectorRecQx[k] = new TProfile3D(Form("fVZEROQVectorRecQx[%d]",k),Form("fVZEROQVectorRecQx[%d]",k),fCRCnRun,0.,1.*fCRCnRun,100,0.,100.,8,0.,8.,"s");
//    fVZEROQVectorRecQx[k]->Sumw2();
//    fVZEROStuffList->Add(fVZEROQVectorRecQx[k]);
//    fVZEROQVectorRecQy[k] = new TProfile3D(Form("fVZEROQVectorRecQy[%d]",k),Form("fVZEROQVectorRecQy[%d]",k),fCRCnRun,0.,1.*fCRCnRun,100,0.,100.,8,0.,8.,"s");
//    fVZEROQVectorRecQy[k]->Sumw2();
//    fVZEROStuffList->Add(fVZEROQVectorRecQy[k]);
//  }
  
  Double_t ptmin[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.8,2.2,3.,4.,6.,8.,12.,20.};
  Double_t phimin[] = {0.,TMath::Pi()/8.,2*TMath::Pi()/8.,3*TMath::Pi()/8.,4*TMath::Pi()/8.,5*TMath::Pi()/8.,6*TMath::Pi()/8.,7*TMath::Pi()/8.,8*TMath::Pi()/8.,9*TMath::Pi()/8.,10*TMath::Pi()/8.,11*TMath::Pi()/8.,12*TMath::Pi()/8.,13*TMath::Pi()/8.,14*TMath::Pi()/8.,15*TMath::Pi()/8.,16*TMath::Pi()/8.};
  Double_t etamin[] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};
  
  for(Int_t r=0;r<fCRCnRun;r++) {
    fCRCQVecListRun[r] = new TList();
    fCRCQVecListRun[r]->SetName(Form("Run %d",fRunList[r]));
    fCRCQVecListRun[r]->SetOwner(kTRUE);
    fOutput->Add(fCRCQVecListRun[r]);
    
    if(!fUseTowerEq) {
      for(Int_t k=0;k<fCRCnTow;k++) {
        fZNCTower[r][k] = new TProfile(Form("fZNCTower[%d][%d]",fRunList[r],k),Form("fZNCTower[%d][%d]",fRunList[r],k),100,0.,100.,"s");
        fZNCTower[r][k]->Sumw2();
        fCRCQVecListRun[r]->Add(fZNCTower[r][k]);
        fZNATower[r][k] = new TProfile(Form("fZNATower[%d][%d]",fRunList[r],k),Form("fZNATower[%d][%d]",fRunList[r],k),100,0.,100.,"s");
        fZNATower[r][k]->Sumw2();
        fCRCQVecListRun[r]->Add(fZNATower[r][k]);
      }
    }
    
//    fhZNSpectraRbR[r] = new TH3D(Form("fhZNSpectraRbR[%d]",fRunList[r]),Form("fhZNSpectraRbR[%d]",fRunList[r]),50,0.,100.,8,0.,8.,100,0.,1.E5);
//    fCRCQVecListRun[r]->Add(fhZNSpectraRbR[r]);
    
    //   for(Int_t i=0;i<fnCen;i++) {
    //     fPtPhiEtaRbRFB128[r][i] = new TH3F(Form("fPtPhiEtaRbRFB128[%d][%d]",r,i),Form("fPtPhiEtaRbRFB128[%d][%d]",r,i),14, ptmin, 16, phimin, 16, etamin);
    //     fCRCQVecListRun[r]->Add(fPtPhiEtaRbRFB128[r][i]);
    //     fPtPhiEtaRbRFB768[r][i] = new TH3F(Form("fPtPhiEtaRbRFB768[%d][%d]",r,i),Form("fPtPhiEtaRbRFB768[%d][%d]",r,i),14, ptmin, 16, phimin, 16, etamin);
    //     fCRCQVecListRun[r]->Add(fPtPhiEtaRbRFB768[r][i]);
    //   }
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
  
  Int_t RunBin=-1, bin=0;
  Int_t RunNum = aod->GetRunNumber();
  for(Int_t c=0;c<fCRCnRun;c++) {
    if(fRunList[c]==RunNum) RunBin=bin;
    else bin++;
  }
  if(RunBin==-1) return;
  if(fDataSet==kAny) RunBin=0;
  
  //DEFAULT - automatically takes care of everything
  if (fAnalysisType == "AUTOMATIC") {
    
    // get centrality
    Double_t centrV0M=300, centrCL1=300, centrCL0=300, centrTRK=300;
    if(fDataSet!=k2015 && fDataSet!=k2015v6) {
      centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
      centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
      centrCL0 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
      centrTRK = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
    } else {
      fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(!fMultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
      } else {
        centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
        centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
        centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
        centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
      }
    }
    
    //check event cuts
    if (InputEvent()) {
      if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
      if(fRejectPileUp) {
        if(fDataSet!=k2015 && fDataSet!=k2015v6) {
          
          Bool_t BisPileup=kFALSE;
          
          // check anyway pileup
          if (plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }
          
          Short_t isPileup = aod->IsPileupFromSPD(3);
          if (isPileup != 0) {
            fPileUpCount->Fill(1.5);
            BisPileup=kTRUE;
          }
          
          if (((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }
          
          if (aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }
          
          // check vertex consistency
          const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
          const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();
          
          if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)  {
            fPileUpCount->Fill(5.5);
            BisPileup=kTRUE;
          }
          
          double covTrc[6], covSPD[6];
          vtTrc->GetCovarianceMatrix(covTrc);
          vtSPD->GetCovarianceMatrix(covSPD);
          
          double dz = vtTrc->GetZ() - vtSPD->GetZ();
          
          double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
          double errTrc = TMath::Sqrt(covTrc[5]);
          double nsigTot = dz/errTot;
          double nsigTrc = dz/errTrc;
          
          if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
            fPileUpCount->Fill(6.5);
            BisPileup=kTRUE;
          }
          
          if (fAnalysisUtil->IsPileUpEvent(InputEvent())) {
            fPileUpCount->Fill(7.5);
            BisPileup=kTRUE;
          }
          
//        if(BisPileup) return;
        } else {
          // pileup from AliMultSelection
          if(!fMultSelection->GetThisEventIsNotPileup()) fPileUpMultSelCount->Fill(0.5);
          if(!fMultSelection->GetThisEventIsNotPileupMV()) fPileUpMultSelCount->Fill(1.5);
          if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) fPileUpMultSelCount->Fill(2.5);
          if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) fPileUpMultSelCount->Fill(3.5);
          if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) fPileUpMultSelCount->Fill(4.5);
          if(!fMultSelection->GetThisEventIsNotAsymmetricInVZERO()) fPileUpMultSelCount->Fill(5.5);
          if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) fPileUpMultSelCount->Fill(6.5);
          if(!fMultSelection->GetThisEventHasGoodVertex2016()) fPileUpMultSelCount->Fill(7.5);
          
          Bool_t BisPileup=kFALSE;
          
          // pile-up a la Dobrin for LHC15o
          if (plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }
          
          Short_t isPileup = aod->IsPileupFromSPD(3);
          if (isPileup != 0) {
            fPileUpCount->Fill(1.5);
            BisPileup=kTRUE;
          }
          
          if (((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }
          
          if (aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }
          
          if(fabs(centrV0M-centrCL1)>7.5)  {
            fPileUpCount->Fill(4.5);
            BisPileup=kTRUE;
          }
          
          // check vertex consistency
          const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
          const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();
          
          if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)  {
            fPileUpCount->Fill(5.5);
            BisPileup=kTRUE;
          }
          
          double covTrc[6], covSPD[6];
          vtTrc->GetCovarianceMatrix(covTrc);
          vtSPD->GetCovarianceMatrix(covSPD);
          
          double dz = vtTrc->GetZ() - vtSPD->GetZ();
          
          double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
          double errTrc = TMath::Sqrt(covTrc[5]);
          double nsigTot = dz/errTot;
          double nsigTrc = dz/errTrc;
          
          if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
            fPileUpCount->Fill(6.5);
            BisPileup=kTRUE;
          }
          
          // cuts on tracks
          const Int_t nTracks = aod->GetNumberOfTracks();
          Int_t multEsd = ((AliAODHeader*)aod->GetHeader())->GetNumberOfESDTracks();

          Int_t multTrk = 0;
          Int_t multTrkBefC = 0;
          Int_t multTrkTOFBefC = 0;
          Int_t multTPC = 0;
          
          for (Int_t it = 0; it < nTracks; it++) {
            
            AliAODTrack* aodTrk = (AliAODTrack*)aod->GetTrack(it);
            if (!aodTrk){
              delete aodTrk;
              continue;
            }
            
//            if (aodTrk->TestFilterBit(32)){
//              multTrkBefC++;
//              
//              if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
//                multTrkTOFBefC++;
//              
//              if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
//                multTrk++;
//            }
            
            if (aodTrk->TestFilterBit(128))
              multTPC++;
            
          } // end of for (Int_t it = 0; it < nTracks; it++)
          
          Double_t multTPCn = multTPC;
          Double_t multEsdn = multEsd;
          Double_t multESDTPCDif = multEsdn - multTPCn*3.38;
          
          if (multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
            fPileUpCount->Fill(7.5);
            BisPileup=kTRUE;
          }
          
          if(fRejectPileUpTight) {
            if(BisPileup==kFALSE) {
              if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
              if(BisPileup) fPileUpCount->Fill(8.5);
            }
          }
          
          if(BisPileup) return;
        }
      }
    }
    
    //first attach all possible information to the cuts
    fCutsRP->SetEvent( InputEvent(), MCEvent() );  //attach event
    fCutsPOI->SetEvent( InputEvent(), MCEvent() );
    
    //then make the event
    fFlowEvent->Fill( fCutsRP, fCutsPOI );
    
    fFlowEvent->SetReferenceMultiplicity(fCutsEvent->GetReferenceMultiplicity(InputEvent(),McEvent));
    
    if(fCentrEstimator==kV0M) fFlowEvent->SetCentrality(centrV0M);
    if(fCentrEstimator==kCL0) fFlowEvent->SetCentrality(centrCL0);
    if(fCentrEstimator==kCL1) fFlowEvent->SetCentrality(centrCL1);
    if(fCentrEstimator==kTRK) fFlowEvent->SetCentrality(centrTRK);
    fFlowEvent->SetCentralityCL1(centrCL1);
    fFlowEvent->SetCentralityTRK(centrTRK);
    //   fFlowEvent->SetNITSCL1(((AliVAODHeader*)aod->GetHeader())->GetNumberOfITSClusters(1));
    Double_t SumV0=0.;
    for(Int_t i=0; i<64; i++) {
      if(std::isfinite(aod->GetVZEROEqMultiplicity(i))) SumV0 += aod->GetVZEROEqMultiplicity(i);
    }
    fFlowEvent->SetNITSCL1(SumV0);
    
    Double_t vtxpos[3]={0.,0.,0.};
    vtxpos[0] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetX();
    vtxpos[1] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetY();
    vtxpos[2] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetZ();
    fFlowEvent->SetVertexPosition(vtxpos);
    
    if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
    
    // run-by-run QA
    //   for(Int_t jTracks = 0; jTracks<aod->GetNumberOfTracks(); jTracks++){
    //     AliAODTrack* track = (AliAODTrack*)aod->GetTrack(jTracks);
    //     if(!track) continue;
    //     // general kinematic & quality cuts
    //     if (track->Pt() < .2 || track->Pt() > 20. || TMath::Abs(track->Eta()) > .8 || track->GetTPCNcls() < 70)  continue;
    //     if (track->TestFilterBit(128)) fPtPhiEtaRbRFB128[RunBin][CenBin]->Fill(track->Pt(),track->Phi(),track->Eta());
    //     if (track->TestFilterBit(768)) fPtPhiEtaRbRFB768[RunBin][CenBin]->Fill(track->Pt(),track->Phi(),track->Eta());
    //   }
    fCenDis->Fill(centrV0M);
    
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
    Double_t centr = 300;
    if(fDataSet==k2015 || fDataSet==k2015v6) {
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
    // centrality bin
    if (centr<fCentrLowLim || centr>=fCentrUpLim ) return;
    Int_t CenBin = -1;
    CenBin = GetCenBin(centr);
    if(CenBin==-1) return;
    
    // reconstructed
    for(Int_t jTracks = 0; jTracks<aod->GetNumberOfTracks(); jTracks++){
      
      AliAODTrack* track = (AliAODTrack*)aod->GetTrack(jTracks);
      if(!track) continue;
      
      // to select primaries
      Int_t lp = TMath::Abs(track->GetLabel());
      
      // general kinematic cuts
      if (track->Pt() < .2 || track->Pt() > 20. || TMath::Abs(track->Eta()) > 0.8) continue;
      
      // cut on DCA
      Double_t DCAxy = track->DCA();
      Double_t DCAz = track->ZAtDCA();
      if(fabs(DCAxy)>2.4 || fabs(DCAz)>3.2) continue;
      
      // various cuts on TPC clusters
      if (track->GetTPCNcls() < 70) continue;
      Double_t chi2_per_tpc = track->Chi2perNDF();
      if (chi2_per_tpc < 0.1 || chi2_per_tpc > 4.) continue;
      Double_t fraction_shared_tpccls = 1.*track->GetTPCnclsS()/track->GetTPCncls();
      if (fraction_shared_tpccls > 0.4) continue;
      
      // test filter bits
      if (((AliAODMCParticle*)fStack->At(lp))->IsPhysicalPrimary()) {
        if (track->TestFilterBit(32))  fPtSpecFB32[0][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(96))  fPtSpecFB96[0][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(128)) fPtSpecFB128[0][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(768)) fPtSpecFB768[0][CenBin]->Fill(track->Pt());
      } else {
        if (track->TestFilterBit(32))  fPtSpecFB32[1][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(96))  fPtSpecFB96[1][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(128)) fPtSpecFB128[1][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(768)) fPtSpecFB768[1][CenBin]->Fill(track->Pt());
      }
      
      fCenDis->Fill(centr);
    }
    
    // generated (physical primaries)
    
    for(Int_t jTracks = 0; jTracks<fStack->GetEntriesFast(); jTracks++) {
      AliAODMCParticle *MCpart = (AliAODMCParticle*)fStack->At(jTracks);
      if (!MCpart) {
        printf("ERROR: Could not receive MC track %d\n", jTracks);
        continue;
      }
      // kinematic cuts
      if ( MCpart->Pt() < .2 || MCpart->Pt() > 20. || TMath::Abs(MCpart->Eta()) > .8 ) continue;
      // select charged primaries
      if ( MCpart->Charge() == 0. || !MCpart->IsPhysicalPrimary()) continue;
      
      fPtSpecGen[0][CenBin]->Fill(MCpart->Pt());
    }
    
    //    fGenHeader = McEvent->GenEventHeader();
    //    if(fGenHeader) fPythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(fGenHeader);
    //  printf("#reconstructed : %d (rejected from cuts %d), #MC primaries : %d (rejected from cuts %d) \n",AODPOIs,AODbads,MCPrims,MCSecos);
    fFlowEvent->SetReferenceMultiplicity(aod->GetNumberOfTracks());
    fFlowEvent->SetCentrality(centr);
    //    if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
    fFlowEvent->SetRun(RunNum);
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
      
      fPtSpecGen[0][CenBin]->Fill(particle->Pt());
      
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
  
  //********************************************************************************************************************************
  
  if(fAnalysisType == "AOD" || fAnalysisType == "AUTOMATIC") {
    
    // PHYSICS SELECTION
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *hdr = (AliInputEventHandler*)am->GetInputEventHandler();
    
    if(hdr->IsEventSelected() && AliVEvent::kAny) {
      
      Double_t centrperc = fFlowEvent->GetCentrality();
      Int_t cenb = (Int_t)centrperc;
      
      AliAODTracklets *trackl = aod->GetTracklets();
      Int_t nTracklets = trackl->GetNumberOfTracklets();
      
      // VZERO
      
      // get VZERO data
      AliAODVZERO *vzeroAOD = aod->GetVZEROData();
      Double_t multV0A = vzeroAOD->GetMTotV0A();
      Double_t multV0C = vzeroAOD->GetMTotV0C();
      Int_t CachednRing = 1;
      Double_t QxTot[fkVZEROnHar] = {0.}, QyTot[fkVZEROnHar] = {0.};
      Double_t denom = 0.;
      Double_t V0TotQC[fkVZEROnHar][2] = {{0.}}, V0TotQA[fkVZEROnHar][2] = {{0.}};
      Double_t MultC[fkVZEROnHar] = {0.}, MultA[fkVZEROnHar] = {0.};
      
      for(Int_t i=0; i<64; i++) {
        
        // correct multiplicity per channel
        Double_t mult = vzeroAOD->GetMultiplicity(i);
        if(fVZEROGainEqHist) {
          Double_t EqFactor = fVZEROGainEqHist->GetBinContent(RunBin+1,i+1);
          if(EqFactor>0.) mult *= EqFactor;
        }
        fVZEROMult->Fill(RunBin+0.5,i+0.5,mult);
        
        // build Q-vector per ring
        Int_t nRing = (Int_t)i/8 + 1;
        Double_t ChPhi = TMath::PiOver4()*(0.5+i%8);
        
        if(i == 63) {
          for (Int_t k=0; k<fkVZEROnHar; k++) {
            QxTot[k] += mult*TMath::Cos((k+1.)*ChPhi);
            QyTot[k] += mult*TMath::Sin((k+1.)*ChPhi);
          }
          denom += mult;
          nRing++;
        }
        
        if(nRing!=CachednRing) {
          for (Int_t k=0; k<fkVZEROnHar; k++) {
            Double_t QxRec = QxTot[k]/denom;
            Double_t QyRec = QyTot[k]/denom;
            // store values for re-centering
//            fVZEROQVectorRecQx[k]->Fill(RunBin+0.5,centrperc,CachednRing-0.5,QxRec);
//            fVZEROQVectorRecQy[k]->Fill(RunBin+0.5,centrperc,CachednRing-0.5,QyRec);
            // do re-centering
            if(fVZEROQVectorRecQxStored[k]) {
              if(!std::isnan(fVZEROQVectorRecQxStored[k]->GetBinContent(fVZEROQVectorRecQxStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5)))) QxRec -= fVZEROQVectorRecQxStored[k]->GetBinContent(fVZEROQVectorRecQxStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5));
              if(!std::isnan(fVZEROQVectorRecQyStored[k]->GetBinContent(fVZEROQVectorRecQyStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5)))) QyRec -= fVZEROQVectorRecQyStored[k]->GetBinContent(fVZEROQVectorRecQyStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5));
            }
            // sum of Q-vectors over all rings (total V0 Q-vector)
            if (CachednRing >= fMinRingVZC && CachednRing <= fMaxRingVZC) {
              V0TotQC[k][0] += QxRec*denom;
              V0TotQC[k][1] += QyRec*denom;
              MultC[k] += denom;
            }
            if (CachednRing >= fMinRingVZA && CachednRing <= fMaxRingVZA) {
              V0TotQA[k][0] += QxRec*denom;
              V0TotQA[k][1] += QyRec*denom;
              MultA[k] += denom;
            }
            QxTot[k] = 0.;
            QyTot[k] = 0.;
          }
          denom = 0.;
          CachednRing = nRing;
        }
        for (Int_t k=0; k<fkVZEROnHar; k++) {
          QxTot[k] += mult*TMath::Cos((k+1.)*ChPhi);
          QyTot[k] += mult*TMath::Sin((k+1.)*ChPhi);
        }
        denom += mult;
      }
      
      for (Int_t k=0; k<fkVZEROnHar; k++) {
        if(MultC[k]>0. && MultA[k]>0.) {
          Double_t QCx = V0TotQC[k][0]/MultC[k], QCy = V0TotQC[k][1]/MultC[k], QAx = V0TotQA[k][0]/MultA[k], QAy = V0TotQA[k][1]/MultA[k];
          if(!std::isnan(QCx) && !std::isnan(QCy) && !std::isnan(QAx) && !std::isnan(QAy)) {
            fFlowEvent->SetV02Qsub(QCx,QCy,MultC[k],QAx,QAy,MultA[k],k+1);
            fVZEROQVectorRecFinal[k][0]->Fill(RunBin+0.5,centrperc,QCx);
            fVZEROQVectorRecFinal[k][1]->Fill(RunBin+0.5,centrperc,QCy);
            fVZEROQVectorRecFinal[k][2]->Fill(RunBin+0.5,centrperc,QAx);
            fVZEROQVectorRecFinal[k][3]->Fill(RunBin+0.5,centrperc,QAy);
            fVZEROQVectorRecFinal[k][4]->Fill(RunBin+0.5,centrperc,QCx*QAx);
            fVZEROQVectorRecFinal[k][5]->Fill(RunBin+0.5,centrperc,QCy*QAy);
            fVZEROQVectorRecFinal[k][6]->Fill(RunBin+0.5,centrperc,QCx*QAy);
            fVZEROQVectorRecFinal[k][7]->Fill(RunBin+0.5,centrperc,QCy*QAx);
          } else {
            fFlowEvent->SetV02Qsub(0.,0.,0.,0.,0.,0.,k+1);
          }
        } else {
          fFlowEvent->SetV02Qsub(0.,0.,0.,0.,0.,0.,k+1);
        }
      }
    
    
    
//      AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(aodEvent->FindListObject("Forward"));
//      const TH2D& d2Ndetadphi = aodForward->GetHistogram();
//      Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
//      Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
//      Double_t ret = 0.;
//      // Loop over eta
//      for (Int_t iEta = 1; iEta <= nEta; iEta++) {
//        Int_t valid = d2Ndetadphi.GetBinContent(iEta, 0);
//        if (!valid) continue; // No data expected for this eta
//        // Loop over phi
//        for (Int_t iPhi = 1; i <= nPhi; i++) {
//          ret = d2Ndetadphi.GetBinContent(iEta, iPhi);
//          printf("eta %e phi %e : %e \n",d2Ndetadphi.GetXaxis()->GetBinCenter(iEta),d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi),ret);
//        }
//      }
      
      AliAODZDC *aodZDC = aod->GetZDCData();
      
      Double_t energyZNC  = (Double_t) (aodZDC->GetZNCEnergy());
      Double_t energyZPC  = (Double_t) (aodZDC->GetZPCEnergy());
      Double_t energyZNA  = (Double_t) (aodZDC->GetZNAEnergy());
      Double_t energyZPA  = (Double_t) (aodZDC->GetZPAEnergy());
      Double_t energyZEM1 = (Double_t) (aodZDC->GetZEM1Energy());
      Double_t energyZEM2 = (Double_t) (aodZDC->GetZEM2Energy());
      
      const Double_t * towZNCraw = aodZDC->GetZNCTowerEnergy();
      const Double_t * towZNAraw = aodZDC->GetZNATowerEnergy();
      
      // Get centroid from ZDCs *******************************************************
      
      Double_t Enucl = (RunNum < 209122 ? 1380. : 2511.);
      Double_t xyZNC[2]={0.,0.}, xyZNA[2]={0.,0.};
      Double_t towZNC[5]={0.}, towZNA[5]={0.};
      
      Double_t ZNCcalib=1., ZNAcalib=1.;
      if(fUseTowerEq) {
        if(RunNum!=fCachedRunNum) {
          for(Int_t i=0; i<5; i++) {
            fTowerGainEq[0][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNCTower[%d][%d]",RunNum,i)));
            fTowerGainEq[1][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNATower[%d][%d]",RunNum,i)));
          }
        }
        for(Int_t i=0; i<5; i++) {
          if(fTowerGainEq[0][i]) towZNC[i] = towZNCraw[i]*fTowerGainEq[0][i]->GetBinContent(fTowerGainEq[0][i]->FindBin(centrperc));
          if(fTowerGainEq[1][i]) towZNA[i] = towZNAraw[i]*fTowerGainEq[1][i]->GetBinContent(fTowerGainEq[1][i]->FindBin(centrperc));
          if(fResetNegativeZDC) {
            if(towZNC[i]<0.) towZNC[i] = 0.;
            if(towZNA[i]<0.) towZNA[i] = 0.;
          }
        }
      } else {
        for(Int_t i=0; i<5; i++) {
          towZNC[i] = towZNCraw[i];
          towZNA[i] = towZNAraw[i];
          if(fResetNegativeZDC) {
            if(towZNC[i]<0.) towZNC[i] = 0.;
            if(towZNA[i]<0.) towZNA[i] = 0.;
          }
          fZNCTower[RunBin][i]->Fill(centrperc,towZNC[i]);
          fZNATower[RunBin][i]->Fill(centrperc,towZNA[i]);
        }
      }
      
      if(RunNum>=245829) towZNA[2] = 0.;
      Double_t zncEnergy=0., znaEnergy=0.;
      for(Int_t i=0; i<5; i++){
        zncEnergy += towZNC[i];
        znaEnergy += towZNA[i];
      }
      if(RunNum>=245829) znaEnergy *= 8./7.;
      fFlowEvent->SetZNCEnergy(towZNC[0]);
      fFlowEvent->SetZNAEnergy(towZNA[0]);
      
      const Double_t x[4] = {-1.75, 1.75, -1.75, 1.75};
      const Double_t y[4] = {-1.75, -1.75, 1.75, 1.75};
      Double_t numXZNC=0., numYZNC=0., denZNC=0., cZNC, wZNC, EZNC, SumEZNC=0.;
      Double_t numXZNA=0., numYZNA=0., denZNA=0., cZNA, wZNA, EZNA, SumEZNA=0., BadChOr;
      Bool_t fAllChONZNC=kTRUE, fAllChONZNA=kTRUE;
      
      if (fUseMCCen) {
        for(Int_t i=0; i<4; i++){
          
          // get energy
          EZNC = towZNC[i+1];
          fhZNSpectra->Fill(centrperc,i+0.5,EZNC);
//          fhZNSpectraRbR[RunBin]->Fill(centrperc,i+0.5,EZNC);
          if(fUseZDCSpectraCorr && EZNC>0.) {
            Double_t mu1 = SpecCorMu1[i]->Interpolate(centrperc);
            Double_t mu2 = SpecCorMu2[i]->Interpolate(centrperc);
            Double_t av = SpecCorAv[i]->Interpolate(centrperc);
            Double_t cor1 = SpecCorSi[i]->Interpolate(centrperc);
            EZNC = exp( (log(EZNC) - mu1 + mu2*cor1)/cor1 ) + av;
            fhZNSpectraCor->Fill(centrperc,i+0.5,EZNC);
          }
          if(fUseZDCSpectraCorr && EZNC<=0.) fAllChONZNC=kFALSE;
          
          SumEZNC += EZNC;
          
          // build centroid
          wZNC = TMath::Power(EZNC, fZDCGainAlpha);
          numXZNC += x[i]*wZNC;
          numYZNC += y[i]*wZNC;
          denZNC += wZNC;
          fhZNSpectraPow->Fill(centrperc,i+0.5,wZNC);
          
          // get energy
          if(i==1) {
            EZNA = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
            if(fUseBadTowerCalib && fBadTowerCalibHist[cenb]) {
              EZNA = GetBadTowerResp(EZNA, fBadTowerCalibHist[cenb]);
            }
          } else {
            EZNA = towZNA[i+1];
          }
          fhZNSpectra->Fill(centrperc,i+4.5,EZNA);
//          fhZNSpectraRbR[RunBin]->Fill(centrperc,i+4.5,EZNA);
          if(fUseZDCSpectraCorr && EZNA>0.) {
            Double_t mu1 = SpecCorMu1[i+4]->Interpolate(centrperc);
            Double_t mu2 = SpecCorMu2[i+4]->Interpolate(centrperc);
            Double_t av = SpecCorAv[i+4]->Interpolate(centrperc);
            Double_t cor1 = SpecCorSi[i+4]->Interpolate(centrperc);
            EZNA = exp( (log(EZNA) - mu1 + mu2*cor1)/cor1 ) + av;
            fhZNSpectraCor->Fill(centrperc,i+4.5,EZNA);
          }
          if(fUseZDCSpectraCorr && EZNA<=0.) fAllChONZNA=kFALSE;
          SumEZNA += EZNA;
          
          // build centroid
          wZNA = TMath::Power(EZNA, fZDCGainAlpha);
          numXZNA += x[i]*wZNA;
          numYZNA += y[i]*wZNA;
          denZNA += wZNA;
          fhZNSpectraPow->Fill(centrperc,i+4.5,wZNA);
        }
        // store distribution for unfolding
        if(RunNum<245829) {
          Double_t recoE = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
          Double_t trueE = towZNA[2];
          fhZNBCCorr->Fill(centrperc,trueE,recoE);
        }
        if(denZNC>0.){
         Double_t nSpecnC = SumEZNC/Enucl;
         cZNC = 1.89358-0.71262/(nSpecnC+0.71789);
         xyZNC[0] = cZNC*numXZNC/denZNC;
         xyZNC[1] = cZNC*numYZNC/denZNC;
          denZNC *= cZNC;
        }
        else{
          xyZNC[0] = xyZNC[1] = 0.;
        }
        if(denZNA>0.){
         Double_t nSpecnA = SumEZNA/Enucl;
         cZNA = 1.89358-0.71262/(nSpecnA+0.71789);
         xyZNA[0] = cZNA*numXZNA/denZNA;
         xyZNA[1] = cZNA*numYZNA/denZNA;
          denZNA *= cZNA;
        }
        else{
          xyZNA[0] = xyZNA[1] = 0.;
        }
      } else {
        for(Int_t i=0; i<4; i++) {
          if(towZNC[i+1]>0.) {
            wZNC = TMath::Power(towZNC[i+1], fZDCGainAlpha);
            numXZNC += x[i]*wZNC;
            numYZNC += y[i]*wZNC;
            denZNC += wZNC;
          }
          if(towZNA[i+1]>0.) {
            wZNA = TMath::Power(towZNA[i+1], fZDCGainAlpha);
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
      
      if(!fAllChONZNC) denZNC=-1.;
      if(!fAllChONZNA) denZNA=-1.;
      
      if(denZNC>0. && pow(xyZNC[0]*xyZNC[0]+xyZNC[1]*xyZNC[1],0.5)>1.E-6) fhZNCenDis[0]->Fill(centrperc,xyZNC[0],xyZNC[1]);
      if(denZNA>0. && pow(xyZNA[0]*xyZNA[0]+xyZNA[1]*xyZNA[1],0.5)>1.E-6) fhZNCenDis[1]->Fill(centrperc,-xyZNA[0], xyZNA[1]);
      
      fFlowEvent->SetZDC2Qsub(xyZNC,denZNC,xyZNA,denZNA);
      
      // ******************************************************************************
      
      Double_t tdcSum = aodZDC->GetZDCTimeSum();
      Double_t tdcDiff = aodZDC->GetZDCTimeDiff();
      fhDebunch->Fill(tdcDiff, tdcSum);
      
      for(int i=0; i<5; i++){
        fhZNCPM[i]->Fill(towZNC[i]);
        if((i<4) && (towZNC[0]>0.)) fhZNCPMQiPMC[i]->Fill(towZNC[i+1]/towZNC[0]);
      }
      for(int i=0; i<5; i++){
        fhZNAPM[i]->Fill(towZNA[i]);
        if(((i<4) && towZNA[0]>0.)) fhZNAPMQiPMC[i]->Fill(towZNA[i+1]/towZNA[0]);
      }
      
      fhZNCvsZNA->Fill(energyZNA, energyZNC);
      fhZDCCvsZDCCA->Fill(energyZNA+energyZPA, energyZNC+energyZPC);
      fhZNCvsZPC->Fill(energyZPC, energyZNC);
      fhZNAvsZPA->Fill(energyZPA, energyZNA);
      fhZNvsZP->Fill(energyZPA+energyZPC, energyZNA+energyZNC);
      fhZNvsVZERO->Fill(multV0A+multV0C, energyZNC+energyZNA);
      fhZDCvsVZERO->Fill(multV0A+multV0C, energyZNA+energyZPA+energyZNC+energyZPC);
      fhZDCvsTracklets->Fill((Double_t) (nTracklets), energyZNA+energyZPA+energyZNC+energyZPC);
      
      Double_t asymmetry = -999.;
      if((energyZNC+energyZNA)>0.) asymmetry = (energyZNC-energyZNA)/(energyZNC+energyZNA);
      fhAsymm->Fill(asymmetry);
      fhZNAvsAsymm->Fill(asymmetry, energyZNA/1000.);
      fhZNCvsAsymm->Fill(asymmetry, energyZNC/1000.);
      
      fhZNCvscentrality->Fill(centrperc, energyZNC/1000.);
      fhZNAvscentrality->Fill(centrperc, energyZNA/1000.);
      
    } // PHYSICS SELECTION
    
  }
  
  // p) cache run number
  fCachedRunNum = fFlowEvent->GetRun();
  
  //  printf("debug: NoRPs %e, NoPOIs %e, RunNum %d, Cen %e \n",fFlowEvent->GetNumberOfRPs(),fFlowEvent->GetNumberOfPOIs(),fCachedRunNum,fFlowEvent->GetCentrality());
  
  PostData(1, fFlowEvent);
  
  PostData(2, fOutput);
}
//________________________________________________________________________

Double_t AliAnalysisTaskCRCZDC::GetBadTowerResp(Double_t Et, TH2D* BadTowerCalibHist)
{
  Double_t EtC = BadTowerCalibHist->ProjectionY("",BadTowerCalibHist->GetXaxis()->FindBin(Et),BadTowerCalibHist->GetXaxis()->FindBin(Et))->GetRandom();
  return EtC;
}

//________________________________________________________________________

Int_t AliAnalysisTaskCRCZDC::GetCenBin(Double_t Centrality)
{
  Int_t CenBin=-1;
  if (Centrality>0. && Centrality<5.) CenBin=0;
  if (Centrality>5. && Centrality<10.) CenBin=1;
  if (Centrality>10. && Centrality<20.) CenBin=2;
  if (Centrality>20. && Centrality<30.) CenBin=3;
  if (Centrality>30. && Centrality<40.) CenBin=4;
  if (Centrality>40. && Centrality<50.) CenBin=5;
  if (Centrality>50. && Centrality<60.) CenBin=6;
  if (Centrality>60. && Centrality<70.) CenBin=7;
  if (Centrality>70. && Centrality<80.) CenBin=8;
  if (Centrality>80. && Centrality<90.) CenBin=9;
  if (CenBin>=fnCen) CenBin=-1;
  if (fnCen==1) CenBin=0;
  return CenBin;
} // end of AliFlowAnalysisCRC::GetCRCCenBin(Double_t Centrality)
//_____________________________________________________________________________

Double_t AliAnalysisTaskCRCZDC::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    printf("One of vertices is not valid\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
  +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}
//________________________________________________________________________

Bool_t AliAnalysisTaskCRCZDC::plpMV(const AliAODEvent* aod)
{
  // check for multi-vertexer pile-up
  
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2 = 5.0;
  const double kMinWDist = 15;
  
  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;
  int nPlp = 0;
  
  if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
  vtPrm = aod->GetPrimaryVertex();
  if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
  
  //int bcPrim = vtPrm->GetBC();
  
  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)aod->GetPileupVertexTracks(ipl);
    //
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
    //  int bcPlp = vtPlp->GetBC();
    //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
    //
    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist) continue;
    //
    return kTRUE; // pile-up: well separated vertices
  }
  
  return kFALSE;
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
