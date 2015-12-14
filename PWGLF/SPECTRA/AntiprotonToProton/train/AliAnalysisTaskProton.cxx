/**************************************************************************
 * Author: Michal Meres.                                               *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//                 AliAnalysisTaskProton class
//            This task is for raw antiproton to proton ratio from ESD only
//-----------------------------------------------------------------
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliAODEvent.h"
#include <AliCFContainer.h>

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TMCProcess.h"

#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"

#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliTPCPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliProdInfo.h"

#include <iostream>
#include "AliLog.h"

using namespace std;

#include "AliAnalysisTaskProton.h"

ClassImp(AliAnalysisTaskProton)

AliAnalysisTaskProton::AliAnalysisTaskProton()
  : AliAnalysisTaskSE(), fPIDMode(kSigma), EventNo(0) ,fProtonContainer(0),fAntiProtonContainer(0),fListAnalysis(0),gHistProtonsDCAxyEtaPt(0),gHistAntiProtonsDCAxyEtaPt(0),fMaxDCAXYFlag(kFALSE), fPtDependentDcaXYFlag(kFALSE), fMaxDCAZFlag(kFALSE),fDebugMode(kFALSE),nbinsPt(6),fLowPt(0.45),fHighPt(1.05), fOADBPath(0), fPIDResponse(0), fRun(0), fOldRun(0), fRecoPass(0), fAnalysisType(0), fCollidingSystems(0), fUsePhysicsSelection(0), fMaxPrimaryVtxPosZ(0), fHistEventStats(0), fGlobalQAList(0), fListQA(0), fQA2DList(0), fHistMultiplicity(0), gHistdEdxP(0), gHistProtonsdEdxP(0), gHistProtonsDCAzEtaPt(0), gHistAntiProtonsDCAzEtaPt(0),gHistProtonsDCAzCentPt (0), gHistAntiProtonsDCAzCentPt(0), gHistFieldProtonsEtaPt(0), gHistFieldAntiProtonsEtaPt(0), gHistFieldProtonsCentPt(0), gHistFieldAntiProtonsCentPt(0), gHistFieldProtonsLengthPt(0), gHistFieldAntiProtonsLengthPt(0), gHistProtonsLengthCentPt(0), gHistAntiProtonsLengthCentPt(0), fPhysicsSelection(0), fMultiplicityMode(0), nbinsY(100), fLowY(0), fHighY(100), fMinTPCClusters(80), fMinITSClusters(2), fMaxChi2PerTPCCluster(3.5), fMaxChi2PerITSCluster(36.), fMaxDCAXY(0.2), fMaxDCAZ(1.), fMaxDCAXYTPCFlag(0), fMaxDCAXYTPC(0),fPtDependentDcaXY(NULL), fNSigmaDCAXY(7), fNBoundP(0.7), fNSigma1(3), fNSigma2(3), fNRatio1(0), fNRatio2(0), fIsMC(kFALSE), fUserDataRecoPass(0), MAXCent(100), MINCent(0) {

}

//________________________________________________________________________
AliAnalysisTaskProton::AliAnalysisTaskProton(const char *name) 
  : AliAnalysisTaskSE(name), fPIDMode(kSigma), EventNo(0) ,fProtonContainer(0),fAntiProtonContainer(0),fListAnalysis(0),gHistProtonsDCAxyEtaPt(0),gHistAntiProtonsDCAxyEtaPt(0),fMaxDCAXYFlag(kFALSE), fPtDependentDcaXYFlag(kFALSE), fMaxDCAZFlag(kFALSE),fDebugMode(kFALSE),nbinsPt(6),fLowPt(0.45),fHighPt(1.05), fOADBPath(0), fPIDResponse(0), fRun(0), fOldRun(0), fRecoPass(0), fAnalysisType(0), fCollidingSystems(0), fUsePhysicsSelection(0), fMaxPrimaryVtxPosZ(0), fHistEventStats(0), fGlobalQAList(0), fListQA(0), fQA2DList(0), fHistMultiplicity(0), gHistdEdxP(0), gHistProtonsdEdxP(0), gHistProtonsDCAzEtaPt(0), gHistAntiProtonsDCAzEtaPt(0),gHistProtonsDCAzCentPt (0), gHistAntiProtonsDCAzCentPt(0), gHistFieldProtonsEtaPt(0), gHistFieldAntiProtonsEtaPt(0), gHistFieldProtonsCentPt(0), gHistFieldAntiProtonsCentPt(0), gHistFieldProtonsLengthPt(0), gHistFieldAntiProtonsLengthPt(0), gHistProtonsLengthCentPt(0), gHistAntiProtonsLengthCentPt(0), fPhysicsSelection(0), fMultiplicityMode(0), nbinsY(100), fLowY(0), fHighY(100), fMinTPCClusters(80), fMinITSClusters(2), fMaxChi2PerTPCCluster(3.5), fMaxChi2PerITSCluster(36.), fMaxDCAXY(0.2), fMaxDCAZ(1.), fMaxDCAXYTPCFlag(0), fMaxDCAXYTPC(0),fPtDependentDcaXY(NULL), fNSigmaDCAXY(7), fNBoundP(0.7), fNSigma1(3), fNSigma2(3), fNRatio1(0), fNRatio2(0), fIsMC(kFALSE), fUserDataRecoPass(0), MAXCent(100), MINCent(0) {
 // Constructor
//fPtDependentDcaXY=NULL;

  // Define output slots only here
   fPhysicsSelection = new AliPhysicsSelection();
  AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if(lmcEvtHandler) fPhysicsSelection->SetAnalyzeMC();

  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (handler) {
    handler->SetEventSelection(fPhysicsSelection);
    AliInfo("Physics Event Selection enabled.");
  } else {
    AliError("No input event handler connected to analysis manager. No Physics Event Selection.");
  }
  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
//  DefineOutput(3, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskProton::~AliAnalysisTaskProton(){
  if (fListAnalysis) { delete fListAnalysis; fListAnalysis = 0x0; }
  if (fPhysicsSelection) {delete fPhysicsSelection;}
  if (fPtDependentDcaXY) {delete fPtDependentDcaXY;}

}
//________________________________________________________________________
void AliAnalysisTaskProton::UserCreateOutputObjects(){
  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");

  //pid response object
  AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if(lmcEvtHandler) inputHandler->CreatePIDResponse(kTRUE);
  else 
  inputHandler->CreatePIDResponse(kFALSE);
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("PIDResponse object was not created");

  fPIDResponse->SetOADBPath(AliAnalysisManager::GetOADBPath());
  if (!fOADBPath.IsNull()) fPIDResponse->SetOADBPath(fOADBPath.Data());

 // fTPCpid = dynamic_cast<AliTPCPIDResponse&>(fPIDResponse->GetTPCResponse());

    if (fMultiplicityMode){
	nbinsY = 100; fLowY = 0; fHighY = 100; //nbinsPt = 6; fLowPt = 0.45; fHighPt = 1.05;
}else
{	nbinsY = 5; fLowY = -0.5; fHighY = 0.5; //nbinsPt = 6; fLowPt = 0.45; fHighPt = 1.05;
}

 const Int_t nBinsdca = 1000;
  Double_t dcamin = -10., dcamax = 10.;

  fListAnalysis = new TList();
  fListAnalysis->SetOwner();
  fListAnalysis->SetName("ListAnalysis");

  TString gCutName[4] = {"Total","Physics Selection",
		       "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
			     "Event statistics;;N_{events}",
			     4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++) 
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());

  fHistMultiplicity = new TH1D("fHistCentrality","Centrality",200,0,200); //110
  fHistMultiplicity->GetXaxis()->SetTitle("Centrality");
  fHistMultiplicity->GetYaxis()->SetTitle("Counts");

  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;
  Double_t *binLimY = new Double_t[nbinsY+1];
  Double_t *binLimPt = new Double_t[nbinsPt+1];
  //values for bin lower bounds
  for(Int_t i = 0; i <= nbinsY; i++) 
    binLimY[i]=(Double_t)fLowY  + (fHighY - fLowY)  /nbinsY*(Double_t)i;
  for(Int_t i = 0; i <= nbinsPt; i++) 
    binLimPt[i]=(Double_t)fLowPt  + (fHighPt - fLowPt)  /nbinsPt*(Double_t)i;

  fProtonContainer = new AliCFContainer("containerProtons",
					"container for protons",
					4,2,iBin);
  fProtonContainer->SetBinLimits(0,binLimY); //rapidity or eta
  fProtonContainer->SetBinLimits(1,binLimPt); //pT
  fProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  fProtonContainer->SetVarTitle(0,"y");
  fAntiProtonContainer = new AliCFContainer("containerAntiProtons",
					    "container for antiprotons",
					    4,2,iBin);
  fAntiProtonContainer->SetBinLimits(0,binLimY); //rapidity or eta
  fAntiProtonContainer->SetBinLimits(1,binLimPt); //pT
  fAntiProtonContainer->SetVarTitle(1,"P_{T} [GeV/c]");
  fAntiProtonContainer->SetVarTitle(0,"y");

fListAnalysis->Add(fProtonContainer);
fListAnalysis->Add(fAntiProtonContainer);
fListAnalysis->Add(fHistEventStats);
fListAnalysis->Add(fHistMultiplicity);

  //========================QA================================//

  fListQA = new TList();
  fListQA->SetOwner();
  fListQA->SetName("fListQA");
  
  fGlobalQAList = new TList();
  fGlobalQAList->SetName("fGlobalQAList");
  fListQA->Add(fGlobalQAList);

  fQA2DList = new TList();
  fQA2DList->SetName("fQA2DList");
  fGlobalQAList->Add(fQA2DList);

  //dEdx plots
  gHistdEdxP = new TH2F("gHistdEdxP","dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
  fQA2DList->Add(gHistdEdxP);
  gHistProtonsdEdxP = new TH2F("gHistProtonsdEdxP","Accepted protons dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
  fQA2DList->Add(gHistProtonsdEdxP);
 //dca vs pT for protons & antiprotons
  gHistProtonsDCAxyEtaPt = new TH3F("gHistProtonsDCAxyEtaPt",
					  ";P_{T} [GeV/c];dca_{xy} [cm]",
					  nbinsY,fLowY,fHighY,
					  nbinsPt,fLowPt,fHighPt,
					  nBinsdca,dcamin, dcamax);

  gHistProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsDCAxyEtaPt);
   gHistAntiProtonsDCAxyEtaPt = new TH3F("gHistAntiProtonsDCAxyEtaPt",
					  ";P_{T} [GeV/c];dca_{xy} [cm]",
					  nbinsY,fLowY,fHighY,
					  nbinsPt,fLowPt,fHighPt,
					  nBinsdca,dcamin, dcamax);

  gHistAntiProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistAntiProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistAntiProtonsDCAxyEtaPt);

/*  gHistFAKEProtonsDCAxyEtaPt = new TH3F("gHistFAKEProtonsDCAxyEtaPt",
					  ";P_{T} [GeV/c];dca_{xy} [cm]",
					  nbinsY,fLowY,fHighY,
					  nbinsPt,fLowPt,fHighPt,
					  nBinsdca,dcamin, dcamax);

  gHistFAKEProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKEProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKEProtonsDCAxyEtaPt);
   gHistFAKEAntiProtonsDCAxyEtaPt = new TH3F("gHistFAKEAntiProtonsDCAxyEtaPt",
					  ";P_{T} [GeV/c];dca_{xy} [cm]",
					  nbinsY,fLowY,fHighY,
					  nbinsPt,fLowPt,fHighPt,
					  nBinsdca,dcamin, dcamax);

  gHistFAKEAntiProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKEAntiProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKEAntiProtonsDCAxyEtaPt);

  gHistFAKESPDProtonsDCAxyEtaPt = new TH3F("gHistFAKESPDProtonsDCAxyEtaPt",
					  ";P_{T} [GeV/c];dca_{xy} [cm]",
					  nbinsY,fLowY,fHighY,
					  nbinsPt,fLowPt,fHighPt,
					  nBinsdca,dcamin, dcamax);

  gHistFAKESPDProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKESPDProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKESPDProtonsDCAxyEtaPt);
   gHistFAKESPDAntiProtonsDCAxyEtaPt = new TH3F("gHistFAKESPDAntiProtonsDCAxyEtaPt",
					  ";P_{T} [GeV/c];dca_{xy} [cm]",
					  nbinsY,fLowY,fHighY,
					  nbinsPt,fLowPt,fHighPt,
					  nBinsdca,dcamin, dcamax);

  gHistFAKESPDAntiProtonsDCAxyEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKESPDAntiProtonsDCAxyEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKESPDAntiProtonsDCAxyEtaPt);*/

//========================DCAz================================//

  gHistProtonsDCAzEtaPt = new TH3F("gHistProtonsDCAzEtaPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsDCAzEtaPt);

  gHistAntiProtonsDCAzEtaPt = new TH3F("gHistAntiProtonsDCAzEtaPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistAntiProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistAntiProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistAntiProtonsDCAzEtaPt);

  gHistProtonsDCAzCentPt = new TH3F("gHistProtonsDCAzCentPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistProtonsDCAzCentPt->GetXaxis()->SetTitle("Centrality");
  gHistProtonsDCAzCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistProtonsDCAzCentPt);

  gHistAntiProtonsDCAzCentPt = new TH3F("gHistAntiProtonsDCAzCentPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistAntiProtonsDCAzCentPt->GetXaxis()->SetTitle("Centrality");
  gHistAntiProtonsDCAzCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistAntiProtonsDCAzCentPt);

//========================DCAz FAKE================================//
/*
  gHistFAKEProtonsDCAzEtaPt = new TH3F("gHistFAKEProtonsDCAzEtaPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKEProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKEProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKEProtonsDCAzEtaPt);

  gHistFAKEAntiProtonsDCAzEtaPt = new TH3F("gHistFAKEAntiProtonsDCAzEtaPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKEAntiProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKEAntiProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKEAntiProtonsDCAzEtaPt);

  gHistFAKEProtonsDCAzCentPt = new TH3F("gHistFAKEProtonsDCAzCentPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKEProtonsDCAzCentPt->GetXaxis()->SetTitle("Centrality");
  gHistFAKEProtonsDCAzCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKEProtonsDCAzCentPt);

  gHistFAKEAntiProtonsDCAzCentPt = new TH3F("gHistFAKEAntiProtonsDCAzCentPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKEAntiProtonsDCAzCentPt->GetXaxis()->SetTitle("Centrality");
  gHistFAKEAntiProtonsDCAzCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKEAntiProtonsDCAzCentPt);

//========================DCAz FAKE SPD================================//

  gHistFAKESPDProtonsDCAzEtaPt = new TH3F("gHistFAKESPDProtonsDCAzEtaPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKESPDProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKESPDProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKESPDProtonsDCAzEtaPt);

  gHistFAKESPDAntiProtonsDCAzEtaPt = new TH3F("gHistFAKESPDAntiProtonsDCAzEtaPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKESPDAntiProtonsDCAzEtaPt->GetXaxis()->SetTitle("y");
  gHistFAKESPDAntiProtonsDCAzEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKESPDAntiProtonsDCAzEtaPt);

  gHistFAKESPDProtonsDCAzCentPt = new TH3F("gHistFAKESPDProtonsDCAzCentPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKESPDProtonsDCAzCentPt->GetXaxis()->SetTitle("Centrality");
  gHistFAKESPDProtonsDCAzCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKESPDProtonsDCAzCentPt);

  gHistFAKESPDAntiProtonsDCAzCentPt = new TH3F("gHistFAKESPDAntiProtonsDCAzCentPt",
					  ";P_{T} [GeV/c];dca_{z} [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  500,0, 10);
  gHistFAKESPDAntiProtonsDCAzCentPt->GetXaxis()->SetTitle("Centrality");
  gHistFAKESPDAntiProtonsDCAzCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFAKESPDAntiProtonsDCAzCentPt);*/

//========================Magnetic field y Pt================================//
  gHistFieldProtonsEtaPt = new TH3F("gHistFieldProtonsEtaPt",
					  ";P_{T} [GeV/c];Mag. field",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  2,-10, 10);
  gHistFieldProtonsEtaPt->GetXaxis()->SetTitle("y");
  gHistFieldProtonsEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFieldProtonsEtaPt);

  gHistFieldAntiProtonsEtaPt = new TH3F("gHistFieldAntiProtonsEtaPt",
					  ";P_{T} [GeV/c];Mag. field",
					  5,-0.5,0.5,
					  nbinsPt,fLowPt,fHighPt,
					  2,-10, 10);
  gHistFieldAntiProtonsEtaPt->GetXaxis()->SetTitle("y");
  gHistFieldAntiProtonsEtaPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFieldAntiProtonsEtaPt);

//========================Magnetic field cent Pt================================//
 
  gHistFieldProtonsCentPt = new TH3F("gHistFieldProtonsCentPt",
					  ";P_{T} [GeV/c];Mag. field",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  2,-10, 10);
  gHistFieldProtonsCentPt->GetXaxis()->SetTitle("Centrality");
  gHistFieldProtonsCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFieldProtonsCentPt);

  gHistFieldAntiProtonsCentPt = new TH3F("gHistFieldAntiProtonsCentPt",
					  ";P_{T} [GeV/c];Mag. field",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  2,-10, 10);
  gHistFieldAntiProtonsCentPt->GetXaxis()->SetTitle("Centrality");
  gHistFieldAntiProtonsCentPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFieldAntiProtonsCentPt);

//========================Magnetic field Length Pt================================//

  gHistFieldProtonsLengthPt = new TH3F("gHistFieldProtonsLengthPt",
					  ";P_{T} [GeV/c];Mag. field",
					  640,0,260,
					  nbinsPt,fLowPt,fHighPt,
					  2,-10, 10);
  gHistFieldProtonsLengthPt->GetXaxis()->SetTitle("track length in TPC (cm)");
  gHistFieldProtonsLengthPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFieldProtonsLengthPt);

  gHistFieldAntiProtonsLengthPt = new TH3F("gHistFieldAntiProtonsLengthPt",
					  ";P_{T} [GeV/c];Mag. field",
					  640,0,260,
					  nbinsPt,fLowPt,fHighPt,
					  2,-10, 10);
  gHistFieldAntiProtonsLengthPt->GetXaxis()->SetTitle("track length in TPC (cm)");
  gHistFieldAntiProtonsLengthPt->SetStats(kTRUE);
  fQA2DList->Add(gHistFieldAntiProtonsLengthPt);

//========================Length Cent Pt================================//
gHistProtonsLengthCentPt = new TH3F("gHistProtonsLengthCentPt",
					  ";P_{T} [GeV/c];Length [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  640,0,260);
gHistProtonsLengthCentPt->GetXaxis()->SetTitle("Centrality");
gHistProtonsLengthCentPt->SetStats(kTRUE);
fQA2DList->Add(gHistProtonsLengthCentPt);

gHistAntiProtonsLengthCentPt = new TH3F("gHistAntiProtonsLengthCentPt",
					  ";P_{T} [GeV/c];Length [cm]",
					  100,0,100,
					  nbinsPt,fLowPt,fHighPt,
					  640,0,260);
gHistAntiProtonsLengthCentPt->GetXaxis()->SetTitle("Centrality");
gHistAntiProtonsLengthCentPt->SetStats(kTRUE);
fQA2DList->Add(gHistAntiProtonsLengthCentPt);

/*
  //==========================Systematics==============================//
  fListSystematics = new TList();
  fListSystematics->SetOwner();
  fListSystematics->SetName("fListSystematics");

//-----------------------LOW--------------------------------------
fSysProtonLow = new TList();
fSysProtonLow->SetName("fListSysProtonLow");

TH2F *gProtonMaxDCAXYLow = new TH2F("gProtonMaxDCAXYLow","Proton-MaxDCAXY-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonLow->Add(gProtonMaxDCAXYLow);
TH2F *gProtonMinTPCClustersLow = new TH2F("gProtonMinTPCClustersLow","Proton-MinTPCClusters-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonLow->Add(gProtonMinTPCClustersLow);
TH2F *gProtonNsigmaLow = new TH2F("gProtonNSigmaLow","Proton-Nsigma-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonLow->Add(gProtonNsigmaLow);
TH2F *gProtonAcceptedVertexDiamondLow = new TH2F("gProtonAcceptedVertexDiamondLow","Proton-Accepted Vertex Diamond-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonLow->Add(gProtonAcceptedVertexDiamondLow);

fSysAntiProtonLow = new TList();
fSysAntiProtonLow->SetName("fListSysAntiProtonLow");

TH2F *gAntiProtonMaxDCAXYLow = new TH2F("gAntiProtonMaxDCAXYLow","AntiProton-MaxDCAXY-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonLow->Add(gAntiProtonMaxDCAXYLow);
TH2F *gAntiProtonMinTPCClustersLow = new TH2F("gAntiProtonMinTPCClustersLow","AntiProton-MinTPCClusters-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonLow->Add(gAntiProtonMinTPCClustersLow);
TH2F *gAntiProtonNsigmaLow = new TH2F("gAntiProtonNSigmaLow","AntiProton-Nsigma-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonLow->Add(gAntiProtonNsigmaLow);
TH2F *gAntiProtonAcceptedVertexDiamondLow = new TH2F("gAntiProtonAcceptedVertexDiamondLow","AntiProton-Accepted Vertex Diamond-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonLow->Add(gAntiProtonAcceptedVertexDiamondLow);
//----------------------END LOW------------------------------------
fSysProtonHigh = new TList();
fSysProtonHigh->SetName("fListSysProtonHigh");

TH2F *gProtonMaxDCAXYHigh = new TH2F("gProtonMaxDCAXYHigh","Proton-MaxDCAXY-High",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonHigh->Add(gProtonMaxDCAXYHigh);
TH2F *gProtonMinTPCClustersHigh = new TH2F("gProtonMinTPCClustersHigh","Proton-MinTPCClusters-High",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonHigh->Add(gProtonMinTPCClustersHigh);
TH2F *gProtonNsigmaHigh = new TH2F("gProtonNSigmaHigh","Proton-Nsigma-High",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonHigh->Add(gProtonNsigmaHigh);
TH2F *gProtonAcceptedVertexDiamondHigh = new TH2F("gProtonAcceptedVertexDiamondHigh","Proton-Accepted Vertex Diamond-High",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysProtonHigh->Add(gProtonAcceptedVertexDiamondHigh);

fSysAntiProtonHigh = new TList();
fSysAntiProtonHigh->SetName("fListSysAntiProtonHigh");

TH2F *gAntiProtonMaxDCAXYHigh = new TH2F("gAntiProtonMaxDCAXYHigh","AntiProton-MaxDCAXY-High",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonHigh->Add(gAntiProtonMaxDCAXYHigh);
TH2F *gAntiProtonMinTPCClustersHigh = new TH2F("gAntiProtonMinTPCClustersHigh","AntiProton-MinTPCClusters-High",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonHigh->Add(gAntiProtonMinTPCClustersHigh);
TH2F *gAntiProtonNsigmaHigh = new TH2F("gAntiProtonNSigmaHigh","AntiProton-Nsigma-High",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonHigh->Add(gAntiProtonNsigmaHigh);
TH2F *gAntiProtonAcceptedVertexDiamondHigh = new TH2F("gAntiProtonAcceptedVertexDiamondHigh","AntiProton-Accepted Vertex Diamond-Low",nbinsY,fLowY,fHighY,nbinsPt,fLowPt,fHighPt);
fSysAntiProtonHigh->Add(gAntiProtonAcceptedVertexDiamondHigh);

fListSystematics->Add(fSysProtonLow);
fListSystematics->Add(fSysAntiProtonLow);
fListSystematics->Add(fSysProtonHigh);
fListSystematics->Add(fSysAntiProtonHigh);*/

  PostData(1, fListAnalysis);
  PostData(2, fListQA);
//  PostData(3, fListSystematics);
}

void AliAnalysisTaskProton::UserExec(Option_t *){

// Called for each event
  AliVEvent* lEvent = InputEvent();
  if (!lEvent) {
    Printf("ERROR: Event not available");
    return;
  }
  AliESDEvent* lESDEvent = (AliESDEvent*)lEvent;
  Double_t B = lESDEvent->GetMagneticField();
//Printf("Magnetic field %lf",B);
  fHistEventStats->Fill(1); //number of analyzed events
  EventNo++;
  //Printf("Event number %i",EventNo);
/*  fRun=lEvent->GetRunNumber();
  if (fRun!=fOldRun){
    SetRecoInfo();
    fOldRun=fRun;
  }


fPIDResponse->InitialiseEvent(lEvent,fRecoPass);*/

Int_t nTracks = 0;

AliCentrality *esdCentrality = lESDEvent->GetCentrality();
Float_t nTracklets = esdCentrality->GetCentralityPercentile("V0M");
  

  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if ((fUsePhysicsSelection)&&(!isSelected)) return;
  fHistEventStats->Fill(2); //number of analyzed events
fHistMultiplicity->Fill(nTracklets);
  if(esdCentrality->GetCentralityPercentile("V0M") > MAXCent || esdCentrality->GetCentralityPercentile("V0M") < MINCent) return;


 fHistEventStats->Fill(3); //number of analyzed events
const AliVVertex *primaryVtx = lEvent->GetPrimaryVertex();
const AliESDVertex *trkVtx = lESDEvent->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0) return;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = trkVtx->GetZ();
const AliESDVertex* spdVtx = lESDEvent->GetPrimaryVertexSPD();
  if (spdVtx->GetNContributors()<=0) return;
  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
  if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
 





  fHistEventStats->Fill(4); //number of analyzed events

  Double_t containerInput[2] ;
  Double_t gPt = 0.0, gP = 0.0;
  Float_t dcaXY = 0.0, dcaZ = 0.0;

 nTracks = lESDEvent->GetNumberOfTracks();
   for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    AliESDtrack* track = lESDEvent->GetTrack(iTracks);
    AliESDtrack trackTPC;

    //Int_t nClustersTPC = track->GetTPCclusters(0x0);
    //Int_t npointsTPCdEdx = track->GetTPCsignalN();
    Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
    Double_t dca3D = 0.0;
    
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) continue;
if (TMath::Abs(zvtx)<fMaxPrimaryVtxPosZ){// event selections
    gPt = tpcTrack->Pt();
    gP = track->GetInnerParam()->P();
    Double_t Length = track->GetLengthInActiveZone(0,3,236,B,0,0);

    AliExternalTrackParam cParam;
    track->RelateToVertex(spdVtx,
			    lESDEvent->GetMagneticField(),
			    100.,&cParam);
    track->GetImpactParameters(dcaXY,dcaZ);
    dca[0] = dcaXY; dca[1] = dcaZ;

	dca3D = TMath::Sqrt(TMath::Power(dca[0],2) + TMath::Power(dca[1],2));

	 if(fMultiplicityMode)
        containerInput[0] = nTracklets;
      else
	containerInput[0] = Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3);
      containerInput[1] = gPt;

      if(IsAccepted(track)) {
	gHistProtonsdEdxP->Fill(gP,track->GetTPCsignal());
		if(IsProton(track)) {
		gHistdEdxP->Fill(gP,track->GetTPCsignal());
			if(tpcTrack->Charge() > 0){

			gHistProtonsDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),
						  tpcTrack->Pt(),
						  dca[1]);

			gHistProtonsDCAzCentPt->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[1]);

			if(fMultiplicityMode){
gHistProtonsDCAxyEtaPt->Fill(nTracklets,tpcTrack->Pt(),dca[0]);
	    }
	    else {
	      gHistProtonsDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),tpcTrack->Pt(),dca[0]);
	    }
	  					}//protons
	else if(tpcTrack->Charge() < 0){




	    gHistAntiProtonsDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),
						  tpcTrack->Pt(),
						  dca[1]);
	    gHistAntiProtonsDCAzCentPt->Fill(nTracklets,
						  tpcTrack->Pt(),
						  dca[1]);

			if(fMultiplicityMode){
	      gHistAntiProtonsDCAxyEtaPt->Fill(nTracklets,tpcTrack->Pt(),dca[0]);
	    }
	    else {

	      gHistAntiProtonsDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),
						  tpcTrack->Pt(),
						  dca[0]);
	    }
	  					}//antiprotons
						if(IsPrimary(lESDEvent,spdVtx,track)) {
	if(((gPt > fLowPt) && (gPt < fHighPt)) && ((tpcTrack->P() > fLowPt) && (tpcTrack->P() < fHighPt)) && ((Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3) > -0.5) && (Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3) < 0.5))) {

		if(tpcTrack->Charge() > 0) {
		gHistFieldProtonsEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),tpcTrack->Pt(),B);
		gHistFieldProtonsCentPt->Fill(nTracklets,tpcTrack->Pt(),B);
		gHistFieldProtonsLengthPt->Fill(Length,tpcTrack->Pt(),B);
		gHistProtonsLengthCentPt->Fill(nTracklets,tpcTrack->Pt(),Length);

		fProtonContainer->Fill(containerInput,3);   
	      }//protons
	      else if(tpcTrack->Charge() < 0) {
		gHistFieldAntiProtonsEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),tpcTrack->Pt(),B);
		gHistFieldAntiProtonsCentPt->Fill(nTracklets,tpcTrack->Pt(),B);
		gHistFieldAntiProtonsLengthPt->Fill(Length,tpcTrack->Pt(),B);
		gHistAntiProtonsLengthCentPt->Fill(nTracklets,tpcTrack->Pt(),Length);
		
		fAntiProtonContainer->Fill(containerInput,3);
	      }//antiprotons

        }//end 	IsInPhaseSpace							
	}//end 	IsPrimary						
	}//end 	IsProton
	}//end 	IsAccepted	
	}//end event selection
//FillSystematics(lESDEvent,spdVtx,track);
		}//Track loop


  PostData(1, fListAnalysis);
  PostData(2, fListQA);
//  PostData(3, fListSystematics);
}

void AliAnalysisTaskProton::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  gStyle->SetPalette(1,0);

  fListAnalysis = dynamic_cast<TList*> (GetOutputData(1));
  if (!fListAnalysis) {
    Printf("ERROR: fListAnalysis not available");
    return;
  }
   
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskProton::IsLabelUsed(TArrayI labelArray, 
					Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskProton::IsProton(AliESDtrack *track) {
  //Function that checks if a track is a proton
  
   Double_t gPt = 0.0, gP = 0.0, gEta = 0.0;
  
  //Ratio of the measured over the theoretical dE/dx a la STAR
  if(fPIDMode == kRatio) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(tpcTrack) {
      gPt = tpcTrack->Pt();
      gP = track->GetInnerParam()->P();
      gEta = tpcTrack->Eta();
    }
    
    Double_t normalizeddEdx = -10.;
    if((track->GetTPCsignal() > 0.0) && (fPIDResponse->GetTPCResponse().GetExpectedSignal(gP,AliPID::kProton) > 0.0))
      normalizeddEdx = TMath::Log(track->GetTPCsignal()/fPIDResponse->GetTPCResponse().GetExpectedSignal(gP,AliPID::kProton));

    if (gP <= fNBoundP) if(normalizeddEdx >= fNRatio1) return kTRUE;
    if (gP >  fNBoundP) if(normalizeddEdx >= fNRatio2) return kTRUE;
  }//kRatio PID mode

  //Definition of an N-sigma area around the dE/dx vs P band
 else if(fPIDMode == kSigma) {
   
    Double_t nsigma = 100.0;
    
    Double_t mom = track->GetP();
    const AliExternalTrackParam *in = track->GetInnerParam();
    if (in)
      mom = in->GetP();

    nsigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
  
      if(nsigma <= fNSigma1) return kTRUE;
    //if (mom <= fNBoundP) if(nsigma <= fNSigma1) return kTRUE;
    //if (mom >  fNBoundP) if(nsigma <= fNSigma2) return kTRUE;
 }//kSigma PID method 

  return kFALSE;
}

Bool_t AliAnalysisTaskProton::IsAccepted(AliESDtrack* track) {
  // Checks if the track is excluded from the cuts
  //Int_t  fIdxInt[200];
  //Int_t nClustersITS = track->GetITSclusters(fIdxInt);
  //Int_t nClustersTPC = track->GetTPCclusters(fIdxInt);
  Int_t nClustersITS = track->GetITSclusters(0x0);
  Int_t nClustersTPC = track->GetTPCclusters(0x0);

  Float_t chi2PerClusterITS = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

 // Double_t extCov[15];
 // track->GetExternalCovariance(extCov);

if(nClustersITS < fMinITSClusters) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has %d ITS points (min. requested: %d)",nClustersITS,fMinITSClusters);
      return kFALSE;
    }
if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has a chi2 per ITS cluster %lf (max. requested: %lf)",chi2PerClusterITS,fMaxChi2PerITSCluster);
      return kFALSE;
    }
if(nClustersTPC < fMinTPCClusters) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has %d TPC clusters (min. requested: %d)",nClustersTPC,fMinTPCClusters);
      return kFALSE;
    }
if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has a chi2 per TPC cluster %lf (max. requested: %lf)",chi2PerClusterTPC,fMaxChi2PerTPCCluster);
      return kFALSE; 
    }
if(track->GetTPCsignalN() < fMinTPCClusters) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has %d TPC points for the calculation of the energy loss (min. requested: %d)",track->GetTPCsignalN(),fMinTPCClusters);
      return kFALSE;
    }
return kTRUE;
}

Bool_t AliAnalysisTaskProton::IsPrimary(AliESDEvent *esd,
					const AliESDVertex *vertex, 
					AliESDtrack* track) {
  // Checks if the track is a primary-like candidate
  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  Float_t dcaXY = 0.0, dcaZ = 0.0;

    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    AliExternalTrackParam cParam;
    if(!tpcTrack) {
      gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
      dca[0] = -100.; dca[1] = -100.; dca3D = -100.;
      cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
    }
    else {
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();

      track->RelateToVertex(vertex,
			    esd->GetMagneticField(),
			    100.,&cParam);
      track->GetImpactParameters(dcaXY,dcaZ);
      dca[0] = dcaXY; dca[1] = dcaZ;
    }
  //standalone TPC or hybrid TPC approaches
  dca3D = TMath::Sqrt(TMath::Power(dca[0],2) +
		      TMath::Power(dca[1],2));

   if(fMaxDCAXYFlag) { 
    if(TMath::Abs(dca[0]) > fMaxDCAXY) {
	if(fDebugMode)
	Printf("IsPrimary: Track rejected because it has a value of dca(xy) of %lf (max. requested: %lf)",TMath::Abs(dca[0]),fMaxDCAXY);
      return kFALSE;
    }
  }

 if(fPtDependentDcaXYFlag) {
    if(TMath::Abs(dca[0]) > fNSigmaDCAXY*fPtDependentDcaXY->Eval(gPt)) { // kMicrometer2Centimeter*
      if(fDebugMode)
	Printf("IsPrimary: Track rejected because it has a value of the dca(xy) higher than the %d sigma pt dependent cut: %lf (max. requested: %lf) and Pt is %lf",fNSigmaDCAXY,TMath::Abs(dca[0]),fNSigmaDCAXY*fPtDependentDcaXY->Eval(gPt),gPt);
      return kFALSE;
    }
  }

   if(fMaxDCAZFlag) { 
    if(TMath::Abs(dca[1]) > fMaxDCAZ) {
	if(fDebugMode)
	Printf("IsPrimary: Track rejected because it has a value of dca(z) of %lf (max. requested: %lf)",TMath::Abs(dca[1]),fMaxDCAZ);
      return kFALSE;
    }
  }

return kTRUE;
  }

//____________________________________________________________________//
Double_t AliAnalysisTaskProton::Rapidity(Double_t gPx, 
					      Double_t gPy, 
					      Double_t gPz, Int_t fType) const {
  //returns the rapidity of the proton - to be removed
  Double_t fMass = 9.38270000000000048e-01;
  if (fType == 0) fMass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  if (fType == 1) fMass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  
  Double_t gP = TMath::Sqrt(TMath::Power(gPx,2) + 
                           TMath::Power(gPy,2) + 
			   TMath::Power(gPz,2));
  Double_t energy = TMath::Sqrt(gP*gP + fMass*fMass);
  Double_t y = -999;
  if(energy != gPz) 
    y = 0.5*TMath::Log((energy + gPz)/(energy - gPz));

  return y;
}

void AliAnalysisTaskProton::SetPtDependentDCAxy(Int_t nSigma, 
						Double_t p0, 
						Double_t p1, 
						Double_t p2) {
  //Pt dependent dca cut in xy
  fNSigmaDCAXY = nSigma;
  if (fPtDependentDcaXY != NULL) delete fPtDependentDcaXY;
  fPtDependentDcaXY = new TF1("fPtDependentDcaXY","[0]+[1]/x^[2]",0.1,10.1);
  fPtDependentDcaXY->SetParameter(0,p0); 
  fPtDependentDcaXY->SetParameter(1,p1);
  fPtDependentDcaXY->SetParameter(2,p2);
  fPtDependentDcaXYFlag = kTRUE; 
}


