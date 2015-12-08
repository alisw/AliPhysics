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
//                 AliAnalysisTaskProtonQA class
//            This task is for QAing the protons and antiprotons from ESD only
//-----------------------------------------------------------------
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"

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

#include "AliAnalysisTaskProtonQA.h"

ClassImp(AliAnalysisTaskProtonQA)

AliAnalysisTaskProtonQA::AliAnalysisTaskProtonQA()
  : AliAnalysisTaskSE(), fPIDResponse(0), fAnalysisType("ESD"), fCollidingSystems(0), fNBinsPt(6) ,fMinPt(0.45), fMaxPt(1.05),
    fUsePhysicsSelection(1), fMaxPrimaryVtxPosZ(100.) ,fPhysicsSelection(0),fMultiplicityMode(0), fListHist(),
    gHistPrimaryProtonsDCAxyEtaPt(0), gHistPrimaryAntiProtonsDCAxyEtaPt(0), gHistSecondaryProtonsFromWeakDCAxyEtaPt(0), gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt(0),
    gHistSecondaryProtonsFromHadronicDCAxyEtaPt(0), gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt(0), gHistdEdxP(0), gHistProtonsdEdxP(0),
    fPIDMode(kSigma), fDebugMode(kFALSE), fNBinsY(100), fMinY(0), fMaxY(100), fMinTPCClusters(80), fMinITSClusters(2), fMaxChi2PerTPCCluster(3.5), fMaxChi2PerITSCluster(36.), fMaxDCAXY(10.), fMaxDCAZ(1.), fNBoundP(0.7), fNSigma1(3), fNSigma2(3), fNRatio1(0), fNRatio2(0), MAXCent(100), MINCent(0), fMaxDCAXYFlag(kTRUE), fMaxDCAZFlag(kFALSE), fListQA(0), gHistPrimaryProtonsDCAzEtaPt(0), gHistPrimaryAntiProtonsDCAzEtaPt(0), gHistSecondaryProtonsFromWeakDCAzEtaPt(0), gHistSecondaryAntiProtonsFromWeakDCAzEtaPt(0), gHistSecondaryProtonsFromHadronicDCAzEtaPt(0), gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt(0), gHistPrimaryProtonsDCAzCentPt(0), gHistPrimaryAntiProtonsDCAzCentPt(0), gHistSecondaryProtonsFromWeakDCAzCentPt(0), gHistSecondaryAntiProtonsFromWeakDCAzCentPt(0), gHistSecondaryProtonsFromHadronicDCAzCentPt(0), gHistSecondaryAntiProtonsFromHadronicDCAzCentPt(0), fOADBPath(0)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskProtonQA::AliAnalysisTaskProtonQA(const char *name) 
  : AliAnalysisTaskSE(name), fPIDResponse(0), fAnalysisType("ESD"), fCollidingSystems(0), fNBinsPt(6) ,fMinPt(0.45), fMaxPt(1.05),
    fUsePhysicsSelection(1), fMaxPrimaryVtxPosZ(100.) ,fPhysicsSelection(0),fMultiplicityMode(0), fListHist(),
    gHistPrimaryProtonsDCAxyEtaPt(0), gHistPrimaryAntiProtonsDCAxyEtaPt(0), gHistSecondaryProtonsFromWeakDCAxyEtaPt(0), gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt(0), gHistdEdxP(0), gHistProtonsdEdxP(0),
    gHistSecondaryProtonsFromHadronicDCAxyEtaPt(0), gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt(0),
    fPIDMode(kSigma), fDebugMode(kFALSE), fNBinsY(100), fMinY(0), fMaxY(100), fMinTPCClusters(80), fMinITSClusters(2), fMaxChi2PerTPCCluster(3.5), fMaxChi2PerITSCluster(36.), fMaxDCAXY(10.), fMaxDCAZ(1.), fNBoundP(0.7), fNSigma1(3), fNSigma2(3), fNRatio1(0), fNRatio2(0), MAXCent(100), MINCent(0), fMaxDCAXYFlag(kTRUE), fMaxDCAZFlag(kFALSE), fListQA(0), gHistPrimaryProtonsDCAzEtaPt(0), gHistPrimaryAntiProtonsDCAzEtaPt(0), gHistSecondaryProtonsFromWeakDCAzEtaPt(0), gHistSecondaryAntiProtonsFromWeakDCAzEtaPt(0), gHistSecondaryProtonsFromHadronicDCAzEtaPt(0), gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt(0), gHistPrimaryProtonsDCAzCentPt(0), gHistPrimaryAntiProtonsDCAzCentPt(0), gHistSecondaryProtonsFromWeakDCAzCentPt(0), gHistSecondaryAntiProtonsFromWeakDCAzCentPt(0), gHistSecondaryProtonsFromHadronicDCAzCentPt(0), gHistSecondaryAntiProtonsFromHadronicDCAzCentPt(0), fOADBPath(0)
{
  // Constructor
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
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskProtonQA::~AliAnalysisTaskProtonQA(){
  if (fListHist) { delete fListHist; fListHist = 0x0; }
  if (fPhysicsSelection) {delete fPhysicsSelection;}
}
//________________________________________________________________________
void AliAnalysisTaskProtonQA::UserCreateOutputObjects()
{
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

//  fTPCpid = dynamic_cast<AliTPCPIDResponse&>(fPIDResponse->GetTPCResponse());

  // Create histograms
  // Called once
  fListHist = new TList();
  fListHist->SetOwner();
  fListHist->SetName("acceptedDCAList");

  //fListQA = new TList();
  //fListQA->SetOwner();
  
 // fListHist->Add(fListQA);

  //fListHist->SetName("acceptedDCAList");

    if (fMultiplicityMode){
	fNBinsY = 100; fMinY = 0; fMaxY = 100;// fNBinsPt = 6; fMinPt = 0.45; fMaxPt = 1.05;
}else
{	fNBinsY = 5; fMinY = -0.5; fMaxY = 0.5;// fNBinsPt = 6; fMinPt = 0.45; fMaxPt = 1.05;
}
	


if (!gHistPrimaryProtonsDCAxyEtaPt){
gHistPrimaryProtonsDCAxyEtaPt = new TH3F("gHistPrimaryProtonsDCAxyEtaPt",
						 ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						 fNBinsY,fMinY,fMaxY,
					  	 fNBinsPt,fMinPt,fMaxPt,
						 1000,-10,10);
  gHistPrimaryProtonsDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistPrimaryProtonsDCAxyEtaPt);
}


gHistPrimaryAntiProtonsDCAxyEtaPt = new TH3F("gHistPrimaryAntiProtonsDCAxyEtaPt",
						     ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						     fNBinsY,fMinY,fMaxY,
					  	     fNBinsPt,fMinPt,fMaxPt,
						     1000,-10,10);
  gHistPrimaryAntiProtonsDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistPrimaryAntiProtonsDCAxyEtaPt);


/*gHistFake1PrimaryProtonsDCAxyEtaPt = new TH3F("gHistFake1PrimaryProtonsDCAxyEtaPt",
						 ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						 fNBinsY,fMinY,fMaxY,
					  	 fNBinsPt,fMinPt,fMaxPt,
						 1000,-10,10);
  gHistFake1PrimaryProtonsDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake1PrimaryProtonsDCAxyEtaPt);


gHistFake1PrimaryAntiProtonsDCAxyEtaPt = new TH3F("gHistFake1PrimaryAntiProtonsDCAxyEtaPt",
						     ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						     fNBinsY,fMinY,fMaxY,
					  	     fNBinsPt,fMinPt,fMaxPt,
						     1000,-10,10);
  gHistFake1PrimaryAntiProtonsDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake1PrimaryAntiProtonsDCAxyEtaPt);

gHistFake2PrimaryProtonsDCAxyEtaPt = new TH3F("gHistFake2PrimaryProtonsDCAxyEtaPt",
						 ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						 fNBinsY,fMinY,fMaxY,
					  	 fNBinsPt,fMinPt,fMaxPt,
						 1000,-10,10);
  gHistFake2PrimaryProtonsDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake2PrimaryProtonsDCAxyEtaPt);


gHistFake2PrimaryAntiProtonsDCAxyEtaPt = new TH3F("gHistFake2PrimaryAntiProtonsDCAxyEtaPt",
						     ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
						     fNBinsY,fMinY,fMaxY,
					  	     fNBinsPt,fMinPt,fMaxPt,
						     1000,-10,10);
  gHistFake2PrimaryAntiProtonsDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake2PrimaryAntiProtonsDCAxyEtaPt);*/


if (!gHistSecondaryProtonsFromWeakDCAxyEtaPt){
gHistSecondaryProtonsFromWeakDCAxyEtaPt = new TH3F("gHistSecondaryProtonsFromWeakDCAxyEtaPt",
							   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							   fNBinsY,fMinY,fMaxY,
					  	     	   fNBinsPt,fMinPt,fMaxPt,
							   1000,-10,10);
  gHistSecondaryProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistSecondaryProtonsFromWeakDCAxyEtaPt);
}

if (!gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt){
gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt = new TH3F("gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       fNBinsY,fMinY,fMaxY,
					  	     	       fNBinsPt,fMinPt,fMaxPt,
							       1000,-10,10);
  gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt);
}


/*gHistFake1SecondaryProtonsFromWeakDCAxyEtaPt = new TH3F("gHistFake1SecondaryProtonsFromWeakDCAxyEtaPt",
							   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							   fNBinsY,fMinY,fMaxY,
					  	     	   fNBinsPt,fMinPt,fMaxPt,
							   1000,-10,10);
  gHistFake1SecondaryProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake1SecondaryProtonsFromWeakDCAxyEtaPt);


gHistFake1SecondaryAntiProtonsFromWeakDCAxyEtaPt = new TH3F("gHistFake1SecondaryAntiProtonsFromWeakDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       fNBinsY,fMinY,fMaxY,
					  	     	       fNBinsPt,fMinPt,fMaxPt,
							       1000,-10,10);
  gHistFake1SecondaryAntiProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake1SecondaryAntiProtonsFromWeakDCAxyEtaPt);

gHistFake2SecondaryProtonsFromWeakDCAxyEtaPt = new TH3F("gHistFake2SecondaryProtonsFromWeakDCAxyEtaPt",
							   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							   fNBinsY,fMinY,fMaxY,
					  	     	   fNBinsPt,fMinPt,fMaxPt,
							   1000,-10,10);
  gHistFake2SecondaryProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake2SecondaryProtonsFromWeakDCAxyEtaPt);


gHistFake2SecondaryAntiProtonsFromWeakDCAxyEtaPt = new TH3F("gHistFake2SecondaryAntiProtonsFromWeakDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       fNBinsY,fMinY,fMaxY,
					  	     	       fNBinsPt,fMinPt,fMaxPt,
							       1000,-10,10);
  gHistFake2SecondaryAntiProtonsFromWeakDCAxyEtaPt->SetStats(kFALSE);
  fListHist->Add(gHistFake2SecondaryAntiProtonsFromWeakDCAxyEtaPt);*/


if (!gHistSecondaryProtonsFromHadronicDCAxyEtaPt){
gHistSecondaryProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryProtonsFromHadronicDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       fNBinsY,fMinY,fMaxY,
					  	               fNBinsPt,fMinPt,fMaxPt,
							       1000,-10,10);
  gHistSecondaryProtonsFromHadronicDCAxyEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryProtonsFromHadronicDCAxyEtaPt);
}

if (!gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt){
gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
								   fNBinsY,fMinY,fMaxY,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   1000,-10,10);
  gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt);
}

//NO FAKES
//SPD1
/*if (!gHistSecondaryFake1ProtonsFromHadronicDCAxyEtaPt){
gHistSecondaryFake1ProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryFake1ProtonsFromHadronicDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       fNBinsY,fMinY,fMaxY,
					  	               fNBinsPt,fMinPt,fMaxPt,
							       1000,-10,10);
  gHistSecondaryFake1ProtonsFromHadronicDCAxyEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryFake1ProtonsFromHadronicDCAxyEtaPt);
}

if (!gHistSecondaryFake1AntiProtonsFromHadronicDCAxyEtaPt){
gHistSecondaryFake1AntiProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryFake1AntiProtonsFromHadronicDCAxyEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
								   fNBinsY,fMinY,fMaxY,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   1000,-10,10);
  gHistSecondaryFake1AntiProtonsFromHadronicDCAxyEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryFake1AntiProtonsFromHadronicDCAxyEtaPt);
}
//SPD2
if (!gHistSecondaryFake2ProtonsFromHadronicDCAxyEtaPt){
gHistSecondaryFake2ProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryFake2ProtonsFromHadronicDCAxyEtaPt",
							       ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
							       fNBinsY,fMinY,fMaxY,
					  	               fNBinsPt,fMinPt,fMaxPt,
							       1000,-10,10);
  gHistSecondaryFake2ProtonsFromHadronicDCAxyEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryFake2ProtonsFromHadronicDCAxyEtaPt);
}

if (!gHistSecondaryFake2AntiProtonsFromHadronicDCAxyEtaPt){
gHistSecondaryFake2AntiProtonsFromHadronicDCAxyEtaPt = new TH3F("gHistSecondaryFake2AntiProtonsFromHadronicDCAxyEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{xy} [cm]",
								   fNBinsY,fMinY,fMaxY,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   1000,-10,10);
  gHistSecondaryFake2AntiProtonsFromHadronicDCAxyEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryFake2AntiProtonsFromHadronicDCAxyEtaPt);
}*/

gHistPrimaryProtonsDCAzEtaPt = new TH3F("gHistPrimaryProtonsDCAzEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
								   5,-0.5,0.5,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistPrimaryProtonsDCAzEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistPrimaryProtonsDCAzEtaPt);

gHistPrimaryAntiProtonsDCAzEtaPt = new TH3F("gHistPrimaryAntiProtonsDCAzEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
								   5,-0.5,0.5,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistPrimaryAntiProtonsDCAzEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistPrimaryAntiProtonsDCAzEtaPt);

gHistSecondaryProtonsFromWeakDCAzEtaPt = new TH3F("gHistSecondaryProtonsFromWeakDCAzEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
								   5,-0.5,0.5,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryProtonsFromWeakDCAzEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryProtonsFromWeakDCAzEtaPt);

gHistSecondaryAntiProtonsFromWeakDCAzEtaPt = new TH3F("gHistSecondaryAntiProtonsFromWeakDCAzEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
								   5,-0.5,0.5,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAzEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryAntiProtonsFromWeakDCAzEtaPt);

gHistSecondaryProtonsFromHadronicDCAzEtaPt = new TH3F("gHistSecondaryProtonsFromHadronicDCAzEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
								   5,-0.5,0.5,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryProtonsFromHadronicDCAzEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryProtonsFromHadronicDCAzEtaPt);

gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt = new TH3F("gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt",
								   ";#eta;P_{T} [GeV/c];dca_{z} [cm]",
								   5,-0.5,0.5,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt);

gHistPrimaryProtonsDCAzCentPt = new TH3F("gHistPrimaryProtonsDCAzCentPt",
								   ";Centrality;P_{T} [GeV/c];dca_{z} [cm]",
								   100,0,100,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistPrimaryProtonsDCAzCentPt->SetStats(kTRUE);
  fListHist->Add(gHistPrimaryProtonsDCAzCentPt);

gHistPrimaryAntiProtonsDCAzCentPt = new TH3F("gHistAntiPrimaryProtonsDCAzCentPt",
								   ";Centrality;P_{T} [GeV/c];dca_{z} [cm]",
								   100,0,100,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistPrimaryAntiProtonsDCAzCentPt->SetStats(kTRUE);
  fListHist->Add(gHistPrimaryAntiProtonsDCAzCentPt);

gHistSecondaryProtonsFromWeakDCAzCentPt = new TH3F("gHistSecondaryProtonsFromWeakDCAzCentPt",
								   ";Centrality;P_{T} [GeV/c];dca_{z} [cm]",
								   100,0,100,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryProtonsFromWeakDCAzCentPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryProtonsFromWeakDCAzCentPt);

gHistSecondaryAntiProtonsFromWeakDCAzCentPt = new TH3F("gHistSecondaryAntiProtonsFromWeakDCAzCentPt",
								   ";Centrality;P_{T} [GeV/c];dca_{z} [cm]",
								   100,0,100,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryAntiProtonsFromWeakDCAzCentPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryAntiProtonsFromWeakDCAzCentPt);

gHistSecondaryProtonsFromHadronicDCAzCentPt = new TH3F("gHistSecondaryProtonsFromHadronicDCAzCentPt",
								   ";Centrality;P_{T} [GeV/c];dca_{z} [cm]",
								   100,0,100,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryProtonsFromHadronicDCAzCentPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryProtonsFromHadronicDCAzCentPt);

gHistSecondaryAntiProtonsFromHadronicDCAzCentPt = new TH3F("gHistSecondaryAntiProtonsFromHadronicDCAzCentPt",
								   ";Centrality;P_{T} [GeV/c];dca_{z} [cm]",
								   100,0,100,
					  	                   fNBinsPt,fMinPt,fMaxPt,
								   500,0,10);
  gHistSecondaryAntiProtonsFromHadronicDCAzCentPt->SetStats(kTRUE);
  fListHist->Add(gHistSecondaryAntiProtonsFromHadronicDCAzCentPt);

  //dEdx plots
if (!gHistdEdxP){
  gHistdEdxP = new TH2F("gHistdEdxP","dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
  gHistdEdxP->SetStats(kFALSE);
  fListHist->Add(gHistdEdxP);
}

if (!gHistProtonsdEdxP){
 gHistProtonsdEdxP = new TH2F("gHistProtonsdEdxP","Accepted protons dE/dx (TPC); P [GeV/c]; dE/dx [a.u]",1000,0.05,20.05,600,0,600);
 gHistProtonsdEdxP->SetStats(kFALSE);
 fListHist->Add(gHistProtonsdEdxP);
}

PostData(1, fListHist);
}

void AliAnalysisTaskProtonQA::UserExec(Option_t *){
// Main loop
  // Called for each event


// Create MC part and stack
           AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
           if (!eventHandler) {
	     Printf("ERROR: Could not retrieve MC event handler");
	     return;
           }
           AliMCEvent* mcEvent = eventHandler->MCEvent();
           if (!mcEvent) {
	     Printf("ERROR: Could not retrieve MC event");
	     return;
           }
	  AliStack* mcStack = mcEvent->Stack(); 

//--------------------------------------------------------------------------

// Create ESD part and stack
  AliVEvent* lEvent = InputEvent();
  if (!lEvent) {
    Printf("ERROR: Event not available");
    return;
  }

  AliESDEvent* lESDEvent = (AliESDEvent*)lEvent;

//--------------------------------------------------------------------------
AliCentrality *esdCentrality = lESDEvent->GetCentrality();
Float_t nTracklets = esdCentrality->GetCentralityPercentile("V0M");


  AliCentrality *centrality = lESDEvent->GetCentrality();

  Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if ((fUsePhysicsSelection)&&(!isSelected)) return;

  if(centrality->GetCentralityPercentile("V0M") > MAXCent || centrality->GetCentralityPercentile("V0M") < MINCent) return;
//Printf("Centrality %f ", centrality->GetCentralityPercentile("V0M"));

// Create vertex
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
//--------------------------------------------------------------------------
if (TMath::Abs(zvtx)<fMaxPrimaryVtxPosZ){// event selections

 Double_t dca[2] = {0.0,0.0};
 //ESD track loop
  Int_t nGoodTracks = lESDEvent->GetNumberOfTracks();
  TArrayI labelArray(nGoodTracks);
  Int_t labelCounter = 0;
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
    AliESDtrack* track = lESDEvent->GetTrack(iTracks);
    if(!track) continue;

    Int_t label = TMath::Abs(track->GetLabel());
    //if (label <= 0) continue; 
    if(IsLabelUsed(labelArray,label)) continue;
    labelArray.AddAt(label,labelCounter);
    labelCounter += 1;
    if(label > mcStack->GetNtrack()) continue;

    TParticle *particle = mcStack->Particle(label);
    if(!particle) continue;
    if(TMath::Abs(particle->Eta()) > 1.0) continue;//acceptance

    AliESDtrack trackTPC;

    Double_t gPt,gPy,gPx,gPz,gP;
    Float_t dcaXY = 0.0, dcaZ = 0.0;
    
    //TPC
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) continue;
    gPt = tpcTrack->Pt();
    gPx = tpcTrack->Px();
    gPy = tpcTrack->Py();
    gPz = tpcTrack->Pz();
    AliExternalTrackParam *in = (AliExternalTrackParam *)track->GetInnerParam();
    if (in) gP = in->GetP();

    //FullHybrid
    AliExternalTrackParam cParam;
    track->RelateToVertex(spdVtx,
			    lESDEvent->GetMagneticField(),
			    100.,&cParam);
    track->GetImpactParameters(dcaXY,dcaZ);
    dca[0] = dcaXY; dca[1] = dcaZ;
    
	if(IsAccepted(track)){
		
						gHistdEdxP->Fill(gP,track->GetTPCsignal());
	if(((gPt > fMinPt) && (gPt < fMaxPt)) && ((gP > fMinPt) && (gP < fMaxPt)) && ((Rapidity(gPx,gPy,gPz,3) > -0.5) && (Rapidity(gPx,gPy,gPz,3) < 0.5))) {
	if(IsPrimary(lESDEvent,spdVtx,track)) {
	if(IsProton(track)) {

			//if ((tpcTrack->P()) > (track->P())) Printf("tpcTrack p P %lf > track P %lf",tpcTrack->P(),track->Pt());
						gHistProtonsdEdxP->Fill(gP,track->GetTPCsignal());
			if(label <= mcStack->GetNprimary()) {
			if(track->Charge() > 0) {
gHistPrimaryProtonsDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[1]);
gHistPrimaryProtonsDCAzCentPt->Fill(nTracklets,gPt,dca[1]);

			//Printf("Phase space primary protons: Pt = %lf, P = %lf, Rapidity = %lf)",gPt,gP,Rapidity(gPx,gPy,gPz,3));
			if (fMultiplicityMode){
			gHistPrimaryProtonsDCAxyEtaPt->Fill(nTracklets,gPt,dca[0]);

			}else gHistPrimaryProtonsDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[0]);
			}//accepted primary protons
			else if(track->Charge() < 0) {
gHistPrimaryAntiProtonsDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[1]);
gHistPrimaryAntiProtonsDCAzCentPt->Fill(nTracklets,gPt,dca[1]);

			//Printf("Phase space primary antiprotons: Pt = %lf, P = %lf, Rapidity = %lf)",gPt,gP,Rapidity(gPx,gPy,gPz,3));
			if (fMultiplicityMode){
			gHistPrimaryAntiProtonsDCAxyEtaPt->Fill(nTracklets,gPt,dca[0]);
			}else gHistPrimaryAntiProtonsDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[0]);
			}//accepted primary antiprotons
			}//accepted primary particles
				else if(label > mcStack->GetNprimary()) {
		      		/*  Int_t lPartMother = -1;
	      			  Int_t motherPDGCode = -1;
	      			if(particle) {
				  lPartMother = particle->GetFirstMother();
				  TParticle *motherParticle = mcStack->Particle(lPartMother);
				if(motherParticle) motherPDGCode = motherParticle->GetPdgCode();   TOTO JE MYSLIM NEPOTREBNE
	      			}
	      
	      			if(fMCProcessIdFlag)
				if(particle->GetUniqueID() != fMCProcessId) continue;
	      			if(fMotherParticlePDGCodeFlag)
				if(TMath::Abs(motherPDGCode) != fMotherParticlePDGCode) continue;*/

				if(track->Charge() > 0) {
				if(particle->GetUniqueID() == 4) {
gHistSecondaryProtonsFromWeakDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[1]);
gHistSecondaryProtonsFromWeakDCAzCentPt->Fill(nTracklets,gPt,dca[1]);
				//Printf("Phase space secondary from weak decay protons: Pt = %lf, P = %lf, Rapidity = %lf)",gPt,gP,Rapidity(gPx,gPy,gPz,3));
				if (fMultiplicityMode){
				gHistSecondaryProtonsFromWeakDCAxyEtaPt->Fill(nTracklets,gPt,dca[0]);
				}else gHistSecondaryProtonsFromWeakDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[0]);
				}//end GetUniqueID() == 4
				if(particle->GetUniqueID() == 13) {
gHistSecondaryProtonsFromHadronicDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[1]);
gHistSecondaryProtonsFromHadronicDCAzCentPt->Fill(nTracklets,gPt,dca[1]);
				//Printf("Phase space secondary from material protons: Pt = %lf, P = %lf, Rapidity = %lf)",gPt,gP,Rapidity(gPx,gPy,gPz,3));
				if (fMultiplicityMode){
				gHistSecondaryProtonsFromHadronicDCAxyEtaPt->Fill(nTracklets,gPt,dca[0]);
				}else gHistSecondaryProtonsFromHadronicDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[0]);
				}//end GetUniqueID() == 13
}//accepted secondary protons
				if(track->Charge() < 0) {
gHistSecondaryAntiProtonsFromWeakDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[1]);
gHistSecondaryAntiProtonsFromWeakDCAzCentPt->Fill(nTracklets,gPt,dca[1]);
				if(particle->GetUniqueID() == 4) {
				//Printf("Phase space secondary from weak decay ANTI protons: Pt = %lf, P = %lf, Rapidity = %lf)",gPt,gP,Rapidity(gPx,gPy,gPz,3));
				if (fMultiplicityMode){
				gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt->Fill(nTracklets,gPt,dca[0]);
				}else gHistSecondaryAntiProtonsFromWeakDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[0]);
				}//end GetUniqueID() == 4
				if(particle->GetUniqueID() == 13) {
gHistSecondaryAntiProtonsFromHadronicDCAzEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[1]);
gHistSecondaryAntiProtonsFromHadronicDCAzCentPt->Fill(nTracklets,gPt,dca[1]);
				//Printf("Phase space secondary from material ANTI protons: Pt = %lf, P = %lf, Rapidity = %lf)",gPt,gP,Rapidity(gPx,gPy,gPz,3));
				if (fMultiplicityMode){
				gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt->Fill(nTracklets,gPt,dca[0]);
				}else gHistSecondaryAntiProtonsFromHadronicDCAxyEtaPt->Fill(Rapidity(tpcTrack->Px(),tpcTrack->Py(),tpcTrack->Pz(),3),gPt,dca[0]);
				}//end GetUniqueID() == 13
}//accepted secondary antiprotons    
	}//accepted secondary particles
		}//end IsAccepted
			}//end IsPrimary	
			 	}//end Phase space
					}//end IsProton
//Printf("****************************************************");
						}//end iTracks loop

							}//end event selection

}


/*void AliAnalysisTaskProtonQA::SetRecoInfo()
{
 //
  // Set reconstruction information
  //
 
  //reset information
  fRecoPass=0;
 
  //Get UserInfo from the current tree
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
  if (!inputHandler) return;
 
  TList *uiList = inputHandler->GetUserInfo();
  AliProdInfo prodInfo(uiList);
  prodInfo.List();
 
  TTree *tree= (TTree*)inputHandler->GetTree();
  TFile *file= (TFile*)tree->GetCurrentFile();
  if (!file) {
    AliError("Current file not found, cannot set reconstruction information");
    return;
  } else {
    TString fileName(file->GetName());
    fPIDResponse->SetCurrentFile(fileName.Data());
  }
 
  fPIDResponse->SetCurrentAliRootRev(prodInfo.GetAlirootSvnVersion());
 
  if (prodInfo.IsMC() == kTRUE) fIsMC=kTRUE;         // protection if user didn't use macro switch
  if ( (prodInfo.IsMC() == kFALSE) && (fIsMC == kFALSE) ) {      // reco pass is needed only for data
 
    if (fUserDataRecoPass > -1) {
      AliInfo(Form("Data reconstruction pass is user specified. Setting pass #: %d",fUserDataRecoPass));
      fRecoPass = fUserDataRecoPass;
    } else {
      fRecoPass = prodInfo.GetRecoPass();
      if (fRecoPass < 0) {   // as last resort we find pass from file name (UGLY, but not stored in ESDs/AODs before LHC12d )
        TString fileName(file->GetName());
        if (fileName.Contains("pass1") ) {
          fRecoPass=1;
        } else if (fileName.Contains("pass2") ) {
          fRecoPass=2;
        } else if (fileName.Contains("pass3") ) {
          fRecoPass=3;
        } else if (fileName.Contains("pass4") ) {
          fRecoPass=4;
        } else if (fileName.Contains("pass5") ) {
          fRecoPass=5;
        }
    }
    }
    if (fRecoPass <= 0) {
      AliError(" ******** Failed to find reconstruction pass number *********");
      AliError(" ******** PID information loaded for 'pass 0': parameters unreliable ******");
      AliError("      --> If these are MC data: please set kTRUE first argument of AddTaskPIDResponse");
      AliError("      --> If these are real data: ");
      AliError("          (a) please insert pass number inside the path of your local file  OR");
      AliError("          (b) specify reconstruction pass number when adding PIDResponse task");
      AliFatal(" no pass number specified for this data file, abort. Asserting AliFatal");
    }
  }
}*/
//________________________________________________________________________
Bool_t AliAnalysisTaskProtonQA::IsLabelUsed(TArrayI labelArray, 
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
Bool_t AliAnalysisTaskProtonQA::IsProton(AliESDtrack *track) {
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

Bool_t AliAnalysisTaskProtonQA::IsAccepted(AliESDtrack* track) {
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

Bool_t AliAnalysisTaskProtonQA::IsPrimary(AliESDEvent *esd,
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
Double_t AliAnalysisTaskProtonQA::Rapidity(Double_t gPx, 
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
