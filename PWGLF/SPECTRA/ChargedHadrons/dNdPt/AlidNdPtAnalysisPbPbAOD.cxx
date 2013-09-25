/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
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
//------------------------------------------------------------------------------
// AlidNdPtAnalysisPbPbAOD class. 
// 
// Author: P. Luettig, 15.05.2013
//------------------------------------------------------------------------------


#include "AlidNdPtAnalysisPbPbAOD.h"


using namespace std;

ClassImp(AlidNdPtAnalysisPbPbAOD)

// dummy constructor
AlidNdPtAnalysisPbPbAOD::AlidNdPtAnalysisPbPbAOD() : AliAnalysisTaskSE(),
fOutputList(0),
// Histograms
hPt(0),
hMCPt(0),
hnZvPtEtaCent(0),
hnMCRecPrimZvPtEtaCent(0),
hnMCGenZvPtEtaCent(0),
hnMCRecSecZvPtEtaCent(0),
hEventStatistics(0),
hEventStatisticsCentrality(0),
hMCEventStatisticsCentrality(0),
hAllEventStatisticsCentrality(0),
hEventStatisticsCentralityTrigger(0),
hnZvMultCent(0),
hTriggerStatistics(0),
hMCTrackPdgCode(0),
hMCTrackStatusCode(0),
hCharge(0),
hMCCharge(0),
hMCPdgPt(0),
hMCHijingPrim(0),
hDCAPtAll(0),
hDCAPtAccepted(0),
hMCDCAPtSecondary(0),
hMCDCAPtPrimary(0),
hnCrossedRowsClustersChiPtEtaPhiAll(0),
hnCrossedRowsClustersChiPtEtaPhiAcc(0),
//global
bIsMonteCarlo(0),
// event cut variables
dCutMaxZVertex(10.),  
// track kinematic cut variables
dCutPtMin(0.15),
dCutPtMax(1000.),
dCutEtaMin(-0.8),
dCutEtaMax(0.8),    
// track quality cut variables
bCutRequireTPCRefit(kTRUE),
dCutMinNumberOfCrossedRows(120.),
dCutMinRatioCrossedRowsOverFindableClustersTPC(0.8),
dCutMaxChi2PerClusterTPC(4.),
dCutMaxFractionSharedTPCClusters(0.4),
dCutMaxDCAToVertexZ(3.0),
dCutMaxDCAToVertexXY(3.0),
bCutRequireITSRefit(kTRUE),
dCutMaxChi2PerClusterITS(36.),
dCutDCAToVertex2D(kFALSE),
dCutRequireSigmaToVertex(kFALSE),
dCutMaxDCAToVertexXYPtDepPar0(0.0182),
dCutMaxDCAToVertexXYPtDepPar1(0.0350),
dCutMaxDCAToVertexXYPtDepPar2(1.01),
bCutAcceptKinkDaughters(kFALSE),
dCutMaxChi2TPCConstrainedGlobal(36.),
// binning for THnSparse
fMultNbins(0),
fPtNbins(0),
fPtCorrNbins(0),
fEtaNbins(0),
fZvNbins(0),
fCentralityNbins(0),
fBinsMult(0),
fBinsPt(0),
fBinsPtCorr(0),
fBinsEta(0),
fBinsZv(0),
fBinsCentrality(0)
{
  
  fMultNbins = 0;
  fPtNbins = 0;
  fPtCorrNbins = 0;
  fEtaNbins = 0;
  fZvNbins = 0;
  fCentralityNbins = 0;
  fBinsMult = 0;
  fBinsPt = 0;
  fBinsPtCorr = 0;
  fBinsEta = 0;
  fBinsEta = 0;
  fBinsZv = 0;
  fBinsCentrality = 0;
  
}

AlidNdPtAnalysisPbPbAOD::AlidNdPtAnalysisPbPbAOD(const char *name) : AliAnalysisTaskSE(name),
fOutputList(0),
// Histograms
hPt(0),
hMCPt(0),
hnZvPtEtaCent(0),
hnMCRecPrimZvPtEtaCent(0),
hnMCGenZvPtEtaCent(0),
hnMCRecSecZvPtEtaCent(0),
hEventStatistics(0),
hEventStatisticsCentrality(0),
hMCEventStatisticsCentrality(0),
hAllEventStatisticsCentrality(0),
hEventStatisticsCentralityTrigger(0),
hnZvMultCent(0),
hTriggerStatistics(0),
hMCTrackPdgCode(0),
hMCTrackStatusCode(0),
hCharge(0),
hMCCharge(0),
hMCPdgPt(0),
hMCHijingPrim(0),
hDCAPtAll(0),
hDCAPtAccepted(0),
hMCDCAPtSecondary(0),
hMCDCAPtPrimary(0),
hnCrossedRowsClustersChiPtEtaPhiAll(0),
hnCrossedRowsClustersChiPtEtaPhiAcc(0),
//global
bIsMonteCarlo(0),
// event cut variables
dCutMaxZVertex(10.),  
// track kinematic cut variables
dCutPtMin(0.15),
dCutPtMax(200.),
dCutEtaMin(-0.8),
dCutEtaMax(0.8),    
// track quality cut variables
bCutRequireTPCRefit(kTRUE),
dCutMinNumberOfCrossedRows(120.),
dCutMinRatioCrossedRowsOverFindableClustersTPC(0.8),
dCutMaxChi2PerClusterTPC(4.),
dCutMaxFractionSharedTPCClusters(0.4),
dCutMaxDCAToVertexZ(3.0),
dCutMaxDCAToVertexXY(3.0),
bCutRequireITSRefit(kTRUE),
dCutMaxChi2PerClusterITS(36.),
dCutDCAToVertex2D(kFALSE),
dCutRequireSigmaToVertex(kFALSE),
dCutMaxDCAToVertexXYPtDepPar0(0.0182),
dCutMaxDCAToVertexXYPtDepPar1(0.0350),
dCutMaxDCAToVertexXYPtDepPar2(1.01),
bCutAcceptKinkDaughters(kFALSE),
dCutMaxChi2TPCConstrainedGlobal(36.),
// binning for THnSparse
fMultNbins(0),
fPtNbins(0),
fPtCorrNbins(0),
fEtaNbins(0),
fZvNbins(0),
fCentralityNbins(0),
fBinsMult(0),
fBinsPt(0),
fBinsPtCorr(0),
fBinsEta(0),
fBinsZv(0),
fBinsCentrality(0)
{
  fMultNbins = 0;
  fPtNbins = 0;
  fPtCorrNbins = 0;
  fEtaNbins = 0;
  fZvNbins = 0;
  fCentralityNbins = 0;
  fBinsMult = 0;
  fBinsPt = 0;
  fBinsPtCorr = 0;
  fBinsEta = 0;
  fBinsEta = 0;
  fBinsZv = 0;
  fBinsCentrality = 0;
  
  DefineOutput(1, TList::Class());
}

// destructor
AlidNdPtAnalysisPbPbAOD::~AlidNdPtAnalysisPbPbAOD()
{
  
  if(hnZvPtEtaCent) delete hnZvPtEtaCent; hnZvPtEtaCent = 0;
  if(hPt) delete hPt; hPt = 0;
  if(hnMCRecPrimZvPtEtaCent) delete hnMCRecPrimZvPtEtaCent; hnMCRecPrimZvPtEtaCent = 0;
  if(hnMCGenZvPtEtaCent) delete hnMCGenZvPtEtaCent; hnMCGenZvPtEtaCent = 0;
  if(hnMCRecSecZvPtEtaCent) delete hnMCRecSecZvPtEtaCent; hnMCRecSecZvPtEtaCent = 0;
  if(hMCPt) delete hMCPt; hMCPt = 0;
  if(hEventStatistics) delete hEventStatistics; hEventStatistics = 0;
  if(hEventStatisticsCentrality) delete hEventStatisticsCentrality; hEventStatisticsCentrality = 0;
  if(hMCEventStatisticsCentrality) delete hMCEventStatisticsCentrality; hMCEventStatisticsCentrality = 0;
  if(hAllEventStatisticsCentrality) delete hAllEventStatisticsCentrality; hAllEventStatisticsCentrality = 0;
  if(hEventStatisticsCentralityTrigger) delete hEventStatisticsCentralityTrigger; hEventStatisticsCentralityTrigger = 0;
  if(hnZvMultCent) delete hnZvMultCent; hnZvMultCent = 0;
  if(hTriggerStatistics) delete hTriggerStatistics; hTriggerStatistics = 0;
  if(hMCTrackPdgCode) delete hMCTrackPdgCode; hMCTrackPdgCode = 0;
  if(hMCTrackStatusCode) delete hMCTrackStatusCode; hMCTrackStatusCode = 0;
  if(hCharge) delete hCharge; hCharge = 0;
  if(hMCCharge) delete hMCCharge; hMCCharge = 0;
  if(hMCPdgPt) delete hMCPdgPt; hMCPdgPt = 0;
  if(hMCHijingPrim) delete hMCHijingPrim; hMCHijingPrim = 0;
  if(hDCAPtAll) delete hDCAPtAll; hDCAPtAll = 0;
  if(hDCAPtAccepted) delete hDCAPtAccepted; hDCAPtAccepted = 0;
  if(hMCDCAPtSecondary) delete hMCDCAPtSecondary; hMCDCAPtSecondary = 0;
  if(hMCDCAPtPrimary) delete hMCDCAPtPrimary; hMCDCAPtPrimary = 0;
  if(hMCDCAPtSecondary) delete hMCDCAPtSecondary; hMCDCAPtSecondary = 0;
  if(hMCDCAPtPrimary) delete hMCDCAPtPrimary; hMCDCAPtPrimary = 0;
  if(hnCrossedRowsClustersChiPtEtaPhiAll) delete hnCrossedRowsClustersChiPtEtaPhiAll; hnCrossedRowsClustersChiPtEtaPhiAll = 0;
  if(hnCrossedRowsClustersChiPtEtaPhiAcc) delete hnCrossedRowsClustersChiPtEtaPhiAcc; hnCrossedRowsClustersChiPtEtaPhiAcc = 0;
}

void AlidNdPtAnalysisPbPbAOD::UserCreateOutputObjects()
{
  // create all output histograms here
  OpenFile(1, "RECREATE");
  
  fOutputList = new TList();
  fOutputList->SetOwner();
  
  //define default binning
  Double_t binsMultDefault[48] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,19.5, 20.5, 30.5, 40.5 , 50.5 , 60.5 , 70.5 , 80.5 , 90.5 , 100.5,200.5, 300.5, 400.5, 500.5, 600.5, 700.5, 800.5, 900.5, 1000.5, 2000.5, 3000.5, 4000.5, 5000.5, 6000.5, 7000.5, 8000.5, 9000.5, 10000.5 }; 
  Double_t binsPtDefault[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
  Double_t binsPtCorrDefault[37] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,3.0,4.0,200.0};    
  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
  Double_t binsCentralityDefault[12] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};  
  
  // if no binning is set, use the default
  if (!fBinsMult)	{ SetBinsMult(48,binsMultDefault); }
  if (!fBinsPt)		{ SetBinsPt(82,binsPtDefault); }
  if (!fBinsPtCorr)	{ SetBinsPtCorr(37,binsPtCorrDefault); }
  if (!fBinsEta)	{ SetBinsEta(31,binsEtaDefault); }
  if (!fBinsZv)		{ SetBinsZv(13,binsZvDefault); }  
  if (!fBinsCentrality)	{ SetBinsCentrality(12,binsCentralityDefault); }
  
  Int_t binsZvPtEtaCent[4]={fZvNbins-1,fPtNbins-1,fEtaNbins-1,fCentralityNbins-1};
  Int_t binsZvMultCent[3]={fZvNbins-1,fMultNbins-1,fCentralityNbins-1};
  
  // define Histograms
  hnZvPtEtaCent = new THnSparseF("hnZvPtEtaCent","Zv:Pt:Eta:Centrality",4,binsZvPtEtaCent);
  hnZvPtEtaCent->SetBinEdges(0,fBinsZv);
  hnZvPtEtaCent->SetBinEdges(1,fBinsPt);
  hnZvPtEtaCent->SetBinEdges(2,fBinsEta);
  hnZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  hnZvPtEtaCent->GetAxis(0)->SetTitle("Zv (cm)");
  hnZvPtEtaCent->GetAxis(1)->SetTitle("Pt (GeV/c)");
  hnZvPtEtaCent->GetAxis(2)->SetTitle("Eta");
  hnZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  hnZvPtEtaCent->Sumw2();
  
  hnMCRecPrimZvPtEtaCent = new THnSparseF("hnMCRecPrimZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtEtaCent);
  hnMCRecPrimZvPtEtaCent->SetBinEdges(0,fBinsZv);
  hnMCRecPrimZvPtEtaCent->SetBinEdges(1,fBinsPt);
  hnMCRecPrimZvPtEtaCent->SetBinEdges(2,fBinsEta);
  hnMCRecPrimZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  hnMCRecPrimZvPtEtaCent->GetAxis(0)->SetTitle("MC Zv (cm)");
  hnMCRecPrimZvPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  hnMCRecPrimZvPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  hnMCRecPrimZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  hnMCRecPrimZvPtEtaCent->Sumw2();
  
  hnMCGenZvPtEtaCent = new THnSparseF("hnMCGenZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtEtaCent);
  hnMCGenZvPtEtaCent->SetBinEdges(0,fBinsZv);
  hnMCGenZvPtEtaCent->SetBinEdges(1,fBinsPt);
  hnMCGenZvPtEtaCent->SetBinEdges(2,fBinsEta);
  hnMCGenZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  hnMCGenZvPtEtaCent->GetAxis(0)->SetTitle("MC Zv (cm)");
  hnMCGenZvPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  hnMCGenZvPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  hnMCGenZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  hnMCGenZvPtEtaCent->Sumw2();
  
  hnMCRecSecZvPtEtaCent = new THnSparseF("hnMCRecSecZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtEtaCent);
  hnMCRecSecZvPtEtaCent->SetBinEdges(0,fBinsZv);
  hnMCRecSecZvPtEtaCent->SetBinEdges(1,fBinsPt);
  hnMCRecSecZvPtEtaCent->SetBinEdges(2,fBinsEta);
  hnMCRecSecZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  hnMCRecSecZvPtEtaCent->GetAxis(0)->SetTitle("MC Sec Zv (cm)");
  hnMCRecSecZvPtEtaCent->GetAxis(1)->SetTitle("MC Sec Pt (GeV/c)");
  hnMCRecSecZvPtEtaCent->GetAxis(2)->SetTitle("MC Sec Eta");
  hnMCRecSecZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  hnMCRecSecZvPtEtaCent->Sumw2();
  
  hPt = new TH1F("hPt","hPt",2000,0,200);
  hPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPt->GetYaxis()->SetTitle("dN/dp_{T}");
  hPt->Sumw2();
  
  hMCPt = new TH1F("hMCPt","hMCPt",2000,0,200);
  hMCPt->GetXaxis()->SetTitle("MC p_{T} (GeV/c)");
  hMCPt->GetYaxis()->SetTitle("dN/dp_{T}");
  hMCPt->Sumw2();
  
  hEventStatistics = new TH1F("hEventStatistics","hEventStatistics",10,0,10);
  hEventStatistics->GetYaxis()->SetTitle("number of events");
  hEventStatistics->SetBit(TH1::kCanRebin);
  
  hEventStatisticsCentrality = new TH1F("hEventStatisticsCentrality","hEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  hEventStatisticsCentrality->GetYaxis()->SetTitle("number of events");
  
  hMCEventStatisticsCentrality = new TH1F("hMCEventStatisticsCentrality","hMCEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  hMCEventStatisticsCentrality->GetYaxis()->SetTitle("number of MC events");
  
  hAllEventStatisticsCentrality = new TH1F("hAllEventStatisticsCentrality","hAllEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  hAllEventStatisticsCentrality->GetYaxis()->SetTitle("number of events");
  
  hEventStatisticsCentralityTrigger = new TH2F("hEventStatisticsCentralityTrigger","hEventStatisticsCentralityTrigger;centrality;trigger",100,0,100,3,0,3);
  hEventStatisticsCentralityTrigger->Sumw2();
  
  hnZvMultCent = new THnSparseF("hnZvMultCent","Zv:mult:Centrality",3,binsZvMultCent);
  hnZvMultCent->SetBinEdges(0,fBinsZv);
  hnZvMultCent->SetBinEdges(1,fBinsMult);
  hnZvMultCent->SetBinEdges(2,fBinsCentrality);
  hnZvMultCent->GetAxis(0)->SetTitle("Zv (cm)");
  hnZvMultCent->GetAxis(1)->SetTitle("N_{acc}");
  hnZvMultCent->GetAxis(2)->SetTitle("Centrality");
  hnZvMultCent->Sumw2();
  
  hTriggerStatistics = new TH1F("hTriggerStatistics","hTriggerStatistics",10,0,10);
  hTriggerStatistics->GetYaxis()->SetTitle("number of events");
  
  hMCTrackPdgCode = new TH1F("hMCTrackPdgCode","hMCTrackPdgCode",100,0,10);
  hMCTrackPdgCode->GetYaxis()->SetTitle("number of tracks");
  hMCTrackPdgCode->SetBit(TH1::kCanRebin);
  
  hMCTrackStatusCode = new TH1F("hMCTrackStatusCode","hMCTrackStatusCode",100,0,10);
  hMCTrackStatusCode->GetYaxis()->SetTitle("number of tracks");
  hMCTrackStatusCode->SetBit(TH1::kCanRebin);
  
  hCharge = new TH1F("hCharge","hCharge",30, -5, 5);
  hCharge->GetXaxis()->SetTitle("Charge");
  hCharge->GetYaxis()->SetTitle("number of tracks");
  
  hMCCharge = new TH1F("hMCCharge","hMCCharge",30, -5, 5);
  hMCCharge->GetXaxis()->SetTitle("MC Charge");
  hMCCharge->GetYaxis()->SetTitle("number of tracks");
  
  hMCPdgPt = new TH2F("hMCPdgPt","hMCPdgPt",fPtNbins-1, fBinsPt, 100,0,100);
  hMCPdgPt->GetYaxis()->SetTitle("particle");
  hMCPdgPt->GetXaxis()->SetTitle("Pt (GeV/c)");
  
  hMCHijingPrim = new TH1F("hMCHijingPrim","hMCHijingPrim",2,0,2);
  hMCPdgPt->GetYaxis()->SetTitle("number of particles");
  
  
  
  Int_t binsDCAxyDCAzPtEtaPhi[6] = { 200,200, fPtNbins-1, fEtaNbins-1, 36, fCentralityNbins-1};
  Double_t minDCAxyDCAzPtEtaPhi[6] = { -5, -5, 0, -1.5, 0., 0, };
  Double_t maxDCAxyDCAzPtEtaPhi[6] = { 5., 5., 100, 1.5, 2.*TMath::Pi(), 100};
  
  hDCAPtAll = new THnSparseF("hDCAPtAll","hDCAPtAll",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  hDCAPtAccepted = new THnSparseF("hDCAPtAccepted","hDCAPtAccepted",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  hMCDCAPtSecondary = new THnSparseF("hMCDCAPtSecondary","hMCDCAPtSecondary",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  hMCDCAPtPrimary = new THnSparseF("hMCDCAPtPrimary","hMCDCAPtPrimary",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  
  hDCAPtAll->SetBinEdges(2, fBinsPt);
  hDCAPtAccepted->SetBinEdges(2, fBinsPt);
  hMCDCAPtSecondary->SetBinEdges(2, fBinsPt);
  hMCDCAPtPrimary->SetBinEdges(2, fBinsPt);
  
  hDCAPtAll->SetBinEdges(3, fBinsEta);
  hDCAPtAccepted->SetBinEdges(3, fBinsEta);
  hMCDCAPtSecondary->SetBinEdges(3, fBinsEta);
  hMCDCAPtPrimary->SetBinEdges(3, fBinsEta);
  
  hDCAPtAll->SetBinEdges(5, fBinsCentrality);
  hDCAPtAccepted->SetBinEdges(5, fBinsCentrality);
  hMCDCAPtSecondary->SetBinEdges(5, fBinsCentrality);
  hMCDCAPtPrimary->SetBinEdges(5, fBinsCentrality);
  
  hDCAPtAll->Sumw2();
  hDCAPtAccepted->Sumw2();
  hMCDCAPtSecondary->Sumw2();
  hMCDCAPtPrimary->Sumw2();
  
  hDCAPtAll->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  hDCAPtAll->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  hDCAPtAll->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  hDCAPtAll->GetAxis(3)->SetTitle("#eta");
  hDCAPtAll->GetAxis(4)->SetTitle("#phi");
  hDCAPtAll->GetAxis(5)->SetTitle("Centrality");
  
  hDCAPtAccepted->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  hDCAPtAccepted->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  hDCAPtAccepted->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  hDCAPtAccepted->GetAxis(3)->SetTitle("#eta");
  hDCAPtAccepted->GetAxis(4)->SetTitle("#phi");
  hDCAPtAccepted->GetAxis(5)->SetTitle("Centrality");
  
  hMCDCAPtSecondary->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  hMCDCAPtSecondary->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  hMCDCAPtSecondary->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  hMCDCAPtSecondary->GetAxis(3)->SetTitle("#eta");
  hMCDCAPtSecondary->GetAxis(4)->SetTitle("#phi");
  hMCDCAPtSecondary->GetAxis(5)->SetTitle("Centrality");
  
  hMCDCAPtPrimary->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  hMCDCAPtPrimary->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  hMCDCAPtPrimary->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  hMCDCAPtPrimary->GetAxis(3)->SetTitle("#eta");
  hMCDCAPtPrimary->GetAxis(4)->SetTitle("#phi");
  hMCDCAPtPrimary->GetAxis(5)->SetTitle("Centrality");
  
  Int_t binsRowsClusterChiPtEtaPhi[7] = { 32,32, 100, fPtNbins-1, fEtaNbins-1, 36, fCentralityNbins-1};
  Double_t minRowsClusterChiPtEtaPhi[7] = { 0, 0, 0, 0, -1.5, 0., 0.};
  Double_t maxRowsClusterChiPtEtaPhi[7] = { 159, 159, 10, 100, 1.5, 2.*TMath::Pi(), 100.};
  
  hnCrossedRowsClustersChiPtEtaPhiAll = new THnSparseF("hnCrossedRowsClustersChiPtEtaPhiAll","hnCrossedRowsClustersChiPtEtaPhiAll",7,binsRowsClusterChiPtEtaPhi, minRowsClusterChiPtEtaPhi, maxRowsClusterChiPtEtaPhi);
  
  hnCrossedRowsClustersChiPtEtaPhiAcc = new THnSparseF("hnCrossedRowsClustersChiPtEtaPhiAcc","hnCrossedRowsClustersChiPtEtaPhiAcc",7,binsRowsClusterChiPtEtaPhi, minRowsClusterChiPtEtaPhi, maxRowsClusterChiPtEtaPhi);
  
  hnCrossedRowsClustersChiPtEtaPhiAll->Sumw2();
  hnCrossedRowsClustersChiPtEtaPhiAll->SetBinEdges(3, fBinsPt);
  hnCrossedRowsClustersChiPtEtaPhiAll->SetBinEdges(4, fBinsEta);
  hnCrossedRowsClustersChiPtEtaPhiAll->SetBinEdges(6, fBinsCentrality);
  hnCrossedRowsClustersChiPtEtaPhiAll->GetAxis(0)->SetTitle("NcrossedRows before Cut");
  hnCrossedRowsClustersChiPtEtaPhiAll->GetAxis(1)->SetTitle("Nclusters before Cut");
  hnCrossedRowsClustersChiPtEtaPhiAll->GetAxis(2)->SetTitle("#Chi^{2}/cluster before Cut");
  hnCrossedRowsClustersChiPtEtaPhiAll->GetAxis(3)->SetTitle("p_{T} (GeV/c)");
  hnCrossedRowsClustersChiPtEtaPhiAll->GetAxis(4)->SetTitle("#eta");
  hnCrossedRowsClustersChiPtEtaPhiAll->GetAxis(5)->SetTitle("#phi");
  hnCrossedRowsClustersChiPtEtaPhiAll->GetAxis(6)->SetTitle("Centrality");
  
  hnCrossedRowsClustersChiPtEtaPhiAcc->Sumw2();
  hnCrossedRowsClustersChiPtEtaPhiAcc->SetBinEdges(3, fBinsPt);
  hnCrossedRowsClustersChiPtEtaPhiAcc->SetBinEdges(4, fBinsEta);
  hnCrossedRowsClustersChiPtEtaPhiAcc->SetBinEdges(6, fBinsCentrality);
  hnCrossedRowsClustersChiPtEtaPhiAcc->GetAxis(0)->SetTitle("NcrossedRows after Cut");
  hnCrossedRowsClustersChiPtEtaPhiAcc->GetAxis(1)->SetTitle("Nclusters after Cut");
  hnCrossedRowsClustersChiPtEtaPhiAcc->GetAxis(2)->SetTitle("#Chi^{2}/cluster after Cut");
  hnCrossedRowsClustersChiPtEtaPhiAcc->GetAxis(3)->SetTitle("p_{T} (GeV/c)");
  hnCrossedRowsClustersChiPtEtaPhiAcc->GetAxis(4)->SetTitle("#eta");
  hnCrossedRowsClustersChiPtEtaPhiAcc->GetAxis(5)->SetTitle("#phi");
  hnCrossedRowsClustersChiPtEtaPhiAcc->GetAxis(6)->SetTitle("Centrality");
  
  // Add Histos, Profiles etc to List
  fOutputList->Add(hnZvPtEtaCent);
  fOutputList->Add(hPt);
  fOutputList->Add(hnMCRecPrimZvPtEtaCent);
  fOutputList->Add(hnMCGenZvPtEtaCent);
  fOutputList->Add(hnMCRecSecZvPtEtaCent);
  fOutputList->Add(hMCPt);
  fOutputList->Add(hEventStatistics);
  fOutputList->Add(hEventStatisticsCentrality);
  fOutputList->Add(hMCEventStatisticsCentrality);
  fOutputList->Add(hAllEventStatisticsCentrality);
  fOutputList->Add(hEventStatisticsCentralityTrigger);
  fOutputList->Add(hnZvMultCent);
  fOutputList->Add(hTriggerStatistics);
  fOutputList->Add(hMCTrackPdgCode);
  fOutputList->Add(hMCTrackStatusCode);
  fOutputList->Add(hCharge);
  fOutputList->Add(hMCCharge);
  fOutputList->Add(hMCPdgPt);
  fOutputList->Add(hMCHijingPrim);
  fOutputList->Add(hDCAPtAll);
  fOutputList->Add(hDCAPtAccepted);
  fOutputList->Add(hMCDCAPtSecondary);
  fOutputList->Add(hMCDCAPtPrimary);
  fOutputList->Add(hnCrossedRowsClustersChiPtEtaPhiAll);
  fOutputList->Add(hnCrossedRowsClustersChiPtEtaPhiAcc);
  
  PostData(1, fOutputList);
}

void AlidNdPtAnalysisPbPbAOD::UserExec(Option_t *option)
{
  
  // Main Loop
  // called for each event
  hEventStatistics->Fill("all events",1);
  
  // set ZERO pointers:
  AliInputEventHandler *inputHandler = NULL;
  AliAODTrack *track = NULL;
  AliAODMCParticle *mcPart = NULL;
  AliAODMCHeader *mcHdr = NULL;
  AliGenHijingEventHeader *genHijingHeader = NULL;
  //AliGenPythiaEventHeader *genPythiaHeader = NULL;
  
  Bool_t bIsEventSelectedMB = kFALSE;
  Bool_t bIsEventSelectedSemi = kFALSE;
  Bool_t bIsEventSelectedCentral = kFALSE;
  Bool_t bIsEventSelected = kFALSE;
  Bool_t bIsPrimary = kFALSE;
  Bool_t bIsHijingParticle = kFALSE;
  Bool_t bMotherIsHijingParticle = kFALSE;
  //Bool_t bIsPythiaParticle = kFALSE;
  Bool_t bEventHasATrack = kFALSE;
  Bool_t bEventHasATrackInRange = kFALSE;
  Int_t nTriggerFired = 0;
  
  
  Double_t dMCTrackZvPtEtaCent[4] = {0};
  Double_t dTrackZvPtEtaCent[4] = {0};
  
  Double_t dDCA[2] = {0};
  
  Double_t dMCEventZv = -100;
  Double_t dEventZv = -100;
  Int_t iAcceptedMultiplicity = 0;
  
  bIsMonteCarlo = kFALSE;
  
  AliAODEvent *eventAOD = 0x0;
  eventAOD = dynamic_cast<AliAODEvent*>( InputEvent() );
  if (!eventAOD) {
    AliWarning("ERROR: eventAOD not available \n");
    return;
  }
  
  // check, which trigger has been fired
  inputHandler = (AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  bIsEventSelectedMB = ( inputHandler->IsEventSelected() & AliVEvent::kMB);
  bIsEventSelectedSemi = ( inputHandler->IsEventSelected() & AliVEvent::kSemiCentral);
  bIsEventSelectedCentral = ( inputHandler->IsEventSelected() & AliVEvent::kCentral);
  
  if(bIsEventSelectedMB || bIsEventSelectedSemi || bIsEventSelectedCentral) hTriggerStatistics->Fill("all triggered events",1);
  if(bIsEventSelectedMB) { hTriggerStatistics->Fill("MB trigger",1); nTriggerFired++; }
  if(bIsEventSelectedSemi) { hTriggerStatistics->Fill("SemiCentral trigger",1); nTriggerFired++; }
  if(bIsEventSelectedCentral) { hTriggerStatistics->Fill("Central trigger",1); nTriggerFired++; }
  if(nTriggerFired == 0) { hTriggerStatistics->Fill("No trigger",1); }
  
  bIsEventSelected = ( inputHandler->IsEventSelected() & GetCollisionCandidates() );
  
  // only take tracks of events, which are triggered
  if(nTriggerFired == 0) { return; } 
  
  //   if( !bIsEventSelected || nTriggerFired>1 ) return;
  
  //   hEventStatistics->Fill("events with only coll. cand.", 1);
  
  
  
  // check if there is a stack, if yes, then do MC loop
  TList *list = eventAOD->GetList();
  TClonesArray *stack = 0x0;
  stack = (TClonesArray*)list->FindObject(AliAODMCParticle::StdBranchName());
  
  if( stack )
  {
    bIsMonteCarlo = kTRUE;
    
    mcHdr = (AliAODMCHeader*)list->FindObject(AliAODMCHeader::StdBranchName());
    
    genHijingHeader = GetHijingEventHeader(mcHdr);
    //     genPythiaHeader = GetPythiaEventHeader(mcHdr);
    
    if(!genHijingHeader) { return; }
    
    //     if(!genPythiaHeader)  { return; }
    
    dMCEventZv = mcHdr->GetVtxZ();
    dMCTrackZvPtEtaCent[0] = dMCEventZv;
    hEventStatistics->Fill("MC all events",1);
  }
  
  AliCentrality* aCentrality = eventAOD->GetCentrality();
  Double_t dCentrality = aCentrality->GetCentralityPercentile("V0M");
  
  if( dCentrality < 0 ) return;
  hEventStatistics->Fill("after centrality selection",1);
  
  
  
  // start with MC truth analysis
  if(bIsMonteCarlo)
  {
    
    if( dMCEventZv > dCutMaxZVertex )  { return; }
    
    dMCTrackZvPtEtaCent[0] = dMCEventZv;
    
    hEventStatistics->Fill("MC afterZv cut",1);
    
    for(Int_t iMCtrack = 0; iMCtrack < stack->GetEntriesFast(); iMCtrack++)
    {
      mcPart =(AliAODMCParticle*)stack->At(iMCtrack);
      
      // check for charge
      if( !(IsMCTrackAccepted(mcPart)) ) continue;
      
      if(!IsHijingParticle(mcPart, genHijingHeader)) { continue; }
      
      if(mcPart->IsPhysicalPrimary() ) 
      {
	hMCHijingPrim->Fill("IsPhysicalPrimary",1);
      }
      else
      {
	hMCHijingPrim->Fill("NOT a primary",1);
	continue;
      }
      
      
      //       
      // ======================== fill histograms ========================
      dMCTrackZvPtEtaCent[1] = mcPart->Pt();
      dMCTrackZvPtEtaCent[2] = mcPart->Eta();
      dMCTrackZvPtEtaCent[3] = dCentrality;
      
      bEventHasATrack = kTRUE;
      
      hnMCGenZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
      
      if( (dMCTrackZvPtEtaCent[1] > dCutPtMin) &&
	(dMCTrackZvPtEtaCent[1] < dCutPtMax) &&
	(dMCTrackZvPtEtaCent[2] > dCutEtaMin) &&
	(dMCTrackZvPtEtaCent[2] < dCutEtaMax) )
      {
	hMCPt->Fill(mcPart->Pt());
	hMCCharge->Fill(mcPart->Charge()/3.);
	bEventHasATrackInRange = kTRUE;
      }
      
    }
  } // isMonteCarlo
  if(bEventHasATrack) { hEventStatistics->Fill("MC events with tracks",1); }
  if(bEventHasATrackInRange) 
  { 
    hEventStatistics->Fill("MC events with tracks in range",1); 
    hMCEventStatisticsCentrality->Fill(dCentrality);
  }
  bEventHasATrack = kFALSE;
  bEventHasATrackInRange = kFALSE;
  
  
  
  // Loop over recontructed tracks
  
  dEventZv = eventAOD->GetPrimaryVertex()->GetZ();
  if( TMath::Abs(dEventZv) > dCutMaxZVertex ) return;
  
  hAllEventStatisticsCentrality->Fill(dCentrality/*, nTriggerFired*/);
  
  hEventStatistics->Fill("after Zv cut",1);
  
  dTrackZvPtEtaCent[0] = dEventZv;
  
  for(Int_t itrack = 0; itrack < eventAOD->GetNumberOfTracks(); itrack++)
  {
    track = eventAOD->GetTrack(itrack);
    if(!track) continue;
    
    mcPart = NULL;
    dMCTrackZvPtEtaCent[1] = 0;
    dMCTrackZvPtEtaCent[2] = 0;
    dMCTrackZvPtEtaCent[3] = 0;
    
    bIsPrimary = kFALSE;
    
    GetDCA(track, eventAOD, dDCA);
    
    Double_t dDCAxyDCAzPt[5] = { dDCA[0], dDCA[1], track->Pt(), track->Eta(), track->Phi() };
    
    hDCAPtAll->Fill(dDCAxyDCAzPt);
    
    if( !(IsTrackAccepted(track)) ) continue;
    
    dTrackZvPtEtaCent[1] = track->Pt();
    dTrackZvPtEtaCent[2] = track->Eta();
    dTrackZvPtEtaCent[3] = dCentrality;
    
    if( bIsMonteCarlo )
    {
      mcPart = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
      if( !mcPart ) { continue; }
      
      // check for charge
      if( !(IsMCTrackAccepted(mcPart)) ) {  continue; } 
      
      bIsHijingParticle = IsHijingParticle(mcPart, genHijingHeader);
      //       bIsPythiaParticle = IsPythiaParticle(mcPart, genPythiaHeader);
      
      //       if(!bIsHijingParticle) continue; // only take real tracks, not injected ones
      
      bIsPrimary = mcPart->IsPhysicalPrimary();
      
      dMCTrackZvPtEtaCent[1] = mcPart->Pt();
      dMCTrackZvPtEtaCent[2] = mcPart->Eta();
      dMCTrackZvPtEtaCent[3] = dCentrality;
      
      if(bIsPrimary && bIsHijingParticle)
      {
	hnMCRecPrimZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
	hMCDCAPtPrimary->Fill(dDCAxyDCAzPt);
      }
      
      if(!bIsPrimary /*&& !bIsHijingParticle*/)
      {
	Int_t indexMoth = mcPart->GetMother(); 
	if(indexMoth >= 0)
	{
	  AliAODMCParticle* moth = (AliAODMCParticle*)stack->At(indexMoth);
	  bMotherIsHijingParticle = IsHijingParticle(moth, genHijingHeader);
	  
	  if(bMotherIsHijingParticle) // only store secondaries, which come from a not embedded signal!
	  {
	    hMCTrackStatusCode->Fill(Form("%d",mcPart->GetStatus()), 1);
	    if(TMath::Abs(mcPart->Eta()) < 0.8) { hMCPdgPt->Fill(mcPart->Pt(), Form("%s",GetParticleName(mcPart->GetPdgCode())), 1); }
	    
	    hnMCRecSecZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
	    hMCDCAPtSecondary->Fill(dDCAxyDCAzPt);
	    hMCTrackPdgCode->Fill(Form("%s_H%i_H%i",GetParticleName(moth->GetPdgCode()),bMotherIsHijingParticle, bIsHijingParticle), 1);
	    // 	  delete moth;
	  }
	} 	
      }
    } // end isMonteCarlo 
    
    // ======================== fill histograms ========================
    
    // only keep prim and sec from not embedded signal
    Bool_t bKeepMCTrack = kFALSE;
    if(bIsMonteCarlo) 
    {
      if( (bIsHijingParticle && bIsPrimary) ^ (bMotherIsHijingParticle && !bIsPrimary) )
      {
	bKeepMCTrack = kTRUE;
      }
      else
      {
	continue;
      }
    }
    
    bEventHasATrack = kTRUE;
    
    hnZvPtEtaCent->Fill(dTrackZvPtEtaCent);
    hDCAPtAccepted->Fill(dDCAxyDCAzPt);
    
    if( (dTrackZvPtEtaCent[1] > dCutPtMin) &&
      (dTrackZvPtEtaCent[1] < dCutPtMax) &&
      (dTrackZvPtEtaCent[2] > dCutEtaMin) &&
      (dTrackZvPtEtaCent[2] < dCutEtaMax) )
    {
      iAcceptedMultiplicity++;
      bEventHasATrackInRange = kTRUE;
      hPt->Fill(track->Pt());
      hCharge->Fill(track->Charge());
    }
  } // end track loop
  
  if(bEventHasATrack) { hEventStatistics->Fill("events with tracks",1); bEventHasATrack = kFALSE; }
  
  if(bEventHasATrackInRange) 
  { 
    hEventStatistics->Fill("events with tracks in range",1); 
    hEventStatisticsCentrality->Fill(dCentrality); 
    bEventHasATrackInRange = kFALSE; 
    
    if(bIsEventSelectedMB) hEventStatisticsCentralityTrigger->Fill(dCentrality, 0);
    if(bIsEventSelectedSemi) hEventStatisticsCentralityTrigger->Fill(dCentrality, 1);
    if(bIsEventSelectedCentral) hEventStatisticsCentralityTrigger->Fill(dCentrality, 2);
  }
  
  Double_t dEventZvMultCent[3] = {dEventZv, iAcceptedMultiplicity, dCentrality};
  hnZvMultCent->Fill(dEventZvMultCent);
  
  PostData(1, fOutputList);
  
}

Bool_t AlidNdPtAnalysisPbPbAOD::IsTrackAccepted(AliAODTrack *tr)
{
  if(!tr) return kFALSE;
  
  if(tr->Charge()==0) { return kFALSE; }
  
  Double_t dNClustersTPC = tr->GetTPCNcls();
  Double_t dCrossedRowsTPC = tr->GetTPCClusterInfo(2,1);
  Double_t dChi2PerClusterTPC = (dNClustersTPC>0)?tr->Chi2perNDF()*(dNClustersTPC-5)/dNClustersTPC:-1.; // see AliDielectronVarManager.h
  
  //   hAllCrossedRowsTPC->Fill(dCrossedRowsTPC);
  
  Double_t dRowClusterChiPtEtaPhi[6] = {dCrossedRowsTPC, dNClustersTPC, dChi2PerClusterTPC, tr->Pt(), tr->Eta(), tr->Phi() };
  hnCrossedRowsClustersChiPtEtaPhiAll->Fill(dRowClusterChiPtEtaPhi);
  
  // filter bit 5
  if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { return kFALSE; } 
  
  // filter bit 4
  //   if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) ) { return kFALSE; }
  
  //   hFilterCrossedRowsTPC->Fill(dCrossedRowsTPC);
  
  if(dCrossedRowsTPC < GetCutMinNCrossedRowsTPC()) { return kFALSE; }
  
  hnCrossedRowsClustersChiPtEtaPhiAcc->Fill(dRowClusterChiPtEtaPhi);
  
  //   hAccNclsTPC->Fill(dNClustersTPC);
  //   hAccCrossedRowsTPC->Fill(dCrossedRowsTPC);
  //   Double_t dFindableClustersTPC = tr->GetTPCNclsF();
  //   Double_t dChi2PerClusterTPC = (dNClustersTPC>0)?tr->Chi2perNDF()*(dNClustersTPC-5)/dNClustersTPC:-1.; // see AliDielectronVarManager.h
  //   
  //   Bool_t bIsFromKink = kFALSE;
  //   if(tr->GetProdVertex()->GetType() == AliAODVertex::kKink) bIsFromKink = kTRUE;
  //   
  //   // from AliAnalysisTaskPIDqa.cxx
  //   ULong_t uStatus = tr->GetStatus();
  //   Bool_t bHasRefitTPC = kFALSE;
  //   Bool_t bHasRefitITS = kFALSE;
  //   
  //   if ((uStatus & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) bHasRefitTPC = kTRUE;
  //   if ((uStatus & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) bHasRefitITS = kTRUE;
  //   
  //   // from AliDielectronVarManager.h
  //   Bool_t bHasHitInSPD = kFALSE;
  //   for (Int_t iC=0; iC<2; iC++) 
  //   {
    //     if (((tr->GetITSClusterMap()) & (1<<(iC))) > 0) {  bHasHitInSPD = kTRUE;  }
    //   }
    //   
    //   Double_t dNClustersITS = tr->GetITSNcls();
    
    // cuts to be done:
    // TPC  
    //   esdTrackCuts->SetMinNCrossedRowsTPC(70);
    //   esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    //   
    //   esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    //   esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    //   esdTrackCuts->SetRequireTPCRefit(kTRUE);
    // ITS
    //   esdTrackCuts->SetRequireITSRefit(kTRUE);
    //   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,	 AliESDtrackCuts::kAny);
    //   
    //   esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
    //   esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
    //   
    //   esdTrackCuts->SetMaxDCAToVertexZ(2);
    //   esdTrackCuts->SetDCAToVertex2D(kFALSE);
    //   esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    //   
    //   esdTrackCuts->SetMaxChi2PerClusterITS(36);
    
    
    return kTRUE;
}

Bool_t AlidNdPtAnalysisPbPbAOD::GetDCA(const AliAODTrack *track, AliAODEvent *evt, Double_t d0z0[2])
{
  // function adapted from AliDielectronVarManager.h
  
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    return kTRUE;
  }
  
  Bool_t ok=kFALSE;
  if(evt) {
    Double_t covd0z0[3];
    //AliAODTrack copy(*track);
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);
    
    Float_t xstart = etp.GetX();
    if(xstart>3.) {
      d0z0[0]=-999.;
      d0z0[1]=-999.;
      //printf("This method can be used only for propagation inside the beam pipe \n");
      return kFALSE;
    }
    
    
    AliAODVertex *vtx =(AliAODVertex*)(evt->GetPrimaryVertex());
    Double_t fBzkG = evt->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
    //ok = copy.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  }
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}


Bool_t AlidNdPtAnalysisPbPbAOD::IsMCTrackAccepted(AliAODMCParticle *part)
{
  if(!part) return kFALSE;
  
  Double_t charge = part->Charge()/3.;
  if (TMath::Abs(charge) < 0.001) return kFALSE;
  
  return kTRUE;
}

const char * AlidNdPtAnalysisPbPbAOD::GetParticleName(Int_t pdg) 
{
  TParticlePDG * p1 = TDatabasePDG::Instance()->GetParticle(pdg);
  if(p1) return p1->GetName();
  return Form("%d", pdg);
}

AliGenHijingEventHeader* AlidNdPtAnalysisPbPbAOD::GetHijingEventHeader(AliAODMCHeader *header)
{
  //
  // inspired by PWGJE/AliPWG4HighPtSpectra.cxx
  //
  
  if(!header) return 0x0;
  AliGenHijingEventHeader* hijingGenHeader = NULL;
  
  TList* headerList = header->GetCocktailHeaders();
  
  for(Int_t i = 0; i < headerList->GetEntries(); i++)
  {
    hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(headerList->At(i));
    if(hijingGenHeader) break;
  }
  
  if(!hijingGenHeader) return 0x0;
  
  return hijingGenHeader;
}

AliGenPythiaEventHeader* AlidNdPtAnalysisPbPbAOD::GetPythiaEventHeader(AliAODMCHeader *header)
{
  //
  // inspired by PWGJE/AliPWG4HighPtSpectra.cxx
  //
  
  if(!header) return 0x0;
  AliGenPythiaEventHeader* PythiaGenHeader = NULL;
  
  TList* headerList = header->GetCocktailHeaders();
  
  for(Int_t i = 0; i < headerList->GetEntries(); i++)
  {
    PythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
    if(PythiaGenHeader) break;
  }
  
  if(!PythiaGenHeader) return 0x0;
  
  return PythiaGenHeader;
}

//________________________________________________________________________
Bool_t AlidNdPtAnalysisPbPbAOD::IsHijingParticle(const AliAODMCParticle *part, AliGenHijingEventHeader* hijingGenHeader){
  
  // Check whether a particle is from Hijing or some injected 
  
  if(part->Label() > (hijingGenHeader->NProduced()-1)) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AlidNdPtAnalysisPbPbAOD::IsPythiaParticle(const AliAODMCParticle *part, AliGenPythiaEventHeader* pythiaGenHeader){
  
  // Check whether a particle is from Pythia or some injected 
  
  if(part->Label() > (pythiaGenHeader->NProduced()-1)) return kFALSE;
  return kTRUE;
}    

Double_t* AlidNdPtAnalysisPbPbAOD::GetArrayClone(Int_t n, Double_t* source)
{
  if (!source || n==0) return 0;
  Double_t* dest = new Double_t[n];
  for (Int_t i=0; i<n ; i++) { dest[i] = source[i]; }
  return dest;
}

void AlidNdPtAnalysisPbPbAOD::Terminate(Option_t *)
{
  
}


