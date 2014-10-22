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
// last modified: 10.06.2014
//------------------------------------------------------------------------------
/*
* This task analysis measured data in PbPb collisions stored in AODs and extract 
* transverse momentum spectra for unidentified charged hadrons vs. centrality.
* Based on MC the efficiency and secondary contamination are determined,
* to correct the measured pT distribution.
* Histograms for the pT resolution correction are also filled.
*
*/ 


#include "AlidNdPtAnalysisPbPbAOD.h"

#include "AliAnalysisTaskSE.h"

using namespace std;

ClassImp(AlidNdPtAnalysisPbPbAOD)

AlidNdPtAnalysisPbPbAOD::AlidNdPtAnalysisPbPbAOD(const char *name) : AliAnalysisTaskSE(name),
fOutputList(0),
// Histograms
fPt(0),
fMCPt(0),
fZvPtEtaCent(0),
fDeltaphiPtEtaCent(0),
fPtResptCent(0),
fMCRecPrimZvPtEtaCent(0),
fMCGenZvPtEtaCent(0),
fMCRecSecZvPtEtaCent(0),
fMCRecPrimDeltaphiPtEtaCent(0),
fMCGenDeltaphiPtEtaCent(0),
fMCRecSecDeltaphiPtEtaCent(0),
fEventStatistics(0),
fEventStatisticsCentrality(0),
fMCEventStatisticsCentrality(0),
fAllEventStatisticsCentrality(0),
fEventStatisticsCentralityTrigger(0),
fZvMultCent(0),
fTriggerStatistics(0),
fCharge(0),
fMCCharge(0),
fDCAPtAll(0),
fDCAPtAccepted(0),
fMCDCAPtSecondary(0),
fMCDCAPtPrimary(0),
fCutPercClusters(0),
fCutPercCrossed(0),
fCrossCheckRowsLength(0),
fCrossCheckClusterLength(0),
fCrossCheckRowsLengthAcc(0),
fCrossCheckClusterLengthAcc(0),
fCutSettings(0),
fEventplaneDist(0),
fEventplaneRunDist(0),
fMCEventplaneDist(0),
fCorrelEventplaneMCDATA(0),
fCorrelEventplaneDefaultCorrected(0),
fEventplaneSubtractedPercentage(0),
// cross check for event plane resolution
fEPDistCent(0),
fPhiCent(0),
fPcosEPCent(0),
fPsinEPCent(0),
fPcosPhiCent(0),
fPsinPhiCent(0),
// cross check for event plane determination
fDeltaPhiCent(0),
//global
fIsMonteCarlo(0),
fEPselector("Q"),
// event cut variables
fCutMaxZVertex(10.),  
// track kinematic cut variables
fCutPtMin(0.15),
fCutPtMax(200.),
fCutEtaMin(-0.8),
fCutEtaMax(0.8),    
// track quality cut variables
fFilterBit(AliAODTrack::kTrkGlobal),
fUseRelativeCuts(kFALSE),
fCutRequireTPCRefit(kTRUE),
fCutRequireITSRefit(kTRUE),
fCutMinNumberOfClusters(60),
fCutPercMinNumberOfClusters(0.2),
fCutMinNumberOfCrossedRows(120.),
fCutPercMinNumberOfCrossedRows(0.2),
fCutMinRatioCrossedRowsOverFindableClustersTPC(0.8),
fCutMaxChi2PerClusterTPC(4.),
fCutMaxFractionSharedTPCClusters(0.4),
fCutMaxDCAToVertexZ(3.0),
fCutMaxDCAToVertexXY(3.0),
fCutMaxChi2PerClusterITS(36.),
fCutDCAToVertex2D(kFALSE),
fCutRequireSigmaToVertex(kFALSE),
fCutMaxDCAToVertexXYPtDepPar0(0.0182),
fCutMaxDCAToVertexXYPtDepPar1(0.0350),
fCutMaxDCAToVertexXYPtDepPar2(1.01),
fCutAcceptKinkDaughters(kFALSE),
fCutMaxChi2TPCConstrainedGlobal(36.),
fCutLengthInTPCPtDependent(kFALSE),
fPrefactorLengthInTPCPtDependent(1),
// binning for THnSparse
fMultNbins(0),
fPtNbins(0),
fPtCorrNbins(0),
fPtCheckNbins(0),
fEtaNbins(0),
fEtaCheckNbins(0),
fZvNbins(0),
fCentralityNbins(0),
fPhiNbins(0),
fRunNumberNbins(0),
fBinsMult(0),
fBinsPt(0),
fBinsPtCorr(0),
fBinsPtCheck(0),
fBinsEta(0),
fBinsEtaCheck(0),
fBinsZv(0),
fBinsCentrality(0),
fBinsPhi(0),
fBinsRunNumber(0)
{
  
  for(Int_t i = 0; i < cqMax; i++)
  {
	fCrossCheckAll[i] = 0;
	fCrossCheckAcc[i] = 0;
  }
  
  fMultNbins = 0;
  fPtNbins = 0;
  fPtCorrNbins = 0;
  fPtCheckNbins = 0;
  fEtaNbins = 0;
  fEtaCheckNbins = 0;
  fZvNbins = 0;
  fCentralityNbins = 0;
  fPhiNbins = 0;
  fRunNumberNbins = 0;
  fBinsMult = 0;
  fBinsPt = 0;
  fBinsPtCorr = 0;
  fBinsPtCheck = 0;
  fBinsEta = 0;
  fBinsEtaCheck = 0;
  fBinsZv = 0;
  fBinsCentrality = 0;
  fBinsPhi = 0;
  fBinsRunNumber = 0;
  
  DefineOutput(1, TList::Class());
}

// destructor
AlidNdPtAnalysisPbPbAOD::~AlidNdPtAnalysisPbPbAOD()
{ 
  //
  //  because task is owner of the output list, all objects are deleted, when list->Clear() is called
  //
  if(fOutputList)
  {
	fOutputList->Clear();
	delete fOutputList;
  }
  fOutputList = 0;
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
  Double_t binsPtCorrDefault[37] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 3.0, 4.0, 200.0}; 
  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZvDefault[7] = {-30.,-10.,-5.,0.,5.,10.,30.};
  Double_t binsCentralityDefault[12] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};  
  
  Double_t binsPhiDefault[37] = { 0., 0.174533, 0.349066, 0.523599, 0.698132, 0.872665, 1.0472, 1.22173, 1.39626, 1.5708, 1.74533, 1.91986, 2.0944, 2.26893, 2.44346, 2.61799, 2.79253, 2.96706, 3.14159, 3.31613, 3.49066, 3.66519, 3.83972, 4.01426, 4.18879, 4.36332, 4.53786, 4.71239, 4.88692, 5.06145, 5.23599, 5.41052, 5.58505, 5.75959, 5.93412, 6.10865, 2.*TMath::Pi()};
  
  Double_t binsPtCheckDefault[20] = {0.,0.15,0.5,1.0,2.0,3.0,4.0, 5.0, 10.0, 13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0};  
  Double_t binsEtaCheckDefault[7] = {-1.0,-0.8,-0.4,0.,0.4,0.8,1.0};
  
  Double_t binsRunNumbers2011[186] = {
	167693, 167706, 167711, 167712, 167713, 167806, 167807, 167808, 167813, 167814, 167818, 167841, 167842, 167844, 167846, 167902, 167903, 167909, 167915, 167920, 167921, 167985, 167986, 167987, 167988, 168066, 168068, 168069, 168076, 168103, 168104, 168105, 168107, 168108, 168115, 168171, 168172, 168173, 168175, 168177, 168181, 168203, 168204, 168205, 168206, 168207, 168208, 168212, 168213, 168310, 168311, 168318, 168322, 168325, 168341, 168342, 168356, 168361, 168362, 168458, 168460, 168461, 168464, 168467, 168511, 168512, 168514, 168644, 168777, 168826, 168984, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169143, 169144, 169145, 169148, 169156, 169160, 169167, 169236, 169238, 169377, 169382, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169584, 169586, 169587, 169588, 169590, 169591, 169628, 169683, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169914, 169918, 169919, 169920, 169922, 169923, 169924, 169926, 169956, 169961, 169965, 169969, 169975, 169981, 170027, 170036, 170038, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170152, 170155, 170159, 170162, 170163, 170193, 170195, 170203, 170204, 170205, 170207, 170208, 170228, 170230, 170264, 170267, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170374, 170387, 170388, 170389, 170390, 170546, 170552, 170556, 170572, 170593, 170593+1 
  };
  
  // if no binning is set, use the default
  if (!fBinsMult)	{ SetBinsMult(48,binsMultDefault); }
  if (!fBinsPt)		{ SetBinsPt(82,binsPtDefault); }
  if (!fBinsPtCorr)	{ SetBinsPtCorr(37,binsPtCorrDefault); }
  if (!fBinsPtCheck)	{ SetBinsPtCheck(20,binsPtCheckDefault); }
  if (!fBinsEta)	{ SetBinsEta(31,binsEtaDefault); }
  if (!fBinsEtaCheck)	{ SetBinsEtaCheck(7,binsEtaCheckDefault); }
  if (!fBinsZv)		{ SetBinsZv(7,binsZvDefault); }  
  if (!fBinsCentrality)	{ SetBinsCentrality(12,binsCentralityDefault); }
  if (!fBinsPhi)	{ SetBinsPhi(37,binsPhiDefault); }
  if (!fBinsRunNumber) {SetBinsRunNumber(186, binsRunNumbers2011); }
  
  Int_t binsZvPtEtaCent[4]={fZvNbins-1,fPtNbins-1,fEtaNbins-1,fCentralityNbins-1};
  Int_t binsPhiPtEtaCent[4]={fPhiNbins-1,fPtNbins-1,fEtaNbins-1,fCentralityNbins-1};
  Int_t binsZvMultCent[3]={fZvNbins-1,fMultNbins-1,fCentralityNbins-1};
  
  Int_t binsOneOverPtPtResCent[3]={400,300,11};
  Double_t minbinsOneOverPtPtResCent[3]={0,0,0}; 
  Double_t maxbinsOneOverPtPtResCent[3]={1,0.015,100};
  
  // define Histograms
  fZvPtEtaCent = new THnSparseF("fZvPtEtaCent","Zv:Pt:Eta:Centrality",4,binsZvPtEtaCent);
  fZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fZvPtEtaCent->SetBinEdges(1,fBinsPt);
  fZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fZvPtEtaCent->GetAxis(0)->SetTitle("Zv (cm)");
  fZvPtEtaCent->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fZvPtEtaCent->GetAxis(2)->SetTitle("Eta");
  fZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fZvPtEtaCent->Sumw2();
  
  fDeltaphiPtEtaCent = new THnSparseF("fDeltaphiPtEtaCent","Deltaphi:Pt:Eta:Centrality",4,binsPhiPtEtaCent);
  fDeltaphiPtEtaCent->SetBinEdges(0,fBinsPhi);
  fDeltaphiPtEtaCent->SetBinEdges(1,fBinsPt);
  fDeltaphiPtEtaCent->SetBinEdges(2,fBinsEta);
  fDeltaphiPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fDeltaphiPtEtaCent->GetAxis(0)->SetTitle("#Delta phi to ep");
  fDeltaphiPtEtaCent->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fDeltaphiPtEtaCent->GetAxis(2)->SetTitle("Eta");
  fDeltaphiPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fDeltaphiPtEtaCent->Sumw2();
  
  fPtResptCent = new THnSparseF("fPtResptCent","OneOverPt:PtRes:Centrality",3,binsOneOverPtPtResCent, minbinsOneOverPtPtResCent, maxbinsOneOverPtPtResCent);
  fPtResptCent->SetBinEdges(2, fBinsCentrality);
  fPtResptCent->GetAxis(0)->SetTitle("1/pT (GeV/c)^{-1}");
  fPtResptCent->GetAxis(1)->SetTitle("#sigma(1/pT)");
  fPtResptCent->GetAxis(2)->SetTitle("centrality");
  fPtResptCent->Sumw2();
  
  fMCRecPrimZvPtEtaCent = new THnSparseF("fMCRecPrimZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtEtaCent);
  fMCRecPrimZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fMCRecPrimZvPtEtaCent->SetBinEdges(1,fBinsPt);
  fMCRecPrimZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCRecPrimZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCRecPrimZvPtEtaCent->GetAxis(0)->SetTitle("MC Zv (cm)");
  fMCRecPrimZvPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  fMCRecPrimZvPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  fMCRecPrimZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecPrimZvPtEtaCent->Sumw2();
  
  fMCGenZvPtEtaCent = new THnSparseF("fMCGenZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtEtaCent);
  fMCGenZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fMCGenZvPtEtaCent->SetBinEdges(1,fBinsPt);
  fMCGenZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCGenZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCGenZvPtEtaCent->GetAxis(0)->SetTitle("MC Zv (cm)");
  fMCGenZvPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  fMCGenZvPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  fMCGenZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCGenZvPtEtaCent->Sumw2();
  
  fMCRecSecZvPtEtaCent = new THnSparseF("fMCRecSecZvPtEtaCent","mcZv:mcPt:mcEta:Centrality",4,binsZvPtEtaCent);
  fMCRecSecZvPtEtaCent->SetBinEdges(0,fBinsZv);
  fMCRecSecZvPtEtaCent->SetBinEdges(1,fBinsPt);
  fMCRecSecZvPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCRecSecZvPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCRecSecZvPtEtaCent->GetAxis(0)->SetTitle("MC Sec Zv (cm)");
  fMCRecSecZvPtEtaCent->GetAxis(1)->SetTitle("MC Sec Pt (GeV/c)");
  fMCRecSecZvPtEtaCent->GetAxis(2)->SetTitle("MC Sec Eta");
  fMCRecSecZvPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecSecZvPtEtaCent->Sumw2();
  
  fMCRecPrimDeltaphiPtEtaCent = new THnSparseF("fMCRecPrimDeltaphiPtEtaCent","mcDeltaphi:mcPt:mcEta:Centrality",4,binsPhiPtEtaCent);
  fMCRecPrimDeltaphiPtEtaCent->SetBinEdges(0,fBinsPhi);
  fMCRecPrimDeltaphiPtEtaCent->SetBinEdges(1,fBinsPt);
  fMCRecPrimDeltaphiPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCRecPrimDeltaphiPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCRecPrimDeltaphiPtEtaCent->GetAxis(0)->SetTitle("MC #Delta phi to rp");
  fMCRecPrimDeltaphiPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  fMCRecPrimDeltaphiPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  fMCRecPrimDeltaphiPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecPrimDeltaphiPtEtaCent->Sumw2();
  
  fMCGenDeltaphiPtEtaCent = new THnSparseF("fMCGenDeltaphiPtEtaCent","mcDeltaphi:mcPt:mcEta:Centrality",4,binsPhiPtEtaCent);
  fMCGenDeltaphiPtEtaCent->SetBinEdges(0,fBinsPhi);
  fMCGenDeltaphiPtEtaCent->SetBinEdges(1,fBinsPt);
  fMCGenDeltaphiPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCGenDeltaphiPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCGenDeltaphiPtEtaCent->GetAxis(0)->SetTitle("MC #Delta phi to rp");
  fMCGenDeltaphiPtEtaCent->GetAxis(1)->SetTitle("MC Pt (GeV/c)");
  fMCGenDeltaphiPtEtaCent->GetAxis(2)->SetTitle("MC Eta");
  fMCGenDeltaphiPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCGenDeltaphiPtEtaCent->Sumw2();
  
  fMCRecSecDeltaphiPtEtaCent = new THnSparseF("fMCRecSecDeltaphiPtEtaCent","mcDeltaphi:mcPt:mcEta:Centrality",4,binsPhiPtEtaCent);
  fMCRecSecDeltaphiPtEtaCent->SetBinEdges(0,fBinsPhi);
  fMCRecSecDeltaphiPtEtaCent->SetBinEdges(1,fBinsPt);
  fMCRecSecDeltaphiPtEtaCent->SetBinEdges(2,fBinsEta);
  fMCRecSecDeltaphiPtEtaCent->SetBinEdges(3,fBinsCentrality);
  fMCRecSecDeltaphiPtEtaCent->GetAxis(0)->SetTitle("MC Sec #Delta phi to rp");
  fMCRecSecDeltaphiPtEtaCent->GetAxis(1)->SetTitle("MC Sec Pt (GeV/c)");
  fMCRecSecDeltaphiPtEtaCent->GetAxis(2)->SetTitle("MC Sec Eta");
  fMCRecSecDeltaphiPtEtaCent->GetAxis(3)->SetTitle("Centrality");
  fMCRecSecDeltaphiPtEtaCent->Sumw2();
  
  fPt = new TH1F("fPt","fPt",2000,0,200);
  fPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fPt->GetYaxis()->SetTitle("dN/dp_{T}");
  fPt->Sumw2();
  
  fMCPt = new TH1F("fMCPt","fMCPt",2000,0,200);
  fMCPt->GetXaxis()->SetTitle("MC p_{T} (GeV/c)");
  fMCPt->GetYaxis()->SetTitle("dN/dp_{T}");
  fMCPt->Sumw2();
  
  fEventStatistics = new TH1F("fEventStatistics","fEventStatistics",10,0,10);
  fEventStatistics->GetYaxis()->SetTitle("number of events");
  fEventStatistics->SetBit(TH1::kCanRebin);
  
  fEventStatisticsCentrality = new TH1F("fEventStatisticsCentrality","fEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  fEventStatisticsCentrality->GetYaxis()->SetTitle("number of events");
  
  fMCEventStatisticsCentrality = new TH1F("fMCEventStatisticsCentrality","fMCEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  fMCEventStatisticsCentrality->GetYaxis()->SetTitle("number of MC events");
  
  fAllEventStatisticsCentrality = new TH1F("fAllEventStatisticsCentrality","fAllEventStatisticsCentrality",fCentralityNbins-1, fBinsCentrality);
  fAllEventStatisticsCentrality->GetYaxis()->SetTitle("number of events");
  
  fEventStatisticsCentralityTrigger = new TH2F("fEventStatisticsCentralityTrigger","fEventStatisticsCentralityTrigger;centrality;trigger",100,0,100,3,0,3);
  fEventStatisticsCentralityTrigger->Sumw2();
  
  fZvMultCent = new THnSparseF("fZvMultCent","Zv:mult:Centrality",3,binsZvMultCent);
  fZvMultCent->SetBinEdges(0,fBinsZv);
  fZvMultCent->SetBinEdges(1,fBinsMult);
  fZvMultCent->SetBinEdges(2,fBinsCentrality);
  fZvMultCent->GetAxis(0)->SetTitle("Zv (cm)");
  fZvMultCent->GetAxis(1)->SetTitle("N_{acc}");
  fZvMultCent->GetAxis(2)->SetTitle("Centrality");
  fZvMultCent->Sumw2();
  
  fTriggerStatistics = new TH1F("fTriggerStatistics","fTriggerStatistics",10,0,10);
  fTriggerStatistics->GetYaxis()->SetTitle("number of events");
  
  fCharge = new TH1F("fCharge","fCharge",30, -5, 5);
  fCharge->GetXaxis()->SetTitle("Charge");
  fCharge->GetYaxis()->SetTitle("number of tracks");
  
  fMCCharge = new TH1F("fMCCharge","fMCCharge",30, -5, 5);
  fMCCharge->GetXaxis()->SetTitle("MC Charge");
  fMCCharge->GetYaxis()->SetTitle("number of tracks");  
  
  Int_t binsDCAxyDCAzPtEtaPhi[6] =   { 10 , 10 , fPtCheckNbins-1, fEtaCheckNbins-1,             18, fCentralityNbins-1 };
  Double_t minDCAxyDCAzPtEtaPhi[6] = { -5 , -5 ,               0,             -1.5,             0.,                  0 };
  Double_t maxDCAxyDCAzPtEtaPhi[6] = {  5.,  5.,             100,              1.5, 2.*TMath::Pi(),                100 };
  
  fDCAPtAll = new THnSparseF("fDCAPtAll","fDCAPtAll",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  fDCAPtAccepted = new THnSparseF("fDCAPtAccepted","fDCAPtAccepted",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  fMCDCAPtSecondary = new THnSparseF("fMCDCAPtSecondary","fMCDCAPtSecondary",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  fMCDCAPtPrimary = new THnSparseF("fMCDCAPtPrimary","fMCDCAPtPrimary",6, binsDCAxyDCAzPtEtaPhi, minDCAxyDCAzPtEtaPhi, maxDCAxyDCAzPtEtaPhi);
  
  fDCAPtAll->SetBinEdges(2, fBinsPtCheck);
  fDCAPtAccepted->SetBinEdges(2, fBinsPtCheck);
  fMCDCAPtSecondary->SetBinEdges(2, fBinsPtCheck);
  fMCDCAPtPrimary->SetBinEdges(2, fBinsPtCheck);
  
  fDCAPtAll->SetBinEdges(3, fBinsEtaCheck);
  fDCAPtAccepted->SetBinEdges(3, fBinsEtaCheck);
  fMCDCAPtSecondary->SetBinEdges(3, fBinsEtaCheck);
  fMCDCAPtPrimary->SetBinEdges(3, fBinsEtaCheck);
  
  fDCAPtAll->SetBinEdges(5, fBinsCentrality);
  fDCAPtAccepted->SetBinEdges(5, fBinsCentrality);
  fMCDCAPtSecondary->SetBinEdges(5, fBinsCentrality);
  fMCDCAPtPrimary->SetBinEdges(5, fBinsCentrality);
  
  fDCAPtAll->Sumw2();
  fDCAPtAccepted->Sumw2();
  fMCDCAPtSecondary->Sumw2();
  fMCDCAPtPrimary->Sumw2();
  
  fDCAPtAll->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fDCAPtAll->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fDCAPtAll->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fDCAPtAll->GetAxis(3)->SetTitle("#eta");
  fDCAPtAll->GetAxis(4)->SetTitle("#phi");
  fDCAPtAll->GetAxis(5)->SetTitle("Centrality");
  
  fDCAPtAccepted->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fDCAPtAccepted->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fDCAPtAccepted->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fDCAPtAccepted->GetAxis(3)->SetTitle("#eta");
  fDCAPtAccepted->GetAxis(4)->SetTitle("#phi");
  fDCAPtAccepted->GetAxis(5)->SetTitle("Centrality");
  
  fMCDCAPtSecondary->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fMCDCAPtSecondary->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fMCDCAPtSecondary->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fMCDCAPtSecondary->GetAxis(3)->SetTitle("#eta");
  fMCDCAPtSecondary->GetAxis(4)->SetTitle("#phi");
  fMCDCAPtSecondary->GetAxis(5)->SetTitle("Centrality");
  
  fMCDCAPtPrimary->GetAxis(0)->SetTitle("DCA_{xy} (cm)");
  fMCDCAPtPrimary->GetAxis(1)->SetTitle("DCA_{z} (cm)");
  fMCDCAPtPrimary->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fMCDCAPtPrimary->GetAxis(3)->SetTitle("#eta");
  fMCDCAPtPrimary->GetAxis(4)->SetTitle("#phi");
  fMCDCAPtPrimary->GetAxis(5)->SetTitle("Centrality");
  
  
  char cFullTempTitle[255];
  char cTempTitleAxis0All[255];
  char cTempTitleAxis0Acc[255];
  //   char cTempTitleAxis1[255];
  char cFullTempName[255];
  char cTempNameAxis0[255];
  //   char cTempNameAxis1[255];
  const Int_t iNbinRowsClusters = 21;
  //   Double_t dBinsRowsClusters[iNbinRowsClusters] = {0, 7.95, 15.9, 23.85, 31.8, 39.75, 47.7, 55.65, 63.6, 71.55, 79.5, 87.45, 95.4, 103.35, 111.3, 119.25, 127.2, 135.15, 143.1, 151.05, 159.};
  
  const Int_t iNbinChi = 51;
  const Int_t iNbinLength = 165;
  const Int_t iNbinRowsOverClusters = 60;
  //   Double_t dBinsChi[iNbinChi] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.2, 5.4, 5.6, 5.8, 6, 6.2, 6.4, 6.6, 6.8, 7, 7.2, 7.4, 7.6, 7.8, 8, 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6, 9.8,10.};
  
  Int_t iNbin = 0;
  //   Double_t *dBins = 0x0;
  Double_t dBinMin = 0;
  Double_t dBinMax = 0;
  
  for(Int_t iCheckQuant = 0; iCheckQuant < cqMax; iCheckQuant++)
  {
	// iCheckQuant: 0 = CrossedRows, 1 = Nclusters, 2 = Chi^2/clusterTPC
	if(iCheckQuant == cqCrossedRows) 
	{
	  snprintf(cTempTitleAxis0All,255, "NcrossedRows before Cut"); 
	  snprintf(cTempTitleAxis0Acc,255, "NcrossedRows after Cut"); 
	  snprintf(cTempNameAxis0,255, "CrossedRows");
	  iNbin = iNbinRowsClusters;
	  dBinMin = 0;
	  dBinMax = 159.;
	}
	else if(iCheckQuant == cqNcluster) 
	{
	  snprintf(cTempTitleAxis0All,255, "Nclusters before Cut"); 
	  snprintf(cTempTitleAxis0Acc,255, "Nclusters after Cut"); 
	  snprintf(cTempNameAxis0,255, "Clusters");
	  iNbin = iNbinRowsClusters;
	  dBinMin = 0;
	  dBinMax = 159.;
	}
	else if(iCheckQuant == cqChi) 
	{
	  snprintf(cTempTitleAxis0All,255, "#Chi^{2}/cluster before Cut"); 
	  snprintf(cTempTitleAxis0Acc,255, "#Chi^{2}/cluster after Cut"); 
	  snprintf(cTempNameAxis0,255, "Chi");
	  iNbin = iNbinChi;
	  dBinMin = 0;
	  dBinMax = 10.;
	}
	else if(iCheckQuant == cqLength) 
	{
	  snprintf(cTempTitleAxis0All,255, "Length in TPC before Cut (cm)"); 
	  snprintf(cTempTitleAxis0Acc,255, "Length in TPC after Cut (cm)"); 
	  snprintf(cTempNameAxis0,255, "Length");
	  iNbin = iNbinLength;
	  dBinMin = 0;
	  dBinMax = 165.;
	}
	else if(iCheckQuant == cqRowsOverFindable) 
	{
	  snprintf(cTempTitleAxis0All,255, "Number of Crossed Rows / Number of Findable Clusters before Cut"); 
	  snprintf(cTempTitleAxis0Acc,255, "Number of Crossed Rows / Number of Findable Clusters before Cut"); 
	  snprintf(cTempNameAxis0,255, "RowsOverFindable");
	  iNbin = iNbinRowsOverClusters;
	  dBinMin = 0.6;
	  dBinMax = 1.2;
	}
	
	
	Int_t binsCheckPtEtaPhi[5] = { iNbin, fPtCheckNbins-1, fEtaCheckNbins-1, 18, fCentralityNbins-1};
	//     Int_t binsCheckPtEtaPhi[5] = { iNbin, fPtNbins-1, fEtaCheckNbins-1, 18, fCentralityNbins-1};
	Double_t minCheckPtEtaPhi[5] = { dBinMin,  0, -1.5, 0., 0, };
	Double_t maxCheckPtEtaPhi[5] = { dBinMax, 100, 1.5, 2.*TMath::Pi(), 100};
	
	snprintf(cFullTempName, 255, "f%sPtEtaPhiAll",cTempNameAxis0);
	snprintf(cFullTempTitle, 255,"%s;%s;p_{T} (GeV/c);#eta;#phi;Centrality", cFullTempName, cTempTitleAxis0All);
	fCrossCheckAll[iCheckQuant] = new THnF(cFullTempName, cFullTempTitle, 5, binsCheckPtEtaPhi, minCheckPtEtaPhi, maxCheckPtEtaPhi);
	fCrossCheckAll[iCheckQuant]->SetBinEdges(1, fBinsPtCheck);
	fCrossCheckAll[iCheckQuant]->SetBinEdges(2, fBinsEtaCheck);
	fCrossCheckAll[iCheckQuant]->Sumw2();
	
	snprintf(cFullTempName, 255, "f%sPtEtaPhiAcc",cTempNameAxis0);
	snprintf(cFullTempTitle, 255,"%s;%s;p_{T} (GeV/c);#eta;#phi;Centrality", cFullTempName, cTempTitleAxis0Acc);
	fCrossCheckAcc[iCheckQuant] = new THnF(cFullTempName, cFullTempTitle, 5, binsCheckPtEtaPhi, minCheckPtEtaPhi, maxCheckPtEtaPhi);
	fCrossCheckAcc[iCheckQuant]->SetBinEdges(1, fBinsPtCheck);
	fCrossCheckAcc[iCheckQuant]->SetBinEdges(2, fBinsEtaCheck);
	fCrossCheckAcc[iCheckQuant]->Sumw2();
  } // end iCheckQuant
  
  fCutPercClusters = new TH1F("fCutPercClusters","fCutPercClusters;NclustersTPC;counts",160,0,160);
  fCutPercClusters->Sumw2();
  fCutPercCrossed = new TH1F("fCutPercCrossed","fCutPercCrossed;NcrossedRowsTPC;counts",160,0,160);
  fCutPercCrossed->Sumw2();
  
  fCrossCheckRowsLength = new TH2F("fCrossCheckRowsLength","fCrossCheckRowsLength;Length in TPC;NcrossedRows",170,0,170,170,0,170);
  fCrossCheckRowsLength->Sumw2();
  
  fCrossCheckClusterLength = new TH2F("fCrossCheckClusterLength","fCrossCheckClusterLength;Length in TPC;NClusters",170,0,170,170,0,170);
  fCrossCheckClusterLength->Sumw2();
  
  fCrossCheckRowsLengthAcc = new TH2F("fCrossCheckRowsLengthAcc","fCrossCheckRowsLengthAcc;Length in TPC;NcrossedRows",170,0,170,170,0,170);
  fCrossCheckRowsLengthAcc->Sumw2();
  
  fCrossCheckClusterLengthAcc = new TH2F("fCrossCheckClusterLengthAcc","fCrossCheckClusterLengthAcc;Length in TPC;NClusters",170,0,170,170,0,170);
  fCrossCheckClusterLengthAcc->Sumw2();
  
  fCutSettings = new TH1F("fCutSettings","fCutSettings",100,0,10);
  fCutSettings->GetYaxis()->SetTitle("cut value");
  fCutSettings->SetBit(TH1::kCanRebin);
  
  fEventplaneDist = new TH1F("fEventplaneDist","fEventplaneDist",200, 0, 2.*TMath::Pi());
  fEventplaneDist->GetXaxis()->SetTitle("#phi (event plane)");
  fEventplaneDist->Sumw2();
  
  fEventplaneRunDist = new TH2F("fEventplaneRunDist","fEventplaneRunDist",200, 0, 2.*TMath::Pi(),fRunNumberNbins-1, fBinsRunNumber );
  fEventplaneRunDist->GetXaxis()->SetTitle("#phi (event plane)");
  fEventplaneRunDist->GetYaxis()->SetTitle("runnumber");
  fEventplaneRunDist->Sumw2();
  
  fMCEventplaneDist = new TH1F("fMCEventplaneDist","fMCEventplaneDist",20, -1.*TMath::Pi(), TMath::Pi());
  fMCEventplaneDist->GetXaxis()->SetTitle("#phi (MC event plane)");
  fMCEventplaneDist->Sumw2();
  
  fCorrelEventplaneMCDATA = new TH2F("fCorrelEventplaneMCDATA","fCorrelEventplaneMCDATA",40, -2.*TMath::Pi(), 2.*TMath::Pi(), 40, -2.*TMath::Pi(), 2.*TMath::Pi());
  fCorrelEventplaneMCDATA->GetXaxis()->SetTitle("#phi (event plane)");
  fCorrelEventplaneMCDATA->GetYaxis()->SetTitle("#phi (MC event plane)");
  fCorrelEventplaneMCDATA->Sumw2();
  
  Int_t binsCorrelPhiPhiCent[3] = { 40, 40, 10};
  Double_t minCorrelPhiPhiCent[3] = { -2.*TMath::Pi(), -2.*TMath::Pi(), 0};
  Double_t maxCorrelPhiPhiCent[3] = { 2.*TMath::Pi(), 2.*TMath::Pi(), 100};
  
  fCorrelEventplaneDefaultCorrected = new THnSparseF("fCorrelEventplaneDefaultCorrected","fCorrelEventplaneDefaultCorrected",3,binsCorrelPhiPhiCent, minCorrelPhiPhiCent, maxCorrelPhiPhiCent);
  fCorrelEventplaneDefaultCorrected->SetBinEdges(2, fBinsCentrality);
  fCorrelEventplaneDefaultCorrected->GetAxis(0)->SetTitle("#phi (event plane)");
  fCorrelEventplaneDefaultCorrected->GetAxis(1)->SetTitle("#phi (corrected event plane)");
  fCorrelEventplaneDefaultCorrected->GetAxis(2)->SetTitle("centrality");
  fCorrelEventplaneDefaultCorrected->Sumw2();
  
  fEventplaneSubtractedPercentage = new TH2F("fEventplaneSubtractedPercentage","fEventplaneSubtractedPercentage",100, 0,1, fCentralityNbins-1, fBinsCentrality);
  fEventplaneSubtractedPercentage->GetXaxis()->SetTitle("percentage of tracks, which have been subtracted during analysis");
  fEventplaneSubtractedPercentage->GetYaxis()->SetTitle("centrality");
  fEventplaneSubtractedPercentage->Sumw2();
  
  // cross check for event plane resolution
  fEPDistCent = new TH2F("fEPDistCent","fEPDistCent",20, -2.*TMath::Pi(), 2.*TMath::Pi(), fCentralityNbins-1, fBinsCentrality);
  fEPDistCent->GetXaxis()->SetTitle("#phi (#Psi_{EP})");
  fEPDistCent->GetYaxis()->SetTitle("Centrality");
  fEPDistCent->Sumw2();
  
  fPhiCent = new TH2F("fPhiCent","fPhiCent",200, -2.*TMath::Pi(), 2.*TMath::Pi(), fCentralityNbins-1, fBinsCentrality);
  fPhiCent->GetXaxis()->SetTitle("#phi");
  fPhiCent->GetYaxis()->SetTitle("Centrality");
  fPhiCent->Sumw2();
  
  fPcosEPCent = new TProfile("fPcosEPCent","fPcosEPCent", 100,0,100);
  fPcosEPCent->GetXaxis()->SetTitle("Centrality");
  fPcosEPCent->GetYaxis()->SetTitle("#LT cos 2 #Psi_{EP} #GT");
  fPcosEPCent->Sumw2();
  
  fPsinEPCent = new TProfile("fPsinEPCent","fPsinEPCent", 100,0,100);
  fPsinEPCent->GetXaxis()->SetTitle("Centrality");
  fPsinEPCent->GetYaxis()->SetTitle("#LT sin 2 #Psi_{EP} #GT");
  fPsinEPCent->Sumw2();
  
  fPcosPhiCent = new TProfile("fPcosPhiCent","fPcosPhiCent", 100,0,100);
  fPcosPhiCent->GetXaxis()->SetTitle("Centrality");
  fPcosPhiCent->GetYaxis()->SetTitle("#LT cos 2 #phi #GT");
  fPcosPhiCent->Sumw2();
  
  fPsinPhiCent = new TProfile("fPsinPhiCent","fPsinPhiCent", 100,0,100);
  fPsinPhiCent->GetXaxis()->SetTitle("Centrality");
  fPsinPhiCent->GetYaxis()->SetTitle("#LT sin 2 #phi #GT");
  fPsinPhiCent->Sumw2();
  
  fDeltaPhiCent = new TH2F("fDeltaPhiCent","fDeltaPhiCent",200, -2.*TMath::Pi(), 2.*TMath::Pi(), fCentralityNbins-1, fBinsCentrality);
  fDeltaPhiCent->GetXaxis()->SetTitle("#Delta #phi");
  fDeltaPhiCent->GetYaxis()->SetTitle("Centrality");
  fDeltaPhiCent->Sumw2();
  
  // Add Histos, Profiles etc to List
  fOutputList->Add(fZvPtEtaCent);
  fOutputList->Add(fDeltaphiPtEtaCent);
  fOutputList->Add(fPtResptCent);
  fOutputList->Add(fPt);
  fOutputList->Add(fMCRecPrimZvPtEtaCent);
  fOutputList->Add(fMCGenZvPtEtaCent);
  fOutputList->Add(fMCRecSecZvPtEtaCent);
  fOutputList->Add(fMCRecPrimDeltaphiPtEtaCent);
  fOutputList->Add(fMCGenDeltaphiPtEtaCent);
  fOutputList->Add(fMCRecSecDeltaphiPtEtaCent);
  fOutputList->Add(fMCPt);
  fOutputList->Add(fEventStatistics);
  fOutputList->Add(fEventStatisticsCentrality);
  fOutputList->Add(fMCEventStatisticsCentrality);
  fOutputList->Add(fAllEventStatisticsCentrality);
  fOutputList->Add(fEventStatisticsCentralityTrigger);
  fOutputList->Add(fZvMultCent);
  fOutputList->Add(fTriggerStatistics);
  fOutputList->Add(fCharge);
  fOutputList->Add(fMCCharge);
  fOutputList->Add(fDCAPtAll);
  fOutputList->Add(fDCAPtAccepted);
  fOutputList->Add(fMCDCAPtSecondary);
  fOutputList->Add(fMCDCAPtPrimary);
  for(Int_t i = 0; i < cqMax; i++)
  {   
	fOutputList->Add(fCrossCheckAll[i]);
	fOutputList->Add(fCrossCheckAcc[i]);
  }
  fOutputList->Add(fCutPercClusters);
  fOutputList->Add(fCutPercCrossed);
  fOutputList->Add(fCrossCheckRowsLength);
  fOutputList->Add(fCrossCheckClusterLength);
  fOutputList->Add(fCrossCheckRowsLengthAcc);
  fOutputList->Add(fCrossCheckClusterLengthAcc);
  fOutputList->Add(fCutSettings);
  fOutputList->Add(fEventplaneDist);
  fOutputList->Add(fEventplaneRunDist);
  fOutputList->Add(fMCEventplaneDist);
  fOutputList->Add(fCorrelEventplaneMCDATA);
  fOutputList->Add(fCorrelEventplaneDefaultCorrected);
  fOutputList->Add(fEventplaneSubtractedPercentage);
  
  fOutputList->Add(fEPDistCent);
  fOutputList->Add(fPhiCent);
  fOutputList->Add(fPcosEPCent);
  fOutputList->Add(fPsinEPCent);
  fOutputList->Add(fPcosPhiCent);
  fOutputList->Add(fPsinPhiCent);
  
  fOutputList->Add(fDeltaPhiCent);
    
  StoreCutSettingsToHistogram();
  
  PostData(1, fOutputList);
}

void AlidNdPtAnalysisPbPbAOD::UserExec(Option_t *option)
{
  //
  // Main Loop
  // called for each event
  //
  
  fEventStatistics->Fill("all events",1);
  
  // set ZERO pointers:
  AliInputEventHandler *inputHandler = NULL;
  AliAODTrack *track = NULL;
  AliAODMCParticle *mcPart = NULL;
  AliAODMCHeader *mcHdr = NULL;
  AliGenHijingEventHeader *genHijingHeader = NULL;
  //AliGenPythiaEventHeader *genPythiaHeader = NULL;
  AliEventplane *ep = NULL;
  
  TVector2 *epQvector = NULL;
  
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
  
  Double_t dMCTrackPhiPtEtaCent[4] = {0};
  Double_t dTrackPhiPtEtaCent[4] = {0};
  
  Double_t dDCA[2] = {0};
  
  Double_t dMCEventZv = -100;
  Double_t dEventZv = -100;
  Int_t iAcceptedMultiplicity = 0;
  Double_t dEventplaneAngle = -10;
  Double_t dEventplaneAngleCorrected = -10; // event plane angle, where tracks contributing to this angle have been subtracted
  Double_t dMCEventplaneAngle = -10;
  
  fIsMonteCarlo = kFALSE;
  
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
  
  if(bIsEventSelectedMB || bIsEventSelectedSemi || bIsEventSelectedCentral) fTriggerStatistics->Fill("all triggered events",1);
  if(bIsEventSelectedMB) { fTriggerStatistics->Fill("MB trigger",1); nTriggerFired++; }
  if(bIsEventSelectedSemi) { fTriggerStatistics->Fill("SemiCentral trigger",1); nTriggerFired++; }
  if(bIsEventSelectedCentral) { fTriggerStatistics->Fill("Central trigger",1); nTriggerFired++; }
  if(nTriggerFired == 0) { fTriggerStatistics->Fill("No trigger",1); }
  
  bIsEventSelected = ( inputHandler->IsEventSelected() & GetCollisionCandidates() );
  
  // only take tracks of events, which are triggered
  if(nTriggerFired == 0) { return; } 
  
  
  //   if( !bIsEventSelected || nTriggerFired>1 ) return;
  
  //   fEventStatistics->Fill("events with only coll. cand.", 1);
  
  
  
  // check if there is a stack, if yes, then do MC loop
  TList *list = eventAOD->GetList();
  TClonesArray *stack = 0x0;
  stack = (TClonesArray*)list->FindObject(AliAODMCParticle::StdBranchName());
  
  if( stack )
  {
	fIsMonteCarlo = kTRUE;
	
	mcHdr = (AliAODMCHeader*)list->FindObject(AliAODMCHeader::StdBranchName());
	
	genHijingHeader = GetHijingEventHeader(mcHdr);
	//     genPythiaHeader = GetPythiaEventHeader(mcHdr);
	
	if(!genHijingHeader) { return; }
	
	//     if(!genPythiaHeader)  { return; }
	
	
	dMCEventZv = mcHdr->GetVtxZ();
	dMCTrackZvPtEtaCent[0] = dMCEventZv;
	dMCEventplaneAngle = MoveEventplane(genHijingHeader->ReactionPlaneAngle());
	fEventStatistics->Fill("MC all events",1);
	fMCEventplaneDist->Fill(dMCEventplaneAngle);
  }
  
  AliCentrality* aCentrality = eventAOD->GetCentrality();
  Double_t dCentrality = aCentrality->GetCentralityPercentile("V0M");
  
  if( dCentrality < 0 ) return;
  
  // protection for bias on pt spectra if all triggers selected
  //   if( (bIsEventSelectedCentral) /*&& (!bIsEventSelectedSemi) && (!bIsEventSelectedMB)*/ && (dCentrality > 10) ) return;
  //   if( /*(!bIsEventSelectedCentral) &&*/ (bIsEventSelectedSemi) /*&& (!bIsEventSelectedMB)*/ && (dCentrality < 20) && (dCentrality > 50)) return;
  if( (bIsEventSelectedCentral)  && (dCentrality > 10) ) return;
  if( (bIsEventSelectedSemi) && ((dCentrality < 20) || (dCentrality > 50))) return;
  
  fEventStatistics->Fill("after centrality selection",1);
  
  // start with track analysis
//   Int_t *iIndexAcceptedTracks = new Int_t[eventAOD->GetNumberOfTracks()]; // maximum number of track indices, this array can have
//   Int_t nTotalNumberAcceptedTracks = 0;
//   for(Int_t i = 0; i < eventAOD->GetNumberOfTracks(); i++) { iIndexAcceptedTracks[i] = 0; }
//   for(Int_t itrack = 0; itrack < eventAOD->GetNumberOfTracks(); itrack++) 
//   { 
// 	track = eventAOD->GetTrack(itrack);
// 	if(!track) continue;
// 	
// 	GetDCA(track, eventAOD, dDCA);
// 	
// 	Double_t dDCAxyDCAzPt[5] = { dDCA[0], dDCA[1], track->Pt(), track->Eta(), track->Phi() };
// 	
// 	fDCAPtAll->Fill(dDCAxyDCAzPt);
// 	
// 	if( !(IsTrackAccepted(track, dCentrality, eventAOD->GetMagneticField())) ) continue;
// 	
// 	iIndexAcceptedTracks[nTotalNumberAcceptedTracks] = itrack;
// 	nTotalNumberAcceptedTracks++;
//   }
  
  // get event plane Angle from AODHeader, default is Q
  ep = const_cast<AliAODEvent*>(eventAOD)->GetEventplane();
  if(ep) {
	dEventplaneAngle = MoveEventplane(ep->GetEventplane(GetEventplaneSelector().Data(),eventAOD));
	if(GetEventplaneSelector().CompareTo("Q") == 0) 
	{
	  epQvector = ep->GetQVector(); 
	  if(epQvector) dEventplaneAngle = epQvector->Phi();//MoveEventplane(epQvector->Phi());
	}
  }
  
  if( (GetEventplaneSelector().CompareTo("Q") == 0) && !epQvector )
  {
	AliWarning("ERROR: epQvector not available \n");
	return;
  }
  
  //   cout << dEventplaneAngle << endl;
  fEventplaneDist->Fill(dEventplaneAngle);
  fEventplaneRunDist->Fill(dEventplaneAngle, (Double_t)eventAOD->GetRunNumber());
  
  // fill crosscheck histos
  fEPDistCent->Fill(dEventplaneAngle, dCentrality);
  fPcosEPCent->Fill(dCentrality, TMath::Cos(2.*dEventplaneAngle));
  fPsinEPCent->Fill(dCentrality, TMath::Sin(2.*dEventplaneAngle));
  
  // start with MC truth analysis
  if(fIsMonteCarlo)
  {
	
	if( dMCEventZv > GetCutMaxZVertex() )  { return; }
	
	dMCTrackZvPtEtaCent[0] = dMCEventZv;
	
	fEventStatistics->Fill("MC afterZv cut",1);
	
	for(Int_t iMCtrack = 0; iMCtrack < stack->GetEntriesFast(); iMCtrack++)
	{
	  mcPart =(AliAODMCParticle*)stack->At(iMCtrack);
	  
	  // check for charge
	  if( !(IsMCTrackAccepted(mcPart)) ) continue;
	  
	  if(!IsHijingParticle(mcPart, genHijingHeader)) { continue; }
	  
	  if(mcPart->IsPhysicalPrimary() ) 
	  {
		// 	fMCHijingPrim->Fill("IsPhysicalPrimary",1);
	  }
	  else
	  {
		// 	fMCHijingPrim->Fill("NOT a primary",1);
		continue;
	  }
	  
	  
	  //       
	  // ======================== fill histograms ========================
	  dMCTrackZvPtEtaCent[1] = mcPart->Pt();
	  dMCTrackZvPtEtaCent[2] = mcPart->Eta();
	  dMCTrackZvPtEtaCent[3] = dCentrality;
	  fMCGenZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
	  
	  dMCTrackPhiPtEtaCent[0] = RotatePhi(mcPart->Phi(), dEventplaneAngle); // use eventplane and not reactionplan, similar to centrality vs impact paramter
	  // 	  if( dMCTrackPhiPtEtaCent[0] < 0) dMCTrackPhiPtEtaCent[0] += 2.*TMath::Pi();
	  // 	  else if( dMCTrackPhiPtEtaCent[0] > 2.*TMath::Pi()) dMCTrackPhiPtEtaCent[0] -= 2.*TMath::Pi();
	  dMCTrackPhiPtEtaCent[1] = mcPart->Pt();
	  dMCTrackPhiPtEtaCent[2] = mcPart->Eta();
	  dMCTrackPhiPtEtaCent[3] = dCentrality;
	  fMCGenDeltaphiPtEtaCent->Fill(dMCTrackPhiPtEtaCent);
	  
	  bEventHasATrack = kTRUE;
	  
	  
	  if( (dMCTrackZvPtEtaCent[1] > GetCutPtMin() ) &&
		(dMCTrackZvPtEtaCent[1] < GetCutPtMax() ) &&
		(dMCTrackZvPtEtaCent[2] > GetCutEtaMin() ) &&
		(dMCTrackZvPtEtaCent[2] < GetCutEtaMax() ) )
	  {
		fMCPt->Fill(mcPart->Pt());
		fMCCharge->Fill(mcPart->Charge()/3.);
		bEventHasATrackInRange = kTRUE;
	  }
	  
	}
  } // isMonteCarlo
  
  if(bEventHasATrack) { fEventStatistics->Fill("MC events with tracks",1); }
  if(bEventHasATrackInRange) 
  { 
	fEventStatistics->Fill("MC events with tracks in range",1); 
	fMCEventStatisticsCentrality->Fill(dCentrality);
  }
  bEventHasATrack = kFALSE;
  bEventHasATrackInRange = kFALSE;
  
  
  //
  // Loop over recontructed tracks
  //
  
  dEventZv = eventAOD->GetPrimaryVertex()->GetZ();
  if( TMath::Abs(dEventZv) > GetCutMaxZVertex() ) return;
  
  // count all events, which are within zv distribution
  fAllEventStatisticsCentrality->Fill(dCentrality/*, nTriggerFired*/);
  
  fEventStatistics->Fill("after Zv cut",1);
  
  dTrackZvPtEtaCent[0] = dEventZv;
  
  
  
  if(AreRelativeCutsEnabled())
  {
	if(!SetRelativeCuts(eventAOD)) return;
  }
  
  Int_t iSubtractedTracks = 0;
  
  for(Int_t itrack = 0; itrack < eventAOD->GetNumberOfTracks(); itrack++)
//   for(Int_t itrack = 0; itrack < nTotalNumberAcceptedTracks; itrack++)
  {
	track = eventAOD->GetTrack(itrack);
// 	track = eventAOD->GetTrack(iIndexAcceptedTracks[itrack]);
	if(!track) continue;
	
	mcPart = NULL;
	dMCTrackZvPtEtaCent[1] = 0;
	dMCTrackZvPtEtaCent[2] = 0;
	dMCTrackZvPtEtaCent[3] = 0;
	
	dMCTrackPhiPtEtaCent[0] = 0;
	dMCTrackPhiPtEtaCent[1] = 0;
	dMCTrackPhiPtEtaCent[2] = 0;
	dMCTrackPhiPtEtaCent[3] = 0;
	
	bIsPrimary = kFALSE;
	
	GetDCA(track, eventAOD, dDCA);
	
	Double_t dDCAxyDCAzPt[5] = { dDCA[0], dDCA[1], track->Pt(), track->Eta(), track->Phi() };
	
	fDCAPtAll->Fill(dDCAxyDCAzPt);
	
	if( !(IsTrackAccepted(track, dCentrality, eventAOD->GetMagneticField())) ) continue;
	
	dTrackZvPtEtaCent[1] = track->Pt();
	dTrackZvPtEtaCent[2] = track->Eta();
	dTrackZvPtEtaCent[3] = dCentrality;
	
	if(GetEventplaneSelector().CompareTo("Q") == 0) 
	{
	  // subtract track contribution from eventplane
	  Double_t dX = -1000;
	  Double_t dY = -1000;
	  
	  dX = epQvector->X();
	  dY = epQvector->Y();
	  if( (dX>-1000) && (dY>-1000) ) // only subtract, if not default!
	  {
		dX -= ep->GetQContributionX(track);
		dY -= ep->GetQContributionY(track);
		iSubtractedTracks++;
	  }
	  TVector2 epCorrected(dX, dY);
	  dEventplaneAngleCorrected = epCorrected.Phi(); // see AlEPSelectionTask.cxx:354 - without dividing by 2!
	}
	else
	{
	  dEventplaneAngleCorrected = dEventplaneAngle; 
	}
	
	Double_t dFillEPCorrectionCheck[] = {dEventplaneAngle, dEventplaneAngleCorrected, dCentrality};
	fCorrelEventplaneDefaultCorrected->Fill(dFillEPCorrectionCheck);
	
	
	dTrackPhiPtEtaCent[0] = RotatePhi(track->Phi(), dEventplaneAngleCorrected); 
	
	// 	if( dTrackPhiPtEtaCent[0] < -1.0*TMath::Pi()) dTrackPhiPtEtaCent[0] += 2.*TMath::Pi();
	// 	else if( dTrackPhiPtEtaCent[0] > TMath::Pi()) dTrackPhiPtEtaCent[0] -= 2.*TMath::Pi();
	dTrackPhiPtEtaCent[1] = track->Pt();
	dTrackPhiPtEtaCent[2] = track->Eta();
	dTrackPhiPtEtaCent[3] = dCentrality;
	
	
	if( fIsMonteCarlo )
	{
	  mcPart = (AliAODMCParticle*)stack->At(TMath::Abs(track->GetLabel()));
	  if( !mcPart ) { continue; }
	  
	  // check for charge
	  // if( !(IsMCTrackAccepted(mcPart)) ) {  continue; } 
	  
	  bIsHijingParticle = IsHijingParticle(mcPart, genHijingHeader);
	  //       bIsPythiaParticle = IsPythiaParticle(mcPart, genPythiaHeader);
	  
	  bIsPrimary = mcPart->IsPhysicalPrimary();
	  
	  dMCTrackZvPtEtaCent[1] = mcPart->Pt();
	  dMCTrackZvPtEtaCent[2] = mcPart->Eta();
	  dMCTrackZvPtEtaCent[3] = dCentrality;
	  
	  dMCTrackPhiPtEtaCent[0] = RotatePhi(mcPart->Phi(), dEventplaneAngle); // use eventplane and not reactionplan, similar to centrality vs impact paramter
	  
	  // 	  if( dMCTrackPhiPtEtaCent[0] < -1.0*TMath::Pi()) dMCTrackPhiPtEtaCent[0] += 2.*TMath::Pi();
	  // 	  else if( dMCTrackPhiPtEtaCent[0] > TMath::Pi()) dMCTrackPhiPtEtaCent[0] -= 2.*TMath::Pi();
	  dMCTrackPhiPtEtaCent[1] = mcPart->Pt();
	  dMCTrackPhiPtEtaCent[2] = mcPart->Eta();
	  dMCTrackPhiPtEtaCent[3] = dCentrality;
	  
	  if(bIsPrimary && bIsHijingParticle)
	  {
		fMCRecPrimZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
		fMCRecPrimDeltaphiPtEtaCent->Fill(dMCTrackPhiPtEtaCent);
		fMCDCAPtPrimary->Fill(dDCAxyDCAzPt);
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
			fMCRecSecZvPtEtaCent->Fill(dMCTrackZvPtEtaCent);
			fMCRecSecDeltaphiPtEtaCent->Fill(dMCTrackPhiPtEtaCent);
			fMCDCAPtSecondary->Fill(dDCAxyDCAzPt);
			// 	  delete moth;
		  }
		} 	
	  }
	} // end isMonteCarlo 
	
	// ======================== fill histograms ========================
	
	// only keep prim and sec from not embedded signal
	Bool_t bKeepMCTrack = kFALSE;
	if(fIsMonteCarlo) 
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
	
	fZvPtEtaCent->Fill(dTrackZvPtEtaCent);
	fDeltaphiPtEtaCent->Fill(dTrackPhiPtEtaCent);
	
	fDCAPtAccepted->Fill(dDCAxyDCAzPt);
	
	if( (dTrackZvPtEtaCent[1] > GetCutPtMin()) &&
	  (dTrackZvPtEtaCent[1] < GetCutPtMax()) &&
	  (dTrackZvPtEtaCent[2] > GetCutEtaMin()) &&
	  (dTrackZvPtEtaCent[2] < GetCutEtaMax()) )
	{
	  iAcceptedMultiplicity++;
	  bEventHasATrackInRange = kTRUE;
	  fPt->Fill(track->Pt());
	  fCharge->Fill(track->Charge());
	  
	  fPhiCent->Fill(track->Phi(), dCentrality);
	  fPcosPhiCent->Fill(dCentrality, TMath::Cos(2.*track->Phi()));
	  fPsinPhiCent->Fill(dCentrality, TMath::Sin(2.*track->Phi()));
	  
	  Double_t deltaphi = track->Phi() - dEventplaneAngleCorrected;
// 	  if(deltaphi > TMath::Pi()) deltaphi -= 2.*TMath::Pi();
	  
	  fDeltaPhiCent->Fill(deltaphi, dCentrality);
	}
  } // end track loop
  
  Int_t iContributorsQVector = ep->GetQContributionXArray()->GetSize();
  if(iContributorsQVector) fEventplaneSubtractedPercentage->Fill((Double_t)iSubtractedTracks/(Double_t)iContributorsQVector, dCentrality);
  
  if(bEventHasATrack) { fEventStatistics->Fill("events with tracks",1); bEventHasATrack = kFALSE; }
  
  if(bEventHasATrackInRange) 
  { 
	fEventStatistics->Fill("events with tracks in range",1); 
	fEventStatisticsCentrality->Fill(dCentrality); 
	
	bEventHasATrackInRange = kFALSE; 
  }
  
  if(bIsEventSelectedMB) fEventStatisticsCentralityTrigger->Fill(dCentrality, 0);
  if(bIsEventSelectedSemi) fEventStatisticsCentralityTrigger->Fill(dCentrality, 1);
  if(bIsEventSelectedCentral) fEventStatisticsCentralityTrigger->Fill(dCentrality, 2);
  
  Double_t dEventZvMultCent[3] = {dEventZv, static_cast<Double_t>(iAcceptedMultiplicity), dCentrality};
  fZvMultCent->Fill(dEventZvMultCent);
  
  // store correlation between data and MC eventplane
  if(fIsMonteCarlo) fCorrelEventplaneMCDATA->Fill(dEventplaneAngle, dMCEventplaneAngle);
  
  PostData(1, fOutputList);
  
  // delete pointers:
//   delete [] iIndexAcceptedTracks;
}

Double_t AlidNdPtAnalysisPbPbAOD::MoveEventplane(Double_t dMCEP)
{
  Double_t retval = 0;
  retval = dMCEP;
  
  if( (dMCEP > 0) && (dMCEP < 1./2.*TMath::Pi()) ) 
  {
	return retval;
  }
  
  if( (dMCEP >= 1./2.*TMath::Pi()) && (dMCEP <= 3./2.*TMath::Pi()) )
  {
	retval -= TMath::Pi();
	return retval;
  }
  
  if(dMCEP > 3./2.*TMath::Pi())
  {
	retval -= 2.*TMath::Pi();
	return retval;
  }
  
  return -9999.;
}

Double_t AlidNdPtAnalysisPbPbAOD::RotatePhi(Double_t phiTrack, Double_t phiEP)
{
  Double_t dPhi = 0;
  dPhi = TMath::Abs(phiTrack - phiEP);
  
  if( dPhi <= TMath::Pi() )
  {
	return dPhi;
  }
  if( dPhi > TMath::Pi() )
  {
	dPhi = TMath::Pi()/2. - dPhi;
	return dPhi;
  }
//   if( (dPhi > TMath::Pi()) && (dPhi <= 3./2.*TMath::Pi()) )
//   {
// 	dPhi = dPhi - TMath::Pi()/2.;
// 	return dPhi;
//   }
//   if( (dPhi > 3./2.*TMath::Pi()) )
//   {
// 	dPhi = dPhi - 3./2.*TMath::Pi();
// 	return dPhi;
//   }
//   if( dPhi < 0 )
//   
//   if ((dPhi >= -1./2. * TMath::Pi() ) && 
// 	(dPhi <= 1./2. * TMath::Pi() ) )
//   {
// 	return dPhi;
//   }
//   
//   if( (dPhi < 0) )
//   {
// 	dPhi += 2.*TMath::Pi();
//   }
//   
//   if ((dPhi > 0) && 
// 	(dPhi > 1./2. * TMath::Pi() ) && 
// 	(dPhi <= 3./2. * TMath::Pi() ) )
//   {
// 	dPhi -= TMath::Pi();
// 	return dPhi;
//   }	
//   
//   if ((dPhi > 0) && 
// 	(dPhi > 3./2. * TMath::Pi() )) 
//   {
// 	dPhi -= 2.*TMath::Pi();
// 	return dPhi;
//   }
  
  //   Printf("[E] dphi = %.4f , phiTrack = %.4f, phiEP = %.4f", dPhi, phiTrack, phiEP);
  
  return -9999.;
}

Bool_t AlidNdPtAnalysisPbPbAOD::SetRelativeCuts(AliAODEvent *event)
{
  //
  // this function determines the absolute cut event-by-event based on the 
  // the percentage given from outside
  //  - cut set on Nclusters and NcrossedRows
  //
  
  if(!event) return kFALSE; 
  
  AliAODTrack *tr = 0x0;
  TH1F *hCluster = new TH1F("hCluster","hCluster",160,0,160);
  TH1F *hCrossed = new TH1F("hCrossed","hCrossed",160,0,160);
  
  for(Int_t itrack = 0; itrack < event->GetNumberOfTracks(); itrack++)
  {
	tr = event->GetTrack(itrack);
	if(!tr) continue;
	
	// do some selection already
	//if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { continue; } 
	
	Double_t dNClustersTPC = tr->GetTPCNcls();
	Double_t dCrossedRowsTPC = tr->GetTPCClusterInfo(2,1);
	
	hCluster->Fill(dNClustersTPC);
	hCrossed->Fill(dCrossedRowsTPC);
  }
  
  // loop trough histogram to check, where percentage is reach
  Double_t dTotIntCluster = hCluster->Integral();
  Double_t dTotIntCrossed = hCrossed->Integral();
  Float_t dIntCluster = 0;
  Float_t dIntCrossed = 0;
  
  if(dTotIntCluster)
  {
	for(Int_t i = 0; i < hCluster->GetNbinsX(); i++)
	{
	  if(hCluster->GetBinCenter(i) < 0) continue;
	  dIntCluster += hCluster->GetBinContent(i);
	  if(dIntCluster/dTotIntCluster > (1-GetCutPercMinNClustersTPC())) 
	  {
		SetCutMinNClustersTPC(hCluster->GetBinCenter(i));
		fCutPercClusters->Fill(hCluster->GetBinCenter(i));
		break;
	  }
	}
  }
  
  if(dTotIntCrossed)
  {
	for(Int_t i = 0; i < hCrossed->GetNbinsX(); i++)
	{
	  if(hCrossed->GetBinCenter(i) < 0) continue;
	  dIntCrossed += hCrossed->GetBinContent(i);
	  if(dIntCrossed/dTotIntCrossed > (1-GetCutPercMinNCrossedRowsTPC())) 
	  {
		SetCutMinNClustersTPC(hCrossed->GetBinCenter(i));
		fCutPercCrossed->Fill(hCrossed->GetBinCenter(i));
		break;
	  }
	}
  }
  
  delete hCrossed;
  delete hCluster;
  return kTRUE;
  
}

Bool_t AlidNdPtAnalysisPbPbAOD::IsTrackAccepted(AliAODTrack *tr, Double_t dCentrality, Double_t bMagZ)
{
  //
  // this function checks the track parameters for quality
  // returns kTRUE if track is accepted
  //
  // - debug histograms (cuts vs pt,eta,phi) are filled in this function
  // - histogram for pt resolution correction are filled here as well
  //
  
  if(!tr) return kFALSE;
   
  if(tr->Charge()==0) { return kFALSE; }
  
  //
  // as done in AliAnalysisTaskFragmentationFunction
  //
  
  Short_t sign = tr->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[21];
  
  for(Int_t i = 0; i < 21; i++) cv[i] = 0;
  for(Int_t i = 0; i < 50; i++) xyz[i] = 0;
  for(Int_t i = 0; i < 50; i++) pxpypz[i] = 0;
  
  tr->GetXYZ(xyz);
  tr->GetPxPyPz(pxpypz);
  tr->GetCovarianceXYZPxPyPz(cv);
  
  // similar error occured as this one:
  // See https://savannah.cern.ch/bugs/?102721
  // which is one of the two 11h re-filtering follow-ups:
  // Andrea Dainese now first does the beam pipe
  // check and then copies from the vtrack (was the other
  // way around) to avoid the crash in the etp::Set()
  
  //   if(xyz[0]*xyz[0]+xyz[1]*xyz[1] > 3.*3.) { return kFALSE; }
  
  AliExternalTrackParam par(xyz, pxpypz, cv, sign);
  //   AliExternalTrackParam *par = new AliExternalTrackParam(xyz, pxpypz, cv, sign); // high mem consumption!!!!
  static AliESDtrack dummy;
  //   Double_t dLength = dummy.GetLengthInActiveZone(par,3,236, -5 ,0,0);
  //   Double_t dLengthInTPC = GetLengthInTPC(tr, 1.8, 220, bMagZ);
  
  Double_t dLengthInTPC = 0;
  if ( DoCutLengthInTPCPtDependent() ) { dLengthInTPC = dummy.GetLengthInActiveZone(&par,3,236, bMagZ ,0,0); }
   
  Double_t dNClustersTPC = tr->GetTPCNcls();
  Double_t dCrossedRowsTPC = tr->GetTPCNCrossedRows();//GetTPCClusterInfo(2,1);
  Double_t dFindableClustersTPC = tr->GetTPCNclsF();
  Double_t dChi2PerClusterTPC = (dNClustersTPC>0)?tr->Chi2perNDF()*(dNClustersTPC-5)/dNClustersTPC:-1.; // see AliDielectronVarManager.h
  Double_t dOneOverPt = tr->OneOverPt();
  Double_t dSigmaOneOverPt = TMath::Sqrt(par.GetSigma1Pt2());
  
  //   hAllCrossedRowsTPC->Fill(dCrossedRowsTPC);
  
  Double_t dCrossedRowsTPCOverFindableClustersTPC = 0;
  if(dFindableClustersTPC) dCrossedRowsTPCOverFindableClustersTPC = dCrossedRowsTPC/dFindableClustersTPC;
  Double_t dCheck[cqMax] = {dCrossedRowsTPC, dNClustersTPC, dChi2PerClusterTPC, dLengthInTPC, dCrossedRowsTPCOverFindableClustersTPC};// = new Double_t[cqMax];
  Double_t dKine[kqMax] = {tr->Pt(), tr->Eta(), tr->Phi()};// = new Double_t[kqMax];
  
  //   dKine[0] = tr->Pt();
  //   dKine[1] = tr->Eta();
  //   dKine[2] = tr->Phi();
  //   
  //   dCheck[0] = dCrossedRowsTPC;
  //   dCheck[1] = dNClustersTPC;
  //   dCheck[2] = dChi2PerClusterTPC;
  
  
  FillDebugHisto(dCheck, dKine, dCentrality, kFALSE);
  
  // first cut on length
  
  if( DoCutLengthInTPCPtDependent() && ( dLengthInTPC < GetPrefactorLengthInTPCPtDependent()*(130-5*TMath::Abs(1./tr->Pt())) )  ) { return kFALSE; }
  
  // filter bit 5
  //   if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobal)) ) { return kFALSE; }
  if(!(tr->TestFilterBit(GetFilterBit())) ) { return kFALSE; } 
  
  // filter bit 4
  //   if(!(tr->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) ) { return kFALSE; }
  
  //   hFilterCrossedRowsTPC->Fill(dCrossedRowsTPC);
  
  
  if(dFindableClustersTPC == 0) {return kFALSE; }
  if(dCrossedRowsTPC < GetCutMinNCrossedRowsTPC()) { return kFALSE; }
  if( (dCrossedRowsTPCOverFindableClustersTPC) < GetCutMinRatioCrossedRowsOverFindableClustersTPC() ) { return kFALSE; }
  if(dNClustersTPC < GetCutMinNClustersTPC()) { return kFALSE; }
  
  if (IsITSRefitRequired() && !(tr->GetStatus() & AliVTrack::kITSrefit)) { return kFALSE; } // no ITS refit
	
	// do a relativ cut in Nclusters, both time at 80% of mean
	//   if(fIsMonteCarlo) 
	//   { 
	  //     if(dNClustersTPC < 88) { return kFALSE; }
	  //   }
	  //   else
	  //   {
		//     if(dNClustersTPC < 76) { return kFALSE; }
		//   }
		
		// fill histogram for pT resolution correction
		Double_t dPtResolutionHisto[3] = { dOneOverPt, dSigmaOneOverPt, dCentrality };
		fPtResptCent->Fill(dPtResolutionHisto);
		
		// fill debug histogram for all accepted tracks
		FillDebugHisto(dCheck, dKine, dCentrality, kTRUE);
		
		// delete pointers
		
		return kTRUE;
}

Bool_t AlidNdPtAnalysisPbPbAOD::FillDebugHisto(Double_t *dCrossCheckVar, Double_t *dKineVar, Double_t dCentrality, Bool_t bIsAccepted)
{
  if(bIsAccepted)
  {
	for(Int_t iCrossCheck = 0; iCrossCheck < cqMax; iCrossCheck++)
	{
	  Double_t dFillIt[5] = {dCrossCheckVar[iCrossCheck], dKineVar[0], dKineVar[1], dKineVar[2], dCentrality};
	  fCrossCheckAcc[iCrossCheck]->Fill(dFillIt);
	}
	
	fCrossCheckRowsLengthAcc->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqCrossedRows]);
	fCrossCheckClusterLengthAcc->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqNcluster]);
  }
  else
  {
	for(Int_t iCrossCheck = 0; iCrossCheck < cqMax; iCrossCheck++)
	{
	  Double_t dFillIt[5] = {dCrossCheckVar[iCrossCheck], dKineVar[0], dKineVar[1], dKineVar[2], dCentrality};
	  fCrossCheckAll[iCrossCheck]->Fill(dFillIt);
	}
	
	fCrossCheckRowsLength->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqCrossedRows]);
	fCrossCheckClusterLength->Fill(dCrossCheckVar[cqLength], dCrossCheckVar[cqNcluster]);
  }
  
  return kTRUE;
  
}

void AlidNdPtAnalysisPbPbAOD::StoreCutSettingsToHistogram()
{
  //
  // this function stores all cut settings to a histograms
  //
  
  fCutSettings->Fill("IsMonteCarlo",fIsMonteCarlo);
  
  fCutSettings->Fill("fCutMaxZVertex", fCutMaxZVertex);
  
  // kinematic cuts
  fCutSettings->Fill("fCutPtMin", fCutPtMin);
  fCutSettings->Fill("fCutPtMax", fCutPtMax);
  fCutSettings->Fill("fCutEtaMin", fCutEtaMin);
  fCutSettings->Fill("fCutEtaMax", fCutEtaMax);
  
  // track quality cut variables
  fCutSettings->Fill("fFilterBit", fFilterBit);
  if(fUseRelativeCuts) fCutSettings->Fill("fUseRelativeCuts", 1);
  if(fCutRequireTPCRefit) fCutSettings->Fill("fCutRequireTPCRefit", 1);
  if(fCutRequireITSRefit) fCutSettings->Fill("fCutRequireITSRefit", 1);
  
  fCutSettings->Fill("fCutMinNumberOfClusters", fCutMinNumberOfClusters);
  fCutSettings->Fill("fCutPercMinNumberOfClusters", fCutPercMinNumberOfClusters);
  fCutSettings->Fill("fCutMinNumberOfCrossedRows", fCutMinNumberOfCrossedRows);
  fCutSettings->Fill("fCutPercMinNumberOfCrossedRows", fCutPercMinNumberOfCrossedRows);
  
  fCutSettings->Fill("fCutMinRatioCrossedRowsOverFindableClustersTPC", fCutMinRatioCrossedRowsOverFindableClustersTPC);
  fCutSettings->Fill("fCutMaxFractionSharedTPCClusters", fCutMaxFractionSharedTPCClusters);
  fCutSettings->Fill("fCutMaxDCAToVertexXY", fCutMaxDCAToVertexXY);
  fCutSettings->Fill("fCutMaxChi2PerClusterITS", fCutMaxChi2PerClusterITS);
  
  if(fCutDCAToVertex2D) fCutSettings->Fill("fCutDCAToVertex2D", 1);
  if(fCutRequireSigmaToVertex) fCutSettings->Fill("fCutRequireSigmaToVertex",1);
  fCutSettings->Fill("fCutMaxDCAToVertexXYPtDepPar0", fCutMaxDCAToVertexXYPtDepPar0);
  fCutSettings->Fill("fCutMaxDCAToVertexXYPtDepPar1", fCutMaxDCAToVertexXYPtDepPar1);
  fCutSettings->Fill("fCutMaxDCAToVertexXYPtDepPar2", fCutMaxDCAToVertexXYPtDepPar2);
  
  if(fCutAcceptKinkDaughters) fCutSettings->Fill("fCutAcceptKinkDaughters", 1);
  fCutSettings->Fill("fCutMaxChi2TPCConstrainedGlobal", fCutMaxChi2TPCConstrainedGlobal);
  if(fCutLengthInTPCPtDependent) fCutSettings->Fill("fCutLengthInTPCPtDependent", 1);
  fCutSettings->Fill("fPrefactorLengthInTPCPtDependent", fPrefactorLengthInTPCPtDependent);
  fCutSettings->Fill(Form("EP selector %s", fEPselector.Data()), 1);
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
  // returns kFALSE if particle is injected
  
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


