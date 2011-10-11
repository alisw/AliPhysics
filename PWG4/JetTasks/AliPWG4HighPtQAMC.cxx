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

//-----------------------------------------------------------------------
// This class compares the global reconstruction with the MC information
// Momentum resolution is stored as function of track cuts and pt.
// Output: Histograms for different set of cuts
//-----------------------------------------------------------------------
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------

#ifndef ALIPWG4HighPtQAMC_CXX
#define ALIPWG4HighPtQAMC_CXX

#include "AliPWG4HighPtQAMC.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TList.h"
#include "TFile.h"
#include "TChain.h"
#include "TH3F.h"
#include "TKey.h"
#include "TSystem.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
//#include "AliAnalysisHelperJetTasks.h"

using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliPWG4HighPtQAMC)

AliPWG4HighPtQAMC::AliPWG4HighPtQAMC()
: AliAnalysisTask("AliPWG4HighPtQAMC", ""), 
  fESD(0), 
  fMC(0),
  fStack(0),
  fVtx(0x0),
  fTrackCuts(0), 
  fTrackCutsReject(),
  fTrackType(0),
  fSigmaConstrainedMax(1e6),
  fPtMax(100.),
  fAvgTrials(1),
  fNEventAll(0),
  fNEventSel(0),
  fNEventReject(0),
  fh1Xsec(0),
  fh1Trials(0),
  fh1PtHard(0),
  fh1PtHardTrials(0),
  fPtAll(0),  
  fPtSel(0),  
  fPtSelFakes(0),
  fNPointTPCFakes(0),
  fPtSelLargeLabel(0),
  fMultRec(0),
  fNPointTPCMultRec(0),
  fDeltaPtMultRec(0),
  fPtAllvsPtMC(0),
  fPtAllminPtMCvsPtAll(0),
  fPtAllvsPtMCvsMult(0),
  fPtAllminPtMCvsPtAllvsMult(0),
  fPtAllminPtMCvsPtAllNPointTPC(0),
  fPtAllminPtMCvsPtAllNPointTPCIter1(0),
  fPtAllminPtMCvsPtAllChi2TPC(0),
  fPtAllminPtMCvsPtAllChi2TPCIter1(0),
  fPtAllminPtMCvsPtAllDCAR(0),
  fPtAllminPtMCvsPtAllDCAZ(0),
  fPtAllminPtMCvsPtAllPhi(0),
  fPtAllminPtMCvsPtAllNPointITS(0),
  fPtAllminPtMCvsPtAllNSigmaToVertex(0),
  fPtAllminPtMCvsPtAllChi2C(0),
  fPtAllminPtMCvsPtAllRel1PtUncertainty(0),
  fPtAllMC(0),
  fPtSelMC(0),
  fHistList(0)
{

  fPtBinEdges[0][0] = 10.;
  fPtBinEdges[0][1] = 1.;
  fPtBinEdges[1][0] = 20.;
  fPtBinEdges[1][1] = 2.;
  fPtBinEdges[2][0] = 100.;
  fPtBinEdges[2][1] = 10.;

}
//________________________________________________________________________
AliPWG4HighPtQAMC::AliPWG4HighPtQAMC(const char *name): 
  AliAnalysisTask(name,""), 
  fESD(0),
  fMC(0),
  fStack(0),
  fVtx(0x0),
  fTrackCuts(),
  fTrackCutsReject(),
  fTrackType(0),
  fSigmaConstrainedMax(1e6),
  fPtMax(100.),
  fAvgTrials(1),
  fNEventAll(0),
  fNEventSel(0),
  fNEventReject(0),
  fh1Xsec(0),
  fh1Trials(0),
  fh1PtHard(0),
  fh1PtHardTrials(0),
  fPtAll(0),
  fPtSel(0),
  fPtSelFakes(0),
  fNPointTPCFakes(0),
  fPtSelLargeLabel(0),
  fMultRec(0),
  fNPointTPCMultRec(0),
  fDeltaPtMultRec(0),
  fPtAllvsPtMC(0),
  fPtAllminPtMCvsPtAll(0),
  fPtAllvsPtMCvsMult(0),
  fPtAllminPtMCvsPtAllvsMult(0),
  fPtAllminPtMCvsPtAllNPointTPC(0),
  fPtAllminPtMCvsPtAllNPointTPCIter1(0),
  fPtAllminPtMCvsPtAllChi2TPC(0),
  fPtAllminPtMCvsPtAllChi2TPCIter1(0),
  fPtAllminPtMCvsPtAllDCAR(0),
  fPtAllminPtMCvsPtAllDCAZ(0),
  fPtAllminPtMCvsPtAllPhi(0),
  fPtAllminPtMCvsPtAllNPointITS(0),
  fPtAllminPtMCvsPtAllNSigmaToVertex(0),
  fPtAllminPtMCvsPtAllChi2C(0),
  fPtAllminPtMCvsPtAllRel1PtUncertainty(0),
  fPtAllMC(0),
  fPtSelMC(0),
  fHistList(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  AliDebug(2,Form("AliPWG4HighPtQAMC Calling Constructor"));

  fPtBinEdges[0][0] = 10.;
  fPtBinEdges[0][1] = 1.;
  fPtBinEdges[1][0] = 20.;
  fPtBinEdges[1][1] = 2.;
  fPtBinEdges[2][0] = 100.;
  fPtBinEdges[2][1] = 10.;

  // Input slot #0 works with a TChain ESD
  DefineInput(0, TChain::Class());
  // Output slot #0, #1 write into a TList
  DefineOutput(0, TList::Class());
}

//________________________________________________________________________
void AliPWG4HighPtQAMC::SetPtBinEdges(Int_t region, Double_t ptmax, Double_t ptBinWidth) {

  if(region<3) {
    fPtBinEdges[region][0] = ptmax;
    fPtBinEdges[region][1] = ptBinWidth;
  }
  else {
    AliError("Only 3 regions alowed. Use region 0/1/2\n");
    return;
  }

}

//________________________________________________________________________
void AliPWG4HighPtQAMC::ConnectInputData(Option_t *) 
{
  // Connect ESD and MC event handler here
  // Called once
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::ConnectInputData \n"));
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    AliDebug(2,Form( "ERROR: Could not read chain from input slot 0 \n"));
    return;
  }

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  if (!esdH) {
    AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
    return;
  } else
    fESD = esdH->GetEvent();
  
  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    AliDebug(2,Form( "ERROR: Could not retrieve MC event handler \n"));
  }
  else
    fMC = eventHandler->MCEvent();

}


//________________________________________________________________________
void AliPWG4HighPtQAMC::CreateOutputObjects() {
  //Create output objects
  AliDebug(2,Form(">> AliPWG4HighPtQATPConly::CreateOutputObjects \n"));

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); 

  OpenFile(0);
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);

  Float_t fgkPtMin = 0.;

  //fPtBinEdges[region][0] = ptmax of region ; fPtBinEdges[region][1] = binWidth of region
  const Float_t ptmin1 =  fgkPtMin;
  const Float_t ptmax1 =  fPtBinEdges[0][0];
  const Float_t ptmin2 =  ptmax1 ;
  const Float_t ptmax2 =  fPtBinEdges[1][0];
  const Float_t ptmin3 =  ptmax2 ;
  const Float_t ptmax3 =  fPtBinEdges[2][0];//fgkPtMax;
  const Int_t nbin11 = (int)((ptmax1-ptmin1)/fPtBinEdges[0][1]);
  const Int_t nbin12 = (int)((ptmax2-ptmin2)/fPtBinEdges[1][1])+nbin11;
  const Int_t nbin13 = (int)((ptmax3-ptmin3)/fPtBinEdges[2][1])+nbin12;
  Int_t fgkNPtBins=nbin13;
  //Create array with low edges of each bin
  Double_t *binsPt=new Double_t[fgkNPtBins+1];
  for(Int_t i=0; i<=fgkNPtBins; i++) {
    if(i<=nbin11) binsPt[i]=(Double_t)ptmin1 + (ptmax1-ptmin1)/nbin11*(Double_t)i ;
    if(i<=nbin12 && i>nbin11) binsPt[i]=(Double_t)ptmin2 + (ptmax2-ptmin2)/(nbin12-nbin11)*((Double_t)i-(Double_t)nbin11) ;
    if(i<=nbin13 && i>nbin12) binsPt[i]=(Double_t)ptmin3 + (ptmax3-ptmin3)/(nbin13-nbin12)*((Double_t)i-(Double_t)nbin12) ;
  }

  Int_t fgkNPhiBins = 18*6;
  Float_t kMinPhi   = 0.;
  Float_t kMaxPhi   = 2.*TMath::Pi();
  Double_t *binsPhi = new Double_t[fgkNPhiBins+1];
  for(Int_t i=0; i<=fgkNPhiBins; i++) binsPhi[i]=(Double_t)kMinPhi + (kMaxPhi-kMinPhi)/fgkNPhiBins*(Double_t)i ;

  Int_t fgkResPtBins=80;
  Float_t kMinResPt = -1.;
  Float_t kMaxResPt = 1.;
  Double_t *binsResPt = new Double_t[fgkResPtBins+1];
  for(Int_t i=0; i<=fgkResPtBins; i++) binsResPt[i]=(Double_t)kMinResPt + (kMaxResPt-kMinResPt)/fgkResPtBins*(Double_t)i ;

  Int_t fgkMultBins=50;
  Float_t kMinMult = 0.;
  Float_t kMaxMult = 3000.;
  Double_t *binsMult = new Double_t[fgkMultBins+1];
  for(Int_t i=0; i<=fgkMultBins; i++) binsMult[i]=(Double_t)kMinMult + (kMaxMult-kMinMult)/fgkMultBins*(Double_t)i ;

  Int_t fgkNEtaBins=20;
  Float_t fgkEtaMin = -1.;
  Float_t fgkEtaMax = 1.;
  Double_t *binsEta=new Double_t[fgkNEtaBins+1];
  for(Int_t i=0; i<=fgkNEtaBins; i++) binsEta[i]=(Double_t)fgkEtaMin + (fgkEtaMax-fgkEtaMin)/fgkNEtaBins*(Double_t)i ;

  Int_t fgkNNClustersTPCBins=80;
  Float_t fgkNClustersTPCMin = 0.5;
  Float_t fgkNClustersTPCMax = 160.5;
  Double_t *binsNClustersTPC=new Double_t[fgkNNClustersTPCBins+1];
  for(Int_t i=0; i<=fgkNNClustersTPCBins; i++) binsNClustersTPC[i]=(Double_t)fgkNClustersTPCMin + (fgkNClustersTPCMax-fgkNClustersTPCMin)/fgkNNClustersTPCBins*(Double_t)i ;

  Int_t fgkNDCA2DBins=80;
  Float_t fgkDCA2DMin = -0.2;
  Float_t fgkDCA2DMax = 0.2;
  if(fTrackType==1 || fTrackType==2 || fTrackType==4) {
    fgkDCA2DMin = -2.;
    fgkDCA2DMax = 2.;
  }
  Double_t *binsDCA2D=new Double_t[fgkNDCA2DBins+1];
  for(Int_t i=0; i<=fgkNDCA2DBins; i++) binsDCA2D[i]=(Double_t)fgkDCA2DMin + (fgkDCA2DMax-fgkDCA2DMin)/fgkNDCA2DBins*(Double_t)i ;

  Int_t fgkNDCAZBins=80;
  Float_t fgkDCAZMin = -2.;
  Float_t fgkDCAZMax = 2.;
  if(fTrackType==1 || fTrackType==2 || fTrackType==4) {
    fgkDCAZMin = -5.;
    fgkDCAZMax = 5.;
  }
  Double_t *binsDCAZ=new Double_t[fgkNDCAZBins+1];
  for(Int_t i=0; i<=fgkNDCAZBins; i++) binsDCAZ[i]=(Double_t)fgkDCAZMin + (fgkDCAZMax-fgkDCAZMin)/fgkNDCAZBins*(Double_t)i ;

 Int_t fgkNNPointITSBins=9;
  Float_t fgkNPointITSMin = -0.5;
  Float_t fgkNPointITSMax = 8.5;
  Double_t *binsNPointITS=new Double_t[fgkNNPointITSBins+1];
  for(Int_t i=0; i<=fgkNNPointITSBins; i++) binsNPointITS[i]=(Double_t)fgkNPointITSMin + (fgkNPointITSMax-fgkNPointITSMin)/fgkNNPointITSBins*(Double_t)i ;

  Int_t fgkNNSigmaToVertexBins=20;
  Float_t fgkNSigmaToVertexMin = 0.;
  Float_t fgkNSigmaToVertexMax = 8.;
  Double_t *binsNSigmaToVertex=new Double_t[fgkNNSigmaToVertexBins+1];
  for(Int_t i=0; i<=fgkNNSigmaToVertexBins; i++) binsNSigmaToVertex[i]=(Double_t)fgkNSigmaToVertexMin + (fgkNSigmaToVertexMax-fgkNSigmaToVertexMin)/fgkNNSigmaToVertexBins*(Double_t)i ;

  Int_t fgkNChi2CBins=20;
  Float_t fgkChi2CMin = 0.;
  Float_t fgkChi2CMax = 100.;
  Double_t *binsChi2C=new Double_t[fgkNChi2CBins+1];
  for(Int_t i=0; i<=fgkNChi2CBins; i++) binsChi2C[i]=(Double_t)fgkChi2CMin + (fgkChi2CMax-fgkChi2CMin)/fgkNChi2CBins*(Double_t)i ;

  Int_t fgkNRel1PtUncertaintyBins=50;
  Float_t fgkRel1PtUncertaintyMin = 0.;
  Float_t fgkRel1PtUncertaintyMax = 1.;
  Double_t *binsRel1PtUncertainty=new Double_t[fgkNRel1PtUncertaintyBins+1];
  for(Int_t i=0; i<=fgkNRel1PtUncertaintyBins; i++) binsRel1PtUncertainty[i]=(Double_t)fgkRel1PtUncertaintyMin + (fgkRel1PtUncertaintyMax-fgkRel1PtUncertaintyMin)/fgkNRel1PtUncertaintyBins*(Double_t)i ;

  Float_t fgkChi2PerClusMin = 0.;
  Float_t fgkChi2PerClusMax = 4.;
  Int_t fgkNChi2PerClusBins = (int)(fgkChi2PerClusMax*10.);
  Double_t *binsChi2PerClus=new Double_t[fgkNChi2PerClusBins+1];
  for(Int_t i=0; i<=fgkNChi2PerClusBins; i++) binsChi2PerClus[i]=(Double_t)fgkChi2PerClusMin + (fgkChi2PerClusMax-fgkChi2PerClusMin)/fgkNChi2PerClusBins*(Double_t)i ;


  fNEventAll = new TH1F("fNEventAll","NEventAll",1,-0.5,0.5);
  fHistList->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel","NEvent Selected for analysis",1,-0.5,0.5);
  fHistList->Add(fNEventSel);
  fNEventReject = new TH1F("fNEventReject","Reason events are rejectected for analysis",20,0,20);
  //Set labels
  fNEventReject->Fill("noESD",0);
  fNEventReject->Fill("Trigger",0);
  fNEventReject->Fill("noMCEvent",0);
  fNEventReject->Fill("noStack",0);
  fNEventReject->Fill("NTracks<2",0);
  fNEventReject->Fill("noVTX",0);
  fNEventReject->Fill("VtxStatus",0);
  fNEventReject->Fill("NCont<2",0);
  fNEventReject->Fill("ZVTX>10",0);
  fHistList->Add(fNEventReject);

  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fHistList->Add(fh1Xsec);

  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fHistList->Add(fh1Trials);

  fh1PtHard       = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",350,-.5,349.5);
  fHistList->Add(fh1PtHard);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",350,-.5,349.5);
  fHistList->Add(fh1PtHardTrials);

  fPtAll = new TH1F("fPtAll","PtAll",fgkNPtBins, binsPt);
  fHistList->Add(fPtAll);
  fPtSel = new TH1F("fPtSel","PtSel",fgkNPtBins, binsPt);
  fHistList->Add(fPtSel);

  fPtSelFakes = new TH1F("fPtSelFakes","PtSelFakes",fgkNPtBins, binsPt);
  fHistList->Add(fPtSelFakes);

  fNPointTPCFakes = new TH1F("fNPointTPCFakes","fNPointTPCFakes",fgkNNClustersTPCBins, binsNClustersTPC);
  fHistList->Add(fNPointTPCFakes);

  fPtSelLargeLabel = new TH1F("fPtSelLargeLabel","PtSelLargeLabel",fgkNPtBins, binsPt);
  fHistList->Add(fPtSelLargeLabel);

  fMultRec = new TH1F("fMultRec","Multiple reconstruction of tracks",fgkMultBins, binsMult);
  fHistList->Add(fMultRec);

  fNPointTPCMultRec = new TH1F("fNPointTPCMultRec","Multiple reconstruction of tracks",fgkNNClustersTPCBins, binsNClustersTPC);
  fHistList->Add(fNPointTPCMultRec);

  fDeltaPtMultRec = new TH2F("fDeltaPtMultRec","Delta pT vs pT for multiple reconstructed tracks",50,0.,50.,40,-20.,20.);
  fHistList->Add(fDeltaPtMultRec);

  fPtAllvsPtMC = new TH2F("fPtAllvsPtMC","fPtAllvsPtMC;p_{T,MC};p_{T,rec}",fgkNPtBins, binsPt,fgkNPtBins, binsPt);
  fHistList->Add(fPtAllvsPtMC);

  fPtAllminPtMCvsPtAll = new TH2F("fPtAllminPtMCvsPtAll","PtAllminPtMCvsPtAll",fgkNPtBins, binsPt,fgkResPtBins,binsResPt);
  fPtAllminPtMCvsPtAll->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAll->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fHistList->Add(fPtAllminPtMCvsPtAll);

  fPtAllvsPtMCvsMult = new TH3F("fPtAllvsPtMCvsMult","fPtAllvsPtMCvsMult;p_{T,MC};p_{T,rec};N_{tracks}",fgkNPtBins, binsPt,fgkNPtBins, binsPt,fgkMultBins,binsMult);
  fHistList->Add(fPtAllvsPtMCvsMult);

  fPtAllminPtMCvsPtAllvsMult = new TH3F("fPtAllminPtMCvsPtAllvsMult","fPtAllminPtMCvsPtAllvsMult",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkMultBins,binsMult);
  fPtAllminPtMCvsPtAllvsMult->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllvsMult->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllvsMult->SetZTitle("N_{tracks}");
  fHistList->Add(fPtAllminPtMCvsPtAllvsMult);
  
  fPtAllminPtMCvsPtAllNPointTPC = new TH3F("fPtAllminPtMCvsPtAllNPointTPC","PtAllminPtMCvsPtAllNPointTPC",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNNClustersTPCBins, binsNClustersTPC);
  fPtAllminPtMCvsPtAllNPointTPC->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllNPointTPC->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllNPointTPC->SetZTitle("N_{cls,TPC}");
  fHistList->Add(fPtAllminPtMCvsPtAllNPointTPC);

  fPtAllminPtMCvsPtAllNPointTPCIter1 = new TH3F("fPtAllminPtMCvsPtAllNPointTPCIter1","PtAllminPtMCvsPtAllNPointTPCIter1",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNNClustersTPCBins, binsNClustersTPC);
  fPtAllminPtMCvsPtAllNPointTPCIter1->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllNPointTPCIter1->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllNPointTPCIter1->SetZTitle("N_{clsIter1,TPC}");
  fHistList->Add(fPtAllminPtMCvsPtAllNPointTPCIter1);

  fPtAllminPtMCvsPtAllChi2TPC = new TH3F("fPtAllminPtMCvsPtAllChi2TPC","PtAllminPtMCvsPtAllChi2TPC",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNChi2PerClusBins, binsChi2PerClus);
  fPtAllminPtMCvsPtAllChi2TPC->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllChi2TPC->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllChi2TPC->SetZTitle("#chi^{2} TPC");
  fHistList->Add(fPtAllminPtMCvsPtAllChi2TPC);

  fPtAllminPtMCvsPtAllChi2TPCIter1 = new TH3F("fPtAllminPtMCvsPtAllChi2TPCIter1","PtAllminPtMCvsPtAllChi2TPCIter1",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNChi2PerClusBins, binsChi2PerClus);
  fPtAllminPtMCvsPtAllChi2TPCIter1->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllChi2TPCIter1->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllChi2TPCIter1->SetZTitle("#chi^{2} TPC Iter1");
  fHistList->Add(fPtAllminPtMCvsPtAllChi2TPCIter1);

  fPtAllminPtMCvsPtAllDCAR = new TH3F("fPtAllminPtMCvsPtAllDCAR","PtAllminPtMCvsPtAllDCAR",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNDCA2DBins,binsDCA2D);
  fPtAllminPtMCvsPtAllDCAR->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllDCAR->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllDCAR->SetZTitle("DCA_{R}");
  fHistList->Add(fPtAllminPtMCvsPtAllDCAR);

  fPtAllminPtMCvsPtAllDCAZ = new TH3F("fPtAllminPtMCvsPtAllDCAZ","PtAllminPtMCvsPtAllDCAZ",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNDCAZBins,binsDCAZ);
  fPtAllminPtMCvsPtAllDCAZ->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllDCAZ->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllDCAZ->SetZTitle("DCA_{Z}");
  fHistList->Add(fPtAllminPtMCvsPtAllDCAZ);

  fPtAllminPtMCvsPtAllPhi = new TH3F("fPtAllminPtMCvsPtAllPhi","PtAllminPtMCvsPtAllPhi",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNPhiBins,binsPhi);
  fPtAllminPtMCvsPtAllPhi->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllPhi->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllPhi->SetZTitle("#phi");
  fHistList->Add(fPtAllminPtMCvsPtAllPhi);

  fPtAllminPtMCvsPtAllNPointITS = new TH3F("fPtAllminPtMCvsPtAllNPointITS","PtAllminPtMCvsPtAllNPointITS",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNNPointITSBins,binsNPointITS);
  fPtAllminPtMCvsPtAllNPointITS->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllNPointITS->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllNPointITS->SetZTitle("N_{point,ITS}}");
  fHistList->Add(fPtAllminPtMCvsPtAllNPointITS);
  
  fPtAllminPtMCvsPtAllNSigmaToVertex = new TH3F("fPtAllminPtMCvsPtAllNSigmaToVertex","PtAllminPtMCvsPtAllNSigmaToVertex",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNNSigmaToVertexBins,binsNSigmaToVertex);
  fPtAllminPtMCvsPtAllNSigmaToVertex->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllNSigmaToVertex->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistList->Add(fPtAllminPtMCvsPtAllNSigmaToVertex);

  fPtAllminPtMCvsPtAllChi2C = new TH3F("fPtAllminPtMCvsPtAllChi2C","PtAllminPtMCvsPtAllChi2C",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNChi2CBins,binsChi2C);
  fPtAllminPtMCvsPtAllChi2C->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllChi2C->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllChi2C->SetZTitle("Constrained #chi^{2}");
  fHistList->Add(fPtAllminPtMCvsPtAllChi2C);

  fPtAllminPtMCvsPtAllRel1PtUncertainty = new TH3F("fPtAllminPtMCvsPtAllRel1PtUncertainty","PtAllminPtMCvsPtAllRel1PtUncertainty",fgkNPtBins, binsPt,fgkResPtBins,binsResPt,fgkNRel1PtUncertaintyBins,binsRel1PtUncertainty);
  fPtAllminPtMCvsPtAllRel1PtUncertainty->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllRel1PtUncertainty->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistList->Add(fPtAllminPtMCvsPtAllRel1PtUncertainty);


  fPtAllMC = new TH1F("fPtAllMC","PtAll",fgkNPtBins, binsPt);
  fHistList->Add(fPtAllMC);
  fPtSelMC = new TH1F("fPtSelMC","PtSel",fgkNPtBins, binsPt);
  fHistList->Add(fPtSelMC);

  TH1::AddDirectory(oldStatus); 

  PostData(0, fHistList);

  if(binsPt)                delete [] binsPt;
  if(binsResPt)             delete [] binsResPt;
  if(binsPhi)               delete [] binsPhi;
  if(binsNClustersTPC)      delete [] binsNClustersTPC;
  if(binsDCA2D)             delete [] binsDCA2D;
  if(binsDCAZ)              delete [] binsDCAZ;
  if(binsNPointITS)         delete [] binsNPointITS;
  if(binsNSigmaToVertex)    delete [] binsNSigmaToVertex;
  if(binsChi2C)             delete [] binsChi2C;
  if(binsEta)               delete [] binsEta;
  if(binsRel1PtUncertainty) delete [] binsRel1PtUncertainty;
  if(binsChi2PerClus)       delete [] binsChi2PerClus;
  if(binsMult)              delete [] binsMult;
  
}

//________________________________________________________________________
Bool_t AliPWG4HighPtQAMC::SelectEvent() {
  //
  // Decide if event should be selected for analysis
  //

  // Checks following requirements:
  // - fESD available
  // - trigger info from AliPhysicsSelection
  // - MCevent available
  // - number of reconstructed tracks > 1
  // - primary vertex reconstructed
  // - z-vertex < 10 cm

  Bool_t selectEvent = kTRUE;

  //fESD object available?
  if (!fESD) {
    AliDebug(2,Form("ERROR: fInputEvent not available\n"));
    fNEventReject->Fill("noESD",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Trigger
  UInt_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(isSelected&AliVEvent::kMB)) { //Select collison candidates
    AliDebug(2,Form(" Trigger Selection: event REJECTED ... "));
    fNEventReject->Fill("Trigger",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //MCEvent available?
  //if yes: get stack
  if(fMC) {
    AliDebug(2,Form("MC particles: %d", fMC->GetNumberOfTracks()));
    fStack = fMC->Stack();                //Particles Stack
    if(fStack) { 
      AliDebug(2,Form("MC particles stack: %d", fStack->GetNtrack()));
    }
    else {
      AliDebug(2,Form("ERROR: Could not retrieve MC eventHandler"));
      fNEventReject->Fill("noStack",1);
      selectEvent = kFALSE;
      return selectEvent;
    }
  } else {
    AliDebug(2,Form("ERROR: Could not retrieve stack"));
    fNEventReject->Fill("noMCEvent",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Check if number of reconstructed tracks is larger than 1
  if(!fESD->GetNumberOfTracks() || fESD->GetNumberOfTracks()<2)  {
    fNEventReject->Fill("NTracks<2",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Check if vertex is reconstructed
  if(fTrackType==1) fVtx = fESD->GetPrimaryVertexTPC();
  else              fVtx = fESD->GetPrimaryVertexSPD();

  if(!fVtx) {
    fNEventReject->Fill("noVTX",1);
    selectEvent = kFALSE;
    return selectEvent;
  }
   
  if(!fVtx->GetStatus()) {
    fNEventReject->Fill("VtxStatus",1);
    selectEvent = kFALSE;
    return selectEvent;
  }
  
  // Need vertex cut
  //  TString vtxName(fVtx->GetName());
  if(fVtx->GetNContributors()<2) {
    fNEventReject->Fill("NCont<2",1);
    selectEvent = kFALSE;
    return selectEvent;
  }
  
  //Check if z-vertex < 10 cm
  double primVtx[3];
  fVtx->GetXYZ(primVtx);
  if(TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2]>10.)){
    fNEventReject->Fill("ZVTX>10",1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  AliDebug(2,Form("Vertex title %s, status %d, nCont %d\n",fVtx->GetTitle(), fVtx->GetStatus(), fVtx->GetNContributors()));

  return selectEvent;

}

//________________________________________________________________________
void AliPWG4HighPtQAMC::Exec(Option_t *) {  
  // Main loop
  // Called for each event
  AliDebug(2,Form(">> AliPWG4HighPtQATPConly::Exec \n"));  
  // All events without selection
  fNEventAll->Fill(0.);

  if(!SelectEvent()) {
    // Post output data
    PostData(0, fHistList);
    return;
  }

  // ---- Get MC Header information (for MC productions in pThard bins) ----
  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data
  
  if(fMC){
    AliGenPythiaEventHeader*  pythiaGenHeader = GetPythiaEventHeader(fMC);
     if(pythiaGenHeader){
       nTrials = pythiaGenHeader->Trials();
       ptHard  = pythiaGenHeader->GetPtHard();
       
       fh1PtHard->Fill(ptHard);
       fh1PtHardTrials->Fill(ptHard,nTrials);
       
       fh1Trials->Fill("#sum{ntrials}",fAvgTrials);
     }
  }

  //Need to keep track of selected events
  fNEventSel->Fill(0.);

  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2,Form("nTracks ESD%d", nTracks));

  int nMCtracks = fStack->GetNtrack();

  Float_t pt      = 0.;
  Float_t ptMC    = 0.;
  Float_t phi     = 0.;
  Float_t dca2D   = 0.;
  Float_t dcaZ    = 0.;
  Int_t nPointITS = 0;
  Float_t chi2C   = 0.;
  Float_t nSigmaToVertex    = 0.;
  Float_t relUncertainty1Pt = 0.;

  int mult = fTrackCuts->CountAcceptedTracks(fESD);

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    AliESDtrack *track = 0;
    AliESDtrack *esdtrack = fESD->GetTrack(iTrack);
    if(!esdtrack) continue;

    if(fTrackType==1) {
      track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
      if(!track) continue;
    }
    else if(fTrackType==2) {
      track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
      if(!track) continue;
      
      AliExternalTrackParam exParam;
      Bool_t relate = track->RelateToVertexTPC(fVtx,fESD->GetMagneticField(),kVeryBig,&exParam);
      if( !relate ) {
	delete track;
	continue;
      }
      track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
    }
    else if(fTrackType==7) {
      //use global constrained track
      track = new AliESDtrack(*esdtrack);
      //      track = esdtrack;
      //      track->Set(esdtrack->GetConstrainedParam()->GetX(),esdtrack->GetConstrainedParam()->GetAlpha(),esdtrack->GetConstrainedParam()->GetParameter(),esdtrack->GetConstrainedParam()->GetCovariance());
    }
    else
      track = esdtrack;
    
    if(fTrackType==2) {
      //Cut on chi2 of constrained fit
      if(track->GetConstrainedChi2TPC() > fSigmaConstrainedMax*fSigmaConstrainedMax) {
	delete track;
	continue;
      }
    }

    Int_t label = TMath::Abs(track->GetLabel());
    if(label>=nMCtracks) {
      if (fTrackCuts->AcceptTrack(track)) { fPtSelLargeLabel->Fill(pt); }
      if(fTrackType==1 || fTrackType==2 || fTrackType==7) delete track;
      continue;
    }

    TParticle *particle = fStack->Particle(label) ;
    if(!particle) {
      if(fTrackType==1 || fTrackType==2 || fTrackType==7) {
	if(track) delete track;
      }
      continue;
    }

    ptMC = particle->Pt();

    if(fTrackType==1 || fTrackType==2) 
      track->GetImpactParametersTPC(dca2D,dcaZ);  //TPConly
    else
      track->GetImpactParameters(dca2D,dcaZ);     //Global
    
    UChar_t itsMap = track->GetITSClusterMap();
    for (Int_t i=0; i < 6; i++) {
      if (itsMap & (1 << i))
	nPointITS ++;
    }
    
    //    fPtAll->Fill(pt);
    fPtAllMC->Fill(ptMC);

    if (fTrackCuts->AcceptTrack(track)) {

      if(fTrackType==7) {
	if(fTrackCutsReject ) {
	  if(fTrackCutsReject->AcceptTrack(track) ) {
	    if(track) delete track;
	    continue;
	  }
	}
	
	if(esdtrack->GetConstrainedParam()) 
	  track->Set(esdtrack->GetConstrainedParam()->GetX(),esdtrack->GetConstrainedParam()->GetAlpha(),esdtrack->GetConstrainedParam()->GetParameter(),esdtrack->GetConstrainedParam()->GetCovariance());
      }

      pt  = track->Pt();
      phi = track->Phi();

      nSigmaToVertex = fTrackCuts->GetSigmaToVertex(track);// Calculates the number of sigma to the vertex for a track.
      chi2C = track->GetConstrainedChi2();
      relUncertainty1Pt = TMath::Sqrt(track->GetSigma1Pt2())*pt;

      fPtSel->Fill(pt);
      if(track->GetLabel()<0) {
        fPtSelFakes->Fill(pt);
        fNPointTPCFakes->Fill(track->GetTPCNcls());
      }
      fPtSelMC->Fill(ptMC);
      fPtAllvsPtMC->Fill(ptMC,pt);
      fPtAllminPtMCvsPtAll->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC) );
      fPtAllminPtMCvsPtAllNPointTPC->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),track->GetTPCNcls());
      fPtAllminPtMCvsPtAllNPointTPCIter1->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),track->GetTPCNclsIter1());
      if(track->GetTPCNcls()>0.) 
	fPtAllminPtMCvsPtAllChi2TPC->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),track->GetTPCchi2()/track->GetTPCNcls());
      if(track->GetTPCNclsIter1()>0.) 
	fPtAllminPtMCvsPtAllChi2TPCIter1->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),track->GetTPCchi2Iter1()/track->GetTPCNclsIter1());
      fPtAllminPtMCvsPtAllDCAR->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dca2D);
      fPtAllminPtMCvsPtAllDCAZ->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dcaZ);
      fPtAllminPtMCvsPtAllPhi->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),phi);
      fPtAllminPtMCvsPtAllNPointITS->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),nPointITS);
      fPtAllminPtMCvsPtAllNSigmaToVertex->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),nSigmaToVertex);
      fPtAllminPtMCvsPtAllChi2C->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),chi2C);
      fPtAllminPtMCvsPtAllRel1PtUncertainty->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),relUncertainty1Pt);

      fPtAllvsPtMCvsMult->Fill(ptMC,pt,mult);
      fPtAllminPtMCvsPtAllvsMult->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC), mult);

      //Check if track is reconstructed multiple times
      /*
      int multCounter = 1;
      for (Int_t iTrack2 = iTrack+1; iTrack2 < nTracks; iTrack2++) {
	//   AliESDtrack *track2 = GetTrackForAnalysis(iTrack2);
	AliESDtrack *track2;
	AliESDtrack *esdtrack2 = fESD->GetTrack(iTrack2);
	if(!esdtrack2) continue;
      
	if(fTrackType==1)
	  track2 = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack2->GetID());
	else if(fTrackType==2) {
	  track2 = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack2->GetID());
	  if(!track2) { 
	    continue;
	  }
	  AliExternalTrackParam exParam2;
	  Bool_t relate = track2->RelateToVertexTPC(fVtx,fESD->GetMagneticField(),kVeryBig,&exParam2);
	  if( !relate ) {
	    delete track2;
	    continue;
	  }
	  track2->Set(exParam2.GetX(),exParam2.GetAlpha(),exParam2.GetParameter(),exParam2.GetCovariance());
	}
	else
	  track2 = esdtrack2;
	if(!track2) {
	  continue;
	}
      
	if (fTrackCuts->AcceptTrack(track2)) {
	  Int_t label2 = TMath::Abs(track2->GetLabel());
	  if(label==label2) {
	    fNPointTPCMultRec->Fill(track->GetTPCNcls());
	    fNPointTPCMultRec->Fill(track2->GetTPCNcls());
	    fDeltaPtMultRec->Fill(track->Pt(),track->Pt()-track2->Pt());
	    multCounter++;
	  }
	}
	if(fTrackType==1 || fTrackType==2) delete track2;
      }//track2 loop
      if(multCounter>1) fMultRec->Fill(multCounter);
      */

    }//fTrackCuts selection


    if(fTrackType==1 || fTrackType==2 || fTrackType==7) {
      if(track) delete track;    
    }

  }//ESD track loop
   
  // Post output data
  PostData(0, fHistList);
  
}
//________________________________________________________________________
Bool_t AliPWG4HighPtQAMC::PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials){
  //
  // get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // This is to called in Notify and should provide the path to the AOD/ESD file
  // Copied from AliAnalysisTaskJetSpectrum2
  //

  TString file(currFile);  
  fXsec = 0;
  fTrials = 1;

  if(file.Contains("root_archive.zip#")){
    Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    file.Replace(pos+1,20,"");
  }
  else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  //  Printf("%s",file.Data());
   

  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec){
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec){
	// not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else{
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if(!key){
	fxsec->Close();
	return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list){
	fxsec->Close();
	return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree){
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliPWG4HighPtQAMC::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // Copied from AliAnalysisTaskJetSpectrum2
  // 

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;

  fAvgTrials = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      //      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    PythiaInfoFromFile(curfile->GetName(),xsection,ftrials);
    fh1Xsec->Fill("<#sigma>",xsection);
    // construct a poor man average trials 
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if(ftrials>=nEntries && nEntries>0.)fAvgTrials = ftrials/nEntries;
  }  
  return kTRUE;
}

//________________________________________________________________________
AliGenPythiaEventHeader*  AliPWG4HighPtQAMC::GetPythiaEventHeader(AliMCEvent *mcEvent){
  
  if(!mcEvent)return 0;
  AliGenEventHeader* genHeader = mcEvent->GenEventHeader();
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
  if(!pythiaGenHeader){
    // cocktail ??
    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(genHeader);
    
    if (!genCocktailHeader) {
      //      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Unknown header type (not Pythia or Cocktail)");
      //      AliWarning(Form("%s %d: Unknown header type (not Pythia or Cocktail)",(char*)__FILE__,__LINE__));
      return 0;
    }
    TList* headerList = genCocktailHeader->GetHeaders();
    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }
    if(!pythiaGenHeader){
      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Pythia event header not found");
      return 0;
    }
  }
  return pythiaGenHeader;

}

//________________________________________________________________________
void AliPWG4HighPtQAMC::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

#endif
