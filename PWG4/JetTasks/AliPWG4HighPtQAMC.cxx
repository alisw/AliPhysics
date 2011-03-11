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
  fTrackCutsITS(0),
  fTrackType(0),
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
  fPtAllminPtMCvsPtAllDCAR(0),
  fPtAllminPtMCvsPtAllDCAZ(0),
  fPtAllminPtMCvsPtAllPhi(0),
  fPtAllminPtMCvsPtAllNPointITS(0),
  fPtAllminPtMCvsPtAllNSigmaToVertex(0),
  fPtAllminPtMCvsPtAllChi2C(0),
  fPtAllminPtMCvsPtAllRel1PtUncertainty(0),
  fPtAllMC(0),
  fPtSelMC(0),
  fPtSelMCITS(0),
  fHistList(0),
  fPtSelITS(0),
  fPtITSminPtMCvsPtITS(0),
  fPtITSminPtMCvsPtITSNPointTPC(0),
  fPtITSminPtMCvsPtITSDCAR(0),
  fPtITSminPtMCvsPtITSDCAZ(0),
  fPtITSminPtMCvsPtITSPhi(0),
  fPtITSminPtMCvsPtITSNPointITS(0),
  fPtITSminPtMCvsPtITSNSigmaToVertex(0),
  fPtITSminPtMCvsPtITSChi2C(0),
  fPtITSminPtMCvsPtITSRel1PtUncertainty(0),
  fHistListITS(0)
{

}
//________________________________________________________________________
AliPWG4HighPtQAMC::AliPWG4HighPtQAMC(const char *name): 
  AliAnalysisTask(name,""), 
  fESD(0),
  fMC(0),
  fStack(0),
  fVtx(0x0),
  fTrackCuts(),
  fTrackCutsITS(),
  fTrackType(0),
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
  fPtAllminPtMCvsPtAllDCAR(0),
  fPtAllminPtMCvsPtAllDCAZ(0),
  fPtAllminPtMCvsPtAllPhi(0),
  fPtAllminPtMCvsPtAllNPointITS(0),
  fPtAllminPtMCvsPtAllNSigmaToVertex(0),
  fPtAllminPtMCvsPtAllChi2C(0),
  fPtAllminPtMCvsPtAllRel1PtUncertainty(0),
  fPtAllMC(0),
  fPtSelMC(0),
  fPtSelMCITS(0),
  fHistList(0),
  fPtSelITS(0),
  fPtITSminPtMCvsPtITS(0),
  fPtITSminPtMCvsPtITSNPointTPC(0),
  fPtITSminPtMCvsPtITSDCAR(0),
  fPtITSminPtMCvsPtITSDCAZ(0),
  fPtITSminPtMCvsPtITSPhi(0),
  fPtITSminPtMCvsPtITSNPointITS(0),
  fPtITSminPtMCvsPtITSNSigmaToVertex(0),
  fPtITSminPtMCvsPtITSChi2C(0),
  fPtITSminPtMCvsPtITSRel1PtUncertainty(0),
  fHistListITS(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  AliDebug(2,Form("AliPWG4HighPtQAMC Calling Constructor"));
  // Input slot #0 works with a TChain ESD
  DefineInput(0, TChain::Class());
  // Output slot #0, #1 write into a TList
  DefineOutput(0, TList::Class());
  DefineOutput(1, TList::Class());
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
  OpenFile(1);
  fHistListITS = new TList();
  fHistListITS->SetOwner(kTRUE);

  Int_t fgkNPhiBins=18;
  Float_t kMinPhi = 0.;
  Float_t kMaxPhi = 2.*TMath::Pi();
  
  Int_t fgkNPtBins=100;
  Float_t fgkPtMin=0.;//2.;
  Float_t fgkPtMax=fPtMax;
  Int_t fgkResPtBins=80;

  Int_t fgkNMultBins = 50;
  Float_t fgkMultMin = 0;
  Float_t fgkMultMax = 500;

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

  fPtAll = new TH1F("fPtAll","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtAll);
  fPtSel = new TH1F("fPtSel","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSel);

  fPtSelFakes = new TH1F("fPtSelFakes","PtSelFakes",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSelFakes);

  fNPointTPCFakes = new TH1F("fNPointTPCFakes","fNPointTPCFakes",160,0.5,160.5);
  fHistList->Add(fNPointTPCFakes);

  fPtSelLargeLabel = new TH1F("fPtSelLargeLabel","PtSelLargeLabel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSelLargeLabel);

  fMultRec = new TH1F("fMultRec","Multiple reconstruction of tracks",20,0,20);
  fHistList->Add(fMultRec);

  fNPointTPCMultRec = new TH1F("fNPointTPCMultRec","Multiple reconstruction of tracks",160,0.5,160.5);
  fHistList->Add(fNPointTPCMultRec);

  fDeltaPtMultRec = new TH2F("fDeltaPtMultRec","Delta pT vs pT for multiple reconstructed tracks",100,0.,50.,100,-20.,20.);
  fHistList->Add(fDeltaPtMultRec);

  fPtAllvsPtMC = new TH2F("fPtAllvsPtMC","fPtAllvsPtMC;p_{T,MC};p_{T,rec}",fgkNPtBins, fgkPtMin,fgkPtMax,fgkNPtBins, fgkPtMin,fgkPtMax);
  fHistList->Add(fPtAllvsPtMC);

  fPtAllminPtMCvsPtAll = new TH2F("fPtAllminPtMCvsPtAll","PtAllminPtMCvsPtAll",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtAllminPtMCvsPtAll->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAll->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fHistList->Add(fPtAllminPtMCvsPtAll);

  fPtAllvsPtMCvsMult = new TH3F("fPtAllvsPtMCvsMult","fPtAllvsPtMCvsMult;p_{T,MC};p_{T,rec};N_{tracks}",fgkNPtBins, fgkPtMin,fgkPtMax,fgkNPtBins, fgkPtMin,fgkPtMax,fgkNMultBins,fgkMultMin,fgkMultMax);
  fHistList->Add(fPtAllvsPtMCvsMult);

  fPtAllminPtMCvsPtAllvsMult = new TH3F("fPtAllminPtMCvsPtAllvsMult","fPtAllminPtMCvsPtAllvsMult",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNMultBins,fgkMultMin,fgkMultMax);
  fPtAllminPtMCvsPtAllvsMult->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllvsMult->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllvsMult->SetZTitle("N_{tracks}");
  fHistList->Add(fPtAllminPtMCvsPtAllvsMult);
  
  fPtAllminPtMCvsPtAllNPointTPC = new TH3F("fPtAllminPtMCvsPtAllNPointTPC","PtAllminPtMCvsPtAllNPointTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,160,0.5,160.5);
  fPtAllminPtMCvsPtAllNPointTPC->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllNPointTPC->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllNPointTPC->SetZTitle("N_{point,TPC}");
  fHistList->Add(fPtAllminPtMCvsPtAllNPointTPC);

  fPtAllminPtMCvsPtAllDCAR = new TH3F("fPtAllminPtMCvsPtAllDCAR","PtAllminPtMCvsPtAllDCAR",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-0.2,0.2);
  fPtAllminPtMCvsPtAllDCAR->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllDCAR->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllDCAR->SetZTitle("DCA_{R}");
  fHistList->Add(fPtAllminPtMCvsPtAllDCAR);

  fPtAllminPtMCvsPtAllDCAZ = new TH3F("fPtAllminPtMCvsPtAllDCAZ","PtAllminPtMCvsPtAllDCAZ",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-2.,2.);
  fPtAllminPtMCvsPtAllDCAZ->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllDCAZ->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllDCAZ->SetZTitle("DCA_{Z}");
  fHistList->Add(fPtAllminPtMCvsPtAllDCAZ);

  fPtAllminPtMCvsPtAllPhi = new TH3F("fPtAllminPtMCvsPtAllPhi","PtAllminPtMCvsPtAllPhi",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtAllminPtMCvsPtAllPhi->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllPhi->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllPhi->SetZTitle("#phi");
  fHistList->Add(fPtAllminPtMCvsPtAllPhi);

  fPtAllminPtMCvsPtAllNPointITS = new TH3F("fPtAllminPtMCvsPtAllNPointITS","PtAllminPtMCvsPtAllNPointITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,9,-0.5,8.5);
  fPtAllminPtMCvsPtAllNPointITS->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllNPointITS->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllNPointITS->SetZTitle("N_{point,ITS}}");
  fHistList->Add(fPtAllminPtMCvsPtAllNPointITS);
  
  fPtAllminPtMCvsPtAllNSigmaToVertex = new TH3F("fPtAllminPtMCvsPtAllNSigmaToVertex","PtAllminPtMCvsPtAllNSigmaToVertex",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,0.,8.);
  fPtAllminPtMCvsPtAllNSigmaToVertex->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllNSigmaToVertex->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistList->Add(fPtAllminPtMCvsPtAllNSigmaToVertex);

  fPtAllminPtMCvsPtAllChi2C = new TH3F("fPtAllminPtMCvsPtAllChi2C","PtAllminPtMCvsPtAllChi2C",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,0.,10.);
  fPtAllminPtMCvsPtAllChi2C->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllChi2C->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllChi2C->SetZTitle("Constrained #chi^{2}");
  fHistList->Add(fPtAllminPtMCvsPtAllChi2C);

  fPtAllminPtMCvsPtAllRel1PtUncertainty = new TH3F("fPtAllminPtMCvsPtAllRel1PtUncertainty","PtAllminPtMCvsPtAllRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,30,0.,0.3);
  fPtAllminPtMCvsPtAllRel1PtUncertainty->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAllRel1PtUncertainty->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtAllminPtMCvsPtAllRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistList->Add(fPtAllminPtMCvsPtAllRel1PtUncertainty);

  //ITSrefit
  fPtSelITS = new TH1F("fPtSelITSrefit","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListITS->Add(fPtSelITS);
  
  fPtITSminPtMCvsPtITS = new TH2F("fPtITSminPtMCvsPtITS","PtITSminPtMCvsPtITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSminPtMCvsPtITS->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITS->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fHistListITS->Add(fPtITSminPtMCvsPtITS);
  
  fPtITSminPtMCvsPtITSNPointTPC = new TH3F("fPtITSminPtMCvsPtITSNPointTPC","PtITSminPtMCvsPtITSNPointTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,160,0.5,160.5);
  fPtITSminPtMCvsPtITSNPointTPC->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSNPointTPC->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtITSminPtMCvsPtITSNPointTPC->SetZTitle("N_{point,TPC}");
  fHistListITS->Add(fPtITSminPtMCvsPtITSNPointTPC);
    
  fPtITSminPtMCvsPtITSDCAR = new TH3F("fPtITSminPtMCvsPtITSDCAR","PtITSminPtMCvsPtITSDCAR",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-0.2,0.2);
  fPtITSminPtMCvsPtITSDCAR->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSDCAR->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtITSminPtMCvsPtITSDCAR->SetZTitle("DCA_{R}");
  fHistListITS->Add(fPtITSminPtMCvsPtITSDCAR);
  
  fPtITSminPtMCvsPtITSDCAZ = new TH3F("fPtITSminPtMCvsPtITSDCAZ","PtITSminPtMCvsPtITSDCAZ",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-2.,2.);
  fPtITSminPtMCvsPtITSDCAZ->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSDCAZ->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtITSminPtMCvsPtITSDCAZ->SetZTitle("DCA_{Z}");
  fHistListITS->Add(fPtITSminPtMCvsPtITSDCAZ);
  
  fPtITSminPtMCvsPtITSPhi = new TH3F("fPtITSminPtMCvsPtITSPhi","PtITSminPtMCvsPtITSPhi",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtITSminPtMCvsPtITSPhi->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSPhi->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtITSminPtMCvsPtITSPhi->SetZTitle("#phi");
  fHistListITS->Add(fPtITSminPtMCvsPtITSPhi);
  
  fPtITSminPtMCvsPtITSNPointITS = new TH3F("fPtITSminPtMCvsPtITSNPointITS","PtITSminPtMCvsPtITSNPointITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,9,-0.5,8.5);
  fPtITSminPtMCvsPtITSNPointITS->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSNPointITS->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})}");
  fPtITSminPtMCvsPtITSNPointITS->SetZTitle("N_{point,ITS}}");
  fHistListITS->Add(fPtITSminPtMCvsPtITSNPointITS); 
  
  fPtITSminPtMCvsPtITSNSigmaToVertex = new TH3F("fPtITSminPtMCvsPtITSNSigmaToVertex","PtITSminPtMCvsPtITSNSigmaToVertex",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,0.,8.);
  fPtITSminPtMCvsPtITSNSigmaToVertex->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSNSigmaToVertex->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtITSminPtMCvsPtITSNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistListITS->Add(fPtITSminPtMCvsPtITSNSigmaToVertex);

  fPtITSminPtMCvsPtITSChi2C = new TH3F("fPtITSminPtMCvsPtITSChi2C","PtITSminPtMCvsPtITSChi2C",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,0.,10.);
  fPtITSminPtMCvsPtITSChi2C->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSChi2C->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtITSminPtMCvsPtITSChi2C->SetZTitle("Constrained #chi^{2}");
  fHistListITS->Add(fPtITSminPtMCvsPtITSChi2C);

  fPtITSminPtMCvsPtITSRel1PtUncertainty = new TH3F("fPtITSminPtMCvsPtITSRel1PtUncertainty","PtITSminPtMCvsPtITSRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,30,0.,0.3);
  fPtITSminPtMCvsPtITSRel1PtUncertainty->SetXTitle("p_{t}^{MC}");
  fPtITSminPtMCvsPtITSRel1PtUncertainty->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fPtITSminPtMCvsPtITSRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistListITS->Add(fPtITSminPtMCvsPtITSRel1PtUncertainty);

  fPtAllMC = new TH1F("fPtAllMC","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtAllMC);
  fPtSelMC = new TH1F("fPtSelMC","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSelMC);
  fPtSelMCITS = new TH1F("fPtSelMCITS","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSelMCITS);
  
  TH1::AddDirectory(oldStatus); 

  PostData(0, fHistList);
  PostData(1, fHistListITS);

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
    PostData(1, fHistListITS);
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

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    AliESDtrack *track;
    AliESDtrack *esdtrack = fESD->GetTrack(iTrack);
    if(!esdtrack) continue;

    if(fTrackType==1 || fTrackType==3)
      track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
    else if(fTrackType==2) {
      track = AliESDtrackCuts::GetTPCOnlyTrack(fESD,esdtrack->GetID());
      if(!track) {
	delete track;
	continue;
      }
      AliExternalTrackParam exParam;
      Bool_t relate = track->RelateToVertexTPC(fVtx,fESD->GetMagneticField(),kVeryBig,&exParam);
      if( !relate ) {
	delete track;
	continue;
      }
      track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),exParam.GetCovariance());
    }
    else
      track = esdtrack;
    
    if(!track) {
      if(fTrackType==1 || fTrackType==2) delete track;
      continue;
    }

    Int_t label = TMath::Abs(track->GetLabel());
    if(label>=nMCtracks) {
      if (fTrackCuts->AcceptTrack(track)) { fPtSelLargeLabel->Fill(pt); }
      if(fTrackType==1 || fTrackType==2) delete track;
      continue;
    }
    TParticle *particle = fStack->Particle(label) ;
    if(!particle) {
      if(fTrackType==1 || fTrackType==2) delete track;
      continue;
    }

    ptMC = particle->Pt();

    if(fTrackType==0) {       //Global
      pt  = track->Pt();
      phi = track->Phi();
      track->GetImpactParameters(dca2D,dcaZ);
    }
    else if(fTrackType==1 || fTrackType==2) {  //TPConly
      pt  = track->Pt();
      phi = track->Phi();
      track->GetImpactParametersTPC(dca2D,dcaZ);
    }
    else {continue;}

    
    UChar_t itsMap = track->GetITSClusterMap();
    for (Int_t i=0; i < 6; i++) {
      if (itsMap & (1 << i))
	nPointITS ++;
    }
    nSigmaToVertex = fTrackCuts->GetSigmaToVertex(track);// Calculates the number of sigma to the vertex for a track.
    chi2C = track->GetConstrainedChi2();
    relUncertainty1Pt = TMath::Sqrt(track->GetSigma1Pt2())*pt;
    
    fPtAll->Fill(pt);
    fPtAllMC->Fill(ptMC);

    if (fTrackCuts->AcceptTrack(track)) {

      fPtSel->Fill(pt);
      if(track->GetLabel()<0) {
        fPtSelFakes->Fill(pt);
        fNPointTPCFakes->Fill(track->GetTPCNcls());
      }
      fPtSelMC->Fill(ptMC);
      fPtAllvsPtMC->Fill(ptMC,pt);
      fPtAllminPtMCvsPtAll->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC) );
      fPtAllminPtMCvsPtAllNPointTPC->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),track->GetTPCNcls());
      fPtAllminPtMCvsPtAllDCAR->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dca2D);
      fPtAllminPtMCvsPtAllDCAZ->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dcaZ);
      fPtAllminPtMCvsPtAllPhi->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),phi);
      fPtAllminPtMCvsPtAllNPointITS->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),nPointITS);
      fPtAllminPtMCvsPtAllNSigmaToVertex->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),nSigmaToVertex);
      fPtAllminPtMCvsPtAllChi2C->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),chi2C);
      fPtAllminPtMCvsPtAllRel1PtUncertainty->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),relUncertainty1Pt);

      int mult = fTrackCuts->CountAcceptedTracks(fESD);
      fPtAllvsPtMCvsMult->Fill(ptMC,pt,mult);
      fPtAllminPtMCvsPtAllvsMult->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC), mult);

      //Check if track is reconstructed multiple times
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
	    delete track2;
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
	  if(fTrackType==1 || fTrackType==2) delete track2;
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

    }//fTrackCuts selection

    //ITSrefit selection
    if (fTrackCutsITS->AcceptTrack(track) && fTrackType==0) {
      
      fPtSelITS->Fill(pt);
      fPtSelMCITS->Fill(ptMC);
      fPtITSminPtMCvsPtITS->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC) );
      fPtITSminPtMCvsPtITSNPointTPC->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),track->GetTPCNcls());
      fPtITSminPtMCvsPtITSDCAR->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dca2D);
      fPtITSminPtMCvsPtITSDCAZ->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dcaZ);
      fPtITSminPtMCvsPtITSPhi->Fill(ptMC,(pt-ptMC)/(pt),phi);
      fPtITSminPtMCvsPtITSNPointITS->Fill(ptMC,(pt-ptMC)/(pt),nPointITS);
      fPtITSminPtMCvsPtITSNSigmaToVertex->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),nSigmaToVertex);
      fPtITSminPtMCvsPtITSChi2C->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),chi2C);
      fPtITSminPtMCvsPtITSRel1PtUncertainty->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),relUncertainty1Pt);
    }//fTrackCutsITS loop

    if(fTrackType==1 || fTrackType==2) delete track;    

  }//ESD track loop
   
  // Post output data
  PostData(0, fHistList);
  PostData(1, fHistListITS);

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
