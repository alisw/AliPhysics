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
#include "TList.h"
#include "TChain.h"
#include "TH3F.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"

using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliPWG4HighPtQAMC)

AliPWG4HighPtQAMC::AliPWG4HighPtQAMC(): AliAnalysisTask("AliPWG4HighPtQAMC", ""), 
  fESD(0), 
  fTrackCuts(0), 
  fTrackCutsITS(0),
  fNEventAll(0),
  fNEventSel(0),
  fPtAll(0),  
  fPtSel(0),  
  fPtAllminPtMCvsPtAll(0),
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
  AliAnalysisTask(name, ""), 
  fESD(0),
  fTrackCuts(),
  fTrackCutsITS(),
  fNEventAll(0),
  fNEventSel(0),
  fPtAll(0),
  fPtSel(0),
  fPtAllminPtMCvsPtAll(0),
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
  // Output slot #0 writes into a TList
  DefineOutput(0, TList::Class());
  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
  // Output slot #2 writes into a TList
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliPWG4HighPtQAMC::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  // Called once
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::ConnectInputData \n"));
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    AliDebug(2,Form("ERROR: Could not read chain from input slot 0"));
  } else {
    
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      AliDebug(2,Form("ERROR: Could not get ESDInputHandler"));
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliPWG4HighPtQAMC::CreateOutputObjects() {
  //Create output objects
  AliDebug(2,Form(">> AliPWG4HighPtQATPConly::CreateOutputObjects \n"));

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); 

  OpenFile(0);
  fHistList = new TList();
  OpenFile(1);
  fHistListITS = new TList();

  Int_t fgkNPhiBins=18;
  Float_t kMinPhi = 0.;
  Float_t kMaxPhi = 2.*TMath::Pi();
  
  Int_t fgkNPtBins=98;
  Float_t fgkPtMin=2.;
  Float_t fgkPtMax=100.;
  Int_t fgkResPtBins=80;

  fNEventAll = new TH1F("fNEventAll","NEventAll",1,-0.5,0.5);
  fHistList->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel","NEvent Selected for analysis",1,-0.5,0.5);
  fHistList->Add(fNEventSel);
  fPtAll = new TH1F("fPtAll","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtAll);
  fPtSel = new TH1F("fPtSel","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSel);
  
  fPtAllminPtMCvsPtAll = new TH2F("fPtAllminPtMCvsPtAll","PtAllminPtMCvsPtAll",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtAllminPtMCvsPtAll->SetXTitle("p_{t}^{MC}");
  fPtAllminPtMCvsPtAll->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{MC})/(1/p_{t}^{MC})");
  fHistList->Add(fPtAllminPtMCvsPtAll);
  
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

}
//________________________________________________________________________
void AliPWG4HighPtQAMC::Exec(Option_t *) {  
  // Main loop
  // Called for each event
  AliDebug(2,Form(">> AliPWG4HighPtQATPConly::Exec \n"));  
  
  // All events without selection
  fNEventAll->Fill(0.);

  if (!fESD) {
    AliDebug(2,Form("ERROR: fESD not available"));
    PostData(0, fHistList);
    PostData(1, fHistListITS);
    return;
  }

  Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!isSelected) { //Select collison candidates
    AliDebug(2,Form(" Trigger Selection: event REJECTED ... "));
    // Post output data
    PostData(0, fHistList);
    PostData(1, fHistListITS);
    return;
  }
  
 AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    AliDebug(2,Form("ERROR: Could not retrieve MC event handler"));
    PostData(0, fHistList);
    PostData(1, fHistListITS);
    return;
  }

  AliMCEvent* mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
    AliDebug(2,Form("ERROR: Could not retrieve MC event"));
    PostData(0, fHistList);
    PostData(1, fHistListITS);
    return;
  }

  AliDebug(2,Form("MC particles: %d", mcEvent->GetNumberOfTracks()));

  if (!fESD) {
    AliDebug(2,Form("ERROR: fESD not available"));
    PostData(0, fHistList);
    PostData(1, fHistListITS);
    return;
  }

  AliStack* stack = mcEvent->Stack();                //Particles Stack

  AliDebug(2,Form("MC particles stack: %d", stack->GetNtrack()));

  const AliESDVertex *vtx = fESD->GetPrimaryVertex();
  // Need vertex cut
  TString vtxName(vtx->GetName());
  if(vtx->GetNContributors() < 2 || (vtxName.Contains("TPCVertex")) ) {
    // SPD vertex
    vtx = fESD->GetPrimaryVertexSPD();
    if(vtx->GetNContributors()<2) {
      vtx = 0x0;
      // Post output data
      PostData(0, fHistList);
      PostData(1, fHistListITS);
      return;
    }
  }

  double primVtx[3];
  vtx->GetXYZ(primVtx);
  //  printf("primVtx: %g  %g  %g \n",primVtx[0],primVtx[1],primVtx[2]);
  if(TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2]>10.)){
    // Post output data
    PostData(0, fHistList);
    PostData(1, fHistListITS);
    return;
  }
  
  AliDebug(2,Form("Vertex title %s, status %d, nCont %d\n",vtx->GetTitle(), vtx->GetStatus(), vtx->GetNContributors()));

  // Need to keep track of evts without vertex
  fNEventSel->Fill(0.);

  if(!fESD->GetNumberOfTracks() || fESD->GetNumberOfTracks()<2)  {
    PostData(0, fHistList);
    PostData(1, fHistListITS);
    return;
  }

  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2,Form("nTracks %d", nTracks));

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    AliESDtrack *track = fESD->GetTrack(iTrack);
    if(!track) continue;
    Int_t label = TMath::Abs(track->GetLabel());
    TParticle *particle = stack->Particle(label) ;
    if(!particle) continue;

    Float_t pt = track->Pt();
    Float_t ptMC = particle->Pt();
    Float_t phi = track->Phi();
    Float_t dca2D, dcaZ;
    track->GetImpactParameters(dca2D,dcaZ);
    UChar_t itsMap = track->GetITSClusterMap();
    Int_t nPointITS = 0;
    for (Int_t i=0; i < 6; i++) {
      if (itsMap & (1 << i))
	nPointITS ++;
    }
    Float_t nSigmaToVertex = fTrackCuts->GetSigmaToVertex(track);// Calculates the number of sigma to the vertex for a track.
    Float_t chi2C = track->GetConstrainedChi2();
    Float_t relUncertainty1Pt = TMath::Sqrt(track->GetSigma1Pt2())*pt;

    fPtAll->Fill(pt);
    fPtAllMC->Fill(ptMC);
    
    if (fTrackCuts->AcceptTrack(track)) {

      fPtSel->Fill(pt);
      
      fPtSelMC->Fill(ptMC);
      fPtAllminPtMCvsPtAll->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC) );
      fPtAllminPtMCvsPtAllNPointTPC->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),track->GetTPCNcls());
      fPtAllminPtMCvsPtAllDCAR->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dca2D);
      fPtAllminPtMCvsPtAllDCAZ->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),dcaZ);
      fPtAllminPtMCvsPtAllPhi->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),phi);
      fPtAllminPtMCvsPtAllNPointITS->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),nPointITS);
      fPtAllminPtMCvsPtAllNSigmaToVertex->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),nSigmaToVertex);
      fPtAllminPtMCvsPtAllChi2C->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),chi2C);
      fPtAllminPtMCvsPtAllRel1PtUncertainty->Fill(ptMC,(1./pt-1./ptMC)/(1./ptMC),relUncertainty1Pt);
    }//fTrackCuts selection
    
    
    //ITSrefit selection
    if (fTrackCutsITS->AcceptTrack(track)) {
      
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
      
  }//ESD track loop
   
  // Post output data
  PostData(0, fHistList);
  PostData(1, fHistListITS);

}
//________________________________________________________________________
void AliPWG4HighPtQAMC::Terminate(Option_t *)
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

#endif
