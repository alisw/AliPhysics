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
// This class compares the global reconstruction with the TPConly 
// reconstruction
// Momentum resolution is stored as function of track cuts and pt.
// Output: Histograms for different set of cuts
//-----------------------------------------------------------------------
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------

#ifndef ALIPWG4HIGHPTQATPCONLY_CXX
#define ALIPWG4HIGHPTQATPCONLY_CXX

#include "AliPWG4HighPtQATPConly.h"

#include "TVector3.h"
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TChain.h"
#include "TH3F.h"
#include <Bytes.h>
#include <TTree.h>

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
//#include "AliAnalysisHelperJetTasks.h"

#include "AliStack.h"
#include "TParticle.h"
#include "TH1I.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliPWG4HighPtQATPConly)

AliPWG4HighPtQATPConly::AliPWG4HighPtQATPConly(): AliAnalysisTask("AliPWG4HighPtQATPConly", ""), 
  fESD(0), 
  fESDfriend(0), 
  fMC(0),
  fTrackCuts(0), 
  fTrackCutsITS(0),
  fMaxCosmicAngle(0.002),
  fNEventAll(0),
  fNEventSel(0),
  fPtAll(0),
  fPtSel(0),
  fPtAllminPtTPCvsPtAll(0),
  fPtAllminPtTPCvsPtAllEtaPos(0),
  fPtAllminPtTPCvsPtAllEtaNeg(0),
  fPtAllminPtTPCvsPtAllNPointTPC(0),
  fPtAllminPtTPCvsPtAllNPointTPCS(0),
  fPtAllminPtTPCvsPtAllDCAR(0),
  fPtAllminPtTPCvsPtAllDCAZ(0),
  fPtAllminPtTPCvsPtAllPhi(0),
  fPtAllminPtTPCvsPtAllNPointITS(0),
  fPtAllminPtTPCvsPtAllNSigmaToVertex(0),
  fPtAllminPtTPCvsPtAllChi2C(0),
  fPtAllminPtTPCvsPtAllRel1PtUncertainty(0),
  fPtAllminPtTPCvsPtAllChi2PerNClusTPC(0),
  fPtAllminPtTPCvsPtAllChi2PerNClusITS(0),
  fPtAllminPtTPCvsNPointTPCPhi(0),
  fPtAllminPtTPCvsNPointITSPhi(0),
  fPtAllminPtTPCvsRel1PtUncertaintyPhi(0),
  fEtaPhiOutliers(0),
  fPtSelITSouter(0),
  fPtITSouterminPtTPCvsPtAll(0),
  fPtITSouterminPtTPCvsPtAllEtaPos(0),
  fPtITSouterminPtTPCvsPtAllEtaNeg(0),
  fPtITSouterminPtTPCvsPtAllNPointTPC(0),
  fPtITSouterminPtTPCvsPtAllNPointTPCS(0),
  fPtITSouterminPtTPCvsPtAllDCAR(0),
  fPtITSouterminPtTPCvsPtAllDCAZ(0),
  fPtITSouterminPtTPCvsPtAllPhi(0),
  fPtITSouterminPtTPCvsPtAllNPointITS(0),
  fPtITSouterminPtTPCvsPtAllNSigmaToVertex(0),
  fPtITSouterminPtTPCvsPtAllChi2C(0),
  fPtITSouterminPtTPCvsPtAllRel1PtUncertainty(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITS(0),
  fPtITSouterminPtTPCvsPtAllITSLayer0(0),
  fPtITSouterminPtTPCvsPtAllITSLayer1(0),
  fPtITSouterminPtTPCvsPtAllITSLayer2(0),
  fPtITSouterminPtTPCvsPtAllITSLayer3(0),
  fPtITSouterminPtTPCvsPtAllITSLayer4(0),
  fPtITSouterminPtTPCvsPtAllITSLayer5(0),
  fPtITSouterminPtTPCvsPtAllNoSPD(0),
  fPtITSouterminPtTPCvsPtAllNoSDD(0),
  fPtITSouterminPtTPCvsPtAllNoSSD(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD(0),
  fHistList(0),
  fPtAllTPC(0),
  fPtSelTPC(0),
  fPtSelTPCITS(0),
  fHistListTPC(0),
  fPtSelITS(0),
  fPtITSminPtTPCvsPtITS(0),
  fPtITSminPtTPCvsPtITSEtaPos(0),
  fPtITSminPtTPCvsPtITSEtaNeg(0),
  fPtITSminPtTPCvsPtITSNPointTPC(0),
  fPtITSminPtTPCvsPtITSNPointTPCS(0),
  fPtITSminPtTPCvsPtITSDCAR(0),
  fPtITSminPtTPCvsPtITSDCAZ(0),
  fPtITSminPtTPCvsPtITSPhi(0),
  fPtITSminPtTPCvsPtITSNPointITS(0),
  fPtITSminPtTPCvsPtITSNSigmaToVertex(0),
  fPtITSminPtTPCvsPtITSChi2C(0),
  fPtITSminPtTPCvsPtITSRel1PtUncertainty(0),
  fPtITSminPtTPCvsPtITSChi2PerNClusTPC(0),
  fPtITSminPtTPCvsPtITSChi2PerNClusITS(0),
  fPtITSminPtTPCvsNPointTPCPhi(0),
  fPtITSminPtTPCvsNPointITSPhi(0),
  fPtITSminPtTPCvsRel1PtUncertaintyPhi(0),
  fPtRel1PtUncertaintyChi2PerClusTPC(0),
  fPtNPointTPCSChi2PerClusTPC(0),
  fPtNPointTPCSRel1PtUncertainty(0),
  fPtRel1PtUncertaintyChi2PerClusTPCSharedSel(0),
  fHistListITS(0),
  fPtCosmicCandidates(0),
  fDeltaPtCosmicCandidates(0),
  fDeltaPhi(0),
  fDeltaEta(0),
  fHistListCosmics(0)
{

}
//________________________________________________________________________
AliPWG4HighPtQATPConly::AliPWG4HighPtQATPConly(const char *name): 
  AliAnalysisTask(name, ""), 
  fESD(0),
  fESDfriend(0), 
  fMC(0),
  fTrackCuts(),
  fTrackCutsITS(),
  fMaxCosmicAngle(0.002),
  fNEventAll(0),
  fNEventSel(0),
  fPtAll(0),
  fPtSel(0),
  fPtAllminPtTPCvsPtAll(0),
  fPtAllminPtTPCvsPtAllEtaPos(0),
  fPtAllminPtTPCvsPtAllEtaNeg(0),
  fPtAllminPtTPCvsPtAllNPointTPC(0),
  fPtAllminPtTPCvsPtAllNPointTPCS(0),
  fPtAllminPtTPCvsPtAllDCAR(0),
  fPtAllminPtTPCvsPtAllDCAZ(0),
  fPtAllminPtTPCvsPtAllPhi(0),
  fPtAllminPtTPCvsPtAllNPointITS(0),
  fPtAllminPtTPCvsPtAllNSigmaToVertex(0),
  fPtAllminPtTPCvsPtAllChi2C(0),
  fPtAllminPtTPCvsPtAllRel1PtUncertainty(0),
  fPtAllminPtTPCvsPtAllChi2PerNClusTPC(0),
  fPtAllminPtTPCvsPtAllChi2PerNClusITS(0),
  fPtAllminPtTPCvsNPointTPCPhi(0),
  fPtAllminPtTPCvsNPointITSPhi(0),
  fPtAllminPtTPCvsRel1PtUncertaintyPhi(0),
  fEtaPhiOutliers(0),
  fPtSelITSouter(0),
  fPtITSouterminPtTPCvsPtAll(0),
  fPtITSouterminPtTPCvsPtAllEtaPos(0),
  fPtITSouterminPtTPCvsPtAllEtaNeg(0),
  fPtITSouterminPtTPCvsPtAllNPointTPC(0),
  fPtITSouterminPtTPCvsPtAllNPointTPCS(0),
  fPtITSouterminPtTPCvsPtAllDCAR(0),
  fPtITSouterminPtTPCvsPtAllDCAZ(0),
  fPtITSouterminPtTPCvsPtAllPhi(0),
  fPtITSouterminPtTPCvsPtAllNPointITS(0),
  fPtITSouterminPtTPCvsPtAllNSigmaToVertex(0),
  fPtITSouterminPtTPCvsPtAllChi2C(0),
  fPtITSouterminPtTPCvsPtAllRel1PtUncertainty(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITS(0),
  fPtITSouterminPtTPCvsPtAllITSLayer0(0),
  fPtITSouterminPtTPCvsPtAllITSLayer1(0),
  fPtITSouterminPtTPCvsPtAllITSLayer2(0),
  fPtITSouterminPtTPCvsPtAllITSLayer3(0),
  fPtITSouterminPtTPCvsPtAllITSLayer4(0),
  fPtITSouterminPtTPCvsPtAllITSLayer5(0),
  fPtITSouterminPtTPCvsPtAllNoSPD(0),
  fPtITSouterminPtTPCvsPtAllNoSDD(0),
  fPtITSouterminPtTPCvsPtAllNoSSD(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD(0),
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD(0),
  fHistList(0),
  fPtAllTPC(0),
  fPtSelTPC(0),
  fPtSelTPCITS(0),
  fHistListTPC(0),
  fPtSelITS(0),
  fPtITSminPtTPCvsPtITS(0),
  fPtITSminPtTPCvsPtITSEtaPos(0),
  fPtITSminPtTPCvsPtITSEtaNeg(0),
  fPtITSminPtTPCvsPtITSNPointTPC(0),
  fPtITSminPtTPCvsPtITSNPointTPCS(0),
  fPtITSminPtTPCvsPtITSDCAR(0),
  fPtITSminPtTPCvsPtITSDCAZ(0),
  fPtITSminPtTPCvsPtITSPhi(0),
  fPtITSminPtTPCvsPtITSNPointITS(0),
  fPtITSminPtTPCvsPtITSNSigmaToVertex(0),
  fPtITSminPtTPCvsPtITSChi2C(0),
  fPtITSminPtTPCvsPtITSRel1PtUncertainty(0),
  fPtITSminPtTPCvsPtITSChi2PerNClusTPC(0),
  fPtITSminPtTPCvsPtITSChi2PerNClusITS(0),
  fPtITSminPtTPCvsNPointTPCPhi(0),
  fPtITSminPtTPCvsNPointITSPhi(0),
  fPtITSminPtTPCvsRel1PtUncertaintyPhi(0),
  fPtRel1PtUncertaintyChi2PerClusTPC(0),
  fPtNPointTPCSChi2PerClusTPC(0),
  fPtNPointTPCSRel1PtUncertainty(0),
  fPtRel1PtUncertaintyChi2PerClusTPCSharedSel(0),
  fHistListITS(0),
  fPtCosmicCandidates(0),
  fDeltaPtCosmicCandidates(0),
  fDeltaPhi(0),
  fDeltaEta(0),
  fHistListCosmics(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliPWG4HighPtQATPConly","Calling Constructor");
  // Input slot #0 works with a TChain ESD
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0, TList::Class());
  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
  // Output slot #2 writes into a TList
  DefineOutput(2, TList::Class());
  // Output slot #3 writes into a TList
  DefineOutput(3, TList::Class());
}

//________________________________________________________________________
void AliPWG4HighPtQATPConly::ConnectInputData(Option_t *) 
{
  // Connect ESD here
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
 // AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*>
 //                        (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    AliDebug(2,Form( "ERROR: Could not retrieve MC event handler \n"));
  }
  else
    fMC = eventHandler->MCEvent();

  //attach the ESD friend
   //  tree->SetBranchStatus("*", kTRUE);
//   tree->SetBranchStatus("Tracks*", kTRUE);
//   tree->SetBranchStatus("ESDfriend*", kTRUE);
  //  fESD->ReadFromTree(tree);

  fESDfriend = (AliESDfriend*)fESD->FindListObject("AliESDfriend");
  if (!fESDfriend)
    {
      // works for both, we just want to avoid setting the branch adress twice
      // in case of the new ESD
      tree->SetBranchAddress("ESDfriend.",&fESDfriend);
    }
  
}

//________________________________________________________________________
void AliPWG4HighPtQATPConly::CreateOutputObjects() {
  //Create output objects
  AliDebug(2,Form(">> AliPWG4HighPtQATPConly::CreateOutputObjects \n")); 

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE); 

  OpenFile(0);
  fHistList = new TList();
  OpenFile(1);
  fHistListTPC = new TList();
  OpenFile(2);
  fHistListITS = new TList();
  OpenFile(3);
  fHistListCosmics = new TList();


  Int_t fgkNPhiBins=36;
  Float_t kMinPhi = 0.;
  Float_t kMaxPhi = 2.*TMath::Pi();
  
  Float_t fgkPtMin=0.;
  Float_t fgkPtMax=100.;
  Int_t fgkNPtBins=(int)(fgkPtMax-fgkPtMin);

  Float_t fgkChi2PerClusMin = 0.;
  Float_t fgkChi2PerClusMax = 4.;
  Int_t fgkChi2PerClusBins = (int)(fgkChi2PerClusMax*10.);
  

  Int_t fgkResPtBins=80;

  fNEventAll = new TH1F("fNEventAll","NEventAll",1,-0.5,0.5);
  fHistList->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel","NEvent Selected for analysis",1,-0.5,0.5);
  fHistList->Add(fNEventSel);
  fPtAll = new TH1F("fPtAll","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtAll);
  fPtSel = new TH1F("fPtSel","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSel);
 
  fPtAllminPtTPCvsPtAll = new TH2F("fPtAllminPtTPCvsPtAll","PtAllminPtTPCvsPtAll",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtAllminPtTPCvsPtAll->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAll->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistList->Add(fPtAllminPtTPCvsPtAll);
  
  fPtAllminPtTPCvsPtAllEtaPos = new TH3F("fPtAllminPtTPCvsPtAllEtaPos","PtAllminPtTPCvsPtAllEtaPos",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtAllminPtTPCvsPtAllEtaPos->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllEtaPos->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistList->Add(fPtAllminPtTPCvsPtAllEtaPos);

  fPtAllminPtTPCvsPtAllEtaNeg = new TH3F("fPtAllminPtTPCvsPtAllEtaNeg","PtAllminPtTPCvsPtAllEtaNeg",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtAllminPtTPCvsPtAllEtaNeg->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllEtaNeg->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistList->Add(fPtAllminPtTPCvsPtAllEtaNeg);

  fPtAllminPtTPCvsPtAllNPointTPC = new TH3F("fPtAllminPtTPCvsPtAllNPointTPC","PtAllminPtTPCvsPtAllNPointTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,0.5,160.5);
  fPtAllminPtTPCvsPtAllNPointTPC->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllNPointTPC->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllNPointTPC->SetZTitle("N_{point,TPC}");
  fHistList->Add(fPtAllminPtTPCvsPtAllNPointTPC);

  fPtAllminPtTPCvsPtAllNPointTPCS = new TH3F("fPtAllminPtTPCvsPtAllNPointTPCS","PtAllminPtTPCvsPtAllNPointTPCS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,100,0.,1.);
  fPtAllminPtTPCvsPtAllNPointTPCS->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllNPointTPCS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllNPointTPCS->SetZTitle("N_{point,TPC}^{Shared}/N_{point,TPC}");
  fHistList->Add(fPtAllminPtTPCvsPtAllNPointTPCS);

  fPtAllminPtTPCvsPtAllDCAR = new TH3F("fPtAllminPtTPCvsPtAllDCAR","PtAllminPtTPCvsPtAllDCAR",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-0.2,0.2);
  fPtAllminPtTPCvsPtAllDCAR->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllDCAR->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllDCAR->SetZTitle("DCA_{R}");
  fHistList->Add(fPtAllminPtTPCvsPtAllDCAR);

  fPtAllminPtTPCvsPtAllDCAZ = new TH3F("fPtAllminPtTPCvsPtAllDCAZ","PtAllminPtTPCvsPtAllDCAZ",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-2.,2.);
  fPtAllminPtTPCvsPtAllDCAZ->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllDCAZ->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllDCAZ->SetZTitle("DCA_{Z}");
  fHistList->Add(fPtAllminPtTPCvsPtAllDCAZ);

  fPtAllminPtTPCvsPtAllPhi = new TH3F("fPtAllminPtTPCvsPtAllPhi","PtAllminPtTPCvsPtAllPhi",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtAllminPtTPCvsPtAllPhi->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllPhi->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllPhi->SetZTitle("#phi");
  fHistList->Add(fPtAllminPtTPCvsPtAllPhi);

  fPtAllminPtTPCvsPtAllNPointITS = new TH3F("fPtAllminPtTPCvsPtAllNPointITS","PtAllminPtTPCvsPtAllNPointITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,9,-0.5,8.5);
  fPtAllminPtTPCvsPtAllNPointITS->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllNPointITS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllNPointITS->SetZTitle("N_{point,ITS}");
  fHistList->Add(fPtAllminPtTPCvsPtAllNPointITS);
  
  fPtAllminPtTPCvsPtAllNSigmaToVertex = new TH3F("fPtAllminPtTPCvsPtAllNSigmaToVertex","PtAllminPtTPCvsPtAllNSigmaToVertex",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,0.,8.);
  fPtAllminPtTPCvsPtAllNSigmaToVertex->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllNSigmaToVertex->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistList->Add(fPtAllminPtTPCvsPtAllNSigmaToVertex);

  fPtAllminPtTPCvsPtAllChi2C = new TH3F("fPtAllminPtTPCvsPtAllChi2C","PtAllminPtTPCvsPtAllChi2C",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,0.,10.);
  fPtAllminPtTPCvsPtAllChi2C->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllChi2C->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllChi2C->SetZTitle("Constrained #chi^{2}");
  fHistList->Add(fPtAllminPtTPCvsPtAllChi2C);

  fPtAllminPtTPCvsPtAllRel1PtUncertainty = new TH3F("fPtAllminPtTPCvsPtAllRel1PtUncertainty","PtAllminPtTPCvsPtAllRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,30,0.,0.3);
  fPtAllminPtTPCvsPtAllRel1PtUncertainty->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllRel1PtUncertainty->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistList->Add(fPtAllminPtTPCvsPtAllRel1PtUncertainty);

  fPtAllminPtTPCvsPtAllChi2PerNClusTPC = new TH3F("fPtAllminPtTPCvsPtAllChi2PerNClusTPC","PtAllminPtTPCvsPtAllChi2PerNClusTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtAllminPtTPCvsPtAllChi2PerNClusTPC->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllChi2PerNClusTPC->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllChi2PerNClusTPC->SetZTitle("#chi^{2}/(2*NClusTPC-5)");
  fHistList->Add(fPtAllminPtTPCvsPtAllChi2PerNClusTPC);

  fPtAllminPtTPCvsPtAllChi2PerNClusITS = new TH3F("fPtAllminPtTPCvsPtAllChi2PerNClusITS","PtAllminPtTPCvsPtAllChi2PerNClusITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtAllminPtTPCvsPtAllChi2PerNClusITS->SetXTitle("p_{t}^{Global}");
  fPtAllminPtTPCvsPtAllChi2PerNClusITS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsPtAllChi2PerNClusITS->SetZTitle("#chi^{2}/(2*NClusITS-5)");
  fHistList->Add(fPtAllminPtTPCvsPtAllChi2PerNClusITS);

  fPtAllminPtTPCvsNPointTPCPhi = new TH3F("fPtAllminPtTPCvsNPointTPCPhi","PtAllminPtTPCvsNPointTPCPhi",fgkResPtBins,-1,1.,80,0.5,160.5,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtAllminPtTPCvsNPointTPCPhi->SetXTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsNPointTPCPhi->SetYTitle("N_{point,TPC}");
  fPtAllminPtTPCvsNPointTPCPhi->SetZTitle("#phi");
  fHistList->Add(fPtAllminPtTPCvsNPointTPCPhi);

  fPtAllminPtTPCvsNPointITSPhi = new TH3F("fPtAllminPtTPCvsNPointITSPhi","PtAllminPtTPCvsNPointITSPhi",fgkResPtBins,-1,1.,9,-0.5,8.5,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtAllminPtTPCvsNPointITSPhi->SetXTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsNPointITSPhi->SetYTitle("N_{point,ITS}");
  fPtAllminPtTPCvsNPointITSPhi->SetZTitle("#phi");
  fHistList->Add(fPtAllminPtTPCvsNPointITSPhi);

  fPtAllminPtTPCvsRel1PtUncertaintyPhi = new TH3F("fPtAllminPtTPCvsRel1PtUncertaintyPhi","PtAllminPtTPCvsRel1PtUncertaintyPhi",fgkResPtBins,-1,1.,30,0.,0.3,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtAllminPtTPCvsRel1PtUncertaintyPhi->SetXTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtAllminPtTPCvsRel1PtUncertaintyPhi->SetYTitle("Rel1PtUncertainty");
  fPtAllminPtTPCvsRel1PtUncertaintyPhi->SetZTitle("#phi");
  fHistList->Add(fPtAllminPtTPCvsRel1PtUncertaintyPhi);

  fEtaPhiOutliers = new TH2F("fEtaPhiOutliers","PtAllminPtTPCvsPtAll",20, -1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fEtaPhiOutliers->SetXTitle("#eta");
  fEtaPhiOutliers->SetYTitle("#phi");
  fHistList->Add(fEtaPhiOutliers);

  //Global vs ITSouter-TPCinner
  fPtSelITSouter = new TH1F("fPtSelITSouter","PtSelITSouter",fgkNPtBins,fgkPtMin,fgkPtMax);
  fHistList->Add(fPtSelITSouter);
  
  fPtITSouterminPtTPCvsPtAll = new TH2F("fPtITSouterminPtTPCvsPtAll","PtITSouterminPtTPCvsPtAll",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAll->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAll->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAll);
  
  fPtITSouterminPtTPCvsPtAllEtaPos = new TH3F("fPtITSouterminPtTPCvsPtAllEtaPos","PtITSouterminPtTPCvsPtAllEtaPos",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtITSouterminPtTPCvsPtAllEtaPos->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllEtaPos->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllEtaPos);

  fPtITSouterminPtTPCvsPtAllEtaNeg = new TH3F("fPtITSouterminPtTPCvsPtAllEtaNeg","PtITSouterminPtTPCvsPtAllEtaNeg",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtITSouterminPtTPCvsPtAllEtaNeg->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllEtaNeg->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllEtaNeg);

  fPtITSouterminPtTPCvsPtAllNPointTPC = new TH3F("fPtITSouterminPtTPCvsPtAllNPointTPC","PtITSouterminPtTPCvsPtAllNPointTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,0.5,160.5);
  fPtITSouterminPtTPCvsPtAllNPointTPC->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllNPointTPC->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllNPointTPC->SetZTitle("N_{point,TPC}");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllNPointTPC);

  fPtITSouterminPtTPCvsPtAllNPointTPCS = new TH3F("fPtITSouterminPtTPCvsPtAllNPointTPCS","PtITSouterminPtTPCvsPtAllNPointTPCS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,100,0.,1.);
  fPtITSouterminPtTPCvsPtAllNPointTPCS->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllNPointTPCS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSouterminPtTPCvsPtAllNPointTPCS->SetZTitle("N_{point,TPC}^{Shared}/N_{point,TPC}");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllNPointTPCS);

  fPtITSouterminPtTPCvsPtAllDCAR = new TH3F("fPtITSouterminPtTPCvsPtAllDCAR","PtITSouterminPtTPCvsPtAllDCAR",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-0.2,0.2);
  fPtITSouterminPtTPCvsPtAllDCAR->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllDCAR->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllDCAR->SetZTitle("DCA_{R}");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllDCAR);

  fPtITSouterminPtTPCvsPtAllDCAZ = new TH3F("fPtITSouterminPtTPCvsPtAllDCAZ","PtITSouterminPtTPCvsPtAllDCAZ",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-2.,2.);
  fPtITSouterminPtTPCvsPtAllDCAZ->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllDCAZ->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllDCAZ->SetZTitle("DCA_{Z}");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllDCAZ);

  fPtITSouterminPtTPCvsPtAllPhi = new TH3F("fPtITSouterminPtTPCvsPtAllPhi","PtITSouterminPtTPCvsPtAllPhi",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtITSouterminPtTPCvsPtAllPhi->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllPhi->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllPhi->SetZTitle("#phi");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllPhi);

  fPtITSouterminPtTPCvsPtAllNPointITS = new TH3F("fPtITSouterminPtTPCvsPtAllNPointITS","PtITSouterminPtTPCvsPtAllNPointITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,9,-0.5,8.5);
  fPtITSouterminPtTPCvsPtAllNPointITS->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllNPointITS->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllNPointITS->SetZTitle("N_{point,ITS}");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllNPointITS);
  
  fPtITSouterminPtTPCvsPtAllNSigmaToVertex = new TH3F("fPtITSouterminPtTPCvsPtAllNSigmaToVertex","PtITSouterminPtTPCvsPtAllNSigmaToVertex",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,0.,8.);
  fPtITSouterminPtTPCvsPtAllNSigmaToVertex->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllNSigmaToVertex->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllNSigmaToVertex);

  fPtITSouterminPtTPCvsPtAllChi2C = new TH3F("fPtITSouterminPtTPCvsPtAllChi2C","PtITSouterminPtTPCvsPtAllChi2C",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,0.,10.);
  fPtITSouterminPtTPCvsPtAllChi2C->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2C->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2C->SetZTitle("Constrained #chi^{2}");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2C);

  fPtITSouterminPtTPCvsPtAllRel1PtUncertainty = new TH3F("fPtITSouterminPtTPCvsPtAllRel1PtUncertainty","PtITSouterminPtTPCvsPtAllRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,30,0.,0.3);
  fPtITSouterminPtTPCvsPtAllRel1PtUncertainty->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllRel1PtUncertainty->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllRel1PtUncertainty);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC","PtITSouterminPtTPCvsPtAllChi2PerNClusTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,35,0.,3.5);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC->SetZTitle("#chi^{2}/(2*NClusTPC-5)");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITS = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITS","PtITSouterminPtTPCvsPtAllChi2PerNClusITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,35,0.,3.5);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITS->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITS->SetZTitle("#chi^{2}/(2*NClusITS-5)");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITS);

  //As function of ITS layers
  fPtITSouterminPtTPCvsPtAllITSLayer0 = new TH2F("fPtITSouterminPtTPCvsPtAllITSLayer0","PtITSouterminPtTPCvsPtAllITSLayer0",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllITSLayer0->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllITSLayer0->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllITSLayer0);

  fPtITSouterminPtTPCvsPtAllITSLayer1 = new TH2F("fPtITSouterminPtTPCvsPtAllITSLayer1","PtITSouterminPtTPCvsPtAllITSLayer1",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllITSLayer1->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllITSLayer1->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllITSLayer1);

  fPtITSouterminPtTPCvsPtAllITSLayer2 = new TH2F("fPtITSouterminPtTPCvsPtAllITSLayer2","PtITSouterminPtTPCvsPtAllITSLayer2",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllITSLayer2->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllITSLayer2->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllITSLayer2);

  fPtITSouterminPtTPCvsPtAllITSLayer3 = new TH2F("fPtITSouterminPtTPCvsPtAllITSLayer3","PtITSouterminPtTPCvsPtAllITSLayer3",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllITSLayer3->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllITSLayer3->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllITSLayer3);

  fPtITSouterminPtTPCvsPtAllITSLayer4 = new TH2F("fPtITSouterminPtTPCvsPtAllITSLayer4","PtITSouterminPtTPCvsPtAllITSLayer4",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllITSLayer4->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllITSLayer4->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllITSLayer4);

  fPtITSouterminPtTPCvsPtAllITSLayer5 = new TH2F("fPtITSouterminPtTPCvsPtAllITSLayer5","PtITSouterminPtTPCvsPtAllITSLayer5",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllITSLayer5->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllITSLayer5->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllITSLayer5);

  fPtITSouterminPtTPCvsPtAllNoSPD = new TH2F("fPtITSouterminPtTPCvsPtAllNoSPD","PtITSouterminPtTPCvsPtAllNoSPD",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllNoSPD->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllNoSPD->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllNoSPD);

  fPtITSouterminPtTPCvsPtAllNoSDD = new TH2F("fPtITSouterminPtTPCvsPtAllNoSDD","PtITSouterminPtTPCvsPtAllNoSDD",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllNoSDD->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllNoSDD->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllNoSDD);

  fPtITSouterminPtTPCvsPtAllNoSSD = new TH2F("fPtITSouterminPtTPCvsPtAllNoSSD","PtITSouterminPtTPCvsPtAllNoSSD",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSouterminPtTPCvsPtAllNoSSD->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllNoSSD->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllNoSSD);

  //
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0 = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0","PtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0->SetZTitle("#chi^{2}/NPointITS");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1 = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1","PtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1->SetZTitle("#chi^{2}/NPointITS");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2 = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2","PtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2->SetZTitle("#chi^{2}/NPointITS");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3 = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3","PtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3->SetZTitle("#chi^{2}/NPointITS");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4 = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4","PtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4->SetZTitle("#chi^{2}/NPointITS");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5 = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5","PtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5->SetZTitle("#chi^{2}/NPointITS");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD","PtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD->SetZTitle("#chi^{2}/(2*NPointITS-5)");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD","PtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD->SetZTitle("#chi^{2}/(2*NPointITS-5)");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD);

  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD = new TH3F("fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD","PtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkChi2PerClusBins,fgkChi2PerClusMin,fgkChi2PerClusMax);
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD->SetXTitle("p_{t}^{Global}");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD->SetYTitle("(1/p_{t}^{ITSouter}-1/p_{t}^{TPCinner})/(1/p_{t}^{ITSouter})");
  fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD->SetZTitle("#chi^{2}/(2*NPointITS-5)");
  fHistList->Add(fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD);


  //ITSrefit
  fPtSelITS = new TH1F("fPtSelITSrefit","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListITS->Add(fPtSelITS);
  
  fPtITSminPtTPCvsPtITS = new TH2F("fPtITSminPtTPCvsPtITS","PtITSminPtTPCvsPtITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSminPtTPCvsPtITS->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistListITS->Add(fPtITSminPtTPCvsPtITS);

  fPtITSminPtTPCvsPtITSEtaPos = new TH3F("fPtITSminPtTPCvsPtITSEtaPos","PtITSminPtTPCvsPtITSEtaPos",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtITSminPtTPCvsPtITSEtaPos->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSEtaPos->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSEtaPos);

  fPtITSminPtTPCvsPtITSEtaNeg = new TH3F("fPtITSminPtTPCvsPtITSEtaNeg","PtITSminPtTPCvsPtITSEtaNeg",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtITSminPtTPCvsPtITSEtaNeg->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSEtaNeg->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSEtaNeg);
  
  fPtITSminPtTPCvsPtITSNPointTPC = new TH3F("fPtITSminPtTPCvsPtITSNPointTPC","PtITSminPtTPCvsPtITSNPointTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,0.5,160.5);
  fPtITSminPtTPCvsPtITSNPointTPC->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSNPointTPC->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSNPointTPC->SetZTitle("N_{point,TPC}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSNPointTPC);

  fPtITSminPtTPCvsPtITSNPointTPCS = new TH3F("fPtITSminPtTPCvsPtITSNPointTPCS","PtITSminPtTPCvsPtITSNPointTPCS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,100,0.,1.);
  fPtITSminPtTPCvsPtITSNPointTPCS->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSNPointTPCS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSNPointTPCS->SetZTitle("N_{point,TPC}^{Shared}/N_{point,TPC}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSNPointTPCS);    

  fPtITSminPtTPCvsPtITSDCAR = new TH3F("fPtITSminPtTPCvsPtITSDCAR","PtITSminPtTPCvsPtITSDCAR",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-0.2,0.2);
  fPtITSminPtTPCvsPtITSDCAR->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSDCAR->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSDCAR->SetZTitle("DCA_{R}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSDCAR);
  
  fPtITSminPtTPCvsPtITSDCAZ = new TH3F("fPtITSminPtTPCvsPtITSDCAZ","PtITSminPtTPCvsPtITSDCAZ",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,80,-2.,2.);
  fPtITSminPtTPCvsPtITSDCAZ->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSDCAZ->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSDCAZ->SetZTitle("DCA_{Z}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSDCAZ);
  
  fPtITSminPtTPCvsPtITSPhi = new TH3F("fPtITSminPtTPCvsPtITSPhi","PtITSminPtTPCvsPtITSPhi",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtITSminPtTPCvsPtITSPhi->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSPhi->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSPhi->SetZTitle("#phi");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSPhi);
  
  fPtITSminPtTPCvsPtITSNPointITS = new TH3F("fPtITSminPtTPCvsPtITSNPointITS","PtITSminPtTPCvsPtITSNPointITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,9,-0.5,8.5);
  fPtITSminPtTPCvsPtITSNPointITS->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSNPointITS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSNPointITS->SetZTitle("N_{point,ITS}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSNPointITS); 
  
  fPtITSminPtTPCvsPtITSNSigmaToVertex = new TH3F("fPtITSminPtTPCvsPtITSNSigmaToVertex","PtITSminPtTPCvsPtITSNSigmaToVertex",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,0.,8.);
  fPtITSminPtTPCvsPtITSNSigmaToVertex->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSNSigmaToVertex->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSNSigmaToVertex);

  fPtITSminPtTPCvsPtITSChi2C = new TH3F("fPtITSminPtTPCvsPtITSChi2C","PtITSminPtTPCvsPtITSChi2C",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,0.,10.);
  fPtITSminPtTPCvsPtITSChi2C->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSChi2C->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSChi2C->SetZTitle("Constrained #chi^{2}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSChi2C);

  fPtITSminPtTPCvsPtITSRel1PtUncertainty = new TH3F("fPtITSminPtTPCvsPtITSRel1PtUncertainty","PtITSminPtTPCvsPtITSRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,30,0.,0.3);
  fPtITSminPtTPCvsPtITSRel1PtUncertainty->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSRel1PtUncertainty->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSRel1PtUncertainty);

  fPtITSminPtTPCvsPtITSChi2PerNClusTPC = new TH3F("fPtITSminPtTPCvsPtITSChi2PerNClusTPC","PtITSminPtTPCvsPtITSChi2PerNClusTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,35,0.,3.5);
  fPtITSminPtTPCvsPtITSChi2PerNClusTPC->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSChi2PerNClusTPC->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSChi2PerNClusTPC->SetZTitle("#chi^{2}/(2*NClusTPC-5)");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSChi2PerNClusTPC);

  fPtITSminPtTPCvsPtITSChi2PerNClusITS = new TH3F("fPtITSminPtTPCvsPtITSChi2PerNClusITS","PtITSminPtTPCvsPtITSChi2PerNClusITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,35,0.,3.5);
  fPtITSminPtTPCvsPtITSChi2PerNClusITS->SetXTitle("p_{t}^{Global}");
  fPtITSminPtTPCvsPtITSChi2PerNClusITS->SetYTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsPtITSChi2PerNClusITS->SetZTitle("#chi^{2}/(2*NClusITS-5)");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSChi2PerNClusITS);

  fPtITSminPtTPCvsNPointTPCPhi = new TH3F("fPtITSminPtTPCvsNPointTPCPhi","PtITSminPtTPCvsNPointTPCPhi",fgkResPtBins,-1,1.,80,0.5,160.5,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtITSminPtTPCvsNPointTPCPhi->SetXTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsNPointTPCPhi->SetYTitle("N_{point,TPC}");
  fPtITSminPtTPCvsNPointTPCPhi->SetZTitle("#phi");
  fHistListITS->Add(fPtITSminPtTPCvsNPointTPCPhi);

  fPtITSminPtTPCvsNPointITSPhi = new TH3F("fPtITSminPtTPCvsNPointITSPhi","PtITSminPtTPCvsNPointITSPhi",fgkResPtBins,-1,1.,9,-0.5,8.5,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtITSminPtTPCvsNPointITSPhi->SetXTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsNPointITSPhi->SetYTitle("N_{point,ITS}");
  fPtITSminPtTPCvsNPointITSPhi->SetZTitle("#phi");
  fHistListITS->Add(fPtITSminPtTPCvsNPointITSPhi);

  fPtITSminPtTPCvsRel1PtUncertaintyPhi = new TH3F("fPtITSminPtTPCvsRel1PtUncertaintyPhi","PtITSminPtTPCvsRel1PtUncertaintyPhi",fgkResPtBins,-1,1.,30,0.,0.3,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtITSminPtTPCvsRel1PtUncertaintyPhi->SetXTitle("(1/p_{t}^{Global}-1/p_{t}^{TPC})/(1/p_{t}^{Global})");
  fPtITSminPtTPCvsRel1PtUncertaintyPhi->SetYTitle("Rel1PtUncertainty");
  fPtITSminPtTPCvsRel1PtUncertaintyPhi->SetZTitle("#phi");
  fHistListITS->Add(fPtITSminPtTPCvsRel1PtUncertaintyPhi);

  fPtRel1PtUncertaintyChi2PerClusTPC = new TH3F("fPtRel1PtUncertaintyChi2PerClusTPC","PtITSminPtTPCvsPtITSRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,30,0.,0.3,35,0.,3.5);
  fPtRel1PtUncertaintyChi2PerClusTPC->SetXTitle("p_{t}^{global}");
  fPtRel1PtUncertaintyChi2PerClusTPC->SetYTitle("Rel1PtUncertainty");
  fPtRel1PtUncertaintyChi2PerClusTPC->SetZTitle("#chi^{2}/(2*N_{clusters}^{TPC}-5)");
  fHistListITS->Add(fPtRel1PtUncertaintyChi2PerClusTPC);

  fPtNPointTPCSChi2PerClusTPC = new TH3F("fPtNPointTPCSChi2PerClusTPC","PtITSminPtTPCvsPtITSNPointTPCS",fgkNPtBins, fgkPtMin,fgkPtMax,100,0.,1.,35,0.,3.5);
  fPtNPointTPCSChi2PerClusTPC->SetXTitle("p_{t}^{global}");
  fPtNPointTPCSChi2PerClusTPC->SetYTitle("N_{Point,TPC}^{Shared}/N_{Point,TPC}");
  fPtNPointTPCSChi2PerClusTPC->SetZTitle("#chi^{2}/(2*N_{clusters}^{TPC}-5)");
  fHistListITS->Add(fPtNPointTPCSChi2PerClusTPC);

  fPtNPointTPCSRel1PtUncertainty = new TH3F("fPtNPointTPCSRel1PtUncertainty","PtITSminPtTPCvsPtITSNPointTPCS",fgkNPtBins, fgkPtMin,fgkPtMax,100,0.,1.,30,0.,0.3);
  fPtNPointTPCSRel1PtUncertainty->SetXTitle("p_{t}^{global}");
  fPtNPointTPCSRel1PtUncertainty->SetYTitle("N_{Point,TPC}^{Shared}/N_{Point,TPC}");
  fPtNPointTPCSRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistListITS->Add(fPtNPointTPCSRel1PtUncertainty);

  fPtRel1PtUncertaintyChi2PerClusTPCSharedSel = new TH3F("fPtRel1PtUncertaintyChi2PerClusTPCSharedSel","PtITSminPtTPCvsPtITSRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,30,0.,0.3,35,0.,3.5);
  fPtRel1PtUncertaintyChi2PerClusTPCSharedSel->SetXTitle("p_{t}^{global}");
  fPtRel1PtUncertaintyChi2PerClusTPCSharedSel->SetYTitle("Rel1PtUncertainty");
  fPtRel1PtUncertaintyChi2PerClusTPCSharedSel->SetZTitle("#chi^{2}/(2*N_{clusters}^{TPC}-5)");
  fHistListITS->Add(fPtRel1PtUncertaintyChi2PerClusTPCSharedSel);

  fPtAllTPC = new TH1F("fPtAllTPC","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListTPC->Add(fPtAllTPC);
  fPtSelTPC = new TH1F("fPtSelTPC","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListTPC->Add(fPtSelTPC);
  fPtSelTPCITS = new TH1F("fPtSelTPCITS","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListTPC->Add(fPtSelTPCITS);

  //****************************************************************************************************************//
  //                                              Cosmic Candidates                                                 //
  //****************************************************************************************************************//
  fPtCosmicCandidates = new TH1F("fPtCosmicCandidates","fPtCosmicCandidates",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListCosmics->Add(fPtCosmicCandidates);  
  fDeltaPtCosmicCandidates = new TH1F("fDeltaPtCosmicCandidates","fDeltaPtCosmicCandidates",fgkNPtBins, -50., 50.);
  fHistListCosmics->Add(fDeltaPtCosmicCandidates);  
  fDeltaPhi = new TH1F("fDeltaPhi","fDeltaPhi",fgkNPhiBins*2,-1*kMaxPhi,kMaxPhi);
  fHistListCosmics->Add(fDeltaPhi);  
  fDeltaEta = new TH1F("fDeltaEta","fDeltaEta",20, -2.,2.);
  fHistListCosmics->Add(fDeltaEta);  

  TH1::AddDirectory(oldStatus);   

}
//________________________________________________________________________
void AliPWG4HighPtQATPConly::Exec(Option_t *) {  
  // Main loop
  // Called for each event
  AliDebug(2,Form(">> AliPWG4HighPtQATPConly::Exec \n"));  

  // All events without selection
  fNEventAll->Fill(0.);

  if (!fESD) {
    AliDebug(2,Form("ERROR: fESD not available"));
    // Post output data
     PostData(0, fHistList);
     PostData(1, fHistListTPC);
     PostData(2, fHistListITS);
     PostData(3, fHistListCosmics);
    return;
  }

  fESD->SetESDfriend(fESDfriend); //Attach the friend to the ESD
  if (!fESDfriend) {
    AliDebug(2,Form("ERROR: fESDfriend not available"));
    // Post output data
     PostData(0, fHistList);
     PostData(1, fHistListTPC);
     PostData(2, fHistListITS);
     PostData(3, fHistListCosmics);
    return;
  }

  Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!isSelected) { //Select collison candidates
    AliDebug(2,Form(" Trigger Selection: event REJECTED ... "));
    // Post output data
     PostData(0, fHistList);
     PostData(1, fHistListTPC);
     PostData(2, fHistListITS);
     PostData(3, fHistListCosmics);
    return;
  }

  //  AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
//   //  AliMCEventHandler* eventHandler = (AliMCEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
  AliStack* stack = 0x0;
  AliMCEvent* mcEvent = 0x0;
  
  if(fMC) {
    mcEvent = fMC;
    if (!mcEvent) {
      AliDebug(2,Form("ERROR: Could not retrieve MC event"));
      PostData(0, fHistList);
      PostData(1, fHistListTPC);
      PostData(2, fHistListITS);
      PostData(3, fHistListCosmics);
      return;
    }
    
    AliDebug(2,Form("MC particles: %d", mcEvent->GetNumberOfTracks()));
    
    stack = mcEvent->Stack();                //Particles Stack
    
    AliDebug(2,Form("MC particles stack: %d", stack->GetNtrack()));
  }
  

  const AliESDVertex *vtx = fESD->GetPrimaryVertexTracks();
  // Need vertex cut
  if (vtx->GetNContributors() < 2) {
    // Post output data
    PostData(0, fHistList);
    PostData(1, fHistListTPC);
    PostData(2, fHistListITS);
    PostData(3, fHistListCosmics);
    return;
  }

  AliDebug(2,Form("Vertex title %s, status %d, nCont %d\n",vtx->GetTitle(), vtx->GetStatus(), vtx->GetNContributors()));
  double primVtx[3];
  vtx->GetXYZ(primVtx);
  //  printf("primVtx: %g  %g  %g \n",primVtx[0],primVtx[1],primVtx[2]);
  if(TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2]>10.)){
    // Post output data
    PostData(0, fHistList);
    PostData(1, fHistListTPC);
    PostData(2, fHistListITS);
    PostData(3, fHistListCosmics);
    return;
  }
  if(!fESD->GetNumberOfTracks() || fESD->GetNumberOfTracks()<2){ 
    // Post output data
    PostData(0, fHistList);
    PostData(1, fHistListTPC);
    PostData(2, fHistListITS);
    PostData(3, fHistListCosmics);
    return;
  }
  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2,Form("nTracks %d\n", nTracks));

  if(!fTrackCuts) return;

  // Selected events for analysis
  fNEventSel->Fill(0.);

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    AliESDtrack *track = fESD->GetTrack(iTrack);
    AliExternalTrackParam *trackTPC = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!track || !trackTPC) continue;
    const AliESDfriendTrack* constfriendtrack = track->GetFriendTrack();
    if (!constfriendtrack) { continue;}
    AliESDfriendTrack friendtrack(*constfriendtrack);
 
    Float_t pt = track->Pt();
    Float_t ptTPC = trackTPC->Pt();
    Float_t phi = track->Phi();
    Float_t dca2D, dcaZ;
    track->GetImpactParameters(dca2D,dcaZ);
    // Float_t dca2DTPC, dcaZTPC;
    //track->GetImpactParametersTPC(dca2DTPC,dcaZTPC);
    UChar_t itsMap = track->GetITSClusterMap();
    Int_t nPointITS = 0;
    for (Int_t i=0; i < 6; i++) {
      if (itsMap & (1 << i))
	nPointITS ++;
    }
    double mom[3];
    track->GetPxPyPz(mom);
    double momTPC[3];
    trackTPC->GetPxPyPz(momTPC);
    Float_t nSigmaToVertex = fTrackCuts->GetSigmaToVertex(track);// Calculates the number of sigma to the vertex for a track.
    Float_t chi2C = track->GetConstrainedChi2();
    Float_t relUncertainty1Pt = TMath::Sqrt(TMath::Abs(track->GetSigma1Pt2()))*pt;
    Float_t chi2PerClusterTPC = -1.;
    Float_t nClustersTPC = track->GetTPCNcls();
    if(nClustersTPC>0.) chi2PerClusterTPC = track->GetTPCchi2()/(2.*nClustersTPC-5.);
    Float_t chi2PerNPointITS = -1.;
    if(nPointITS>3) chi2PerNPointITS = track->GetITSchi2()/(2.*(float)nPointITS-5.);

    fPtAll->Fill(pt);
    fPtAllTPC->Fill(ptTPC);


    if (fTrackCuts->AcceptTrack(track)) {

      Bool_t cosmic = kFALSE;    
      if(pt>6.) { cosmic = IsCosmic(track,iTrack,6.); }
      //    if(cosmic) continue;

      fPtSel->Fill(pt);
      fPtSelTPC->Fill(ptTPC);
      if(ptTPC==0. || pt==0.) continue;
      fPtAllminPtTPCvsPtAll->Fill(pt,(1./pt-1./ptTPC)/(1./pt) );

      if(track->GetSign()>0.) fPtAllminPtTPCvsPtAllEtaPos->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->Eta());
      if(track->GetSign()<0.) fPtAllminPtTPCvsPtAllEtaNeg->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->Eta());

      fPtAllminPtTPCvsPtAllNPointTPC->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nClustersTPC);
      if(nClustersTPC>0.) fPtAllminPtTPCvsPtAllNPointTPCS->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->GetTPCnclsS()/nClustersTPC);
      fPtAllminPtTPCvsPtAllDCAR->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dca2D);
      fPtAllminPtTPCvsPtAllDCAZ->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dcaZ);
      fPtAllminPtTPCvsPtAllPhi->Fill(pt,(1./pt-1./ptTPC)/(1./pt),phi);
      fPtAllminPtTPCvsPtAllNPointITS->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nPointITS);
      fPtAllminPtTPCvsPtAllNSigmaToVertex->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nSigmaToVertex);
      fPtAllminPtTPCvsPtAllChi2C->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2C);
      fPtAllminPtTPCvsPtAllRel1PtUncertainty->Fill(pt,(1./pt-1./ptTPC)/(1./pt),relUncertainty1Pt);
      fPtAllminPtTPCvsPtAllChi2PerNClusTPC->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2PerClusterTPC);
      if(nPointITS>3) fPtAllminPtTPCvsPtAllChi2PerNClusITS->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2PerNPointITS);
      fPtAllminPtTPCvsNPointTPCPhi->Fill((1./pt-1./ptTPC)/(1./pt),nClustersTPC,phi);
      fPtAllminPtTPCvsNPointITSPhi->Fill((1./pt-1./ptTPC)/(1./pt),nPointITS,phi);
      fPtAllminPtTPCvsRel1PtUncertaintyPhi->Fill((1./pt-1./ptTPC)/(1./pt),relUncertainty1Pt,phi);

      if(TMath::Abs((1./pt-1./ptTPC)/(1./pt))>0.8) fEtaPhiOutliers->Fill(track->Eta(),phi);
      if (friendtrack.GetITSOut()) {
	AliExternalTrackParam trackITSouter(*(friendtrack.GetITSOut())); 
	Float_t ptITSouter = trackITSouter.Pt();
	if(ptITSouter==0.) continue;
	fPtSelITSouter->Fill(ptITSouter);
	fPtITSouterminPtTPCvsPtAll->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter) );

	if(trackITSouter.GetSign()>0.) fPtAllminPtTPCvsPtAllEtaPos->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),trackITSouter.Eta());
	if(trackITSouter.GetSign()<0.) fPtAllminPtTPCvsPtAllEtaNeg->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),trackITSouter.Eta());

	fPtITSouterminPtTPCvsPtAllNPointTPC->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),nClustersTPC);
	if(nClustersTPC>0.) fPtITSouterminPtTPCvsPtAllNPointTPCS->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),track->GetTPCnclsS()/nClustersTPC);
	fPtITSouterminPtTPCvsPtAllDCAR->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),dca2D);
	fPtITSouterminPtTPCvsPtAllDCAZ->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),dcaZ);
	fPtITSouterminPtTPCvsPtAllPhi->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),phi);
	fPtITSouterminPtTPCvsPtAllNPointITS->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),nPointITS);
	fPtITSouterminPtTPCvsPtAllNSigmaToVertex->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),nSigmaToVertex);
	fPtITSouterminPtTPCvsPtAllChi2C->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2C);
	fPtITSouterminPtTPCvsPtAllRel1PtUncertainty->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),relUncertainty1Pt);
	fPtITSouterminPtTPCvsPtAllChi2PerNClusTPC->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerClusterTPC);
	if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITS->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	if(track->HasPointOnITSLayer(0)) {
	  fPtITSouterminPtTPCvsPtAllITSLayer0->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer0->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
	if(!track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1)) {
	  fPtITSouterminPtTPCvsPtAllITSLayer1->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer1->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1) && track->HasPointOnITSLayer(2)) {
	  fPtITSouterminPtTPCvsPtAllITSLayer2->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer2->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1) && !track->HasPointOnITSLayer(2) && track->HasPointOnITSLayer(3)) {
	  fPtITSouterminPtTPCvsPtAllITSLayer3->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer3->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1) && !track->HasPointOnITSLayer(2) && !track->HasPointOnITSLayer(3) && track->HasPointOnITSLayer(4)) {
	  fPtITSouterminPtTPCvsPtAllITSLayer4->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer4->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1) && !track->HasPointOnITSLayer(2) && !track->HasPointOnITSLayer(3) && !track->HasPointOnITSLayer(4) && track->HasPointOnITSLayer(5)) {
	  fPtITSouterminPtTPCvsPtAllITSLayer5->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSLayer5->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}

	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) {
	  fPtITSouterminPtTPCvsPtAllNoSPD->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSPD->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
	if(!track->HasPointOnITSLayer(2) && !track->HasPointOnITSLayer(3)) {
	  fPtITSouterminPtTPCvsPtAllNoSDD->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSDD->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
	if(!track->HasPointOnITSLayer(4) && !track->HasPointOnITSLayer(5)) {
	  fPtITSouterminPtTPCvsPtAllNoSSD->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter));
	  if(nPointITS>3) fPtITSouterminPtTPCvsPtAllChi2PerNClusITSNoSSD->Fill(pt,(1./ptITSouter-1./ptTPC)/(1./ptITSouter),chi2PerNPointITS);
	}
      }

    }//fTrackCuts selection
    
    
    //ITSrefit selection
    if (fTrackCutsITS->AcceptTrack(track)) {
      
      fPtSelITS->Fill(pt);
      fPtSelTPCITS->Fill(ptTPC);
      fPtITSminPtTPCvsPtITS->Fill(pt,(1./pt-1./ptTPC)/(1./pt) );
      if(track->GetSign()>0.) fPtITSminPtTPCvsPtITSEtaPos->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->Eta());
      if(track->GetSign()<0.) fPtITSminPtTPCvsPtITSEtaNeg->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->Eta());
      fPtITSminPtTPCvsPtITSNPointTPC->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nClustersTPC);
      if(nClustersTPC>0.) fPtITSminPtTPCvsPtITSNPointTPCS->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->GetTPCnclsS()/nClustersTPC);
      fPtITSminPtTPCvsPtITSDCAR->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dca2D);
      fPtITSminPtTPCvsPtITSDCAZ->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dcaZ);
      fPtITSminPtTPCvsPtITSPhi->Fill(pt,(1./pt-1./ptTPC)/(1./pt),phi);
      fPtITSminPtTPCvsPtITSNPointITS->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nPointITS);
      fPtITSminPtTPCvsPtITSNSigmaToVertex->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nSigmaToVertex);
      fPtITSminPtTPCvsPtITSChi2C->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2C);
      fPtITSminPtTPCvsPtITSRel1PtUncertainty->Fill(pt,(1./pt-1./ptTPC)/(1./pt),relUncertainty1Pt);
      fPtITSminPtTPCvsPtITSChi2PerNClusTPC->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2PerClusterTPC);
      if(nPointITS>3) fPtITSminPtTPCvsPtITSChi2PerNClusITS->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2PerNPointITS);
      fPtITSminPtTPCvsNPointTPCPhi->Fill((1./pt-1./ptTPC)/(1./pt),nClustersTPC,phi);
      fPtITSminPtTPCvsNPointITSPhi->Fill((1./pt-1./ptTPC)/(1./pt),nPointITS,phi);
      fPtITSminPtTPCvsRel1PtUncertaintyPhi->Fill((1./pt-1./ptTPC)/(1./pt),relUncertainty1Pt,phi);

      fPtRel1PtUncertaintyChi2PerClusTPC->Fill(pt,relUncertainty1Pt,chi2PerClusterTPC);
      fPtNPointTPCSChi2PerClusTPC->Fill(pt,track->GetTPCnclsS()/nClustersTPC,chi2PerClusterTPC);
      fPtNPointTPCSRel1PtUncertainty->Fill(pt,track->GetTPCnclsS()/nClustersTPC,relUncertainty1Pt);
      if(track->GetTPCnclsS()/nClustersTPC>0.05) fPtRel1PtUncertaintyChi2PerClusTPCSharedSel->Fill(pt,relUncertainty1Pt,chi2PerClusterTPC);
    }//fTrackCutsITS loop
      
  }//ESD track loop
   
  // Post output data
  PostData(0, fHistList);
  PostData(1, fHistListTPC);
  PostData(2, fHistListITS);
  PostData(3, fHistListCosmics);

}
//________________________________________________________________________
Bool_t AliPWG4HighPtQATPConly::IsCosmic(const AliESDtrack *track1 , Int_t trackNumber, Double_t ptMin)
{
  Bool_t candidate1 = kFALSE;
  Bool_t candidate2 = kFALSE;
  if(!track1) return candidate1;

  Int_t nTracks = fESD->GetNumberOfTracks();
 
  for (Int_t iTrack = trackNumber+1; iTrack < nTracks; iTrack++) {
    candidate2 = kFALSE;
    AliESDtrack *track2 = fESD->GetTrack(iTrack);
    if(!track2) continue;
    if(!(fTrackCuts->AcceptTrack(track2))) continue;
    if(track2->Pt()<ptMin) continue;
    
    //Check if same charge. If not same charge, pair is cosmic candidate
    //Removed condition for time being. Not so clear how well we can measure curvature at high momenta
    // if( (track1->GetSign()*track2->GetSign()) > 0. ) continue;
    
    //Check if back-to-back
    Double_t mom1[3],mom2[3];
    track1->GetPxPyPz(mom1);
    track2->GetPxPyPz(mom2);
    Double_t cosTheta = (mom1[0]*mom2[0]+mom1[1]*mom2[1]+mom1[2]*mom2[2])/( TMath::Sqrt(mom1[0]*mom1[0]+mom1[1]*mom1[1]+mom1[2]*mom1[2])*TMath::Sqrt(mom2[0]*mom2[0]+mom2[1]*mom2[1]+mom2[2]*mom2[2]) );
    //Double_t theta = TMath::ACos(cosTheta);
    //if(TMath::Abs(TMath::Pi()-theta)<fMaxCosmicAngle) { candidate1 = kTRUE; candidate2 = kTRUE;}
    Double_t cosMaxCosmicAngle[2] = {TMath::Cos(TMath::Pi()-fMaxCosmicAngle),TMath::Cos(TMath::Pi()+fMaxCosmicAngle)};
    if(cosTheta >= cosMaxCosmicAngle[0] && cosTheta <= cosMaxCosmicAngle[1]) { candidate1 = kTRUE; candidate2 = kTRUE;}
    if(candidate2) {
      fDeltaPtCosmicCandidates->Fill(track1->Pt()-track2->Pt());
      fDeltaPhi->Fill(track1->Phi()-track2->Phi());
      fDeltaEta->Fill(track1->Eta()-track2->Eta());
    }

  }

  if(candidate1) {
    fPtCosmicCandidates->Fill(track1->Pt());
  }

   return candidate1;
}

//________________________________________________________________________
void AliPWG4HighPtQATPConly::Terminate(Option_t *)
{

}

#endif
