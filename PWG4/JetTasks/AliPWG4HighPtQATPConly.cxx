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

#ifndef ALIPWG4HighPtQATPCONLY_CXX
#define ALIPWG4HighPtQATPCONLY_CXX

#include "AliPWG4HighPtQATPConly.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TChain.h"
#include "TH3F.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"

using namespace std; //required for resolving the 'cout' symbol

ClassImp(AliPWG4HighPtQATPConly)

AliPWG4HighPtQATPConly::AliPWG4HighPtQATPConly(): AliAnalysisTask("AliPWG4HighPtQATPConly", ""), 
  fESD(0), 
  fTrackCuts(0), 
  fTrackCutsITS(0),
  fNEvent(0),
  fPtAll(0),
  fPtSel(0),
  fPtAllminPtTPCvsPtAll(0),
  fPtAllminPtTPCvsPtAllNPointTPC(0),
  fPtAllminPtTPCvsPtAllDCAR(0),
  fPtAllminPtTPCvsPtAllDCAZ(0),
  fPtAllminPtTPCvsPtAllPhi(0),
  fPtAllminPtTPCvsPtAllNPointITS(0),
  fPtAllminPtTPCvsPtAllNSigmaToVertex(0),
  fPtAllminPtTPCvsPtAllChi2C(0),
  fPtAllminPtTPCvsPtAllRel1PtUncertainty(0),
  fHistList(0),
  fPtAllTPC(0),
  fPtSelTPC(0),
  fPtSelTPCITS(0),
  fHistListTPC(0),
  fPtSelITS(0),
  fPtITSminPtTPCvsPtITS(0),
  fPtITSminPtTPCvsPtITSNPointTPC(0),
  fPtITSminPtTPCvsPtITSDCAR(0),
  fPtITSminPtTPCvsPtITSDCAZ(0),
  fPtITSminPtTPCvsPtITSPhi(0),
  fPtITSminPtTPCvsPtITSNPointITS(0),
  fPtITSminPtTPCvsPtITSNSigmaToVertex(0),
  fPtITSminPtTPCvsPtITSChi2C(0),
  fPtITSminPtTPCvsPtITSRel1PtUncertainty(0),
  fHistListITS(0)
{

}
//________________________________________________________________________
AliPWG4HighPtQATPConly::AliPWG4HighPtQATPConly(const char *name): 
  AliAnalysisTask(name, ""), 
  fESD(0),
  fTrackCuts(),
  fTrackCutsITS(),
  fNEvent(0),
  fPtAll(0),
  fPtSel(0),
  fPtAllminPtTPCvsPtAll(0),
  fPtAllminPtTPCvsPtAllNPointTPC(0),
  fPtAllminPtTPCvsPtAllDCAR(0),
  fPtAllminPtTPCvsPtAllDCAZ(0),
  fPtAllminPtTPCvsPtAllPhi(0),
  fPtAllminPtTPCvsPtAllNPointITS(0),
  fPtAllminPtTPCvsPtAllNSigmaToVertex(0),
  fPtAllminPtTPCvsPtAllChi2C(0),
  fPtAllminPtTPCvsPtAllRel1PtUncertainty(0),
  fHistList(0),
  fPtAllTPC(0),
  fPtSelTPC(0),
  fPtSelTPCITS(0),
  fHistListTPC(0),
  fPtSelITS(0),
  fPtITSminPtTPCvsPtITS(0),
  fPtITSminPtTPCvsPtITSNPointTPC(0),
  fPtITSminPtTPCvsPtITSDCAR(0),
  fPtITSminPtTPCvsPtITSDCAZ(0),
  fPtITSminPtTPCvsPtITSPhi(0),
  fPtITSminPtTPCvsPtITSNPointITS(0),
  fPtITSminPtTPCvsPtITSNSigmaToVertex(0),
  fPtITSminPtTPCvsPtITSChi2C(0),
  fPtITSminPtTPCvsPtITSRel1PtUncertainty(0),
  fHistListITS(0)
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
}

//________________________________________________________________________
void AliPWG4HighPtQATPConly::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  // Called once
  AliDebug(2,Form(">> AliPWG4HighPtSpectra::ConnectInputData \n"));
  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      AliDebug(2,Form("ERROR: Could not get ESDInputHandler")); 
    } else
      fESD = esdH->GetEvent();
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

  Int_t fgkNPhiBins=18;
  Float_t kMinPhi = 0.;
  Float_t kMaxPhi = 2.*TMath::Pi();
  
  Int_t fgkNPtBins=98;
  Float_t fgkPtMin=2.;
  Float_t fgkPtMax=100.;
  Int_t fgkResPtBins=40;

  fNEvent = new TH1F("fNEvent","NEvent",1,-0.5,0.5);
  fHistList->Add(fNEvent);
  fPtAll = new TH1F("fPtAll","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtAll);
  fPtSel = new TH1F("fPtSel","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSel);
  
  fPtAllminPtTPCvsPtAll = new TH2F("fPtAllminPtTPCvsPtAll","PtAllminPtTPCvsPtAll",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtAllminPtTPCvsPtAll->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAll->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fHistList->Add(fPtAllminPtTPCvsPtAll);
  
  fPtAllminPtTPCvsPtAllNPointTPC = new TH3F("fPtAllminPtTPCvsPtAllNPointTPC","PtAllminPtTPCvsPtAllNPointTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,160,0.5,160.5);
  fPtAllminPtTPCvsPtAllNPointTPC->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllNPointTPC->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllNPointTPC->SetZTitle("N_{point,TPC}");
  fHistList->Add(fPtAllminPtTPCvsPtAllNPointTPC);

  fPtAllminPtTPCvsPtAllDCAR = new TH3F("fPtAllminPtTPCvsPtAllDCAR","PtAllminPtTPCvsPtAllDCAR",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtAllminPtTPCvsPtAllDCAR->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllDCAR->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllDCAR->SetZTitle("DCA_{R}");
  fHistList->Add(fPtAllminPtTPCvsPtAllDCAR);

  fPtAllminPtTPCvsPtAllDCAZ = new TH3F("fPtAllminPtTPCvsPtAllDCAZ","PtAllminPtTPCvsPtAllDCAZ",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,-2.,2.);
  fPtAllminPtTPCvsPtAllDCAZ->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllDCAZ->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllDCAZ->SetZTitle("DCA_{Z}");
  fHistList->Add(fPtAllminPtTPCvsPtAllDCAZ);

  fPtAllminPtTPCvsPtAllPhi = new TH3F("fPtAllminPtTPCvsPtAllPhi","PtAllminPtTPCvsPtAllPhi",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtAllminPtTPCvsPtAllPhi->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllPhi->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllPhi->SetZTitle("#phi");
  fHistList->Add(fPtAllminPtTPCvsPtAllPhi);

  fPtAllminPtTPCvsPtAllNPointITS = new TH3F("fPtAllminPtTPCvsPtAllNPointITS","PtAllminPtTPCvsPtAllNPointITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,9,-0.5,8.5);
  fPtAllminPtTPCvsPtAllNPointITS->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllNPointITS->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllNPointITS->SetZTitle("N_{point,ITS}}");
  fHistList->Add(fPtAllminPtTPCvsPtAllNPointITS);
  
  fPtAllminPtTPCvsPtAllNSigmaToVertex = new TH3F("fPtAllminPtTPCvsPtAllNSigmaToVertex","PtAllminPtTPCvsPtAllNSigmaToVertex",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,0.,8.);
  fPtAllminPtTPCvsPtAllNSigmaToVertex->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllNSigmaToVertex->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistList->Add(fPtAllminPtTPCvsPtAllNSigmaToVertex);

  fPtAllminPtTPCvsPtAllChi2C = new TH3F("fPtAllminPtTPCvsPtAllChi2C","PtAllminPtTPCvsPtAllChi2C",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,0.,10.);
  fPtAllminPtTPCvsPtAllChi2C->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllChi2C->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllChi2C->SetZTitle("Constrained #chi^{2}");
  fHistList->Add(fPtAllminPtTPCvsPtAllChi2C);

  fPtAllminPtTPCvsPtAllRel1PtUncertainty = new TH3F("fPtAllminPtTPCvsPtAllRel1PtUncertainty","PtAllminPtTPCvsPtAllRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,30,0.,0.3);
  fPtAllminPtTPCvsPtAllRel1PtUncertainty->SetXTitle("p_{t}^{All}");
  fPtAllminPtTPCvsPtAllRel1PtUncertainty->SetYTitle("(1/p_{t}^{All}-1/p_{t}^{TPC})/(1/p_{t}^{All})");
  fPtAllminPtTPCvsPtAllRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistList->Add(fPtAllminPtTPCvsPtAllRel1PtUncertainty);

  //ITSrefit
  fPtSelITS = new TH1F("fPtSelITSrefit","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListITS->Add(fPtSelITS);
  
  fPtITSminPtTPCvsPtITS = new TH2F("fPtITSminPtTPCvsPtITS","PtITSminPtTPCvsPtITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.);
  fPtITSminPtTPCvsPtITS->SetXTitle("p_{t}^{ITS}");
  fPtITSminPtTPCvsPtITS->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{TPC})/(1/p_{t}^{ITS})");
  fHistListITS->Add(fPtITSminPtTPCvsPtITS);
  
  fPtITSminPtTPCvsPtITSNPointTPC = new TH3F("fPtITSminPtTPCvsPtITSNPointTPC","PtITSminPtTPCvsPtITSNPointTPC",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,160,0.5,160.5);
  fPtITSminPtTPCvsPtITSNPointTPC->SetXTitle("p_{t}^{ITSrefit}");
  fPtITSminPtTPCvsPtITSNPointTPC->SetYTitle("(1/p_{t}^{ITSrefit}-1/p_{t}^{TPC})/(1/p_{t}^{ITSrefit})");
  fPtITSminPtTPCvsPtITSNPointTPC->SetZTitle("N_{point,TPC}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSNPointTPC);
    
  fPtITSminPtTPCvsPtITSDCAR = new TH3F("fPtITSminPtTPCvsPtITSDCAR","PtITSminPtTPCvsPtITSDCAR",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,-1.,1.);
  fPtITSminPtTPCvsPtITSDCAR->SetXTitle("p_{t}^{ITSrefit}");
  fPtITSminPtTPCvsPtITSDCAR->SetYTitle("(1/p_{t}^{ITSrefit}-1/p_{t}^{TPC})/(1/p_{t}^{ITSrefit})");
  fPtITSminPtTPCvsPtITSDCAR->SetZTitle("DCA_{R}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSDCAR);
  
  fPtITSminPtTPCvsPtITSDCAZ = new TH3F("fPtITSminPtTPCvsPtITSDCAZ","PtITSminPtTPCvsPtITSDCAZ",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,-2.,2.);
  fPtITSminPtTPCvsPtITSDCAZ->SetXTitle("p_{t}^{ITSrefit}");
  fPtITSminPtTPCvsPtITSDCAZ->SetYTitle("(1/p_{t}^{ITSrefit}-1/p_{t}^{TPC})/(1/p_{t}^{ITSrefit})");
  fPtITSminPtTPCvsPtITSDCAZ->SetZTitle("DCA_{Z}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSDCAZ);
  
  fPtITSminPtTPCvsPtITSPhi = new TH3F("fPtITSminPtTPCvsPtITSPhi","PtITSminPtTPCvsPtITSPhi",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,fgkNPhiBins,kMinPhi,kMaxPhi);
  fPtITSminPtTPCvsPtITSPhi->SetXTitle("p_{t}^{ITSrefit}");
  fPtITSminPtTPCvsPtITSPhi->SetYTitle("(1/p_{t}^{ITSrefit}-1/p_{t}^{TPC})/(1/p_{t}^{ITSrefit})");
  fPtITSminPtTPCvsPtITSPhi->SetZTitle("#phi");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSPhi);
  
  fPtITSminPtTPCvsPtITSNPointITS = new TH3F("fPtITSminPtTPCvsPtITSNPointITS","PtITSminPtTPCvsPtITSNPointITS",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,9,-0.5,8.5);
  fPtITSminPtTPCvsPtITSNPointITS->SetXTitle("p_{t}^{ITSrefit}");
  fPtITSminPtTPCvsPtITSNPointITS->SetYTitle("(1/p_{t}^{ITSrefit}-1/p_{t}^{TPC})/(1/p_{t}^{ITSrefit})");
  fPtITSminPtTPCvsPtITSNPointITS->SetZTitle("N_{point,ITS}}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSNPointITS); 
  
  fPtITSminPtTPCvsPtITSNSigmaToVertex = new TH3F("fPtITSminPtTPCvsPtITSNSigmaToVertex","PtITSminPtTPCvsPtITSNSigmaToVertex",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,40,0.,8.);
  fPtITSminPtTPCvsPtITSNSigmaToVertex->SetXTitle("p_{t}^{ITS}");
  fPtITSminPtTPCvsPtITSNSigmaToVertex->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{TPC})/(1/p_{t}^{ITS})");
  fPtITSminPtTPCvsPtITSNSigmaToVertex->SetZTitle("N#sigma to vertex");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSNSigmaToVertex);

  fPtITSminPtTPCvsPtITSChi2C = new TH3F("fPtITSminPtTPCvsPtITSChi2C","PtITSminPtTPCvsPtITSChi2C",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,20,0.,10.);
  fPtITSminPtTPCvsPtITSChi2C->SetXTitle("p_{t}^{ITS}");
  fPtITSminPtTPCvsPtITSChi2C->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{TPC})/(1/p_{t}^{ITS})");
  fPtITSminPtTPCvsPtITSChi2C->SetZTitle("Constrained #chi^{2}");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSChi2C);

  fPtITSminPtTPCvsPtITSRel1PtUncertainty = new TH3F("fPtITSminPtTPCvsPtITSRel1PtUncertainty","PtITSminPtTPCvsPtITSRel1PtUncertainty",fgkNPtBins, fgkPtMin,fgkPtMax,fgkResPtBins,-1.,1.,30,0.,0.3);
  fPtITSminPtTPCvsPtITSRel1PtUncertainty->SetXTitle("p_{t}^{ITS}");
  fPtITSminPtTPCvsPtITSRel1PtUncertainty->SetYTitle("(1/p_{t}^{ITS}-1/p_{t}^{TPC})/(1/p_{t}^{ITS})");
  fPtITSminPtTPCvsPtITSRel1PtUncertainty->SetZTitle("Rel1PtUncertainty");
  fHistListITS->Add(fPtITSminPtTPCvsPtITSRel1PtUncertainty);

  fPtAllTPC = new TH1F("fPtAllTPC","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListTPC->Add(fPtAllTPC);
  fPtSelTPC = new TH1F("fPtSelTPC","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListTPC->Add(fPtSelTPC);
  fPtSelTPCITS = new TH1F("fPtSelTPCITS","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistListTPC->Add(fPtSelTPCITS);

  TH1::AddDirectory(oldStatus);   

}
//________________________________________________________________________
void AliPWG4HighPtQATPConly::Exec(Option_t *) {  
  // Main loop
  // Called for each event
  AliDebug(2,Form(">> AliPWG4HighPtQATPConly::Exec \n"));  

  if (!fESD) {
    AliDebug(2,Form("ERROR: fESD not available"));
    return;
  }

  const AliESDVertex *vtx = fESD->GetPrimaryVertex();

  // Need vertex cut
  if (vtx->GetNContributors() < 2)
    return;

  AliDebug(2,Form("Vertex title %s, status %d, nCont %d\n",vtx->GetTitle(), vtx->GetStatus(), vtx->GetNContributors()));

  // Need to keep track of evts without vertex
  fNEvent->Fill(0.);

  if(!fESD->GetNumberOfTracks() || fESD->GetNumberOfTracks()<2) return;
  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2,Form("nTracks %d", nTracks));

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    
    AliESDtrack *track = fESD->GetTrack(iTrack);
    AliExternalTrackParam *trackTPC = (AliExternalTrackParam *)track->GetTPCInnerParam();
    AliESDtrack *trackTPConly = fTrackCuts->GetTPCOnlyTrack(fESD,iTrack); 
    if(!track || !trackTPC || !trackTPConly) continue;

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
    Float_t nSigmaToVertex = fTrackCuts->GetSigmaToVertex(track);// Calculates the number of sigma to the vertex for a track.
    Float_t chi2C = track->GetConstrainedChi2();
    // Float_t relUncertainty1Pt = TMath::Sqrt(extCov[14])*pt;
    Float_t relUncertainty1Pt = TMath::Sqrt(track->GetSigma1Pt2())*pt;

    fPtAll->Fill(pt);
    fPtAllTPC->Fill(ptTPC);
    
    if (fTrackCuts->AcceptTrack(track)) {

      fPtSel->Fill(pt);
      
      fPtSelTPC->Fill(ptTPC);
      fPtAllminPtTPCvsPtAll->Fill(pt,(1./pt-1./ptTPC)/(1./pt) );
      fPtAllminPtTPCvsPtAllNPointTPC->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->GetTPCNcls());
      fPtAllminPtTPCvsPtAllDCAR->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dca2D);
      fPtAllminPtTPCvsPtAllDCAZ->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dcaZ);
      fPtAllminPtTPCvsPtAllPhi->Fill(pt,(1./pt-1./ptTPC)/(1./pt),phi);
      fPtAllminPtTPCvsPtAllNPointITS->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nPointITS);
      fPtAllminPtTPCvsPtAllNSigmaToVertex->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nSigmaToVertex);
      fPtAllminPtTPCvsPtAllChi2C->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2C);
      fPtAllminPtTPCvsPtAllRel1PtUncertainty->Fill(pt,(1./pt-1./ptTPC)/(1./pt),relUncertainty1Pt);
    }//fTrackCuts selection
    
    
    //ITSrefit selection
    if (fTrackCutsITS->AcceptTrack(track)) {
      
      fPtSelITS->Fill(pt);
      fPtSelTPCITS->Fill(ptTPC);
      fPtITSminPtTPCvsPtITS->Fill(pt,(1./pt-1./ptTPC)/(1./pt) );
      fPtITSminPtTPCvsPtITSNPointTPC->Fill(pt,(1./pt-1./ptTPC)/(1./pt),track->GetTPCNcls());
      fPtITSminPtTPCvsPtITSDCAR->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dca2D);
      fPtITSminPtTPCvsPtITSDCAZ->Fill(pt,(1./pt-1./ptTPC)/(1./pt),dcaZ);
      fPtITSminPtTPCvsPtITSPhi->Fill(pt,(pt-ptTPC)/(pt),phi);
      fPtITSminPtTPCvsPtITSNPointITS->Fill(pt,(pt-ptTPC)/(pt),nPointITS);
      fPtITSminPtTPCvsPtITSNSigmaToVertex->Fill(pt,(1./pt-1./ptTPC)/(1./pt),nSigmaToVertex);
      fPtITSminPtTPCvsPtITSChi2C->Fill(pt,(1./pt-1./ptTPC)/(1./pt),chi2C);
      fPtITSminPtTPCvsPtITSRel1PtUncertainty->Fill(pt,(1./pt-1./ptTPC)/(1./pt),relUncertainty1Pt);
    }//fTrackCutsITS loop
      
  }//ESD track loop
   
  // Post output data
  PostData(0, fHistList);
  PostData(1, fHistListTPC);
  PostData(2, fHistListITS);

}
//________________________________________________________________________
void AliPWG4HighPtQATPConly::Terminate(Option_t *)
{

}

#endif
