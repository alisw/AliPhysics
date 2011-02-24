
/* $Id$ */

#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TString.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDMuonTrack.h"

#include "AliAnalysisTaskSingleMuESD.h"

// analysis task for single muon analysis of the ESD events
// Authors: Bogdan Vulpescu, Nicole Bastid, LPC Clermont-Ferrand

ClassImp(AliAnalysisTaskSingleMuESD)

//________________________________________________________________________
AliAnalysisTaskSingleMuESD::AliAnalysisTaskSingleMuESD(const char *name) 
  : AliAnalysisTask(name, ""), 
    fESD(0), 
    fNtuple(0), 
    fTrigger(0), 
    fMaskTrig1MuL(1),
    fMaskTrig1MuH(2),
    fMaskTrigUSL(4),
    fMaskTrigLSL(16),
    fMaskTrigUSH(8),
    fMaskTrigLSH(32),
    fTriggerType(kFALSE)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TNtuple container
  DefineOutput(0, TNtuple::Class());
  // Output slot #1 writes into a TH1F container
  DefineOutput(1, TH1F::Class());
}

//________________________________________________________________________
void AliAnalysisTaskSingleMuESD::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    tree->SetBranchStatus("*", kFALSE);
    tree->SetBranchStatus("MuonTracks.*", kTRUE);
    tree->SetBranchStatus("SPDVertex.*", kTRUE);
    tree->SetBranchStatus("AliESDHeader.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskSingleMuESD::CreateOutputObjects()
{
  // Create output ntuple

  fNtuple = new TNtuple("ntEsd1mu","ntEsd1mu","VtZ:Mult:TMask:Z:Pt:P:Eta:Theta:Phi:Zpos:Hits:Xi2:Eff:TrigMatch:TrigXi2");
  // histogram with trigger classes
  fTrigger = new TH1F("TriggerMask","L/H/USL/USH/LSL/LSH",6,0,6);

}

//________________________________________________________________________
void AliAnalysisTaskSingleMuESD::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  // Fill the ntuple with variables and the trigger histogram

  Float_t  ntvar[15], pt, p, eta, thetad, phid, charge, chi2, zcoor, chi2MatchTrig;
  ULong64_t mask = 0;
  Int_t masksave = 0, hits, matchTrig;
  Double_t rEff[2];
  Double_t fZVertex = 0;
  Double_t fYVertex = 0;
  Double_t fXVertex = 0;
  Double_t fErXVertex = 0;
  Double_t fErYVertex = 0;
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  // Trigger checks
  // mask values: L/H/USL/USH/LSL/LSH
  // Bit position    Histogram bin     Trigger class
  // 1 (less)        1 ("L")           MULow          single low pt
  // 2               2 ("H")           MUHigh         single high pt
  // 3               3 ("USL")         MULU           unlike-sign low pt
  // 4               4 ("USH")         MUHU           unlike-sign high pt
  // 5               5 ("LSL")         MULL           like-sign low pt
  // 6 (most)        6 ("LSH")         MUHL           like-sign high pt
  if (fTriggerType) {
    mask = fESD->GetTriggerMask();
    if (mask & fMaskTrig1MuL) {
      fTrigger->Fill(0.5);
      masksave |= 0x01;
    }
    if (mask & fMaskTrig1MuH) {
      fTrigger->Fill(1.5);
      masksave |= 0x02;
    }
    if (mask & fMaskTrigUSL)  {
      fTrigger->Fill(2.5);
      masksave |= 0x04;
    }
    if (mask & fMaskTrigUSH)  {
      fTrigger->Fill(3.5);
      masksave |= 0x08;
    }
    if (mask & fMaskTrigLSL)  {
      fTrigger->Fill(4.5);
      masksave |= 0x10;
    }
    if (mask & fMaskTrigLSH)  {
      fTrigger->Fill(5.5);
      masksave |= 0x20;
    }
  } else {
    TString firedTC = fESD->GetFiredTriggerClasses();
    if (firedTC.Contains("MULow")) {
      fTrigger->Fill(0.5);
      masksave |= 0x01;
    }
    if (firedTC.Contains("MUHigh")) {
      fTrigger->Fill(1.5);
      masksave |= 0x02;
    }
    if (firedTC.Contains("MULU")) {
      fTrigger->Fill(2.5);
      masksave |= 0x04;
    }
    if (firedTC.Contains("MUHU")) {
      fTrigger->Fill(3.5);
      masksave |= 0x08;
    }
    if (firedTC.Contains("MULL")) {
      fTrigger->Fill(4.5);
      masksave |= 0x10;
    }
    if (firedTC.Contains("MUHL")) {
      fTrigger->Fill(5.5);
      masksave |= 0x20;
    }
  }

  // get the SPD reconstructed vertex (vertexer) 
  AliESDVertex* vertex = (AliESDVertex*) fESD->GetVertex();
  if (vertex->GetNContributors()) {
    fZVertex = vertex->GetZv();
    fYVertex = vertex->GetYv();
    fXVertex = vertex->GetXv();
    fErXVertex = vertex->GetXRes();
    fErYVertex = vertex->GetYRes();
  }

  Int_t nTracks = (Int_t)fESD->GetNumberOfMuonTracks() ;
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

    AliESDMuonTrack* muonTrack = new AliESDMuonTrack(*(fESD->GetMuonTrack(iTrack)));
    if (!muonTrack) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }

    pt = muonTrack->Pt();
    p = muonTrack->P();
    eta = muonTrack->Eta();
    thetad =  muonTrack->Theta()*180./TMath::Pi();	
    phid = muonTrack->Phi()*180./TMath::Pi();
    charge = muonTrack->Charge();
    hits = muonTrack->GetNHit();	
    chi2 = muonTrack->GetChi2() / (2.0 * hits - 5);
    zcoor = muonTrack->GetZ();
    matchTrig = muonTrack->GetMatchTrigger();
    chi2MatchTrig = muonTrack->GetChi2MatchTrigger();
    
    rEff[0] = rEff[1] = 0.0;
    GetEffFitted(eta,pt,rEff);
        
    ntvar[0] = fZVertex;
    ntvar[1] = nTracks;
    ntvar[2] = masksave;
    ntvar[3] = charge;
    ntvar[4] = pt;
    ntvar[5] = p;
    ntvar[6] = eta;	
    ntvar[7] = thetad;
    ntvar[8] = phid;
    ntvar[9] = zcoor;
    ntvar[10] = hits;	
    ntvar[11] = chi2;
    ntvar[12] = rEff[0];
    ntvar[13] = matchTrig;
    ntvar[14] = chi2MatchTrig;	

    fNtuple->Fill(ntvar);


  } //track loop 

  // Post output data.
  PostData(0, fNtuple);
  PostData(1, fTrigger);

}      

//________________________________________________________________________
void AliAnalysisTaskSingleMuESD::Terminate(Option_t *) 
{
  // the end
}

//________________________________________________________________________
void AliAnalysisTaskSingleMuESD::GetEffFitted(Double_t eta, Double_t pt, Double_t rEff[2])
{
  // tracking efficiency

  Int_t nbinEta =15;
  Float_t etaMin = -4.;
  Float_t etaBin = 0.1;
  Int_t numBin = -1;
  rEff[0]=1.;rEff[1]=0.;
  
  // nbinEta = number of eta bins startig at etaMin = -4 and binning of 0.1 (etaBin)
  Float_t p0[15] = {-0.9338,-0.9441,-0.9397,-0.9586,-0.9713,-0.9752,-0.9711,-0.9751,-0.9746,-0.9729,-0.9754,-0.9722,-0.9722,-0.971,-0.9709};
  Float_t p1[15] = {-0.6383,-1.4729,-1.0109,-1.4316,-1.5926,-1.367,-1.1895,-1.2834,-1.3289,-1.5916,-1.4258,-1.0983,-1.0812,-0.7179, -0.4613};      
  
  Float_t etaInf, etaSup;
  for (Int_t i = 1; i < nbinEta+1; i++){
    etaInf = etaMin + (i-1)*etaBin;
    etaSup = etaMin + i*etaBin;
    if(eta > etaInf && eta <= etaSup) numBin = i-1;
  } 
  
  if(numBin>=0 && pt >=1.) rEff[0] = p0[numBin] * TMath::TanH(p1[numBin]*pt);
  rEff[1] = 0;

}

//________________________________________________________________________
void AliAnalysisTaskSingleMuESD::SetTriggerType(const Char_t *trig)
{
  // set the trigger masks according to the trigger type used in the
  // simulations
  
  TString triggerType(trig);

  if (!triggerType.CompareTo("MUON")) {
    fMaskTrig1MuL =  1;  // MULow
    fMaskTrig1MuH =  2;  // MUHigh
    fMaskTrigUSL  =  4;  // MULU
    fMaskTrigLSL  = 16;  // MULL
    fMaskTrigUSH  =  8;  // MUHU
    fMaskTrigLSH  = 32;  // MUHL
  } else if (!triggerType.CompareTo("p-p")) {
    fMaskTrig1MuL =    128;  // MULow
    fMaskTrig1MuH =    256;  // MUHigh
    fMaskTrigUSL  =  65536;  // MULU
    fMaskTrigLSL  = 131072;  // MULL
    fMaskTrigUSH  = 262144;  // MUHU
    fMaskTrigLSH  = 524288;  // MUHL
  } else {
    // MUON trigger, values by default in the constructor
  }
  
  fTriggerType = kTRUE;

}

