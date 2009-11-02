/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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

//--------------------------------------------------------------------------
//    Implementation of the beam interaction spot location and size estimate
//    via VertexerTracksNoConstraint resolution extrapolation
//    to the infinit multiplicity alias zero resolution
//
// Origin: AliAnalysisTaskVtXY
//         A. Jacholkowski, Catania 
//         adam.jacholkowski@cern.ch
//-----------------------------------------------------------------
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskVtXY.h"

#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
ClassImp(AliAnalysisTaskVtXY)

//________________________________________________________________________
AliAnalysisTaskVtXY::AliAnalysisTaskVtXY(const char *name) 
  : AliAnalysisTask(name, ""),
  fESD(0), 
  fList(0),
  fHistVtx(0),
  fHistVty(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(0, TList::Class());
  // Output slot #0 writes into a TList container
}

//________________________________________________________________________
void AliAnalysisTaskVtXY::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {
    // Disable all branches and enable only the needed ones
    // The next two lines are different when data produced as AliESDEvent is read
    //tree->SetBranchStatus("*", kFALSE);
    //tree->SetBranchStatus("Tracks.*", kTRUE);

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
      fESD = esdH->GetEvent();
  }
}

//________________________________________________________________________
void AliAnalysisTaskVtXY::CreateOutputObjects()
{
  // Create histograms
  // Called once

   fHistVtx = new TProfile("fHistVtx","Xvert-RMS vs sigmaX(NC)",100,0.,0.05,-1.,1.,"s");
   fHistVty = new TProfile("fHistVty","Yvert-RMS vs sigmaY(NC)",100,0.,0.05,-1.,1.,"s");
   fHistVtx->GetXaxis()->SetTitle("sigmaX(NContributors)");
   fHistVtx->GetYaxis()->SetTitle("Xv-av&spread");
   fHistVty->GetXaxis()->SetTitle("sigmaY(NContributors)");
   fHistVty->GetYaxis()->SetTitle("Yv-av&spread");
   fHistVtx->SetMarkerStyle(kFullCircle);
   fHistVty->SetMarkerStyle(kOpenCircle);
   //
   fList = new TList();
   fList->Add(fHistVtx);
   fList->Add(fHistVty);
}

//________________________________________________________________________
void AliAnalysisTaskVtXY::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  // 0 condition / existence of ESD
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  Int_t Ngood = 0;
  // Track loop kept in case of need
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    Ngood++;
  } //end track loop 
  // redo the ESD vertex withour constraint
  AliVertexerTracks vertexer(fESD->GetMagneticField());
  vertexer.SetITSMode();
  AliESDVertex *vertex = vertexer.FindPrimaryVertex(fESD);
  // ***************************************************
  //  const AliESDVertex *vtxESD = 0;
  // if(fESD) {vtxESD = fESD->GetPrimaryVertexTracks();}
  Double_t ESDvtx[3] ={999.,999.,9999.};
  Double_t XRes,YRes ;
  Double_t sig[3];
  // Cond 1.vertex has to exist 
   if (!vertex) {
     Printf("ERROR: vertex not available");
     return;
    }
   // Cond 2 - vertexer has not failed
   if(!vertex->GetStatus()){
     Printf("WARNING: vertexer has failed");
       return;
   }
     TString vtitle = vertex->GetTitle();
   // Cond 3 --> vertexer type check
     if(!vtitle.Contains("VertexerTracksNoConstraint")){
       Printf("WARNING: not VertexerTracksNoConstraint");
      return;
   }
   ULong64_t TrigMask=0;
   TrigMask = fESD->GetTriggerMask();
   // Cond.4 Trigger mask eventually - not activated
   if(TrigMask==0){
     Printf("Trigger Mask = 0");
     // return;
    }
    if(vertex){
    vertex->GetXYZ(ESDvtx);
    XRes = vertex->GetXRes();
    YRes = vertex->GetYRes();
    vertex->GetSigmaXYZ(sig);
    fHistVtx->Fill(XRes,ESDvtx[0]);
    fHistVty->Fill(YRes,ESDvtx[1]);
    PostData(0,fList);
    }
  //
}      

//________________________________________________________________________
void AliAnalysisTaskVtXY::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fList = dynamic_cast<TList*>(GetOutputData(0));
  if(!fList){
    Printf("ERROR: fList not available");
    return;
  }
  fHistVtx = (TProfile *)fList->At(0);
  fHistVty = (TProfile *)fList->At(1);
  TCanvas *c1 = new TCanvas("AliAnalysisTaskVtXY","Vtx analysis",10,10,800,800);
  c1->SetFillColor(10); c1->SetHighLightColor(10);
  c1->Divide(1,2);
  c1->cd(1);
  fHistVtx->DrawCopy("");
  c1->cd(2);
  fHistVty->DrawCopy("");
  // derive histos with fits from the above profiles
  Double_t spread = 0.;
  Double_t content = 0.;
  Double_t entries = 0.;
  Double_t error = 0.;
  TH1D * hsx = new TH1D("hsx","Xvtx spreads",100,0.,0.050);
  //fill
  for(Int_t i=0; i<100; i++){
    spread = fHistVtx->GetBinError(i);
    content = fHistVtx->GetBinContent(i);
    entries = fHistVtx->GetBinEntries(i);
    error = 0.;
    if(entries<10){content = 0.; spread = 0.;}
    if(entries>0){error = 0.0005 + spread/TMath::Sqrt(entries);}
    hsx->SetBinContent(i,spread);
    hsx->SetBinError(i,error);
  }
  hsx->GetXaxis()->SetTitle("sigma-Xvert");
  hsx->GetYaxis()->SetTitle("RMS(Xv) [cm]");
  hsx->GetYaxis()->SetTitleOffset(-0.3);
  hsx->GetYaxis()->SetRangeUser(0.,0.15);
  hsx->SetMarkerColor(3);
  hsx->SetMarkerStyle(20);
  /*TCanvas *cx = */new TCanvas("cx","",50,50,800,800);
  gStyle->SetOptFit(111);
  TPad *px = new TPad("px","",0,0,1,1);
  px->Draw();
  px->cd();
  px->SetFillColor(42);
  px->SetFrameFillColor(10);
  hsx->Fit("pol3","R","",0.001,0.02);
  hsx->Draw();
  // VtY
  TH1D * hsy = new TH1D("hsy","Yvtx spreads",100,0.,0.050);
  //fill
  for(Int_t i=0; i<100; i++){
    spread = fHistVty->GetBinError(i);
    content = fHistVty->GetBinContent(i);
    entries = fHistVty->GetBinEntries(i);
    error = 0.;
    if(entries<10){content = 0.; spread = 0.;}
    if(entries>0){error = 0.0005 + spread/TMath::Sqrt(entries);}
    hsy->SetBinContent(i,spread);
    hsy->SetBinError(i,error);
  }
  hsy->GetXaxis()->SetTitle("sigma-Yvert");
  hsy->GetYaxis()->SetTitle("RMS(Yv) [cm]");
  hsy->GetYaxis()->SetTitleOffset(-0.3);
  hsy->GetYaxis()->SetRangeUser(0.,0.15);
  hsy->SetMarkerColor(2);
  hsy->SetMarkerStyle(20);
  /*TCanvas *cy = */new TCanvas("cy","",100,100,800,800);
  TPad *py = new TPad("py","",0,0,1,1);
  py->Draw();
  py->cd();
  py->SetFillColor(42);
  py->SetFrameFillColor(10);
  hsy->Fit("pol3","R","",0.001,0.02);
  hsy->Draw();
}
