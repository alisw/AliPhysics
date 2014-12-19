/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena.                                                 *
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
/*
  .
  .                  A template class to read tracks (PMD Cluster)
  .		         Runs in Local and Grid Modes
  .			Can be used for PbPb PMD analysis
  .		     Auther: Satyajit Jena <sjena@cern.ch>
  .                               12/04/2012  
  .
  . You are free to ask me anything you want about this code, I will add
  . a readme shortly to this folder about the running and adjusting the
  . code   
**************************************************************************/


#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDPmdTrack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliPMDAnalysisTaskPbPb.h"
#include "AliCentrality.h"


ClassImp(AliPMDAnalysisTaskPbPb)

//________________________________________________________________________
AliPMDAnalysisTaskPbPb::AliPMDAnalysisTaskPbPb(const char *name) 
: AliAnalysisTaskSE(name), fOutputList(0), fTrackCuts(0), fESD(0), fHistPt(0),fHistEta(0),
  fhEsdXYP(0),
  fhEsdXYC(0),
  fIsMC(kFALSE){

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliPMDAnalysisTaskPbPb::CreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner();

  //  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);

  // Track Cut Setup
  // Tune it accoringly
  fTrackCuts = new  AliESDtrackCuts("MyTrack Cuts","tracks");
  
  fTrackCuts->SetMinNClustersTPC(70);
  fTrackCuts->SetMaxChi2PerClusterTPC(4);
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetRequireITSRefit(kTRUE);
  fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fTrackCuts->SetMaxDCAToVertexZ(2);

  fTrackCuts->SetDCAToVertex2D(kFALSE);
  fTrackCuts->SetRequireSigmaToVertex(kFALSE);
    
  // Create the histograms as you like. 
  Int_t etabins = 40;
  Float_t etalow = -2.0, etaup = 2.0;
  fHistEta = new TH1F("fHistEta","#eta distribution for reconstructed",etabins, etalow, etaup);
  fHistEta->GetXaxis()->SetTitle("#eta");
  fHistEta->GetYaxis()->SetTitle("counts");
  fHistEta->SetMarkerStyle(kFullCircle);
  fHistEta->SetMarkerColor(kRed);
  fOutputList->Add(fHistEta);

  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 100, 0.1, 10.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  fHistPt->SetMarkerColor(kRed);
  fOutputList->Add(fHistPt);
  
  fhEsdXYC = new TH2F("hEsdXYC"," Scattered Plot for CPV (ESD) ",2000,-100.,100.,2000,-100.,100.);
  fhEsdXYC->GetXaxis()->SetTitle("X-axis");
  fhEsdXYC->GetYaxis()->SetTitle("Y-axis");
  fOutputList->Add(fhEsdXYC);
  
  fhEsdXYP = new TH2F("hEsdXYP"," Scattered Plot for PRE (ESD) ",2000,-100.,100.,2000,-100.,100.);
  fhEsdXYP->GetXaxis()->SetTitle("X-axis");
  fhEsdXYP->GetYaxis()->SetTitle("Y-axis");
  fOutputList->Add(fhEsdXYP);


  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliPMDAnalysisTaskPbPb::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if(!isSelected) return;

   
  if(!(fESD->GetPrimaryVertex()->GetStatus())) return;
  // if vertex is from spd vertexZ, require more stringent cut
  if (fESD->GetPrimaryVertex()->IsFromVertexerZ()) {
    if (fESD->GetPrimaryVertex()->GetDispersion()>0.02 ||  fESD->GetPrimaryVertex()->GetZRes()>0.25 ) return; // bad vertex from VertexerZ
  }
  

 AliCentrality *centrality = fESD->GetCentrality();

 /* If you want to make centrality bining use following */
 // Double_t nCentrality = -1;
 // nCentrality = centrality->GetCentralityPercentile("V0M");
 // printf(" %f \n", nCentrality);
 // if (nCentrality < 0) return;

 /* otherwise 
    if you want to cuttoff other centrality and keep a window 
    for example this is for 0 - 10% centrality with V0M estimator
*/
 if(!centrality->IsEventInCentralityClass(0,10,"V0M")) return;  
  
  
  // Track loop for reconstructed event - General purpose
  // any central analysis can be done.

  Int_t   nGoodTracks = 0;
  Int_t ntracks = fESD->GetNumberOfTracks();
  for(Int_t i = 0; i < ntracks; i++) {
    AliESDtrack* esdtrack = fESD->GetTrack(i); // pointer to reconstructed to track          
    if(!esdtrack) {
      AliError(Form("ERROR: Could not retrieve esdtrack %d",i));
      continue;
    }
    
    if(!fTrackCuts->AcceptTrack(esdtrack)) continue;
    nGoodTracks++;
    fHistPt->Fill(esdtrack->Pt());
    fHistEta->Fill(esdtrack->Eta());
  } 
  
   
  //_______________________________________________________
  // Only for PMD part
  
  Int_t npmdcl = fESD->GetNumberOfPmdTracks();
  printf ("Number of PMD Clusters: %8d \n",npmdcl) ;  

  for (Int_t i = 0; i < npmdcl; i++) {
    AliESDPmdTrack *pmdtr = fESD->GetPmdTrack(i);
    Int_t   det     = pmdtr->GetDetector();
    Float_t clsX    = pmdtr->GetClusterX();
    Float_t clsY    = pmdtr->GetClusterY();
    
    Float_t adc   = pmdtr->GetClusterADC();
    if(adc>1200)  continue;
        
    if ( det == 0 ) {
      fhEsdXYP->Fill(clsX,clsY);
    }
    else if ( det == 1 ) {
      fhEsdXYC->Fill(clsX,clsY);
    }
  }


  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliPMDAnalysisTaskPbPb::Terminate(Option_t *) {
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: fHistPt not available");
    return;
  }

  fHistPt = dynamic_cast<TH1F*>(fOutputList->FindObject("fHistPt"));
 TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,745,412);
 c1->SetFillColor(10);
 c1->SetBorderMode(0);
 c1->SetBorderSize(0);
 c1->SetFrameBorderMode(0);
 c1->cd(1)->SetLogy();
 fHistPt->DrawCopy("E");


  fHistEta = dynamic_cast<TH1F*>(fOutputList->FindObject("fHistEta"));
  TCanvas *c2 = new TCanvas("eta","Pt",10,10,745,412);
  c2->SetFillColor(10);
  c2->SetBorderMode(0);
  c2->SetBorderSize(0);
  c2->SetFrameBorderMode(0);
  c2->cd();
  fHistEta->DrawCopy("E");


  fhEsdXYP = dynamic_cast<TH2F*>(fOutputList->FindObject("hEsdXYP"));
  TCanvas *c3 = new TCanvas("XYP","Pt",10,10,745,745);
  c3->SetFillColor(10);
  c3->SetBorderMode(0);
  c3->SetBorderSize(0);
  c3->SetFrameBorderMode(0);
  c3->cd();
  fhEsdXYP->DrawCopy();


fhEsdXYC = dynamic_cast<TH2F*>(fOutputList->FindObject("hEsdXYC"));
  TCanvas *c4 = new TCanvas("XYC","Pt",10,10,745,745);
  c4->SetFillColor(10);
  c4->SetBorderMode(0);
  c4->SetBorderSize(0);
  c4->SetFrameBorderMode(0);
  c4->cd();
  fhEsdXYC->DrawCopy();

 

   
}

