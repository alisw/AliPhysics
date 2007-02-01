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

/* $Id$ */

//_________________________________________________________________________
// A test analysis task to check the pt of tracks distribution in simulated data
//
//*-- Panos
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>

#include "AliAnalysisTaskPt.h"
#include "AliESD.h"
#include "AliLog.h"

//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt(const char *name) :
  AliAnalysisTask(name,""),  
  fESD(0), 
  fhPt(0) 
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TObjArray::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPt::Init(Option_t *) 
{
  // Initialisation of branch container and histograms 
  
  AliInfo(Form("*** Initialization of %s", GetName())) ; 
  
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }
  
  if (!fESD) {
    // One should first check if the branch address was taken by some other task
    char ** address = (char **)GetBranchAddress(0, "ESD") ;
    if (address) 
      fESD = (AliESD *)(*address) ; 
    if (!fESD) 
      fChain->SetBranchAddress("ESD", &fESD) ;  
  }
  // The output objects will be written to 
  TDirectory * cdir = gDirectory ; 
  // Open a file for output #0
  char outputName[1024] ; 
  sprintf(outputName, "%s.root", GetName() ) ; 
  OpenFile(0, outputName , "RECREATE") ; 
  if (cdir) 
    cdir->cd() ; 
  
  // create histograms 

  fhPt = new TH1F("fhPt","This is the Pt distribution",15,0.1,3.1);
  fhPt->SetStats(kTRUE);
  fhPt->GetXaxis()->SetTitle("P_{T} [GeV]");
  fhPt->GetYaxis()->SetTitle("#frac{dN}{dP_{T}}");
  fhPt->GetXaxis()->SetTitleColor(1);
  fhPt->SetMarkerStyle(kFullCircle);
 
  // create output container
  
  fOutputContainer = new TObjArray(1) ; 
  fOutputContainer->SetName(GetName()) ;
  
  fOutputContainer->AddAt(fhPt, 0) ;
}

//________________________________________________________________________
void AliAnalysisTaskPt::Exec(Option_t *) 
{
  // Processing of one event
    
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  //************************  Pt tracks  *************************************
 
  for(Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack * ESDTrack = fESD->GetTrack(iTracks);
    Double_t momentum[3];
    ESDTrack->GetPxPyPz(momentum);
    Double_t Pt = sqrt(pow(momentum[0],2) + pow(momentum[1],2));
    fhPt->Fill(Pt);
  }//track loop 

  PostData(0, fOutputContainer);
}      

//________________________________________________________________________
void AliAnalysisTaskPt::Terminate(Option_t *) 
{
  // Processing when the event loop is ended
  
  TCanvas *c1 = new TCanvas("c1","Pt",10,10,310,310);
  c1->SetFillColor(10);
  c1->SetHighLightColor(10);
  
  c1->cd(1)->SetLeftMargin(0.15);
  c1->cd(1)->SetBottomMargin(0.15);  
  c1->cd(1)->SetLogy();
  fhPt->DrawCopy("E");

  c1->Print("TracksPt.eps");

  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);

  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
}
