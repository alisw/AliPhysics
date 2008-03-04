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

/// An analysis task to check the MUON data in simulated data
/// This class checks out the ESD tree, providing the matching with
/// the trigger,trigger responses for Low and High Pt cuts
/// (in Single, Unlike Sign and Like Sign) and gives Pt, Y, ITS vertex
/// and multiplicity distributions. All results are in histogram form.
/// The output is a root file and eps files in MUON.tar.gz. 

//*-- Frederic Yermia, yermia@to.infn.it
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TString.h> 

#include "AliMUONQATask.h" 
#include "AliESD.h" 
#include "AliLog.h"
#include "AliESDVertex.h" 
#include "AliESDMuonTrack.h"
#include <TLorentzVector.h>
//______________________________________________________________________________
AliMUONQATask::AliMUONQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fOutputContainer(0),
  fV1(),
  fnTrackTrig(0), 
  ftracktot(0),
  fnevents(0),
  fSLowpt(0),
  fUSLowpt(0),
  fUSHighpt(0),
  fLSLowpt(0),
  fLSHighpt(0),
  fmuonMass(0.105658389),
  fthetaX(0),
  fthetaY(0),
  fpYZ(0),
  fPxRec1(0),
  fPyRec1(0),
  fPzRec1(0),
  fE1(0),
  fZ1(0),
  fhMUONVertex(0),
  fhMUONMult(0),
  fhPt(0),
  fhY(0),
  fheffMatchT(0),
  fhSLowpt(0),
  fhUSLowpt(0),
  fhUSHighpt(0),
  fhLSLowpt(0),
  fhLSHighpt(0),
  fhChi2(0),  
  fhChi2match(0)  
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 

}


//______________________________________________________________________________
AliMUONQATask::~AliMUONQATask()
{ 
  // dtor
   fOutputContainer->Clear() ; 
   delete fOutputContainer ; 
  
   delete fhMUONVertex ; 
   delete fhMUONMult ; 
   delete fhPt ; 
   delete fhY ;
   delete fheffMatchT ;
   delete fhSLowpt ;
   delete fhUSLowpt ;
   delete fhUSHighpt;
   delete fhLSLowpt ;
   delete fhLSHighpt;
   delete fhChi2   ;  
   delete fhChi2match ;
}

//______________________________________________________________________________
void AliMUONQATask::ConnectInputData(const Option_t*)
{
  // Initialisation of branch container and histograms 
    
  AliInfo(Form("*** Initialization of %s", GetName())) ; 
  
  // Get input data
  fChain = dynamic_cast<TChain *>(GetInputData(0)) ;
  if (!fChain) {
    AliError(Form("Input 0 for %s not found\n", GetName()));
    return ;
  }
  
  // One should first check if the branch address was taken by some other task
  char ** address = (char **)GetBranchAddress(0, "ESD");
  if (address) {
    fESD = (AliESD*)(*address);
  } else {
    fESD = new AliESD();
    SetBranchAddress(0, "ESD", &fESD);
  }
}

//________________________________________________________________________
void AliMUONQATask::CreateOutputObjects()
{  
  // create histograms 

  OpenFile(0) ; 

  fhMUONVertex = new TH1F("hMUONVertex","ITS Vertex"                ,100, -25., 25.);
  fhMUONMult   = new TH1F("hMUONMult"  ,"Multiplicity of ESD tracks",10,  -0.5, 9.5);
  fhPt = new TH1F("hPt","Pt",100, 0.,20.);
  fhY = new TH1F("hY","Rapidity",100,-5.,-1.);
  fheffMatchT =  new TH1F("heff_matchT","Trigger Matching Efficiency",100, 0.,100.);
  fhSLowpt = new TH1F("hSLowpt","Single Low Pt Response (%)",101, 0.,101.);
  fhUSLowpt = new TH1F("hUSLowpt","Unlike Sign Low Pt Response (%)",101, 0.,101.);
  fhUSHighpt = new TH1F("hUSHighpt","Unlike Sign High Pt Response (%)",101, 0.,101.);
  fhLSLowpt = new TH1F("hLSLowpt","Like Sign Low Pt Response (%)",101, 0.,101.);
  fhLSHighpt = new TH1F("hLSHighpt","Like Sign High Pt Response (%)",101, 0.,101.);
  fhChi2 = new TH1F("hChi2","Chi2 by d.o.f.",100, 0.,20.);
  fhChi2match = new TH1F("hChi2match","Chi2 of trig/track matching",100, 0.,20.);
  // create output container
  
  fOutputContainer = new TObjArray(12) ; 
  fOutputContainer->SetName(GetName()) ; 
  fOutputContainer->AddAt(fhMUONVertex,             0) ; 
  fOutputContainer->AddAt(fhMUONMult,               1) ; 
  fOutputContainer->AddAt(fhPt,                     2) ; 
  fOutputContainer->AddAt(fhY,                      3) ; 
  fOutputContainer->AddAt(fheffMatchT,              4) ;
  fOutputContainer->AddAt(fhSLowpt,                 5) ;
  fOutputContainer->AddAt(fhUSLowpt,                6) ;
  fOutputContainer->AddAt(fhUSHighpt,               7) ;
  fOutputContainer->AddAt(fhLSLowpt,                8) ;
  fOutputContainer->AddAt(fhLSHighpt,               9) ;
  fOutputContainer->AddAt(fhChi2,                  10) ;
  fOutputContainer->AddAt( fhChi2match,            11) ;
}

//______________________________________________________________________________
void AliMUONQATask::Exec(Option_t *) 
{
  // Processing of one event
    
  fnevents++ ; 

  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // ************************  MUON *************************************
    
  const AliESDVertex* vertex = dynamic_cast<const AliESDVertex*>(fESD->GetVertex()) ;

  Double_t zVertex = 0. ;
  if (vertex) 
    zVertex = vertex->GetZv() ;
  
  Int_t nTracks = fESD->GetNumberOfMuonTracks() ;
  
  ULong64_t trigWord = fESD->GetTriggerMask() ;

  if (trigWord & 0x80) {
        fSLowpt++;
  }
  if (trigWord & 0x100){
    fLSLowpt++;
  } 
  if (trigWord & 0x200){
    fLSHighpt++;
  }  
  if (trigWord & 0x400){
    fUSLowpt++;
  }
  if (trigWord & 0x800){
    fUSHighpt++;
  }

  Int_t tracktrig  = 0 ;
  Int_t iTrack1 ; 
  
  for (iTrack1 = 0 ; iTrack1 < nTracks ; iTrack1++) { //1st loop
    AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack1) ;
    ftracktot++ ;
      fthetaX = muonTrack->GetThetaX();
      fthetaY = muonTrack->GetThetaY();
      fpYZ     =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
      fPzRec1  = - fpYZ / TMath::Sqrt(1.0 + TMath::Tan(fthetaY)*TMath::Tan(fthetaY));
      fPxRec1  = fPzRec1 * TMath::Tan(fthetaX);
      fPyRec1  = fPzRec1 * TMath::Tan(fthetaY);
      fZ1 = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));
      fE1 = TMath::Sqrt(fmuonMass * fmuonMass + fPxRec1 * fPxRec1 + fPyRec1 * fPyRec1 + fPzRec1 * fPzRec1);
      fV1.SetPxPyPzE(fPxRec1, fPyRec1, fPzRec1, fE1);

      // -----------> transverse momentum
      Float_t pt1 = fV1.Pt();
      // ----------->Rapidity
      Float_t y1 = fV1.Rapidity();
  
    if(muonTrack->GetMatchTrigger()) {
      fnTrackTrig++ ;
      tracktrig++ ;
      Float_t  Chi2match = muonTrack->GetChi2MatchTrigger();
      fhChi2match->Fill(Chi2match);
    }

     Float_t   fitfmin  = muonTrack->GetChi2();
     Int_t    ntrackhits = muonTrack->GetNHit();
     Float_t Chi2= fitfmin  / (2.0 * ntrackhits - 5);
    
     fhChi2->Fill(Chi2);
     fhPt->Fill(pt1);
     fhY->Fill(y1);
  }
  
  fhMUONVertex->Fill(zVertex) ;
  fhMUONMult->Fill(Float_t(nTracks)) ;
  
  PostData(0, fOutputContainer);  
}

//______________________________________________________________________________
void AliMUONQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
    Int_t fSLowPt = fSLowpt;
   if(fnevents){
     fSLowPt = 100 * fSLowpt / fnevents ;
     fhSLowpt->Fill(fSLowPt); }
    Int_t fUSLowPt = fUSLowpt;
   if(fnevents){
     fUSLowPt = 100 * fUSLowpt / fnevents ;
     fhUSLowpt->Fill(fUSLowPt); }
    Int_t fUSHighPt = fUSHighpt;
    if(fnevents){
      fUSHighPt = 100 * fUSHighpt / fnevents ;
      fhUSHighpt->Fill(fUSHighPt); }
    Int_t fLSLowPt = fLSLowpt;
    if(fnevents){
      fLSLowPt = 100 * fLSLowpt / fnevents ;
      fhLSLowpt->Fill(fLSLowPt); }
    Int_t fLSHighPt = fLSHighpt;
    if(fnevents){
      fLSHighPt = 100 * fLSHighpt / fnevents ;
      fhLSHighpt->Fill(fLSHighPt); }
  
    Int_t effMatch = -1 ; 
    if (ftracktot){ 
      effMatch = 100 * fnTrackTrig / ftracktot ;
      fheffMatchT->Fill(effMatch);}

    Bool_t problem = kFALSE ; 
    AliInfo(Form(" *** %s Report:", GetName())) ; 

    fOutputContainer = (TObjArray*)GetOutputData(0);
    fhMUONVertex = (TH1F*)fOutputContainer->At(0);
    fhMUONMult   = (TH1F*)fOutputContainer->At(1); 
    fhPt  = (TH1F*)fOutputContainer->At(2); 
    fhY  = (TH1F*)fOutputContainer->At(3); 
    fheffMatchT=(TH1F*)fOutputContainer->At(4); 
    fhSLowpt=(TH1F*)fOutputContainer->At(5); 
    fhUSLowpt=(TH1F*)fOutputContainer->At(6); 
    fhUSHighpt=(TH1F*)fOutputContainer->At(7);
    fhLSLowpt=(TH1F*)fOutputContainer->At(8); 
    fhLSHighpt=(TH1F*)fOutputContainer->At(9);
    
    printf("         Total number of processed events  %d      \n", fnevents) ;
    printf("     \n")  ;
    printf("     \n")  ;
    printf("     Table 1:                                         \n") ;
    printf("    ===================================================\n") ;
    printf("      Global Trigger output       Low pt  High pt \n") ;
    printf("     number of Single      :\t");
    printf("     %i\t", fSLowpt) ;
    printf("\n");
    printf("     number of UnlikeSign pair  :\t"); 
    printf("     %i\t%i\t", fUSLowpt, fUSHighpt) ;
    printf("\n");
    printf("     number of LikeSign pair    :\t");  
    printf("     %i\t%i\t", fLSLowpt, fLSHighpt) ;
    printf("\n");
    printf("     matching efficiency with the trigger for single tracks = %2d %% \n", effMatch);
    printf("\n") ;
    
    TCanvas * cMUON1 = new TCanvas("cMUON1", "MUON ESD Vert & Mult", 400, 10, 600, 700) ;
    cMUON1->Divide(1,2) ;
    cMUON1->cd(1) ;
    fhMUONVertex->SetXTitle("Vz (cm)");
    fhMUONVertex->Draw() ;
    cMUON1->cd(2) ;
    fhMUONMult->SetXTitle(" Track Multiplicity");
    fhMUONMult->Draw() ;
    cMUON1->Print("MUON1.eps") ; 
    
    TCanvas * cMUON2 = new TCanvas("cMUON2", "MUON ESD Pt & Y", 400, 10, 600, 700) ;
    cMUON2->Divide(1,2) ;
    cMUON2->cd(1) ;
    fhPt->SetXTitle("Pt (GeV)");
    fhPt->Draw() ;
    cMUON2->cd(2) ;
    fhY->SetXTitle("Y");
    fhY->Draw() ;
    cMUON2->Print("MUON2.eps") ;
    
    TCanvas * cMUON3 = new TCanvas("cMUON3", "Track Chi2 by dof and Chi2 of trigger/track matching ", 400, 10, 600, 700) ;
    cMUON3->Divide(1,2) ;
    cMUON3->cd(1) ;
    fhChi2->SetXTitle("Chi2 by d.o.f.");
    fhChi2->Draw();
    cMUON3->cd(2) ;
    fhChi2match->SetXTitle("Chi2 of trig/track matching");
    fhChi2match->Draw();
    cMUON3->Print("MUON3.eps") ;
    
    TCanvas * cMUON4 = new TCanvas("cMUON4", "Trigger Matching and Trigger Response (%)", 400, 10, 600, 700) ;
    cMUON4->Divide(2,3) ;
    cMUON4->cd(1) ;
    fheffMatchT->SetXTitle("%");
    fheffMatchT->Draw() ;
    cMUON4->cd(2) ;
    fhSLowpt->SetXTitle("%");
    fhSLowpt->Draw() ;
    cMUON4->cd(3) ;
    fhUSLowpt->SetXTitle("%");
    fhUSLowpt->Draw() ;
    cMUON4->cd(4) ;
    fhUSHighpt->SetXTitle("%");
    fhUSHighpt->Draw() ;
    cMUON4->cd(5) ;
    fhLSLowpt->SetXTitle("%");
    fhLSLowpt->Draw() ;
    cMUON4->cd(6) ;
    fhLSHighpt->SetXTitle("%");
    fhLSHighpt->Draw() ;
    cMUON4->Print("MUON4.eps") ;
    
    char line[1024] ; 
    sprintf(line, ".!tar -zcf %s.tar.gz *.eps", GetName()) ; 
    gROOT->ProcessLine(line);
    sprintf(line, ".!rm -fR *.eps"); 
    gROOT->ProcessLine(line);
    
    AliInfo(Form("!!! All the eps files are in %s.tar.gz !!!", GetName())) ;
 
   TString report ; 
   if(problem)
      report="Problems found, please check!!!";  
    else 
      report="OK";
    
   AliInfo(Form("*** %s Summary Report: %s \n",GetName(), report.Data())) ; 
}
