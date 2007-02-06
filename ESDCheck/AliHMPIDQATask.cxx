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
// An analysis task to check the HMPID data in simulated data
//
//*-- Annalisa Mastroserio
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h> 
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h> 
#include <TROOT.h>
#include <TVector3.h> 

#include "AliHMPIDQATask.h" 
#include "AliESD.h" 
#include "AliLog.h"
#include "AliPID.h"

//______________________________________________________________________________
AliHMPIDQATask::AliHMPIDQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fhHMPIDCkovP(0),
  fhHMPIDMipXY(0),
  fhHMPIDDifXY(0),
  fhHMPIDSigP(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 

  Int_t i ; 
  for(i = 0 ; i < 5 ; i++) 
    fhHMPIDProb[i]=0;
}

//______________________________________________________________________________
AliHMPIDQATask::~AliHMPIDQATask()
{
  // dtor
  fOutputContainer->Clear() ; 
  delete fOutputContainer ; 
  
  delete fhHMPIDCkovP ;  
  delete fhHMPIDMipXY ;  
  delete fhHMPIDDifXY ;  
  delete fhHMPIDSigP ;   
  delete [] fhHMPIDProb ;
}

//______________________________________________________________________________
void AliHMPIDQATask::ConnectInputData(const Option_t*)
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
void AliHMPIDQATask::CreateOutputObjects()
{  
  // create histograms 
  fhHMPIDCkovP    = new TH2F("CkovP" , "#theta_{c}, [rad];P, [GeV]", 150,   0,  7  ,100, -3, 1); 
  fhHMPIDSigP     = new TH2F("SigP"  ,"#sigma_{#theta_c}"          , 150,   0,  7  ,100, 0, 1e20);
  fhHMPIDMipXY    = new TH2F("MipXY" ,"mip position"               , 260,   0,130  ,252,0,126); 
  fhHMPIDDifXY    = new TH2F("DifXY" ,"diff"                       , 260, -10, 10  ,252,-10,10); 
  
  fhHMPIDProb[0] = new TH1F("PidE" ,"PID: e yellow #mu magenta"  ,100,0,1); 
  fhHMPIDProb[0]->SetLineColor(kYellow);
  fhHMPIDProb[1] = new TH1F("PidMu","pid of #mu"                 ,100,0,1); 
  fhHMPIDProb[1]->SetLineColor(kMagenta);
  fhHMPIDProb[2] = new TH1F("PidPi","PID: #pi red K green p blue",100,0,1); 
  fhHMPIDProb[2]->SetLineColor(kRed);
  fhHMPIDProb[3] = new TH1F("PidK" ,"pid of K"                   ,100,0,1); 
  fhHMPIDProb[3]->SetLineColor(kGreen);
  fhHMPIDProb[4] = new TH1F("PidP" ,"pid of p"                   ,100,0,1); 
  fhHMPIDProb[4]->SetLineColor(kBlue);
 

  
  // create output container
  
  fOutputContainer = new TObjArray(9) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhHMPIDCkovP,      0) ; 
  fOutputContainer->AddAt(fhHMPIDSigP,       1) ; 
  fOutputContainer->AddAt(fhHMPIDMipXY,      2) ; 
  fOutputContainer->AddAt(fhHMPIDDifXY,      3) ; 
  fOutputContainer->AddAt(fhHMPIDProb[0],    4) ; 
  fOutputContainer->AddAt(fhHMPIDProb[1],    5) ; 
  fOutputContainer->AddAt(fhHMPIDProb[2],    6) ; 
  fOutputContainer->AddAt(fhHMPIDProb[3],    7) ; 
  fOutputContainer->AddAt(fhHMPIDProb[4],    8) ; 
}

//______________________________________________________________________________
void AliHMPIDQATask::Exec(Option_t *) 
{
  // Processing of one event
    
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // ************************  HMPID *************************************
  Int_t iTrk ; 
  for(iTrk = 0 ; iTrk < fESD->GetNumberOfTracks() ; iTrk++){
    AliESDtrack *pTrk = fESD->GetTrack(iTrk) ;

    fhHMPIDCkovP->Fill( pTrk->GetP(), pTrk->GetHMPIDsignal() ) ; 
    fhHMPIDSigP ->Fill( pTrk->GetP(), TMath::Sqrt(pTrk->GetHMPIDchi2()) ) ;
     
//     Float_t xm,ym; Int_t q,np;  pTrk->GetHMPIDmip(xm,ym,q,np);  fMipXY->Fill(xm,ym); //mip info
//     Float_t xd,yd,th,ph;        pTrk->GetHMPIDtrk(xd,yd,th,ph); fDifXY->Fill(xd,yd); //track info 
     
    Double_t pid[5] ;  
    pTrk->GetHMPIDpid(pid) ; 
    Int_t i ; 
    for(i = 0 ; i < 5 ; i++) 
      fhHMPIDProb[i]->Fill(pid[i]) ;
  }//tracks loop 
       
  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliHMPIDQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  fOutputContainer = (TObjArray*)GetOutputData(0);
  fhHMPIDCkovP   = (TH2F*)fOutputContainer->At(0);
  fhHMPIDSigP    = (TH2F*)fOutputContainer->At(1);
  fhHMPIDMipXY   = (TH2F*)fOutputContainer->At(2);
  fhHMPIDDifXY   = (TH2F*)fOutputContainer->At(3);
  fhHMPIDProb[0] = (TH1F*)fOutputContainer->At(4);
  fhHMPIDProb[1] = (TH1F*)fOutputContainer->At(5);
  fhHMPIDProb[2] = (TH1F*)fOutputContainer->At(6);
  fhHMPIDProb[3] = (TH1F*)fOutputContainer->At(7);
  fhHMPIDProb[4] = (TH1F*)fOutputContainer->At(8);
  
  Float_t n = 1.292 ; //mean freon ref idx 
  TF1 * hHMPIDpPi = new TF1("RiPiTheo", "acos(sqrt(x*x+[0]*[0])/(x*[1]))", 1.2, 7) ; 
  hHMPIDpPi->SetLineWidth(1) ; 
  hHMPIDpPi->SetParameter(1,n) ; 

  AliPID ppp ;                 
  hHMPIDpPi->SetLineColor(kRed);   
  hHMPIDpPi->SetParameter(0,AliPID::ParticleMass(AliPID::kPion));    //mass

  TF1 * hHMPIDK = static_cast<TF1*>(hHMPIDpPi->Clone()) ; 
  hHMPIDK ->SetLineColor(kGreen) ; 
  hHMPIDK ->SetParameter(0, AliPID::ParticleMass(AliPID::kKaon)) ; 

  TF1 * hHMPIDP=static_cast<TF1*>(hHMPIDpPi->Clone()) ; 
  hHMPIDP ->SetLineColor(kBlue) ;  
  hHMPIDP ->SetParameter(0,AliPID::ParticleMass(AliPID::kProton)) ; 

  TCanvas * cHMPID = new TCanvas("cHMPID","HMPID ESD Test") ;
  cHMPID->SetFillColor(10) ; 
  cHMPID->SetHighLightColor(10) ; 
  cHMPID->Divide(3,2) ;

  cHMPID->cd(1); 
  fhHMPIDCkovP->Draw() ; 
  hHMPIDpPi->Draw("same") ; 
  hHMPIDK->Draw("same") ; 
  hHMPIDP->Draw("same") ;   

  cHMPID->cd(2) ; 
  fhHMPIDMipXY->Draw() ;   

  cHMPID->cd(3) ; 
  fhHMPIDProb[0]->Draw() ; 
  fhHMPIDProb[1]->Draw("same") ;
  
  cHMPID->cd(4) ; 
  fhHMPIDSigP ->Draw() ;                                                          

  cHMPID->cd(5) ; 
  fhHMPIDDifXY->Draw() ;   

  cHMPID->cd(6) ; 
  fhHMPIDProb[2]->Draw() ; 
  fhHMPIDProb[3]->Draw("same") ; 
  fhHMPIDProb[4]->Draw("same") ; 

  cHMPID->Print("HMPID.eps");
  
  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
  
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
}
