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

// An analysis task to check the MUON data in simulated data
//
//*-- Ivana Hrivnacova
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h> 
#include <TH1F.h>
#include <TROOT.h>

#include "AliMUONQATask.h" 
#include "AliESD.h" 
#include "AliLog.h"
#include "AliESDVertex.h" 
#include "AliESDMuonTrack.h"

//______________________________________________________________________________
AliMUONQATask::AliMUONQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fnTrackTrig(0), 
  ftracktot(0),
  fnevents(0),
  fSPLowpt(0),
  fSPHighpt(0),
  fSPAllpt(0),
  fSMLowpt(0),
  fSMHighpt(0),
  fSMAllpt(0),
  fSULowpt(0),
  fSUHighpt(0),
  fSUAllpt(0),
  fUSLowpt(0),
  fUSHighpt(0),
  fUSAllpt(0), 
  fLSLowpt(0),
  fLSHighpt(0),
  fLSAllpt(0),
  fhMUONVertex(0),
  fhMUONMult(0)
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
  fhMUONVertex = new TH1F("hMUONVertex","ITS Vertex"                ,100, -25., 25.);
  fhMUONMult   = new TH1F("hMUONMult"  ,"Multiplicity of ESD tracks",10,  -0.5, 9.5);

  
  // create output container
  
  fOutputContainer = new TObjArray(2) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhMUONVertex,             0) ; 
  fOutputContainer->AddAt(fhMUONMult,               1) ; 
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

  if (trigWord & 0x01) 
    fSPLowpt++;
  if (trigWord & 0x02)
    fSPHighpt++;
  if (trigWord & 0x04)
    fSPAllpt++;
  if (trigWord & 0x08)
    fSMLowpt++;
  if (trigWord & 0x010)
    fSMHighpt++;
  if (trigWord & 0x020)
    fSMAllpt++;
  if (trigWord & 0x040)
    fSULowpt++;
  if (trigWord & 0x080)
    fSUHighpt++; 
  if (trigWord & 0x100)
    fSUAllpt++;
  if (trigWord & 0x200)
    fUSLowpt++;     
  if (trigWord & 0x400)
    fUSHighpt++;
  if (trigWord & 0x800)
    fUSAllpt++;
  if (trigWord & 0x1000)
    fLSLowpt++;
  if (trigWord & 0x2000)
    fLSHighpt++;
  if (trigWord & 0x4000)
    fLSAllpt++;

  Int_t tracktrig  = 0 ;
  Int_t iTrack1 ; 
  
  for (iTrack1 = 0 ; iTrack1 < nTracks ; iTrack1++) { //1st loop
    AliESDMuonTrack* muonTrack = fESD->GetMuonTrack(iTrack1) ;
    ftracktot++ ;
    if(muonTrack->GetMatchTrigger()) {
      fnTrackTrig++ ;
      tracktrig++ ;
    }
  }

  fhMUONVertex->Fill(zVertex) ;
  fhMUONMult->Fill(Float_t(nTracks)) ;

  PostData(0, fOutputContainer);  
}

//______________________________________________________________________________
void AliMUONQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended

  AliInfo(Form("Terminate %s:", GetName())) ;
  fOutputContainer = (TObjArray*)GetOutputData(0);
  fhMUONVertex = (TH1F*)fOutputContainer->At(0);
  fhMUONMult   = (TH1F*)fOutputContainer->At(1); 
  
  Int_t eff_match = -1 ; 
  if (ftracktot) 
    eff_match = 100 * fnTrackTrig / ftracktot ;

  printf("===================================================\n") ;
  printf("================  %s ESD SUMMARY    ==============\n", GetName()) ;
  printf("                                                   \n") ;
  printf("         Total number of processed events  %d      \n", fnevents) ;
  printf("\n")  ;
  printf("\n")  ;
  printf("Table 4:                                         \n") ;
  printf(" Global Trigger output       Low pt  High pt   All\n") ;
  printf(" number of Single Plus      :\t");
  printf("%i\t%i\t%i\t", fSPLowpt, fSPHighpt, fSPAllpt) ;
  printf("\n");
  printf(" number of Single Minus     :\t");
  printf("%i\t%i\t%i\t", fSMLowpt, fSMHighpt, fSMAllpt) ;
  printf("\n");
  printf(" number of Single Undefined :\t"); 
  printf("%i\t%i\t%i\t", fSULowpt, fSUHighpt, fSUAllpt) ;
  printf("\n");
  printf(" number of UnlikeSign pair  :\t"); 
  printf("%i\t%i\t%i\t", fUSLowpt, fUSHighpt, fUSAllpt) ;
  printf("\n");
  printf(" number of LikeSign pair    :\t");  
  printf("%i\t%i\t%i\t", fLSLowpt, fLSHighpt, fLSAllpt) ;
  printf("\n");
  printf("===================================================\n") ;
  printf("\n") ;
  printf("matching efficiency with the trigger for single tracks = %2d %% \n", eff_match);
  
  TCanvas * cMUON = new TCanvas("cMUON", "MUON ESD Test", 400, 10, 600, 700) ;
  cMUON->Divide(1,2) ;
  cMUON->cd(1) ;
  fhMUONVertex->Draw() ;
  cMUON->cd(2) ;
  fhMUONMult->Draw() ;  
  cMUON->Print("MUON.eps") ; 

  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
 
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
}
