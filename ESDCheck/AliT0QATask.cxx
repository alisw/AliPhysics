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

//______________________________________________________________________________
// An analysis task to check the T0 data in simulated data
//
//*-- Alla Maevskaya
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h> 
#include <TH1F.h>
#include <TLegend.h> 
#include <TROOT.h>

#include "AliT0QATask.h" 
#include "AliESD.h" 
#include "AliLog.h"
#include "AliESDVertex.h" 

//______________________________________________________________________________
AliT0QATask::AliT0QATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fhT01(0),
  fhT02(0),
  fhT03(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
AliT0QATask::~AliT0QATask()
{
  // dtor
  fOutputContainer->Clear() ; 
  delete fOutputContainer ;
 
  delete fhT01 ; 
  delete fhT02 ;
  delete fhT03 ; 
}

//______________________________________________________________________________
void AliT0QATask::Init(const Option_t*)
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
  
  fhT01 = new TH1F("hRealVertex", "Primary vertex", 100,   -20,    20);
  fhT02 = new TH1F("hT0start",    "T0 start time",  100, 12400, 12600);
  fhT03 = new TH1F("hT0vertex",   "T0vertex",       100,   -20,    20);


  // create output container
  
  fOutputContainer = new TObjArray(3) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhT01,             0) ; 
  fOutputContainer->AddAt(fhT02,             1) ; 
  fOutputContainer->AddAt(fhT03,             2) ; 
}

//______________________________________________________________________________
void AliT0QATask::Exec(Option_t *) 
{
  // Processing of one event
    
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // ************************  T0 *************************************
  
  const AliESDVertex * vertex = fESD->GetPrimaryVertex() ; 
  Double_t posZ = vertex->GetZv() ; 
  fhT01->Fill( posZ ) ;

  fhT02->Fill( fESD->GetT0() ) ;
  
  fhT03->Fill( fESD->GetT0zVertex() / 2. ) ;
  
  PostData(0, fOutputContainer);
  
}

//______________________________________________________________________________
void AliT0QATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  
  Float_t mean = fhT02->GetMean();

  printf ("mean time T0 ps %f\n", mean) ;

  if ( mean > 12600 || mean < 12400 ) 
   AliWarning (" !!!!!!!!!!-----events sample is WRONG - T0 unreal -------");  

  TCanvas  * cTO1 = new TCanvas("cT01", "T0 ESD Test", 400, 10, 600, 700) ;
  cTO1->Divide(2, 2) ;

  cTO1->cd(1) ; 
  fhT01->Draw() ; 
    
  cTO1->cd(2) ; 
  fhT02->Draw() ;

  cTO1->cd(3) ; 
  fhT03->Draw() ; 


  cTO1->Print("T0.eps");

  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
 
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
}
