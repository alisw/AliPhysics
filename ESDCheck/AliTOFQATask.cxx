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
//_________________________________________________________________________
// An analysis task to check the TOF data in simulated data
//
//*-- Silvia Arcelli
//////////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h> 

#include "AliTOFQATask.h" 
#include "AliESD.h" 
#inclued "AliESDtrack.h" 
#include "AliLog.h"

//______________________________________________________________________________
AliTOFQATask::AliTOFQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fhTOF(0),
  fhTOFEnergy(0),
  fhTOFDigits(0),
  fhTOFRecParticles(0),
  fhTOFPhotons(0),
  fhTOFInvariantMass(0),
  fhTOFDigitsEvent(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
void AliTOFQATask::Init(const Option_t*)
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
  fhTOFSector   = new TH1F("hSector",   " TOF TrackRefs, sector # ",  18, 0., 18.) ;
  fhTOFSectorM  = new TH1F("hSectorM",  " TOF Matched, sector # ",    18, 0., 18.) ;
  fhTOFSectorMF = new TH1F("hSectorMF", " TOF Matched G, sector # ",  18, 0., 18.) ;
  fhTOFSectorMG = new TH1F("hSectorMG", " TOF Matched F , sector # ", 18, 0., 18.) ;
  fhTOFprimP    = new TH1F("hprimP",    " TPC mom  tracks",           20, 0., 4.) ;
  fhTOFprimPpi  = new TH1F("hprimPpi",  " TPC mom  tracks",           20, 0., 4.) ;
  fhTOFprimPka  = new TH1F("hprimPka",  " TPC mom  tracks",           20, 0., 4.) ;
  fhTOFprimPpr  = new TH1F("hprimPpr",  " TPC mom  tracks",           20, 0., 4.) ;


   // Reaching TOF,prim
  fhTOFprimPTOF   = new TH1F("hprimPTOF",   " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpi = new TH1F("hprimPTOFpi", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFka = new TH1F("hprimPTOFka", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpr = new TH1F("hprimPTOFpr", " TPC mom  tracks", 20, 0., 4.) ;


   // Well matched,prim
  fhTOFprimPTOF3   = new TH1F("hprimPTOF3",   " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpi3 = new TH1F("hprimPTOFpi3", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFka3 = new TH1F("hprimPTOFka3", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpr3 = new TH1F("hprimPTOFpr3", " TPC mom  tracks", 20, 0., 4.) ;

   // bad matched,prim
  fhTOFprimPTOF4   = new TH1F("hprimPTOF4",   " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpi4 = new TH1F("hprimPTOFpi4", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFka4 = new TH1F("hprimPTOFka4", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpr4 = new TH1F("hprimPTOFpr4", " TPC mom  tracks", 20, 0., 4.) ;

   // matched,prim
  fhTOFprimPTOF34   = new TH1F("hprimPTOF34",   " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpi34 = new TH1F("hprimPTOFpi34", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFka34 = new TH1F("hprimPTOFka34", " TPC mom  tracks", 20, 0., 4.) ;
  fhTOFprimPTOFpr34 = new TH1F("hprimPTOFpr34", " TPC mom  tracks", 20,0 ., 4.) ;
  
  // create output container
  
  fOutputContainer = new TObjArray(24) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhTOFSector,            0) ; 
  fOutputContainer->AddAt(fhTOFSectorM,           1) ; 
  fOutputContainer->AddAt(fhTOFSectorMF,          2) ; 
  fOutputContainer->AddAt(fhTOFSectorMG,          3) ; 
  fOutputContainer->AddAt(fhTOFprimP,             4) ; 
  fOutputContainer->AddAt(fhTOFprimPpi,           5) ; 
  fOutputContainer->AddAt(fhTOFprimPka,           6) ; 
  fOutputContainer->AddAt(fhTOFprimPpr,           7) ; 
  fOutputContainer->AddAt(fhTOFprimPTOF,          8) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpi,        9) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFka,       10) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpr,       11) ; 
  fOutputContainer->AddAt(fhTOFprimPTOF3,        12) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpi3,      13) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFka3,      14) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpr3,      15) ; 
  fOutputContainer->AddAt(fhTOFprimPTOF4,        16) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpi4,      17) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFka4,      18) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpr4,      19) ; 
  fOutputContainer->AddAt(fhTOFprimPTOF34,       20) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpi34,     21) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFka34,     22) ; 
  fOutputContainer->AddAt(fhTOFprimPTOFpr34,     23) ; 
 
}

//______________________________________________________________________________
void AliTOFQATask::Exec(Option_t *) 
{
  // Processing of one event
    
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld",  (dynamic_cast<TChain *>(fChain))->GetFile()->GetName(), entry)) ; 
  
  // ************************  TOF *************************************
  const Int_t knCalinSec = 8736 ;
  
  Int_t ntrk = fESD->GetNumberOfTracks() ;
  while ( ntrk-- ) {
    AliESDtrack * t = fESD->GetTrack(ntrk) ;
    if ( (t->GetStatus() & AliESDtrack::kTIME)==0 )
      continue;
    Int_t label               = TMath::Abs(t->GetLabel()) ;
    Double_t p                = t->GetP() ; 
    UInt_t assignedTOFcluster = t->GetTOFcluster() ;    //index of the assigned TOF cluster, >0 ?
    Int_t detid               = t->GetTOFCalChannel() ; //index of the assigned TOF cluster, >0 ?
    
    Int_t sector = detid / knCalinSec ;
    
    if(assignedTOFcluster){ //matched
      hSectorM->Fill(sector);
    }
  }
 
  PostData(0, fOutputContainer);

  
}

//______________________________________________________________________________
void AliTOFQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
  
  // some plots

  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
 
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
}
