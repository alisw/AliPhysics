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
// An analysis task to check the FMD data in simulated data
//
//*-- Hans Hjersing Dalsgaard 
//////////////////////////////////////////////////////////////////////////////

#include <TCanvas.h> 
#include <TChain.h>
#include <TF1.h> 
#include <TFile.h> 
#include <TH1D.h> 
#include <TROOT.h>

#include "AliFMDQATask.h" 
#include "AliESD.h" 
#include "AliLog.h"

//______________________________________________________________________________
AliFMDQATask::AliFMDQATask(const char *name) : 
  AliAnalysisTask(name,""),  
  fChain(0),
  fESD(0), 
  fhFMD1i(0),
  fhFMD2i(0), 
  fhFMD2o(0), 
  fhFMD3i(0), 
  fhFMD3o(0) 
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
AliFMDQATask::~AliFMDQATask()
{
  // dtor

  fOutputContainer->Clear() ; 
  delete fOutputContainer ; 

  delete fhFMD1i ;
  delete fhFMD2i ; 
  delete fhFMD2o ; 
  delete fhFMD3i ; 
  delete fhFMD3o ;

}

//______________________________________________________________________________
void AliFMDQATask::Init(const Option_t*)
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
  fhFMD1i = new TH1D("FMD1i", "FMD1i", 100, -0.5, 3);
  fhFMD2i = new TH1D("FMD2i", "FMD2i", 100, -0.5, 3);
  fhFMD2o = new TH1D("FMD2o", "FMD2o", 100, -0.5, 3);
  fhFMD3i = new TH1D("FMD3i", "FMD3i", 100, -0.5, 3);
  fhFMD3o = new TH1D("FMD3o", "FMD3o", 100, -0.5, 3);

  // create output container
  
  fOutputContainer = new TObjArray(5) ; 
  fOutputContainer->SetName(GetName()) ; 

  fOutputContainer->AddAt(fhFMD1i,             0) ; 
  fOutputContainer->AddAt(fhFMD2i,             1) ; 
  fOutputContainer->AddAt(fhFMD2o,             2) ; 
  fOutputContainer->AddAt(fhFMD3i,             3) ; 
  fOutputContainer->AddAt(fhFMD3o,             4) ; 
}

//______________________________________________________________________________
void AliFMDQATask::Exec(Option_t *) 
{
  // Processing of one event
    
  Long64_t entry = fChain->GetReadEntry() ;
  
  if (!fESD) {
    AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  TFile * currentFile = (dynamic_cast<TChain *>(fChain))->GetFile() ; 
  if ( !((entry-1)%100) ) 
    AliInfo(Form("%s ----> Processing event # %lld", currentFile->GetName(), entry)) ; 
	    
  // ************************  FMD *************************************
  
  AliESDFMD * fmd = fESD->GetFMDData() ;
  
  fmd->CheckNeedUShort(currentFile);

  Int_t nFMD1 = 0, nFMD2i = 0, nFMD2o = 0, nFMD3i = 0, nFMD3o = 0 ;
  
  UShort_t detector = 1 ;
  for(detector = 1 ; detector <= fmd->MaxDetectors() ; detector++) {
    Char_t ring = 'O' ;
    UShort_t sector ; 
    for(sector = 0 ;sector < fmd->MaxSectors() ; sector++) {
      UShort_t strip ;
      for(strip = 0 ; strip < fmd->MaxStrips(); strip++) {
	if(fmd->Multiplicity(detector, ring, sector, strip) != AliESDFMD::kInvalidMult)
	  RingSelector(detector, ring, fmd->Multiplicity(detector, ring, sector, strip)) ;
	if( (fmd->Multiplicity(detector, ring, sector, strip) == AliESDFMD::kInvalidMult) && detector == 2 )
	  nFMD2o++ ;
	if( (fmd->Multiplicity(detector, ring, sector, strip)==AliESDFMD::kInvalidMult) && detector == 3 )
	  nFMD3o++ ;
      }
    }
    ring='I';
    for(sector = 0; sector < fmd->MaxSectors(); sector++) {
      UShort_t strip ;
      for(strip = 0 ; strip < fmd->MaxStrips() ; strip++) {
	if(fmd->Multiplicity(detector, ring, sector, strip) != AliESDFMD::kInvalidMult)
	  RingSelector(detector, ring, fmd->Multiplicity(detector, ring, sector, strip));
	if( (fmd->Multiplicity(detector, ring, sector, strip) == AliESDFMD::kInvalidMult) && detector == 1 )
	  nFMD1++;
	if( (fmd->Multiplicity(detector, ring, sector, strip) == AliESDFMD::kInvalidMult) && detector == 2 )
	  nFMD2i++;
	if( (fmd->Multiplicity(detector, ring, sector, strip) == AliESDFMD::kInvalidMult) && detector == 3 )
	  nFMD3i++;
      }     
    }
  }

  if(nFMD1>100+10240)
    AliWarning(Form("number of missing strips in FMD1i too high in event number %lld in file", entry, fChain->GetCurrentFile()->GetName())) ;
  if(nFMD2i>100+10240)
    AliWarning(Form("number of missing strips in FMD2i too high in event number %lld in file", entry, fChain->GetCurrentFile()->GetName())) ;
  if(nFMD2o>100+10240)
    AliWarning(Form("number of missing strips in FMD2o too high in event number %lld in file", entry, fChain->GetCurrentFile()->GetName())) ;
  if(nFMD3i>100+10240)
    AliWarning(Form("number of missing strips in FMD3i too high in event number %lld in file", entry, fChain->GetCurrentFile()->GetName())) ;
  if(nFMD3o>100+10240)
    AliWarning(Form("number of missing strips in FMD3o too high in event number %lld in file", entry, fChain->GetCurrentFile()->GetName())) ;
  
  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliFMDQATask::Terminate(Option_t *)
{
  // Processing when the event loop is ended
 
  TCanvas * cFMD1 = new TCanvas("cFMD1", "FMD ESD Test", 400, 10, 600, 700);
  cFMD1->Divide(3, 2) ; 

  cFMD1->cd(1) ;; 
  fhFMD1i->Draw() ; 
  
  cFMD1->cd(2) ;; 
  fhFMD2i->Draw() ; 

  cFMD1->cd(3) ;; 
  fhFMD2o->Draw() ; 

  cFMD1->cd(4) ;; 
  fhFMD3i->Draw() ; 

  cFMD1->cd(5) ;; 
  fhFMD3o->Draw() ; 
  
  cFMD1->Print("FMD.eps") ;

  Bool_t rv1i = TestHisto(fhFMD1i) ;
  Bool_t rv2i = TestHisto(fhFMD2i) ;
  Bool_t rv2o = TestHisto(fhFMD2o) ;
  Bool_t rv3i = TestHisto(fhFMD2i) ;
  Bool_t rv3o = TestHisto(fhFMD2o) ;
  
  if (  !(rv1i * rv2i * rv2o * rv3i * rv3o) )
    AliWarning("Possible problem in file !!! Check output!") ;



  char line[1024] ; 
  sprintf(line, ".!tar -zcvf %s.tar.gz *.eps", GetName()) ; 
  gROOT->ProcessLine(line);
  sprintf(line, ".!rm -fR *.eps"); 
  gROOT->ProcessLine(line);
 
  AliInfo(Form("!!! All the eps files are in %s.tar.gz !!! \n", GetName())) ;
}

//______________________________________________________________________________
void AliFMDQATask::RingSelector(const UShort_t detector, const Char_t ring, const Float_t mult) const 
{ 
  // fill the histograms for each ring in each detector layer

  if(ring == 'I' && detector == 1)
    fhFMD1i->Fill(mult) ;
  if(ring == 'I' && detector == 2)
    fhFMD2i->Fill(mult) ;
  if(ring == 'O' && detector == 2)
    fhFMD2o ->Fill(mult) ;
  if(ring == 'I' && detector == 3)
    fhFMD3i ->Fill(mult) ;
  if(ring == 'O' && detector == 3)
    fhFMD3o ->Fill(mult) ;
}
//______________________________________________________________________________
Bool_t AliFMDQATask::TestHisto(TH1D * hTest) const 
{  
  // analyses the histogram with a Landau function

  Float_t chiMax = 3, chiLow=0.5 ;
  Float_t chiSq ;
  Float_t mpv ;
  Int_t   ndf ;
  
  FitAll(hTest, chiSq, ndf, mpv, chiMax, chiLow);

  if( (chiSq > chiMax || chiSq < chiLow) || mpv < 0.6 || mpv > 1 ) {
    hTest->Rebin(2) ;
    FitAll(hTest, chiSq, ndf, mpv, chiMax, chiLow) ;
  }

  Bool_t   ret   = kFALSE ;
  Char_t * test  = "not OK";
  Char_t * test2 = "not OK";
  
  if(chiSq < chiMax && chiSq > chiLow)
    test = "OK" ;
  if(mpv > 0.6 && mpv < 1)
    test2 = "OK" ;
 
  if(test == "OK" && test2 == "OK")
    ret = kTRUE;
  
  if(test == "not OK" || test2 == "not OK") {
    AliWarning("Bad fit results") ; 
    printf("Detector : %s\n", hTest->GetName()) ;
    printf("Landau fit Chi Square / NDF = %f / %d which is %s\n", chiSq*ndf, ndf, test) ; 
    printf("Landau fit MPV is: %f which is %s\n", mpv, test2) ;
  }
  return ret;  
}
//______________________________________________________________________________
void AliFMDQATask::FitAll(TH1D* hTest, Float_t &chiSq, Int_t &ndf, Float_t &mpv, Float_t chiMax, Float_t chiLow ) const 
{
  // fit a histogram with a Landau distribution and returns chi square

  Float_t fitMax = hTest->GetXaxis()->GetXmax() ;

  TH1D hTmp = *hTest ;
  hTmp.SetAxisRange(0.4,fitMax) ;
  Int_t   maxBin = hTmp.GetMaximumBin();
  Float_t max    = hTmp.GetBinCenter(maxBin);
 
  hTest->Fit("landau", "QOI", "", max-0.3, fitMax) ;
  TF1* fitfunc = hTest->GetFunction("landau") ;
  chiSq = fitfunc->GetChisquare() / fitfunc->GetNDF() ;
  mpv   = fitfunc->GetParameter(1) ;
  ndf   = fitfunc->GetNDF() ;

  if( ( chiSq > chiMax || chiSq < chiLow ) || ( mpv < 0.6 || mpv > 1 ) ) {
    hTest->Fit("landau", "QOI", "", max-0.2, fitMax) ;
    fitfunc = hTest->GetFunction("landau") ;
    chiSq   = fitfunc->GetChisquare() / fitfunc->GetNDF() ;
    mpv     = fitfunc->GetParameter(1) ;
    ndf     = fitfunc->GetNDF() ;
  }
  if( ( chiSq >chiMax || chiSq < chiLow ) || ( mpv < 0.6 || mpv > 1 ) ) {
    hTest->Fit("landau", "QOI", "", max-0.1, fitMax) ;
    fitfunc = hTest->GetFunction("landau") ;
    chiSq   = fitfunc->GetChisquare() / fitfunc->GetNDF() ;
    mpv     = fitfunc->GetParameter(1) ;
    ndf     = fitfunc->GetNDF() ;
  }
  if( ( chiSq > chiMax || chiSq <chiLow ) || ( mpv < 0.6 || mpv > 1 ) ) {
    hTest->Fit("landau", "QOI", "", max, fitMax) ;
    fitfunc = hTest->GetFunction("landau") ;
    chiSq   = fitfunc->GetChisquare() / fitfunc->GetNDF() ;
    mpv     = fitfunc->GetParameter(1) ;
    ndf     = fitfunc->GetNDF(); 
  }
  if( ( chiSq > chiMax || chiSq < chiLow ) || ( mpv < 0.6 ||mpv > 1 ) ) {
    hTest->Fit("landau", "QOI", "", max+0.1, fitMax) ;
    fitfunc = hTest->GetFunction("landau") ;
    chiSq   = fitfunc->GetChisquare() / fitfunc->GetNDF() ;
    mpv     = fitfunc->GetParameter(1) ;
    ndf     = fitfunc->GetNDF() ;
  } 
}
