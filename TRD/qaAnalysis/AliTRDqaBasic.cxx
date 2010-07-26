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

/* $Id: AliTRDqaBasic.cxx 23387 2008-01-17 17:25:16Z cblume $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliTRDqaBasic.h"
#include "AliTRDqaAT.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TChain.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

//______________________________________________________________________________

AliTRDqaBasic::AliTRDqaBasic() 
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fTrackCuts(0),
    fStatus(0),
    fnTracks(0),
    fPtIn(0),
    fPtOut(0),
    fPtVtx(0),
    fPtVtxSec(0),
    fPtPt(0)
{
  //
  // default constructor
  //
 
}
//______________________________________________________________________________

AliTRDqaBasic:: AliTRDqaBasic(const AliTRDqaBasic & /*trd*/)
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fTrackCuts(0),
    fStatus(0),
    fnTracks(0),
    fPtIn(0),
    fPtOut(0),
    fPtVtx(0),
    fPtVtxSec(0),
    fPtPt(0)
{
  //
  // Copy constructor
  //

}

//______________________________________________________________________________
AliTRDqaBasic::AliTRDqaBasic(const char *name) 
  : AliAnalysisTask(name,""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fTrackCuts(0),
    fStatus(0),
    fnTracks(0),
    fPtIn(0),
    fPtOut(0),
    fPtVtx(0),
    fPtVtxSec(0),
    fPtPt(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
void AliTRDqaBasic::ConnectInputData(const Option_t *)
{
  // Initialisation of branch container and histograms 

  //AliInfo(Form("*** Initialization of %s", GetName())) ; 

  fChain = (TChain*)GetInputData(0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
}

//________________________________________________________________________
void AliTRDqaBasic::CreateOutputObjects()
{
  // build histograms
 
  fStatus   = new TH1D("status", ";status bit", 32, -0.5, 31.5);
  fPtIn     = new TH1D("ptIn", ";p_{T}, GeV/c", 100, 0, 15);
  fPtOut    = new TH1D("ptOut", ";p_{T}, GeV/c", 100, 0, 15);
  fPtVtx    = new TH1D("ptVtx", ";p_{T}, GeV/c", 100, 0, 15);
  fPtVtxSec = new TH1D("ptVtxSec", ";p_{T}, GeV/c", 100, 0, 15);
  fnTracks  = new TH1D("nTracks", "", 200, -0.5, 199.5);
  
  fPtPt     = new TH2D("ptpt", ";p_{T} from inner plane, GeV/c;p_{T} at vertex, GeV/c", 
		       100, 0, 10, 100, 0, 10);

  Int_t c = 0;
  fOutputContainer = new TObjArray(50);

  fOutputContainer->AddAt(fStatus, c++);
  fOutputContainer->AddAt(fPtIn, c++);
  fOutputContainer->AddAt(fPtOut, c++);
  fOutputContainer->AddAt(fPtVtx, c++);
  fOutputContainer->AddAt(fPtVtxSec, c++);
  fOutputContainer->AddAt(fnTracks, c++);
  fOutputContainer->AddAt(fPtPt, c++);  

  printf("n hist = %d\n", c);
  
  // initialize cuts
  fTrackCuts =  new AliESDtrackCuts("AliESDtrackCuts");
  fTrackCuts->DefineHistograms(0);

  // default cuts for ITS+TPC
  Double_t cov1 = 2;
  Double_t cov2 = 2;
  Double_t cov3 = 0.5;
  Double_t cov4 = 0.5;
  Double_t cov5 = 2;
  Double_t nSigma = 3;

  TString tag("Global tracking");

  fTrackCuts->SetMaxCovDiagonalElements(cov1, cov2, cov3, cov4, cov5);

  fTrackCuts->SetMaxNsigmaToVertex(nSigma);
  fTrackCuts->SetRequireSigmaToVertex(kTRUE);

  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);

  fTrackCuts->SetMinNClustersTPC(50);
  fTrackCuts->SetMaxChi2PerClusterTPC(3.5);
}
//______________________________________________________________________________
void AliTRDqaBasic::Exec(Option_t *) 
{
  // Process one event
  Long64_t entry = fChain->GetReadEntry() ;
  if (!(entry%100)) Info("Exec", "Entry = %lld", entry);

  // Processing of one event 
   
  if (!fESD) {
    //AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  Int_t nTracks = fESD->GetNumberOfTracks();
  fnTracks->Fill(nTracks); 

  // track loop
  for(Int_t i=0; i<nTracks; i++) {
    
    //
    // track selection 
    //
    // param in and Out
    // TRDrefit and TRDPid bit
    //
 
    AliESDtrack *track = fESD->GetTrack(i);
    const AliExternalTrackParam *paramOut = track->GetOuterParam();
    const AliExternalTrackParam *paramIn = track->GetInnerParam();

    // long track ..
    if (!paramIn) continue;
    if (!paramOut) continue;
    
    UInt_t status = track->GetStatus();
    if (!(status & AliESDtrack::kTPCrefit)) continue;

    //if (!fTrackCuts->AcceptTrack(track)) continue;
    
    AliTRDqaAT::FillStatus(fStatus, status);

    fPtIn->Fill(paramIn->Pt());
    fPtOut->Fill(paramOut->Pt());
    fPtVtx->Fill(track->Pt());
    if (track->GetConstrainedParam())
      fPtVtxSec->Fill(track->GetConstrainedParam()->Pt());
    
    fPtPt->Fill(paramIn->Pt(), track->Pt());
  }

  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliTRDqaBasic::Terminate(Option_t *)
{
  // save histograms
  fOutputContainer = (TObjArray*)GetOutputData(0);
  
  TFile *file = new TFile("outBasic.root", "RECREATE");
  fOutputContainer->Write();
 
  file->Flush();
  file->Close();
  delete file;

  //for(Int_t i=0; i<fOutputContainer->GetEntries(); i++) {
  //  TObject *obj = fOu
  // }
}

//______________________________________________________________________________
