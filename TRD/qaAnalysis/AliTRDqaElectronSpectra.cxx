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

/* $Id: AliTRDqaElectronSpectra.cxx  $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
//
// The transverse momentum spectrum is analyzed stack-by-stack
// for all tracks, and for electron tracks. 
// Tracks have to pass quality cuts. 
// Electrons are waighted with the PID LQ
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliTRDqaElectronSpectra.h"
#include "AliTRDqaAT.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"

//______________________________________________________________________________

AliTRDqaElectronSpectra::AliTRDqaElectronSpectra() 
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fStatus(0),
    fSector(0),
    fTheta(0),  
    fStack(0),  
    fnTracks(0),     
    fnElTracks(0),   
    fTracksRatio(0), 
    fPt(0),         
    fPtElectron(0), 
    fMeanPt(0),         
    fMeanPtElectron(0), 
    fPtStack(0),        
    fPtStackElectron(0),
    fElectronLQ(0)
{
  //
  // default dummy constructor
  //
 
}
//______________________________________________________________________________

AliTRDqaElectronSpectra:: AliTRDqaElectronSpectra(AliTRDqaElectronSpectra& /*trd*/)
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fStatus(0),
    fSector(0), 
    fTheta(0),  
    fStack(0),  
    fnTracks(0),     
    fnElTracks(0),   
    fTracksRatio(0), 
    fPt(0),         
    fPtElectron(0), 
    fMeanPt(0),         
    fMeanPtElectron(0), 
    fPtStack(0),        
    fPtStackElectron(0),
    fElectronLQ(0)
{
  //
  // Dummy copy constructor
  //

  //return *this;
}


//______________________________________________________________________________
AliTRDqaElectronSpectra::AliTRDqaElectronSpectra(const char *name) 
  : AliAnalysisTask(name,""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fStatus(0),
    fSector(0), 
    fTheta(0),  
    fStack(0),  
    fnTracks(0),     
    fnElTracks(0),   
    fTracksRatio(0), 
    fPt(0),         
    fPtElectron(0), 
    fMeanPt(0),         
    fMeanPtElectron(0), 
    fPtStack(0),
    fPtStackElectron(0),
    fElectronLQ(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
void AliTRDqaElectronSpectra::ConnectInputData(const Option_t *)
{
  // Initialisation of branch container and histograms 

  //AliInfo(Form("*** Initialization of %s", GetName())) ; 

  fChain = (TChain*)GetInputData(0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
}

//________________________________________________________________________
void AliTRDqaElectronSpectra::CreateOutputObjects()
{
  // build histograms
 
  fStatus = new TH1D("status", ";status bit", 32, -0.5, 31.5);
  fSector = new TH1D("sector", ";sector", 18, -0.5, 17.5);
  fTheta  = new TH1D("theta", ";theta (rad)", 100, -1, 1);
  fStack  = new TH1D("stack", ";stack", 90, -0.5, 89.5);
  
  fnTracks     = new TH1D("tracks", ";stack;number of tracks", 90, -0.5, 89.5);
  fnElTracks   = new TH1D("elTracks", ";stack;number of electron tracks", 90, -0.5, 89.5); 
  fTracksRatio = new TH1D("fractionElectrons", ";stack;fraction of electron tracks", 90, -0.5, 89.5);
  
  fPt         = new TH1D("pt", "p_{T} (GeV/c)", 50, 0, 10);
  fPtElectron = new TH1D("ptElectron", "p_{T} (GeV/c)", 50, 0, 10);
  
  fMeanPt         = new TH1D("meanPt", ";<p_{T}> (GeV/c)", 100, 0, 5);
  fMeanPtElectron = new TH1D("meanPtElectron", ";<P_{T}> (GeV/c)", 100, 0, 5);
    
  fPtStack         = new TH2D("stackPt", ";stack;p_{T} (GeV/c)", 90, -0.5, 89.5, 50, 0, 10);
  fPtStackElectron = new TH2D("stackPtEl", ";stack;p_{T} (GeV/c)", 90, -0.5, 89.5, 50, 0, 10);
  
  fElectronLQ = new TH1D("elLQ", ";likelyhood", 100, 0, 1);

  Int_t c = 0;
  fOutputContainer = new TObjArray(50);

  fOutputContainer->AddAt(fStatus, c++);
  fOutputContainer->AddAt(fSector, c++);
  fOutputContainer->AddAt(fTheta, c++);
  fOutputContainer->AddAt(fStack, c++);
  
  fOutputContainer->AddAt(fnTracks, c++);
  fOutputContainer->AddAt(fnElTracks, c++);
  fOutputContainer->AddAt(fTracksRatio, c++);
  
  fOutputContainer->AddAt(fPt, c++);
  fOutputContainer->AddAt(fPtElectron, c++);

  fOutputContainer->AddAt(fMeanPt, c++);
  fOutputContainer->AddAt(fMeanPtElectron, c++);

  fOutputContainer->AddAt(fPtStack, c++);
  fOutputContainer->AddAt(fPtStackElectron, c++);
  
  fOutputContainer->AddAt(fElectronLQ, c++);

  printf("n hist = %d\n", c);
}
//______________________________________________________________________________
void AliTRDqaElectronSpectra::Exec(Option_t *) 
{
  // Process one event
  Long64_t entry = fChain->GetReadEntry() ;
  if (!(entry%100)) Info("Exec", "Entry = %ld", entry);

  // Processing of one event 
   
  if (!fESD) {
    //AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  Int_t nTracks = fESD->GetNumberOfTracks();
  //fNTracks->Fill(nTracks); 

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
    if (!(status & AliESDtrack::kTRDrefit)) continue;
    if (!(status & AliESDtrack::kTRDpid)) continue;
    if (track->GetTRDntrackletsPID() < 6) continue;

    Int_t sm = AliTRDqaAT::GetSector(paramOut->GetAlpha());
    Int_t stack = 5*sm + AliTRDqaAT::GetStack(paramOut);
    Double_t lq = track->GetTRDpid(AliPID::kElectron);
    Double_t pt = paramOut->Pt();

    //TH1D *fStatus;  // track status
    fSector->Fill(sm);  
    fStack->Fill(stack); 
    fElectronLQ->Fill(lq);

    fTheta->Fill(paramOut->GetZ() / paramOut->GetX());

    fnTracks->Fill(stack);
    fnElTracks->Fill(stack, lq);

    fPt->Fill(pt);
    fPtElectron->Fill(pt, lq);
        
    fPtStack->Fill(stack, pt);
    fPtStackElectron->Fill(stack, pt, lq);
  }

  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliTRDqaElectronSpectra::Terminate(Option_t *)
{
  // save histograms
  fOutputContainer = (TObjArray*)GetOutputData(0);
  
  // build ratios
  fnTracks     = (TH1D*)fOutputContainer->FindObject("tracks");
  fnElTracks   = (TH1D*)fOutputContainer->FindObject("elTracks");
  fTracksRatio = (TH1D*)fOutputContainer->FindObject("fractionElectrons");

  AliTRDqaAT::BuildRatio(fTracksRatio, fnElTracks, fnTracks);

  // save the results

  TFile *file = new TFile("outElSpectra.root", "RECREATE");
  fOutputContainer->Write();
 
  file->Flush();
  file->Close();
  delete file;

  //for(Int_t i=0; i<fOutputContainer->GetEntries(); i++) {
  //  TObject *obj = fOu
  // }
}

//______________________________________________________________________________
