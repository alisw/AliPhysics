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

/* $Id: AliTRDqaESDFriends.cxx $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
// The residuals of cluster with respect to tracklets are analyzed 
// in this class. This class needs ESDfriends.root
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliTRDqaESDFriends.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliKalmanTrack.h"
#include "AliESDfriend.h"

//______________________________________________________________________________

AliTRDqaESDFriends::AliTRDqaESDFriends() 
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fResiduals(0),
    fResidualsAngle(0) 
{
  // Dummy default constructor
}
//______________________________________________________________________________

AliTRDqaESDFriends:: AliTRDqaESDFriends(AliTRDqaESDFriends& /*trd*/)
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fResiduals(0),
    fResidualsAngle(0) 
{
  // dummy copy constructor

  //return *this;
}


//______________________________________________________________________________
AliTRDqaESDFriends::AliTRDqaESDFriends(const char *name) 
  : AliAnalysisTask(name,""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fResiduals(0),
    fResidualsAngle(0) 
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
void AliTRDqaESDFriends::ConnectInputData(const Option_t *)
{
  // Initialisation of branch container and histograms 

  //AliInfo(Form("*** Initialization of %s", GetName())) ; 

  fChain = (TChain*)GetInputData(0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
}

//________________________________________________________________________
void AliTRDqaESDFriends::CreateOutputObjects()
{
  // build histograms
 
  fResiduals = new TH1D("residuals", ";residuals (cm)", 1000, -1, 1);
  fResidualsAngle = new TH2D("residualsAngle", ";angle (rad);residuals (cm)", 100, -1, 1, 100, -1, 1);

  Int_t c = 0;
  fOutputContainer = new TObjArray(50);

  fOutputContainer->AddAt(fResiduals, c++);
  fOutputContainer->AddAt(fResidualsAngle, c++);

  printf("n hist = %d\n", c);
}
//______________________________________________________________________________
void AliTRDqaESDFriends::Exec(Option_t *) 
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

    // standard selection
    AliESDfriend *fr = (AliESDfriend*)fESD->FindListObject("AliESDfriend");
    if (fr) fESD->SetESDfriend(fr);

    AliESDfriendTrack *f = (AliESDfriendTrack*)track->GetFriendTrack();
    
    if (!f) continue;

    //AliKalmanTrack *trdTrack = 0;
    //if (f) trdTrack = f->GetTRDtrack();
    //if (trdTrack) trdTrack->Print();

    //if (f) f->Dump();
    //if (f) f->Print();
      
    //  fESD->GetList()->Print();
    //AliESDfriend *f = (AliESDfriend*)fESD->FindListObject("ESDfriend");
    //if (f) f->Print();
    // AliKalmanTrack *trdTrack = track->GetTRDtrack();
    //if (trdTrack) trdTrack->Print();
  }

  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliTRDqaESDFriends::Terminate(Option_t *)
{
  // retrieve histograms
  fOutputContainer = (TObjArray*)GetOutputData(0);
  
  // post processing


  // save the results
  TFile *file = new TFile("outESDFriends.root", "RECREATE");
  fOutputContainer->Write();
 
  file->Flush();
  file->Close();
  delete file;

  //for(Int_t i=0; i<fOutputContainer->GetEntries(); i++) {
  //  TObject *obj = fOu
  // }
}

//______________________________________________________________________________
