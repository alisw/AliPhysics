
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

/* $Id: AliTRDqaEnergyDeposit.cxx $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
// In this class the dEdX is analyzed as a function of the particle type,
// momentum and chamber.
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliTRDqaEnergyDeposit.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliKalmanTrack.h"

//______________________________________________________________________________

AliTRDqaEnergyDeposit::AliTRDqaEnergyDeposit() 
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0)
{
  // Dummy default constructor
}
//______________________________________________________________________________

AliTRDqaEnergyDeposit:: AliTRDqaEnergyDeposit(AliTRDqaEnergyDeposit& /*trd*/)
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0)
{
  // Dummy copy constructor
  
  //return *this;
}

//______________________________________________________________________________
AliTRDqaEnergyDeposit::AliTRDqaEnergyDeposit(const char *name) 
  : AliAnalysisTask(name,""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0)
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 
}

//______________________________________________________________________________
void AliTRDqaEnergyDeposit::ConnectInputData(const Option_t *)
{
  // Initialisation of branch container and histograms 

  //AliInfo(Form("*** Initialization of %s", GetName())) ; 

  fChain = (TChain*)GetInputData(0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
}

//________________________________________________________________________
void AliTRDqaEnergyDeposit::CreateOutputObjects()
{
  // build histograms
  
  // prepare the scale from 0.5 to 10 GeV
  const Int_t knbinsx = 50;
  Double_t scalex[knbinsx+1];
  Double_t dd = (TMath::Log(10) - TMath::Log(0.5)) / knbinsx;
  for(Int_t ix=0; ix<knbinsx+1; ix++) {
    scalex[ix] = 0.5 * TMath::Exp(dd * ix);
  }

  const Int_t knbinsy = 50;
  Double_t scaley[knbinsy+1];
  for(Int_t iy=0; iy<knbinsy+1; iy++) {
    scaley[iy] = iy * (3e3/100.);
  }
  
  const char *title = ";p_{T};dEdX (a. u.)";
  const char *charge[2] = {"Pos", "Neg"};

  // build histograms
  fOutputContainer = new TObjArray(50);
  Int_t c=0;

  for(Int_t i=0; i<2; i++) {

    fSignalPtSum[i] = new TH2D(Form("ptSig%s", charge[i]), title, knbinsx, scalex, knbinsy, scaley);
    fOutputContainer->AddAt(fSignalPtSum[i], c++);

    for(Int_t j=0; j<AliPID::kSPECIES; j++) {
      Int_t idx = AliPID::kSPECIES*i+j;

      //
      fSignalPtType[idx] = 
	new TH2D(Form("ptSig%s%d", charge[i], j), title, knbinsx, scalex, knbinsy, scaley);
      fOutputContainer->AddAt(fSignalPtType[idx], c++);

      //
      fSignalPtPure[idx] = 
	new TH2D(Form("ptSigPure%s%d", charge[i], j), title, knbinsx, scalex, knbinsy, scaley);
      fOutputContainer->AddAt(fSignalPtPure[idx], c++);
      
      fProb[idx] = new TH1D(Form("prob%s%d", charge[i], j), ";LQ", 100, 0, 1);
      fOutputContainer->AddAt(fProb[idx], c++);
    }  
  }

  printf("n hist = %d\n", c);
}
//______________________________________________________________________________
void AliTRDqaEnergyDeposit::Exec(Option_t *) 
{
  // Process one event
  
  //  Long64_t entry = fChain->GetReadEntry() ;
  
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
    
    Int_t idx = (track->GetSign() > 0) ? 0 : 1;
    Double_t pt = paramOut->Pt();

    Double_t signal = 0;
    for(Int_t k=0; k<6; ++k)
      signal += track->GetTRDslice(k, -1);
    signal /= 6;

    fSignalPtSum[idx]->Fill(pt, signal);
    
    for(Int_t k=0; k<AliPID::kSPECIES; ++k) {
      
      Double_t lq = track->GetTRDpid(k);
      fProb[AliPID::kSPECIES*idx+k]->Fill(lq);
      if (lq > 0.8) fSignalPtType[AliPID::kSPECIES*idx+k]->Fill(pt, signal);
    }
  }

  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliTRDqaEnergyDeposit::FillElectrons() {

  

}

//______________________________________________________________________________
void AliTRDqaEnergyDeposit::Terminate(Option_t *)
{
  // retrieve histograms
  fOutputContainer = (TObjArray*)GetOutputData(0);
  
  // post processing


  // save the results
  TFile *file = new TFile("outEnergyDeposit.root", "RECREATE");
  fOutputContainer->Write();
 
  file->Flush();
  file->Close();
  delete file;

  //for(Int_t i=0; i<fOutputContainer->GetEntries(); i++) {
  //  TObject *obj = fOu
  // }
}

//______________________________________________________________________________
