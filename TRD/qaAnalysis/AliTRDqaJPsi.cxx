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
 * provided "as iLs" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
// This class is a part of a package of high level QA monitoring for TRD.
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliTRDqaJPsi.h"
#include "AliTRDqaAT.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"

#include "AliESDEvent.h" 
#include "AliESDtrack.h"
#include "AliKFParticle.h"

#include "TLorentzVector.h"

//______________________________________________________________________________

AliTRDqaJPsi::AliTRDqaJPsi() 
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fnKFtracks(0)
{
  //
  // default dummy constructor
  //

  for (Int_t k = 0; k < fgknSteps; k++) {
    fStatus[k]      = 0x0;
    fAngleSM[k]     = 0x0;
    fInvMass[k]     = 0x0;
    fInvMassVec[k]  = 0x0;
    fInvMassDiff[k] = 0x0;
    fPtAngle[k]     = 0x0;
  }

  for (Int_t l = 0; l < 2*fgknSteps; l++) {
    fnTracks[l]     = 0x0;
    fPt[l]          = 0x0;
    fPID[l]         = 0x0;
  }

  for (Int_t i = 0; i < 1000; i++) {
    fSM[i]          = 0;
    fTracks[i]      = 0x0;
    fVec[i]         = 0x0;
    for (Int_t j = 0; j < fgknSteps; j++) {
      fInSample[i][j] = 0;
    }
  }

}

//______________________________________________________________________________
AliTRDqaJPsi:: AliTRDqaJPsi(const AliTRDqaJPsi & /*trd*/)
  : AliAnalysisTask("",""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fnKFtracks(0)
{
  //
  // Copy constructor
  //

  for (Int_t k = 0; k < fgknSteps; k++) {
    fStatus[k]      = 0x0;
    fAngleSM[k]     = 0x0;
    fInvMass[k]     = 0x0;
    fInvMassVec[k]  = 0x0;
    fInvMassDiff[k] = 0x0;
    fPtAngle[k]     = 0x0;
  }

  for (Int_t l = 0; l < 2*fgknSteps; l++) {
    fnTracks[l]     = 0x0;
    fPt[l]          = 0x0;
    fPID[l]         = 0x0;
  }

  for (Int_t i = 0; i < 1000; i++) {
    fSM[i]          = 0;
    fTracks[i]      = 0x0;
    fVec[i]         = 0x0;
    for (Int_t j = 0; j < fgknSteps; j++) {
      fInSample[i][j] = 0;
    }
  }

}

//______________________________________________________________________________
AliTRDqaJPsi::AliTRDqaJPsi(const char *name) 
  : AliAnalysisTask(name,""),  
    fChain(0),
    fESD(0),
    fOutputContainer(0),
    fnKFtracks(0)
    
{
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0,  TObjArray::Class()) ; 

  for (Int_t k = 0; k < fgknSteps; k++) {
    fStatus[k]      = 0x0;
    fAngleSM[k]     = 0x0;
    fInvMass[k]     = 0x0;
    fInvMassVec[k]  = 0x0;
    fInvMassDiff[k] = 0x0;
    fPtAngle[k]     = 0x0;
  }

  for (Int_t l = 0; l < 2*fgknSteps; l++) {
    fnTracks[l]     = 0x0;
    fPt[l]          = 0x0;
    fPID[l]         = 0x0;
  }

  for (Int_t i = 0; i < 1000; i++) {
    fSM[i]          = 0;
    fTracks[i]      = 0x0;
    fVec[i]         = 0x0;
    for (Int_t j = 0; j < fgknSteps; j++) {
      fInSample[i][j] = 0;
    }
  }

}

//______________________________________________________________________________
void AliTRDqaJPsi::ConnectInputData(const Option_t *)
{
  // Initialisation of branch container and histograms 

  //AliInfo(Form("*** Initialization of %s", GetName())) ; 

  fChain = (TChain*)GetInputData(0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
}

//________________________________________________________________________
void AliTRDqaJPsi::CreateOutputObjects()
{
  //
  // Create the objects that contain the analysis output
  //
  
  Int_t c = 0;
  fOutputContainer = new TObjArray(100);
  
  const char *charge[2] = {"Neg", "Pos"};
  
  // build histograms
  
  for(Int_t i=0; i<fgknSteps; i++) {

    fStatus[i] = new TH1D(Form("status_%d", i), "status", 32, -0.5, 31.5);
    fOutputContainer->AddAt(fStatus[i], c++);
    
    fInvMass[i] = new TH1D(Form("mass_%d", i), ";m_{inv} (GeV);", 100, 0, 5);
    fOutputContainer->AddAt(fInvMass[i], c++);
    
    fInvMassVec[i] = new TH1D(Form("massVec_%d", i), ";m_{inv} (GeV);", 100, 0, 5);
    fOutputContainer->AddAt(fInvMassVec[i], c++);
    
    fInvMassDiff[i] = new TH1D(Form("massDiff_%d", i), ";m_{inv} (GeV);", 100, -1, 1);
    fOutputContainer->AddAt(fInvMassDiff[i], c++);
    
    fAngleSM[i] = new TH1D(Form("angleSM_%d", i), ";#delta SM", 19, -0.5, 18.5); 	
    fOutputContainer->AddAt(fAngleSM[i], c++);
    
    fPtAngle[i] = new TH2D(Form("ptAngle_%d", i), ";p_{T} (GeV/c);#delta SM",
			   20, 0, 5, 10, -0.5, 9.5);
    fOutputContainer->AddAt(fPtAngle[i], c++);
    

    for(Int_t j=0; j<2; j++) {
      fnTracks[j*fgknSteps+i] = 
	new TH1D(Form("nTracks%s_%d", charge[j],i), Form("%s;number of tracks",charge[j]), 100, -0.5, 99.5);
      fPt[j*fgknSteps+i] = new TH1D(Form("pt%s_%d", charge[j], i), Form("%s;p_{T} (GeV/c)", charge[j]), 100, 0, 5);
      fPID[j*fgknSteps+i] = new TH1D(Form("pid%s_%d", charge[j], i), ";electron LQ", 100, 0, 1);

      fOutputContainer->AddAt(fnTracks[j*fgknSteps+i], c++);
      fOutputContainer->AddAt(fPt[j*fgknSteps+i], c++); 
      fOutputContainer->AddAt(fPID[j*fgknSteps+i], c++); 
    }
  }

  //TH2D *fnGoodTracks;  

  printf("n hist = %d\n", c);
}
//______________________________________________________________________________
void AliTRDqaJPsi::Exec(Option_t *) 
{
  /*
    Selection steps:
    - Parameters In and Out
    - TRDrefit bit
    - TRDpid and quality
    - pt > 0.8 GeV/c
    - PID > 0.9
   */


  // Process one event
  Long64_t entry = fChain->GetReadEntry() ;
  if (!(entry%100)) Info("Exec", "Entry = %lld", entry);

  // Processing of one event 
   
  if (!fESD) {
    //AliError("fESD is not connected to the input!") ; 
    return ; 
  }
  
  
  Int_t nTracks = fESD->GetNumberOfTracks();
  Int_t cTracks[2*fgknSteps] = {0,0,0,0,0,0,0,0,0,0};
  fnKFtracks = 0;

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

    Int_t step    = 0;
    Int_t charge  = (track->Charge() > 0) ? 1 : 0;
    UInt_t status = track->GetStatus();
    Double_t pt   = track->Pt();
    Double_t pid  = track->GetTRDpid(AliPID::kElectron);
    
    Double_t esdPid[5];
    track->GetESDpid(esdPid); 

    // create a kalman particle
    Int_t pdg = (charge == 0)? -11 : 11;
    for(Int_t k=0; k<fgknSteps; k++) fInSample[fnKFtracks][k] = 0;  
 
    fVec[fnKFtracks] = CreateVector(track);
    fTracks[fnKFtracks] = new AliKFParticle(*track, pdg);
    fSM[fnKFtracks] = AliTRDqaAT::GetSector(paramOut->GetAlpha());
    fnKFtracks++;

    //AliTRDqaAT::PrintPID(track);

    // apply the cuts

    cTracks[fgknSteps *charge + step]++;
    FillHist(track, step++);

    if (!(status & AliESDtrack::kTRDrefit)) continue;

    cTracks[fgknSteps *charge + step]++;    
    FillHist(track, step++);    

    if (!(status & AliESDtrack::kTRDpid)) continue;
    if (track->GetTRDntracklets() < 6) continue;

    cTracks[fgknSteps *charge + step]++;
    FillHist(track, step++);  

    if (pt < 0.8) continue;

    cTracks[fgknSteps *charge + step]++;
    FillHist(track, step++);      

    if (pid < 0.3) continue; // 
    //if (esdPid[AliPID::kElectron] < 0.5) continue;
    
    cTracks[fgknSteps *charge + step]++;
    FillHist(track, step);  

    for(Int_t k=0; k<2*fgknSteps; k++) fnTracks[k]->Fill(cTracks[k]);
  }

  // calculate invariant mass

  for(Int_t k=0; k<fgknSteps; k++) {
    for(Int_t i=0; i<fnKFtracks; i++) {
      if (!fInSample[i][k]) continue;
      for(Int_t j=i+1; j<fnKFtracks; j++) {
	if (!fInSample[j][k]) continue;
	AliKFParticle jpsi(*(fTracks[i]), *(fTracks[j]));
	TLorentzVector jpsiVec = (*(fVec[i])) + (*fVec[j]);
	fInvMass[k]->Fill(jpsi.GetMass());
	fInvMassVec[k]->Fill(jpsiVec.M());
	fInvMassDiff[k]->Fill(jpsiVec.M() - jpsi.GetMass());
	
	if (jpsi.GetMass() > 2.5 && jpsi.GetMass() < 3.5) {
	  Int_t dSM = TMath::Abs(fSM[i] - fSM[j]);
	  if (dSM > 9) dSM = (18-dSM);
	  fAngleSM[k]->Fill(dSM);
	  fPtAngle[k]->Fill(TMath::Hypot(jpsi.GetPx(), jpsi.GetPy()), dSM);
	}
	
      }
    }
  }

  PostData(0, fOutputContainer);
}

//______________________________________________________________________________
void AliTRDqaJPsi::Terminate(Option_t *)
{
  // save histograms
  fOutputContainer = (TObjArray*)GetOutputData(0);
  
  TFile *file = new TFile("outJPsi.root", "RECREATE");
  fOutputContainer->Write();
  
  file->Flush();
  file->Close();
  delete file;  

  //for(Int_t i=0; i<fOutputContainer->GetEntries(); i++) {
  //  TObject *obj = fOu
  // }
}

//______________________________________________________________________________

void AliTRDqaJPsi::FillHist(AliESDtrack *track, Int_t step) {
  //
  // Fill the histograms
  //

  Int_t charge  = (track->Charge() > 0) ? 1 : 0;
  UInt_t status = track->GetStatus();
  Double_t pt   = track->Pt();
  Double_t pid  = track->GetTRDpid(AliPID::kElectron);
  
  Double_t esdPid[5];
  track->GetESDpid(esdPid);

  Int_t id = charge * fgknSteps + step;
  AliTRDqaAT::FillStatus(fStatus[step], status);
  fPt[id]->Fill(pt);
  fPID[id]->Fill(pid);
  //fPID[id]->Fill(esdPid[AliPID::kElectron]);

  fInSample[fnKFtracks-1][step] = 1;
}

//______________________________________________________________________________

TLorentzVector *AliTRDqaJPsi::CreateVector(AliESDtrack *track) {

  TLorentzVector *vec = new TLorentzVector();
  vec->SetXYZM(track->Px(), track->Py(), track->Pz(), 0);
  return vec;
}

//______________________________________________________________________________
