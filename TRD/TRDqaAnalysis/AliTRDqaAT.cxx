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

/* $Id: AliTRDqaAT.cxx  $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
// In this class provides a commonly used tools as static functions.
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliTRDqaAT.h"

#include "TMath.h"
#include "TH1D.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"

//______________________________________________________________________________

AliTRDqaAT::AliTRDqaAT() {
  //
  // Dummy contructor
  //

}

//______________________________________________________________________________
Int_t AliTRDqaAT::GetSector(const Double_t alpha) 
{
  // Gets the sector number 

  Double_t size = TMath::DegToRad() * 20.;
  Int_t sector = (Int_t)((alpha + TMath::Pi())/size);
  return sector;
}

//______________________________________________________________________________

Int_t AliTRDqaAT::GetStack(const AliExternalTrackParam *paramOut) 
{
  //
  // calculates the stack the track is in
  //
  
  const Double_t l = -0.9;
  const Double_t w = (2*l)/5;

  Double_t tan = paramOut->GetZ() / paramOut->GetX();
  Double_t pos = (tan - l) / w;
  return (Int_t) pos;
}

//______________________________________________________________________________

void AliTRDqaAT::BuildRatio(TH1D *ratio, TH1D *histN, TH1D*histD) {
  //
  // Calculate the ratio of two histograms 
  // error are calculated assuming the histos have the same counts
  //

  // calclate

  Int_t nbins = histN->GetXaxis()->GetNbins();
  for(Int_t i=1; i<nbins+2; i++) {
    
    Double_t valueN = histN->GetBinContent(i);
    Double_t valueD = histD->GetBinContent(i);
    
    if (valueD < 1) {
      ratio->SetBinContent(i, 0);
      ratio->SetBinError(i, 0);
      continue;
    }

    Double_t eps = (valueN < valueD-valueN)? valueN : valueD-valueN;
    
    ratio->SetBinContent(i, valueN/valueD);
    ratio->SetBinError(i, TMath::Sqrt(eps)/valueD);
  }

  // style
  ratio->SetMinimum(-0.1);
  ratio->SetMaximum(1.1);
  ratio->SetMarkerStyle(20);
}
//__________________________________________________________________________

void AliTRDqaAT::FillStatus(TH1D *fStatusHist, UInt_t status) {

  UInt_t u = 1;
  //UInt_t status = track->GetStatus();
  for(Int_t bit=0; bit<32; bit++) 
    if (u<<bit & status) fStatusHist->Fill(bit);
}

//__________________________________________________________________________

void AliTRDqaAT::PrintPID(const AliESDtrack *track) {

  Int_t id = AliPID::kElectron;
  Double_t pidESD[5], pidITS[5], pidTPC[5], pidTRD[5];
  
  track->GetESDpid(pidESD);
  track->GetITSpid(pidITS);
  track->GetTPCpid(pidTPC);
  track->GetTRDpid(pidTRD);
  //track->GetTOFpid(pidTOF);
  
  Double_t comb = pidITS[id] * pidTPC[id] * pidTRD[id];

  Double_t normTot = 0;
  Double_t norm[5] = {0,0,0,0,0};
  for(Int_t i=0; i<5; i++) {
    norm[i] = pidITS[i] * pidTPC[i] * pidTRD[i];
    normTot += norm[i];
  }
  
  Double_t trdTot = 0;
  for(Int_t i=0; i<AliPID::kSPECIES; i++) trdTot += pidTRD[i];
    

  //  comb /= norm;
  
  printf("%.3f | %.3f %.3f | %.3f %.3f %.3f | ",
	 pidESD[id], comb, normTot, pidITS[id], pidTPC[id], pidTRD[id]);

  for(Int_t i=0; i<5; i++) printf("%.2f:%.2f:%.2f | ", pidITS[i], pidTPC[i], pidTRD[i]);
  printf("| %.3f |\n", trdTot);
}

//__________________________________________________________________________
