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

#include "AliTRDPartID.h"
#include "AliESDtrack.h"
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <Riostream.h>

ClassImp(AliTRDPartID)


AliTRDPartID::AliTRDPartID()
{
// create a TRD particle identifier

  fBetheBloch = NULL;
  fRes = 0.2;
  fRange = 3.;
}

AliTRDPartID::AliTRDPartID(TF1* betheBloch, Double_t res, Double_t range)
{
// create a TRD particle identifier with custom settings

  fBetheBloch = betheBloch;
  fRes = res;
  fRange = range;
}

AliTRDPartID::~AliTRDPartID()
{
  if (fBetheBloch) delete fBetheBloch;
}


Bool_t AliTRDPartID::MakePID(AliESDtrack* track)
{
// This function calculates the "detector response" PID probabilities 

  if (((track->GetStatus()&AliESDtrack::kTRDin) == 0) &&
      ((track->GetStatus()&AliESDtrack::kTRDout) == 0)) return kFALSE;
  Double_t momentum = track->GetP();
  if (momentum < 0.001) return kFALSE;

  // get the probability densities
  Double_t pSum = 0;
  for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) {
    Double_t expectedSignal = fBetheBloch->Eval(momentum/AliPID::ParticleMass(iSpecies));
    Double_t expectedError = fRes * expectedSignal;
    Double_t measuredSignal = track->GetTRDsignal();
    if (TMath::Abs(measuredSignal - expectedSignal) > fRange * expectedError) {
      track->SetTRDpid(iSpecies, 0.);
    } else {
      Double_t delta = (measuredSignal-expectedSignal) / expectedError;
      const Double_t kInvSqr2Pi = 0.398942280401432703;
      Double_t p = kInvSqr2Pi / expectedError * TMath::Exp(-delta*delta / 2.);
      pSum += p;
      track->SetTRDpid(iSpecies, p);
    }
  }

  // "normalize" the probability densities
  if (pSum <= 0) return kFALSE;
  for (Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++) {
    track->SetTRDpid(iSpecies, track->GetTRDpid(iSpecies) / pSum);
  }

  return kTRUE;
}


void AliTRDPartID::FitBetheBloch(TProfile* dEdxVsBetaGamma)
{
// determine the parameters of the bethe bloch function

  if (fBetheBloch) delete fBetheBloch;
  fBetheBloch = new TF1("fBetheBlochTRD", fcnBetheBloch, 0.001, 100000., 5);
  fBetheBloch->SetParameters(1, 10, 0.00002, 2, 2);
  fBetheBloch->SetParLimits(2, 0., 0.01);
  fBetheBloch->SetParLimits(3, 0., 10.);
  fBetheBloch->SetParLimits(4, 0., 10.);
  fBetheBloch->SetFillStyle(0);
  fBetheBloch->SetLineColor(kRed);
  dEdxVsBetaGamma->Fit(fBetheBloch, "NIR", "goff", 0.6, dEdxVsBetaGamma->GetXaxis()->GetXmax());
}

TF1* AliTRDPartID::CreateBetheBloch(Double_t mass)
{
// create a function for expected dE/dx vs p

  TF1* result = new TF1("betheBlochMass", fcnBetheBlochMass, 
			0.001, 100000., 6);
  result->SetParameter(0, mass);
  if (fBetheBloch) {
    for (Int_t iPar = 0; iPar < 5; iPar++) {
      result->SetParameter(iPar+1, fBetheBloch->GetParameter(iPar));
      result->SetParError(iPar+1, fBetheBloch->GetParError(iPar));
    }
  }
  result->SetFillStyle(0);
  return result;
}

AliTRDPartID* AliTRDPartID::GetFromFile(const char* fileName)
{
// read an AliTRDPartID object from a file

  TFile* pidFile = (TFile*) gROOT->GetListOfFiles()->FindObject(fileName);
  Bool_t fileOpened = kFALSE;
  if (!pidFile) {
    pidFile = TFile::Open(fileName);
    fileOpened = kTRUE;
  }
  if (!pidFile->IsOpen()) {
    cerr << "Can't open " << fileName << " !\n";
    if (fileOpened) delete pidFile;
    return NULL;
  }
  gROOT->cd();

  AliTRDPartID* trdPID = (AliTRDPartID*) pidFile->Get("AliTRDPartID")->Clone();
  if (!trdPID) {
    cerr << "Can't get PID object !\n";
  } else {
    trdPID->GetBetheBloch()->SetFunction(fcnBetheBloch);
  }

  if (fileOpened) {
    pidFile->Close();
    delete pidFile;
  }

  return trdPID;
}

Double_t AliTRDPartID::fcnBetheBloch(Double_t* xx, Double_t* par)
{
// parametrized bethe bloch function (xx[0] = beta*gamma = p/m):
//
//   p0/beta^p3 * [ p1 - log(p2 + 1/(beta*gamma)^p4) - beta^p3 ]
//

  Double_t betaGamma2 = xx[0] * xx[0];
  Double_t beta2 = betaGamma2 / (1. + betaGamma2);
  Double_t betaPar3 = TMath::Power(beta2, par[3]/2.);
  return par[0]/betaPar3 * (par[1] - TMath::Log(TMath::Abs(par[2] + TMath::Power(betaGamma2, -par[4]/2.))) - betaPar3);
}

Double_t AliTRDPartID::fcnBetheBlochMass(Double_t* xx, Double_t* par)
{
// parametrized bethe bloch function:  xx[0] = p,  p0 = mass

  Double_t betaGamma = xx[0] / par[0];
  return fcnBetheBloch(&betaGamma, &par[1]);
}
