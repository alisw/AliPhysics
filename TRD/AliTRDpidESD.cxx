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

//-----------------------------------------------------------------
//           Implementation of the TRD PID class
// Assigns the electron and pion liklihoods for each ESD track.
// The AliTRDprobdist class is instantiated here.
// The function MakePID(AliESD *event) calculates the probability
// of having dedx and the probability of having timbin at a given 
// momentum (mom) and particle type k (0 for e) and (2 for pi)
// from the precalculated timbin distributions. 
// Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>
//-----------------------------------------------------------------

#include "AliTRDpidESD.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliTRDprobdist.h"

ClassImp(AliTRDpidESD)

//_________________________________________________________________________
AliTRDpidESD::AliTRDpidESD(Double_t *param)
{
  //
  //  The main constructor
  //
  fMIP=param[0];   // MIP signal
  fRes=param[1];   // relative resolution
  fRange=param[2]; // PID "range" (in sigmas)
}

Double_t AliTRDpidESD::Bethe(Double_t bg) 
{
  //
  // Parametrization of the Bethe-Bloch-curve
  // The parametrization is the same as for the TPC and is taken from Lehrhaus.
  //

  // This parameters have been adjusted to averaged values from GEANT
  const Double_t kP1 = 7.17960e-02;
  const Double_t kP2 = 8.54196;
  const Double_t kP3 = 1.38065e-06;
  const Double_t kP4 = 5.30972;
  const Double_t kP5 = 2.83798;

  // This parameters have been adjusted to Xe-data found in:
  // Allison & Cobb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kP1 = 0.76176E-1;
  //const Double_t kP2 = 10.632;
  //const Double_t kP3 = 3.17983E-6;
  //const Double_t kP4 = 1.8631;
  //const Double_t kP5 = 1.9479;

  // Lower cutoff of the Bethe-Bloch-curve to limit step sizes
  const Double_t kBgMin = 0.8;
  const Double_t kBBMax = 6.83298;
  //const Double_t kBgMin = 0.6;
  //const Double_t kBBMax = 17.2809;
  //const Double_t kBgMin = 0.4;
  //const Double_t kBBMax = 82.0;

  if (bg > kBgMin) {
    Double_t yy = bg / TMath::Sqrt(1. + bg*bg);
    Double_t aa = TMath::Power(yy,kP4);
    Double_t bb = TMath::Power((1./bg),kP5);
             bb = TMath::Log(kP3 + bb);
    return ((kP2 - aa - bb)*kP1 / aa);
  }
  else {
    return kBBMax;
  }

}

//_________________________________________________________________________
Int_t AliTRDpidESD::MakePID(AliESD *event)
{
  //
  //  This function calculates the "detector response" PID probabilities 
  //
  // The class AliTRDprobdist contains precalculated prob dis.
  AliTRDprobdist *pd = new AliTRDprobdist();
  pd->SetADCNorm(1.0); // The factor is the ratio of Mean of pi charge dist.
                    // for the New TRD code divided by the Mean of pi charge
                    // dist. given in AliTRDprobdist object

  //  Example to get mean for particle 2 (pi) and momentum number 4 (2 GeV)
  //  printf("%.2f \n", pd->GetMean(2, 4));
  //  Example of use of Copy Constructor 
  //  AliTRDprobdist *pd1 = new AliTRDprobdist(*pd);

  Int_t ntrk=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    if ((t->GetStatus()&AliESDtrack::kTRDin)==0)
      if ((t->GetStatus()&AliESDtrack::kTRDout)==0)
	if ((t->GetStatus()&AliESDtrack::kTRDrefit)==0) continue;
    if(t->GetTRDsignal()==0) continue;
    //    Int_t ns=AliESDtrack::kSPECIES;
    Int_t ns=AliPID::kSPECIES;
    Double_t p[10];
    Double_t mom=t->GetP();
    Double_t probTotal=0.0;
    for (Int_t j=0; j<ns; j++) {
      p[j]=1.;
      for (Int_t ilayer=0; ilayer <6; ilayer++) {
        Double_t dedx=t->GetTRDsignals(ilayer);
        Int_t timbin=t->GetTRDTimBin(ilayer);
	p[j]*= pd->GetProbability(j,mom,dedx);
	p[j]*= pd->GetProbabilityT(j,mom,timbin);
	p[j]*= 100;
      } // loop over layers
      probTotal+=p[j];
    } //loop over particle species
    //  printf(" %f  %d  %f  %f  %f \n", mom, timbin, p[0], p[1], p[2]);
    for (Int_t j=0; j<ns; j++) {
      if(probTotal) p[j]/= probTotal;
      else p[j]=1.0;
      //      p[j]=1.;
    } //loop over particle species
    t->SetTRDpid(p);
  } //loop over tracks
  delete pd;
  return 0;
}
