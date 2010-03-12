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


#include "AliOmegaDalitz.h"
#include <TMath.h>
#include <AliLog.h>
#include <TH1.h>
#include <TPDGCode.h>
#include <TRandom.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

ClassImp(AliOmegaDalitz)


//-----------------------------------------------------------------------------
//
// Generate lepton-pair mass distributions for Dalitz decays according
// to the Kroll-Wada parametrization: N. Kroll, W. Wada: Phys. Rev 98(1955)1355
//
// For the electromagnetic form factor the parameterization from
// Lepton-G is used: L.G. Landsberg et al.: Phys. Rep. 128(1985)301
//
//-----------------------------------------------------------------------------
//
// Adaption for AliRoot
//
// R. Averbeck
// A. Morsch
//
AliOmegaDalitz::AliOmegaDalitz():
	AliDecayer(),
	fLPMass(0)
{
    // Constructor
}

void AliOmegaDalitz::Init()
{
  // Initialisation
  Int_t    ibin, nbins = 1000;
  Double_t min, max, binwidth;
  Double_t pmass, lmass, omass;
  Double_t epsilon, delta, mLL, q, kwHelp, krollWada, formFactor, weight;

  // Get the particle masses
  // electron
  lmass = (TDatabasePDG::Instance()->GetParticle(kElectron))->Mass();
  // omega
  pmass = (TDatabasePDG::Instance()->GetParticle(223))      ->Mass();
  // pi0
  omass = (TDatabasePDG::Instance()->GetParticle(kPi0))     ->Mass();

  min = 2.0 * lmass;
  max = pmass - omass;
  binwidth = (max - min) / (Double_t)nbins;
  fLPMass = new TH1F("fLPMass", "Dalitz", nbins, min, max);

  epsilon = (lmass / pmass) * (lmass / pmass);
  delta   = (omass / pmass) * (omass / pmass);

  for ( ibin = 1; ibin <= nbins; ibin++ )
  {
    mLL = min + (Double_t)(ibin - 1) * binwidth + binwidth / 2.0;
    q    = (mLL / pmass) * (mLL / pmass);
    if ( q <= 4.0 * epsilon )
    {
	AliFatal("Error in calculating Dalitz mass histogram binning!");
	fLPMass = 0;
    }	

    kwHelp = (1.0 + q /  (1.0 - delta)) * (1.0 + q / (1.0 - delta))
            - 4.0 * q / ((1.0 - delta) * (1.0 - delta));    
    if ( kwHelp <= 0.0 )
    {
	AliFatal("Error in calculating Dalitz mass histogram binning!");
	fLPMass = 0;
    }	
    krollWada = (2.0 / mLL) * TMath::Exp(1.5 * TMath::Log(kwHelp))
                            * TMath::Sqrt(1.0 - 4.0 * epsilon / q)   
                            * (1.0 + 2.0 * epsilon / q);
    formFactor = 
	(TMath::Power(TMath::Power(0.6519,2),2)) 
	/ (TMath::Power(TMath::Power(0.6519,2)-TMath::Power(mLL, 2), 2) 
	   + TMath::Power(0.04198, 2)*TMath::Power(0.6519, 2));
    weight = krollWada * formFactor;   
    printf(" %5d %13.3f \n", ibin, weight);
    
    fLPMass->AddBinContent(ibin, weight);
  }
}


void AliOmegaDalitz::Decay(Int_t /*idpart*/, TLorentzVector* pparent)
{
//-----------------------------------------------------------------------------
//
//  Generate Omega Dalitz decay
//
//-----------------------------------------------------------------------------

  Double_t pmass, lmass, omass, lpmass;
  Double_t e1, p1, e3, p3;
  Double_t betaSquare, lambda;
  Double_t costheta, sintheta, cosphi, sinphi, phi;
  
  // Get the particle masses
  // electron
  lmass = (TDatabasePDG::Instance()->GetParticle(kElectron))->Mass();
  // omega
  pmass = (TDatabasePDG::Instance()->GetParticle(223))      ->Mass();
  // pi0
  omass = (TDatabasePDG::Instance()->GetParticle(kPi0))     ->Mass();

  // Sample the lepton pair mass from a histogram
  for( ;; ) 
  {
      lpmass = fLPMass->GetRandom();
      if ( pmass - omass > lpmass && lpmass / 2. > lmass ) break;
  }
  
  // lepton pair kinematics in virtual photon rest frame
  e1 = lpmass / 2.;
  p1 = TMath::Sqrt((e1 + lmass) * (e1 - lmass));
  betaSquare = 1.0 - 4.0 * (lmass * lmass) / (lpmass * lpmass);
  lambda      = betaSquare / (2.0 - betaSquare);
  costheta = (2.0 * gRandom->Rndm()) - 1.;
  sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
  phi      = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
  sinphi   = TMath::Sin(phi);
  cosphi   = TMath::Cos(phi);
  // momentum vectors of leptons in virtual photon rest frame
  Double_t pProd1[3] = {p1 * sintheta * cosphi, 
		       p1 * sintheta * sinphi, 
		       p1 * costheta};
  Double_t pProd2[3] = {-1.0 * p1 * sintheta * cosphi, 
		       -1.0 * p1 * sintheta * sinphi, 
		       -1.0 * p1 * costheta};
  
  // pizero kinematics in omega rest frame
  e3       = (pmass * pmass + omass * omass - lpmass * lpmass)/(2. * pmass);
  p3       = TMath::Sqrt((e3 + omass)  * (e3 - omass));
  costheta = (2.0 * gRandom->Rndm()) - 1.;
  sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
  phi      = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
  sinphi   = TMath::Sin(phi);
  cosphi   = TMath::Cos(phi);	    
  // pizero 4-vector in omega rest frame
  fProducts[2].SetPx(p3 * sintheta * cosphi);
  fProducts[2].SetPy(p3 * sintheta * sinphi);
  fProducts[2].SetPz(p3 * costheta);
  fProducts[2].SetE(e3);

  // lepton 4-vectors in properly rotated virtual photon rest frame
  Double_t pRot1[3] = {0.};
  Rot(pProd1, pRot1, costheta, -sintheta, -cosphi, -sinphi);
  Double_t pRot2[3] = {0.};
  Rot(pProd2, pRot2, costheta, -sintheta, -cosphi, -sinphi); 
  fProducts[0].SetPx(pRot1[0]);
  fProducts[0].SetPy(pRot1[1]);
  fProducts[0].SetPz(pRot1[2]);
  fProducts[0].SetE(e1);
  fProducts[1].SetPx(pRot2[0]);
  fProducts[1].SetPy(pRot2[1]);
  fProducts[1].SetPz(pRot2[2]);
  fProducts[1].SetE(e1);

  // boost the dilepton into the omega's rest frame 
  Double_t eLPparent = TMath::Sqrt(p3 * p3 + lpmass * lpmass);
  TVector3 boostPair( -1.0 * fProducts[2].Px() / eLPparent, 
		      -1.0 * fProducts[2].Py() / eLPparent,
		      -1.0 * fProducts[2].Pz() / eLPparent);
  fProducts[0].Boost(boostPair);
  fProducts[1].Boost(boostPair);

  // boost all decay products into the lab frame 
  TVector3 boostLab( pparent->Px() / pparent->E(), 
		     pparent->Py() / pparent->E(),
		     pparent->Pz() / pparent->E());
  fProducts[0].Boost(boostLab);
  fProducts[1].Boost(boostLab);
  fProducts[2].Boost(boostLab);
  
  return;

}

Int_t AliOmegaDalitz::ImportParticles(TClonesArray *particles)
{
    // Import TParticles in array particles
    // e+ e- pi0
    
    TClonesArray &clonesParticles = *particles;

    Int_t pdg   [3] = {kPi0, kElectron, -kElectron};
    Int_t parent[3] = {-1, 0, 0};
    Int_t d1    [3] = {1, -1, -1};
    Int_t d2    [3] = {2, -1, -1};
    for (Int_t i = 2; i > -1; i--) {
	Double_t px = fProducts[i].Px();
	Double_t py = fProducts[i].Py();
	Double_t pz = fProducts[i].Pz();
	Double_t e  = fProducts[i].E();
	//
	new(clonesParticles[2 - i]) TParticle(pdg[i], 1, parent[i], -1, d1[i], d2[i], px, py, pz, e, 0., 0., 0., 0.);
    }
    return (3);
}



void AliOmegaDalitz::
Rot(Double_t pin[3], Double_t pout[3], Double_t costheta, Double_t sintheta,
    Double_t cosphi, Double_t sinphi)
{
  pout[0] = pin[0]*costheta*cosphi-pin[1]*sinphi+pin[2]*sintheta*cosphi;
  pout[1] = pin[0]*costheta*sinphi+pin[1]*cosphi+pin[2]*sintheta*sinphi;
  pout[2] = -1.0  * pin[0] * sintheta + pin[2] * costheta;
  return;
}

