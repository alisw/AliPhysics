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
//____________________________________________________________________________
//*-- 
//*-- Author: Andreas Morsch (CERN)
//*--
//*--
////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include <TRandom.h>

#include "AliEMCALFast.h"


ClassImp(AliEMCALFast)

//____________________________________________________________________________


Float_t AliEMCALFast::SmearMomentum(Int_t ind, Float_t p)
{
//
//  The relative momentum error, i.e. 
//  (delta p)/p = sqrt (a**2 + (b*p)**2) * 10**-2,
//  where typically a = 0.75 and b = 0.16 - 0.24 depending on multiplicity
//  (the lower value is for dn/d(eta) about 2000, and the higher one for 8000)
//
//  If we include information from TRD b will be by a factor 2/3 smaller.
//
//  ind = 1: high multiplicity
//  ind = 2: low  multiplicity
//  ind = 3: high multiplicity + TRD
//  ind = 4: low  multiplicity + TRD

    Float_t pSmeared;
    Float_t a = 0.75;
    Float_t b = 0.24;

    if (ind == 1) b = 0.24;
    if (ind == 2) b = 0.16;
    if (ind == 3) b = 0.16;    
    if (ind == 4) b = 0.11;    
    
    Float_t sigma = p*TMath::Sqrt(a*a+b*b*p*p)*0.01;
    pSmeared = p + gRandom->Gaus(0., sigma);
    return pSmeared;
}


Float_t AliEMCALFast::Efficiency(Int_t ind, Float_t p)
{
// Tracking efficiency:
// above pt 0.5 GeV practically constant, between 90 and 95 % (again,
// depending on multplicity)
// below 0.5 GeV goes down to about 70% at 0.2 GeV.
// On top of that there is 90% geometrical acceptance for tracking due
// to TPC (dead zones between readout chambers).  
// Tracking 
//      
    Float_t eff = 0.9;
    if (ind == 2) eff  = 0.95;
    if (p < 0.5) eff -= (0.5-p)*0.2/0.3;
// Geometry
    eff *= 0.9;
// Acceptance    
    return eff;
}

Bool_t AliEMCALFast::EmcalAcceptance(Float_t eta, Float_t phi)
{
// 27-oct-05 -> should be change
// EMCAL eta-phi acceptance
    Bool_t acc = kFALSE;
    if (TMath::Abs(eta) < 0.7 &&
	phi > 60.*TMath::Pi()/180. &&
	phi < TMath::Pi())
	acc = kTRUE;
    return acc;
}


Bool_t AliEMCALFast::RandomReject(Float_t eff)
{
//
// Random rejection 
    Bool_t rej = kFALSE;
    Float_t ran = gRandom->Rndm();
    if (ran > eff) rej = kTRUE;
    return rej;
}







