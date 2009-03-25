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

/* $Id:$ */

/* History of cvs commits:
 * $Log$
 */

//______________________________________________________
// Author : Aleksei Pavlinov; IHEP, Protvino, Russia
// Feb 17, 2009
// Implementation of fit procedure from ALICE-INT-2008-026:
// "Time and amplitude reconstruction from sampling 
//  measurements of the PHOS signal profile"
//  M.Yu.Bogolyubsky and ..
//
//  Fit by function en*x*x*exp(-2.*x): x = (t-t0)/tau.
//  Tme main goal is fast estimation of amplitude.
//

// --- AliRoot header files ---
#include "AliPHOSFastAltroFit.h"

#include <math.h>

ClassImp(AliPHOSFastAltroFit)

//____________________________________________________________________________
  AliPHOSFastAltroFit:: AliPHOSFastAltroFit() 
: TNamed(), 
  fSig(0),fTau(0),fN(0), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0) 
{
}

//____________________________________________________________________________
AliPHOSFastAltroFit::AliPHOSFastAltroFit(const char* name, const char* title, const Double_t tau)
  : TNamed(name, title), 
  fSig(0),fTau(0),fN(2), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0) 
{
  if(strlen(name)==0) SetName("FastAltroFit");
}

//____________________________________________________________________________
AliPHOSFastAltroFit::~AliPHOSFastAltroFit() 
{  
}

void AliPHOSFastAltroFit::FastFit(Int_t* t, Int_t* y, Int_t n, Double_t sig, Double_t tau, Double_t ped)
{
  fSig = sig;
  fTau = tau;

  Double_t* tn = new Double_t[n]; 
  Double_t* yn = new Double_t[n]; 
  Int_t nn=0;

  DeductPedestal(t,y,n,  tau,ped,  tn,yn,nn);
  //  printf(" n %i : nn %i : ped %f \n", n, nn, ped);
  // for(int i=0; i<nn; i++) 
  // printf(" i %i : yn %7.2f : tn %7.2f \n", i, yn[i], tn[i]); 

  FastFit(tn,yn,nn,sig,tau, fAmp,fAmpErr, fT0,fT0Err,fChi2);

  if(fChi2> 0.0) fNDF = nn - 2;
  else           fNDF = 0; 

  delete [] tn;
  delete [] yn;
}

void AliPHOSFastAltroFit::DeductPedestal(Int_t* t, Int_t* y, Int_t n, Double_t tau, Double_t ped, 
  Double_t* tn, Double_t* yn, Int_t &nn)
{
  Int_t ymax=0, nmax=0;
  for(Int_t i=0; i<n; i++){
    if(y[i]>ymax) {
      ymax = y[i];
      nmax = i;
    }
  }
  Int_t i1 = nmax - Int_t(tau);
  Int_t i2 = n; // No rejection of points from right - should correct for EMCAL
  i1 = i1<0?0:i1;

  nn = 0;
  Double_t yd=0.0;
  for(Int_t i=i1; i<i2; i++) {
    yd = Double_t(y[i]) - ped;
    if(yd>1.) {
      yn[nn] = yd;
      tn[nn] = Double_t(t[i]);
      nn++;
    }
  }
}

void AliPHOSFastAltroFit::FastFit(Double_t* t, Double_t* y, Int_t n, Double_t sig, Double_t tau,
Double_t &amp, Double_t &eamp, Double_t &t0, Double_t &et0, Double_t &chi2)
{
  // It is case of n=k=2 : fnn = x*x*exp(2 - 2*x)
  // Input: 
  //     n  - number of points 
  //   t[n] - array of time bins
  //   y[n] - array of amplitudes after pedestal subtractions;
  //   sig   - error of amplitude measurement (one value for all channels)
  //   tau   - filter time response (in timebin units)
  // Output:
  //       amp - amplitude at t0;
  //        t0 - timeof max amplitude; 
  static Double_t xx; // t/tau
  static Double_t a, b, c;
  static Double_t f02, f12, f22;    // functions
  static Double_t f02d, f12d, f22d; // functions derivations
  chi2 = -1;
  if(n<=0) {
    printf("<I> FastFit : n<=0 \n"); 
    return;
  }

  a = b = c = 0.0;
  for(Int_t i=0; i<n; i++){
    xx  = t[i]/tau;
    f02 = exp(-2.*xx);
    f12 = xx*f02;
    f22 = xx*f12;
    // Derivations
    f02d = -2.*f02;
    f12d = f02 - 2.*f12;
    f22d = 2.*(f12 - f22); 
    //
    a += f02d * y[i];
    b -= 2.*f12d * y[i];
    c += f22d * y[i];
  }
  Double_t t0_1=0.0, t0_2=0.0;
  Double_t amp_1=0.0, amp_2=0.0, chi2_1=0.0, chi2_2=0.0;
  if(QuadraticRoots(a,b,c, t0_1,t0_2)) {
    t0_1 *= tau;
    t0_2 *= tau;
    Amplitude(t,y,n, sig, tau, t0_1, amp_1, chi2_1);
    Amplitude(t,y,n, sig, tau, t0_2, amp_2, chi2_2);
    if(0) {
      printf(" t0_1 %f : t0_2 %f \n", t0_1, t0_2);
      printf(" amp_1 %f : amp_2 %f \n", amp_1, amp_2);
      printf(" chi2_1 %f : chi2_2 %f \n", chi2_1, chi2_2);
    }
    // t0 less on one tau with comparing with value from  "canonical equation"
    amp  = amp_1;
    t0   = t0_1;
    chi2 = chi2_1; 
    if(chi2_1 > chi2_2) {
      amp  = amp_2;
      t0   = t0_2; 
      chi2 = chi2_2; 
    }
    CalculateParsErrors(t, y, n, sig, tau, amp, t0, eamp, et0);

    // Fill1();
    
    // DrawFastFunction(amp, t0, fUtils->GetPedestalValue(), "1");
    //    DrawFastFunction(amp_1, t0_1, fUtils->GetPedestalValue(), "1");
    // DrawFastFunction(amp_2, t0_2, fUtils->GetPedestalValue(), "2");
  } else {
    chi2 = -1.; // bad fit - negative chi2
  }
}

Bool_t AliPHOSFastAltroFit::QuadraticRoots(Double_t a, Double_t b, Double_t c, Double_t &x1, Double_t &x2)
{
  // Resolve quadratic equations a*x**2 + b*x + c
  //printf(" a %12.5e b %12.5e c %12.5e \n", a, b, c);
  static Double_t dtmp = 0.0;
  dtmp = b*b - 4.*a*c;

  if(dtmp>=-1.0e-7 && dtmp<0.0) {
    printf("QuadraticRoots : small negative square : dtmp %f ", dtmp);
    dtmp = 0.0;
  }
  if(dtmp>=0.0) {
    dtmp = sqrt(dtmp);
    x1   = (-b + dtmp) / (2.*a);
    x2   = (-b - dtmp) / (2.*a);

    //    printf(" x1 %f : x2 %f \n", x1, x2);
    return kTRUE;
  } else {
    printf("QuadraticRoots : negative square : dtmp %f ", dtmp);
    return kFALSE;
  }
}

void AliPHOSFastAltroFit::Amplitude(Double_t* t,Double_t* y,Int_t n, Double_t sig, Double_t tau, 
Double_t t0, Double_t &amp, Double_t &chi2)
{  
  // Calculate parameters error too - Mar 24,09
  // sig is independent from points
  amp = 0.;
  Double_t x=0.0, f=0.0, den=0.0, f02;
  for(Int_t i=0; i<n; i++){
    x    = (t[i] - t0)/tau;
    f02  = exp(-2.*x);
    f    = x*x*f02;     
    amp += f*y[i];
    den += f*f;
  }
  amp /= den;
  //
  // chi2 calculation
  //
  Double_t dy=0.0;
  chi2=0.;
  for(Int_t i=0; i<n; i++){
    x    = (t[i] - t0)/tau;
    f02  = exp(-2.*x);
    f    = amp*x*x*f02;
    dy   = y[i]-f;
    chi2 += dy*dy;
    //printf(" %i : y %f -> f %f : dy %f \n", i, y[i], f, dy); 
  }
  chi2 /= (sig*sig);
}

void AliPHOSFastAltroFit::CalculateParsErrors(Double_t* t, Double_t* y, Int_t n, Double_t sig, Double_t tau,
Double_t &amp, Double_t &t0, Double_t &eamp, Double_t &et0)
{
  // fmax_nk = (n/k)**n*exp(-n) => n=k=2 => exp(-n) = exp(-2.)
  static Double_t cc = exp(-2.);

  Double_t sumf2=0.0, sumfd2=0.0, x, f02, f12, f22, f22d;

  for(Int_t i=0; i<n; i++){
    x    = (t[i] - t0)/tau;
    f02  = amp*exp(-2.*x);
    f12  = x*f02;
    f22  = x*f12;
    sumf2 += f22 * f22;
    //
    f22d = 2.*(f12 - f22); 
    sumfd2 += f22d * f22d;
  }

  et0  = sig/amp/sqrt(sumfd2);
  amp *= cc;
  eamp = sig*cc/sqrt(sumf2); 
}
