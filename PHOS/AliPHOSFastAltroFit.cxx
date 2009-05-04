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

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <math.h>

ClassImp(AliPHOSFastAltroFit)

//____________________________________________________________________________
  AliPHOSFastAltroFit:: AliPHOSFastAltroFit() 
: TNamed(), 
  fSig(0),fTau(0),fN(0),fPed(0), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0)
,fNfit(0),fTfit(0),fAmpfit(0), fStdFun(0)
{
}

//____________________________________________________________________________
AliPHOSFastAltroFit::AliPHOSFastAltroFit(const char* name, const char* title, const Double_t tau)
  : TNamed(name, title), 
  fSig(0),fTau(tau),fN(2),fPed(0), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0) 
 ,fNfit(0),fTfit(0),fAmpfit(0), fStdFun(0)
{
  if(strlen(name)==0) SetName("FastAltroFit");
}

//____________________________________________________________________________
AliPHOSFastAltroFit::AliPHOSFastAltroFit(const AliPHOSFastAltroFit &obj)
  : TNamed(obj), 
  fSig(0),fTau(0),fN(2),fPed(0), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0) 
 ,fNfit(0),fTfit(0),fAmpfit(0), fStdFun(0)
{
}

//____________________________________________________________________________
AliPHOSFastAltroFit::~AliPHOSFastAltroFit() 
{  
  if(fTfit) delete [] fTfit;
  if(fAmpfit) delete [] fAmpfit;
}

//____________________________________________________________________________
AliPHOSFastAltroFit& AliPHOSFastAltroFit::operator= (const AliPHOSFastAltroFit &/*obj*/)
{
  // Not implemented yet
  return *this;
}

//____________________________________________________________________________
void AliPHOSFastAltroFit::FastFit(TH1F* h, Double_t sig, Double_t tau, Double_t ped)
{
  // Service method for convinience only  
  Reset();

  if(h==0) return;
  Int_t n = h->GetNbinsX();
  if(n<=0) return;

  Int_t* t = new Int_t[n];
  Int_t* y = new Int_t[n];

  for(Int_t i=0; i<n; i++) {
    t[i] = Int_t(h->GetBinCenter(i+1));
    y[i] = Int_t(h->GetBinContent(i+1));
  }
  FastFit(t,y,n, sig,tau,ped);

  delete [] t;
  delete [] y;
}

void AliPHOSFastAltroFit::FastFit(Int_t* t, Int_t* y, Int_t n, Double_t sig, Double_t tau, Double_t ped)
{
  Reset();

  fSig = sig;
  fTau = tau;
  fPed = ped;

  if(fTfit) delete [] fTfit;
  if(fAmpfit) delete [] fAmpfit;

  fNfit   = 0;
  fTfit   = new Double_t[n]; 
  fAmpfit = new Double_t[n]; 

  DeductPedestal(t,y,n,  tau,ped,  fTfit,fAmpfit,fNfit);
  //  printf(" n %i : fNfit %i : ped %f \n", n, fNfit, ped);
  // for(int i=0; i<fNfit; i++) 
  // printf(" i %i : fAmpfit %7.2f : fTfit %7.2f \n", i, fAmpfit[i], fTfit[i]); 

  if(fNfit>=2) {
    FastFit(fTfit,fAmpfit,fNfit,sig,tau, fAmp,fAmpErr, fT0,fT0Err,fChi2);

    if(fChi2> 0.0) fNDF = fNfit - 2;
    else           fNDF = 0; 
  } else if(fNfit==1){
    Reset(); // What to do here => fT0 = fTfit[0]; fAmp = fAmpFit[0] ??
  } else {
    Reset();
  }
}

void AliPHOSFastAltroFit::Reset()
{
  fSig  = fTau = 0.0;
  fAmp  = fAmpErr = fT0 = fT0Err = 0.0;
  fChi2 = -1.;
  fNDF  = fNfit = 0;
  fTfit = fAmpfit = 0;
}


void AliPHOSFastAltroFit::GetFitResult(Double_t &amp,Double_t &eamp,Double_t &t0,Double_t &et0, 
Double_t &chi2, Int_t &ndf) 
{
  amp  = fAmp;
  eamp = fAmpErr;
  t0   = fT0;
  et0  = fT0Err;
  chi2 = fChi2;
  ndf  = fNDF;
}

void AliPHOSFastAltroFit::GetFittedPoints(Int_t &nfit, Double_t* ar[2])
{
  nfit  = fNfit;
  ar[0] = fTfit;
  ar[1] = fAmpfit;
}

void AliPHOSFastAltroFit::DeductPedestal(Int_t* t, Int_t* y, Int_t n, Double_t tau, Double_t ped, 
  Double_t* tn, Double_t* yn, Int_t &nn)
{
  static Double_t yMinUnderPed=2.; // should be tune
  Int_t ymax=0, nmax=0;
  for(Int_t i=0; i<n; i++){
    if(y[i]>ymax) {
      ymax = y[i];
      nmax = i;
    }
  }
  Int_t i1 = nmax - Int_t(tau);
  //i1 = 0;
  i1 = i1<0?0:i1;
  Int_t i2 = n;

  nn = 0;
  Double_t yd=0.0, tdiff=0.0;;
  for(Int_t i=i1; i<i2; i++) {
    if(ped>0.0) {
      yd = Double_t(y[i]) - ped;
    } else {
      yd = Double_t(y[i]);
    }
    if(yd < yMinUnderPed) continue;

    if(i>i1 && nn>0){
      tdiff = t[i] - tn[nn-1];
      //      printf(" i %i : nn %i : tdiff %6.2f : tn[nn] %6.2f \n", i,nn, tdiff, tn[nn-1]);
      if(tdiff>1.) {
     // discard previous points if its are before maximum point and with gap>1
        if(i<nmax ) {
          nn = 0; // nn--;
     // if point with gap after maximum - finish selection
        } else if(i>=nmax ) {
          break;
        }
      }
     // Far away from maximum
     //if(i-nmax > Int_t(5*tau))              break;
    }
    tn[nn] = Double_t(t[i]);
    yn[nn] = yd;
    //printf("i %i : nn %i : tn %6.2f : yn %6.2f \n", i, nn, tn[nn], yn[nn]); 
    nn++;
  }
  //printf(" nmax %i : n %i : nn %i i1 %i \n", nmax, n, nn, i1);
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
  //        t0 - time of max amplitude; 
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
    if(tau<3.) { // EMCAL case : small tau 
      t0 += -0.03; 
      Amplitude(t,y,n, sig, tau, t0, amp, chi2);
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
    printf("QuadraticRoots : small negative square : dtmp %f \n", dtmp);
    dtmp = 0.0;
  }
  if(dtmp>=0.0) {
    dtmp = sqrt(dtmp);
    x1   = (-b + dtmp) / (2.*a);
    x2   = (-b - dtmp) / (2.*a);

    //    printf(" x1 %f : x2 %f \n", x1, x2);
    return kTRUE;
  } else {
    printf("QuadraticRoots : negative square : dtmp %f \n", dtmp);
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
  if(den>0.0) amp /= den;
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
    //    printf(" %i : y %f -> f %f : dy %f \n", i, y[i], f, dy); 
  }
  chi2 /= (sig*sig);
}

void AliPHOSFastAltroFit::CalculateParsErrors(Double_t* t, Double_t* /*y*/, Int_t n, Double_t sig, 
					      Double_t tau, Double_t &amp, Double_t &t0, 
					      Double_t &eamp, Double_t &et0)
{
  // fmax_nk = (n/k)**n*exp(-n) => n=k=2 => exp(-n) = exp(-2.)
  static Double_t cc = exp(-2.);
  //   static Double_t cc = exp(-fN); // mean(N)~1.5 ??

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

//
// Drawing
//
TCanvas* AliPHOSFastAltroFit::DrawFastFunction()
{
  // QA of fitting
  if(fNfit<=0) return 0; // no points

  static TCanvas *c = 0;
  if(c==0) {
    c =  new TCanvas("fastFun","fastFun",800,600);
  }

  c->cd();
  
  Double_t* eamp = new Double_t[fNfit];
  Double_t* et   = new Double_t[fNfit];

  for(Int_t i=0; i<fNfit; i++) {
    eamp[i] = fSig;
    et[i]   = 0.0;
  }

  TGraphErrors *gr = new TGraphErrors(fNfit, fTfit,fAmpfit, et,eamp);
  gr->Draw("Ap");
  gr->SetTitle(Form("Fast Fit : #chi^{2}/ndf = %8.2f / %i", GetChi2(), GetNDF()));
  gr->GetHistogram()->SetXTitle(" time bin ");
  gr->GetHistogram()->SetYTitle(" amplitude ");

  if(fStdFun==0) {
     fStdFun = new TF1("stdFun", StdResponseFunction, 0., fTfit[fNfit-1]+2., 5);    
     fStdFun->SetParNames("amp","t0","tau","N","ped");
  }
  fStdFun->SetParameter(0, GetEnergy());
  fStdFun->SetParameter(1, GetTime() + GetTau());
  fStdFun->SetParameter(2, GetTau());  // 
  fStdFun->SetParameter(3, GetN());    // 2
  fStdFun->SetParameter(4, 0.);  // 
  
  fStdFun->SetLineColor(kBlue);
  fStdFun->SetLineWidth(1);

  fStdFun->Draw("same");

  delete [] eamp;
  delete [] et;

  c->Update();

  return c;
}

Double_t AliPHOSFastAltroFit::StdResponseFunction(Double_t *x, Double_t *par)
{
  // Standard Response Function : 
  // look to Double_t AliEMCALRawUtils::RawResponseFunction(Double_t *x, Double_t *par)
  // Using for drawing only.
  // 
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function (n=2) of the general form
  //
  // t' = (t - t0 + tau) / tau
  // F = A * t**N * exp( N * ( 1 - t) )   for t >= 0
  // F = 0                                for t < 0 
  //
  // parameters:
  // A:   par[0]   // Amplitude = peak value
  // t0:  par[1]
  // tau: par[2]
  // N:   par[3]
  // ped: par[4]
  //
  static Double_t signal , tau, n, ped, xx;

  tau = par[2];
  n   = par[3];
  ped = par[4];
  xx = ( x[0] - par[1] + tau ) / tau ;

  if (xx <= 0) 
    signal = ped ;  
  else {  
    signal = ped + par[0] * TMath::Power(xx , n) * TMath::Exp(n * (1 - xx )) ; 
  }
  return signal ;  
}
