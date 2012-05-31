/**************************************************************************
 * Copyright(c) 1998-2010 ALICE Experiment at CERN, All rights reserved. *
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
//  The main goal is fast estimation of amplitude and t0.
//

// --- AliRoot header files ---
#include "AliCaloFastAltroFitv0.h"

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>

#include <math.h>

ClassImp(AliCaloFastAltroFitv0)

//____________________________________________________________________________
  AliCaloFastAltroFitv0::AliCaloFastAltroFitv0() 
: TNamed(), 
  fSig(0),fTau(0),fN(0),fPed(0), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0)
,fNoFit(0),fNfit(0),fTfit(0),fAmpfit(0), fStdFun(0)
{
}

//____________________________________________________________________________
AliCaloFastAltroFitv0::AliCaloFastAltroFitv0(const char* name, const char* title, 
  const Double_t sig, const Double_t tau, const Double_t n)
  : TNamed(name, title), 
  fSig(sig),fTau(tau),fN(n),fPed(0), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0) 
 ,fNoFit(0),fNfit(0),fTfit(0),fAmpfit(0), fStdFun(0)
{
  if(strlen(name)==0) SetName("CaloFastAltroFitv0");
}

//____________________________________________________________________________
AliCaloFastAltroFitv0::AliCaloFastAltroFitv0(const AliCaloFastAltroFitv0 &obj)
  : TNamed(obj), 
  fSig(0),fTau(0),fN(2.),fPed(0), fAmp(0),fAmpErr(0),fT0(0),fT0Err(0),fChi2(0.),fNDF(0) 
 ,fNoFit(0),fNfit(0),fTfit(0),fAmpfit(0), fStdFun(0)
{
}

//____________________________________________________________________________
AliCaloFastAltroFitv0::~AliCaloFastAltroFitv0() 
{  
  if(fTfit) delete [] fTfit;
  if(fAmpfit) delete [] fAmpfit;
}

//____________________________________________________________________________
AliCaloFastAltroFitv0& AliCaloFastAltroFitv0::operator= (const AliCaloFastAltroFitv0 &/*obj*/)
{
  // Not implemented yet
  return (*this);
}

void AliCaloFastAltroFitv0::FastFit(Int_t* t, Int_t* y, Int_t nPoints, Double_t sig, Double_t tau, 
				    Double_t /*n*/, Double_t ped, Double_t tMax)
{
  Reset();

  fSig = sig;
  fTau = tau;
  fPed = ped;

  Int_t ii=0;
  CutRightPart(t,y,nPoints,tMax, ii);
  nPoints = ii;

  fNfit   = 0;
  fTfit   = new Double_t[nPoints]; 
  fAmpfit = new Double_t[nPoints]; 

  
  DeductPedestal(t,y,nPoints,  tau,ped,  fTfit,fAmpfit,fNfit);
  //  printf(" n %i : fNfit %i : ped %f \n", n, fNfit, ped);
  // for(int i=0; i<fNfit; i++) 
  // printf(" i %i : fAmpfit %7.2f : fTfit %7.2f \n", i, fAmpfit[i], fTfit[i]); 

  if(fNfit>=2) {
    FastFit(fTfit,fAmpfit,fNfit,sig,tau, fAmp,fAmpErr, fT0,fT0Err,fChi2);

    if(fChi2> 0.0) {
      fNDF = fNfit - 2;
    } else {
      fNDF = 0;
      fNoFit++;
    }
  } else if(fNfit==1){
    Reset(); // What to do here => fT0 = fTfit[0]; fAmp = fAmpFit[0] ??
  } else {
    Reset();
  }
}

//____________________________________________________________________________
void AliCaloFastAltroFitv0::FastFit(TH1F* h, Double_t sig, Double_t tau, Double_t n,
Double_t ped, Double_t tMax)
{
  // Service method for convinience only
  // h - hist with altro response, could have empty bin
  // and center of bin could be different from bin number  
  Reset();

  if(h==0) return;
  Int_t nPoints = h->GetNbinsX();
  if(nPoints<2) return; // Sep 07, 09

  Int_t* t = new Int_t[nPoints];
  Int_t* y = new Int_t[nPoints];

  Int_t nPositive=0;
  for(Int_t i=0; i<nPoints; i++) {
    if(h->GetBinContent(i+1) > 0.0){ // Get only positive
      y[nPositive] = Int_t(h->GetBinContent(i+1));
      t[nPositive] = Int_t(h->GetBinCenter(i+1)+0.0001);
      nPositive++;
    }
  }

  if(nPositive >= 2) {
    FastFit(t,y,nPoints, sig,tau,n,ped, tMax);
  }
  if(fChi2<=0.0) fNoFit++;

  delete [] t;
  delete [] y;
}

void AliCaloFastAltroFitv0::Reset()
{
  // Reset variables
  fSig  = fTau = 0.0;
  fAmp  = fAmpErr = fT0 = fT0Err = 0.0;
  fChi2 = -1.;
  fNDF  = fNfit = 0;

  if(fTfit)   delete [] fTfit;
  if(fAmpfit) delete [] fAmpfit;
  fTfit = fAmpfit = 0;
}


void AliCaloFastAltroFitv0::GetFitResult(Double_t &amp,Double_t &eamp,Double_t &t0,Double_t &et0, 
Double_t &chi2, Int_t &ndf) const
{
  // Return results of fitting procedure
  amp  = fAmp;
  eamp = fAmpErr;
  t0   = fT0;
  et0  = fT0Err;
  chi2 = fChi2;
  ndf  = fNDF;
}

void AliCaloFastAltroFitv0::GetFittedPoints(Int_t &nfit, Double_t* ar[2]) const
{
  nfit  = fNfit;
  ar[0] = fTfit;
  ar[1] = fAmpfit;
}
//
// static functions
// 
void  AliCaloFastAltroFitv0::CutRightPart(Int_t *t,Int_t *y,Int_t nPoints,Double_t tMax, Int_t &ii)
{
  // Cut right part of altro sample : static function
  Int_t tt=0;
  for(Int_t i=0; i<nPoints; i++) {
    tt = t[i];
    if(tMax && tt <= Int_t(tMax)) {
      t[ii] = tt;
      y[ii] = y[i];
      ii++;
    }
  }
  if(0) printf(" ii %i -> ii %i : tMax %7.2f \n", nPoints, ii, tMax);
}

void AliCaloFastAltroFitv0::DeductPedestal(Int_t* t, Int_t* y, Int_t nPoints, Double_t tau, Double_t ped, 
  Double_t* tn, Double_t* yn, Int_t &nPointsOut)
{
  // Pedestal deduction if ped is positive : static function
  // Discard part od samle if it is not compact.
  static Double_t yMinUnderPed=2.; // should be tune
  Int_t ymax=0, nmax=0;
  for(Int_t i=0; i<nPoints; i++){
    if(y[i]>ymax) {
      ymax = y[i];
      nmax = i;
    }
  }
  Int_t i1 = nmax - Int_t(tau);
  //i1 = 0;
  i1 = i1<0?0:i1;
  Int_t i2 = nPoints;

  nPointsOut = 0;
  Double_t yd=0.0, tdiff=0.0;;
  for(Int_t i=i1; i<i2; i++) {
    if(ped>0.0) {
      yd = Double_t(y[i]) - ped;
    } else {
      yd = Double_t(y[i]);
    }
    if(yd < yMinUnderPed) continue;

    if(i>i1 && nPointsOut>0){
      tdiff = t[i] - tn[nPointsOut-1];
      //      printf(" i %i : nPointsOut %i : tdiff %6.2f : tn[nPointsOut] %6.2f \n", i,nPointsOut, tdiff, tn[nPointsOut-1]);
      if(tdiff>1.) {
     // discard previous points if its are before maximum point and with gap>1
        if(i<nmax ) {
          nPointsOut = 0; // nPointsOut--;
     // if point with gap after maximum - finish selection
        } else if(i>=nmax ) {
          break;
        }
      }
     // Far away from maximum
     //if(i-nmax > Int_t(5*tau))              break;
    }
    tn[nPointsOut] = Double_t(t[i]);
    yn[nPointsOut] = yd;
    //printf("i %i : nPointsOut %i : tn %6.2f : yn %6.2f \n", i, nPointsOut, tn[nPointsOut], yn[nPointsOut]); 
    nPointsOut++;
  }
  //printf(" nmax %i : nPointsIn %i :nPointsOut %i i1 %i \n", nmax, nPointsIn, nPointsOut, i1);
}

void AliCaloFastAltroFitv0::FastFit(const Double_t* t, const Double_t* y, const Int_t nPoints, 
                                    const Double_t sig, const Double_t tau,
                                    Double_t &amp, Double_t &eamp, Double_t &t0, Double_t &et0, Double_t &chi2)
{
  // Static function
  // It is case of n=k=2 : fnn = x*x*exp(2 - 2*x)
  // Input: 
  //nPoints  - number of points 
  //   t[]   - array of time bins
  //   y[]   - array of amplitudes after pedestal subtractions;
  //   sig   - error of amplitude measurement (one value for all channels)
  //   tau   - filter time response (in timebin units)
  // Output:
  //       amp - amplitude at t0;
  //      eamp - error of amplitude; 
  //        t0 - time of max amplitude; 
  //       et0 - error of t0;
  //      chi2 - chi2
  static Double_t xx; // t/tau
  static Double_t a, b, c;
  static Double_t f02, f12, f22;    // functions
  static Double_t f02d, f12d, f22d; // functions derivations

  chi2 = -1.;

  if(nPoints<2) {
    printf(" AliCaloFastAltroFitv0::FastFit : nPoints<=%i \n", nPoints); 
    return;
  }

  a = b = c = 0.0;
  for(Int_t i=0; i<nPoints; i++){
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
  Double_t t01=0.0, t02=0.0;
  Double_t amp1=0.0, amp2=0.0, chi21=0.0, chi22=0.0;
  if(QuadraticRoots(a,b,c, t01,t02)) {
    t01 *= tau;
    t02 *= tau;
    Amplitude(t,y,nPoints, sig, tau, t01, amp1, chi21);
    Amplitude(t,y,nPoints, sig, tau, t02, amp2, chi22);
    if(0) {
      printf(" t01 %f : t02 %f \n", t01, t02);
      printf(" amp1 %f : amp2 %f \n", amp1, amp2);
      printf(" chi21 %f : chi22 %f \n", chi21, chi22);
    }
    // t0 less on one tau with comparing with value from  "canonical equation"
    amp  = amp1;
    t0   = t01;
    chi2 = chi21; 
    if(chi21 > chi22) {
      amp  = amp2;
      t0   = t02; 
      chi2 = chi22; 
    }
    if(tau<3.) { // EMCAL case : small tau 
      t0 += -0.03; // Discard bias in t0
      Amplitude(t,y,nPoints, sig, tau, t0, amp, chi2);
    }
    CalculateParsErrors(t, y, nPoints, sig, tau, amp, t0, eamp, et0);

    // Fill1();
    
    // DrawFastFunction(amp, t0, fUtils->GetPedestalValue(), "1");
    //    DrawFastFunction(amp1, t01, fUtils->GetPedestalValue(), "1");
    // DrawFastFunction(amp2, t02, fUtils->GetPedestalValue(), "2");
  } else {
    chi2 = t01; // no roots, bad fit - negative chi2
  }
}

Bool_t AliCaloFastAltroFitv0::QuadraticRoots(const Double_t a, const Double_t b, const Double_t c, 
                                             Double_t &x1, Double_t &x2)
{
  // Resolve quadratic equations a*x**2 + b*x + c
  //printf(" a %12.5e b %12.5e c %12.5e \n", a, b, c);
  static Double_t dtmp = 0.0, dtmpCut = -1.e-6;
  static Int_t iWarning=0, iNoSolution=0;
  dtmp = b*b - 4.*a*c;

  if(dtmp>=dtmpCut && dtmp<0.0) {
    if(iWarning<5 || iWarning%1000==0)
      printf("<W> %i small neg. sq. : dtmp %12.5e \n", iWarning, dtmp);
    iWarning++;
    dtmp = 0.0;
  }
  if(dtmp>=0.0) {
    dtmp = sqrt(dtmp);
    x1   = (-b + dtmp) / (2.*a);
    x2   = (-b - dtmp) / (2.*a);

    //    printf(" x1 %f : x2 %f \n", x1, x2);
    return kTRUE;
  } else {
    x1 = dtmp;
    if(iNoSolution<5 || iNoSolution%1000==0)
      printf("<No solution> %i neg. sq. : dtmp %12.5e \n", iNoSolution, dtmp);
    iNoSolution++;
    return kFALSE;
  }
}

void AliCaloFastAltroFitv0::Amplitude(const Double_t* t,const Double_t* y,const Int_t nPoints, 
                                      const Double_t sig, const Double_t tau, const Double_t t0, 
                                      Double_t &amp, Double_t &chi2)
{  
  // Calculate parameters error too - Mar 24,09
  // sig is independent from points
  amp = 0.;
  Double_t x=0.0, f=0.0, den=0.0, f02;
  for(Int_t i=0; i<nPoints; i++){
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
  for(Int_t i=0; i<nPoints; i++){
    x    = (t[i] - t0)/tau;
    f02  = exp(-2.*x);
    f    = amp*x*x*f02;
    dy   = y[i]-f;
    chi2 += dy*dy;
    //    printf(" %i : y %f -> f %f : dy %f \n", i, y[i], f, dy); 
  }
  chi2 /= (sig*sig);
}

void AliCaloFastAltroFitv0::CalculateParsErrors(const Double_t* t, const Double_t* /*y*/, const Int_t nPoints, 
                                                const Double_t sig, const Double_t tau, 
                                                Double_t &amp, Double_t &t0, Double_t &eamp, Double_t &et0)
{
  // Remember that fmax = exp(-n);
  // fmax_nk = (n/k)**n*exp(-n) => n=k=2 => exp(-n) = exp(-2.)
  static Double_t cc = exp(-2.);
  //   static Double_t cc = exp(-fN); // mean(N)~1.5 ??

  Double_t sumf2=0.0, sumfd2=0.0, x, f02, f12, f22, f22d;

  for(Int_t i=0; i<nPoints; i++){
    x    = (t[i] - t0)/tau;
    f02  = exp(-2.*x);
    f12  = x*f02;
    f22  = x*f12;
    sumf2 += f22 * f22;
    //
    f22d = 2.*(f12 - f22); 
    sumfd2 += f22d * f22d;
  }
  et0  = (sig/amp)/sqrt(sumfd2);
  eamp = sig/sqrt(sumf2);

  amp  *= cc;
  eamp *= cc;
}

//
// Drawing
//
TCanvas* AliCaloFastAltroFitv0::DrawFastFunction()
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

Double_t AliCaloFastAltroFitv0::StdResponseFunction(const Double_t *x, const Double_t *par)
{
  // Static Standard Response Function : 
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
