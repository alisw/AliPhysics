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

/*
$Log$
Revision 1.18  2002/11/07 09:13:42  morsch
Use "Vogt" to label new distributions.

Revision 1.17  2002/11/07 09:06:10  morsch
J/Psi and Upsilon pt and y distributions from R. Vogt 2002 added.

Revision 1.16  2002/10/14 14:55:35  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.14.6.1  2002/06/10 14:57:41  hristov
Merged with v3-08-02

Revision 1.15  2002/04/17 10:11:51  morsch
Coding Rule violations corrected.

Revision 1.14  2002/02/22 17:26:43  morsch
Eta and omega added.

Revision 1.13  2001/03/27 11:01:04  morsch
Charm pt-distribution corrected. More realistic y-distribution for pi and K.

Revision 1.12  2001/03/09 13:01:41  morsch
- enum constants for paramterisation type (particle family) moved to AliGen*lib.h
- use AliGenGSIlib::kUpsilon, AliGenPHOSlib::kEtaPrime to access the constants

Revision 1.11  2000/11/30 07:12:50  alibrary
Introducing new Rndm and QA classes

Revision 1.10  2000/06/29 21:08:27  morsch
All paramatrisation libraries derive from the pure virtual base class AliGenLib.
This allows to pass a pointer to a library directly to AliGenParam and avoids the
use of function pointers in Config.C.

Revision 1.9  2000/06/14 15:20:56  morsch
Include clean-up (IH)

Revision 1.8  2000/06/09 20:32:11  morsch
All coding rule violations except RS3 corrected

Revision 1.7  2000/05/02 08:12:13  morsch
Coding rule violations corrected.

Revision 1.6  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/

// Library class for particle pt and y distributions used for 
// muon spectrometer simulations.
// To be used with AliGenParam.
// The following particle typed can be simulated:
// pi, K, phi, omega, eta, J/Psi, Upsilon, charm and beauty mesons. 
//
// andreas.morsch@cern.ch
//

#include "TMath.h"
#include "TRandom.h"

#include "AliGenMUONlib.h"

ClassImp(AliGenMUONlib)
//
//  Pions
Double_t AliGenMUONlib::PtPion(Double_t *px, Double_t *dummy)
{
//
//     PT-PARAMETERIZATION CDF, PRL 61(88) 1819
//     POWER LAW FOR PT > 500 MEV
//     MT SCALING BELOW (T=160 MEV)
//
  const Double_t kp0 = 1.3;
  const Double_t kxn = 8.28;
  const Double_t kxlim=0.5;
  const Double_t kt=0.160;
  const Double_t kxmpi=0.139;
  const Double_t kb=1.;
  Double_t y, y1, xmpi2, ynorm, a;
  Double_t x=*px;
  //
  y1=TMath::Power(kp0/(kp0+kxlim),kxn);
  xmpi2=kxmpi*kxmpi;
  ynorm=kb*(TMath::Exp(-sqrt(kxlim*kxlim+xmpi2)/kt));
  a=ynorm/y1;
  if (x > kxlim)
    y=a*TMath::Power(kp0/(kp0+x),kxn);
  else
    y=kb*TMath::Exp(-sqrt(x*x+xmpi2)/kt);
  return y*x;
}
//
// y-distribution
//
Double_t AliGenMUONlib::YPion( Double_t *py, Double_t *dummy)
{
// Pion y
  Double_t y=TMath::Abs(*py);
/*
  const Double_t ka    = 7000.;
  const Double_t kdy   = 4.;
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
*/
  return 1.16526e+04+y*-3.79886e+03+y*y*4.31130e+02;
  
}
//                 particle composition
//
Int_t AliGenMUONlib::IpPion(TRandom *ran)
{
// Pion composition 
    if (ran->Rndm() < 0.5) {
	return  211;
    } else {
	return -211;
    }
}

//____________________________________________________________
//
// Mt-scaling

Double_t AliGenMUONlib::PtScal(Double_t pt, Int_t np)
{
  //    SCALING EN MASSE PAR RAPPORT A PTPI
  //    MASS PI,K,ETA,RHO,OMEGA,ETA',PHI
  const Double_t khm[10] = {.13957,.493,.5488,.769,.7826,.958,1.02,0,0,0};
  //     VALUE MESON/PI AT 5 GEV
  const Double_t kfmax[10]={1.,0.3,0.55,1.0,1.0,1.0,1.0,0,0,0};
  np--;
  Double_t f5=TMath::Power(((sqrt(100.018215)+2.)/(sqrt(100.+khm[np]*khm[np])+2.0)),12.3);
  Double_t fmax2=f5/kfmax[np];
  // PIONS
  Double_t ptpion=100.*PtPion(&pt, (Double_t*) 0);
  Double_t fmtscal=TMath::Power(((sqrt(pt*pt+0.018215)+2.)/
				 (sqrt(pt*pt+khm[np]*khm[np])+2.0)),12.3)/ fmax2;
  return fmtscal*ptpion;
}
//
// kaon
//
//                pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtKaon( Double_t *px, Double_t *dummy)
{
// Kaon pT
  return PtScal(*px,2);
}

// y-distribution
//____________________________________________________________
Double_t AliGenMUONlib::YKaon( Double_t *py, Double_t *dummy)
{
// Kaon y
  Double_t y=TMath::Abs(*py);
/*
  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
*/

  return 1.16526e+04+y*-3.79886e+03+y*y*4.31130e+02;
}

//                 particle composition
//
Int_t AliGenMUONlib::IpKaon(TRandom *ran)
{
// Kaon composition
    if (ran->Rndm() < 0.5) {
	return  321;
    } else {
	return -321;
    }
}

//                    J/Psi 
//
//
//                pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtJpsi( Double_t *px, Double_t *dummy)
{
// J/Psi pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtJpsiPbPb( Double_t *px, Double_t *dummy)
{
// J/Psi pT spectrum
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Float_t x = px[0];
    Float_t c[8] = {
	-2.13098e+00, 9.46552e+00, -5.06799e+00, 1.27260e+00, 
	-1.83806e-01, 1.55853e-02, -7.23241e-04, 1.42105e-05
    };
    
    Int_t j;
    Double_t y = c[j = 7];
    while (j > 0) y  = y * x +c[--j];
    y = x * TMath::Exp(y);
    return y;
}
//
//               y-distribution
//____________________________________________________________
Double_t AliGenMUONlib::YJpsi(Double_t *py, Double_t *dummy)
{
// J/psi y
  const Double_t ky0 = 4.;
  const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;
}


Double_t AliGenMUONlib::YJpsiPbPb( Double_t *px, Double_t *dummy)
{

//
// J/Psi y
//
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Double_t c[5] = {-6.03425e+02, 4.98257e+02, -1.38794e+02, 1.62209e+01, -6.85955e-01};
    Double_t x = TMath::Abs(px[0]);
    Double_t y;
    
    if (x < 4.) {
	y = 31.754;
    } else if (x < 6) {
	Int_t j;
	y = c[j = 4];
	while (j > 0) y  = y * x + c[--j];
    } else {
	y =0.;
    }
    
    return y;
}

//                 particle composition
//
Int_t AliGenMUONlib::IpJpsi(TRandom *)
{
// J/Psi composition
    return 443;
}

//                      Upsilon
//
//
//                  pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtUpsilon( Double_t *px, Double_t *dummy )
{
// Upsilon pT
  const Double_t kpt0 = 5.3;
  const Double_t kxn  = 2.5;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

Double_t AliGenMUONlib::PtUpsilonPbPb( Double_t *px, Double_t *dummy)
{

//
// Upsilon pT
//
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//
    Float_t x = px[0];
    Double_t c[8] = {
	-1.03488e+01, 1.28065e+01, -6.60500e+00, 1.66140e+00,       
	-2.34293e-01, 1.86925e-02, -7.80708e-04, 1.30610e-05
    };
    Int_t j;
    Double_t y = c[j = 7];
    while (j > 0) y  = y * x +c[--j];
    y = x * TMath::Exp(y);
    return y;
}

//
//                    y-distribution
//
//____________________________________________________________
Double_t AliGenMUONlib::YUpsilon(Double_t *py, Double_t *dummy)
{
// Upsilon y
  const Double_t ky0 = 3.;
  const Double_t kb=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yu=kb;
  else
    yu=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yu;
}


Double_t AliGenMUONlib::YUpsilonPbPb( Double_t *px, Double_t *dummy)
{

//
// Upsilon y
//
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//

    Double_t c[7] = {3.40036e-01, -3.98882e-07, -4.48398e-03, 8.46411e-08, -6.10854e-04,
		     -2.99753e-09, 1.28895e-05};
        
    Double_t x = px[0];
    if (TMath::Abs(x) > 5.55) return 0.;
    Int_t j;
    Double_t y = c[j = 6];
    while (j > 0) y  = y * x +c[--j];
    return y;
}

//                 particle composition
//
Int_t AliGenMUONlib::IpUpsilon(TRandom *)
{
// y composition
    return 553;
}

//
//                        Phi
//
//
//    pt-distribution (by scaling of pion distribution)
//____________________________________________________________
Double_t AliGenMUONlib::PtPhi( Double_t *px, Double_t *dummy)
{
// Phi pT
  return PtScal(*px,7);
}
//    y-distribution
Double_t AliGenMUONlib::YPhi( Double_t *px, Double_t *dummy)
{
// Phi y
    Double_t *dum=0;
    return YJpsi(px,dum);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpPhi(TRandom *)
{
// Phi composition
    return 333;
}

//
//                        omega
//
//
//    pt-distribution (by scaling of pion distribution)
//____________________________________________________________
Double_t AliGenMUONlib::PtOmega( Double_t *px, Double_t *dummy)
{
// Omega pT
  return PtScal(*px,5);
}
//    y-distribution
Double_t AliGenMUONlib::YOmega( Double_t *px, Double_t *dummy)
{
// Omega y
    Double_t *dum=0;
    return YJpsi(px,dum);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpOmega(TRandom *)
{
// Omega composition
    return 223;
}


//
//                        Eta
//
//
//    pt-distribution (by scaling of pion distribution)
//____________________________________________________________
Double_t AliGenMUONlib::PtEta( Double_t *px, Double_t *dummy)
{
// Eta pT
  return PtScal(*px,3);
}
//    y-distribution
Double_t AliGenMUONlib::YEta( Double_t *px, Double_t *dummy)
{
// Eta y
    Double_t *dum=0;
    return YJpsi(px,dum);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpEta(TRandom *)
{
// Eta composition
    return 221;
}

//
//                        Charm
//
//
//                    pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtCharm( Double_t *px, Double_t *dummy)
{
// Charm pT
  const Double_t kpt0 = 4.08;
  const Double_t kxn  = 9.40;

  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
//                  y-distribution
Double_t AliGenMUONlib::YCharm( Double_t *px, Double_t *dummy)
{
// Charm y
    Double_t *dum=0;
    return YJpsi(px,dum);
}

Int_t AliGenMUONlib::IpCharm(TRandom *ran)
{  
// Charm composition
    Float_t random;
    Int_t ip;
//    411,421,431,4122
    random = ran->Rndm();
    if (random < 0.5) {
	ip=411;
    } else if (random < 0.75) {
	ip=421;
    } else if (random < 0.90) {
	ip=431;
    } else {
	ip=4122;
    }
    if (ran->Rndm() < 0.5) {ip=-ip;}
    
    return ip;
}


//
//                        Beauty
//
//
//                    pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtBeauty( Double_t *px, Double_t *dummy)
{
// Beauty pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
//                     y-distribution
Double_t AliGenMUONlib::YBeauty( Double_t *px, Double_t *dummy)
{
// Beauty y
    Double_t *dum=0;
    return YJpsi(px,dum);
}

Int_t AliGenMUONlib::IpBeauty(TRandom *ran)
{  
// Beauty Composition
    Float_t random;
    Int_t ip;
    random = ran->Rndm();
    if (random < 0.5) {
	ip=511;
    } else if (random < 0.75) {
	ip=521;
    } else if (random < 0.90) {
	ip=531;
    } else {
	ip=5122;
    }
    if (ran->Rndm() < 0.5) {ip=-ip;}
    
    return ip;
}

typedef Double_t (*GenFunc) (Double_t*,  Double_t*);
GenFunc AliGenMUONlib::GetPt(Int_t param,  const char* tname) const
{
// Return pointer to pT parameterisation
    TString sname = TString(tname);
    GenFunc func;
    switch (param) 
    {
    case kPhi:
	func=PtPhi;
	break;
    case kOmega:
	func=PtOmega;
	break;
    case kEta:
	func=PtEta;
	break;
    case kJpsi:
	if (sname == "Vogt") {
	    func=PtJpsiPbPb;
	} else {
	    func=PtJpsi;
	}
	break;
    case kUpsilon:
	if (sname == "Vogt") {
	    func=PtUpsilonPbPb;
	} else {
	    func=PtUpsilon;
	}
	break;
    case kCharm:
	func=PtCharm;
	break;
    case kBeauty:
	func=PtBeauty;
	break;
    case kPion:
	func=PtPion;
	break;
    case kKaon:
	func=PtKaon;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetPt> unknown parametrisation\n");
    }
    return func;
}

GenFunc AliGenMUONlib::GetY(Int_t param, const char* tname) const
{
    TString sname = TString(tname);
    
// Return pointer to y- parameterisation
    GenFunc func;
    switch (param) 
    {
    case kPhi:
	func=YPhi;
	break;
    case kEta:
	func=YEta;
	break;
    case kOmega:
	func=YOmega;
	break;
    case kJpsi:
	if (sname == "Vogt") {
	    func=YJpsiPbPb;
	} else {
	    func=YJpsi;
	}
	break;
    case kUpsilon:
	if (sname == "Vogt") {
	    func=YUpsilonPbPb;
	} else {
	    func=YUpsilon;
	}
	break;
    case kCharm:
	func=YCharm;
	break;
    case kBeauty:
	func=YBeauty;
	break;
    case kPion:
	func=YPion;
	break;
    case kKaon:
	func=YKaon;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetY> unknown parametrisation\n");
    }
    return func;
}
typedef Int_t (*GenFuncIp) (TRandom *);
GenFuncIp AliGenMUONlib::GetIp(Int_t param,  const char* tname) const
{
// Return pointer to particle type parameterisation
    GenFuncIp func;
    switch (param) 
    {
    case kPhi:
	func=IpPhi;
	break;
    case kEta:
	func=IpEta;
	break;
    case kOmega:
	func=IpOmega;
	break;
    case kJpsi:
	func=IpJpsi;
	break;
    case kUpsilon:
	func=IpUpsilon;
	break;
    case kCharm:
	func=IpCharm;
	break;
    case kBeauty:
	func=IpBeauty;
	break;
    case kPion:
	func=IpPion;
	break;
    case kKaon:
	func=IpKaon;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetIp> unknown parametrisation\n");
    }
    return func;
}



Float_t AliGenMUONlib::Interpolate(Float_t x, Float_t* y, Float_t x0, 
				   Float_t dx,
				   Int_t n, Int_t no)
{
//
// Neville's alorithm for interpolation
//
// x:  x-value
// y:  Input array
// x0: minimum x 
// dx: step size
//  n: number of data points
// no: order of polynom 
//
    Float_t*  c = new Float_t[n];
    Float_t*  d = new Float_t[n];
    Int_t m, i;
    for (i = 0; i < n; i++) {
	c[i] = y[i];
	d[i] = y[i];
    }
    
    Int_t   ns  = int((x - x0)/dx);
    
    Float_t y1  = y[ns];
    ns--;    
    for (m = 0; m < no; m++) {
	for (i = 0; i < n-m; i++) {	
	    Float_t ho = x0 + Float_t(i) * dx - x;
	    Float_t hp = x0 + Float_t(i+m+1) * dx - x;
	    Float_t w  = c[i+1] - d[i];
	    Float_t den = ho-hp;
	    den = w/den;
	    d[i] = hp * den;
	    c[i] = ho * den;
	}
	Float_t dy;
	
	if (2*ns < (n-m-1)) {
	    dy  = c[ns+1];
	} else {
	    dy  = d[ns--];
	}
	y1 += dy;}
    delete[] c;
    delete[] d;

    return y1;
}


