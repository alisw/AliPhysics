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
// J/Psi pT
//
// R. Vogt 2002
// PbPb 5.5 TeV
// MRST HO
// mc = 1.4 GeV, pt-kick 1 GeV
//

    Float_t ptJpsi[100] = {
        0.0000e-01,  4.5870e+01,  6.5200e+01,  7.1740e+01,  6.5090e+01,
	5.5070e+01,  4.9420e+01,  3.9780e+01,  3.2390e+01,  2.8120e+01,
	2.3870e+01,  1.9540e+01,  1.6510e+01,  1.4180e+01,  1.2050e+01,
	1.0390e+01,  8.7970e+00,  7.8680e+00,  6.7710e+00,  5.9360e+00,
	5.3460e+00,  4.5670e+00,  4.6500e+00,  3.9360e+00,  3.5070e+00,
	3.2070e+00,  2.8310e+00,  2.6340e+00,  2.4900e+00,  2.2410e+00,
	2.1090e+00,  1.9070e+00,  1.7360e+00,  1.6120e+00,  1.5450e+00,
	1.4350e+00,  1.3890e+00,  1.2610e+00,  1.0880e+00,  1.0930e+00,
	1.0680e+00,  9.2500e-01,  8.6790e-01,  8.1790e-01,  7.9770e-01,
	7.4660e-01,  7.3110e-01,  6.5120e-01,  6.8140e-01,  5.7960e-01,
	5.8210e-01,  5.4640e-01,  5.1700e-01,  5.0760e-01,  4.8280e-01,
	4.5360e-01,  4.4910e-01,  4.2410e-01,  4.2100e-01,  3.9530e-01,
	3.7220e-01,  3.4840e-01,  3.4550e-01,  3.3000e-01,  3.1670e-01,
	3.1470e-01,  2.8920e-01,  2.7650e-01,  2.6860e-01,  2.5390e-01,
	2.4190e-01,  2.5200e-01,  2.2960e-01,  2.2540e-01,  2.0950e-01,
	2.0250e-01,  1.8720e-01,  1.8200e-01,  1.7860e-01,  1.8290e-01,
	1.6970e-01,  1.7130e-01,  1.6310e-01,  1.5500e-01,  1.5100e-01,
	1.5770e-01,  1.4240e-01,  1.4560e-01,  1.3330e-01,  1.4190e-01,
	1.2010e-01,  1.2430e-01,  1.2430e-01,  1.1340e-01,  1.1840e-01,
	1.1380e-01,  1.0330e-01,  1.0130e-01,  1.0390e-01,  9.5810e-02
    };
    Float_t x = px[0] * px[0];
    if (x  < 1.5 || x > 100) {
	return 0.0;
    } else {
	Float_t y =  Interpolate(x, ptJpsi, 0.5, 1., 100, 4);
	return px[0] * y;
    }
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

    Float_t yJpsi[62] = 
	{
	  0.3981E-03, 0.1169E-01, 0.6143E-01, 0.3554E+00, 0.1249E+01 
	, 0.1677E+01, 0.3634E+01, 0.5414E+01, 0.8242E+01, 0.1102E+02
	, 0.1353E+02, 0.1964E+02, 0.2357E+02, 0.2662E+02, 0.3023E+02 
	, 0.3250E+02, 0.3137E+02, 0.3243E+02, 0.3120E+02, 0.3249E+02 
	, 0.3166E+02, 0.3104E+02, 0.3203E+02, 0.3149E+02, 0.3117E+02 
	, 0.3210E+02, 0.3170E+02, 0.3279E+02, 0.3079E+02, 0.3208E+02 
	, 0.3218E+02, 0.3218E+02, 0.3208E+02, 0.3079E+02, 0.3279E+02 
	, 0.3170E+02, 0.3210E+02, 0.3118E+02, 0.3149E+02, 0.3203E+02 
	, 0.3104E+02, 0.3167E+02, 0.3250E+02, 0.3120E+02, 0.3243E+02 
	, 0.3137E+02, 0.3250E+02, 0.3022E+02, 0.2662E+02, 0.2357E+02 
	, 0.1964E+02, 0.1353E+02, 0.1102E+02, 0.8242E+01, 0.5414E+01 
	, 0.3634E+01, 0.1677E+01, 0.1249E+01, 0.3554E+00, 0.6142E-01
        , 0.1169E-01, 0.3981E-03};

    return Interpolate(px[0], yJpsi, -7.625, 0.25, 62, 2);
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

    Float_t ptUps[100] = 
	{   
	    0.0000e-01, -1.5290e-02,  2.6020e-01,  3.4220e-01,  3.3710e-01,
	    3.1880e-01,  3.3420e-01,  2.7740e-01,  2.3730e-01,  2.0640e-01,
	    1.7690e-01,  1.6190e-01,  1.4500e-01,  1.3310e-01,  1.1440e-01,
	    1.0800e-01,  1.0210e-01,  8.4690e-02,  8.0050e-02,  7.0710e-02,
	    6.4160e-02,  6.5200e-02,  6.6890e-02,  6.0600e-02,  5.4030e-02,
	    5.1140e-02,  4.6120e-02,  4.4800e-02,  4.2490e-02,  4.1440e-02,
	    4.0310e-02,  3.7110e-02,  3.5890e-02,  3.5420e-02,  3.0370e-02,
	    2.9970e-02,  3.0770e-02,  2.6380e-02,  2.7740e-02,  2.6690e-02,
	    2.4210e-02,  2.5200e-02,  2.3760e-02,  2.1370e-02,  2.2290e-02,
	    2.2700e-02,  2.0110e-02,  1.9320e-02,  1.8830e-02,  1.9910e-02,
	    1.9740e-02,  1.8460e-02,  1.8240e-02,  1.6740e-02,  1.6140e-02,
	    1.7340e-02,  1.5950e-02,  1.5430e-02,  1.4780e-02,  1.2750e-02,
	    1.4370e-02,  1.2810e-02,  1.2900e-02,  1.1070e-02,  1.1830e-02,
	    1.1150e-02,  1.1260e-02,  1.1610e-02,  1.0700e-02,  1.1600e-02,
	    1.0390e-02,  1.0280e-02,  1.0180e-02,  1.0030e-02,  9.6050e-03,
	    8.8050e-03,  8.9680e-03,  9.0120e-03,  8.4110e-03,  8.6660e-03,
	    8.3060e-03,  8.5850e-03,  8.2600e-03,  8.3800e-03,  8.4200e-03,
	    7.5690e-03,  7.2100e-03,  7.1230e-03,  7.3350e-03,  7.1980e-03,
	    6.7500e-03,  6.6190e-03,  6.3370e-03,  6.6270e-03,  6.8290e-03,
	    6.0880e-03,  6.6310e-03,  6.0490e-03,  5.8900e-03,  5.6100e-03
	};
    Float_t x = px[0] * px[0];
    if (x  < 1.5 || x > 100) {
	return 0.0;
    } else {
	Float_t y =  Interpolate(x, ptUps, 0.5, 1., 100, 4);
	return px[0] * y;
    }
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

    Float_t yUps[52] = 
	{
	    0.000000, 0.000065, 0.000767, 0.004332, 0.014790,
	    0.029370, 0.052060, 0.077930, 0.122500, 0.135800,
	    0.184000, 0.207900, 0.228300, 0.258500, 0.269500,
	    0.288500, 0.316200, 0.304100, 0.315800, 0.323300,
	    0.322400, 0.322600, 0.345500, 0.338100, 0.331900,
	    0.343700, 0.343700, 0.331900, 0.338100, 0.345500,
	    0.322600, 0.322400, 0.323300, 0.315800, 0.304100,
	    0.316200, 0.288500, 0.269500, 0.258500, 0.228300,
	    0.207900, 0.184000, 0.135800, 0.122500, 0.077930,
	    0.052060, 0.029380, 0.014780, 0.004332, 0.000767,
	    0.6479E-04, 0.1013E-06       
	};
    


    return Interpolate(px[0], yUps, -6.5, 0.25, 52, 2);
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
	if (sname == "PbPb") {
	    func=PtJpsiPbPb;
	} else {
	    func=PtJpsi;
	}
	break;
    case kUpsilon:
	if (sname == "PbPb") {
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
	if (sname == "PbPb") {
	    func=YJpsiPbPb;
	} else {
	    func=YJpsi;
	}
	break;
    case kUpsilon:
	if (sname == "PbPb") {
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


