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
Revision 1.7  2000/05/02 08:12:13  morsch
Coding rule violations corrected.

Revision 1.6  1999/09/29 09:24:14  fca
Introduction of the Copyright and cvs Log

*/

#include "AliGenMUONlib.h"
#include "AliMC.h"
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
  const Double_t ka    = 7000.;
  const Double_t kdy   = 4.;

  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpPion()
{
// Pion composition 
    Float_t random[1];
    gMC->Rndm(random,1);
    if (random[0] < 0.5) {
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
  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;
  

  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

//                 particle composition
//
Int_t AliGenMUONlib::IpKaon()
{
// Kaon composition
    Float_t random[1];
    gMC->Rndm(random,1);
    if (random[0] < 0.5) {
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
//                 particle composition
//
Int_t AliGenMUONlib::IpJpsi()
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
//                 particle composition
//
Int_t AliGenMUONlib::IpUpsilon()
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
Int_t AliGenMUONlib::IpPhi()
{
// Phi composition
    return 41;
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
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}
//                  y-distribution
Double_t AliGenMUONlib::YCharm( Double_t *px, Double_t *dummy)
{
// Charm y
    Double_t *dum=0;
    return YJpsi(px,dum);
}

Int_t AliGenMUONlib::IpCharm()
{  
// Charm composition
    Float_t random[2];
    Int_t ip;
//    411,421,431,4122
    gMC->Rndm(random,2);
    if (random[0] < 0.5) {
	ip=411;
    } else if (random[0] < 0.75) {
	ip=421;
    } else if (random[0] < 0.90) {
	ip=431;
    } else {
	ip=4122;
    }
    if (random[1] < 0.5) {ip=-ip;}
    
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

Int_t AliGenMUONlib::IpBeauty()
{  
// Beauty Composition
    Float_t random[2];
    Int_t ip;
    gMC->Rndm(random,2);
    if (random[0] < 0.5) {
	ip=511;
    } else if (random[0] < 0.75) {
	ip=521;
    } else if (random[0] < 0.90) {
	ip=531;
    } else {
	ip=5122;
    }
    if (random[1] < 0.5) {ip=-ip;}
    
    return ip;
}

typedef Double_t (*GenFunc) (Double_t*,  Double_t*);
GenFunc AliGenMUONlib::GetPt(Param_t param)
{
// Return pointer to pT parameterisation
    GenFunc func;
    switch (param) 
    {
    case phi_p:
	func=PtPhi;
	break;
    case jpsi_p:
	func=PtJpsi;
	break;
    case upsilon_p:
	func=PtUpsilon;
	break;
    case charm_p:
	func=PtCharm;
	break;
    case beauty_p:
	func=PtBeauty;
	break;
    case pion_p:
	func=PtPion;
	break;
    case kaon_p:
	func=PtKaon;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetPt> unknown parametrisation\n");
    }
    return func;
}

GenFunc AliGenMUONlib::GetY(Param_t param)
{
// Return pointer to y- parameterisation
    GenFunc func;
    switch (param) 
    {
    case phi_p:
	func=YPhi;
	break;
    case jpsi_p:
	func=YJpsi;
	break;
    case upsilon_p:
	func=YUpsilon;
	break;
    case charm_p:
	func=YCharm;
	break;
    case beauty_p:
	func=YBeauty;
	break;
    case pion_p:
	func=YPion;
	break;
    case kaon_p:
	func=YKaon;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetY> unknown parametrisation\n");
    }
    return func;
}
typedef Int_t (*GenFuncIp) ();
GenFuncIp AliGenMUONlib::GetIp(Param_t param)
{
// Return pointer to particle type parameterisation
    GenFuncIp func;
    switch (param) 
    {
    case phi_p:
	func=IpPhi;
	break;
    case jpsi_p:
	func=IpJpsi;
	break;
    case upsilon_p:
	func=IpUpsilon;
	break;
    case charm_p:
	func=IpCharm;
	break;
    case beauty_p:
	func=IpBeauty;
	break;
    case pion_p:
	func=IpPion;
	break;
    case kaon_p:
	func=IpKaon;
	break;
    default:
        func=0;
        printf("<AliGenMUONlib::GetIp> unknown parametrisation\n");
    }
    return func;
}




