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
*/

#include "AliGenMUONlib.h"
#include "AliRun.h"
ClassImp(AliGenMUONlib)
//
//  Pions
Double_t AliGenMUONlib::PtPion(Double_t *px, Double_t *)
{
//
//     PT-PARAMETERIZATION CDF, PRL 61(88) 1819
//     POWER LAW FOR PT > 500 MEV
//     MT SCALING BELOW (T=160 MEV)
//
  const Double_t p0 = 1.3;
  const Double_t xn = 8.28;
  const Double_t xlim=0.5;
  const Double_t t=0.160;
  const Double_t xmpi=0.139;
  const Double_t b=1.;
  Double_t y, y1, xmpi2, ynorm, a;
  Double_t x=*px;
  //
  y1=TMath::Power(p0/(p0+xlim),xn);
  xmpi2=xmpi*xmpi;
  ynorm=b*(TMath::Exp(-sqrt(xlim*xlim+xmpi2)/t));
  a=ynorm/y1;
  if (x > xlim)
    y=a*TMath::Power(p0/(p0+x),xn);
  else
    y=b*TMath::Exp(-sqrt(x*x+xmpi2)/t);
  return y*x;
}
//
// y-distribution
//
Double_t AliGenMUONlib::YPion( Double_t *py, Double_t *)
{
  const Double_t a    = 7000.;
  const Double_t dy   = 4.;

  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*dy*dy);
  return a*TMath::Exp(-ex);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpPion()
{
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
  const Double_t hm[10] = {.13957,.493,.5488,.769,.7826,.958,1.02,0,0,0};
  //     VALUE MESON/PI AT 5 GEV
  const Double_t fmax[10]={1.,0.3,0.55,1.0,1.0,1.0,1.0,0,0,0};
  np--;
  Double_t f5=TMath::Power(((sqrt(100.018215)+2.)/(sqrt(100.+hm[np]*hm[np])+2.0)),12.3);
  Double_t fmax2=f5/fmax[np];
  // PIONS
  Double_t ptpion=100.*PtPion(&pt, (Double_t*) 0);
  Double_t fmtscal=TMath::Power(((sqrt(pt*pt+0.018215)+2.)/
				 (sqrt(pt*pt+hm[np]*hm[np])+2.0)),12.3)/ fmax2;
  return fmtscal*ptpion;
}
//
// kaon
//
//                pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtKaon( Double_t *px, Double_t *)
{
  return PtScal(*px,2);
}

// y-distribution
//____________________________________________________________
Double_t AliGenMUONlib::YKaon( Double_t *py, Double_t *)
{
  const Double_t a    = 1000.;
  const Double_t dy   = 4.;
  

  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*dy*dy);
  return a*TMath::Exp(-ex);
}

//                 particle composition
//
Int_t AliGenMUONlib::IpKaon()
{
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
Double_t AliGenMUONlib::PtJpsi( Double_t *px, Double_t *)
{
  const Double_t pt0 = 4.;
  const Double_t xn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/pt0)*(x/pt0);
  return x/TMath::Power(pass1,xn);
}
//
//               y-distribution
//____________________________________________________________
Double_t AliGenMUONlib::YJpsi(Double_t *py, Double_t *)
{
  const Double_t y0 = 4.;
  const Double_t b=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < y0)
    yj=b;
  else
    yj=b*TMath::Exp(-(y-y0)*(y-y0)/2);
  return yj;
}
//                 particle composition
//
Int_t AliGenMUONlib::IpJpsi()
{
    return 443;
}

//                      Upsilon
//
//
//                  pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtUpsilon( Double_t *px, Double_t * )
{
  const Double_t pt0 = 5.3;
  const Double_t xn  = 2.5;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/pt0)*(x/pt0);
  return x/TMath::Power(pass1,xn);
}
//
//                    y-distribution
//
//____________________________________________________________
Double_t AliGenMUONlib::YUpsilon(Double_t *py, Double_t *)
{
  const Double_t y0 = 3.;
  const Double_t b=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);
  //
  if (y < y0)
    yu=b;
  else
    yu=b*TMath::Exp(-(y-y0)*(y-y0)/2);
  return yu;
}
//                 particle composition
//
Int_t AliGenMUONlib::IpUpsilon()
{
    return 553;
}

//
//                        Phi
//
//
//    pt-distribution (by scaling of pion distribution)
//____________________________________________________________
Double_t AliGenMUONlib::PtPhi( Double_t *px, Double_t *)
{
  return PtScal(*px,7);
}
//    y-distribution
Double_t AliGenMUONlib::YPhi( Double_t *px, Double_t *)
{
    Double_t *dummy=0;
    return YJpsi(px,dummy);
}
//                 particle composition
//
Int_t AliGenMUONlib::IpPhi()
{
    return 41;
}

//
//                        Charm
//
//
//                    pt-distribution
//____________________________________________________________
Double_t AliGenMUONlib::PtCharm( Double_t *px, Double_t *)
{
  const Double_t pt0 = 4.08;
  const Double_t xn  = 9.40;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/pt0)*(x/pt0);
  return x/TMath::Power(pass1,xn);
}
//                  y-distribution
Double_t AliGenMUONlib::YCharm( Double_t *px, Double_t *)
{
    Double_t *dummy=0;
    return YJpsi(px,dummy);
}

Int_t AliGenMUONlib::IpCharm()
{  
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
Double_t AliGenMUONlib::PtBeauty( Double_t *px, Double_t *)
{
  const Double_t pt0 = 4.;
  const Double_t xn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/pt0)*(x/pt0);
  return x/TMath::Power(pass1,xn);
}
//                     y-distribution
Double_t AliGenMUONlib::YBeauty( Double_t *px, Double_t *)
{
    Double_t *dummy=0;
    return YJpsi(px,dummy);
}

Int_t AliGenMUONlib::IpBeauty()
{  
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




