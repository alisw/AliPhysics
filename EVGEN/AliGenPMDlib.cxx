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

//
// Library of generators for PMD
// providing y and pt parametrisation
// for generated tracks
// Specific for PMD needs
// Author: PMD Offline Group
//

#include <TMath.h>
#include <TPDGCode.h>

#include "AliGenPMDlib.h"

ClassImp(AliGenPMDlib)
//
//  Neutral Pions

Double_t AliGenPMDlib::PtPi0(const Double_t *px, const Double_t */*dummy*/)
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
Double_t AliGenPMDlib::YPi0( const Double_t *py, const Double_t */*dummy*/)
{
  //
  // y parametrisation for pi0
  //
    const Double_t ka1    = 4913.;
    const Double_t ka2    = 1819.;
    const Double_t keta1  = 0.22;
    const Double_t keta2  = 3.66;
    const Double_t kdeta1 = 1.47;
    const Double_t kdeta2 = 1.51;
    Double_t y=TMath::Abs(*py);
    //
    Double_t ex1 = (y-keta1)*(y-keta1)/(2*kdeta1*kdeta1);
    Double_t ex2 = (y-keta2)*(y-keta2)/(2*kdeta2*kdeta2);
    return ka1*TMath::Exp(-ex1)+ka2*TMath::Exp(-ex2);
}

//                 particle composition
//
Int_t AliGenPMDlib::IpPi0(TRandom *)
{
// Pi0
    return kPi0;
}

//____________________________________________________________
//
// Mt-scaling

Double_t AliGenPMDlib::PtScal(Double_t pt, Int_t np)
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
  Double_t ptpion=100.*PtPi0(&pt, (Double_t*) 0);
  Double_t fmtscal=TMath::Power(((sqrt(pt*pt+0.018215)+2.)/
				 (sqrt(pt*pt+khm[np]*khm[np])+2.0)),12.3)/ fmax2;
  return fmtscal*ptpion;
}
//
// kaon
//
//                pt-distribution
//____________________________________________________________

Double_t AliGenPMDlib::PtEta( const Double_t *px, const Double_t */*dummy*/)
{
// Kaon pT
  return PtScal(*px,3);
}

// y-distribution
//____________________________________________________________
Double_t AliGenPMDlib::YEta( const Double_t *py, const Double_t */*dummy*/)
{
    //
    // y parametrisation for etas
    //
    const Double_t ka1    = 4913.;
    const Double_t ka2    = 1819.;
    const Double_t keta1  = 0.22;
    const Double_t keta2  = 3.66;
    const Double_t kdeta1 = 1.47;
    const Double_t kdeta2 = 1.51;
    Double_t y=TMath::Abs(*py);
    //
    Double_t ex1 = (y-keta1)*(y-keta1)/(2*kdeta1*kdeta1);
    Double_t ex2 = (y-keta2)*(y-keta2)/(2*kdeta2*kdeta2);
    return ka1*TMath::Exp(-ex1)+ka2*TMath::Exp(-ex2);
}

//                 particle composition
//
Int_t AliGenPMDlib::IpEta(TRandom *)
{
    return 221;
}


typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);
GenFunc AliGenPMDlib::GetPt(Int_t param,  const char* /*tname*/) const
{
// Return pointer to pT parameterisation
    GenFunc func=NULL;
    switch (param) 
    {
    case kPion:
	func=PtPi0;
	break;
    case kEta:
	func=PtEta;
	break;
    default:
        func=0;
        printf("<AliGenPMDlib::GetPt> unknown parametrisation\n");
    }
    return func;
}

GenFunc AliGenPMDlib::GetY(Int_t param, const char* /*tname*/) const
{
// Return pointer to y- parameterisation
    GenFunc func=NULL;
    switch (param) 
    {
    case kPion:
	func=YPi0;
	break;
    case kEta:
	func=YEta;
	break;
    default:
        func=0;
        printf("<AliGenPMDlib::GetY> unknown parametrisation\n");
    }
    return func;

}
typedef Int_t (*GenFuncIp) (TRandom *);
GenFuncIp AliGenPMDlib::GetIp(Int_t param,  const char* /*tname*/) const
{
// Return pointer to particle type parameterisation
    GenFuncIp func=NULL;
    switch (param) 
    {
    case kPion:
	func=IpPi0;
	break;
    case kEta:
	func=IpEta;
	break;
    default:
        printf("<AliGenPMDlib::GetIp> unknown parametrisation\n");
    }
    return func;
}




