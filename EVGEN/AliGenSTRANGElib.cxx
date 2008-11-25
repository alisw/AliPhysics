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


//======================================================================
//  AliGenSTRANGElib class contains parameterizations of the
//  kaon, phi and hyperon (Lambda, Anti-Lambda, Xi, anti-Xi, Omega,
//  anti-Omega)  for the PPR study of the strange particle production. 
//  These parameterizations are used by the 
//  AliGenParam  class:
//  AliGenParam(npar, param,  AliGenSTRANGElib::GetPt(param),
//                            AliGenSTRANGElib::GetY(param),
//                            AliGenSTRANGElib::GetIp(param) )
//  param represents the particle to be simulated. 
//  ?????????
//  Pt distributions are calculated from the transverse mass scaling 
//  with Pions, using the PtScal function taken from AliGenMUONlib 
//  version aliroot 3.01
//
//     Rocco CALIANDRO. Rosa Anna FINI, Tiziano VIRGILI
//     Rocco.Caliandro@cern.ch Rosanna.Fini@ba.infn.it, 
//     Tiziano.Virgili@roma1.infn.it
//======================================================================

/* $Id$ */

#include "TMath.h"
#include "TRandom.h"

#include "AliGenSTRANGElib.h"

ClassImp(AliGenSTRANGElib)

//============================================================= 
//
 Double_t AliGenSTRANGElib::PtScal(Double_t pt, Int_t np)
{
// Mt-scaling
// Function for the calculation of the Pt distribution for a 
// given particle np, from the pion Pt distribution using the 
// mt scaling. This function was taken from AliGenMUONlib 
// aliroot version 3.01, and was extended for hyperons.
// np = 1=>Pions 2=>Kaons 3=>Etas 4=>Omegas 5=>ETA' 6=>PHI
//      7=>BARYONS-BARYONBARS
//      8=>Lambda-antiLambda
//      9=>Xi-antiXi
//     10=>Omega-antiOmega

  //    MASS SCALING RESPECT TO PIONS
  //    MASS                1=>PI,  2=>K, 3=>ETA,4=>OMEGA,5=>ETA',6=>PHI 
  const Double_t khm[11] = {0.1396, 0.494,0.547, 0.782,   0.957,  1.02, 
  //    MASS               7=>BARYON-BARYONBAR  
                                 0.938, 
  //    MASS               8=>Lambda-antiLambda  
                                  1.1157,
  //    MASS               9=>Xi-antiXi  
                                  1.3213, 
  //    MASS              10=>Omega-antiOmega  
                                  1.6725,
  //    MASS              11=>Lambda(1520)
                                  1.5195};
  //     VALUE MESON/PI AT 5 GEV
  const Double_t kfmax[11]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  np--;
  Double_t f5=TMath::Power(((sqrt(100.018215)+2.)/(sqrt(100.+khm[np]*khm[np])+2.0)),12.3);
  Double_t kfmax2=f5/kfmax[np];
  // PIONS
  Double_t ptpion=100.*PtPion(&pt, (Double_t*) 0);
  Double_t fmtscal=TMath::Power(((sqrt(pt*pt+0.018215)+2.)/
                                 (sqrt(pt*pt+khm[np]*khm[np])+2.0)),12.3)/ kfmax2;
  return fmtscal*ptpion;
}
//============================================================= 
//
 Double_t AliGenSTRANGElib::PtPion(const Double_t *px, const Double_t *)
{
//     Pion transverse momentum distribtuion taken 
//     from AliGenMUONlib class, version 3.01 of aliroot
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
  Double_t y, y1, kxmpi2, ynorm, a;
  Double_t x=*px;
  //
  y1=TMath::Power(kp0/(kp0+kxlim),kxn);
  kxmpi2=kxmpi*kxmpi;
  ynorm=kb*(TMath::Exp(-sqrt(kxlim*kxlim+kxmpi2)/kt));
  a=ynorm/y1;
  if (x > kxlim)
    y=a*TMath::Power(kp0/(kp0+x),kxn);
  else
    y=kb*TMath::Exp(-sqrt(x*x+kxmpi2)/kt);
  return y*x;
}
// End Scaling
//============================================================================
//    K  A  O  N  
 Double_t AliGenSTRANGElib::PtKaon( const Double_t *px, const Double_t *)
{
//                kaon
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,2);  //  2==> Kaon in the PtScal function
}

 Double_t AliGenSTRANGElib::YKaon( const Double_t *py, const Double_t *)
{
// y-distribution
//____________________________________________________________

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

 Int_t AliGenSTRANGElib::IpKaon(TRandom *ran)
{
//                 particle composition
//

    Float_t random = ran->Rndm();
    Float_t random2 = ran->Rndm();
    if (random2 < 0.5) 
    {
      if (random < 0.5) {       
        return  321;   //   K+
      } else {
        return -321;   // K-
      }
    }
    else
    {  
      if (random < 0.5) {       
        return  130;   // K^0 short
      } else {  
        return  310;   // K^0 long
      }
    }
}
// End Kaons
//============================================================================
//============================================================================
//    P  H  I   
 Double_t AliGenSTRANGElib::PtPhi( const Double_t *px, const Double_t *)
{
// phi
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,6);  //  6==> Phi in the PtScal function
}

 Double_t AliGenSTRANGElib::YPhi( const Double_t *py, const Double_t *)
{
// y-distribution
//____________________________________________________________

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

 Int_t AliGenSTRANGElib::IpPhi(TRandom *)
{
//                 particle composition
//
    
        return  333;   //   Phi      
}
// End Phis
//===================================================================
//============================================================================
//    Lambda
 Double_t AliGenSTRANGElib::PtLambda( const Double_t *px, const Double_t *)
{
// Lambda
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,8);  //  8==> Lambda-antiLambda in the PtScal function
}

 Double_t AliGenSTRANGElib::YLambda( const Double_t *py, const Double_t *)
{
// y-distribution
//____________________________________________________________

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

 Int_t AliGenSTRANGElib::IpLambda(TRandom *ran)
{
//                 particle composition
//                 generation of fixed type of particle
//
    Float_t random = ran->Rndm();
    if (random < 0.5) {       
      return  3122;   //   Lambda 
    } else {  
      return -3122;   //   Anti-Lambda
    }
}
// End Lambda
//============================================================================
//    XIminus
 Double_t AliGenSTRANGElib::PtXiMinus( const Double_t *px, const Double_t *)
{
// Xi
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,9);  //  9==> Xi-antiXi in the PtScal function
}

 Double_t AliGenSTRANGElib::YXiMinus( const Double_t *py, const Double_t *)
{
// y-distribution
//____________________________________________________________

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

 Int_t AliGenSTRANGElib::IpXiMinus(TRandom *ran)
{
//                 particle composition
//                 generation of fixed type of particle
//
    Float_t random = ran->Rndm();
    if (random < 0.5) {       
      return  3312;   //   Xi- 
    } else {  
      return -3312;   //   Xi+
    }
}
// End Ximinus
//============================================================================
//    Omegaminus
 Double_t AliGenSTRANGElib::PtOmegaMinus( const Double_t *px, const Double_t *)
{
// Omega
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,10);  //  10==> Omega-antiOmega in the PtScal function
}

 Double_t AliGenSTRANGElib::YOmegaMinus( const Double_t *py, const Double_t *)
{
// y-distribution
//____________________________________________________________

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

 Int_t AliGenSTRANGElib::IpOmegaMinus(TRandom * ran)
{
//                 particle composition
//                 generation of fixed type of particle
//

    Float_t random = ran->Rndm();
    if (random < 0.5) {       
      return  3334;   //   Omega- 
    } else {  
      return -3334;   //   Omega+
    }
}
// End Omegaminus
//============================================================================
//     Lambda(1520)
Double_t AliGenSTRANGElib::PtLambda1520( const Double_t *px, const Double_t *)
{
// Lambda(1520)
//                  pt-distribution
//____________________________________________________________

  return PtScal(*px,11);   //   11=> Lambda(1520) in the PtScal function
}

Double_t AliGenSTRANGElib::YLambda1520( const Double_t *py, const Double_t *)
{
// y-distribution
//____________________________________________________________

  const Double_t ka   = 1000.;
  const Double_t kdy  = 4.;

  
  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

Int_t AliGenSTRANGElib::IpLambda1520(TRandom * ran)
{
//                 particle composition
//                 generation of fixed type of particle
//

   Float_t random = ran->Rndm();
   if (random < 0.5) {       
     return  3124;   //   Lambda(1520) 
   } else {  
     return -3124;   //   antiLambda(1520)
   }
}
// End Lambda(1520)
//============================================================================

typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);
 GenFunc AliGenSTRANGElib::GetPt(Int_t param, const char* /*tname*/) const
{
// Return pinter to pT parameterisation
    GenFunc func;
    
    switch (param)
    {
    case kKaon:
        func=PtKaon;
        break;
    case kPhi:
        func=PtPhi;
        break;
    case kLambda:
        func=PtLambda;
        break;
    case kXiMinus:
        func=PtXiMinus;
        break;
    case kOmegaMinus:
        func=PtOmegaMinus;
        break;
    case kLambda1520:
        func=PtLambda1520;
        break;
    default:
        func=0;
        printf("<AliGenSTRANGElib::GetPt> unknown parametrisationn");
    }
    return func;
}

 GenFunc AliGenSTRANGElib::GetY(Int_t param, const char* /*tname*/) const
{
// Return pointer to Y parameterisation
    GenFunc func;
    switch (param)
    {
    case kKaon:
        func=YKaon;
        break;
    case kPhi:
        func=YPhi;
        break;
    case kLambda:
        func=YLambda;
        break;
    case kXiMinus:
        func=YXiMinus;
        break;
    case kOmegaMinus:
        func=YOmegaMinus;
        break;
    case kLambda1520:
        func=YLambda1520;
        break;
    default:
        func=0;
        printf("<AliGenSTRANGElib::GetY> unknown parametrisationn");
    }
    return func;
}
typedef Int_t (*GenFuncIp) (TRandom *);
 GenFuncIp AliGenSTRANGElib::GetIp(Int_t param,  const char* /*tname*/) const
{
// Return pointer to particle composition
    GenFuncIp func;
    switch (param)
    {
    case kKaon:
        func=IpKaon;
        break;
    case kPhi:
        func=IpPhi;
        break;
    case kLambda:
        func=IpLambda;
        break;
    case kXiMinus:
        func=IpXiMinus;
        break;
    case kOmegaMinus:
        func=IpOmegaMinus;
        break;
    case kLambda1520:
        func=IpLambda1520;
        break;
    default:
        func=0;
        printf("<AliGenSTRANGElib::GetIp> unknown parametrisationn");
    }
    return func;
}

