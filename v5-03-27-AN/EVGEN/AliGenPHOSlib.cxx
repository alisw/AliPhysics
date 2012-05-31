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

//======================================================================
//  AliGenPHOSlib class contains parameterizations of the
//  pion, kaon, eta, omega, etaprime, phi and baryon (proton, 
//  antiproton, neutron and anti-neutron) particles for the 
//  study of the neutral background in PHOS detector. 
//  These parameterizations are used by the 
//  AliGenParam  class:
//  AliGenParam(npar, param,  AliGenPHOSlib::GetPt(param),
//                            AliGenPHOSlib::GetY(param),
//                            AliGenPHOSlib::GetIp(param) )
//  param represents the particle to be simulated : 
//  Pion, Kaon, Eta, Omega, Etaprime, Phi or Baryon    
//  Pt distributions are calculated from the transverse mass scaling 
//  with Pions, using the PtScal function taken from AliGenMUONlib 
//  version aliroot 3.01
//
//     Gines MARTINEZ. Laurent APHECETCHE and Yves SCHUTZ
//      GPS @ SUBATECH,  Nantes , France  (October 1999)
//     http://www-subatech.in2p3.fr/~photons/subatech
//     martinez@subatech.in2p3.fr
//  Additional particle species simulation options has been added: 
//  Charged Pion, Charged Kaons, KLong Proton, Anti-Proton, Neutron, 
//  Anti-Neutron --> Changes made by Gustavo Conesa in November 2004
//  Add flat Omega(782) distribution in Nov. 2010 by Renzhuo WAN
//======================================================================

#include "TMath.h"
#include "TRandom.h"

#include "AliGenPHOSlib.h"

ClassImp(AliGenPHOSlib)

//======================================================================
//    P  I  O  N  S
//    (From GetPt, GetY and GetIp as param = Pion)
//    Transverse momentum distribution" PtPion 
//    Rapidity distribution YPion
//    Particle distribution IdPion  111, 211 and -211 (pi0, pi+ and pi-)
//
 Double_t AliGenPHOSlib::PtPion(const Double_t *px, const Double_t *)
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
 Double_t AliGenPHOSlib::YPion( const Double_t *py, const Double_t *)
{
//
// pion y-distribution
//

  const Double_t ka    = 7000.;   
  const Double_t kdy   = 4.;

  Double_t y=TMath::Abs(*py);
  //
  Double_t ex = y*y/(2*kdy*kdy);
  return ka*TMath::Exp(-ex);
}

Int_t AliGenPHOSlib::IpPion(TRandom */*ran*/)
{
//                 particle composition pi+, pi0, pi-
//

        return 111  ;
}
 Int_t AliGenPHOSlib::IpChargedPion(TRandom *ran)
{
//                 particle composition pi+, pi0, pi-
//

     Float_t random = ran->Rndm();

     if ( (2.*random)  < 1. ) 
       {
	 return 211 ;
       } 
     else
       {  
	 return -211 ;
       }
}

//End Pions
//======================================================================
//    Pi 0 Flat Distribution
//    Transverse momentum distribution PtPi0Flat
//    Rapidity distribution YPi0Flat
//    Particle distribution IdPi0Flat  111 (pi0)
//

Double_t AliGenPHOSlib::PtPi0(const Double_t * px, const Double_t *)
{
//     Pion transverse momentum 
    const Double_t kp0 =1.35; 
    const Double_t kxn= 6.18;
    return TMath::Power(kp0 /(kp0 + px[0]), kxn);
}

Double_t AliGenPHOSlib::PtPi0Flat(const Double_t */*px*/, const Double_t *)
{
//     Pion transverse momentum flat distribution 

return 1;

}

Double_t AliGenPHOSlib::YPi0Flat( const Double_t */*py*/, const Double_t *)
{

// pion y-distribution
//
  return 1.;
}

 Int_t AliGenPHOSlib::IpPi0Flat(TRandom *)
{

//                 particle composition pi0
//
        return 111 ;
}
// End Pi0Flat
//============================================================= 
//
 Double_t AliGenPHOSlib::PtScal(Double_t pt, Int_t np)
{
// Mt-scaling
// Fonction for the calculation of the Pt distribution for a 
// given particle np, from the pion Pt distribution using the 
// mt scaling. This function was taken from AliGenMUONlib 
// aliroot version 3.01, and was extended for baryons
// np = 1=>Pions 2=>Kaons 3=>Etas 4=>Omegas 5=>ETA' 6=>PHI
//      7=>BARYONS-BARYONBARS

  //    SCALING EN MASSE PAR RAPPORT A PTPI
  //    MASS                0=>PI,  1=>K, 2=>ETA, 3=>OMEGA,  4=>ETA',5=>PHI
  const Double_t khm[10] = {0.1396, 0.494,  0.547,    0.782,   0.957,   1.02, 
  //    MASS               6=>BARYON-BARYONBAR  
                                         0.938, 0. , 0., 0.};
  //     VALUE MESON/PI AT 5 GEV
  const Double_t kfmax[10]={1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  Double_t f5=TMath::Power(((sqrt(100.018215)+2.)/(sqrt(100.+khm[np]*khm[np])+2.0)),12.3);
  Double_t kfmax2=f5/kfmax[np];
  // PIONS
  Double_t ptpion=100.*PtPion(&pt, (Double_t*) 0);
  Double_t fmtscal=TMath::Power(((sqrt(pt*pt+0.018215)+2.)/
                                 (sqrt(pt*pt+khm[np]*khm[np])+2.0)),12.3)/ kfmax2;
  return fmtscal*ptpion;

}
// End Scaling
//============================================================================
//    K  A  O  N  S
 Double_t AliGenPHOSlib::PtKaon( const Double_t *px, const Double_t *)
{
//                kaon
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,1);  //  1==> Kaon in the PtScal function
}

 Double_t AliGenPHOSlib::YKaon( const Double_t *py, const Double_t *)
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

 Int_t AliGenPHOSlib::IpKaon(TRandom *ran)
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

 Int_t AliGenPHOSlib::IpChargedKaon(TRandom *ran)
{
//                 particle composition
//
  
  Float_t random = ran->Rndm();
  
  if (random < 0.5) {       
    return  321;   //   K+
  } else {
    return -321;   // K-
  }
  
  
}
Int_t AliGenPHOSlib::IpKaon0L(TRandom *)
{
  //                 particle composition
  //
  
	return  130;   // K^0 long
}
// End Kaons
//============================================================================
//============================================================================
//   E  T  A  S
 Double_t AliGenPHOSlib::PtEta( const Double_t *px, const Double_t *)
{
//                etas
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,2);  //  2==> Eta in the PtScal function
}

 Double_t AliGenPHOSlib::YEta( const Double_t *py, const Double_t *)
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

 Int_t AliGenPHOSlib::IpEta(TRandom *)
{
//                 particle composition
//

        return  221;   //   eta
}
// End Etas

//======================================================================
//    Eta Flat Distribution
//    Transverse momentum distribution PtEtaFlat
//    Rapidity distribution YEtaFlat
//    Particle distribution IdEtaFlat  111 (pi0)
//

Double_t AliGenPHOSlib::PtEtaFlat(const Double_t */*px*/, const Double_t *)
{
//     Eta transverse momentum flat distribution 

  return 1;

}

Double_t AliGenPHOSlib::YEtaFlat( const Double_t */*py*/, const Double_t *)
{
//
// pion y-distribution
//
  return 1.;
}

 Int_t AliGenPHOSlib::IpEtaFlat(TRandom *)
{
//
//                 particle composition eta
//
        return 221 ;
}
// End EtaFlat
//============================================================================
//============================================================================
//    O  M  E  G  A  S
 Double_t AliGenPHOSlib::PtOmega( const Double_t *px, const Double_t *)
{
// omegas
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,3);  //  3==> Omega in the PtScal function
}

 Double_t AliGenPHOSlib::YOmega( const Double_t *py, const Double_t *)
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

 Int_t AliGenPHOSlib::IpOmega(TRandom *)
{
//                 particle composition
//

        return  223;   // Omega
}
// End Omega
//============================================================================
//======================================================================
//    Omega(782) Flat Distribution
//    Transverse momentum distribution PtOmegaFlat
//    Rapidity distribution YOmegaFlat
//    Particle distribution IdOmegaFlat  223(0mega)
//

Double_t AliGenPHOSlib::PtOmegaFlat(const Double_t */*px*/, const Double_t *)
{
//     omega transverse momentum flat distribution 

return 1;

}

Double_t AliGenPHOSlib::YOmegaFlat( const Double_t */*py*/, const Double_t *)
{

// omega y-distribution
//
  return 1.;
}

 Int_t AliGenPHOSlib::IpOmegaFlat(TRandom *)
{

//                 particle composition omega
//
        return 223 ;
}
// End OmegaFlat


//============================================================================
//    E  T  A  P  R  I  M  E
 Double_t AliGenPHOSlib::PtEtaprime( const Double_t *px, const Double_t *)
{
// etaprime
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,4);  //  4==> Etaprime in the PtScal function
}

 Double_t AliGenPHOSlib::YEtaprime( const Double_t *py, const Double_t *)
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

 Int_t AliGenPHOSlib::IpEtaprime(TRandom *)
{
//                 particle composition
//

        return  331;   //   Etaprime
}
// End EtaPrime
//===================================================================
//============================================================================
//    P  H  I   S
 Double_t AliGenPHOSlib::PtPhi( const Double_t *px, const Double_t *)
{
// phi
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,5);  //  5==> Phi in the PtScal function
}

 Double_t AliGenPHOSlib::YPhi( const Double_t *py, const Double_t *)
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

 Int_t AliGenPHOSlib::IpPhi(TRandom *)
{
//                 particle composition
//
    
        return  333;   //   Phi      
}
// End Phis
//===================================================================
//============================================================================
//    B  A  R  Y  O  N  S  == protons, protonsbar, neutrons, and neutronsbars
 Double_t AliGenPHOSlib::PtBaryon( const Double_t *px, const Double_t *)
{
// baryons
//                pt-distribution
//____________________________________________________________

  return PtScal(*px,6);  //  6==> Baryon in the PtScal function
}

 Double_t AliGenPHOSlib::YBaryon( const Double_t *py, const Double_t *)
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

 Int_t AliGenPHOSlib::IpBaryon(TRandom *ran)
{
//                 particle composition
//

    Float_t random = ran->Rndm();
    Float_t random2 = ran->Rndm();
    if (random2 < 0.5) 
    {
      if (random < 0.5) {       
        return  2212;   //   p
      } else {
        return -2212;   // pbar
      }
    }
    else
    {  
      if (random < 0.5) {       
        return  2112;   // n
      } else {  
        return -2112;   // n bar
      }
    }
}

 Int_t AliGenPHOSlib::IpProton(TRandom *)
{
//                 particle composition
//  
        return  2212;   //   p

}
 Int_t AliGenPHOSlib::IpAProton(TRandom *)
{
//                 particle composition
//  
        return  -2212;   //   p bar

}

 Int_t AliGenPHOSlib::IpNeutron(TRandom *)
{
//                 particle composition
//  
        return  2112;   //   n

}
 Int_t AliGenPHOSlib::IpANeutron(TRandom *)
{
//                 particle composition
//  
        return  -2112;   //   n

}
// End Baryons
//===================================================================

typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);
GenFunc AliGenPHOSlib::GetPt(Int_t param, const char* /*tname*/) const
{
// Return pinter to pT parameterisation
    GenFunc func;
    
    switch (param)
      {
      case kPion:     
        func=PtPion;
        break;
      case kPi0:     
        func=PtPi0;
        break;
      case kPi0Flat:     
        func=PtPi0Flat;
        break;
      case kKaon:
        func=PtKaon;
        break;
      case kEta:
        func=PtEta;
        break;
      case kEtaFlat:
        func=PtEtaFlat;
        break;
      case kOmega:
        func=PtOmega;
        break;
      case kOmegaFlat:
        func=PtOmegaFlat;
      break;
      case kEtaPrime:
        func=PtEtaprime;
        break;
      case kBaryon:
        func=PtBaryon;
        break;
      default:
        func=0;
        printf("<AliGenPHOSlib::GetPt> unknown parametrisationn");
      }
    return func;
}

GenFunc AliGenPHOSlib::GetY(Int_t param, const char* /*tname*/) const
{
  // Return pointer to Y parameterisation
  GenFunc func;
  switch (param)
    {
    case kPion:
      func=YPion;
      break;
    case kPi0:     
    case kPi0Flat:
      func=YPi0Flat;
      break;
    case kKaon:
      func=YKaon;
      break;
    case kEta:
      func=YEta;
      break;
    case kEtaFlat:
      func=YEtaFlat;
      break;
    case kOmega:
      func=YOmega;
      break;
    case kOmegaFlat:
      func=YOmegaFlat;
      break;
    case kEtaPrime:
      func=YEtaprime;
      break;
    case kPhi:
      func=YPhi;
      break;
    case kBaryon:
      func=YBaryon;
      break;
    default:
      func=0;
      printf("<AliGenPHOSlib::GetY> unknown parametrisationn");
    }
  return func;
}
typedef Int_t (*GenFuncIp) (TRandom *);
GenFuncIp AliGenPHOSlib::GetIp(Int_t param,  const char* /*tname*/) const
{
  // Return pointer to particle composition
  GenFuncIp func;
  switch (param)
    {
    case kPion:       
      func=IpPion;
      break;
    case kChargedPion:       
      func=IpChargedPion;
      break;
    case kPi0:     
    case kPi0Flat:       
      func=IpPi0Flat;
      break;
    case kKaon:
      func=IpKaon;
      break;
    case kChargedKaon:
      func=IpChargedKaon;
      break;
    case kKaon0L:
      func=IpKaon0L;
      break;
    case kEta:
      func=IpEta;
      break;
    case kEtaFlat:
      func=IpEtaFlat;
      break;
      
    case kOmega:
      func=IpOmega;
      break;
    case kOmegaFlat:
      func=IpOmegaFlat;
      break;
    case kEtaPrime:
      func=IpEtaprime;
      break;
    case kPhi:
      func=IpPhi;
      break;
    case kBaryon:
      func=IpBaryon;
      break;
    case kProton:
      func=IpProton;
      break;
    case kAProton:
      func=IpAProton;
      break;
    case kNeutron:
      func=IpNeutron;
      break;
    case kANeutron:
      func=IpANeutron;
      break;
      
    default:
      func=0;
      printf("<AliGenPHOSlib::GetIp> unknown parametrisationn");
    }
  return func;
}

