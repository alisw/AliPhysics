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

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Implementation of AliGenlib to collect parametrisations used for        //
// GSI simulations.                                                        //
// It is an extension of AliMUONLib providing in addition the option       //
// for different parametrisations of pt, y and ip for every particle type  //
//                                                                         //
// Responsible: Andres.Sandoval@cern.ch                                    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "AliGenGSIlib.h"


ClassImp(AliGenGSIlib)

//==========================================================================
//
//              Definition of Particle Distributions
//                    
//==========================================================================
//
//                         Upsilon 
//
//--------------------------------------------------------------------------
//
//                  upsilon particle composition
//
//--------------------------------------------------------------------------
Int_t AliGenGSIlib::IpUpsilon(TRandom *)
{
// Return upsilon pdg code

  return 553;     

}
Double_t AliGenGSIlib::PtUpsilonFlat( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                   upsilon pt-distribution FLAT
//
//____________________________________________________________--------------
    
  const Double_t kptmin = 0.0;
  const Double_t kptmax  = 15.0;
  Double_t x=*px;
  Double_t weight = 0.;

  if ((x > kptmin) &&  (x < kptmax)) weight = 1.;

  return weight;
   
}
Double_t AliGenGSIlib::YUpsilonFlat(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    upsilon y-distribution FLAT
//
//--------------------------------------------------------------------------

  const Double_t ky0 = 1.5;
  const Double_t kb=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);

  if (y < ky0)
    yu=kb;
  else
    yu = 0.;

  return yu;

}
Double_t AliGenGSIlib::PtUpsilonRitman( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                 upsilon pt-distribution RITMAN
//
//--------------------------------------------------------------------------

  const Double_t kpt0 = 4.7;
  const Double_t kxn  = 3.5;
  Double_t x=*px;

  Double_t pass1 = 1.+((x*x)/(kpt0*kpt0));

  return x/TMath::Power(pass1,kxn);
   
}
Double_t AliGenGSIlib::YUpsilonRitman(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                  upsilon y-distribution RITMAN
//
//--------------------------------------------------------------------------

  const Double_t ky0 = 3.;
  const Double_t kb=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);

  if (y < ky0)
    yu=kb;
  else
    yu=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);

  return yu;
   
}
Double_t AliGenGSIlib::PtUpsilonKarel( const Double_t */*px*/, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                 upsilon pt-distribution kAREL
//
//--------------------------------------------------------------------------
// to implement

  return 0.1;   

}
Double_t AliGenGSIlib::YUpsilonKarel(const Double_t */*py*/, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                  upsilon y-distribution KAREL
//
//--------------------------------------------------------------------------
  
  //to implement

  return 0.2;  

}
Double_t AliGenGSIlib::PtUpsilonMUON( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                 upsilon pt-distribution MUONlib
//
//--------------------------------------------------------------------------

  const Double_t kpt0 = 5.3;
  const Double_t kxn  = 2.5;
  Double_t x=*px;

  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);

  return x/TMath::Power(pass1,kxn);

}
Double_t AliGenGSIlib::YUpsilonMUON(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                   upsilon y-distribution MUONlib
//
//--------------------------------------------------------------------------

  const Double_t ky0 = 3.;
  const Double_t kb=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);

  if (y < ky0)
    yu=kb;
  else
    yu=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);

  return yu;

}
//--------------------------------------------------------------------------
//
//                             J/Psi
//
Int_t AliGenGSIlib::IpJpsi(TRandom *)
{
//--------------------------------------------------------------------------
//
//                    J/Psi particle composition
//
//--------------------------------------------------------------------------

  return 443;     

}
Double_t AliGenGSIlib::PtJpsiFlat( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                   J/Psi pt-distribution FLAT
//
//--------------------------------------------------------------------------

  const Double_t kptmin = 0.0;
  const Double_t kptmax  = 15.0;
  Double_t x=*px;
  Double_t weight = 0.;

  if ((x > kptmin) && (x < kptmax)) weight = 1.;

  return weight;
   
}
Double_t AliGenGSIlib::YJpsiFlat(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    J/Psi y-distribution FLAT
//
//--------------------------------------------------------------------------

  const Double_t ky0 = 1.5;
  const Double_t kb=1.;
  Double_t yu;
  Double_t y=TMath::Abs(*py);

  if (y < ky0)
    yu=kb;
  else
    yu = 0.;

  return yu;

}
Double_t AliGenGSIlib::PtJpsiMUON( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                   J/Psi pt-distribution MUONlib
//
//--------------------------------------------------------------------------

  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;

  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
   
}
Double_t AliGenGSIlib::PtJpsiRitman( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                   J/Psi pt-distribution Ritman
//
//--------------------------------------------------------------------------

  const Double_t kpt0 = 2.3;
  const Double_t kxn  = 3.5;
  Double_t x=*px;

  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);

  return x/TMath::Power(pass1,kxn);
   
}
Double_t AliGenGSIlib::YJpsiMUON(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    J/Psi y-distribution MUONlib
//
//--------------------------------------------------------------------------

  const Double_t ky0 = 4.;
  const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);

  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;

}
//--------------------------------------------------------------------------
//
//                  J/Psi pt-distribution by Sergei
//
//--------------------------------------------------------------------------
//Double_t AliGenGSIlib::PtJpsi( Double_t *px, Double_t */*dummy*/ )
//{

//  return x = gRandom->Rndm()*10.;

//}
//--------------------------------------------------------------------------
//
//                  J/Psi y-distribution by Sergei
//
//--------------------------------------------------------------------------
/*Double_t AliGenGSIlib::YJpsi(Double_t *py, Double_t *dummy)
{

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
*/
//--------------------------------------------------------------------------
//
//                        Charm
//
//--------------------------------------------------------------------------
Int_t AliGenGSIlib::IpCharm(TRandom *ran)
{
//
//                    charm particle composition
//
//--------------------------------------------------------------------------
  
    Float_t random;
    Int_t ip;
    // 411,421,431,4122
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
Double_t AliGenGSIlib::PtCharmFlat( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    charm pt-distribution, FLAT
//
//--------------------------------------------------------------------------

  Double_t x=*px;

  if (x>10.)  x = 0.;
  else x=1.;
  return x ;

}
Double_t AliGenGSIlib::PtCharmGSI( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    charm pt-distribution, from Dariuzs Miskowiec
//
//--------------------------------------------------------------------------

  //Taken from PYTHIA with MRS D-' (3031 from PDFLIB), K=3.0
  const Double_t kp1 = 1.3;
  const Double_t kp2  = 0.39;
  const Double_t kp3  = 0.018;
  const Double_t kp4  = 0.91;
  Double_t x=*px;

  Double_t pass1 =TMath::Exp(-x/kp2)  ;
  Double_t pass2 =TMath::Exp(-x/kp4)  ;
  return TMath::Power(x,kp1) * (pass1 + kp3 * pass2);

}
Double_t AliGenGSIlib::PtCharmMUON( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    charm pt-distribution, from MUONlib
//
//--------------------------------------------------------------------------

  const Double_t kpt0 = 4.08;
  const Double_t kxn  = 9.40;
  Double_t x=*px;

  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);

  return x/TMath::Power(pass1,kxn);

}
Double_t AliGenGSIlib::YCharm( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    charm y-distribution
//
//--------------------------------------------------------------------------

    Double_t *dum=0;

    return YJpsiMUON(px,dum);

}
//--------------------------------------------------------------------------
//
//                        Beauty
//
//--------------------------------------------------------------------------
Int_t AliGenGSIlib::IpBeauty(TRandom *ran)
{  
//
//                    beauty particle composition
//
//--------------------------------------------------------------------------

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
Double_t AliGenGSIlib::PtBeautyFlat( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    beauty pt-distribution, FLAT
//
//--------------------------------------------------------------------------

  Double_t x=*px;

  if (x>10.) x=0.;
  else x = 1.;
  return x ;

}
Double_t AliGenGSIlib::PtBeautyGSI( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//
//                    beauty pt-distribution, from D. Miskowiec
//
//--------------------------------------------------------------------------

  //Taken from PYTHIA with MRS D-' (3031 from PDFLIB), K=3.0
  const Double_t kp1 = 1.3;
  const Double_t kp2  = 1.78;
  const Double_t kp3  = 0.0096;
  const Double_t kp4  = 4.16;
  Double_t x=*px;

  Double_t pass1 =TMath::Exp(-x/kp2)  ;
  Double_t pass2 =TMath::Exp(-x/kp4)  ;

  return TMath::Power(x,kp1) * (pass1 + kp3 * pass2);

}
Double_t AliGenGSIlib::PtBeautyMUON( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    beauty pt-distribution, from MUONlib
//
//--------------------------------------------------------------------------

  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;

  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);

  return x/TMath::Power(pass1,kxn);

}
Double_t AliGenGSIlib::YBeauty( const Double_t *px, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    beauty y-distribution
//
//--------------------------------------------------------------------------

    Double_t *dum=0;

    return YJpsiMUON(px,dum);

}
//--------------------------------------------------------------------------
//
//                          Eta
//
//--------------------------------------------------------------------------
Int_t AliGenGSIlib::IpEta(TRandom *)
{
//
//                 eta particle composition 
//
//--------------------------------------------------------------------------

  return 221;     

}
Double_t AliGenGSIlib::PtEtaPHOS( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                  eta pt-distribution
//
//____________________________________________________________--------------

  return PtScal(*px,3);  //  3==> Eta in the PtScal function
   
}
Double_t AliGenGSIlib::YEtaPHOS(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                   eta y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;

  Double_t y=TMath::Abs(*py);

  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
//--------------------------------------------------------------------------
//
//                       Etaprime
//
//--------------------------------------------------------------------------
Int_t AliGenGSIlib::IpEtaprime(TRandom *)
{
//
//                 etaprime particle composition 
//
//--------------------------------------------------------------------------

  return 331;     

}
Double_t AliGenGSIlib::PtEtaprimePHOS( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                 etaprime pt-distribution
//
//____________________________________________________________--------------

  return PtScal(*px,5);  //  5==> Etaprime in the PtScal function
   
}
Double_t AliGenGSIlib::YEtaprimePHOS(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                  etaprime y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;

  Double_t y=TMath::Abs(*py);

  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
//--------------------------------------------------------------------------
//
//                       omega 
//
//--------------------------------------------------------------------------
Int_t AliGenGSIlib::IpOmega(TRandom *)
{
//
//                 omega particle composition 
//
//--------------------------------------------------------------------------

  return 223;     

}
Double_t AliGenGSIlib::PtOmega( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                  omega pt-distribution
//
//____________________________________________________________--------------

  return PtScal(*px,4);  //  4==> Omega in the PtScal function
   
}
Double_t AliGenGSIlib::YOmega(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                   omega y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);

  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
//--------------------------------------------------------------------------
//
//                       Rho 
//
//--------------------------------------------------------------------------

Int_t AliGenGSIlib::IpRho(TRandom *)
{
//
//                 rho particle composition 
//
//--------------------------------------------------------------------------

  return 113;     

}
Double_t AliGenGSIlib::PtRho( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                  rho pt-distribution
//
//____________________________________________________________--------------

  return PtScal(*px,11);  //  11==> Rho in the PtScal function
   
}
Double_t AliGenGSIlib::YRho(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                   rho y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);

  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
//--------------------------------------------------------------------------
//
//                              Pion
//
//--------------------------------------------------------------------------
Int_t AliGenGSIlib::IpPionPHOS(TRandom *ran)
{
//
//                 particle composition  pi+, pi0, pi-
//
//--------------------------------------------------------------------------

    Float_t random = ran->Rndm();

    if ( (3.*random)  < 1. ) 
    {
          return 211 ;
    } 
    else
    {  
      if ( (3.*random) >= 2.)
      {
         return -211 ;
      }
      else 
      {
        return 111  ;
      }
    }
}
Double_t AliGenGSIlib::PtPion( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                  pion pt-distribution
//
//       Pion transverse momentum distribtuion as in AliGenMUONlib class, 
//       version 3.01 of aliroot:
//       PT-PARAMETERIZATION CDF, PRL 61(88) 1819
//       POWER LAW FOR PT > 500 MEV
//       MT SCALING BELOW (T=160 MEV)
//
//____________________________________________________________--------------

  const Double_t kp0 = 1.3;
  const Double_t kxn = 8.28;
  const Double_t kxlim=0.5;
  const Double_t kt=0.160;
  const Double_t kxmpi=0.139;
  const Double_t kb=1.;
  Double_t y, y1, kxmpi2, ynorm, a;
  Double_t x=*px;

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
Double_t AliGenGSIlib::YPion(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    pion y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 7000.;   
  const Double_t kdy   = 4.;

  Double_t y=TMath::Abs(*py);

  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
Int_t AliGenGSIlib::IpKaonPHOS(TRandom *ran)
{
//--------------------------------------------------------------------------
//
//
//                        Kaon
//--------------------------------------------------------------------------
//
//                kaon particle composition K+, K-, Ko_short, Ko_long
//
//--------------------------------------------------------------------------

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
Double_t AliGenGSIlib::PtKaonPHOS( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                   kaon pt-distribution
//
//____________________________________________________________--------------

  return PtScal(*px,2);  //  2==> Kaon in the PtScal function
   
}
Double_t AliGenGSIlib::YKaonPHOS(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    kaon y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;

  Double_t y=TMath::Abs(*py);

  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
//--------------------------------------------------------------------------
//
//                        Phi  
//
Int_t AliGenGSIlib::IpPhi(TRandom *)
{
//--------------------------------------------------------------------------
//
//                 particle composition 
//
//--------------------------------------------------------------------------

  return 333;     

}
Double_t AliGenGSIlib::PtPhiPHOS( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                   phi pt-distribution
//
//____________________________________________________________--------------

  return PtScal(*px,6);  //  6==> Phi in the PtScal function
   
}
Double_t AliGenGSIlib::YPhiPHOS(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                    phi y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;


  Double_t y=TMath::Abs(*py);
 
  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
Int_t AliGenGSIlib::IpBaryons(TRandom *ran)
{
//--------------------------------------------------------------------------
//
//                          Baryons
//
//--------------------------------------------------------------------------
//
//                 baryons particle composition p, pbar, n, nbar
//
//--------------------------------------------------------------------------

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
Double_t AliGenGSIlib::PtBaryons( const Double_t *px, const Double_t */*dummy*/ )
{
//--------------------------------------------------------------------------
//
//                  baryons pt-distribution
//
//____________________________________________________________--------------

  return PtScal(*px,7);  //  7==> Baryon in the PtScal function
   
}
Double_t AliGenGSIlib::YBaryons(const Double_t *py, const Double_t */*dummy*/)
{
//--------------------------------------------------------------------------
//
//                   baryons y-distribution 
//
//--------------------------------------------------------------------------

  const Double_t ka    = 1000.;
  const Double_t kdy   = 4.;

  Double_t y=TMath::Abs(*py);

  Double_t ex = y*y/(2*kdy*kdy);

  return ka*TMath::Exp(-ex);

}
//============================================================= 
//
//                    Mt-scaling as in AliGenPHOSlib  
//
//============================================================= 
//
 Double_t AliGenGSIlib::PtScal(Double_t pt, Int_t np)
{
// Function for the calculation of the Pt distribution for a 
// given particle np, from the pion Pt distribution using the 
// mt scaling. 

// It was taken from AliGenPHOSlib aliroot version 3.04, which 
// is an update of the one in AliGenMUONlib aliroot version 3.01
// with an extension for Baryons but supressing Rhos
// np = 1=>Pions 2=>Kaons 3=>Etas 4=>Omegas 5=>ETA' 6=>PHI
//      7=>BARYONS-BARYONBARS

// The present adds the Rhos

// MASS   1=>PI, 2=>K, 3=>ETA, 4=>OMEGA, 5=>ETA', 6=>PHI 
//        7=>BARYON-BARYONBAR, 11==>RHO

  const Double_t khm[11] = {0.1396, 0.494,  0.547,    0.782,   0.957,   1.02, 
                                         0.938, 0. , 0., 0., 0.769};

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

//==========================================================================
//
//                     Set Getters 
//    
//==========================================================================

typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);

typedef Int_t (*GenFuncIp) (TRandom *);

GenFunc AliGenGSIlib::GetPt(Int_t param, const char * tname) const
{
// Return pointer to pT parameterisation
   GenFunc func=0;
   TString sname(tname);

   switch (param) 
    {
    case kUpsilon:
      if (sname=="FLAT"){
        func= PtUpsilonFlat;
        break;
      }
      if (sname=="MUON"){
        func= PtUpsilonMUON;
        break;
      }
      if (sname=="RITMAN"){
        func=PtUpsilonRitman;
        break;
      }
      if (sname=="KAREL"){
        func=PtUpsilonKarel;
        break;
      }
      func=0;
      printf("<AliGenGSIlib::GetPt> unknown parametrisation\n");
      break;

    case kJPsi:
      if (sname=="FLAT"){
        func= PtJpsiFlat;
        break;
      }
      if (sname=="MUON"){
        func= PtJpsiMUON;
        break;
      }
      //      if (sname=="SERGEI"){
      //	func= PtJpsi;
      //	break;
      //     }
      func=0;
      printf("<AliGenGSIlib::GetPt> unknown parametrisation\n");
      break;

    case kCharm: 
      if (sname=="FLAT"){
        func= PtCharmFlat;
        break;
      }

      if (sname=="MUON"){
        func= PtCharmMUON;
        break;
      }

      if (sname=="GSI"){
        func= PtCharmGSI;
        break;
      }
      func=0;
      printf("<AliGenGSIlib::GetPt> unknown parametrisation\n");
      break;

    case kBeauty: 
      if (sname=="FLAT"){
        func= PtBeautyFlat;
        break;
      }
      if (sname=="MUON"){
        func= PtBeautyMUON;
        break;
      }
      if (sname=="GSI"){
        func= PtBeautyGSI;
        break;
      }
      func=0;
      printf("<AliGenGSIlib::GetPt> unknown parametrisation\n");
      break;


    case kEta:
      func=PtEtaPHOS;
      break;

    case kEtaprime:
      func=PtEtaprimePHOS;
      break;

    case kOmega:
      func=PtOmega;
      break;

    case kRho:
      func=PtRho;
      break;

    case kKaon:
      func=PtKaonPHOS;
      break;

    case kPion:
      func=PtPion;
      break;

    case kPhi:
      func=PtPhiPHOS;
      break;

	 //    case kLambda:
	 //         func=PtLambda;
	 //         break;

    case kBaryons:
      func=PtBaryons;
      break;

    default:
      func=0;
      printf("<AliGenGSIlib::GetPt> unknown parametrisation\n");
    }
   return func;
}

GenFunc AliGenGSIlib::GetY(Int_t param, const char * tname) const
{
// Return pointer to y- parameterisation
   GenFunc func=0;
   TString sname(tname);

   switch (param) 
    {
    case kUpsilon:
      if (sname=="FLAT"){
	func= YUpsilonFlat;
	break;
      }
      if (sname=="MUON"){
	func= YUpsilonMUON;
	break;
      }
      if (sname=="RITMAN"){
	func=YUpsilonRitman;
	break;
      }
      if (sname=="KAREL"){
	func=YUpsilonKarel;
	break;
      }
      func=0;
      printf("<AliGenGSIlib::GetY> unknown parametrisation\n");
      break;

    case kJPsi:
      if (sname=="FLAT"){
	func= YJpsiFlat;
	break;
      }
      if (sname=="MUON"){
	func= YJpsiMUON;
	break;
      }
      //      if (sname=="SERGEI"){
      //	func= YJpsi;
      //	break;
      //      }
   
      func=0;
      printf("<AliGenGSIlib::GetY> unknown parametrisation\n");
      break;
     
    case kCharm: 
	func= YCharm;
	break;

    case kBeauty: 
	func= YBeauty;
	break;

    case kEta:
         func=YEtaPHOS;
         break;

    case kEtaprime:
         func=YEtaprimePHOS;
         break;

    case kOmega:
         func=YOmega;
         break;

    case kRho:
	 func=YRho;
	 break;

    case kKaon:
         func=YKaonPHOS;
         break;

    case kPion:
         func=YPion;
         break;

    case kPhi:
         func=YPhiPHOS;
         break;

	 //    case kLambda:
	 //         func=YLambda;
	 //         break;

    case kBaryons:
         func=YBaryons;
         break;

    default:
        func=0;
        printf("<AliGenGSIlib::GetY> unknown parametrisation\n");
    }
    return func;
}

GenFuncIp AliGenGSIlib::GetIp(Int_t param, const char * tname) const
{
// Return pointer to particle type parameterisation
   GenFuncIp func=0;
   TString sname(tname);

   switch (param) 
    {
    case kUpsilon:
	func= IpUpsilon;
	break;

    case kJPsi:
	func= IpJpsi;
	break;

    case kCharm: 
	func= IpCharm;
	break;

    case kBeauty: 
	func= IpBeauty;
	break;

    case kEta:
         func=IpEta;
         break;

    case kEtaprime:
         func=IpEtaprime;
         break;

    case kOmega:
         func=IpOmega;
         break;

    case kRho:
	 func=IpRho;
	 break;

    case kKaon:
         func=IpKaonPHOS;
         break;

    case kPion:
         func=IpPionPHOS;
         break;

    case kPhi:
         func=IpPhi;
         break;

	 //    case kLambda:
	 //         func=IpLambda;
	 //         break;

    case kBaryons:
         func=IpBaryons;
         break;

    default:
        func=0;
        printf("<AliGenGSIlib::GetIp> unknown parametrisation\n");
    }
    return func;
}












