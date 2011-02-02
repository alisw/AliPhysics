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
  
/* $Id: AliGenEMlib.cxx 30052 2008-11-25 14:54:18Z morsch $ */

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Implementation of AliGenEMlib for electron, di-electron, and photon     //
// cocktail calculations.                                                  //
// It is based on AliGenGSIlib.                                            //
//                                                                         //
// Responsible: R.Averbeck@gsi.de                                          //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "AliGenEMlib.h"


ClassImp(AliGenEMlib)

//==========================================================================
//
//              Definition of Particle Distributions
//                    
//==========================================================================

//--------------------------------------------------------------------------
//
//                              Pizero
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpPizero(TRandom *)
{
// Return pizero pdg code
  return 111;     
}

Double_t AliGenEMlib::PtPizero( const Double_t *px, const Double_t */*dummy*/ )
{
// Generate pizero pT distribution from modified Hagedorn parameterization
// taken from fit to unidentified hadrons in pp at 7 TeV
  const Double_t kc=0.000565;
  const Double_t kp0=0.2472;
  const Double_t kp1=4.354;
  const Double_t kn=7.007;
  Double_t invYield;
  Double_t x=*px;

  invYield = kc/TMath::Power(kp0+x/kp1,kn);

  return invYield*(2*TMath::Pi()*x);
   
}

Double_t AliGenEMlib::YPizero( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

//--------------------------------------------------------------------------
//
//                              Eta
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpEta(TRandom *)
{
// Return eta pdg code
  return 221;     
}

Double_t AliGenEMlib::PtEta( const Double_t *px, const Double_t */*dummy*/ )
{
// Eta pT
  return MtScal(*px,1);
}

Double_t AliGenEMlib::YEta( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

//--------------------------------------------------------------------------
//
//                              Rho
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpRho(TRandom *)
{
// Return rho pdg code
  return 113;     
}

Double_t AliGenEMlib::PtRho( const Double_t *px, const Double_t */*dummy*/ )
{
// Rho pT
  return MtScal(*px,2);
}

Double_t AliGenEMlib::YRho( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

//--------------------------------------------------------------------------
//
//                              Omega
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpOmega(TRandom *)
{
// Return omega pdg code
  return 223;     
}

Double_t AliGenEMlib::PtOmega( const Double_t *px, const Double_t */*dummy*/ )
{
// Omega pT
  return MtScal(*px,3);
}

Double_t AliGenEMlib::YOmega( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

//--------------------------------------------------------------------------
//
//                              Etaprime
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpEtaprime(TRandom *)
{
// Return etaprime pdg code
  return 331;     
}

Double_t AliGenEMlib::PtEtaprime( const Double_t *px, const Double_t */*dummy*/ )
{
// Eta pT
  return MtScal(*px,4);
}

Double_t AliGenEMlib::YEtaprime( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

//--------------------------------------------------------------------------
//
//                              Phi
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpPhi(TRandom *)
{
// Return phi pdg code
  return 333;     
}

Double_t AliGenEMlib::PtPhi( const Double_t *px, const Double_t */*dummy*/ )
{
// Phi pT
  return MtScal(*px,5);
}

Double_t AliGenEMlib::YPhi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::YFlat(Double_t /*y*/)
{
//--------------------------------------------------------------------------
//
//                    flat rapidity distribution 
//
//--------------------------------------------------------------------------

  Double_t dNdy = 1.;   

  return dNdy;

}

//============================================================= 
//
//                    Mt-scaling  
//
//============================================================= 
//
 Double_t AliGenEMlib::MtScal(Double_t pt, Int_t np)
{
// Function for the calculation of the Pt distribution for a 
// given particle np, from the pizero Pt distribution using  
// mt scaling. 

// MASS   0=>PIZERO, 1=>ETA, 2=>RHO, 3=>OMEGA, 4=>ETAPRIME, 5=>PHI

  const Double_t khm[6] = {0.13498, 0.54751, 0.7755, 0.78265, 0.95778, 1.01946};

  Double_t scaledPt = sqrt(pt*pt + khm[np]*khm[np] - khm[0]*khm[0]);
  Double_t scaledYield = PtPizero(&scaledPt, (Double_t*) 0);

  //     VALUE MESON/PI AT 5 GEV

  Double_t normPt = 5.;
  Double_t scaledNormPt = sqrt(normPt*normPt + khm[np]*khm[np] - khm[0]*khm[0]);
  const Double_t kfmax[6]={1., 0.48, 1.0, 0.9, 0.25, 0.4};

  Double_t norm = kfmax[np] * (PtPizero(&normPt, (Double_t*) 0) / PtPizero(&scaledNormPt, (Double_t*) 0));

  return norm*(pt/scaledPt)*scaledYield;
}

//==========================================================================
//
//                     Set Getters 
//    
//==========================================================================

typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);

typedef Int_t (*GenFuncIp) (TRandom *);

GenFunc AliGenEMlib::GetPt(Int_t param, const char * tname) const
{
// Return pointer to pT parameterisation
   GenFunc func=0;
   TString sname(tname);

   switch (param) 
    {
    case kPizero:
      func=PtPizero;
      break;
    case kEta:
      func=PtEta;
      break;
    case kRho:
      func=PtRho;
      break;
    case kOmega:
      func=PtOmega;
      break;
    case kEtaprime:
      func=PtEtaprime;
      break;
    case kPhi:
      func=PtPhi;
      break;

    default:
      func=0;
      printf("<AliGenEMlib::GetPt> unknown parametrisation\n");
    }
   return func;
}

GenFunc AliGenEMlib::GetY(Int_t param, const char * tname) const
{
// Return pointer to y- parameterisation
   GenFunc func=0;
   TString sname(tname);

   switch (param) 
    {
    case kPizero:
         func=YPizero;
         break;
    case kEta:
         func=YEta;
         break;
    case kRho:
         func=YRho;
         break;
    case kOmega:
         func=YOmega;
         break;
    case kEtaprime:
         func=YEtaprime;
         break;
    case kPhi:
         func=YPhi;
         break;

    default:
        func=0;
        printf("<AliGenEMlib::GetY> unknown parametrisation\n");
    }
    return func;
}

GenFuncIp AliGenEMlib::GetIp(Int_t param, const char * tname) const
{
// Return pointer to particle type parameterisation
   GenFuncIp func=0;
   TString sname(tname);

   switch (param) 
    {
    case kPizero:
         func=IpPizero;
         break;
    case kEta:
         func=IpEta;
         break;
    case kRho:
         func=IpRho;
         break;
    case kOmega:
         func=IpOmega;
         break;
    case kEtaprime:
         func=IpEtaprime;
         break;
    case kPhi:
         func=IpPhi;
         break;

    default:
        func=0;
        printf("<AliGenEMlib::GetIp> unknown parametrisation\n");
    }
    return func;
}












