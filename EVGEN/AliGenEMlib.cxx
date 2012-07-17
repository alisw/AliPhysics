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

AliGenEMlib::Centbins_t AliGenEMlib::fCentbin=k2030;

Double_t AliGenEMlib::CrossOverLc(const double a, const double b, const double x){
if(x<b-a/2) return 1.0;
  else if(x>b+a/2) return 0.0;
  else return cos(((x-b)/a+0.5)*TMath::Pi())/2+0.5;
}
Double_t AliGenEMlib::CrossOverRc(const double a, const double b, const double x){
  return 1-CrossOverLc(a,b,x);
}

const Double_t AliGenEMlib::fpTparam[8][14] = {
  // charged pion
  { 0.0005650000, 4.3540000000, 7.0070000000, 0.2472, 0.0000000000, 1.0000000000, 0.0005650000, 4.3540000000, 7.0070000000, 0.2472, 9999.0000000, 1.0000000000, 0.0000000000, 0.0000000000 },   //
  { 0.0005650000, 4.3540000000, 7.0070000000, 0.2472, 0.0000000000, 1.0000000000, 0.0005650000, 4.3540000000, 7.0070000000, 0.2472, 9999.0000000, 1.0000000000, 0.0000000000, 0.0000000000 },   //
  { 3.217746e+03, 1.632570e+00, 9.340162e+00, 1.0000, 2.650000e+00, 4.200635e+00, 9.039589e+08, 1.191548e-01, 6.907293e+00, 1.0000, 6.750000e+00, 4.099970e+00, 6.036772e+01, 5.928279e+00 },   //10-20
  { 2.125480e+03, 1.711882e+00, 9.585665e+00, 1.0000, 3.100000e+00, 5.041511e+00, 2.431808e+08, 1.155071e-01, 6.574663e+00, 1.0000, 6.250000e+00, 1.842070e+00, 3.928902e+01, 5.860970e+00 },   //20-30
  { 1.577897e+03, 1.411948e+00, 8.638815e+00, 1.0000, 2.550000e+00, 4.066432e+00, 3.774439e+08, 1.000330e-01, 6.535971e+00, 1.0000, 6.750000e+00, 2.482514e+00, 3.495494e+01, 5.954321e+00 },   //30-40
  { 7.061859e+02, 1.223810e+00, 8.532113e+00, 1.0000, 1.350000e+00, 1.956213e+00, 1.318169e+04, 5.658401e-01, 7.157575e+00, 1.0000, 5.250000e+00, 3.900000e+00, 1.958841e+01, 6.114787e+00 },   //40-50
  { 7.061859e+02, 1.223810e+00, 8.532113e+00, 1.0000, 1.350000e+00, 1.956213e+00, 1.318169e+04, 5.658401e-01, 7.157575e+00, 1.0000, 5.250000e+00, 3.900000e+00, 1.958841e+01, 6.114787e+00 },   //
  { 0.0005650000, 4.3540000000, 7.0070000000, 0.2472, 0.0000000000, 1.0000000000, 0.0005650000, 4.3540000000, 7.0070000000, 0.2472, 9999.0000000, 1.0000000000, 0.0000000000, 0.0000000000 } }; //pp

const Double_t AliGenEMlib::fv2param[2][8][9] = { { 
    // charged pion
    {  3.971545e-02, 2.111261e+00, 6.499885e-02, 3.659836e+00, 7.812234e+00, 1.479471e-02, 1.241708e+00, 2.817950e+00, 1.910663e-01 }   //0-5
    ,{ 7.424082e-02, 2.417588e+00, 7.891084e-02, 2.232841e+00, 3.398147e+00, 3.410206e-02, 4.138276e-01,-1.152430e+00, 5.729093e-01 }   //5-10
    ,{ 1.094608e-01, 2.357420e+00, 7.614904e-02, 2.294245e+00, 3.538364e+00, 4.932739e-02, 4.557926e-01, 9.218702e-03, 6.428540e-01 }   //10-20
    ,{ 1.456377e-01, 2.322408e+00, 7.747166e-02, 2.148424e+00, 3.238113e+00, 5.414722e-02, 4.042938e-01, 1.903040e-01, 6.816790e-01 }   //20-30
    ,{ 1.745154e-01, 2.234489e+00, 7.962225e-02, 2.007075e+00, 2.925536e+00, 5.499138e-02, 3.958957e-01, 1.780793e+00, 4.852930e-01 }   //30-40
    ,{ 2.384302e-01, 1.578935e+00, 9.234643e-02, 4.328926e+00, 1.020669e+01, 1.001340e-01, 2.433905e+00, 2.966673e+00, 2.646807e-01 }   //40-50
    ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000 }   //
    ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000 }   //pp no v2
  },{
    // charged kaon
    {  0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000 }  //
    ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000 }  //
    ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000 }  //
    ,{ 1.432314e-01, 1.607297e+00, 2.296868e-01, 1.821156e+00, 2.951652e+00,-2.201385e-01, 1.110419e-01, 8.228701e+00, 4.469488e-01 }  //20-30
    ,{ 1.691586e-01, 1.617165e+00, 2.159350e-01, 1.649338e+00, 2.607635e+00,-9.536279e+00, 3.860086e-01, 1.802996e+01, 2.509780e-01 }  //30-40
    ,{ 1.733831e-01, 1.712705e+00, 1.935862e-01, 1.745523e+00, 2.845436e+00,-5.772356e-01, 6.861616e-02, 5.292878e+00, 9.058815e-01 }  //40-50
    ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000 }  //
    ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000 }  //pp no v2
  } };

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
  // const Double_t kc=0.000565;
  // const Double_t kp0=0.2472;
  // const Double_t kp1=4.354;
  // const Double_t kn=7.007;
  Double_t invYield;
  Double_t x=*px;

  // invYield = kc/TMath::Power(kp0+x/kp1,kn);

  if(!x)return 0;
  const Double_t *c=fpTparam[fCentbin];
   
  invYield = CrossOverLc(c[5],c[4],x)*c[0]/TMath::Power(c[3]+x/c[1],c[2]);
  invYield +=CrossOverRc(c[5],c[4],x)*c[6]/TMath::Power(c[9]+x/c[7],c[8])*CrossOverLc(c[11],c[10],x);
  invYield +=CrossOverRc(c[11],c[10],x)*c[12]/TMath::Power(x,c[13]);

  return invYield*(2*TMath::Pi()*x);
}

Double_t AliGenEMlib::YPizero( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2Pizero( const Double_t *px, const Double_t */*dummy*/ )
{
  return V2Param(px,fv2param[0][fCentbin]);
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

Double_t AliGenEMlib::V2Eta( const Double_t *px, const Double_t */*dummy*/ )
{
  return V2Param(px,fv2param[1][fCentbin]);
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

Double_t AliGenEMlib::V2Param(const Double_t *px, const Double_t *par)
{
  // Very general parametrization of the v2

  return TMath::Max(CrossOverLc(par[4],par[3],px[0])*(2*par[0]/(1+TMath::Exp(par[1]*(par[2]-px[0])))-par[0])+CrossOverRc(par[4],par[3],px[0])*((par[8]-par[5])/(1+TMath::Exp(par[6]*(px[0]-par[7])))+par[5]),0.0);
}

Double_t AliGenEMlib::V2Flat(const Double_t */*px*/, const Double_t */*param*/)
{
  // Flat v2

  return 1.0;
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

GenFunc AliGenEMlib::GetV2(Int_t param, const char * tname) const
{
  // Return pointer to v2-parameterisation
  GenFunc func=0;
  TString sname(tname);

  switch (param) 
    {
    case kPizero:
      func=V2Pizero;
      break;
    case kEta:
      func=V2Eta;
      break;
    case kRho:
      func=V2Pizero;
      break;
    case kOmega:
      func=V2Pizero;
      break;
    case kEtaprime:
      func=V2Pizero;
      break;
    case kPhi:
      func=V2Pizero;
      break;

    default:
      func=0;
      printf("<AliGenEMlib::GetV2> unknown parametrisation\n");
    }
  return func;
}
