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


#include <Riostream.h>
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "AliGenEMlib.h"


ClassImp(AliGenEMlib)

//Initializers for static members
Int_t AliGenEMlib::fgSelectedPtParam=AliGenEMlib::kPizero7TeVpp;
Int_t AliGenEMlib::fgSelectedCentrality=AliGenEMlib::kpp;
Int_t AliGenEMlib::fgSelectedV2Systematic=AliGenEMlib::kNoV2Sys;

Double_t AliGenEMlib::CrossOverLc(const double a, const double b, const double x){
  if(x<b-a/2) return 1.0;
  else if(x>b+a/2) return 0.0;
  else return cos(((x-b)/a+0.5)*TMath::Pi())/2+0.5;
}
Double_t AliGenEMlib::CrossOverRc(const double a, const double b, const double x){
  return 1-CrossOverLc(a,b,x);
}

const Double_t AliGenEMlib::fgkV2param[16][15] = {
  // charged pion                                                                                                                        cent, based on
  {  0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // pp no V2
  ,{ 6.792545e-02, 1.324227e+00, 2.071866e-02, 3.385157e+00, 5.183718e+00, 2.783082e-02, 5.000000e+00, 4.078812e+00, 5.444857e-02, 1.555569e-03, 1.099253e+00, 0, 1, 8.823293e-04, 1.310387e+00 }   // 0-5,  Francesco
  ,{ 1.219655e-01, 1.237698e+00, 2.595132e-02, 3.157085e+00, 3.697252e+00, 4.576262e-02, 1.103851e+00, 5.102592e+00, 1.072299e-01, 2.338602e-03, 8.475128e-01, 0, 1, 1.631554e-03, 9.432292e-01 }   // 5-10, Francesco
  ,{ 1.758981e-01, 1.253238e+00, 2.823197e-02, 3.455476e+00, 4.283204e+00, 6.053191e-02, 6.737159e-01, 4.420205e+00, 1.733567e-01, 3.292614e-03, 6.389341e-01, 0, 1, 3.234254e-03, 5.304084e-01 }   // 10-20,Francesco
  ,{ 2.231108e-01, 1.316781e+00, 3.537289e-02, 3.228243e+00, 3.697226e+00, 7.131631e-02, 5.902192e-01, 4.938186e+00, 2.123481e-01, 4.775233e-03, 5.707533e-01, 0, 1, 3.906720e-03, 5.261247e-01 }   // 20-30,Francesco
  ,{ 2.045341e-01, 1.324054e+00, 1.125395e-01, 2.179936e+00, 5.614603e+00, 8.946014e-02, 5.601938e-01, 2.959228e+00, 3.669442e-01, 5.318697e-03, 6.032741e-01, 0, 1, 4.231316e-03, 5.563485e-01 }   // 30-40,Francesco
  ,{ 2.074778e-01, 1.894322e+00, 6.019224e-02, 2.414500e+00, 3.822721e+00, 9.891310e-02, 4.976126e-01, 1.256500e+00, 5.210237e-01, 5.687290e-03, 6.785355e-01, 0, 1, 3.954051e-03, 6.913428e-01 }   // 40-50,Francesco
  ,{ 1.961099e-01, 2.122719e+00, 7.633636e-02, 1.731213e+00, 2.427536e+00, 1.026013e-01, 7.932261e-01, 4.082226e+00, 2.495494e-01, 5.677386e-03, 8.405024e-01, 0, 1, 2.961753e-03, 1.008892e+00 }   // 50-60,Francesco
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 0-10
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 20-40
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 40-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 60-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 0-20
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 0-40
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 20-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000 }   // 40-80
};

// MASS   0=>PIZERO, 1=>ETA, 2=>RHO, 3=>OMEGA, 4=>ETAPRIME, 5=>PHI, 6=>JPSI
const Double_t AliGenEMlib::fgkHM[7] = {0.13498, 0.54751, 0.7755, 0.78265, 0.95778, 1.01946, 3.0969};

const Double_t AliGenEMlib::fgkMtFactor[2][7] = { 
  // {1.0, 0.5, 1.0, 0.9, 0.4, 0.23, 0.054},  // factor for pp from arXiv:1110.3929
  // {1.0, 0.55, 1.0, 0.9, 0.4, 0.25, 0.004}    // factor for PbPb from arXiv:1110.3929
  //{1., 0.48, 1.0, 0.9, 0.25, 0.4}, (old values)
  //{1., 0.48, 1.0, 0.9, 0.4, 0.25}, (nlo values)
  //{1., 0.48, 1.0, 0.8, 0.4, 0.2, 0.06} (combination of nlo and LHC measurements)
  //https://aliceinfo.cern.ch/Figure/node/2634
  //https://aliceinfo.cern.ch/Figure/node/2788
  //https://aliceinfo.cern.ch/Figure/node/4403
  //J/Psi PbPb from Comparison with Julian Books J/Psi -> e+e-, might be contradicting with https://aliceinfo.cern.ch/Figure/node/3457
  //best guess:
  {1., 0.4, 1.0, 0.8, 0.4, 0.2, 0.0054}, //pp
  {1., 0.4, 1.0, 0.8, 0.4, 0.25, 0.01}  //PbPb
};

//==========================================================================
//
//              Definition of Particle Distributions
//                    
//==========================================================================

//--------------------------------------------------------------------------
//
//                              General functions
//
//--------------------------------------------------------------------------
Double_t AliGenEMlib::PtModifiedHagedornThermal(const Double_t pt,
                                                const Double_t c,
                                                const Double_t p0,
                                                const Double_t p1,
                                                const Double_t n,
                                                const Double_t cT,
                                                const Double_t T)
{
  // Modified Hagedorn Thermal fit to Picharged for PbPb:
  Double_t invYield;
  invYield = c/TMath::Power(p0+pt/p1,n) + cT*exp(-1.0*pt/T);

  return invYield*(2*TMath::Pi()*pt);
}



Double_t AliGenEMlib::PtModifiedHagedornExp(const Double_t pt,
					    const Double_t c,
					    const Double_t p1,
					    const Double_t p2,
					    const Double_t p0,
					    const Double_t n)
{
  // Modified Hagedorn exponentiel fit to Pizero for PbPb:
  Double_t invYield;
  invYield = c*TMath::Power(exp(-1*(p1*pt-p2*pt*pt))+pt/p0,-n);

  return invYield*(2*TMath::Pi()*pt);
}


Double_t AliGenEMlib::PtModifiedHagedornExp2(const Double_t pt,
                                             const Double_t c,
                                             const Double_t a,
                                             const Double_t b,
                                             const Double_t p0,
                                             const Double_t p1,
                                             const Double_t d,
                                             const Double_t n)
{
  // Modified Hagedorn exponential fit to charged pions for pPb:
  Double_t invYield;
  invYield = c*TMath::Power(exp(-a*pt-b*pt*pt)+pt/p0+TMath::Power(pt/p1,d),-n);

  return invYield*(2*TMath::Pi()*pt);
}

Double_t AliGenEMlib::PtTsallis(const Double_t pt,
                                const Double_t m,
                                const Double_t c,
                                const Double_t T,
                                const Double_t n)
{
  // Tsallis fit to Pizero for pp:
  Double_t mt;
  Double_t invYield;
 
  mt = sqrt(m*m + pt*pt);
  invYield = c*((n-1.)*(n-2.))/(n*T*(n*T+m*(n-2.)))*pow(1.+(mt-m)/(n*T),-n);

  return invYield*(2*TMath::Pi()*pt);
}



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
  // fit functions and corresponding parameter of Pizero pT for pp @ 2.76 TeV and @ 7 TeV and for PbPb @ 2.76 TeV 

  Double_t km=0.;
  Double_t kc=0.;
  Double_t kn=0.;
  Double_t kcT=0.;
  Double_t kT=0.;
  Double_t kp0=0.;
  Double_t kp1=0.;
  Double_t kp2=0.;
  Double_t ka=0.;
  Double_t kb=0.;
  Double_t kd=0.;

  switch(fgSelectedPtParam|fgSelectedCentrality) {
    // fit to pi charged v1
    // charged pion from ToF, unidentified hadrons scaled with pion from TPC
    // for Pb-Pb @ 2.76 TeV
  case kPichargedPbPb|k0005:
    kc=1347.5; kp0=0.9393; kp1=2.254; kn=11.294; kcT=0.002537; kT=2.414;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k0510:
    kc=1256.1; kp0=0.9545; kp1=2.248; kn=11.291; kcT=0.002662; kT=2.326;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k2030:
    kc=7421.6; kp0=1.2059; kp1=1.520; kn=10.220; kcT=0.002150; kT=2.196;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k3040:
    kc=1183.2; kp0=1.0478; kp1=1.623; kn=9.8073; kcT=0.00198333; kT=2.073;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
    // the following is what went into the Pb-Pb preliminary approval (0-10%)
  case kPichargedPbPb|k0010:
    kc=1296.0; kp0=0.968; kp1=2.567; kn=12.27; kcT=0.004219; kT=2.207;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k1020:
    kc=986.0; kp0=0.9752; kp1=2.376; kn=11.62; kcT=0.003116; kT=2.213;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k2040:
    kc=17337.0; kp0=1.337; kp1=1.507; kn=10.629; kcT=0.00184; kT=2.234;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k4050:
    kc=6220.0; kp0=1.322; kp1=1.224; kn=9.378; kcT=0.000595; kT=2.383;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k5060:
    kc=2319.0; kp0=1.267; kp1=1.188; kn=9.044; kcT=0.000437; kT=2.276;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k4060:
    kc=4724.0; kp0=1.319; kp1=1.195; kn=9.255; kcT=0.000511; kT=2.344;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPichargedPbPb|k6080:
    kc=2842.0; kp0=1.465; kp1=0.8324; kn=8.167; kcT=0.0001049; kT=2.29;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
   

    // fit to pizero from conversion analysis
    // for PbPb @ 2.76 TeV
    // Pi0 spectra --> not final results 
  case kPizeroPbPb|k0005:
       kc=1952.832; kp1=0.264; kp2=0.069; kp0=1.206; kn=9.732;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;
  case kPizeroPbPb|k0010:
       kc=1810.029; kp1=0.291; kp2=0.059; kp0=1.170; kn=9.447;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;
  case kPizeroPbPb|k0020:
       kc=856.241; kp1=-0.409; kp2=-0.127; kp0=1.219; kn=9.030;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;     
  case kPizeroPbPb|k1020:
       kc=509.781; kp1=-0.784; kp2=-0.120; kp0=0.931; kn=7.299;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;
  case kPizeroPbPb|k2040:
       kc=541.049; kp1=0.542; kp2=-0.069; kp0=0.972; kn=7.866;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;
  case kPizeroPbPb|k2080:
       kc=222.577; kp1=0.634; kp2=0.009; kp0=0.915; kn=7.431;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;
  case kPizeroPbPb|k4080:
       kc=120.153; kp1=0.7; kp2=-0.14; kp0=0.835; kn=6.980;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;
  case kPizeroPbPb|k0040:
       kc=560.532; kp1=0.548; kp2=-0.048; kp0=1.055; kn=8.132;
       return PtModifiedHagedornExp(*px,kc,kp1,kp2,kp0,kn);
       break;  
  
  
    // fit to charged pions for p-Pb @ 5.02TeV     
  case kPichargedPPb:
       kc=235.5; ka=0.6903; kb=0.06864; kp0=2.289; kp1=0.5872; kd=0.6474; kn=7.842; 
       return PtModifiedHagedornExp2(*px,kc,ka,kb,kp0,kp1,kd,kn);
       break;


    // Tsallis fit to final pizero (PHOS+PCM) -> used for publication
    // for pp @ 7 TeV
  case kPizero7TeVpp:
  case kPizeroEta7TeVpp:
    km=0.13498; kc=28.01; kT=0.139; kn=6.875;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
  case kPizero7TeVpplow:
  case kPizeroEta7TeVpplow: 
    km=0.13498; kc=23.84; kT=0.147; kn=7.025;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
  case kPizero7TeVpphigh:
  case kPizeroEta7TeVpphigh:
    km=0.13498; kc=32.47; kT=0.132; kn=6.749;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
    // Tsallis fit to pizero: preliminary result from PCM and PHOS (QM'11)
    // for pp @ 2.76 TeV
  case kPizero2760GeVpp:
  case kPizeroEta2760GeVpp:     
    km = 0.13498; kc = 19.75; kT = 0.130; kn = 7.051;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
  case kPizero2760GeVpplow:
  case kPizeroEta2760GeVpplow:
    km = 0.13498; kc = 16.12; kT = 0.142; kn = 7.327;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
  case kPizero2760GeVpphigh:
  case kPizeroEta2760GeVpphigh:
    km = 0.13498; kc = 25.18; kT = 0.118; kn = 6.782;
    return PtTsallis(*px,km,kc,kT,kn);
    break;

  default:
    return NULL;
  }

}

Double_t AliGenEMlib::YPizero( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2Pizero( const Double_t *px, const Double_t */*dummy*/ )
{
  int selectedCentrality;
  double n1,n2,v1,v2;
  switch(fgSelectedCentrality) {
  case k0010:
    selectedCentrality=fgSelectedCentrality;
    fgSelectedCentrality=k0005;
    n1=PtPizero(px,(Double_t*) 0);
    v1=V2Param(px,fgkV2param[fgSelectedCentrality]);
    fgSelectedCentrality=k0510;
    n2=PtPizero(px,(Double_t*) 0);
    v2=V2Param(px,fgkV2param[fgSelectedCentrality]);
    fgSelectedCentrality=selectedCentrality;
    return (n1*v1+n2*v2)/(n1+n2);
    break;
  case k2040:
    selectedCentrality=fgSelectedCentrality;
    fgSelectedCentrality=k2030;
    n1=PtPizero(px,(Double_t*) 0);
    v1=V2Param(px,fgkV2param[fgSelectedCentrality]);
    fgSelectedCentrality=k3040;
    n2=PtPizero(px,(Double_t*) 0);
    v2=V2Param(px,fgkV2param[fgSelectedCentrality]);
    fgSelectedCentrality=selectedCentrality;
    return (n1*v1+n2*v2)/(n1+n2);
    break;

  default:
    return V2Param(px,fgkV2param[fgSelectedCentrality]);
  }
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

  // fit functions and corresponding parameter of Eta pT for pp @ 2.76 TeV and @ 7 TeV
  // and mtscaled pT 

  Double_t km = 0.;
  Double_t kc = 0.;
  Double_t kT = 0.;
  Double_t kn = 0.;

  switch(fgSelectedPtParam){ 
    // Tsallis fit to final eta (PHOS+PCM) -> used for final publication
    // for pp @ 7 TeV
  case kPizeroEta7TeVpp:
    km = 0.547853; kc = 2.496; kT = 0.229; kn = 6.985;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
  case kPizeroEta7TeVpplow:
    km = 0.547853; kc = 1.970; kT = 0.253; kn = 7.591;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
  case kPizeroEta7TeVpphigh:
    km = 0.547853; kc = 3.060; kT = 0.212; kn = 6.578;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
    // Tsallis fit to preliminary eta (QM'11)
    // for pp @ 2.76 TeV
  case kPizeroEta2760GeVpp:
    km = 0.547853; kc = 1.971; kT = 0.188; kn = 6.308;
    return PtTsallis(*px,km,kc,kT,kn);
  case kPizeroEta2760GeVpplow:
    km = 0.547853; kc = 1.228; kT = 0.220; kn = 7.030;
    return PtTsallis(*px,km,kc,kT,kn);
    break;
  case kPizeroEta2760GeVpphigh:
    km = 0.547853; kc = 2.802; kT = 0.164; kn = 5.815;
    return PtTsallis(*px,km,kc,kT,kn);
    break;

  default:
  return MtScal(*px,1);
    break;

  }

}

Double_t AliGenEMlib::YEta( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlib::V2Eta( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,1); //V2Param(px,fgkV2param[1][fgSelectedV2Param]);
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

Double_t AliGenEMlib::V2Rho( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,2);
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

Double_t AliGenEMlib::V2Omega( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,3);
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

Double_t AliGenEMlib::V2Etaprime( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,4);
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

Double_t AliGenEMlib::V2Phi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,5);
}

//--------------------------------------------------------------------------
//
//                              Jpsi
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpJpsi(TRandom *)
{
  // Return phi pdg code
  return 443;
}

Double_t AliGenEMlib::PtJpsi( const Double_t *px, const Double_t */*dummy*/ )
{
  // Jpsi pT
  return MtScal(*px,6);
}

Double_t AliGenEMlib::YJpsi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlib::V2Jpsi( const Double_t *px, const Double_t */*dummy*/ )
{
  const int oldSys=fgSelectedV2Systematic;
  fgSelectedV2Systematic=kNoV2Sys;
  double ret=0;

  switch(oldSys){
  case kLoV2Sys: ret=0; break;
  case kNoV2Sys: ret=KEtScal(*px,6)/2; break;
  case kUpV2Sys: ret=KEtScal(*px,6); break;
  }

  fgSelectedV2Systematic=oldSys;
  return ret;
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

  Double_t scaledPt = sqrt(pt*pt + fgkHM[np]*fgkHM[np] - fgkHM[0]*fgkHM[0]);
  Double_t scaledYield = PtPizero(&scaledPt, (Double_t*) 0);

  //     VALUE MESON/PI AT 5 GeV/c
  Double_t NormPt = 5.;
  Double_t scaledNormPt = sqrt(NormPt*NormPt + fgkHM[np]*fgkHM[np] - fgkHM[0]*fgkHM[0]);

  Double_t norm = fgkMtFactor[int(bool(fgSelectedCentrality))][np] * (PtPizero(&NormPt, (Double_t*) 0) / PtPizero(&scaledNormPt, (Double_t*) 0));

  return norm*(pt/scaledPt)*scaledYield;
}

Double_t AliGenEMlib::KEtScal(const Double_t pt, Int_t np)
{
  const int nq=2; //number of quarks for particle np, here always 2
  Double_t scaledPt = sqrt(pow(2/nq*(sqrt(pt*pt+fgkHM[np]*fgkHM[np])-fgkHM[np])+fgkHM[0],2)-fgkHM[0]*fgkHM[0]);
  double val=V2Pizero(&scaledPt, (Double_t*) 0);
  static const double syserr[12]={0., 0.09, 0.07, 0.06, 0.04, 0.04, 0.04, 0.05, 0., 0., 0., 0.}; //based on pi vs kaon
  double sys=fgSelectedV2Systematic*min(fgkV2param[fgSelectedCentrality][0],fgkV2param[fgSelectedCentrality][8])*syserr[fgSelectedCentrality];
  return TMath::Max(val+sys,0.0);
}

Double_t AliGenEMlib::V2Param(const Double_t *px, const Double_t *par)
{
  // Very general parametrization of the v2

  const double &pt=px[0];
  double val=CrossOverLc(par[4],par[3],pt)*(2*par[0]/(1+TMath::Exp(par[1]*(par[2]-pt)))-par[0])+CrossOverRc(par[4],par[3],pt)*((par[8]-par[5])/(1+TMath::Exp(par[6]*(pt-par[7])))+par[5]);
  double sys=fgSelectedV2Systematic*par[11+fgSelectedV2Systematic*2]*pow(pt,par[12+fgSelectedV2Systematic*2]);
  return TMath::Max(val+sys,0.0);
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
    case kJpsi:
      func=PtJpsi;
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
    case kJpsi:
      func=YJpsi;
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
    case kJpsi:
      func=IpJpsi;
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
    case kJpsi:
      func=V2Jpsi;
      break;

    default:
      func=0;
      printf("<AliGenEMlib::GetV2> unknown parametrisation\n");
    }
  return func;
}

