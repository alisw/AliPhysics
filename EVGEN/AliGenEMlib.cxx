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

Double_t AliGenEMlib::CrossOverLc(double a, double b, double x){
  if(x<b-a/2) return 1.0;
  else if(x>b+a/2) return 0.0;
  else return cos(((x-b)/a+0.5)*TMath::Pi())/2+0.5;
}
Double_t AliGenEMlib::CrossOverRc(double a, double b, double x){
  return 1-CrossOverLc(a,b,x);
}

const Double_t AliGenEMlib::fgkV2param[kCentralities][16] = {
  // charged pion                                                                                                                        cent, based on: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/FlowPAGQM2012talkIdentified
  {  0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // pp no V2
  ,{ 6.554571e-02, 1.436915e+00, 4.610598e-02, 2.554090e+00, 1.300948e+00, 2.970850e-02, 4.767877e+00, 4.228885e+00, 6.025959e-02, 1.570851e-03, 1.108941e+00, 0, 1, 1.715434e-03, 4.088070e-01, 25 }  // 0-5
  ,{ 1.171348e-01, 1.333067e+00, 4.537086e-02, 3.046348e+00, 3.903416e+00, 4.407152e-02, 9.123846e-01, 4.834531e+00, 1.186227e-01, 2.259463e-03, 8.916458e-01, 0, 1, 2.300647e-03, 4.531172e-01, 25 }  // 5-10
  ,{ 1.748434e-01, 1.285199e+00, 4.219881e-02, 4.018858e+00, 4.255082e+00, 7.955896e-03, 1.183264e-01,-9.329627e+00, 5.826570e-01, 3.368057e-03, 5.437668e-01, 0, 1, 3.178663e-03, 3.617552e-01, 25 }  // 10-20
  ,{ 2.149526e-01, 1.408792e+00, 5.062101e-02, 3.206279e+00, 3.988517e+00, 3.724655e-02, 1.995791e-01,-1.571536e+01, 6.494227e+00, 4.957874e-03, 4.903140e-01, 0, 1, 4.214626e-03, 3.457922e-01, 25 }  // 20-30
  ,{ 2.408942e-01, 1.477541e+00, 5.768983e-02, 3.333347e+00, 3.648508e+00,-2.044309e-02, 1.004145e-01,-2.386625e+01, 3.301913e+00, 5.666750e-03, 5.118686e-01, 0, 1, 4.626802e-03, 3.188974e-01, 25 }  // 30-40
  ,{ 2.495109e-01, 1.543680e+00, 6.217835e-02, 3.518863e+00, 4.557145e+00, 6.014553e-02, 1.491814e-01,-5.443647e+00, 5.403300e-01, 6.217285e-03, 5.580218e-01, 0, 1, 4.620486e-03, 3.792879e-01, 25 }  // 40-50
  ,{ 2.166399e-01, 1.931033e+00, 8.195390e-02, 2.226640e+00, 3.106649e+00, 1.058755e-01, 8.557791e-01, 4.006501e+00, 2.476449e-01, 6.720714e-03, 6.342966e-01, 0, 1, 4.449839e-03, 4.968750e-01, 25 }  // 50-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-10 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 20-40
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 40-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 60-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-20 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-40 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 20-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 40-80
};

const Double_t AliGenEMlib::fgkRawPtOfV2Param[kCentralities][10] = {
  {  0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // pp no V2
  ,{ 2.181446e+08, 9.412925e-01, 1.158774e-01, 3.020303e+01, 6.790828e+00, 9.999996e+01, 2.616827e+00, 3.980492e+00, 1.225169e+07, 5.575243e+00 } // 0-5
  ,{ 3.006215e+08, 9.511881e-01, 1.192788e-01, 2.981931e+01, 5.068175e+01, 9.999993e+01, 2.650635e+00, 4.073982e+00, 2.508045e+07, 5.621039e+00 } // 5-10
  ,{ 1.643438e+09, 9.604242e-01, 1.218512e-01, 2.912684e+01, 1.164242e+00, 9.999709e+01, 2.662326e+00, 4.027795e+00, 7.020810e+07, 5.696860e+00 } // 10-20
  ,{ 8.109985e+08, 9.421935e-01, 1.328020e-01, 2.655910e+01, 1.053677e+00, 9.999812e+01, 2.722949e+00, 3.964547e+00, 6.104096e+07, 5.694703e+00 } // 20-30
  ,{ 5.219789e+08, 9.417339e-01, 1.417541e-01, 2.518080e+01, 7.430803e-02, 9.303295e+01, 2.780227e+00, 3.909570e+00, 4.723116e+07, 5.778375e+00 } // 30-40
  ,{ 2.547159e+08, 9.481459e-01, 2.364858e-01, 1.689288e+01, 3.858883e+00, 6.352619e+00, 2.742270e+00, 3.855226e+00, 3.120535e+07, 5.878677e+00 } // 40-50
  ,{ 9.396097e+07, 9.304491e-01, 3.244940e-01, 1.291696e+01, 2.854367e+00, 6.325908e+00, 2.828258e+00, 4.133699e+00, 1.302739e+07, 5.977896e+00 } // 50-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-10 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 20-40
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 40-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 60-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-20 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-40 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 20-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 40-80
};

const Double_t AliGenEMlib::fgkPtParam[kCentralities][10] = {
  {  0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // pp no V2
  ,{ 7.641493e+01, 7.203468e-01, 3.651383e-01, 1.047542e+01, 3.494331e+00, 5.129019e+00, 3.081716e+00, 5.154525e+00, 3.065719e+01, 5.823718e+00 } // 0-5  
  ,{ 1.704676e+02, 7.852682e-01, 4.177172e-01, 1.014652e+01, 3.324409e+00, 4.825894e+00, 2.889738e+00, 4.772249e+00, 3.502627e+01, 5.938773e+00 } // 5-10 
  ,{ 1.823377e+02, 8.144309e-01, 4.291562e-01, 1.022767e+01, 3.585469e+00, 5.275078e+00, 3.144351e+00, 5.259097e+00, 2.675708e+01, 5.892506e+00 } // 10-20
  ,{ 4.851407e+02, 9.341151e-01, 4.716673e-01, 1.058090e+01, 4.681218e+00, 7.261284e+00, 3.883227e+00, 6.638627e+00, 1.562806e+01, 5.772127e+00 } // 20-30
  ,{ 3.157060e+01, 6.849451e-01, 4.868669e-01, 8.394558e+00, 3.539142e+00, 5.495280e+00, 4.102638e+00, 3.722991e+00, 1.638622e+01, 5.935963e+00 } // 30-40
  ,{ 1.069397e+01, 5.816587e-01, 6.542961e-01, 6.472858e+00, 2.643870e+00, 3.929020e+00, 3.339224e+00, 2.410371e+00, 9.606748e+00, 6.116685e+00 } // 40-50
  ,{ 1.857919e+01, 6.185989e-01, 5.878869e-01, 7.035064e+00, 2.892415e+00, 4.339383e+00, 3.549679e+00, 2.821061e+00, 1.529318e+01, 6.091388e+00 } // 50-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-10 
  ,{ 1.271594e+02, 7.790165e-01, 5.793214e-01, 8.050008e+00, 3.211312e+00, 4.825258e+00, 3.840509e+00, 3.046231e+00, 2.172177e+01, 5.983496e+00 } // 20-40
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 40-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 60-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-20 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-40 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 20-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 40-80
};

const Double_t AliGenEMlib::fgkThermPtParam[kCentralities][2] = {
  {  0.0000000000, 0.0000000000 } // pp no V2
  ,{ 0.0000000000, 0.0000000000 } // 0-5
  ,{ 0.0000000000, 0.0000000000 } // 5-10
  ,{ 3.335661e+01, 3.400585e+00 } // 10-20 //based on: https://aliceinfo.cern.ch/Notes/node/249
  ,{ 0.0000000000, 0.0000000000 } // 20-30
  ,{ 0.0000000000, 0.0000000000 } // 30-40
  ,{ 0.0000000000, 0.0000000000 } // 40-50
  ,{ 0.0000000000, 0.0000000000 } // 50-60
  ,{ 3.648327e+02, 4.477749e+00 } // 0-10  //based on: https://aliceinfo.cern.ch/Notes/node/249
  ,{ 1.696223e+00, 2.429660e+00 } // 20-40 //based on: https://twiki.cern.ch/twiki/pub/ALICE/ALICEDirectPhotonSpectrumPaper/directPbPb.pdf
  ,{ 0.0000000000, 0.0000000000 } // 40-60
  ,{ 0.0000000000, 0.0000000000 } // 60-80
  ,{ 1.492160e+01, 2.805213e+00 } // 0-20  //based on: https://twiki.cern.ch/twiki/pub/ALICE/ALICEDirectPhotonSpectrumPaper/directPbPb.pdf
  ,{ 4.215110e+01, 3.242719e+00 } // 0-40  //based on: https://aliceinfo.cern.ch/Figure/node/2866
  ,{ 0.0000000000, 0.0000000000 } // 20-80
  ,{ 0.0000000000, 0.0000000000 } // 40-80
};

// MASS   0=>PIZERO, 1=>ETA, 2=>RHO, 3=>OMEGA, 4=>ETAPRIME, 5=>PHI, 6=>JPSI
const Double_t AliGenEMlib::fgkHM[8] = {0.13498, 0.54751, 0.7755, 0.78265, 0.95778, 1.01946, 3.0969, 0.0};

const Double_t AliGenEMlib::fgkMtFactor[2][8] = { 
  // {1.0, 0.5, 1.0, 0.9, 0.4, 0.23, 0.054},  // factor for pp from arXiv:1110.3929
  // {1.0, 0.55, 1.0, 0.9, 0.4, 0.25, 0.004}    // factor for PbPb from arXiv:1110.3929
  //{1., 0.48, 1.0, 0.9, 0.25, 0.4}, (old values)
  //{1., 0.48, 1.0, 0.9, 0.4, 0.25}, (nlo values)
  //{1., 0.48, 1.0, 0.8, 0.4, 0.2, 0.06} (combination of nlo and LHC measurements)
  //https://aliceinfo.cern.ch/Figure/node/2634
  //https://aliceinfo.cern.ch/Figure/node/2788
  //https://aliceinfo.cern.ch/Figure/node/4403
  //https://aliceinfo.cern.ch/Notes/node/87
  //best guess:
  {1., 0.48, 1.0, 0.9, 0.4, 0.25, 0., 0.}, //pp
  {1., 0.48, 1.0, 0.9, 0.4, 0.25, 0., 0.}  //PbPb
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
Double_t AliGenEMlib::PtModifiedHagedornThermal(Double_t pt,
                                                Double_t c,
                                                Double_t p0,
                                                Double_t p1,
                                                Double_t n,
                                                Double_t cT,
                                                Double_t T)
{
  // Modified Hagedorn Thermal fit to Picharged for PbPb:
  Double_t invYield;
  invYield = c/TMath::Power(p0+pt/p1,n) + cT*exp(-1.0*pt/T);

  return invYield*(2*TMath::Pi()*pt);
}



Double_t AliGenEMlib::PtModifiedHagedornExp(Double_t pt,
					    Double_t c,
					    Double_t p1,
					    Double_t p2,
					    Double_t p0,
					    Double_t n)
{
  // Modified Hagedorn exponentiel fit to Pizero for PbPb:
  Double_t invYield;
  invYield = c*TMath::Power(exp(-1*(p1*pt-p2*pt*pt))+pt/p0,-n);

  return invYield*(2*TMath::Pi()*pt);
}


Double_t AliGenEMlib::PtModifiedHagedornExp2(Double_t pt,
                                             Double_t c,
                                             Double_t a,
                                             Double_t b,
                                             Double_t p0,
                                             Double_t p1,
                                             Double_t d,
                                             Double_t n)
{
  // Modified Hagedorn exponential fit to charged pions for pPb:
  Double_t invYield;
  invYield = c*TMath::Power(exp(-a*pt-b*pt*pt)+pt/p0+TMath::Power(pt/p1,d),-n);

  return invYield*(2*TMath::Pi()*pt);
}

Double_t AliGenEMlib::PtTsallis(Double_t pt,
                                Double_t m,
                                Double_t c,
                                Double_t T,
                                Double_t n)
{
  // Tsallis fit to Pizero for pp:
  Double_t mt;
  Double_t invYield;
 
  mt = sqrt(m*m + pt*pt);
  invYield = c*((n-1.)*(n-2.))/(n*T*(n*T+m*(n-2.)))*pow(1.+(mt-m)/(n*T),-n);

  return invYield*(2*TMath::Pi()*pt);
}

// Exponential
Double_t AliGenEMlib::PtExponential(const Double_t *px, const Double_t *c){
  const double &pt=px[0];
  Double_t invYield = c[0]*exp(-pt*c[1]);
  
  return invYield*(2*TMath::Pi()*pt);
}

// Hagedorn with additional Powerlaw
Double_t AliGenEMlib::PtModifiedHagedornPowerlaw(const Double_t *px, const Double_t *c){
  const double &pt=px[0];
  Double_t invYield = c[0]*pow(c[1]+pt*c[2],-c[3])*CrossOverLc(c[5],c[4],pt)+CrossOverRc(c[7],c[6],pt)*c[8]*pow(pt+0.001,-c[9]); //pt+0.001: prevent powerlaw from exploding for pt->0
  
  return invYield*(2*TMath::Pi()*pt+0.001); //pt+0.001: be sure to be > 0
}

// double powerlaw for J/Psi yield
Double_t AliGenEMlib::PtDoublePowerlaw(const Double_t *px, const Double_t *c){
  const double &pt=px[0];
  Double_t yield = c[0]*pt*pow(1+pow(pt*c[1],2),-c[2]);
  
  return yield;
}

// integral over krollwada with S=1^2*(1-mee^2/mh^2)^3 from mee=0 up to mee=mh
// approximation is perfect for mh>20MeV
Double_t AliGenEMlib::IntegratedKrollWada(const Double_t *mh, const Double_t *){
  if(*mh<0.002941) return 0;
  return 2*log(*mh/0.000511/exp(1.75))/411.11/TMath::Pi();
}

//--------------------------------------------------------------------------
//
//                             DirectRealGamma
//
//--------------------------------------------------------------------------
Double_t AliGenEMlib::PtPromptRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  const static Double_t promptGammaPtParam[10] = { 8.715017e-02, 4.439243e-01, 1.011650e+00, 5.193789e+00, 2.194442e+01, 1.062124e+01, 2.469876e+01, 6.052479e-02, 5.611410e-02, 5.169743e+00 };
 
  return PtModifiedHagedornPowerlaw(px,promptGammaPtParam)*GetTAA(fgSelectedCentrality);
}

Double_t AliGenEMlib::PtThermalRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return PtExponential(px,fgkThermPtParam[fgSelectedCentrality]);
}

Double_t AliGenEMlib::PtDirectRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return PtPromptRealGamma(px,px)+PtThermalRealGamma(px,px);
}

Int_t AliGenEMlib::IpDirectRealGamma(TRandom *)
{
  return 22;
}

Double_t AliGenEMlib::YDirectRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return YFlat(*px);
}

Double_t AliGenEMlib::V2DirectRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  const static Double_t v2Param[3][16] = {
    { 1.004245e-01, 1.057645e+00, 0.000000e+00, 2.836492e+00, 2.819767e+00, -6.231529e-02, 1.173054e+00, 2.836492e+00, 1.881590e-01, 1.183293e-02, 1.252249e+00, 0, 1, 4.876263e-03, 1.518526e+00, 4.5 } // 00-20, based on: https://aliceinfo.cern.ch/Notes/node/249
    ,{ 1.619000e-01, 1.868201e+00, 6.983303e-15, 2.242170e+00, 4.484339e+00, -1.695734e-02, 2.301359e+00, 2.871469e+00, 1.619000e-01, 2.264320e-02, 1.028641e+00, 0, 1, 8.172203e-03, 1.271637e+00, 4.5 } // 20-40
    ,{ 1.335000e-01, 1.076916e+00, 1.462605e-08, 2.785732e+00, 5.571464e+00, -2.356156e-02, 2.745437e+00, 2.785732e+00, 1.335000e-01, 1.571589e-02, 1.001131e+00, 0, 1, 5.179715e-03, 1.329344e+00, 4.5 } // 00-40
  };
  switch(fgSelectedCentrality){
  case k0020: return V2Param(px,v2Param[0]); break;
  case k2040: return V2Param(px,v2Param[1]); break;
  case k0040: return V2Param(px,v2Param[2]); break;
    // case k0010: return 0.43*V2Param(px,v2Param[1]); break;  //V2Pizero(0010)/V2Pizero(2040)=0.43 +-0.025
    // case k1020: return 0.75*V2Param(px,v2Param[1]); break;  //V2Pizero(1020)/V2Pizero(2040)=0.75 +-0.04
  case k0010: return 0.53*V2Param(px,v2Param[2]); break;  //V2Pizero(0010)/V2Pizero(0040)=0.53 +-0.03
  case k1020: return 0.91*V2Param(px,v2Param[2]); break;  //V2Pizero(1020)/V2Pizero(0040)=0.91 +-0.04
  }
  return 0;
}


//--------------------------------------------------------------------------
//
//                             DirectVirtGamma
//
//--------------------------------------------------------------------------
Double_t AliGenEMlib::PtPromptVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return IntegratedKrollWada(px,px)*PtPromptRealGamma(px,px);
}

Double_t AliGenEMlib::PtThermalVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return IntegratedKrollWada(px,px)*PtThermalRealGamma(px,px);
}

Double_t AliGenEMlib::PtDirectVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return IntegratedKrollWada(px,px)*PtDirectRealGamma(px,px);
}

Int_t AliGenEMlib::IpDirectVirtGamma(TRandom *)
{
  return 220000;
}

Double_t AliGenEMlib::YDirectVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return YFlat(*px);
}

Double_t AliGenEMlib::V2DirectVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return V2DirectRealGamma(px,px);
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
  // double pigammacorr=1; //misuse pion for direct gammas, tuned for 0040, iteration 0
  // pigammacorr*=2.258900e-01*log(*px+0.001)+1.591291e+00;  //iteration 1
  // pigammacorr*=6.601943e-03*log(*px+0.001)+9.593698e-01;  //iteration 2
  // pigammacorr*=4.019933e-03*log(*px+0.001)+9.843412e-01;  //iteration 3
  // pigammacorr*=-4.543991e-03*log(*px+0.001)+1.010886e+00; //iteration 4
  // return pigammacorr*PtPromptRealGamma(px,px); //now the gammas from the pi->gg decay have the pt spectrum of prompt real gammas
  
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

  double n1,n2,n3,n4;
  int oldCent;

  switch(fgSelectedPtParam|fgSelectedCentrality) {
    // fit to pi charged, same data like in kPiOldChargedPbPb,
    // but tested and compared against newest (2014) neutral pi measurement
  case kPichargedPbPb|k0005:
  case kPichargedPbPb|k0510:
  case kPichargedPbPb|k1020:
  case kPichargedPbPb|k2030:
  case kPichargedPbPb|k3040:
  case kPichargedPbPb|k4050:
  case kPichargedPbPb|k5060:
  case kPichargedPbPb|k2040:
    return PtModifiedHagedornPowerlaw(px,fgkPtParam[fgSelectedCentrality]);
    break;
  case kPichargedPbPb|k0010:
    n1=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0005]);
    n2=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0510]);
    return (n1+n2)/2;
    break;
  case kPichargedPbPb|k0020:
    n1=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0005]);
    n2=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0510]);
    n3=PtModifiedHagedornPowerlaw(px,fgkPtParam[k1020]);
    return (n1+n2+2*n3)/4;
    break;
  case kPichargedPbPb|k0040:
    n1=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0005]);
    n2=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0510]);
    n3=PtModifiedHagedornPowerlaw(px,fgkPtParam[k1020]);
    n4=PtModifiedHagedornPowerlaw(px,fgkPtParam[k2040]);
    return (n1+n2+2*n3+4*n4)/8;
    break;
  case kPichargedPbPb|k4060:
    n1=PtModifiedHagedornPowerlaw(px,fgkPtParam[k4050]);
    n2=PtModifiedHagedornPowerlaw(px,fgkPtParam[k5060]);
    return (n1+n2)/2;
    break;


    // fit to pi charged v1
    // charged pion from ToF, unidentified hadrons scaled with pion from TPC
    // for Pb-Pb @ 2.76 TeV
  case kPiOldChargedPbPb|k0005:
    kc=1347.5; kp0=0.9393; kp1=2.254; kn=11.294; kcT=0.002537; kT=2.414;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k0510:
    kc=1256.1; kp0=0.9545; kp1=2.248; kn=11.291; kcT=0.002662; kT=2.326;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k2030:
    kc=7421.6; kp0=1.2059; kp1=1.520; kn=10.220; kcT=0.002150; kT=2.196;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k3040:
    kc=1183.2; kp0=1.0478; kp1=1.623; kn=9.8073; kcT=0.00198333; kT=2.073;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
    // the following is what went into the Pb-Pb preliminary approval (0-10%)
  case kPiOldChargedPbPb|k0010:
    kc=1296.0; kp0=0.968; kp1=2.567; kn=12.27; kcT=0.004219; kT=2.207;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k1020:
    kc=986.0; kp0=0.9752; kp1=2.376; kn=11.62; kcT=0.003116; kT=2.213;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k2040:
    kc=17337.0; kp0=1.337; kp1=1.507; kn=10.629; kcT=0.00184; kT=2.234;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k4050:
    kc=6220.0; kp0=1.322; kp1=1.224; kn=9.378; kcT=0.000595; kT=2.383;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k5060:
    kc=2319.0; kp0=1.267; kp1=1.188; kn=9.044; kcT=0.000437; kT=2.276;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k4060:
    kc=4724.0; kp0=1.319; kp1=1.195; kn=9.255; kcT=0.000511; kT=2.344;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k6080:
    kc=2842.0; kp0=1.465; kp1=0.8324; kn=8.167; kcT=0.0001049; kT=2.29;
    return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
    break;
  case kPiOldChargedPbPb|k0020:
    oldCent=fgSelectedCentrality;
    fgSelectedCentrality=k0010;
    n1=PtPizero(px,px);
    fgSelectedCentrality=k1020;
    n2=PtPizero(px,px);
    fgSelectedCentrality=oldCent;
    return (n1+n2)/2;
    break;
  case kPiOldChargedPbPb|k0040:
    oldCent=fgSelectedCentrality;
    fgSelectedCentrality=k0010;
    n1=PtPizero(px,px);
    fgSelectedCentrality=k1020;
    n2=PtPizero(px,px);
    fgSelectedCentrality=k2040;
    n3=PtPizero(px,px);
    fgSelectedCentrality=oldCent;
    return (n1+n2+2*n3)/4;
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
  double n1,n2,n3,n4,n5;
  double v1,v2,v3,v4,v5;
  switch(fgSelectedCentrality) {
  case k0010:
    n1=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0005]);
    v1=V2Param(px,fgkV2param[k0005]);
    n2=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0510]);
    v2=V2Param(px,fgkV2param[k0510]);
    return (n1*v1+n2*v2)/(n1+n2);
    break;
  case k0020:
    n1=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0005]);
    v1=V2Param(px,fgkV2param[k0005]);
    n2=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0510]);
    v2=V2Param(px,fgkV2param[k0510]);
    n3=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k1020]);
    v3=V2Param(px,fgkV2param[k1020]);
    return (n1*v1+n2*v2+2*n3*v3)/(n1+n2+2*n3);
    break;
  case k2040:
    n1=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k2030]);
    v1=V2Param(px,fgkV2param[k2030]);
    n2=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k3040]);
    v2=V2Param(px,fgkV2param[k3040]);
    return (n1*v1+n2*v2)/(n1+n2);
    break;
  case k0040:
    n1=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0005]);
    v1=V2Param(px,fgkV2param[k0005]);
    n2=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0510]);
    v2=V2Param(px,fgkV2param[k0510]);
    n3=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k1020]);
    v3=V2Param(px,fgkV2param[k1020]);
    n4=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k2030]);
    v4=V2Param(px,fgkV2param[k2030]);
    n5=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k3040]);
    v5=V2Param(px,fgkV2param[k3040]);
    return (n1*v1+n2*v2+2*n3*v3+2*n4*v4+2*n5*v5)/(n1+n2+2*n3+2*n4+2*n5);
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
  // based on: //https://aliceinfo.cern.ch/Notes/node/242, https://aliceinfo.cern.ch/Figure/node/3457, www.sciencedirect.com/science/article/pii/S0370269312011446
  const static Double_t jpsiPtParam[2][3] = {
    {  9.686337e-03, 2.629441e-01, 4.552044e+00 }
    ,{ 3.403549e-03, 2.897061e-01, 3.644278e+00 }
  };
  const double pt=px[0]*2.28/2.613;
  switch(fgSelectedCentrality) {
  case k0020: return 2.405*PtDoublePowerlaw(&pt,jpsiPtParam[0]); break;
  case k2040: return 2.405*PtDoublePowerlaw(&pt,jpsiPtParam[1]); break;
  case k0040: return 0.5*2.405*(PtDoublePowerlaw(&pt,jpsiPtParam[0])+PtDoublePowerlaw(&pt,jpsiPtParam[1])); break;
  }
  return 0;
}

Double_t AliGenEMlib::YJpsi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlib::V2Jpsi( const Double_t *px, const Double_t */*dummy*/ )
{
  const static Double_t v2Param[16] = { 1.156000e-01, 8.936854e-01, 0.000000e+00, 4.000000e+00, 6.222375e+00, -1.600314e-01, 8.766676e-01, 7.824143e+00, 1.156000e-01, 3.484503e-02, 4.413685e-01, 0, 1, 3.484503e-02, 4.413685e-01, 7.2 };
  switch(fgSelectedCentrality){
  case k2040: return V2Param(px,v2Param); break;
  case k0010: return 0.43*V2Param(px,v2Param); break;  //V2Pizero(0010)/V2Pizero(2040)=0.43 +-0.025
  case k1020: return 0.75*V2Param(px,v2Param); break;  //V2Pizero(1020)/V2Pizero(2040)=0.75 +-0.04
  case k0020: return 0.66*V2Param(px,v2Param); break;  //V2Pizero(0020)/V2Pizero(2040)=0.66 +-0.035
  case k0040: return 0.82*V2Param(px,v2Param); break;  //V2Pizero(0040)/V2Pizero(2040)=0.82 +-0.05
  }
  return 0;
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

Double_t AliGenEMlib::KEtScal(Double_t pt, Int_t np)
{
  const int nq=2; //number of quarks for particle np, here always 2
  Double_t scaledPt = sqrt(pow(2.0/nq*(sqrt(pt*pt+fgkHM[np]*fgkHM[np])-fgkHM[np])+fgkHM[0],2)-fgkHM[0]*fgkHM[0]);
  return V2Pizero(&scaledPt, (Double_t*) 0);
}

Double_t AliGenEMlib::V2Param(const Double_t *px, const Double_t *par)
{
  // Very general parametrization of the v2

  const double &pt=px[0];
  double val=CrossOverLc(par[4],par[3],pt)*(2*par[0]/(1+TMath::Exp(par[1]*(par[2]-pt)))-par[0])+CrossOverRc(par[4],par[3],pt)*((par[8]-par[5])/(1+TMath::Exp(par[6]*(pt-par[7])))+par[5]);
  double sys=0;
  if(fgSelectedV2Systematic){
    double syspt=pt>par[15]?par[15]:pt;
    sys=fgSelectedV2Systematic*par[11+fgSelectedV2Systematic*2]*pow(syspt,par[12+fgSelectedV2Systematic*2]);
  }
  return std::max(val+sys,0.0);
}

Double_t AliGenEMlib::V2Flat(const Double_t */*px*/, const Double_t */*param*/)
{
  // Flat v2

  return 0.0;
}

Double_t AliGenEMlib::GetTAA(Int_t cent){
  const static Double_t taa[16] = { 1.0,    // pp
				    26.32,  // 0-5
				    20.56,  // 5-10
				    14.39,  // 10-20
				    8.70,   // 20-30
				    5.001,  // 30-40
				    2.675,  // 40-50
				    1.317,  // 50-60
				    23.44,  // 0-10
				    6.85,   // 20-40
				    1.996,  // 40-60
				    0.4174, // 60-80
				    18.91,  // 0-20
				    12.88,  // 0-40
				    3.088,  // 20-80
				    1.207}; // 40-80
  return taa[cent];  
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
     case kDirectRealGamma:
       func=PtDirectRealGamma;
      break;
     case kDirectVirtGamma:
       func=PtDirectVirtGamma;
      break;
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
    case kDirectRealGamma:
      func=YDirectRealGamma;
      break;
    case kDirectVirtGamma:
      func=YDirectVirtGamma;
      break;
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
    case kDirectRealGamma:
      func=IpDirectRealGamma;
      break;
    case kDirectVirtGamma:
      func=IpDirectVirtGamma;
      break;
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
    case kDirectRealGamma:
      func=V2DirectRealGamma;
      break;
    case kDirectVirtGamma:
      func=V2DirectVirtGamma;
      break;
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

