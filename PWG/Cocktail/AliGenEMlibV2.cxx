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

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Implementation of AliGenEMlibV2 for electron, di-electron, and photon   //
// cocktail calculations.                                                  //
// It is based on AliGenEMlib                                              //
//                                                                         //
// Responsible: Friederike Bock (friederike.bock@cern.ch)                  //
//              Lucas Altenkaemper (lucas.altenkamper@cern.ch)             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TFile.h"
#include "TFormula.h"
#include "AliLog.h"
#include "AliGenEMlibV2.h"
#include "TH1D.h"

using std::cout;
using std::endl;

ClassImp(AliGenEMlibV2)

//Initializers for static members
TF1*  AliGenEMlibV2::fPtParametrization[]       = {0x0};
TF1*  AliGenEMlibV2::fPtParametrizationProton   = NULL;
TH1D* AliGenEMlibV2::fMtFactorHisto             = NULL;
TH2F* AliGenEMlibV2::fPtYDistribution[]         = {0x0};
Int_t AliGenEMlibV2::fgSelectedCollisionsSystem = AliGenEMlibV2::kpp7TeV;
Int_t AliGenEMlibV2::fgSelectedCentrality       = AliGenEMlibV2::kpp;
Int_t AliGenEMlibV2::fgSelectedV2Systematic     = AliGenEMlibV2::kNoV2Sys;

Double_t AliGenEMlibV2::CrossOverLc(double a, double b, double x){
  if(x<b-a/2) return 1.0;
  else if(x>b+a/2) return 0.0;
  else return cos(((x-b)/a+0.5)*TMath::Pi())/2+0.5;
}

Double_t AliGenEMlibV2::CrossOverRc(double a, double b, double x){
  return 1-CrossOverLc(a,b,x);
}

const Double_t AliGenEMlibV2::fgkV2param[kCentralities][16] = {
  // charged pion                                                                                                                        cent, based on: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/FlowPAGQM2012talkIdentified
  {  0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10  }  // pp no V2
  ,{ 6.516808e-02, 1.449877e+00, 4.705013e-02, 3.025555e+00, 3.391947e+00, 2.364841e-02, 2.303211e+00, 4.443068e+00, 5.989572e-02, 2.079546e-03, 8.365175e-01, 0, 1, 2.079546e-03, 8.365175e-01, 6.0 }  // 0-5
  ,{ 1.199628e-01, 1.288148e+00, 4.230109e-02, 3.071488e+00, 3.649381e+00, 4.643994e-02, 1.210442e+00, 5.055715e+00, 1.089290e-01, 2.561783e-03, 7.907385e-01, 0, 1, 2.561783e-03, 7.907385e-01, 7.0 }  // 5-10
  ,{ 1.672604e-01, 1.367064e+00, 4.704060e-02, 3.005428e+00, 3.816287e+00, 5.845636e-02, 8.660308e-01, 4.587477e+00, 1.819530e-01, 2.945173e-03, 7.256487e-01, 0, 1, 2.945173e-03, 7.256487e-01, 8.0 }  // 10-20
  ,{ 2.079093e-01, 1.480510e+00, 5.518112e-02, 2.845102e+00, 3.485697e+00, 3.363945e-02, 3.599748e-01, 3.207685e+00, 3.494331e-01, 4.398452e-03, 6.399492e-01, 0, 1, 4.398452e-03, 6.399492e-01, 16.0 }  // 20-30
  ,{ 2.275701e-01, 1.622512e+00, 6.793718e-02, 2.742464e+00, 3.229573e+00, 3.316905e-02, 3.128005e-01, 2.465896e+00, 4.296590e-01, 4.837696e-03, 6.855293e-01, 0, 1, 4.837696e-03, 6.855293e-01, 17.0 }  // 30-40
  ,{ 2.394412e-01, 1.663730e+00, 7.163620e-02, 2.406896e+00, 2.343323e+00, 1.033637e-01, 1.028906e+00, 4.512137e+00, 2.383232e-01, 5.013671e-03, 7.776760e-01, 0, 1, 5.013671e-03, 7.776760e-01, 8.0 }  // 40-50
  ,{ 2.115576e-01, 1.981979e+00, 8.221619e-02, 2.203373e+00, 3.213309e+00, 1.113934e-01, 7.620465e-01, 3.575663e+00, 2.673220e-01, 4.555673e-03, 9.927123e-01, 0, 1, 4.555673e-03, 9.927123e-01, 7.0 }  // 50-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-10
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 20-40
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 40-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 60-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-20
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-40
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 20-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 40-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 20-50
};

const Double_t AliGenEMlibV2::fgkRawPtOfV2Param[kCentralities][10] = {
   { 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // pp no V2
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
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 20-50
};

const Double_t AliGenEMlibV2::fgkThermPtParam[kCentralities][2] = {
  // might also be interesting: https://aliceinfo.cern.ch/Notes/node/226
   { 0.0000000000, 0.0000000000 } // pp no V2
  ,{ 0.0000000000, 0.0000000000 } // 0-5
  ,{ 0.0000000000, 0.0000000000 } // 5-10
  ,{ 1.944330e+01, 3.047106e+00 } // 10-20 //based on: https://aliceinfo.cern.ch/Notes/node/87
  ,{ 0.0000000000, 0.0000000000 } // 20-30
  ,{ 0.0000000000, 0.0000000000 } // 30-40
  ,{ 0.0000000000, 0.0000000000 } // 40-50
  ,{ 0.0000000000, 0.0000000000 } // 50-60
  ,{ 3.557511e+02, 4.459934e+00 } // 0-10  //based on: https://aliceinfo.cern.ch/Notes/node/87
  ,{ 8.088291e-01, 2.013231e+00 } // 20-40 //based on: https://twiki.cern.ch/twiki/bin/view/ALICE/ALICEDirectPhotonSpectrumPaper
  ,{ 0.0000000000, 0.0000000000 } // 40-60
  ,{ 0.0000000000, 0.0000000000 } // 60-80
  ,{ 1.363037e+01, 2.696863e+00 } // 0-20  //based on: https://twiki.cern.ch/twiki/bin/view/ALICE/ALICEDirectPhotonSpectrumPaper
  ,{ 4.351422e+01, 3.267624e+00 } // 0-40  //based on: https://aliceinfo.cern.ch/Figure/node/2866
  ,{ 0.0000000000, 0.0000000000 } // 20-80
  ,{ 0.0000000000, 0.0000000000 } // 40-80
  ,{ 8.088291e-01, 2.013231e+00 } // 20-50
};

// MASS   0=>PIZERO, 1=>ETA, 2=>RHO0, 3=>OMEGA, 4=>ETAPRIME, 5=>PHI, 6=>JPSI, 7=>SIGMA, 8=>K0s, 9=>DELTA++, 10=>DELTA+, 11=>DELTA-, 12=>DELTA0, 13=>Rho+, 14=>Rho-, 15=>K0*, 16=>K0l, 17=>Lambda, 18=>K+, 19=>K-, 20=>Omega+, 21=>Omega-, 22=>Xi+, 23=>Xi-, 24=>Sigma+, 25=>Sigma-
const Double_t AliGenEMlibV2::fgkHM[26] = {0.1349766, 0.547853, 0.77549, 0.78265, 0.95778, 1.019455, 3.096916, 1.192642, 0.497614, 1.2311, 1.2349, 1.2349, 1.23340, 0.77549, 0.77549, 0.896, 0.497614, 1.115683, 0.493677, 0.493677, 1.67245, 1.67245, 1.32171, 1.32171, 1.3828, 1.3872};

const Double_t AliGenEMlibV2::fgkMtFactor[3][26] = {
  // {1.0, 0.5, 1.0, 0.9, 0.4, 0.23, 0.054},  // factor for pp from arXiv:1110.3929
  // {1.0, 0.55, 1.0, 0.9, 0.4, 0.25, 0.004}    // factor for PbPb from arXiv:1110.3929
  //{1., 0.48, 1.0, 0.9, 0.25, 0.4}, (old values)
  //{1., 0.48, 1.0, 0.9, 0.4, 0.25}, (nlo values)
  //{1., 0.48, 1.0, 0.8, 0.4, 0.2, 0.06} (combination of nlo and LHC measurements)
  //https://aliceinfo.cern.ch/Figure/node/2634
  //https://aliceinfo.cern.ch/Figure/node/2788
  //https://aliceinfo.cern.ch/Figure/node/4403
  //https://aliceinfo.cern.ch/Figure/node/5842
  //https://aliceinfo.cern.ch/Notes/node/87
  /*best guess:
   - pp values for eta/pi0 [arXiv:1205.5724], omega/pi0 [arXiv:1210.5749], phi/(pi+/-) [arXiv:1208.5717], K+-/pi+- [arXiv:1504.00024v2] from measured 7 Tev data
   */
  {1., 0.476, 1.0, 0.85, 0.4, 0.13, 1., 0.49, 0.575, 1, 1, 1, 1, 1.0, 1.0, 1.0, 0.575, 0.18, 0.41, 0.41, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, //pp
  {1., 0.476, 1.0, 0.85, 0.4, 0.25, 1., 0.49, 0.575, 1, 1, 1, 1, 1.0, 1.0, 1.0, 0.575, 0.18, 0.41, 0.41, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, //pPb
  {1., 0.476, 1.0, 0.85, 0.4, 0.25, 1., 0.49, 0.575, 1, 1, 1, 1, 1.0, 1.0, 1.0, 0.575, 0.18, 0.41, 0.41, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}  //PbPb
};

// Exponential
Double_t AliGenEMlibV2::PtExponential(const Double_t *px, const Double_t *c){
  const double &pt=px[0];
  Double_t invYield = c[0]*exp(-pt*c[1]);
  
  return invYield*(2*TMath::Pi()*pt);
}

// Hagedorn with additional Powerlaw
Double_t AliGenEMlibV2::PtModifiedHagedornPowerlaw(const Double_t *px, const Double_t *c){
  const double &pt=px[0];
  Double_t invYield = c[0]*pow(c[1]+pt*c[2],-c[3])*CrossOverLc(c[5],c[4],pt)+CrossOverRc(c[7],c[6],pt)*c[8]*pow(pt+0.001,-c[9]); //pt+0.001: prevent powerlaw from exploding for pt->0
  
  return invYield*(2*TMath::Pi()*pt+0.001); //+0.001: be sure to be > 0
}

// integral over krollwada with S=1^2*(1-mee^2/mh^2)^3 from mee=0 up to mee=mh
// approximation is perfect for mh>20MeV
Double_t AliGenEMlibV2::IntegratedKrollWada(const Double_t *mh, const Double_t *){
  if(*mh<0.002941) return 0;
  return 2*log(*mh/0.000511/exp(1.75))/411.11/TMath::Pi();
}

//--------------------------------------------------------------------------
//
//                             DirectRealGamma
//
//--------------------------------------------------------------------------
Double_t AliGenEMlibV2::PtPromptRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  const static Double_t promptGammaPtParam[10] = { 1.908746e-02, 3.326402e-01, 7.525743e-01, 5.251425e+00, 9.275261e+00, 1.855052e+01, 9.855216e+00, 1.867316e+01, 1.198770e-01, 5.407858e+00 };
  //{ 2.146541e-02, 5.540414e-01, 6.664706e-01, 5.739829e+00, 1.816496e+01, 2.561591e+01, 1.121881e+01, 3.569223e+01, 6.624561e-02, 5.234547e+00 };
  
  return PtModifiedHagedornPowerlaw(px,promptGammaPtParam)*GetTAA(fgSelectedCentrality)*1.15;
  //correction factor 1.15 comes from the global fit of all ALICE direct gamma measurements (see definition of fgkThermPtParam), showing that direct gamma is about 15% above the NLO expectation
}

Double_t AliGenEMlibV2::PtThermalRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return PtExponential(px,fgkThermPtParam[fgSelectedCentrality]);
}

Double_t AliGenEMlibV2::PtDirectRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return PtPromptRealGamma(px,px)+PtThermalRealGamma(px,px);
}

Int_t AliGenEMlibV2::IpDirectRealGamma(TRandom *)
{
  return 220000;
}

Double_t AliGenEMlibV2::YDirectRealGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return YFlat(*px);
}

Double_t AliGenEMlibV2::V2DirectRealGamma( const Double_t *px, const Double_t */*dummy*/ )
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
Double_t AliGenEMlibV2::PtPromptVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return IntegratedKrollWada(px,px)*PtPromptRealGamma(px,px);
}

Double_t AliGenEMlibV2::PtThermalVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return IntegratedKrollWada(px,px)*PtThermalRealGamma(px,px);
}

Double_t AliGenEMlibV2::PtDirectVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return IntegratedKrollWada(px,px)*PtDirectRealGamma(px,px);
}

Int_t AliGenEMlibV2::IpDirectVirtGamma(TRandom *)
{
  return 220001;
}

Double_t AliGenEMlibV2::YDirectVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return YFlat(*px);
}

Double_t AliGenEMlibV2::V2DirectVirtGamma( const Double_t *px, const Double_t */*dummy*/ )
{
  return V2DirectRealGamma(px,px);
}


//--------------------------------------------------------------------------
//
//                              Pizero
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpPizero(TRandom *)
{
  // Return pizero pdg code
  return 111;
}

Double_t AliGenEMlibV2::PtPizero( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kPizero]->Eval(pt);
}

Double_t AliGenEMlibV2::YPizero( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2Pizero( const Double_t *px, const Double_t */*dummy*/ )
{
  double n1,n2,n3,n4,n5;
  double v1,v2,v3,v4,v5;
  switch(fgSelectedCollisionsSystem|fgSelectedCentrality) {
    case kPbPb|k0010:
      n1=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0005]);
      v1=V2Param(px,fgkV2param[k0005]);
      n2=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0510]);
      v2=V2Param(px,fgkV2param[k0510]);
      return (n1*v1+n2*v2)/(n1+n2);
      break;
    case kPbPb|k0020:
      n1=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0005]);
      v1=V2Param(px,fgkV2param[k0005]);
      n2=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k0510]);
      v2=V2Param(px,fgkV2param[k0510]);
      n3=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k1020]);
      v3=V2Param(px,fgkV2param[k1020]);
      // raw yeilds are not normalized per event
      return (n1*v1+n2*v2+n3*v3)/(n1+n2+n3);
      break;
    case kPbPb|k2040:
      n1=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k2030]);
      v1=V2Param(px,fgkV2param[k2030]);
      n2=PtModifiedHagedornPowerlaw(px,fgkRawPtOfV2Param[k3040]);
      v2=V2Param(px,fgkV2param[k3040]);
      return (n1*v1+n2*v2)/(n1+n2);
      break;
    case kPbPb|k0040:
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
      // raw yeilds are not normalized per event
      return (n1*v1+n2*v2+n3*v3+n4*v4+n5*v5)/(n1+n2+n3+n4+n5);
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
Int_t AliGenEMlibV2::IpEta(TRandom *)
{
  // Return eta pdg code
  return 221;
}

Double_t AliGenEMlibV2::PtEta( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kEta]->Eval(pt);
}

Double_t AliGenEMlibV2::YEta( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2Eta( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kEta); //V2Param(px,fgkV2param[1][fgSelectedV2Param]);
}


//--------------------------------------------------------------------------
//
//                              Rho
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpRho0(TRandom *)
{
  // Return rho pdg code
  return 113;
}

Double_t AliGenEMlibV2::PtRho0( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kRho0]->Eval(pt);
}

Double_t AliGenEMlibV2::YRho0( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2Rho0( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kRho0);
}


//--------------------------------------------------------------------------
//
//                              Omega
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpOmega(TRandom *)
{
  // Return omega pdg code
  return 223;
}

Double_t AliGenEMlibV2::PtOmega( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kOmega]->Eval(pt);
}

Double_t AliGenEMlibV2::YOmega( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2Omega( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kOmega);
  
}


//--------------------------------------------------------------------------
//
//                              Etaprime
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpEtaprime(TRandom *)
{
  // Return etaprime pdg code
  return 331;
}

Double_t AliGenEMlibV2::PtEtaprime( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kEtaprime]->Eval(pt);
}

Double_t AliGenEMlibV2::YEtaprime( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2Etaprime( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kEtaprime);
}


//--------------------------------------------------------------------------
//
//                              Phi
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpPhi(TRandom *)
{
  // Return phi pdg code
  return 333;
}

Double_t AliGenEMlibV2::PtPhi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kPhi]->Eval(pt);
}

Double_t AliGenEMlibV2::YPhi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2Phi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kPhi);
}


//--------------------------------------------------------------------------
//
//                              Jpsi
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpJpsi(TRandom *)
{
  // Return jpsi pdg code
  return 443;
}

Double_t AliGenEMlibV2::PtJpsi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kJpsi]->Eval(pt);
}

Double_t AliGenEMlibV2::YJpsi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2Jpsi( const Double_t *px, const Double_t */*dummy*/ )
{
  const static Double_t v2Param[16] = { 1.156000e-01, 8.936854e-01, 0.000000e+00, 4.000000e+00, 6.222375e+00, -1.600314e-01, 8.766676e-01, 7.824143e+00, 1.156000e-01, 3.484503e-02, 4.413685e-01, 0, 1, 3.484503e-02, 4.413685e-01, 7.2 };
  switch(fgSelectedCollisionsSystem|fgSelectedCentrality){
    case kPbPb|k2040: return V2Param(px,v2Param); break;
    case kPbPb|k0010: return 0.43*V2Param(px,v2Param); break;  //V2Pizero(0010)/V2Pizero(2040)=0.43 +-0.025
    case kPbPb|k1020: return 0.75*V2Param(px,v2Param); break;  //V2Pizero(1020)/V2Pizero(2040)=0.75 +-0.04
    case kPbPb|k0020: return 0.66*V2Param(px,v2Param); break;  //V2Pizero(0020)/V2Pizero(2040)=0.66 +-0.035
    case kPbPb|k0040: return 0.82*V2Param(px,v2Param); break;  //V2Pizero(0040)/V2Pizero(2040)=0.82 +-0.05
    default:
      return KEtScal(*px,kJpsi);
  }
  return 0;
}


//--------------------------------------------------------------------------
//
//                              Sigma0
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpSigma0(TRandom *)
{
  // Return Sigma pdg code
  return 3212;
}

Double_t AliGenEMlibV2::PtSigma0( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kSigma0]->Eval(pt);
}

Double_t AliGenEMlibV2::YSigma0( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2Sigma0( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kSigma0,3);
}


//--------------------------------------------------------------------------
//
//                              K0short
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpK0short(TRandom *)
{
  // Return kzeroshort pdg code
  return 310;
}

Double_t AliGenEMlibV2::PtK0short( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kK0s]->Eval(pt);
}

Double_t AliGenEMlibV2::YK0short( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2K0short( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kK0s);
}


//--------------------------------------------------------------------------
//
//                              K0long
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpK0long(TRandom *)
{
  // Return kzerolong pdg code
  return 130;
}

Double_t AliGenEMlibV2::PtK0long( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kK0l]->Eval(pt);
}

Double_t AliGenEMlibV2::YK0long( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2K0long( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kK0l);
}


//--------------------------------------------------------------------------
//
//                              Lambda
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpLambda(TRandom *)
{
  // Return kzerolong pdg code
  return 3122;
}

Double_t AliGenEMlibV2::PtLambda( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kLambda]->Eval(pt);
}

Double_t AliGenEMlibV2::YLambda( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2Lambda( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kLambda);
}


//--------------------------------------------------------------------------
//
//                              Delta ++
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpDeltaPlPl(TRandom *)
{
  // Return Delta++ pdg code
  return 2224;
}

Double_t AliGenEMlibV2::PtDeltaPlPl( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kDeltaPlPl]->Eval(pt);
}

Double_t AliGenEMlibV2::YDeltaPlPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2DeltaPlPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaPlPl,3);
}


//--------------------------------------------------------------------------
//
//                              Delta +
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpDeltaPl(TRandom *)
{
  // Return Delta+ pdg code
  return 2214;
}

Double_t AliGenEMlibV2::PtDeltaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kDeltaPl]->Eval(pt);
}

Double_t AliGenEMlibV2::YDeltaPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2DeltaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaPl,3);
}


//--------------------------------------------------------------------------
//
//                              Delta -
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpDeltaMi(TRandom *)
{
  // Return Delta- pdg code
  return 1114;
}

Double_t AliGenEMlibV2::PtDeltaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kDeltaMi]->Eval(pt);
}

Double_t AliGenEMlibV2::YDeltaMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2DeltaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaMi,3);
}


//--------------------------------------------------------------------------
//
//                              Delta 0
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpDeltaZero(TRandom *)
{
  // Return Delta0 pdg code
  return 2114;
}

Double_t AliGenEMlibV2::PtDeltaZero( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kDeltaZero]->Eval(pt);
}

Double_t AliGenEMlibV2::YDeltaZero( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2DeltaZero( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaZero,3);
}


//--------------------------------------------------------------------------
//
//                              rho +
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpRhoPl(TRandom *)
{
  // Return rho+ pdg code
  return 213;
}

Double_t AliGenEMlibV2::PtRhoPl( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kRhoPl]->Eval(pt);
}

Double_t AliGenEMlibV2::YRhoPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2RhoPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kRhoPl);
}


//--------------------------------------------------------------------------
//
//                              rho -
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpRhoMi(TRandom *)
{
  // Return rho- pdg code
  return -213;
}

Double_t AliGenEMlibV2::PtRhoMi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kRhoMi]->Eval(pt);
}

Double_t AliGenEMlibV2::YRhoMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2RhoMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kRhoMi);
}


//--------------------------------------------------------------------------
//
//                             K0 *
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpK0star(TRandom *)
{
  // Return K0 * pdg code
  return 313;
}

Double_t AliGenEMlibV2::PtK0star( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kK0star]->Eval(pt);
}

Double_t AliGenEMlibV2::YK0star( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
  
}

Double_t AliGenEMlibV2::V2K0star( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kK0star);
}


//--------------------------------------------------------------------------
//
//                             K+
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpKPl(TRandom *)
{
  // Return K+ pdg code
  return 321;
}

Double_t AliGenEMlibV2::PtKPl( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kKPl]->Eval(pt);
}

Double_t AliGenEMlibV2::YKPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2KPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kKPl);
}


//--------------------------------------------------------------------------
//
//                             K-
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpKMi(TRandom *)
{
  // Return K- pdg code
  return -321;
}

Double_t AliGenEMlibV2::PtKMi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kKMi]->Eval(pt);
}

Double_t AliGenEMlibV2::YKMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2KMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kKMi);
}


//--------------------------------------------------------------------------
//
//                             Omega+
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpOmegaPl(TRandom *)
{
  // Return Omega+ pdg code
  return -3334;
}

Double_t AliGenEMlibV2::PtOmegaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kOmegaPl]->Eval(pt);
}

Double_t AliGenEMlibV2::YOmegaPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2OmegaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kOmegaPl);
}


//--------------------------------------------------------------------------
//
//                             Omega-
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpOmegaMi(TRandom *)
{
  // Return Omega- pdg code
  return 3334;
}

Double_t AliGenEMlibV2::PtOmegaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kOmegaMi]->Eval(pt);
}

Double_t AliGenEMlibV2::YOmegaMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2OmegaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kOmegaMi);
}


//--------------------------------------------------------------------------
//
//                             Xi+
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpXiPl(TRandom *)
{
  // Return Xi+ pdg code
  return -3312;
}

Double_t AliGenEMlibV2::PtXiPl( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kXiPl]->Eval(pt);
}

Double_t AliGenEMlibV2::YXiPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2XiPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kXiPl);
}


//--------------------------------------------------------------------------
//
//                             Xi-
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpXiMi(TRandom *)
{
  // Return Xi- pdg code
  return 3312;
}

Double_t AliGenEMlibV2::PtXiMi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kXiMi]->Eval(pt);
}

Double_t AliGenEMlibV2::YXiMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2XiMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kXiMi);
}


//--------------------------------------------------------------------------
//
//                             Simga(1385)+
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpSigmaPl(TRandom *)
{
  // Return Simga(1385)+ pdg code (called Simga*+ in http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
  return 3224;
}

Double_t AliGenEMlibV2::PtSigmaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kSigmaPl]->Eval(pt);
}

Double_t AliGenEMlibV2::YSigmaPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2SigmaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kSigmaPl);
}


//--------------------------------------------------------------------------
//
//                             Simga(1385)-
//
//--------------------------------------------------------------------------
Int_t AliGenEMlibV2::IpSigmaMi(TRandom *)
{
  // Return Simga(1385)- pdg code (called Simga*- in http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf)
  return 3114;
}

Double_t AliGenEMlibV2::PtSigmaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  const double &pt=px[0];
  return fPtParametrization[kSigmaMi]->Eval(pt);
}

Double_t AliGenEMlibV2::YSigmaMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlibV2::V2SigmaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kSigmaMi);
}


//--------------------------------------------------------------------------
//
//                    flat rapidity distribution
//
//--------------------------------------------------------------------------
Double_t AliGenEMlibV2::YFlat(Double_t /*y*/)
{
  Double_t dNdy = 1.;
  return dNdy;
}


//--------------------------------------------------------------------------
//
//                             Mt-scaling
//
//--------------------------------------------------------------------------
TF1* AliGenEMlibV2::MtScal(Int_t np, TString name, Bool_t isMeson)
{
  // function that calculates the pt distribution of a given particle np
  // by mt scaling the pi0 pt distribution for mesons and the proton pt
  // distribution for baryons

  Double_t xmin, xmax;
  Int_t nPar;
  TString formulaBaseScaled, scaledPt;
  
  // value meson/pi0 (baryon/p) at 5 GeV/c
  Double_t NormPt       = 5.;
  Double_t scaledNormPt, norm;
  
  if (!isMeson && fPtParametrizationProton) {
    // scale baryons from protons
    fPtParametrizationProton->GetRange(xmin, xmax);
    nPar                  = fPtParametrizationProton->GetNpar();
    formulaBaseScaled     = fPtParametrizationProton->GetExpFormula();
    scaledPt              = Form("(TMath::Sqrt(x*x + %.7f*%.7f - %.7f*%.7f))",fgkHM[np],fgkHM[np],0.9382720,0.9382720);
    scaledNormPt          = TMath::Sqrt(NormPt*NormPt + fgkHM[np]*fgkHM[np] - 0.9382720*0.9382720);
    norm                  = fMtFactorHisto->GetBinContent(np+1) * fPtParametrizationProton->Eval(NormPt) / fPtParametrizationProton->Eval(scaledNormPt);
  } else {
    // scale mesons from pi0 (also baryons if proton is not provided)
    fPtParametrization[0]->GetRange(xmin, xmax);
    nPar                  = fPtParametrization[0]->GetNpar();
    formulaBaseScaled     = fPtParametrization[0]->GetExpFormula();
    scaledPt              = Form("(TMath::Sqrt(x*x + %.7f*%.7f - %.7f*%.7f))",fgkHM[np],fgkHM[np],fgkHM[0],fgkHM[0]);
    scaledNormPt          = TMath::Sqrt(NormPt*NormPt + fgkHM[np]*fgkHM[np] - fgkHM[0]*fgkHM[0]);
    norm                  = fMtFactorHisto->GetBinContent(np+1) * fPtParametrization[0]->Eval(NormPt) / fPtParametrization[0]->Eval(scaledNormPt);
  }
  
  TString formulaBaseScaledTemp = "";
  TString sub1                  = "";
  TString sub2                  = "";
  TString sub3                  = "";
  for (Int_t i=0; i<formulaBaseScaled.Length(); i++) {
    if (i>0) sub1               = formulaBaseScaled(i-1, 1);
    else sub1                   = "";
    sub2                        = formulaBaseScaled(i, 1);
    if (i<formulaBaseScaled.Length()-1)
      sub3                      = formulaBaseScaled(i+1, 1);
    else sub3                   = "";
    
    if (sub2.CompareTo("x")!=0) {
      formulaBaseScaledTemp += sub2;
    } else if (sub2.CompareTo("x")==0) {
      if (i==0) {
        formulaBaseScaledTemp += scaledPt;
      } else if (sub1.CompareTo("e")!=0 && sub3.CompareTo("p")!=0) {
        formulaBaseScaledTemp += scaledPt;
      } else {
        formulaBaseScaledTemp += sub2;
      }
    } else {
      formulaBaseScaledTemp += sub2;
    }
  }
  formulaBaseScaled = formulaBaseScaledTemp;
  
  TF1* result = new TF1(name.Data(), Form("%.10f * (x/%s) * (%s)", norm, scaledPt.Data(), formulaBaseScaled.Data()), xmin, xmax);
  if (!isMeson && fPtParametrizationProton) {
    for (Int_t i=0; i<nPar; i++) {
      result->SetParameter(i, fPtParametrizationProton->GetParameter(i));
    }
  } else {
    for (Int_t i=0; i<nPar; i++) {
      result->SetParameter(i, fPtParametrization[0]->GetParameter(i));
    }
  }
  return result;
}


//--------------------------------------------------------------------------
//
//                             Et-scaling
//
//--------------------------------------------------------------------------
Double_t AliGenEMlibV2::KEtScal(Double_t pt, Int_t np, Int_t nq)
{
  Double_t scaledPt = sqrt(pow(2.0/nq*(sqrt(pt*pt+fgkHM[np]*fgkHM[np])-fgkHM[np])+fgkHM[0],2)-fgkHM[0]*fgkHM[0]);
  return V2Pizero(&scaledPt, (Double_t*) 0);
}


//--------------------------------------------------------------------------
//
//                                  V2
//
//--------------------------------------------------------------------------
Double_t AliGenEMlibV2::V2Param(const Double_t *px, const Double_t *par)
{
  // Very general parametrization of the v2
  const double &pt=px[0];
  double val=CrossOverLc(par[4],par[3],pt)*(2*par[0]/(1+TMath::Exp(par[1]*(par[2]-pt)))-par[0])+CrossOverRc(par[4],par[3],pt)*((par[8]-par[5])/(1+TMath::Exp(par[6]*(pt-par[7])))+par[5]);
  double sys=0;
  if(fgSelectedV2Systematic){
    double syspt=((pt>par[15])&(fgSelectedV2Systematic>0))?par[15]:pt;
    sys=fgSelectedV2Systematic*par[11+fgSelectedV2Systematic*2]*pow(syspt,par[12+fgSelectedV2Systematic*2]);
  }
  return std::max(val+sys,0.0);
}

Double_t AliGenEMlibV2::V2Flat(const Double_t */*px*/, const Double_t */*param*/)
{
  // Flat v2
  return 0.0;
}


//--------------------------------------------------------------------------
//
//                                  TAA
//
//--------------------------------------------------------------------------
Double_t AliGenEMlibV2::GetTAA(Int_t cent){
  const static Double_t taa[17] = { 1.0,    // pp
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
    1.207, // 40-80
    6.85};   // 20-50
  return taa[cent];
}


//--------------------------------------------------------------------------
//
//                        set pt parametrizations
//
//--------------------------------------------------------------------------
Bool_t AliGenEMlibV2::SetPtParametrizations(TString fileName, TString dirName) {
  
  // open parametrizations file
  TFile* fParametrizationFile = TFile::Open(fileName.Data());
  if (!fParametrizationFile) AliFatalClass(Form("File %s not found",fileName.Data()));
  TDirectory* fParametrizationDir = (TDirectory*)fParametrizationFile->Get(dirName.Data());
  if (!fParametrizationDir) AliFatalClass(Form("Directory %s not found",dirName.Data()));
  
  // check for pi0 parametrization
  TF1* fPtParametrizationTemp = (TF1*)fParametrizationDir->Get("111_pt");
  if (!fPtParametrizationTemp) AliFatalClass(Form("File %s doesn't contain pi0 parametrization",fileName.Data()));
  fPtParametrization[0] = new TF1(*fPtParametrizationTemp);
  fPtParametrization[0]->SetName("111_pt");

  // check for proton parametrization (base for baryon mt scaling)
  TF1* fPtParametrizationProtonTemp = (TF1*)fParametrizationDir->Get("2212_pt");
  if (!fPtParametrizationProtonTemp) {
    AliWarningClass(Form("File %s doesn't contain proton parametrization, scaling baryons from pi0.", fileName.Data()));
    fPtParametrizationProton = NULL;
  } else {
    fPtParametrizationProton = new TF1(*fPtParametrizationProtonTemp);
    fPtParametrizationProton->SetName("2212_pt");
  }
  
  AliGenEMlibV2 lib;
  TRandom* rndm;

  // get parametrizations from file
  for (Int_t i=1; i<26; i++) {
    Int_t ip = (Int_t)(lib.GetIp(i, ""))(rndm);
    fPtParametrizationTemp = (TF1*)fParametrizationDir->Get(Form("%d_pt", ip));
    if (fPtParametrizationTemp) {
      fPtParametrization[i] = new TF1(*fPtParametrizationTemp);
      fPtParametrization[i]->SetName(Form("%d_pt", ip));
    } else {
      if (i==7 || i==9 || i==10 || i==11 || i==12 || i==17 || (i>=20 && i<=25))
        fPtParametrization[i] = (TF1*)MtScal(i, Form("%d_pt_mtScaled", ip), 0);
      else
        fPtParametrization[i] = (TF1*)MtScal(i, Form("%d_pt_mtScaled", ip), 1);
    }
  }
  
  fParametrizationFile->Close();
  delete fParametrizationFile;
  
  return kTRUE;
}


//--------------------------------------------------------------------------
//
//                     return pt parametrization
//
//--------------------------------------------------------------------------
TF1* AliGenEMlibV2::GetPtParametrization(Int_t np) {
  if (np<26)
    return fPtParametrization[np];
  else if (np==26)
    return fPtParametrizationProton;
  else
    return NULL;
}


//--------------------------------------------------------------------------
//
//                     set mt scaling factor histo
//
//--------------------------------------------------------------------------
void AliGenEMlibV2::SetMtScalingFactors(TString fileName, TString dirName) {
  
  // set collision system
  Int_t selectedCol;
  switch (fgSelectedCollisionsSystem){
    case kpp900GeV:
      selectedCol=0;
      break;
    case kpp2760GeV:
      selectedCol=0;
      break;
    case kpp7TeV:
      selectedCol=0;
      break;
    case kpPb:
      selectedCol=1;
      break;
    case kPbPb:
      selectedCol=2;
      break;
    default:
      selectedCol=0;
      printf("<AliGenEMlibV2::SetMtScalingFactors> no collision system has been given\n");
  }
  
  // open file
  TFile*        fMtFactorFile = TFile::Open(fileName.Data());
  TDirectory*   fMtFactorDir  = (TDirectory*)fMtFactorFile->Get(dirName.Data());

  // set bin labels
  fMtFactorHisto = new TH1D("histoMtScaleFactor", "", 26, 0.5, 26.5);
  fMtFactorHisto->GetYaxis()->SetTitle("mt scaling factor");
  fMtFactorHisto->GetXaxis()->SetBinLabel(1,"111");
  fMtFactorHisto->GetXaxis()->SetBinLabel(2,"221");
  fMtFactorHisto->GetXaxis()->SetBinLabel(3,"113");
  fMtFactorHisto->GetXaxis()->SetBinLabel(4,"223");
  fMtFactorHisto->GetXaxis()->SetBinLabel(5,"331");
  fMtFactorHisto->GetXaxis()->SetBinLabel(6,"333");
  fMtFactorHisto->GetXaxis()->SetBinLabel(7,"443");
  fMtFactorHisto->GetXaxis()->SetBinLabel(8,"3212");
  fMtFactorHisto->GetXaxis()->SetBinLabel(9,"310");
  fMtFactorHisto->GetXaxis()->SetBinLabel(10,"2224");
  fMtFactorHisto->GetXaxis()->SetBinLabel(11,"2214");
  fMtFactorHisto->GetXaxis()->SetBinLabel(12,"1114");
  fMtFactorHisto->GetXaxis()->SetBinLabel(13,"2114");
  fMtFactorHisto->GetXaxis()->SetBinLabel(14,"213");
  fMtFactorHisto->GetXaxis()->SetBinLabel(15,"-213");
  fMtFactorHisto->GetXaxis()->SetBinLabel(16,"313");
  fMtFactorHisto->GetXaxis()->SetBinLabel(17,"130");
  fMtFactorHisto->GetXaxis()->SetBinLabel(18,"3122");
  fMtFactorHisto->GetXaxis()->SetBinLabel(19,"321");
  fMtFactorHisto->GetXaxis()->SetBinLabel(20,"-321");
  fMtFactorHisto->GetXaxis()->SetBinLabel(21,"-3334");
  fMtFactorHisto->GetXaxis()->SetBinLabel(22,"3334");
  fMtFactorHisto->GetXaxis()->SetBinLabel(23,"-3312");
  fMtFactorHisto->GetXaxis()->SetBinLabel(24,"3312");
  fMtFactorHisto->GetXaxis()->SetBinLabel(25,"3224");
  fMtFactorHisto->GetXaxis()->SetBinLabel(26,"3114");
  fMtFactorHisto->SetDirectory(0);

  // check for mt scaling factor histo
  TH1D*             fMtFactorHistoTemp = NULL;
  if (fMtFactorDir) fMtFactorHistoTemp = (TH1D*)fMtFactorDir->Get("histoMtScaleFactor");
  if (fMtFactorHistoTemp) {
    AliGenEMlibV2 lib;
    TRandom* rndm;
    for (Int_t i=0; i<26; i++) {
      Int_t ip = (Int_t)(lib.GetIp(i, ""))(rndm);
      Double_t factor = 0.;
      for (Int_t j=1; j<fMtFactorHistoTemp->GetNbinsX()+1; j++) {
        factor = 0.;
        TString tempLabel = Form("%s", fMtFactorHistoTemp->GetXaxis()->GetBinLabel(j));
        if (tempLabel.Atoi()==ip) {
          factor = fMtFactorHistoTemp->GetBinContent(j);
          break;
        }
      }
      if (factor>0) fMtFactorHisto->SetBinContent(i+1, factor);
      else          fMtFactorHisto->SetBinContent(i+1, fgkMtFactor[selectedCol][i]);
    }
  } else {
    for (Int_t i=1; i<27; i++)
      fMtFactorHisto->SetBinContent(i, fgkMtFactor[selectedCol][i-1]);
  }

  fMtFactorFile->Close();
  delete fMtFactorFile;
}


//--------------------------------------------------------------------------
//
//                 return mt scaling factor histo
//
//--------------------------------------------------------------------------
TH1D* AliGenEMlibV2::GetMtScalingFactors() {
  return fMtFactorHisto;
}


//--------------------------------------------------------------------------
//
//                        set pt-y distributions
//
//--------------------------------------------------------------------------
Bool_t AliGenEMlibV2::SetPtYDistributions(TString fileName, TString dirName) {

  // open parametrizations file
  TFile* fPtYDistributionFile = TFile::Open(fileName.Data());
  if (!fPtYDistributionFile) AliFatalClass(Form("File %s not found",fileName.Data()));
  TDirectory* fPtYDistributionDir = (TDirectory*)fPtYDistributionFile->Get(dirName.Data());
  if (!fPtYDistributionDir) AliFatalClass(Form("Directory %s not found",dirName.Data()));

  // check for pt-y parametrizations
  AliGenEMlibV2 lib;
  TRandom* rndm;
  TH2F* ptYTemp = NULL;
  for (Int_t i=0; i<26; i++) {
    Int_t ip = (Int_t)(lib.GetIp(i, ""))(rndm);
    ptYTemp = (TH2F*)fPtYDistributionDir->Get(Form("%d_pt_y", ip));
    if (ptYTemp) {
      fPtYDistribution[i] = new TH2F(*ptYTemp);
      fPtYDistribution[i]->SetName(Form("%d_pt_y", ip));
      fPtYDistribution[i]->SetDirectory(0);
    } else {
      fPtYDistribution[i] = NULL;
    }
  }

  if (!fPtYDistribution[0]) AliFatalClass(Form("File %s doesn't contain pi0 pt-y distribution",fileName.Data()));

  fPtYDistributionFile->Close();
  delete fPtYDistributionFile;

  return kTRUE;
}


//--------------------------------------------------------------------------
//
//                     return pt-y distribution
//
//--------------------------------------------------------------------------
TH2F* AliGenEMlibV2::GetPtYDistribution(Int_t np) {
  if (np<26 && fPtYDistribution[np])
    return fPtYDistribution[np];
  else
    return NULL;
}


//==========================================================================
//
//                     Set Getters
//
//==========================================================================
typedef Double_t (*GenFunc) (const Double_t*,  const Double_t*);
typedef Int_t (*GenFuncIp) (TRandom *);

GenFunc AliGenEMlibV2::GetPt(Int_t param, const char * tname) const
{
  // Return pointer to pT parameterisation
  GenFunc func=0;
  TString sname(tname);
  
  switch (param) {
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
    case kRho0:
      func=PtRho0;
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
    case kSigma0:
      func=PtSigma0;
      break;
    case kK0s:
      func=PtK0short;
      break;
    case kK0l:
      func=PtK0long;
      break;
    case kLambda:
      func=PtLambda;
      break;
    case kDeltaPlPl:
      func=PtDeltaPlPl;
      break;
    case kDeltaPl:
      func=PtDeltaPl;
      break;
    case kDeltaMi:
      func=PtDeltaMi;
      break;
    case kDeltaZero:
      func=PtDeltaZero;
      break;
    case kRhoPl:
      func=PtRhoPl;
      break;
    case kRhoMi:
      func=PtRhoMi;
      break;
    case kK0star:
      func=PtK0star;
      break;
    case kKPl:
      func=PtKPl;
      break;
    case kKMi:
      func=PtKMi;
      break;
    case kOmegaPl:
      func=PtOmegaPl;
      break;
    case kOmegaMi:
      func=PtOmegaMi;
      break;
    case kXiPl:
      func=PtXiPl;
      break;
    case kXiMi:
      func=PtXiMi;
      break;
    case kSigmaPl:
      func=PtSigmaPl;
      break;
    case kSigmaMi:
      func=PtSigmaMi;
      break;
    default:
      func=0;
      printf("<AliGenEMlibV2::GetPt> unknown parametrisation\n");
  }
  return func;
}

GenFunc AliGenEMlibV2::GetY(Int_t param, const char * tname) const
{
  // Return pointer to y- parameterisation
  GenFunc func=0;
  TString sname(tname);
  
  switch (param) {
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
    case kRho0:
      func=YRho0;
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
    case kSigma0:
      func=YSigma0;
      break;
    case kK0s:
      func=YK0short;
      break;
    case kK0l:
      func=YK0long;
      break;
    case kLambda:
      func=YLambda;
      break;
    case kDeltaPlPl:
      func=YDeltaPlPl;
      break;
    case kDeltaPl:
      func=YDeltaPl;
      break;
    case kDeltaMi:
      func=YDeltaMi;
      break;
    case kDeltaZero:
      func=YDeltaZero;
      break;
    case kRhoPl:
      func=YRhoPl;
      break;
    case kRhoMi:
      func=YRhoMi;
      break;
    case kK0star:
      func=YK0star;
      break;
    case kKPl:
      func=YKPl;
      break;
    case kKMi:
      func=YKMi;
      break;
    case kOmegaPl:
      func=YOmegaPl;
      break;
    case kOmegaMi:
      func=YOmegaMi;
      break;
    case kXiPl:
      func=YXiPl;
      break;
    case kXiMi:
      func=YXiMi;
      break;
    case kSigmaPl:
      func=YSigmaPl;
      break;
    case kSigmaMi:
      func=YSigmaMi;
      break;
    default:
      func=0;
      printf("<AliGenEMlibV2::GetY> unknown parametrisation\n");
  }
  return func;
}

GenFuncIp AliGenEMlibV2::GetIp(Int_t param, const char * tname) const
{
  // Return pointer to particle type parameterisation
  GenFuncIp func=0;
  TString sname(tname);
  
  switch (param) {
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
    case kRho0:
      func=IpRho0;
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
    case kSigma0:
      func=IpSigma0;
      break;
    case kK0s:
      func=IpK0short;
      break;
    case kK0l:
      func=IpK0long;
      break;
    case kLambda:
      func=IpLambda;
      break;
    case kDeltaPlPl:
      func=IpDeltaPlPl;
      break;
    case kDeltaPl:
      func=IpDeltaPl;
      break;
    case kDeltaMi:
      func=IpDeltaMi;
      break;
    case kDeltaZero:
      func=IpDeltaZero;
      break;
    case kRhoPl:
      func=IpRhoPl;
      break;
    case kRhoMi:
      func=IpRhoMi;
      break;
    case kK0star:
      func=IpK0star;
      break;
    case kKPl:
      func=IpKPl;
      break;
    case kKMi:
      func=IpKMi;
      break;
    case kOmegaPl:
      func=IpOmegaPl;
      break;
    case kOmegaMi:
      func=IpOmegaMi;
      break;
    case kXiPl:
      func=IpXiPl;
      break;
    case kXiMi:
      func=IpXiMi;
      break;
    case kSigmaPl:
      func=IpSigmaPl;
      break;
    case kSigmaMi:
      func=IpSigmaMi;
      break;
    default:
      func=0;
      printf("<AliGenEMlibV2::GetIp> unknown parametrisation\n");
  }
  return func;
}

GenFunc AliGenEMlibV2::GetV2(Int_t param, const char * tname) const
{
  // Return pointer to v2-parameterisation
  GenFunc func=0;
  TString sname(tname);
  
  switch (param) {
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
    case kRho0:
      func=V2Rho0;
      break;
    case kOmega:
      func=V2Omega;
      break;
    case kEtaprime:
      func=V2Etaprime;
      break;
    case kPhi:
      func=V2Phi;
      break;
    case kJpsi:
      func=V2Jpsi;
      break;
    case kSigma0:
      func=V2Sigma0;
      break;
    case kK0s:
      func=V2K0short;
      break;
    case kK0l:
      func=V2K0long;
      break;
    case kLambda:
      func=V2Lambda;
      break;
    case kDeltaPlPl:
      func=V2DeltaPlPl;
      break;
    case kDeltaPl:
      func=V2DeltaPl;
      break;
    case kDeltaMi:
      func=V2DeltaMi;
      break;
    case kDeltaZero:
      func=V2DeltaZero;
      break;
    case kRhoPl:
      func=V2RhoPl;
      break;
    case kRhoMi:
      func=V2RhoMi;
      break;
    case kK0star:
      func=V2K0star;
      break;
    case kKPl:
      func=V2KPl;
      break;
    case kKMi:
      func=V2KMi;
      break;
    case kOmegaPl:
      func=V2OmegaPl;
      break;
    case kOmegaMi:
      func=V2OmegaMi;
      break;
    case kXiPl:
      func=V2XiPl;
      break;
    case kXiMi:
      func=V2XiMi;
      break;
    case kSigmaPl:
      func=V2SigmaPl;
      break;
    case kSigmaMi:
      func=V2SigmaMi;
      break;
    default:
      func=0;
      printf("<AliGenEMlibV2::GetV2> unknown parametrisation\n");
  }
  return func;
}
