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

using std::cout;
using std::endl;

ClassImp(AliGenEMlib)

//Initializers for static members
Int_t AliGenEMlib::fgSelectedCollisionsSystem=AliGenEMlib::kpp7TeV; 
Int_t AliGenEMlib::fgSelectedPtParamPi0=AliGenEMlib::kPizeroParam; 
Int_t AliGenEMlib::fgSelectedPtParamEta=AliGenEMlib::kEtaParamRatiopp; 
Int_t AliGenEMlib::fgSelectedPtParamOmega=AliGenEMlib::kOmegaParampp; 
Int_t AliGenEMlib::fgSelectedPtParamPhi=AliGenEMlib::kPhiParampp;
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

// const Double_t AliGenEMlib::fgkV2param[kCentralities][16] = {
//   // charged pion                                                                                                                        cent, based on: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/FlowPAGQM2012talkIdentified
//   {  0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // pp no V2
//   ,{ 6.554571e-02, 1.436915e+00, 4.610598e-02, 2.554090e+00, 1.300948e+00, 2.970850e-02, 4.767877e+00, 4.228885e+00, 6.025959e-02, 1.570851e-03, 1.108941e+00, 0, 1, 1.715434e-03, 4.088070e-01, 25 }  // 0-5
//   ,{ 1.171348e-01, 1.333067e+00, 4.537086e-02, 3.046348e+00, 3.903416e+00, 4.407152e-02, 9.123846e-01, 4.834531e+00, 1.186227e-01, 2.259463e-03, 8.916458e-01, 0, 1, 2.300647e-03, 4.531172e-01, 25 }  // 5-10
//   ,{ 1.748434e-01, 1.285199e+00, 4.219881e-02, 4.018858e+00, 4.255082e+00, 7.955896e-03, 1.183264e-01,-9.329627e+00, 5.826570e-01, 3.368057e-03, 5.437668e-01, 0, 1, 3.178663e-03, 3.617552e-01, 25 }  // 10-20
//   ,{ 2.149526e-01, 1.408792e+00, 5.062101e-02, 3.206279e+00, 3.988517e+00, 3.724655e-02, 1.995791e-01,-1.571536e+01, 6.494227e+00, 4.957874e-03, 4.903140e-01, 0, 1, 4.214626e-03, 3.457922e-01, 25 }  // 20-30
//   ,{ 2.408942e-01, 1.477541e+00, 5.768983e-02, 3.333347e+00, 3.648508e+00,-2.044309e-02, 1.004145e-01,-2.386625e+01, 3.301913e+00, 5.666750e-03, 5.118686e-01, 0, 1, 4.626802e-03, 3.188974e-01, 25 }  // 30-40
//   ,{ 2.495109e-01, 1.543680e+00, 6.217835e-02, 3.518863e+00, 4.557145e+00, 6.014553e-02, 1.491814e-01,-5.443647e+00, 5.403300e-01, 6.217285e-03, 5.580218e-01, 0, 1, 4.620486e-03, 3.792879e-01, 25 }  // 40-50
//   ,{ 2.166399e-01, 1.931033e+00, 8.195390e-02, 2.226640e+00, 3.106649e+00, 1.058755e-01, 8.557791e-01, 4.006501e+00, 2.476449e-01, 6.720714e-03, 6.342966e-01, 0, 1, 4.449839e-03, 4.968750e-01, 25 }  // 50-60
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-10 
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 20-40
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 40-60
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 60-80
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-20 
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 0-40 
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 20-80
//   ,{ 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 2.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10 }  // 40-80
// };

const Double_t AliGenEMlib::fgkV2param[kCentralities][16] = {
  // charged pion                                                                                                                        cent, based on: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/FlowPAGQM2012talkIdentified
  {  0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.0000000000, 0, 1, 0.0000000000, 1.0000000000, 10  }  // pp no V2
  ,{ 6.516808e-02, 1.449877e+00, 4.705013e-02, 3.025555e+00, 3.391947e+00, 2.364841e-02, 2.303211e+00, 4.443068e+00, 5.989572e-02, 2.079546e-03, 8.365175e-01, 0, 1, 2.079546e-03, 8.365175e-01, 6.0 }  // 0-5
  ,{ 1.199628e-01, 1.288148e+00, 4.230109e-02, 3.071488e+00, 3.649381e+00, 4.643994e-02, 1.210442e+00, 5.055715e+00, 1.089290e-01, 2.561783e-03, 7.907385e-01, 0, 1, 2.561783e-03, 7.907385e-01, 7.0 }  // 5-10
  ,{ 1.672604e-01, 1.367064e+00, 4.704060e-02, 3.005428e+00, 3.816287e+00, 5.845636e-02, 8.660308e-01, 4.587477e+00, 1.819530e-01, 2.945173e-03, 7.256487e-01, 0, 1, 2.945173e-03, 7.256487e-01, 8.0 }  // 10-20
  ,{ 2.079093e-01, 1.480510e+00, 5.518112e-02, 2.845102e+00, 3.485697e+00, 3.363945e-02, 3.599748e-01, 3.207685e+00, 3.494331e-01, 4.398452e-03, 6.399492e-01, 0, 1, 4.398452e-03, 6.399492e-01, 16.0}  // 20-30
  ,{ 2.275701e-01, 1.622512e+00, 6.793718e-02, 2.742464e+00, 3.229573e+00, 3.316905e-02, 3.128005e-01, 2.465896e+00, 4.296590e-01, 4.837696e-03, 6.855293e-01, 0, 1, 4.837696e-03, 6.855293e-01, 17.0}  // 30-40
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

//based on: arXiv:1401.1250 [nucl-ex], https://twiki.cern.ch/twiki/bin/view/ALICE/PWGLFPAGSPECTRALowToHighPt
const Double_t AliGenEMlib::fgkPtParam[kCentralities][10] = {
  {  0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // pp no V2
  ,{ 1.240475e+02, 7.238972e-01, 4.297020e-01, 9.315793e+00, 3.304777e+00, 4.821807e+00, 3.729274e+00, 2.400991e+00, 2.487974e+01, 5.728662e+00 } // 0-5
  ,{ 7.912895e+01, 7.185643e-01, 3.995308e-01, 9.887192e+00, 3.191115e+00, 4.495417e+00, 2.724148e+00, 4.382795e+00, 3.024189e+01, 5.860432e+00 } // 5-10
  ,{ 6.956488e+01, 6.911720e-01, 4.646870e-01, 8.550547e+00, 3.172843e+00, 4.624597e+00, 3.679825e+00, 2.529246e+00, 2.454819e+01, 5.846622e+00 } // 10-20
  ,{ 4.851407e+02, 9.341151e-01, 4.716673e-01, 1.058090e+01, 4.681218e+00, 7.261284e+00, 3.883227e+00, 6.638627e+00, 1.562806e+01, 5.772127e+00 } // 20-30 (based on older unpublished)
  ,{ 3.157060e+01, 6.849451e-01, 4.868669e-01, 8.394558e+00, 3.539142e+00, 5.495280e+00, 4.102638e+00, 3.722991e+00, 1.638622e+01, 5.935963e+00 } // 30-40 (based on older unpublished)
  ,{ 1.857919e+01, 6.185989e-01, 5.878869e-01, 7.035064e+00, 2.892415e+00, 4.339383e+00, 3.549679e+00, 2.821061e+00, 1.529318e+01, 6.091388e+00 } // 40-50 (based on older unpublished)
  ,{ 1.069397e+01, 5.816587e-01, 6.542961e-01, 6.472858e+00, 2.643870e+00, 3.929020e+00, 3.339224e+00, 2.410371e+00, 9.606748e+00, 6.116685e+00 } // 50-60 (based on older unpublished)
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-10 
  ,{ 3.936495e+01, 7.043984e-01, 4.552076e-01, 9.142929e+00, 3.036362e+00, 4.248826e+00, 2.477718e+00, 3.911314e+00, 2.045442e+01, 5.946635e+00 } // 20-40
  ,{ 1.944300e+01, 6.102861e-01, 6.715858e-01, 6.491069e+00, 2.578329e+00, 3.671677e+00, 3.162338e+00, 2.317503e+00, 1.224814e+01, 6.089513e+00 } // 40-60
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 60-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-20 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 0-40 
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 20-80
  ,{ 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,-1.0000000000, 1.0000000000,-1.0000000000, 1.0000000000, 0.0000000000, 0.0000000000 } // 40-80
};

const Double_t AliGenEMlib::fgkThermPtParam[kCentralities][2] = {
  // might also be interesting: https://aliceinfo.cern.ch/Notes/node/226
  {  0.0000000000, 0.0000000000 } // pp no V2
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
};

const Double_t AliGenEMlib::fgkModTsallisParamPi0PbPb[kCentralities][7] = { 
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // pp
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 0-5
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 5-10
  {2.09114,  2.7482,  6.911,  51.1383, -10.6896, 1.30818, -1.59137 }, // 10-20
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 20-30
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 30-40
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 40-50
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 50-60
  {970.684,  0.46451, 5.52064, 21.7707, -4.2495, 4.62292, -3.47699 }, // 00-10
  {3.22534,  2.83717, 7.50505, 38.8015, -8.93242, 0.9842,  -0.821508 }, // 20-40
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 40-60
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 60-80
  {44.3972,  1.0644,  5.92254, 40.9254, -8.07759, 3.30333, -2.25078 }, // 00-20
  {38.187,  1.58985, 6.81705, 31.2526, -8.35251, 3.39482, -1.56172 }, // 00-40
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 20-80
  {89.2412,  0.150509, 4.97424, 112.986, 643.257, 3.41969, -2.46034 }, // 40-80 
};

const Double_t AliGenEMlib::fgkModTsallisParamPiChargedPbPb[kCentralities][7] = { 
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // pp
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 0-5
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 5-10
  {117.693,  1.20567, 6.57362, 29.2275, -7.71065, 4.50046, -1.91093 }, // 10-20
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 20-30
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 30-40
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 40-50
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 50-60
  {229.035,  1.11283, 6.53404, 28.9739, -8.15381, 5.05703, -2.32825 }, // 00-10
  {45.5686,  1.39268, 6.72519, 29.4513, -7.23335, 3.80987, -1.18599 }, // 20-40
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 40-60
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 60-80
  {171.812,  1.14849, 6.54742, 29.0473, -7.96679, 4.82915, -2.15967 }, // 00-20
  {108.083,  1.21125, 6.59362, 28.9651, -7.70072, 4.56583, -1.89318 }, // 00-40
  {0.,   0.,   0.,   0.,   0.,   0.,   0.   }, // 20-80
  {3.14057,  0.329224, 4.8235,  491.145, 186.041, 2.91138, -2.20281 }, // 40-80 
};

const Double_t AliGenEMlib::fgkParamSetPi07TeV[kNPi0Param][7] = { 
  {0.134977,  2.31335/(2*TMath::Pi()), 0.1433,   7.003,   0,   0,   0   }, //kPizeroParam
  {-0.0580977, 0.580804,     5.0964,   -669.913,  -160.2,  4.45906, -0.465125 }, //kPizeroParamlow
  {-7.09454,  0.218554,     5.2278,   77.9877,  -359.457, 4.67862, -0.996938 }, //kPizeroParamhigh
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPichargedParam
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPichargedParamlow
  {-0.000632204, 0.371249,     4.56778,  -111488,  -15573.4, 4.33064, -1.5506  }, //kPichargedParamhigh
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPizeroParamAlter
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPizeroParamAlterlow
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPizeroParamAlterhigh
};

const Double_t AliGenEMlib::fgkParamSetPi02760GeV[kNPi0Param][7] = { 
  {0.134977,  1.7/(2*TMath::Pi()),  0.135,   7.1,   0,   0,   0   }, //kPizeroParam
  {-1.4583,  0.188108,     5.34499,  546.328,  -2142.93, 4.55572, -1.82455 }, //kPizeroParamlow
  {-0.0648943, 0.231029,     4.39238,  -3705.4,  -35.9761, 4.87127, -2.00983 }, //kPizeroParamhigh
  {0.755554,  0.20772,     5.24167,  11.8279,  1492.53, 3.93863, -1.72404 }, //kPichargedParam
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPichargedParamlow
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPichargedParamhigh
  {0.0610774,  5.86289,     -7.83516,  1.29301,  2.62416, 0.,   0.   }, //kPizeroParamAlter
  {0.065338,  5.86705,     -9.05494,  1.38435,  3.58848, 0.,   0.   }, //kPizeroParamAlterlow
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPizeroParamAlterhigh
};

const Double_t AliGenEMlib::fgkParamSetPi0900GeV[kNPi0Param][7] = { 
  {0.134977,  1.5/(2*TMath::Pi()),  0.132,   7.8,   0,   0,   0   }, //kPizeroParam
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPizeroParamlow
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPizeroParamhigh
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPichargedParam
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPichargedParamlow
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPichargedParamhigh
  {0.0467285,  59.6495,     -212.865,  0.0584012,  2.83816, 0.,   0.   }, //kPizeroParamAlter
  {0.0468173,  66.9511,     -136.287,  0.0298099,  1.17147, 0.,   0.   }, //kPizeroParamAlterlow
  {0.,   0.,       0.,    0.,    0.,   0.,   0.   }, //kPizeroParamAlterhigh
};


// MASS   0=>PIZERO, 1=>ETA, 2=>RHO0, 3=>OMEGA, 4=>ETAPRIME, 5=>PHI, 6=>JPSI, 7=>SIGMA, 8=>K0s, 9=>DELTA++, 10=>DELTA+, 11=>DELTA-, 12=>DELTA0, 13=>Rho+, 14=>Rho-, 15=>K0*
const Double_t AliGenEMlib::fgkHM[16] = {0.1349766, 0.547853, 0.77549, 0.78265, 0.95778, 1.019455, 3.096916, 1.192642, 0.497614, 1.2311, 1.2349, 1.2349, 1.23340, 0.77549, 0.77549, 0.896};

const Double_t AliGenEMlib::fgkMtFactor[3][16] = { 
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
    - pp values for eta/pi0 [arXiv:1205.5724], omega/pi0 [arXiv:1210.5749], phi/(pi+/-) [arXiv:1208.5717] from measured 7 Tev data
  */
  {1., 0.476, 1.0, 0.85, 0.4, 0.13, 1., 0.49, 0.575, 1, 1, 1, 1, 1.0, 1.0, 1.0}, //pp
  {1., 0.476, 1.0, 0.85, 0.4, 0.25, 1., 0.49, 0.575, 1, 1, 1, 1, 1.0, 1.0, 1.0}, //pPb
  {1., 0.476, 1.0, 0.85, 0.4, 0.25, 1., 0.49, 0.575, 1, 1, 1, 1, 1.0, 1.0, 1.0}  //PbPb
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


Double_t AliGenEMlib::PtXQCD( Double_t pt,
                              Double_t a,
                              Double_t b,
                              Double_t c,
                              Double_t d,
                              Double_t e,
                              Double_t f)
{
  // QCD inspired function by Martin Wilde
  // DISCLAIMER: Please be careful outside of the measured pT range
  Double_t invYield = 0;
  if(pt>0.05){
    invYield = a*pow(pt,-1*(b+c/(pow(pt,d)+e*pow(pt,f))));
  } else invYield = 100;
  return invYield*(2*TMath::Pi()*pt);
}

Double_t AliGenEMlib::PtQCD(  Double_t pt,
                              Double_t a,
                              Double_t b,
                              Double_t c,
                              Double_t d,
                              Double_t e)
{
  // QCD inspired function by Martin Wilde
  // DISCLAIMER: Please be careful outside of the measured pT range
  Double_t invYield = 0;
  if(pt>0.05){
    invYield = a*pow(pt,-1*(b+c/(pow(pt,d)+e)));
  } else invYield = 100;
  return invYield*(2*TMath::Pi()*pt);
}

Double_t AliGenEMlib::PtModTsallis( Double_t pt,
                                    Double_t a,
                                    Double_t b,
                                    Double_t c,
                                    Double_t d,
                                    Double_t e,
                                    Double_t f,
                                    Double_t g,
                                    Double_t mass)
{

  Double_t invYield = 0;
  Double_t mt = sqrt(mass*mass + pt*pt);
  Double_t pt2 = pt*pt;
  if(pt>0.05){
    invYield = a*TMath::Power((1.+(mt-mass)/(b)),-c)*(d+e*pt+pt2)/(f+g*pt+pt2);
  } else invYield = 100;
  return invYield*(2*TMath::Pi()*pt);
}

Double_t AliGenEMlib::PtParticleRatiopp(Double_t pt,
                                        Double_t m1,
                                        Double_t m2,
                                        Double_t c1,
                                        Double_t c2,
                                        Double_t T1,
                                        Double_t T2,
                                        Double_t n)
{
 
  Double_t ratio = 0;
  if (PtTsallis (pt, m2, c2, T2, n)>0) ratio = PtTsallis (pt, m1, c1, T1, n)/ PtTsallis (pt, m2, c2, T2, n);
  return ratio;
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
 
  return invYield*(2*TMath::Pi()*pt+0.001); //+0.001: be sure to be > 0
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
  const static Double_t promptGammaPtParam[10] = { 1.908746e-02, 3.326402e-01, 7.525743e-01, 5.251425e+00, 9.275261e+00, 1.855052e+01, 9.855216e+00, 1.867316e+01, 1.198770e-01, 5.407858e+00 };
  //{ 2.146541e-02, 5.540414e-01, 6.664706e-01, 5.739829e+00, 1.816496e+01, 2.561591e+01, 1.121881e+01, 3.569223e+01, 6.624561e-02, 5.234547e+00 };
 
  return PtModifiedHagedornPowerlaw(px,promptGammaPtParam)*GetTAA(fgSelectedCentrality)*1.15;
  //correction factor 1.15 comes from the global fit of all ALICE direct gamma measurements (see definition of fgkThermPtParam), showing that direct gamma is about 15% above the NLO expectation
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
  return 220000;
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
  return 220001;
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
//   std::cout << "intitializing collision system: " << fgSelectedCollisionsSystem <<"\t" << kpp900GeV <<"\t" << kpp2760GeV <<"\t" << kpp7TeV <<"\t" << kpPb <<"\t" << kPbPb << std::endl;
//   std::cout << "centrality: " << fgSelectedCentrality <<"\t"<< kpp <<"\t"<<  k0005<<"\t"<< k0510<<"\t"<< k1020<<"\t"<< k2030<<"\t"<< std::endl
//                                                           << k3040<<"\t"<< k4050<<"\t"<< k5060<<"\t"<< k0010<<"\t"<< k2040<<"\t"<< std::endl
//                                                           << k4060<<"\t"<< k6080<<"\t"<< k0020<<"\t"<< k0040<<"\t"<< k2080<<"\t"<< std::endl
//                                                           << k4080<< std::endl;
//   std::cout << "parametrisation: " << fgSelectedPtParamPi0 << std::endl;

  
  Double_t kc=0.;
  Double_t kn=0.;
  Double_t kcT=0.;
  Double_t kT=0.;
  Double_t kp0=0.;
  Double_t kp1=0.;
  Double_t ka=0.;
  Double_t kb=0.;
  Double_t kd=0.;
  
  double n1,n2,n3,n4;
  int oldCent;
  
  switch(fgSelectedCollisionsSystem) {
    case kPbPb:
      switch (fgSelectedPtParamPi0){
        case kPichargedParamNew:
          // fit to pi charged, same data like in kPiOldChargedPbPb,
          // but tested and compared against newest (2014) neutral pi measurement
          switch (fgSelectedCentrality){
            case k0005:
            case k0510:
            case k1020:
            case k2030:
            case k3040:
            case k4050:
            case k5060:
            case k2040:
            case k4060:
              return PtModifiedHagedornPowerlaw(px,fgkPtParam[fgSelectedCentrality]);
              break;
            case k0010:
              n1=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0005]);
              n2=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0510]);
              return (n1+n2)/2;
              break;
            case k0020:
              n1=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0005]);
              n2=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0510]);
              n3=PtModifiedHagedornPowerlaw(px,fgkPtParam[k1020]);
              return (n1+n2+2*n3)/4;
              break;
            case k0040:
              n1=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0005]);
              n2=PtModifiedHagedornPowerlaw(px,fgkPtParam[k0510]);
              n3=PtModifiedHagedornPowerlaw(px,fgkPtParam[k1020]);
              n4=PtModifiedHagedornPowerlaw(px,fgkPtParam[k2040]);
              return (n1+n2+2*n3+4*n4)/8;
              break;
            default:
              return 0; 
          }

        case kPichargedParamOld:
          switch (fgSelectedCentrality){
      // fit to pi charged v1
      // charged pion from ToF, unidentified hadrons scaled with pion from TPC
      // for Pb-Pb @ 2.76 TeV
            case k0005:
              kc=1347.5; kp0=0.9393; kp1=2.254; kn=11.294; kcT=0.002537; kT=2.414;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT); 
              break;
            case k0510:
              kc=1256.1; kp0=0.9545; kp1=2.248; kn=11.291; kcT=0.002662; kT=2.326;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k2030:
              kc=7421.6; kp0=1.2059; kp1=1.520; kn=10.220; kcT=0.002150; kT=2.196;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k3040:
              kc=1183.2; kp0=1.0478; kp1=1.623; kn=9.8073; kcT=0.00198333; kT=2.073;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            // the following is what went into the Pb-Pb preliminary approval (0-10%)
            case k0010:
              kc=1296.0; kp0=0.968; kp1=2.567; kn=12.27; kcT=0.004219; kT=2.207;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k1020:
              kc=986.0; kp0=0.9752; kp1=2.376; kn=11.62; kcT=0.003116; kT=2.213;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k2040:
              kc=17337.0; kp0=1.337; kp1=1.507; kn=10.629; kcT=0.00184; kT=2.234;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k4050:
              kc=6220.0; kp0=1.322; kp1=1.224; kn=9.378; kcT=0.000595; kT=2.383;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k5060:
              kc=2319.0; kp0=1.267; kp1=1.188; kn=9.044; kcT=0.000437; kT=2.276;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k4060:
              kc=4724.0; kp0=1.319; kp1=1.195; kn=9.255; kcT=0.000511; kT=2.344;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k6080:
              kc=2842.0; kp0=1.465; kp1=0.8324; kn=8.167; kcT=0.0001049; kT=2.29;
              return PtModifiedHagedornThermal(*px,kc,kp0,kp1,kn,kcT,kT);
              break;
            case k0020:
              oldCent=fgSelectedCentrality;
              fgSelectedCentrality=k0010;
              n1=PtPizero(px,px);
              fgSelectedCentrality=k1020;
              n2=PtPizero(px,px);
              fgSelectedCentrality=oldCent;
              return (n1+n2)/2;
              break;
            case k0040:
              oldCent=fgSelectedCentrality;
              fgSelectedCentrality=k0010;
              n1=PtPizero(px,px);
              fgSelectedCentrality=k1020;
              n2=PtPizero(px,px);
              fgSelectedCentrality=k2040;
              n3=PtPizero(px,px);
              fgSelectedCentrality=oldCent;
              return (n1+n2+2*n3)/4;
            default:
              return 0; 
          }
          
        case kPichargedParam:
          switch (fgSelectedCentrality){
            case k0010:
            case k1020:
            case k2040:
            case k0020:
            case k0040:
            case k4080:
              return PtModTsallis( *px,
                                  fgkModTsallisParamPiChargedPbPb[fgSelectedCentrality][0],
                                  fgkModTsallisParamPiChargedPbPb[fgSelectedCentrality][1],
                                  fgkModTsallisParamPiChargedPbPb[fgSelectedCentrality][2],
                                  fgkModTsallisParamPiChargedPbPb[fgSelectedCentrality][3],
                                  fgkModTsallisParamPiChargedPbPb[fgSelectedCentrality][4],
                                  fgkModTsallisParamPiChargedPbPb[fgSelectedCentrality][5],
                                  fgkModTsallisParamPiChargedPbPb[fgSelectedCentrality][6],
                                  0.135);
              break;
            default:
              return 0;       
          }
          
        case kPizeroParam:
          switch (fgSelectedCentrality){
            case k0010:
            case k1020:
            case k2040:
            case k0020:
            case k0040:
            case k4080:
              return PtModTsallis( *px,
                                    fgkModTsallisParamPi0PbPb[fgSelectedCentrality][0],
                                    fgkModTsallisParamPi0PbPb[fgSelectedCentrality][1],
                                    fgkModTsallisParamPi0PbPb[fgSelectedCentrality][2],
                                    fgkModTsallisParamPi0PbPb[fgSelectedCentrality][3],
                                    fgkModTsallisParamPi0PbPb[fgSelectedCentrality][4],
                                    fgkModTsallisParamPi0PbPb[fgSelectedCentrality][5],
                                    fgkModTsallisParamPi0PbPb[fgSelectedCentrality][6],
                                    0.135);
                break;
            default:
              return 0;       
          }
        
        default:
          return 0;
      }
      
    case kpPb:
      // fit to charged pions for p-Pb @ 5.02TeV     
      switch (fgSelectedPtParamPi0){
        case kPichargedParam:
          //kc=235.5; ka=0.6903; kb=0.06864; kp0=2.289; kp1=0.5872; kd=0.6474; kn=7.842; 
          kc = 80.718314; ka = 0.510550; kb = 0.081444; kp0 = 3.415777; kp1 = 0.722887; kd = 0.820271; kn = 7.140332;
          return PtModifiedHagedornExp2(*px,kc,ka,kb,kp0,kp1,kd,kn);
          break;
        case kPichargedParamlow:
          kc = 62.548591; ka = 0.410163; kb = 0.111024; kp0 = 3.643849; kp1 = 0.734388; kd = 0.889554; kn = 6.741050;
          return PtModifiedHagedornExp2(*px,kc,ka,kb,kp0,kp1,kd,kn);
          break;
        case kPichargedParamhigh:
          kc = 105.785183; ka = 0.598377; kb = 0.055882; kp0 = 3.751795; kp1 = 0.679347; kd = 0.767044; kn = 7.516741;
          return PtModifiedHagedornExp2(*px,kc,ka,kb,kp0,kp1,kd,kn);
          break;  
        default:
          return 0;
  
      }
    case kpp7TeV:
      switch (fgSelectedPtParamPi0){
          // Tsallis fit to final pizero (PHOS+PCM) -> used for publication
          // for pp @ 7 TeV    
        case kPizeroParam: // fit to combined spectrum with stat errors only
          return PtTsallis(*px,fgkParamSetPi07TeV[kPizeroParam][0],fgkParamSetPi07TeV[kPizeroParam][1],fgkParamSetPi07TeV[kPizeroParam][2],fgkParamSetPi07TeV[kPizeroParam][3]);
          break;   
        case kPizeroParamlow:
          return PtModTsallis(    *px, 
                fgkParamSetPi07TeV[kPizeroParamlow][0],
                fgkParamSetPi07TeV[kPizeroParamlow][1],
                fgkParamSetPi07TeV[kPizeroParamlow][2],
                fgkParamSetPi07TeV[kPizeroParamlow][3],
                fgkParamSetPi07TeV[kPizeroParamlow][4],
                fgkParamSetPi07TeV[kPizeroParamlow][5],
                fgkParamSetPi07TeV[kPizeroParamlow][6],
                0.135);
          break;
        case kPizeroParamhigh:
          return PtModTsallis(    *px, 
                fgkParamSetPi07TeV[kPizeroParamhigh][0],
                fgkParamSetPi07TeV[kPizeroParamhigh][1],
                fgkParamSetPi07TeV[kPizeroParamhigh][2],
                fgkParamSetPi07TeV[kPizeroParamhigh][3],
                fgkParamSetPi07TeV[kPizeroParamhigh][4],
                fgkParamSetPi07TeV[kPizeroParamhigh][5],
                fgkParamSetPi07TeV[kPizeroParamhigh][6],
                0.135);
          break;
        case kPichargedParamhigh: 
          return PtModTsallis(    *px, 
                fgkParamSetPi07TeV[kPichargedParamhigh][0],
                fgkParamSetPi07TeV[kPichargedParamhigh][1],
                fgkParamSetPi07TeV[kPichargedParamhigh][2],
                fgkParamSetPi07TeV[kPichargedParamhigh][3],
                fgkParamSetPi07TeV[kPichargedParamhigh][4],
                fgkParamSetPi07TeV[kPichargedParamhigh][5],
                fgkParamSetPi07TeV[kPichargedParamhigh][6],
                0.135);
          break;
        
        default:
          return 0;
      }

      
    case kpp2760GeV:
      switch (fgSelectedPtParamPi0){
          // Tsallis fit to pizero: published pi0
          // for pp @ 2.76 TeV
        case kPizeroParam: //published fit parameters
          return PtTsallis(*px,fgkParamSetPi02760GeV[kPizeroParam][0],fgkParamSetPi02760GeV[kPizeroParam][1],fgkParamSetPi02760GeV[kPizeroParam][2],fgkParamSetPi02760GeV[kPizeroParam][3]);
          break;
        case kPizeroParamlow:
          return PtModTsallis(  *px, 
              fgkParamSetPi02760GeV[kPizeroParamlow][0], 
              fgkParamSetPi02760GeV[kPizeroParamlow][1], 
              fgkParamSetPi02760GeV[kPizeroParamlow][2], 
              fgkParamSetPi02760GeV[kPizeroParamlow][3], 
              fgkParamSetPi02760GeV[kPizeroParamlow][4], 
              fgkParamSetPi02760GeV[kPizeroParamlow][5], 
              fgkParamSetPi02760GeV[kPizeroParamlow][6], 
              0.135);
          break;
        case kPizeroParamhigh:
          return PtModTsallis(  *px, 
              fgkParamSetPi02760GeV[kPizeroParamhigh][0], 
              fgkParamSetPi02760GeV[kPizeroParamhigh][1], 
              fgkParamSetPi02760GeV[kPizeroParamhigh][2], 
              fgkParamSetPi02760GeV[kPizeroParamhigh][3], 
              fgkParamSetPi02760GeV[kPizeroParamhigh][4], 
              fgkParamSetPi02760GeV[kPizeroParamhigh][5], 
              fgkParamSetPi02760GeV[kPizeroParamhigh][6], 
              0.135);
          break;
        case kPichargedParam: 
          return PtModTsallis(    *px, 
                fgkParamSetPi02760GeV[kPichargedParam][0],
                fgkParamSetPi02760GeV[kPichargedParam][1],
                fgkParamSetPi02760GeV[kPichargedParam][2],
                fgkParamSetPi02760GeV[kPichargedParam][3],
                fgkParamSetPi02760GeV[kPichargedParam][4],
                fgkParamSetPi02760GeV[kPichargedParam][5],
                fgkParamSetPi02760GeV[kPichargedParam][6],
                0.135);
          break;
        case kPizeroParamAlter: 
          return PtQCD(  *px, 
            fgkParamSetPi02760GeV[kPizeroParamAlter][0], 
            fgkParamSetPi02760GeV[kPizeroParamAlter][1], 
            fgkParamSetPi02760GeV[kPizeroParamAlter][2], 
            fgkParamSetPi02760GeV[kPizeroParamAlter][3], 
            fgkParamSetPi02760GeV[kPizeroParamAlter][4]);
          break;
        case kPizeroParamAlterlow:
          return PtQCD(  *px, 
            fgkParamSetPi02760GeV[kPizeroParamAlter][0], 
            fgkParamSetPi02760GeV[kPizeroParamAlter][1], 
            fgkParamSetPi02760GeV[kPizeroParamAlter][2], 
            fgkParamSetPi02760GeV[kPizeroParamAlter][3],
            fgkParamSetPi02760GeV[kPizeroParamAlter][4]);
          break;
        default:
          return 0;   
      }  
    case kpp900GeV:
      switch (fgSelectedPtParamPi0){
          // Tsallis fit to pizero: published pi0
          // for pp @ 0.9 TeV
        case kPizeroParam: //published fit parameters
          return PtTsallis( *px,
          fgkParamSetPi0900GeV[kPizeroParam][0],
          fgkParamSetPi0900GeV[kPizeroParam][1],
          fgkParamSetPi0900GeV[kPizeroParam][2],
          fgkParamSetPi0900GeV[kPizeroParam][3]);
          break;
        case kPizeroParamAlter:
          return PtQCD(  *px, 
            fgkParamSetPi0900GeV[kPizeroParamAlter][0], 
            fgkParamSetPi0900GeV[kPizeroParamAlter][1], 
            fgkParamSetPi0900GeV[kPizeroParamAlter][2], 
            fgkParamSetPi0900GeV[kPizeroParamAlter][3], 
            fgkParamSetPi0900GeV[kPizeroParamAlter][4]);
          break;
        case kPizeroParamhigh:
          return PtQCD(  *px, 
            fgkParamSetPi0900GeV[kPizeroParamAlterlow][0], 
            fgkParamSetPi0900GeV[kPizeroParamAlterlow][1], 
            fgkParamSetPi0900GeV[kPizeroParamAlterlow][2], 
            fgkParamSetPi0900GeV[kPizeroParamAlterlow][3], 
            fgkParamSetPi0900GeV[kPizeroParamAlterlow][4]);
          break;
        default:
          return 0;   
      }  
    
    default:
      cout << "ERROR:: No valid collision system defined: "<< fgSelectedCollisionsSystem << endl;;
      return 0;
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
Int_t AliGenEMlib::IpEta(TRandom *)
{
  // Return eta pdg code
  return 221;     
}

Double_t AliGenEMlib::PtEta( const Double_t *px, const Double_t */*dummy*/ )
{

  // fit functions and corresponding parameter of Eta pT for pp @ 2.76 TeV and @ 7 TeV
  // and mtscaled pT 

  // parameters for Tsallis fit to eta
  Double_t km = 0.;
  Double_t kc = 0.;
  Double_t kT = 0.;
  Double_t kn = 0.;

  // parameters for Tsallis fit to pi0
  Double_t kmPi0 = 0.;
  Double_t kcPi0 = 0.;
  Double_t kTPi0 = 0.;
  Double_t knPi0 = 0.;

  // parameters for fit to eta/pi0
  Double_t krm1 = 0.;
  Double_t krm2 = 0.;
  Double_t krc1 = 0.;
  Double_t krc2 = 0.;
  Double_t krT1 = 0.;
  Double_t krT2 = 0.;
  Double_t krn = 0.;
 
  switch(fgSelectedCollisionsSystem){
    case kpp7TeV:
      switch(fgSelectedPtParamEta){ 
        // Tsallis fit to final eta (PHOS+PCM) -> used stat errors only for final publication
        // for pp @ 7 TeV
        case kEtaParamRatiopp:
          krm1 = 0.547853; krm2 = 0.134977; krc1 = 1.44198e+11; krc2 = 2.06751e+12 ; krT1 = 0.154567 ; krT2 = 0.139634 ; krn=32.0715; 
          kmPi0=0.134977; kcPi0=2.31335/(2*TMath::Pi()); kTPi0=0.1433; knPi0=7.003;
          return PtParticleRatiopp(*px, krm1, krm2, krc1, krc2, krT1, krT2, krn) * PtTsallis(*px,kmPi0,kcPi0,kTPi0,knPi0);
          break;
        case kEtaParampp:
          km = 0.547853; kc = 0.290164/(2*TMath::Pi()); kT = 0.212; kn = 7.352;
          return PtTsallis(*px,km,kc,kT,kn);
          break;
          // NOTE: None of these parametrisations look right - no idea where they come from 
        case kEtaParampplow:
          km = 0.547853; kc = 1.970; kT = 0.253; kn = 7.591;
          return PtTsallis(*px,km,kc,kT,kn);
          break;
        case kEtaParampphigh:
          km = 0.547853; kc = 3.060; kT = 0.212; kn = 6.578;
          return PtTsallis(*px,km,kc,kT,kn);
          break;
        case kEtaMtScal:
        default:
          return MtScal(*px,kEta);
      }
    case kpp2760GeV: 
      switch(fgSelectedPtParamEta){ 
        // Tsallis fit to preliminary eta (QM'11)
        // for pp @ 2.76 TeV
        // NOTE: None of these parametrisations look right - no idea where they come from
        case kEtaParampp:
          km = 0.547853; kc = 1.971; kT = 0.188; kn = 6.308;
          return PtTsallis(*px,km,kc,kT,kn);
        case kEtaParampplow:
          km = 0.547853; kc = 1.228; kT = 0.220; kn = 7.030;
          return PtTsallis(*px,km,kc,kT,kn);
          break;
        case kEtaParampphigh:
          km = 0.547853; kc = 2.802; kT = 0.164; kn = 5.815;
          return PtTsallis(*px,km,kc,kT,kn);
          break;
        case kEtaMtScal:
        default:
          return MtScal(*px,kEta);
          break;
      }
    default:
      return MtScal(*px,kEta);
      break; 
  }
}

Double_t AliGenEMlib::YEta( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlib::V2Eta( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kEta); //V2Param(px,fgkV2param[1][fgSelectedV2Param]);
}

//--------------------------------------------------------------------------
//
//                              Rho
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpRho0(TRandom *)
{
  // Return rho pdg code
  return 113;     
}

Double_t AliGenEMlib::PtRho0( const Double_t *px, const Double_t */*dummy*/ )
{
  // Rho pT
  return MtScal(*px,kRho0);
}

Double_t AliGenEMlib::YRho0( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlib::V2Rho0( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kRho0);
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
  // fit functions and corresponding parameter of Omega pT for preliminary pp @ 7 TeV data
  // and mtscaled pT 

  // parameters for Tsallis fit to omega
  Double_t km = 0.;
  Double_t kc = 0.;
  Double_t kT = 0.;
  Double_t kn = 0.;

  // parameters for Tsallis fit to pi0
  Double_t kmPi0 = 0.;
  Double_t kcPi0 = 0.;
  Double_t kTPi0 = 0.;
  Double_t knPi0 = 0.;

  // parameters for fit to omega/pi0
  Double_t krm1 = 0.;
  Double_t krm2 = 0.;
  Double_t krc1 = 0.;
  Double_t krc2 = 0.;
  Double_t krT1 = 0.;
  Double_t krT2 = 0.;
  Double_t krn = 0.;
 
  switch(fgSelectedCollisionsSystem){
    case kpp7TeV:
      switch(fgSelectedPtParamOmega){ 
        // Tsallis fit to final omega (PHOS) -> stat errors only, preliminary QM12
        // for pp @ 7 TeV
      case kOmegaParamRatiopp:
        krm1 = 0.78265; krm2 = 0.134977; krc1 = 21240028553.4600143433; krc2 = 168266377865.0805969238 ; krT1 = 0.21175 ; krT2 = 0.14328 ; krn=12.8831831756; 
        kmPi0=0.134977; kcPi0=2.31335/(2*TMath::Pi()); kTPi0=0.1433; knPi0=7.003;
        return PtParticleRatiopp(*px, krm1, krm2, krc1, krc2, krT1, krT2, krn) * PtTsallis(*px,kmPi0,kcPi0,kTPi0,knPi0);
        break;
      case kOmegaParampp:
        km = 0.78265; kc = 0.340051/(2*TMath::Pi()); kT = 0.206; kn = 6.31422;
        return PtTsallis(*px,km,kc,kT,kn);
        break;
      case kOmegaMtScal:
      default:
        return MtScal(*px,kOmega);
      }
    default:
      return MtScal(*px,kOmega);
      break; 
  }
 
  return MtScal(*px,kOmega);

}

Double_t AliGenEMlib::YOmega( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlib::V2Omega( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kOmega);

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
  return MtScal(*px,kEtaprime);
}

Double_t AliGenEMlib::YEtaprime( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2Etaprime( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kEtaprime);
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
  // fit functions and corresponding parameter of Phi pT for preliminary pp @ 7 TeV data
  // and PbPb collisions
  // and mtscaled pT 

  // parameters for Tsallis fit to phi
  Double_t km = 0.;
  Double_t kc = 0.;
  Double_t kT = 0.;
  Double_t kn = 0.;

 
  switch(fgSelectedCollisionsSystem){
    //   case kPbPb:
    //    switch(fgSelectedCentrality){
    //     // Tsallis fit to final phi->K+K- (TPC, ITS) -> stat+syst
    //     case k0010:
    //      switch(fgSelectedPtParamPhi){ 
    //       case kPhiParamPbPb:
    //        km = 0.78265; kc = 0.340051/(2*TMath::Pi()); kT = 0.206; kn = 6.31422;
    //        return PtTsallis(*px,km,kc,kT,kn);
    //        break;
    //     case kPhiMtScal:
    //       default:
    //        return MtScal(*px,kPhi);
    //      } 
    //     default:
    //      return MtScal(*px,kPhi);
    //    }
    case kpp7TeV:
      switch(fgSelectedPtParamPhi){ 
        // Tsallis fit to final phi->K+K- (TPC, ITS) -> stat+syst
        // for pp @ 7 TeV
      case kPhiParampp:
        km = 1.01946; kc = 0.0269578/(2*TMath::Pi()); kT = 0.2718119311; kn = 6.6755739295;
        return PtTsallis(*px,km,kc,kT,kn);
        break;
      case kPhiMtScal:
      default:
        return MtScal(*px,kPhi);
      }
    case kpPb: 
      switch(fgSelectedPtParamPhi){ 
        // for pPb @ 5.02 TeV
      case kPhiParamPPb:
        km = 1.01946; kc = 0.13484191317 / (2*TMath::Pi()); kT = 0.44560169252; kn = 12.78215005772;
        return PtTsallis(*px,km,kc,kT,kn);
        break;
      case kPhiParamPPblow:
        km = 1.01946; kc = 0.12420303835 / (2*TMath::Pi()); kT = 0.45097611367; kn = 13.18917228630;
        return PtTsallis(*px,km,kc,kT,kn);
        break;
      case kPhiParamPPbhigh:
        km = 1.01946; kc = 0.14548654894 / (2*TMath::Pi()); kT = 0.44101775738; kn = 12.45359262330; 
        return PtTsallis(*px,km,kc,kT,kn);
        break;    
      case kPhiMtScal:
      default:
        return MtScal(*px,kPhi);
      } 
    default:
      return MtScal(*px,kPhi);
      break; 
  }
  return MtScal(*px,kPhi);

}

Double_t AliGenEMlib::YPhi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);
}

Double_t AliGenEMlib::V2Phi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kPhi);
}

//--------------------------------------------------------------------------
//
//                              Jpsi
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpJpsi(TRandom *)
{
  // Return jpsi pdg code
  return 443;
}

Double_t AliGenEMlib::PtJpsi( const Double_t *px, const Double_t */*dummy*/ )
{
  // Jpsi pT
  // based on: //https://aliceinfo.cern.ch/Notes/node/242, https://aliceinfo.cern.ch/Figure/node/3457, http://arxiv.org/abs/1203.3641
  const static Double_t jpsiPtParam[2][3] = {
    {  9.686337e-03, 2.629441e-01, 4.552044e+00 }
    ,{ 3.403549e-03, 2.897061e-01, 3.644278e+00 }
  };
  const double pt=px[0]*2.28/2.613;
  switch(fgSelectedCollisionsSystem|fgSelectedCentrality) {
    case kPbPb|k0020: return 2.405*PtDoublePowerlaw(&pt,jpsiPtParam[0]); break;
    case kPbPb|k2040: return 2.405*PtDoublePowerlaw(&pt,jpsiPtParam[1]); break;
    case kPbPb|k0040: return 0.5*2.405*(PtDoublePowerlaw(&pt,jpsiPtParam[0])+PtDoublePowerlaw(&pt,jpsiPtParam[1])); break;
    default:
      return MtScal(*px,kJpsi);
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
//                              Sigma
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpSigma(TRandom *)
{
  // Return Sigma pdg code
  return 3212;     
}

Double_t AliGenEMlib::PtSigma( const Double_t *px, const Double_t */*dummy*/ )
{
  // Sigma pT
  return MtScal(*px,kSigma0);
}

Double_t AliGenEMlib::YSigma( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2Sigma0( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kSigma0,3);
}


//--------------------------------------------------------------------------
//
//                              K0short
//
//--------------------------------------------------------------------------

Int_t AliGenEMlib::IpK0short(TRandom *)
{
  // Return kzeroshort pdg code
  return 310;
}

Double_t AliGenEMlib::PtK0short( const Double_t *px, const Double_t */*dummy*/ )
{
  // K0short pT

  Double_t ka = 0; 
  Double_t kb = 0; 
  Double_t kc = 0; 
  Double_t kd = 0; 
  Double_t ke = 0; 
  Double_t kf = 0;
  
  switch (fgSelectedCentrality){
    case k0010:
      ka =9.21859; kb=5.71299; kc=-3.34251; kd=0.48796; ke=0.0192272; kf=3.82224;
      return PtXQCD( *px, ka, kb, kc, kd, ke, kf);
      break;
    case k1020:
      ka=6.2377; kb=5.6133; kc=-117.295; kd=3.51154; ke=36.3047; kf=0.456243;
      return PtXQCD( *px, ka, kb, kc, kd, ke, kf);
      break;
    case k0020:
      ka=7.7278; kb=5.6686; kc=-3.29259; kd=0.475403; ke=0.0223951; kf=3.69326;
      return PtXQCD( *px, ka, kb, kc, kd, ke, kf);
      break;
    case k2040:
      ka=3.38301; kb= 5.5323; kc=-96.078; kd=3.30782; ke=31.085; kf=0.466908;
      return PtXQCD( *px, ka, kb, kc, kd, ke, kf);
      break;
    case k0040:
      ka=5.55478; kb=5.61919; kc=-125.635; kd=3.5637; ke=38.9668; kf=0.47068;
      return PtXQCD( *px, ka, kb, kc, kd, ke, kf);
      break;
    case k4080:
      ka=0.731606; kb=5.49931; kc=-25.3106; kd=2.2439; ke=8.25063; kf= 0.289288;
      return PtXQCD( *px, ka, kb, kc, kd, ke, kf);
      break;
    default:
      return MtScal(*px,kK0s);
      break;     
  }
}
Double_t AliGenEMlib::YK0short( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2K0sshort( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kK0s);
}


//--------------------------------------------------------------------------
//
//                              Delta ++
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpDeltaPlPl(TRandom *)
{
  // Return Delta++ pdg code
  return 2224;     
}

Double_t AliGenEMlib::PtDeltaPlPl( const Double_t *px, const Double_t */*dummy*/ )
{
  // Delta++ pT
  return MtScal(*px,kDeltaPlPl);
}

Double_t AliGenEMlib::YDeltaPlPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2DeltaPlPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaPlPl,3);
}


//--------------------------------------------------------------------------
//
//                              Delta +
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpDeltaPl(TRandom *)
{
  // Return Delta+ pdg code
  return 2214;     
}

Double_t AliGenEMlib::PtDeltaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  // Delta+ pT
  return MtScal(*px,kDeltaPl);
}

Double_t AliGenEMlib::YDeltaPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2DeltaPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaPl,3);
}


//--------------------------------------------------------------------------
//
//                              Delta -
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpDeltaMi(TRandom *)
{
  // Return Delta- pdg code
  return 1114;     
}

Double_t AliGenEMlib::PtDeltaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  // Delta- pT
  return MtScal(*px,kDeltaMi);
}

Double_t AliGenEMlib::YDeltaMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2DeltaMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaMi,3);
}



//--------------------------------------------------------------------------
//
//                              Delta 0
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpDeltaZero(TRandom *)
{
  // Return Delta0 pdg code
  return 2114;     
}

Double_t AliGenEMlib::PtDeltaZero( const Double_t *px, const Double_t */*dummy*/ )
{
  // Delta0 pT
  return MtScal(*px,kDeltaZero);
}

Double_t AliGenEMlib::YDeltaZero( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2DeltaZero( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kDeltaZero,3);
}


//--------------------------------------------------------------------------
//
//                              rho +
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpRhoPl(TRandom *)
{
  // Return rho+ pdg code
  return 213;     
}

Double_t AliGenEMlib::PtRhoPl( const Double_t *px, const Double_t */*dummy*/ )
{
  // rho + pT
  return MtScal(*px,kRhoPl);
}

Double_t AliGenEMlib::YRhoPl( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2RhoPl( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kRhoPl);
}


//--------------------------------------------------------------------------
//
//                              rho -
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpRhoMi(TRandom *)
{
  // Return rho- pdg code
  return -213;     
}

Double_t AliGenEMlib::PtRhoMi( const Double_t *px, const Double_t */*dummy*/ )
{
  // rho- pT
  return MtScal(*px,kRhoMi);
}

Double_t AliGenEMlib::YRhoMi( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2RhoMi( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kRhoMi);
}


//--------------------------------------------------------------------------
//
//                             K0 *
//
//--------------------------------------------------------------------------
Int_t AliGenEMlib::IpK0star(TRandom *)
{
  // Return K0 * pdg code
  return 313;     
}

Double_t AliGenEMlib::PtK0star( const Double_t *px, const Double_t */*dummy*/ )
{
  // K0 * pT
  return MtScal(*px,kK0star);
}

Double_t AliGenEMlib::YK0star( const Double_t *py, const Double_t */*dummy*/ )
{
  return YFlat(*py);

}

Double_t AliGenEMlib::V2K0star( const Double_t *px, const Double_t */*dummy*/ )
{
  return KEtScal(*px,kK0star);
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
      printf("<AliGenEMlib::MtScal> no collision system has been given\n");
  }
 
  Double_t norm = fgkMtFactor[selectedCol][np] * (PtPizero(&NormPt, (Double_t*) 0) / PtPizero(&scaledNormPt, (Double_t*) 0));

  return norm*(pt/scaledPt)*scaledYield;
}

Double_t AliGenEMlib::KEtScal(Double_t pt, Int_t np, Int_t nq)
{
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
    double syspt=((pt>par[15])&(fgSelectedV2Systematic>0))?par[15]:pt;
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
      func= PtSigma;
      break;
    case kK0s:
      func= PtK0short;
      break;
    case kDeltaPlPl:
      func= PtDeltaPlPl;
      break;
    case kDeltaPl:
      func= PtDeltaPlPl;
      break;
    case kDeltaMi:
      func= PtDeltaMi;
      break;
    case kDeltaZero:
      func= PtDeltaZero;
      break;
    case kRhoPl:
      func= PtRhoPl;
      break;
    case kRhoMi:
      func= PtRhoMi;
      break;
    case kK0star:
      func= PtK0star;
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
      func=YSigma;
      break;   
    case kK0s:
      func=YK0short;
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
      func=IpSigma;
      break; 
    case kK0s:
      func=IpK0short;
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
      func=V2K0sshort;
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

    default:
      func=0;
      printf("<AliGenEMlib::GetV2> unknown parametrisation\n");
  }
  return func;
}

