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

/*  
  
  Unlinearities fitting:

  Model for Outer field cage:
  Unlinearities at the edge aproximated using two exponential decays.
 
  dz = dz0(r,z) +dr(r,z)*tan(theta) 
  dy =          +dr(r,z)*tan(phi)

   
  
    
*/


#include "TLinearFitter.h"
#include "Riostream.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "AliLog.h"
#include "AliTPCROC.h"
#include "AliTPCPointCorrection.h"

AliTPCPointCorrection* AliTPCPointCorrection::fgInstance = 0;

ClassImp(AliTPCPointCorrection)

AliTPCPointCorrection::AliTPCPointCorrection():
  TNamed(),
  fParamsOutR(38),
  fParamsOutZ(38),
  fParamOutRVersion(0),
  fErrorsOutR(38),
  fErrorsOutZ(38),
  fParamOutZVersion(0)
{
  //
  // Default constructor
  //
}

AliTPCPointCorrection::AliTPCPointCorrection(const Text_t *name, const Text_t *title):
  TNamed(name,title),
  fParamsOutR(38),
  fParamsOutZ(38),
  fParamOutRVersion(0),
  fErrorsOutR(38),
  fErrorsOutZ(38),
  fParamOutZVersion(0)
{
  //
  // Non default constructor
  //
}

AliTPCPointCorrection::~AliTPCPointCorrection(){
  //
  //
  //
}


AliTPCPointCorrection* AliTPCPointCorrection::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliTPCPointCorrection();
  }
  return fgInstance;
}



Double_t AliTPCPointCorrection::GetDrOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  //  return radial correction
  //
  if (fParamOutRVersion==0) return CorrectionOutR0(isGlobal, type,cx,cy,cz,sector);
  return 0;
}

Double_t      AliTPCPointCorrection::SGetDrOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  // return radial correction - static function
  // 
  return fgInstance->GetDrOut(isGlobal, type,cx,cy,cz,sector);
}




Double_t AliTPCPointCorrection::GetDzOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  //
  //
  if (fParamOutZVersion==0) return CorrectionOutZ0(isGlobal, type,cx,cy,cz,sector);
  return 0;
}


Double_t      AliTPCPointCorrection::SGetDzOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  //
  //
  return fgInstance->GetDzOut(isGlobal, type,cx,cy,cz,sector);
}




Double_t AliTPCPointCorrection::CorrectionOutR0(Bool_t isGlobal, Bool_t type,  Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  // return dR correction - for correction version 0 
  // Parameters:
  // isGlobal   - is the x in global frame
  // type       - kTRUE - use sectors - kFALSE - only Side param
  // cx, cy,cz  - cluster position
  // sector     - sector number
  if (isGlobal){
    // recalculate sector if in global frame
    Double_t alpha    = TMath::ATan2(cy,cx);
    if (alpha<0)  alpha+=TMath::Pi()*2;
    sector = Int_t(18*alpha/(TMath::Pi()*2));
  }

  if (type==kFALSE) sector=36+(sector%36>=18);  // side level parameters
  // dout
  Double_t radius = TMath::Sqrt(cx*cx+cy*cy);  
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xout = roc->GetPadRowRadiiUp(roc->GetNRows(37)-1);
  Double_t dout = xout-radius;
  //drift
  Double_t dr   = 0.5 - TMath::Abs(cz/250.);
  //
  //
  TVectorD * vec = GetParamOutR(sector);
  if (!vec) return 0;
  Double_t eout10 = TMath::Exp(-dout/10.);
  Double_t eout5 = TMath::Exp(-dout/5.);		    
  Double_t eoutp  = (eout10+eout5)*0.5;
  Double_t eoutm  = (eout10-eout5)*0.5;
  //
  Double_t result=0;
  result+=(*vec)[1]*eoutp;
  result+=(*vec)[2]*eoutm;
  result+=(*vec)[3]*eoutp*dr;
  result+=(*vec)[4]*eoutm*dr;
  result+=(*vec)[5]*eoutp*dr*dr;
  result+=(*vec)[6]*eoutm*dr*dr;
  return result;

}

Double_t AliTPCPointCorrection::CorrectionOutZ0(Bool_t isGlobal, Bool_t type,  Double_t cx, Double_t cy, Double_t cz, Int_t sector){
  //
  // return dR correction - for correction version 0 
  // Parameters:
  // isGlobal   - is the x in global frame
  // type       - kTRUE - use sectors - kFALSE - only Side param
  // cx, cy,cz  - cluster position
  // sector     - sector number
  if (isGlobal){
    // recalculate sector if in global frame
    Double_t alpha    = TMath::ATan2(cy,cx);
    if (alpha<0)  alpha+=TMath::Pi()*2;
    sector = Int_t(18*alpha/(TMath::Pi()*2));
  }

  if (type==kFALSE) sector=36+(sector%36>=18);  // side level parameters
  // dout
  Double_t radius = TMath::Sqrt(cx*cx+cy*cy);  
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t xout = roc->GetPadRowRadiiUp(roc->GetNRows(37)-1);
  Double_t dout = xout-radius;
  //drift
  Double_t dr   = 0.5 - TMath::Abs(cz/250.);
  //
  //
  TVectorD * vec = GetParamOutR(sector);
  if (!vec) return 0;
  Double_t eout10 = TMath::Exp(-dout/10.);
  Double_t eout5 = TMath::Exp(-dout/5.);		    
  Double_t eoutp  = (eout10+eout5)*0.5;
  Double_t eoutm  = (eout10-eout5)*0.5;
  //
  Double_t result=0;
  result+=(*vec)[1]*eoutp;
  result+=(*vec)[2]*eoutm;
  result+=(*vec)[3]*eoutp*dr;
  result+=(*vec)[4]*eoutm*dr;
  result+=(*vec)[5]*eoutp*dr*dr;
  result+=(*vec)[6]*eoutm*dr*dr;
  return result;

}



