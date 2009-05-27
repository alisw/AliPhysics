/*
  
 Class AliTPCtransformation:
 Should represent general non linear transformation. Currently tune for TPConly.
 To be used:
 1. Simulation-Digitization
 2. Reconstruction - AliTPCTransform
 3. Calibration/Alignment (KalmanFilter, Milipedde)
 4. Set of transformation to be stored/retrieved as OCDB entry 

 Base functionality:
 
 1. Double_t GetDeltaXYZ(Int_t coord, Int_t volID, Double_t param, Double_t x, Double_t y, Double_t z)
 Get correction - return the delta of coordinate coord dx or dy or dz for given volID for point at point (x,y,z)
 All coordinates are global

 2. The transformation should work only for given volIDs and detector IDs
    Currently Bitmask is used for filtering 



*/

#include <string.h>
#include "TRandom.h"
#include "TMath.h"
#include "TBits.h"
#include "TFormula.h"
#include "TF1.h"
#include "TLinearFitter.h"
#include "TFile.h"
#include "TObjString.h"

#include "TTreeStream.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
#include "AliTPCTransformation.h"

ClassImp(AliTPCTransformation)


AliTPCTransformation::GenFuncG    AliTPCTransformation::fgFormulas[10000];
TObjArray*  AliTPCTransformation::fgFormulasName = new TObjArray(10000);


void AliTPCTransformation::RegisterFormula(const char * name, GenFuncG formula){
  //
  // Add Formula to the list of formulas
  //
  Int_t last= fgFormulasName->GetEntries();
  fgFormulasName->AddLast(new TObjString(name));
  fgFormulas[last]=formula;
}

Int_t  AliTPCTransformation::BuildBasicFormulas(){
  //
  //build list of basic formulas
  //
  RegisterFormula("TPCscalingRPol",(GenFuncG)(AliTPCTransformation::TPCscalingRPol));
  RegisterFormula("TPCscalingZDr",(GenFuncG)(AliTPCTransformation::TPCscalingZDr));
  RegisterFormula("TPCscalingPhiLocal",(GenFuncG)(AliTPCTransformation::TPCscalingPhiLocal));
  return 0;
}

AliTPCTransformation::GenFuncG  AliTPCTransformation::FindFormula(const char * name){
  //
  // find formula - if registered
  //
  if (fgFormulasName->FindObject(name)==0) return 0;
  Int_t entries = fgFormulasName->GetEntries();
  for (Int_t i=0;i<entries;i++){
    if (strcmp(fgFormulasName->At(i)->GetName(), name)==0){
      return fgFormulas[i];
    }
  }
  return 0;
}

Double_t AliTPCTransformation::Eval(const char * name, const Double_t*x,const Double_t*par){
  //
  // Only for test purposes - very slow
  //
  GenFuncG fun = FindFormula(name);
  if (!fun) return 0;
  return fun(x,par);
}


AliTPCTransformation::AliTPCTransformation():
  TNamed(),
  fNameX(0),          // x formula name
  fNameY(0),          // y formula name
  fNameZ(0),          // z formula name 
  //
  fBitMask(0),        // bitmaps - transformation only for specified volID
  fCoordSystem(0),    // coord system of  output deltas 
  fParam(0),          // free parameter of transformation 
  fSigma(0),          // uncertainty of the parameter
  fFixedParam(0),     // fixed parameters of tranformation  
  //
  fInit(kFALSE),      // initialization flag - set to kTRUE if corresponding formulas found
  fFormulaX(0),       // x formula - pointer to the function
  fFormulaY(0),       // y formula - pointer to the function
  fFormulaZ(0)       // z formula - pointer to the function
  //
{
  //
  // default constructor
  //
}



AliTPCTransformation::AliTPCTransformation(const char *name, TBits *mask, const char *fx, const char *fy, const char *fz, Int_t coordSystem,Double_t param, Double_t sigma, TVectorD *fixedParams):
  TNamed(name,name),
  fNameX(0),       // x formula name
  fNameY(0),       // y formula name
  fNameZ(0),       // z formula name
  fBitMask(mask),   // bitmaps - transformation only for specified volID
  fCoordSystem(coordSystem), // coordinate system of output deltas
  fParam(param),          // free parameter of transformation 
  fSigma(sigma),
  fFixedParam(0),     // fixed parameters of tranformation  
  //
  fInit(kFALSE),      // initialization flag - set to kTRUE if corresponding formulas found
  fFormulaX(0),       // x formula - pointer to the function
  fFormulaY(0),       // y formula - pointer to the function
  fFormulaZ(0)       // z formula - pointer to the function
{
  //
  // non default constructor
  //
  if (fx) fNameX= new TString(fx);
  if (fy) fNameY= new TString(fy);
  if (fz) fNameZ= new TString(fz);
  if (fixedParams) fFixedParam = new TVectorD(*fixedParams);
  if (!fFixedParam) fFixedParam = new TVectorD(1);
  Init();
}



Bool_t AliTPCTransformation::Init(){
  //
  // associate formulas with pointer to the function
  //
  Bool_t isOK=kTRUE;
  if (fNameX) {
    fFormulaX=FindFormula(fNameX->Data());
    if (fFormulaX==0) isOK=kFALSE;
  }
  if (fNameY) {
    fFormulaY=FindFormula(fNameY->Data());
    if (fFormulaY==0) isOK=kFALSE;
  }
  if (fNameZ) {
    fFormulaZ=FindFormula(fNameZ->Data());
    if (!fFormulaZ) isOK=kFALSE;
  }
  return isOK;
}

TBits * AliTPCTransformation::BitsSide(Bool_t aside){
  //
  TBits * bits = new TBits(72);
  for (Int_t i=0; i<72;i++){
    if (i%36<18 && aside) (*bits)[i]=kTRUE;
    if (i%36<18 && (!aside)) (*bits)[i]=kFALSE;
    if (i%36>=18 && aside) (*bits)[i]=kFALSE;
    if (i%36>=18 && (!aside)) (*bits)[i]=kTRUE;
  }
  return bits;
}

TBits * AliTPCTransformation::BitsAll(){
  //
  //
  //
  TBits * bits = new TBits(72);
  for (Int_t i=0; i<72;i++){
    (*bits)[i]=kTRUE;
  }
  return bits;
}

Double_t AliTPCTransformation::GetDeltaXYZ(Int_t coord, Int_t volID, Double_t param, Double_t x, Double_t y, Double_t z){
  //
  //
  //
  if (fBitMask && (!(*fBitMask)[volID])) return 0;
  Double_t xyz[4]={x,y,z, param};
  if (fCoordSystem==0){
    // cartezian system
    if (coord==0 && fFormulaX) return fFormulaX(xyz,fFixedParam->GetMatrixArray()); 
    if (coord==1 && fFormulaY) return fFormulaY(xyz,fFixedParam->GetMatrixArray()); 
    if (coord==2 && fFormulaY) return fFormulaY(xyz,fFixedParam->GetMatrixArray());
  }
  if (fCoordSystem==1){  
    // cylindrical system
    if (coord==2) {
      if (fFormulaZ==0) return 0;
      return fFormulaZ(xyz,fFixedParam->GetMatrixArray());
    }
    Double_t rrphiz[3]={0,0,0};
    if (fFormulaX) rrphiz[0] = fFormulaX(xyz,fFixedParam->GetMatrixArray());
    if (fFormulaY) rrphiz[1] = fFormulaY(xyz,fFixedParam->GetMatrixArray());
    Double_t alpha = TMath::ATan2(y,x);
    Double_t ca    = TMath::Cos(alpha);
    Double_t sa    = TMath::Sin(alpha);
    if (coord==0) return ca*rrphiz[0]-sa*rrphiz[1];
    if (coord==1) return sa*rrphiz[0]+ca*rrphiz[1];
  }
  return 0;
}

Double_t  AliTPCTransformation::TPCscalingRPol(Double_t *xyz, Double_t * param){
  //
  // Scaling and shift of TPC radius
  // xyz[0..2] - global xyz of point 
  // xyz[3]    - scale parameter
  Double_t radius  = 0.5 - TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])/250.;
  Double_t driftM  = 0.5 - TMath::Abs(xyz[2]/250.);
  Double_t deltaR  = TMath::Power(radius,param[0])*TMath::Power(driftM,param[1]);
  return deltaR*xyz[3];
}


Double_t  AliTPCTransformation::TPCscalingZDr(Double_t *xyz, Double_t * param){
  //
  //
  // Scaling and shift of TPC radius
  // xyz[0..2] - global xyz of point 
  // xyz[3]    - scale parameter
  Double_t driftP  = TMath::Power(1. - TMath::Abs(xyz[2]/250.), param[0]);
  Double_t deltaZ  = (xyz[2]>0) ? -driftP : driftP;
  return deltaZ*xyz[3];
}

Double_t  AliTPCTransformation::TPCscalingPhiLocal(Double_t *xyz, Double_t * param){
  //
  //
  // Scaling if the local y -phi
  // xyz[0..2] - global xyz of point 
  // xyz[3]    - scale parameter
  Double_t alpha       = TMath::ATan2(xyz[1],xyz[0]);
  Double_t sector      = TMath::Nint(9*alpha/TMath::Pi()-0.5);
  Double_t localAlpha  = (alpha-(sector+0.5)*TMath::Pi()/9.);
  Double_t radius      = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1])/250.;
  Double_t deltaAlpha  = TMath::Power(9*localAlpha*radius/TMath::Pi(),param[0]);
  return deltaAlpha*xyz[3];
}


