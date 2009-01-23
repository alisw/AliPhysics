#ifndef ALITPCPOINTCORRECTION_H
#define ALITPCPOINTCORRECTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TNamed.h"
#include "TArrayD.h"
#include "TObjArray.h"
#include "TVectorD.h"


 
class AliTPCPointCorrection:public TNamed {
public:
  AliTPCPointCorrection(); 
  AliTPCPointCorrection(const Text_t *name, const Text_t *title);
  virtual ~AliTPCPointCorrection();
  //  
  TVectorD * GetParamOutR(Int_t sector) {return (TVectorD*)fParamsOutR.At(sector);}
  TVectorD * GetParamOutZ(Int_t sector) {return (TVectorD*)fParamsOutZ.At(sector);}
  //
  Double_t      GetDrOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector);
  Double_t      GetDzOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector);

  static Double_t      SGetDrOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector);
  static Double_t      SGetDzOut(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector);
  
  static AliTPCPointCorrection* Instance();
  void SetInstance(AliTPCPointCorrection*param){fgInstance = param;}
  //
  Double_t CorrectionOutR0(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector);
  Double_t CorrectionOutZ0(Bool_t isGlobal, Bool_t type, Double_t cx, Double_t cy, Double_t cz, Int_t sector);
  
public: 
  //
  // Correction out
  //
  TObjArray   fParamsOutR;       // Parameters  for radial distortion  - outer field cage
  TObjArray   fParamsOutZ;       // Parameters  for z      distortion  - outer field cage
  Int_t       fParamOutRVersion;  // version of the parameterization
  TObjArray   fErrorsOutR;       // Parameters  for radial distortion  - outer field cage
  TObjArray   fErrorsOutZ;       // Parameters  for z      distortion  - outer field cage
  Int_t       fParamOutZVersion;  // version of the parameterization
  //
private:
  AliTPCPointCorrection(const AliTPCPointCorrection&); 
  AliTPCPointCorrection& operator=(const AliTPCPointCorrection&); 
  static AliTPCPointCorrection*   fgInstance; //! Instance of this class (singleton implementation)
  ClassDef(AliTPCPointCorrection, 1); 
};

#endif


