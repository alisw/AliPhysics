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
  
  Double_t RPhiCOGCorrection(Int_t sector, Int_t padrow, Float_t pad, Float_t cy, Float_t y, Float_t z, Float_t ky, Float_t qMax, Float_t threhsold);
  Double_t SRPhiCOGCorrection(Int_t sector, Int_t padrow, Float_t pad, Float_t cy, Float_t y, Float_t z, Float_t ky,  Float_t qMax, Float_t threhsold);
  //
  Double_t GetEdgeQ0(Int_t sector, Int_t padrow, Float_t y);
  static   Double_t SGetEdgeQ0(Int_t sector, Int_t padrow, Float_t y);
  //
  // IROC -OROC+Quadrant alignment
  //
  void     AddCorrectionSector(TObjArray & sideAPar, TObjArray &sideCPar, TObjArray & sideACov, TObjArray &sideCCov, Bool_t reset);
  Double_t GetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz, Int_t quadrant =-1); 
  static Double_t SGetCorrectionSector(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz, Int_t quadrant=-1); 
  //
  // Global alignment
  //
  Double_t GetCorrection(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz);
  static Double_t SGetCorrection(Int_t coord, Int_t sector, Double_t lx, Double_t ly, Double_t lz);
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
  //  Edge rfi
  // 
  //
  // Alignment part
  //
  Double_t fXIO;               // OROC-IROC boundary
  Double_t fXmiddle;           // center of the TPC sector local X
  Double_t fXquadrant;         // x quadrant
  //
  // IROC OROC alignment
  //
  TObjArray fArraySectorIntParam; // array of sector alignment parameters
  TObjArray fArraySectorIntCovar; // array of sector alignment covariances 
  //
  // Kalman filter for global alignment
  //
  TMatrixD  *fSectorParam;     // Kalman parameter   
  TMatrixD  *fSectorCovar;     // Kalman covariance  
  //
  //
  //
private:

  AliTPCPointCorrection(const AliTPCPointCorrection&); 
  AliTPCPointCorrection& operator=(const AliTPCPointCorrection&); 
  static AliTPCPointCorrection*   fgInstance; //! Instance of this class (singleton implementation)
  ClassDef(AliTPCPointCorrection, 3); 
};

#endif


