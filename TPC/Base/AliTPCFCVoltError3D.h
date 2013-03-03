#ifndef ALITPCFCVOLTERROR3D_H
#define ALITPCFCVOLTERROR3D_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
// AliTPCFCVoltError3D class                                              //
// Authors: Jim Thomas, Stefan Rossegger                                  //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"


class AliTPCFCVoltError3D : public AliTPCCorrection {
public:
  AliTPCFCVoltError3D();
  virtual ~AliTPCFCVoltError3D();

  // initialization and update functions
  virtual void Init();
  virtual void Update(const TTimeStamp &timeStamp);

  // common setters and getters for tangled ExB effect
  virtual void SetOmegaTauT1T2(Float_t omegaTau,Float_t t1,Float_t t2) {
    fT1=t1; fT2=t2;
    const Double_t wt0=t2*omegaTau;     fC0=1./(1.+wt0*wt0);
    const Double_t wt1=t1*omegaTau;     fC1=wt1/(1.+wt1*wt1);
  };
  void SetC0C1(Float_t c0,Float_t c1) {fC0=c0;fC1=c1;} // CAUTION: USE WITH CARE
  Float_t GetC0() const {return fC0;}
  Float_t GetC1() const {return fC1;}

  // setters and getters 

  // Set rod shift in Voltage equivalents (40V ~ 1mm)
  // rod numbers: 0-17 (IFC), 18-35 (OFC)
  // note: strips move accordingly
  void SetRodVoltShiftA(Int_t rod, Float_t voltOffset, Bool_t doInit=kTRUE) {fRodVoltShiftA[rod]=voltOffset; fInitLookUp=doInit;}
  void SetRodVoltShiftC(Int_t rod, Float_t voltOffset, Bool_t doInit=kTRUE) {fRodVoltShiftC[rod]=voltOffset; fInitLookUp=doInit;}
  Float_t GetRodVoltShiftA(Int_t i) const {return fRodVoltShiftA[i];}// 0-17: IFC, 18-35; OFC
  Float_t GetRodVoltShiftC(Int_t i) const {return fRodVoltShiftC[i];}// 0-17: IFC, 18-35; OFC

  // Set rotated clip (just at High Voltage RODs) in Voltage equivalents (40V ~ 1mm)
  // rod number: 0 (IFC), 1 (OFC)
  void SetRotatedClipVoltA(Int_t rod, Float_t voltOffset, Bool_t doInit=kTRUE) {fRotatedClipVoltA[rod]=voltOffset; fInitLookUp=doInit;}
  void SetRotatedClipVoltC(Int_t rod, Float_t voltOffset, Bool_t doInit=kTRUE) {fRotatedClipVoltC[rod]=voltOffset; fInitLookUp=doInit;}
  Float_t GetRotatedClipVoltA(Int_t i) const {return fRotatedClipVoltA[i];}// (0,1):(IFC,OFC)
  Float_t GetRotatedClipVoltC(Int_t i) const {return fRotatedClipVoltC[i];}// (0,1):(IFC,OFC)

  // Set rod shift in Voltage equivalents (40V ~ 1mm)
  // rod numbers: 0-17 (OFC)
  // note: strips DO NOT move, only the copper rods do ...
  void SetCopperRodShiftA(Int_t rod, Float_t voltOffset, Bool_t doInit=kTRUE) {fCopperRodShiftA[rod]=voltOffset; fInitLookUp=doInit;}
  void SetCopperRodShiftC(Int_t rod, Float_t voltOffset, Bool_t doInit=kTRUE) {fCopperRodShiftC[rod]=voltOffset; fInitLookUp=doInit;}
  Float_t GetCopperRodShiftA(Int_t i) const {return fCopperRodShiftA[i];}// 0-17: IFC, 18-35; OFC
  Float_t GetCopperRodShiftC(Int_t i) const {return fCopperRodShiftC[i];}// 0-17: IFC, 18-35; OFC


  void InitFCVoltError3D(); // Fill the lookup tables
  void ForceInitFCVoltError3D() { fInitLookUp=kFALSE; InitFCVoltError3D(); }

  virtual void Print(const Option_t* option="") const;



protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

private:

  AliTPCFCVoltError3D(const AliTPCFCVoltError3D &);               // not implemented
  AliTPCFCVoltError3D &operator=(const AliTPCFCVoltError3D &);    // not implemented

  Float_t fC0; // coefficient C0           (compare Jim Thomas's notes for definitions)
  Float_t fC1; // coefficient C1           (compare Jim Thomas's notes for definitions)
  Float_t fRodVoltShiftA[36];      // Rod (&strips) shift in Volt (40V~1mm) 
  Float_t fRodVoltShiftC[36];      // Rod (&strips) shift in Volt (40V~1mm) 
  Float_t fRotatedClipVoltA[2];    // rotated clips at HV rod
  Float_t fRotatedClipVoltC[2];    // rotated clips at HV rod
  Float_t fCopperRodShiftA[36];    // only Rod shift 
  Float_t fCopperRodShiftC[36];    // only Rod shift 

  Bool_t fInitLookUp;           // flag to check if the Look Up table was created (SUM)
  Bool_t fInitLookUpBasic[6];   // ! flag if the basic lookup was created (shifted Rod (IFC,OFC) or rotated clip (IFC,OFC))


  TMatrixF *fLookUpErOverEz[kNPhi];   // Array to store electric field integral (int Er/Ez)
  TMatrixF *fLookUpEphiOverEz[kNPhi]; // Array to store electric field integral (int Er/Ez)
  TMatrixF *fLookUpDeltaEz[kNPhi];    // Array to store electric field integral (int Er/Ez)

  // basic numbers for the poisson relaxation //can be set individually in each class
  enum {kRows   =257}; // grid size in r direction used in the poisson relaxation // ( 2**n + 1 ) eg. 65, 129, 257 etc.
  enum {kColumns=129}; // grid size in z direction used in the poisson relaxation // ( 2**m + 1 ) eg. 65, 129, 257 etc.
  enum {kPhiSlicesPerSector = 10 }; // number of points in phi slices
  enum {kPhiSlices = 1+kPhiSlicesPerSector*3 };      // number of points in phi for the basic lookup tables
  enum {kIterations=100}; // Number of iterations within the poisson relaxation 

  // ugly way to store "partial" look up tables
  // needed for the faster calculation of the final distortion table

  // for Rod and Strip shift
  TMatrixD *fLookUpBasic1ErOverEz[kPhiSlices];   // ! Array to store electric field integral (int Er/Ez)
  TMatrixD *fLookUpBasic1EphiOverEz[kPhiSlices]; // ! Array to store electric field integral (int Ephi/Ez)
  TMatrixD *fLookUpBasic1DeltaEz[kPhiSlices];    // ! Array to store electric field integral (int Ez)

  TMatrixD *fLookUpBasic2ErOverEz[kPhiSlices];   // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic2EphiOverEz[kPhiSlices]; // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic2DeltaEz[kPhiSlices];    // ! Array to store electric field integral 

  // for rotated clips
  TMatrixD *fLookUpBasic3ErOverEz[kPhiSlices];   // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic3EphiOverEz[kPhiSlices]; // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic3DeltaEz[kPhiSlices];    // ! Array to store electric field integral 

  TMatrixD *fLookUpBasic4ErOverEz[kPhiSlices];   // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic4EphiOverEz[kPhiSlices]; // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic4DeltaEz[kPhiSlices];    // ! Array to store electric field integral 

  // for (only rod) shift (copper rods)
  TMatrixD *fLookUpBasic5ErOverEz[kPhiSlices];   // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic5EphiOverEz[kPhiSlices]; // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic5DeltaEz[kPhiSlices];    // ! Array to store electric field integral 

  TMatrixD *fLookUpBasic6ErOverEz[kPhiSlices];   // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic6EphiOverEz[kPhiSlices]; // ! Array to store electric field integral 
  TMatrixD *fLookUpBasic6DeltaEz[kPhiSlices];    // ! Array to store electric field integral 


  ClassDef(AliTPCFCVoltError3D,3); //
};

#endif
