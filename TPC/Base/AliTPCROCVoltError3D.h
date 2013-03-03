#ifndef ALITPCROCVOLTERROR3D_H
#define ALITPCROCVOLTERROR3D_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
// AliTPCROCVoltError3D class                                             //
// Authors: Jim Thomas, Stefan Rossegger                                  //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"
#include "TH2F.h"


class AliTPCROCVoltError3D : public AliTPCCorrection {
public:
  AliTPCROCVoltError3D();
  virtual ~AliTPCROCVoltError3D();

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
  void SetROCData(TMatrixD * matrix);
  // setters and getters 
  void SetROCDataFileName(const char * fname);
  const Char_t* GetROCDataFileName() const {return fROCDataFileName.Data();}

  // flag to wheter or not include the z aligment in the dz calculation 
  // if FALSE, the dz offset is purely due to the electric field change
  void SetROCDisplacement(Bool_t flag) { 
    if (flag!=fROCdisplacement) { fROCdisplacement = flag; fInitLookUp=kFALSE; }
  }
  Bool_t GetROCDisplacement() const { return fROCdisplacement; }
  
  // flag on wheter to consider the difference in the electron arrival between IROC and OROC
  // due to the different position of the Anode wires
  void SetElectronArrivalCorrection(Bool_t flag) { 
    if (flag!=fElectronArrivalCorrection) { fElectronArrivalCorrection = flag; fInitLookUp=kFALSE; }
  }
  Bool_t GetElectronArrivalCorrection() const { return fElectronArrivalCorrection; }


  void InitROCVoltError3D(); // Fill the lookup tables
  void ForceInitROCVoltError3D() { fInitLookUp=kFALSE; InitROCVoltError3D(); }; 

  Float_t GetROCVoltOffset(Int_t side, Float_t r0, Float_t phi0);
  TH2F* CreateHistoOfZAlignment(Int_t side, Int_t nx=250, Int_t ny=250);

  virtual void Print(const Option_t* option="") const;
  TMatrixD *GetMatrix() const {return fdzDataLinFit;}  // Linear fits of dz survey points (each sector=72) (z0,slopeX,slopeY)         

protected:
  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);

private:

  AliTPCROCVoltError3D(const AliTPCROCVoltError3D &);               // not implemented
  AliTPCROCVoltError3D &operator=(const AliTPCROCVoltError3D &);    // not implemented

  Float_t fC0; // coefficient C0           (compare Jim Thomas's notes for definitions)
  Float_t fC1; // coefficient C1           (compare Jim Thomas's notes for definitions)

  Bool_t fROCdisplacement;      // flag on wheter to consider the ROC displacement 
                                // when calculating the z distortions
  Bool_t fElectronArrivalCorrection; // flag on wheter to consider the difference 
                                      // in the electron arrival between IROC and OROC
                                      // due to the different position of the Anode wires

  Bool_t fInitLookUp;           // flag to check it the Look Up table was created (SUM)

  TMatrixF *fLookUpErOverEz[kNPhi];   // Array to store electric field integral (int Er/Ez)
  TMatrixF *fLookUpEphiOverEz[kNPhi]; // Array to store electric field integral (int Er/Ez)
  TMatrixF *fLookUpDeltaEz[kNPhi];    // Array to store electric field integral (int Er/Ez)

  TString  fROCDataFileName;         // filename of the survey data containing the lin Fit values
  TMatrixD *fdzDataLinFit;  // Linear fits of dz survey points (each sector=72) (z0,slopeX,slopeY)         

  // basic numbers for the poisson relaxation //can be set individually in each class
  enum {kRows   =257}; // grid size in r direction used in the poisson relaxation // ( 2**n + 1 ) eg. 65, 129, 257 etc.
  enum {kColumns=129}; // grid size in z direction used in the poisson relaxation // ( 2**m + 1 ) eg. 65, 129, 257 etc.
  enum {kPhiSlicesPerSector=10};  // phi slices per sector
  enum {kPhiSlices = 18*kPhiSlicesPerSector };    // number of points in phi for the basic lookup tables
  enum {kIterations=100}; // Number of iterations within the poisson relaxation 

  ClassDef(AliTPCROCVoltError3D,2); 
};

#endif
