/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////
// AliTPCCorrectionLookupTable class                                              //
// Authors: Jens Wiechula                                            //
////////////////////////////////////////////////////////////////////////////

#include "AliTPCCorrection.h"
#include <TVectorD.h>
#include <TMatrixFfwd.h>

class AliTPCCorrectionLookupTable : public AliTPCCorrection {

public:
  AliTPCCorrectionLookupTable();
  virtual ~AliTPCCorrectionLookupTable();

  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);
  virtual void GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]);

  void CreateLookupTable(AliTPCCorrection &tpcCorr, Float_t stepSize=5.);


private:

  // sizes of lookup tables
  // TODO: Remove, since it will be stored in the TVectorD anyhow?
  Int_t     fNR;                   // number of rows (r) used for lookup table
  Int_t     fNPhi;                 // number of phi slices used for lookup table
  Int_t     fNZ;                   // number of columns (z) used for lookup table
  //
  TVectorD  fLimitsR;              // bin limits in row direction
  TVectorD  fLimitsPhi;            // bin limits in phi direction
  TVectorD  fLimitsZ;              // bin limits in z direction
  // for distortion
  TMatrixF **fLookUpDxDist;        //[fNPhi] Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDyDist;        //[fNPhi] Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDzDist;        //[fNPhi] Array to store electric field integral (int Er/Ez)

  // for correction
  TMatrixF **fLookUpDxCorr;        //[fNPhi] Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDyCorr;        //[fNPhi] Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDzCorr;        //[fNPhi] Array to store electric field integral (int Er/Ez)

  void InitTables();
  void ResetTables();
  void SetupDefaultLimits();

  void GetInterpolation(const Float_t x[],const Short_t roc,Float_t dx[],
                        TMatrixF **mR, TMatrixF **mPhi, TMatrixF **mZ);

  AliTPCCorrectionLookupTable(const AliTPCCorrectionLookupTable &corr);
  AliTPCCorrectionLookupTable& operator= (const AliTPCCorrectionLookupTable &corr);
  ClassDef(AliTPCCorrectionLookupTable,1);  // TPC corrections dumped into a lookup table
};


