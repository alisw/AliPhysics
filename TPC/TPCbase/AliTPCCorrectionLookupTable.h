/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliTPCCorrectionLookupTable
/// \brief AliTPCCorrectionLookupTable class
///
/// \author Jens Wiechula

#include "AliTPCCorrection.h"
#include <TVectorD.h>
#include <TMatrixFfwd.h>
#include <THn.h>

class AliTPCCorrectionLookupTable : public AliTPCCorrection {

public:
  AliTPCCorrectionLookupTable();
  virtual ~AliTPCCorrectionLookupTable();

  virtual void GetCorrection(const Float_t x[],const Short_t roc,Float_t dx[]);
  virtual void GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]);

  void SetupDefaultLimits();
  void CreateLookupTable(AliTPCCorrection &tpcCorr, Float_t stepSize=5.);
  void CreateLookupTableSinglePhi(AliTPCCorrection &tpcCorr, Int_t iPhi, Float_t stepSize=5.);

  void CreateLookupTableFromResidualDistortion(THn &resDist);
  void CreateResidual(AliTPCCorrection *distortion, AliTPCCorrection* correction);

  void MergePhiTables(const char* files);

  void   SetFillCorrection(Bool_t fill) { fFillCorrection=fill;   }
  Bool_t GetFillCorrection() const      { return fFillCorrection; }
  void BuildExactInverse();

  Int_t GetNR()   const { return fNR;   }
  Int_t GetNPhi() const { return fNPhi; }
  Int_t GetNZ()   const { return fNZ;   }

  const TVectorD& GetLimitsR()   const { return fLimitsR; }
  const TVectorD& GetLimitsPhi() const { return fLimitsPhi; }
  const TVectorD& GetLimitsZ()   const { return fLimitsZ; }

  void SetCorrScaleFactor(Float_t    val) { fCorrScaleFactor = val; }
  Float_t    GetCorrScaleFactor() const { return fCorrScaleFactor; }

private:

  // sizes of lookup tables
  // TODO: Remove, since it will be stored in the TVectorD anyhow?
  Int_t     fNR;                   ///< number of rows (r) used for lookup table
  Int_t     fNPhi;                 ///< number of phi slices used for lookup table
  Int_t     fNZ;                   ///< number of columns (z) used for lookup table

  Float_t   fCorrScaleFactor;      ///< overall scaling factor for the correction

  Bool_t    fFillCorrection;       ///< whether to also fill the correction tables
  //
  TVectorD  fLimitsR;              ///< bin limits in row direction
  TVectorD  fLimitsPhi;            ///< bin limits in phi direction
  TVectorD  fLimitsZ;              ///< bin limits in z direction
  // for distortion
  /// Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDxDist;        //[fNPhi]
  /// Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDyDist;        //[fNPhi]
  /// Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDzDist;        //[fNPhi]

  // for correction
  /// Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDxCorr;        //[fNPhi]
  /// Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDyCorr;        //[fNPhi]
  /// Array to store electric field integral (int Er/Ez)
  TMatrixF **fLookUpDzCorr;        //[fNPhi]

  void InitTables();
  void InitTableArrays();
  void InitTablesPhiBin(Int_t iPhi);

  void ResetTables();
  void ResetLimits();

  void GetInterpolation(const Float_t x[],const Short_t roc,Float_t dx[],
                        TMatrixF **mR, TMatrixF **mPhi, TMatrixF **mZ);

  void CreateLookupTablePhiBin(AliTPCCorrection &tpcCorr, Int_t iPhi, Float_t stepSize);

  void FindClosestPosition(const Int_t binR, const Int_t binZ, const Int_t binPhi,
                           const Float_t xref[3], Float_t xret[3]);
  AliTPCCorrectionLookupTable(const AliTPCCorrectionLookupTable &corr);
  AliTPCCorrectionLookupTable& operator= (const AliTPCCorrectionLookupTable &corr);

  /// \cond CLASSIMP
  ClassDef(AliTPCCorrectionLookupTable,3);  // TPC corrections dumped into a lookup table
  /// \endcond
};


