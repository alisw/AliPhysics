#ifndef ALI_TPC_SPACECHARGE3D_DRIFTLINE_H
#define ALI_TPC_SPACECHARGE3D_DRIFTLINE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
 * See cxx source for full Copyright notice                                */

/// \class AliTPCSpaceCharge3DDriftLine
/// \brief This class provides correction-distortion following the drift line          
///  
/// Usage:
///
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date November 11, 2015///

#include "TVectorD.h"
#include "TFormula.h"
#include "AliTPCCorrection.h"
#include "TStopwatch.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliLog.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector.h"
#include "TVector3.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "AliSysInfo.h"
#include "AliTPCROC.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCPoissonSolver.h"
#include "AliTPCLookUpTable3DInterpolatorD.h"
#include "AliTPC3DCylindricalInterpolator.h"
#include "AliTPCLookUpTable3DInterpolatorDFull.h"
#include "AliTPC3DCylindricalInterpolatorFull.h"
#include "AliTPC3DCylindricalInterpolator.h"


class TCollection;

class TTimeStamp;

class TFormula;

class TH3F;

class TH3;

class TH2F;

class TH2;

class TF1;


class AliTPCSpaceCharge3DDriftLine : public AliTPCCorrection {
public:
  AliTPCSpaceCharge3DDriftLine();

  AliTPCSpaceCharge3DDriftLine(const char * name, const char *title);

  AliTPCSpaceCharge3DDriftLine(const char * name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice);

  AliTPCSpaceCharge3DDriftLine(const char * name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice, Int_t interpolationOrder,
                               Int_t irregularGridSize, Int_t strategy, Int_t rbfKernelType);

  virtual ~AliTPCSpaceCharge3DDriftLine();

  void InitSpaceCharge3DPoissonIntegralDz(Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration,
                                          Double_t stoppingConv);

  void InitSpaceCharge3DPoisson(Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration, Double_t stoppingConv);

  void ForceInitSpaceCharge3DPoissonIntegralDz(Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration,
                                               Double_t stoppingConv);

  void GetDistortionCyl(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetCorrectionCyl(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetCorrectionCylAC(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetCorrectionCylACIrregular(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetDistortion(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetCorrection(const Float_t x[], Short_t roc, Float_t dx[]);

  Double_t GetChargeCylAC(const Float_t x[], Short_t roc);

  Double_t GetInverseChargeCylAC(const Float_t x[], Short_t roc);

  void SetStrategyType(Int_t strategy) {
    fStrategy = strategy;
  }

  void SetCorrectionType(Int_t correctionType) {
    fCorrectionType = correctionType;
  }

  TH2F *CreateHistoDistDRinXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistoDistDRPhiinXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistoDistDZinXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistoCorrDRinXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistoCorrDRPhiinXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistoCorrDZinXY(Float_t z, Int_t nx, Int_t ny);

  enum {
    kNumSector = 18
  };

  enum StrategyType {
    kNaive = 0, ///< Naive, calculate each drift line point
    kUseInterpolator = 1,  ///< Use interpolation for
  };

  enum CorrectionType {
    kRegularInterpolator = 0,     ///< use interpolation with regular interpolator for correction look up table
    kIrregularInterpolator = 1,   ///< use irregular interpolator for correction look up table
  };


  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Double_t norm);

  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D) { SetInputSpaceCharge(hisSpaceCharge3D, 1); }

  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Double_t norm, Int_t side);

  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Int_t side) { SetInputSpaceCharge(hisSpaceCharge3D, 1, side); }


  void SetInputSpaceChargeA(TMatrixD **matricesLookUpCharge) {
    fInterpolatorChargeA->SetValue(matricesLookUpCharge);
    fInterpolatorChargeA->InitCubicSpline();
  }

  void SetInputSpaceChargeC(TMatrixD **matricesLookUpCharge) {
    fInterpolatorChargeC->SetValue(matricesLookUpCharge);
    fInterpolatorChargeC->InitCubicSpline();
  }


  TTree *CreateDistortionTree(Double_t step);

  TH2F *CreateHistoSCinXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistoSCinZR(Float_t phi, Int_t nz, Int_t nr);

  void SetNRRows(Int_t nRRow) { fNRRows = nRRow; }

  void SetNPhiSlices(Int_t nPhiSlice) { fNPhiSlices = nPhiSlice; }

  void SetNZColumns(Int_t nZColumn) { fNZColumns = nZColumn; }

  Int_t GetNRRows() { return fNRRows; }

  Int_t GetNPhiSlices() { return fNPhiSlices; }

  Int_t GetNZColumns() { return fNZColumns; }

  void SetPoissonSolver(AliTPCPoissonSolver *poissonSolver) { fPoissonSolver = poissonSolver; }

  AliTPCPoissonSolver *GetPoissonSolver() { return fPoissonSolver; }

  void SetInterpolationOrder(Int_t order) { fInterpolationOrder = order; }

  Int_t GetInterpolationOrder() { return fInterpolationOrder; }

  void SetOmegaTauT1T2(Float_t omegaTau, Float_t t1, Float_t t2) {
    fT1 = t1;
    fT2 = t2;
    const Double_t wt0 = t2 * omegaTau;
    fC0 = 1. / (1. + wt0 * wt0);
    const Double_t wt1 = t1 * omegaTau;
    fC1 = wt1 / (1. + wt1 * wt1);
  };

  void SetC0C1(Float_t c0, Float_t c1) {
    fC0 = c0;
    fC1 = c1;
  }

  Float_t GetC0() const { return fC0; }

  Float_t GetC1() const { return fC1; }

  void SetCorrectionFactor(Float_t correctionFactor) { fCorrectionFactor = correctionFactor; }

  Float_t GetCorrectionFactor() const { return fCorrectionFactor; }

  void InverseDistortionMaps(TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEphi,
                             TMatrixD **matricesEz, TMatrixD **matricesInvLocalIntErDz,
                             TMatrixD **matricesInvLocalIntEphiDz, TMatrixD **matricesInvLocalEz,
                             TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz, TMatrixD **matricesDistDz,
                             const Int_t nRRow, const Int_t nZColumn, const Int_t phiSliceconst, const Int_t nStep,
                             const Bool_t useCylAC, Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t interpType,
                             Int_t inverseType

  );

  void InverseDistortionMapsNoDrift(TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEphi,
                                    TMatrixD **matricesEz, TMatrixD **matricesInvLocalIntErDz,
                                    TMatrixD **matricesInvLocalIntEphiDz, TMatrixD **matricesInvLocalEz,
                                    TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz,
                                    TMatrixD **matricesDistDz, const Int_t nRRow, const Int_t nZColumn,
                                    const Int_t phiSlice);

  void GetCorrectionCylNoDrift(const Float_t x[], const Short_t roc, Float_t dx[]);

  void GetDistortionCylNoDrift(const Float_t x[], Short_t roc, Float_t dx[]);

  void InverseGlobalToLocalDistortionNoDrift(TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz,
                                             TMatrixD **matricesDistDz, Double_t *rList, Double_t *zList,
                                             Double_t *phiList, const Int_t nRRow, const Int_t nZColumn,
                                             const Int_t phiSlice);

  void GetChargeDensity(TMatrixD **matricesChargeA, TMatrixD **matricesChargeC, TH3 *spaceChargeHistogram3D,
                        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice);

  void GetInverseLocalDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetLocalDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]);

  void SetIrregularGridSize(Int_t size) { fIrregularGridSize = size; }

  Int_t GetIrregularGridSize() { return fIrregularGridSize; }

  Int_t GetRBFKernelType() { return fRBFKernelType; }

  void SetPotentialBoundaryAndCharge(TFormula *vTestFunction, TFormula *rhoTestFunction);

  void SetBoundaryIFCA(TF1 *f1) {
    fFormulaBoundaryIFCA = new TF1(*f1);
  }

  void SetBoundaryIFCC(TF1 *f1) {
    fFormulaBoundaryIFCC = new TF1(*f1);
  }

  void SetBoundaryOFCA(TF1 *f1) {
    fFormulaBoundaryOFCA = new TF1(*f1);
  }

  void SetBoundaryOFCC(TF1 *f1) {
    fFormulaBoundaryOFCC = new TF1(*f1);
  }

  void SetBoundaryROCA(TF1 *f1) {
    fFormulaBoundaryROCA = new TF1(*f1);
  }

  void SetBoundaryROCC(TF1 *f1) {
    fFormulaBoundaryROCC = new TF1(*f1);
  }

  void SetBoundaryCE(TF1 *f1) {
    fFormulaBoundaryCE = new TF1(*f1);
  }

  Float_t GetSpaceChargeDensity(Float_t r, Float_t phi, Float_t z);


  void Init();
private:
  static const Int_t kNMaxPhi = 360;

  Int_t fNRRows;     ///< the maximum on row-slices so far ~ 2cm slicing
  Int_t fNPhiSlices; ///< the maximum of phi-slices so far = (8 per sector)
  Int_t fNZColumns;  ///< the maximum on column-slices so  ~ 2cm slicing

  Float_t fC0; ///< coefficient C0 (compare Jim Thomas's notes for definitions)
  Float_t fC1; ///< coefficient C1 (compare Jim Thomas's notes for definitions)
  Float_t fCorrectionFactor; ///< Space Charge Correction factor in comparison to initialized

  Bool_t fInitLookUp; ///< flag to check if the Look Up table was created

  Double_t *fListR; //[fNRRows] list of r-coordinate of grids
  Double_t *fListPhi; //[fNPhiSlices] list of \f$ \phi\f$ -coordinate of grids
  Double_t *fListZ; //[fNZColumns]
  Double_t *fListZA; //[fNZColumns]  list of z-coordinate of grids
  Double_t *fListZC; //[fNZColumns] list of z-coordinate of grids
  Double_t *fListPotentialBoundaryA; //[fNRRows + fNNColumns] * 2 * fNPhiSlices
  Double_t *fListPotentialBoundaryC; //[fNRRows + fNNColumns] * 2 * fNPhiSlices

  Int_t fStrategy; ///>  Strategy for building  (0-> naive algorithm (O(n^4)), 1->use interpolation (O(n^3))
  Int_t fCorrectionType; ///> C
  Int_t fInterpolationOrder; ///>  Order of interpolation (1-> tri linear, 2->Lagrange interpolation order 2)
  Int_t fIrregularGridSize; ///>  Size of irregular grid cubes for interpolation (min 3)
  Int_t fRBFKernelType; ///>  RBF kernel type

  TMatrixD *fMatrixIntDistDrEz[kNMaxPhi];   //[kNPhi] Matrices for storing global distortion  \f$ R \f$ direction
  TMatrixD *fMatrixIntDistDphiREz[kNMaxPhi];   //[kNPhi] Matrices for storing global distortion  \f$ \phi\f$ direction
  TMatrixD *fMatrixIntDistDz[kNMaxPhi];  //[kNPhi] Matrices for storing Global distortion \f$ z \f$ Distortion

  TMatrixD *fMatrixIntDistDrEzA[kNMaxPhi];  //[kNPhi] Matrices for storing Global distortion  \f$ R \f$ direction for A Sise
  TMatrixD *fMatrixIntDistDphiREzA[kNMaxPhi]; //[kNPhi] Matrices for storing Global \f$ \phi R \f$ Distortion for A Side
  TMatrixD *fMatrixIntDistDzA[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ z \f$ Distortion for A Side

  TMatrixD *fMatrixIntDistDrEzC[kNMaxPhi]; //[kNPhi] Matrices for storing Global  \f$ R \f$ direction for C side
  TMatrixD *fMatrixIntDistDphiREzC[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ \phi R \f$ Distortion for C side
  TMatrixD *fMatrixIntDistDzC[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ z \f$ Distortion for C side

  TMatrixD *fMatrixErOverEz[kNMaxPhi]; //[kNPhi] Matrices for storing Er Over Ez for intermediate value
  TMatrixD *fMatrixEphiOverEz[kNMaxPhi]; //[kNPhi] Matrices for storing Er Over Ez for intermediate value
  TMatrixD *fMatrixDeltaEz[kNMaxPhi]; //[kNPhi] Matrices for storing Er Over Ez for intermediate value

  TMatrixD *fMatrixErOverEzA[kNMaxPhi]; //[kNPhi] Matrices for storing Er Over Ez for intermediate value for A side
  TMatrixD *fMatrixEphiOverEzA[kNMaxPhi];
  TMatrixD *fMatrixDeltaEzA[kNMaxPhi];

  TMatrixD *fMatrixErOverEzC[kNMaxPhi]; //[kNPhi] Matrices for storing Er Over Ez for intermediate value for C side
  TMatrixD *fMatrixEphiOverEzC[kNMaxPhi];
  TMatrixD *fMatrixDeltaEzC[kNMaxPhi];

  TMatrixD *fMatrixIntCorrDrEz[kNMaxPhi]; //[kNPhi] Matrices for storing Global  \f$  R \f$ correction
  TMatrixD *fMatrixIntCorrDphiREz[kNMaxPhi];   //[kNPhi] Matrices for storing Global  \f$ \phi R \f$  correction
  TMatrixD *fMatrixIntCorrDz[kNMaxPhi];  //[kNPhi]  Matrices for storing Global  \f$ X \f$ correctiona

  TMatrixD *fMatrixIntCorrDrEzA[kNMaxPhi]; //[kNPhi] Matrices for storing Global  \f$  R \f$ correction
  TMatrixD *fMatrixIntCorrDphiREzA[kNMaxPhi];   //[kNPhi] Matrices for storing Global  \f$ \phi R \f$  correction
  TMatrixD *fMatrixIntCorrDzA[kNMaxPhi]; //[kNPhi] Matrices for storing Global  \f$ X \f$ correctiona

  TMatrixD *fMatrixIntCorrDrEzC[kNMaxPhi]; //[kNPhi]  Matrices for storing Global  \f$  R \f$ correction
  TMatrixD *fMatrixIntCorrDphiREzC[kNMaxPhi];   //[kNPhi] Matrices for storing Global  \f$ \phi R \f$  correction
  TMatrixD *fMatrixIntCorrDzC[kNMaxPhi];  //[kNPhi] Matrices for storing Global  \f$ X \f$ correctiona

  // Correction with irregular interpolation side A
  TMatrixD *fMatrixIntCorrDrEzIrregularA[kNMaxPhi]; //[kNPhi] Matrices for storing Global  \f$ R \f$ direction
  TMatrixD *fMatrixIntCorrDphiREzIrregularA[kNMaxPhi];   //[kNPhi] Matrices for storing Globar \f$ \phi R \f$ Distortion
  TMatrixD *fMatrixIntCorrDzIrregularA[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ z \f$ Distortion

  TMatrixD *fMatrixRListIrregularA[kNMaxPhi]; //[kNPhi]
  TMatrixD *fMatrixPhiListIrregularA[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ \phi R \f$ Distortion
  TMatrixD *fMatrixZListIrregularA[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ z \f$ Distortion

  // Correction with irregular interpolation side C
  TMatrixD *fMatrixIntCorrDrEzIrregularC[kNMaxPhi]; //[kNPhi] Matrices for storing Global  \f$ R \f$ direction
  TMatrixD *fMatrixIntCorrDphiREzIrregularC[kNMaxPhi];   //[kNPhi] Matrices for storing Globar \f$ \phi R \f$ Distortion
  TMatrixD *fMatrixIntCorrDzIrregularC[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ z \f$ Distortion

  TMatrixD *fMatrixRListIrregularC[kNMaxPhi]; //[kNPhi]
  TMatrixD *fMatrixPhiListIrregularC[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ \phi R \f$ Distortion
  TMatrixD *fMatrixZListIrregularC[kNMaxPhi]; //[kNPhi] Matrices for storing Globar \f$ z \f$ Distortion

  // look up for charge densities
  TMatrixD *fMatrixChargeA[kNMaxPhi]; //[kNPhi] Matrices for storing input charge densities side A
  TMatrixD *fMatrixChargeC[kNMaxPhi];  //[kNPhi] Matrices for storing input charge densities side C
  TMatrixD *fMatrixChargeInverseA[kNMaxPhi];  //[kNPhi] Matrices for storing charge densities from backward algorithm side A
  TMatrixD *fMatrixChargeInverseC[kNMaxPhi]; //[kNPhi] Matrices for storing charge densities from backward algorithm side A

  AliTPC3DCylindricalInterpolator *fInterpolatorChargeA; //-> interpolator for charge densities side A
  AliTPC3DCylindricalInterpolator *fInterpolatorChargeC; //-> interpolator for charge densities side C
  AliTPC3DCylindricalInterpolator *fInterpolatorInverseChargeA; //-> interpolator for inverse charge densities side A
  AliTPC3DCylindricalInterpolator *fInterpolatorInverseChargeC; //-> interpolator for inverse charge densities side C

  AliTPCLookUpTable3DInterpolatorD *fLookupIntDist; //-> interpolator for global distortion
  AliTPCLookUpTable3DInterpolatorD *fLookupIntCorr; //-> interpolator for global correction


  AliTPCLookUpTable3DInterpolatorD *fLookupIntDistA; //-> interpolator for global distortion side A
  AliTPCLookUpTable3DInterpolatorD *fLookupIntCorrA; //-> interpolator for global correction side A
  AliTPCLookUpTable3DInterpolatorD *fLookupIntDistC; //-> interpolator for global distortion side C
  AliTPCLookUpTable3DInterpolatorD *fLookupIntCorrC; //-> interpolator for global correction side C
  AliTPCLookUpTable3DInterpolatorDFull *fLookupIntCorrIrregularA; //-> interpolator for global correction side A (method irregular)
  AliTPCLookUpTable3DInterpolatorDFull *fLookupIntCorrIrregularC; //-> interpolator for global correction side C (method irregular)
  AliTPCLookUpTable3DInterpolatorD *fLookupIntENoDriftA; //-> interpolator for no drift integration side A
  AliTPCLookUpTable3DInterpolatorD *fLookupIntENoDriftC; //-> interpolator for no drift integration side C
  AliTPCLookUpTable3DInterpolatorD *fLookupIntENoDrift; //-> interpolator for no drift integration
  AliTPCLookUpTable3DInterpolatorD *fLookupDistA; //-> interpolator for local distortion side A
  AliTPCLookUpTable3DInterpolatorD *fLookupDistC; //-> interpolator for local distortion side C
  AliTPCLookUpTable3DInterpolatorD *fLookupInverseDistA; //-> interpolator for local distortion (from inverse) side A
  AliTPCLookUpTable3DInterpolatorD *fLookupInverseDistC; //-> interpolator for local distortion (from inverse) side C

  TH3 *fHistogram3DSpaceCharge;  ///<- Histogram with the input space charge histogram - used as an optional input
  TH3 *fHistogram3DSpaceChargeA;  ///<- Histogram with the input space charge histogram - used as an optional input side A
  TH3 *fHistogram3DSpaceChargeC;  ///<- Histogram with the input space charge histogram - used as an optional input side C
  TF1 *fFormulaBoundaryIFCA = NULL; ///<- function define boundary values for V(z) assuming uniformity in phi and r.
  TF1 *fFormulaBoundaryIFCC = NULL; ///<- function define boundary values for V(z) assuming uniformity in phi and r.
  TF1 *fFormulaBoundaryOFCA = NULL; ///<- function define boundary values for V(z) assuming uniformity in phi and r.
  TF1 *fFormulaBoundaryOFCC = NULL; ///<- function define boundary values for V(z) assuming uniformity in phi and r.
  TF1 *fFormulaBoundaryROCA = NULL; ///<- function define boundary values for V(z) assuming uniformity in phi and r.
  TF1 *fFormulaBoundaryROCC = NULL; ///<- function define boundary values for V(z) assuming uniformity in phi and r.
  TF1 *fFormulaBoundaryCE = NULL; ///<- function define boundary values for V(z) assuming uniformity in phi and r.

  TFormula * fFormulaPotentialV = NULL; ///<- potential V(r,rho,z) function
  TFormula * fFormulaChargeRho = NULL; ///<- charge density Rho(r,rho,z) function

  AliTPCPoissonSolver *fPoissonSolver; //-> Pointer to a poisson solver

  void ElectricField(TMatrixD **matricesV, TMatrixD **matricesEr, TMatrixD **matricesEphi, TMatrixD **matricesEz,
                     const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlices, const Float_t gridSizeR,
                     const Float_t gridSizePhi, const Float_t gridSizeZ, const Int_t symmetry,
                     const Float_t innerRadius);

  void
  LocalDistCorrDz(TMatrixD **matricesEr, TMatrixD **matricesEphi, TMatrixD **matricesEz, TMatrixD **matricesDistDrDz,
                  TMatrixD **matricesDistDphiRDz, TMatrixD **matricesDistDz, TMatrixD **matricesCorrDrDz,
                  TMatrixD **matricesCorrDphiRDz, TMatrixD **matricesCorrDz, const Int_t nRRow, const Int_t nZColumn,
                  const Int_t phiSlice, const Float_t gridSizeZ, const Double_t ezField);

  void IntegrateDistCorrDriftLineDz(AliTPCLookUpTable3DInterpolatorD *lookupLocalDist, TMatrixD **matricesGDistDrDz,
                                    TMatrixD **matricesGDistDphiRDz, TMatrixD **matricesGDistDz,
                                    AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr, TMatrixD **matricesGCorrDrDz,
                                    TMatrixD **matricesGCorrDphiRDz, TMatrixD **matricesGCorrDz,


                                    TMatrixD **matricesGCorrIrregularDrDz, TMatrixD **matricesGCorrIrregularDphiRDz,
                                    TMatrixD **matricesGCorrIrregularDz,

                                    TMatrixD **matricesRIrregular, TMatrixD **matricesPhiIrregular,
                                    TMatrixD **matricesZIrregular,

                                    const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice,
                                    const Double_t *rlist, const Double_t *phiList, const Double_t *zList);


  void FillLookUpTable(AliTPCLookUpTable3DInterpolatorD *lookupGlobal, TMatrixD **lookupRDz, TMatrixD **lookupPhiRDz,
                       TMatrixD **lookupDz, const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice,
                       const Double_t *rlist, const Double_t *phiList, const Double_t *zList, Int_t side);


  void FillLookUpTableA(AliTPCLookUpTable3DInterpolatorD *lookupGlobal, TMatrixD **lookupRDz, TMatrixD **lookupPhiRDz,
                        TMatrixD **lookupDz, const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice,
                        const Double_t *rlist, const Double_t *phiList, const Double_t *zList);

  void FillLookUpTableC(AliTPCLookUpTable3DInterpolatorD *lookupGlobal, TMatrixD **lookupRDz, TMatrixD **lookupPhiRDz,
                        TMatrixD **lookupDz, const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice,
                        const Double_t *rlist, const Double_t *phiList, const Double_t *zList);

  Double_t Interpolate3DTableCyl(Int_t order, Double_t r, Double_t z, Double_t phi, Int_t nr, Int_t nz, Int_t nphi,
                                 const Double_t rlist[], const Double_t zList[], const Double_t phiList[],
                                 TMatrixD **matricesArrays, const Int_t zlow);

  Double_t
  InterpolatePhi(const Double_t xArray[], const Int_t ilow, const Int_t nx, const Float_t yArray[], Int_t order,
                 Double_t x);


  void
  IntegrateDistCorrDriftLineDzOpt2(AliTPCLookUpTable3DInterpolatorD *lookupLocalDist, TMatrixD **matricesGDistDrDz,
                                   TMatrixD **matricesGDistDphiRDz, TMatrixD **matricesGDistDz,
                                   AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr, TMatrixD **matricesGCorrDrDz,
                                   TMatrixD **matricesGCorrDphiRDz, TMatrixD **matricesGCorrDz,

                                   const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice, const Double_t *rlist,
                                   const Double_t *phiList, const Double_t *zList);

  Double_t InterpolatePhi(TH3 *h3, const Double_t r, const Double_t phi, const Double_t z);

  void InverseGlobalToLocalDistortion(TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz,
                                      TMatrixD **matricesDistDz, Double_t *rList, Double_t *zList, Double_t *phiList,
                                      const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice, const Int_t nStep,
                                      const Bool_t useCylAC, Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t type);

  void InverseGlobalToLocalDistortionTwoStages(TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz,
                                               TMatrixD **matricesDistDz, Double_t *rList, Double_t *zList,
                                               Double_t *phiList, const Int_t nRRow, const Int_t nZColumn,
                                               const Int_t phiSlice, const Int_t nStep, const Bool_t useCylAC,
                                               Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t type);

  void InverseGlobalToLocalDistortionGlobalInvTable(TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz,
                                                    TMatrixD **matricesDistDz, Double_t *rList, Double_t *zList,
                                                    Double_t *phiList, const Int_t nRRow, const Int_t nZColumn,
                                                    const Int_t phiSlice, const Int_t nStep, const Bool_t useCylAC,
                                                    Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t type);

  void InverseLocalDistortionToElectricField(TMatrixD **matricesEr, TMatrixD **matricesEphi, TMatrixD **matricesEz,
                                             TMatrixD **matricesInvLocalIntErDz, TMatrixD **matricesInvLocalIntEphiDz,
                                             TMatrixD **matricesInvLocalIntEz, TMatrixD **matricesDistDrDz,
                                             TMatrixD **matricesDistDphiRDz, TMatrixD **matricesDistDz,
                                             Double_t *rList, Double_t *zList, Double_t *phiList, const Int_t nRRow,
                                             const Int_t nZColumn, const Int_t phiSlice);

  void InverseElectricFieldToCharge(TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEphi,
                                    TMatrixD **matricesEz, Double_t *rList, Double_t *zList, Double_t *phiList,
                                    const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice);

  void CalculateEField(TMatrixD **matricesArrayV, TMatrixD **matricesEroverEz, TMatrixD **matricesEPhioverEz,
                       TMatrixD **matricesDeltaEz, const Int_t nRRow, const Int_t nZColumn, const Int_t nPhiSlice,
                       const Int_t symmetry, Bool_t rocDisplacement = kFALSE);

  void
  IntegrateEz(TMatrixD **matricesArrayExoverEz, TMatrixD **matricesArrayEx, const Int_t nRRow, const Int_t nZColumn,
              const Int_t nPhiSlice, const Double_t ezField);

  void InitAllocateMemory();

/// \cond CLASSIMP
  ClassDef(AliTPCSpaceCharge3DDriftLine,
  1);
/// \endcond
};

#endif
