#ifndef ALI_TPC_SPACECHARGE3D_DRIFTLINE_H
#define ALI_TPC_SPACECHARGE3D_DRIFTLINE_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \class AliTPCSpaceCharge3DDriftLine
/// \brief This class provides distortion and correction map with integration following electron drift
/// TODO: validate distortion z by comparing with exisiting classes
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Nov 20, 2017

#include "TF1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMatrixD.h"
#include "AliTPCCorrection.h"
#include "AliTPCSpaceCharge3DCalc.h"

class TFormula;
class AliTPCPoissonSolver;

class AliTPCSpaceCharge3DDriftLine : public AliTPCCorrection {
public:
  AliTPCSpaceCharge3DDriftLine();
  AliTPCSpaceCharge3DDriftLine(const char *name, const char *title);
  AliTPCSpaceCharge3DDriftLine(const char *name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice);
  AliTPCSpaceCharge3DDriftLine(const char *name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice,
                               Int_t interpolationOrder, Int_t irregularGridSize, Int_t rbfKernelType);
  virtual ~AliTPCSpaceCharge3DDriftLine();
  void InitSpaceCharge3DPoissonIntegralDz(Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration,
                                          Double_t stopConvergence);
  void InitSpaceCharge3DPoissonIntegralDz(
    Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration, Double_t stopConvergence,
    TMatrixD **matricesErA, TMatrixD **matricesEphiA, TMatrixD **matricesEzA,
    TMatrixD **matricesErC, TMatrixD **matricesEphiC, TMatrixD **matricesEzC,
    TMatrixD **matricesDistDrDzA, TMatrixD **matricesDistDPhiRDzA, TMatrixD **matricesDistDzA,
    TMatrixD **matricesCorrDrDzA, TMatrixD **matricesCorrDPhiRDzA, TMatrixD **matricesCorrDzA,
    TMatrixD **matricesDistDrDzC, TMatrixD **matricesDistDPhiRDzC, TMatrixD **matricesDistDzC,
    TMatrixD **matricesCorrDrDzC, TMatrixD **matricesCorrDPhiRDzC, TMatrixD **matricesCorrDzC,
    TFormula *intErDzTestFunction, TFormula *intEPhiRDzTestFunction, TFormula *intDzTestFunction);

  void
  InitSpaceCharge3DPoisson(Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration, Double_t stopConvergence);
  void ForceInitSpaceCharge3DPoissonIntegralDz(Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration,
                                               Double_t stopConvergence);
  void GetDistortionCyl(const Float_t x[], Short_t roc, Float_t dx[]);
  void GetDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]);
  void GetCorrectionCyl(const Float_t x[], Short_t roc, Float_t dx[]);
  void GetCorrectionCylAC(const Float_t x[], Short_t roc, Float_t dx[]);
  void GetCorrectionCylACIrregular(const Float_t x[], Short_t roc, Float_t dx[]);
  void GetDistortion(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetCorrection(const Float_t x[], Short_t roc, Float_t dx[]);

  Double_t GetChargeCylAC(const Float_t x[], Short_t roc);
  Double_t GetPotentialCylAC(const Float_t x[], Short_t roc);

  Double_t GetInverseChargeCylAC(const Float_t x[], Short_t roc);

  void SetCorrectionType(Int_t correctionType) {
    fSpaceCharge3DCalc.SetCorrectionType(correctionType);
  }

  TH2F *CreateHistogramDistDRInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramDistDRPhiInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramDistDZInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramCorrDRInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramCorrDRPhiInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramCorrDZInXY(Float_t z, Int_t nx, Int_t ny);

  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Double_t norm);
  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D) { SetInputSpaceCharge(hisSpaceCharge3D, 1); }
  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Double_t norm, Int_t side);
  void SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Int_t side) { SetInputSpaceCharge(hisSpaceCharge3D, 1, side); }

  void SetInputSpaceChargeA(TMatrixD **matricesLookUpCharge) {
    fSpaceCharge3DCalc.SetInputSpaceChargeA(matricesLookUpCharge);
  }

  void SetInputSpaceChargeC(TMatrixD **matricesLookUpCharge) {
    fSpaceCharge3DCalc.SetInputSpaceChargeC(matricesLookUpCharge);
  }

  TTree *CreateDistortionTree(Double_t step);

  TTree *CreateDistortionTree(const Int_t nRRowTest, const Int_t nZColTest, const Int_t nPhiSliceTest);

  TH2F *CreateHistogramSCInXY(Float_t z, Int_t nx, Int_t ny);

  TH2F *CreateHistogramSCInZR(Float_t phi, Int_t nz, Int_t nr);

  void SetNRRows(Int_t nRRow) { fSpaceCharge3DCalc.SetNRRows(nRRow); }

  void SetNPhiSlices(Int_t nPhiSlice) { fSpaceCharge3DCalc.SetNPhiSlices(nPhiSlice); }

  void SetNZColumns(Int_t nZColumn) { fSpaceCharge3DCalc.SetNZColumns(nZColumn); }

  Int_t GetNRRows() { return fSpaceCharge3DCalc.GetNRRows(); }

  Int_t GetNPhiSlices() { return fSpaceCharge3DCalc.GetNPhiSlices(); }

  Int_t GetNZColumns() { return fSpaceCharge3DCalc.GetNZColumns(); }

  void SetPoissonSolver(AliTPCPoissonSolver *poissonSolver) {
    fSpaceCharge3DCalc.SetPoissonSolver(poissonSolver);
  }

  AliTPCPoissonSolver *GetPoissonSolver() { return fSpaceCharge3DCalc.GetPoissonSolver(); }

  void SetInterpolationOrder(Int_t order) { fSpaceCharge3DCalc.SetInterpolationOrder(order); }

  Int_t GetInterpolationOrder() { return fSpaceCharge3DCalc.GetInterpolationOrder(); }

  void SetOmegaTauT1T2(Float_t omegaTau, Float_t t1, Float_t t2) {
    fT1 = t1;
    fT2 = t2;
    fSpaceCharge3DCalc.SetOmegaTauT1T2(omegaTau, t1, t2);
  }

  void SetC0C1(Float_t c0, Float_t c1) {
    fSpaceCharge3DCalc.SetC0C1(c0, c1);
  }

  Float_t GetC0() const { return fSpaceCharge3DCalc.GetC0(); }

  Float_t GetC1() const { return fSpaceCharge3DCalc.GetC1(); }

  void SetCorrectionFactor(Float_t correctionFactor) { fSpaceCharge3DCalc.SetCorrectionFactor(correctionFactor); }

  Float_t GetCorrectionFactor() const { return fSpaceCharge3DCalc.GetCorrectionFactor(); }

  void InverseDistortionMaps(TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEPhi,
                             TMatrixD **matricesEz, TMatrixD **matricesInvLocalIntErDz,
                             TMatrixD **, TMatrixD **matricesInvLocalEz,
                             TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz, TMatrixD **matricesDistDz,
                             const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice, const Int_t nStep,
                             const Bool_t useCylAC, Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t interpType);

  void InverseDistortionMapsNoDrift(TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEPhi,
                                    TMatrixD **matricesEz, TMatrixD **matricesInvLocalIntErDz,
                                    TMatrixD **matricesInvLocalIntEPhiDz, TMatrixD **matricesInvLocalEz,
                                    TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz,
                                    TMatrixD **matricesDistDz, const Int_t nRRow, const Int_t nZColumn,
                                    const Int_t phiSlice);

  void GetCorrectionCylNoDrift(const Float_t x[], const Short_t roc, Float_t dx[]);

  void GetDistortionCylNoDrift(const Float_t x[], Short_t roc, Float_t dx[]);

  void InverseGlobalToLocalDistortionNoDrift(TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz,
                                             TMatrixD **matricesDistDz, Double_t *rList, Double_t *zList,
                                             Double_t *phiList, const Int_t nRRow, const Int_t nZColumn,
                                             const Int_t phiSlice);

  void GetChargeDensity(TMatrixD **matricesChargeA, TMatrixD **matricesChargeC, TH3 *spaceChargeHistogram3D,
                        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice);

  void GetInverseLocalDistortionCyl(const Float_t x[], Short_t roc, Float_t dx[]);

  void GetLocalDistortionCyl(const Float_t x[], Short_t roc, Float_t dx[]);

  void SetIrregularGridSize(Int_t size) { fSpaceCharge3DCalc.SetIrregularGridSize(size); }

  Int_t GetIrregularGridSize() { return fSpaceCharge3DCalc.GetIrregularGridSize(); }

  Int_t GetRBFKernelType() { return fSpaceCharge3DCalc.GetRBFKernelType(); }

  void SetPotentialBoundaryAndChargeFormula(TFormula *vTestFunction, TFormula *rhoTestFunction);

  void SetBoundaryIFCA(TF1 *f1) {
    fSpaceCharge3DCalc.SetBoundaryIFCA(f1);
  }

  void SetBoundaryIFCC(TF1 *f1) {
    fSpaceCharge3DCalc.SetBoundaryIFCC(f1);
  }

  void SetBoundaryOFCA(TF1 *f1) {
    fSpaceCharge3DCalc.SetBoundaryOFCA(f1);
  }

  void SetBoundaryOFCC(TF1 *f1) {
    fSpaceCharge3DCalc.SetBoundaryOFCC(f1);
  }

  void SetBoundaryROCA(TF1 *f1) {
    fSpaceCharge3DCalc.SetBoundaryROCA(f1);
  }

  void SetBoundaryROCC(TF1 *f1) {
    fSpaceCharge3DCalc.SetBoundaryROCC(f1);
  }

  void SetBoundaryCE(TF1 *f1) {
    fSpaceCharge3DCalc.SetBoundaryCE(f1);
  }

  void SetElectricFieldFormula(TFormula *formulaEr, TFormula *formulaEPhi, TFormula *formulaEz) {
    fSpaceCharge3DCalc.SetElectricFieldFormula(formulaEr, formulaEPhi, formulaEz);
  }

  Float_t GetSpaceChargeDensity(Float_t r, Float_t phi, Float_t z);
  Float_t GetPotential(Float_t r, Float_t phi, Float_t z);
  void GetElectricFieldCyl(const Float_t x[], Short_t roc, Double_t dx[]);
  void Init();

private:
AliTPCSpaceCharge3DCalc fSpaceCharge3DCalc; // Lookup table calculator

/// \cond CLASSIMP
  ClassDef(AliTPCSpaceCharge3DDriftLine,
  1);
/// \endcond
};


#endif //ALIROOT_ALITPCSPACECHARGE3DDRIFTLINE_H
