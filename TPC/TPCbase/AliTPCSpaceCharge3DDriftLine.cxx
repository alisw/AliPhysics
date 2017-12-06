/*************************************************************************
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

#include "AliTPCSpaceCharge3DDriftLine.h"

/// \cond CLASSIMP
ClassImp(AliTPCSpaceCharge3DDriftLine)
/// \endcond

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Default values
/// ~~~
/// fInterpolationOrder = 5; // interpolation cubic spline with 5 points
/// fNRRows = 129;
/// fNPhiSlices = 180; // the maximum of phi-slices so far = (8 per sector)
/// fNZColumns = 129; // the maximum on column-slices so  ~ 2cm slicing
/// ~~~
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine()
        : AliTPCCorrection(), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE), fInterpolationOrder(5),
          fIrregularGridSize(3), fRBFKernelType(0), fNRRows(129), fNZColumns(129), fNPhiSlices(180),
          fCorrectionType(0) {
  InitAllocateMemory();
}

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Default values
/// ~~~
/// fInterpolationOrder = 5; interpolation cubic spline with 5 points
/// fNRRows = 129;
/// fNPhiSlices = 180; // the maximum of phi-slices so far = (8 per sector)
/// fNZColumns = 129; // the maximum on column-slices so  ~ 2cm slicing
/// ~~~
///
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(const char *name, const char *title)
        : AliTPCCorrection(name, title), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE),
          fInterpolationOrder(5),
          fIrregularGridSize(3), fRBFKernelType(0), fNRRows(129), fNZColumns(129), fNPhiSlices(144),
          fCorrectionType(0) {
  InitAllocateMemory();
}

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Member values from params
///
/// \param nRRow Int_t number of grid in r direction
/// \param nZColumn Int_t number of grid in z direction
/// \param nPhiSlice Int_t number of grid in \f$ \phi \f$ direction
///
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(const char *name, const char *title, Int_t nRRow,
                                                           Int_t nZColumn, Int_t nPhiSlice) :
        AliTPCCorrection(name, title), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE),
        fInterpolationOrder(5),
        fIrregularGridSize(3), fRBFKernelType(0), fCorrectionType(0) {
  fNRRows = nRRow;
  fNPhiSlices = nPhiSlice; // the maximum of phi-slices so far = (8 per sector)
  fNZColumns = nZColumn; // the maximum on column-slices so  ~ 2cm slicing

  InitAllocateMemory();
}

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Member values from params
///
/// \param nRRow Int_t number of grid in r direction
/// \param nZColumn Int_t number of grid in z direction
/// \param nPhiSlice Int_t number of grid in \f$ \phi \f$ direction
/// \param interpolationOrder Int_t order of interpolation
/// \param strategy Int_t strategy for global distortion
/// \param rbfKernelType Int_t strategy for global distortion
///
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(
        const char *name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice, Int_t interpolationOrder,
        Int_t irregularGridSize, Int_t rbfKernelType)
        : AliTPCCorrection(name, title), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE),
          fCorrectionType(0) {
  fInterpolationOrder = interpolationOrder;
  fIrregularGridSize = irregularGridSize;

  fNRRows = nRRow;
  fNPhiSlices = nPhiSlice;
  fNZColumns = nZColumn;
  fRBFKernelType = rbfKernelType;

  InitAllocateMemory();
}

/// Memory allocation for working/output memory
///
void AliTPCSpaceCharge3DDriftLine::InitAllocateMemory() {
  fListR = new Double_t[fNRRows];
  fListPhi = new Double_t[fNPhiSlices];
  fListZ = new Double_t[fNZColumns];
  fListZA = new Double_t[fNZColumns];
  fListZC = new Double_t[fNZColumns];

  // allocate for boundary
  Int_t len = 2 * fNPhiSlices * (fNZColumns + fNRRows) - (4 * fNPhiSlices);
  fListPotentialBoundaryA = new Double_t[len];
  fListPotentialBoundaryC = new Double_t[len];

  Int_t phiSlicesPerSector = fNPhiSlices / kNumSector;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (fNRRows - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (fNZColumns - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / fNPhiSlices;

  for (Int_t k = 0; k < fNPhiSlices; k++) fListPhi[k] = gridSizePhi * k;
  for (Int_t i = 0; i < fNRRows; i++) fListR[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < fNZColumns; j++) fListZ[j] = (j * gridSizeZ) - fgkTPCZ0;

  for (Int_t j = 0; j < fNZColumns; j++) {
    fListZA[j] = (j * gridSizeZ);
  }

  for (Int_t j = 0; j < fNZColumns; j++) {
    fListZC[j] = (j * gridSizeZ);
  }

  for (Int_t k = 0; k < fNPhiSlices; k++) {

    fMatrixIntDistDrEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDPhiREzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDzA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntDistDrEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDPhiREzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDzC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDPhiREzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDzA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDPhiREzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDzC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixErOverEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixEPhiOverEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixDeltaEzA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixErOverEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixEPhiOverEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixDeltaEzC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEzIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);        //[kNPhi]
    fMatrixIntCorrDPhiREzIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);   //[kNPhi]
    fMatrixIntCorrDzIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixRListIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixPhiListIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixZListIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEzIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);        //[kNPhi]
    fMatrixIntCorrDPhiREzIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);   //[kNPhi]
    fMatrixIntCorrDzIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixRListIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixPhiListIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixZListIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixChargeA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixChargeC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixChargeInverseA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixChargeInverseC[k] = new TMatrixD(fNRRows, fNZColumns);
  }

  fLookupIntDistA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntDistDrEzA, fListR, fNPhiSlices, fMatrixIntDistDPhiREzA, fListPhi,
                  fNZColumns, fMatrixIntDistDzA, fListZA, fInterpolationOrder);
  fLookupIntDistC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntDistDrEzC, fListR, fNPhiSlices, fMatrixIntDistDPhiREzC, fListPhi,
                  fNZColumns, fMatrixIntDistDzC, fListZC, fInterpolationOrder);
  fLookupIntCorrA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntCorrDrEzA, fListR, fNPhiSlices, fMatrixIntCorrDPhiREzA, fListPhi,
                  fNZColumns, fMatrixIntCorrDzA, fListZA, fInterpolationOrder);
  fLookupIntCorrC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntCorrDrEzC, fListR, fNPhiSlices, fMatrixIntCorrDPhiREzC, fListPhi,
                  fNZColumns, fMatrixIntCorrDzC, fListZC, fInterpolationOrder);

  fLookupIntENoDriftA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixErOverEzA, fListR, fNPhiSlices, fMatrixEPhiOverEzA, fListPhi,
                  fNZColumns, fMatrixDeltaEzA, fListZA, fInterpolationOrder);
  fLookupIntENoDriftC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixErOverEzC, fListR, fNPhiSlices, fMatrixEPhiOverEzC, fListPhi,
                  fNZColumns, fMatrixDeltaEzC, fListZC, fInterpolationOrder);
  fLookupIntCorrIrregularA =
          new AliTPCLookUpTable3DInterpolatorIrregularD(
                  fNRRows, fMatrixIntCorrDrEzIrregularA, fMatrixRListIrregularA, fNPhiSlices,
                  fMatrixIntCorrDPhiREzIrregularA, fMatrixPhiListIrregularA,  fNZColumns,
                  fMatrixIntCorrDzIrregularA, fMatrixZListIrregularA, 2, GetIrregularGridSize(),
                  GetIrregularGridSize(), GetIrregularGridSize(), 1);

  fLookupIntCorrIrregularC =
          new AliTPCLookUpTable3DInterpolatorIrregularD(
                  fNRRows, fMatrixIntCorrDrEzIrregularC, fMatrixRListIrregularC,  fNPhiSlices,
                  fMatrixIntCorrDPhiREzIrregularC, fMatrixPhiListIrregularC,  fNZColumns,
                  fMatrixIntCorrDzIrregularC, fMatrixZListIrregularC,  2, GetIrregularGridSize(),
                  GetIrregularGridSize(), GetIrregularGridSize(), 1);

  fInterpolatorChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorChargeC = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeC = new AliTPC3DCylindricalInterpolator();

  fInterpolatorChargeA->SetNR(fNRRows);
  fInterpolatorChargeA->SetNZ(fNZColumns);
  fInterpolatorChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorChargeA->SetRList(fListR);
  fInterpolatorChargeA->SetZList(fListZA);
  fInterpolatorChargeA->SetPhiList(fListPhi);
  fInterpolatorChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorChargeC->SetNR(fNRRows);
  fInterpolatorChargeC->SetNZ(fNZColumns);
  fInterpolatorChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorChargeC->SetRList(fListR);
  fInterpolatorChargeC->SetZList(fListZC);
  fInterpolatorChargeC->SetPhiList(fListPhi);
  fInterpolatorChargeC->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeA->SetNR(fNRRows);
  fInterpolatorInverseChargeA->SetNZ(fNZColumns);
  fInterpolatorInverseChargeA->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeA->SetRList(fListR);
  fInterpolatorInverseChargeA->SetZList(fListZA);
  fInterpolatorInverseChargeA->SetPhiList(fListPhi);
  fInterpolatorInverseChargeA->SetOrder(fInterpolationOrder);

  fInterpolatorInverseChargeC->SetNR(fNRRows);
  fInterpolatorInverseChargeC->SetNZ(fNZColumns);
  fInterpolatorInverseChargeC->SetNPhi(fNPhiSlices);
  fInterpolatorInverseChargeC->SetRList(fListR);
  fInterpolatorInverseChargeC->SetZList(fListZC);
  fInterpolatorInverseChargeC->SetPhiList(fListPhi);
  fInterpolatorInverseChargeC->SetOrder(fInterpolationOrder);

  fLookupDistA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, NULL, fListR, fNPhiSlices, NULL, fListPhi, fNZColumns, NULL, fListZA,
                  fInterpolationOrder);

  fLookupDistC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, NULL, fListR, fNPhiSlices, NULL, fListPhi, fNZColumns, NULL, fListZA,
                  fInterpolationOrder);

  fLookupInverseDistA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, NULL, fListR, fNPhiSlices, NULL, fListPhi, fNZColumns, NULL, fListZA,
                  fInterpolationOrder);

  fLookupInverseDistC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, NULL, fListR, fNPhiSlices, NULL, fListPhi, fNZColumns, NULL, fListZA,
                  fInterpolationOrder);

  fLookupIntCorrIrregularA->SetKernelType(fRBFKernelType);
  fLookupIntCorrIrregularC->SetKernelType(fRBFKernelType);
}

/// Destruction for AliTPCSpaceCharge3DDriftLine
/// Deallocate memory for lookup table and charge distribution
///
AliTPCSpaceCharge3DDriftLine::~AliTPCSpaceCharge3DDriftLine() {
  for (Int_t k = 0; k < fNPhiSlices; k++) {
    delete fMatrixIntDistDrEzA[k];
    delete fMatrixIntDistDPhiREzA[k];
    delete fMatrixIntDistDzA[k];
    delete fMatrixIntDistDrEzC[k];
    delete fMatrixIntDistDPhiREzC[k];
    delete fMatrixIntDistDzC[k];
    delete fMatrixIntCorrDrEzA[k];
    delete fMatrixIntCorrDPhiREzA[k];
    delete fMatrixIntCorrDzA[k];
    delete fMatrixIntCorrDrEzC[k];
    delete fMatrixIntCorrDPhiREzC[k];
    delete fMatrixIntCorrDzC[k];
    delete fMatrixErOverEzA[k];
    delete fMatrixEPhiOverEzA[k];
    delete fMatrixDeltaEzA[k];
    delete fMatrixErOverEzC[k];
    delete fMatrixEPhiOverEzC[k];
    delete fMatrixDeltaEzC[k];
    delete fMatrixIntCorrDrEzIrregularA[k];
    delete fMatrixIntCorrDPhiREzIrregularA[k];
    delete fMatrixIntCorrDzIrregularA[k];
    delete fMatrixRListIrregularA[k];
    delete fMatrixPhiListIrregularA[k];
    delete fMatrixZListIrregularA[k];
    delete fMatrixIntCorrDrEzIrregularC[k];
    delete fMatrixIntCorrDPhiREzIrregularC[k];
    delete fMatrixIntCorrDzIrregularC[k];
    delete fMatrixRListIrregularC[k];
    delete fMatrixPhiListIrregularC[k];
    delete fMatrixZListIrregularC[k];
    delete fMatrixChargeA[k];
    delete fMatrixChargeC[k];
    delete fMatrixChargeInverseA[k];
    delete fMatrixChargeInverseC[k];
  }
  delete[] fListR;
  delete[] fListPhi;
  delete[] fListZ;
  delete[] fListZA;
  delete[] fListZC;

  delete fLookupIntDistA;
  delete fLookupIntDistC;
  delete fLookupIntENoDriftA;
  delete fLookupIntENoDriftC;
  delete fLookupIntCorrA;
  delete fLookupIntCorrC;
  delete fLookupIntCorrIrregularA;
  delete fLookupIntCorrIrregularC;
  delete fLookupDistA;
  delete fLookupDistC;
  delete fLookupInverseDistA;
  delete fLookupInverseDistC;

  delete fInterpolatorChargeA;
  delete fInterpolatorChargeC;
  delete fInterpolatorInverseChargeA;
  delete fInterpolatorInverseChargeC;

  delete[] fListPotentialBoundaryA;
  delete[] fListPotentialBoundaryC;
}

/// Init copy from AliTPCSpaceCharge3D
void AliTPCSpaceCharge3DDriftLine::Init() {
  AliMagF *magF = (AliMagF *) TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magnetic field - not initialized");
  Double_t bzField = magF->SolenoidField() / 10.; //field in T
  AliTPCParam *param = AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vDrift = param->GetDriftV() / 1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField * 10) * vDrift / ezField;
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt, fT1, fT2);
}

/// Creating look-up tables of Correction/Distortion by integration following
/// drift line, input from space charge 3d histogram (fSpaceCharge3D) and boundary values are filled with zeroes
///
/// TODO: provide an interface for setting boundary values
///
/// The algorithm and implementations of this function is the following:
///
/// Do for each side A,C
///
/// 1) Solving \f$ \nabla^2 \Phi(r,\phi,z) = -  \rho(r,\phi,z)\f$
/// ~~~ Calling poisson solver
/// fPoissonSolver->PoissonSolver3D( matricesV, matricesCharge, nRRow, nZColumn, phiSlice, maxIteration, symmetry ) ;
/// ~~~
///
/// 2) Get the electric field \f$ \vec{E} = - \nabla \Phi(r,\phi,z) \f$
/// ~~~
/// ElectricField( matricesV, matricesEr,  matricesEPhi, matricesEz, nRRow, nZColumn, phiSlice,
/// gridSizeR, gridSizePhi ,gridSizeZ,symmetry, fgkIFCRadius);
/// ~~~
///
/// 3) Calculate local distortion and correction, using Langevin formula
/// ~~~ cxx
/// LocalDistCorrDz (matricesEr, matricesEPhi, 	matricesEz,
///	matricesDistDrDz,  matricesDistDPhiRDz, matricesDistDz,
///	matricesCorrDrDz,  matricesCorrDPhiRDz, matricesCorrDz,
///	nRRow,  nZColumn, phiSlice, gridSizeZ, ezField);
/// ~~~
///
/// 4) Integrate distortion by following the drift line
///
/// 5) Fill look up table for Integral distortion
///
/// 6) Fill look up table for Integral correction
///
/// \param nRRow Int_t Number of nRRow in r-direction
/// \param nZColumn Int_t Number of nZColumn in z-direction
/// \param phiSlice Int_t Number of phi slice in \f$ phi \f$ direction
/// \param maxIteration Int_t Maximum iteration for poisson solver
/// \param stoppingConvergence Convergence error stopping condition for poisson solver
///
/// \post Lookup tables for distortion:
/// ~~~
/// fLookUpIntDistDrEz,fLookUpIntDistDPhiREz,fLookUpIntDistDz
/// ~~~
/// and correction:
/// ~~~
/// fLookUpIntCorrDrEz,fLookUpIntCorrDPhiREz,fLookUpIntCorrDz
/// ~~~
/// are initialized
///
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoissonIntegralDz(Int_t nRRow, Int_t nZColumn, Int_t phiSlice,
                                                                      Int_t maxIteration,
                                                                      Double_t stoppingConvergence) {
  Int_t phiSlicesPerSector = phiSlice / kNumSector;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;
  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;

  // local variables
  Float_t radius0, phi0, z0;

  // memory allocation for temporary matrices:
  // potential (boundary values), charge distribution
  TMatrixD *matricesV[phiSlice], *matricesCharge[phiSlice];
  TMatrixD *matricesEr[phiSlice], *matricesEPhi[phiSlice], *matricesEz[phiSlice];
  TMatrixD *matricesDistDrDz[phiSlice], *matricesDistDPhiRDz[phiSlice], *matricesDistDz[phiSlice];
  TMatrixD *matricesCorrDrDz[phiSlice], *matricesCorrDPhiRDz[phiSlice], *matricesCorrDz[phiSlice];
  TMatrixD *matricesGDistDrDz[phiSlice], *matricesGDistDPhiRDz[phiSlice], *matricesGDistDz[phiSlice];
  TMatrixD *matricesGCorrDrDz[phiSlice], *matricesGCorrDPhiRDz[phiSlice], *matricesGCorrDz[phiSlice];

  for (Int_t k = 0; k < phiSlice; k++) {
    matricesV[k] = new TMatrixD(nRRow, nZColumn);
    matricesCharge[k] = new TMatrixD(nRRow, nZColumn);
    matricesEr[k] = new TMatrixD(nRRow, nZColumn);
    matricesEPhi[k] = new TMatrixD(nRRow, nZColumn);
    matricesEz[k] = new TMatrixD(nRRow, nZColumn);
    matricesDistDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesDistDPhiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesDistDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDPhiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGDistDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGDistDPhiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGDistDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGCorrDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGCorrDPhiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGCorrDz[k] = new TMatrixD(nRRow, nZColumn);

  }

  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t rList[nRRow], zList[nZColumn], phiList[phiSlice];

  // pointer to current TF1 for potential boundary values
  TF1 *f1BoundaryIFC = NULL;
  TF1 *f1BoundaryOFC = NULL;
  TF1 *f1BoundaryROC = NULL;

  TStopwatch w;

  for (Int_t k = 0; k < phiSlice; k++) phiList[k] = gridSizePhi * k;
  for (Int_t i = 0; i < nRRow; i++) rList[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < nZColumn; j++) zList[j] = j * gridSizeZ;


  // allocate look up local distortion
  AliTPCLookUpTable3DInterpolatorD *lookupLocalDist =
          new AliTPCLookUpTable3DInterpolatorD(
                  nRRow, matricesDistDrDz, rList, phiSlice, matricesDistDPhiRDz, phiList, nZColumn, matricesDistDz,
                  zList, fInterpolationOrder);

  // allocate look up local correction
  AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr =
          new AliTPCLookUpTable3DInterpolatorD(
                  nRRow, matricesCorrDrDz, rList, phiSlice, matricesCorrDPhiRDz, phiList, nZColumn, matricesCorrDz,
                  zList, fInterpolationOrder);

  // allocate look up for global distortion
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalDist =
          new AliTPCLookUpTable3DInterpolatorD(
                  nRRow, matricesGDistDrDz, rList, phiSlice, matricesGDistDPhiRDz, phiList, nZColumn, matricesGDistDz,
                  zList, fInterpolationOrder);
  // allocate look up for global distortion
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalCorr =
          new AliTPCLookUpTable3DInterpolatorD(
                  nRRow, matricesGCorrDrDz, rList, phiSlice, matricesGCorrDPhiRDz, phiList, nZColumn, matricesGCorrDz,
                  zList, fInterpolationOrder);

  // should be set, in another place
  const Int_t symmetry = 0; // fSymmetry

  // for irregular
  TMatrixD **matricesIrregularDrDz = NULL;
  TMatrixD **matricesIrregularDPhiRDz = NULL;
  TMatrixD **matricesIrregularDz = NULL;
  TMatrixD **matricesPhiIrregular = NULL;
  TMatrixD **matricesRIrregular = NULL;
  TMatrixD **matricesZIrregular = NULL;

  // for charge
  TMatrixD **matricesLookUpCharge = NULL;
  AliTPC3DCylindricalInterpolator *chargeInterpolator = NULL;

  Double_t *potentialBoundary = NULL;

  TMatrixD *matrixV;
  TMatrixD *matrixCharge;

  Int_t pIndex = 0;

  // do if look up table haven't be initialized
  if (!fInitLookUp) {
    // initialize for working memory
    for (Int_t side = 0; side < 2; side++) {
      // zeroing global distortion/correction
      for (Int_t k = 0; k < phiSlice; k++) {
        matricesDistDrDz[k]->Zero();
        matricesDistDPhiRDz[k]->Zero();
        matricesDistDz[k]->Zero();
        matricesCorrDrDz[k]->Zero();
        matricesCorrDPhiRDz[k]->Zero();
        matricesCorrDz[k]->Zero();

        matricesGDistDrDz[k]->Zero();
        matricesGDistDPhiRDz[k]->Zero();
        matricesGDistDz[k]->Zero();
        matricesGCorrDrDz[k]->Zero();
        matricesGCorrDPhiRDz[k]->Zero();
        matricesGCorrDz[k]->Zero();

        if (side == 0) {
          matricesIrregularDrDz = fMatrixIntCorrDrEzIrregularA;
          matricesIrregularDPhiRDz = fMatrixIntCorrDPhiREzIrregularA;
          matricesIrregularDz = fMatrixIntCorrDzIrregularA;
          matricesPhiIrregular = fMatrixPhiListIrregularA;
          matricesRIrregular = fMatrixRListIrregularA;
          matricesZIrregular = fMatrixZListIrregularA;
          matricesLookUpCharge = fMatrixChargeA;
          chargeInterpolator = fInterpolatorChargeA;
          fLookupDistA->SetLookUpR(matricesDistDrDz);
          fLookupDistA->SetLookUpPhi(matricesDistDPhiRDz);
          fLookupDistA->SetLookUpZ(matricesDistDz);
          potentialBoundary = fListPotentialBoundaryA;
          f1BoundaryIFC = fFormulaBoundaryIFCA;
          f1BoundaryOFC = fFormulaBoundaryOFCA;
          f1BoundaryROC = fFormulaBoundaryROCA;
        } else {
          matricesIrregularDrDz = fMatrixIntCorrDrEzIrregularC;
          matricesIrregularDPhiRDz = fMatrixIntCorrDPhiREzIrregularC;
          matricesIrregularDz = fMatrixIntCorrDzIrregularC;
          matricesPhiIrregular = fMatrixPhiListIrregularC;
          matricesRIrregular = fMatrixRListIrregularC;
          matricesZIrregular = fMatrixZListIrregularC;
          matricesLookUpCharge = fMatrixChargeC;
          chargeInterpolator = fInterpolatorChargeC;
          fLookupDistC->SetLookUpR(matricesDistDrDz);
          fLookupDistC->SetLookUpPhi(matricesDistDPhiRDz);
          fLookupDistC->SetLookUpZ(matricesDistDz);
          potentialBoundary = fListPotentialBoundaryC;
          f1BoundaryIFC = fFormulaBoundaryIFCC;
          f1BoundaryOFC = fFormulaBoundaryOFCC;
          f1BoundaryROC = fFormulaBoundaryROCC;
        }
      }

      // fill the potential boundary
      // guess the initial potential
      // fill also charge
      //pIndex = 0;

      AliInfo(Form("Step = 0: Fill Boundary and Charge Densities"));
      for (Int_t k = 0; k < phiSlice; k++) {
        phi0 = k * gridSizePhi;
        matrixV = matricesV[k];
        matrixCharge = matricesCharge[k];
        for (Int_t i = 0; i < nRRow; i++) {
          radius0 = fgkIFCRadius + i * gridSizeR;
          for (Int_t j = 0; j < nZColumn; j++) {
            z0 = j * gridSizeZ;
            (*matrixCharge)(i, j) = chargeInterpolator->GetValue(rList[i], phiList[k], zList[j]);
            (*matrixV)(i, j) = 0.0; // fill zeros

            if (fFormulaPotentialV == NULL) {
              // boundary IFC
              if (i == 0) {
                if (f1BoundaryIFC != NULL) {
                  (*matrixV)(i, j) = f1BoundaryIFC->Eval(z0);
                }
              }
              if (i == (nRRow - 1)) {
                if (f1BoundaryOFC != NULL)
                  (*matrixV)(i, j) = f1BoundaryOFC->Eval(z0);
              }
              if (j == 0) {
                if (fFormulaBoundaryCE) {
                  (*matrixV)(i, j) = fFormulaBoundaryCE->Eval(radius0);
                }
              }
              if (j == (nZColumn - 1)) {
                if (f1BoundaryROC != NULL)
                  (*matrixV)(i, j) = f1BoundaryROC->Eval(radius0);
              }
            } else {
              if ((i == 0) || (i == (nRRow - 1)) || (j == 0) || (j == (nZColumn - 1))) {
                (*matrixV)(i, j) = fFormulaPotentialV->Eval(radius0, phi0, z0);
              }
            }
          }
        }
      }
      AliInfo(Form("Step 0: Preparing Charge interpolator: %f\n", w.CpuTime()));
      AliTPCPoissonSolver::fgConvergenceError = stoppingConvergence;

      fPoissonSolver->SetStrategy(AliTPCPoissonSolver::kMultiGrid);
      (fPoissonSolver->fMgParameters).cycleType = AliTPCPoissonSolver::kFCycle;
      (fPoissonSolver->fMgParameters).isFull3D = kFALSE;
      (fPoissonSolver->fMgParameters).nMGCycle = maxIteration;
      (fPoissonSolver->fMgParameters).maxLoop = 6;

      w.Start();
      fPoissonSolver->PoissonSolver3D(matricesV, matricesCharge, nRRow, nZColumn, phiSlice, maxIteration, symmetry);
      w.Stop();

      AliInfo(Form("Step 1: Poisson solver: %f\n", w.CpuTime()));
      w.Start();
      ElectricField(matricesV,
                    matricesEr, matricesEPhi, matricesEz, nRRow, nZColumn, phiSlice,
                    gridSizeR, gridSizePhi, gridSizeZ, symmetry, fgkIFCRadius);
      w.Stop();

      AliInfo(Form("Step 2: Electric Field Calculation: %f\n", w.CpuTime()));
      w.Start();
      LocalDistCorrDz(matricesEr, matricesEPhi, matricesEz,
                      matricesDistDrDz, matricesDistDPhiRDz, matricesDistDz,
                      matricesCorrDrDz, matricesCorrDPhiRDz, matricesCorrDz,
                      nRRow, nZColumn, phiSlice, gridSizeZ, ezField);
      w.Stop();

      // copy to interpolator
      if (side == 0) {
        lookupLocalDist->CopyFromMatricesToInterpolator();
        lookupLocalCorr->CopyFromMatricesToInterpolator();
        fLookupDistA->CopyFromMatricesToInterpolator();
      } else {
        lookupLocalDist->CopyFromMatricesToInterpolator();
        lookupLocalCorr->CopyFromMatricesToInterpolator();
        fLookupDistC->CopyFromMatricesToInterpolator();
      }

      AliInfo(Form("Step 3: Local distortion and correction: %f\n", w.CpuTime()));
      w.Start();

      IntegrateDistCorrDriftLineDz(
              lookupLocalDist,
              matricesGDistDrDz, matricesGDistDPhiRDz, matricesGDistDz,
              lookupLocalCorr,
              matricesGCorrDrDz, matricesGCorrDPhiRDz, matricesGCorrDz,
              matricesIrregularDrDz, matricesIrregularDPhiRDz, matricesIrregularDz,
              matricesRIrregular, matricesPhiIrregular, matricesZIrregular,
              nRRow, nZColumn, phiSlice, rList, phiList, zList
      );

      w.Stop();
      AliInfo(Form("Step 4: Global correction/distortion: %f\n", w.CpuTime()));
      w.Start();

      //// copy to 1D interpolator /////
      lookupGlobalDist->CopyFromMatricesToInterpolator();
      lookupGlobalCorr->CopyFromMatricesToInterpolator();
      ////


      w.Stop();
      AliInfo(Form("Step 5: Filling up the look up: %f\n", w.CpuTime()));

      if (side == 0) {
        FillLookUpTable(lookupGlobalDist,
                        fMatrixIntDistDrEzA, fMatrixIntDistDPhiREzA, fMatrixIntDistDzA,
                        nRRow, nZColumn, phiSlice, rList, phiList, zList);

        FillLookUpTable(lookupGlobalCorr,
                        fMatrixIntCorrDrEzA, fMatrixIntCorrDPhiREzA, fMatrixIntCorrDzA,
                        nRRow, nZColumn, phiSlice, rList, phiList, zList);

        fLookupIntDistA->CopyFromMatricesToInterpolator();
        fLookupIntCorrA->CopyFromMatricesToInterpolator();

        fLookupIntCorrIrregularA->CopyFromMatricesToInterpolator();

        AliInfo(" A side done");
      }
      if (side == 1) {
        FillLookUpTable(lookupGlobalDist,
                        fMatrixIntDistDrEzC, fMatrixIntDistDPhiREzC, fMatrixIntDistDzC,
                        nRRow, nZColumn, phiSlice, rList, phiList, zList);

        FillLookUpTable(lookupGlobalCorr,
                        fMatrixIntCorrDrEzC, fMatrixIntCorrDPhiREzC, fMatrixIntCorrDzC,
                        nRRow, nZColumn, phiSlice, rList, phiList, zList);

        fLookupIntDistC->CopyFromMatricesToInterpolator();
        fLookupIntCorrC->CopyFromMatricesToInterpolator();
        fLookupIntCorrIrregularC->CopyFromMatricesToInterpolator();
        AliInfo(" C side done");
      }

    }

    fInitLookUp = kTRUE;
  }



  // memory de-allocation for temporary matrices
  for (Int_t k = 0; k < phiSlice; k++) {
    delete matricesV[k];
    delete matricesCharge[k];
    delete matricesEr[k];
    delete matricesEPhi[k];
    delete matricesEz[k];
    delete matricesDistDrDz[k];
    delete matricesDistDPhiRDz[k];
    delete matricesDistDz[k];

    delete matricesCorrDrDz[k];
    delete matricesCorrDPhiRDz[k];
    delete matricesCorrDz[k];
    delete matricesGDistDrDz[k];
    delete matricesGDistDPhiRDz[k];
    delete matricesGDistDz[k];

    delete matricesGCorrDrDz[k];
    delete matricesGCorrDPhiRDz[k];
    delete matricesGCorrDz[k];

  }
  delete lookupLocalDist;
  delete lookupLocalCorr;
  delete lookupGlobalDist;
  delete lookupGlobalCorr;
}

/// Creating look-up tables of Correction/Distortion by linear integration
/// on z line line with known distributions for potential and spacecharge.
///
/// \param nRRow	Int_t  number of grid in row direction
///	\param nZColumn Int_t number of grid in z direction
/// \param phiSlice 	Int_t number of slices in phi direction
/// \param maxIteration Int_t max iteration for convergence
/// \param stoppingConvergence Double_t stopping criteria for convergence
/// \post Lookup tables for distortion:
/// ~~~
/// fLookUpIntDistDrEz,fLookUpIntDistDPhiREz,fLookUpIntDistDz
/// ~~~ fo
/// and correction:
/// ~~~
/// fLookUpIntCorrDrEz,fLookUpIntCorrDPhiREz,fLookUpIntCorrDz
/// ~~~
/// are initialized
///
void
AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoisson(Int_t nRRow, Int_t nZColumn, Int_t phiSlice, Int_t maxIteration,
                                                       Double_t stoppingConvergence) {
  // Compute grid size for all direction
  Int_t phiSlicesPerSector = phiSlice / kNumSector;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;
  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;
  printf("gridSizeZ in init: %f\n", gridSizeZ);

  // local variables
  Float_t radius0, phi0, z0;


  // memory allocation for temporary matrices:
  // potential (boundary values), charge distribution

  TMatrixD *matricesV[phiSlice], *matricesCharge[phiSlice];
  TMatrixD *matricesEr[phiSlice], *matricesEPhi[phiSlice], *matricesEz[phiSlice];

  for (Int_t k = 0; k < phiSlice; k++) {
    matricesEr[k] = new TMatrixD(nRRow, nZColumn);
    matricesEPhi[k] = new TMatrixD(nRRow, nZColumn);
    matricesEz[k] = new TMatrixD(nRRow, nZColumn);
    matricesV[k] = new TMatrixD(nRRow, nZColumn);
    matricesCharge[k] = new TMatrixD(nRRow, nZColumn);
  }

  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t rList[nRRow], zList[nZColumn], phiList[phiSlice];

  for (Int_t k = 0; k < phiSlice; k++) phiList[k] = gridSizePhi * k;
  for (Int_t i = 0; i < nRRow; i++) rList[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < nZColumn; j++) zList[j] = j * gridSizeZ;
  // should be set, in another place
  const Int_t symmetry = 0;
  // do if look up table haven't be initialized
  if (!fInitLookUp) {
    for (Int_t side = 0; side < 2; side++) {
      for (Int_t k = 0; k < phiSlice; k++) {
        TMatrixD *mV = matricesV[k];
        TMatrixD *mCharge = matricesCharge[k];
        phi0 = phiList[k];
        for (Int_t i = 0; i < nRRow; i++) {
          radius0 = rList[i];
          for (Int_t j = 0; j < nZColumn; j++) {
            z0 = zList[j];
            if (side == 1) z0 = -TMath::Abs(zList[j]);
            if (fHistogram3DSpaceCharge != NULL) {
              // * Boundary values and charge distribution setup
              (*mV)(i, j) = 0.0;
              (*mCharge)(i, j) = -1 * InterpolatePhi(fHistogram3DSpaceCharge, phi0, radius0, z0);
            }
          }
        }
      }
      AliTPCLookUpTable3DInterpolatorD *lookupEField =
              new AliTPCLookUpTable3DInterpolatorD(
                      nRRow,
                      matricesEr,
                      rList, phiSlice,
                      matricesEPhi,
                      phiList, nZColumn,
                      matricesEz,
                      zList,
                      fInterpolationOrder
              );

      AliInfo("Step 1: Solving poisson solver");
      fPoissonSolver->PoissonSolver3D(matricesV, matricesCharge, nRRow, nZColumn, phiSlice, maxIteration, symmetry);
      AliInfo("Step 2: Calculate electric field");
      CalculateEField(
              matricesV,
              matricesEr,
              matricesEPhi,
              matricesEz,
              nRRow,
              nZColumn,
              phiSlice,
              maxIteration,
              symmetry
      );
      lookupEField->CopyFromMatricesToInterpolator();
      AliInfo("Step 3: Fill the ");

      if (side == 0) {
        FillLookUpTable(lookupEField,
                        fMatrixErOverEzA, fMatrixEPhiOverEzA, fMatrixDeltaEzA,
                        nRRow, nZColumn, phiSlice, rList, phiList, zList);
        fLookupIntENoDriftA->CopyFromMatricesToInterpolator();
        AliInfo(" A side done");
      }
      if (side == 1) {
        FillLookUpTable(lookupEField,
                        fMatrixErOverEzC, fMatrixEPhiOverEzC, fMatrixDeltaEzC,
                        nRRow, nZColumn, phiSlice, rList, phiList, zList);
        fLookupIntENoDriftC->CopyFromMatricesToInterpolator();

        AliInfo(" C side done");
      }
      delete lookupEField;
    }
    fInitLookUp = kTRUE;
  }

  for (Int_t k = 0; k < phiSlice; k++) {
    delete matricesV[k];
    delete matricesCharge[k];
    delete matricesEr[k];
    delete matricesEPhi[k];
    delete matricesEz[k];
  }
}

/// Force creating look-up table of Correction/Distortion by integration following
/// drift line.
///
/// \param nRRow Int_t Number of nRRow in r-direction
/// \param nZColumn Int_t Number of nZColumn in z-direction
/// \param phiSlice Int_t Number of phi slices in \f$ phi \f$ direction
/// \param maxIteration Int_t Maximum iteration for poisson solver
/// \param stoppingConvergence Convergence error stopping condition for poisson solver
///
void AliTPCSpaceCharge3DDriftLine::ForceInitSpaceCharge3DPoissonIntegralDz
        (Int_t nRRow, Int_t nZColumn, Int_t phiSlice,
         Int_t maxIteration, Double_t stoppingConvergence) {
  fInitLookUp = kFALSE;
  InitSpaceCharge3DPoissonIntegralDz(nRRow, nZColumn, phiSlice, maxIteration, stoppingConvergence);
}

/// Electric field Calculation:
///
///
/// \param matricesV 
/// \param matricesEr 
/// \param matricesEPhi 
/// \param matricesEz 
/// \param nRRow 
/// \param nZColumn 
/// \param phiSlice 
/// \param gridSizeR 
/// \param gridSizePhi 
/// \param gridSizeZ 
/// \param symmetry 
/// \param innerRadius
///
/// \pre   Matrix matricesV is assumed had been calculated  by Poisson solver
/// \post  Results of  E-fields are calculated by measuring gradient at potential distribution
///
///
///	* Differentiate potential on all direction (r,z and phi)
/// * Non-boundary -> Central difference (3 stencil) TODO: 5 Stencil
///
///   \f$  \nabla_{r} V(r_{i},\phi_{j},z_{k}) \approx -( V_{i+1,j,k} - V_{i-1,j,k}) / (2* h_{r}) \f$
///
///   \f$ -\nabla_{\phi} V(r_{i},\phi_{j},z_{k}) \approx -( V_{i,j-1,k} - V_{i,j+1,k}) / (2* r_{j} * h_{\phi}) \f$
///
///   \f$ -\nabla_{z} V(r_{i},\phi_{j},z_{k}) \approx -( V_{i,j,k+1} - V_{i,j,k-1}) / (2* h_{z}) \f$
///
///   ~~~ cxx
///   matrixEr(i,j) = -1 * ( arrayV(i+1,j) - arrayV(i-1,j) ) / (2*gridSizeR); // r direction
///		matrixEz(i,j) = -1 * ( arrayV(i,j+1) - arrayV(i,j-1) ) / (2*gridSizeZ) ; // z direction
///		matrixEPhi(i,j) = -1 * (signPlus * arrayVP(i,j) - signMinus * arrayVM(i,j) ) / (2*radius*gridSizePhi)
///   ~~~
///
/// * Boundary -> Forward/Backward difference (3 stencil) TODO: 5 Stencil
///
///   \f$ -\nabla_{r} V(r_{0},\phi_{j},z_{k}) \approx -( -0.5 V_{2,j,k} + 2 V_{1,j,k} - 1.5 * V_{0,j,k}) /  h_{r} \f$
///
///   \f$ -\nabla_{r} V(r_{nRRow - 1},\phi_{j},z_{k}) \approx -( 1.5 V_{nRRow-1,j,k} - 2.0 V_{nRRow-2,j,k} + 0.5 V_{nRRow -3,j,k}) / h_{\phi} \f$
///  
void AliTPCSpaceCharge3DDriftLine::ElectricField(TMatrixD **matricesV, TMatrixD **matricesEr, TMatrixD **matricesEPhi,
                                                 TMatrixD **matricesEz, const Int_t nRRow, const Int_t nZColumn,
                                                 const Int_t phiSlice,
                                                 const Float_t gridSizeR, const Float_t gridSizePhi,
                                                 const Float_t gridSizeZ,
                                                 const Int_t symmetry, const Float_t innerRadius) {
  Float_t radius;
  Int_t mPlus, mMinus, signPlus, signMinus;
  for (Int_t m = 0; m < phiSlice; m++) {
    mPlus = m + 1;
    signPlus = 1;
    mMinus = m - 1;
    signMinus = 1;
    if (symmetry == 1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (mPlus > phiSlice - 1) mPlus = phiSlice - 2;
      if (mMinus < 0) mMinus = 1;
    } else if (symmetry == -1) {       // Anti-symmetry in phi
      if (mPlus > phiSlice - 1) {
        mPlus = phiSlice - 2;
        signPlus = -1;
      }
      if (mMinus < 0) {
        mMinus = 1;
        signMinus = -1;
      }
    } else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if (mPlus > phiSlice - 1) mPlus = m + 1 - phiSlice;
      if (mMinus < 0) mMinus = m - 1 + phiSlice;
    }

    TMatrixD &arrayVP = *matricesV[mPlus];
    TMatrixD &arrayVM = *matricesV[mMinus];
    TMatrixD &arrayV = *matricesV[m];
    TMatrixD &matrixEr = *matricesEr[m];
    TMatrixD &matrixEz = *matricesEz[m];
    TMatrixD &matrixEPhi = *matricesEPhi[m];

    // for non-boundary V
    for (Int_t i = 1; i < nRRow - 1; i++) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 1; j < nZColumn - 1; j++) {
        matrixEr(i, j) = -1 * (arrayV(i + 1, j) - arrayV(i - 1, j)) / (2 * gridSizeR); // r direction
        matrixEz(i, j) = -1 * (arrayV(i, j + 1) - arrayV(i, j - 1)) / (2 * gridSizeZ); // z direction
        matrixEPhi(i, j) = -1 * (signPlus * arrayVP(i, j) - signMinus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi direction
      }
    }

    // for boundary-r
    for (Int_t j = 0; j < nZColumn; j++) {
      matrixEr(0, j) = -1 * (-0.5 * arrayV(2, j) + 2.0 * arrayV(1, j) - 1.5 * arrayV(0, j)) /
                      gridSizeR; // forward difference
      matrixEr(nRRow - 1, j) =
              -1 * (1.5 * arrayV(nRRow - 1, j) - 2.0 * arrayV(nRRow - 2, j) + 0.5 * arrayV(nRRow - 3, j)) /
              gridSizeR; // backward difference
    }

    for (Int_t i = 0; i < nRRow; i += nRRow - 1) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 1; j < nZColumn - 1; j++) {
        matrixEz(i, j) = -1 * (arrayV(i, j + 1) - arrayV(i, j - 1)) / (2 * gridSizeZ); // z direction
        matrixEPhi(i, j) = -1 * (signPlus * arrayVP(i, j) - signMinus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi direction
      }
    }

    // for boundary-z
    for (Int_t i = 0; i < nRRow; i++) {
      matrixEz(i, 0) = -1 * (-0.5 * arrayV(i, 2) + 2.0 * arrayV(i, 1) - 1.5 * arrayV(i, 0)) / gridSizeZ;
      matrixEz(i, nZColumn - 1) =
              -1 * (1.5 * arrayV(i, nZColumn - 1) - 2.0 * arrayV(i, nZColumn - 2) + 0.5 * arrayV(i, nZColumn - 3)) /
              gridSizeZ;
    }

    for (Int_t i = 1; i < nRRow - 1; i++) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 0; j < nZColumn; j += nZColumn - 1) {
        matrixEr(i, j) = -1 * (arrayV(i + 1, j) - arrayV(i - 1, j)) / (2 * gridSizeR); // r direction
        matrixEPhi(i, j) = -1 * (signPlus * arrayVP(i, j) - signMinus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi direction
      }
    }

    // corner points for EPhi
    for (Int_t i = 0; i < nRRow; i += nRRow - 1) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 0; j < nZColumn; j += nZColumn - 1) {
        matrixEPhi(i, j) = -1 * (signPlus * arrayVP(i, j) - signMinus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi direction
      }
    }
  }
}

///
/// Local distortion and correction, calculate local distortion/correction
/// based on simplified langevin equation, see internal note ALICE-INT-2010-016.
///
/// <b> Local Distortion </b>
///
/// Local distortion is calculated based on formulation in ALICE-INT-2010-016, this function assume that
/// electric field \f$\vec{E}(r_{i},z_{j},\phi_{m})\f$ is provided.
///
/// First, we calculate integration of the Electric field in z-direction for all direction.
/// Assumption: \f$ z_{0} \f$ is location of CE (Central Electrode) and \f$ z_{nZColumn - 1} \f$ is location of End Plate.
///
/// This integration is in \f$z\f$ direction we can only use trapezoidal rule.
///
/// Let suppose we want to calculate local distortion at \f$(r_{i},z_{j},\phi_{m})\f$.
/// Assume \f$\vec{E}(r_{i},z_{j+1},\phi_{m}) \f$ and \f$\vec{E}(r_{i},z_{j},\phi_{m}) \f$ are known, see Figure \ref fig1 (a),
///
/// \anchor fig1
/// ![Local Distortion](localdist.png)
///
/// Than we can calculate definite integrations for each directions in respect of $z$  from \f$ z_{j} \f$  to \f$ z_{j + 1} \f$   as follows:
///
/// \f$  \int^{z_{j+1}}_{z_{j}} \frac{E_{r}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{r}(r_{i},z_{j},\phi_{m}) + E_{r}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j+1}}_{z_{j}} \frac{E_{\phi}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx  \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{\phi}(r_{i},z_{j},\phi_{m}) + E_{\phi}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j+1}}_{z_{j}} E_{z}(r_{i},z_{j},\phi_{m}) dz  \approx   \frac{h_{z}}{2.0} \left( E_{z}(r_{i},z_{j},\phi_{m}) + E_{z}(r_{i},z_{j+1},\phi_{m}) \right) \f$
///
/// Code sample at \ref impllocaldist is an implementation of the local integration of electric field.
///
/// \anchor impllocaldist
/// ~~~
/// Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = Electric Field (V/cm) Magnitude ~ -400 V/cm;
///
/// localIntErOverEz = (gridSizeZ/2.0)*((*eR)(i,j)+(*eR)(i,j+1))/(-1*ezField) ;
/// localIntEPhiOverEz = (gridSizeZ/2.0)*((*ePhi)(i,j)+(*ePhi)(i,j+1))/(-1*ezField) ;
/// localIntDeltaEz = (gridSizeZ/2.0)*((*eZ)(i,j)+(*eZ)(i,j+1)) ;
/// ~~~
///
///
/// After we have local integrations for electric fields in each direction,
/// local distortion \f$\hat{\delta}(r_{i},z_{j},\phi_{m})\f$ is calculated by simplified Langevin equation (see Figure \ref1 (b) for illustration):
///
/// \f$ \hat{\delta}_{rE}(r_{i},z_{j},\phi_{m}) = c_{0} \int^{z_{j+1}}_{z_{j}} \frac{E_{r}}{E_{z}} dz   + c_{1}  \int^{z_{j+1}}_{z_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEPhiOverEz;
/// ~~~
///
/// \f$ r\hat{\delta}_{\phi E}(r_{i},z_{j},\phi_{m})  = - c_{1} \int^{z_{j+1}}_{z_{j}} \frac{E_{j}}{E_{j}} dz  + c_{0} \int^{z_{j+1}}_{j_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDPhiRDz)(i,j) = fC0*localIntEPhiOverEz - fC1*localIntErOverEz ;
/// ~~~
///
/// \f$ \hat{\delta}_{z}(r_{i},z_{j},\phi_{m})  = \int_{z_{j}}^{z_{j+1}} \frac{v^{\prime}(E)}{v_{0}} (E - E_{0}) dz\f$
///
/// ~~~
/// (*distDz)(i,j) = localIntDeltaEz*fgkdvdE;
/// ~~~
///
/// Where \f$c_{0}\f$ and \f$c_{1}\f$ are constants (see the ALICE-INT-2010-016 for further details).
///
/// <b> Local correction </b>
///
/// Local correction is computed as local distortion where the electric fields are in opposite direction (see Figure \ref fig2 (a)).
///
/// \anchor fig2
/// ![Local Correction](localcorr.png)
///
/// Let suppose we want to calculate local correction at \f$(r_{i},\mathbf{z_{j+1}},\phi_{m})\f$.
/// Assume \f$\vec{E}(r_{i},z_{j+1},\phi_{m}) \f$ and \f$\vec{E}(r_{i},z_{j},\phi_{m}) \f$ are known.
///
/// Than we can calculate definite integrations for each directions in respect of \f$z\f$  from \f$ z_{j+1} \f$  to \f$ z_{j} \f$   as follows:
///
/// \f$  \int^{z_{j}}_{z_{j+1}} \frac{E_{r}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx -1 * \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{r}(r_{i},z_{j},\phi_{m}) + E_{r}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j}}_{z_{j+1}} \frac{E_{\phi}}{E_{z}}(r_{i},z_{j},\phi_{m}) dz  \approx  -1 *  \frac{-1}{\mathrm{ezField}} \frac{h_{z}}{2.0} \left( E_{\phi}(r_{i},z_{j},\phi_{m}) + E_{\phi}(r_{i},z_{j+1},\phi_{m}) \right)\f$
///
/// \f$  \int^{z_{j}}_{z_{j+1}} E_{z}(r_{i},z_{j},\phi_{m}) dz  \approx  -1 *   \frac{h_{z}}{2.0} \left( E_{z}(r_{i},z_{j},\phi_{m}) + E_{z}(r_{i},z_{j+1},\phi_{m}) \right) \f$
///
/// Local correction at \f$\hat{\delta'}(r_{i},\mathbf{z_{j+1}},\phi_{m})\f$ is calculated by simplified Langevin equation (see Figure \ref fig2 (b) for illustration):
///
/// \f$ \hat{\delta'}_{rE}(r_{i},z_{j+1},\phi_{m}) = c_{0} \int^{z_{j}}_{z_{j+1}} \frac{E_{r}}{E_{z}} dz   + c_{1}  \int^{z_{j-1}}_{z_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// \f$ r\hat{\delta'}_{\phi E}(r_{i},z_{j+1},\phi_{m})  = - c_{1} \int^{z_{j}}_{z_{j+1}} \frac{E_{j}}{E_{j}} dz  + c_{0} \int^{z_{j-1}}_{j_{k}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// \f$ \hat{\delta'}_{z}(r_{i},z_{j+1},\phi_{m})  = \int_{z_{j}}^{z_{j+1}} \frac{v^{\prime}(E)}{v_{0}} (E - E_{0}) dz\f$
///
/// For implementation, we use the fact that
///
/// \f$ \hat{\delta'}_{rE}(r_{i},z_{j+1},\phi_{m}) = -1 * \hat{\delta}_{rE}(r_{i},z_{j},\phi_{m}) \f$
///
/// \f$ r\hat{\delta'}_{\phi E}(r_{i},z_{j+1},\phi_{m}) =  -1 *  r\hat{\delta}_{\phi E}(r_{i},z_{j},\phi_{m}) \f$
///
/// \f$ \hat{\delta'}_{z}(r_{i},z_{j+1},\phi_{m}) =  -1 *  \hat{\delta}_{z}(r_{i},z_{j},\phi_{m}) \f$
///
/// ~~~
///	(*corrDrDz)(i,j+1) 		= -1* (*distDrDz)(i,j) ;
/// (*corrDPhiRDz)(i,j+1) = -1* (*distDPhiRDz)(i,j);
/// (*corrDz)(i,j+1)      = -1* (*distDz)(i,j);
/// ~~~
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEPhi TMatrixD** electric field for \f$\phi\f$ component
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component
///	\param matricesDistDrDz TMatrixD**  local distortion \f$\hat{\delta}_{r}\f$
///	\param matricesDistDPhiRDz TMatrixD** local distortion \f$r \hat{\delta}_{\phi}\f$
///	\param matricesDistDz TMatrixD**   local distortion \f$ \hat{\delta}_{z}\f$
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDPhiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param nRRow Int_t Number of nRRow in r-direction
/// \param nZColumn Int_t Number of nZColumn in z-direction
/// \param phiSlice Int_t Number of phi slices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction
/// \param ezField const Double_t ezField calculated from the invoking operation
///
/// \pre matricesEr, matricesEPhi, matrices Ez assume already been calculated
/// \post Local distortion and correction are computed according simplified Langevin equation
/// ~~~
/// matricesDistDrDz,matricesDistDPhiRDz,matricesDistDz
/// ~~~
/// and correction:
/// ~~~
/// matricesCorrDrDz,matricesCorrDPhiRDz,matricesCorrDz
/// ~~~
///
void
AliTPCSpaceCharge3DDriftLine::LocalDistCorrDz(TMatrixD **matricesEr, TMatrixD **matricesEPhi, TMatrixD **matricesEz,
                                              TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz,
                                              TMatrixD **matricesDistDz,
                                              TMatrixD **matricesCorrDrDz, TMatrixD **matricesCorrDPhiRDz,
                                              TMatrixD **matricesCorrDz,
                                              const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice,
                                              const Float_t gridSizeZ,
                                              const Double_t ezField) {
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEPhiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;
  TMatrixD *eR;
  TMatrixD *ePhi;
  TMatrixD *eZ;
  TMatrixD *distDrDz;
  TMatrixD *distDPhiRDz;
  TMatrixD *distDz;
  TMatrixD *corrDrDz;
  TMatrixD *corrDPhiRDz;
  TMatrixD *corrDz;


  // Initialization for j == column-1 integration is 0.0
  for (Int_t m = 0; m < phiSlice; m++) {
    distDrDz = matricesDistDrDz[m];
    distDPhiRDz = matricesDistDPhiRDz[m];
    distDz = matricesDistDz[m];

    corrDrDz = matricesCorrDrDz[m];
    corrDPhiRDz = matricesCorrDPhiRDz[m];
    corrDz = matricesCorrDz[m];

    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDPhiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

      (*corrDrDz)(i, 0) = 0.0;
      (*corrDPhiRDz)(i, 0) = 0.0;
      (*corrDz)(i, 0) = 0.0;
    }
  }

  // for this case
  // use trapezoidal rule assume no ROC displacement
  for (Int_t m = 0; m < phiSlice; m++) {
    eR = matricesEr[m];
    ePhi = matricesEPhi[m];
    eZ = matricesEz[m];
    distDrDz = matricesDistDrDz[m];
    distDPhiRDz = matricesDistDPhiRDz[m];
    distDz = matricesDistDz[m];

    corrDrDz = matricesCorrDrDz[m];
    corrDPhiRDz = matricesCorrDPhiRDz[m];
    corrDz = matricesCorrDz[m];

    for (Int_t j = 0; j < nZColumn - 1; j++) {
      for (Int_t i = 0; i < nRRow; i++) {
        localIntErOverEz = (gridSizeZ / 2.0) * ((*eR)(i, j) + (*eR)(i, j + 1)) / (-1 * ezField);
        localIntEPhiOverEz = (gridSizeZ / 2.0) * ((*ePhi)(i, j) + (*ePhi)(i, j + 1)) / (-1 * ezField);
        localIntDeltaEz = (gridSizeZ / 2.0) * ((*eZ)(i, j) + (*eZ)(i, j + 1));

        (*distDrDz)(i, j) = fC0 * localIntErOverEz + fC1 * localIntEPhiOverEz;
        (*distDPhiRDz)(i, j) = fC0 * localIntEPhiOverEz - fC1 * localIntErOverEz;
        (*distDz)(i, j) = localIntDeltaEz * fgkdvdE * fgkdvdE;// two times?


        (*corrDrDz)(i, j + 1) = -1 * (*distDrDz)(i, j);
        (*corrDPhiRDz)(i, j + 1) = -1 * (*distDPhiRDz)(i, j);
        (*corrDz)(i, j + 1) = -1 * (*distDz)(i, j);

      }
    }
  }
}


/// IntegrateDistCorrDriftLineDz, integration of local distortion by following electron drift

/// See explanation at LocalDistCorrDz
///
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEPhi TMatrixD** electric field for \f$\phi\f$ component
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDPhiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param nRRow Int_t Number of nRRow in r-direction
/// \param nZColumn Int_t Number of nZColumn in z-direction
/// \param phiSlice Int_t Number of phi slices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction
/// \param ezField const Double_t ezField calculate from the invoking operation
///
/// \pre matricesEr, matricesEPhi, matrices Ez are provided
/// \post Local correction are computed according simplified Langevin equation
/// ~~~
/// matricesCorrDz,matricesCorrDPhiRDz,matricesDistDz
/// ~~~
///
void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDz (
        AliTPCLookUpTable3DInterpolatorD *lookupLocalDist, TMatrixD **matricesGDistDrDz,
        TMatrixD **matricesGDistDPhiRDz, TMatrixD **matricesGDistDz, AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr,
        TMatrixD **matricesGCorrDrDz, TMatrixD **matricesGCorrDPhiRDz, TMatrixD **matricesGCorrDz,
        TMatrixD **matricesGCorrIrregularDrDz, TMatrixD **matricesGCorrIrregularDPhiRDz, TMatrixD **matricesGCorrIrregularDz,
        TMatrixD **matricesRIrregular, TMatrixD **matricesPhiIrregular, TMatrixD **matricesZIrregular,
        const Int_t nRRow, const Int_t nZColumn,     const Int_t phiSlice,
        const Double_t *rList, const Double_t *phiList, const Double_t *zList) {
  Float_t dr, dRPhi, dz, ddR, ddRPhi, ddZ;
  Float_t radius0, phi0, z0, radius, phi, z, radiusCorrection;
  radiusCorrection = 0.0;
  radius = 0.0;
  TMatrixD *mDistDrDz;
  TMatrixD *mDistDPhiRDz;
  TMatrixD *mDistDz;
  TMatrixD *mCorrDrDz;
  TMatrixD *mCorrDPhiRDz;
  TMatrixD *mCorrDz;

  TMatrixD *mCorrIrregularDrDz;
  TMatrixD *mCorrIrregularDPhiRDz;
  TMatrixD *mCorrIrregularDz;

  TMatrixD *mRIrregular;
  TMatrixD *mPhiIrregular;
  TMatrixD *mZIrregular;

  Int_t j = nZColumn - 1;
  z0 = zList[j];


  for (Int_t m = 0; m < phiSlice; m++) {
    phi0 = phiList[m];

    mDistDrDz = matricesGDistDrDz[m];
    mDistDPhiRDz = matricesGDistDPhiRDz[m];
    mDistDz = matricesGDistDz[m];

    //
    mCorrDrDz = matricesGCorrDrDz[m];
    mCorrDPhiRDz = matricesGCorrDPhiRDz[m];
    mCorrDz = matricesGCorrDz[m];

    mCorrIrregularDrDz = matricesGCorrIrregularDrDz[m];
    mCorrIrregularDPhiRDz = matricesGCorrIrregularDPhiRDz[m];
    mCorrIrregularDz = matricesGCorrIrregularDz[m];

    mRIrregular = matricesRIrregular[m];
    mPhiIrregular = matricesPhiIrregular[m];
    mZIrregular = matricesZIrregular[m];

    for (Int_t i = 0; i < nRRow; i++) {
      // do from j to 0
      // follow the drift
      radius0 = rList[i];
      phi = phi0;
      radius = radius0;

      dr = 0.0;
      dRPhi = 0.0;
      dz = 0.0;
      ddRPhi = 0.0;

      ///
      (*mDistDrDz)(i, j) = dr;
      (*mDistDPhiRDz)(i, j) = dRPhi;
      (*mDistDz)(i, j) = dz;

//////////////// use irregular grid look up table for correction
      // set
      (*mCorrIrregularDrDz)(i, j) = -dr;
      (*mCorrIrregularDPhiRDz)(i, j) = -dRPhi;
      (*mCorrIrregularDz)(i, j) = -dz;


      // distorted point
      (*mRIrregular)(i, j) = radius0 + dr;
      (*mPhiIrregular)(i, j) = phi0 + (dRPhi / radius0);
      (*mZIrregular)(i, j) = z0 + dz;
///////////////

    }
  }

  // from j one column near end cap
  for (j = nZColumn - 2; j >= 0; j--) {

    z0 = zList[j];
    for (Int_t m = 0; m < phiSlice; m++) {
      phi0 = phiList[m];

      mDistDrDz = matricesGDistDrDz[m];
      mDistDPhiRDz = matricesGDistDPhiRDz[m];
      mDistDz = matricesGDistDz[m];

      //
      mCorrDrDz = matricesGCorrDrDz[m];
      mCorrDPhiRDz = matricesGCorrDPhiRDz[m];
      mCorrDz = matricesGCorrDz[m];

      mCorrIrregularDrDz = matricesGCorrIrregularDrDz[m];
      mCorrIrregularDPhiRDz = matricesGCorrIrregularDPhiRDz[m];
      mCorrIrregularDz = matricesGCorrIrregularDz[m];

      mRIrregular = matricesRIrregular[m];
      mPhiIrregular = matricesPhiIrregular[m];
      mZIrregular = matricesZIrregular[m];

      for (Int_t i = 0; i < nRRow; i++) {
        // do from j to 0
        // follow the drift
        radius0 = rList[i];
        phi = phi0;
        radius = radius0;

        dr = 0.0;
        dRPhi = 0.0;
        dz = 0.0;
        ddRPhi = 0.0;


        // follow the drift line from z=j --> nZColumn - 1
        for (Int_t jj = j; jj < nZColumn; jj++) {
          // interpolation the local distortion for current position
          phi += ddRPhi / radius;
          radius = radius0 + dr;
          z = zList[jj] + dz;

          // regulate phi
          while (phi < 0.0) phi = TMath::TwoPi() + phi;
          while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();

          
          lookupLocalDist->GetValue(radius, phi, z, ddR, ddRPhi, ddZ);

          
          // add local distortion
          dr += ddR;
          dRPhi += ddRPhi;
          dz += ddZ;

        }
        // set the global distortion after following the electron drift
        (*mDistDrDz)(i, j) = dr;
        (*mDistDPhiRDz)(i, j) = dRPhi;
        (*mDistDz)(i, j) = dz;
/////////////// use irregular grid look up table for correction
        // set
        (*mCorrIrregularDrDz)(i, j) = -dr;
        (*mCorrIrregularDPhiRDz)(i, j) = -dRPhi;
        (*mCorrIrregularDz)(i, j) = -dz;


        // distorted point
        (*mRIrregular)(i, j) = radius0 + dr;
        (*mPhiIrregular)(i, j) = phi0 + (dRPhi / radius0);
        (*mZIrregular)(i, j) = z0 + dz;
///////////////

        // put the radius to the original value
        if (j == nZColumn - 2) radiusCorrection = radius0;

        // get global correction from j+1
        dr = (*mCorrDrDz)(i, j + 1);
        dRPhi = (*mCorrDPhiRDz)(i, j + 1);
        dz = (*mCorrDz)(i, j + 1);
        
        radiusCorrection = radius0 + dr;
        phi = phi0 + dRPhi / radiusCorrection;
        z = zList[j + 1] + dz;

        while (phi < 0.0) phi = TMath::TwoPi() + phi;
        while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();

        lookupLocalCorr->GetValue(radiusCorrection, phi, z, ddR, ddRPhi, ddZ);

        dr += ddR;
        dz += ddZ;
        dRPhi += ddRPhi;

        (*mCorrDrDz)(i, j) = dr;
        (*mCorrDPhiRDz)(i, j) = dRPhi;
        (*mCorrDz)(i, j) = dz;

      }
    }
  }
}

///
/// \param lookupGlobal 
/// \param lookupRDz 
/// \param lookupPhiRDz 
/// \param lookupDz 
/// \param nRRow 
/// \param nZColumn 
/// \param phiSlice 
/// \param rList 
/// \param phiList 
/// \param zList 
void AliTPCSpaceCharge3DDriftLine::FillLookUpTable(AliTPCLookUpTable3DInterpolatorD *lookupGlobal, TMatrixD **lookupRDz,
                                                   TMatrixD **lookupPhiRDz, TMatrixD **lookupDz, const Int_t nRRow,
                                                   const Int_t nZColumn, const Int_t phiSlice, const Double_t *rList,
                                                   const Double_t *phiList, const Double_t *zList) {
  Double_t r, phi, z;
  TMatrixD *mR;
  TMatrixD *mPhiR;
  TMatrixD *mDz;

  /// * Interpolate basicLookup tables; once for each rod, then sum the results
  for (Int_t k = 0; k < fNPhiSlices; k++) {
    phi = fListPhi[k];

    mR = lookupRDz[k];
    mPhiR = lookupPhiRDz[k];
    mDz = lookupDz[k];
    for (Int_t j = 0; j < fNZColumns; j++) {
      z = fListZA[j];  // Symmetric solution in Z that depends only on ABS(Z)

      for (Int_t i = 0; i < fNRRows; i++) {
        r = fListR[i];

        lookupGlobal->GetValue(r, phi, z, (*mR)(i, j), (*mPhiR)(i, j), (*mDz)(i, j));

      }
    }
  }
}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetDistortionCyl(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  GetDistortionCylAC(x, roc, dx);

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  Float_t dR, dRPhi, dZ;
  Double_t r, phi, z;
  Int_t sign;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE

  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  if (z > 0)
    fLookupIntDistA->GetValue(r, phi, z, dR, dRPhi, dZ);
  else {
    fLookupIntDistC->GetValue(r, phi, -1 * z, dR, dRPhi, dZ);
    dZ = -1 * dZ;
  }
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dRPhi;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with drift velocity dependency on the Ez field and the overall scaling factor)
}

/// Get Correction from irregular table
///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylACIrregular(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  Double_t dR, dRPhi, dZ;
  Double_t r, phi, z;
  Int_t sign;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (fNRRows - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (fNZColumns - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / fNPhiSlices;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");


  // get distortion from irregular table


  Int_t iAnchor = TMath::FloorNint((r - fgkIFCRadius) / gridSizeR);
  Int_t kAnchor = TMath::FloorNint(phi / gridSizePhi);
  Int_t zAnchor = TMath::FloorNint(z / gridSizeZ);

  if (z > 0)
    fLookupIntCorrIrregularA->GetValue(r, phi, z, dR, dRPhi, dZ, iAnchor, kAnchor, zAnchor, fNRRows / 8 + 1,
                                       fNPhiSlices / 8 + 1, fNZColumns / 8 + 1, 0);
  else {
    fLookupIntCorrIrregularC->GetValue(r, phi, -z, dR, dRPhi, dZ, iAnchor, kAnchor, -zAnchor, fNRRows / 8 + 1,
                                       fNPhiSlices / 8 + 1, fNZColumns / 8 + 1, 0);
    dZ = -1 * dZ;
  }

  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dRPhi;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with drift velocity dependency on the Ez field and the overall scaling factor)

}

/// Get correction regular grid by following electron
/// 
/// \param x 
/// \param roc 
/// \param dx 
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylAC(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  Float_t dR, dRPhi, dZ;
  Double_t r, phi, z;
  Int_t sign;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  if (z > 0)
    fLookupIntCorrA->GetValue(r, phi, z, dR, dRPhi, dZ);
  else {
    fLookupIntCorrC->GetValue(r, phi, -z, dR, dRPhi, dZ);
    dZ = -1 * dZ;
  }
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dRPhi;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with drift velocity dependency on the Ez field and the overall scaling factor)
}

void AliTPCSpaceCharge3DDriftLine::GetDistortion(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  Float_t pCyl[3]; // a point in cylindrical coordinate
  Float_t dCyl[3]; // distortion

  pCyl[0] = TMath::Sqrt(x[0] * x[0] + x[1] * x[1]);
  pCyl[1] = TMath::ATan2(x[1], x[0]);

  // normalize phi
  while (pCyl[1] > TMath::Pi()) pCyl[1] -= TMath::TwoPi();
  while (pCyl[1] < -TMath::Pi()) pCyl[1] += TMath::TwoPi();

  pCyl[2] = x[2];                                         // Create temporary copy of x[2]

  GetDistortionCylAC(pCyl, roc, dCyl);


  // Calculate distorted position
  if (pCyl[0] > 0.0) {
    pCyl[0] = pCyl[0] + fCorrectionFactor * dCyl[0];
    pCyl[1] = pCyl[1] + fCorrectionFactor * dCyl[1] / pCyl[0];
  }

  dCyl[2] = fCorrectionFactor * dCyl[2];

  // distortion in x,y and z
  dx[0] = (pCyl[0] * TMath::Cos(pCyl[1]) - x[0]);
  dx[1] = (pCyl[0] * TMath::Sin(pCyl[1]) - x[1]);
  dx[2] = dCyl[2];

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCyl(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }
  if (fCorrectionType == kRegularInterpolator)
    GetCorrectionCylAC(x, roc, dx);
  else
    GetCorrectionCylACIrregular(x, roc, dx);
}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetCorrection(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  Float_t pCyl[3]; // a point in cylindrical coordinate
  Float_t dCyl[3]; // distortion

  pCyl[0] = TMath::Sqrt(x[0] * x[0] + x[1] * x[1]);
  pCyl[1] = TMath::ATan2(x[1], x[0]);

  // normalize phi
  while (pCyl[1] > TMath::Pi()) pCyl[1] -= TMath::TwoPi();
  while (pCyl[1] < -TMath::Pi()) pCyl[1] += TMath::TwoPi();

  pCyl[2] = x[2];                                         // Create temporary copy of x[2]


  if (fCorrectionType == kRegularInterpolator)
    GetCorrectionCylAC(pCyl, roc, dCyl);
  else
    GetCorrectionCylACIrregular(pCyl, roc, dCyl);


  // Calculate distorted position
  if (pCyl[0] > 0.0) {
    pCyl[0] = pCyl[0] + fCorrectionFactor * dCyl[0];
    pCyl[1] = pCyl[1] + fCorrectionFactor * dCyl[1] / pCyl[0];
  }

  dCyl[2] = fCorrectionFactor * dCyl[2];

  // distortion in x,y and z
  dx[0] = (pCyl[0] * TMath::Cos(pCyl[1]) - x[0]);
  dx[1] = (pCyl[0] * TMath::Sin(pCyl[1]) - x[1]);
  dx[2] = dCyl[2];

}

///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramDistDRInXY(Float_t z, Int_t nx, Int_t ny) {
  /// Simple plot functionality.
  /// Returns a 2d histogram which represents the corrections in radial direction (dr)
  /// in respect to position z within the XY plane.
  /// The histogram nx times ny entries.

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dr_xy", TString::Format("%s: DRInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dr [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3], xCyl[3];
  Float_t r0, phi0;

  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      phi0 = TMath::ATan2(x[1], x[0]);

      while (phi0 > TMath::Pi()) phi0 -= TMath::TwoPi();
      while (phi0 < -TMath::Pi()) phi0 += TMath::TwoPi();
      xCyl[0] = r0;
      xCyl[1] = phi0;
      GetDistortionCylAC(xCyl, roc, dx);

      //Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));

      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        //Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
        h->SetBinContent(ix, iy, dx[0]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in radial direction (dr)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramCorrDRInXY
        (
                Float_t z,
                Int_t nx,
                Int_t ny
        ) {

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dr_xy", TString::Format("%s: DRInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dr [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3], xCyl[3];
  Float_t r0, phi0;

  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      phi0 = TMath::ATan2(x[1], x[0]);

      while (phi0 > TMath::Pi()) phi0 -= TMath::TwoPi();
      while (phi0 < -TMath::Pi()) phi0 += TMath::TwoPi();
      xCyl[0] = r0;
      xCyl[1] = phi0;
      GetCorrectionCylAC(xCyl, roc, dx);

      //Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));

      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        //Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
        h->SetBinContent(ix, iy, dx[0]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in r Phi direction (drPhi)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramDistDRPhiInXY(Float_t z, Int_t nx, Int_t ny) {

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("drPhi_xy", TString::Format("%s: DRPhiInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drPhi [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetDistortion(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        Float_t phi0 = TMath::ATan2(x[1], x[0]);
        Float_t phi1 = TMath::ATan2(x[1] + dx[1], x[0] + dx[0]);

        Float_t dPhi = phi1 - phi0;
        if (dPhi < TMath::Pi()) dPhi += TMath::TwoPi();
        if (dPhi > TMath::Pi()) dPhi -= TMath::TwoPi();

        h->SetBinContent(ix, iy, r0 * dPhi);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }

  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in r phi direction (drPhi)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramCorrDRPhiInXY(Float_t z, Int_t nx, Int_t ny) {

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("drPhi_xy", TString::Format("%s: DRPhiInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drPhi [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetCorrection(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        Float_t phi0 = TMath::ATan2(x[1], x[0]);
        Float_t phi1 = TMath::ATan2(x[1] + dx[1], x[0] + dx[0]);

        Float_t dPhi = phi1 - phi0;
        if (dPhi < TMath::Pi()) dPhi += TMath::TwoPi();
        if (dPhi > TMath::Pi()) dPhi -= TMath::TwoPi();

        h->SetBinContent(ix, iy, r0 * dPhi);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }

  delete tpcParam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in longitudinal direction (dz)
/// in respect to position z within the XY plane.
/// The histogram nx times ny entries.
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramDistDZInXY(Float_t z, Int_t nx, Int_t ny) {

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dz_xy", TString::Format("%s: DZInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dz [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetDistortion(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        h->SetBinContent(ix, iy, dx[2]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramCorrDZInXY(Float_t z, Int_t nx, Int_t ny) {
  /// Simple plot functionality.
  /// Returns a 2d histogram which represents the corrections in longitudinal direction (dz)
  /// in respect to position z within the XY plane.
  /// The histogram nx times ny entries.

  AliTPCParam *tpcParam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dz_xy", TString::Format("%s: DZInXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dz [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3];
  x[2] = z;
  Int_t roc = z > 0. ? 0 : 18; // FIXME
  for (Int_t iy = 1; iy <= ny; ++iy) {
    x[1] = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      x[0] = h->GetXaxis()->GetBinCenter(ix);
      GetCorrection(x, roc, dx);
      Float_t r0 = TMath::Sqrt((x[0]) * (x[0]) + (x[1]) * (x[1]));
      if (tpcParam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcParam->GetPadRowRadii(36, 95)) {
        h->SetBinContent(ix, iy, dx[2]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcParam;
  return h;
}

/// Use 3D space charge map as an optional input
/// The layout of the input histogram is assumed to be: (phi,r,z)
/// Density histogram is  expected to bin in  C/m^3
///
/// Standard histogram interpolation is used in order to use the density at center of bin
///
/// \param hisSpaceCharge3D
/// \param norm
void AliTPCSpaceCharge3DDriftLine::SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Double_t norm) {
  fHistogram3DSpaceCharge = hisSpaceCharge3D;
  fInitLookUp = kFALSE;
}

/// SetInputCharge
///
/// \param hisSpaceCharge3D TH3* histogram for spacecharge
/// \param norm Double_t norm/weight
/// \param side Int_t side = 0 => side A, side = 1 => side C
///
/// side effects: create Charge interpolator
void AliTPCSpaceCharge3DDriftLine::SetInputSpaceCharge(TH3 *hisSpaceCharge3D, Double_t norm, Int_t side) {
  if (side == 0)
    fHistogram3DSpaceChargeA = hisSpaceCharge3D;
  else
    fHistogram3DSpaceChargeC = hisSpaceCharge3D;

  Double_t rMin = hisSpaceCharge3D->GetYaxis()->GetBinCenter(0);
  Double_t rMax = hisSpaceCharge3D->GetYaxis()->GetBinUpEdge(hisSpaceCharge3D->GetYaxis()->GetNbins());
  Double_t zMin = hisSpaceCharge3D->GetZaxis()->GetBinCenter(0);
  Double_t zMax = hisSpaceCharge3D->GetZaxis()->GetBinCenter(hisSpaceCharge3D->GetZaxis()->GetNbins());
  Double_t radius0, z0, phi0;
  TMatrixD *charge;

  for (Int_t k = 0; k < fNPhiSlices; k++) {
    if (side == 0)
      charge = fMatrixChargeA[k];
    else
      charge = fMatrixChargeC[k];

    phi0 = fListPhi[k];
    for (Int_t i = 0; i < fNRRows; i++) {
      radius0 = fListR[i];
      for (Int_t j = 0; j < fNZColumns; j++) {
        z0 = fListZ[j];

        if (radius0 > rMin && radius0 < rMax && z0 > zMin && z0 < zMax) {
          (*charge)(i, j) = norm * InterpolatePhi(hisSpaceCharge3D, phi0, radius0, z0);
        }
      } // end j
    } // end i
  } // end phi

  if (side == 0) {
    fInterpolatorChargeA->SetValue(fMatrixChargeA);
    fInterpolatorChargeA->InitCubicSpline();
  } else {
    fInterpolatorChargeC->SetValue(fMatrixChargeC);
    fInterpolatorChargeC->InitCubicSpline();
  }

  fInitLookUp = kFALSE;
}

/// create the distortion tree on a mesh with granularity given by step
/// return the tree with distortions at given position
/// Map is created on the mesh with given step size
/// type - 0: Call GetDistortion()
///        1: Call GetDistortionIntegralDz()
///
/// \param step
/// \return
TTree *AliTPCSpaceCharge3DDriftLine::CreateDistortionTree(Double_t step) {
  TTreeSRedirector *pcStream = new TTreeSRedirector(Form("distortion%s.root", GetName()));
  Float_t xyz[3];     // current point
  Float_t dist[3];    // distortion
  Float_t localDist[3];    // distortion
  Float_t corr[3];    // correction
  Float_t xyzDist[3]; // distorted point
  Float_t xyzCorr[3]; // corrected point

  //AliMagF* mag= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  //if (!mag) AliError("Magnetic field - not initialized");

  Int_t roc;
  AliTPCParam *tpcParam = new AliTPCParamSR;
  Double_t r, phi, rDist, phiDist, drDist, drPhiDist, rCorr, phiCorr, drCorr, drPhiCorr;
  for (Double_t x = -250; x < 250; x += step) {
    for (Double_t y = -250; y < 250; y += step) {

      r = TMath::Sqrt(x * x + y * y);

      if (tpcParam->GetPadRowRadii(0, 0) > r || r > tpcParam->GetPadRowRadii(36, 95)) continue;
      
      phi = TMath::ATan2(y, x);

      for (Double_t z = -250; z < 250; z += step) {
        roc = (z > 0) ? 0 : 18;
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;

        GetDistortion(xyz, roc, dist);

        for (Int_t i = 0; i < 3; ++i) {
          xyzDist[i] = xyz[i] + dist[i];
        }

        xyz[0] = r;
        xyz[1] = phi;
        xyz[2] = z;

        GetLocalDistortionCylAC(xyz, roc, localDist);

        GetCorrection(xyzDist, roc, corr);

        for (Int_t i = 0; i < 3; ++i) {
          xyzCorr[i] = xyzDist[i] + corr[i];
        }

        // === r, rPhi + residuals for the distorted point =========================
        rDist = TMath::Sqrt(xyzDist[0] * xyzDist[0] + xyzDist[1] * xyzDist[1]);
        phiDist = TMath::ATan2(xyzDist[1], xyzDist[0]);
        //rDist = xyzDist[0];
        //phiDist = xyzDist[1];

        while ((phiDist - phi) > TMath::Pi()) phiDist -= TMath::TwoPi();
        while ((phiDist - phi) < -TMath::Pi()) phiDist += TMath::TwoPi();

        drDist = rDist - r;
        drPhiDist = (phiDist - phi) * r;

        // === r, rPhi + residuals for the corrected point =========================
        rCorr = TMath::Sqrt(xyzCorr[0] * xyzCorr[0] + xyzCorr[1] * xyzCorr[1]);
        phiCorr = TMath::ATan2(xyzCorr[1], xyzCorr[0]);
        //rCorr = xyzCorr[0];
        //phiCorr = xyzCorr[1];

        while ((phiCorr - phiDist) > TMath::Pi()) phiCorr -= TMath::TwoPi();
        while ((phiCorr - phiDist) < -TMath::Pi()) phiCorr += TMath::TwoPi();

        drCorr = rCorr - rDist;
        drPhiCorr = (phiCorr - phiDist) * r;
        (*pcStream) << "distortion" <<
                    "x=" << x <<           // original position
                    "y=" << y <<
                    "z=" << z <<
                    "r=" << r <<
                    "phi=" << phi <<
                    "xDist=" << xyzDist[0] <<      // distorted position
                    "yDist=" << xyzDist[1] <<
                    "zDist=" << xyzDist[2] <<
                    "rDist=" << rDist <<
                    "phiDist=" << phiDist <<
                    "dxDist=" << dist[0] <<     // distortion
                    "dyDist=" << dist[1] <<
                    "dzDist=" << dist[2] <<
                    "drDist=" << drDist <<
                    "drPhiDist=" << drPhiDist <<
                    "drLocalDist=" << localDist[0] <<
                    "drPhiLocalDist=" << localDist[1] <<
                    "dzLocalDist=" << localDist[2] <<
                    "xCorr=" << xyzCorr[0] <<      // corrected position
                    "yCorr=" << xyzCorr[1] <<
                    "zCorr=" << xyzCorr[2] <<
                    "rCorr=" << rCorr <<
                    "phiCorr=" << phiCorr <<
                    //
                    "dxCorr=" << corr[0] <<     // correction
                    "dyCorr=" << corr[1] <<
                    "dzCorr=" << corr[2] <<
                    "drCorr=" << drCorr <<
                    "drPhiCorr=" << drPhiCorr <<
                    "\n";
      }
    }
  }
  delete pcStream;
  TFile f(Form("distortion%s.root", GetName()));
  TTree *tree = (TTree *) f.Get("distortion");

  return tree;
}

/// InterpolationPhi is only used for reading from TH3F (since it is not cylindrical)
///
/// \param r
/// \param z
/// \return
Double_t AliTPCSpaceCharge3DDriftLine::InterpolatePhi
        (
                TH3 *h3,
                const Double_t phi,
                const Double_t r,
                const Double_t z
        ) {

  Int_t ubx = h3->GetXaxis()->FindBin(phi);
  if (phi < h3->GetXaxis()->GetBinCenter(ubx)) ubx -= 1;
  Int_t obx = ubx + 1;

  Int_t uby = h3->GetYaxis()->FindBin(r);
  if (r < h3->GetYaxis()->GetBinCenter(uby)) uby -= 1;
  Int_t oby = uby + 1;

  Int_t ubz = h3->GetZaxis()->FindBin(z);
  if (z < h3->GetZaxis()->GetBinCenter(ubz)) ubz -= 1;
  Int_t obz = ubz + 1;

  if (uby <= 0 || ubz <= 0 ||
      oby > h3->GetYaxis()->GetNbins() || obz > h3->GetZaxis()->GetNbins()) {
    return 0;
  }

  if (ubx <= 0) ubx = h3->GetXaxis()->GetNbins();

  if (obx > h3->GetXaxis()->GetNbins()) obx = 1;

  Double_t xw = h3->GetXaxis()->GetBinCenter(obx) - h3->GetXaxis()->GetBinCenter(ubx);
  Double_t yw = h3->GetYaxis()->GetBinCenter(oby) - h3->GetYaxis()->GetBinCenter(uby);
  Double_t zw = h3->GetZaxis()->GetBinCenter(obz) - h3->GetZaxis()->GetBinCenter(ubz);

  Double_t xd = (phi - h3->GetXaxis()->GetBinCenter(ubx)) / xw;
  Double_t yd = (r - h3->GetYaxis()->GetBinCenter(uby)) / yw;
  Double_t zd = (z - h3->GetZaxis()->GetBinCenter(ubz)) / zw;

  Double_t v[] = {h3->GetBinContent(ubx, uby, ubz), h3->GetBinContent(ubx, uby, obz),
                  h3->GetBinContent(ubx, oby, ubz), h3->GetBinContent(ubx, oby, obz),
                  h3->GetBinContent(obx, uby, ubz), h3->GetBinContent(obx, uby, obz),
                  h3->GetBinContent(obx, oby, ubz), h3->GetBinContent(obx, oby, obz)};

  Double_t i1 = v[0] * (1 - zd) + v[1] * zd;
  Double_t i2 = v[2] * (1 - zd) + v[3] * zd;
  Double_t j1 = v[4] * (1 - zd) + v[5] * zd;
  Double_t j2 = v[6] * (1 - zd) + v[7] * zd;
  Double_t w1 = i1 * (1 - yd) + i2 * yd;
  Double_t w2 = j1 * (1 - yd) + j2 * yd;
  Double_t result = w1 * (1 - xd) + w2 * xd;
  return result;

}

///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramSCInXY(Float_t z, Int_t nx, Int_t ny) {
  TH2F *h = CreateTH2F("spaceCharge", GetTitle(), "x [cm]", "y [cm]", "#rho_{sc} [C/m^{3}/e_{0}]",
                       nx, -250., 250., ny, -250., 250.);

  for (Int_t iy = 1; iy <= ny; ++iy) {
    Double_t yp = h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix = 1; ix <= nx; ++ix) {
      Double_t xp = h->GetXaxis()->GetBinCenter(ix);

      Float_t r = TMath::Sqrt(xp * xp + yp * yp);
      Float_t phi = TMath::ATan2(yp, xp);

      if (85. <= r && r <= 250.) {

        Float_t sc = GetSpaceChargeDensity(r, phi, z) / fgke0; // in [C/m^3/e0]
        h->SetBinContent(ix, iy, sc);
      } else {
        h->SetBinContent(ix, iy, 0.);
      }
    }
  }

  return h;
}

/// return a simple histogram containing the space charge distribution (input for the calculation)
///
/// \param phi
/// \param nz
/// \param nr
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistogramSCInZR(Float_t phi, Int_t nz, Int_t nr) {

  TH2F *h = CreateTH2F("spaceCharge ", GetTitle(), "z [cm]", "r [cm]", "#rho_{sc} [C/m^{3}/e_{0}]",
                       nz, -250., 250., nr, 85., 250.);
  for (Int_t ir = 1; ir <= nr; ++ir) {
    Float_t r = h->GetYaxis()->GetBinCenter(ir);
    for (Int_t iz = 1; iz <= nz; ++iz) {
      Float_t z = h->GetXaxis()->GetBinCenter(iz);
      if (85. <= r && r <= 250.) {
        Float_t sc = GetSpaceChargeDensity(r, phi, z) / fgke0; // in [C/m^3/e0]
        h->SetBinContent(iz, ir, sc);
      } else {
        h->SetBinContent(iz, ir, 0.);
      }
    }
  }
  return h;
}

/// returns the (input) space charge density at a given point according
/// Note: input in [cm], output in [C/m^3/e0] !!
Float_t AliTPCSpaceCharge3DDriftLine::GetSpaceChargeDensity(Float_t r, Float_t phi, Float_t z) {
  while (phi < 0) phi += TMath::TwoPi();
  while (phi > TMath::TwoPi()) phi -= TMath::TwoPi();
  
  const Int_t order = 1; //

  const Float_t x[] = {r, phi, z};
  Float_t sc = 0;
  if (z > 0)
    sc = GetChargeCylAC(x, 0);
  else
    sc = GetChargeCylAC(x, 18);

  return sc;
}

///
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDPhiRDz TMatrixD **  matrix of global distortion dRPhi (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rList Double_t * points of r in the grid (ascending mode)
/// \param zList Double_t * points of z in the grid (ascending mode)
/// \param phiList Double_t * points of phi in the grid (ascending mode)
/// \param nRRow Int_t number of grid in r direction
/// \param nZColumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///	\param nStep Int_t number of step to calculate local dist
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionGlobalInvTable(
        TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz, TMatrixD **matricesDistDz, Double_t *rList,
        Double_t *zList,  Double_t *phiList, const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice,
        const Int_t nStep,  const Bool_t useCylAC,    Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t type) {
  Double_t z, phi, r, zAfter, zPrevious,  ddR, ddRPhi, ddZ, zl, dr, dRPhi, dz, ddPhi, dPhi, deltaZ, r0, z0, phi0;
  Float_t x[3], dx[3], pdx[3];

  Int_t roc;

  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;

  TMatrixD *distDrDz;
  TMatrixD *distDPhiRDz;
  TMatrixD *distDz;

  // correction build up for inverse flow
  TMatrixD *corrDrDz;
  TMatrixD *corrDPhiRDz;
  TMatrixD *corrDz;

  TMatrixD *listR;
  TMatrixD *listPhi;
  TMatrixD *listZ;

  TMatrixD *matricesCorrDrDz[phiSlice];
  TMatrixD *matricesCorrDPhiRDz[phiSlice];
  TMatrixD *matricesCorrDz[phiSlice];

  TMatrixD *matricesRList[phiSlice];
  TMatrixD *matricesPhiList[phiSlice];
  TMatrixD *matricesZList[phiSlice];

  for (Int_t m = 0; m < phiSlice; m++) {
    matricesCorrDrDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDPhiRDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDz[m] = new TMatrixD(nRRow, nZColumn);

    matricesRList[m] = new TMatrixD(nRRow, nZColumn);
    matricesPhiList[m] = new TMatrixD(nRRow, nZColumn);
    matricesZList[m] = new TMatrixD(nRRow, nZColumn);
  }

  AliTPCLookUpTable3DInterpolatorIrregularD *lookupInverseCorr = new AliTPCLookUpTable3DInterpolatorIrregularD(
          nRRow, matricesCorrDrDz, matricesRList,  phiSlice,  matricesCorrDPhiRDz,
          matricesPhiList,   nZColumn,  matricesCorrDz, matricesZList,  2,
          stepR, stepZ, stepPhi, type  );
  
  lookupInverseCorr->SetKernelType(GetRBFKernelType());

  for (Int_t k = 0; k < phiSlice; k++) {
    distDrDz = matricesDistDrDz[k];
    distDPhiRDz = matricesDistDPhiRDz[k];
    distDz = matricesDistDz[k];

    listR = matricesRList[k];
    listPhi = matricesPhiList[k];
    listZ = matricesZList[k];

    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDPhiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

      for (Int_t j = 0; j < nZColumn; j++) {
        (*listR)(i, j) = rList[i];
        (*listPhi)(i, j) = phiList[k];
        (*listZ)(i, j) = zList[j];

      }

    }
  }



  // 1) create global correction
  deltaZ = (zList[1] - zList[0]);
  Int_t iAnchor, kAnchor, zAnchor;

  for (Int_t j = nZColumn - 2; j >= 0; j--) {

    roc = 0; // FIXME
    for (Int_t k = 0; k < phiSlice; k++) {

      corrDrDz = matricesCorrDrDz[k];
      corrDPhiRDz = matricesCorrDPhiRDz[k];
      corrDz = matricesCorrDz[k];

      listR = matricesRList[k];
      listPhi = matricesPhiList[k];
      listZ = matricesZList[k];

      for (Int_t i = 0; i < nRRow; i++) {
        // get global distortion

        r = rList[i];
        phi = phiList[k];
        z = zList[j];

        dr = 0.0;
        dz = 0.0;
        dRPhi = 0.0;


        x[0] = r;
        x[1] = phi;
        x[2] = z;

        if (useCylAC == kTRUE)
          GetDistortionCylAC(x, roc, dx);
        else
          GetDistortionCyl(x, roc, dx);

        dr = dx[0];
        dz = dx[2];
        dRPhi = dx[1];

        r = rList[i];
        phi = phiList[k];
        z = zList[j];

        (*corrDrDz)(i, j + 1) = -dr;
        (*corrDz)(i, j + 1) = -dz;
        (*corrDPhiRDz)(i, j + 1) = -dRPhi;

        (*listR)(i, j + 1) = r + dr;
        (*listPhi)(i, j + 1) = phi + dRPhi / r;
        (*listZ)(i, j + 1) = z + dz;
      }
    }
    lookupInverseCorr->CopyFromMatricesToInterpolator(j + 1);
  }
  // 2) calculate local distortion
  for (Int_t j = nZColumn - 2; j >= 0; j--) {
    roc = 0; // FIXME
    for (Int_t k = 0; k < phiSlice; k++) {
      distDrDz = matricesDistDrDz[k];
      distDPhiRDz = matricesDistDPhiRDz[k];
      distDz = matricesDistDz[k];
      for (Int_t i = 0; i < nRRow; i++) {
        // get global distortion
        r = rList[i];
        phi = phiList[k];
        z = zList[j];
        dr = 0.0;
        dz = 0.0;
        dRPhi = 0.0;

        if (j < nZColumn - 2) {
          // get global distortion of this point
          x[0] = r;
          x[1] = phi;
          x[2] = z;
          if (useCylAC == kTRUE)
            GetDistortionCylAC(x, roc, dx);
          else
            GetDistortionCyl(x, roc, dx);

          r0 = r + dx[0];
          z0 = zList[j + 1] + dx[2];
          phi0 = phi + (dx[1] / r);
          iAnchor = TMath::FloorNint((r0 - fgkIFCRadius) / gridSizeR);
          kAnchor = TMath::FloorNint(phi0 / gridSizePhi);
          zAnchor = TMath::FloorNint(z0 / gridSizeZ);

          if (j > nZColumn - (GetIrregularGridSize() + 2))
            lookupInverseCorr->GetValue(r0, phi0, z0, dr, dRPhi, dz, iAnchor, kAnchor, zAnchor,
                                        nRRow / 4 + 1, phiSlice / 4 + 1, 1, 0);
          else
            lookupInverseCorr->GetValue(r0, phi0, z0, dr, dRPhi, dz, iAnchor, kAnchor, zAnchor,
                                        nRRow / 4 + 1, phiSlice / 4 + 1, GetIrregularGridSize(), 0);

          phi0 = phi0 + ((dRPhi) / r0);
          r0 = r0 + (dr);
          z0 += dz;

          x[0] = r0;
          x[1] = phi0;
          x[2] = z0;

          if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
          if (phi0 > TMath::TwoPi()) phi0 = phi0 - TMath::TwoPi();

          if (useCylAC == kTRUE)
            GetDistortionCylAC(x, roc, pdx);
          else
            GetDistortionCyl(x, roc, pdx);

          dr = (dx[0] - pdx[0]);
          dz = (dx[2] - pdx[2]);
          dRPhi = (dx[1] - pdx[1]);

        } else if (j == (nZColumn - 2)) {

          x[0] = r;
          x[1] = phi;
          x[2] = zList[j];
          if (useCylAC == kTRUE)
            GetDistortionCylAC(x, roc, dx);
          else
            GetDistortionCyl(x, roc, dx);

          x[2] = zList[j + 1];
          if (useCylAC == kTRUE)
            GetDistortionCylAC(x, roc, pdx);
          else
            GetDistortionCyl(x, roc, pdx);
          dr = (dx[0] - pdx[0]);
          dz = (dx[2] - pdx[2]);
          dRPhi = (dx[1] - pdx[1]);
        }

        (*distDrDz)(i, j) = dr;
        (*distDz)(i, j) = dz;
        (*distDPhiRDz)(i, j) = dRPhi;
      }
    }
  }

  for (Int_t m = 0; m < phiSlice; m++) {
    delete matricesCorrDrDz[m];
    delete matricesCorrDPhiRDz[m];
    delete matricesCorrDz[m];
    delete matricesRList[m];
    delete matricesPhiList[m];
    delete matricesZList[m];
  }
  delete lookupInverseCorr;
}

///
/// \param matricesEr
/// \param matricesEPhi
/// \param matricesEz
/// \param matricesInvLocalIntErDz
/// \param matricesInvLocalIntEPhiDz
/// \param matricesInvLocalIntEz
/// \param matricesDistDrDz
/// \param matricesDistDPhiRDz
/// \param matricesDistDz
/// \param rList
/// \param zList
/// \param phiList
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseLocalDistortionToElectricField (
        TMatrixD **matricesEr, TMatrixD **matricesEPhi, TMatrixD **matricesEz,
        TMatrixD **matricesInvLocalIntErDz, TMatrixD **matricesInvLocalIntEPhiDz,
        TMatrixD **matricesInvLocalIntEz, TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz, 
        TMatrixD **matricesDistDz,  Double_t *rList, Double_t *zList, Double_t *phiList,
        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice ) {
  // calculate integral
  Float_t localIntErOverEz, localIntEPhiOverEz, localIntDeltaEz, z2;
  Double_t r;
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;

  TMatrixD *distDrDz;
  TMatrixD *distDz;
  TMatrixD *distDPhiRDz;
  TMatrixD *tDistDz;
  TMatrixD *tDistDPhiRDz;
  TMatrixD *tDistDrDz;

  Float_t c02c12 = fC0 * fC0 + fC1 * fC1;

  // solve local integration
  for (Int_t j = 0; j < nZColumn; j++) {
    for (Int_t k = 0; k < phiSlice; k++) {
      distDrDz = matricesDistDrDz[k];
      distDz = matricesDistDz[k];
      distDPhiRDz = matricesDistDPhiRDz[k];

      tDistDrDz = matricesInvLocalIntErDz[k];
      tDistDz = matricesInvLocalIntEz[k];
      tDistDPhiRDz = matricesInvLocalIntEPhiDz[k];
      for (Int_t i = 0; i < nRRow; i++) {
        localIntErOverEz = fC0 * (*distDrDz)(i, j) - fC1 * (*distDPhiRDz)(i, j);
        localIntErOverEz = localIntErOverEz / (fC0 * fC0 + fC1 * fC1);
        localIntEPhiOverEz = ((*distDrDz)(i, j) - (fC0 * localIntErOverEz)) / fC1;
        localIntDeltaEz = (*distDz)(i, j) / (fgkdvdE * fgkdvdE); // two times?
        (*tDistDrDz)(i, j) = localIntErOverEz;
        (*tDistDPhiRDz)(i, j) = localIntEPhiOverEz;
        (*tDistDz)(i, j) = localIntDeltaEz;
      }
    }
  }
  TMatrixD *mEPhi;
  TMatrixD *mEr;
  TMatrixD *mEz;

  // use central-backward-forward difference for calculating Electric field component
  for (Int_t m = 0; m < phiSlice; m++) {
    mEPhi = matricesEPhi[m];
    mEr = matricesEr[m];
    mEz = matricesEz[m];
    distDrDz = matricesInvLocalIntErDz[m];
    distDPhiRDz = matricesInvLocalIntEPhiDz[m];
    distDz = matricesInvLocalIntEz[m];
    for (Int_t i = 0; i < nRRow; i++) {
      (*mEr)(i, 0) = ((*distDrDz)(i, 0) / gridSizeZ) * -1 * ezField;
      (*mEPhi)(i, 0) = ((*distDPhiRDz)(i, 0) / gridSizeZ) * -1 * ezField;
      (*mEz)(i, 0) = ((*distDz)(i, 0) / gridSizeZ);
      (*mEr)(i, nZColumn - 1) =
              ((-0.5 * (*distDrDz)(i, nZColumn - 3) + 1.5 * (*distDrDz)(i, nZColumn - 2)) / gridSizeZ) * -1 *
              ezField;
      (*mEPhi)(i, nZColumn - 1) =
              ((-0.5 * (*distDPhiRDz)(i, nZColumn - 3) + 1.5 * (*distDPhiRDz)(i, nZColumn - 2)) / gridSizeZ) * -1 *
              ezField;
      (*mEz)(i, nZColumn - 1) = (-0.5 * (*distDz)(i, nZColumn - 3) + 1.5 * (*distDz)(i, nZColumn - 2)) / gridSizeZ;
    }

    for (Int_t i = 0; i < nRRow; i++) {
      for (Int_t j = 1; j < nZColumn - 1; j++) {
        (*mEr)(i, j) = (((*distDrDz)(i, j) + (*distDrDz)(i, j - 1)) / (2 * gridSizeZ)) * -1 *
                       ezField; // z direction
        (*mEPhi)(i, j) = (((*distDPhiRDz)(i, j) + (*distDPhiRDz)(i, j - 1)) / (2 * gridSizeZ)) * -1 *
                         ezField; // z direction
        (*mEz)(i, j) = ((*distDz)(i, j) + (*distDz)(i, j - 1)) / (2 * gridSizeZ); // z direction
      }
    }
  }
}

/// Inverse Electric Field to Charge
/// using partial differential
///
/// \param matricesCharge
/// \param matricesEr
/// \param matricesEPhi
/// \param matricesEz
/// \param rList
/// \param zList
/// \param phiList
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseElectricFieldToCharge (
        TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEPhi, TMatrixD **matricesEz, 
        Double_t *rList, Double_t *zList, Double_t *phiList, const Int_t nRRow, 
        const Int_t nZColumn, const Int_t phiSlice) {

  Float_t radius;
  Double_t dr, dz, dPhi;

  Int_t mPlus, mMinus, mPlus2, mMinus2, signPlus, signMinus;
  Int_t symmetry = 0;

  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;

  for (Int_t m = 0; m < phiSlice; m++) {
    mPlus = m + 1;
    signPlus = 1;
    mMinus = m - 1;
    signMinus = 1;
    mPlus2 = m + 2;
    mMinus2 = m - 2;
    if (symmetry == 1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (mPlus > phiSlice - 1) mPlus = phiSlice - 2;
      if (mMinus < 0) mMinus = 1;
    } else if (symmetry == -1) {       // Anti-symmetry in phi
      if (mPlus > phiSlice - 1) {
        mPlus = phiSlice - 2;
        signPlus = -1;
      }
      if (mMinus < 0) {
        mMinus = 1;
        signMinus = -1;
      }
    } else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if (mPlus > phiSlice - 1) mPlus = m + 1 - phiSlice;
      if (mMinus < 0) mMinus = m - 1 + phiSlice;
      if (mPlus2 > phiSlice - 1) mPlus2 = m + 2 - phiSlice;
      if (mMinus2 < 0) mMinus2 = m - 2 + phiSlice;
    }

    TMatrixD &matrixCharge = *matricesCharge[m];
    TMatrixD &matrixEr = *matricesEr[m];
    TMatrixD &matrixEz = *matricesEz[m];
    TMatrixD &matrixEPhi = *matricesEPhi[m];
    TMatrixD &matrixEPhiM = *matricesEPhi[mMinus];
    TMatrixD &matrixEPhiP = *matricesEPhi[mPlus];
    TMatrixD &matrixEPhiM2 = *matricesEPhi[mMinus2];
    TMatrixD &matrixEPhiP2 = *matricesEPhi[mPlus2];


    // for non-boundary V
    for (Int_t i = 2; i < nRRow - 2; i++) {
      radius = fgkIFCRadius + i * gridSizeR;
      for (Int_t j = 2; j < nZColumn - 2; j++) {
        dr = (-matrixEr(i + 2, j) + 8 * matrixEr(i + 1, j) - 8 * matrixEr(i - 1, j) + matrixEr(i - 2, j)) /
             (12 * gridSizeR); // r direction
        dz = (-matrixEz(i, j + 2) + 8 * matrixEz(i, j + 1) - 8 * matrixEz(i, j - 1) + matrixEz(i, j - 2)) /
             (12 * gridSizeZ); // r direction
        dPhi = (-matrixEPhiP2(i, j) + 8 * matrixEPhiP(i, j) - 8 * matrixEPhiM(i, j) + matrixEPhiM2(i, j)) /
               (12 * gridSizePhi); // phi

        matrixCharge(i, j) = -1 * (matrixEr(i, j) / radius + dr + dPhi / radius + dz);
      }
    }

    // for boundary in r
    for (Int_t j = 2; j < nZColumn - 2; j++) {

      // r near inner radius
      // for index r[0]
      radius = fgkIFCRadius;
      dr = (-(11.0 / 6.0) * matrixEr(0, j) + (3.0 * matrixEr(1, j)) - (1.5 * matrixEr(2, j)) +
            ((1.0 / 3.0) * matrixEr(3, j))) / gridSizeR; // forward difference

      //	dr 	=  ( -(1.5)*matrixEr(0,j) + (2.0*matrixEr(1,j)) - (0.5*matrixEr(2,j)) )  / gridSizeR;

      dz = (-matrixEz(0, j + 2) + 8 * matrixEz(0, j + 1) - 8 * matrixEz(0, j - 1) + matrixEz(0, j - 2)) /
           (12.0 * gridSizeZ);; // z direction
      dPhi = (-matrixEPhiP2(0, j) + 8 * matrixEPhiP(0, j) - 8 * matrixEPhiM(0, j) + matrixEPhiM2(0, j)) /
             (12.0 * gridSizePhi);

      matrixCharge(0, j) = -1 * (matrixEr(0, j) / radius + dr + dPhi / radius + dz);


      // index use central difference 3-point center
      radius = fgkIFCRadius + gridSizeR;
      //	dr 	=  (-matrixEr(3,j)  +6.0*matrixEr(2,j) - 3.0*matrixEr(1,j) - 2*matrixEr(0,j) ) / (6.0*gridSizeR) ; // forward difference
      dr = (matrixEr(2, j) - matrixEr(0, j)) / (2.0 * gridSizeR);

      dz = (-matrixEz(1, j + 2) + 8 * matrixEz(1, j + 1) - 8 * matrixEz(1, j - 1) + matrixEz(1, j - 2)) /
           (12 * gridSizeZ);   // z direction
      dPhi = (-matrixEPhiP2(1, j) + 8 * matrixEPhiP(1, j) - 8 * matrixEPhiM(1, j) + matrixEPhiM2(1, j)) /
             (12 * gridSizePhi);
      matrixCharge(1, j) = -1 * (matrixEr(1, j) / radius + dr + dPhi / radius + dz);


      // index use central difference 3-point center
      radius = fgkIFCRadius + (nRRow - 2) * gridSizeR;
      //	dr =   (2.0 * matrixEr(nRRow - 1,j)  + 3.0*matrixEr(nRRow - 2,j) - 6.0*matrixEr(nRRow -3,j) + matrixEr(nRRow-4,j) ) / (6.0*gridSizeR) ;
      dr = (matrixEr(nRRow - 1, j) - matrixEr(nRRow - 3, j)) / (2.0 * gridSizeR);

      dz = (-matrixEz(nRRow - 2, j + 2) + 8 * matrixEz(nRRow - 2, j + 1) - 8 * matrixEz(nRRow - 2, j - 1) +
            matrixEz(nRRow - 2, j - 2)) / (12 * gridSizeZ);
      dPhi = (-matrixEPhiP2(nRRow - 2, j) + 8 * matrixEPhiP(nRRow - 2, j) - 8 * matrixEPhiM(nRRow - 2, j) +
              matrixEPhiM2(nRRow - 2, j)) / (12.0 * gridSizePhi);
      matrixCharge(nRRow - 2, j) = -1 * (matrixEr(nRRow - 2, j) / radius + dr + dPhi / radius + dz);

      // index r[nRRow -1] backward difference
      radius = fgkIFCRadius + (nRRow - 1) * gridSizeR;
      //dr =  ( 1.5*matrixEr(nRRow-1,j) - 2.0*matrixEr(nRRow-2,j) + 0.5*matrixEr(nRRow-3,j) ) / gridSizeR ; // backward difference
      dr = (-(11.0 / 6.0) * matrixEr(nRRow - 1, j) + (3.0 * matrixEr(nRRow - 2, j)) - (1.5 * matrixEr(nRRow - 3, j)) +
            ((1.0 / 3.0) * matrixEr(nRRow - 4, j))) / (-1 * gridSizeR);

      //dz    =  ( matrixEz(nRRow-1,j+1) - matrixEz(nRRow-1,j-1) ) / (2*gridSizeZ) ; // z direction
      dz = (-matrixEz(nRRow - 1, j + 2) + 8 * matrixEz(nRRow - 1, j + 1) - 8 * matrixEz(nRRow - 1, j - 1) +
            matrixEz(nRRow - 1, j - 2)) / (12 * gridSizeZ);

      dPhi = (-matrixEPhiP2(nRRow - 1, j) + 8 * matrixEPhiP(nRRow - 1, j) - 8 * matrixEPhiM(nRRow - 1, j) +
              matrixEPhiM2(nRRow - 1, j)) / (12 * gridSizePhi);
      matrixCharge(nRRow - 1, j) = -1 * (matrixEr(nRRow - 1, j) / radius + dr + dPhi / radius + dz);
    }

    // boundary z
    for (Int_t i = 2; i < nRRow - 2; i++) {
      // z[0]
      radius = fgkIFCRadius + i * gridSizeR;
      dz = (-(11.0 / 6.0) * matrixEz(i, 0) + (3.0 * matrixEz(i, 1)) - (1.5 * matrixEz(i, 2)) +
            ((1.0 / 3.0) * matrixEz(i, 3))) / (1 * gridSizeZ); // forward difference
      dr = (-matrixEr(i + 2, 0) + 8 * matrixEr(i + 1, 0) - 8 * matrixEr(i - 1, 0) + matrixEr(i - 2, 0)) /
           (12 * gridSizeR);; // z direction
      dPhi = (-matrixEPhiP2(i, 0) + 8 * matrixEPhiP(i, 0) - 8 * matrixEPhiM(i, 0) + matrixEPhiM2(i, 0)) /
             (12 * gridSizePhi);
      matrixCharge(i, 0) = -1 * (matrixEr(i, 0) / radius + dr + dPhi / radius + dz);

      dz = (matrixEz(i, 2) - matrixEz(i, 0)) / (2.0 * gridSizeZ); // forward difference

      dr = (-matrixEr(i + 2, 1) + 8 * matrixEr(i + 1, 1) - 8 * matrixEr(i - 1, 1) + matrixEr(i - 2, 1)) /
           (12 * gridSizeR);; // z direction
      dPhi = (-matrixEPhiP2(i, 1) + 8 * matrixEPhiP(i, 1) - 8 * matrixEPhiM(i, 1) + matrixEPhiM2(i, 1)) /
             (12 * gridSizePhi);
      matrixCharge(i, 1) = -1 * (matrixEr(i, 1) / radius + dr + dPhi / radius + dz);

      dz = (matrixEz(i, nZColumn - 1) - matrixEz(i, nZColumn - 3)) / (2.0 * gridSizeZ); // forward difference

      dr = (-matrixEr(i + 2, nZColumn - 2) + 8 * matrixEr(i + 1, nZColumn - 2) - 8 * matrixEr(i - 1, nZColumn - 2) +
            matrixEr(i - 2, nZColumn - 2)) / (12 * gridSizeR);; // z direction
      dPhi = (-matrixEPhiP2(i, nZColumn - 2) + 8 * matrixEPhiP(i, nZColumn - 2) - 8 * matrixEPhiM(i, nZColumn - 2) +
              matrixEPhiM2(i, nZColumn - 2)) / (12 * gridSizePhi);
      matrixCharge(i, nZColumn - 2) = -1 * (matrixEr(i, nZColumn - 2) / radius + dr + dPhi / radius + dz);

      dz = (-(11.0 / 6.0) * matrixEz(i, nZColumn - 1) + (3.0 * matrixEz(i, nZColumn - 2)) -
            (1.5 * matrixEz(i, nZColumn - 3)) + ((1.0 / 3.0) * matrixEz(i, nZColumn - 4))) /
           (-gridSizeZ); // backward difference
      dr = (-matrixEr(i + 2, nZColumn - 1) + 8 * matrixEr(i + 1, nZColumn - 1) - 8 * matrixEr(i - 1, nZColumn - 1) +
            matrixEr(i - 2, nZColumn - 1)) / (12 * gridSizeR);; // z direction
      dPhi = (-matrixEPhiP2(i, nZColumn - 1) + 8 * matrixEPhiP(i, nZColumn - 1) - 8 * matrixEPhiM(i, nZColumn - 1) +
              matrixEPhiM2(i, nZColumn - 1)) / (12 * gridSizePhi);
      matrixCharge(i, nZColumn - 1) = -1 * (matrixEr(i, nZColumn - 1) / radius + dr + dPhi / radius + dz);
    }
    // for corner points
    // corner points for EPhi
    radius = fgkIFCRadius;
    dr = (-0.5 * matrixEr(2, 0) + 2.0 * matrixEr(1, 0) - 1.5 * matrixEr(0, 0)) / gridSizeR; // forward difference
    dz = (-0.5 * matrixEz(0, 2) + 2.0 * matrixEz(0, 1) - 1.5 * matrixEz(0, 0)) / gridSizeZ; // forward difference
    dPhi = (-matrixEPhiP2(0, 0) + 8 * matrixEPhiP(0, 0) - 8 * matrixEPhiM(0, 0) + matrixEPhiM2(0, 0)) /
           (12 * gridSizePhi);
    matrixCharge(0, 0) = -1 * (matrixEr(0, 0) / radius + dr + dPhi / radius + dz);
    dr = (-0.5 * matrixEr(2, 1) + 2.0 * matrixEr(1, 1) - 1.5 * matrixEr(0, 1)) / gridSizeR; // forward difference
    dz = (matrixEz(0, 2) - matrixEz(0, 0)) / (2.0 * gridSizeZ); // forward difference
    dPhi = (-matrixEPhiP2(0, 1) + 8 * matrixEPhiP(0, 1) - 8 * matrixEPhiM(0, 1) + matrixEPhiM2(0, 1)) /
           (12 * gridSizePhi);
    matrixCharge(0, 1) = -1 * (matrixEr(0, 1) / radius + dr + dPhi / radius + dz);
    dr = (-0.5 * matrixEr(2, nZColumn - 2) + 2.0 * matrixEr(1, nZColumn - 2) - 1.5 * matrixEr(0, nZColumn - 2)) /
         gridSizeR; // forward difference
    dz = (2.0 * matrixEz(0, nZColumn - 1) + 3.0 * matrixEz(0, nZColumn - 2) - 6.0 * matrixEz(0, nZColumn - 3) +
          matrixEz(0, nZColumn - 4)) / (6.0 * gridSizeZ); // backward difference
    dPhi = (-matrixEPhiP2(0, nZColumn - 2) + 8 * matrixEPhiP(0, nZColumn - 2) - 8 * matrixEPhiM(0, nZColumn - 2) +
            matrixEPhiM2(0, nZColumn - 2)) / (12 * gridSizePhi);
    matrixCharge(0, nZColumn - 2) = -1 * (matrixEr(0, nZColumn - 2) / radius + dr + dPhi / radius + dz);
    dr = (-0.5 * matrixEr(2, nZColumn - 1) + 2.0 * matrixEr(1, nZColumn - 1) - 1.5 * matrixEr(0, nZColumn - 1)) /
         gridSizeR; // forward difference
    dz = (1.5 * matrixEz(0, nZColumn - 1) - 2.0 * matrixEz(0, nZColumn - 2) + 0.5 * matrixEz(0, nZColumn - 3)) /
         gridSizeZ; // backward difference
    dPhi = (-matrixEPhiP2(0, nZColumn - 1) + 8 * matrixEPhiP(0, nZColumn - 1) - 8 * matrixEPhiM(0, nZColumn - 1) +
            matrixEPhiM2(0, nZColumn - 1)) / (12 * gridSizePhi);
    matrixCharge(0, nZColumn - 1) = -1 * (matrixEr(0, nZColumn - 1) / radius + dr + dPhi / radius + dz);

    radius = fgkIFCRadius + gridSizeR;
    dr = (-matrixEr(3, 0) + 6.0 * matrixEr(2, 0) - 3.0 * matrixEr(1, 0) - 2 * matrixEr(0, 0)) /
         (6.0 * gridSizeR); // forward difference
    dz = (-0.5 * matrixEz(1, 2) + 2.0 * matrixEz(1, 1) - 1.5 * matrixEz(1, 0)) / gridSizeZ; // forward difference
    dPhi = (-matrixEPhiP2(1, 0) + 8 * matrixEPhiP(1, 0) - 8 * matrixEPhiM(1, 0) + matrixEPhiM2(1, 0)) /
           (12 * gridSizePhi);
    matrixCharge(1, 0) = -1 * (matrixEr(1, 0) / radius + dr + dPhi / radius + dz);
    dr = (-matrixEr(3, 1) + 6.0 * matrixEr(2, 1) - 3.0 * matrixEr(1, 1) - 2 * matrixEr(0, 1)) /
         (6.0 * gridSizeR); // forward difference
    dz = (-matrixEz(1, 3) + 6.0 * matrixEz(1, 2) - 3.0 * matrixEz(1, 1) - 2 * matrixEz(1, 0)) /
         (6.0 * gridSizeZ); // forward difference
    dPhi = (-matrixEPhiP2(1, 1) + 8 * matrixEPhiP(1, 1) - 8 * matrixEPhiM(1, 1) + matrixEPhiM2(1, 1)) /
           (12 * gridSizePhi);
    matrixCharge(1, 1) = -1 * (matrixEr(1, 1) / radius + dr + dPhi / radius + dz);
    dr = (-matrixEr(3, nZColumn - 2) + 6.0 * matrixEr(2, nZColumn - 2) - 3.0 * matrixEr(1, nZColumn - 2) -
          2 * matrixEr(0, nZColumn - 2)) / (6.0 * gridSizeR); // forward difference
    dz = (2.0 * matrixEz(1, nZColumn - 1) + 3.0 * matrixEz(1, nZColumn - 2) - 6.0 * matrixEz(1, nZColumn - 3) +
          matrixEz(1, nZColumn - 4)) / (6.0 * gridSizeZ); // backward difference
    dPhi = (-matrixEPhiP2(1, nZColumn - 2) + 8 * matrixEPhiP(1, nZColumn - 2) - 8 * matrixEPhiM(1, nZColumn - 2) +
            matrixEPhiM2(1, nZColumn - 2)) / (12 * gridSizePhi);
    matrixCharge(1, nZColumn - 2) = -1 * (matrixEr(1, nZColumn - 2) / radius + dr + dPhi / radius + dz);

    dr = (-matrixEr(3, nZColumn - 1) + 6.0 * matrixEr(2, nZColumn - 1) - 3.0 * matrixEr(1, nZColumn - 1) -
          2 * matrixEr(0, nZColumn - 1)) / (6.0 * gridSizeR); // forward difference
    dz = (1.5 * matrixEz(1, nZColumn - 1) - 2.0 * matrixEz(1, nZColumn - 2) + 0.5 * matrixEz(1, nZColumn - 3)) /
         gridSizeZ; // backward difference
    dPhi = (-matrixEPhiP2(1, nZColumn - 1) + 8 * matrixEPhiP(1, nZColumn - 1) - 8 * matrixEPhiM(1, nZColumn - 1) +
            matrixEPhiM2(1, nZColumn - 1)) / (12 * gridSizePhi);
    matrixCharge(1, nZColumn - 1) = -1 * (matrixEr(1, nZColumn - 1) / radius + dr + dPhi / radius + dz);

    radius = fgkIFCRadius + (nRRow - 2) * gridSizeR;
    dr = (2.0 * matrixEr(nRRow - 1, 0) + 3.0 * matrixEr(nRRow - 2, 0) - 6.0 * matrixEr(nRRow - 3, 0) +
          matrixEr(nRRow - 4, 0)) / (6.0 * gridSizeR); // backward difference
    dz = (-0.5 * matrixEz(nRRow - 2, 2) + 2.0 * matrixEz(nRRow - 2, 1) - 1.5 * matrixEz(nRRow - 2, 0)) /
         gridSizeZ; // forward difference
    dPhi = (-matrixEPhiP2(nRRow - 2, 0) + 8 * matrixEPhiP(nRRow - 2, 0) - 8 * matrixEPhiM(nRRow - 2, 0) +
            matrixEPhiM2(nRRow - 2, 0)) / (12 * gridSizePhi);

    matrixCharge(nRRow - 2, 0) = -1 * (matrixEr(nRRow - 2, 0) / radius + dr + dPhi / radius + dz);
    dr = (2.0 * matrixEr(nRRow - 1, 1) + 3.0 * matrixEr(nRRow - 2, 1) - 6.0 * matrixEr(nRRow - 3, 1) +
          matrixEr(nRRow - 4, 1)) / (6.0 * gridSizeR); // backward difference
    dz = (-matrixEz(nRRow - 2, 3) + 6.0 * matrixEz(nRRow - 2, 2) - 3.0 * matrixEz(nRRow - 2, 1) -
          2 * matrixEz(nRRow - 2, 0)) / (6.0 * gridSizeZ); // forward difference
    dPhi = (-matrixEPhiP2(nRRow - 2, 1) + 8 * matrixEPhiP(nRRow - 2, 1) - 8 * matrixEPhiM(nRRow - 2, 1) +
            matrixEPhiM2(nRRow - 2, 1)) / (12 * gridSizePhi);
    matrixCharge(nRRow - 2, 1) = -1 * (matrixEr(nRRow - 2, 1) / radius + dr + dPhi / radius + dz);
    dr = (2.0 * matrixEr(nRRow - 1, nZColumn - 2) + 3.0 * matrixEr(nRRow - 2, nZColumn - 2) -
          6.0 * matrixEr(nRRow - 3, nZColumn - 2) + matrixEr(nRRow - 4, nZColumn - 2)) /
         (6.0 * gridSizeR); // backward difference
    dz = (2.0 * matrixEz(nRRow - 2, nZColumn - 1) + 3.0 * matrixEz(nRRow - 2, nZColumn - 2) -
          6.0 * matrixEz(nRRow - 2, nZColumn - 3) + matrixEz(nRRow - 2, nZColumn - 4)) /
         (6.0 * gridSizeZ); // backward difference
    dPhi = (-matrixEPhiP2(nRRow - 2, nZColumn - 2) + 8 * matrixEPhiP(nRRow - 2, nZColumn - 2) -
            8 * matrixEPhiM(nRRow - 2, nZColumn - 2) + matrixEPhiM2(nRRow - 2, nZColumn - 2)) / (12 * gridSizePhi);
    matrixCharge(nRRow - 2, nZColumn - 2) = -1 * (matrixEr(nRRow - 2, nZColumn - 2) / radius + dr + dPhi / radius + dz);
    dr = (2.0 * matrixEr(nRRow - 1, nZColumn - 1) + 3.0 * matrixEr(nRRow - 2, nZColumn - 1) -
          6.0 * matrixEr(nRRow - 3, nZColumn - 1) + matrixEr(nRRow - 4, nZColumn - 1)) /
         (6.0 * gridSizeR); // backward difference
    dz = (1.5 * matrixEz(0, nZColumn - 1) - 2.0 * matrixEz(0, nZColumn - 2) + 0.5 * matrixEz(0, nZColumn - 3)) /
         gridSizeZ; // backward difference
    dPhi = (-matrixEPhiP2(nRRow - 2, nZColumn - 1) + 8 * matrixEPhiP(nRRow - 2, nZColumn - 1) -
            8 * matrixEPhiM(nRRow - 2, nZColumn - 1) + matrixEPhiM2(nRRow - 2, nZColumn - 1)) / (12 * gridSizePhi);

    matrixCharge(nRRow - 2, nZColumn - 1) = -1 * (matrixEr(nRRow - 2, nZColumn - 1) / radius + dr + dPhi / radius + dz);
    radius = fgkIFCRadius + (nRRow - 1) * gridSizeR;
    dr = (1.5 * matrixEr(nRRow - 1, 0) - 2.0 * matrixEr(nRRow - 2, 0) + 0.5 * matrixEr(nRRow - 3, 0)) /
         gridSizeR; // backward difference
    dz = (-0.5 * matrixEz(nRRow - 1, 2) + 2.0 * matrixEz(nRRow - 1, 1) - 1.5 * matrixEz(nRRow - 1, 0)) /
         gridSizeZ; // forward difference
    dPhi = (-matrixEPhiP2(nRRow - 1, 0) + 8 * matrixEPhiP(nRRow - 1, 0) - 8 * matrixEPhiM(nRRow - 1, 0) +
            matrixEPhiM2(nRRow - 1, 0)) / (12 * gridSizePhi);
    matrixCharge(nRRow - 1, 0) = -1 * (matrixEr(nRRow - 1, 0) / radius + dr + dPhi / radius + dz);
    dr = (1.5 * matrixEr(nRRow - 1, 1) - 2.0 * matrixEr(nRRow - 2, 1) + 0.5 * matrixEr(nRRow - 3, 1)) /
         gridSizeR; // backward difference
    dz = (-matrixEz(nRRow - 1, 3) + 6.0 * matrixEz(nRRow - 1, 2) - 3.0 * matrixEz(nRRow - 1, 1) -
          2 * matrixEz(nRRow - 1, 0)) / (6.0 * gridSizeZ); // forward difference
    dPhi = (-matrixEPhiP2(nRRow - 1, 1) + 8 * matrixEPhiP(nRRow - 1, 1) - 8 * matrixEPhiM(nRRow - 1, 1) +
            matrixEPhiM2(nRRow - 1, 1)) / (12 * gridSizePhi);
    matrixCharge(nRRow - 1, 1) = -1 * (matrixEr(nRRow - 1, 1) / radius + dr + dPhi / radius + dz);

    dr = (1.5 * matrixEr(nRRow - 1, nZColumn - 2) - 2.0 * matrixEr(nRRow - 2, nZColumn - 2) +
          0.5 * matrixEr(nRRow - 3, nZColumn - 2)) / gridSizeR; // backward difference
    dz = (2.0 * matrixEz(nRRow - 1, nZColumn - 1) + 3.0 * matrixEz(nRRow - 1, nZColumn - 2) -
          6.0 * matrixEz(nRRow - 1, nZColumn - 3) + matrixEz(nRRow - 1, nZColumn - 4)) /
         (6.0 * gridSizeZ); // backward difference
    dPhi = (-matrixEPhiP2(nRRow - 1, nZColumn - 2) + 8 * matrixEPhiP(nRRow - 1, nZColumn - 2) -
            8 * matrixEPhiM(nRRow - 1, nZColumn - 2) + matrixEPhiM2(nRRow - 1, nZColumn - 2)) / (12 * gridSizePhi);
    matrixCharge(nRRow - 1, nZColumn - 2) = -1 * (matrixEr(nRRow - 1, nZColumn - 2) / radius + dr + dPhi / radius + dz);

    dr = (1.5 * matrixEr(nRRow - 1, nZColumn - 1) - 2.0 * matrixEr(nRRow - 2, nZColumn - 1) +
          0.5 * matrixEr(nRRow - 3, nZColumn - 1)) / gridSizeR; // backward difference
    dz = (1.5 * matrixEz(nRRow - 1, nZColumn - 1) - 2.0 * matrixEz(nRRow - 1, nZColumn - 2) +
          0.5 * matrixEz(nRRow - 1, nZColumn - 3)) / gridSizeZ; // backward difference

    dPhi = (-matrixEPhiP2(nRRow - 1, nZColumn - 1) + 8 * matrixEPhiP(nRRow - 1, nZColumn - 1) -
            8 * matrixEPhiM(nRRow - 1, nZColumn - 1) + matrixEPhiM2(nRRow - 1, nZColumn - 1)) / (12 * gridSizePhi);

    matrixCharge(nRRow - 1, nZColumn - 1) = -1 * (matrixEr(nRRow - 1, nZColumn - 1) / radius + dr + dPhi / radius + dz);
  }
}

///
/// \param matricesCharge
/// \param matricesEr
/// \param matricesEPhi
/// \param matricesEz
/// \param matricesInvLocalIntErDz
/// \param matricesInvLocalIntEPhiDz
/// \param matricesInvLocalEz
/// \param matricesDistDrDz
/// \param matricesDistDPhiRDz
/// \param matricesDistDz
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
/// \param nSize
/// \param useCylAC
/// \param stepR
/// \param stepZ
/// \param stepPhi
/// \param interpType
/// \param inverseType
void AliTPCSpaceCharge3DDriftLine::InverseDistortionMaps(
        TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEPhi, TMatrixD **matricesEz,
        TMatrixD **matricesInvLocalIntErDz, TMatrixD **matricesInvLocalIntEPhiDz, TMatrixD **matricesInvLocalEz,
        TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz, TMatrixD **matricesDistDz,
        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice, const Int_t nSize,
        const Bool_t useCylAC, Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t interpType) {
  // can inverse after lookup table for global distortion been calculated
  Double_t *rList = new Double_t[nRRow];
  Double_t *zList = new Double_t[nZColumn];
  Double_t *phiList = new Double_t[phiSlice];

  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;

  for (Int_t k = 0; k < phiSlice; k++) phiList[k] = gridSizePhi * k;
  for (Int_t i = 0; i < nRRow; i++) rList[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < nZColumn; j++) zList[j] = (j * gridSizeZ);
  // memory allocation
  if (fInitLookUp) {
    // 1)  get local distortion
    InverseGlobalToLocalDistortionGlobalInvTable(matricesDistDrDz, matricesDistDPhiRDz, matricesDistDz, rList,
                                                 zList, phiList, nRRow, nZColumn, phiSlice, nSize, useCylAC,
                                                 stepR, stepZ, stepPhi, interpType);

    fLookupInverseDistA->SetLookUpR(matricesDistDrDz);
    fLookupInverseDistA->SetLookUpPhi(matricesDistDPhiRDz);
    fLookupInverseDistA->SetLookUpZ(matricesDistDz);
    fLookupInverseDistA->CopyFromMatricesToInterpolator();

    // 2)  calculate local integral
    InverseLocalDistortionToElectricField(matricesEr, matricesEPhi, matricesEz, matricesInvLocalIntErDz,
                                          matricesInvLocalIntEPhiDz, matricesInvLocalEz,
                                          matricesDistDrDz, matricesDistDPhiRDz, matricesDistDz, rList, zList,
                                          phiList, nRRow, nZColumn, phiSlice);
    // 3)  get potential from electric field assuming zero boundaries
    InverseElectricFieldToCharge(matricesCharge, matricesEr, matricesEPhi, matricesEz, rList, zList, phiList, nRRow,
                                 nZColumn, phiSlice);
  }

  // copy charge inverse here just for side A (TODO: do for side C)
  for (Int_t k = 0; k < phiSlice; k++) *(fMatrixChargeInverseA[k]) = *(matricesCharge[k]);
  fInterpolatorInverseChargeA->SetValue(fMatrixChargeInverseA);
  fInterpolatorInverseChargeA->InitCubicSpline();

  delete[] zList;
  delete[] rList;
  delete[] phiList;
}

/// CalculateEField (New Version: with reorganization of modules)
/// Calculate E field based on look-up table created by Poisson Solver
/// * Differentiate V(r) and solve for E(r) using special equations for the first and last row
/// * Integrate E(r)/E(z) from point of origin to pad plane
/// * Differentiate V(r) and solve for E(phi)
/// * Integrate E(phi)/E(z) from point of origin to pad plane
/// * Differentiate V(r) and solve for E(z) using special equations for the first and last row
/// * Integrate (E(z)-Ez(ROC)) from point of origin to pad plane
///
/// \param matricesV TMatrixD** 3D matrix representing calculated potential
/// \param matricesErOverEz TMatrix** 3D matrix representing e-field at Er/Ez
/// \param matricesEPhiOverEz TMatrix** 3D matrix representing e-field at EPhi/Ez
/// \param matricesDeltaZ TMatrix** 3D matrix representing e-field at DeltaZ
/// \param nRRow Int_t number of nRRow (in R direction)
/// \param nZColumn Int_t number of nZColumn (in Z direction)
/// \param phiSlice Int_t number of (phi slices in phi direction)
/// \param symmetry Int_t symmetry?
/// \param rocDisplace rocDisplacement
///
/// \pre   Matrix matricesV is assumed had been calculated  by Poisson solver
/// \post  Results of Integration and Derivations for E-field calculation are stored in matricesErOverEz, matricesEPhiOverEz, matricesDeltaZ
///
void AliTPCSpaceCharge3DDriftLine::CalculateEField(
        TMatrixD **matricesV, TMatrixD **matricesErOverEz, TMatrixD **matricesEPhiOverEz,
        TMatrixD **matricesDeltaEz, const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice,
        const Int_t symmetry, Bool_t rocDisplacement) {

  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;

  TMatrixD *matricesEr[phiSlice], *matricesEz[phiSlice], *matricesEPhi[phiSlice];

  AliSysInfo::AddStamp("CalcField", 100, 0, 0);

  //Allocate memory for electric field r,z, phi direction
  for (Int_t k = 0; k < phiSlice; k++) {
    matricesEr[k] = new TMatrixD(nRRow, nZColumn);
    matricesEz[k] = new TMatrixD(nRRow, nZColumn);
    matricesEPhi[k] = new TMatrixD(nRRow, nZColumn);
  }

  //Differentiate V(r) and solve for E(r) using special equations for the first and last row
  TStopwatch w;
  w.Start();

  ElectricField(matricesV, matricesEr, matricesEPhi, matricesEz, nRRow, nZColumn,
                phiSlice, gridSizeR, gridSizePhi, gridSizeZ, symmetry, fgkIFCRadius);

  w.Stop();
  AliInfo(Form("Time for calculation E-field CPU = %f s\n", w.CpuTime()));

  AliSysInfo::AddStamp("Electron Drift Calc", 120, 0, 0);

  //Integrate E(r)/E(z) from point of origin to pad plane

  IntegrateEz(matricesErOverEz, matricesEr, nRRow, nZColumn, phiSlice, ezField);
  IntegrateEz(matricesEPhiOverEz, matricesEPhi, nRRow, nZColumn, phiSlice, ezField);
  IntegrateEz(matricesDeltaEz, matricesEz, nRRow, nZColumn, phiSlice, -1.0);

  // calculate z distortion from the integrated Delta Ez residuals
  // and include the equivalence (Volt to cm) of the ROC shift !!
  for (Int_t m = 0; m < phiSlice; m++) {
    TMatrixD &arrayV = *matricesV[m];
    TMatrixD &deltaEz = *matricesDeltaEz[m];

    for (Int_t j = 0; j < nZColumn; j++) {
      for (Int_t i = 0; i < nRRow; i++) {
        // Scale the Ez distortions with the drift velocity  -> delivers cm
        deltaEz(i, j) = deltaEz(i, j) * fgkdvdE;
        // ROC Potential in cm equivalent
        Double_t dzROCShift = arrayV(i, nZColumn - 1) / ezField;
        if (rocDisplacement) deltaEz(i, j) = deltaEz(i, j) + dzROCShift;  // add the ROC mis alignment
      }
    }
  }
  // clear the temporary arrays lists

  for (Int_t k = 0; k < phiSlice; k++) {
    delete matricesEr[k];
    delete matricesEz[k];
    delete matricesEPhi[k];
  }
}

///
/// Integrate at z direction Ez for electron drift calculation
///
///
/// \param matricesExOverEz TMatrixD** 3D matrix representing ExOverEz
/// \param matricesEx TMatrix** 3D matrix representing e-field at x direction
/// \param nRRow const Int_t number of nRRow  (in R direction)
/// \param nZColumn const Int_t number of nZColumn  (in Z direction)
/// \param phiSlice const Int_t number of (phiSlice in phi direction)
/// \param ezField const Double_t Electric field in z direction
///
/// \pre   matricesEx is assumed already been calculated by ElectricFieldCalculation
/// \post  Matrix matricesExOverEz is calculated by integration of matricesEx
///
void AliTPCSpaceCharge3DDriftLine::IntegrateEz (
        TMatrixD **matricesExOverEz, TMatrixD **matricesEx, const Int_t nRRow, const Int_t nZColumn,
        const Int_t phiSlice, const Double_t ezField) {
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  for (Int_t m = 0; m < phiSlice; m++) {
    TMatrixD &eXoverEz = *matricesExOverEz[m];
    TMatrixD &arrayEx = *matricesEx[m];

    for (Int_t j = nZColumn - 1; j >= 0; j--) {
      for (Int_t i = 0; i < nRRow; i++) {

        /// Calculate integration from int^{0}_{j} (TODO: Split the integration)
        if (j < nZColumn - 3) {
          eXoverEz(i, j) = eXoverEz(i, j + 2) +
                           (gridSizeZ / 3.0) * (arrayEx(i, j) + 4 * arrayEx(i, j + 1) + arrayEx(i, j + 2)) /
                           (-1 * ezField);
        } else {
          if (j == nZColumn - 3) {
            eXoverEz(i, j) = (gridSizeZ / 3.0) * (arrayEx(i, nZColumn - 3) + 4 * arrayEx(i, nZColumn - 2) +
                                                  arrayEx(i, nZColumn - 1)) / (-1 * ezField);
          }
          if (j == nZColumn - 2) {
            eXoverEz(i, j) =
                    (gridSizeZ / 3.0) * (1.5 * arrayEx(i, nZColumn - 2) + 1.5 * arrayEx(i, nZColumn - 1)) /
                    (-1 * ezField);
          }
          if (j == nZColumn - 1) {
            eXoverEz(i, j) = 0.0;
          }
        }
      }
    }
  }
}

/// GetCorrection from no-drift
///
/// \param x Float_t point origin
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylNoDrift(const Float_t x[], const Short_t roc, Float_t dx[]) {
  /// Calculates the correction due the Space Charge effect within the TPC drift volume

  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
//    InitSpaceCharge3DDistortion();
    return;
  }

  Float_t intEr, intEPhi, intDEz;
  Double_t r, phi, z;
  Int_t sign;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi -= TMath::TwoPi();

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  if (sign == -1 && z < 0.0) {
    printf("call C side\n");
    fLookupIntENoDriftC->GetValue(r, phi, z, intEr, intEPhi, intDEz);
  } else
    fLookupIntENoDriftA->GetValue(r, phi, z, intEr, intEPhi, intDEz);

  // Calculate distorted position
  if (r > 0.0) {
    phi = phi + fCorrectionFactor * (fC0 * intEPhi - fC1 * intEr) / r;
    r = r + fCorrectionFactor * (fC0 * intEr + fC1 * intEPhi);
  }
  Double_t dz = intDEz * fCorrectionFactor * fgkdvdE;

  // Calculate correction in cartesian coordinates
  dx[0] = -(r - x[0]);
  dx[1] = -(phi - x[1]);
  dx[2] = -dz;  // z distortion - (scaled with drift velocity dependency on the Ez field and the overall scaling factor)

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetDistortionCylNoDrift(const Float_t x[], Short_t roc, Float_t dx[]) {
  /// This function delivers the distortion values dx in respect to the initial coordinates x
  /// roc represents the TPC read out chamber (offline numbering convention)

  GetCorrectionCylNoDrift(x, roc, dx);
  for (Int_t j = 0; j < 3; ++j) dx[j] = -dx[j];
}

/// inverse for no drift
/// inverse from global distortion to local distortion
///
/// \param matricesDistDrDz
/// \param matricesDistDPhiRDz
/// \param matricesDistDz
/// \param rList
/// \param zList
/// \param phiList
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionNoDrift (
        TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz, TMatrixD **matricesDistDz,   
        Double_t *rList, Double_t *zList, Double_t *phiList,
        const Int_t nRRow, const Int_t nZColumn,  const Int_t phiSlice) {
  
  Double_t z, phi, r, zAfter, zPrevious,  ddR, ddRPhi, ddZ, dr, dRPhi, dz;
  Float_t x[3], dx[3], pdx[3], dxp1[3], dxp2[3];

  Int_t roc;
  
  TMatrixD *distDrDz;
  TMatrixD *distDPhiRDz;
  TMatrixD *distDz;

  for (Int_t k = 0; k < phiSlice; k++) {
    distDrDz = matricesDistDrDz[k];
    distDPhiRDz = matricesDistDPhiRDz[k];
    distDz = matricesDistDz[k];
    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDPhiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

    }
  }


  for (Int_t j = nZColumn - 2; j >= 0; j--) {
    roc = 0; // FIXME
    for (Int_t k = 0; k < phiSlice; k++) {

      distDrDz = matricesDistDrDz[k];
      distDPhiRDz = matricesDistDPhiRDz[k];
      distDz = matricesDistDz[k];
      for (Int_t i = 0; i < nRRow; i++) {
        // get global distortion

        r = rList[i];
        phi = phiList[k];
        z = zList[j];
        zPrevious = zList[j + 1];
        //zAfter = zList[j-1];

        (*distDrDz)(i, j) = 0.0;
        (*distDPhiRDz)(i, j) = 0.0;
        (*distDz)(i, j) = 0.0;
        dr = 0.0;
        dRPhi = 0.0;
        dz = 0.0;

        r = rList[i];
        phi = phiList[k];
        z = zList[j];

        x[0] = r;
        x[1] = phi;
        x[2] = z;

        GetDistortionCylNoDrift(x, roc, dx);


        //x[0] = x[0] + dr;
        //x[1] = x[1] + dRPhi/r;
        x[2] = zPrevious;

        GetDistortionCylNoDrift(x, roc, pdx);

        (*distDrDz)(i, j) = (dx[0] - pdx[0]);
        (*distDPhiRDz)(i, j) = (dx[1] - pdx[1]) * r;
        (*distDz)(i, j) = (dx[2] - pdx[2]);

      }
    }
  }
}

///
/// \param matricesCharge
/// \param matricesEr
/// \param matricesEPhi
/// \param matricesEz
/// \param matricesInvLocalIntErDz
/// \param matricesInvLocalIntEPhiDz
/// \param matricesInvLocalEz
/// \param matricesDistDrDz
/// \param matricesDistDPhiRDz
/// \param matricesDistDz
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseDistortionMapsNoDrift(
        TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEPhi, TMatrixD **matricesEz,
        TMatrixD **matricesInvLocalIntErDz, TMatrixD **matricesInvLocalIntEPhiDz, TMatrixD **matricesInvLocalEz,
        TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz, TMatrixD **matricesDistDz,
        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice) {
  // can inverse after lookup table for global distortion been calculated
  Double_t *rList = new Double_t[nRRow];
  Double_t *zList = new Double_t[nZColumn];
  Double_t *phiList = new Double_t[phiSlice];

  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;

  for (Int_t k = 0; k < phiSlice; k++) phiList[k] = gridSizePhi * k;
  for (Int_t i = 0; i < nRRow; i++) rList[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < nZColumn; j++) zList[j] = (j * gridSizeZ);
  // memory allocation
  if (fInitLookUp) {
    // 1)  get local distortion
    InverseGlobalToLocalDistortionNoDrift(matricesDistDrDz, matricesDistDPhiRDz, matricesDistDz, rList, zList,
                                          phiList, nRRow, nZColumn, phiSlice);
    // 2)  calculate local integral
    InverseLocalDistortionToElectricField(matricesEr, matricesEPhi, matricesEz, matricesInvLocalIntErDz,
                                          matricesInvLocalIntEPhiDz, matricesInvLocalEz,
                                          matricesDistDrDz, matricesDistDPhiRDz, matricesDistDz, rList, zList,
                                          phiList, nRRow, nZColumn, phiSlice);
    // 3)  get potential from electric field assuming zero boundaries
    InverseElectricFieldToCharge(matricesCharge, matricesEr, matricesEPhi, matricesEz, rList, zList, phiList, nRRow,
                                 nZColumn, phiSlice);
  }
  delete[] zList;
  delete[] rList;
  delete[] phiList;
}

///
/// \param matricesChargeA
/// \param matricesChargeC
/// \param spaceChargeHistogram3D
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::GetChargeDensity (
        TMatrixD **matricesChargeA, TMatrixD **matricesChargeC, TH3 *spaceChargeHistogram3D,
        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice ) {
  Int_t phiSlicesPerSector = phiSlice / kNumSector;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;
  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;
  // local variables
  Float_t radius0, phi0, z0;
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t rList[nRRow], zList[nZColumn], phiList[phiSlice];
  for (Int_t k = 0; k < phiSlice; k++) phiList[k] = gridSizePhi * k;
  for (Int_t i = 0; i < nRRow; i++) rList[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < nZColumn; j++) zList[j] = j * gridSizeZ;

  TMatrixD *mCharge;
  for (Int_t side = 0; side < 2; side++) {
    for (Int_t k = 0; k < phiSlice; k++) {
      if (side == 0)
        mCharge = matricesChargeA[k];
      else
        mCharge = matricesChargeC[k];

      phi0 = phiList[k];
      for (Int_t i = 0; i < nRRow; i++) {
        radius0 = rList[i];
        for (Int_t j = 0; j < nZColumn; j++) {
          z0 = zList[j];
          if (side == 1) z0 = -TMath::Abs(zList[j]);
          if (spaceChargeHistogram3D != NULL) {
            (*mCharge)(i, j) = InterpolatePhi(spaceChargeHistogram3D, phi0, radius0, z0);
            //InterpolatePhi(spaceChargeHistogram3D,phi0,radius0,z0);
          }
        }
      }
    }
  }
}

///
/// \param x
/// \param roc
/// \return
Double_t AliTPCSpaceCharge3DDriftLine::GetChargeCylAC(const Float_t x[], Short_t roc) {
  Double_t r, phi, z;
  Int_t sign;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  if (z > 0)
    return fInterpolatorChargeA->GetValue(r, phi, z);
  else
    return fInterpolatorChargeC->GetValue(r, phi, -z);
}

/// chargeInverse
///
/// \param x
/// \param roc
/// \return
Double_t AliTPCSpaceCharge3DDriftLine::GetInverseChargeCylAC(const Float_t x[], Short_t roc) {
  Double_t r, phi, z;
  Int_t sign;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  if (z > -1e-6)
    return fInterpolatorInverseChargeA->GetValue(r, phi, z);
  else
    return fInterpolatorInverseChargeC->GetValue(r, phi, z);
}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetLocalDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]) {
  Float_t dR, dRPhi, dZ;
  Double_t r, phi, z;
  Int_t sign;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  if (z > -1e-6)
    fLookupDistA->GetValue(r, phi, z, dR, dRPhi, dZ);
  else
    fLookupDistC->GetValue(r, phi, -z, dR, dRPhi, dZ);

  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dRPhi;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with drift velocity dependency on the Ez field and the overall scaling factor)

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetInverseLocalDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the initialization now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  Float_t dR, dRPhi, dZ;
  Double_t r, phi, z;
  Int_t sign;

  r = x[0];
  phi = x[1];
  if (phi < 0) phi += TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();                   // Table uses phi from 0 to 2*Pi

  z = x[2];                                         // Create temporary copy of x[2]

  if ((roc % 36) < 18) {
    sign = 1;       // (TPC A side)
  } else {
    sign = -1;       // (TPC C side)
  }

  if (sign == 1 && z < fgkZOffSet) z = fgkZOffSet;    // Protect against discontinuity at CE
  if (sign == -1 && z > -fgkZOffSet) z = -fgkZOffSet;    // Protect against discontinuity at CE


  if ((sign == 1 && z < 0) || (sign == -1 && z > 0)) // just a consistency check
    AliError("ROC number does not correspond to z coordinate! Calculation of distortions is most likely wrong!");

  if (z > -1e-6)
    fLookupInverseDistA->GetValue(r, phi, z, dR, dRPhi, dZ);
  else
    fLookupInverseDistC->GetValue(r, phi, -z, dR, dRPhi, dZ);

  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dRPhi;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with drift velocity dependency on the Ez field and the overall scaling factor)

}

/// Function for setting Potential Boundary Values and Charge distribution input TFormula
///
/// \param vTestFunction
/// \param rhoTestFunction
///
void AliTPCSpaceCharge3DDriftLine::SetPotentialBoundaryAndCharge (TFormula *vTestFunction,  TFormula *rhoTestFunction) {
  /**** allocate memory for charge ***/
  // we allocate pointer to TMatrixD array to picture 3D (slices), this representation should be revised
  // since it is easier for GPU implementation to run for 1D memory

  fFormulaPotentialV = new TFormula(*vTestFunction);
  fFormulaChargeRho = new TFormula(*rhoTestFunction);

  // grid size for one side
  TMatrixD *chargeA;
  TMatrixD *chargeC;

  Double_t radius0, z0, phi0, z0neg;

  Int_t indexB = 0;
  for (Int_t k = 0; k < fNPhiSlices; k++) {
    chargeA = fMatrixChargeA[k];
    chargeC = fMatrixChargeC[k];

    phi0 = fListPhi[k];

    /// Fill the non-boundary values
    for (Int_t i = 0; i < fNRRows; i++) {
      radius0 = fListR[i];
      for (Int_t j = 0; j < fNZColumns; j++) {
        z0 = fListZ[j];
        z0neg = -z0;

        (*chargeA)(i, j) = -1.0 * rhoTestFunction->Eval(radius0, phi0, z0);
        (*chargeC)(i, j) = -1.0 * rhoTestFunction->Eval(radius0, phi0, z0neg);

        if ((i == 0) || (i == fNRRows - 1) || (j == 0) || (j == fNZColumns - 1)) {
          fListPotentialBoundaryA[indexB] = vTestFunction->Eval(radius0, phi0, z0);
          fListPotentialBoundaryC[indexB] = vTestFunction->Eval(radius0, phi0, z0neg);
          indexB++;
        }

      } // end j
    } // end i
  } // end phi

  fInterpolatorChargeA->SetValue(fMatrixChargeA);
  fInterpolatorChargeA->InitCubicSpline();
  fInterpolatorChargeC->SetValue(fMatrixChargeC);
  fInterpolatorChargeC->InitCubicSpline();
}

