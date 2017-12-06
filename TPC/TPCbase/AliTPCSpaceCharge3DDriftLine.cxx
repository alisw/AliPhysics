/*************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyrigh notice and this permission notice   *
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
/// fInterpolationOrder = 1; //triliniear 3d //2 there is an error
/// fStrategy = kUseInterpolator;	
/// fNRRows = 72;
/// fNPhiSlices = 180; // the maximum of phi-slices so far = (8 per sector)
/// fNZColumns = 166; // the maximum on column-slices so  ~ 2cm slicing
/// ~~~
///  
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine()
        : AliTPCCorrection(), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE), fInterpolationOrder(5),
          fIrregularGridSize(3), fStrategy(0), fRBFKernelType(0), fNRRows(129), fNZColumns(129), fNPhiSlices(144),
          fCorrectionType(0) {
  InitAllocateMemory();
}

/// Construction for AliTPCSpaceCharge3DDriftLine class
/// Default values
/// ~~~
/// fInterpolationOrder = 1; //triliniear 3d //2 there is an error
/// fStrategy = kUseInterpolator;
/// fNRRows = 72;
/// fNPhiSlices = 180; // the maximum of phi-slices so far = (8 per sector)
/// fNZColumns = 166; // the maximum on column-slices so  ~ 2cm slicing
/// ~~~
///
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(const char * name, const char *title)
        : AliTPCCorrection(name,title), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE), fInterpolationOrder(5),
          fIrregularGridSize(3), fStrategy(0), fRBFKernelType(0), fNRRows(129), fNZColumns(129), fNPhiSlices(144),
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
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine(const char * name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice ):
        AliTPCCorrection(name,title), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE), fInterpolationOrder(5),
        fIrregularGridSize(3), fStrategy(0), fRBFKernelType(0), fCorrectionType(0) {
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
AliTPCSpaceCharge3DDriftLine::AliTPCSpaceCharge3DDriftLine (
        const char * name, const char *title, Int_t nRRow, Int_t nZColumn, Int_t nPhiSlice, Int_t interpolationOrder,
        Int_t irregularGridSize, Int_t strategy, Int_t rbfKernelType )
        : AliTPCCorrection(name,title), fC0(0.), fC1(0.), fCorrectionFactor(1.), fInitLookUp(kFALSE), fCorrectionType(0) {
  fInterpolationOrder = interpolationOrder; //triliniear 3d //2 there is an error
  fIrregularGridSize = irregularGridSize; // default 3
  fStrategy = strategy;

  fNRRows = nRRow;
  fNPhiSlices = nPhiSlice; // the maximum of phi-slices so far = (8 per sector)
  fNZColumns = nZColumn; // the maximum on column-slices so  ~ 2cm slicing
  fRBFKernelType = rbfKernelType;

  InitAllocateMemory();
}

/// Memory allocation for working/output memory
///
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
    fMatrixIntDistDrEz[k] = new TMatrixD(fNRRows, fNZColumns * 2);
    fMatrixIntDistDphiREz[k] = new TMatrixD(fNRRows, fNZColumns * 2);
    fMatrixIntDistDz[k] = new TMatrixD(fNRRows, fNZColumns * 2);

    fMatrixIntDistDrEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDphiREzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDzA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntDistDrEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDphiREzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntDistDzC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEz[k] = new TMatrixD(fNRRows, fNZColumns * 2);
    fMatrixIntCorrDphiREz[k] = new TMatrixD(fNRRows, fNZColumns * 2);
    fMatrixIntCorrDz[k] = new TMatrixD(fNRRows, fNZColumns * 2);

    fMatrixIntCorrDrEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDphiREzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDzA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDphiREzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixIntCorrDzC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixErOverEz[k] = new TMatrixD(fNRRows, fNZColumns * 2);
    fMatrixEphiOverEz[k] = new TMatrixD(fNRRows, fNZColumns * 2);
    fMatrixDeltaEz[k] = new TMatrixD(fNRRows, fNZColumns * 2);

    fMatrixErOverEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixEphiOverEzA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixDeltaEzA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixErOverEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixEphiOverEzC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixDeltaEzC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEzIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);        //[kNPhi]
    fMatrixIntCorrDphiREzIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);   //[kNPhi]
    fMatrixIntCorrDzIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixRListIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixPhiListIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixZListIrregularA[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixIntCorrDrEzIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);        //[kNPhi]
    fMatrixIntCorrDphiREzIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);   //[kNPhi]
    fMatrixIntCorrDzIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixRListIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixPhiListIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixZListIrregularC[k] = new TMatrixD(fNRRows, fNZColumns);

    fMatrixChargeA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixChargeC[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixChargeInverseA[k] = new TMatrixD(fNRRows, fNZColumns);
    fMatrixChargeInverseC[k] = new TMatrixD(fNRRows, fNZColumns);
  }

  fLookupIntDist =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntDistDrEz, fListR, fNPhiSlices, fMatrixIntDistDphiREz, fListPhi,
                  fNZColumns * 2, fMatrixIntDistDz, fListZ, fInterpolationOrder);
  fLookupIntCorr =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntCorrDrEz, fListR, fNPhiSlices, fMatrixIntCorrDphiREz, fListPhi,
                  fNZColumns * 2, fMatrixIntCorrDz, fListZ, fInterpolationOrder);
  fLookupIntDistA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntDistDrEzA, fListR, fNPhiSlices, fMatrixIntDistDphiREzA, fListPhi,
                  fNZColumns, fMatrixIntDistDzA, fListZA, fInterpolationOrder);
  fLookupIntDistC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntDistDrEzC, fListR, fNPhiSlices, fMatrixIntDistDphiREzC, fListPhi,
                  fNZColumns, fMatrixIntDistDzC, fListZC, fInterpolationOrder);
  fLookupIntCorrA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntCorrDrEzA, fListR, fNPhiSlices, fMatrixIntCorrDphiREzA, fListPhi,
                  fNZColumns, fMatrixIntCorrDzA, fListZA, fInterpolationOrder);
  fLookupIntCorrC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixIntCorrDrEzC, fListR, fNPhiSlices, fMatrixIntCorrDphiREzC, fListPhi,
                  fNZColumns, fMatrixIntCorrDzC, fListZC, fInterpolationOrder);
  fLookupIntENoDrift =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixErOverEz, fListR, fNPhiSlices, fMatrixEphiOverEz, fListPhi,
                  fNZColumns * 2, fMatrixDeltaEz, fListZ, fInterpolationOrder);
  fLookupIntENoDriftA =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixErOverEzA, fListR, fNPhiSlices, fMatrixEphiOverEzA, fListPhi,
                  fNZColumns, fMatrixDeltaEzA, fListZA, fInterpolationOrder);
  fLookupIntENoDriftC =
          new AliTPCLookUpTable3DInterpolatorD(
                  fNRRows, fMatrixErOverEzC, fListR, fNPhiSlices, fMatrixEphiOverEzC, fListPhi,
                  fNZColumns, fMatrixDeltaEzC, fListZC, fInterpolationOrder);
  fLookupIntCorrIrregularA =
          new AliTPCLookUpTable3DInterpolatorDFull(
                  fNRRows, fMatrixIntCorrDrEzIrregularA, fMatrixRListIrregularA, fListR, fNPhiSlices,
                  fMatrixIntCorrDphiREzIrregularA, fMatrixPhiListIrregularA, fListPhi, fNZColumns,
                  fMatrixIntCorrDzIrregularA, fMatrixZListIrregularA, fListZA, 2, GetIrregularGridSize(),
                  GetIrregularGridSize(), GetIrregularGridSize(), 1);

  fLookupIntCorrIrregularC =
          new AliTPCLookUpTable3DInterpolatorDFull(
                  fNRRows, fMatrixIntCorrDrEzIrregularC, fMatrixRListIrregularC, fListR, fNPhiSlices,
                  fMatrixIntCorrDphiREzIrregularC, fMatrixPhiListIrregularC, fListPhi, fNZColumns,
                  fMatrixIntCorrDzIrregularC, fMatrixZListIrregularC, fListZC, 2, GetIrregularGridSize(),
                  GetIrregularGridSize(), GetIrregularGridSize(), 1);

  fInterpolatorChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorChargeC = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeA = new AliTPC3DCylindricalInterpolator();
  fInterpolatorInverseChargeC = new AliTPC3DCylindricalInterpolator();

  // should be in contructor
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
    delete fMatrixIntDistDrEz[k];
    delete fMatrixIntDistDphiREz[k];
    delete fMatrixIntDistDz[k];
    delete fMatrixIntDistDrEzA[k];
    delete fMatrixIntDistDphiREzA[k];
    delete fMatrixIntDistDzA[k];
    delete fMatrixIntDistDrEzC[k];
    delete fMatrixIntDistDphiREzC[k];
    delete fMatrixIntDistDzC[k];
    delete fMatrixIntCorrDrEz[k];
    delete fMatrixIntCorrDphiREz[k];
    delete fMatrixIntCorrDz[k];
    delete fMatrixIntCorrDrEzA[k];
    delete fMatrixIntCorrDphiREzA[k];
    delete fMatrixIntCorrDzA[k];
    delete fMatrixIntCorrDrEzC[k];
    delete fMatrixIntCorrDphiREzC[k];
    delete fMatrixIntCorrDzC[k];
    delete fMatrixErOverEz[k];
    delete fMatrixEphiOverEz[k];
    delete fMatrixDeltaEz[k];
    delete fMatrixErOverEzA[k];
    delete fMatrixEphiOverEzA[k];
    delete fMatrixDeltaEzA[k];
    delete fMatrixErOverEzC[k];
    delete fMatrixEphiOverEzC[k];
    delete fMatrixDeltaEzC[k];
    delete fMatrixIntCorrDrEzIrregularA[k];
    delete fMatrixIntCorrDphiREzIrregularA[k];
    delete fMatrixIntCorrDzIrregularA[k];
    delete fMatrixRListIrregularA[k];
    delete fMatrixPhiListIrregularA[k];
    delete fMatrixZListIrregularA[k];
    delete fMatrixIntCorrDrEzIrregularC[k];
    delete fMatrixIntCorrDphiREzIrregularC[k];
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

  delete fLookupIntDist;
  delete fLookupIntDistA;
  delete fLookupIntDistC;
  delete fLookupIntENoDriftA;
  delete fLookupIntENoDrift;
  delete fLookupIntENoDriftC;
  delete fLookupIntCorr;
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
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magnetic field - not initialized");
  Double_t bzField = magF->SolenoidField()/10.; //field in T
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  if (!param) AliError("Parameters - not initialized");
  Double_t vdrift = param->GetDriftV()/1000000.; // [cm/us]   // From dataBase: to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ;
  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  SetOmegaTauT1T2(wt,fT1,fT2);
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
/// ElectricField( matricesV, matricesEr,  matricesEphi, matricesEz, nRRow, nZColumn, phiSlice,
/// gridSizeR, gridSizePhi ,gridSizeZ,symmetry, fgkIFCRadius);
/// ~~~
///
/// 3) Calculate local distortion and correction, useing Langevin formula
/// ~~~ cxx
/// LocalDistCorrDz (matricesEr, matricesEphi, 	matricesEz,
///	matricesDistDrDz,  matricesDistDphiRDz, matricesDistDz,
///	matricesCorrDrDz,  matricesCorrDphiRDz, matricesCorrDz,
///	nRRow,  nZColumn, phiSlice, gridSizeZ, ezField);
/// ~~~
///
/// 4) Integrate distortion by following the drift line
///
/// 5) Fill look up table for Integral distortion
///
/// 6) Fill look up table for Integral correction
///
/// \param nRRow Int_t Number of rows in r-direction
/// \param nZColumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
/// \param maxIteration Int_t Maximum iteration for poisson solver
/// \param stoppingConv Convergence error stopping conditioin for poisson solver
///
/// \post Lookup tables for distortion:
/// ~~~
/// fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz
/// ~~~
/// and correction:
/// ~~~
/// fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz
/// ~~~
/// are initialized
///
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoissonIntegralDz
        (
                Int_t nRRow,
                Int_t nZColumn,
                Int_t phiSlice,
                Int_t maxIteration,
                Double_t stoppingConv
        ) {
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
  TMatrixD *matricesEr[phiSlice], *matricesEphi[phiSlice], *matricesEz[phiSlice];
  TMatrixD *matricesDistDrDz[phiSlice], *matricesDistDphiRDz[phiSlice], *matricesDistDz[phiSlice];
  TMatrixD *matricesCorrDrDz[phiSlice], *matricesCorrDphiRDz[phiSlice], *matricesCorrDz[phiSlice];
  TMatrixD *matricesGDistDrDz[phiSlice], *matricesGDistDphiRDz[phiSlice], *matricesGDistDz[phiSlice];
  TMatrixD *matricesGCorrDrDz[phiSlice], *matricesGCorrDphiRDz[phiSlice], *matricesGCorrDz[phiSlice];


  for (Int_t k = 0; k < phiSlice; k++) {
    matricesV[k] = new TMatrixD(nRRow, nZColumn);
    matricesCharge[k] = new TMatrixD(nRRow, nZColumn);
    matricesEr[k] = new TMatrixD(nRRow, nZColumn);
    matricesEphi[k] = new TMatrixD(nRRow, nZColumn);
    matricesEz[k] = new TMatrixD(nRRow, nZColumn);
    matricesDistDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesDistDphiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesDistDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDphiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGDistDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGDistDphiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGDistDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGCorrDrDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGCorrDphiRDz[k] = new TMatrixD(nRRow, nZColumn);
    matricesGCorrDz[k] = new TMatrixD(nRRow, nZColumn);

  }

  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t rList[nRRow], zList[nZColumn], phiList[phiSlice];

  // pointer to current TF1 for potentital boundary values
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
                  nRRow, matricesDistDrDz, rList, phiSlice, matricesDistDphiRDz, phiList, nZColumn, matricesDistDz,
                  zList, fInterpolationOrder);

  // allocate look up local corection
  AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr =
          new AliTPCLookUpTable3DInterpolatorD(
                  nRRow, matricesCorrDrDz, rList, phiSlice, matricesCorrDphiRDz, phiList, nZColumn, matricesCorrDz,
                  zList, fInterpolationOrder);

  // allocate look up for global distortion
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalDist =
          new AliTPCLookUpTable3DInterpolatorD(
                  nRRow, matricesGDistDrDz, rList, phiSlice, matricesGDistDphiRDz, phiList, nZColumn, matricesGDistDz,
                  zList, fInterpolationOrder);
  // allocate look up for global distortion
  AliTPCLookUpTable3DInterpolatorD *lookupGlobalCorr =
          new AliTPCLookUpTable3DInterpolatorD(
                  nRRow, matricesGCorrDrDz, rList, phiSlice, matricesGCorrDphiRDz, phiList, nZColumn, matricesGCorrDz,
                  zList, fInterpolationOrder);

  // should be set, in another place
  const Int_t symmetry = 0; // fSymmetry

  // for irregular
  TMatrixD **matricesIrregularDrDz = NULL;
  TMatrixD **matricesIrregularDphiRDz = NULL;
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
        matricesDistDphiRDz[k]->Zero();
        matricesDistDz[k]->Zero();
        matricesCorrDrDz[k]->Zero();
        matricesCorrDphiRDz[k]->Zero();
        matricesCorrDz[k]->Zero();


        matricesGDistDrDz[k]->Zero();
        matricesGDistDphiRDz[k]->Zero();
        matricesGDistDz[k]->Zero();
        matricesGCorrDrDz[k]->Zero();
        matricesGCorrDphiRDz[k]->Zero();
        matricesGCorrDz[k]->Zero();

        if (side == 0) {
          matricesIrregularDrDz = fMatrixIntCorrDrEzIrregularA;
          matricesIrregularDphiRDz = fMatrixIntCorrDphiREzIrregularA;
          matricesIrregularDz = fMatrixIntCorrDzIrregularA;
          matricesPhiIrregular = fMatrixPhiListIrregularA;
          matricesRIrregular = fMatrixRListIrregularA;
          matricesZIrregular = fMatrixZListIrregularA;
          matricesLookUpCharge = fMatrixChargeA;
          chargeInterpolator = fInterpolatorChargeA;
          fLookupDistA->SetLookUpR(matricesDistDrDz);
          fLookupDistA->SetLookUpPhi(matricesDistDphiRDz);
          fLookupDistA->SetLookUpZ(matricesDistDz);
          potentialBoundary = fListPotentialBoundaryA;
          f1BoundaryIFC = fFormulaBoundaryIFCA;
          f1BoundaryOFC = fFormulaBoundaryOFCA;
          f1BoundaryROC = fFormulaBoundaryROCA;
        } else {
          matricesIrregularDrDz = fMatrixIntCorrDrEzIrregularC;
          matricesIrregularDphiRDz = fMatrixIntCorrDphiREzIrregularC;
          matricesIrregularDz = fMatrixIntCorrDzIrregularC;
          matricesPhiIrregular = fMatrixPhiListIrregularC;
          matricesRIrregular = fMatrixRListIrregularC;
          matricesZIrregular = fMatrixZListIrregularC;
          matricesLookUpCharge = fMatrixChargeC;
          chargeInterpolator = fInterpolatorChargeC;
          fLookupDistC->SetLookUpR(matricesDistDrDz);
          fLookupDistC->SetLookUpPhi(matricesDistDphiRDz);
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
              if (i==0) {
                if (f1BoundaryIFC != NULL) {
                  (*matrixV)(i, j) = f1BoundaryIFC->Eval(z0);
                }
              }
              if (i== (nRRow - 1)) {
                if (f1BoundaryOFC != NULL)
                  (*matrixV)(i, j) = f1BoundaryOFC->Eval(z0);
              }
              if (j== 0) {
                if (fFormulaBoundaryCE) {
                  (*matrixV)(i, j) = fFormulaBoundaryCE->Eval(radius0);
                }
              }
              if (j== (nZColumn - 1)) {
                if (f1BoundaryROC != NULL)
                  (*matrixV)(i, j) = f1BoundaryROC->Eval(radius0);
              }
            } else {
              if ((i==0) || (i == (nRRow -1)) || (j==0) || (j == (nZColumn-1)))
              {
                (*matrixV)(i, j) = fFormulaPotentialV->Eval(radius0, phi0, z0);
              }
            }
          }
        }
      }
      AliInfo(Form("Step 0: Preparing Charge interpolator: %f\n", w.CpuTime()));
      AliTPCPoissonSolver::fgConvErr = stoppingConv;

      fPoissonSolver->SetStrategy(AliTPCPoissonSolver::kMultigrid);
      (fPoissonSolver->fMgParameters).cycleType = AliTPCPoissonSolver::kFCycle;
      (fPoissonSolver->fMgParameters).isFull3D = kFALSE;
      (fPoissonSolver->fMgParameters).NMGCYCLE = maxIteration;
      (fPoissonSolver->fMgParameters).MAXLOOP = 6;

      w.Start();
      fPoissonSolver->PoissonSolver3D(matricesV, matricesCharge, nRRow, nZColumn, phiSlice, maxIteration, symmetry);
      w.Stop();

      AliInfo(Form("Step 1: Poisson solver: %f\n", w.CpuTime()));
      w.Start();
      ElectricField(matricesV,
                    matricesEr, matricesEphi, matricesEz, nRRow, nZColumn, phiSlice,
                    gridSizeR, gridSizePhi, gridSizeZ, symmetry, fgkIFCRadius);
      w.Stop();


      AliInfo(Form("Step 2: Electric Field Calculation: %f\n", w.CpuTime()));
      w.Start();
      LocalDistCorrDz(matricesEr, matricesEphi, matricesEz,
                      matricesDistDrDz, matricesDistDphiRDz, matricesDistDz,
                      matricesCorrDrDz, matricesCorrDphiRDz, matricesCorrDz,
                      nRRow, nZColumn, phiSlice, gridSizeZ, ezField);
      w.Stop();

      // copy to interpolator
      if (side == 0) {
        lookupLocalDist->CopyVals();
        lookupLocalCorr->CopyVals();
        fLookupDistA->CopyVals();
      } else {
        lookupLocalDist->CopyVals();
        lookupLocalCorr->CopyVals();
        fLookupDistC->CopyVals();
      }

      AliInfo(Form("Step 3: Local distortion and correction: %f\n", w.CpuTime()));
      w.Start();

      if (fStrategy == kNaive)
        IntegrateDistCorrDriftLineDz(
                lookupLocalDist,
                matricesGDistDrDz, matricesGDistDphiRDz, matricesGDistDz,
                lookupLocalCorr,
                matricesGCorrDrDz, matricesGCorrDphiRDz, matricesGCorrDz,
                matricesIrregularDrDz, matricesIrregularDphiRDz, matricesIrregularDz,
                matricesRIrregular, matricesPhiIrregular, matricesZIrregular,
                nRRow, nZColumn, phiSlice, rList, phiList, zList
        );
      else
        IntegrateDistCorrDriftLineDzOpt2(
                lookupLocalDist,
                matricesGDistDrDz, matricesGDistDphiRDz, matricesGDistDz,
                lookupLocalCorr,
                matricesGCorrDrDz, matricesGCorrDphiRDz, matricesGCorrDz,
                nRRow, nZColumn, phiSlice, rList, phiList, zList
        );


      w.Stop();
      AliInfo(Form("Step 4: Global correction/distortion: %f\n", w.CpuTime()));
      w.Start();

      //// copy to 1D interpolator /////
      lookupGlobalDist->CopyVals();
      lookupGlobalCorr->CopyVals();
      ////


      w.Stop();
      AliInfo(Form("Step 5: Filling up the look up: %f\n", w.CpuTime()));


      if (side == 0) {
        FillLookUpTableA(lookupGlobalDist,
                         fMatrixIntDistDrEzA, fMatrixIntDistDphiREzA, fMatrixIntDistDzA,
                         nRRow, nZColumn, phiSlice, rList, phiList, zList);

        FillLookUpTableA(lookupGlobalCorr,
                         fMatrixIntCorrDrEzA, fMatrixIntCorrDphiREzA, fMatrixIntCorrDzA,
                         nRRow, nZColumn, phiSlice, rList, phiList, zList);


        fLookupIntDistA->CopyVals();
        fLookupIntCorrA->CopyVals();

        fLookupIntCorrIrregularA->CopyVals();


        AliInfo(" A side done");
      }
      if (side == 1) {
        FillLookUpTableC(lookupGlobalDist,
                       fMatrixIntDistDrEzC, fMatrixIntDistDphiREzC, fMatrixIntDistDzC,
                         nRRow, nZColumn, phiSlice, rList, phiList, zList);

        FillLookUpTableC(lookupGlobalCorr,
                        fMatrixIntCorrDrEzC, fMatrixIntCorrDphiREzC, fMatrixIntCorrDzC,
                        nRRow, nZColumn, phiSlice, rList, phiList, zList);

        fLookupIntDistC->CopyVals();
        fLookupIntCorrC->CopyVals();
        fLookupIntCorrIrregularC->CopyVals();
        AliInfo(" C side done");
      }

    }

    fInitLookUp = kTRUE;
  }



  // memory deallocation for temporary matrices
  for (Int_t k = 0; k < phiSlice; k++) {
    delete matricesV[k];
    delete matricesCharge[k];
    delete matricesEr[k];
    delete matricesEphi[k];
    delete matricesEz[k];
    delete matricesDistDrDz[k];
    delete matricesDistDphiRDz[k];
    delete matricesDistDz[k];

    delete matricesCorrDrDz[k];
    delete matricesCorrDphiRDz[k];
    delete matricesCorrDz[k];
    delete matricesGDistDrDz[k];
    delete matricesGDistDphiRDz[k];
    delete matricesGDistDz[k];

    delete matricesGCorrDrDz[k];
    delete matricesGCorrDphiRDz[k];
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
/// \param phiSlice 	Int_t number of slicees in phi direction
/// \param maxIteration Int_t max iteration for convergence
/// \param stoppingConv Double_t stopping criteria for convergence
/// \post Lookup tables for distortion:
/// ~~~
/// fLookUpIntDistDrEz,fLookUpIntDistDphiREz,fLookUpIntDistDz
/// ~~~ fo
/// and correction:
/// ~~~
/// fLookUpIntCorrDrEz,fLookUpIntCorrDphiREz,fLookUpIntCorrDz
/// ~~~
/// are initialized
///
void AliTPCSpaceCharge3DDriftLine::InitSpaceCharge3DPoisson
        (
                Int_t nRRow,
                Int_t nZColumn,
                Int_t phiSlice,
                Int_t maxIteration,
                Double_t stoppingConv
        ) {
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
  TMatrixD *matricesEr[phiSlice], *matricesEphi[phiSlice], *matricesEz[phiSlice];

  for (Int_t k = 0; k < phiSlice; k++) {
    matricesEr[k] = new TMatrixD(nRRow, nZColumn);
    matricesEphi[k] = new TMatrixD(nRRow, nZColumn);
    matricesEz[k] = new TMatrixD(nRRow, nZColumn);
    matricesV[k] = new TMatrixD(nRRow, nZColumn);
    matricesCharge[k] = new TMatrixD(nRRow, nZColumn);
  }

  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t rlist[nRRow], zedlist[nZColumn], philist[phiSlice];

  for (Int_t k = 0; k < phiSlice; k++) philist[k] = gridSizePhi * k;
  for (Int_t i = 0; i < nRRow; i++) rlist[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < nZColumn; j++) zedlist[j] = j * gridSizeZ;
  // should be set, in another place
  const Int_t symmetry = 0;
  // do if look up table haven't be initialized
  if (!fInitLookUp) {
    for (Int_t side = 0; side < 2; side++) {
      for (Int_t k = 0; k < phiSlice; k++) {
        TMatrixD *mV = matricesV[k];
        TMatrixD *mCharge = matricesCharge[k];
        phi0 = philist[k];
        for (Int_t i = 0; i < nRRow; i++) {
          radius0 = rlist[i];
          for (Int_t j = 0; j < nZColumn; j++) {
            z0 = zedlist[j];
            if (side == 1) z0 = -TMath::Abs(zedlist[j]);
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
                      rlist, phiSlice,
                      matricesEphi,
                      philist, nZColumn,
                      matricesEz,
                      zedlist,
                      fInterpolationOrder
              );


      AliInfo("Step 1: Solving poisson solver");
      fPoissonSolver->PoissonSolver3D(matricesV, matricesCharge, nRRow, nZColumn, phiSlice, maxIteration, symmetry);
      AliInfo("Step 2: Calculate electric field");
      CalculateEField(
              matricesV,
              matricesEr,
              matricesEphi,
              matricesEz,
              nRRow,
              nZColumn,
              phiSlice,
              maxIteration,
              symmetry
      );
      lookupEField->CopyVals();
      AliInfo("Step 3: Fill the ");
      FillLookUpTable(lookupEField,
                      fMatrixErOverEz, fMatrixEphiOverEz, fMatrixDeltaEz,
                      nRRow, nZColumn, phiSlice, rlist, philist, zedlist, side);

      if (side == 0) {
        FillLookUpTableA(lookupEField,
                         fMatrixErOverEzA, fMatrixEphiOverEzA, fMatrixDeltaEzA,
                         nRRow, nZColumn, phiSlice, rlist, philist, zedlist);
        fLookupIntENoDriftA->CopyVals();
        AliInfo(" A side done");
      }
      if (side == 1) {
        FillLookUpTableC(lookupEField,
                         fMatrixErOverEzC, fMatrixEphiOverEzC, fMatrixDeltaEzC,
                         nRRow, nZColumn, phiSlice, rlist, philist, zedlist);
        fLookupIntENoDriftC->CopyVals();

        AliInfo(" C side done");
      }
      delete lookupEField;
    }
    fInitLookUp = kTRUE;
    fLookupIntENoDrift->CopyVals();
  }

  for (Int_t k = 0; k < phiSlice; k++) {
    delete matricesV[k];
    delete matricesCharge[k];
    delete matricesEr[k];
    delete matricesEphi[k];
    delete matricesEz[k];
  }
}


/// Force creating look-up table of Correction/Distortion by integration following
/// drift line.
///
/// \param nRRow Int_t Number of rows in r-direction
/// \param nZColumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
/// \param maxIteration Int_t Maximum iteration for poisson solver
/// \param stoppingConv Convergence error stopping conditioin for poisson solver
///
void AliTPCSpaceCharge3DDriftLine::ForceInitSpaceCharge3DPoissonIntegralDz
        (
                Int_t nRRow,
                Int_t nZColumn,
                Int_t phiSlice,
                Int_t maxIteration,
                Double_t stoppingConv
        ) {
  fInitLookUp = kFALSE;
  InitSpaceCharge3DPoissonIntegralDz(
          nRRow,
          nZColumn,
          phiSlice,
          maxIteration,
          stoppingConv
  );
}

///
/// Electricfield Calculation:
///
/// \param arrayofArrayV TMatrixD** 3D matrix representing calculated potential
/// \param arrayofArrayEr TMatrix** 3D matrix representing e-field at Er
/// \param arrayofArrayEz TMatrix** 3D matrix representing e-field at Ez
/// \param arrayofArrayEphi TMatrix** 3D matrix representing e-field at Ephi
/// \param rows Int_t number of rows of discritization (in R direction)
/// \param columns Int_t number of columns  of discritization (in Z direction)
/// \param phislices Int_t number of (phislices in phi direction)
/// \param symmetry Int_t symmetry?
///
/// \pre   Matrix arrayofArrayV is assumed had been calculated  by Poisson solver
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
///   arrayEr(i,j) = -1 * ( arrayV(i+1,j) - arrayV(i-1,j) ) / (2*gridSizeR); // r direction
///		arrayEz(i,j) = -1 * ( arrayV(i,j+1) - arrayV(i,j-1) ) / (2*gridSizeZ) ; // z direction
///		arrayEphi(i,j) = -1 * (signplus * arrayVP(i,j) - signminus * arrayVM(i,j) ) / (2*radius*gridSizePhi) ; // phi didrection
///   ~~~
///
/// * Boundary -> Forward/Backward difference (3 stencil) TODO: 5 Stencil
///
///   \f$ -\nabla_{r} V(r_{0},\phi_{j},z_{k}) \approx -( -0.5 V_{2,j,k} + 2 V_{1,j,k} - 1.5 * V_{0,j,k}) /  h_{r} \f$
///
///   \f$ -\nabla_{r} V(r_{nRRow - 1},\phi_{j},z_{k}) \approx -( 1.5 V_{nRRow-1,j,k} - 2.0 V_{nRRow-2,j,k} + 0.5 V_{nRRow -3,j,k}) / h_{\phi} \f$
///
void AliTPCSpaceCharge3DDriftLine::ElectricField
        (
                TMatrixD **matricesV,
                TMatrixD **matricesEr,
                TMatrixD **matricesEphi,
                TMatrixD **matricesEz,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Float_t gridSizeR,
                const Float_t gridSizePhi,
                const Float_t gridSizeZ,
                const Int_t symmetry,
                const Float_t innerRadius
        ) {


  Float_t radius;
  Int_t mplus, mminus, signplus, signminus;
  // iterate over phislices


  for (Int_t m = 0; m < phiSlice; m++) {
    mplus = m + 1;
    signplus = 1;
    mminus = m - 1;
    signminus = 1;
    if (symmetry == 1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (mplus > phiSlice - 1) mplus = phiSlice - 2;
      if (mminus < 0) mminus = 1;
    } else if (symmetry == -1) {       // Anti-symmetry in phi
      if (mplus > phiSlice - 1) {
        mplus = phiSlice - 2;
        signplus = -1;
      }
      if (mminus < 0) {
        mminus = 1;
        signminus = -1;
      }
    } else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if (mplus > phiSlice - 1) mplus = m + 1 - phiSlice;
      if (mminus < 0) mminus = m - 1 + phiSlice;
    }

    TMatrixD &arrayVP = *matricesV[mplus];
    TMatrixD &arrayVM = *matricesV[mminus];
    TMatrixD &arrayV = *matricesV[m];
    TMatrixD &arrayEr = *matricesEr[m];
    TMatrixD &arrayEz = *matricesEz[m];
    TMatrixD &arrayEphi = *matricesEphi[m];

    // for non-boundary V
    for (Int_t i = 1; i < nRRow - 1; i++) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 1; j < nZColumn - 1; j++) {
        arrayEr(i, j) = -1 * (arrayV(i + 1, j) - arrayV(i - 1, j)) / (2 * gridSizeR); // r direction
        arrayEz(i, j) = -1 * (arrayV(i, j + 1) - arrayV(i, j - 1)) / (2 * gridSizeZ); // z direction
        arrayEphi(i, j) = -1 * (signplus * arrayVP(i, j) - signminus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi didrection
      }
    }

    // for boundary-r
    for (Int_t j = 0; j < nZColumn; j++) {
      arrayEr(0, j) = -1 * (-0.5 * arrayV(2, j) + 2.0 * arrayV(1, j) - 1.5 * arrayV(0, j)) /
                      gridSizeR; // forward difference
      arrayEr(nRRow - 1, j) =
              -1 * (1.5 * arrayV(nRRow - 1, j) - 2.0 * arrayV(nRRow - 2, j) + 0.5 * arrayV(nRRow - 3, j)) /
              gridSizeR; // backward difference
    }

    for (Int_t i = 0; i < nRRow; i += nRRow - 1) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 1; j < nZColumn - 1; j++) {
        arrayEz(i, j) = -1 * (arrayV(i, j + 1) - arrayV(i, j - 1)) / (2 * gridSizeZ); // z direction
        arrayEphi(i, j) = -1 * (signplus * arrayVP(i, j) - signminus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi didrection
      }
    }

    // for boundary-z
    for (Int_t i = 0; i < nRRow; i++) {
      arrayEz(i, 0) = -1 * (-0.5 * arrayV(i, 2) + 2.0 * arrayV(i, 1) - 1.5 * arrayV(i, 0)) / gridSizeZ;
      arrayEz(i, nZColumn - 1) =
              -1 * (1.5 * arrayV(i, nZColumn - 1) - 2.0 * arrayV(i, nZColumn - 2) + 0.5 * arrayV(i, nZColumn - 3)) /
              gridSizeZ;
    }

    for (Int_t i = 1; i < nRRow - 1; i++) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 0; j < nZColumn; j += nZColumn - 1) {
        arrayEr(i, j) = -1 * (arrayV(i + 1, j) - arrayV(i - 1, j)) / (2 * gridSizeR); // r direction
        arrayEphi(i, j) = -1 * (signplus * arrayVP(i, j) - signminus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi didrection
      }
    }

    // corner points for Ephi
    for (Int_t i = 0; i < nRRow; i += nRRow - 1) {
      radius = innerRadius + i * gridSizeR;
      for (Int_t j = 0; j < nZColumn; j += nZColumn - 1) {
        arrayEphi(i, j) = -1 * (signplus * arrayVP(i, j) - signminus * arrayVM(i, j)) /
                          (2 * radius * gridSizePhi); // phi didrection
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
/// This integration is in \f$z\f$ direction. If no change in granurality, then we can only use trapezoidal rule.
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
/// The code sniplet at \ref impllocaldist is an implementation of the local integration of electric field.
///
/// \anchor impllocaldist
/// ~~~
/// Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = Electric Field (V/cm) Magnitude ~ -400 V/cm;
///
/// localIntErOverEz = (gridSizeZ/2.0)*((*eR)(i,j)+(*eR)(i,j+1))/(-1*ezField) ;
/// localIntEphiOverEz = (gridSizeZ/2.0)*((*ePhi)(i,j)+(*ePhi)(i,j+1))/(-1*ezField) ;
/// localIntDeltaEz = (gridSizeZ/2.0)*((*eZ)(i,j)+(*eZ)(i,j+1)) ;
/// ~~~
///
///
/// After we have local integrations for ellectric fields in each direction,
/// local distortion \f$\hat{\delta}(r_{i},z_{j},\phi_{m})\f$ is calculated by simplified Langevin equation (see Figure \ref1 (b) for illustration):
///
/// \f$ \hat{\delta}_{rE}(r_{i},z_{j},\phi_{m}) = c_{0} \int^{z_{j+1}}_{z_{j}} \frac{E_{r}}{E_{z}} dz   + c_{1}  \int^{z_{j+1}}_{z_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDrDz)(i,j) 		= fC0*localIntErOverEz   + fC1*localIntEphiOverEz;
/// ~~~
///
/// \f$ r\hat{\delta}_{\phi E}(r_{i},z_{j},\phi_{m})  = - c_{1} \int^{z_{j+1}}_{z_{j}} \frac{E_{j}}{E_{j}} dz  + c_{0} \int^{z_{j+1}}_{j_{j}} \frac{E_{\phi}}{E_{z}} dz \f$
///
/// ~~~
///	(*distDphiRDz)(i,j) = fC0*localIntEphiOverEz - fC1*localIntErOverEz ;
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
/// (*corrDphiRDz)(i,j+1) = -1* (*distDphiRDz)(i,j);
/// (*corrDz)(i,j+1)      = -1* (*distDz)(i,j);
/// ~~~
///
/// \param matricesEr TMatrixD**  electric field for \f$r\f$ component
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component
///	\param matricesDistDrDz TMatrixD**  local distortion \f$\hat{\delta}_{r}\f$
///	\param matricesDistDphiRDz TMatrixD** local distortion \f$r \hat{\delta}_{\phi}\f$
///	\param matricesDistDz TMatrixD**   local distortion \f$ \hat{\delta}_{z}\f$
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param nRRow Int_t Number of rows in r-direction
/// \param nZColumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez assummed already been calculated
/// \post Local distortion and correction are computed according simplified Langevin equation
/// ~~~
/// matricesDistDrDz,matricesDistDphiRDz,matricesDistDz
/// ~~~
/// and correction:
/// ~~~
/// matricesCorrDrDz,matricesCorrDphiRDz,matricesCorrDz
/// ~~~
///
void AliTPCSpaceCharge3DDriftLine::LocalDistCorrDz
        (
                TMatrixD **matricesEr,
                TMatrixD **matricesEphi,
                TMatrixD **matricesEz,
                TMatrixD **matricesDistDrDz,
                TMatrixD **matricesDistDphiRDz,
                TMatrixD **matricesDistDz,
                TMatrixD **matricesCorrDrDz,
                TMatrixD **matricesCorrDphiRDz,
                TMatrixD **matricesCorrDz,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Float_t gridSizeZ,
                const Double_t ezField
        ) {
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEphiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;

  // pointer declaration
  TMatrixD *eR;
  TMatrixD *ePhi;
  TMatrixD *eZ;
  TMatrixD *distDrDz;
  TMatrixD *distDphiRDz;
  TMatrixD *distDz;
  TMatrixD *corrDrDz;
  TMatrixD *corrDphiRDz;
  TMatrixD *corrDz;


  // Initialization for j == colomn-1 integration is 0.0
  for (Int_t m = 0; m < phiSlice; m++) {
    distDrDz = matricesDistDrDz[m];
    distDphiRDz = matricesDistDphiRDz[m];
    distDz = matricesDistDz[m];

    corrDrDz = matricesCorrDrDz[m];
    corrDphiRDz = matricesCorrDphiRDz[m];
    corrDz = matricesCorrDz[m];

    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDphiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

      (*corrDrDz)(i, 0) = 0.0;
      (*corrDphiRDz)(i, 0) = 0.0;
      (*corrDz)(i, 0) = 0.0;
    }
  }

  // for this case
  // use trapezoidal rule assume no ROC displacement
  for (Int_t m = 0; m < phiSlice; m++) {
    eR = matricesEr[m];
    ePhi = matricesEphi[m];
    eZ = matricesEz[m];
    distDrDz = matricesDistDrDz[m];
    distDphiRDz = matricesDistDphiRDz[m];
    distDz = matricesDistDz[m];

    corrDrDz = matricesCorrDrDz[m];
    corrDphiRDz = matricesCorrDphiRDz[m];
    corrDz = matricesCorrDz[m];

    for (Int_t j = 0; j < nZColumn - 1; j++) {
      for (Int_t i = 0; i < nRRow; i++) {
        localIntErOverEz = (gridSizeZ / 2.0) * ((*eR)(i, j) + (*eR)(i, j + 1)) / (-1 * ezField);
        localIntEphiOverEz = (gridSizeZ / 2.0) * ((*ePhi)(i, j) + (*ePhi)(i, j + 1)) / (-1 * ezField);
        localIntDeltaEz = (gridSizeZ / 2.0) * ((*eZ)(i, j) + (*eZ)(i, j + 1));


        (*distDrDz)(i, j) = fC0 * localIntErOverEz + fC1 * localIntEphiOverEz;
        (*distDphiRDz)(i, j) = fC0 * localIntEphiOverEz - fC1 * localIntErOverEz;
        (*distDz)(i, j) = localIntDeltaEz * fgkdvdE * fgkdvdE;// two times?


        (*corrDrDz)(i, j + 1) = -1 * (*distDrDz)(i, j);
        (*corrDphiRDz)(i, j + 1) = -1 * (*distDphiRDz)(i, j);
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
///	\param matricesEphi TMatrixD** electric field for \f$\phi\f$ component
///	\param matricesEz TMatrixD** electric field for \f$z\f$ component
///	\param matricesCorrDrDz TMatrixD** local correction \f$\hat{\delta}_{r}\f$
///	\param matricesCorrDphiRDz TMatrixD** local correction \f$r \hat{\delta}_{\phi}\f$
///	\param matricesCorrDz TMatrixD** local correction \f$ \hat{\delta}_{z}\f$
/// \param nRRow Int_t Number of rows in r-direction
/// \param nZColumn Int_t Number of columns in z-direction
/// \param phiSlice Int_t Number of phislices in \f$ phi \f$ direction
///	\param gridSizeZ const Float_t grid size in z direction
/// \param ezField const Double_t ezField calulated from the invoking operation
///
/// \pre matricesEr, matricesEphi, matrices Ez are provided
/// \post Local correction are computed according simplified Langevin equation
/// ~~~
/// matricesCorrDz,matricesCorrDphiRDz,matricesDistDz
/// ~~~
///
void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDz
        (
                AliTPCLookUpTable3DInterpolatorD *lookupLocalDist, TMatrixD **matricesGDistDrDz, TMatrixD **matricesGDistDphiRDz,
                TMatrixD **matricesGDistDz, AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr,  TMatrixD **matricesGCorrDrDz,
                TMatrixD **matricesGCorrDphiRDz, TMatrixD **matricesGCorrDz, TMatrixD **matricesGCorrIrregularDrDz,
                TMatrixD **matricesGCorrIrregularDphiRDz, TMatrixD **matricesGCorrIrregularDz, TMatrixD **matricesRIrregular,
                TMatrixD **matricesPhiIrregular,  TMatrixD **matricesZIrregular, const Int_t nRRow,  const Int_t nZColumn,
                const Int_t phiSlice,  const Double_t *rlist,  const Double_t *philist,  const Double_t *zlist ) {

  Float_t dr, dphir, dz, ddr, ddrphi, ddz;
  Float_t radius0, phi0, z0, radius, phi, z, radiusc;
  radiusc = 0.0;
  radius = 0.0;
  TMatrixD *mDistDrDz;
  TMatrixD *mDistDphiRDz;
  TMatrixD *mDistDz;
  TMatrixD *mCorrDrDz;
  TMatrixD *mCorrDphiRDz;
  TMatrixD *mCorrDz;

  TMatrixD *mCorrIrregularDrDz;
  TMatrixD *mCorrIrregularDphiRDz;
  TMatrixD *mCorrIrregularDz;

  TMatrixD *mRIrregular;
  TMatrixD *mPhiIrregular;
  TMatrixD *mZIrregular;

  // initialized for nZColumn - 1
  Int_t j = nZColumn - 1;
  z0 = zlist[j];

  // initialized for mdist, mcorr and irregular
  // for j = zn -1 (this is end cap)
  for (Int_t m = 0; m < phiSlice; m++) {
    phi0 = philist[m];

    mDistDrDz = matricesGDistDrDz[m];
    mDistDphiRDz = matricesGDistDphiRDz[m];
    mDistDz = matricesGDistDz[m];

    //
    mCorrDrDz = matricesGCorrDrDz[m];
    mCorrDphiRDz = matricesGCorrDphiRDz[m];
    mCorrDz = matricesGCorrDz[m];


    mCorrIrregularDrDz = matricesGCorrIrregularDrDz[m];
    mCorrIrregularDphiRDz = matricesGCorrIrregularDphiRDz[m];
    mCorrIrregularDz = matricesGCorrIrregularDz[m];

    mRIrregular = matricesRIrregular[m];
    mPhiIrregular = matricesPhiIrregular[m];
    mZIrregular = matricesZIrregular[m];


    for (Int_t i = 0; i < nRRow; i++) {
      // do from j to 0
      // follow the drift
      radius0 = rlist[i];
      phi = phi0;
      radius = radius0;

      dr = 0.0;
      dphir = 0.0;
      dz = 0.0;
      ddrphi = 0.0;

      ///
      (*mDistDrDz)(i, j) = dr;
      (*mDistDphiRDz)(i, j) = dphir;
      (*mDistDz)(i, j) = dz;

//////////////// use irregular grid look up table for correction
      // set
      (*mCorrIrregularDrDz)(i, j) = -dr;
      (*mCorrIrregularDphiRDz)(i, j) = -dphir;
      (*mCorrIrregularDz)(i, j) = -dz;


      // distorted point
      (*mRIrregular)(i, j) = radius0 + dr;
      (*mPhiIrregular)(i, j) = phi0 + (dphir / radius0);
      (*mZIrregular)(i, j) = z0 + dz;
///////////////

    }
  }

  // from j one column near end cap
  for (j = nZColumn - 2; j >= 0; j--) {


    z0 = zlist[j];
    //printf("global dist j:%d\n",j);


    for (Int_t m = 0; m < phiSlice; m++) {
      phi0 = philist[m];

      mDistDrDz = matricesGDistDrDz[m];
      mDistDphiRDz = matricesGDistDphiRDz[m];
      mDistDz = matricesGDistDz[m];

      //
      mCorrDrDz = matricesGCorrDrDz[m];
      mCorrDphiRDz = matricesGCorrDphiRDz[m];
      mCorrDz = matricesGCorrDz[m];


      mCorrIrregularDrDz = matricesGCorrIrregularDrDz[m];
      mCorrIrregularDphiRDz = matricesGCorrIrregularDphiRDz[m];
      mCorrIrregularDz = matricesGCorrIrregularDz[m];

      mRIrregular = matricesRIrregular[m];
      mPhiIrregular = matricesPhiIrregular[m];
      mZIrregular = matricesZIrregular[m];


      for (Int_t i = 0; i < nRRow; i++) {
        // do from j to 0
        // follow the drift
        radius0 = rlist[i];
        phi = phi0;
        radius = radius0;

        dr = 0.0;
        dphir = 0.0;
        dz = 0.0;
        ddrphi = 0.0;


        // follow the drift line from z=j --> nZColumn - 1
        for (Int_t jj = j; jj < nZColumn; jj++) {
          // interpolation the local distortion for current position
          phi += ddrphi / radius;
          radius = radius0 + dr;
          z = zlist[jj] + dz;

          // regulate phi
          while (phi < 0.0) phi = TMath::TwoPi() + phi;
          while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();

          //  printf("lookupdist (%f,%f,%f)\n",radius,phi,z);

          lookupLocalDist->GetValue(radius, phi, z, ddr, ddrphi, ddz);

          //printf("lookdist     (%f,%f,%f)\n",ddr,ddrphi,dz);


          // add local distortion
          dr += ddr;
          dphir += ddrphi;
          dz += ddz;


        }
        // set the global distortion after following the electron drift
        (*mDistDrDz)(i, j) = dr;
        (*mDistDphiRDz)(i, j) = dphir;
        (*mDistDz)(i, j) = dz;
/////////////// use irregular grid look up table for correction
        // set
        (*mCorrIrregularDrDz)(i, j) = -dr;
        (*mCorrIrregularDphiRDz)(i, j) = -dphir;
        (*mCorrIrregularDz)(i, j) = -dz;


        // distorted point
        (*mRIrregular)(i, j) = radius0 + dr;
        (*mPhiIrregular)(i, j) = phi0 + (dphir / radius0);
        (*mZIrregular)(i, j) = z0 + dz;
///////////////

        // put the radius to the original value
        if (j == nZColumn - 2) radiusc = radius0;

        // get global correction from j+1
        dr = (*mCorrDrDz)(i, j + 1);
        dphir = (*mCorrDphiRDz)(i, j + 1);
        dz = (*mCorrDz)(i, j + 1);
        //phi = phi0 + dphir/radiusc;
        radiusc = radius0 + dr;
        phi = phi0 + dphir / radiusc;
        z = zlist[j + 1] + dz;

        //if (phi < 0.0) printf("lookupcorr (%f,%f,%f)\n",radiusc,phi,z);
        while (phi < 0.0) phi = TMath::TwoPi() + phi;
        while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();

        lookupLocalCorr->GetValue(radiusc, phi, z, ddr, ddrphi, ddz);

        dr += ddr;
        dz += ddz;
        dphir += ddrphi;


        (*mCorrDrDz)(i, j) = dr;
        (*mCorrDphiRDz)(i, j) = dphir;
        (*mCorrDz)(i, j) = dz;


      }
    }
  }
}

/// 
/// \param lookupLocalDist 
/// \param matricesGDistDrDz 
/// \param matricesGDistDphiRDz 
/// \param matricesGDistDz 
/// \param lookupLocalCorr 
/// \param matricesGCorrDrDz 
/// \param matricesGCorrDphiRDz 
/// \param matricesGCorrDz 
/// \param nRRow 
/// \param nZColumn 
/// \param phiSlice 
/// \param rlist 
/// \param philist 
/// \param zlist 
void AliTPCSpaceCharge3DDriftLine::IntegrateDistCorrDriftLineDzOpt2
        (
                AliTPCLookUpTable3DInterpolatorD *lookupLocalDist,
                TMatrixD **matricesGDistDrDz,
                TMatrixD **matricesGDistDphiRDz,
                TMatrixD **matricesGDistDz,
                AliTPCLookUpTable3DInterpolatorD *lookupLocalCorr,
                TMatrixD **matricesGCorrDrDz,
                TMatrixD **matricesGCorrDphiRDz,
                TMatrixD **matricesGCorrDz,

                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Double_t *rlist,
                const Double_t *philist,
                const Double_t *zlist
        ) {
  Float_t dr, dphir, dz, ddr, ddrphi, ddz;
  Float_t radius0, phi0, z0, radius, phi, z;

  radius = 0.0;

  TMatrixD *mDistDrDz;
  TMatrixD *mDistDphiRDz;
  TMatrixD *mDistDz;

  TMatrixD *mCorrDrDz;
  TMatrixD *mCorrDphiRDz;
  TMatrixD *mCorrDz;


  dr = 0.0;
  dz = 0.0;
  dphir = 0.0;

  for (Int_t j = nZColumn - 2; j >= 0; j--) {
    z0 = zlist[j];


    for (Int_t m = 0; m < phiSlice; m++) {
      phi0 = philist[m];

      mCorrDrDz = matricesGCorrDrDz[m];
      mCorrDphiRDz = matricesGCorrDphiRDz[m];
      mCorrDz = matricesGCorrDz[m];;

      for (Int_t i = 0; i < nRRow; i++) {
        // do from j to 0
        // follow the drift
        radius0 = rlist[i];
        //if (j == nZColumn-2) radius = radius0;

        dr = (*mCorrDrDz)(i, j + 1);
        dphir = (*mCorrDphiRDz)(i, j + 1);
        dz = (*mCorrDz)(i, j + 1);


        radius = radius0 + dr;
        phi = phi0 + dphir / radius;
        z = zlist[j + 1] + dz;


        lookupLocalCorr->GetValue(radius, phi, z, ddr, ddrphi, ddz);

        dr += ddr;
        dz += ddz;
        dphir += ddrphi;
        (*mCorrDrDz)(i, j) = dr;
        (*mCorrDphiRDz)(i, j) = dphir;
        (*mCorrDz)(i, j) = dz;
      }
    }
  }


  dr = 0.0;
  dz = 0.0;
  dphir = 0.0;


  for (Int_t j = nZColumn - 2; j >= 0; j--) {
    z0 = zlist[j];

    for (Int_t m = 0; m < phiSlice; m++) {
      phi0 = philist[m];
      mDistDrDz = matricesGDistDrDz[m];
      mDistDphiRDz = matricesGDistDphiRDz[m];
      mDistDz = matricesGDistDz[m];

      for (Int_t i = 0; i < nRRow; i++) {
        // do from j to 0
        // follow the drift
        radius0 = rlist[i];
        phi = phi0;
        radius = radius0;
        z = z0;

        lookupLocalDist->GetValue(radius, phi, z, ddr, ddrphi, ddz);

        phi += ddrphi / radius;
        radius = radius0 + ddr;
        z = zlist[j + 1] + ddz;

        dr = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi,
                                   nRRow, nZColumn, phiSlice, rlist, zlist, philist, matricesGDistDrDz, j + 1);
        dphir = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi,
                                      nRRow, nZColumn, phiSlice, rlist, zlist, philist, matricesGDistDphiRDz,
                                      j + 1);
        dz = Interpolate3DTableCyl(fInterpolationOrder, radius, z, phi,
                                   nRRow, nZColumn, phiSlice, rlist, zlist, philist, matricesGDistDz, j + 1);


        (*mDistDrDz)(i, j) = dr + ddr;
        (*mDistDphiRDz)(i, j) = dphir + ddrphi;
        (*mDistDz)(i, j) = dz + ddz;


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
/// \param rlist 
/// \param philist 
/// \param zlist 
/// \param side 
void AliTPCSpaceCharge3DDriftLine::FillLookUpTable
        (
                AliTPCLookUpTable3DInterpolatorD *lookupGlobal,
                TMatrixD **lookupRDz,
                TMatrixD **lookupPhiRDz,
                TMatrixD **lookupDz,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Double_t *rlist,
                const Double_t *philist,
                const Double_t *zlist,
                Int_t side
        ) {
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
      z = TMath::Abs(fListZ[j]);  // Symmetric solution in Z that depends only on ABS(Z)
      if (side == 0 && fListZ[j] < 0) continue; // Skip rest of this loop if on the wrong side
      if (side == 1 && fListZ[j] > 0) continue; // Skip rest of this loop if on the wrong side
      for (Int_t i = 0; i < fNRRows; i++) {
        r = fListR[i];

        lookupGlobal->GetValue(r, phi, z, (*mR)(i, j), (*mPhiR)(i, j), (*mDz)(i, j));

        if (side == 1) (*mDz)(i, j) = -(*mDz)(i, j); // negative coordinate system on C side
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
/// \param rlist 
/// \param philist 
/// \param zlist 
void AliTPCSpaceCharge3DDriftLine::FillLookUpTableA
        (
                AliTPCLookUpTable3DInterpolatorD *lookupGlobal,
                TMatrixD **lookupRDz,
                TMatrixD **lookupPhiRDz,
                TMatrixD **lookupDz,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Double_t *rlist,
                const Double_t *philist,
                const Double_t *zlist
        ) {
  Double_t r, phi, z, rl, phil, zl;
  TMatrixD *mR;
  TMatrixD *mPhiR;
  TMatrixD *mDz;

  /// * Interpolate basicLookup tables; once for each rod, then sum the results
  for (Int_t k = 0; k < fNPhiSlices; k++) {
    phi = fListPhi[k];
    phil = philist[k];

    mR = lookupRDz[k];
    mPhiR = lookupPhiRDz[k];
    mDz = lookupDz[k];
    for (Int_t j = 0; j < fNZColumns; j++) {
      z = fListZA[j];  // Symmetric solution in Z that depends only on ABS(Z)

      zl = zlist[j];
      //			printf("(%f,%f)(%f,%f)\n",z,zl);

      for (Int_t i = 0; i < fNRRows; i++) {
        r = fListR[i];
        rl = rlist[i];


        lookupGlobal->GetValue(r, phi, z, (*mR)(i, j), (*mPhiR)(i, j), (*mDz)(i, j));

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
/// \param rlist 
/// \param philist 
/// \param zlist 
void AliTPCSpaceCharge3DDriftLine::FillLookUpTableC
        (
                AliTPCLookUpTable3DInterpolatorD *lookupGlobal,
                TMatrixD **lookupRDz,
                TMatrixD **lookupPhiRDz,
                TMatrixD **lookupDz,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Double_t *rlist,
                const Double_t *philist,
                const Double_t *zlist
        ) {
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
      z = TMath::Abs(fListZC[j]);  // Symmetric solution in Z that depends only on ABS(Z)

      for (Int_t i = 0; i < fNRRows; i++) {
        r = fListR[i];

        lookupGlobal->GetValue(r, phi, z, (*mR)(i, j), (*mPhiR)(i, j), (*mDz)(i, j));
        //(*mDz)(i, j) = -(*mDz)(i, j);
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
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  GetDistortionCylAC(x, roc, dx);

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetDistortionCylAC
        (
                const Float_t x[],
                Short_t roc,
                Float_t dx[]
        ) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }


  Float_t dR, dPhiR, dZ;
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
    fLookupIntDistA->GetValue(r, phi, z, dR, dPhiR, dZ);
  else {
    fLookupIntDistC->GetValue(r, phi, -1 * z, dR, dPhiR, dZ);
    dZ = -1 * dZ;
  }


  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}


// Get Correction from irreggular table
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylACIrregular(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }


  Double_t dR, dPhiR, dZ;
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


  // get distorion from irregular table


  Int_t ianchor = TMath::FloorNint((r - fgkIFCRadius) / gridSizeR);
  Int_t kanchor = TMath::FloorNint(phi / gridSizePhi);
  Int_t zanchor = TMath::FloorNint(z / gridSizeZ);

  if (z > 0)
    fLookupIntCorrIrregularA->GetValue(r, phi, z, dR, dPhiR, dZ, ianchor, kanchor, zanchor, fNRRows / 8 + 1,
                                       fNPhiSlices / 8 + 1, fNZColumns / 8 + 1, 0);
  else {
    fLookupIntCorrIrregularC->GetValue(r, phi, -z, dR, dPhiR, dZ, ianchor, kanchor, -zanchor, fNRRows / 8 + 1,
                                       fNPhiSlices / 8 + 1, fNZColumns / 8 + 1, 0);
    dZ = -1 * dZ;
  }


  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}

/// Get correction regular grid by following electron
/// 
/// \param x 
/// \param roc 
/// \param dx 
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylAC(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }


  Float_t dR, dPhiR, dZ;
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
    fLookupIntCorrA->GetValue(r, phi, z, dR, dPhiR, dZ);
  else {
    fLookupIntCorrC->GetValue(r, phi, -z, dR, dPhiR, dZ);
    dZ = -1 * dZ;
  }
  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)
}


void AliTPCSpaceCharge3DDriftLine::GetDistortion(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }


  Float_t pCyl[3]; // a point in cylindrical coordintae
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
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }
  //printf("calling get correction spacecharge3dintegraldz\n");
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
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }

  Float_t pCyl[3]; // a point in cylindrical coordintae
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
/// \param order
/// \param r
/// \param z
/// \param phi
/// \param nr
/// \param nz
/// \param nphi
/// \param rlist
/// \param zlist
/// \param philist
/// \param arrayofArrays
/// \param zlow
/// \return
Double_t AliTPCSpaceCharge3DDriftLine::Interpolate3DTableCyl
        (
                Int_t order,
                Double_t r,
                Double_t z,
                Double_t phi,
                Int_t nr,
                Int_t nz,
                Int_t nphi,
                const Double_t rlist[],
                const Double_t zlist[],
                const Double_t philist[],

                TMatrixD **arrayofArrays,
                const Int_t zlow
        ) {
  /// Interpolate table (TMatrix format) - 3D interpolation
  /// Float version (in order to decrease the OCDB size)

  static Int_t ilow = 0, jlow = 0, klow = 0, m = 0;
  Bool_t bExtrapolate = kFALSE;
  Float_t saveArray[5] = {0., 0., 0., 0., 0.};
  Float_t savedArray[5] = {0., 0., 0., 0., 0.};

  if (phi < 0.0) phi = TMath::TwoPi() + phi;
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();

  Search(nr, rlist, r, ilow);
  Search(nz, zlist, z, jlow);
  Search(nphi, philist, phi, klow);

  //if (jlow + order < jhigh)
  //printf("z=%f,jlow=%d,zlow=%d\n",z,jlow,zlow);
  jlow = zlow;


  if (ilow < 0) {
    ilow = 0;
  }
  if (jlow < 0) {
    jlow = 0;
  }

  if (klow < 0) klow = nphi + klow;

  if (ilow + order >= nr - 1) {
    ilow = nr - 1 - order;
  }
  if (jlow + order >= nz - 1) {
    jlow = nz - 1 - order;
  }

  for (Int_t k = 0; k < order + 1; k++) {
    m = (klow + k) % nphi;
    TMatrixD &table = *arrayofArrays[m];

    for (Int_t i = ilow; i < ilow + order + 1; i++) {

      saveArray[i - ilow] = Interpolate(&zlist[jlow], &table(i, jlow), order, z);

    }

    savedArray[k] = Interpolate(&rlist[ilow], saveArray, order, r);

  }
  return (InterpolatePhi(&philist[0], klow, nphi, savedArray, order, phi));
}

///
/// \param xArray
/// \param ilow
/// \param nx
/// \param yArray
/// \param order
/// \param x
/// \return
Double_t AliTPCSpaceCharge3DDriftLine::InterpolatePhi ( const Double_t xArray[],  const Int_t ilow,
                const Int_t nx,  const Float_t yArray[],  Int_t order, Double_t x ) {
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.

  Int_t i0 = ilow;
  Double_t xi0 = xArray[ilow];

  Int_t i1 = (ilow + 1) % nx;
  Double_t xi1 = xArray[i1];
  Int_t i2 = (ilow + 2) % nx;
  Double_t xi2 = xArray[i2];

  if (x < 0) x += TMath::TwoPi();
  if (xi0 < 0) xi0 += TMath::TwoPi();
  if (xi1 < 0) xi1 += TMath::TwoPi();
  if (xi2 < 0) xi2 += TMath::TwoPi();

  if (xi0 > xi1)
    xi1 += TMath::TwoPi();

  if (xi1 > xi2)
    xi2 += TMath::TwoPi();

  if (xi0 > x)
    x += TMath::TwoPi();

  Double_t y;
  if (order == 2) {                // Quadratic Interpolation = 2

    y = (x - xi1) * (x - xi2) * yArray[0] / ((xi0 - xi1) * (xi0 - xi2));
    y += (x - xi2) * (x - xi0) * yArray[1] / ((xi1 - xi2) * (xi1 - xi0));
    y += (x - xi0) * (x - xi1) * yArray[2] / ((xi2 - xi0) * (xi2 - xi1));
  } else {                           // Linear Interpolation = 1
    y = yArray[0] + (yArray[1] - yArray[0]) * (x - xi0) / (xi1 - xi0);
  }

  return (y);
}


///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoDistDRinXY( Float_t z, Int_t nx, Int_t ny ) {
  /// Simple plot functionality.
  /// Returns a 2d hisogram which represents the corrections in radial direction (dr)
  /// in respect to position z within the XY plane.
  /// The histogramm has nx times ny entries.

  AliTPCParam *tpcparam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dr_xy", TString::Format("%s: DRinXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dr [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3], xcyl[3];
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
      xcyl[0] = r0;
      xcyl[1] = phi0;
      GetDistortionCylAC(xcyl, roc, dx);

      //Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));

      if (tpcparam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcparam->GetPadRowRadii(36, 95)) {
        //Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
        h->SetBinContent(ix, iy, dx[0]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcparam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d histogram which represents the corrections in radial direction (dr)
/// in respect to position z within the XY plane.
/// The histogramm has nx times ny entries.
///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoCorrDRinXY
        (
                Float_t z,
                Int_t nx,
                Int_t ny
        ) {

  AliTPCParam *tpcparam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dr_xy", TString::Format("%s: DRinXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "dr [cm]",
                       nx, -250., 250., ny, -250., 250.);
  Float_t x[3], dx[3], xcyl[3];
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
      xcyl[0] = r0;
      xcyl[1] = phi0;
      GetCorrectionCylAC(xcyl, roc, dx);

      //Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));

      if (tpcparam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcparam->GetPadRowRadii(36, 95)) {
        //Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
        h->SetBinContent(ix, iy, dx[0]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcparam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d hisogram which represents the corrections in rphi direction (drphi)
/// in respect to position z within the XY plane.
/// The histogramm has nx times ny entries.
///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoDistDRPhiinXY(Float_t z, Int_t nx, Int_t ny) {

  AliTPCParam *tpcparam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("drphi_xy", TString::Format("%s: DRPhiinXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drphi [cm]",
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
      if (tpcparam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcparam->GetPadRowRadii(36, 95)) {
        Float_t phi0 = TMath::ATan2(x[1], x[0]);
        Float_t phi1 = TMath::ATan2(x[1] + dx[1], x[0] + dx[0]);

        Float_t dphi = phi1 - phi0;
        if (dphi < TMath::Pi()) dphi += TMath::TwoPi();
        if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();

        h->SetBinContent(ix, iy, r0 * dphi);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }

  delete tpcparam;
  return h;
}

/// Simple plot functionality.
/// Returns a 2d hisogram which represents the corrections in rphi direction (drphi)
/// in respect to position z within the XY plane.
/// The histogramm has nx times ny entries.
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoCorrDRPhiinXY(Float_t z, Int_t nx, Int_t ny) {

  AliTPCParam *tpcparam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("drphi_xy", TString::Format("%s: DRPhiinXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
                       "drphi [cm]",
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
      if (tpcparam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcparam->GetPadRowRadii(36, 95)) {
        Float_t phi0 = TMath::ATan2(x[1], x[0]);
        Float_t phi1 = TMath::ATan2(x[1] + dx[1], x[0] + dx[0]);

        Float_t dphi = phi1 - phi0;
        if (dphi < TMath::Pi()) dphi += TMath::TwoPi();
        if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();

        h->SetBinContent(ix, iy, r0 * dphi);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }

  delete tpcparam;
  return h;
}


/// Simple plot functionality.
/// Returns a 2d hisogram which represents the corrections in longitudinal direction (dz)
/// in respect to position z within the XY plane.
/// The histogramm has nx times ny entries.
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoDistDZinXY(Float_t z, Int_t nx, Int_t ny) {

  AliTPCParam *tpcparam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dz_xy", TString::Format("%s: DZinXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
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
      if (tpcparam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcparam->GetPadRowRadii(36, 95)) {
        h->SetBinContent(ix, iy, dx[2]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcparam;
  return h;
}

///
/// \param z
/// \param nx
/// \param ny
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoCorrDZinXY(Float_t z, Int_t nx, Int_t ny) {
  /// Simple plot functionality.
  /// Returns a 2d hisogram which represents the corrections in longitudinal direction (dz)
  /// in respect to position z within the XY plane.
  /// The histogramm has nx times ny entries.

  AliTPCParam *tpcparam = new AliTPCParamSR;

  TH2F *h = CreateTH2F("dz_xy", TString::Format("%s: DZinXY Z=%2.0f", GetTitle(), z).Data(), "x [cm]", "y [cm]",
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
      if (tpcparam->GetPadRowRadii(0, 0) <= r0 && r0 <= tpcparam->GetPadRowRadii(36, 95)) {
        h->SetBinContent(ix, iy, dx[2]);
      } else
        h->SetBinContent(ix, iy, 0.);
    }
  }
  delete tpcparam;
  return h;
}


/// Use 3D space charge map as an optional input
/// The layout of the input histogram is assumed to be: (phi,r,z)
/// Density histogram is expreseed is expected to bin in  C/m^3
///
/// Standard histogram interpolation is used in order to use the density at center of voxel
/// Warning: Since  value at phi=0, is discontinued
///
/// \param hisSpaceCharge3D
/// \param norm
void AliTPCSpaceCharge3DDriftLine::SetInputSpaceCharge
        (
                TH3 *hisSpaceCharge3D,
                Double_t norm
        ) {
  fHistogram3DSpaceCharge = hisSpaceCharge3D;
  fInitLookUp = kFALSE;
}

///
/// \param hisSpaceCharge3D
/// \param norm
/// \param side
void AliTPCSpaceCharge3DDriftLine::SetInputSpaceCharge
        (
                TH3 *hisSpaceCharge3D,
                Double_t norm,
                Int_t side
        ) {

  // grid size for one side

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

    /// Fill the non-boundary values
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
TTree *AliTPCSpaceCharge3DDriftLine::CreateDistortionTree
        (
                Double_t step
        ) {


  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("distortion%s.root", GetName()));
  Float_t xyz[3];     // current point
  Float_t dist[3];    // distorion
  Float_t localDist[3];    // distorion
  Float_t corr[3];    // correction
  Float_t xyzdist[3]; // distorted point
  Float_t xyzcorr[3]; // corrected point

  //AliMagF* mag= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  //if (!mag) AliError("Magnetic field - not initialized");

  Int_t roc;
  AliTPCParam *tpcparam = new AliTPCParamSR;
  Double_t r, phi, rdist, phidist, drdist, drphidist, rcorr, phicorr, drcorr, drphicorr;
  for (Double_t x = -250; x < 250; x += step) {
    for (Double_t y = -250; y < 250; y += step) {

      r = TMath::Sqrt(x * x + y * y);

      if (tpcparam->GetPadRowRadii(0, 0) > r || r > tpcparam->GetPadRowRadii(36, 95)) continue;

      //printf("(%f,%f)\n",x,y);
      phi = TMath::ATan2(y, x);


      for (Double_t z = -250; z < 250; z += step) {
        roc = (z > 0) ? 0 : 18;
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;

        GetDistortion(xyz, roc, dist);

        for (Int_t i = 0; i < 3; ++i) {
          xyzdist[i] = xyz[i] + dist[i];
        }

        xyz[0] = r;
        xyz[1] = phi;
        xyz[2] = z;

        GetLocalDistortionCylAC(xyz, roc, localDist);

        GetCorrection(xyzdist, roc, corr);

        for (Int_t i = 0; i < 3; ++i) {
          xyzcorr[i] = xyzdist[i] + corr[i];
        }

        // === r, rphi + residuals for the distorted point =========================
        rdist = TMath::Sqrt(xyzdist[0] * xyzdist[0] + xyzdist[1] * xyzdist[1]);
        phidist = TMath::ATan2(xyzdist[1], xyzdist[0]);
        //rdist = xyzdist[0];
        //phidist = xyzdist[1];

        while ((phidist - phi) > TMath::Pi()) phidist -= TMath::TwoPi();
        while ((phidist - phi) < -TMath::Pi()) phidist += TMath::TwoPi();

        drdist = rdist - r;
        drphidist = (phidist - phi) * r;

        // === r, rphi + residuals for the corrected point =========================
        rcorr = TMath::Sqrt(xyzcorr[0] * xyzcorr[0] + xyzcorr[1] * xyzcorr[1]);
        phicorr = TMath::ATan2(xyzcorr[1], xyzcorr[0]);
        //rcorr = xyzcorr[0];
        //phicorr = xyzcorr[1];

        while ((phicorr - phidist) > TMath::Pi()) phicorr -= TMath::TwoPi();
        while ((phicorr - phidist) < -TMath::Pi()) phicorr += TMath::TwoPi();


        drcorr = rcorr - rdist;
        drphicorr = (phicorr - phidist) * r;
        (*pcstream) << "distortion" <<
                    "x=" << x <<           // original position
                    "y=" << y <<
                    "z=" << z <<
                    "r=" << r <<
                    "phi=" << phi <<
                    "xDist=" << xyzdist[0] <<      // distorted position
                    "yDist=" << xyzdist[1] <<
                    "zDist=" << xyzdist[2] <<
                    "rDist=" << rdist <<
                    "phiDist=" << phidist <<
                    "dxDist=" << dist[0] <<     // distortion
                    "dyDist=" << dist[1] <<
                    "dzDist=" << dist[2] <<
                    "drDist=" << drdist <<
                    "drPhiDist=" << drphidist <<
                    "drLocalDist=" << localDist[0] <<
                    "drPhiLocalDist=" << localDist[1] <<
                    "dzLocalDist=" << localDist[2] <<
                    "xCorr=" << xyzcorr[0] <<      // corrected position
                    "yCorr=" << xyzcorr[1] <<
                    "zCorr=" << xyzcorr[2] <<
                    "rCorr=" << rcorr <<
                    "phiCorr=" << phicorr <<
                    //
                    "dxCorr=" << corr[0] <<     // correction
                    "dyCorr=" << corr[1] <<
                    "dzCorr=" << corr[2] <<
                    "drCorr=" << drcorr <<
                    "drphiCorr=" << drphicorr <<
                    "\n";
      }
    }
  }
  delete pcstream;
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
    // 	printf("Phi --> Interpolate Cannot interpolate outside histogram domain. (%f,%f,%f)\n",phi,r,z);
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
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoSCinXY
        (
                Float_t z,
                Int_t nx,
                Int_t ny
        ) {
  /// return a simple histogramm containing the space charge distribution (input for the calculation)

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


/// return a simple histogramm containing the space charge distribution (input for the calculation)
///
/// \param phi
/// \param nz
/// \param nr
/// \return
TH2F *AliTPCSpaceCharge3DDriftLine::CreateHistoSCinZR
        (
                Float_t phi,
                Int_t nz,
                Int_t nr
        ) {


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
  // Float_t sc =fSCdensityDistribution->Interpolate(r0,phi0,z0);
  const Int_t order = 1; //

  const Float_t x[] = {r, phi, z};
  Float_t sc = 0;
  if (z > 0)
    sc = GetChargeCylAC(x, 0);
  else
    sc = GetChargeCylAC(x, 18);

  return sc;
}


/// Inverse from Global to Local Distortion
///
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDphiRDz TMatrixD **  matrix of global distortion dphir (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rlist Double_t * points of r in the grid (ascending mode)
/// \param zlist Double_t * points of z in the grid (ascending mode)
/// \param philist Double_t * points of phi in the grid (ascending mode)
/// \param nRRow Int_t number of grid in r direction
/// \param nZColumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///	\param nstep Int_t number of step to calculate local dist
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortion(TMatrixD **matricesDistDrDz,
                                                                  TMatrixD **matricesDistDphiRDz,
                                                                  TMatrixD **matricesDistDz, Double_t *rList,
                                                                  Double_t *zList, Double_t *phiList, const Int_t nRRow,
                                                                  const Int_t nZColumn, const Int_t phiSlice,
                                                                  const Int_t nstep, const Bool_t useCylAC, Int_t stepR,
                                                                  Int_t stepZ, Int_t stepPhi, Int_t type) {
  Double_t z, phi, r, zaft, zprev, zstep, ddr, ddphir, ddz, zl, dr, dphir, dz, ddphi, dphi, deltaz, zm1, zm2, zp1, zp2, gradient_r, gradient_phir, gradient_z, r0, z0, phi0;
  Float_t x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];
  Int_t roc;

  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;


  TMatrixD *distDrDz;
  TMatrixD *distDphiRDz;
  TMatrixD *distDz;

  // correction build up for inverse flow
  TMatrixD *corrDrDz;
  TMatrixD *corrDphiRDz;
  TMatrixD *corrDz;

  TMatrixD *listR;
  TMatrixD *listPhi;
  TMatrixD *listZ;

  /// make interpolator for inverse flow
  TMatrixD *matricesCorrDrDz[phiSlice];
  TMatrixD *matricesCorrDphiRDz[phiSlice];
  TMatrixD *matricesCorrDz[phiSlice];

  TMatrixD *matricesRList[phiSlice];
  TMatrixD *matricesPhiList[phiSlice];
  TMatrixD *matricesZList[phiSlice];

  for (Int_t m = 0; m < phiSlice; m++) {
    matricesCorrDrDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDphiRDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDz[m] = new TMatrixD(nRRow, nZColumn);

    matricesRList[m] = new TMatrixD(nRRow, nZColumn);
    matricesPhiList[m] = new TMatrixD(nRRow, nZColumn);
    matricesZList[m] = new TMatrixD(nRRow, nZColumn);
  }

  AliTPCLookUpTable3DInterpolatorDFull *lookupInverseCorr =
          new AliTPCLookUpTable3DInterpolatorDFull
                  (
                          nRRow,
                          matricesCorrDrDz,
                          matricesRList,
                          rList,
                          phiSlice,
                          matricesCorrDphiRDz,
                          matricesPhiList,
                          phiList,
                          nZColumn,
                          matricesCorrDz,
                          matricesZList,
                          zList,
                          2,
                          stepR,
                          stepZ,
                          stepPhi,
                          type
                  );


  for (Int_t k = 0; k < phiSlice; k++) {
    distDrDz = matricesDistDrDz[k];
    distDphiRDz = matricesDistDphiRDz[k];
    distDz = matricesDistDz[k];

    listR = matricesRList[k];
    listPhi = matricesPhiList[k];
    listZ = matricesZList[k];


    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDphiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

      for (Int_t j = 0; j < nZColumn; j++) {
        (*listR)(i, j) = rList[i];
        (*listPhi)(i, j) = phiList[k];
        (*listZ)(i, j) = zList[j];

      }

    }
  }


  deltaz = (zList[1] - zList[0]);
  Int_t ianchor, kanchor, zanchor;

  for (Int_t j = nZColumn - 2; j >= 0; j--) {

    printf("inversion global on z column  = %d\n", j);

    roc = 0; // FIXME
    for (Int_t k = 0; k < phiSlice; k++) {

      distDrDz = matricesDistDrDz[k];
      distDphiRDz = matricesDistDphiRDz[k];
      distDz = matricesDistDz[k];


      corrDrDz = matricesCorrDrDz[k];
      corrDphiRDz = matricesCorrDphiRDz[k];
      corrDz = matricesCorrDz[k];


      listR = matricesRList[k];
      listPhi = matricesPhiList[k];
      listZ = matricesZList[k];


      for (Int_t i = 0; i < nRRow; i++) {
        // get global distortion

        r = rList[i];
        phi = phiList[k];
        z = zList[j];
        zp1 = zList[j + 1];

        dr = 0.0;
        dz = 0.0;
        dphir = 0.0;

        if (j < nZColumn - 2) {
          // get global distortion of this point

          //printf("inversion global on z column  = %d\n",j);
          x[0] = r;
          x[1] = phi;
          x[2] = zList[j];
          if (useCylAC == kTRUE)
            GetDistortionCylAC(x, roc, dx);
          else
            GetDistortionCyl(x, roc, dx);



          // get position on 0 based on global distortion
          r0 = r + dx[0];
          z0 = zList[nZColumn - 1] + dx[2];
          phi0 = phi + (dx[1] / r);

          // follow electron path from correction table
          for (Int_t jj = nZColumn - 1; jj > j + 1; jj--) {

            //if (r0 < fgkIFCRadius)
            //  r0 = fgkIFCRadius;
            //if (r0 > fgkOFCRadius)
            //    printf("(%d,%d,%d,%f)\n",i,j,k,r0);
            //  r0 = fgkOFCRadius;

            if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
            if (phi0 > TMath::TwoPi()) phi0 = phi0 - TMath::TwoPi();

            ianchor = TMath::FloorNint((r0 - fgkIFCRadius) / gridSizeR);
            kanchor = TMath::FloorNint(phi0 / gridSizePhi);
            zanchor = TMath::FloorNint(z0 / gridSizeZ);

            if (j > nZColumn - 5)
              lookupInverseCorr->GetValue(r0, phi0, z0, dr, dphir, dz, ianchor, kanchor, zanchor,
                                          nRRow / 4 + 1, phiSlice / 4 + 1, 1, j + 2);
            else
              lookupInverseCorr->GetValue(r0, phi0, z0, dr, dphir, dz, ianchor, kanchor, zanchor,
                                          nRRow / 4 + 1, phiSlice / 4 + 1, 3, j + 2);

            phi0 = phi0 + ((dphir) / r0);
            r0 = r0 + (dr);

            z0 = z0 - deltaz + (dz);

          }

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
          dphir = (dx[1] - pdx[1]);

        } else if (j == (nZColumn - 2)) {

          // in case of nZColumn-2, the global distortion is the local distortion
          //printf("inversion global on z column  = %d\n",j);

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
          dphir = (dx[1] - pdx[1]);

        }

        (*distDrDz)(i, j) = dr;
        (*distDz)(i, j) = dz;
        (*distDphiRDz)(i, j) = dphir;

        r = rList[i];
        phi = phiList[k];
        z = zList[j];

        (*corrDrDz)(i, j + 1) = -dr;
        (*corrDz)(i, j + 1) = -dz;
        (*corrDphiRDz)(i, j + 1) = -dphir;

        (*listR)(i, j + 1) = r + dr;
        (*listPhi)(i, j + 1) = phi + dphir / r;
        (*listZ)(i, j + 1) = zp1 + dz;


      }


    }

    // copy Vals to correction look up table
    // in case RBF should compute weighted to the interpolant
    lookupInverseCorr->CopyVals(j + 1);


  }

  for (Int_t m = 0; m < phiSlice; m++) {
    delete matricesCorrDrDz[m];
    delete matricesCorrDphiRDz[m];
    delete matricesCorrDz[m];
    delete matricesRList[m];
    delete matricesPhiList[m];
    delete matricesZList[m];
  }
  delete lookupInverseCorr;

}


///
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDphiRDz TMatrixD **  matrix of global distortion dphir (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rlist Double_t * points of r in the grid (ascending mode)
/// \param zlist Double_t * points of z in the grid (ascending mode)
/// \param philist Double_t * points of phi in the grid (ascending mode)
/// \param nRRow Int_t number of grid in r direction
/// \param nZColumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///	\param nstep Int_t number of step to calculate local dist
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionTwoStages
        (
                TMatrixD **matricesDistDrDz,
                TMatrixD **matricesDistDphiRDz,
                TMatrixD **matricesDistDz,
                Double_t *rList,
                Double_t *zList,
                Double_t *phiList,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Int_t nstep,
                const Bool_t useCylAC,
                Int_t stepR,
                Int_t stepZ,
                Int_t stepPhi,
                Int_t type // 0 1 grid, 1 two stages
        ) {
  Double_t z, phi, r, zaft, zprev, zstep, ddr, ddphir, ddz, zl, dr, dphir, dz, ddphi, dphi, deltaz, zm1, zm2, zp1, zp2, gradient_r, gradient_phir, gradient_z, r0, z0, phi0;
  Float_t x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];


  Int_t roc;
  Int_t ianchor, kanchor, zanchor;


  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;


  TMatrixD *distDrDz;
  TMatrixD *distDphiRDz;
  TMatrixD *distDz;
  /// get





  // correction build up for inverse flow
  TMatrixD *corrDrDz;
  TMatrixD *corrDphiRDz;
  TMatrixD *corrDz;

  TMatrixD *listR;
  TMatrixD *listPhi;
  TMatrixD *listZ;


  /// make interpolator for inverse flow
  ///
  TMatrixD *matricesCorrDrDz[phiSlice];
  TMatrixD *matricesCorrDphiRDz[phiSlice];
  TMatrixD *matricesCorrDz[phiSlice];

  TMatrixD *matricesRList[phiSlice];
  TMatrixD *matricesPhiList[phiSlice];
  TMatrixD *matricesZList[phiSlice];

  for (Int_t m = 0; m < phiSlice; m++) {
    matricesCorrDrDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDphiRDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDz[m] = new TMatrixD(nRRow, nZColumn);

    matricesRList[m] = new TMatrixD(nRRow, nZColumn);
    matricesPhiList[m] = new TMatrixD(nRRow, nZColumn);
    matricesZList[m] = new TMatrixD(nRRow, nZColumn);
  }

  AliTPCLookUpTable3DInterpolatorDFull *lookupInverseCorr =
          new AliTPCLookUpTable3DInterpolatorDFull
                  (
                          nRRow,
                          matricesCorrDrDz,
                          matricesRList,
                          rList,
                          phiSlice,
                          matricesCorrDphiRDz,
                          matricesPhiList,
                          phiList,
                          nZColumn,
                          matricesCorrDz,
                          matricesZList,
                          zList,
                          2,
                          stepR,
                          stepZ,
                          stepPhi,
                          type
                  );


  for (Int_t k = 0; k < phiSlice; k++) {
    distDrDz = matricesDistDrDz[k];
    distDphiRDz = matricesDistDphiRDz[k];
    distDz = matricesDistDz[k];

    listR = matricesRList[k];
    listPhi = matricesPhiList[k];
    listZ = matricesZList[k];


    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDphiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

      for (Int_t j = 0; j < nZColumn; j++) {
        (*listR)(i, j) = rList[i];
        (*listPhi)(i, j) = phiList[k];
        (*listZ)(i, j) = zList[j];

      }

    }
  }


// Do for the lowest (17,17,18)

/////////////////// get local distortion table for coarser grid /////
  Int_t nRRowCoarser = (nRRow / 2) + 1;
  Int_t nZColumnCoarser = (nZColumn / 2) + 1;
  Int_t phiSliceCoarser = phiSlice / 2;
  Double_t *rListCoarser = new Double_t[nRRowCoarser];
  Double_t *zListCoarser = new Double_t[nZColumnCoarser];
  Double_t *phiListCoarser = new Double_t[phiSliceCoarser];

  // copy rList
  for (int i = 0; i < nRRowCoarser; i++) rListCoarser[i] = rList[i * 2];
  for (int i = 0; i < nZColumnCoarser; i++) zListCoarser[i] = zList[i * 2];
  for (int i = 0; i < phiSliceCoarser; i++) phiListCoarser[i] = phiList[i * 2];


  TMatrixD *matricesDistDrDzCoarser[phiSliceCoarser];
  TMatrixD *matricesDistDphiRDzCoarser[phiSliceCoarser];
  TMatrixD *matricesDistDzCoarser[phiSliceCoarser];


  for (Int_t m = 0; m < phiSliceCoarser; m++) {
    matricesDistDrDzCoarser[m] = new TMatrixD(nRRowCoarser, nZColumnCoarser);
    matricesDistDphiRDzCoarser[m] = new TMatrixD(nRRowCoarser, nZColumnCoarser);
    matricesDistDzCoarser[m] = new TMatrixD(nRRowCoarser, nZColumnCoarser);
  }


  // prepare interpolator
  AliTPCLookUpTable3DInterpolatorD *lookupLocalDistCoarser =
          new AliTPCLookUpTable3DInterpolatorD
                  (
                          nRRowCoarser,
                          matricesDistDrDzCoarser,
                          rListCoarser,
                          phiSliceCoarser,
                          matricesDistDphiRDzCoarser,
                          phiListCoarser,
                          nZColumnCoarser,
                          matricesDistDzCoarser,
                          zListCoarser,
                          fInterpolationOrder
                  );

  InverseGlobalToLocalDistortion
          (
                  matricesDistDrDzCoarser,
                  matricesDistDphiRDzCoarser,
                  matricesDistDzCoarser,
                  rListCoarser,
                  zListCoarser,
                  phiListCoarser,
                  nRRowCoarser,
                  nZColumnCoarser,
                  phiSliceCoarser,
                  nstep,
                  useCylAC,
                  stepR,
                  stepZ,
                  stepPhi,
                  type
          );

  // copy vals from mat** to mat* FIX
  lookupLocalDistCoarser->CopyVals();
  //

  /** For current grid we have two cases
	**/

  deltaz = (zList[1] - zList[0]);

  for (Int_t j = nZColumn - 2; j >= 0; j--) {

    printf("inversion global on z column  = %d\n", j);

    roc = 0; // FIXME
    for (Int_t k = 0; k < phiSlice; k++) {

      distDrDz = matricesDistDrDz[k];
      distDphiRDz = matricesDistDphiRDz[k];
      distDz = matricesDistDz[k];


      corrDrDz = matricesCorrDrDz[k];
      corrDphiRDz = matricesCorrDphiRDz[k];
      corrDz = matricesCorrDz[k];


      listR = matricesRList[k];
      listPhi = matricesPhiList[k];
      listZ = matricesZList[k];


      for (Int_t i = 0; i < nRRow; i++) {
        // get global distortion

        r = rList[i];
        phi = phiList[k];
        z = zList[j];
        zp1 = zList[j + 1];

        dr = 0.0;
        dz = 0.0;
        dphir = 0.0;

        if (j < nZColumn - 2) {
          // get global distortion of this point

          //printf("inversion global on z column  = %d\n",j);
          //				if (j  % 2 == 0) {
          // in this case we know local distortion from coarser grid
          // find local distortion at


          x[0] = r;
          x[1] = phi;
          x[2] = zList[j];
          if (useCylAC == kTRUE)
            GetDistortionCylAC(x, roc, dx);
          else
            GetDistortionCyl(x, roc, dx);


//						printf("original point (%f,%f,%f)\n",r,phi,z);

          lookupLocalDistCoarser->GetValue(r, phi, z, dr, dphir, dz);

          // distorted point
          r0 = r + dr;
          z0 = zList[j + 2] + dz;
          phi0 = phi + (dphir / r);

          // find local correction at this point
//						printf("distorted point (%f,%f,%f)\n",r0,phi0,z0);

          if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
          if (phi0 > TMath::TwoPi()) phi0 = phi0 - TMath::TwoPi();

          ianchor = TMath::FloorNint((r0 - fgkIFCRadius) / gridSizeR);
          kanchor = TMath::FloorNint(phi0 / gridSizePhi);
          zanchor = TMath::FloorNint(z0 / gridSizeZ);


          if (j > nZColumn - 5)
            lookupInverseCorr->GetValue(r0, phi0, z0, dr, dphir, dz, ianchor, kanchor, zanchor,
                                        nRRow / 4 + 1, phiSlice / 4 + 1, 1, j + 2);
          else
            lookupInverseCorr->GetValue(r0, phi0, z0, dr, dphir, dz, ianchor, kanchor, zanchor,
                                        nRRow / 4 + 1, phiSlice / 4 + 1, 3, j + 2);

          // get the last point to calculate local distortion
          phi0 = phi0 + ((dphir) / r0);
          r0 = r0 + (dr);
          z0 = z0 - deltaz + (dz);


          // find local correction at this point
          //					printf("last distorted point (%f,%f,%f)\n",r0,phi0,z0);

          // get final distortion
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
          dphir = (dx[1] - pdx[1]);

/**

					} else {


					    x[0] = r ;
					    x[1] = phi;
					    x[2] = zList[j];
					    if (useCylAC == kTRUE)
					      GetDistortionCylAC(x,roc,dx);
					    else
					      GetDistortionCyl(x,roc,dx);



					    // get position on 0 based on global distortion
					    r0    = r + dx[0];
					    z0    = zList[nZColumn - 1] + dx[2];
					    phi0  = phi + (dx[1]/r);

					    // follow electron path from correction table
					    for (Int_t jj=nZColumn-1;jj > j+1; jj--)
					    {


							if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
							if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

							ianchor = TMath::FloorNint((r0 - fgkIFCRadius)/gridSizeR);
							kanchor = TMath::FloorNint(phi0/gridSizePhi);
							zanchor = TMath::FloorNint(z0/gridSizeZ);

							//printf("call %f,%f,%f,%d,%d,%d\n",r0,phi0,z0,ianchor,kanchor,zanchor);

							if (j > nZColumn - 5)
							  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,nRRow/4  + 1, phiSlice/4 + 1,1,j+2);
							else
							  lookupInverseCorr->GetValue(r0,phi0,z0,dr,dphir,dz,ianchor,kanchor,zanchor,nRRow/4  + 1, phiSlice/4 + 1,3,j+2);

							phi0 = phi0 + ((dphir)/r0);
							r0 = r0 + ( dr );

							z0 = z0 - deltaz + (dz);

					  	}

				    	x[0] = r0 ;
				    	x[1] = phi0;
				    	x[2] = z0;

				    	if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
				    	if (phi0 > TMath::TwoPi())  phi0 = phi0 - TMath::TwoPi();

				    	if (useCylAC == kTRUE)
				      		GetDistortionCylAC(x,roc,pdx);
				    	else
				      		GetDistortionCyl(x,roc,pdx);

				    	dr   = (dx[0] - pdx[0]);
				    	dz   = (dx[2] - pdx[2]);
				    	dphir = (dx[1] - pdx[1]);
					}
**/
        } else if (j == (nZColumn - 2)) {

          // in case of nZColumn-2, the global distortion is the local distortion
          //printf("inversion global on z column  = %d\n",j);

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
          dphir = (dx[1] - pdx[1]);

        }

        (*distDrDz)(i, j) = dr;
        (*distDz)(i, j) = dz;
        (*distDphiRDz)(i, j) = dphir;

        r = rList[i];
        phi = phiList[k];
        z = zList[j];

        (*corrDrDz)(i, j + 1) = -dr;
        (*corrDz)(i, j + 1) = -dz;
        (*corrDphiRDz)(i, j + 1) = -dphir;

        (*listR)(i, j + 1) = r + dr;
        (*listPhi)(i, j + 1) = phi + dphir / r;
        (*listZ)(i, j + 1) = zp1 + dz;


      }


    }

    // copy Vals to correction look up table
    // in case RBF should compute weighted to the interpolant
    lookupInverseCorr->CopyVals(j + 1);


  }



  // dealocate memory
  delete lookupLocalDistCoarser;

  for (Int_t m = 0; m < phiSliceCoarser; m++) {
    delete matricesDistDrDzCoarser[m];
    delete matricesDistDphiRDzCoarser[m];
    delete matricesDistDzCoarser[m];
  }


  delete[] rListCoarser;
  delete[] zListCoarser;
  delete[] phiListCoarser;


  for (Int_t m = 0; m < phiSlice; m++) {
    delete matricesCorrDrDz[m];
    delete matricesCorrDphiRDz[m];
    delete matricesCorrDz[m];
    delete matricesRList[m];
    delete matricesPhiList[m];
    delete matricesZList[m];
  }
  delete lookupInverseCorr;

  /**
	delete[] matricesCorrDrDz;
	delete[] matricesCorrDphiRDz;
	delete[] matricesCorrDz;

	delete[] matricesRList;
	delete[] matricesPhiList;
	delete[] matricesZList;
	**/
}


///
///
/// \param matricesDistDrDz TMatrixD **  matrix of global distortion dr (r direction)
/// \param matricesDistDphiRDz TMatrixD **  matrix of global distortion dphir (phi r direction)
/// \param matricesDistDz TMatrixD **  matrix of global distortion dz (z direction)
/// \param rlist Double_t * points of r in the grid (ascending mode)
/// \param zlist Double_t * points of z in the grid (ascending mode)
/// \param philist Double_t * points of phi in the grid (ascending mode)
/// \param nRRow Int_t number of grid in r direction
/// \param nZColumn Int_t number of grid in z direction
/// \param phiSlice Int_t number of grid in phi direction
///	\param nstep Int_t number of step to calculate local dist
///
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionGlobalInvTable
        (
                TMatrixD **matricesDistDrDz,
                TMatrixD **matricesDistDphiRDz,
                TMatrixD **matricesDistDz,
                Double_t *rList,
                Double_t *zList,
                Double_t *phiList,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice,
                const Int_t nstep,
                const Bool_t useCylAC,
                Int_t stepR,
                Int_t stepZ,
                Int_t stepPhi,
                Int_t type // 0 1 grid, 1 two stages
        ) {

  Double_t z, phi, r, zaft, zprev, zstep, ddr, ddphir, ddz, zl, dr, dphir, dz, ddphi, dphi, deltaz, zm1, zm2, zp1, zp2, gradient_r, gradient_phir, gradient_z, r0, z0, phi0;
  Float_t x[3], dx[3], pdx[3], dxm2[3], dxm1[3], dxp1[3], dxp2[3];


  Int_t roc;

  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;


  TMatrixD *distDrDz;
  TMatrixD *distDphiRDz;
  TMatrixD *distDz;

  // correction build up for inverse flow
  TMatrixD *corrDrDz;
  TMatrixD *corrDphiRDz;
  TMatrixD *corrDz;

  TMatrixD *listR;
  TMatrixD *listPhi;
  TMatrixD *listZ;


  /// make interpolator for inverse flow
  ///
  TMatrixD *matricesCorrDrDz[phiSlice];
  TMatrixD *matricesCorrDphiRDz[phiSlice];
  TMatrixD *matricesCorrDz[phiSlice];

  TMatrixD *matricesRList[phiSlice];
  TMatrixD *matricesPhiList[phiSlice];
  TMatrixD *matricesZList[phiSlice];

  for (Int_t m = 0; m < phiSlice; m++) {
    matricesCorrDrDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDphiRDz[m] = new TMatrixD(nRRow, nZColumn);
    matricesCorrDz[m] = new TMatrixD(nRRow, nZColumn);

    matricesRList[m] = new TMatrixD(nRRow, nZColumn);
    matricesPhiList[m] = new TMatrixD(nRRow, nZColumn);
    matricesZList[m] = new TMatrixD(nRRow, nZColumn);
  }

  AliTPCLookUpTable3DInterpolatorDFull *lookupInverseCorr =
          new AliTPCLookUpTable3DInterpolatorDFull
                  (
                          nRRow,
                          matricesCorrDrDz,
                          matricesRList,
                          rList,
                          phiSlice,
                          matricesCorrDphiRDz,
                          matricesPhiList,
                          phiList,
                          nZColumn,
                          matricesCorrDz,
                          matricesZList,
                          zList,
                          2,
                          stepR,
                          stepZ,
                          stepPhi,
                          type
                  );
  lookupInverseCorr->SetKernelType(GetRBFKernelType());

  for (Int_t k = 0; k < phiSlice; k++) {
    distDrDz = matricesDistDrDz[k];
    distDphiRDz = matricesDistDphiRDz[k];
    distDz = matricesDistDz[k];

    listR = matricesRList[k];
    listPhi = matricesPhiList[k];
    listZ = matricesZList[k];


    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDphiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

      for (Int_t j = 0; j < nZColumn; j++) {
        (*listR)(i, j) = rList[i];
        (*listPhi)(i, j) = phiList[k];
        (*listZ)(i, j) = zList[j];

      }

    }
  }



  // 1) create global correction
  deltaz = (zList[1] - zList[0]);
  Int_t ianchor, kanchor, zanchor;

  for (Int_t j = nZColumn - 2; j >= 0; j--) {

//		printf("create inversion global inversion on z column  = %d\n",j);

    roc = 0; // FIXME
    for (Int_t k = 0; k < phiSlice; k++) {


      corrDrDz = matricesCorrDrDz[k];
      corrDphiRDz = matricesCorrDphiRDz[k];
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
        dphir = 0.0;


        // in case of nZColumn-2, the global distortion is the local distortion
        //printf("inversion global on z column  = %d\n",j);

        x[0] = r;
        x[1] = phi;
        x[2] = z;

        if (useCylAC == kTRUE)
          GetDistortionCylAC(x, roc, dx);
        else
          GetDistortionCyl(x, roc, dx);


        dr = dx[0];
        dz = dx[2];
        dphir = dx[1];
        //	printf("global distortion (%f,%f,%f) => (%f,%f,%f)\n",r,phi,z,dr,dphir,dz);

        r = rList[i];
        phi = phiList[k];
        z = zList[j];

        (*corrDrDz)(i, j + 1) = -dr;
        (*corrDz)(i, j + 1) = -dz;
        (*corrDphiRDz)(i, j + 1) = -dphir;

        (*listR)(i, j + 1) = r + dr;
        (*listPhi)(i, j + 1) = phi + dphir / r;
        (*listZ)(i, j + 1) = z + dz;
      }
    }
    lookupInverseCorr->CopyVals(j + 1);
  }
  // 2) calculate local distortion
  for (Int_t j = nZColumn - 2; j >= 0; j--) {
    roc = 0; // FIXME
    for (Int_t k = 0; k < phiSlice; k++) {
      distDrDz = matricesDistDrDz[k];
      distDphiRDz = matricesDistDphiRDz[k];
      distDz = matricesDistDz[k];
      for (Int_t i = 0; i < nRRow; i++) {
        // get global distortion
        r = rList[i];
        phi = phiList[k];
        z = zList[j];
        dr = 0.0;
        dz = 0.0;
        dphir = 0.0;

        if (j < nZColumn - 2) {
          // get global distortion of this point

          //printf("inversion global on z column  = %d\n",j);

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
          ianchor = TMath::FloorNint((r0 - fgkIFCRadius) / gridSizeR);
          kanchor = TMath::FloorNint(phi0 / gridSizePhi);
          zanchor = TMath::FloorNint(z0 / gridSizeZ);

          if (j > nZColumn - (GetIrregularGridSize() + 2))
            lookupInverseCorr->GetValue(r0, phi0, z0, dr, dphir, dz, ianchor, kanchor, zanchor,
                                        nRRow / 4 + 1, phiSlice / 4 + 1, 1, 0);
          else
            lookupInverseCorr->GetValue(r0, phi0, z0, dr, dphir, dz, ianchor, kanchor, zanchor,
                                        nRRow / 4 + 1, phiSlice / 4 + 1, GetIrregularGridSize(), 0);

          phi0 = phi0 + ((dphir) / r0);
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
          dphir = (dx[1] - pdx[1]);

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
          dphir = (dx[1] - pdx[1]);
        }

        (*distDrDz)(i, j) = dr;
        (*distDz)(i, j) = dz;
        (*distDphiRDz)(i, j) = dphir;
      }
    }
  }

  for (Int_t m = 0; m < phiSlice; m++) {
    delete matricesCorrDrDz[m];
    delete matricesCorrDphiRDz[m];
    delete matricesCorrDz[m];
    delete matricesRList[m];
    delete matricesPhiList[m];
    delete matricesZList[m];
  }
  delete lookupInverseCorr;
}

///
/// \param matricesEr
/// \param matricesEphi
/// \param matricesEz
/// \param matricesInvLocalIntErDz
/// \param matricesInvLocalIntEphiDz
/// \param matricesInvLocalIntEz
/// \param matricesDistDrDz
/// \param matricesDistDphiRDz
/// \param matricesDistDz
/// \param rList
/// \param zList
/// \param phiList
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseLocalDistortionToElectricField
        (
                TMatrixD **matricesEr, TMatrixD **matricesEphi,  TMatrixD **matricesEz,
                TMatrixD **matricesInvLocalIntErDz, TMatrixD **matricesInvLocalIntEphiDz, TMatrixD **matricesInvLocalIntEz,
                TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz, TMatrixD **matricesDistDz,
                Double_t *rList, Double_t *zList, Double_t *phiList,
                const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice
        ) {
  // calculate integral
  Float_t localIntErOverEz, localIntEphiOverEz, localIntDeltaEz, z2, gridzsize;
  Double_t r;
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;

  TMatrixD *distDrDz;
  TMatrixD *distDz;
  TMatrixD *distDphiRDz;
  TMatrixD *tdistDz;
  TMatrixD *tdistDphiRDz;
  TMatrixD *tdistDrDz;

  Float_t c02c12 = fC0 * fC0 + fC1 * fC1;

  // solve local integration
  for (Int_t j = 0; j < nZColumn; j++) {
    for (Int_t k = 0; k < phiSlice; k++) {
      distDrDz = matricesDistDrDz[k];
      distDz = matricesDistDz[k];
      distDphiRDz = matricesDistDphiRDz[k];

      tdistDrDz = matricesInvLocalIntErDz[k];
      tdistDz = matricesInvLocalIntEz[k];
      tdistDphiRDz = matricesInvLocalIntEphiDz[k];
      for (Int_t i = 0; i < nRRow; i++) {
        localIntErOverEz = fC0 * (*distDrDz)(i, j) - fC1 * (*distDphiRDz)(i, j);
        localIntErOverEz = localIntErOverEz / (fC0 * fC0 + fC1 * fC1);
        localIntEphiOverEz = ((*distDrDz)(i, j) - (fC0 * localIntErOverEz)) / fC1;
        localIntDeltaEz = (*distDz)(i, j) / (fgkdvdE * fgkdvdE); // two times?
        (*tdistDrDz)(i, j) = localIntErOverEz;
        (*tdistDphiRDz)(i, j) = localIntEphiOverEz;
        (*tdistDz)(i, j) = localIntDeltaEz;
      }
    }
  }
  TMatrixD *mEphi;
  TMatrixD *mEr;
  TMatrixD *mEz;

  // use central-backward-forward difference for calculating Electric field component
  for (Int_t m = 0; m < phiSlice; m++) {
    mEphi = matricesEphi[m];
    mEr = matricesEr[m];
    mEz = matricesEz[m];
    distDrDz = matricesInvLocalIntErDz[m];
    distDphiRDz = matricesInvLocalIntEphiDz[m];
    distDz = matricesInvLocalIntEz[m];
    for (Int_t i = 0; i < nRRow; i++) {
      (*mEr)(i, 0) = ((*distDrDz)(i, 0) / gridSizeZ) * -1 * ezField;
      (*mEphi)(i, 0) = ((*distDphiRDz)(i, 0) / gridSizeZ) * -1 * ezField;
      (*mEz)(i, 0) = ((*distDz)(i, 0) / gridSizeZ);
      (*mEr)(i, nZColumn - 1) =
              ((-0.5 * (*distDrDz)(i, nZColumn - 3) + 1.5 * (*distDrDz)(i, nZColumn - 2)) / gridSizeZ) * -1 *
              ezField;
      (*mEphi)(i, nZColumn - 1) =
              ((-0.5 * (*distDphiRDz)(i, nZColumn - 3) + 1.5 * (*distDphiRDz)(i, nZColumn - 2)) / gridSizeZ) * -1 *
              ezField;
      (*mEz)(i, nZColumn - 1) = (-0.5 * (*distDz)(i, nZColumn - 3) + 1.5 * (*distDz)(i, nZColumn - 2)) / gridSizeZ;
    }

    for (Int_t i = 0; i < nRRow; i++) {
      for (Int_t j = 1; j < nZColumn - 1; j++) {
        (*mEr)(i, j) = (((*distDrDz)(i, j) + (*distDrDz)(i, j - 1)) / (2 * gridSizeZ)) * -1 *
                       ezField; // z direction
        (*mEphi)(i, j) = (((*distDphiRDz)(i, j) + (*distDphiRDz)(i, j - 1)) / (2 * gridSizeZ)) * -1 *
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
/// \param matricesEphi
/// \param matricesEz
/// \param rList
/// \param zList
/// \param phiList
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseElectricFieldToCharge
        ( TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEphi, TMatrixD **matricesEz,
          Double_t *rList, Double_t *zList, Double_t *phiList, const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice) {

  Float_t radius;
  Double_t dr, dz, dphi;

  Int_t mplus, mminus, mplus2, mminus2, signplus, signminus;
  Int_t symmetry = 0;
  // iterate over phislices

  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;

  for (Int_t m = 0; m < phiSlice; m++) {
    mplus = m + 1;
    signplus = 1;
    mminus = m - 1;
    signminus = 1;
    mplus2 = m + 2;
    mminus2 = m - 2;
    if (symmetry == 1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if (mplus > phiSlice - 1) mplus = phiSlice - 2;
      if (mminus < 0) mminus = 1;
    } else if (symmetry == -1) {       // Anti-symmetry in phi
      if (mplus > phiSlice - 1) {
        mplus = phiSlice - 2;
        signplus = -1;
      }
      if (mminus < 0) {
        mminus = 1;
        signminus = -1;
      }
    } else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if (mplus > phiSlice - 1) mplus = m + 1 - phiSlice;
      if (mminus < 0) mminus = m - 1 + phiSlice;
      if (mplus2 > phiSlice - 1) mplus2 = m + 2 - phiSlice;
      if (mminus2 < 0) mminus2 = m - 2 + phiSlice;
    }

    TMatrixD &arrayCharge = *matricesCharge[m];
    TMatrixD &arrayEr = *matricesEr[m];
    TMatrixD &arrayEz = *matricesEz[m];
    TMatrixD &arrayEphi = *matricesEphi[m];
    TMatrixD &arrayEphiM = *matricesEphi[mminus];
    TMatrixD &arrayEphiP = *matricesEphi[mplus];
    TMatrixD &arrayEphiM2 = *matricesEphi[mminus2];
    TMatrixD &arrayEphiP2 = *matricesEphi[mplus2];


    // for non-boundary V
    for (Int_t i = 2; i < nRRow - 2; i++) {
      radius = fgkIFCRadius + i * gridSizeR;
      for (Int_t j = 2; j < nZColumn - 2; j++) {
        dr = (-arrayEr(i + 2, j) + 8 * arrayEr(i + 1, j) - 8 * arrayEr(i - 1, j) + arrayEr(i - 2, j)) /
             (12 * gridSizeR); // r direction
        dz = (-arrayEz(i, j + 2) + 8 * arrayEz(i, j + 1) - 8 * arrayEz(i, j - 1) + arrayEz(i, j - 2)) /
             (12 * gridSizeZ); // r direction
        dphi = (-arrayEphiP2(i, j) + 8 * arrayEphiP(i, j) - 8 * arrayEphiM(i, j) + arrayEphiM2(i, j)) /
               (12 * gridSizePhi); // phi// didrection

        arrayCharge(i, j) = -1 * (arrayEr(i, j) / radius + dr + dphi / radius + dz);
      }
    }

    // for boundary in r
    for (Int_t j = 2; j < nZColumn - 2; j++) {

      // r near inner radius
      // for index r[0]
      radius = fgkIFCRadius;
      dr = (-(11.0 / 6.0) * arrayEr(0, j) + (3.0 * arrayEr(1, j)) - (1.5 * arrayEr(2, j)) +
            ((1.0 / 3.0) * arrayEr(3, j))) / gridSizeR; // forward difference

      //	dr 	=  ( -(1.5)*arrayEr(0,j) + (2.0*arrayEr(1,j)) - (0.5*arrayEr(2,j)) )  / gridSizeR;

      dz = (-arrayEz(0, j + 2) + 8 * arrayEz(0, j + 1) - 8 * arrayEz(0, j - 1) + arrayEz(0, j - 2)) /
           (12.0 * gridSizeZ);; // z direction
      dphi = (-arrayEphiP2(0, j) + 8 * arrayEphiP(0, j) - 8 * arrayEphiM(0, j) + arrayEphiM2(0, j)) /
             (12.0 * gridSizePhi);

      arrayCharge(0, j) = -1 * (arrayEr(0, j) / radius + dr + dphi / radius + dz);


      // index use central difference 3-point center
      radius = fgkIFCRadius + gridSizeR;
      //	dr 	=  (-arrayEr(3,j)  +6.0*arrayEr(2,j) - 3.0*arrayEr(1,j) - 2*arrayEr(0,j) ) / (6.0*gridSizeR) ; // forward difference
      dr = (arrayEr(2, j) - arrayEr(0, j)) / (2.0 * gridSizeR);

      dz = (-arrayEz(1, j + 2) + 8 * arrayEz(1, j + 1) - 8 * arrayEz(1, j - 1) + arrayEz(1, j - 2)) /
           (12 * gridSizeZ);   // z direction
      dphi = (-arrayEphiP2(1, j) + 8 * arrayEphiP(1, j) - 8 * arrayEphiM(1, j) + arrayEphiM2(1, j)) /
             (12 * gridSizePhi);
      arrayCharge(1, j) = -1 * (arrayEr(1, j) / radius + dr + dphi / radius + dz);


      // index use centra difference 3-pont center
      radius = fgkIFCRadius + (nRRow - 2) * gridSizeR;
      //	dr =   (2.0 * arrayEr(nRRow - 1,j)  + 3.0*arrayEr(nRRow - 2,j) - 6.0*arrayEr(nRRow -3,j) + arrayEr(nRRow-4,j) ) / (6.0*gridSizeR) ;
      dr = (arrayEr(nRRow - 1, j) - arrayEr(nRRow - 3, j)) / (2.0 * gridSizeR);

      dz = (-arrayEz(nRRow - 2, j + 2) + 8 * arrayEz(nRRow - 2, j + 1) - 8 * arrayEz(nRRow - 2, j - 1) +
            arrayEz(nRRow - 2, j - 2)) / (12 * gridSizeZ);
      dphi = (-arrayEphiP2(nRRow - 2, j) + 8 * arrayEphiP(nRRow - 2, j) - 8 * arrayEphiM(nRRow - 2, j) +
              arrayEphiM2(nRRow - 2, j)) / (12.0 * gridSizePhi);
      arrayCharge(nRRow - 2, j) = -1 * (arrayEr(nRRow - 2, j) / radius + dr + dphi / radius + dz);

      // index r[nRRow -1] backward difference
      radius = fgkIFCRadius + (nRRow - 1) * gridSizeR;
      //dr =  ( 1.5*arrayEr(nRRow-1,j) - 2.0*arrayEr(nRRow-2,j) + 0.5*arrayEr(nRRow-3,j) ) / gridSizeR ; // backward difference
      dr = (-(11.0 / 6.0) * arrayEr(nRRow - 1, j) + (3.0 * arrayEr(nRRow - 2, j)) - (1.5 * arrayEr(nRRow - 3, j)) +
            ((1.0 / 3.0) * arrayEr(nRRow - 4, j))) / (-1 * gridSizeR);

      //dz    =  ( arrayEz(nRRow-1,j+1) - arrayEz(nRRow-1,j-1) ) / (2*gridSizeZ) ; // z direction
      dz = (-arrayEz(nRRow - 1, j + 2) + 8 * arrayEz(nRRow - 1, j + 1) - 8 * arrayEz(nRRow - 1, j - 1) +
            arrayEz(nRRow - 1, j - 2)) / (12 * gridSizeZ);


      dphi = (-arrayEphiP2(nRRow - 1, j) + 8 * arrayEphiP(nRRow - 1, j) - 8 * arrayEphiM(nRRow - 1, j) +
              arrayEphiM2(nRRow - 1, j)) / (12 * gridSizePhi);
      arrayCharge(nRRow - 1, j) = -1 * (arrayEr(nRRow - 1, j) / radius + dr + dphi / radius + dz);
    }

    // boundary z
    for (Int_t i = 2; i < nRRow - 2; i++) {
      // z[0]
      radius = fgkIFCRadius + i * gridSizeR;
      dz = (-(11.0 / 6.0) * arrayEz(i, 0) + (3.0 * arrayEz(i, 1)) - (1.5 * arrayEz(i, 2)) +
            ((1.0 / 3.0) * arrayEz(i, 3))) / (1 * gridSizeZ); // forward difference
      dr = (-arrayEr(i + 2, 0) + 8 * arrayEr(i + 1, 0) - 8 * arrayEr(i - 1, 0) + arrayEr(i - 2, 0)) /
           (12 * gridSizeR);; // z direction
      dphi = (-arrayEphiP2(i, 0) + 8 * arrayEphiP(i, 0) - 8 * arrayEphiM(i, 0) + arrayEphiM2(i, 0)) /
             (12 * gridSizePhi);
      arrayCharge(i, 0) = -1 * (arrayEr(i, 0) / radius + dr + dphi / radius + dz);

      dz = (arrayEz(i, 2) - arrayEz(i, 0)) / (2.0 * gridSizeZ); // forward difference

      dr = (-arrayEr(i + 2, 1) + 8 * arrayEr(i + 1, 1) - 8 * arrayEr(i - 1, 1) + arrayEr(i - 2, 1)) /
           (12 * gridSizeR);; // z direction
      dphi = (-arrayEphiP2(i, 1) + 8 * arrayEphiP(i, 1) - 8 * arrayEphiM(i, 1) + arrayEphiM2(i, 1)) /
             (12 * gridSizePhi);
      arrayCharge(i, 1) = -1 * (arrayEr(i, 1) / radius + dr + dphi / radius + dz);

      dz = (arrayEz(i, nZColumn - 1) - arrayEz(i, nZColumn - 3)) / (2.0 * gridSizeZ); // forward difference

      dr = (-arrayEr(i + 2, nZColumn - 2) + 8 * arrayEr(i + 1, nZColumn - 2) - 8 * arrayEr(i - 1, nZColumn - 2) +
            arrayEr(i - 2, nZColumn - 2)) / (12 * gridSizeR);; // z direction
      dphi = (-arrayEphiP2(i, nZColumn - 2) + 8 * arrayEphiP(i, nZColumn - 2) - 8 * arrayEphiM(i, nZColumn - 2) +
              arrayEphiM2(i, nZColumn - 2)) / (12 * gridSizePhi);
      arrayCharge(i, nZColumn - 2) = -1 * (arrayEr(i, nZColumn - 2) / radius + dr + dphi / radius + dz);

      dz = (-(11.0 / 6.0) * arrayEz(i, nZColumn - 1) + (3.0 * arrayEz(i, nZColumn - 2)) -
            (1.5 * arrayEz(i, nZColumn - 3)) + ((1.0 / 3.0) * arrayEz(i, nZColumn - 4))) /
           (-gridSizeZ); // backward difference
      dr = (-arrayEr(i + 2, nZColumn - 1) + 8 * arrayEr(i + 1, nZColumn - 1) - 8 * arrayEr(i - 1, nZColumn - 1) +
            arrayEr(i - 2, nZColumn - 1)) / (12 * gridSizeR);; // z direction
      dphi = (-arrayEphiP2(i, nZColumn - 1) + 8 * arrayEphiP(i, nZColumn - 1) - 8 * arrayEphiM(i, nZColumn - 1) +
              arrayEphiM2(i, nZColumn - 1)) / (12 * gridSizePhi);
      arrayCharge(i, nZColumn - 1) = -1 * (arrayEr(i, nZColumn - 1) / radius + dr + dphi / radius + dz);
    }
    // for corner points
    // corner points for Ephi
    radius = fgkIFCRadius;
    dr = (-0.5 * arrayEr(2, 0) + 2.0 * arrayEr(1, 0) - 1.5 * arrayEr(0, 0)) / gridSizeR; // forward difference
    dz = (-0.5 * arrayEz(0, 2) + 2.0 * arrayEz(0, 1) - 1.5 * arrayEz(0, 0)) / gridSizeZ; // forward difference
    dphi = (-arrayEphiP2(0, 0) + 8 * arrayEphiP(0, 0) - 8 * arrayEphiM(0, 0) + arrayEphiM2(0, 0)) /
           (12 * gridSizePhi);
    arrayCharge(0, 0) = -1 * (arrayEr(0, 0) / radius + dr + dphi / radius + dz);
    dr = (-0.5 * arrayEr(2, 1) + 2.0 * arrayEr(1, 1) - 1.5 * arrayEr(0, 1)) / gridSizeR; // forward difference
    dz = (arrayEz(0, 2) - arrayEz(0, 0)) / (2.0 * gridSizeZ); // forward difference
    dphi = (-arrayEphiP2(0, 1) + 8 * arrayEphiP(0, 1) - 8 * arrayEphiM(0, 1) + arrayEphiM2(0, 1)) /
           (12 * gridSizePhi);
    arrayCharge(0, 1) = -1 * (arrayEr(0, 1) / radius + dr + dphi / radius + dz);
    dr = (-0.5 * arrayEr(2, nZColumn - 2) + 2.0 * arrayEr(1, nZColumn - 2) - 1.5 * arrayEr(0, nZColumn - 2)) /
         gridSizeR; // forward difference
    dz = (2.0 * arrayEz(0, nZColumn - 1) + 3.0 * arrayEz(0, nZColumn - 2) - 6.0 * arrayEz(0, nZColumn - 3) +
          arrayEz(0, nZColumn - 4)) / (6.0 * gridSizeZ); // backward difference
    dphi = (-arrayEphiP2(0, nZColumn - 2) + 8 * arrayEphiP(0, nZColumn - 2) - 8 * arrayEphiM(0, nZColumn - 2) +
            arrayEphiM2(0, nZColumn - 2)) / (12 * gridSizePhi);
    arrayCharge(0, nZColumn - 2) = -1 * (arrayEr(0, nZColumn - 2) / radius + dr + dphi / radius + dz);
    dr = (-0.5 * arrayEr(2, nZColumn - 1) + 2.0 * arrayEr(1, nZColumn - 1) - 1.5 * arrayEr(0, nZColumn - 1)) /
         gridSizeR; // forward difference
    dz = (1.5 * arrayEz(0, nZColumn - 1) - 2.0 * arrayEz(0, nZColumn - 2) + 0.5 * arrayEz(0, nZColumn - 3)) /
         gridSizeZ; // backward difference
    dphi = (-arrayEphiP2(0, nZColumn - 1) + 8 * arrayEphiP(0, nZColumn - 1) - 8 * arrayEphiM(0, nZColumn - 1) +
            arrayEphiM2(0, nZColumn - 1)) / (12 * gridSizePhi);
    arrayCharge(0, nZColumn - 1) = -1 * (arrayEr(0, nZColumn - 1) / radius + dr + dphi / radius + dz);

    radius = fgkIFCRadius + gridSizeR;
    dr = (-arrayEr(3, 0) + 6.0 * arrayEr(2, 0) - 3.0 * arrayEr(1, 0) - 2 * arrayEr(0, 0)) /
         (6.0 * gridSizeR); // forward difference
    dz = (-0.5 * arrayEz(1, 2) + 2.0 * arrayEz(1, 1) - 1.5 * arrayEz(1, 0)) / gridSizeZ; // forward difference
    dphi = (-arrayEphiP2(1, 0) + 8 * arrayEphiP(1, 0) - 8 * arrayEphiM(1, 0) + arrayEphiM2(1, 0)) /
           (12 * gridSizePhi);
    arrayCharge(1, 0) = -1 * (arrayEr(1, 0) / radius + dr + dphi / radius + dz);
    dr = (-arrayEr(3, 1) + 6.0 * arrayEr(2, 1) - 3.0 * arrayEr(1, 1) - 2 * arrayEr(0, 1)) /
         (6.0 * gridSizeR); // forward difference
    dz = (-arrayEz(1, 3) + 6.0 * arrayEz(1, 2) - 3.0 * arrayEz(1, 1) - 2 * arrayEz(1, 0)) /
         (6.0 * gridSizeZ); // forward difference
    dphi = (-arrayEphiP2(1, 1) + 8 * arrayEphiP(1, 1) - 8 * arrayEphiM(1, 1) + arrayEphiM2(1, 1)) /
           (12 * gridSizePhi);
    arrayCharge(1, 1) = -1 * (arrayEr(1, 1) / radius + dr + dphi / radius + dz);
    dr = (-arrayEr(3, nZColumn - 2) + 6.0 * arrayEr(2, nZColumn - 2) - 3.0 * arrayEr(1, nZColumn - 2) -
          2 * arrayEr(0, nZColumn - 2)) / (6.0 * gridSizeR); // forward difference
    dz = (2.0 * arrayEz(1, nZColumn - 1) + 3.0 * arrayEz(1, nZColumn - 2) - 6.0 * arrayEz(1, nZColumn - 3) +
          arrayEz(1, nZColumn - 4)) / (6.0 * gridSizeZ); // backward difference
    dphi = (-arrayEphiP2(1, nZColumn - 2) + 8 * arrayEphiP(1, nZColumn - 2) - 8 * arrayEphiM(1, nZColumn - 2) +
            arrayEphiM2(1, nZColumn - 2)) / (12 * gridSizePhi);
    arrayCharge(1, nZColumn - 2) = -1 * (arrayEr(1, nZColumn - 2) / radius + dr + dphi / radius + dz);

    dr = (-arrayEr(3, nZColumn - 1) + 6.0 * arrayEr(2, nZColumn - 1) - 3.0 * arrayEr(1, nZColumn - 1) -
          2 * arrayEr(0, nZColumn - 1)) / (6.0 * gridSizeR); // forward difference
    dz = (1.5 * arrayEz(1, nZColumn - 1) - 2.0 * arrayEz(1, nZColumn - 2) + 0.5 * arrayEz(1, nZColumn - 3)) /
         gridSizeZ; // backward difference
    dphi = (-arrayEphiP2(1, nZColumn - 1) + 8 * arrayEphiP(1, nZColumn - 1) - 8 * arrayEphiM(1, nZColumn - 1) +
            arrayEphiM2(1, nZColumn - 1)) / (12 * gridSizePhi);
    arrayCharge(1, nZColumn - 1) = -1 * (arrayEr(1, nZColumn - 1) / radius + dr + dphi / radius + dz);

    radius = fgkIFCRadius + (nRRow - 2) * gridSizeR;
    dr = (2.0 * arrayEr(nRRow - 1, 0) + 3.0 * arrayEr(nRRow - 2, 0) - 6.0 * arrayEr(nRRow - 3, 0) +
          arrayEr(nRRow - 4, 0)) / (6.0 * gridSizeR); // backward difference
    dz = (-0.5 * arrayEz(nRRow - 2, 2) + 2.0 * arrayEz(nRRow - 2, 1) - 1.5 * arrayEz(nRRow - 2, 0)) /
         gridSizeZ; // forward difference
    dphi = (-arrayEphiP2(nRRow - 2, 0) + 8 * arrayEphiP(nRRow - 2, 0) - 8 * arrayEphiM(nRRow - 2, 0) +
            arrayEphiM2(nRRow - 2, 0)) / (12 * gridSizePhi);

    arrayCharge(nRRow - 2, 0) = -1 * (arrayEr(nRRow - 2, 0) / radius + dr + dphi / radius + dz);
    dr = (2.0 * arrayEr(nRRow - 1, 1) + 3.0 * arrayEr(nRRow - 2, 1) - 6.0 * arrayEr(nRRow - 3, 1) +
          arrayEr(nRRow - 4, 1)) / (6.0 * gridSizeR); // backward difference
    dz = (-arrayEz(nRRow - 2, 3) + 6.0 * arrayEz(nRRow - 2, 2) - 3.0 * arrayEz(nRRow - 2, 1) -
          2 * arrayEz(nRRow - 2, 0)) / (6.0 * gridSizeZ); // forward difference
    dphi = (-arrayEphiP2(nRRow - 2, 1) + 8 * arrayEphiP(nRRow - 2, 1) - 8 * arrayEphiM(nRRow - 2, 1) +
            arrayEphiM2(nRRow - 2, 1)) / (12 * gridSizePhi);
    arrayCharge(nRRow - 2, 1) = -1 * (arrayEr(nRRow - 2, 1) / radius + dr + dphi / radius + dz);
    dr = (2.0 * arrayEr(nRRow - 1, nZColumn - 2) + 3.0 * arrayEr(nRRow - 2, nZColumn - 2) -
          6.0 * arrayEr(nRRow - 3, nZColumn - 2) + arrayEr(nRRow - 4, nZColumn - 2)) /
         (6.0 * gridSizeR); // backward difference
    dz = (2.0 * arrayEz(nRRow - 2, nZColumn - 1) + 3.0 * arrayEz(nRRow - 2, nZColumn - 2) -
          6.0 * arrayEz(nRRow - 2, nZColumn - 3) + arrayEz(nRRow - 2, nZColumn - 4)) /
         (6.0 * gridSizeZ); // backward difference
    dphi = (-arrayEphiP2(nRRow - 2, nZColumn - 2) + 8 * arrayEphiP(nRRow - 2, nZColumn - 2) -
            8 * arrayEphiM(nRRow - 2, nZColumn - 2) + arrayEphiM2(nRRow - 2, nZColumn - 2)) / (12 * gridSizePhi);
    arrayCharge(nRRow - 2, nZColumn - 2) = -1 * (arrayEr(nRRow - 2, nZColumn - 2) / radius + dr + dphi / radius + dz);
    dr = (2.0 * arrayEr(nRRow - 1, nZColumn - 1) + 3.0 * arrayEr(nRRow - 2, nZColumn - 1) -
          6.0 * arrayEr(nRRow - 3, nZColumn - 1) + arrayEr(nRRow - 4, nZColumn - 1)) /
         (6.0 * gridSizeR); // backward difference
    dz = (1.5 * arrayEz(0, nZColumn - 1) - 2.0 * arrayEz(0, nZColumn - 2) + 0.5 * arrayEz(0, nZColumn - 3)) /
         gridSizeZ; // backward difference
    dphi = (-arrayEphiP2(nRRow - 2, nZColumn - 1) + 8 * arrayEphiP(nRRow - 2, nZColumn - 1) -
            8 * arrayEphiM(nRRow - 2, nZColumn - 1) + arrayEphiM2(nRRow - 2, nZColumn - 1)) / (12 * gridSizePhi);

    arrayCharge(nRRow - 2, nZColumn - 1) = -1 * (arrayEr(nRRow - 2, nZColumn - 1) / radius + dr + dphi / radius + dz);
    radius = fgkIFCRadius + (nRRow - 1) * gridSizeR;
    dr = (1.5 * arrayEr(nRRow - 1, 0) - 2.0 * arrayEr(nRRow - 2, 0) + 0.5 * arrayEr(nRRow - 3, 0)) /
         gridSizeR; // backward difference
    dz = (-0.5 * arrayEz(nRRow - 1, 2) + 2.0 * arrayEz(nRRow - 1, 1) - 1.5 * arrayEz(nRRow - 1, 0)) /
         gridSizeZ; // forward difference
    dphi = (-arrayEphiP2(nRRow - 1, 0) + 8 * arrayEphiP(nRRow - 1, 0) - 8 * arrayEphiM(nRRow - 1, 0) +
            arrayEphiM2(nRRow - 1, 0)) / (12 * gridSizePhi);
    arrayCharge(nRRow - 1, 0) = -1 * (arrayEr(nRRow - 1, 0) / radius + dr + dphi / radius + dz);
    dr = (1.5 * arrayEr(nRRow - 1, 1) - 2.0 * arrayEr(nRRow - 2, 1) + 0.5 * arrayEr(nRRow - 3, 1)) /
         gridSizeR; // backward difference
    dz = (-arrayEz(nRRow - 1, 3) + 6.0 * arrayEz(nRRow - 1, 2) - 3.0 * arrayEz(nRRow - 1, 1) -
          2 * arrayEz(nRRow - 1, 0)) / (6.0 * gridSizeZ); // forward difference
    dphi = (-arrayEphiP2(nRRow - 1, 1) + 8 * arrayEphiP(nRRow - 1, 1) - 8 * arrayEphiM(nRRow - 1, 1) +
            arrayEphiM2(nRRow - 1, 1)) / (12 * gridSizePhi);
    arrayCharge(nRRow - 1, 1) = -1 * (arrayEr(nRRow - 1, 1) / radius + dr + dphi / radius + dz);

    dr = (1.5 * arrayEr(nRRow - 1, nZColumn - 2) - 2.0 * arrayEr(nRRow - 2, nZColumn - 2) +
          0.5 * arrayEr(nRRow - 3, nZColumn - 2)) / gridSizeR; // backward difference
    dz = (2.0 * arrayEz(nRRow - 1, nZColumn - 1) + 3.0 * arrayEz(nRRow - 1, nZColumn - 2) -
          6.0 * arrayEz(nRRow - 1, nZColumn - 3) + arrayEz(nRRow - 1, nZColumn - 4)) /
         (6.0 * gridSizeZ); // backward difference
    dphi = (-arrayEphiP2(nRRow - 1, nZColumn - 2) + 8 * arrayEphiP(nRRow - 1, nZColumn - 2) -
            8 * arrayEphiM(nRRow - 1, nZColumn - 2) + arrayEphiM2(nRRow - 1, nZColumn - 2)) / (12 * gridSizePhi);
    arrayCharge(nRRow - 1, nZColumn - 2) = -1 * (arrayEr(nRRow - 1, nZColumn - 2) / radius + dr + dphi / radius + dz);


    dr = (1.5 * arrayEr(nRRow - 1, nZColumn - 1) - 2.0 * arrayEr(nRRow - 2, nZColumn - 1) +
          0.5 * arrayEr(nRRow - 3, nZColumn - 1)) / gridSizeR; // backward difference
    dz = (1.5 * arrayEz(nRRow - 1, nZColumn - 1) - 2.0 * arrayEz(nRRow - 1, nZColumn - 2) +
          0.5 * arrayEz(nRRow - 1, nZColumn - 3)) / gridSizeZ; // backward difference

    dphi = (-arrayEphiP2(nRRow - 1, nZColumn - 1) + 8 * arrayEphiP(nRRow - 1, nZColumn - 1) -
            8 * arrayEphiM(nRRow - 1, nZColumn - 1) + arrayEphiM2(nRRow - 1, nZColumn - 1)) / (12 * gridSizePhi);

    arrayCharge(nRRow - 1, nZColumn - 1) = -1 * (arrayEr(nRRow - 1, nZColumn - 1) / radius + dr + dphi / radius + dz);
  }
}

///
/// \param matricesCharge
/// \param matricesEr
/// \param matricesEphi
/// \param matricesEz
/// \param matricesInvLocalIntErDz
/// \param matricesInvLocalIntEphiDz
/// \param matricesInvLocalEz
/// \param matricesDistDrDz
/// \param matricesDistDphiRDz
/// \param matricesDistDz
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
/// \param nsize
/// \param useCylAC
/// \param stepR
/// \param stepZ
/// \param stepPhi
/// \param interpType
/// \param inverseType
void AliTPCSpaceCharge3DDriftLine::InverseDistortionMaps(
        TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEphi, TMatrixD **matricesEz,
        TMatrixD **matricesInvLocalIntErDz,TMatrixD **matricesInvLocalIntEphiDz, TMatrixD **matricesInvLocalEz,
        TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz, TMatrixD **matricesDistDz,
        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice, const Int_t nsize,
        const Bool_t useCylAC, Int_t stepR, Int_t stepZ, Int_t stepPhi, Int_t interpType, Int_t inverseType
) {
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
  // memory alocation
  if (fInitLookUp) {
    // 1)  get local distortion
    if (inverseType == 2)
      InverseGlobalToLocalDistortionGlobalInvTable(matricesDistDrDz, matricesDistDphiRDz, matricesDistDz, rList,
                                                   zList, phiList, nRRow, nZColumn, phiSlice, nsize, useCylAC,
                                                   stepR, stepZ, stepPhi, interpType);
    else if (inverseType == 1)
      InverseGlobalToLocalDistortionTwoStages(matricesDistDrDz, matricesDistDphiRDz, matricesDistDz, rList, zList,
                                              phiList, nRRow, nZColumn, phiSlice, nsize, useCylAC, stepR, stepZ,
                                              stepPhi, interpType);
    else
      InverseGlobalToLocalDistortion(matricesDistDrDz, matricesDistDphiRDz, matricesDistDz, rList, zList, phiList,
                                     nRRow, nZColumn, phiSlice, nsize, useCylAC, stepR, stepZ, stepPhi, interpType);


    fLookupInverseDistA->SetLookUpR(matricesDistDrDz);
    fLookupInverseDistA->SetLookUpPhi(matricesDistDphiRDz);
    fLookupInverseDistA->SetLookUpZ(matricesDistDz);
    fLookupInverseDistA->CopyVals();

    // 2)  calculate local integral
    InverseLocalDistortionToElectricField(matricesEr, matricesEphi, matricesEz, matricesInvLocalIntErDz,
                                          matricesInvLocalIntEphiDz, matricesInvLocalEz,
                                          matricesDistDrDz, matricesDistDphiRDz, matricesDistDz, rList, zList,
                                          phiList, nRRow, nZColumn, phiSlice);
    // 3)  get potential from electric field assuming zero boundaries
    InverseElectricFieldToCharge(matricesCharge, matricesEr, matricesEphi, matricesEz, rList, zList, phiList, nRRow,
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
/// * Integrate (E(z)-Ezstd) from point of origin to pad plane
///
/// \param arrayofArrayV TMatrixD** 3D matrix representing calculated potential
/// \param arrayofEroverEz TMatrix** 3D matrix representing e-field at Er/Ez
/// \param arrayofEPhioverEz TMatrix** 3D matrix representing e-field at Ephi/Ez
/// \param arrayofDeltaZ TMatrix** 3D matrix representing e-field at DeltaZ
/// \param rows Int_t number of rows of discritization (in R direction)
/// \param columns Int_t number of columns  of discritization (in Z direction)
/// \param phislices Int_t number of (phislices in phi direction)
/// \param symmetry Int_t symmetry?
/// \param rocDisplace rocDisplacement
///
/// \pre   Matrix arrayofArrayV is assumed had been calculated  by Poisson solver
/// \post  Results of Integration and Derivations for E-field calculation are stored in arrayofEroverEz, arrayofEPhioverEz, arrayofDeltaZ
///
void AliTPCSpaceCharge3DDriftLine::CalculateEField (
                TMatrixD **arrayofArrayV, TMatrixD **arrayofEroverEz, TMatrixD **arrayofEPhioverEz,
                TMatrixD **arrayofDeltaEz, const Int_t rows, const Int_t columns, const Int_t phislices,
                const Int_t symmetry, Bool_t rocDisplacement
        ) {

  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (rows - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (columns - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phislices;


  TMatrixD *arrayOfArrayEr[phislices], *arrayOfArrayEz[phislices], *arrayOfArrayEphi[phislices];


  AliSysInfo::AddStamp("CalcField", 100, 0, 0);

  //Allocate memory for electric field r,z, phi direction
  for (Int_t k = 0; k < phislices; k++) {
    arrayOfArrayEr[k] = new TMatrixD(rows, columns);
    arrayOfArrayEz[k] = new TMatrixD(rows, columns);
    arrayOfArrayEphi[k] = new TMatrixD(rows, columns);
  }

  //Differentiate V(r) and solve for E(r) using special equations for the first and last row
  TStopwatch w;
  w.Start();


  ElectricField ( arrayofArrayV,  arrayOfArrayEr, arrayOfArrayEphi, arrayOfArrayEz,  rows,  columns,
                  phislices, gridSizeR, gridSizePhi, gridSizeZ, symmetry, fgkIFCRadius);

  w.Stop();
  AliInfo(Form("Time for calculation E-field CPU = %f s\n", w.CpuTime()));

  AliSysInfo::AddStamp("Electron Drift Calc", 120, 0, 0);

  //Integrate E(r)/E(z) from point of origin to pad plane

  IntegrateEz(arrayofEroverEz, arrayOfArrayEr, rows, columns, phislices, ezField);
  IntegrateEz(arrayofEPhioverEz, arrayOfArrayEphi, rows, columns, phislices, ezField);
  IntegrateEz(arrayofDeltaEz, arrayOfArrayEz, rows, columns, phislices, -1.0);

  // calculate z distortion from the integrated Delta Ez residuals
  // and include the aquivalence (Volt to cm) of the ROC shift !!
  for (Int_t m = 0; m < phislices; m++) {
    TMatrixD &arrayV = *arrayofArrayV[m];
    TMatrixD &deltaEz = *arrayofDeltaEz[m];

    for (Int_t j = 0; j < columns; j++) {
      for (Int_t i = 0; i < rows; i++) {
        // Scale the Ez distortions with the drift velocity pertubation -> delivers cm
        deltaEz(i, j) = deltaEz(i, j) * fgkdvdE;
        // ROC Potential in cm aquivalent
        Double_t dzROCShift = arrayV(i, columns - 1) / ezField;
        if (rocDisplacement) deltaEz(i, j) = deltaEz(i, j) + dzROCShift;  // add the ROC misaligment
      }
    }
  }
  // clear the temporary arrays lists

  for (Int_t k = 0; k < phislices; k++) {
    delete arrayOfArrayEr[k];
    delete arrayOfArrayEz[k];
    delete arrayOfArrayEphi[k];
  }
}


///
/// Integrate at z direction Ez for electron drift calculation
///
///
/// \param arrayofArrayExoverEz TMatrixD** 3D matrix representing ExoverEz
/// \param arrayofArrayEx TMatrix** 3D matrix representing e-field at x direction
/// \param rows const Int_t number of rows of discritization (in R direction)
/// \param columns const Int_t number of columns  of discritization (in Z direction)
/// \param phislices const Int_t number of (phislices in phi direction)
/// \param ezField const Double_t Electric field in z direction
///
/// \pre   arrayofArrayEx is assumed already been calculated by ElectricFieldCalculation
/// \post  Matrix arrayofArrayExoverEz is calculated by integration of arrayofArrayEx
///
void AliTPCSpaceCharge3DDriftLine::IntegrateEz
        (
                TMatrixD **arrayofArrayExoverEz,
                TMatrixD **arrayofArrayEx,
                const Int_t rows,
                const Int_t columns,
                const Int_t phislices,
                const Double_t ezField
        ) {

  const Float_t gridSizeZ = fgkTPCZ0 / (columns - 1);

  for (Int_t m = 0; m < phislices; m++) {
    TMatrixD &eXoverEz = *arrayofArrayExoverEz[m];
    TMatrixD &arrayEx = *arrayofArrayEx[m];

    for (Int_t j = columns - 1; j >= 0; j--) {
      for (Int_t i = 0; i < rows; i++) {

        /// Calculate integration from int^{0}_{j} (TODO: Split the integration)
        if (j < columns - 3) {
          eXoverEz(i, j) = eXoverEz(i, j + 2) +
                           (gridSizeZ / 3.0) * (arrayEx(i, j) + 4 * arrayEx(i, j + 1) + arrayEx(i, j + 2)) /
                           (-1 * ezField);
        } else {
          if (j == columns - 3) {
            eXoverEz(i, j) = (gridSizeZ / 3.0) * (arrayEx(i, columns - 3) + 4 * arrayEx(i, columns - 2) +
                                                  arrayEx(i, columns - 1)) / (-1 * ezField);
          }
          if (j == columns - 2) {
            eXoverEz(i, j) =
                    (gridSizeZ / 3.0) * (1.5 * arrayEx(i, columns - 2) + 1.5 * arrayEx(i, columns - 1)) /
                    (-1 * ezField);
          }
          if (j == columns - 1) {
            eXoverEz(i, j) = 0.0;
          }
        }
      }
    }
  }

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetCorrectionCylNoDrift(const Float_t x[], const Short_t roc, Float_t dx[]) {
  /// Calculates the correction due the Space Charge effect within the TPC drift volume

  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
//    InitSpaceCharge3DDistortion();
    return;
  }

  Float_t intEr, intEphi, intdEz;
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

  // Get the Er and Ephi field integrals plus the integral over DeltaEz
  if (sign == -1 && z < 0.0) {
    printf("call C side\n");
    fLookupIntENoDriftC->GetValue(r, phi, z, intEr, intEphi, intdEz);
  } else
    fLookupIntENoDriftA->GetValue(r, phi, z, intEr, intEphi, intdEz);

  // Calculate distorted position
  if (r > 0.0) {
    phi = phi + fCorrectionFactor * (fC0 * intEphi - fC1 * intEr) / r;
    r = r + fCorrectionFactor * (fC0 * intEr + fC1 * intEphi);
  }
  Double_t dz = intdEz * fCorrectionFactor * fgkdvdE;

  // Calculate correction in cartesian coordinates
  dx[0] = -(r - x[0]);
  dx[1] = -(phi - x[1]);
  dx[2] = -dz;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetDistortionCylNoDrift(const Float_t x[], Short_t roc, Float_t dx[]) {
  /// This function delivers the distortion values dx in respect to the inital coordinates x
  /// roc represents the TPC read out chamber (offline numbering convention)

  GetCorrectionCylNoDrift(x, roc, dx);
  for (Int_t j = 0; j < 3; ++j) dx[j] = -dx[j];
}

/// inverse for no drift
/// inverse from global distortion to local distortion
///
/// \param matricesDistDrDz
/// \param matricesDistDphiRDz
/// \param matricesDistDz
/// \param rList
/// \param zList
/// \param phiList
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseGlobalToLocalDistortionNoDrift
        (
                TMatrixD **matricesDistDrDz,
                TMatrixD **matricesDistDphiRDz,
                TMatrixD **matricesDistDz,
                Double_t *rList,
                Double_t *zList,
                Double_t *phiList,
                const Int_t nRRow,
                const Int_t nZColumn,
                const Int_t phiSlice
        ) {
  Double_t z, phi, r, zaft, zprev, zstep, ddr, ddphir, ddz, zl, dr, dphir, dz;
  Float_t x[3], dx[3], pdx[3], dxp1[3], dxp2[3];

  Int_t roc;
  Int_t maxzstep = 1;
  //Int_t izstep;


  TMatrixD *distDrDz;
  TMatrixD *distDphiRDz;
  TMatrixD *distDz;


  for (Int_t k = 0; k < phiSlice; k++) {
    distDrDz = matricesDistDrDz[k];
    distDphiRDz = matricesDistDphiRDz[k];
    distDz = matricesDistDz[k];

    for (Int_t i = 0; i < nRRow; i++) {
      (*distDrDz)(i, nZColumn - 1) = 0.0;
      (*distDphiRDz)(i, nZColumn - 1) = 0.0;
      (*distDz)(i, nZColumn - 1) = 0.0;

    }
  }


  zstep = (zList[1] - zList[0]) / maxzstep;

  for (Int_t j = nZColumn - 2; j >= 0; j--) {

    roc = 0; // FIXME

    for (Int_t k = 0; k < phiSlice; k++) {

      distDrDz = matricesDistDrDz[k];
      distDphiRDz = matricesDistDphiRDz[k];
      distDz = matricesDistDz[k];
      for (Int_t i = 0; i < nRRow; i++) {
        // get global distortion

        r = rList[i];
        phi = phiList[k];
        z = zList[j];
        zprev = zList[j + 1];
        //zaft = zList[j-1];
        zl = fListZA[j];  // Symmetric solution in Z that depends only on ABS(Z)

        //printf("(%f,%f)\n",z,zl);



        (*distDrDz)(i, j) = 0.0;
        (*distDphiRDz)(i, j) = 0.0;
        (*distDz)(i, j) = 0.0;
        dr = 0.0;
        dphir = 0.0;
        dz = 0.0;

        r = rList[i];
        phi = phiList[k];
        z = zList[j];

        x[0] = r;
        x[1] = phi;
        x[2] = z;

        GetDistortionCylNoDrift(x, roc, dx);


        //x[0] = x[0] + dr;
        //x[1] = x[1] + dphir/r;
        x[2] = zprev;


        GetDistortionCylNoDrift(x, roc, pdx);


        (*distDrDz)(i, j) = (dx[0] - pdx[0]);
        (*distDphiRDz)(i, j) = (dx[1] - pdx[1]) * r;
        (*distDz)(i, j) = (dx[2] - pdx[2]);

      }
    }
  }
}

///
/// \param matricesCharge
/// \param matricesEr
/// \param matricesEphi
/// \param matricesEz
/// \param matricesInvLocalIntErDz
/// \param matricesInvLocalIntEphiDz
/// \param matricesInvLocalEz
/// \param matricesDistDrDz
/// \param matricesDistDphiRDz
/// \param matricesDistDz
/// \param nRRow
/// \param nZColumn
/// \param phiSlice
void AliTPCSpaceCharge3DDriftLine::InverseDistortionMapsNoDrift(
        TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEphi, TMatrixD **matricesEz,
        TMatrixD **matricesInvLocalIntErDz, TMatrixD **matricesInvLocalIntEphiDz, TMatrixD **matricesInvLocalEz,
        TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDphiRDz, TMatrixD **matricesDistDz,
        const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice
) {
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
  // memory alocation
  if (fInitLookUp) {
    // 1)  get local distortion
    InverseGlobalToLocalDistortionNoDrift(matricesDistDrDz, matricesDistDphiRDz, matricesDistDz, rList, zList,
                                          phiList, nRRow, nZColumn, phiSlice);
    // 2)  calculate local integral
    InverseLocalDistortionToElectricField(matricesEr, matricesEphi, matricesEz, matricesInvLocalIntErDz,
                                          matricesInvLocalIntEphiDz, matricesInvLocalEz,
                                          matricesDistDrDz, matricesDistDphiRDz, matricesDistDz, rList, zList,
                                          phiList, nRRow, nZColumn, phiSlice);
    // 3)  get potential from electric field assuming zero boundaries
    InverseElectricFieldToCharge(matricesCharge, matricesEr, matricesEphi, matricesEz, rList, zList, phiList, nRRow,
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
void AliTPCSpaceCharge3DDriftLine::GetChargeDensity
        (
                TMatrixD **matricesChargeA, TMatrixD **matricesChargeC, TH3 *spaceChargeHistogram3D,
                const Int_t nRRow, const Int_t nZColumn, const Int_t phiSlice
        ) {
  Int_t phiSlicesPerSector = phiSlice / kNumSector;
  const Float_t gridSizeR = (fgkOFCRadius - fgkIFCRadius) / (nRRow - 1);
  const Float_t gridSizeZ = fgkTPCZ0 / (nZColumn - 1);
  const Float_t gridSizePhi = TMath::TwoPi() / phiSlice;
  const Double_t ezField = (fgkCathodeV - fgkGG) / fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;
  // local variables
  Float_t radius0, phi0, z0;
  // list of point as used in the poisson relaxation and the interpolation (for interpolation)
  Double_t rlist[nRRow], zedlist[nZColumn], philist[phiSlice];
  for (Int_t k = 0; k < phiSlice; k++) philist[k] = gridSizePhi * k;
  for (Int_t i = 0; i < nRRow; i++) rlist[i] = fgkIFCRadius + i * gridSizeR;
  for (Int_t j = 0; j < nZColumn; j++) zedlist[j] = j * gridSizeZ;

  TMatrixD *mCharge;
  for (Int_t side = 0; side < 2; side++) {
    for (Int_t k = 0; k < phiSlice; k++) {
      if (side == 0)
        mCharge = matricesChargeA[k];
      else
        mCharge = matricesChargeC[k];

      phi0 = philist[k];
      for (Int_t i = 0; i < nRRow; i++) {
        radius0 = rlist[i];
        for (Int_t j = 0; j < nZColumn; j++) {
          z0 = zedlist[j];
          if (side == 1) z0 = -TMath::Abs(zedlist[j]);
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
Double_t AliTPCSpaceCharge3DDriftLine::GetChargeCylAC ( const Float_t x[], Short_t roc) {
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
Double_t AliTPCSpaceCharge3DDriftLine::GetInverseChargeCylAC ( const Float_t x[], Short_t roc ) {
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
  Float_t dR, dPhiR, dZ;
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
    fLookupDistA->GetValue(r, phi, z, dR, dPhiR, dZ);
  else
    fLookupDistC->GetValue(r, phi, -z, dR, dPhiR, dZ);


  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor *
          dZ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}

///
/// \param x
/// \param roc
/// \param dx
void AliTPCSpaceCharge3DDriftLine::GetInverseLocalDistortionCylAC(const Float_t x[], Short_t roc, Float_t dx[]) {
  if (!fInitLookUp) {
    AliInfo("Lookup table was not initialized! Performing the inizialisation now ...");
    InitSpaceCharge3DPoissonIntegralDz(129, 129, 144, 100, 1e-8);
  }


  Float_t dR, dPhiR, dZ;
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
    fLookupInverseDistA->GetValue(r, phi, z, dR, dPhiR, dZ);
  else
    fLookupInverseDistC->GetValue(r, phi, -z, dR, dPhiR, dZ);

  dx[0] = fCorrectionFactor * dR;
  dx[1] = fCorrectionFactor * dPhiR;
  dx[2] = fCorrectionFactor * dZ;  // z distortion - (scaled with driftvelocity dependency on the Ez field and the overall scaling factor)

}

/// Function for setting Potential Boundary Values and Charge distribution input TFormula
///
/// \param vTestFunction
/// \param rhoTestFunction
///
void AliTPCSpaceCharge3DDriftLine::SetPotentialBoundaryAndCharge
        (
                TFormula *vTestFunction,
                TFormula *rhoTestFunction
        ) {
  /**** allocate memory for charge ***/
  // we allocate pointer to TMatrixD array to picture 3D (slices), this representation should be revised
  // since it is easier for GPU implementation to run for 1D memory

  fFormulaPotentialV = new TFormula(*vTestFunction);
  fFormulaChargeRho = new TFormula(*rhoTestFunction);

  // grid size for one side
  TMatrixD *chargeA;
  TMatrixD *chargeC;

  Double_t radius0, z0, phi0, z0neg;

  Int_t indexb = 0;
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
          fListPotentialBoundaryA[indexb] = vTestFunction->Eval(radius0, phi0, z0);
          fListPotentialBoundaryC[indexb] = vTestFunction->Eval(radius0, phi0, z0neg);
          indexb++;
        }

      } // end j
    } // end i
  } // end phi

  fInterpolatorChargeA->SetValue(fMatrixChargeA);
  fInterpolatorChargeA->InitCubicSpline();
  fInterpolatorChargeC->SetValue(fMatrixChargeC);
  fInterpolatorChargeC->InitCubicSpline();
}

