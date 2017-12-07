/*
  Unit test for the  AliTPCSpaceCharge3DDriftLine class:

  Tests:
  1.) Check consistency for distortion and correction

  usage:

  .x AliTPCSpaceCharge3DDriftLineTest.C+

*/

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TLegend.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"
#include "TH3.h"
#include "TH2.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TMap.h"
#include "TFile.h"
#include "TGraph.h"
#include "AliTPCSpaceCharge3DDriftLine.h"
#include "AliTPCPoissonSolver.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include "TGrid.h"
#include "TStatToolkit.h"

#endif

void UnitTestConsistencyRDistortion();
void AliTPCSpaceCharge3DDriftLineTest() {
  UnitTestConsistencyRDistortion();
}

/// unit test for R distortion
void UnitTestConsistencyRDistortion() {
  //
  // Unit Test Check the Consistency of correction --
  // generate nPoints of point, calculate its distortion
  // calculate correction for distorted points (R direction)
  // return OK if  distance(corrected distorted point, original point)  < epsilon
  // constructor with parameters: number of grids

  const Int_t phiSlice = 18;
  const Int_t rRow = 17;
  const Int_t zColumn = 17;
  const Int_t maxIter = 200;
  const Float_t convergenceError = 1e-8;

  const Int_t nPhiSliceTest = 36;
  const Int_t nRRowTest = 33;
  const Int_t nZColumnTest = 33;
  const Int_t nPoints = nPhiSliceTest * nRRowTest * nZColumnTest;
  const Double_t maxEpsilon = 1e-2;

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Begin");
  AliTPCSpaceCharge3DDriftLine *sc = new AliTPCSpaceCharge3DDriftLine("unitTestDrDist", "unitTestDrDist", rRow, zColumn,
                                                                      phiSlice);

  TFormula vTestFunction("v", " [0]*(-x^4 + 338.0 *x^3 - 21250.75 * x^2)*cos([1]* y)^2*exp(-1* [2] * z^2)");
  TFormula rhoTestFunction("rho",
                           "[0]*(((-16.0 * x^2 + 9.0 * 338.0 * x - 4.0*21250.75) *cos([1] * y)^2 * exp(-1 *[2]*z^2)) - ((-1* x^2 + 338.0 * x - 21250.75) * 2 * [1]^2 * cos(2 * [1] * y) * exp(-1 *[2]*z^2)) + ((-1* x^4 +  338.0 * x^3 - 21250.75 * x^2) * cos([1] * y)^2 * (4*[2]^2*z^2 - 2 * [2]) * exp(-1 *[2]*z^2)))");

  // set parameters for function
  Double_t a = AliTPCPoissonSolver::fgkOFCRadius * AliTPCPoissonSolver::fgkOFCRadius;
  a *= (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius);
  a *= (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius);
  a = -1 * (1000.0 / a);
  Double_t b = 0.5;
  Double_t c = 1.0 / (((AliTPCPoissonSolver::fgkTPCZ0) / 2.0) * ((AliTPCPoissonSolver::fgkTPCZ0) / 2.0));

  vTestFunction.SetParameters(a, b, c);
  rhoTestFunction.SetParameters(a, b, c);


  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Setting up Poisson Solver and Boundary Charge");
  // create Poisson Solver
  AliTPCPoissonSolver *poissonSolver = new AliTPCPoissonSolver();

  sc->SetPoissonSolver(poissonSolver);

  //	Set Potential Boundary And Charge function
  sc->SetPotentialBoundaryAndCharge(&vTestFunction, &rhoTestFunction);

  // space charge parameters
  Double_t c0 = 0.904466;
  Double_t c1 = 0.904466;// set constant parameters
  sc->SetC0C1(c0, c1);
  sc->SetCorrectionType(1);

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Create Distortion/Correction maps");
  sc->InitSpaceCharge3DPoissonIntegralDz(rRow, zColumn, phiSlice, maxIter, convergenceError);

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Creating test points");

  Double_t residueMean = 0.0;
  Double_t residueRMS = 0.0;
  Double_t residueMax = 0.0;

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion",
         TString::Format("Residue mean=%f,rms=%f, max=%f", residueMean, residueRMS, residueMax).Data());

  if (residueMean < maxEpsilon)
    ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Test OK: Mean Residue=%.2E < %.2E",
           residueMean, maxEpsilon);
  else
    ::Error("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Test FAILED: Mean Residue=%.2E > %.2E",
            residueMean, maxEpsilon);
}
