/*
  Unit test for the  AliTPCSpaceCharge3DDriftLine class:

  Tests:
  1.) Check consistency for distortion and correction

  usage:

  .L loadClasses.C
  loadClasses
  .L AliTPCSpaceCharge3DDriftLineTest.C
  UnitTestConsistency()

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

void UnitTestConsistencyPotentialBoundaryInZ();

void CreateTestDump(AliTPCSpaceCharge3DDriftLine *sc, const Int_t nPoint);


/// unit test for R distortion
void UnitTestConsistencyRDistortion() {
  //
  // Unit Test Check the Consistency of correction --
  // generate npoints of point, calculate its distortion
  // calculate correction for distorted points (R direction)
  // return OK if  distance(corrected distorted point, original point)  < epsilon
  // contructor with parameters: number of grids

  const Int_t phiSlice = 36;
  const Int_t rRow = 33;
  const Int_t zColumn = 33;
  const Int_t maxIter = 200;
  const Float_t convError = 1e-8;

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

  sc->InitSpaceCharge3DPoissonIntegralDz(rRow, zColumn, phiSlice, maxIter, convError);

  // create tree fber or many random points
  //
  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Creating test points");
  CreateTestDump(sc, nRRowTest, nZColumnTest, nPhiSliceTest);
  // calculate error

  TFile f("TpcDumpDistortionDriftInverseAnalytic.root");
  TTree *t = f.Get("testDump");
  // get the biggest distortion r
  t->Draw("abs(dr)", "", "goff");
  Double_t maxVar = TMath::MaxElement(nPoints, t->GetV1());
  t->Draw(TString::Format("abs(dr+cr)/%f", maxVar).Data(), "", "goff");

  Double_t residueMean = TMath::Mean(nPoints, t->GetV1());
  Double_t residueRMS = TMath::RMS(nPoints, t->GetV1());
  Double_t residueMax = TMath::MaxElement(nPoints, t->GetV1());

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion",
         TString::Format("Residue mean=%f,rms=%f, max=%f", residueMean, residueRMS, residueMax).Data());

  // result
  if (residueMean < maxEpsilon)
    ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Mean Residue=%.2E < %.2E - OK",
           residueMean, maxEpsilon);
  else
    ::Error("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyRDistortion", "Mean Residue=%.2E > %.2E - FAILED",
            residueMean, maxEpsilon);
}


/// create a tree for storing testing num
void CreateTestDump(AliTPCSpaceCharge3DDriftLine *sc, const Int_t nRRowTest, const Int_t nZColTest,
                    const Int_t nPhiSliceTest) {
  //	create a tree file
  // 	create a tree

  TFile f("TpcDumpDistortionDriftInverseAnalytic.root", "RECREATE");
  TTree *testDump = new TTree("testDump", "Collections of AliTPCSpaceCharge3DDriftLine");


  Int_t nrrow;
  Int_t nzcol;
  Int_t nphislice;
  Int_t orderinterp;
  Int_t rbfkernel;
  Int_t ireggridsize;


  const Int_t nPoint = nRRowTest * nZColTest * nPhiSliceTest;

//	printf("nPoint = %d\n",nPoint);

  Float_t ofcRadius = AliTPCPoissonSolver::fgkOFCRadius;
  Float_t ifcRadius = AliTPCPoissonSolver::fgkIFCRadius;
  Float_t tpcZ0 = AliTPCPoissonSolver::fgkTPCZ0;

  Float_t dRadius = (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius) / (nRRowTest - 1);
  Float_t dZ = AliTPCPoissonSolver::fgkTPCZ0 / (nZColTest - 1);
  Float_t dphi = TMath::TwoPi() / nPhiSliceTest;


  testDump->Branch("nPoint", &nPoint, "nPoint/I");
  testDump->Branch("nrrow", &nrrow, "nrrow/I");
  testDump->Branch("nzcol", &nzcol, "nzcol/I");
  testDump->Branch("nphislice", &nphislice, "nphislice/I");
  testDump->Branch("rbfkernel", &rbfkernel, "rbfkernel/I");
  testDump->Branch("orderinterp", &orderinterp, "orderinterp/I");
  testDump->Branch("ireggridsize", &ireggridsize, "ireggridsize/I");


  Double_t phi0[nPoint];
  Double_t z0[nPoint];
  Double_t r0[nPoint];

  Float_t dr[nPoint];
  Float_t dphir[nPoint];
  Float_t dz[nPoint];

  Float_t cr[nPoint];
  Float_t cphir[nPoint];
  Float_t cz[nPoint];


  Double_t charge[nPoint];


  testDump->Branch("phi0", phi0, "phi0[nPoint]/D");
  testDump->Branch("z0", z0, "z0[nPoint]/D");
  testDump->Branch("r0", r0, "r0[nPoint]/D");

  testDump->Branch("dr", dr, "dr[nPoint]/F");
  testDump->Branch("dphir", dphir, "dphir[nPoint]/F");
  testDump->Branch("dz", dz, "dz[nPoint]/F");

  testDump->Branch("cr", cr, "cr[nPoint]/F");
  testDump->Branch("cphir", cphir, "cphir[nPoint]/F");
  testDump->Branch("cz", cz, "cz[nPoint]/F");

  testDump->Branch("charge", charge, "charge[nPoint]/D");

  // store sc configuration
  nrrow = sc->GetNRRows();
  nzcol = sc->GetNZColumns();
  nphislice = sc->GetNPhiSlices();
  orderinterp = sc->GetInterpolationOrder();
  ireggridsize = sc->GetIrregularGridSize();
  rbfkernel = sc->GetRBFKernelType();

  Int_t rocnum;

  Float_t point0[] = {0.0, 0.0, 0.0};
  Float_t point1[] = {0.0, 0.0, 0.0};
  Float_t dist[] = {0.0, 0.0, 0.0};
  Float_t corr[] = {0.0, 0.0, 0.0};

  Float_t r1;
  Float_t phi1;
  Float_t z1;


  Int_t index = 0;

  for (Int_t m = 0; m < nPhiSliceTest; m++) {
    printf("slice: %d\n", m);
    for (Int_t i = 0; i < nRRowTest; i++) {
      for (Int_t j = 0; j < nZColTest; j++) {
        index = m * nRRowTest * nZColTest + i * nZColTest + j;
//				index = m * (nRRowTest - 2) * (nZColTest - 2) + (i - 1) * (nZColTest - 2) + (j - 1);
        printf("index = %d\n", index);

        phi0[index] = m * dphi;
        r0[index] = ifcRadius + (dRadius * i);
        z0[index] = dZ * j;

        //	printf("(%f,%f,%f)\t",phi0[index],r0[index],z0[index]);


        // get charge density at this point (interpolation)
        charge[index] = sc->GetSpaceChargeDensity(r0[index], phi0[index], z0[index]);

        point0[0] = r0[index];
        point0[1] = phi0[index];
        point0[2] = z0[index];

        rocnum = (z0 > 0.0) ? 0 : 18;

        // lookup distortion on the selected point
        sc->GetDistortionCyl(point0, rocnum, dist);

        dr[index] = dist[0];
        dphir[index] = dist[1];
        dz[index] = dist[2];

        //	printf("(%f,%f,%f)\n",dist[0],dist[1],dist[2]);


        // create distorted point
        r1 = r0[index] + dist[0];
        phi1 = phi0[index] + dist[1] / r0[index];
        z1 = z0[index] + dist[2];

        point1[0] = r1;
        point1[1] = phi1;
        point1[2] = z1;

        // lookup distortion on the selected point
        sc->GetCorrectionCyl(point1, rocnum, corr);

        cr[index] = corr[0];
        cphir[index] = corr[1];
        cz[index] = corr[2];
      }
    }
  }

  //	AddMetaData(testDump);
  testDump->Fill();
  f.Write();
  f.Close();
}


/// unit test consistency
/// knowing charge densities and potential boundary
void UnitTestConsistencyPotentialBoundaryInZ() {

  const Int_t phiSlice = 18;
  const Int_t rRow = 17;
  const Int_t zColumn = 17;
  const Int_t maxIter = 200;
  const Float_t convError = 1e-8;
  const Double_t maxEpsilon = 1e-2;

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyPotentialBoundaryInZ", "Begin");
  AliTPCSpaceCharge3DDriftLine *sc = new AliTPCSpaceCharge3DDriftLine(rRow, zColumn, phiSlice);


  // set potential boundary in V
  TF1 *potentialBoundaryFunctionInZ = new TF1("fa1", "TMath::DiLog(x)", -250, 250);
  // charge in histogram
  TH3 *chargeA = new TH3F("charge A", "charge A", 250, 0.0, TMath::TwoPi(), 250, 85.0, 250.0, 250, 0.0, 250.0);
  TH3 *chargeC = new TH3F("charge C", "charge C", 250, 0.0, TMath::TwoPi(), 250, 85.0, 250.0, 250, 0.0, 250.0);


  // set the input boundary potential and charge densities
  sc->SetInputSpaceCharge(chargeA, 0);
  sc->SetInputSpaceCharge(chargeC, 1);
  sc->SetPotentialBoundaryInZ(potentialBoundaryFunctionInZ, 0);
  sc->SetPotentialBoundaryInZ(potentialBoundaryFunctionInZ, 1);

  // create Poisson Solver
  AliTPCPoissonSolver *poissonSolver = new AliTPCPoissonSolver();
  sc->SetPoissonSolver(poissonSolver);

  // space charge parameters
  Double_t c0 = 0.904466;
  Double_t c1 = 0.904466;// set constant parameters
  sc->SetC0C1(c0, c1);


  sc->InitSpaceCharge3DPoissonIntegralDz(rRow, zColumn, phiSlice, maxIter, convError);

  // create tree fber or many random points
  //
  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyPotentialBoundaryInZ", "Creating distortion tree");
  sc->CreateDistortionTree(1.0);


  // calculate error

  // draw the input space charge
  //TCanvas * canvas1 = new TCanvas();
  //sc->CreateHistoSCinXY(5,200,200)->Draw("colz");
  //TCanvas * canvas2 = new TCanvas();
  //sc->CreateHistoSCinZR(1,200,200)->Draw("colz");

  // draw distortion
  //TCanvas * canvas3 = new TCanvas("canvas distortion");
  //sc->CreateHistoDistDRinZR(0,200,200)->Draw("colz");
}

