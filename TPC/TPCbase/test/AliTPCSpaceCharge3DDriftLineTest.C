/*
  Unit test for the  AliTPCSpaceCharge3DDriftLine class:

  Tests:
  1.) Check consistency for distortion and correction

  usage:

  .x AliTPCSpaceCharge3DDriftLineTest.C+

*/

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TROOT.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "AliTPCSpaceCharge3DDriftLine.h"
#include "TTreeStream.h"

#endif


void UnitTestCorrectnessDistortion(const Int_t rRow,const Int_t zColumn,const  Int_t phiSlice, Int_t rRowTest, Int_t zColumnTest, const Int_t phiSliceTest, Int_t correctionType,TTreeSRedirector *pcStream);
void UnitTestConsistencyDistortionZShort(const Int_t rRow, const Int_t zColumn, const Int_t phiSlice, const Int_t rRowTest, const Int_t zColumnTest, const Int_t phiSliceTest, const Int_t correctionType, TTreeSRedirector *pcStream);
void GetResidueFromDistortionTree(Int_t unitTestId, const char * filename,const char *numericName, const char *analyticName, Double_t *errorList, Int_t &indexErrorList,TTreeSRedirector *pcStream);
void GetResidueFromDistortionAnalyticTree(Int_t unitTestId, const char * numericFileName, const char * analyticFileName, const char *varName, Double_t *errorList, Int_t &indexErrorList,TTreeSRedirector *pcStream);
void WriteErrorToPCStream(Int_t unitTestId, Int_t rRow, Int_t zColumn, Int_t phiSlice, Int_t rRowTest, Int_t zColumnTest, Int_t phiSliceTest, Int_t correctionType, Int_t varNameId, Double_t *errorList, Int_t  indexErrorList,TTreeSRedirector *pcStream);
void PrintErrorStatus();

Double_t dFunctionVZ(Double_t *x, Double_t *par); // function test for boundary
const Double_t gkVCE  = -50000;
const Double_t gkVROC = -70;
const Int_t gkMaxIter = 200;
const Float_t gkConvError = 1e-8;
const Double_t gkMaxEpsilon = 1e-2;


// TODO: Move the helper classes to AliTPCSpaceCharge3DDriftLine
// helper class
void
LocalDistCorrDzExact(TFormula *intDrDzF, TFormula *intDPhiDzF, TFormula *intDzDzF, Double_t *rList, Double_t *phiList,
                     Double_t *zList, TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz,
                     TMatrixD **matricesDistDz, TMatrixD **matricesCorrDrDz, TMatrixD **matricesCorrDPhiRDz,
                     TMatrixD **matricesCorrDz, const Int_t rRow, const Int_t zColumn, const Int_t phiSlice,
                     const Double_t c0, const Double_t c1, const Double_t ezField, const Double_t dvdE,
                     const Int_t side);

void InitPotentialAndCharge3D(TFormula *vTestFunction, TFormula *rhoTestFunction, TFormula *erTestFunction,
                         TFormula *ePhiTestFunction, TFormula *ezTestFunction, TFormula *intErDzTestFunction,
                         TFormula *intEPhiRDzTestFunction, TFormula *intDzTestFunction, TMatrixD **matricesV,
                         TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEPhi,
                         TMatrixD **matricesEz, TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz,
                         TMatrixD **matricesDistDz, TMatrixD **matricesCorrDrDz, TMatrixD **matricesCorrDPhiRDz,
                         TMatrixD **matricesCorrDz, const Int_t rRow, const Int_t zColumn, const Int_t phiSlice,
                         const Int_t side, const Double_t c0, const Double_t c1, const Double_t ezField,
                         const Double_t dvdE);


/// main macro function, call unit test
void AliTPCSpaceCharge3DDriftLineTest() {

  Int_t rRowList[] = {17,33,65,129,257};
  Int_t zColumnList[] = {17,33,65,129,257};
  Int_t phiSliceList[] = {18,36,72,144,288};
  Int_t rRow, zColumn, phiSlice;
  const Int_t numberExperiment = 1;
  Int_t rRowTest = 50;
  Int_t zColumnTest = 50;
  Int_t phiSliceTest = 50;
  Int_t nPoint =rRowTest * zColumnTest * phiSliceTest;

  TTreeSRedirector *pcStream = new TTreeSRedirector("spaceChargeDriftLinePerformance.root","RECREATE");


  for (Int_t iExperiment = 0;iExperiment < numberExperiment;iExperiment++) {
    rRow = rRowList[iExperiment];
    zColumn = zColumnList[iExperiment];
    phiSlice = phiSliceList[iExperiment];
    UnitTestCorrectnessDistortion(rRow, zColumn, phiSlice, rRowTest, zColumnTest, phiSliceTest, 0,  pcStream); // regular
    UnitTestCorrectnessDistortion(rRow, zColumn, phiSlice, rRowTest, zColumnTest, phiSliceTest, 1,  pcStream); // irregular
  }

  delete pcStream;

  PrintErrorStatus();
}



/// unit test for Correctness
/// - a set of known functions and stored it in scExact
/// - generate numerical calculation in scNumeric
/// - calculate relative error (accepted if less than 10^-2)
void UnitTestCorrectnessDistortion(const Int_t rRow,const  Int_t zColumn, const Int_t phiSlice,   Int_t rRowTest, Int_t zColumnTest, Int_t phiSliceTest, Int_t correctionType, TTreeSRedirector *pcStream) {

  const Int_t maxIter = 100;
  const Float_t convergenceError = 1e-8;

  const Int_t nPoints = phiSliceTest * rRowTest * zColumnTest;
  const Double_t maxEpsilon = 1e-2;

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestCorrectnessDistortion", "Begin");
  AliTPCSpaceCharge3DDriftLine *sc = new AliTPCSpaceCharge3DDriftLine(TString::Format("unitTestNumeric%d-%d-%d-%d",correctionType,rRow,zColumn,phiSlice).Data(), "unitTestNumeric", rRow, zColumn, phiSlice);
  AliTPCSpaceCharge3DDriftLine *scExact = new AliTPCSpaceCharge3DDriftLine(TString::Format("unitTestExact%d-%d-%d-%d",correctionType,rRow,zColumn,phiSlice).Data(), "unitTestExact", rRow, zColumn, phiSlice);

  
  TFormula vTestFunction1("f1", "[0]*(x^4 - 338.0 *x^3 + 21250.75 * x^2)*cos([1]* y)^2*exp(-1* [2] * z^2)");
  TFormula rhoTestFunction1("ff1", "[0]*(((16.0 * x^2 - 9.0 * 338.0 * x + 4.0*21250.75) *cos([1] * y)^2 * exp(-1 *[2]*z^2)) - ((x^2 -  338.0 * x + 21250.75) * 2 * [1]^2 * cos(2 * [1] * y) * exp(-1 *[2]*z^2)) + ((x^4 -  338.0 * x^3 + 21250.75 * x^2) * cos([1] * y)^2 * (4*[2]^2*z^2 - 2 * [2]) * exp(-1 *[2]*z^2)))");

  TFormula erTestFunction1("er", " [0]*(4*x^3 - 3 * 338.0 *x^2 + 2 * 21250.75 * x)*cos([1]* y)^2*exp(-1* [2] * z^2)");
  TFormula ePhiTestFunction1("ePhi",
                            "  [0]*(x^3 - 338.0 *x^2 +  21250.75 * x)* -1  * [1] * sin(2 * [1]* y)*exp(-1* [2] * z^2)");
  TFormula ezTestFunction1("ez",
                          " [0]*(x^4 - 338.0 *x^3 + 21250.75 * x^2)*cos([1]* y)^2*-1*2*[2]*z*exp(-1* [2] * z^2)");

  TFormula intErDzTestFunction1("intErDz",
                               " [0]*(4*x^3 - 3 * 338.0 *x^2 + 2 * 21250.75 * x)*cos([1]* y)^2*((sqrt(pi)*TMath::Erf(sqrt([2]) * z))/(2 * sqrt([2]))) ");
  TFormula intEPhiRDzTestFunction1("intEPhiDz",
                                  "[0]* (x^3 - 338.0 *x^2 +  21250.75 * x)* -1  * [1] * sin(2 * [1]* y)*((sqrt(pi)*TMath::Erf(sqrt([2]) * z))/(2 * sqrt([2])))");
  TFormula intDzTestFunction1("intEzDz",
                             "[0]* (x^4 - 338.0 *x^3 + 21250.75 * x^2)*cos([1]* y)^2*exp(-1* [2] * z^2)");


  // set parameters for function
  Double_t a = AliTPCPoissonSolver::fgkOFCRadius * AliTPCPoissonSolver::fgkOFCRadius;
  a *= (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius);
  a *= (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius);
  a =  (1000.0 / a);
  a = 1e-5;
  Double_t b = 0.5;
  Double_t c = 1.0 / (((AliTPCPoissonSolver::fgkTPCZ0) / 2.0) * ((AliTPCPoissonSolver::fgkTPCZ0) / 2.0));
  c = 1e-4;

  vTestFunction1.SetParameters(a, b, c);
  rhoTestFunction1.SetParameters(a, b, c);

  erTestFunction1.SetParameters(-a, b, c);
  ePhiTestFunction1.SetParameters(-a, b, c);
  ezTestFunction1.SetParameters(-a, b, c);
  intErDzTestFunction1.SetParameters(-a, b, c);
  intEPhiRDzTestFunction1.SetParameters(-a, b, c);
  intDzTestFunction1.SetParameters(-a, b, c);

  TFormula *vTestFunction = &vTestFunction1;
  TFormula *rhoTestFunction = &rhoTestFunction1;


  TFormula *erTestFunction = &erTestFunction1;
  TFormula *ePhiTestFunction = &ePhiTestFunction1;
  TFormula *ezTestFunction = &ezTestFunction1;

  TFormula *intErDzTestFunction = &intErDzTestFunction1;
  TFormula *intEPhiRDzTestFunction = &intEPhiRDzTestFunction1;
  TFormula *intDzTestFunction = &intDzTestFunction1;



  Float_t ofcRadius = AliTPCPoissonSolver::fgkOFCRadius;
  Float_t ifcRadius = AliTPCPoissonSolver::fgkIFCRadius;
  Float_t tpcZ0 = AliTPCPoissonSolver::fgkTPCZ0;

  Double_t ezField = (AliTPCPoissonSolver::fgkCathodeV - AliTPCPoissonSolver::fgkGG) /
                     AliTPCPoissonSolver::fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm;
  Double_t dvdE = AliTPCPoissonSolver::fgkdvdE;

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestCorrectnessDistortion","Setting up Poisson Solver and Boundary Charge");
  // create Poisson Solver
  if (correctionType == 0)
    ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestCorrectnessDistortion","Case regular grid interpolation for Correction");
  else
    ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestCorrectnessDistortion","Case irregular grid interpolation for Correction");

  sc->SetPotentialBoundaryAndChargeFormula(vTestFunction, rhoTestFunction);
  sc->SetElectricFieldFormula(erTestFunction,ePhiTestFunction,ezTestFunction);
  sc->SetCorrectionType(correctionType);
  sc->SetOmegaTauT1T2(-0.35, 1., 1.);

  if (correctionType == 0) {
	  scExact->SetPotentialBoundaryAndChargeFormula(vTestFunction, rhoTestFunction);
	  scExact->SetElectricFieldFormula(erTestFunction,ePhiTestFunction,ezTestFunction);
	  scExact->SetCorrectionType(correctionType);
	  scExact->SetOmegaTauT1T2(-0.35, 1., 1.);
  
  	 // 	Float_t c0 = scExact->GetC0();
	 //      Float_t c1 = scExact->GetC1();
  }
  Float_t c0 = sc->GetC0();
  Float_t c1 = sc->GetC1();

  // for numeric
  sc->InitSpaceCharge3DPoissonIntegralDz(rRow, zColumn, phiSlice, maxIter, convergenceError);


  TMatrixD *matricesDistDrDzExactA[phiSlice];
  TMatrixD *matricesDistDrDzExactC[phiSlice];
  TMatrixD *matricesDistDPhiRDzExactA[phiSlice];
  TMatrixD *matricesDistDPhiRDzExactC[phiSlice];
  TMatrixD *matricesDistDzExactA[phiSlice];
  TMatrixD *matricesDistDzExactC[phiSlice];
  TMatrixD *matricesCorrDrDzExactA[phiSlice];
  TMatrixD *matricesCorrDrDzExactC[phiSlice];
  TMatrixD *matricesCorrDPhiRDzExactA[phiSlice];
  TMatrixD *matricesCorrDPhiRDzExactC[phiSlice];
  TMatrixD *matricesCorrDzExactA[phiSlice];
  TMatrixD *matricesCorrDzExactC[phiSlice];

  TMatrixD *matricesChargeA[phiSlice];
  TMatrixD *matricesChargeC[phiSlice];
  TMatrixD *matricesVExactA[phiSlice];
  TMatrixD *matricesVExactC[phiSlice];
  TMatrixD *matricesErExactA[phiSlice];
  TMatrixD *matricesErExactC[phiSlice];
  TMatrixD *matricesEPhiExactA[phiSlice];
  TMatrixD *matricesEPhiExactC[phiSlice];
  TMatrixD *matricesEzExactA[phiSlice];
  TMatrixD *matricesEzExactC[phiSlice];

  for (Int_t m = 0; m < phiSlice; m++) {
    matricesDistDrDzExactA[m] = new TMatrixD(rRow, zColumn);
    matricesDistDPhiRDzExactA[m] = new TMatrixD(rRow, zColumn);
    matricesDistDzExactA[m] = new TMatrixD(rRow, zColumn);
    matricesCorrDrDzExactA[m] = new TMatrixD(rRow, zColumn);
    matricesCorrDPhiRDzExactA[m] = new TMatrixD(rRow, zColumn);
    matricesCorrDzExactA[m] = new TMatrixD(rRow, zColumn);

    matricesDistDrDzExactC[m] = new TMatrixD(rRow, zColumn);
    matricesDistDPhiRDzExactC[m] = new TMatrixD(rRow, zColumn);
    matricesDistDzExactC[m] = new TMatrixD(rRow, zColumn);
    matricesCorrDrDzExactC[m] = new TMatrixD(rRow, zColumn);
    matricesCorrDPhiRDzExactC[m] = new TMatrixD(rRow, zColumn);
    matricesCorrDzExactC[m] = new TMatrixD(rRow, zColumn);

    matricesVExactA[m] = new TMatrixD(rRow, zColumn);
    matricesChargeA[m] = new TMatrixD(rRow, zColumn);
    matricesErExactA[m] = new TMatrixD(rRow, zColumn);
    matricesEPhiExactA[m] = new TMatrixD(rRow, zColumn);
    matricesEzExactA[m] = new TMatrixD(rRow, zColumn);

    matricesVExactC[m] = new TMatrixD(rRow, zColumn);
    matricesChargeC[m] = new TMatrixD(rRow, zColumn);
    matricesErExactC[m] = new TMatrixD(rRow, zColumn);
    matricesEPhiExactC[m] = new TMatrixD(rRow, zColumn);
    matricesEzExactC[m] = new TMatrixD(rRow, zColumn);
  }

  if (correctionType == 0) {
	  InitPotentialAndCharge3D(vTestFunction, rhoTestFunction, erTestFunction, ePhiTestFunction, ezTestFunction,
                      intErDzTestFunction, intEPhiRDzTestFunction, intDzTestFunction, matricesVExactA,
                      matricesChargeA, matricesErExactA, matricesEPhiExactA, matricesEzExactA, matricesDistDrDzExactA,
                      matricesDistDPhiRDzExactA, matricesDistDzExactA, matricesCorrDrDzExactA,
                      matricesCorrDPhiRDzExactA, matricesCorrDzExactA, rRow, zColumn, phiSlice, 0, c0, c1, ezField,
                      dvdE);

	  InitPotentialAndCharge3D(vTestFunction, rhoTestFunction, erTestFunction, ePhiTestFunction, ezTestFunction,
                      intErDzTestFunction, intEPhiRDzTestFunction, intDzTestFunction, matricesVExactC,
                      matricesChargeC, matricesErExactC, matricesEPhiExactC, matricesEzExactC, matricesDistDrDzExactC,
                      matricesDistDPhiRDzExactC, matricesDistDzExactC, matricesCorrDrDzExactC,
                      matricesCorrDPhiRDzExactC, matricesCorrDzExactC, rRow, zColumn, phiSlice, 1, c0, c1, ezField,
                      dvdE);

	  scExact->InitSpaceCharge3DPoissonIntegralDz(rRow, zColumn, phiSlice, maxIter, convergenceError,
                                              matricesErExactA, matricesEPhiExactA, matricesEzExactA,
                                              matricesErExactC, matricesEPhiExactC, matricesEzExactC,
                                              matricesDistDrDzExactA, matricesDistDPhiRDzExactA, matricesDistDzExactA,
                                                                                                                                                                                                                                                                                        matricesCorrDrDzExactA, matricesCorrDPhiRDzExactA, matricesCorrDzExactA,
                                              matricesDistDrDzExactC, matricesDistDPhiRDzExactC, matricesDistDzExactC,
                                              matricesCorrDrDzExactC, matricesCorrDPhiRDzExactC, matricesCorrDzExactC,
                                              intErDzTestFunction, intEPhiRDzTestFunction, intDzTestFunction);

  }
  printf("creating distortion tree for sc\n");
  sc->CreateDistortionTree(rRowTest,zColumnTest,phiSliceTest);

  if (correctionType == 0) {
	  printf("creating distortion tree for scExact\n");
	  scExact->CreateDistortionTree(rRowTest,zColumnTest,phiSliceTest);
  }




  printf("================= Correctness Testing (numeric vs analytic solution) ========================\n");
  printf("v(r,phi,z)\t\t= %-100s\n",  (vTestFunction->GetExpFormula()).Data());
  printf("rho(r,phi,z)\t= %-100s\n", (rhoTestFunction->GetExpFormula()).Data());
  printf("---------------------------------------------------------------------------------------------\n");
  printf("%-50.50s%-15s%-15s%-15s\n","var","mean","rms","max");
  printf("---------------------------------------------------------------------------------------------\n");


  Int_t indexErrorList = 0;
  Int_t varNameId = 0;
  Double_t errorList[200];
  Int_t oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"rho","rhoAnalytic",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);


  oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"v","vAnalytic",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"ER","ERAnalytic",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"EPhi","EPhiAnalytic",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"EZ","EZAnalytic",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;

  if (correctionType == 0) {
	  GetResidueFromDistortionAnalyticTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),TString::Format("distortion%s.root",scExact->GetName()),"drLocalDist",errorList,indexErrorList, pcStream);
	  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
	  oldIndexErrorList = indexErrorList;

	  GetResidueFromDistortionAnalyticTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),TString::Format("distortion%s.root",scExact->GetName()),"drPhiLocalDist",errorList,indexErrorList, pcStream);
	  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
	  oldIndexErrorList = indexErrorList;

	  GetResidueFromDistortionAnalyticTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),TString::Format("distortion%s.root",scExact->GetName()),"dzLocalDist",errorList,indexErrorList, pcStream);
	  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
	  oldIndexErrorList = indexErrorList;

	  GetResidueFromDistortionAnalyticTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),TString::Format("distortion%s.root",scExact->GetName()),"drDist",errorList,indexErrorList, pcStream);
	  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
	  oldIndexErrorList = indexErrorList;

	  GetResidueFromDistortionAnalyticTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),TString::Format("distortion%s.root",scExact->GetName()),"drPhiDist",errorList,indexErrorList, pcStream);
	  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
	  oldIndexErrorList = indexErrorList;

	  GetResidueFromDistortionAnalyticTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),TString::Format("distortion%s.root",scExact->GetName()),"dzDist",errorList,indexErrorList, pcStream);
	  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
	  oldIndexErrorList = indexErrorList;
  } else {
          varNameId += 6;
  }

  printf("\n\n");
  printf("================= Consistency testing  (distortion - correction) ============================\n");
  if (correctionType == 0)
    printf("Correction type: case regular grid interpolation\n");
  else
    printf("Correction type: case irregular grid interpolation\n");
  printf("---------------------------------------------------------------------------------------------\n");
  printf("%-50.50s%-15s%-15s%-15s\n","var","mean","rms","max");
  printf("---------------------------------------------------------------------------------------------\n");
  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"-drCorr","drDist",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;
  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"-drPhiCorr","drPhiDist",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(0,TString::Format("distortion%s.root",sc->GetName()).Data(),"-dzCorr","dzDist",errorList,indexErrorList, pcStream);
  WriteErrorToPCStream(0,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;


  //(*pcStream) << "\n";

  for (Int_t m = 0;m <phiSlice;m++) {
    delete matricesVExactA[m];
    delete matricesChargeA[m];
    delete matricesVExactC[m];
    delete matricesChargeC[m];
    delete matricesErExactA[m];
    delete matricesEPhiExactA[m];
    delete matricesEzExactA[m];
    delete matricesErExactC[m];
    delete matricesEPhiExactC[m];
    delete matricesEzExactC[m];
    delete matricesDistDrDzExactA[m];
    delete matricesDistDPhiRDzExactA[m];
    delete matricesDistDzExactA[m];
    delete matricesCorrDrDzExactA[m];
    delete matricesCorrDPhiRDzExactA[m];
    delete matricesCorrDzExactA[m];
    delete matricesDistDrDzExactC[m];
    delete matricesDistDPhiRDzExactC[m];
    delete matricesDistDzExactC[m];
    delete matricesCorrDrDzExactC[m];
    delete matricesCorrDPhiRDzExactC[m];
    delete matricesCorrDzExactC[m];
  }

  delete scExact;
  delete sc;

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestCorrectnessDistortion", "End");

}




/// unit test for Consistency
/// - a known potential function at boundary, charge distribution
/// - generate numerical calculation of distortions and corrections in sc
/// - calc
/// \param correctionType Int_t 0-> use regular interpolation, 1->use irregularinterpolation
void UnitTestConsistencyDistortionZShort(const Int_t rRow, const Int_t zColumn, const Int_t phiSlice, const Int_t rRowTest, const Int_t zColumnTest, const Int_t phiSliceTest, const Int_t correctionType, TTreeSRedirector *pcStream) {
  const Int_t maxIter = 100;
  const Float_t convergenceError = 1e-8;

  const Int_t nPoints = phiSliceTest * rRowTest * zColumnTest;
  const Double_t maxEpsilon = 1e-2;

  const Double_t zShort = 40;
  const Double_t amplitude = -1000;

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyDistortionZShort", "Begin");
  if (correctionType == 0)
    ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyDistortionZShort", "Regular interpolation");
  else
    ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyDistortionZShort", "Irregular interpolation");

  // set potential boundary in V
  TF1 * potentialBoundaryFunctionInZ = new TF1("dFunctionVZ", dFunctionVZ, -250.0, 250.0,2);
  potentialBoundaryFunctionInZ->SetParName(0,"zShort");
  potentialBoundaryFunctionInZ->SetParName(1,"amplitude");
  potentialBoundaryFunctionInZ->SetParameter(0,zShort);
  potentialBoundaryFunctionInZ->SetParameter(1,amplitude);
  // charge in histogram
  TH3 * chargeA = new TH3F ("charge A","charge A",250, 0.0, TMath::TwoPi(), 250, 85.0, 250.0, 250,0.0,250.0);
  TH3 * chargeC = new TH3F ("charge C","charge C",250, 0.0, TMath::TwoPi(), 250, 85.0, 250.0, 250,0.0,250.0);

  AliTPCSpaceCharge3DDriftLine * sc = new AliTPCSpaceCharge3DDriftLine(TString::Format("distortionIFCDeltaPotentialRRow%dZCol%dPhiSlice%d",rRow,zColumn,phiSlice).Data(),TString::Format("distortionIFCDeltaPotentialRRow%dZCol%dPhiSlice%d",rRow,zColumn,phiSlice).Data(), rRow,zColumn,phiSlice);

  // set the input boundary potential and charge densities
  sc->SetInputSpaceCharge(chargeA,0);
  sc->SetInputSpaceCharge(chargeC,1);
  sc->SetBoundaryIFCA(potentialBoundaryFunctionInZ);
  sc->SetOmegaTauT1T2(-0.35,1.,1.);

  sc->SetCorrectionType(correctionType);

  sc->InitSpaceCharge3DPoissonIntegralDz(rRow,zColumn,phiSlice,gkMaxIter,gkConvError);

  ::Info("AliTPCSpaceCharge3DDriftLineTest::UnitTestConsistencyDistortionZShort", "Creating test points");
  sc->CreateDistortionTree(rRowTest,zColumnTest,phiSliceTest);


  printf("================= Consistency testing  (distortion - correction) ============================\n");
  printf("v(r,phi,z)\t\t= %-100s\n", "(z < z0) ?-k0 * (TMath::Abs(z) / TMath::Abs(z0)): -k0 + k1 * (TMath::Abs(z) - z0);");
  printf("zShort = %.3f cm, V (max)=%.f V\n",zShort,amplitude);

  if (correctionType == 0)
    printf("Correction type: case regular grid interpolation\n");
  else
    printf("Correction type: case irregular grid interpolation\n");
  printf("---------------------------------------------------------------------------------------------\n");
  printf("%-50.50s%-15s%-15s%-15s\n","var","mean","rms","max");
  printf("---------------------------------------------------------------------------------------------\n");

  Int_t indexErrorList = 90;
  Double_t errorList[200];
  Int_t varNameId = 11;
  Int_t oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(1,TString::Format("distortion%s.root",sc->GetName()).Data(),"-drCorr","drDist",errorList,indexErrorList,pcStream);
  WriteErrorToPCStream(1,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;

  GetResidueFromDistortionTree(1,TString::Format("distortion%s.root",sc->GetName()).Data(),"-drPhiCorr","drPhiDist",errorList,indexErrorList,pcStream);
  WriteErrorToPCStream(1,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  oldIndexErrorList = indexErrorList;
  GetResidueFromDistortionTree(1,TString::Format("distortion%s.root",sc->GetName()).Data(),"-dzCorr","dzDist",errorList,indexErrorList,pcStream);
  WriteErrorToPCStream(1,rRow,zColumn,phiSlice,rRowTest,zColumnTest,phiSliceTest,correctionType,varNameId++,errorList,oldIndexErrorList,pcStream);
  //printf("varNameId = %d\n",varNameId);
  //oldIndexErrorList = indexErrorList;






  delete chargeA;
  delete chargeC;
  delete sc;
}


/// for exact functions
///
/// \param vTestFunction
/// \param rhoTestFunction
/// \param erTestFunction
/// \param ePhiTestFunction
/// \param ezTestFunction
/// \param intErDzTestFunction
/// \param intEPhiRDzTestFunction
/// \param intDzTestFunction
/// \param matricesV
/// \param matricesCharge
/// \param matricesEr
/// \param matricesEPhi
/// \param matricesEz
/// \param matricesDistDrDz
/// \param matricesDistDPhiRDz
/// \param matricesDistDz
/// \param matricesCorrDrDz
/// \param matricesCorrDPhiRDz
/// \param matricesCorrDz
/// \param rRow
/// \param zColumn
/// \param phiSlice
/// \param side
/// \param c0
/// \param c1
/// \param ezField
/// \param dvdE
void InitPotentialAndCharge3D(TFormula *vTestFunction, TFormula *rhoTestFunction, TFormula *erTestFunction,
                         TFormula *ePhiTestFunction, TFormula *ezTestFunction, TFormula *intErDzTestFunction,
                         TFormula *intEPhiRDzTestFunction, TFormula *intDzTestFunction, TMatrixD **matricesV,
                         TMatrixD **matricesCharge, TMatrixD **matricesEr, TMatrixD **matricesEPhi,
                         TMatrixD **matricesEz, TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz,
                         TMatrixD **matricesDistDz, TMatrixD **matricesCorrDrDz, TMatrixD **matricesCorrDPhiRDz,
                         TMatrixD **matricesCorrDz, const Int_t rRow, const Int_t zColumn, const Int_t phiSlice,
                         const Int_t side, const Double_t c0, const Double_t c1, const Double_t ezField,
                         const Double_t dvdE) {
  Double_t rList[rRow], zedList[zColumn], phiList[phiSlice];
  Float_t ofcRadius = AliTPCPoissonSolver::fgkOFCRadius;
  Float_t ifcRadius = AliTPCPoissonSolver::fgkIFCRadius;
  Float_t tpcZ0 = AliTPCPoissonSolver::fgkTPCZ0;

  Float_t dRadius = (AliTPCPoissonSolver::fgkOFCRadius - AliTPCPoissonSolver::fgkIFCRadius) / (rRow - 1);
  Float_t dZ = AliTPCPoissonSolver::fgkTPCZ0 / (zColumn - 1);
  Float_t dPhi = TMath::TwoPi() / phiSlice;



  // list points on grid in cm
  for (Int_t k = 0; k < phiSlice; k++)
    phiList[k] = dPhi * k;
  for (Int_t i = 0; i < rRow; i++)
    rList[i] = AliTPCPoissonSolver::fgkIFCRadius + i * dRadius;
  for (Int_t j = 0; j < zColumn; j++)
    zedList[j] = j * dZ;

  TMatrixD *arrayV;
  TMatrixD *charge;
  TMatrixD *er, *ePhi, *ez;

  Double_t radius0, z0, phi0;

  for (Int_t k = 0; k < phiSlice; k++) {
    arrayV = matricesV[k];
    charge = matricesCharge[k];
    er = matricesEr[k];
    ePhi = matricesEPhi[k];
    ez = matricesEz[k];
    phi0 = phiList[k];

    /// Fill the non-boundary values
    for (Int_t i = 0; i < rRow; i++) {
      radius0 = rList[i];
      for (Int_t j = 0; j < zColumn; j++) {
        z0 = zedList[j];

        // Exact solution
        (*arrayV)(i, j) = vTestFunction->Eval(radius0, phi0, z0);
        (*charge)(i, j) = -1.0 * rhoTestFunction->Eval(radius0, phi0, z0);


        (*er)(i, j) = -1.0 * erTestFunction->Eval(radius0, phi0, z0);
        (*ePhi)(i, j) = -1.0 * ePhiTestFunction->Eval(radius0, phi0, z0);
        (*ez)(i, j) = -1.0 * ezTestFunction->Eval(radius0, phi0, z0);
        //if (side == 1) (*ez)(i, j) = -(*ez)(i, j);
      } // end j
    } // end i
  } // end phi


  LocalDistCorrDzExact(intErDzTestFunction, intEPhiRDzTestFunction, intDzTestFunction, rList, phiList, zedList,
                       matricesDistDrDz, matricesDistDPhiRDz, matricesDistDz, matricesCorrDrDz, matricesCorrDPhiRDz,
                       matricesCorrDz, rRow, zColumn, phiSlice, c0, c1, ezField, dvdE, side);
}

/// for correctness analysis
///
/// \param intDrDzF
/// \param intDPhiDzF
/// \param intDzDzF
/// \param rList
/// \param phiList
/// \param zList
/// \param matricesDistDrDz
/// \param matricesDistDPhiRDz
/// \param matricesDistDz
/// \param matricesCorrDrDz
/// \param matricesCorrDPhiRDz
/// \param matricesCorrDz
/// \param rRow
/// \param zColumn
/// \param phiSlice
/// \param c0
/// \param c1
/// \param ezField
/// \param dvdE
/// \param side
void
LocalDistCorrDzExact(TFormula *intDrDzF, TFormula *intDPhiDzF, TFormula *intDzDzF, Double_t *rList, Double_t *phiList,
                     Double_t *zList, TMatrixD **matricesDistDrDz, TMatrixD **matricesDistDPhiRDz,
                     TMatrixD **matricesDistDz, TMatrixD **matricesCorrDrDz, TMatrixD **matricesCorrDPhiRDz,
                     TMatrixD **matricesCorrDz, const Int_t rRow, const Int_t zColumn, const Int_t phiSlice,
                     const Double_t c0, const Double_t c1, const Double_t ezField, const Double_t dvdE,
                     const Int_t side) {
  Double_t r0, z0, phi0, z1;
  Float_t localIntErOverEz = 0.0;
  Float_t localIntEPhiOverEz = 0.0;
  Float_t localIntDeltaEz = 0.0;

  // pointer declaration
  TMatrixD *distDrDz;
  TMatrixD *distDPhiRDz;
  TMatrixD *distDz;
  TMatrixD *corrDrDz;
  TMatrixD *corrDPhiRDz;
  TMatrixD *corrDz;


  for (Int_t m = 0; m < phiSlice; m++) {
    distDrDz = matricesDistDrDz[m];
    distDPhiRDz = matricesDistDPhiRDz[m];
    distDz = matricesDistDz[m];

    corrDrDz = matricesCorrDrDz[m];
    corrDPhiRDz = matricesCorrDPhiRDz[m];
    corrDz = matricesCorrDz[m];

    for (Int_t i = 0; i < rRow; i++) {
      (*distDrDz)(i, zColumn - 1) = 0.0;
      (*distDPhiRDz)(i, zColumn - 1) = 0.0;
      (*distDz)(i, zColumn - 1) = 0.0;

      (*corrDrDz)(i, 0) = 0.0;
      (*corrDPhiRDz)(i, 0) = 0.0;
      (*corrDz)(i, 0) = 0.0;
    }

  }

  // for this case
  // use trapezoidal rule assume no ROC displacement
  for (Int_t m = 0; m < phiSlice; m++) {
    phi0 = phiList[m];

    distDrDz = matricesDistDrDz[m];
    distDPhiRDz = matricesDistDPhiRDz[m];
    distDz = matricesDistDz[m];

    corrDrDz = matricesCorrDrDz[m];
    corrDPhiRDz = matricesCorrDPhiRDz[m];
    corrDz = matricesCorrDz[m];

    for (Int_t j = 0; j < zColumn - 1; j++) {
      z0 = zList[j];
      z1 = zList[j + 1];
      for (Int_t i = 0; i < rRow; i++) {
        r0 = rList[i];

        if (side == 0) {
          localIntErOverEz = (intDrDzF->Eval(r0, phi0, z1) - intDrDzF->Eval(r0, phi0, z0)) / (-1 * ezField);
          localIntEPhiOverEz = (intDPhiDzF->Eval(r0, phi0, z1) - intDPhiDzF->Eval(r0, phi0, z0)) / (-1 * ezField);
          localIntDeltaEz = intDzDzF->Eval(r0, phi0, z1) - intDzDzF->Eval(r0, phi0, z0);
        } else {
          localIntErOverEz = (intDrDzF->Eval(r0, phi0, z1) - intDrDzF->Eval(r0, phi0, z0)) / (-1 * ezField);
          localIntEPhiOverEz = (intDPhiDzF->Eval(r0, phi0, z1) - intDPhiDzF->Eval(r0, phi0, z0)) / (-1 * ezField);
          localIntDeltaEz = intDzDzF->Eval(r0, phi0, z1) - intDzDzF->Eval(r0, phi0, z0);
        }

        (*distDrDz)(i, j) = c0 * localIntErOverEz + c1 * localIntEPhiOverEz;
        (*distDPhiRDz)(i, j) = c0 * localIntEPhiOverEz - c1 * localIntErOverEz;
        (*distDz)(i, j) = localIntDeltaEz * dvdE * dvdE; // two times?


        (*corrDrDz)(i, j + 1) = -1 * (*distDrDz)(i, j);
        (*corrDPhiRDz)(i, j + 1) = -1 * (*distDPhiRDz)(i, j);
        (*corrDz)(i, j + 1) = -1 * (*distDz)(i, j);
      }
    }
  }
}




// function for dV(z)
Double_t dFunctionVZ(Double_t *x, Double_t *par)
{
  /// dU=0 at (0 and at 250)
  /// -k0*z     at z(0,z0)
  /// -k0*z0+k1*z   at z(z0,250)
  /// -k0*z0+k1*250=0
  ///  k1=k0*z0/250.

  Float_t xx =x[0];
  const Float_t z0 =par[0]; // U at central electrode
  const Float_t k0 =par[1]; // z0
  //const Float_t k1 = k0 * z0 / AliTPCPoissonSolver::fgkTPCZ0; // k1

  const Float_t k1 = k0  / (AliTPCPoissonSolver::fgkTPCZ0 - z0);
  Double_t f = 0.0;

  // z location near
  if (TMath::Abs(xx) < TMath::Abs(z0)) {
    // ideal case
    f = -k0 * (TMath::Abs(xx) / TMath::Abs(z0)) ;
  } else {
    f = -k0 + k1 * (TMath::Abs(xx) - z0);
  }

  return f;
}


// open distortion tree and print result
//
void GetResidueFromDistortionTree(Int_t unitTestId, const char * fileName,const char *numericName, const char *analyticName, Double_t *errorList, Int_t &indexErrorList,TTreeSRedirector *pcStream) {
  TFile fileNumeric(fileName);
  TTree *treeNumeric = (TTree *)fileNumeric.Get("distortion");


  Int_t nPoint = treeNumeric->GetEntries();
  // Potential
  treeNumeric->SetEstimate(nPoint);
  treeNumeric->Draw(TString::Format("abs(%s)",analyticName),"", "goff");
  Double_t *numericList = treeNumeric->GetV1();
  Double_t maxVar = TMath::MaxElement(nPoint, numericList);
  Double_t meanVar = TMath::Mean(nPoint, numericList);
  Double_t rmsVar = TMath::RMS(nPoint, numericList);

  errorList[indexErrorList++] = meanVar;
  errorList[indexErrorList++] = rmsVar;
  errorList[indexErrorList++] = maxVar;

  printf("%-50.50s%-15.3E%-15.3E%-15.3E\n",TString::Format("abs(%s)",analyticName).Data(),meanVar,rmsVar,maxVar);

  Double_t maxVarAnalytic = maxVar;
  if (numericName[0] == 'E')
    maxVarAnalytic = 400;


  treeNumeric->SetEstimate(nPoint);
  treeNumeric->Draw(TString::Format("abs(%s - %s)/%f",analyticName,numericName,maxVarAnalytic),"", "goff");
  numericList = treeNumeric->GetV1();
  maxVar = TMath::MaxElement(nPoint, numericList);
  meanVar = TMath::Mean(nPoint, numericList);
  rmsVar = TMath::RMS(nPoint, numericList);
  if (numericName[0] == 'E')
    printf("%-50.50s%-15.3E%-15.3E%-15.3E\n",TString::Format("abs(%s - %s)/%f",analyticName,numericName,maxVarAnalytic).Data(),meanVar,rmsVar,maxVar);
  else
    printf("%-50.50s%-15.3E%-15.3E%-15.3E\n",TString::Format("abs(%s - %s)/max(abs(%s))",analyticName,numericName,analyticName).Data(),meanVar,rmsVar,maxVar);

  errorList[indexErrorList++] = meanVar;
  errorList[indexErrorList++] = rmsVar;
  errorList[indexErrorList++] = maxVar;


  treeNumeric->SetEstimate(nPoint);
  treeNumeric->Draw(TString::Format("abs(%s - %s)",analyticName,numericName),"", "goff");
  numericList = treeNumeric->GetV1();
  maxVar = TMath::MaxElement(nPoint, numericList);
  meanVar = TMath::Mean(nPoint, numericList);
  rmsVar = TMath::RMS(nPoint, numericList);

  errorList[indexErrorList++] = meanVar;
  errorList[indexErrorList++] = rmsVar;
  errorList[indexErrorList++] = maxVar;

  printf("%-50.50s%-15.3E%-15.3E%-15.3E\n",TString::Format("abs(%s - %s)",analyticName,numericName).Data(),meanVar,rmsVar,maxVar);
  treeNumeric->SetEstimate(nPoint);



  printf("---------------------------------------------------------------------------------------------\n");
  fileNumeric.Close();
}



// open distortion tree from numeric and analytic
//
void GetResidueFromDistortionAnalyticTree(Int_t unitTestId, const char * numericFileName, const char * analyticFileName, const char *varName, Double_t *errorList, Int_t &indexErrorList,TTreeSRedirector *pcStream) {
  TFile fileNumeric(numericFileName);
  TTree *treeNumeric = (TTree *)fileNumeric.Get("distortion");
  Int_t nPoint = treeNumeric->GetEntries();


  treeNumeric->AddFriend("analytic = distortion",analyticFileName);

  treeNumeric->Draw(TString::Format("abs(analytic.%s)",varName).Data(), "", "goff");
  Double_t *numericList = treeNumeric->GetV1();
  Double_t maxVar = TMath::MaxElement(nPoint, numericList);
  Double_t meanVar = TMath::Mean(nPoint, numericList);
  Double_t rmsVar = TMath::RMS(nPoint, numericList);
  printf("%-50.50s%-15.3E%-15.3E%-15.3E\n",TString::Format("abs(analytic.%s)",varName).Data(),meanVar,rmsVar,maxVar);

  errorList[indexErrorList++] = meanVar;
  errorList[indexErrorList++] = rmsVar;
  errorList[indexErrorList++] = maxVar;

  Double_t maxVarAnalytic = maxVar;

  treeNumeric->SetEstimate(nPoint);
  treeNumeric->Draw(TString::Format("abs(analytic.%s - %s)/%f",varName,varName,maxVarAnalytic),"", "goff");
  numericList = treeNumeric->GetV1();
  maxVar = TMath::MaxElement(nPoint, numericList);
  meanVar = TMath::Mean(nPoint, numericList);
  rmsVar = TMath::RMS(nPoint, numericList);

  errorList[indexErrorList++] = meanVar;
  errorList[indexErrorList++] = rmsVar;
  errorList[indexErrorList++] = maxVar;


  printf("%-50.50s%-15.3E%-15.3E%-15.3E\n",TString::Format("abs(analytic.%s - %s)/",varName,varName).Data(),meanVar,rmsVar,maxVar);
  printf("%-50.50s\n",TString::Format("max(abs(analytic.%s))",varName).Data());



  treeNumeric->SetEstimate(nPoint);
  treeNumeric->Draw(TString::Format("abs(analytic.%s - %s)",varName,varName),"", "goff");
  numericList = treeNumeric->GetV1();
  maxVar = TMath::MaxElement(nPoint, numericList);
  meanVar = TMath::Mean(nPoint, numericList);
  rmsVar = TMath::RMS(nPoint, numericList);
  errorList[indexErrorList++] = meanVar;
  errorList[indexErrorList++] = rmsVar;
  errorList[indexErrorList++] = maxVar;

  printf("%-50.50s%-15.3E%-15.3E%-15.3E\n",TString::Format("abs(analytic.%s - %s)",varName,varName).Data(),meanVar,rmsVar,maxVar);

  printf("---------------------------------------------------------------------------------------------\n");
  fileNumeric.Close();
}



/// write error to a tree
///
void WriteErrorToPCStream(Int_t unitTestId, Int_t rRow, Int_t zColumn, Int_t phiSlice, Int_t rRowTest, Int_t zColumnTest, Int_t phiSliceTest, Int_t correctionType, Int_t varNameId, Double_t *errorList, Int_t  indexErrorList,TTreeSRedirector *pcStream) {

  (*pcStream) << "residue" <<
              "unitTestId=" << unitTestId <<
              "rRow=" << rRow   << "zColumn=" << zColumn  << "phiSlice=" << phiSlice <<
              "rRowTest=" << rRowTest << "zColumnTest=" << zColumnTest << "phiSliceTest=" << phiSliceTest <<
              "correctionType=" << correctionType << "varNameId=" << varNameId <<
              "varMean=" << errorList[indexErrorList]  <<
              "varRms=" << errorList[indexErrorList + 1] <<
              "varMax=" << errorList[indexErrorList + 2] <<
              "errorRelativeMean=" << errorList[indexErrorList + 3]  <<
              "errorRelativeRms=" << errorList[indexErrorList + 4] <<
              "errorRelativeMax=" << errorList[indexErrorList + 5] <<
              "errorAbsMean=" << errorList[indexErrorList + 6]  <<
              "errorAbsRms=" << errorList[indexErrorList + 7] <<
              "errorAbsMax=" << errorList[indexErrorList + 8] <<
              "\n";

}


/// print error status
///
void PrintErrorStatus() {
  const Double_t maxEpsilon = 1e-2;

  TFile file("spaceChargeDriftLinePerformance.root");
  TTree *t = (TTree *)file.Get("residue");


  const Int_t numberOfUnitTest = 2;
  const Int_t numberOfCorrectionType = 2;

  const char * unitTestName[] = {"UnitTestCorrectnessDistortion","UnitTestCorrectnessDistortionZShort"};
  const char * correctionTypeName[] = {"Regular interpolation","Irregular interpolation"};
  const char * varName[] = {"rho","v","ER","EPhi","EZ","drLocalDist","drPhiLocalDist","dzLocalDist","drDist","drPhiDist","dzDist","drDist+drCorr","drPhiDist+drPhiCorr","dzDist+dzCorr"};

  Int_t unitTestId, correctionType, varNameId;
  Double_t varMean,errorRelativeMean,errorAbsMean;
  Double_t varMax,errorRelativeMax,errorAbsMax;
  Double_t varRms,errorRelativeRms,errorAbsRms;
  t->SetBranchAddress("unitTestId",&unitTestId);
  t->SetBranchAddress("varNameId",&varNameId);
  t->SetBranchAddress("correctionType",&correctionType);

  t->SetBranchAddress("varMean",&varMean);
  t->SetBranchAddress("varRms",&varRms);
  t->SetBranchAddress("varMax",&varMax);
  t->SetBranchAddress("errorRelativeMean",&errorRelativeMean);
  t->SetBranchAddress("errorRelativeRms",&errorRelativeRms);
  t->SetBranchAddress("errorRelativeMax",&errorRelativeMax);
  t->SetBranchAddress("errorAbsMean",&errorAbsMean);
  t->SetBranchAddress("errorAbsRms",&errorAbsRms);
  t->SetBranchAddress("errorAbsMax",&errorAbsMax);


  ::Info(TString::Format("AliTPCSpaceCharge3DDriftLineTest::* (%d)",0).Data(),"%s",correctionTypeName[0]);
  ::Info(TString::Format("AliTPCSpaceCharge3DDriftLineTest::* (%d)",1).Data(),"%s",correctionTypeName[1]);

  for (Int_t entry = 0;  entry < t->GetEntries(); entry++) {
    t->GetEntry(entry);
    if (varNameId  == 0) continue;
    if (varNameId < 2 || varNameId > 10) {	
	    if ( errorRelativeMean < maxEpsilon)
	      ::Info(TString::Format("AliTPCSpaceCharge3DDriftLineTest::%-30.30s (%d) %-20.20s",unitTestName[unitTestId],correctionType,varName[varNameId]).Data(),      "Test OK: Mean Relative Error=%.2E < %.2E", errorRelativeMean, maxEpsilon);
	    else
	      ::Error(TString::Format("AliTPCSpaceCharge3DDriftLineTest::%-30.30s (%d) %-20.20s",unitTestName[unitTestId],correctionType,varName[varNameId]).Data(),"Test FAILED: Mean Relative Error=%.2E > %.2E",errorAbsMean, maxEpsilon);
   } else if (varNameId < 5) {
	    if ( errorRelativeMean < maxEpsilon)
	      ::Info(TString::Format("AliTPCSpaceCharge3DDriftLineTest::%-30.30s (%d) %-20.20s",unitTestName[unitTestId],correctionType,varName[varNameId]).Data(),      "Test OK: Mean Absolute Error/400.0=%.2E < %.2E", errorRelativeMean, maxEpsilon);
	    else
	      ::Error(TString::Format("AliTPCSpaceCharge3DDriftLineTest::%-30.30s (%d) %-20.20s",unitTestName[unitTestId],correctionType,varName[varNameId]).Data(),"Test FAILED: Mean Absolute Error/400.0=%.2E > %.2E",errorRelativeMean, maxEpsilon);
   }  else {
	    if ( errorAbsMean < maxEpsilon)
	      ::Info(TString::Format("AliTPCSpaceCharge3DDriftLineTest::%-30.30s (%d) %-20.20s",unitTestName[unitTestId],correctionType,varName[varNameId]).Data(),      "Test OK: Mean Absolute Error=%.2E < %.2E", errorAbsMean, maxEpsilon);
	    else
	      ::Error(TString::Format("AliTPCSpaceCharge3DDriftLineTest::%-30.30s (%d) %-20.20s",unitTestName[unitTestId],correctionType,varName[varNameId]).Data(),"Test FAILED: Mean Absolute Error=%.2E > %.2E",errorAbsMean, maxEpsilon);
   }  	   
  }

  ::Info(TString::Format("AliTPCSpaceCharge3DDriftLineTest").Data(),"Finer grid test report can be found in https://alice.its.cern.ch/jira/browse/ATO-433");
}
