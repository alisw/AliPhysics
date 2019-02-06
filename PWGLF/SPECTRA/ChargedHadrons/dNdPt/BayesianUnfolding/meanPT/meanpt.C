#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THn.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TObjArray.h"
#include "TGFrame.h"
#include "TGFileBrowser.h"



#include <string>
#include <iostream>


using namespace std;


void meanpt(string fileName, string inputPath = "Input/Input_", string outputPath = "Histograms/"){

  // Settings
  const Int_t nIter = 10;
//  string nIterstr = "10";
//  string iterCheckHistName = "momentUnfolded1_nIter_" + nIterstr;
  TStopwatch* stopwatch = new TStopwatch();

  Bool_t useSecondMC = kFALSE;
  Bool_t useExtrapolation = kFALSE;
  Bool_t useModulatedResponse = kFALSE;
  Bool_t useSmoothing = kFALSE;
  Bool_t includeSystematics = kTRUE; //swich that off when calculating track systematix
  Bool_t includeSimulation = kFALSE;
//  if(outputPath == "Histograms/") includeSystematics = kTRUE;

  /// Define input and output files
//  string inputFileResponseName = inputPath + fileName + "_Pythia6" + ".root";
  string inputFileResponseName = inputPath + fileName + "_EPOS" + ".root";
  TFile* inputFileResponse = NULL;

  string inputFileName = inputPath + fileName + ".root";
  string outputFileName = outputPath + fileName + ".root";
  string inputFileNameExtrapol = "Histograms/" + fileName + "_extrapol.root";
  string inputFileNameSyst = outputPath+ fileName + "_Sys.root";

  /// Include RooUnfold Software
  //  (first make library... ; TUnfold Interface does not work with newest root version, needs to be disabled in makefile)
  gSystem->Load("$MEANPT/RooUnfold-1.1.1/libRooUnfold");
  stopwatch->Start();

  ///----------------------------------------------------------------------------------------------------
  /// Retrieve Histograms
  ///----------------------------------------------------------------------------------------------------

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");

  /// :::::::::: Data Input Histograms ::::::::::

  TH1D* multDistMeasured = (TH1D*) inputFile->FindObjectAny("multDistMeasured");
  if(!multDistMeasured) cout << "ERROR: Could not find multDistMeasured" << endl;

  TH2D* multPtMeasured = (TH2D*) inputFile->FindObjectAny("multPtMeasured");
  if(!multPtMeasured) cout << "ERROR: Could not find multPtMeasured" << endl;

  /// ::::::::::: Detector Response Histograms :::::::::::

  if(useSecondMC) {inputFileResponse = TFile::Open(inputFileResponseName.c_str(),"READ"); cout << endl << "-> using detector response from: " << inputFileResponse->GetName() << endl;}
  else {inputFileResponse = inputFile;}

  TH2D* responseMatrixOrig = (TH2D*) inputFileResponse->FindObjectAny("responseMatrixOrig");
  if(!responseMatrixOrig) cout << "ERROR: Could not find responseMatrixOrig" << endl;

  TH2D* responseMatrix = (TH2D*) inputFileResponse->FindObjectAny("responseMatrix");
  if(!responseMatrix) cout << "ERROR: Could not find responseMatrix" << endl;

  TH2D* responseMatrixTracksOrig = (TH2D*) inputFileResponse->FindObjectAny("responseMatrixTracksOrig");
  if(!responseMatrixTracksOrig) cout << "ERROR: Could not find responseMatrixTracksOrig" << endl;

  TH2D* responseMatrixTracks = (TH2D*) inputFileResponse->FindObjectAny("responseMatrixTracks");
  if(!responseMatrixTracks) cout << "ERROR: Could not find responseMatrixTracks" << endl;

  /// ::::::::::: MC Input Histograms :::::::::::

  TH2D* multPtMeasuredMC = (TH2D*) inputFile->FindObjectAny("multPtMeasuredMC");
  TH2D* multPtMeasuredPtResUnfoldedMC = (TH2D*) inputFile->FindObjectAny("multPtMeasuredPtResUnfoldedMC");
//  multPtMeasuredMC->SetName("multPtMeasuredMC");
  if(!multPtMeasuredMC) cout << "ERROR: Could not find multPtMeasuredMC" << endl;

  TH2D* multPtGeneratedMC = (TH2D*) inputFile->FindObjectAny("multPtGeneratedMC");
  if(!multPtGeneratedMC) cout << "ERROR: Could not find multPtGeneratedMC" << endl;

//  TH2D* multPtMultGen = (TH2D*) inputFile->FindObjectAny("multPtMultGen");
//  if(!multPtMultGen) cout << "ERROR: Could not find multPtMultGen" << endl;

  TH2D* ptResolution = (TH2D*) inputFile->FindObjectAny("ptResolution");
  if(!ptResolution) cout << "ERROR: Could not find ptResolution" << endl;
  TH2D* ptResponseMatrix = (TH2D*) inputFile->FindObjectAny("ptResponseMatrix");
  if(!ptResponseMatrix) cout << "ERROR: Could not find ptResponseMatrix" << endl;

  ///----------------------------------------------------------------------------------------------------
  /// Unfolding for Data
  ///----------------------------------------------------------------------------------------------------

  // n(N_ch) : Unfolded Multiplicity Distribution after nIter iterations
  cout  << endl << "-> Unfolding measured mutliplicity distribution (data)...." << endl;
  TH1D* multDistUnfolded = getUnfoldedMultDist(multDistMeasured, responseMatrixOrig, nIter, useSmoothing);

  // P(N_ch|N_acc) : Unfolding Matrix based on the unfolded multiplicity distribution
  TH2D* unfoldingMatrix = getUnfoldingMatrix(responseMatrix, multDistUnfolded);

  // Unfolded pT Spectra (pT vs. N_ch)
  cout  << endl << "-> Unfolding measured pT spectra (data)...." << endl;
  TH2D* multPtUnfolded = unfoldPtSpectra(multPtMeasured, responseMatrixTracksOrig, nIter, useSmoothing);

  // Measured Moments <(pt)^n>(N_acc)
  TH1D* momentMeasured1 = getMomentFromSpectra(multPtMeasured, 1);
  TH1D* momentMeasured2 = getMomentFromSpectra(multPtMeasured, 2);
  TH1D* momentMeasured3 = getMomentFromSpectra(multPtMeasured, 3);
  TH1D* momentMeasuredRMS = calcRMS(momentMeasured2);

  // Reweighted Moments <(pt)^n>(N_ch) from <(pt)^n>(N_acc) using response matrix
  TH1D* momentReweighted1 = reweightMoment(momentMeasured1, responseMatrix, 1);
  TH1D* momentReweighted2 = reweightMoment(momentMeasured2, responseMatrix, 2);
  TH1D* momentReweighted3 = reweightMoment(momentMeasured3, responseMatrix, 3);
  TH1D* momentReweightedRMS = reweightMoment(momentMeasuredRMS, responseMatrix, 0);

  // Moments from unfolded Spectra <(pt)^n>(N_ch)
  TH1D* momentUnfolded1 = getMomentFromSpectra(multPtUnfolded, 1);
  TH1D* momentUnfolded2 = getMomentFromSpectra(multPtUnfolded, 2);
  TH1D* momentUnfolded3 = getMomentFromSpectra(multPtUnfolded, 3);
  TH1D* momentUnfoldedRMS = calcRMS(momentUnfolded2);



  ///----------------------------------------------------------------------------------------------------
  /// Unfolding for MC
  ///----------------------------------------------------------------------------------------------------

  // n_0(N_ch) : Initial multiplicity distribution (starting point for RooUnfold)
  TH1D* multDistGeneratedMC = responseMatrixOrig->ProjectionY("multDistGeneratedMC");
  multDistGeneratedMC->GetYaxis()->SetTitle("#it{n}(#it{N}_{ch})");

  // n(N_acc)_MC
  TH1D* multDistMeasuredMC = responseMatrixOrig->ProjectionX("multDistMeasuredMC");

  // n(N_ch)_MC : Unfolded Multiplicity Distribution from MC 'measurement' after nIter iterations
  cout  << endl << "-> Unfolding measured mutliplicity distribution (MC)...." << endl;
  TH1D* multDistUnfoldedMC = getUnfoldedMultDist(multDistMeasuredMC, responseMatrixOrig, nIter, useSmoothing);

  // P(N_ch|N_acc)_MC : Unfolding Matrix_MC based on the unfolded multiplicity distribution
  TH2D* unfoldingMatrixMC = getUnfoldingMatrix(responseMatrix, multDistUnfoldedMC);

  // MC Closure Test for Unfolding of Multiplicity Distribution
  TH1D* multDistMeasuredClosureTest = responseMatrixOrig->ProjectionX("multDistMeasuredClosureTest");
  TH1D* multDistGeneratedClosureTest = responseMatrixOrig->ProjectionY("multDistGeneratedClosureTest");
  TH1D* multDistUnfoldedClosureTest = getUnfoldedMultDist(multDistMeasuredClosureTest, responseMatrixOrig, nIter, useSmoothing);
  multDistUnfoldedClosureTest->SetName("multDistUnfoldedClosureTest");

  // MC Closure Test for Unfolding of Multiplicity Distribution with flat start
  TH1D* multDistInitialClosureTestFlat = responseMatrix->ProjectionY("multDistInitialClosureTestFlat");
  TH1D* multDistUnfoldedClosureTestFlat = getUnfoldedMultDist(multDistMeasuredClosureTest, responseMatrix, nIter, useSmoothing);
  multDistUnfoldedClosureTestFlat->SetName("multDistUnfoldedClosureTestFlat");

  // MC Closure Test for Unfolding of Integrated Spectra
  TH1D* multDistTracksMeasuredClosureTest = multPtMeasuredMC->ProjectionX("multDistTracksMeasuredClosureTest");
  TH1D* multDistTracksGeneratedClosureTest = multPtGeneratedMC->ProjectionX("multDistTracksGeneratedClosureTest");
  TH1D* multDistTracksUnfoldedClosureTest = getUnfoldedMultDist(multDistTracksMeasuredClosureTest, responseMatrixTracksOrig, nIter, useSmoothing);
  multDistTracksUnfoldedClosureTest->SetName("multDistTracksUnfoldedClosureTest");
  TH1D* multDistTracksInitialClosureTest = responseMatrixTracksOrig->ProjectionY("multDistTracksInitialClosureTest");

  TH1D* ptDistMeasuredMC = multPtMeasuredMC->ProjectionY("ptDistMeasuredMC");
  TH1D* ptDistGeneratedMC = multPtGeneratedMC->ProjectionY("ptDistGeneratedMC");


  // Unfolded pT Spectra (pT vs. N_ch)
  cout  << endl << "-> Unfolding measured pT spectra (MC)...." << endl;
  TH2D* multPtUnfoldedMC = unfoldPtSpectra(multPtMeasuredMC, responseMatrixTracksOrig, nIter, useSmoothing);
  multPtUnfoldedMC->SetName("multPtUnfoldedMC");

  // Moments from MC truth Spectra <(pT)^n>(N_ch)
  TH1D* momentGeneratedMC1 = getMomentFromSpectra(multPtGeneratedMC, 1);
  TH1D* momentGeneratedMC2 = getMomentFromSpectra(multPtGeneratedMC, 2);
  TH1D* momentGeneratedMC3 = getMomentFromSpectra(multPtGeneratedMC, 3);
  TH1D* momentGeneratedMCRMS = calcRMS(momentGeneratedMC2);

  // Measured Moments MC <(pT)^n>(N_acc)
  TH1D* momentMeasuredMC1 = getMomentFromSpectra(multPtMeasuredMC, 1);
  TH1D* momentMeasuredMC2 = getMomentFromSpectra(multPtMeasuredMC, 2);
  TH1D* momentMeasuredMC3 = getMomentFromSpectra(multPtMeasuredMC, 3);
  TH1D* momentMeasuredMCRMS = calcRMS(momentMeasuredMC2);

  // Reweighted Moments MC <(pT)^n>(N_ch) from <(pT)^n>(N_acc) using response matrix
  TH1D* momentReweightedMC1 = reweightMoment(momentMeasuredMC1, responseMatrix, 1);
  TH1D* momentReweightedMC2 = reweightMoment(momentMeasuredMC2, responseMatrix, 2);
  TH1D* momentReweightedMC3 = reweightMoment(momentMeasuredMC3, responseMatrix, 3);
  TH1D* momentReweightedMCRMS = reweightMoment(momentMeasuredMCRMS, responseMatrix, 0);

  // Moments from unfolded Spectra <(pT)^n>(N_ch)
  TH1D* momentUnfoldedMC1 = getMomentFromSpectra(multPtUnfoldedMC, 1);
  TH1D* momentUnfoldedMC2 = getMomentFromSpectra(multPtUnfoldedMC, 2);
  TH1D* momentUnfoldedMC3 = getMomentFromSpectra(multPtUnfoldedMC, 3);
  TH1D* momentUnfoldedMCRMS = calcRMS(momentUnfoldedMC2);

  if(includeSimulation){

    TH2D* multPtSimulatedMC = (TH2D*) inputFile->FindObjectAny("multPtSimulatedMC");
    if(!multPtSimulatedMC) cout << "ERROR: Could not find multPtSimulatedMC" << endl;

    TH2D* multPtSimulatedMPI1MC = (TH2D*) inputFile->FindObjectAny("multPtSimulatedMPI1MC");
    if(!multPtSimulatedMPI1MC) cout << "ERROR: Could not find multPtSimulatedMPI1MC" << endl;

    TH2D* multPtSimulatedMPI2MC = (TH2D*) inputFile->FindObjectAny("multPtSimulatedMPI2MC");
    if(!multPtSimulatedMPI2MC) cout << "ERROR: Could not find multPtSimulatedMPI2MC" << endl;

    TH2D* multPtSimulatedMPI3MC = (TH2D*) inputFile->FindObjectAny("multPtSimulatedMPI3MC");
    if(!multPtSimulatedMPI3MC) cout << "ERROR: Could not find multPtSimulatedMPI3MC" << endl;

    TH1D* momentSimulatedMPI1MC1 = getMomentFromSpectra(multPtSimulatedMPI1MC, 1);
    momentSimulatedMPI1MC1->SetName("momentSimulatedMPI1MC1");
    TH1D* momentSimulatedMPI2MC1 = getMomentFromSpectra(multPtSimulatedMPI2MC, 1);
    momentSimulatedMPI2MC1->SetName("momentSimulatedMPI2MC1");
    TH1D* momentSimulatedMPI3MC1 = getMomentFromSpectra(multPtSimulatedMPI3MC, 1);
    momentSimulatedMPI3MC1->SetName("momentSimulatedMPI3MC1");

    // Moments from MC truth Spectra <(pT)^n>(N_ch)
    TH1D* momentSimulatedMC1 = getMomentFromSpectra(multPtSimulatedMC, 1);
    TH1D* momentSimulatedMC2 = getMomentFromSpectra(multPtSimulatedMC, 2);
    TH1D* momentSimulatedMC3 = getMomentFromSpectra(multPtSimulatedMC, 3);
    TH1D* momentSimulatedMCRMS = calcRMS(momentSimulatedMC2);
  }

  ///----------------------------------------------------------------------------------------------------
  /// Calculate Systematics
  ///----------------------------------------------------------------------------------------------------

  if(includeSystematics){
    TFile* systFile = TFile::Open(inputFileNameSyst.c_str(),"READ");
    if(!systFile) {cout << "Error: Could not open Systematics file!" << endl;}
    else{
      TH1D* momentReweighted1SystErr = getSystematix(momentReweighted1, systFile, momentReweightedMC1, momentGeneratedMC1);
      TH1D* momentReweighted2SystErr = getSystematix(momentReweighted2, systFile, momentReweightedMC2, momentGeneratedMC2);
      TH1D* momentReweighted3SystErr = getSystematix(momentReweighted3, systFile, momentReweightedMC3, momentGeneratedMC3);
      TH1D* momentReweightedRMSSystErr = getSystematix(momentReweightedRMS, systFile, momentReweightedMCRMS, momentGeneratedMCRMS);

      TH1D* momentUnfolded1SystErr = getSystematix(momentUnfolded1, systFile, momentUnfoldedMC1, momentGeneratedMC1);
      TH1D* momentUnfolded2SystErr = getSystematix(momentUnfolded2, systFile, momentUnfoldedMC2, momentGeneratedMC2);
      TH1D* momentUnfolded3SystErr = getSystematix(momentUnfolded3, systFile, momentUnfoldedMC3, momentGeneratedMC3);
      TH1D* momentUnfoldedRMSSystErr = getSystematix(momentUnfoldedRMS, systFile, momentUnfoldedMCRMS, momentGeneratedMCRMS);
    }
  }

  ///----------------------------------------------------------------------------------------------------------------------
  /// Output
  ///----------------------------------------------------------------------------------------------------------------------

//tempstuff
  TH1D* efficiencyPt = (TH1D*) inputFile->FindObjectAny("efficiencyPt");
  TH1D* secContPt = (TH1D*) inputFile->FindObjectAny("secContPt");


//tempstuf

  /// Container to store histograms
  TObjArray* outputHistos = new TObjArray(1);
  outputHistos->Add(efficiencyPt);
  outputHistos->Add(secContPt);

  if(includeSimulation){
    outputHistos->Add(momentSimulatedMC1);
    outputHistos->Add(momentSimulatedMC2);
    outputHistos->Add(momentSimulatedMC3);
    outputHistos->Add(momentSimulatedMCRMS);
    outputHistos->Add(momentSimulatedMPI1MC1);
    outputHistos->Add(momentSimulatedMPI2MC1);
    outputHistos->Add(momentSimulatedMPI3MC1);

  }

  outputHistos->Add(multPtMeasuredPtResUnfoldedMC);
  outputHistos->Add(ptResolution);
  outputHistos->Add(ptResponseMatrix);

  outputHistos->Add(responseMatrixOrig);
  outputHistos->Add(responseMatrix);
  outputHistos->Add(unfoldingMatrix);
  outputHistos->Add(responseMatrixTracksOrig);
  outputHistos->Add(responseMatrixTracks);

  outputHistos->Add(multDistMeasured);
  outputHistos->Add(multDistUnfolded);
  outputHistos->Add(multDistGeneratedMC);
  outputHistos->Add(multDistMeasuredMC);
  outputHistos->Add(multDistUnfoldedMC);

  outputHistos->Add(multPtMeasured);
  outputHistos->Add(multPtUnfolded);

  outputHistos->Add(multPtMeasuredMC);
  outputHistos->Add(multPtUnfoldedMC);
  outputHistos->Add(multPtGeneratedMC);

  outputHistos->Add(momentMeasuredMC1);
  outputHistos->Add(momentMeasuredMC2);
  outputHistos->Add(momentMeasuredMC3);
  outputHistos->Add(momentMeasuredMCRMS);
  outputHistos->Add(momentReweightedMC1);
  outputHistos->Add(momentReweightedMC2);
  outputHistos->Add(momentReweightedMC3);
  outputHistos->Add(momentReweightedMCRMS);
  outputHistos->Add(momentUnfoldedMC1);
  outputHistos->Add(momentUnfoldedMC2);
  outputHistos->Add(momentUnfoldedMC3);
  outputHistos->Add(momentUnfoldedMCRMS);
  outputHistos->Add(momentGeneratedMC1);
  outputHistos->Add(momentGeneratedMC2);
  outputHistos->Add(momentGeneratedMC3);
  outputHistos->Add(momentGeneratedMCRMS);

  outputHistos->Add(momentMeasured1);
  outputHistos->Add(momentMeasured2);
  outputHistos->Add(momentMeasured3);
  outputHistos->Add(momentMeasuredRMS);

  outputHistos->Add(momentReweighted1);
  outputHistos->Add(momentReweighted2);
  outputHistos->Add(momentReweighted3);
  outputHistos->Add(momentReweightedRMS);

  outputHistos->Add(momentUnfolded1);
  outputHistos->Add(momentUnfolded2);
  outputHistos->Add(momentUnfolded3);
  outputHistos->Add(momentUnfoldedRMS);

  outputHistos->Add(unfoldingMatrixMC);
  outputHistos->Add(multDistGeneratedClosureTest);
  outputHistos->Add(multDistMeasuredClosureTest);
  outputHistos->Add(multDistUnfoldedClosureTest);

  outputHistos->Add(multDistInitialClosureTestFlat);
  outputHistos->Add(multDistUnfoldedClosureTestFlat);


  outputHistos->Add(multDistTracksGeneratedClosureTest);
  outputHistos->Add(multDistTracksMeasuredClosureTest);
  outputHistos->Add(multDistTracksUnfoldedClosureTest);
  outputHistos->Add(multDistTracksInitialClosureTest);

  outputHistos->Add(ptDistMeasuredMC);
  outputHistos->Add(ptDistGeneratedMC);

//  outputHistos->Add(multPtMultGen);


  if(includeSystematics){
    if(systFile){
      outputHistos->Add(momentReweighted1SystErr);
      outputHistos->Add(momentReweighted2SystErr);
      outputHistos->Add(momentReweighted3SystErr);
      outputHistos->Add(momentReweightedRMSSystErr);

      outputHistos->Add(momentUnfolded1SystErr);
      outputHistos->Add(momentUnfolded2SystErr);
      outputHistos->Add(momentUnfolded3SystErr);
      outputHistos->Add(momentUnfoldedRMSSystErr);
    }
  }


  // Write output to file
  TFile* outputFile =  new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  outputHistos->Write();

//  TFile* iterationCheckFile = TFile::Open((outputPath + "iterationCheck.root").c_str(),"update");
//  iterationCheckFile->cd();
//  TH1D* iterCheckHist = momentUnfolded1->Clone(iterCheckHistName.c_str());
//  iterCheckHist->Write();

  // Print some useful information
  cout << endl << "-----------------------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------------------" << endl << endl;

  cout << "    Input File:                  " 	<< inputFileName.c_str() << endl;
  cout << "    Output File:                 " 	<< outputFileName.c_str() << endl;
  cout << "    Histograms written to file:  " 	<< outputHistos->GetEntries() << endl;
  cout << "    Iterations of Unfolding:     " 	<< nIter << endl << endl;

  cout << "    Data - Events measured:      " 	<< multDistMeasured->Integral() << endl;
  cout << "    Data - Events unfolded:      " 	<< multDistUnfolded->Integral() << endl;
  cout << "    MC   - Events generated:     " 	<< multDistGeneratedMC->Integral() << endl;
  cout << "    MC   - Events measured:      " 	<< multDistMeasuredMC->Integral() << endl;
  cout << "    MC   - Events unfolded:      " 	<< multDistUnfoldedMC->Integral() << endl << endl;

  cout << "    MC   - Tracks measured:      " 	<< multPtMeasuredMC->Integral() << endl;
  cout << "    MC   - Tracks generated:     " 	<< multPtGeneratedMC->Integral() << endl;
  cout << "    MC   - Tracks unfolded:      " 	<< multPtUnfoldedMC->Integral() << endl << endl;

  cout << endl << "-----------------------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------------------" << endl << endl;
  cout << "Running time: " << stopwatch->RealTime() / 60 << " minutes" << endl << endl;

//  TFile* respmatFile = TFile::Open("detectorResponse_pp5TeV.root","RECREATE");
//  respmatFile->cd();
//  responseMatrix->Write();

  TBrowser* outputInspector = new TBrowser();

}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TH2D* matrixMutliplication(TH2D* matrix1, TH2D* matrix2){

  TH2D* product = matrix1->Clone("Unity");
  product->Reset();

  for(Int_t k = 1; k <= matrix1->GetNbinsY(); k++){
    for(Int_t l = 1; l <= matrix1->GetNbinsX(); l++){

      Double_t newValue = 0;
      for(Int_t i = 1; i <= matrix1->GetNbinsX(); i++){
	      newValue += matrix1->GetBinContent(i, k) * matrix2->GetBinContent(l, i);
      }

      product->SetBinContent(l,k, newValue);
    }
  }
  return product;
}


TH1D* getSystematix(TH1D* histogram, TFile* systFile, TH1D* mcUnfRew, TH1D* mcGen, Bool_t linearizeClosureSyst = kTRUE){
  string origName = histogram->GetName();
  string name = origName + "SystErr";

  TH1D* histogramSys = histogram->Clone(name.c_str());

  TList* list = (TList*)systFile->Get(origName.c_str());
  TH1D* trackSyst = (TH1D*)list->FindObject("total");

  TH1D* closureTestSyst = mcUnfRew->Clone();
  closureTestSyst->Divide(mcGen);

  for(Int_t i = 2; i <= closureTestSyst->GetNbinsX(); i++){
    Double_t value = closureTestSyst->GetBinContent(i);
    if(value) value = TMath::Abs(value - 1);
    closureTestSyst->SetBinContent(i, value);
//    closureTestSyst->SetBinError(i, 0);
  }

  TF1* flinear = new TF1("flinear", "[0] + [1]*x", 0, 100);
  closureTestSyst->Fit(flinear, "Q0", "", 20., 50.);


  for(Int_t i = 2; i <= histogramSys->GetNbinsX(); i++){

    Double_t sysError = trackSyst->GetBinContent(i) * trackSyst->GetBinContent(i);
    Double_t closureSyst = closureTestSyst->GetBinContent(i);
//    cout << i << ": " << closureSyst*100 << " | ";
    if(linearizeClosureSyst && closureTestSyst->GetBinCenter(i) > 50) closureSyst = flinear->Eval(closureTestSyst->GetBinCenter(i));
//    cout << closureSyst*100 << endl;
    sysError = sysError + closureSyst * closureSyst;
    sysError = TMath::Sqrt(sysError);
    sysError = sysError * histogramSys->GetBinContent(i);
    if(histogramSys->GetBinContent(i)) histogramSys->SetBinError(i, sysError);
  }


  return histogramSys;
}




TH1D* divideByMult(TH1D* trackHist){

  for(Int_t i = trackHist->GetXaxis()->FindBin(1); i <= trackHist->GetNbinsX(); i++){

    Double_t nTracks = trackHist->GetBinContent(i);
    Double_t nTracksErr = trackHist->GetBinError(i);
    Double_t multEvent = trackHist->GetBinCenter(i);
    trackHist->SetBinContent(i, nTracks/multEvent);
    trackHist->SetBinError(i, nTracksErr/multEvent);

  }
  return trackHist;
}

/// ---------------------------------------------------------------------------
/// Functions to normalize Multiplicity distributions and mult dependandt spectra
/// ---------------------------------------------------------------------------

void normalizeBinWidth(TH1D* multDist){
  for(Int_t multBin = 1; multBin <=  multDist->GetNbinsX(); multBin++){
    Double_t binWidth =  multDist->GetXaxis()->GetBinWidth(multBin);
    Double_t binContent =  multDist->GetBinContent(multBin);
    Double_t binError =  multDist->GetBinWidth(multBin);
    multDist->SetBinContent(multBin, binContent/binWidth);
    multDist->SetBinError(multBin, binError/binWidth);
  }
}

void normalizeBinWidth(TH2D* multPt){
  for(Int_t multBin = 1; multBin <=  multPt->GetNbinsX(); multBin++){
    Double_t binWidth =  multPt->GetXaxis()->GetBinWidth(multBin);
    for(Int_t ptBin = 1; ptBin <=  multPt->GetNbinsY(); ptBin++){
      Double_t binContent = multPt->GetBinContent(multBin, ptBin);
      Double_t binError = multPt->GetBinError(multBin, ptBin);
      multPt->SetBinContent(multBin, ptBin, binContent/binWidth);
      multPt->SetBinError(multBin, ptBin, binError/binWidth);
    }
  }
}

/// ---------------------------------------------------------------------------
/// Function to clean responseMatrix before normalizing
/// ---------------------------------------------------------------------------

void deleteMatrixZeroes(TH2D* responseMatrixOrig){
  for(Int_t i = 1; i <=  responseMatrixOrig->GetNbinsX(); i++){
    responseMatrixOrig->SetBinContent(i, 1, 0);
    responseMatrixOrig->SetBinContent(1, i, 0);
    responseMatrixOrig->SetBinError(i, 1, 0);
    responseMatrixOrig->SetBinError(1, i, 0);
  }
}

void deleteMatrixZeroes(TH3D* responseMatrixOrig){
  for(Int_t pt = 1; pt <=  responseMatrixOrig->GetNbinsY(); pt++){
    for(Int_t i = 1; i <=  responseMatrixOrig->GetNbinsX(); i++){
      responseMatrixOrig->SetBinContent(i,pt, 1, 0);
      responseMatrixOrig->SetBinContent(1,pt, i, 0);
      responseMatrixOrig->SetBinError(i,pt, 1, 0);
      responseMatrixOrig->SetBinError(1,pt, i, 0);
    }
  }
}

/// ---------------------------------------------------------------------------
/// Function to clean multPt measured MC from Nacc=0 tracks and those in RespMat overflow
/// ---------------------------------------------------------------------------

void cleanMultPtMeasuredMC(TH2D* responseMatrixOrig, TH2D* multPtMeasuredMC){

  // first remove entries for Nacc = 0;
//  cout << "Removing tracks from multPtMeasuredMC for Nacc = 0..." << endl;
  for (Int_t j = 1; j <= multPtMeasuredMC->GetNbinsY(); j++){
    multPtMeasuredMC->SetBinContent(1,j,0);
    multPtMeasuredMC->SetBinError(1,j,0);
  }


  Int_t overflowBin = 1 + responseMatrixOrig->GetNbinsY();

  for (Int_t i = 1; i <= multPtMeasuredMC->GetNbinsX(); i++){

    if(responseMatrixOrig->GetBinContent(i, overflowBin)){
//      cout << "Nacc = " << i-1 << " has overflow! Removing tracks from multPtMeasuredMC..." << endl;

      for (Int_t j = 1; j <= multPtMeasuredMC->GetNbinsY(); j++){
      	multPtMeasuredMC->SetBinContent(i,j,0);
      	multPtMeasuredMC->SetBinError(i,j,0);
      }
    }
  }

  for (Int_t j = 1; j <= multPtMeasuredMC->GetNbinsX(); j++){
    Int_t pt150 = multPtMeasuredMC->GetYaxis()->FindBin(0.15-0.01);
//    if(multPtMeasuredMC->GetBinContent(j,pt150))  cout << "Removing tracks from multPtMeasuredMC for pt = 150MeV, Nch = " << j-1 << endl;
    multPtMeasuredMC->SetBinContent(j,pt150,0);
    multPtMeasuredMC->SetBinError(j,pt150,0);
  }


}


void cleanMultPtGeneratedMC(TH2D* multPtGeneratedMC){
  for (Int_t j = 1; j <= multPtGeneratedMC->GetNbinsX(); j++){
    Int_t pt150 = multPtGeneratedMC->GetYaxis()->FindBin(0.15-0.01);
//    if(multPtGeneratedMC->GetBinContent(j,pt150))  cout << "Removing tracks from multPtGeneratedMC for pt = 150MeV, Nch = " << j-1 << endl;
    multPtGeneratedMC->SetBinContent(j,pt150,0);
    multPtGeneratedMC->SetBinError(j,pt150,0);
  }
}



/// ---------------------------------------------------------------------------
/// Function to obtain track unfolding Matrix from Event unfolding Matrix
/// ---------------------------------------------------------------------------

TH2D* getTrackUnfoldingMatrix(TH2D* unfoldingMatrix){

  TH2D* unfoldingMatrixTracks = unfoldingMatrix->Clone("unfoldingMatrixTracks");
  unfoldingMatrixTracks->Reset();

  TString originalName = unfoldingMatrix->GetName();
  if (originalName.Contains("MC")) unfoldingMatrixTracks->SetName("unfoldingMatrixTracksMC");


  for(Int_t Nacc = 1; Nacc <= unfoldingMatrixTracks->GetNbinsX(); Nacc++){
  /* old stuff
    Double_t meanNch = 0;

    for(Int_t Nch = 1; Nch <= unfoldingMatrixTracks->GetNbinsY(); Nch++){
      meanNch += (Nch-1) * unfoldingMatrix->GetBinContent(Nacc, Nch);
    }

*/
    Double_t denominator = (Nacc-1);
    if(denominator){
      for(Int_t Nch = 1; Nch <= unfoldingMatrixTracks->GetNbinsY(); Nch++){
	unfoldingMatrixTracks->SetBinContent(Nacc, Nch, (Nch-1)/denominator * unfoldingMatrix->GetBinContent(Nacc, Nch));
      }
    }
  }
  return unfoldingMatrixTracks;
}








TH2D* getTrackResponseMatrix(TH2D* responseMatrixOrig){

  TH2D* responseMatrixTracksOrig = responseMatrixOrig->Clone("responseMatrixTracksOrig");
  responseMatrixTracksOrig->Reset();

  for(Int_t NaccBin = 1; NaccBin <= responseMatrixTracksOrig->GetNbinsX(); NaccBin++){
    Double_t Nacc = responseMatrixTracksOrig->GetXaxis()->GetBinCenter(NaccBin);
    for(Int_t NchBin = 1; NchBin <= responseMatrixTracksOrig->GetNbinsY(); NchBin++){
      Double_t content = responseMatrixOrig->GetBinContent(NaccBin, NchBin);
      Double_t error = responseMatrixOrig->GetBinError(NaccBin, NchBin);
      responseMatrixTracksOrig->SetBinContent(NaccBin, NchBin, Nacc * content);
      responseMatrixTracksOrig->SetBinError(NaccBin, NchBin, Nacc * error);
    }
  }
  return responseMatrixTracksOrig;
}







/// ---------------------------------------------------------------------------
/// Function to scale unfolded pt spectra according to the unfolded Multiplicity
/// Distribution
/// ---------------------------------------------------------------------------

TH2D* normalizeSpectra(TH2D* multPtUnfoldedYield, TH1D* multDistUnfolded){

  string originalName = multPtUnfoldedYield->GetName();
  originalName += "Norm";

  TH2D* normSpectra = multPtUnfoldedYield->Clone(originalName.c_str());
  normSpectra->Reset();


//  Double_t nEvents = multDistUnfolded->Integral();

  for(Int_t multBin = 1; multBin <= multPtUnfoldedYield->GetNbinsX(); multBin++){
    if(!multDistUnfolded->GetBinContent(multBin)) continue;
    Double_t scalingFactor = 1 / multDistUnfolded->GetBinContent(multBin);
    for(Int_t ptBin = 1; ptBin <= multPtUnfoldedYield->GetNbinsY(); ptBin++){
      Double_t currentContent = multPtUnfoldedYield->GetBinContent(multBin, ptBin);
      normSpectra->SetBinContent(multBin, ptBin, currentContent * scalingFactor);
    }
  }

  return normSpectra;
}



TH2D* undoInvYield(TH2D* multPt, TH1D* multDist){

  string originalName = multPt->GetName();
  originalName += "Yield";
  TH2D* yieldSpectra = multPt->Clone(originalName.c_str());
  yieldSpectra->Reset();

  Double_t nEvents = multDist->Integral();
  Double_t scalingFactor = nEvents * 2* TMath::Pi();

  for(Int_t multBin = 1; multBin <= yieldSpectra->GetNbinsX(); multBin++){
    for(Int_t ptBin = 1; ptBin <= yieldSpectra->GetNbinsY(); ptBin++){
      Double_t currentContent = multPt->GetBinContent(multBin, ptBin);
      Double_t pt = yieldSpectra->GetYaxis()->GetBinCenter(ptBin);
      yieldSpectra->SetBinContent(multBin, ptBin, currentContent * scalingFactor * pt);
    }
  }
//  yieldSpectra->GetXaxis()->SetRangeUser(0.5,82.5);
  yieldSpectra->GetYaxis()->SetRangeUser(0.15,10);
  yieldSpectra->GetZaxis()->SetRangeUser(1e-8,1e2);
//  yieldSpectra->GetZaxis()->SetTitle("");

  return yieldSpectra;
}







TH2D* dividePt(TH2D* multPt){

  string originalName = multPt->GetName();

  TH2D* yieldSpectra = multPt->Clone("dummyName");
  yieldSpectra->Reset();

  for(Int_t multBin = 1; multBin <= yieldSpectra->GetNbinsX(); multBin++){
    for(Int_t ptBin = 1; ptBin <= yieldSpectra->GetNbinsY(); ptBin++){
      Double_t currentContent = multPt->GetBinContent(multBin, ptBin);
      Double_t currentError = multPt->GetBinError(multBin, ptBin);
      Double_t pt = yieldSpectra->GetYaxis()->GetBinCenter(ptBin);
      Double_t width_pt = yieldSpectra->GetYaxis()->GetBinWidth(ptBin);

      yieldSpectra->SetBinContent(multBin, ptBin, currentContent / (width_pt * pt));
      yieldSpectra->SetBinError(multBin, ptBin, currentError / (width_pt * pt));
    }
  }

  yieldSpectra->SetName(originalName.c_str());
  return yieldSpectra;
}


TH3D* dividePt(TH3D* multPt){

  string originalName = multPt->GetName();

  TH3D* yieldSpectra = multPt->Clone("dummyName");
  yieldSpectra->Reset();

  for(Int_t multBin = 1; multBin <= yieldSpectra->GetNbinsX(); multBin++){
    for(Int_t Nch = 1; Nch <= yieldSpectra->GetNbinsZ(); Nch++){
      for(Int_t ptBin = 1; ptBin <= yieldSpectra->GetNbinsY(); ptBin++){
	Double_t currentContent = multPt->GetBinContent(multBin, ptBin, Nch);
	Double_t currentError = multPt->GetBinError(multBin, ptBin, Nch);
	Double_t pt = yieldSpectra->GetYaxis()->GetBinCenter(ptBin);
	Double_t width_pt = yieldSpectra->GetYaxis()->GetBinWidth(ptBin);

	yieldSpectra->SetBinContent(multBin, ptBin, Nch, currentContent / (width_pt * pt));
	yieldSpectra->SetBinError(multBin, ptBin, Nch, currentError / (width_pt * pt));
      }
    }
  }

  yieldSpectra->SetName(originalName.c_str());
  return yieldSpectra;
}










/// ---------------------------------------------------------------------------
/// Function to obtain moments from momentum distributions
/// Assumes Pt distribution to be the invariant yield.
/// ---------------------------------------------------------------------------
TH1D* getMomentFromSpectra(TH2D* multPt, Int_t moment){

  TH1D* multMoment = (TH1D*) multPt->ProjectionX("dummyName");
  multMoment->Reset();

  TString originalName = multPt->GetName();
  if (originalName.Contains("Generated")) {multMoment->SetName(getSpectrumName(moment, "GeneratedMC").c_str());}
  else if (originalName.Contains("Simulated")) {multMoment->SetName(getSpectrumName(moment, "SimulatedMC").c_str());}
  else{
    if (originalName.Contains("Unfolded")){
      if (originalName.Contains("MC")) multMoment->SetName(getSpectrumName(moment, "UnfoldedMC").c_str());
      else multMoment->SetName(getSpectrumName(moment, "Unfolded").c_str());
    }else{
      if (originalName.Contains("MC")) multMoment->SetName(getSpectrumName(moment, "MeasuredMC").c_str());
      else multMoment->SetName(getSpectrumName(moment, "Measured").c_str());
    }
  }

  multMoment->GetYaxis()->SetTitle(getAxisTitle(moment).c_str());

  for (Int_t multBin = multPt->GetXaxis()->FindBin(1); multBin <= multPt->GetNbinsX(); multBin++){

    Double_t numerator = 0;
    Double_t denominator = 0;
    Double_t ptBinCenter = 0;
    Double_t ptBinWidth = 0;
    Double_t yield = 0;
    Double_t yieldError = 0;
    Double_t weightedQuantity = 0;

    // calculate moment (startbin at Pt > 150MeV)
    for (Int_t ptBin = multPt->GetYaxis()->FindBin(0.1501); ptBin <= multPt->GetNbinsY(); ptBin++){

      ptBinCenter = multPt->GetYaxis()->GetBinCenter(ptBin);
      ptBinWidth = multPt->GetYaxis()->GetBinWidth(ptBin);
      yield = multPt->GetBinContent(multBin, ptBin) * ptBinCenter ;

      weightedQuantity = 1.;
      for(Int_t i = 0; i < moment; i++) weightedQuantity *= ptBinCenter;

      numerator += yield * weightedQuantity * ptBinWidth;
      denominator += yield * ptBinWidth;
    }
    Double_t mean =  numerator / denominator;

    // calculate errors
    Double_t errSum = 0;
    for (Int_t ptBin = multPt->GetYaxis()->FindBin(0.1501); ptBin <= multPt->GetNbinsY(); ptBin++){

      ptBinCenter = multPt->GetYaxis()->GetBinCenter(ptBin);
      ptBinWidth = multPt->GetYaxis()->GetBinWidth(ptBin);
      yieldError = multPt->GetBinError(multBin, ptBin) * ptBinCenter;

      weightedQuantity = 1.;
      for(Int_t i = 0; i < moment; i++) weightedQuantity *= ptBinCenter;

      errSum += TMath::Power((weightedQuantity - mean) * yieldError * ptBinWidth, 2);
    }

    // fill results in histogram
    if(denominator){
      multMoment->SetBinContent(multBin, mean);
      Double_t error = 1/denominator * TMath::Sqrt(errSum);
      multMoment->SetBinError(multBin, error);
    }
  }
  return multMoment;
}


/// ---------------------------------------------------------------------------
/// Function to obtain n(N_ch) and the
/// response Matrix P(N_acc|N_ch).
/// ---------------------------------------------------------------------------
TH1D* getUnfoldedMultDist(TH1D* multDistMeasured, TH2D* responseMatrix, Int_t nIter, Bool_t useSmoothing = kFALSE, Int_t verbose = 0){

  TH1D* meas = responseMatrix->ProjectionX("projX"); 	// includes measurements with no corresponding truth (fakes/background)
  TH1D* truth = responseMatrix->ProjectionY("projY"); 	// includes unmeasured events (inefficiency)

  RooUnfoldResponse response(meas, truth, responseMatrix);
  RooUnfoldBayes unfoldObject(&response, multDistMeasured, nIter, useSmoothing);
  unfoldObject.SetVerbose(verbose);


  TH1D* multDistUnfolded = (TH1D*) unfoldObject.Hreco(RooUnfold::kCovariance);
  multDistUnfolded->SetTitle("");
  multDistUnfolded->GetYaxis()->SetTitle("#it{n}(#it{N}_{ch})");


  TString originalName = multDistMeasured->GetName();
  if (originalName.Contains("MC"))   multDistUnfolded->SetName("multDistUnfoldedMC");
  else   multDistUnfolded->SetName("multDistUnfolded");

  delete meas;
  delete truth;

  return multDistUnfolded;
}



/// ------------------------------------------------------------------------------------------------------------------------------
/// Function to calculate unfolding matrix P(N_ch|N_acc) using response Matrix P(N_acc|N_ch) and n(N_ch).
/// ------------------------------------------------------------------------------------------------------------------------------
TH2D* getUnfoldingMatrix(TH2D* responseMatrix, TH1D* multDistUnfolded){

  TH2D* unfoldingMatrix = responseMatrix->Clone("unfoldingMatrix");
  unfoldingMatrix->Reset();

  TString originalName = multDistUnfolded->GetName();
  if (originalName.Contains("MC")) unfoldingMatrix->SetName("unfoldingMatrixMC");

  // normalize response matrix to ensure it contains the reconstruction prbabilities and not measured abundances
  responseMatrix = normalize(responseMatrix);

  for (Int_t Nacc = 1 ; Nacc <= unfoldingMatrix->GetNbinsX() ; Nacc++){

    Double_t denominator = 0;
    for (Int_t Nch = 1 ; Nch <= unfoldingMatrix->GetNbinsY() ; Nch++){
      denominator += (responseMatrix->GetBinContent(Nacc, Nch)) * (multDistUnfolded->GetBinContent(Nch));
    }

    if (denominator){
      for (Int_t Nch = 1 ; Nch <= unfoldingMatrix->GetNbinsY() ; Nch++){
	Double_t value = (responseMatrix->GetBinContent(Nacc, Nch) * multDistUnfolded->GetBinContent(Nch))/denominator;
	if (value) unfoldingMatrix->SetBinContent(Nacc, Nch, value);
      }
    }
  }

  delete responseMatrix;
  return unfoldingMatrix;
}

/// --------------------------------------------------------------------------------------
/// Helper-Function to unfold pt(Nacc)-Spectra using the unfolding Matrix  P(N_ch|N_acc).
/// --------------------------------------------------------------------------------------

TH1D* unfoldSpectrum(TH1D* spectrum, TH2D* unfoldingMatrix){

  TH1D* unfoldedSpectrum = spectrum->Clone("dummyName");
  unfoldedSpectrum->Reset();

  for (Int_t Nch = 1; Nch <= unfoldedSpectrum->GetNbinsX(); Nch++){

    Double_t value = 0;
    Double_t errSum = 0;

    for (Int_t Nacc = 1; Nacc <= spectrum->GetNbinsX(); Nacc++){
      value += (unfoldingMatrix->GetBinContent(Nacc, Nch))*(spectrum->GetBinContent(Nacc));
      errSum += TMath::Power(unfoldingMatrix->GetBinContent(Nacc, Nch) * spectrum->GetBinError(Nacc), 2);
    }

    unfoldedSpectrum->SetBinContent(Nch, value);
    unfoldedSpectrum->SetBinError(Nch, TMath::Sqrt(errSum));
  }

  return unfoldedSpectrum;
}


/// --------------------------------------------------------------------------------------
/// Function to reweight measured moments using the response Matrix  P(N_acc|N_ch).
/// --------------------------------------------------------------------------------------
TH1D* reweightMoment(TH1D* spectrum, TH2D* responseMatrix, Int_t moment){

  TH1D* reweightedSpectrum = spectrum->Clone("dummyName");
  reweightedSpectrum->Reset();

  TString originalName = spectrum->GetName();
  if (originalName.Contains("MC")) reweightedSpectrum->SetName(getSpectrumName(moment, "ReweightedMC").c_str());
  else reweightedSpectrum->SetName(getSpectrumName(moment, "Reweighted").c_str());

  reweightedSpectrum->GetXaxis()->SetTitle("#it{N}_{ch}");
  reweightedSpectrum->GetYaxis()->SetTitle(getAxisTitle(moment).c_str());

  // only count N_acc up to maximum of measurement (otherwise the normalization sum(P(N_ch|N_acc)) falsifies the reweighted value)
  Int_t highestMeasuredMultBin = spectrum->GetNbinsX();
  while ((highestMeasuredMultBin > 0) && spectrum->GetBinContent(highestMeasuredMultBin) == 0) highestMeasuredMultBin--;
  if(highestMeasuredMultBin == 0) return reweightedSpectrum;


  for (Int_t Nch = spectrum->GetXaxis()->FindBin(1); Nch <= reweightedSpectrum->GetNbinsX(); Nch++){

    Double_t value = 0;
    Double_t errSum = 0;


    for (Int_t Nacc = spectrum->GetXaxis()->FindBin(1); Nacc <= highestMeasuredMultBin; Nacc++){
      value += (responseMatrix->GetBinContent(Nacc, Nch))*(spectrum->GetBinContent(Nacc));
      errSum += TMath::Power(responseMatrix->GetBinContent(Nacc, Nch) * spectrum->GetBinError(Nacc), 2);
    }

    // normalization (only nessecary because of tempstuff)
//    Double_t norm = 0;
//    for (Int_t Nacc = spectrum->GetXaxis()->FindBin(1); Nacc <= highestMeasuredMultBin; Nacc++){
//      norm += (responseMatrix->GetBinContent(Nacc, Nch));
//    }
//    if(norm) value = value/norm;

    reweightedSpectrum->SetBinContent(Nch, value);
    reweightedSpectrum->SetBinError(Nch, TMath::Sqrt(errSum));
  }
//TODO nachrechnen + was ist mit fehler von respmat???
  return reweightedSpectrum;
}


/// ---------------------------------------------------------------------------
/// Function to calculate root mean squared from second moment of distribution
/// TODO with correct propagation of uncertainties.
/// ---------------------------------------------------------------------------
TH1D* calcRMS(TH1D* meanSquared){

  TH1D* rms = meanSquared->Clone("dummyName");
  rms->Reset();

  TString originalName = meanSquared->GetName();
  if (originalName.Contains("Measured")) rms->SetName(getSpectrumName(0, "Measured").c_str());
  if (originalName.Contains("MeasuredMC")) rms->SetName(getSpectrumName(0, "MeasuredMC").c_str());
  if (originalName.Contains("Unfolded")) rms->SetName(getSpectrumName(0, "Unfolded").c_str());
  if (originalName.Contains("UnfoldedMC")) rms->SetName(getSpectrumName(0, "UnfoldedMC").c_str());
  if (originalName.Contains("GeneratedMC")) rms->SetName(getSpectrumName(0, "GeneratedMC").c_str());
  if (originalName.Contains("SimulatedMC")) rms->SetName(getSpectrumName(0, "SimulatedMC").c_str());


  rms->GetYaxis()->SetTitle("#sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})");

  for (Int_t i = 1 ; i <= rms->GetNbinsX() ; i++){

    Double_t value = TMath::Sqrt(meanSquared->GetBinContent(i));

    if(value){
      Double_t error = 0.5 / value * meanSquared->GetBinError(i);
      rms->SetBinContent(i, value);
      rms->SetBinError(i, error);
    }
  }
  return rms;
}


/// ---------------------------------------------------------------------------
/// Function to normalize response Matrix
/// ---------------------------------------------------------------------------
TH2D* normalize(TH2D* matrix){

  TH2D* normalizedMatrix = matrix->Clone("responseMatrix");
  normalizedMatrix->Reset();

  for (Int_t y = 2 ; y <= matrix->GetNbinsY() ; y++){

    Double_t yIntegral = 0;
    for (Int_t x = 2 ; x <= matrix->GetNbinsX() ; x++) yIntegral += matrix->GetBinContent(x, y);

    for (Int_t x = 2 ; x <= matrix->GetNbinsX() ; x++){

      Double_t value = matrix->GetBinContent(x, y);
      if (value) normalizedMatrix->SetBinContent(x, y, value/yIntegral);

    }
  }
  return normalizedMatrix;
}


/// ---------------------------------------------------------------------------
/// Function to calculate covariance Matrix
/// ---------------------------------------------------------------------------

TMatrixD* getCovariance(const TMatrixD& jacobi, const TMatrixD& inputCov)
{
  // outCov = J * Cov * J^T.
  TMatrixD temp (inputCov, TMatrixD::kMultTranspose, jacobi);
  TMatrixD* outputCov = new TMatrixD(jacobi.GetNrows(), jacobi.GetNrows());
  outputCov->Mult(jacobi, temp);
  return outputCov;
}


/// ---------------------------------------------------------------------------
/// Helper-Functions to get the naming of spectra and axes right
///
/// ---------------------------------------------------------------------------
string getSpectrumName(Int_t moment, char* type){

  char name[30] = "";

  sprintf(name, "moment%s%d", type, moment);
  if(moment == 0) sprintf(name, "moment%sRMS", type);

  return string(name);
}

string getAxisTitle(Int_t moment){

  char axisTitle[50] = "#LT#it{p}_{T}#GT (GeV/#it{c})";
  if (moment > 1) sprintf(axisTitle, "#LT#it{p}^{%d}_{T}#GT (GeV/#it{c})^{%d}", moment, moment);
  if (moment == 0) sprintf(axisTitle, "#sqrt{#LT#it{p}^{2}_{T}#GT} (GeV/#it{c})"); // for rms
  if (moment == -1) sprintf(axisTitle, "x (GeV/#it{c})"); // for mean and rms
  return string(axisTitle);
}





/// ---------------------------------------------------------------------------
/// Function to obtain event composition Matrix from Event unfolding Matrix
/// ---------------------------------------------------------------------------
TH2D* getCompositionMatrix(TH2D* unfoldingMatrix, TH1D* multDistMeasured){

  TH2D* compositionMatrix = unfoldingMatrix->Clone("compositionMatrix");
  compositionMatrix->Reset();


  for(Int_t Nch = 1; Nch <= compositionMatrix->GetNbinsY(); Nch++){

    Double_t nNch = 0;

    for(Int_t Nacc = 1; Nacc <= compositionMatrix->GetNbinsX(); Nacc++){
      nNch += unfoldingMatrix->GetBinContent(Nacc, Nch) * multDistMeasured->GetBinContent(Nacc);
    }

    if(nNch){
      for(Int_t Nacc = 1; Nacc <= compositionMatrix->GetNbinsX(); Nacc++){
	compositionMatrix->SetBinContent(Nacc, Nch, unfoldingMatrix->GetBinContent(Nacc, Nch) * multDistMeasured->GetBinContent(Nacc));// / nNch);
      }
    }
  }
  return compositionMatrix;

}



/// ---------------------------------------------------------------------------
/// Quick calculation of mean pt
/// ---------------------------------------------------------------------------

Double_t getMeanFromHist(TH1D* hist){

  Double_t numerator = 0;
  Double_t denominator = 0;

  for (Int_t ptBin = hist->GetXaxis()->FindBin(0.1501); ptBin <= hist->GetNbinsX() ; ptBin++){

    Double_t ptBinCenter = hist->GetBinCenter(ptBin);
    Double_t ptBinWidth = hist->GetBinWidth(ptBin);
    Double_t yield = hist->GetBinContent(ptBin) * ptBinCenter;

    numerator += yield * ptBinCenter * ptBinWidth;
    denominator += yield * ptBinWidth;
  }

  Double_t mean =  0;
  if(denominator) mean = numerator / denominator;

  return mean;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Artefacts
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// ---------------------------------------------------------------------------
/// Function to get Unfolding Matrix P(N_ch|N_acc) from RooUnfold
/// for crosschecking with the one calculated by hand.
/// ---------------------------------------------------------------------------
TH2D* getUnfoldingMatrixRooUnfold(TH1D* multDistMeasured, TH2D* responseMatrix, Int_t nIter, Bool_t useSmoothing = kFALSE){

  RooUnfoldResponse response(0, 0, responseMatrix);
  RooUnfoldBayes unfold(&response, multDistMeasured, nIter, useSmoothing);
  unfold.SetVerbose(0);

  TH1D* multDistUnfolded = (TH1D*) unfold.Hreco(RooUnfold::kCovariance);
  TMatrixD& M = unfold.UnfoldingMatrix();

  TH2D* unfoldingMatrix = responseMatrix->Clone("dummyName");
  unfoldingMatrix->Reset();

  // put values in the right bins
  for(Int_t x = 1; x <= unfoldingMatrix->GetNbinsX(); x++){
    for(Int_t y = 1; y <= unfoldingMatrix->GetNbinsY(); y++){

      Double_t value = M(y-1,x-1);
      if(value) unfoldingMatrix->SetBinContent(x,y, value);
//      cout << value;

    }
  }

  return unfoldingMatrix;
}


/// raw estimate of systematix am ende sind nullen im spektrum???
TH1D* getSystematix (TH2D* multPt, TH1D* meanPt, Double_t sysUnc){

  TH1D* meanPtSys = meanPt->Clone("Sys");

  TH2D* multPtUp = multPt->Clone("Up");
  TH2D* multPtDown = multPt->Clone("Down");

  Double_t factorUp = 1 + sysUnc;
  Double_t factorDown = 1 - sysUnc;

  for(Int_t multBin = 1; multBin <= multPt->GetNbinsX(); multBin++){

    Double_t mean = meanPt->GetBinContent(multBin);
    if(!mean) continue;

    for(Int_t ptBin = 1; ptBin < multPt->GetYaxis()->FindBin(mean); ptBin++){

      Double_t content = multPt->GetBinContent(multBin, ptBin);
      if(!content) continue;

      multPtUp->SetBinContent(multBin, ptBin, factorUp * content);
      multPtDown->SetBinContent(multBin, ptBin, factorDown * content);

    }
    for(Int_t ptBin = multPt->GetYaxis()->FindBin(mean) + 1; ptBin <= multPt->GetNbinsY(); ptBin++){

      Double_t content = multPt->GetBinContent(multBin, ptBin);
      if(!content) continue;

      multPtUp->SetBinContent(multBin, ptBin, factorDown * content);
      multPtDown->SetBinContent(multBin, ptBin, factorUp * content);

    }
  }

  TH1D* meanPtUp = getMomentFromSpectra(multPtUp, 1);
  TH1D* meanPtDown = getMomentFromSpectra(multPtDown, 1);

  for(Int_t multBin = 1; multBin <= meanPt->GetNbinsX(); multBin++){

//    cout << multBin << ": (" << meanPtUp->GetBinContent(multBin) <<" | " << meanPt->GetBinContent(multBin) << " | " << meanPtDown->GetBinContent(multBin) << ")" << endl;

    Double_t mean = meanPt->GetBinContent(multBin);
    Double_t diffUp = TMath::Abs(meanPtUp->GetBinContent(multBin) - mean);
    Double_t diffDown = TMath::Abs(meanPtDown->GetBinContent(multBin) - mean);

//    if(mean) cout << "MultBin " << multBin << ": " << setprecision(1) << "(" << 100* diffUp/mean << "% | " << 100 * diffDown/mean << "%)" << endl;

    if(diffUp > diffDown){
      meanPtSys->SetBinError(multBin, diffUp);
    }else{
      meanPtSys->SetBinError(multBin, diffDown);
    }


  }

  return meanPtSys;

}



/// ---------------------------------------------------------------------------
/// Function to normalize frequency distribution
/// ---------------------------------------------------------------------------
TH1D* normalize(TH1D* hist){

    TH1D* normalizedHist = hist->Clone("Normalized Hist");
    normalizedHist->Reset();
    normalizedHist->SetTitle("Normalized Distribution");

    Double_t integral = hist->Integral();

    for (Int_t x = 1; x <= hist->GetNbinsX(); x++){

      Double_t value = hist->GetBinContent(x);
      Double_t error = hist->GetBinError(x);
      if (value) {
	normalizedHist->SetBinContent(x, value/integral);
	normalizedHist->SetBinError(x, error/integral);
      }
    }

  return normalizedHist;
}




/// -------------------------------------------------------------------------------
/// Function for transposing the response matrix (assumes same binning in x and y)
/// -------------------------------------------------------------------------------
TH2D* transpose(TH2D* matrix){

  TH2D* transposedMatrix = matrix->Clone("Transposed Matrix");
  transposedMatrix->Reset();
  transposedMatrix->SetTitle("Transposed Matrix");

  for (Int_t x = 1 ; x <= matrix->GetNbinsX() ; x++){
    for (Int_t y = 1 ; y <= matrix->GetNbinsY() ; y++){

      Double_t value = matrix->GetBinContent(x, y);
      Double_t error = matrix->GetBinError(x, y);
      if (value){
	transposedMatrix->SetBinContent(y, x, value);
	transposedMatrix->SetBinError(y, x, error);
      }
    }
  }
  const char* xTitle = matrix->GetXaxis()->GetTitle();
  transposedMatrix->GetYaxis()->SetTitle(xTitle);

  const char* yTitle = matrix->GetYaxis()->GetTitle();
  transposedMatrix->GetXaxis()->SetTitle(yTitle);

  Double_t entries = matrix->GetEntries();
  transposedMatrix->SetEntries(entries);

  return transposedMatrix;
}



/// ---------------------------------------------------------------------------
/// Function to calculate chi^2 of Multiplicity between Iterations
///
/// ---------------------------------------------------------------------------
Double_t getChi2(TH1D* multDist1, TH1D* multDist2){

  Int_t nBins = multDist1->GetNbinsX();
  if(nBins != multDist2->GetNbinsX()){cout << "Dimensions do not match!" << endl; return -1;}

  Double_t  value1 = 0.;
  Double_t  value2 = 0.;

  Double_t integral = multDist1->Integral();

  Double_t chi2 = 0.;

  for(Int_t i = 0; i <= nBins; i++){

    value1 = multDist1->GetBinContent(i);
    value2 = multDist2->GetBinContent(i);

    chi2 += (value2 - value1) / value2;

  }

  return chi2;

}



void cutHisto(TH1D* inputHist, Int_t lastBin){

  for(Int_t i = lastBin; i <= inputHist->GetNbinsX(); i++){

    inputHist->SetBinContent(i, 0.);
    inputHist->SetBinError(i, 0.);

  }

}


/// ---------------------------------------------------------------------------
/// Function to convert Nch on xAxis to dNch/dEta
///
/// ---------------------------------------------------------------------------
TH1D* compressAxis(TH1D* inputSpectrum, Double_t factor){

  ostringstream stringStream;
  stringStream << inputSpectrum->GetName() << " dnchdeta";
  string name = stringStream.str();

  Double_t lowerBound = factor * inputSpectrum->GetBinLowEdge(1);
  Double_t lastBin = inputSpectrum->GetNbinsX();
  Double_t upperBound = factor * inputSpectrum->GetBinLowEdge(lastBin + inputSpectrum->GetBinWidth(lastBin));

  TH1D* outputSpectrum = new TH1D(name.c_str(), name.c_str(), inputSpectrum->GetNbinsX(), lowerBound, upperBound);
  outputSpectrum->GetXaxis()->SetTitle("d#it{N}_{ch} / d#it{#eta}");
  outputSpectrum->GetYaxis()->SetTitle(inputSpectrum->GetYaxis()->GetTitle());

  for(Int_t i = 1; i <= lastBin; i++){
    outputSpectrum->SetBinContent(i, inputSpectrum->GetBinContent(i));
    outputSpectrum->SetBinError(i, inputSpectrum->GetBinError(i));
  }


  return outputSpectrum;

}

/// ---------------------------------------------------------------------------
/// Function to cut Parts of a Matrix
///
/// ---------------------------------------------------------------------------
  TH2D* cutMatrix(TH2D* inputMatrix, Int_t xCut, Int_t yCut){

    TH2D* matrix = inputMatrix->Clone();

    for (Int_t x = 1 ; x <= matrix->GetNbinsX() ; x++){
      for (Int_t y = 1 ; y <= matrix->GetNbinsY() ; y++){
	if (x > xCut || y > yCut){
	    matrix->SetBinContent(x, y, 0);
	}
      }
    }
    return matrix;
  }


TH1D* addSysError(TH1D* inputHist, Double_t sysUnc){

  TH1D* inputHistSys = inputHist->Clone();

  for(Int_t i = 1; i <= inputHistSys->GetNbinsX(); i++){

    inputHistSys->SetBinError(i, sysUnc * inputHist->GetBinContent(i));

  }
  return inputHistSys;

}





/// bla meins startet mit -0.5 aber das aus paper mit 0.5 daher mit 2. bin anfangen...
TH1D* reduceRange(TH1D* largeHist, TH1D* smallHist){

  TH1D* shrinkedHist = smallHist->Clone("new");
  shrinkedHist->Reset();

  for(Int_t i = 1; i <= shrinkedHist->GetNbinsX(); i++){

    shrinkedHist->SetBinContent(i, largeHist->GetBinContent(i+1));
    shrinkedHist->SetBinError(i, largeHist->GetBinError(i+1));

  }
  return shrinkedHist;

}

/// --------------------------------------------------------------------------------------
/// Funktion to unfold pT-Spectra.
/// Assumes same binning in Nacc and Nch!
/// --------------------------------------------------------------------------------------

TH2D* unfoldPtSpectra(TH2D* multPtMeasured, TH2D* responseMatrixTracksOrig, Int_t nIter, Bool_t useSmoothing = kFALSE){

  TH2D* multPtUnfolded = multPtMeasured->Clone("dummyName");

  TString originalName = multPtMeasured->GetName();
  if (originalName.Contains("MC")) multPtUnfolded->SetName("multPtUnfoldedMC");
  else multPtUnfolded->SetName("multPtUnfolded");

  multPtUnfolded->Reset();
  multPtUnfolded->GetXaxis()->SetTitle("#it{N}_{ch}");

  TH1D* currentSpectrum = multPtMeasured->ProjectionX("dummyName");
  cout << "   |";

  for (Int_t ptBin = 1; ptBin <= multPtUnfolded->GetNbinsY(); ptBin++){

    currentSpectrum->Reset();
    currentSpectrum = multPtMeasured->ProjectionX("dummy", ptBin, ptBin);
    if(currentSpectrum->GetEntries() < 0.1) continue;
    cout << "=" << flush;
    TH1D* currentUnfoldedSpectrum = getUnfoldedMultDist(currentSpectrum, responseMatrixTracksOrig, nIter, useSmoothing);
    for (Int_t multBin = 1; multBin <= multPtUnfolded->GetNbinsX(); multBin++){
      multPtUnfolded->SetBinContent(multBin, ptBin, currentUnfoldedSpectrum->GetBinContent(multBin));
      multPtUnfolded->SetBinError(multBin, ptBin, currentUnfoldedSpectrum->GetBinError(multBin));
    }
    delete currentUnfoldedSpectrum;
  }
  cout << "|" << endl;
  delete currentSpectrum;
  return multPtUnfolded;
}



void applyEfficiencyCorrection(TH2D* multPtUncorrMC, TH2D* multPtGeneratedMC, TH2D* multPt){

  TH1D* corr = multPtGeneratedMC->ProjectionY("eff");
  corr->Divide(multPtUncorrMC->ProjectionY("eff2"));

  for(Int_t ptBin = 1; ptBin <= multPt->GetNbinsY(); ptBin++){

    Double_t corrFactor = corr->GetBinContent(ptBin);
    for(Int_t multBin = 1; multBin <= multPt->GetNbinsX(); multBin++){
      Double_t value = corrFactor * multPt->GetBinContent(multBin, ptBin);
      Double_t error = corrFactor * multPt->GetBinError(multBin, ptBin);
      multPt->SetBinContent(multBin, ptBin, value);
      multPt->SetBinError(multBin, ptBin, error);
    }
  }


}
TH1D* getEfficiencyCorrection(TH2D* multPtUncorrMC, TH2D* multPtGeneratedMC){

  TH1D* corr = multPtGeneratedMC->ProjectionY("corr");
  corr->Divide(multPtUncorrMC->ProjectionY("corr2"));
  return corr;

}

TH1D* getEfficiency(TH2D* multPtUncorrMC, TH2D* multPtGeneratedMC){

  TH1D* corr = multPtUncorrMC->ProjectionY("eff");
  corr->Divide(multPtGeneratedMC->ProjectionY("eff2"));
  return corr;
}
