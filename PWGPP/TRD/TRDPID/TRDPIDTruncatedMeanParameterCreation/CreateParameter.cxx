
#include <iostream>
#include <fstream>      
#include <ctime>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TList.h>
#include <TF1.h>
#include <TStyle.h>
#include "V0Tracks.cxx"
#include "constants.h"
#include "utils.h"
#include "signal.cxx"
#include "resolutionIncl.cxx"
#include "etaAsymmetryMap.cxx"


void CreateParameter() {
  TString InputFile = "/media/ikp/DATA/PID/TruncatedMean/constants.h";
  TString InputOutput= Form("output/%s_InputFile.txt", FileSpec.Data());
    
  copyFile(InputFile.Data(), InputOutput.Data());
  for (Int_t i=3; i<4; i++) { // loop only if different periods are fittet at the same time
    RunPeriod = RunPeriodArray[i];
    inputSpec = inputSpecArray[i];
    FileSpec = FileSpecArray[i];
    NclsCorrection=kFALSE; // no correction in first round
    fEtaCorrection=kFALSE;

    input = inputPath + inputSpec + inputFile;
    signal();
    resolutionIncl();
    etaAsymmetryMap();

    //cluster correction is applied in second round
    NclsCorrection = kTRUE;  
    signal();
    resolutionIncl(); 
    etaAsymmetryMap();

    //eta correction in third round
    fEtaCorrection = kTRUE;
    signal();
    resolutionIncl();
  }
}
