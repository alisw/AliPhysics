#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TList.h>
#include <TF1.h>
#include <TObject.h>
#include <TMath.h>

using std::map;
using std::string;
using std::vector;

const char   kLetter[2] = {'M','A'}; // M -> Matter, A -> Anti-matter
const string kNames[2] = {"deuterons","antideuterons"};

const string kBaseOutputDir = "../results/";
const string kBaseInputDir = "../";

const string kDataFilename = kBaseInputDir + "data/PbPb5TeV/treno20191123/LHC18qr.root";
const string kNormDataFilename = kBaseInputDir + "data/PbPb5TeV/treno20191109/LHC18qr_PIDqa.root";
const string kMCfilename = kBaseInputDir + "mc/PbPb5TeV/treno20191123/LHC19d2.root";

const string kFilterListNames = "mpuccio_deuterons_";
const string kNormalisationList = "mpuccio_deuterons_";

const string kEfficiencyOutput = kBaseOutputDir + "efficiency.root";
const string kSignalOutput = kBaseOutputDir + "signal.root";
const string kSecondariesOutput = kBaseOutputDir + "secondaries.root";
const string kSecondariesOutputRooFit = kBaseOutputDir + "RooSec.root";
const string kSecondariesTPCoutput = kBaseOutputDir + "secondaries_TPC.root";
const string kSecondariesTPCoutputRooFit = kBaseOutputDir + "RooSecTPC.root";
const string kMaterialOutput = kBaseOutputDir + "materialbudget.root";
const string kSpectraOutput = kBaseOutputDir + "spectra.root";
const string kFitSystematicsOutput = kBaseOutputDir + "fitsystematics.root";
const string kSystematicsOutput = kBaseOutputDir + "systematics.root";
const string kSystematicsOutputTPC = kBaseOutputDir + "systematics_TPC.root";
const string kFinalOutput = kBaseOutputDir + "final.root";

const bool   kPrintFigures{true};
const string kFiguresFolder = "../results/images/";
const string kMacrosFolder = "../results/images/macros/";

const int    kNPtBins = 30;
const int    kNCentralityBins = 11;
const float  kPtBins[kNPtBins] = {0.2f,0.3f,0.4f,0.5f,0.6f,0.7f,0.8f,0.9f,1.0f,1.1f,
    1.2f,1.4f,1.6f,1.8f,2.0f,2.2f,2.4f,2.6f,2.8f,3.0f,
    3.2f,3.4f,3.6f,3.8f,4.0f,4.2f,4.4f,5.0f,6.0f,8.0f};
const float  kCentralityBins[kNCentralityBins] = {0.f,5.f,10.f,20.f,30.f,40.f,50.f,60.f,70.f,80.f,90.f};

// Full cent bins
// const int    kCentLength = 10;
// const int    kCentBinsArray[kCentLength][2] = {{0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9},{9,10}};
// const float  kCentPtLimits[kCentLength] = {7,7,7,7,7,7,7,7,7,7};
// const float  kCentLabels[kCentLength][2] = {{0.,5.},{1.,5.},{5.,10.},{10.,20.},{20.,30.},{30.,40.},{40.,50.},{50.,60.},{60.,70.},{70.,80.}};
// const float  kPtRebin[kCentLength] = {2.6,2.6,2.2,2.2,2.2,2.,0.6,3.4};

// for specific bins
const int    kCentLength = 10;
const int    kCentBinsArray[kCentLength][2] = {{0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9},{9,10}};
const float  kCentPtLimits[kCentLength] = {7,7,7,7,7,7,7,7,7,7};
const float  kCentLabels[kCentLength][2] = {{0.,5.},{5.,10.},{10.,20.},{20.,30.},{30.,40.},{40.,50.},{50.,60.},{60.,70.},{70.,80.},{80.,90.}};
//const float  kPtRebin[kCentLength] = {2.6,2.6,2.2,2.2,2.2,2.,0.6,3.4};

// Fit setup
const float  kFitminPt = -2.0f;
const float  kFitmaxPt = 3.0f;
const float  kFitmaxNBkg = 1.0e+9f;

const float  kTPCmaxPt = 1.4f;
const float  kTOFminPt = 0.6f;
const float  kPtRange[2] = {0.6,8};
const float  kPtRangeMatCorrection[2] = {0.65,1.55};
const float  kPtRangeMatCorrectionTPC[2] = {0.65,1.55};

const bool   kUseBarlow{true};
const bool   kSmoothSystematics{true};
const float  kAbsSyst[2] = {0.08,0.1f};

// Plotter
const int kPtBinLimit[kCentLength] = {28,28,28,28,28,28,28,28,28,28};

// Secondaries
const int knDCAbins = 52;
const double kDcabins[53] = {
    -1.30,-1.20,-1.10,-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,
    -0.35,-0.30,-0.25,-0.20,-0.15,-0.12,-0.10,-0.09,-0.08,-0.07,
    -0.06,-0.05,-0.04,-0.03,-0.02,-0.01, 0.00, 0.01, 0.02, 0.03,
     0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.15, 0.20,
     0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00,
     1.10, 1.20, 1.30
  };
  

const map<string,vector<float> > kCutNames {
	{"dcaz",{0.75f,1.25f,1.5f,2.f}},
	{"tpc",{60.f,65.f,75.f,80.f}},
	{"pid",{25.f,35.f}}
};


#endif
