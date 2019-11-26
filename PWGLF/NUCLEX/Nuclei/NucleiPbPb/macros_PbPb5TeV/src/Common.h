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
const string kMCfilename = kBaseInputDir + "mc/PbPb5TeV/treno20191123/LHC19d2.root";

const string kFilterListNames = "mpuccio_deuterons_";
const string kNormalisationList = "mpuccio_NucleiPIDqa";

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

const float  kPtBins[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};
const int    kNPtBins = 15;
const float  kCentralityBins[13] = {0.f,1.f,5.f,10.f,20.f,30.f,40.f,50.f,60.f,70.f,80.f,90.f,100.f};
const int    kNCentralityBins = 12;

const int    kCentLength = 8;
const int    kCentBinsArray[kCentLength][2] = {{2,2},{3,3},{4,4},{5,5},{6,6},{7,8},{9,10},{2,10}};
const float  kCentPtLimits[kCentLength] = {3.4,3.,3.,2.6,2.6,2.2,2.,3.8};
const float  kCentLabels[kCentLength][2] = {{0.,5.},{5.,10.},{10.,20.},{20.,30.},{30.,40.},{40.,60.},{60.,100.},{0.,100.}};
const float  kPtRebin[kCentLength] = {2.6,2.6,2.2,2.2,2.2,2.,0.6,3.4};

const float  kTPCmaxPt = 1.4f;
const float  kTOFminPt = 1.f;
const float  kPtRange[2] = {0.6,3.8};
const float  kPtRangeMatCorrection[2] = {1.05,1.55};
const float  kPtRangeMatCorrectionTPC[2] = {0.65,1.35};

const bool   kUseBarlow{true};
const bool   kSmoothSystematics{true};
const float  kAbsSyst[2] = {0.08,0.1f};

const map<string,vector<float> > kCutNames {
	{"dcaz",{0.75f,1.25f,1.5f,2.f}},
	{"tpc",{60.f,65.f,75.f,80.f}},
	{"pid",{25.f,35.f}}
};


#endif
