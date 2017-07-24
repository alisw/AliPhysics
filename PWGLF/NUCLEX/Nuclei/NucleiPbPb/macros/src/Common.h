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

const string kBaseOutputDir = "/Users/lbariogl/cernbox/Deuterons13TeV/macros/results/";
const string kBaseInputDir = "/Users/lbariogl/cernbox/Deuterons13TeV/";

const string kDataFilename = kBaseInputDir + "Dati/data.root";
const string kMCfilename = kBaseInputDir + "mc.root";

const string kFilterListNames = "mpuccio_deuterons_";
const string kNormalisationList = "mpuccio_NucleiPIDqa";

const string kEfficiencyOutput = kBaseOutputDir + "efficiency.root";
const string kSignalOutput = kBaseOutputDir + "signal.root";
const string kSecondariesOutput = kBaseOutputDir + "secondaries.root";
const string kSecondariesTPCoutput = kBaseOutputDir + "secondaries_TPC.root";
const string kMaterialOutput = kBaseOutputDir + "materialbudget.root";
const string kSpectraOutput = kBaseOutputDir + "spectra.root";
const string kFitSystematicsOutput = kBaseOutputDir + "fitsystematics.root";
const string kSystematicsOutput = kBaseOutputDir + "systematics.root";
const string kSystematicsOutputTPC = kBaseOutputDir + "systematics_TPC.root";
const string kFinalOutput = kBaseOutputDir + "final.root";

const bool   kPrintFigures{true};
const string kFiguresFolder = "/Users/lbariogl/cernbox/Deuterons13TeV/macros/results/images/";
const string kMacrosFolder = "/Users/lbariogl/cernbox/Deuterons13TeV/macros/results/images/macros/";

const float  kTPCmaxPt = 1.4f;
const float  kTOFminPt = 1.f;
const float  kPtRange[2] = {0.6,3.8};
const float  kPtRangeMatCorrection[2] = {0.7,1.3};

const bool   kUseBarlow{true};
const bool   kSmoothSystematics{true};
const float  kAbsSyst[2] = {0.08,0.1f};

const map<string,vector<float> > kCutNames {
	{"dcaz",{0.75f,1.25f,1.5f,2.f}},
	{"tpc",{60.f,65.f,75.f,80.f}},
	{"pid",{25.f,35.f}}
};


#endif
