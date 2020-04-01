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
const string kNamePlot[2] = {"deuterons", "anti-deuterons"};
const char* kSymbols[2] = {"d","#bar{d}"};

const string kBaseOutputDir = "../results/pp13TeV/";
//const string kBaseOutputDir = "../results/pp5TeV/";
const string kBaseInputDir = "../";

const string kDataFilename = kBaseInputDir + "data/pp13TeV/multiplicity/treno20190427/AnalysisResults.root";//treno20190227/data.root";//"data/latest_pp13TeV/data.root";
const string kDataFilenameMB = kBaseInputDir + "data/pp13TeV/MB/treno20190510/AnalysisResults.root";
const string kMCfilename = kBaseInputDir + "mc/pp13TeV/multiplicity/treno20190427/AnalysisResults.root";//"mc/latest_pp13TeV/mc_enriched.root";
const string kMCfilenameMB = kBaseInputDir + "mc/pp13TeV/MB/treno20190515/AnalysisResults.root";
//const string kDataFilename = kBaseInputDir + "data/latest_pp5TeV/data.root";
//const string kMCfilename = kBaseInputDir + "mc/latest_pp5TeV/mc.root";
const string kMCMaterialBudget1 = kBaseInputDir + "mc/MaterialBudget/LHC16h7c.root";
const string kMCMaterialBudget2 = kBaseInputDir + "mc/MaterialBudget/LHC17d5a.root";
const string kMCMaterialBudget3 = kBaseInputDir + "mc/MaterialBudget/LHC17d5b.root";

const string kMCgeant3 = kBaseInputDir + "mc/geant_stuff/treno20190518/geant3.root";
const string kMCgeant4 = kBaseInputDir + "mc/geant_stuff/treno20190518/geant4.root";

const string kFilterListNames = "nuclei_deuterons_";
const string kFilterListNamesMCasData = "nuclei_deuteronsMCasData_";
const string kNormalisationList = "mpuccio_NucleiPIDqa_default";

const string kEfficiencyOutput = kBaseOutputDir + "efficiency.root";
const string kSignalOutput = kBaseOutputDir + "signal.root";
const string kSignalMCOutput = kBaseOutputDir + "signalMC.root";
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
const string kFinalCheckOutput = kBaseOutputDir + "final_check.root";
const string kJoinSystematicsOutput = kBaseOutputDir + "joinsystematics.root";
const string kBWfitsOutput = kBaseOutputDir + "BWfits.root";
const string kSignalLossInput = kBaseOutputDir + "SignalLoss.root";
const string kSignalLossOutput = kBaseOutputDir + "SignalLoss_Correction.root";
const string kRatioOutput = kBaseOutputDir + "ratio.root";
const string kG3G4output = kBaseOutputDir + "g3g4syst.root";

// const string kEfficiencyOutput = kBaseOutputDir + "efficiency_2017.root";
// const string kSignalOutput = kBaseOutputDir + "signal_2017.root";
// const string kSecondariesOutput = kBaseOutputDir + "secondaries.root";
// const string kSecondariesOutputRooFit = kBaseOutputDir + "RooSec_2017.root";
// const string kSecondariesTPCoutput = kBaseOutputDir + "secondaries_TPC.root";
// const string kSecondariesTPCoutputRooFit = kBaseOutputDir + "RooSecTPC_2017.root";
// const string kMaterialOutput = kBaseOutputDir + "materialbudget.root";
// const string kSpectraOutput = kBaseOutputDir + "spectra_2017.root";
// const string kFitSystematicsOutput = kBaseOutputDir + "fitsystematics_2017.root";
// const string kSystematicsOutput = kBaseOutputDir + "systematics_2017.root";
// const string kSystematicsOutputTPC = kBaseOutputDir + "systematics_TPC_2017.root";
// const string kFinalOutput = kBaseOutputDir + "final_2017.root";
// const string kJoinSystematicsOutput = kBaseOutputDir + "joinsystematics_2017.root";
// const string kBWfitsOutput = kBaseOutputDir + "BWfits_2017.root";
// const string kSignalLossInput = kBaseOutputDir + "SignalLoss.root";
// const string kSignalLossOutput = kBaseOutputDir + "SignalLoss_Correction.root";

const bool kUseIntegratedForMB = true;

const bool   kPrintFigures{true};
const string kFiguresFolder = kBaseOutputDir + "images/";
const string kMacrosFolder = kBaseOutputDir + "images/macros/";

const int    kNPtBins = 17;
const double kPtBins[kNPtBins+1] = {0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8,4.4f,5.0f,6.0f};
const int    kNCentralityBins = 12;
const float  kCentralityBins[kNCentralityBins+1] = {0.f,1.f,5.f,10.f,20.f,30.f,40.f,50.f,60.f,70.f,80.f,90.f,100.f};

const int    kCentLength = 10;

const float kTOFlowPt = 1.0;
const float kBinCountingCut[kCentLength] = {1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,2.0,1.8};
const float kNoBkgTOF[kCentLength] = {1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,2.0,1.8};
const float kSingleExpBkg = 2.2;
const float kSingleExpSideBkg = 2.2;
const float kFixSigma = 2.0;

const float kNoBkgTPC = 0.9;

const float kLimitsTPC[2] = {-5.6,5.6};
const float kLimitsTOF[2] = {-1.20, 1.76};
const float kLimitsTOFbkg[2] = {-1.20, 1.76};
const int    kCentBinsArray[kCentLength][2] = {{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,10},{11,13},{2,13}};
const float  kCentPtLimits[kCentLength] = {3.8,3.4,3.4,3.4,3.,3.,2.6,2.6,2.6,3.8};
//const float  kCentPtLimits[kCentLength] = {3.4,3.4,3.4,3.4,3.,3.,3.,3.,2.6,3.8};
const float  kCentLabels[kCentLength][2] = {{0.,1.},{1.,5.},{5.,10.},{10.,20.},{20.,30.},{30.,40.},{40.,50.},{50,70},{70.,100.},{0.,100.}};
const float  kPtRebin[kCentLength] = {2.6,2.6,2.6,2.6,2.6,2.6,2.2,2.,1.8,2.6};

const float  kTPCmaxPt = 1.2f;
const float  kTOFminPt = 0.9f;
const float  kPtRange[2] = {0.7,6.};
const float  kPtHadronicInteractionRange[2] = {0.9,1.4};
const float  kPtRangeMatCorrectionTOF[2] = {0.8,1.6};
const float  kPtRangeMatCorrectionTPC[2] = {0.6,1.2};

const int kNDCAbins = 28;
const double kDCAbins[kNDCAbins+1] = {-1.20,-1.10,-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20};

const bool   kUseBarlow{true};
const bool   kSmoothSystematics{true};
const float  kAbsSyst[2] = {0.04,0.075};
const float  kAbsCorr[2] = {1.118,1.118};
const float  kAbsCorrMax[2] = {1.076,1.153};
const float  kCorrG3G4tpc[2] = {1.,1.1};
const float  kCorrG3G4tof[2] = {1.11,1.26};

//const float  kAbsCorr[2] = {0.803091,0.841291};
//const float  kAbsSyst[2] = {0.049459, 0.020119};
//const float  kAbsCorr[2] = {0.694284, 0.886095};

const map<string,vector<float> > kCutNames {
 {"dcaxy",{0.10f,0.14f}},
 {"dcaz",{0.5f,0.75f,1.5f,2.f}},
 {"tpc",{60.f,65.f,75.f,80.f}},
 {"pid",{3.25f,3.50f}}
};

const char* kRomanLabels[10] = {"I","II","III","IV + V","VI","VII","VIII","IX","X","MB"};
const int kScaleFactor[10] = {256,128,64,32,16,8,4,2,1,8192};
const int kExponent[10] = {8,7,6,5,4,3,2,1,0,13};

const int kNsigmaVar = 5; 
//const int kNshiftVar = 4;
const int kNshiftVar = 5;
int vSigmaWidth[kNsigmaVar] = {0,1,-1,2,-2};
//float vSigmaShift[kNsigmaVar] = {-0.1,-0.05,0.05,0.1};
int vSigmaShift[kNshiftVar] = {0,-2,-1,1,2};

Double_t sigmoid_function(Double_t *x, Double_t *par) {
    Float_t xx =x[0];
    Double_t f = 1./(1.+par[0]*TMath::Exp(par[1]*xx));
    return f;
}

#endif
