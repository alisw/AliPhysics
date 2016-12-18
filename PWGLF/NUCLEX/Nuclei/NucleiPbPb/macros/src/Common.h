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

const string kBaseOutputDir = "/tmp/";
const string kBaseInputDir = "~/Desktop/Repositories/root-files/deuterons/";

const string kDataFilename = kBaseInputDir + "data.root";
const string kMCfilename = kBaseInputDir + "mc.root";

const string kFilterListNames = "mpuccio_deuterons6cent";

const string kEfficiencyOutput = kBaseOutputDir + "efficiency.root";
const string kSignalOutput = kBaseOutputDir + "signal.root";
const string kSecondariesOutput = kBaseOutputDir + "secondaries.root";
const string kSpectraOutput = kBaseOutputDir + "spectra.root";
const string kFitSystematicsOutput = kBaseOutputDir + "fitsystematics.root";
const string kSystematicsOutput = kBaseOutputDir + "systematics.root";
const string kFinalOutput = kBaseOutputDir + "final.root";

const bool   kPrintFigures{true};
const string kFiguresFolder = "/tmp/images/";
const string kMacrosFolder = "/tmp/images/macros/";

const float  kTPCmaxPt = 1.2f;
const float  kPtRange[2] = {1.,6.};
const float  kPtRangeMatCorrection[2] = {0.8,1.3};

const bool   kUseBarlow{false};
const bool   kUseBarlowFit{false};
const bool   kSmoothSystematics{true};
const float  kMatSyst = 0.03f;
const float  kAbsSyst[2] = {0.08,0.1f};

const map<string,vector<float> > kCutNames {
	{"chisquare",{3.5f,4.5f,5.f,5.5f,6.f}},
	{"dcaz",{0.75f,1.25f,1.5f,2.f}},
	{"tpc",{60.f,65.f,75.f,80.f}},
	{"pid",{25.f,35.f}}
};

TF1 G3G4corr[2] = {
  {"g3g4M","0.9438455-0.006023916*x+0.2508776*TMath::Exp(-0.7569186*x)",0.5,6},
  {"g3g4A","0.8451751+0.009187905*x+2.289737*TMath::Exp(-3.02719*x)",0.5,6}};

template<class T> void Requires(T* obj, string msg = "") {
  if (!obj || obj == nullptr) {
    std::cout << "Missing object. " << msg.data() << std::endl;
    abort();
  }
}

class TTList : public TList {
  public:
    TObject* Get(std::string name) {
      TObject* obj = this->FindObject(name.data());
      Requires(obj);
      return obj;
    }
};

template<typename F> F RobustRMS(vector<F> &data) {
  F mean = TMath::Mean(data.size(),data.data());
  F rms = TMath::RMS(data.size(),data.data());
  F robust = 0.f;
  F norm = 0.f;
  for (const auto& val : data) {
    const F delta = TMath::Abs(val - mean);
    if (delta < 3 * rms) {
      robust += delta * delta;
      norm++;
    }
  }
  return norm > 0 ? TMath::Sqrt(robust) / norm : 0.f;
}

#endif
