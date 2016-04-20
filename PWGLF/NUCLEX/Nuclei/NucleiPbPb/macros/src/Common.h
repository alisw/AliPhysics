#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TList.h>
#include <TObject.h>

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
const string kSystematicsOutput = kBaseOutputDir + "systematics.root";
const string kFinalOutput = kBaseOutputDir + "final.root";

const float  kTPCmaxPt = 1.2f;
const float  kPtRange[2] = {4.,6.};
const float  kPtRangeMatCorrection[2] = {0.6,1.3};

const bool   kUseBarlow{false};
const float  kMatSyst = 0.04f;
const float  kAbsSyst = 0.1f;

const map<string,vector<float> > kCutNames {
	{"chisquare",{3.5f,4.5f,5.f,5.5f,6.f}},
	{"dcaz",{0.75f,1.25f,1.5f,2.f}},
	{"tpc",{60.f,65.f,75.f,80.f}},
	{"pid",{25.f,35.f}}
};

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

#endif
