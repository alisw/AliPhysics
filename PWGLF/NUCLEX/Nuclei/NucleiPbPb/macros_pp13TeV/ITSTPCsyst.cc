#include "src/Common.h"
#include "src/FitModules.h"
#include "src/Utils.h"
using namespace utils;

#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>

const char* periods[21] = {"LHC16d", "LHC16e", "LHC16g", "LHC16h", "LHC16i", "LHC16j", "LHC16k", "LHC16l", "LHC16o", "LHC16p",
"LHC17c", "LHC17e", "LHC17f", "LHC17g", "LHC17i", "LHC17j", "LHC17k", "LHC17l", "LHC17m", "LHC17o", "LHC17r"};

double events[21] = {};

double systematics_input[21][4] = {{0.7,2.6,4.9,4.8},{0.62,1.92,3.16,3.1},{0.83,2.85,4.83,4.46},{1.00,2.65,4.20,3.88},{0.99,0.85,0.91,0.41},{0.88,0.77,0.73,0.54},{0.74,0.58,0.88,1.12},{0.80,0.82,1.35,1.97},{0.64,0.47,0.74,1.13},{0.60,0.54,0.87,1.18},{0.68,1.63,2.30,2.44},{0.85,1.72,2.48,2.36},{0.88,1.77,2.55,2.42},{0.91,1.77,2.57,2.77},{-1,-1,-1,-1},{0.93,1.82,2.54,2.68},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}};

double systematics_output[4] = {0.,0.,0.,0.,};

double weights[21]{0.};

double total_events = 0.;

void ITSTPCsyst(){

  for(int iPeriod = 0; iPeriod<21; iPeriod++){
    TFile* fInput = new TFile(Form("../data/pp13TeV/multiplicity/treno20190427/%s/AnalysisResults.root",periods[iPeriod]));
    TTList* list = (TTList*)fInput->Get("NucleiPIDqa_default_default/mpuccio_NucleiPIDqa_default");
    TH1D* fCutStats = (TH1D*)list->Get("fCutStats");
    double counts = fCutStats->GetBinContent(fCutStats->GetNbinsX());
    events[iPeriod] = counts;
    if(systematics_input[iPeriod][0]<0) continue;
    total_events += events[iPeriod];
  }
  for(int iPeriod=0; iPeriod<21; iPeriod++){
    if(systematics_input[iPeriod][0]<0){
      weights[iPeriod] = 0.;
    } else {
      weights[iPeriod] = events[iPeriod]/total_events;
    }
    for(int iPt=0; iPt<4; iPt++){
      systematics_output[iPt] += systematics_input[iPeriod][iPt]*weights[iPeriod];
    }
  }
  std::cout << "ITS-TPC matching mean systematics" << std::endl;
  std::cout << "[0.8,1.] GeV/c : " << systematics_output[0] << std::endl;
  std::cout << "[1.,2.] GeV/c : " << systematics_output[1] << std::endl;
  std::cout << "[2.,3.] GeV/c : " << systematics_output[2] << std::endl;
  std::cout << "[3.,4.] GeV/c : " << systematics_output[3] << std::endl;
}