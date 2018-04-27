#include "src/Plotting.h"
using namespace plotting;

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <Riostream.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include "src/Common.h"
#include "YieldMean.C"
#include "AliPWGFunc.h"
#include "AdditionalFunctions.h"

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;

constexpr double kParticleMass = 1.87561;
constexpr int kNfitFunctions = 3;
const string kFitFunctionNames[kNfitFunctions] = {"LevyTsallis", "Boltzmann", "Mt-exp"};

//Levy-Tsallis
enum e_param_tsallis {e_mass, e_n, e_C, e_norm, e_chi2};
const char * param_names_tsallis[5] = {"mass","n","C","norm","chi2"};

//Levi-Tsallis parameters
const Float_t normal=2e-4, normMin=1e-6, normMax=1.;
const Float_t n=10, nMin=2, nMax=100;
const Float_t C=0.2, CMin=0.01, CMax=0.4;

//Boltzmann and MT-scaling
enum e_param_bolmt {b_norm,b_T,b_chi2};
const char* param_names_bolmt[3] = {"norm","T","chi2"};

void Denormalize(TH1 * h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    h->SetBinContent(i,h->GetBinContent(i) * TMath::TwoPi() * h->GetBinCenter(i));
    h->SetBinError(i,h->GetBinError(i) * TMath::TwoPi() * h->GetBinCenter(i));
  }
}

void BWFits(bool antimatter_analysys = false) {
  const char* kind_of_particle = (antimatter_analysys) ? "anti" : "";
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.);
  TFile *mineF = TFile::Open(kFinalOutput.data());
  TCanvas *mResults[kCentLength];
  TCanvas *mCanvM[kCentLength],*mCanvA[kCentLength];
  TH1D * mysystM[kCentLength], *mysystA[kCentLength],*speM[kCentLength],*speA[kCentLength];
  TF1 *bw[3];

  AliPWGFunc pwgfunc;
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  TF1* fit_functions[kNfitFunctions] = {
    LevyTsallis(kFitFunctionNames[0].data(), kParticleMass),
    pwgfunc.GetBoltzmann(kParticleMass, 0.1, 1, kFitFunctionNames[1].data()),
    pwgfunc.GetMTExp(kParticleMass, 0.1, 1, kFitFunctionNames[2].data())
  };
  //levi-tsallis parameters
  double tsallis_param_val[kCentLength][5];
  double tsallis_param_err[kCentLength][5];
  const char* tsallis_param_name[5] = {"$m$ ($Gev/$c^{2}$)","$n$","$C$ (GeV)" ,"\\mathrm{d}N/\\mathrm{d}y","\\chi^2 ndf"};
  //Boltzmann parameters
  double boltzmann_param_val[kCentLength][3];
  double boltzmann_param_err[kCentLength][3];
  const char* boltzmann_param_name[3] = {"$A (GeV^{-2}c^{3}$)","$T$ (GeV)","\\chi^2/ndf"};
  //mt exponential
  double mtexp_param_val[kCentLength][3];
  double mtexp_param_err[kCentLength][3];
  const char* mtexp_param_name[3] = {"$A (GeV^{-1}c$)","$T$ (GeV)","\\chi^2/ndf"};
  //
  TFile bwfile(Form("%s/%sfits.root",kBaseOutputDir.data(),kind_of_particle),"recreate");
  TDirectory* datadir = bwfile.mkdir("data");
  TDirectory* function_dir[4]{nullptr};
  for (int iF = 0; iF < kNfitFunctions; ++iF)
    function_dir[iF] = bwfile.mkdir(kFitFunctionNames[iF].data());
  for (int iC = 0; iC < kCentLength; ++iC) {
    speM[iC] = (TH1D*)mineF->Get(Form("%sdeuterons/%i/stat_all",kind_of_particle,iC));
    mysystM[iC] = (TH1D*)mineF->Get(Form("%sdeuterons/%i/syst_all",kind_of_particle,iC));
    if (!mysystM[iC]) cout << "Missing " << Form("syst%d",iC) << endl;
    if (!speM[iC]) cout << "Missing " << Form("stat%d",iC) << endl;
    if (!speM[iC] || !mysystM[iC]) return;
    TH1D *m = new TH1D(Form("m%d",iC),";#it{p}_{T} (GeV / #it{c}); 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",speM[iC]->GetNbinsX(), speM[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *sm = new TH1D(Form("sm%d",iC),";#it{p}_{T} (GeV/#it{c});1/(2#pi#it{p}_{T}) 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-2}",speM[iC]->GetNbinsX(), speM[iC]->GetXaxis()->GetXbins()->GetArray());

    for (int j = 1; j <= speM[iC]->GetNbinsX(); ++j) {
      double x = speM[iC]->GetBinCenter(j);
      if (x < 0.6 || x > kCentPtLimits[iC]) continue;
      float stat = speM[iC]->GetBinError(speM[iC]->FindBin(m->GetBinCenter(j)));
      float syst = mysystM[iC]->GetBinError(speM[iC]->FindBin(m->GetBinCenter(j)));
      m->SetBinContent(j,speM[iC]->GetBinContent(speM[iC]->FindBin(m->GetBinCenter(j))));
      m->SetBinError(j,stat);
      sm->SetBinContent(j,speM[iC]->GetBinContent(speM[iC]->FindBin(m->GetBinCenter(j))));
      sm->SetBinError(j,syst);
    }

    for (int iF = 0; iF < kNfitFunctions; ++iF) {
      if(iF==0){
        fit_functions[iF]->SetParameters(kParticleMass, n, C, normal);
        fit_functions[iF]->SetParLimits(1, nMin, nMax);
        fit_functions[iF]->SetParLimits(2, CMin, CMax);
        fit_functions[iF]->SetParLimits(3, normMin, normMax);
      }
      function_dir[iF]->cd();
      TH1* h = YieldMean(m,sm,fit_functions[iF],0,10.1);
      cout << "\n*****************************" << endl;
      printf("Function: %s\n", kFitFunctionNames[iF].data());
      printf("iC: %d\n", iC);
      for (int iP = 0; iP < fit_functions[iF]->GetNpar(); ++iP) {
        cout << fit_functions[iF]->GetParName(iP) << ": " << fit_functions[iF]->GetParameter(iP);
        cout << " +/- " << fit_functions[iF]->GetParError(iP) << endl;
      }
      cout << "*****************************" << endl << endl;
      fit_functions[iF]->Write(Form("%s%d",kFitFunctionNames[iF].data(),iC));
      if(iF==0){
        for(int iParam=0; iParam<4; iParam++){
          tsallis_param_val[iC][iParam] = fit_functions[iF]->GetParameter(iParam);
          tsallis_param_err[iC][iParam] = fit_functions[iF]->GetParError(iParam);
        }
        tsallis_param_val[iC][4] = fit_functions[iF]->GetChisquare()/fit_functions[iF]->GetNDF();
        tsallis_param_err[iC][4] = 0;
      }
      else if(iF==1){
        for(int iParam=0; iParam<2; iParam++){
          boltzmann_param_val[iC][iParam] = fit_functions[iF]->GetParameter(iParam);
          boltzmann_param_err[iC][iParam] = fit_functions[iF]->GetParError(iParam);
        }
        boltzmann_param_val[iC][2] = fit_functions[iF]->GetChisquare()/fit_functions[iF]->GetNDF();
        boltzmann_param_err[iC][2] = 0;
      }
      else{
        for(int iParam=0; iParam<2; iParam++){
          mtexp_param_val[iC][iParam] = fit_functions[iF]->GetParameter(iParam);
          mtexp_param_err[iC][iParam] = fit_functions[iF]->GetParError(iParam);
        }
        mtexp_param_val[iC][2] = fit_functions[iF]->GetChisquare()/fit_functions[iF]->GetNDF();
        mtexp_param_err[iC][2] = 0;
      }
      h->Write(Form("result%d",iC));
    }
    datadir->cd();
    m->Write(Form("stat%d",iC));
    sm->Write(Form("syst%d",iC));
  }
  //Writign Tsallis parameters
  printf("****************************************\n");
  printf("\t\tLevy-Tsallis\n");
  printf("****************************************\n");
  std::cout << "Multiplicity (%) & ";
  for(int iCol=0; iCol<4; iCol++){
    std::cout << tsallis_param_name[iCol] << "\t& ";
  }
  std::cout << "\\chi^{2} / ndf \\\\" << endl;
  for(int iC=0; iC<kCentLength; iC++){
    cout << Form("%4.0f - %2.0f",kCentLabels[iC][0],kCentLabels[iC][1]) << "\t& ";
    for(int iParam=0; iParam<4; iParam++){
      std::cout << tsallis_param_val[iC][iParam] << " $\\pm$ " << tsallis_param_err[iC][iParam] << "\t& ";
    }
    std::cout << tsallis_param_val[iC][4] << " \\\\"<<endl;;
  }
  //Writign Boltzmann parameters
  printf("****************************************\n");
  printf("\t\tBoltzmann\n");
  printf("****************************************\n");
  std::cout << "Multiplicity (%) & ";
  for(int iCol=0; iCol<2; iCol++){
    std::cout << boltzmann_param_name[iCol] << "\t& ";
  }
  std::cout << "\\chi^{2} / ndf \\\\" << endl;
  for(int iC=0; iC<kCentLength; iC++){
    cout << Form("%4.0f - %2.0f",kCentLabels[iC][0],kCentLabels[iC][1]) << "\t& ";
    for(int iParam=0; iParam<2; iParam++){
      std::cout << boltzmann_param_val[iC][iParam] << " $\\pm$ " << boltzmann_param_err[iC][iParam] << "\t& ";
    }
    std::cout << boltzmann_param_val[iC][2] << " \\\\"<<endl;;
  }
  //Writign Mt-exponential parameters
  printf("****************************************\n");
  printf("\t\tMt-exponential\n");
  printf("****************************************\n");
  std::cout << "Multiplicity (%) & ";
  for(int iCol=0; iCol<2; iCol++){
    std::cout << mtexp_param_name[iCol] << "\t& ";
  }
  std::cout << "\\chi^{2} / ndf \\\\" << endl;
  for(int iC=0; iC<kCentLength; iC++){
    cout << Form("%4.0f - %2.0f",kCentLabels[iC][0],kCentLabels[iC][1]) << "\t& ";
    for(int iParam=0; iParam<2; iParam++){
      std::cout << mtexp_param_val[iC][iParam] << " $\\pm$ " << mtexp_param_err[iC][iParam] << "\t& ";
    }
    std::cout << mtexp_param_val[iC][2] << " \\\\"<<endl;;
  }

  bwfile.Close();
}
