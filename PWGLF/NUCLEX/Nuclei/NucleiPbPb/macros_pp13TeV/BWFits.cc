#include "src/Plotting.h"
using namespace plotting;

#include "src/Utils.h"
using namespace utils;

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
#include "src/YieldMean.C"
#include "src/YieldMeanNew.cc"
#include "AliPWGFunc.h"
#include "AdditionalFunctions.h"

#include <cstdio>
#include <fstream>

static TF1 *fBGBlastWave_Integrand = NULL;
static TF1 *fBGBlastWave_Integrand_num = NULL;
static TF1 *fBGBlastWave_Integrand_den = NULL;

constexpr double kParticleMass = 1.87561;
constexpr int kNfitFunctions = 4;
const string kFitFunctionNames[kNfitFunctions] = {"LevyTsallis", "Boltzmann", "Mt-exp", "Pt-exp"};

//Levy-Tsallis
enum e_param_tsallis {e_mass, e_n, e_C, e_norm, e_chi2};
const char * param_names_tsallis[5] = {"mass","n","C","norm","chi2"};

//Levi-Tsallis parameters
const float normal=2e-4, normMin=1e-6, normMax=1.;
const float n=10, nMin=2, nMax=100;
const float C=0.2, CMin=0.01, CMax=0.4;

//Boltzmann and MT-scaling and Pt-scaling
enum e_param_bolmt {b_norm,b_T,b_chi2};
const char* param_names_bolmt[3] = {"norm","T","chi2"};

void Denormalize(TH1 * h) {
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    h->SetBinContent(i,h->GetBinContent(i) * TMath::TwoPi() * h->GetBinCenter(i));
    h->SetBinError(i,h->GetBinError(i) * TMath::TwoPi() * h->GetBinCenter(i));
  }
}

void BWFits(bool antimatter_analysys = false) {
  const char* kind_of_particle = (antimatter_analysys) ? "antideuterons" : "deuterons";
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.);
  TFile *input_file = TFile::Open(kFinalOutput.data());
  TH1D *stat[kCentLength],*syst[kCentLength],*syst_pt_uncorr[kCentLength],*syst_pt_corr[kCentLength],*syst_mult_corr[kCentLength],*syst_mult_uncorr[kCentLength];

  remove(Form("%s/%s_log.root",kBaseOutputDir.data(),kind_of_particle));
  std::ofstream output_param(Form("%s/%s_fits_param.txt",kBaseOutputDir.data(),kind_of_particle));

  AliPWGFunc pwgfunc;
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  TF1* fit_functions[kNfitFunctions] = {
    LevyTsallis(kFitFunctionNames[0].data(), kParticleMass),
    pwgfunc.GetBoltzmann(kParticleMass, 0.1, 1, kFitFunctionNames[1].data()),
    pwgfunc.GetMTExp(kParticleMass, 0.1, 1, kFitFunctionNames[2].data()),
    pwgfunc.GetPTExp(0.1, 1, kFitFunctionNames[3].data()),
  };
  //levi-tsallis parameters'
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
  //pt exponential
  double ptexp_param_val[kCentLength][3];
  double ptexp_param_err[kCentLength][3];
  const char* ptexp_param_name[3] = {"$A (GeV^{-1}c$)","$T$ (GeV)","\\chi^2/ndf"};
  //
  TFile bwfile(Form("%s/%s_fits.root",kBaseOutputDir.data(),kind_of_particle),"recreate");
  TDirectory* datadir = bwfile.mkdir("data");
  TDirectory* function_dir[4]{nullptr};
  for (int iF = 0; iF < kNfitFunctions; ++iF)
    function_dir[iF] = bwfile.mkdir(kFitFunctionNames[iF].data());
  for (int iC = 0; iC < kCentLength; ++iC) {
    stat[iC] = (TH1D*)input_file->Get(Form("%s/%i/stat",kind_of_particle,iC));
    Requires(stat[iC],Form("%s/%i/stat",kind_of_particle,iC));
    syst[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst",kind_of_particle,iC));
    Requires(syst[iC],Form("%s/%i/syst",kind_of_particle,iC));
    syst_pt_uncorr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_pt_uncorr",kind_of_particle,iC));
    Requires(syst_pt_uncorr[iC],Form("%s/%i/syst_pt_uncorr",kind_of_particle,iC));
    syst_pt_corr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_pt_corr",kind_of_particle,iC));
    Requires(syst_pt_corr[iC],Form("%s/%i/syst_pt_corr",kind_of_particle,iC));
    syst_mult_corr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_mult_corr",kind_of_particle,iC));
    Requires(syst_mult_corr[iC],Form("%s/%i/syst_mult_corr",kind_of_particle,iC));
    syst_mult_uncorr[iC] = (TH1D*)input_file->Get(Form("%s/%i/syst_mult_uncorr",kind_of_particle,iC));
    Requires(syst_mult_uncorr[iC],Form("%s/%i/syst_mult_uncorr",kind_of_particle,iC));
    TH1D *mStat = new TH1D(Form("mStat%d",iC),";#it{p}_{T} (GeV / #it{c}); 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-1}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystTot = new TH1D(Form("mSystTot%d",iC),";#it{p}_{T} (GeV/#it{c});1/(2#pi#it{p}_{T}) 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-2}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystPtUncorr = new TH1D(Form("mSystPtUncorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/(2#pi#it{p}_{T}) 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-2}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystPtCorr = new TH1D(Form("mSystPtCorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/(2#pi#it{p}_{T}) 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-2}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystMultCorr = new TH1D(Form("mSystMultCorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/(2#pi#it{p}_{T}) 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-2}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());
    TH1D *mSystMultUncorr = new TH1D(Form("mSystMultUncorr%d",iC),";#it{p}_{T} (GeV/#it{c});1/(2#pi#it{p}_{T}) 1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} (GeV/#it{c})^{-2}",stat[iC]->GetNbinsX(), stat[iC]->GetXaxis()->GetXbins()->GetArray());

    for (int j = 1; j <= stat[iC]->GetNbinsX(); ++j) {
      double x = stat[iC]->GetBinCenter(j);
      if (x < 0.6 || x > kCentPtLimits[iC]) continue;
      float val = stat[iC]->GetBinContent(stat[iC]->FindBin(mStat->GetBinCenter(j)));
      float stat_err = stat[iC]->GetBinError(stat[iC]->FindBin(mStat->GetBinCenter(j)));
      float syst_err = syst[iC]->GetBinError(stat[iC]->FindBin(mStat->GetBinCenter(j)));
      float syst_pt_corr_err = syst_pt_corr[iC]->GetBinError(stat[iC]->FindBin(mStat->GetBinCenter(j)));
      float syst_pt_uncorr_err = syst_pt_uncorr[iC]->GetBinError(stat[iC]->FindBin(mStat->GetBinCenter(j)));
      float syst_mult_corr_err = syst_mult_corr[iC]->GetBinError(stat[iC]->FindBin(mStat->GetBinCenter(j)));
      float syst_mult_uncorr_err = syst_mult_uncorr[iC]->GetBinError(stat[iC]->FindBin(mStat->GetBinCenter(j)));
      mStat->SetBinContent(j,val);
      mStat->SetBinError(j,stat_err);
      mSystTot->SetBinContent(j,val);
      mSystTot->SetBinError(j,syst_err);
      mSystPtUncorr->SetBinContent(j,val);
      mSystPtUncorr->SetBinError(j,syst_pt_uncorr_err);
      mSystPtCorr->SetBinContent(j,val);
      mSystPtCorr->SetBinError(j,syst_pt_corr_err);
      mSystMultCorr->SetBinContent(j,val);
      mSystMultCorr->SetBinError(j,syst_mult_corr_err);
      mSystMultUncorr->SetBinContent(j,val);
      mSystMultUncorr->SetBinError(j,syst_mult_uncorr_err);
      
    }
    for (int iF = 0; iF < kNfitFunctions; ++iF) {
      if(iF==0){
        fit_functions[iF]->SetParameters(kParticleMass, n, C, normal);
        fit_functions[iF]->SetParLimits(1, nMin, nMax);
        fit_functions[iF]->SetParLimits(2, CMin, CMax);
        fit_functions[iF]->SetParLimits(3, normMin, normMax);
      }
      TF1* fout;
      //TH1* h = yieldmean::YieldMean(mStat,mSystTot,fout,fit_functions[iF],0,10.,0.01,0.1,true,Form("%s/%slog.root",kBaseOutputDir.data(),kind_of_particle),Form("%s/%d", kFitFunctionNames[iF].data(),iC));
      TH1* h = yieldmeannew::YieldMeanNew(mStat,mSystTot,mSystPtUncorr,mSystPtCorr,mSystMultCorr,fout,fit_functions[iF],0,10.,0.01,0.1,true,Form("%s%s_log.root",kBaseOutputDir.data(),kind_of_particle),Form("%s/%d", kFitFunctionNames[iF].data(),iC));
      cout << "\n*****************************" << endl;
      printf("Function: %s\n", kFitFunctionNames[iF].data());
      printf("iC: %d\n", iC);
      for (int iP = 0; iP < fout->GetNpar(); ++iP) {
        cout << fout->GetParName(iP) << ": " << fout->GetParameter(iP);
        cout << " +/- " << fout->GetParError(iP) << endl;
      }
      cout << "*****************************" << endl << endl;
      function_dir[iF]->mkdir(Form("%d",iC));
      function_dir[iF]->cd(Form("%d",iC));
      fout->Write(Form("%s%d",kFitFunctionNames[iF].data(),iC));
      if(iF==0){
        for(int iParam=0; iParam<4; iParam++){
          tsallis_param_val[iC][iParam] = fout->GetParameter(iParam);
          tsallis_param_err[iC][iParam] = fout->GetParError(iParam);
        }
        tsallis_param_val[iC][4] = fout->GetChisquare()/fout->GetNDF();
        tsallis_param_err[iC][4] = 0;
      }
      else if(iF==1){
        for(int iParam=0; iParam<2; iParam++){
          boltzmann_param_val[iC][iParam] = fout->GetParameter(iParam);
          boltzmann_param_err[iC][iParam] = fout->GetParError(iParam);
        }
        boltzmann_param_val[iC][2] = fout->GetChisquare()/fout->GetNDF();
        boltzmann_param_err[iC][2] = 0;
      }
      else if(iF==2){
        for(int iParam=0; iParam<2; iParam++){
          mtexp_param_val[iC][iParam] = fout->GetParameter(iParam);
          mtexp_param_err[iC][iParam] = fout->GetParError(iParam);
        }
        mtexp_param_val[iC][2] = fout->GetChisquare()/fout->GetNDF();
        mtexp_param_err[iC][2] = 0;
      }
      else{
        for(int iParam=0; iParam<2; iParam++){
          ptexp_param_val[iC][iParam] = fout->GetParameter(iParam);
          ptexp_param_err[iC][iParam] = fout->GetParError(iParam);
        }
        ptexp_param_val[iC][2] = fout->GetChisquare()/fout->GetNDF();
        ptexp_param_err[iC][2] = 0;
      }
      h->Write(Form("result%d",iC));
    }
    datadir->mkdir(Form("%d",iC));
    datadir->cd(Form("%d",iC));
    mStat->Write(Form("stat%d",iC));
    mSystTot->Write(Form("syst%d",iC));
    mSystPtUncorr->Write(Form("syst_mult_uncorr%d",iC));
    mSystMultCorr->Write(Form("syst_mult_corr%d",iC));
    mSystMultUncorr->Write(Form("syst_mult_uncorr%d",iC));
  }
  //Writign Tsallis parameters
  output_param << "****************************************\n";
  output_param << "\t\tLevy-Tsallis\n" ;
  output_param << "****************************************\n";
  output_param << "Multiplicity (%) & ";
  for(int iCol=0; iCol<4; iCol++){
    output_param << tsallis_param_name[iCol] << "\t& ";
  }
  output_param << "\\chi^{2} / ndf \\\\" << endl;
  for(int iC=0; iC<kCentLength; iC++){
    output_param << Form("%4.0f - %2.0f",kCentLabels[iC][0],kCentLabels[iC][1]) << "\t& ";
    for(int iParam=0; iParam<4; iParam++){
      output_param << tsallis_param_val[iC][iParam] << " $\\pm$ " << tsallis_param_err[iC][iParam] << "\t& ";
    }
    output_param << tsallis_param_val[iC][4] << " \\\\"<<endl;;
  }
  //Writign Boltzmann parameters
  output_param<<"****************************************\n";
  output_param<<"\t\tBoltzmann\n";
  output_param<<"****************************************\n";
  output_param<< "Multiplicity (%) & ";
  for(int iCol=0; iCol<2; iCol++){
    output_param<< boltzmann_param_name[iCol] << "\t& ";
  }
  output_param<< "\\chi^{2} / ndf \\\\" << endl;
  for(int iC=0; iC<kCentLength; iC++){
    output_param<< Form("%4.0f - %2.0f",kCentLabels[iC][0],kCentLabels[iC][1]) << "\t& ";
    for(int iParam=0; iParam<2; iParam++){
      output_param<< boltzmann_param_val[iC][iParam] << " $\\pm$ " << boltzmann_param_err[iC][iParam] << "\t& ";
    }
    output_param<< boltzmann_param_val[iC][2] << " \\\\"<<endl;;
  }
  //Writign Mt-exponential parameters
  output_param<<"****************************************\n";
  output_param<<"\t\tMt-exponential\n";
  output_param<<"****************************************\n";
  output_param<< "Multiplicity (%) & ";
  for(int iCol=0; iCol<2; iCol++){
    output_param<< mtexp_param_name[iCol] << "\t& ";
  }
  output_param<< "\\chi^{2} / ndf \\\\" << endl;
  for(int iC=0; iC<kCentLength; iC++){
    output_param<< Form("%4.0f - %2.0f",kCentLabels[iC][0],kCentLabels[iC][1]) << "\t& ";
    for(int iParam=0; iParam<2; iParam++){
      output_param<< mtexp_param_val[iC][iParam] << " $\\pm$ " << mtexp_param_err[iC][iParam] << "\t& ";
    }
    output_param<< mtexp_param_val[iC][2] << " \\\\"<<endl;;
  }

  //Writign Pt-exponential parameters
  output_param<<"****************************************\n";
  output_param<<"\t\tPt-exponential\n";
  output_param<<"****************************************\n";
  output_param<< "Multiplicity (%) & ";
  for(int iCol=0; iCol<2; iCol++){
    output_param<< ptexp_param_name[iCol] << "\t& ";
  }
  output_param<< "\\chi^{2} / ndf \\\\" << endl;
  for(int iC=0; iC<kCentLength; iC++){
    output_param<< Form("%4.0f - %2.0f",kCentLabels[iC][0],kCentLabels[iC][1]) << "\t& ";
    for(int iParam=0; iParam<2; iParam++){
      output_param<< ptexp_param_val[iC][iParam] << " $\\pm$ " << ptexp_param_err[iC][iParam] << "\t& ";
    }
    output_param<< ptexp_param_val[iC][2] << " \\\\"<<endl;;
  }
  output_param.close();

  bwfile.Close();
}
