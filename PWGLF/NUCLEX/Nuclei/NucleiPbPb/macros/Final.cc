#include "src/Common.h"
#include "src/Plotting.h"
#include "src/Utils.h"
using namespace utils;

#include <map>
#include <vector>
#include <array>
using std::array;
#include <memory>

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>

#include <AliPID.h>
#include <AliPWGFunc.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <TVirtualFitter.h>

//Defining Levi-Tsallis and parameters
const Float_t normal=2e-4, normMin=1e-5, normMax=1.;
const Float_t n=7, nMin=2, nMax=100;
const Float_t C=0.2, CMin=0.01, CMax=0.6;
const int n_parameters = 5; //parameters of the function
enum e_parameters { e_mass, e_n, e_C, e_norm, e_chi2};
const char * param_names[5] ={"mass","n","C","norm","chi2"};

TF1* LevyTsallis(const Char_t *name, Double_t mass, Double_t n, Double_t nMin, Double_t nMax, Double_t C, Double_t CMin, Double_t CMax, Double_t normal, Double_t normMin, Double_t normMax);
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p);

void Final() {
  TFile spectra_file(kSpectraOutput.data());
  TFile TOF_systematics_file(kSystematicsOutput.data());
  TFile TPC_systematics_file(kSystematicsOutputTPC.data());
  TFile final_file(kFinalOutput.data(),"recreate");

  TVirtualFitter::SetDefaultFitter("Minuit");

  const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
  const int n_pt_bins = 15;

  const int n_centralities = kCentLength;

  const char* syst_names[5] = {"cutsyst","countsyst","abssyst","matsyst","shiftsyst"};

  enum syst_enum { cutsyst, countsyst, abssyst, matsyst, shiftsyst};

  //histograms fot TOF analysis
  TH1F* stat_tof[2][n_centralities];
  TH1F* syst_tof[2][n_centralities];
  TH1F* vec_syst_tof[2][n_centralities][5];
  //histograms for TPC analysis
  TH1F* stat_tpc[2][n_centralities];
  TH1F* syst_tpc[2][n_centralities];
  TH1F* vec_syst_tpc[2][n_centralities][5];
  //inclusive histogram with both tpc and TOF
  TH1F* stat_all[2][n_centralities];
  TH1F* syst_all[2][n_centralities];

  TH1F* all[2][n_centralities];

  //fit function and parameters;
  TF1* fit_function[2][n_centralities];
  TH1F* fit_parameters[2][n_parameters];
  TH1F* ratio_fit_data[2][n_centralities];

  //AliPWGFunc funcMaker; // function for fit to the Spectra
  //funcMaker.SetVarType(AliPWGFunc::kdNdpt);

  TH1F* comp[2][n_centralities];

  RooRealVar x_pull("x_pull","#frac{#mu - 1}{#sigma}",-4.,4.);
  std::unique_ptr<RooDataSet> ratio_pull[n_centralities];
  for(int iC = 0; iC < n_centralities; iC++){
    ratio_pull[iC] = make_unique<RooDataSet>(Form("pull_ratio_%d",iC),Form("pull_ratio_%d",iC),RooArgSet(x_pull));
  }
  RooRealVar pull_mean("#mu","#mu",0,-2,2);
  RooRealVar pull_sigma("#sigma","#sigma",1.,0.1,3.);
  RooGaussian pull_gauss("pull_gauss","pull_gauss",x_pull,pull_mean,pull_sigma);

  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = final_file.mkdir(kNames[iS].data());
    //mass parameter
    fit_parameters[iS][e_mass] = new TH1F(Form("fit_param_%d_%s",iS,param_names[e_mass]),";iC;mass (GeV/#it{c}^{2})",n_centralities,0.,float(n_centralities));
    fit_parameters[iS][e_n] = new TH1F(Form("fit_param_%d_%s",iS,param_names[e_n]),";iC;n",n_centralities,0.,float(n_centralities));
    fit_parameters[iS][e_C] = new TH1F(Form("fit_param_%d_%s",iS,param_names[e_C]),";iC;C (GeV)",n_centralities,0.,float(n_centralities));
    fit_parameters[iS][e_norm] = new TH1F(Form("fit_param_%d_%s",iS,param_names[e_norm]),";iC;#frac{d#it{N}}{d#it{y}}",n_centralities,0.,float(n_centralities));
    fit_parameters[iS][e_chi2] = new TH1F(Form("fit_param_%d_%s",iS,param_names[e_chi2]),";iC;#chi^{2}/NDF",n_centralities,0.,float(n_centralities));
    for (int iC = 0; iC < n_centralities; ++iC) {
      fit_parameters[iS][e_mass]->GetXaxis()->SetBinLabel(iC+1,Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]));
      fit_parameters[iS][e_norm]->GetXaxis()->SetBinLabel(iC+1,Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]));
      fit_parameters[iS][e_n]->GetXaxis()->SetBinLabel(iC+1,Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]));
      fit_parameters[iS][e_C]->GetXaxis()->SetBinLabel(iC+1,Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]));
      fit_parameters[iS][e_chi2]->GetXaxis()->SetBinLabel(iC+1,Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]));
    }
    for (int iC = 0; iC < n_centralities; ++iC) {
      TDirectory *c_dir = s_dir->mkdir(to_string(iC).data());
      c_dir->cd();
      ratio_fit_data[iS][iC] = new TH1F(Form("ratio_fit_data_%d_%d",iS,iC),";#it{p}_{T}(GeV/#it{c});data / fit",n_pt_bins,pt_bin_limits);
      ratio_fit_data[iS][iC]->GetYaxis()->SetRangeUser(0.,2.);
      TH1F* totsyst_tmp = (TH1F*)TOF_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst").data());
      Requires(totsyst_tmp, "Missing totsyt");
      TH1F* totsyst = (TH1F*)totsyst_tmp->Rebin(n_pt_bins,Form("totsyst_%d_%d",iS,iC),pt_bin_limits);
      TH1F* totsyst_tpc_tmp = (TH1F*)TPC_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst_tpc").data());
      Requires(totsyst_tpc_tmp, "Missing totsyt_tpc");
      TH1F* totsyst_tpc = (TH1F*)totsyst_tpc_tmp->Rebin(n_pt_bins,Form("totsyst_tpc_%d_%d",iS,iC),pt_bin_limits);
      string tof_basepath = kFilterListNames + "/" + kNames[iS] + "/TOFspectra" + to_string(iC);
      TH1F* spectra_tof_tmp  = (TH1F*)spectra_file.Get(tof_basepath.data());
      Requires(spectra_tof_tmp,tof_basepath.data());
      stat_tof[iS][iC] =(TH1F*)spectra_tof_tmp->Rebin(n_pt_bins,Form("stat_%d_%d",iS,iC),pt_bin_limits);
      auto ptAxis = stat_tof[iS][iC]->GetXaxis();
      string tpc_basepath = kFilterListNames + "/" + kNames[iS] + "/TPCspectra" + to_string(iC);
      TH1F* spectra_tpc_tmp  = (TH1F*)spectra_file.Get(tpc_basepath.data());
      Requires(spectra_tpc_tmp,tpc_basepath.data());
      stat_tpc[iS][iC] = (TH1F*)spectra_tpc_tmp->Rebin(n_pt_bins,Form("stat_tpc_%d_%d",iS,iC),pt_bin_limits);
      syst_tof[iS][iC]  = (TH1F*)totsyst->Clone(("syst_tof" + to_string(iC)).data());
      syst_tof[iS][iC]->Reset();
      syst_tpc[iS][iC]  = (TH1F*)totsyst->Clone(("syst_tpc" + to_string(iC)).data());
      syst_tpc[iS][iC]->Reset();
      stat_all[iS][iC] = (TH1F*)totsyst->Clone(("stat_all" + to_string(iC)).data());
      stat_all[iS][iC]->Reset();
      syst_all[iS][iC] = (TH1F*)totsyst->Clone(("syst_all" + to_string(iC)).data());
      syst_all[iS][iC]->Reset();
      all[iS][iC] = (TH1F*)totsyst->Clone(("all" + to_string(iC)).data());
      all[iS][iC]->Reset();
      //syst_all[iS][iC] = (TH1F*)totsyst->Clone(("syst_all" + to_string(iC)).data());
      //syst_all[iS][iC]->Reset();
      for(int iSyst=0; iSyst<5; iSyst++){
        TH1F* partial_syst_tmp = (TH1F*) TOF_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst]).data());
        vec_syst_tof[iS][iC][iSyst] = (TH1F*) partial_syst_tmp->Rebin(n_pt_bins,Form("%s_%d_%d",syst_names[iSyst],iS,iC),pt_bin_limits);
        TH1F* partial_syst_tpc_tmp = (TH1F*) TPC_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst] + "_tpc").data());
        if(iSyst==shiftsyst) vec_syst_tpc[iS][iC][iSyst] = nullptr;
        else{
          vec_syst_tpc[iS][iC][iSyst] = (TH1F*) partial_syst_tpc_tmp->Rebin(n_pt_bins,Form("%s_tpc_%d_%d",syst_names[iSyst],iS,iC),pt_bin_limits);
        }
      }
      // Compute the minimum of the systematic uncertainty related to the material budget
      float min_material_tof = 1.,min_material_tpc = 1.;
      float tof_tmp = 0., tpc_tmp = 0.;
      for (int iB = 1; iB <= n_pt_bins; ++iB) {
        tof_tmp = vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB);
        tpc_tmp = vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB);
        if(tof_tmp>1e-8){
          if(tof_tmp<min_material_tof) min_material_tof = tof_tmp;
        }
        if(tpc_tmp>1e-8){
          if(tpc_tmp<min_material_tpc) min_material_tpc = tpc_tmp;
        }
      }
      comp[iS][iC] = new TH1F(Form("comp_%d_%d",iS,iC),";#it{p}_{T} (GeV/#it{c}); #frac{N_{TPC} - N_{TOF}}{#sigma}",2,1.,1.2);
      plotting::SetHistStyle(comp[iS][iC],plotting::kSpectraColors[iC]);
      // Fit function
      fit_function[iS][iC]=LevyTsallis(Form("Levi-Tsallis_%d_%d",iS,iC), AliPID::ParticleMass(AliPID::kDeuteron), n, nMin, nMax, C, CMin, CMax, normal, normMin, normMax);
      //funcMaker.GetTsallis(AliPID::ParticleMass(AliPID::kDeuteron),.1,1./0.9,1.,Form("func_%d_%d",iS,iC));
      fit_function[iS][iC]->SetLineColor(kBlack);
      for (int iB = 1; iB <= n_pt_bins; ++iB) {
        if(ptAxis->GetBinCenter(iB)<1.){
          stat_tof[iS][iC]->SetBinContent(iB,0.);
          stat_tof[iS][iC]->SetBinError(iB,0.);
          syst_tof[iS][iC]->SetBinContent(iB,0.);
          syst_tof[iS][iC]->SetBinError(iB,0.);
          //
          syst_tpc[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          syst_tpc[iS][iC]->SetBinError(iB,totsyst_tpc->GetBinContent(iB) * stat_tpc[iS][iC]->GetBinContent(iB));
          //
          stat_all[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          stat_all[iS][iC]->SetBinError(iB,stat_tpc[iS][iC]->GetBinError(iB));
          syst_all[iS][iC]->SetBinContent(iB,syst_tpc[iS][iC]->GetBinContent(iB));
          syst_all[iS][iC]->SetBinError(iB,syst_tpc[iS][iC]->GetBinError(iB));
          //
          all[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
          all[iS][iC]->SetBinError(iB,TMath::Sqrt(stat_tpc[iS][iC]->GetBinError(iB)*stat_tpc[iS][iC]->GetBinError(iB) +                      syst_tpc[iS][iC]->GetBinError(iB)*syst_tpc[iS][iC]->GetBinError(iB)-kAbsSyst[iS]*kAbsSyst[iS]*stat_tpc[iS][iC]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)-min_material_tpc*min_material_tpc*stat_tpc[iS][iC]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)));
        }
        else{
          stat_tpc[iS][iC]->SetBinContent(iB,0.);
          stat_tpc[iS][iC]->SetBinError(iB,0.);
          syst_tpc[iS][iC]->SetBinContent(iB,0.);
          syst_tpc[iS][iC]->SetBinError(iB,0.);
          //
          syst_tof[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          syst_tof[iS][iC]->SetBinError(iB,totsyst->GetBinContent(iB) * stat_tof[iS][iC]->GetBinContent(iB));
          //
          stat_all[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          stat_all[iS][iC]->SetBinError(iB,stat_tof[iS][iC]->GetBinError(iB));
          syst_all[iS][iC]->SetBinContent(iB,syst_tof[iS][iC]->GetBinContent(iB));
          syst_all[iS][iC]->SetBinError(iB,syst_tof[iS][iC]->GetBinError(iB));
          //
          all[iS][iC]->SetBinContent(iB,stat_tof[iS][iC]->GetBinContent(iB));
          all[iS][iC]->SetBinError(iB,TMath::Sqrt(stat_tof[iS][iC]->GetBinError(iB)*stat_tof[iS][iC]->GetBinError(iB) +            syst_tof[iS][iC]->GetBinError(iB)*syst_tof[iS][iC]->GetBinError(iB)-kAbsSyst[iS]*kAbsSyst[iS]*stat_tof[iS][iC]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)-min_material_tof*min_material_tof*stat_tof[iS][iC]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)));
        }
        if(ptAxis->GetBinCenter(iB)>1. && ptAxis->GetBinCenter(iB)<1.2){ //comparing TPC and TOF where both
          float z = stat_tpc[iS][iC]->GetBinContent(iB) - stat_tof[iS][iC]->GetBinContent(iB);
          float err2 = (stat_tof[iS][iC]->GetBinError(iB)-stat_tpc[iS][iC]->GetBinError(iB))*
                       (stat_tof[iS][iC]->GetBinError(iB)-stat_tpc[iS][iC]->GetBinError(iB)) +

                      (vec_syst_tof[iS][iC][cutsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][cutsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
                      (vec_syst_tof[iS][iC][cutsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][cutsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

                      (vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][abssyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
                      (vec_syst_tof[iS][iC][abssyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][abssyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

                      (vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB))*
                      (vec_syst_tof[iS][iC][matsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) - vec_syst_tpc[iS][iC][matsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)) +

                      vec_syst_tof[iS][iC][countsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)*
                      vec_syst_tof[iS][iC][countsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) +

                      vec_syst_tpc[iS][iC][countsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB)*
                      vec_syst_tpc[iS][iC][countsyst]->GetBinContent(iB)*stat_tpc[iS][iC]->GetBinContent(iB) +

                      vec_syst_tof[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB)*
                      vec_syst_tof[iS][iC][shiftsyst]->GetBinContent(iB)*stat_tof[iS][iC]->GetBinContent(iB) +

                      0.03*stat_tof[iS][iC]->GetBinContent(iB)*0.03*stat_tof[iS][iC]->GetBinContent(iB);
          z /= TMath::Sqrt(err2);
          comp[iS][iC]->SetBinError(comp[iS][iC]->FindBin(ptAxis->GetBinCenter(iB)),0.);
          comp[iS][iC]->SetBinContent(comp[iS][iC]->FindBin(ptAxis->GetBinCenter(iB)),z);

        }
      }
      plotting::SetHistStyle(stat_tof[iS][iC],plotting::kSpectraColors[iC]);
      plotting::SetHistStyle(syst_tof[iS][iC],plotting::kSpectraColors[iC]);
      stat_tof[iS][iC]->Scale(1<<(kCentLength-iC-1));
      stat_tof[iS][iC]->Write("stat_tof");
      syst_tof[iS][iC]->Scale(1<<(kCentLength-iC-1));
      syst_tof[iS][iC]->Write("syst_tof");
      plotting::SetHistStyle(stat_tpc[iS][iC],plotting::kSpectraColors[iC],25);
      plotting::SetHistStyle(syst_tpc[iS][iC],plotting::kSpectraColors[iC],25);
      stat_tpc[iS][iC]->Scale(1<<(kCentLength-iC-1));
      stat_tpc[iS][iC]->Write("stat_tpc");
      syst_tpc[iS][iC]->Scale(1<<(kCentLength-iC-1));
      syst_tpc[iS][iC]->Write("syst_tpc");
      all[iS][iC]->Scale(1<<(kCentLength-iC-1));
      all[iS][iC]->Fit(Form("Levi-Tsallis_%d_%d",iS,iC),"IQ");
      plotting::SetHistStyle(stat_all[iS][iC],plotting::kSpectraColors[iC],25);
      plotting::SetHistStyle(syst_all[iS][iC],plotting::kSpectraColors[iC],25);
      stat_all[iS][iC]->Write("stat_all");
      syst_all[iS][iC]->Write("syst_all");
      //filling ratio_fit_data[iS][iC]
      for (int i = 1; i <= n_pt_bins; ++i) {
        ratio_fit_data[iS][iC]->SetBinContent(i, all[iS][iC]->GetBinContent(i) / fit_function[iS][iC]->Eval(all[iS][iC]->GetBinCenter(i)));
        ratio_fit_data[iS][iC]->SetBinError(i, all[iS][iC]->GetBinError(i) / fit_function[iS][iC]->Eval(all[iS][iC]->GetBinCenter(i)));
      }
      ratio_fit_data[iS][iC]->Write();
      fit_parameters[iS][e_mass]->SetBinContent(iC+1,fit_function[iS][iC]->GetParameter(e_mass));
      fit_parameters[iS][e_mass]->SetBinError(iC+1,fit_function[iS][iC]->GetParError(e_mass));
      fit_parameters[iS][e_n]->SetBinContent(iC+1,fit_function[iS][iC]->GetParameter(e_n));
      fit_parameters[iS][e_n]->SetBinError(iC+1,fit_function[iS][iC]->GetParError(e_n));
      fit_parameters[iS][e_C]->SetBinContent(iC+1,fit_function[iS][iC]->GetParameter(e_C));
      fit_parameters[iS][e_C]->SetBinError(iC+1,fit_function[iS][iC]->GetParError(e_C));
      fit_parameters[iS][e_norm]->SetBinContent(iC+1,fit_function[iS][iC]->GetParameter(e_norm));
      fit_parameters[iS][e_norm]->SetBinError(iC+1,fit_function[iS][iC]->GetParError(e_norm));
      fit_parameters[iS][e_chi2]->SetBinContent(iC+1,fit_function[iS][iC]->GetChisquare()/fit_function[iS][iC]->GetNDF());
      //syst_all[iS][iC]->Scale(1<<(kCentLength-iC-1));

    }

    s_dir->cd();
    fit_parameters[iS][e_mass]->Write();
    fit_parameters[iS][e_n]->Write();
    fit_parameters[iS][e_C]->Write();
    fit_parameters[iS][e_norm]->Write();
    fit_parameters[iS][e_chi2]->Write();
    TCanvas spectra("spectra","spectra",3200,2400);
    spectra.DrawFrame(
        0.4,
        1e-6,
        3.8,
        1,
        ";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}}"
        );
    spectra.SetLeftMargin(0.1425591);
    spectra.SetRightMargin(0.027121);
    spectra.SetTopMargin(0.06053269);
    spectra.SetBottomMargin(0.1598063);
    TLegend final_leg(0.771116,0.202655,0.927529,0.5451330);
    final_leg.SetBorderSize(0);
    final_leg.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kNames[iS].data()));
    final_leg.SetNColumns(2);
    final_leg.AddEntry((TObject*)nullptr,"TPC","");
    final_leg.AddEntry((TObject*)nullptr,"TPC + TOF","");
    for (int iC = 0; iC < n_centralities -1; ++iC) {
      stat_tof[iS][iC]->Draw("esamex0");
      syst_tof[iS][iC]->Draw("e2same");
      stat_tpc[iS][iC]->Draw("esamex0");
      syst_tpc[iS][iC]->Draw("e2same");
      //all[iS][iC]->Draw("esamex0");
      fit_function[iS][iC]->SetRange(0.5,1.05*kCentPtLimits[iC]);
      fit_function[iS][iC]->SetLineStyle(2);
      fit_function[iS][iC]->Draw("lsame");
      final_leg.AddEntry(syst_tpc[iS][iC],Form("%4.0f - %2.0f %% (#times %d)",kCentLabels[iC][0],kCentLabels[iC][1],1<<(kCentLength-iC-1)),"fp");
      final_leg.AddEntry(syst_tof[iS][iC],Form("%4.0f - %2.0f %% (#times %d)",kCentLabels[iC][0],kCentLabels[iC][1],1<<(kCentLength-iC-1)),"fp");
    }
    final_leg.Draw();
    spectra.SetLogy();
    spectra.Write();
    if (kPrintFigures) {
      spectra.SaveAs((kFiguresFolder + "spectraTOF" + kLetter[iS] + ".eps").data());
      spectra.SaveAs((kMacrosFolder + "spectraTOF" + kLetter[iS] + ".C").data());
    }
    // comparing tpc and tof where both present
    TCanvas cComp("cComp","cComp",3200,2400);
    cComp.DrawFrame(0.95,-4.,1.25,4.,";#it{p}_{T} (GeV/#it{c});#frac{N_{TPC} - N_{TOF}}{#sigma}");
    TLegend comp_leg(0.76,0.61,0.90,0.89);
    comp_leg.SetBorderSize(0);
    comp_leg.SetHeader(Form("%s, pp #sqrt{s} = 13 TeV",kNames[iS].data()));
    for (int iC = 0; iC < n_centralities -1; ++iC) {
      comp[iS][iC]->Draw("psame");
      comp_leg.AddEntry(comp[iS][iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"p");
    }
    for(int i=-3; i<=3; i++){
      TLine *line = new TLine(0.95,i,1.25,i);
      int absi = (i<0) ? -1*i : i;
      switch(absi){
        case(3):
          line->SetLineColor(kRed);
          break;
        case(2):
          line->SetLineColor(kOrange);
          break;
        case(1):
          line->SetLineColor(kGreen+3);
          break;
        default:
          line->SetLineColor(kBlack);
          break;
      }
      line->SetLineStyle(2);
      line->Draw();
    }
    comp_leg.Draw();
    cComp.Write();
  }

  TDirectory* r_dir = final_file.mkdir("ratio");
  TF1* funcpol[n_centralities];

  r_dir->cd();
  TPad* pad[9] = {nullptr};
  TCanvas ratio("ratio","ratio",3200,3200);
  plotting::CanvasPartition(&ratio,pad,3,3);

  for (int iC = 0; iC < n_centralities -1; ++iC) {
    r_dir->mkdir(to_string(iC).data())->cd();
    stat_tof[1][iC]->Divide(stat_tof[0][iC]);
    syst_tof[1][iC]->Divide(syst_tof[0][iC]);
    all[1][iC]->Divide(all[0][iC]);
    funcpol[iC] = new TF1(Form("ratiopol_%d",iC),"pol0",0.6,kCentPtLimits[iC]);
    int nx = iC/3;
    int ny = iC%3;
    pad[iC]->cd();
    all[1][iC]->Fit(Form("ratiopol_%d",iC),"Q");
    for(int iB=1; iB<=n_pt_bins; iB++){
      if(stat_tof[1][iC]->GetBinCenter(iB)<1.){
        stat_tof[1][iC]->SetBinContent(iB, 0.);
        stat_tof[1][iC]->SetBinError(iB, 0.);
        syst_tof[1][iC]->SetBinContent(iB, 0.);
        syst_tof[1][iC]->SetBinError(iB, 0.);
      }
      else{
        if(TMath::Abs(stat_tof[1][iC]->GetBinContent(iB)) < 1e-8) break;
        x_pull.setVal( (stat_tof[1][iC]->GetBinContent(iB) - 1) / stat_tof[1][iC]->GetBinError(iB));// TMath::Sqrt(stat_tof[1][iC]->GetBinError(iB)*stat_tof[1][iC]->GetBinError(iB) + syst_tof[1][iC]->GetBinError(iB)*syst_tof[1][iC]->GetBinError(iB)) );
        ratio_pull[iC]->add(RooArgSet(x_pull));
      }
    }
    stat_tof[1][iC]->Write("stat_tof");
    syst_tof[1][iC]->Write("syst_tof");

    stat_tpc[1][iC]->Divide(stat_tpc[0][iC]);
    syst_tpc[1][iC]->Divide(syst_tpc[0][iC]);
    for(int iB=1; iB<=n_pt_bins; iB++){
      if(stat_tpc[1][iC]->GetBinCenter(iB)>1.2){
        stat_tpc[1][iC]->SetBinContent(iB, 0.);
        stat_tpc[1][iC]->SetBinError(iB, 0.);
        syst_tpc[1][iC]->SetBinContent(iB, 0.);
        syst_tpc[1][iC]->SetBinError(iB, 0.);
      }
      else{
        x_pull.setVal( (stat_tpc[1][iC]->GetBinContent(iB) - 1) / TMath::Sqrt(stat_tpc[1][iC]->GetBinError(iB)*stat_tpc[1][iC]->GetBinError(iB) + syst_tpc[1][iC]->GetBinError(iB)*syst_tpc[1][iC]->GetBinError(iB)) );
        ratio_pull[iC]->add(RooArgSet(x_pull));
      }
    }

    TCanvas cPull("pull","pull",3200,2400);
    RooMsgService::instance().setSilentMode(true);
    pull_gauss.fitTo(*ratio_pull[iC],RooFit::PrintEvalErrors(-1));
    RooPlot* pull_frame = x_pull.frame();
    pull_frame->SetTitle("");
    ratio_pull[iC]->plotOn(pull_frame,RooFit::Name("data"));
    pull_gauss.plotOn(pull_frame,RooFit::Name("model"));
    float chi2 = pull_frame->chiSquare("model","data");
    pull_gauss.paramOn(pull_frame, RooFit::Label(Form("#chi^{2}/NDF = %2.4f",chi2)),RooFit::Layout(0.71,0.90,0.85));
    pull_frame->getAttLine()->SetLineWidth(0);
    pull_frame->SetMinimum(1e-5);
    pull_frame->Draw();
    cPull.Write();

    //Fill the unique canvas with all the manti-matter/matter ratios
    pad[iC]->cd();
    double XaxisEdge = 0.;
    switch (nx) {
      case 0:
        XaxisEdge=3.9;
        break;
      case 1:
        XaxisEdge=3.2;
        break;
      default:
        XaxisEdge=2.4;
        break;
    }
    gPad->DrawFrame(
        0.35,
        0.1,
        XaxisEdge,
        1.9,
    ";#it{p}_{T} (GeV/#it{c});#bar{d}/d"
    );
    stat_tof[1][iC]->Draw("esamex0");
    syst_tof[1][iC]->Draw("e2same");
    stat_tpc[1][iC]->Draw("esamex0");
    syst_tpc[1][iC]->Draw("e2same");
    //all[1][iC]->Draw("esamex0");
    funcpol[iC]->Draw("same");
    TLine *line_one = new TLine(0.35,1.,XaxisEdge,1.);
    line_one->SetLineColor(kBlack);
    line_one->SetLineStyle(2);
    line_one->Draw();
    TLegend* ratio_leg_one = new TLegend(0.719499,0.128318,0.920602,0.337758);
    ratio_leg_one->SetHeader(Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]));
    ratio_leg_one->SetBorderSize(0);
    ratio_leg_one->AddEntry(syst_tof[0][iC],"TPC + TOF","p");
    ratio_leg_one->AddEntry(syst_tpc[0][iC],"TPC","p");
    ratio_leg_one->AddEntry((TObject*)0, Form("p0: %.2f #pm %.2f", funcpol[iC]->GetParameter(0),funcpol[iC]->GetParError(0)), "");
    ratio_leg_one->AddEntry((TObject*)0, Form("#chi^{2}/NDF: %.2f / %d",funcpol[iC]->GetChisquare(),funcpol[iC]->GetNDF()),"");
    ratio_leg_one->Draw();
    // ratio->Update();

  }
  r_dir->cd();
  if (kPrintFigures) ratio.SaveAs((kFiguresFolder + "ratio.eps").data());
  ratio.Write();


}

//________________________________________________________________________________________________
TF1* LevyTsallis(const Char_t *name, Double_t mass,
          Double_t n, Double_t nMin, Double_t nMax,
          Double_t C, Double_t CMin, Double_t CMax,
          Double_t norm, Double_t normMin, Double_t normMax)
{

  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->SetParameters(mass, n, C, norm);

  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, nMin, nMax);
  fLevyTsallis->SetParLimits(2, CMin, CMax);
  fLevyTsallis->SetParLimits(3, normMin, normMax);

  return fLevyTsallis;
}
//________________________________________________________________________________________________
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{

  /* dN/dpt */

  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;

}
