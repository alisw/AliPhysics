#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"
using namespace plotting;

#include <TFile.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TLegendEntry.h>
#include "Math/QuantFuncMathCore.h"

#include <string>
using std::string;
#include <fstream>

const char* kPrefix[2] = {"","anti"};
const int kScaleFactorB2[10] = {1,2,4,8,16,32,64,128,256,512};

int contatore = 0;

const char* kParticleName[2] = {"Deuterons", "Antideuterons"};

double GetErrorY(const TGraphErrors* gr, double x0);
void MakeItInvariant(TH1* d_stat);
TGraphErrors* GetBAPt(const TGraphErrors* grPrtInvDYieldPt,
    const TGraphErrors* grNucInvDYieldPt, const TString& name, bool syst=true);

const char* kFixedPtLabel[3] = {"75","55","105"};

constexpr int kInputCent = 10;
constexpr int kMultInput[11]{0,1,5,10,15,20,30,40,50,70,100};
constexpr int myCent[kCentLength][2]{
  {0,0}, /// 0-1%
  {1,1}, /// 1-5%
  {2,2}, /// 5-10%
  {3,4}, /// 10-20%
  {5,5}, /// 20-30%
  {6,6}, /// 30-40%
  {7,7}, /// 40-50%
  {8,8}, /// 50-70%
  {9,9}, /// 70-100%
  {0,9}  /// 0-100%
};
void myB2(){

  double dNdEta_tmp, dNdEtaErr_tmp, proton_yields_tmp, proton_yields_stat_tmp, proton_yields_syst_tmp, proton_yields_syst_uncorr_tmp, proton_yields_syst_corr_tmp;

  std::vector<double> dNdEta_vec, dNdEtaErr_vec, proton_yields_vec, proton_yields_stat_vec, proton_yields_syst_vec;

  std::ifstream protonFile(Form("%sproton_yields.txt",kBaseOutputDir.data()));
  if(!protonFile.is_open()){
    printf("The file %s could not be opened\n", Form("%sproton_yields.txt",kBaseOutputDir.data()));
    exit(1);
  }
  else{
    while(protonFile >> dNdEta_tmp >> dNdEtaErr_tmp >> proton_yields_tmp >> proton_yields_stat_tmp >> proton_yields_syst_tmp >> proton_yields_syst_uncorr_tmp >> proton_yields_syst_corr_tmp){
      dNdEta_vec.push_back(dNdEta_tmp);
      dNdEtaErr_vec.push_back(dNdEtaErr_tmp);
      proton_yields_vec.push_back(proton_yields_tmp);
      proton_yields_stat_vec.push_back(proton_yields_stat_tmp);
      proton_yields_syst_vec.push_back(proton_yields_syst_tmp);
    }
  }

  double dNdEta[kCentLength-1], dNdEtaErr[kCentLength-1], proton_yields[kCentLength-1], proton_yields_stat[kCentLength-1], proton_yields_syst[kCentLength-1];

  int iC=0;
  for(int i=0; i< (int) dNdEta_vec.size();i++){
    if(i!=4){
      dNdEta[iC] = dNdEta_vec[i];
      dNdEtaErr[iC] = dNdEtaErr_vec[i];
      proton_yields[iC] = proton_yields_vec[i];
      proton_yields_stat[iC] = proton_yields_stat_vec[i];
      proton_yields_syst[iC] = proton_yields_syst_vec[i];
      if (i!=3) iC++;
    }
    else{
      dNdEta[iC] = (dNdEta_vec[i]+dNdEta[iC])/2;
      dNdEtaErr[iC] = TMath::Sqrt(Sq(dNdEtaErr_vec[i]) + Sq(dNdEtaErr[iC]))/2;
      proton_yields[iC] = (proton_yields_vec[i]+proton_yields[iC])/2;
      proton_yields_stat[iC] = TMath::Sqrt(Sq(proton_yields_stat_vec[i]) + Sq(proton_yields_stat[iC]))/2;
      proton_yields_syst[iC] = TMath::Sqrt(Sq(proton_yields_syst_vec[i]) + Sq(proton_yields_syst[iC]))/2;
      iC++;
    }
  }

  TFile deuteron_file(kFinalOutput.data());
  TFile output_file(Form("%s/B2.root",kBaseOutputDir.data()),"recreate");

  TFile input_file(Form("%sFinal_combined_spectra_TPCTOFTOFonlyrTPCKinksITSsa_pp13TeV.root",kBaseOutputDir.data()));

  TGraphErrors* grB2atPtFixedStat[2][3];
  TGraphErrors* grB2atPtFixedSyst[2][3];
  //
  TGraphErrors* grB2atPtFixedStatMB[2][3];
  TGraphErrors* grB2atPtFixedSystMB[2][3];
  for(int iS=0; iS<2; iS++){
    for(int iL=0; iL<3; iL++){
      grB2atPtFixedStat[iS][iL] = new TGraphErrors(kCentLength-1);
      grB2atPtFixedStat[iS][iL]->SetName(Form("grB2atPtFixedStat_%c_%s",kLetter[iS],kFixedPtLabel[iL]));
      grB2atPtFixedSyst[iS][iL] = new TGraphErrors(kCentLength-1);
      grB2atPtFixedSyst[iS][iL]->SetName(Form("grB2atPtFixedSyst_%c_%s",kLetter[iS],kFixedPtLabel[iL]));
      //
      grB2atPtFixedStatMB[iS][iL] = new TGraphErrors(1);
      grB2atPtFixedStatMB[iS][iL]->SetName(Form("grB2atPtFixedStatMB_%c_%s",kLetter[iS],kFixedPtLabel[iL]));
      grB2atPtFixedSystMB[iS][iL] = new TGraphErrors(1);
      grB2atPtFixedSystMB[iS][iL]->SetName(Form("grB2atPtFixedSystMB_%c_%s",kLetter[iS],kFixedPtLabel[iL]));
    }
  }

  TH1F* hProtSpectraStat[kCentLength];
  TH1F* hProtSpectraSyst[kCentLength];
  TGraphErrors* grProtSpectraStat[kCentLength];
  TGraphErrors* grProtSpectraSyst[kCentLength];
  for(int iInput=0; iInput < kCentLength-1; iInput++) {
    /// Use an additional factor 2 to get the proton spectra from the p+pbar spectrum
    double totalW = 2. * (kMultInput[myCent[iInput][1] + 1] - kMultInput[myCent[iInput][0]]);

    for (int iC = myCent[iInput][0]; iC <= myCent[iInput][1]; ++iC) {
      double currentW = (kMultInput[iC + 1] - kMultInput[iC]) / totalW;
      auto stat = (TH1F*)input_file.Get(Form("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_%ito%i_stat",kMultInput[iC], kMultInput[iC+1]));
      Requires(stat,Form("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_%ito%i_stat",kMultInput[iC], kMultInput[iC+1]));
      auto syst = (TH1F*)input_file.Get(Form("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_%ito%i_syst",kMultInput[iC], kMultInput[iC+1]));
      Requires(syst,Form("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_%ito%i_syst",kMultInput[iC], kMultInput[iC+1]));
      stat->Scale(currentW);
      syst->Scale(currentW);
      if (iC == myCent[iInput][0]) {
        hProtSpectraStat[iInput] = (TH1F*)stat->Clone(Form("hProton_stat_%d",iInput));
        hProtSpectraSyst[iInput] = (TH1F*)syst->Clone(Form("hProton_syst_%d",iInput));
      } else {
        hProtSpectraStat[iInput]->Add(stat);
        hProtSpectraSyst[iInput]->Add(syst);
      }
    }

    MakeItInvariant(hProtSpectraStat[iInput]);
    MakeItInvariant(hProtSpectraSyst[iInput]);
    grProtSpectraStat[iInput] = new TGraphErrors(hProtSpectraStat[iInput]);
    grProtSpectraSyst[iInput] = new TGraphErrors(hProtSpectraSyst[iInput]);
  }

  // MB
  if(kUseIntegratedForMB){
    auto stat = (TH1F*)input_file.Get("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_0to100_stat");
    Requires(stat,"hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_0to100_stat");
    auto syst = (TH1F*)input_file.Get("hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_0to100_syst");
    Requires(syst,"hCombinedTPCTOFTOFonlyrTPCKinksITSsa_Pr_0to100_syst");
    stat->Scale(0.5);
    syst->Scale(0.5);
    hProtSpectraStat[9] = (TH1F*)stat->Clone("hProton_stat_9");
    hProtSpectraSyst[9] = (TH1F*)syst->Clone("hProton_syst_9");
    MakeItInvariant(hProtSpectraStat[9]);
    MakeItInvariant(hProtSpectraSyst[9]);
    grProtSpectraStat[9] = new TGraphErrors(hProtSpectraStat[9]);
    grProtSpectraSyst[9] = new TGraphErrors(hProtSpectraSyst[9]);
  }
  else{
    TFile input_file_MB(Form("%spp13TeV.mb.fullpT.INEL.FINAL-2019-05-30.root",kBaseOutputDir.data()));
    auto stat = (TH1F*)input_file_MB.Get("hstat_pp13_mb_proton_sum");
    Requires(stat,"hstat_pp13_mb_proton_sum");
    auto syst = (TH1F*)input_file_MB.Get("hsys_pp13_mb_proton_sum");
    Requires(syst,"hsys_pp13_mb_proton_sum");
    stat->Scale(0.5);
    syst->Scale(0.5);
    hProtSpectraStat[9] = (TH1F*)stat->Clone("hProton_stat_9");
    hProtSpectraSyst[9] = (TH1F*)syst->Clone("hProton_syst_9");
    //MakeItInvariant(hProtSpectraStat[9]);
    //MakeItInvariant(hProtSpectraSyst[9]);
    grProtSpectraStat[9] = new TGraphErrors(hProtSpectraStat[9]);
    grProtSpectraSyst[9] = new TGraphErrors(hProtSpectraSyst[9]);
  }
  //

  TH1F* hDeutSpectraStat[2][kCentLength];
  TH1F* hDeutSpectraSyst[2][kCentLength];

  TGraphErrors* grDeutSpectraStat[2][kCentLength];
  TGraphErrors* grDeutSpectraSyst[2][kCentLength];

  TGraphErrors* grB2PtStat[2][kCentLength];
  TGraphErrors* grB2PtSyst[2][kCentLength];
  TGraphErrors* grB2PtTot[2][kCentLength];
  TGraphErrors* grB2PtTot2[2][kCentLength];

  TGraphErrors* grB2PtStatClone[2][kCentLength];
  TGraphErrors* grB2PtSystClone[2][kCentLength];

  for(int iC=0; iC<kCentLength; iC++){
    //deuteron spectra: normalising by pt
    for(int iS = 0; iS<2; iS++){
      hDeutSpectraStat[iS][iC] = (TH1F*)deuteron_file.Get(Form("%sdeuterons/%d/stat",kPrefix[iS],iC));
      Requires(hDeutSpectraStat[iS][iC],Form("%sdeuterons/%d/stat",kPrefix[iS],iC));
      hDeutSpectraStat[iS][iC]->SetDirectory(0);
      MakeItInvariant(hDeutSpectraStat[iS][iC]);
      hDeutSpectraSyst[iS][iC]=(TH1F*)deuteron_file.Get(Form("%sdeuterons/%d/syst",kPrefix[iS],iC));
      Requires(hDeutSpectraSyst[iS][iC],Form("%sdeuterons/%d/syst",kPrefix[iS],iC));
      hDeutSpectraSyst[iS][iC]->SetDirectory(0);
      MakeItInvariant(hDeutSpectraSyst[iS][iC]);
      grDeutSpectraStat[iS][iC] = new TGraphErrors(hDeutSpectraStat[iS][iC]);
      grDeutSpectraSyst[iS][iC] = new TGraphErrors(hDeutSpectraSyst[iS][iC]);
      //
      grB2PtStat[iS][iC] = GetBAPt(grProtSpectraStat[iC],grDeutSpectraStat[iS][iC],Form("grStat_%c_%d",kLetter[iS],iC),false);
      grB2PtSyst[iS][iC] = GetBAPt(grProtSpectraSyst[iC],grDeutSpectraSyst[iS][iC],Form("grSyst_%c_%d",kLetter[iS],iC));
      // //
      // // Normalising with event loss
      // //
      // for (int i=0;i<grB2PtStat[iS][iC]->GetN();i++){
      //   grB2PtStat[iS][iC]->GetY()[i] /= kNormalisation[iC];
      //   grB2PtStat[iS][iC]->GetEY()[i] /= kNormalisation[iC];
      //   grB2PtSyst[iS][iC]->GetY()[i] /= kNormalisation[iC];
      //   grB2PtSyst[iS][iC]->GetEY()[i] /= kNormalisation[iC];
      // }
      // //
      grB2PtTot[iS][iC] = (TGraphErrors*)grB2PtStat[iS][iC]->Clone(Form("grTot_%c_%d_",kLetter[iS],iC));
      grB2PtTot2[iS][iC] = (TGraphErrors*)grB2PtStat[iS][iC]->Clone(Form("grTot_%c_%d_",kLetter[iS],iC));
      grB2PtStatClone[iS][iC] = (TGraphErrors*)grB2PtStat[iS][iC]->Clone(Form("grStat_%c_%d_clone",kLetter[iS],iC));
      grB2PtSystClone[iS][iC] = (TGraphErrors*)grB2PtSyst[iS][iC]->Clone(Form("grSyst_%c_%d_clone",kLetter[iS],iC));
      //
      grB2PtStat[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/#it{A} (GeV/#it{c})");
      grB2PtStat[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtStat[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtStat[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtStat[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtStat[iS][iC]->SetMarkerStyle(21);
      grB2PtStat[iS][iC]->SetFillStyle(0);
      //
      grB2PtSyst[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/#it{A} (GeV/#it{c})");
      grB2PtSyst[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtSyst[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtSyst[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtSyst[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtSyst[iS][iC]->SetMarkerStyle(21);
      grB2PtSyst[iS][iC]->SetFillStyle(0);
      //
      grB2PtTot[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/#it{A} (GeV/#it{c})");
      grB2PtTot[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtTot[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtTot[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtTot[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtTot[iS][iC]->SetMarkerStyle(21);
      grB2PtTot[iS][iC]->SetFillStyle(0);
      //
      grB2PtTot2[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/#it{A} (GeV/#it{c})");
      grB2PtTot2[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtTot2[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtTot2[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtTot2[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtTot2[iS][iC]->SetMarkerStyle(21);
      grB2PtTot2[iS][iC]->SetFillStyle(0);
      //
      grB2PtStatClone[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/#it{A} (GeV/#it{c})");
      grB2PtStatClone[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtStatClone[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtStatClone[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtStatClone[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtStatClone[iS][iC]->SetMarkerStyle(21);
      grB2PtStatClone[iS][iC]->SetFillStyle(0);
      //
      grB2PtSystClone[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/#it{A} (GeV/#it{c})");
      grB2PtSystClone[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtSystClone[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtSystClone[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtSystClone[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtSystClone[iS][iC]->SetMarkerStyle(21);
      grB2PtSystClone[iS][iC]->SetFillStyle(0);

      if(iC==kCentLength-1){
        //Fixed pT = 0.75 Gev/c
        double stat_tmp, syst_tmp, x_tmp, dummy, val_tmp;
        x_tmp = grB2PtStat[iS][iC]->GetX()[6];
        val_tmp = grB2PtStat[iS][iC]->GetY()[6];
        stat_tmp = grB2PtStat[iS][iC]->GetErrorY(6);
        syst_tmp = grB2PtSyst[iS][iC]->GetErrorY(6);
        //
        printf("*********************************\n");
        printf("expected: 0.75  real: %f\n",x_tmp);
        printf("*********************************\n");
        grB2atPtFixedStatMB[iS][0]->SetPoint(0,dNdEta[iC],val_tmp);
        grB2atPtFixedSystMB[iS][0]->SetPoint(0,dNdEta[iC],val_tmp);
        grB2atPtFixedStatMB[iS][0]->SetPointError(0,0,stat_tmp);
        grB2atPtFixedSystMB[iS][0]->SetPointError(0,dNdEtaErr[iC],syst_tmp);
        // Fixed pT = 0.55 GeV/c
        double stat1_tmp, syst1_tmp, x_tmp1, dummy1, val1_tmp;
        x_tmp = grB2PtStat[iS][iC]->GetX()[3];
        val_tmp = grB2PtStat[iS][iC]->GetY()[3];//5,6,10
        x_tmp1  = grB2PtStat[iS][iC]->GetX()[4];
        val1_tmp = grB2PtStat[iS][iC]->GetY()[4];
        //
        stat_tmp = grB2PtStat[iS][iC]->GetErrorY(3);
        syst_tmp = grB2PtSyst[iS][iC]->GetErrorY(3);
        stat1_tmp = grB2PtStat[iS][iC]->GetErrorY(4);
        syst1_tmp = grB2PtSyst[iS][iC]->GetErrorY(4);
        //
        grB2atPtFixedStatMB[iS][1]->SetPoint(iC,dNdEta[iC],0.5*val_tmp+0.5*val1_tmp);
        grB2atPtFixedSystMB[iS][1]->SetPoint(iC,dNdEta[iC],0.5*val_tmp+0.5*val1_tmp);
        grB2atPtFixedStatMB[iS][1]->SetPointError(iC,0,TMath::Sqrt(stat_tmp*stat_tmp+stat1_tmp*stat1_tmp)/2);
        grB2atPtFixedSystMB[iS][1]->SetPointError(iC,dNdEtaErr[iC],TMath::Sqrt(syst_tmp*syst_tmp+syst1_tmp*syst1_tmp)/2);

        printf("*********************************\n");
        printf("expected: 0.55  real: %f\n",(x_tmp+x_tmp1)/2);
        printf("*********************************\n");

        if(iC>=kCentLength-3) continue;
        //Fixed pT = 1.05 GeV/c
        x_tmp = grB2PtStat[iS][iC]->GetX()[9];
        val_tmp = grB2PtStat[iS][iC]->GetY()[9];
        stat_tmp = grB2PtStat[iS][iC]->GetErrorY(9);
        syst_tmp = grB2PtSyst[iS][iC]->GetErrorY(9);
        // printf("%s\n",kFixedPtLabel[1]);
        // printf("x_val: %f\n", x_tmp);
        // printf("val:  %f stat: %f syst %f dNdEta:  %f dNdEtaErr: %f\n", val_tmp, stat_tmp, syst_tmp,dNdEta_vec[iC], dNdEtaErr_vec[iC]);
        grB2atPtFixedStatMB[iS][2]->SetPoint(iC,dNdEta[iC],val_tmp);
        grB2atPtFixedSystMB[iS][2]->SetPoint(iC,dNdEta[iC],val_tmp);
        grB2atPtFixedStatMB[iS][2]->SetPointError(iC,0,stat_tmp);
        grB2atPtFixedSystMB[iS][2]->SetPointError(iC,dNdEtaErr[iC],syst_tmp);
        printf("*********************************\n");
        printf("expected: 1.05  real: %f\n",x_tmp);
        printf("*********************************\n");
        continue;
      }
      //Fixed pT = 0.75 Gev/c
      double stat_tmp, syst_tmp, x_tmp, dummy, val_tmp;
      x_tmp = grB2PtStat[iS][iC]->GetX()[6];
      val_tmp = grB2PtStat[iS][iC]->GetY()[6];
      stat_tmp = grB2PtStat[iS][iC]->GetErrorY(6);
      syst_tmp = grB2PtSyst[iS][iC]->GetErrorY(6);
      printf("*********************************\n");
      printf("expected: 0.75  real: %f\n",x_tmp);
      printf("*********************************\n");
      //printf("%s\n",kFixedPtLabel[0]);
      //printf("x_val: %f\n", x_tmp);
      //printf("val:  %f stat: %f syst %f dNdEta:  %f dNdEtaErr: %f\n", val_tmp, stat_tmp, syst_tmp,dNdEta_vec[iC], dNdEtaErr_vec[iC]);
      grB2atPtFixedStat[iS][0]->SetPoint(iC,dNdEta[iC],val_tmp);
      grB2atPtFixedSyst[iS][0]->SetPoint(iC,dNdEta[iC],val_tmp);
      grB2atPtFixedStat[iS][0]->SetPointError(iC,0,stat_tmp);
      grB2atPtFixedSyst[iS][0]->SetPointError(iC,dNdEtaErr[iC],syst_tmp);

      // Fixed pT = 0.55 GeV/c
      double stat1_tmp, syst1_tmp, x_tmp1, dummy1, val1_tmp;

      //double stat1_tmp, syst1_tmp, x_tmp1, dummy1;
      //grB2PtStat[iS][iC]->GetPoint(6,x_tmp,val_tmp);
      x_tmp = grB2PtStat[iS][iC]->GetX()[3];
      val_tmp = grB2PtStat[iS][iC]->GetY()[3];//5,6,10
      x_tmp1  = grB2PtStat[iS][iC]->GetX()[4];
      val1_tmp = grB2PtStat[iS][iC]->GetY()[4];

      stat_tmp = grB2PtStat[iS][iC]->GetErrorY(3);
      syst_tmp = grB2PtSyst[iS][iC]->GetErrorY(3);
      stat1_tmp = grB2PtStat[iS][iC]->GetErrorY(4);
      syst1_tmp = grB2PtSyst[iS][iC]->GetErrorY(4);
      // printf("%s\n",kFixedPtLabel[2]);
      // printf("x_val: %f\n", x_tmp);
      // printf("val:  %f stat: %f syst %f dNdEta:  %f dNdEtaErr: %f\n", val_tmp, stat_tmp, syst_tmp, dNdEta_vec[iC], dNdEtaErr_vec[iC]);
      // printf("x_val: %f\n", x_tmp1);
      // printf("val:  %f stat: %f syst %f dNdEta:  %f dNdEtaErr: %f\n", val1_tmp, stat1_tmp, syst1_tmp, dNdEta_vec[iC], dNdEtaErr_vec[iC]);
      grB2atPtFixedStat[iS][1]->SetPoint(iC,dNdEta[iC],0.5*val_tmp+0.5*val1_tmp);
      grB2atPtFixedSyst[iS][1]->SetPoint(iC,dNdEta[iC],0.5*val_tmp+0.5*val1_tmp);
      grB2atPtFixedStat[iS][1]->SetPointError(iC,0,TMath::Sqrt(stat_tmp*stat_tmp+stat1_tmp*stat1_tmp)/2);
      grB2atPtFixedSyst[iS][1]->SetPointError(iC,dNdEtaErr[iC],TMath::Sqrt(syst_tmp*syst_tmp+syst1_tmp*syst1_tmp)/2);

      printf("*********************************\n");
      printf("expected: 0.55  real: %f\n",(x_tmp+x_tmp1)/2);
      printf("*********************************\n");

      if(iC>=kCentLength-3) continue;
      //Fixed pT = 1.05 GeV/c
      x_tmp = grB2PtStat[iS][iC]->GetX()[9];
      val_tmp = grB2PtStat[iS][iC]->GetY()[9];
      stat_tmp = grB2PtStat[iS][iC]->GetErrorY(9);
      syst_tmp = grB2PtSyst[iS][iC]->GetErrorY(9);
      // printf("%s\n",kFixedPtLabel[1]);
      // printf("x_val: %f\n", x_tmp);
      // printf("val:  %f stat: %f syst %f dNdEta:  %f dNdEtaErr: %f\n", val_tmp, stat_tmp, syst_tmp,dNdEta_vec[iC], dNdEtaErr_vec[iC]);
      grB2atPtFixedStat[iS][2]->SetPoint(iC,dNdEta[iC],val_tmp);
      grB2atPtFixedSyst[iS][2]->SetPoint(iC,dNdEta[iC],val_tmp);
      grB2atPtFixedStat[iS][2]->SetPointError(iC,0,stat_tmp);
      grB2atPtFixedSyst[iS][2]->SetPointError(iC,dNdEtaErr[iC],syst_tmp);
      printf("*********************************\n");
      printf("expected: 1.05  real: %f\n",x_tmp);
      printf("*********************************\n");
    }
  }

  for(int iS=0; iS<2; iS++){
    for(int iC=0; iC<kCentLength-1; iC++){
      for(int iBin=0; iBin<grB2PtStat[iS][iC]->GetN(); iBin++){
        float x_val = grB2PtStat[iS][iC]->GetX()[iBin];
        float x_err = grB2PtStat[iS][iC]->GetErrorX(iBin);
        float val = grB2PtStat[iS][iC]->GetY()[iBin];
        float stat_err = grB2PtStat[iS][iC]->GetErrorY(iBin);
        float syst_err = grB2PtSyst[iS][iC]->GetErrorY(iBin);
        float tot_err2 = stat_err* stat_err + syst_err*syst_err;
        float tot_err = TMath::Sqrt(tot_err2);
        grB2PtTot[iS][iC]->SetPoint(iBin,x_val,val);
        grB2PtTot[iS][iC]->SetPointError(iBin,x_err,tot_err);
        grB2PtTot2[iS][iC]->SetPoint(iBin,x_val,val);
        grB2PtTot2[iS][iC]->SetPointError(iBin,x_err,tot_err);
      }
      grB2PtTot[iS][iC]->Fit("pol0");
      grB2PtTot2[iS][iC]->Fit("pol1");
    }
  }

  // fit diversi

  float meanB2[kCentLength-1][2] = {{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.},{0.,0.}};
  float chi2[2] = {0.,0.};
  int vNDF[2] = {0,0};

  for (int iS=0; iS<2; iS++){
    //computing mean B2 for each multiplicity class
    for(int iC=0; iC<kCentLength-1; iC++){
      float val = 0.;
      float stat_err = 0.;
      float syst_err = 0.;
      float tot_err2 = 0.;
      float cumulative_err2 = 0.;
      for(int iBin=0; iBin<grB2PtStat[iS][iC]->GetN(); iBin++){
        val = grB2PtStat[iS][iC]->GetY()[iBin];
        stat_err = grB2PtStat[iS][iC]->GetErrorY(iBin);
        syst_err = grB2PtStat[iS][iC]->GetErrorY(iBin);
        tot_err2 = stat_err* stat_err + syst_err*syst_err;
        meanB2[iC][iS] += val/tot_err2;
        cumulative_err2 += 1/tot_err2;
        vNDF[iS]++;
      }
      meanB2[iC][iS] /= cumulative_err2;
    }
    //computing global chi2
    for(int iC=0; iC<kCentLength-1; iC++){
      float val = 0.;
      float stat_err = 0.;
      float syst_err = 0.;
      float tot_err2 = 0.;
      float cumulative_err2 = 0.;
      for(int iBin=0; iBin<grB2PtStat[iS][iC]->GetN(); iBin++){
        val = grB2PtStat[iS][iC]->GetY()[iBin];
        stat_err = grB2PtStat[iS][iC]->GetErrorY(iBin);
        syst_err = grB2PtStat[iS][iC]->GetErrorY(iBin);
        tot_err2 = stat_err* stat_err + syst_err*syst_err;
        chi2[iS] += (val-meanB2[iC][iS])*(val-meanB2[iC][iS])/tot_err2;
      }
      vNDF[iS]--;
    }
  }

  std::cout << "*************  LOCAL   ******************" << std::endl;
  std::cout << kNames[0] << " chi2 / NDF : " << chi2[0] << " / " << vNDF[0] << " = " << chi2[0] / vNDF[0] <<
  " sigma : " << ROOT::Math::gaussian_quantile_c(TMath::Prob(chi2[0],vNDF[0]),1) << std::endl;
  std::cout << kNames[1] << " chi2 / NDF : " << chi2[1] << " / " << vNDF[1] << " = " << chi2[1] / vNDF[1] <<
  " sigma : " << ROOT::Math::gaussian_quantile_c(TMath::Prob(chi2[1],vNDF[1]),1)<< std::endl;
  std::cout << "*******************************" << std::endl;


  // simultaneo

  float meanB2sim[2] = {0.,0.};
  float chi2sim[2] = {0.,0.};
  int vNDFsim[2] = {0,0};

  for (int iS=0; iS<2; iS++){
    //computing mean B2 for each multiplicity class
    float val = 0.;
    float stat_err = 0.;
    float syst_err = 0.;
    float tot_err2 = 0.;
    float cumulative_err2 = 0.;
    for(int iC=0; iC<kCentLength-1; iC++){
      for(int iBin=0; iBin<grB2PtStat[iS][iC]->GetN(); iBin++){
        val = grB2PtStat[iS][iC]->GetY()[iBin];
        stat_err = grB2PtStat[iS][iC]->GetErrorY(iBin);
        syst_err = grB2PtSyst[iS][iC]->GetErrorY(iBin);
        tot_err2 = stat_err* stat_err + syst_err*syst_err;
        meanB2sim[iS] += val/tot_err2;
        cumulative_err2 += 1/tot_err2;
        vNDFsim[iS]++;
      }
    }
    meanB2sim[iS] /= cumulative_err2;
    //computing global chi2
    for(int iC=0; iC<kCentLength-1; iC++){
      for(int iBin=0; iBin<grB2PtStat[iS][iC]->GetN(); iBin++){
        val = grB2PtStat[iS][iC]->GetY()[iBin];
        stat_err = grB2PtStat[iS][iC]->GetErrorY(iBin);
        syst_err = grB2PtSyst[iS][iC]->GetErrorY(iBin);
        tot_err2 = stat_err* stat_err + syst_err*syst_err;
        chi2sim[iS] += (val-meanB2sim[iS])*(val-meanB2sim[iS])/tot_err2;
      }
    }
    vNDFsim[iS]--;
  }

  std::cout << "**************  GLOBAL   *****************" << std::endl;
  std::cout << kNames[0] << " chi2sim / NDFsim : " << chi2sim[0] << " / " << vNDFsim[0] << " = " << chi2sim[0] / vNDFsim[0] <<
  "sigma : " << ROOT::Math::gaussian_quantile_c(TMath::Prob(chi2sim[0],vNDFsim[0]),1)<< std::endl;
  std::cout << kNames[1] << " chi2sim / NDFsim : " << chi2sim[1] << " / " << vNDFsim[1] << " = " << chi2sim[1] / vNDFsim[1] <<
  "sigma : " << ROOT::Math::gaussian_quantile_c(TMath::Prob(chi2sim[1],vNDFsim[1]),1)<< std::endl;
  std::cout << "*******************************" << std::endl;

  TCanvas* cFit[2];
  TCanvas* cB2[2];
  TCanvas* cB2notScaled[2];
  for(int iS=0; iS<2; iS++){
    TDirectory* s_dir = output_file.mkdir(kNames[iS].data());
    TLatex text_one;
    s_dir->cd();
    //
    cFit[iS] = new TCanvas(Form("B2_Fit_%s",kNames[iS].data()),Form("B2_Fit_%s",kNames[iS].data()),800,600);
    cFit[iS]->Divide(3,3);
    for (int iC = 0; iC < 9; ++iC) {
      cFit[iS]->cd(iC+1);
      text_one.DrawLatex(0.5,1.7,Form("#bf{%s}",kRomanLabels[iC]));
      grB2PtTot[iS][iC]->Draw();
    }
    s_dir->cd();
    cFit[iS]->Write();
    //
    cB2[iS] = new TCanvas(Form("B2_%s",kNames[iS].data()),Form("B2_%s",kNames[iS].data()),800,600);
    cB2[iS]->DrawFrame(
      0.2,
      1e-3,
      2.6,
      80,
      ";#it{p}_{T}/#it{A} (GeV/#it{c});#it{B}_{2} (GeV^{2}/#it{c}^{3})"
    );
    cB2[iS]->SetLeftMargin(0.1425591);
    cB2[iS]->SetRightMargin(0.027121);
    cB2[iS]->SetTopMargin(0.06053269);
    cB2[iS]->SetBottomMargin(0.1598063);
    cB2[iS]->SetLogy();
    TLatex text;
    text.SetTextFont(63);
    text.SetTextSize(22);
    text.DrawLatex(0.22,97.4,Form("#bf{%s, pp, #sqrt{#it{s}} = 13 TeV}",kNamePlot[iS].data()));
    text.DrawLatex(1.03464,0.00267768,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}#GT = 26.02}}",plotting::kSpectraColors[0]));
    text.DrawLatex(0.324014,27.6452,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}#GT = 2.55}}",plotting::kSpectraColors[8]));
    TLegend leg(0.690476,0.393043,0.889724,0.803478);
    leg.SetBorderSize(0);
    leg.SetHeader("V0M Multiplicity Classes");
    TLegendEntry *header = (TLegendEntry*)leg.GetListOfPrimitives()->First();
    header->SetTextAlign(13);
    //leg.SetNColumns(2);
    leg.SetTextSize(0.03);
    // MB
    cB2notScaled[iS] = new TCanvas(Form("B2notScaled%s",kNames[iS].data()),Form("B2_%s",kNames[iS].data()),800,600);
    cB2notScaled[iS]->DrawFrame(
      0.2,
      0.,
      2.6,
      0.1,
      ";#it{p}_{T}/#it{A} (GeV/#it{c});#it{B}_{2} (GeV^{2}/#it{c}^{3})"
    );
    cB2notScaled[iS]->SetLeftMargin(0.1425591);
    cB2notScaled[iS]->SetRightMargin(0.027121);
    cB2notScaled[iS]->SetTopMargin(0.06053269);
    cB2notScaled[iS]->SetBottomMargin(0.1598063);
    cB2notScaled[iS]->cd();
    TLatex text_nS;
    text_nS.SetTextFont(63);
    text_nS.SetTextSize(22);
    text_nS.DrawLatex(0.22,0.101,Form("#bf{%s, pp, #sqrt{#it{s}} = 13 TeV}",kNames[iS].data()));
    TLegend leg_nS(0.25188,0.518261,0.642857,0.834783);
    leg_nS.SetBorderSize(0);
    leg_nS.SetHeader("V0M Multiplicity Classes");
    TLegendEntry *header_nS = (TLegendEntry*)leg_nS.GetListOfPrimitives()->First();
    header_nS->SetTextAlign(13);
    leg_nS.SetNColumns(2);
    leg_nS.SetTextSize(0.03);
    for(int iC=0; iC<kCentLength; iC++){
      TDirectory *c_dir = s_dir->mkdir(std::to_string(iC).data());
      c_dir->cd();
      //write to file
      grProtSpectraStat[iC]->Write(Form("grProtStat_%c_%d",kLetter[iS],iC));
      grProtSpectraSyst[iC]->Write(Form("grProtSyst_%c_%d",kLetter[iS],iC));
      grDeutSpectraStat[iS][iC]->Write(Form("grDeutStat_%c_%d",kLetter[iS],iC));
      grDeutSpectraSyst[iS][iC]->Write(Form("grDeutSyst_%c_%d",kLetter[iS],iC));
      grB2PtStat[iS][iC]->Write(Form("grB2Stat_%c_%d",kLetter[iS],iC));
      grB2PtSyst[iS][iC]->Write(Form("grB2Syst_%c_%d",kLetter[iS],iC));
      grB2PtTot[iS][iC]->Write(Form("grB2Tot_%c_%d",kLetter[iS],iC));
      grB2PtTot2[iS][iC]->Write(Form("grB2Tot2_%c_%d",kLetter[iS],iC));
      cB2[iS]->cd();
      for (int i=0;i<grB2PtStatClone[iS][iC]->GetN();i++){
        grB2PtStatClone[iS][iC]->GetY()[i] *= kScaleFactorB2[iC];
        grB2PtStatClone[iS][iC]->GetEY()[i] *= kScaleFactorB2[iC];
        grB2PtSystClone[iS][iC]->GetY()[i] *= kScaleFactorB2[iC];
        grB2PtSystClone[iS][iC]->GetEY()[i] *= kScaleFactorB2[iC];
      }
      grB2PtStatClone[iS][iC]->Draw("samepz");
      grB2PtSystClone[iS][iC]->Draw("e2same");
      leg.AddEntry(grB2PtSystClone[iS][iC],Form("%s (#times %d)",kRomanLabels[iC],kScaleFactorB2[iC]),"pf");
      // not scaled
      cB2notScaled[iS]->cd();
      grB2PtStat[iS][iC]->Draw("samepz");
      grB2PtSyst[iS][iC]->Draw("e2same");
      leg_nS.AddEntry(grB2PtSystClone[iS][iC],Form("%s",kRomanLabels[iC]),"pf");
    }
    s_dir->cd();
    cB2[iS]->cd();
    leg.Draw();
    cB2[iS]->Write();
    cB2notScaled[iS]->cd();
    leg_nS.Draw();
    cB2notScaled[iS]->Write();
  }

  for(int iS=0; iS<2; iS++){
    for(int iL=0; iL<3; iL++){
      output_file.cd();
      grB2atPtFixedStat[iS][iL]->SetLineColor(kOrange-3);
      grB2atPtFixedStat[iS][iL]->SetMarkerColor(kOrange-3);
      grB2atPtFixedStat[iS][iL]->SetMarkerStyle(21);
      grB2atPtFixedStat[iS][iL]->SetFillStyle(0);
      grB2atPtFixedSyst[iS][iL]->SetLineColor(kOrange-3);
      grB2atPtFixedSyst[iS][iL]->SetMarkerColor(kOrange-3);
      grB2atPtFixedSyst[iS][iL]->SetMarkerStyle(21);
      grB2atPtFixedSyst[iS][iL]->SetFillStyle(0);
      grB2atPtFixedStat[iS][iL]->Write();
      grB2atPtFixedSyst[iS][iL]->Write();
      // MB
      grB2atPtFixedStatMB[iS][iL]->SetLineColor(kViolet-1);
      grB2atPtFixedStatMB[iS][iL]->SetMarkerColor(kViolet-1);
      grB2atPtFixedStatMB[iS][iL]->SetMarkerStyle(21);
      grB2atPtFixedStatMB[iS][iL]->SetFillStyle(0);
      grB2atPtFixedSystMB[iS][iL]->SetLineColor(kViolet-1);
      grB2atPtFixedSystMB[iS][iL]->SetMarkerColor(kViolet-1);
      grB2atPtFixedSystMB[iS][iL]->SetMarkerStyle(21);
      grB2atPtFixedSystMB[iS][iL]->SetFillStyle(0);
      grB2atPtFixedStatMB[iS][iL]->Write();
      grB2atPtFixedSystMB[iS][iL]->Write();
    }
  }

  TCanvas* cB2atPtFixed[2] = {new TCanvas(Form("cB2atPtFixed_%c",kLetter[0]),Form("cB2atPtFixed_%c",kLetter[0])),new TCanvas(Form("cB2atPtFixed_%c",kLetter[1]),Form("cB2atPtFixed_%c",kLetter[1]))};

  TFile old_b2_results(Form("%sB2vsCollsAll_pToA_105.root",kBaseOutputDir.data()));
  TCanvas* cB2atPtFixedInput = (TCanvas*)old_b2_results.Get("c1_n5");
  Requires(cB2atPtFixedInput,"old b2 at fixed pt results");
  TPad* pInput = (TPad*)cB2atPtFixedInput->GetPrimitive("c1_n5_1");
  for(int iS=0; iS<2; iS++){
    if(iS==1){
      pInput->cd();
      grB2atPtFixedStat[iS][2]->Draw("pz");
      grB2atPtFixedSyst[iS][2]->Draw("e2");
      pInput->Modified();
      pInput->Update();
      //pInput->GetListOfPrimitives()->ls();
      cB2atPtFixed[iS]->cd();
      pInput->Draw();
      output_file.cd();
      cB2atPtFixed[iS]->Write();
    }
  }
}

double GetErrorY(const TGraphErrors* gr, double x0){
//
// estimate error of gr(x0) with the closest point to x0
//
  const double kEpsilon  = 1.e-6;
  const double kUnknownError = 1.e+6;

  for(int i=0; i<gr->GetN(); ++i){
    double x = gr->GetX()[i];
    if( TMath::Abs(x-x0) < kEpsilon) return gr->GetErrorY(i);
    if( ((i == 0) && (x > x0)) || ((i == (gr->GetN()-1)) && (x < x0)) ){
      std::cout << "GetErrorY: " << x << "is out of bounds" << std::endl;
      return kUnknownError;
    }
    if( x > x0 ){
      std::cout << "GetErrorY: Interpolating error at " << x0 << std::endl;
      return (gr->GetErrorY(i)+gr->GetErrorY(i-1))/2;
    }
  }
  return 0;
}

void MakeItInvariant(TH1* d_stat) {
    for (int iB = 1; iB <= d_stat->GetNbinsX(); ++iB) {
      double norm = 1. / (d_stat->GetBinCenter(iB) * TMath::TwoPi());
      d_stat->SetBinContent(iB,d_stat->GetBinContent(iB) * norm);
      d_stat->SetBinError(iB,d_stat->GetBinError(iB) * norm);
    }
}

TGraphErrors* GetBAPt(const TGraphErrors* grPrtInvDYieldPt,
    const TGraphErrors* grNucInvDYieldPt, const TString& name, bool syst) {
  //
  // coalescence parameter
  //
  TGraphErrors* grBAPt = new TGraphErrors();
  grBAPt->SetName(name.Data());

  for(Int_t i=0, j=0; i < grNucInvDYieldPt->GetN(); ++i)
  {
    Double_t ptNuc, yNuc;

    grNucInvDYieldPt->GetPoint(i, ptNuc, yNuc);

    if(ptNuc/2 < 0.3) continue; // acceptance

    Double_t yPrt = grPrtInvDYieldPt->Eval(ptNuc/2); // interpolate

    if(yPrt == 0 || yNuc == 0 ) continue;

    Double_t bA = yNuc/TMath::Power(yPrt,2);

    // error
    Double_t ePrt = GetErrorY(grPrtInvDYieldPt, ptNuc/2);
    Double_t eNuc = grNucInvDYieldPt->GetErrorY(i);

    Double_t errPt = (syst) ? grNucInvDYieldPt->GetErrorX(i)/2 : 0.;
    Double_t errBA = bA*TMath::Sqrt(TMath::Power(eNuc/yNuc,2) + TMath::Power(2*ePrt/yPrt,2));

    grBAPt->SetPoint(j, ptNuc/2, bA);
    grBAPt->SetPointError(j++, errPt, errBA);
  }

  return grBAPt;
}
