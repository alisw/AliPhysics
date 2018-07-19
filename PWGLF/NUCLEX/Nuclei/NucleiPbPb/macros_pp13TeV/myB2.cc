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

#include <string>
using std::string;
#include <fstream>

const char* kPrefix[2] = {"","anti"};
const int kScaleFactor[10] = {256,128,64,32,16,8,4,2,1,1};

double GetErrorY(const TGraphErrors* gr, double x0);
void MakeItInvariant(TH1* d_stat);
TGraphErrors* GetBAPt(const TGraphErrors* grPrtInvDYieldPt,
    const TGraphErrors* grNucInvDYieldPt, const TString& name, bool syst=true);

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
  {9,9}  /// 70-100%
};
void myB2(){

  double dNdEta_tmp, dNdEtaErr_tmp, proton_yields_tmp, proton_yields_stat_tmp, proton_yields_syst_tmp;

  vector<double> dNdEta_vec, dNdEtaErr_vec, proton_yields_vec, proton_yields_stat_vec, proton_yields_syst_vec;

  std::ifstream protonFile(Form("%s/proton_yields.txt",kBaseOutputDir.data()));
  if(!protonFile.is_open()){
    printf("The file %s could not be opened\n", Form("%s/proton_yields.txt",kBaseOutputDir.data()));
    exit(1);
  }
  else{
    while(protonFile >> dNdEta_tmp >> dNdEtaErr_tmp >> proton_yields_tmp >> proton_yields_stat_tmp >> proton_yields_syst_tmp){
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

  TFile input_file(Form("%sFinal_combined_spectra_pp13TeV.root",kBaseOutputDir.data()));

  TGraphErrors* grB2atPtFixedStat[2];
  TGraphErrors* grB2atPtFixedSyst[2];
  for(int iS=0; iS<2; iS++){
    grB2atPtFixedStat[iS] = new TGraphErrors(kCentLength-1);
    grB2atPtFixedStat[iS]->SetName(Form("grB2atPtFixedStat_%c",kLetter[iS]));
    grB2atPtFixedSyst[iS] = new TGraphErrors(kCentLength-1);
    grB2atPtFixedSyst[iS]->SetName(Form("grB2atPtFixedSyst_%c",kLetter[iS]));
  }

  TH1F* hProtSpectraStat[kCentLength];
  TH1F* hProtSpectraSyst[kCentLength];
  TGraphErrors* grProtSpectraStat[kCentLength];
  TGraphErrors* grProtSpectraSyst[kCentLength];
  for(int iInput=0; iInput < kCentLength; iInput++) {
    /// Use an additional factor 2 to get the proton spectra from the p+pbar spectrum
    double totalW = 2. * (kMultInput[myCent[iInput][1] + 1] - kMultInput[myCent[iInput][0]]);

    for (int iC = myCent[iInput][0]; iC <= myCent[iInput][1]; ++iC) {
      double currentW = (kMultInput[iC + 1] - kMultInput[iC]) / totalW;
      auto stat = (TH1F*)input_file.Get(Form("hCombinedTPCTOFTOF_ITSsa_Pr_%ito%i_stat",kMultInput[iC], kMultInput[iC+1]));
      Requires(stat,Form("hCombinedTPCTOFTOF_ITSsa_Pr_%ito%i_stat",kMultInput[iC], kMultInput[iC+1]));
      auto syst = (TH1F*)input_file.Get(Form("hCombinedTPCTOFTOF_ITSsa_Pr_%ito%i_syst",kMultInput[iC], kMultInput[iC+1]));
      Requires(syst,Form("hCombinedTPCTOFTOF_ITSsa_Pr_%ito%i_syst",kMultInput[iC], kMultInput[iC+1]));
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

  TH1F* hDeutSpectraStat[2][kCentLength];
  TH1F* hDeutSpectraSyst[2][kCentLength];

  TGraphErrors* grDeutSpectraStat[2][kCentLength];
  TGraphErrors* grDeutSpectraSyst[2][kCentLength];

  TGraphErrors* grB2PtStat[2][kCentLength];
  TGraphErrors* grB2PtSyst[2][kCentLength];

  TGraphErrors* grB2PtStatClone[2][kCentLength];
  TGraphErrors* grB2PtSystClone[2][kCentLength];

  for(int iC=0; iC<kCentLength; iC++){
    //deuteron spectra: normalising by pt
    for(int iS = 0; iS<2; iS++){
      hDeutSpectraStat[iS][iC] = (TH1F*)deuteron_file.Get(Form("%sdeuterons/%d/stat_all",kPrefix[iS],iC));
      hDeutSpectraStat[iS][iC]->SetDirectory(0);
      MakeItInvariant(hDeutSpectraStat[iS][iC]);
      hDeutSpectraSyst[iS][iC]=(TH1F*)deuteron_file.Get(Form("%sdeuterons/%d/syst_all",kPrefix[iS],iC));
      hDeutSpectraSyst[iS][iC]->SetDirectory(0);
      MakeItInvariant(hDeutSpectraSyst[iS][iC]);
      grDeutSpectraStat[iS][iC] = new TGraphErrors(hDeutSpectraStat[iS][iC]);
      grDeutSpectraSyst[iS][iC] = new TGraphErrors(hDeutSpectraSyst[iS][iC]);
      //
      grB2PtStat[iS][iC] = GetBAPt(grProtSpectraStat[iC],grDeutSpectraStat[iS][iC],Form("grStat_%c_%d",kLetter[iS],iC),false);
      grB2PtStatClone[iS][iC] = (TGraphErrors*)grB2PtStat[iS][iC]->Clone(Form("grStat_%c_%d_clone",kLetter[iS],iC));
      grB2PtSyst[iS][iC] = GetBAPt(grProtSpectraSyst[iC],grDeutSpectraSyst[iS][iC],Form("grSyst_%c_%d",kLetter[iS],iC));
      grB2PtSystClone[iS][iC] = (TGraphErrors*)grB2PtSyst[iS][iC]->Clone(Form("grSyst_%c_%d_clone",kLetter[iS],iC));
      //
      grB2PtStat[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
      grB2PtStat[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtStat[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtStat[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtStat[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtStat[iS][iC]->SetMarkerStyle(21);
      grB2PtStat[iS][iC]->SetFillStyle(0);
      //
      grB2PtSyst[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
      grB2PtSyst[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtSyst[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtSyst[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtSyst[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtSyst[iS][iC]->SetMarkerStyle(21);
      grB2PtSyst[iS][iC]->SetFillStyle(0);
      //
      grB2PtStatClone[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
      grB2PtStatClone[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtStatClone[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtStatClone[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtStatClone[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtStatClone[iS][iC]->SetMarkerStyle(21);
      grB2PtStatClone[iS][iC]->SetFillStyle(0);
      //
      grB2PtSystClone[iS][iC]->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
      grB2PtSystClone[iS][iC]->GetXaxis()->SetRangeUser(0.1,2.);
      grB2PtSystClone[iS][iC]->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
      grB2PtSystClone[iS][iC]->SetLineColor(kSpectraColors[iC]);
      grB2PtSystClone[iS][iC]->SetMarkerColor(kSpectraColors[iC]);
      grB2PtSystClone[iS][iC]->SetMarkerStyle(21);
      grB2PtSystClone[iS][iC]->SetFillStyle(0);

      if(iC==kCentLength-1) continue;
      double stat_tmp, syst_tmp, x_tmp, dummy;
      double stat1_tmp, syst1_tmp, x_tmp1, dummy1;
      //grB2PtStat[iS][iC]->GetPoint(6,x_tmp,val_tmp);
      double val_tmp = grB2PtStat[iS][iC]->GetY()[5];
      double val1_tmp = grB2PtStat[iS][iC]->GetY()[4];

      stat_tmp = grB2PtStat[iS][iC]->GetErrorY(5);
      syst_tmp = grB2PtSyst[iS][iC]->GetErrorY(5);
      stat1_tmp = grB2PtStat[iS][iC]->GetErrorY(4);
      syst1_tmp = grB2PtSyst[iS][iC]->GetErrorY(4);
      printf("x_val: %f\n", x_tmp);
      printf("val:  %f stat: %f syst %f dNdEta:  %f dNdEtaErr: %f\n", val_tmp, stat_tmp, syst_tmp, dNdEta_vec[iC], dNdEtaErr_vec[iC]);
      grB2atPtFixedStat[iS]->SetPoint(iC,dNdEta[iC],0.5*val_tmp+0.5*val1_tmp);
      grB2atPtFixedSyst[iS]->SetPoint(iC,dNdEta[iC],0.5*val_tmp+0.5*val1_tmp);
      grB2atPtFixedStat[iS]->SetPointError(iC,0,TMath::Sqrt(stat_tmp*stat_tmp+stat1_tmp*stat1_tmp)/2);
      grB2atPtFixedSyst[iS]->SetPointError(iC,dNdEtaErr[iC],TMath::Sqrt(syst_tmp*syst_tmp+syst1_tmp*syst1_tmp)/2);
    }
  }

  TCanvas* cB2[2];
  for(int iS=0; iS<2; iS++){
    TDirectory* s_dir = output_file.mkdir(kNames[iS].data());
    s_dir->cd();
    cB2[iS] = new TCanvas(kNames[iS].data(),kNames[iS].data());
    cB2[iS]->DrawFrame(
      0.,
      0,
      2.,
      0.02,
      ";#it{p}_{T}/A (GeV/#it{c});#it{B}_{2} (GeV^{2}/#it{c}^{3})"
    );
    cB2[iS]->SetLeftMargin(0.1425591);
    cB2[iS]->SetRightMargin(0.027121);
    cB2[iS]->SetTopMargin(0.06053269);
    cB2[iS]->SetBottomMargin(0.1598063);
    //cB2[iS]->SetLogy();
    TLegend leg(0.64,0.20,0.94,0.38);
    leg.SetBorderSize(0);
    leg.SetNColumns(2);
    for(int iC=0; iC<kCentLength; iC++){
      TDirectory *c_dir = s_dir->mkdir(to_string(iC).data());
      c_dir->cd();
      //write to file
      grProtSpectraStat[iC]->Write(Form("grProtStat_%c_%d",kLetter[iS],iC));
      grProtSpectraSyst[iC]->Write(Form("grProtSyst_%c_%d",kLetter[iS],iC));
      grDeutSpectraStat[iS][iC]->Write(Form("grDeutStat_%c_%d",kLetter[iS],iC));
      grDeutSpectraSyst[iS][iC]->Write(Form("grDeutSyst_%c_%d",kLetter[iS],iC));
      grB2PtStat[iS][iC]->Write(Form("grB2Stat_%c_%d",kLetter[iS],iC));
      grB2PtSyst[iS][iC]->Write(Form("grB2Syst_%c_%d",kLetter[iS],iC));
      if(iC==kCentLength-1) continue;
      cB2[iS]->cd();
      // for (int i=0;i<grB2PtStatClone[iS][iC]->GetN();i++){
      //   // double x_tmp, val_tmp;
      //   // grB2PtStat[iS][iC]->GetPoint(i,x_tmp,val_tmp);
      //   // double errx_tmp = grB2PtStat[iS][iC]->GetEX()[i];
      //   // double stat_tmp = grB2PtStat[iS][iC]->GetEY()[i];
      //   // double syst_tmp = grB2PtSyst[iS][iC]->GetEY()[i];
      //   // grB2PtStatClone[iS][iC]->SetPoint(i,x_tmp,val_tmp*kScaleFactor[iC]);
      //   // grB2PtStatClone[iS][iC]->SetPointError(i,errx_tmp,stat_tmp*kScaleFactor[iC]);
      //   // grB2PtSystClone[iS][iC]->SetPoint(i,x_tmp,val_tmp*kScaleFactor[iC]);
      //   // grB2PtSystClone[iS][iC]->SetPointError(i,errx_tmp,syst_tmp*kScaleFactor[iC]);
      //
      //   grB2PtStatClone[iS][iC]->GetY()[i] *= kScaleFactor[iC];
      //   grB2PtStatClone[iS][iC]->GetEY()[i] *= kScaleFactor[iC];
      //   grB2PtSystClone[iS][iC]->GetY()[i] *= kScaleFactor[iC];
      //   grB2PtSystClone[iS][iC]->GetEY()[i] *= kScaleFactor[iC];
      // }
      grB2PtStatClone[iS][iC]->Draw("samepz");
      grB2PtSystClone[iS][iC]->Draw("e2same");
      leg.AddEntry(grB2PtSystClone[iS][iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"fp");
    }
    s_dir->cd();
    leg.Draw();
    cB2[iS]->Write();
  }

  TCanvas* cB2atPtFixed[2] = {new TCanvas(Form("cB2atPtFixed_%c",kLetter[0]),Form("cB2atPtFixed_%c",kLetter[0])),new TCanvas(Form("cB2atPtFixed_%c",kLetter[1]),Form("cB2atPtFixed_%c",kLetter[1]))};

  TFile old_b2_results(Form("%sB2vsColls_pToA_55.root",kBaseOutputDir.data()));
  TCanvas* cB2atPtFixedInput = (TCanvas*)old_b2_results.Get("c1_n5");
  Requires(cB2atPtFixedInput,"old b2 at fixed pt results");
  TPad* pInput = (TPad*)cB2atPtFixedInput->GetPrimitive("c1_n5_1");
  for(int iS=0; iS<2; iS++){
    grB2atPtFixedStat[iS]->SetLineColor(kOrange-3);
    grB2atPtFixedStat[iS]->SetMarkerColor(kOrange-3);
    grB2atPtFixedStat[iS]->SetMarkerStyle(21);
    grB2atPtFixedStat[iS]->SetFillStyle(0);
    grB2atPtFixedSyst[iS]->SetLineColor(kOrange-3);
    grB2atPtFixedSyst[iS]->SetMarkerColor(kOrange-3);
    grB2atPtFixedSyst[iS]->SetMarkerStyle(21);
    grB2atPtFixedSyst[iS]->SetFillStyle(0);
    grB2atPtFixedSyst[iS]->SetName(Form("grB2atPtFixedSyst_%c",kLetter[iS]));
    grB2atPtFixedStat[iS]->SetName(Form("grB2atPtFixedStat_%c",kLetter[iS]));

    if(iS==1){
      pInput->cd();
      grB2atPtFixedStat[iS]->Draw("pz");
      grB2atPtFixedSyst[iS]->Draw("e2");
      pInput->Modified();
      pInput->Update();
      pInput->GetListOfPrimitives()->ls();
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
      cout << "GetErrorY: " << x << "is out of bounds" << endl;
      return kUnknownError;
    }
    if( x > x0 ){
      cout << "GetErrorY: Interpolating error at " << x0 << endl;
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
