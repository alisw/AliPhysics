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

const char* kPrefix[2] = {"","anti"};

double GetErrorY(const TGraphErrors* gr, double x0);
void MakeItInvariant(TH1* d_stat);
TGraphErrors* GetBAPt(const TGraphErrors* grPrtInvDYieldPt,
    const TGraphErrors* grNucInvDYieldPt, const TString& name);

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

  TFile deuteron_file(kFinalOutput.data());
  TFile output_file(Form("%s/B2.root",kBaseOutputDir.data()),"recreate");

  TFile input_file(Form("%sFinal_combined_spectra_pp13TeV.root",kBaseOutputDir.data()));

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
      grB2PtStat[iS][iC] = GetBAPt(grProtSpectraStat[iC],grDeutSpectraStat[iS][iC],Form("grStat_%c_%d",kLetter[iS],iC));
      grB2PtSyst[iS][iC] = GetBAPt(grProtSpectraSyst[iC],grDeutSpectraSyst[iS][iC],Form("grSyst_%c_%d",kLetter[iS],iC));
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
    }
  }

  TCanvas* cB2[2];
  for(int iS=0; iS<2; iS++){
    TDirectory* s_dir = output_file.mkdir(kNames[iS].data());
    s_dir->cd();
    cB2[iS] = new TCanvas(kNames[iS].data(),kNames[iS].data());
    cB2[iS]->DrawFrame(
      0.,
      1e-8,
      2.,
      0.02,
      ";#it{p}_{T}/A (GeV/#it{c});#it{B}_{2} (GeV^{2}/#it{c}^{3})"
    );
    cB2[iS]->SetLeftMargin(0.1425591);
    cB2[iS]->SetRightMargin(0.027121);
    cB2[iS]->SetTopMargin(0.06053269);
    cB2[iS]->SetBottomMargin(0.1598063);
    TLegend leg(0.19,0.65,0.65,0.90);
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
      grB2PtStat[iS][iC]->Draw("samepz");
      grB2PtSyst[iS][iC]->Draw("e2same");
      leg.AddEntry(grB2PtSyst[iS][iC],Form("%4.0f - %2.0f %%",kCentLabels[iC][0],kCentLabels[iC][1]),"fp");
    }
    s_dir->cd();
    leg.Draw();
    cB2[iS]->Write();
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
    const TGraphErrors* grNucInvDYieldPt, const TString& name) {
  //
  // coalescence parameter
  //
  TGraphErrors* grBAPt = new TGraphErrors();
  grBAPt->SetName(name.Data());

  for(Int_t i=0, j=0; i < grNucInvDYieldPt->GetN(); ++i)
  {
    Double_t ptNuc, yNuc;

    grNucInvDYieldPt->GetPoint(i, ptNuc, yNuc);

    if(ptNuc/2 < 0.4) continue; // acceptance

    Double_t yPrt = grPrtInvDYieldPt->Eval(ptNuc/2); // interpolate

    if(yPrt == 0 || yNuc == 0 ) continue;

    Double_t bA = yNuc/TMath::Power(yPrt,2);

    // error
    Double_t ePrt = GetErrorY(grPrtInvDYieldPt, ptNuc/2);
    Double_t eNuc = grNucInvDYieldPt->GetErrorY(i);

    Double_t errPt = grNucInvDYieldPt->GetErrorX(i)/2;
    Double_t errBA = bA*TMath::Sqrt(TMath::Power(eNuc/yNuc,2) + TMath::Power(2*ePrt/yPrt,2));

    grBAPt->SetPoint(j, ptNuc/2, bA);
    grBAPt->SetPointError(j++, errPt, errBA);
  }

  return grBAPt;
}
