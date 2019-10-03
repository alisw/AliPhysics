#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;

#include <TDirectory.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFractionFitter.h>
#include <TH1D.h>
#include <TH3F.h>
#include <TList.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TError.h>
#include <TROOT.h>
#include <TLine.h>

void SecondariesTPC() {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gROOT->ForceStyle();

  /// Taking all the histograms from the MC and data files
  TFile mc_file(kMCfilename.data());
  TFile data_file(kDataFilename.data());
  TFile output_file(kSecondariesTPCoutput.data(),"recreate");

  /// Building the function used to fit the primary fraction distribution
  TF1 fitModel("fitFrac","1/(1+[0]*exp([1]*x))",0.6,4.);
  fitModel.SetParLimits(0, 0, 100000);
  fitModel.SetParLimits(1, -30, 30);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gErrorIgnoreLevel=kError;

  const int nDCAbins = 34;
  const double dcabins[35] = {
    -1.30,-1.20,-1.10,-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,
    -0.35,-0.30,-0.25,-0.20,-0.15,-0.10,-0.05, 0.00, 0.05, 0.10,
     0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80,
     0.90, 1.00, 1.10, 1.20, 1.30
  };

  TObjArray obj(2);
  for (auto&& list_key : *data_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find("_MV") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;
    std::cout << "list_key : " << list_key->GetName() << std::endl;
    string mc_list_name = list_key->GetName();
    //replace(mc_list_name,"nuclei","mpuccio");
    TTList* mcList = (TTList*)mc_file.Get(mc_list_name.data());
    TTList* dtList = (TTList*)data_file.Get(list_key->GetName());

    TH3F *primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTPC");
    TH3F *secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTPC");
    TH3F *data = (TH3F*)dtList->FindObject("fMDCAxyTPC");

    string out_list = list_key->GetName();
    //replace(out_list,"mpuccio","nuclei");

    TDirectory *root = output_file.mkdir(out_list.data());
    TDirectory *datadir = root->mkdir("Data");
    TDirectory *primdir = root->mkdir("Primaries");
    TDirectory *secodir = root->mkdir("Secondaries");
    TDirectory *tffdir = root->mkdir("TFractionFitter");
    TDirectory *resdir = root->mkdir("Results");

    TAxis *pt = data->GetYaxis();
    TAxis *cen = data->GetXaxis();

    int n_cent_bins = data->GetNbinsX();

    TH1D* hResTFF[kCentLength] = {nullptr};
    TH1F* ratio_fit_data[kCentLength] = {nullptr};
    for(int iC=0; iC<kCentLength; iC++){
      hResTFF[iC]= new TH1D(Form("hResTFF_%i",iC),";p_{T} GeV/c;Fraction",pt->GetNbins(),pt->GetXbins()->GetArray());
      hResTFF[iC]->GetXaxis()->SetRangeUser(0.6,1.4);
      ratio_fit_data[iC] = new TH1F(Form("ratio_fit_data_%d",iC),Form("%2.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});data/fit",kCentLabels[iC][0],kCentLabels[iC][1]),pt->GetNbins(),pt->GetXbins()->GetArray());
    }

    for (int iB = pt->FindBin(kPtRangeMatCorrectionTPC[0]); iB <= pt->FindBin(kPtRangeMatCorrectionTPC[1]); ++iB) {

      TH1D *pr_tmp = primaries->ProjectionZ(Form("pr_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* pr = (TH1D*)pr_tmp->Rebin(nDCAbins,Form("pr_%i",iB),dcabins);
      TH1D *sc_tmp = secondaries->ProjectionZ(Form("sc_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* sc = (TH1D*)sc_tmp->Rebin(nDCAbins,Form("sc_%i",iB),dcabins);
      pr->SetTitle(Form("%4.1f < p_{T} #leq %4.1f GeV/#it{c}",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
      sc->SetTitle(Form("%4.1f < p_{T} #leq %4.1f GeV/#it{c}",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
      pr->Scale(1.,"width");
      sc->Scale(1.,"width");
      primdir->cd();
      //pr_tmp->Write();
      pr->Write();
      secodir->cd();
      //sc_tmp->Write();
      sc->Write();

      datadir->cd();
      TH1D* dt[kCentLength] = {nullptr};
      for(int iC=0; iC<kCentLength; iC++){
        TH1D *dt_tmp = data->ProjectionZ(Form("dt_tmp%i",iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iB,iB);
        dt[iC] = (TH1D*)dt_tmp->Rebin(nDCAbins,Form("dt_%i_%i",iB,iC),dcabins);
        dt[iC]->SetTitle(Form("%1.0f - %1.0f %%  %4.1f < p_{T} #leq %4.1f GeV/#it{c}",cen->GetBinLowEdge(kCentBinsArray[iC][0]),cen->GetBinUpEdge(kCentBinsArray[iC][1]),pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
        dt[iC]->Scale(1.,"width");
        dt[iC]->Write();
      }

      obj.Add(pr);
      obj.Add(sc);
      for(int iC=0; iC<kCentLength; iC++){
        TFractionFitter fitter(dt[iC],&obj,"Q");
        fitter.Constrain(0, 0., 1.);
        fitter.Constrain(1, 0., 1.);
        Double_t yieldSec = 0., yieldPri = 1., errorPri = 0., errorSec = 0.;
        char input = 'n';
        TVirtualFitter::SetMaxIterations(10000);
        Int_t result = fitter.Fit();
        if (result == 0) {
          TH1F* hp = (TH1F*)fitter.GetMCPrediction(0);
          TH1F* hs = (TH1F*)fitter.GetMCPrediction(1);
          fitter.GetResult(0, yieldPri, errorPri);
          fitter.GetResult(1, yieldSec, errorSec);
          TH1F* hfit = (TH1F*)fitter.GetPlot();
          Float_t dataIntegral = dt[iC]->Integral();
          hfit->SetLineColor(kGreen + 1);
          hfit->SetLineWidth(3);
          hs->Scale(dataIntegral * yieldSec / hs->Integral());
          hp->Scale(dataIntegral * yieldPri / hp->Integral());
          dt[iC]->SetMarkerStyle(20);
          dt[iC]->SetMarkerSize(0.5);
          dt[iC]->SetMarkerColor(kBlack);


          //hs->Scale(1.,"width");
          //hp->Scale(1.,"width");
          //hfit->Scale(1.,"width");

          TCanvas cv(Form("tff_%i_%i",iB,iC),Form("TFractionFitter_%i_%i",iB,iC));
          cv.cd();
          dt[iC]->Draw("pe");

          float tot_integral = hfit->Integral(hfit->GetXaxis()->FindBin(-0.12),hfit->GetXaxis()->FindBin(0.12));
          float prim_integral = hp->Integral(hp->GetXaxis()->FindBin(-0.12),hp->GetXaxis()->FindBin(0.12));

          float ratio = prim_integral/tot_integral;

          hfit->DrawCopy("same");
          hs->SetLineColor(kRed);
          hp->SetLineColor(kBlue);
          hs->DrawCopy("hist same");
          hp->DrawCopy("hist same");
          TLegend leg (0.6,0.56,0.89,0.84);
          leg.SetBorderSize(0.);
          leg.AddEntry(dt[iC],"Data","pe");
          leg.AddEntry(hp,"Primaries","l");
          leg.AddEntry(hs,"Secondaries","l");
          leg.AddEntry(hfit,"TFF","l");
          leg.Draw();
          tffdir->cd();
          cv.Write();
          hResTFF[iC]->SetBinContent(iB, ratio);
          hResTFF[iC]->SetBinError(iB, std::sqrt(ratio * (1. - ratio) / tot_integral));
        } else {
          cout << "In bin ";
          cout << Form("%4.1f < p_{T} #leq %4.1f",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1));
          cout << " the TFF does not converge." << endl;
        }
        delete dt[iC];
      }
      obj.Remove(pr);
      obj.Remove(sc);
      delete pr;
      delete sc;
    }
    resdir->cd();
    for(int iC=0; iC<kCentLength; iC++){
      hResTFF[iC]->Fit(&fitModel,"R");
      hResTFF[iC]->Write();
      for(int i=2; i<=9; i++){
        float data_val = hResTFF[iC]->GetBinContent(i);
        float fit_val = fitModel.Eval(hResTFF[iC]->GetBinCenter(i));
        float data_err = hResTFF[iC]->GetBinError(i);
        ratio_fit_data[iC]->SetBinContent(i,data_val/fit_val);
        ratio_fit_data[iC]->SetBinError(i,data_err/fit_val);
      }
      ratio_fit_data[iC]->GetXaxis()->SetRangeUser(0.5,1.4);
      ratio_fit_data[iC]->GetYaxis()->SetRangeUser(0.8,1.2);
      ratio_fit_data[iC]->Write();
      TCanvas* cV = new TCanvas(Form("RatioFitData_%d",iC),Form("cRatioFitData_%d",iC));
      cV->Divide(2);
      cV->cd(1);
      hResTFF[iC]->GetYaxis()->SetRangeUser(0.4,1.2);
      hResTFF[iC]->Draw();
      cV->cd(2);
      ratio_fit_data[iC]->GetYaxis()->SetTitleOffset(1.2);
      ratio_fit_data[iC]->Draw();
      TLine *line_one = new TLine(0.5,1.,1.6,1.);
      line_one->SetLineColor(kBlack);
      line_one->SetLineStyle(2);
      line_one->Draw();
      cV->Write();
    }
    root->Write();
  }
}
