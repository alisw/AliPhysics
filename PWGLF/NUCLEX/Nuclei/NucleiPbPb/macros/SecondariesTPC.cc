#include "src/Common.h"

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

void SecondariesTPC() {
  /// Taking all the histograms from the MC and data files
  TFile mc_file(kMCfilename.data());
  TFile data_file(kDataFilename.data());
  TFile output_file(kSecondariesTPCoutput.data(),"recreate");

  /// Building the function used to fit the primary fraction distribution
  TF1 fitModel("fitFrac","1/(1-[0]*exp([1]*x))",0.4,6.);
  fitModel.SetParLimits(0, -100000, 0);
  fitModel.SetParLimits(1, -30, 30);

  /// fit function for systematics evaluation
  TF1* vFitModel[3];
  for(int i=0; i<3; i++){
    vFitModel[i] = new TF1(Form("fitFrac_%d",i),"1/(1-[0]*exp([1]*x))",0.4,6.);
    vFitModel[i]->SetParLimits(0, -100000, 0);
    vFitModel[i]->SetParLimits(1, -30, 30);
  }
  float integral_limits[3] = {0.07,.12,.19};
  TH1D* vResTFFsyst[3];

  gStyle->SetOptFit(1111);

  const int nDCAbins = 34;
  const double dcabins[35] = {
    -1.30,-1.20,-1.10,-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,
    -0.35,-0.30,-0.25,-0.20,-0.15,-0.10,-0.05, 0.00, 0.05, 0.10,
     0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80,
     0.90, 1.00, 1.10, 1.20, 1.30
  };


  TObjArray obj(2);
  for (auto&& list_key : *mc_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    TTList* mcList = (TTList*)mc_file.Get(list_key->GetName());
    TTList* dtList = (TTList*)data_file.Get(list_key->GetName());

    TH3F *primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTPC");
    TH3F *secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTPC");
    TH3F *data = (TH3F*)dtList->FindObject("fMDCAxyTPC");

    TDirectory *root = output_file.mkdir(list_key->GetName());
    TDirectory *datadir = root->mkdir("Data");
    TDirectory *primdir = root->mkdir("Primaries");
    TDirectory *secodir = root->mkdir("Secondaries");
    TDirectory *secodir_rebin = root->mkdir("Secondaries_rebin");
    TDirectory *tffdir = root->mkdir("TFractionFitter");
    TDirectory *resdir = root->mkdir("Results");

    TAxis *pt = data->GetYaxis();

    int n_cent_bins = secondaries->GetNbinsX();

    TH1D hResTFF("hResTFF_TPC",";p_{T} GeV/c;Fraction",pt->GetNbins(),pt->GetXbins()->GetArray());

    if(string(list_key->GetName())==kFilterListNames.data()){
      for(int i=0; i<3; i++){
        vResTFFsyst[i] = (TH1D*) hResTFF.Clone(Form("hResTFFsyst_%d",i));
      }
    }

    for (int iB = pt->FindBin(kPtRangeMatCorrection[0]); iB <= pt->FindBin(1.15); ++iB) {

      TH1D *pr_tmp = primaries->ProjectionZ(Form("pr_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* pr = (TH1D*)pr_tmp->Rebin(nDCAbins,Form("pr_%i",iB),dcabins);
      TH1D *sc_tmp = secondaries->ProjectionZ(Form("sc_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* sc = (TH1D*)sc_tmp->Rebin(nDCAbins,Form("sc_%i",iB),dcabins);
      TH1D* sc_rebin = (TH1D*)sc->Clone(Form("sc_%i",iB));
      TH1D *dt_tmp = data->ProjectionZ(Form("dt_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* dt = (TH1D*)dt_tmp->Rebin(nDCAbins,Form("dt_%i",iB),dcabins);

      dt->SetTitle(Form("%4.1f < p_{T} #leq %4.1f",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
      pr->SetTitle(dt->GetTitle());
      sc->SetTitle(dt->GetTitle());
      sc_rebin->SetTitle(dt->GetTitle());
      sc_rebin->Scale(1.,"width");

      primdir->cd();
      pr->Write();
      secodir->cd();
      sc->Write();
      secodir_rebin->cd();
      sc_rebin->Write();

      datadir->cd();
      dt->Write();

      obj.Add(pr);
      obj.Add(sc);
      TFractionFitter fitter(dt,&obj,"Q");
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
        Float_t dataIntegral = dt->Integral();
        hfit->SetLineColor(kGreen + 1);
        hfit->SetLineWidth(3);
        hs->Scale(dataIntegral * yieldSec / hs->Integral());
        hp->Scale(dataIntegral * yieldPri / hp->Integral());
        dt->SetMarkerStyle(20);
        dt->SetMarkerSize(0.5);
        dt->SetMarkerColor(kBlack);

        dt->Scale(1.,"width");
        hs->Scale(1.,"width");
        hp->Scale(1.,"width");
        hfit->Scale(1.,"width");

        TCanvas cv(Form("tff_%i",iB),Form("TFractionFitter_%i",iB));
        cv.cd();
        dt->DrawCopy("e");

        float tot_integral = hfit->Integral(hfit->GetXaxis()->FindBin(-0.12),hfit->GetXaxis()->FindBin(0.12));
        float prim_integral = hp->Integral(hp->GetXaxis()->FindBin(-0.12),hp->GetXaxis()->FindBin(0.12));

        float ratio = prim_integral/tot_integral;

        hfit->DrawCopy("same");
        hs->SetLineColor(kRed);
        hp->SetLineColor(kBlue);
        hs->DrawCopy("same");
        hp->DrawCopy("same");
        tffdir->cd();
        cv.Write();
        hResTFF.SetBinContent(iB, ratio);
        hResTFF.SetBinError(iB, std::sqrt(ratio * (1. - ratio) / tot_integral));

        if(string(list_key->GetName())==kFilterListNames.data()){
          cout << "sono entrato subito " << endl;
          for(int i=0; i<3; i++){
            tot_integral = hfit->Integral(hfit->GetXaxis()->FindBin(-1*integral_limits[i]),hfit->GetXaxis()->FindBin(integral_limits[i]));
            prim_integral = hp->Integral(hp->GetXaxis()->FindBin(-1*integral_limits[i]),hp->GetXaxis()->FindBin(integral_limits[i]));
            ratio = prim_integral/tot_integral;
            vResTFFsyst[i]->SetBinContent(iB, ratio);
            vResTFFsyst[i]->SetBinError(iB, std::sqrt(ratio * (1. - ratio) / tot_integral));
          }
        }
      } else {
        cout << "In bin ";
        cout << Form("%4.1f < p_{T} #leq %4.1f",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1));
        cout << " the TFF does not converge." << endl;
      }
      obj.Remove(pr);
      obj.Remove(sc);
      delete pr;
      delete sc;
      delete dt;
    }
    resdir->cd();
    hResTFF.Fit(&fitModel);
    hResTFF.Write();
    root->Write();
    if(string(list_key->GetName())==kFilterListNames.data()){
      TDirectory *sysdir = root->mkdir("Systematics");
      sysdir->cd();
      for(int i=0; i<3;i++){
        vResTFFsyst[i]->Fit(vFitModel[i]);
        vResTFFsyst[i]->Write();
      }
      float vec_values[3];
      TH1D hSecondSyst("hSecondSyst",";p_{T} GeV/c;relative error (%)",pt->GetNbins(),pt->GetXbins()->GetArray());
      for(int iB = pt->FindBin(0.65); iB <= pt->FindBin(3.6); ++iB){
        for(int i=0; i<3; i++){
          vec_values[i] = vFitModel[i]->Eval(pt->GetBinCenter(iB));
        }
        float value = (TMath::MaxElement(3,vec_values)-TMath::MinElement(3,vec_values))/2;
        hSecondSyst.SetBinContent(iB,value);
        hSecondSyst.SetBinError(iB,0);
      }
      hSecondSyst.Write();
    }
  }
}
