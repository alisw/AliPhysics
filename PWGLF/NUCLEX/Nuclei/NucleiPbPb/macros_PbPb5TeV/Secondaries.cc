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

void Secondaries() {
  /// Taking all the histograms from the MC and data files
  TFile mc_file(kMCfilename.data());
  TFile data_file(kDataFilename.data());
  TFile output_file(kSecondariesOutput.data(),"recreate");

  /// Building the function used to fit the primary fraction distribution
  TF1 fitModel("fitFrac","1/(1-[0]*exp([1]*x))",0.4,6.);
  fitModel.SetParLimits(0, -100000, 0);
  fitModel.SetParLimits(1, -30, 30);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  

  TObjArray obj(2);
  for (auto&& list_key : *mc_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    TTList* mcList = (TTList*)mc_file.Get(list_key->GetName());
    TTList* dtList = (TTList*)data_file.Get(list_key->GetName());

    TH3F *primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTOF");
    TH3F *secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTOF");
    TH3F *data = (TH3F*)dtList->FindObject("fMDCAxyTOF");

    TDirectory *root = output_file.mkdir(list_key->GetName());
    TDirectory *datadir = root->mkdir("Data");
    TDirectory *primdir = root->mkdir("Primaries");
    TDirectory *secodir = root->mkdir("Secondaries");
    TDirectory *tffdir = root->mkdir("TFractionFitter");
    TDirectory *resdir = root->mkdir("Results");

    TAxis *pt = data->GetYaxis();
    TAxis *cen = data->GetXaxis();

    int n_cent_bins = secondaries->GetNbinsX();

    TH1D* hResTFF[kCentLength] = {nullptr};
    for(int iC=0; iC<kCentLength; iC++){
      hResTFF[iC]= new TH1D(Form("hResTFF_%i",iC),";p_{T} GeV/c;Fraction",pt->GetNbins(),pt->GetXbins()->GetArray());
    }

    for (int iB = pt->FindBin(kPtRangeMatCorrection[0]); iB <= pt->FindBin(kPtRangeMatCorrection[1]); ++iB) {

      TH1D *pr_tmp = primaries->ProjectionZ(Form("pr_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* pr = (TH1D*)pr_tmp->Rebin(knDCAbins,Form("pr_%i",iB),kDcabins);
      TH1D *sc_tmp = secondaries->ProjectionZ(Form("sc_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* sc = (TH1D*)sc_tmp->Rebin(knDCAbins,Form("sc_%i",iB),kDcabins);
      pr->SetTitle(Form("%4.1f < p_{T} #leq %4.1f (GeV/#it{c})",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
      sc->SetTitle(Form("%4.1f < p_{T} #leq %4.1f (GeV/#it{c})",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));

      primdir->cd();
      pr->Write();
      secodir->cd();
      sc->Write();

      datadir->cd();
      TH1D* dt[kCentLength] = {nullptr};
      for(int iC=0; iC<kCentLength; iC++){
        TH1D *dt_tmp = data->ProjectionZ(Form("dt_tmp%i",iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iB,iB);
        dt[iC] = (TH1D*)dt_tmp->Rebin(knDCAbins,Form("dt_%i_%i",iB,iC),kDcabins);
        dt[iC]->SetTitle(Form("%1.0f - %1.0f %%  %4.1f < p_{T} #leq %4.1f",cen->GetBinLowEdge(kCentBinsArray[iC][0]),cen->GetBinUpEdge(kCentBinsArray[iC][1]),pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
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

          dt[iC]->Scale(1.,"width");
          hs->Scale(1.,"width");
          hp->Scale(1.,"width");
          hfit->Scale(1.,"width");

          TCanvas cv(Form("tff_%i_%i",iB,iC),Form("TFractionFitter_%i_%i",iB,iC));
          cv.cd();
          dt[iC]->DrawCopy("e");

          float tot_integral = hfit->Integral(hfit->GetXaxis()->FindBin(-0.12),hfit->GetXaxis()->FindBin(0.12));
          float prim_integral = hp->Integral(hp->GetXaxis()->FindBin(-0.12),hp->GetXaxis()->FindBin(0.12));

          float ratio = prim_integral/tot_integral;

          hfit->DrawCopy("same");
          hs->SetLineColor(kRed);
          hp->SetLineColor(kBlue);
          hs->DrawCopy("same");
          hp->DrawCopy("same");
          tffdir->cd();
          TLegend leg (0.6,0.56,0.89,0.84);
          leg.SetBorderSize(0.);
          leg.AddEntry(dt[iC],"Data","pe");
          leg.AddEntry(hp,"Primaries","l");
          leg.AddEntry(hs,"Secondaries","l");
          leg.AddEntry(hfit,"TFF","l");
          leg.Draw();
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
      hResTFF[iC]->Fit(&fitModel);
      hResTFF[iC]->Write();
    }
    root->Write();
    break;
  }
}
