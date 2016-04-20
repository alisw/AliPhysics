#include "src/Common.h"

#include <TDirectory.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFractionFitter.h>
#include <TH1D.h>
#include <TH3F.h>
#include <TList.h>

void Secondaries() {
  /// Taking all the histograms from the MC and data files
  TFile mc_file(kMCfilename.data());
  TFile data_file(kDataFilename.data());
  TFile output_file(kSecondariesOutput.data(),"recreate");

  /// Building the function used to fit the primary fraction distribution
  TF1 fitModel("fitFrac","1/(1-[0]*exp([1]*x))",0.4,6.);
  fitModel.SetParLimits(0, -100000, 0);
  fitModel.SetParLimits(1, -30, 30);

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


    for (int iC = 1; iC <= secondaries->GetNbinsX(); ++iC) {

      TH1D hResTFF(Form("hResTFF%i",iC),";p_{T} GeV/c;Fraction",pt->GetNbins(),pt->GetXbins()->GetArray());

      for (int iB = pt->FindBin(kPtRangeMatCorrection[0]); iB <= pt->FindBin(kPtRangeMatCorrection[1]); ++iB) {

        TH1D *pr = primaries->ProjectionZ(Form("pr%i_%i",iC,iB),iC,iC,iB,iB);
        TH1D *sc = secondaries->ProjectionZ(Form("sc%i_%i",iC,iB),iC,iC,iB,iB);
        TH1D *dt = data->ProjectionZ(Form("dt%i_%i",iC,iB),iC,iC,iB,iB);

        dt->SetTitle(Form("%4.1f < p_{T} #leq %4.1f",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
        pr->SetTitle(dt->GetTitle());
        sc->SetTitle(dt->GetTitle());

        primdir->cd();
        pr->Write();
        secodir->cd();
        sc->Write();
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

          for (int i = 1; i <= dt->GetNbinsX(); ++i) {
            dt->SetBinContent(i, dt->GetBinContent(i) / dt->GetBinWidth(i));
            hs->SetBinContent(i, hs->GetBinContent(i) / hs->GetBinWidth(i));
            hp->SetBinContent(i, hp->GetBinContent(i) / hp->GetBinWidth(i));
            hfit->SetBinContent(i, hfit->GetBinContent(i) / hfit->GetBinWidth(i));
          }

          TCanvas cv(Form("tff%i_%i",iC,iB),Form("TFractionFitter%i_%i",iC,iB));
          cv.cd();
          dt->DrawCopy("e");

          hfit->DrawCopy("same");
          hs->SetLineColor(kRed);
          hp->SetLineColor(kBlue);
          hs->DrawCopy("same");
          hp->DrawCopy("same");
          tffdir->cd();
          cv.Write();
          hResTFF.SetBinContent(iB, yieldPri);
          hResTFF.SetBinError(iB, errorPri);
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
    }
    root->Write();
  }
}
