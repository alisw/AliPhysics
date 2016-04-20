#include "src/Common.h"
#include "src/Plotting.h"
#include <map>
#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <limits.h>
using std::array;
using std::vector;
using std::string;
using std::to_string;

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TMath.h>

float zTest(const float mu0, const float sig0, const float mu1, const float sig1) {
  const float sigma = sqrt(sig0 * sig0 + sig1 * sig1);
  if (sigma < FLT_MIN * 10.f) return FLT_MAX;
  else return (mu0 - mu1) / sqrt(sig1 * sig1 + sig0 * sig0);
}

void Systematics() {
  TFile input_file(kSpectraOutput.data());
  TFile signal_file(kSignalOutput.data());
  TFile output_file(kSystematicsOutput.data(),"recreate");

  TAxis* centAxis = (TAxis*)input_file.Get("centrality");
  TAxis* ptAxis = (TAxis*)input_file.Get("pt");

  vector<TH1F*> references(centAxis->GetNbins(),nullptr);
  vector<TH1F*> cutsyst(centAxis->GetNbins(),nullptr);
  vector<TH1F*> fitsyst(centAxis->GetNbins(),nullptr);
  vector<TH1F*> matsyst(centAxis->GetNbins(),nullptr);
  vector<TH1F*> abssyst(centAxis->GetNbins(),nullptr);
  vector<TH1F*> totsyst(centAxis->GetNbins(),nullptr);

  for (int iC = 0; iC < centAxis->GetNbins(); ++iC) {
    TDirectory* cent_dir = output_file.mkdir(to_string(iC).data());
    for (int iS = 0; iS < 2; ++iS) {
      TDirectory* species_dir = cent_dir->mkdir(kNames[iS].data());

      string basepath = kFilterListNames + "/" + kNames[iS] + "/TOFspectra" + to_string(iC);
      references[iC] = (TH1F*)input_file.Get(basepath.data());
      if (!references[iC]) {
        cout << basepath.data() << " is missing." << endl;
        return;
      }

      cutsyst[iC] = (TH1F*)references[iC]->Clone(("cutsyst" + to_string(iC)).data());
      cutsyst[iC]->Reset();

      totsyst[iC] = (TH1F*)references[iC]->Clone(("totsyst" + to_string(iC)).data());
      totsyst[iC]->Reset();

      matsyst[iC] = (TH1F*)references[iC]->Clone(("matsyst" + to_string(iC)).data());
      matsyst[iC]->Reset();

      abssyst[iC] = (TH1F*)references[iC]->Clone(("abssyst" + to_string(iC)).data());
      abssyst[iC]->Reset();

      basepath = kFilterListNames + "/" + kNames[iS] + "/Systematics/hSystFit" + kLetter[iS] + to_string(iC);
      fitsyst[iC] = (TH1F*)signal_file.Get(basepath.data());
      if (!fitsyst[iC]) {
        cout << basepath.data() << " is missing." << endl;
        return;
      }
      for (auto& syst : kCutNames) {
        basepath = kFilterListNames + "_" + syst.first.data() + "%i/" + kNames[iS] + "/TOFspectra" + to_string(iC);
        TDirectory* cut_dir = species_dir->mkdir(syst.first.data());
        vector<TH1F*> variations(syst.second.size(),nullptr);
        vector<TH1F*> sigmas(syst.second.size(),nullptr);
        for (size_t iV = 0; iV < syst.second.size(); ++iV) {
          variations[iV] = (TH1F*)input_file.Get(Form(basepath.data(),iV));
          if (!variations[iV]) {
            cout << basepath.data() << " is missing." << endl;
            return;
          }
          variations[iV]->SetName(("cut" + to_string(iV)).data());
          SetHistStyle(variations[iV],iV);
          sigmas[iV] = (TH1F*)variations[iV]->Clone(("sigma" + to_string(iV)).data());
          sigmas[iV]->Reset();
          sigmas[iV]->SetDrawOption("e");
        }

        vector<float> rms(ptAxis->GetNbins(),0.f);
        for (int iB = 1; iB <= ptAxis->GetNbins(); ++iB) {
          if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
              ptAxis->GetBinCenter(iB) > kPtRange[1])
            continue;
          matsyst[iC]->SetBinContent(iB,kMatSyst);
          abssyst[iC]->SetBinContent(iB,kAbsSyst);
          const float m0 = references[iC]->GetBinContent(iB);
          const float s0 = references[iC]->GetBinError(iB);

          vector<float> values{m0};
          for (size_t iV = 0; iV < syst.second.size(); ++iV) {
            const float m1 = variations[iV]->GetBinContent(iB);
            const float s1 = variations[iV]->GetBinError(iB);
            const float z = zTest(m0,s0,m1,s1);
            sigmas[iV]->SetBinContent(iB,z);
            if (z < 1. && kUseBarlow) {
              variations[iV]->SetBinContent(iB,0.f);
              variations[iV]->SetBinError(iB,0.f);
            } else {
              values.push_back(m1);
            }
          }
          rms[iB - 1] = TMath::RMS(values.size(),values.data()) / m0;
          cutsyst[iC]->SetBinContent(iB, cutsyst[iC]->GetBinContent(iB) + rms[iB-1] * rms[iB-1]);
        }

        cut_dir->cd();
        TCanvas cv_variations("cv_variations","cv_variations");
        cv_variations.cd();
        references[iC]->Draw();
        for (auto& var : variations)
          var->Draw("same");
        cv_variations.Write();

        TCanvas cv_ratios("cv_ratios","cv_ratios");
        cv_ratios.DrawFrame(0.01,0.01,6.41,1.99,";#it{p}_{T} (GeV/#it{c});Ratio");
        for (auto& var : variations) {
          var->Divide(references[iC]);
          var->Draw("same");
        }
        cv_ratios.Write();

        TCanvas cv_sigmas("cv_sigmas","cv_sigmas");
        cv_sigmas.DrawFrame(0.01,0.01,6.41,4.99,";#it{p}_{T} (GeV/#it{c});n#sigma");
        for (auto& var : sigmas) {
          var->Draw("same");
        }
        cv_sigmas.Write();

        TH1F* h_rms = (TH1F*)references[iC]->Clone("rms");
        h_rms->GetYaxis()->SetTitle("RMS");
        h_rms->Reset();
        for (int iB = 1; iB <= ptAxis->GetNbins(); ++iB) {
          if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
              ptAxis->GetBinCenter(iB) > kPtRange[1])
            continue;
          h_rms->SetBinContent(iB,rms[iB-1]);
        }
        h_rms->Write();
      }
      for (int iB = 1; iB <= cutsyst[iC]->GetNbinsX(); ++iB) {
        cutsyst[iC]->SetBinContent(iB,sqrt(cutsyst[iC]->GetBinContent(iB)));
        float tot = sqrt(
            cutsyst[iC]->GetBinContent(iB) * cutsyst[iC]->GetBinContent(iB) +
            matsyst[iC]->GetBinContent(iB) * matsyst[iC]->GetBinContent(iB) +
            abssyst[iC]->GetBinContent(iB) * abssyst[iC]->GetBinContent(iB) +
            fitsyst[iC]->GetBinContent(iB) * fitsyst[iC]->GetBinContent(iB)
            );
        totsyst[iC]->SetBinContent(iB,tot);
      }
      species_dir->cd();
      cutsyst[iC]->Write("cutsyst");
      fitsyst[iC]->Write("fitsyst");
      abssyst[iC]->Write("abssyst");
      matsyst[iC]->Write("matsyst");
      totsyst[iC]->Write("totsyst");
    }
  }
  output_file.Close();

}
