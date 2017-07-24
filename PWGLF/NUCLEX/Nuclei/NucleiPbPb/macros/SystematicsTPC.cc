#include "src/Common.h"
#include "src/Utils.h"
#include "src/Plotting.h"
using namespace utils;

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
#include <TLegend.h>

float zTest(const float mu0, const float sig0, const float mu1, const float sig1) {
  const float sigma = sqrt(sig0 * sig0 + sig1 * sig1);
  if (sigma < FLT_MIN * 10.f) return FLT_MAX;
  else return (mu0 - mu1) / sqrt(sig1 * sig1 + sig0 * sig0);
}

void SystematicsTPC() {
  TFile input_file(kSpectraOutput.data());
  TFile countsyst_file(kSignalOutput.data());
  TFile matsyst_file(kMaterialOutput.data());
  //TFile secsyst_tpc_file(kSecondariesTPCoutput.data());
  TFile output_file(kSystematicsOutputTPC.data(),"recreate");

  int n_centralities = 1;

  TAxis* centAxis = (TAxis*)input_file.Get("centrality");

  const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
  const int n_pt_bins = 15;

  vector<TH1F*> references_tpc(n_centralities,nullptr);
  vector<TH1F*> cutsyst_tpc(n_centralities,nullptr);
  vector<TH1F*> matsyst_tpc(n_centralities,nullptr);
  vector<TH1F*> abssyst_tpc(n_centralities,nullptr);
  vector<TH1F*> countsyst_tpc(n_centralities,nullptr);
  vector<TH1F*> totsyst_tpc(n_centralities,nullptr);

  for (int iC = 0; iC < n_centralities; ++iC) {
    TDirectory* cent_dir = output_file.mkdir(to_string(iC).data());
    for (int iS = 0; iS < 2; ++iS) {

      TDirectory* species_dir = cent_dir->mkdir(kNames[iS].data());

      string basepath = kFilterListNames + "/" + kNames[iS] + "/TPCspectra" + to_string(iC);
      TH1F* references_tpc_tmp = (TH1F*)input_file.Get(basepath.data());
      Requires(references_tpc_tmp,"Missing reference");
      references_tpc[iC] = (TH1F*) references_tpc_tmp->Rebin(n_pt_bins,Form("TPCspectra_%d",iC),pt_bin_limits);
      auto ptAxis = references_tpc[iC]->GetXaxis();

      string count_sys_path = kFilterListNames + "/" + kNames[iS] + "/Systematic/hWidenRangeSystTPC" + kLetter[iS] + to_string(iC);
      TH1F* countsyst_tpc_tmp = (TH1F*)countsyst_file.Get(count_sys_path.data());
      Requires(countsyst_tpc_tmp,"Missing systematic");
      countsyst_tpc[iC] = (TH1F*)countsyst_tpc_tmp->Rebin(n_pt_bins,Form("countsyst_tpc_%d",iC),pt_bin_limits);

      string mat_sys_path = Form("deuterons%ctpc",kLetter[iS]);
      TH1F* matsyst_tmp = (TH1F*)matsyst_file.Get(mat_sys_path.data());
      Requires(matsyst_tmp,"Missing matsysttpc");
      matsyst_tpc[iC] = (TH1F*)matsyst_tmp->Rebin(n_pt_bins,Form("hMatSyst_%d",iC),pt_bin_limits);
      matsyst_tpc[iC]->Smooth(1,"R");

      cutsyst_tpc[iC] = (TH1F*)countsyst_tpc[iC]->Clone(("cutsyst_tpc" + to_string(iC)).data());
      cutsyst_tpc[iC]->Reset();

      abssyst_tpc[iC] = (TH1F*)countsyst_tpc[iC]->Clone(("abssyst_tpc" + to_string(iC)).data());
      abssyst_tpc[iC]->Reset();

      totsyst_tpc[iC] = (TH1F*)countsyst_tpc[iC]->Clone(("totsyst_tpc" + to_string(iC)).data());
      totsyst_tpc[iC]->Reset();

      for (auto& syst : kCutNames) {
        basepath = kFilterListNames + syst.first.data() + "%i/" + kNames[iS] + "/TPCspectra" + to_string(iC);
        TDirectory* cut_dir = species_dir->mkdir(syst.first.data());
        vector<TH1F*> variations(syst.second.size(),nullptr);
        vector<TH1F*> sigmas(syst.second.size(),nullptr);
        for (size_t iV = 0; iV < syst.second.size(); ++iV) {
          TH1F* variations_tmp  = (TH1F*)input_file.Get(Form(basepath.data(),iV));
          if (!variations_tmp) {
            cout << basepath.data() << " is missing." << endl;
            return;
          }
          variations[iV] = (TH1F*)variations_tmp->Rebin(n_pt_bins,Form("variation_%d",iV),pt_bin_limits);
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
          if (ptAxis->GetBinCenter(iB) > 1.4) continue;
          abssyst_tpc[iC]->SetBinContent(iB,kAbsSyst[iS]);
          const float m0 = references_tpc[iC]->GetBinContent(iB);
          const float s0 = references_tpc[iC]->GetBinError(iB);

          vector<float> values{m0};
          vector<float> weigths{s0};
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
              weigths.push_back(s1);
            }
          }
          rms[iB - 1] = TMath::RMS(values.begin(),values.end()) / m0;
          cutsyst_tpc[iC]->SetBinContent(iB, cutsyst_tpc[iC]->GetBinContent(iB) + rms[iB-1] * rms[iB-1]);
        }

        cut_dir->cd();
        TCanvas cv_variations("cv_variations","cv_variations");
        cv_variations.cd();
        references_tpc[iC]->Draw();
        for (auto& var : variations)
          var->Draw("same");
        cv_variations.Write();

        TCanvas cv_ratios("cv_ratios","cv_ratios");
        cv_ratios.DrawFrame(0.01,0.01,6.41,1.99,";#it{p}_{T} (GeV/#it{c});Ratio");
        for (auto& var : variations) {
          var->Divide(references_tpc[iC]);
          var->Draw("same");
        }
        cv_ratios.Write();

        TCanvas cv_sigmas("cv_sigmas","cv_sigmas");
        cv_sigmas.DrawFrame(0.01,-5.,6.41,5.,";#it{p}_{T} (GeV/#it{c});n#sigma");
        for (auto& var : sigmas) {
          var->Draw("pesame");
        }
        cv_sigmas.Write();

        TH1F* h_rms = (TH1F*)references_tpc[iC]->Clone("rms");
        h_rms->GetYaxis()->SetTitle("RMS");
        h_rms->Reset();
        for (int iB = 1; iB <= ptAxis->GetNbins(); ++iB) {
          if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
              ptAxis->GetBinCenter(iB) > kPtRange[1])
            continue;
          if (ptAxis->GetBinCenter(iB) > 1.4) continue;
          h_rms->SetBinContent(iB,rms[iB-1]);
        }
        h_rms->Write();
      }
      for (int iB = 1; iB <= cutsyst_tpc[iC]->GetNbinsX(); ++iB) {
        if (ptAxis->GetBinCenter(iB) > 1.4) continue;
        cutsyst_tpc[iC]->SetBinContent(iB,sqrt(cutsyst_tpc[iC]->GetBinContent(iB)));
      }

      if (kSmoothSystematics) {
        cutsyst_tpc[iC]->GetXaxis()->SetRange(cutsyst_tpc[iC]->FindBin(kPtRange[0]+0.01),cutsyst_tpc[iC]->FindBin(kPtRange[1]-0.01));
        cutsyst_tpc[iC]->Smooth(1,"R");
        countsyst_tpc[iC]->Smooth(1,"R");
      }

      for (int iB = 1; iB <= cutsyst_tpc[iC]->GetNbinsX(); ++iB) {
        if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
            ptAxis->GetBinCenter(iB) > kPtRange[1])
          continue;
        if (ptAxis->GetBinCenter(iB) > 1.4) continue;
        float tot = sqrt(
              cutsyst_tpc[iC]->GetBinContent(iB) * cutsyst_tpc[iC]->GetBinContent(iB) +
              matsyst_tpc[iC]->GetBinContent(iB) * matsyst_tpc[iC]->GetBinContent(iB) +
              abssyst_tpc[iC]->GetBinContent(iB) * abssyst_tpc[iC]->GetBinContent(iB) +
              countsyst_tpc[iC]->GetBinContent(iB) * countsyst_tpc[iC]->GetBinContent(iB)
              );
        totsyst_tpc[iC]->SetBinContent(iB,tot);
      }

      TCanvas summary("summary","Summary");
      summary.DrawFrame(0.3,0.,1.7,0.2,";#it{p}_{T} (GeV/#it{c}); Systematics uncertainties");
      TLegend leg (0.6,0.56,0.89,0.84);
      leg.SetBorderSize(0);
      cutsyst_tpc[iC]->SetLineColor(plotting::kHighContrastColors[0]);
      cutsyst_tpc[iC]->Draw("same");
      leg.AddEntry(cutsyst_tpc[iC],"PID and cuts","l");
      countsyst_tpc[iC]->SetLineColor(plotting::kHighContrastColors[1]);
      countsyst_tpc[iC]->Draw("same");
      leg.AddEntry(countsyst_tpc[iC],"Range broadening","l");
      matsyst_tpc[iC]->SetLineColor(plotting::kHighContrastColors[2]);
      matsyst_tpc[iC]->Draw("same");
      leg.AddEntry(matsyst_tpc[iC],"Material budget","l");
      abssyst_tpc[iC]->SetLineColor(plotting::kHighContrastColors[3]);
      abssyst_tpc[iC]->Draw("same");
      leg.AddEntry(abssyst_tpc[iC],"Hadronic interaction","l");
      totsyst_tpc[iC]->SetLineColor(plotting::kHighContrastColors[4]);
      totsyst_tpc[iC]->Draw("same");
      leg.AddEntry(totsyst_tpc[iC],"Total","l");
      totsyst_tpc[iC]->SetLineWidth(2);
      totsyst_tpc[iC]->Draw("same");
      leg.Draw();

      species_dir->cd();
      cutsyst_tpc[iC]->Write("cutsyst_tpc");
      countsyst_tpc[iC]->Write("countsyst_tpc");
      abssyst_tpc[iC]->Write("abssyst_tpc");
      matsyst_tpc[iC]->Write("matsyst_tpc");
      totsyst_tpc[iC]->Write("totsyst_tpc");
      summary.Write();

    }
  }
  output_file.Close();
}
