#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
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
#include <TLegend.h>

float zTest(const float mu0, const float sig0, const float mu1, const float sig1) {
  const float sigma = sqrt(sig0 * sig0 + sig1 * sig1);
  if (sigma < FLT_MIN * 10.f) return FLT_MAX;
  else return (mu0 - mu1) / sqrt(sig1 * sig1 + sig0 * sig0);
}

void Systematics() {
  TFile input_file(kSpectraOutput.data());
  TFile countsyst_file(kSignalOutput.data());
  TFile shiftsyst_file(kSignalOutput.data());
  TFile matsyst_file(kMaterialOutput.data());
  //TFile secsyst_file(kSecondariesOutput.data());
  TFile output_file(kSystematicsOutput.data(),"recreate");

  int n_centralities = 1;

  TAxis* centAxis = (TAxis*)input_file.Get("centrality");

  const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
  const int n_pt_bins = 15;

  vector<TH1F*> references(n_centralities,nullptr);
  vector<TH1F*> cutsyst(n_centralities,nullptr);
  vector<TH1F*> matsyst(n_centralities,nullptr);
  vector<TH1F*> abssyst(n_centralities,nullptr);
  vector<TH1F*> countsyst(n_centralities,nullptr);
  vector<TH1F*> shiftsyst(n_centralities,nullptr);
  //vector<TH1F*> secsyst(n_centralities,nullptr);
  vector<TH1F*> totsyst(n_centralities,nullptr);

  for (int iC = 0; iC < n_centralities; ++iC) {
    TDirectory* cent_dir = output_file.mkdir(to_string(iC).data());
    for (int iS = 0; iS < 2; ++iS) {

      TDirectory* species_dir = cent_dir->mkdir(kNames[iS].data());

      string basepath = kFilterListNames + "/" + kNames[iS] + "/TOFspectra" + to_string(iC);
      TH1F* references_tmp = (TH1F*)input_file.Get(basepath.data());
      Requires(references_tmp,"Missing reference");
      references[iC] = (TH1F*) references_tmp->Rebin(n_pt_bins,Form("TOFspectra_%d",iC),pt_bin_limits);
      auto ptAxis = references[iC]->GetXaxis();

      string count_sys_path = kFilterListNames + "/" + kNames[iS] + "/Systematic/hWidenRangeSyst" + kLetter[iS] + to_string(iC);
      TH1F* countsyst_tmp = (TH1F*)countsyst_file.Get(count_sys_path.data());
      Requires(countsyst_tmp,"Missing systematic");
      countsyst[iC] = (TH1F*)countsyst_tmp->Rebin(n_pt_bins,Form("countsyst_%d",iC),pt_bin_limits);

      string shift_sys_path = kFilterListNames + "/" + kNames[iS] + "/Systematic/hShiftRangeSyst" + kLetter[iS] + to_string(iC);
      TH1F* shiftsyst_tmp = (TH1F*)shiftsyst_file.Get(shift_sys_path.data());
      Requires(shiftsyst_tmp,"Missing systematic");
      shiftsyst[iC] = (TH1F*)shiftsyst_tmp->Rebin(n_pt_bins,Form("shiftsyst_%d",iC),pt_bin_limits);

      string mat_sys_path = Form("deuterons%ctof",kLetter[iS]);
      TH1F* matsyst_tmp = (TH1F*)matsyst_file.Get(mat_sys_path.data());
      Requires(matsyst_tmp,Form("deuterons%ctof",kLetter[iS]));
      matsyst[iC] = (TH1F*)matsyst_tmp->Rebin(n_pt_bins,Form("hMatSyst_%d",iC),pt_bin_limits);
      matsyst[iC]->Smooth(1,"R");

      cutsyst[iC] = (TH1F*)countsyst[iC]->Clone(("cutsyst" + to_string(iC)).data());
      cutsyst[iC]->Reset();

      abssyst[iC] = (TH1F*)countsyst[iC]->Clone(("abssyst" + to_string(iC)).data());
      abssyst[iC]->Reset();

      totsyst[iC] = (TH1F*)countsyst[iC]->Clone(("totsyst" + to_string(iC)).data());
      totsyst[iC]->Reset();

      //basepath = kFilterListNames + "/" + kNames[iS] + "/Systematics/hSystFit" + kLetter[iS] + to_string(iC);

      for (auto& syst : kCutNames) {
        basepath = kFilterListNames + syst.first.data() + "%i/" + kNames[iS] + "/TOFspectra" + to_string(iC);
        TDirectory* cut_dir = species_dir->mkdir(syst.first.data());
        vector<TH1F*> variations(syst.second.size(),nullptr);
        vector<TH1F*> sigmas(syst.second.size(),nullptr);
        for (size_t iV = 0; iV < syst.second.size(); ++iV) {
          TH1F* variations_tmp = (TH1F*)input_file.Get(Form(basepath.data(),iV));
          if (!variations_tmp) {
            cout << basepath.data() << " is missing." << endl;
            return;
          }
          variations[iV] = (TH1F*)variations_tmp->Rebin(n_pt_bins,Form("variation_%zu",iV),pt_bin_limits);
          variations[iV]->SetName(("cut" + to_string(iV)).data());
          plotting::SetHistStyle(variations[iV],iV);
          sigmas[iV] = (TH1F*)variations[iV]->Clone(("sigma" + to_string(iV)).data());
          sigmas[iV]->Reset();
          sigmas[iV]->SetDrawOption("e");
        }

        vector<float> rms(n_pt_bins,0.f);
        for (int iB = 1; iB <= n_pt_bins; ++iB) {
          if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
              ptAxis->GetBinCenter(iB) > kPtRange[1])
            continue;
          if (ptAxis->GetBinCenter(iB)<1) continue;
          abssyst[iC]->SetBinContent(iB,kAbsSyst[iS]);
          const float m0 = references[iC]->GetBinContent(iB);
          const float s0 = references[iC]->GetBinError(iB);

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
        cv_sigmas.DrawFrame(0.01,-5.,6.41,5.,";#it{p}_{T} (GeV/#it{c});n#sigma");
        for (auto& var : sigmas) {
          var->Draw("pesame");
        }
        cv_sigmas.Write();

        TH1F* h_rms = (TH1F*)references[iC]->Clone("rms");
        h_rms->GetYaxis()->SetTitle("RMS");
        h_rms->Reset();
        for (int iB = 1; iB <= n_pt_bins; ++iB) {
          if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
              ptAxis->GetBinCenter(iB) > kPtRange[1])
            continue;
          if(ptAxis->GetBinCenter(iB) <1) continue;
          h_rms->SetBinContent(iB,rms[iB-1]);
        }
        h_rms->Write();
      }
      for (int iB = 1; iB <= n_pt_bins; ++iB) {
        cutsyst[iC]->SetBinContent(iB,sqrt(cutsyst[iC]->GetBinContent(iB)));
      }

      if (kSmoothSystematics) {
        cutsyst[iC]->GetXaxis()->SetRange(cutsyst[iC]->FindBin(kPtRange[0]+0.01),cutsyst[iC]->FindBin(kPtRange[1]-0.01));
        //countsyst[iC]->GetXaxis()->SetRange(countsyst[iC]->FindBin(kPtRange[0]+0.01),countsyst[iC]->FindBin(kPtRange[1]-0.01));
        cutsyst[iC]->Smooth(1,"R");
        countsyst[iC]->Smooth(1,"R");
        shiftsyst[iC]->Smooth(1,"R");
      }

      for (int iB = 1; iB <= n_pt_bins; ++iB) {
        if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
            ptAxis->GetBinCenter(iB) > kPtRange[1])
          continue;
        if (ptAxis->GetBinCenter(iB) < 1.){
          cutsyst[iC]->SetBinContent(iB,0.);
          matsyst[iC]->SetBinContent(iB,0.);
          abssyst[iC]->SetBinContent(iB,0.);
          countsyst[iC]->SetBinContent(iB,0.);
          shiftsyst[iC]->SetBinContent(iB,0.);
          totsyst[iC]->SetBinContent(iB,0.);
        }
        else{
                    std::cout << ptAxis->GetBinCenter(iB) << std::endl;
          float tot = sqrt(
              cutsyst[iC]->GetBinContent(iB) * cutsyst[iC]->GetBinContent(iB) +
              matsyst[iC]->GetBinContent(iB) * matsyst[iC]->GetBinContent(iB) +
              abssyst[iC]->GetBinContent(iB) * abssyst[iC]->GetBinContent(iB) +
              countsyst[iC]->GetBinContent(iB) * countsyst[iC]->GetBinContent(iB)+
              shiftsyst[iC]->GetBinContent(iB) * shiftsyst[iC]->GetBinContent(iB)
              );
          totsyst[iC]->SetBinContent(iB,tot);
        }
      }

      TCanvas summary("summary","Summary");
      summary.DrawFrame(0.7,0.,4.1,0.3,";#it{p}_{T} (GeV/#it{c}); Systematics uncertainties");
      TLegend leg (0.6,0.56,0.89,0.84);
      leg.SetBorderSize(0);
      cutsyst[iC]->SetLineColor(plotting::kHighContrastColors[0]);
      cutsyst[iC]->Draw("same");
      leg.AddEntry(cutsyst[iC],"PID and cuts","l");
      countsyst[iC]->SetLineColor(plotting::kHighContrastColors[1]);
      countsyst[iC]->Draw("same");
      leg.AddEntry(countsyst[iC],"Range broadening","l");
      shiftsyst[iC]->SetLineColor(plotting::kHighContrastColors[5]);
      shiftsyst[iC]->Draw("same");
      leg.AddEntry(shiftsyst[iC],"Range shifting","l");
      matsyst[iC]->SetLineColor(plotting::kHighContrastColors[2]);
      matsyst[iC]->Draw("same");
      leg.AddEntry(matsyst[iC],"Material budget","l");
      abssyst[iC]->SetLineColor(plotting::kHighContrastColors[3]);
      abssyst[iC]->Draw("same");
      leg.AddEntry(abssyst[iC],"Hadronic interaction","l");
      // if(iS==0){
      //   leg.AddEntry(secsyst[iC],"Secondary fraction","l");
      //   secsyst[iC]->SetLineColor(plotting::kHighContrastColors[5]);
      //   secsyst[iC]->Draw("same");
      // }
      totsyst[iC]->SetLineColor(plotting::kHighContrastColors[4]);
      totsyst[iC]->Draw("same");
      leg.AddEntry(totsyst[iC],"Total","l");
      totsyst[iC]->SetLineWidth(2);
      totsyst[iC]->Draw("same");
      leg.Draw();

      species_dir->cd();
      cutsyst[iC]->Write("cutsyst");
      countsyst[iC]->Write("countsyst");
      abssyst[iC]->Write("abssyst");
      matsyst[iC]->Write("matsyst");
      shiftsyst[iC]->Write("shiftsyst");
      //secsyst[iC]->Write("secsyst");
      totsyst[iC]->Write("totsyst");
      summary.Write();

      // if (kPrintFigures) {
      //   summary.SaveAs((kFiguresFolder + "syst" + kLetter[iS] + to_string(iC) + ".eps").data());
      //   summary.SaveAs((kMacrosFolder + "syst" + kLetter[iS] + to_string(iC) + ".C").data());
      // }
    }
  }
  output_file.Close();

}
