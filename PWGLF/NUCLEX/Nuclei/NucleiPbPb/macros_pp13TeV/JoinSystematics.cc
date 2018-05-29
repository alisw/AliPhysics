#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TLegend.h"

void JoinSystematics(){
  TFile TOF_systematics_file(kSystematicsOutput.data());
  TFile TPC_systematics_file(kSystematicsOutputTPC.data());
  TFile joinsystematics_file(kJoinSystematicsOutput.data(),"recreate");

  const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
  const int n_pt_bins = 15;

  const char* syst_names[6] = {"cutsyst","countsyst","shiftsyst","abssyst","matsyst","totsyst"};

  const char* syst_labels[6] = {"PID and cuts", "Range broadening", "Range shifting", "Hadronic interaction", "Material budget", "Total" };

  const int syst_colors[6] = {plotting::kHighContrastColors[0],plotting::kHighContrastColors[1],plotting::kHighContrastColors[5],plotting::kHighContrastColors[3],plotting::kHighContrastColors[2],plotting::kHighContrastColors[4]};

  enum syst_enum {cutsyst, countsyst, shiftsyst, abssyst, matsyst, totsyst};

  const int n_centralities = kCentLength;

  TH1F* vec_syst_tof[2][n_centralities][6];
  TH1F* vec_syst_tpc[2][n_centralities][6];
  TH1F* vec_syst_all[2][n_centralities][6];

  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = joinsystematics_file.mkdir(kNames[iS].data());
    for (int iC = 0; iC < n_centralities; ++iC) {
      TDirectory *c_dir = s_dir->mkdir(to_string(iC).data());
      c_dir->cd();
      for(int iSyst=0; iSyst<6; iSyst++){
        TH1F* partial_syst_tmp = (TH1F*) TOF_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst]).data());
        vec_syst_tof[iS][iC][iSyst] = (TH1F*) partial_syst_tmp->Rebin(n_pt_bins,Form("%s_%d_%d",syst_names[iSyst],iS,iC),pt_bin_limits);
        //
        TH1F* partial_syst_tpc_tmp = (TH1F*) TPC_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/" + syst_names[iSyst] + "_tpc").data());
        vec_syst_tpc[iS][iC][iSyst] = (TH1F*) partial_syst_tpc_tmp->Rebin(n_pt_bins,Form("%s_tpc_%d_%d",syst_names[iSyst],iS,iC),pt_bin_limits);
        //
        vec_syst_all[iS][iC][iSyst] = new TH1F(Form("systematics_%c_%d_%d",kLetter[iS],iC,iSyst), Form("%4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); Systematics uncertainties",kCentLabels[iC][0],kCentLabels[iC][1]),n_pt_bins,pt_bin_limits);
        for(int iB=1; iB<=n_pt_bins; iB++){
          if (vec_syst_all[iS][iC][iSyst]->GetXaxis()->GetBinCenter(iB) >kCentPtLimits[iC]) continue;
          if(iB<=4){
            vec_syst_all[iS][iC][iSyst]->SetBinContent(iB,vec_syst_tpc[iS][iC][iSyst]->GetBinContent(iB));
          }
          else{
            vec_syst_all[iS][iC][iSyst]->SetBinContent(iB,vec_syst_tof[iS][iC][iSyst]->GetBinContent(iB));
          }
          vec_syst_all[iS][iC][iSyst]->SetBinError(iB,0.);
        }
      }
      TCanvas cSyst("Systematics","Systematics",3200,2400);
      cSyst.DrawFrame(
          0.4,
          0.,
          1.1*kCentPtLimits[iC],
          0.3,
          Form("%4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});Systematic uncertainties",kCentLabels[iC][0],kCentLabels[iC][1])
          );
      TLegend leg(0.6,0.56,0.89,0.84);
      for(int iSyst=0; iSyst<6; iSyst++){
        vec_syst_all[iS][iC][iSyst]->SetLineColor(syst_colors[iSyst]);
        vec_syst_all[iS][iC][iSyst]->Draw("same");
        leg.AddEntry(vec_syst_all[iS][iC][iSyst],syst_labels[iSyst],"l");
      }
      leg.Draw();
      cSyst.Write();
    }
  }
}
