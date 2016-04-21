#include "src/Common.h"
#include "src/Plotting.h"

#include <map>
#include <vector>
#include <array>
using std::array;

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>


void Final() {
  TFile spectra_file(kSpectraOutput.data());
  TFile systematics_file(kSystematicsOutput.data());
  TFile final_file(kFinalOutput.data(),"recreate");

  TAxis* centAxis = (TAxis*)spectra_file.Get("centrality");
  TAxis* ptAxis = (TAxis*)spectra_file.Get("pt");

  array<vector<TH1F*>,2> stat{vector<TH1F*>(centAxis->GetNbins(),nullptr),vector<TH1F*>(centAxis->GetNbins(),nullptr)};
  array<vector<TH1F*>,2> syst{vector<TH1F*>(centAxis->GetNbins(),nullptr),vector<TH1F*>(centAxis->GetNbins(),nullptr)};
  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = final_file.mkdir(kNames[iS].data());
    for (int iC = 0; iC < centAxis->GetNbins(); ++iC) {
      TDirectory *c_dir = s_dir->mkdir(to_string(iC).data());
      c_dir->cd();
      TH1F* totsyst = (TH1F*)systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst").data());

      string basepath = kFilterListNames + "/" + kNames[iS] + "/TOFspectra" + to_string(iC);
      stat[iS][iC]  = (TH1F*)spectra_file.Get(basepath.data());
      Requires(stat[iS][iC],basepath.data());
      syst[iS][iC]  = (TH1F*)totsyst->Clone(("syst" + to_string(iC)).data());
      for (int iB = 1; iB <= stat[iS][iC]->GetNbinsX(); ++iB) {
        syst[iS][iC]->SetBinError(iB,syst[iS][iC]->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
        syst[iS][iC]->SetBinContent(iB,stat[iS][iC]->GetBinContent(iB));
      }
      SetHistStyle(stat[iS][iC],iC);
      SetHistStyle(syst[iS][iC],iC);
      stat[iS][iC]->Write("stat");
      syst[iS][iC]->Write("syst");

      TCanvas syst_symmary("syst_summary","syst_symmary");
      TH1F* fitsyst = (TH1F*)systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/fitsyst").data());
      TH1F* cutsyst = (TH1F*)systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/cutsyst").data());
      TH1F* matsyst = (TH1F*)systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/matsyst").data());
      TH1F* abssyst = (TH1F*)systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/abssyst").data());

      SetHistStyle(fitsyst,3,20,2);
      SetHistStyle(cutsyst,2,20,2);
      SetHistStyle(matsyst,0,20,2);
      SetHistStyle(abssyst,5,20,2);
      SetHistStyle(totsyst,6,20,3);

      syst_symmary.DrawFrame(0.9 * kPtRange[0],0.001,1.1 * kPtRange[1],0.51);
      totsyst->Draw("same");
      fitsyst->Draw("same");
      cutsyst->Draw("same");
      matsyst->Draw("same");
      abssyst->Draw("same");
      syst_symmary.Write();

    }
    s_dir->cd();
    TCanvas spectra("spectra","spectra");
    spectra.DrawFrame(
        0.9 * kPtRange[0],
        0.5 * syst[iS][centAxis->GetNbins()-1]->GetBinContent(ptAxis->GetNbins()),
        1.1 * kPtRange[1],
        1.5 * syst[iS][0]->GetMaximum()
        );
    for (int iC = 0; iC < centAxis->GetNbins(); ++iC) {
      stat[iS][iC]->Draw("esamex0");
      syst[iS][iC]->Draw("e2same");
    }
    spectra.Write();
  }

  TDirectory* r_dir = final_file.mkdir("ratio");
  for (int iC = 0; iC < centAxis->GetNbins(); ++iC) {
    r_dir->mkdir(to_string(iC).data())->cd();
    stat[1][iC]->Divide(stat[0][iC]);
    syst[1][iC]->Divide(syst[0][iC]);
    stat[1][iC]->Write("stat");
    syst[1][iC]->Write("syst");
  }
}
