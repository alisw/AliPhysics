#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>

#include <algorithm>

const int kCentBins[12][2] = {{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,9},{10,10},{11,11},{1,11}};
const int kNcentBins = 12;
const int kColors[3]={kBlack,kBlue,kRed};
const char* kLabels[3] = {"Standard material budget", "Decreased material budget (-4.5%)", "Increased material budget (+4.5%)"};

TH1* DoEff(TH1* tof, TH1* tot, string name) {
  TH1* efftof = (TH1*)tof->Clone(name.data());
  efftof->GetYaxis()->SetRangeUser(0.f,1.1f);
  efftof->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  efftof->GetYaxis()->SetTitle("Efficiency x Acceptance");
  efftof->SetMarkerStyle(24);
  efftof->SetMarkerSize(0.7);
  efftof->SetLineColor(kBlack);
  efftof->SetMarkerColor(kBlack);
  for (int iBin = 1; iBin <= efftof->GetNbinsX(); ++iBin) {
    double num = tof->GetBinContent(iBin);
    double den = tot->GetBinContent(iBin);
    if (std::abs(den) < 1.e-24) continue;
    double eff = num/den;
    efftof->SetBinContent(iBin, eff);
    efftof->SetBinError(iBin, std::sqrt(eff * (1. - eff) / den));
  }
  efftof->Write();
  return efftof;

}

void MaterialBudget(bool bRMS = false) {

  gStyle->SetOptStat(0);
  /// Taking all the histograms from the MC file
  unique_ptr<TFile> input_file[3];
  input_file[0] = make_unique<TFile>("/Users/lbariogl/cernbox/Deuterons13TeV/MC/20180309_2314_LHC16h7c/AnalysisResults.root"); //reference
  input_file[1] = make_unique<TFile>("/Users/lbariogl/cernbox/Deuterons13TeV/MC/20180309_2331_LHC17d5a/AnalysisResults.root");
  input_file[2] = make_unique<TFile>("/Users/lbariogl/cernbox/Deuterons13TeV/MC/20180309_2352_LHC17d5b/AnalysisResults.root");
  TFile output_file(kMaterialOutput.data(),"recreate");

  /// Getting all the histograms
  TH2F  *fITS_TPC[2][3],*fTotal[2][3],*fITS_TPC_TOF[2][3];
  TH1D  *tpcMB[2][3],*tofMB[2][3],*totMB[2][3];
  TH1D  *effTofMB[2][3], *effTpcMB[2][3];

  TCanvas* cEff[2] = {new TCanvas("cEff_M","Efficiencies matter"), new TCanvas("cEff_A","Efficiencies anti-matter")};
  TLegend* leg = new TLegend(0.552,0.177,0.923,0.403);
  leg->SetNColumns(2);
  leg->AddEntry((TObject*)nullptr,"TPC","");
  leg->AddEntry((TObject*)nullptr,"TPC + TOF","");

  for(int iFile=0; iFile<3; iFile++){
    TTList* list = (TTList*)input_file[iFile]->Get(kFilterListNames.data());
    for (int iS = 0; iS < 2; ++iS) {
      fTotal[iS][iFile] = dynamic_cast<TH2F*>(list->Get(Form("f%cTotal",kLetter[iS])));
      Requires(fTotal[iS][iFile],Form("fTotal[%d][%d]",iS,iFile));
      fTotal[iS][iFile]->SetDirectory(0);
      fITS_TPC[iS][iFile] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC",kLetter[iS])));
      Requires(fITS_TPC[iS][iFile],Form("fITS_TPC[%d][%d]",iS,iFile));
      fITS_TPC[iS][iFile]->SetDirectory(0);
      fITS_TPC_TOF[iS][iFile] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC_TOF",kLetter[iS])));
      Requires(fITS_TPC_TOF[iS][iFile],Form("fITS_TPC_TOF[%d][%d]",iS,iFile));
      fITS_TPC_TOF[iS][iFile]->SetDirectory(0);

      tpcMB[iS][iFile] = fITS_TPC[iS][iFile]->ProjectionY(Form("tpcMB_%c_%d",kLetter[iS],iFile));
      tofMB[iS][iFile] = fITS_TPC_TOF[iS][iFile]->ProjectionY(Form("tofMB_%c_%d",kLetter[iS],iFile));
      totMB[iS][iFile] = fTotal[iS][iFile]->ProjectionY(Form("totMB_%c_%d",kLetter[iS],iFile));

      effTofMB[iS][iFile] = (TH1D*)DoEff(tofMB[iS][iFile],totMB[iS][iFile],Form("TofMB_%c_%d",kLetter[iS],iFile));
      effTofMB[iS][iFile]->SetDirectory(0);
      effTofMB[iS][iFile]->SetMarkerColor(kColors[iFile]);
      effTofMB[iS][iFile]->SetMarkerSize(1.);
      effTofMB[iS][iFile]->SetLineColor(kColors[iFile]);
      effTofMB[iS][iFile]->SetMarkerStyle(21);

      cEff[iS]->cd();

      if(iFile==0) {
        effTofMB[iS][iFile]->Draw();
      }
      else{
        effTofMB[iS][iFile]->Draw("same");
      }

      effTpcMB[iS][iFile] = (TH1D*)DoEff(tpcMB[iS][iFile],totMB[iS][iFile],Form("TpcMB_%c_%d",kLetter[iS],iFile));
      effTpcMB[iS][iFile]->SetDirectory(0);
      effTpcMB[iS][iFile]->SetMarkerColor(kColors[iFile]);
      effTpcMB[iS][iFile]->SetMarkerSize(1.2);
      effTpcMB[iS][iFile]->SetLineColor(kColors[iFile]);
      effTpcMB[iS][iFile]->SetMarkerStyle(20);
      effTpcMB[iS][iFile]->Draw("same");

      if(iS==0){
        leg->AddEntry(effTpcMB[iS][iFile],kLabels[iFile],"PL");
        leg->AddEntry(effTofMB[iS][iFile],kLabels[iFile],"PL");
      }

      if(iFile==2) leg->Draw();

    }
  }
  int n_pt_bins = effTofMB[0][0]->GetNbinsX();
  TH1D *correctTPC[2], *correctTOF[2];
  for(int iS=0; iS<2;iS++){
    correctTPC[iS] = (TH1D*) effTpcMB[0][0]->Clone(Form("deuterons%ctpc",kLetter[iS]));
    correctTPC[iS]->Reset();
    correctTPC[iS]->GetYaxis()->SetTitle("relative error");
    correctTPC[iS]->GetYaxis()->SetRangeUser(0.,0.1);
    correctTOF[iS] = (TH1D*) effTofMB[0][0]->Clone(Form("deuterons%ctof",kLetter[iS]));
    correctTOF[iS]->Reset();
    correctTOF[iS]->GetYaxis()->SetTitle("relative error");
    correctTOF[iS]->GetYaxis()->SetRangeUser(0.,1.);
  }
  float file_vector_tpc[3],file_vector_tof[3];
  float tof_tmp, tpc_tmp;
  for(int iB=1; iB<=n_pt_bins; iB++){
    for(int iS=0; iS<2; iS++){
      for(int iFile=0; iFile<3; iFile++){
        file_vector_tpc[iFile] = effTpcMB[iS][iFile]->GetBinContent(iB);
        file_vector_tof[iFile] = effTofMB[iS][iFile]->GetBinContent(iB);
      }
      if(bRMS){
        tpc_tmp = (effTpcMB[iS][0]->GetBinContent(iB)>1.e-8) ? TMath::RMS(3,file_vector_tpc)/effTpcMB[iS][0]->GetBinContent(iB) : 0;
      }
      else{
        tpc_tmp = (effTpcMB[iS][0]->GetBinContent(iB)>1.e-8) ? TMath::Abs(*std::max_element(file_vector_tpc,file_vector_tpc+2)-*std::min_element(file_vector_tpc,file_vector_tpc+2))/TMath::Sqrt(12)/effTpcMB[iS][0]->GetBinContent(iB) : 0;
      }
      correctTPC[iS]->SetBinContent(iB,tpc_tmp);
      if(bRMS){
        tof_tmp = (effTofMB[iS][0]->GetBinContent(iB)>1.e-8) ? TMath::RMS(3,file_vector_tof)/effTofMB[iS][0]->GetBinContent(iB) : 0;
      }
      else{
        tof_tmp = (effTofMB[iS][0]->GetBinContent(iB)>1.e-8) ? TMath::Abs(*std::max_element(file_vector_tof,file_vector_tof+2)-*std::min_element(file_vector_tof,file_vector_tof+2))/TMath::Sqrt(12)/effTofMB[iS][0]->GetBinContent(iB) : 0;
      }
      correctTOF[iS]->SetBinContent(iB,tof_tmp);
    }
  }
  for(int iS=0; iS<2; iS++){
    cEff[iS]->Write();
    correctTPC[iS]->Write();
    correctTOF[iS]->Write();
  }
  for(int iFile=0; iFile<3; iFile++){
    input_file[iFile]->Close();
  }
  output_file.Close();
}
