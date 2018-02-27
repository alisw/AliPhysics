#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>

int kCentBins[12][2] = {{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,9},{10,10},{11,11},{1,11}};
int kNcentBins = 12;

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
  /// Taking all the histograms from the MC file
  unique_ptr<TFile> input_file[3];
  input_file[0] = make_unique<TFile>("/Users/lbariogl/cernbox/Deuterons13TeV/MC/20180208_1111_LHC16h7c/AnalysisResults.root"); //reference
  input_file[1] = make_unique<TFile>("/Users/lbariogl/cernbox/Deuterons13TeV/MC/20180208_1124_LHC17d5a/AnalysisResults.root");
  input_file[2] = make_unique<TFile>("/Users/lbariogl/cernbox/Deuterons13TeV/MC/20180208_1144_LHC17d5b/AnalysisResults.root");
  TFile output_file(kMaterialOutput.data(),"recreate");

  /// Getting all the histograms
  TH2F  *fITS_TPC[2][3],*fTotal[2][3],*fITS_TPC_TOF[2][3];
  TH1D  *tpcMB[2][3],*tofMB[2][3],*totMB[2][3];
  TH1D  *effTofMB[2][3], *effTpcMB[2][3];
  for(int iFile=0; iFile<3; iFile++){
    TTList* list = (TTList*)input_file[iFile]->Get("mpuccio_deuterons_");
    for (int iS = 0; iS < 2; ++iS) {
      fTotal[iS][iFile] = dynamic_cast<TH2F*>(list->Get(Form("f%cTotal",kLetter[iS])));
      fITS_TPC[iS][iFile] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC",kLetter[iS])));
      fITS_TPC_TOF[iS][iFile] = dynamic_cast<TH2F*>(list->Get(Form("f%cITS_TPC_TOF",kLetter[iS])));

      tpcMB[iS][iFile] = fITS_TPC[iS][iFile]->ProjectionY(Form("tpcMB_%c_%d",kLetter[iS],iFile));
      tofMB[iS][iFile] = fITS_TPC_TOF[iS][iFile]->ProjectionY(Form("tofMB_%c_%d",kLetter[iS],iFile));
      totMB[iS][iFile] = fTotal[iS][iFile]->ProjectionY(Form("totMB_%c_%d",kLetter[iS],iFile));

      effTofMB[iS][iFile] = (TH1D*)DoEff(tofMB[iS][iFile],totMB[iS][iFile],Form("TofMB_%c_%d",kLetter[iS],iFile));
      effTpcMB[iS][iFile] = (TH1D*)DoEff(tpcMB[iS][iFile],totMB[iS][iFile],Form("TpcMB_%c_%d",kLetter[iS],iFile));

    }
  }
  int n_pt_bins = effTofMB[0][0]->GetNbinsX();
  TH1D *correctTPC[2], *correctTOF[2];
  for(int iS=0; iS<2;iS++){
    correctTPC[iS] = (TH1D*) effTpcMB[0][0]->Clone(Form("deuterons%ctpc",kLetter[iS]));
    correctTPC[iS]->Reset();
    correctTPC[iS]->GetYaxis()->SetTitle("RMS error");
    correctTPC[iS]->GetYaxis()->SetRangeUser(0.,0.1);
    correctTOF[iS] = (TH1D*) effTofMB[0][0]->Clone(Form("deuterons%ctof",kLetter[iS]));
    correctTOF[iS]->Reset();
    correctTOF[iS]->GetYaxis()->SetTitle("RMS error");
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
        tpc_tmp = (effTpcMB[iS][0]->GetBinContent(iB)>1.e-8) ? TMath::Abs(file_vector_tpc[2]-file_vector_tpc[1])/TMath::Sqrt(12)/effTpcMB[iS][0]->GetBinContent(iB) : 0;
      }
      correctTPC[iS]->SetBinContent(iB,tpc_tmp);
      if(bRMS){
        tof_tmp = (effTofMB[iS][0]->GetBinContent(iB)>1.e-8) ? TMath::RMS(3,file_vector_tof)/effTofMB[iS][0]->GetBinContent(iB) : 0;
      }
      else{
        tof_tmp = (effTofMB[iS][0]->GetBinContent(iB)>1.e-8) ? TMath::Abs(file_vector_tof[2]-file_vector_tof[1])/TMath::Sqrt(12)/effTofMB[iS][0]->GetBinContent(iB) : 0;
      }
      correctTOF[iS]->SetBinContent(iB,tof_tmp);
    }
  }
  for(int iS=0; iS<2; iS++){
    correctTPC[iS]->Write();
    correctTOF[iS]->Write();
  }
  for(int iFile=0; iFile<3; iFile++){
    input_file[iFile]->Close();
  }
  output_file.Close();
}
