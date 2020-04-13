#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TMath.h>

#include <algorithm>
using std::unique_ptr;

const int kColors[2]={kBlue,kRed};
const char* kLabels[2] = {"GEANT3", "GEANT4"};

TH1* DoEff(TH1* hRec, TH1* hTot, string name) {
  TH1* hEff = (TH1*)hRec->Clone(name.data());
  hEff->GetYaxis()->SetRangeUser(0.f,1.1f);
  hEff->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hEff->GetYaxis()->SetTitle("Efficiency x Acceptance");
  hEff->SetMarkerStyle(24);
  hEff->SetMarkerSize(0.7);
  hEff->SetLineColor(kBlack);
  hEff->SetMarkerColor(kBlack);
  for (int iBin = 1; iBin <= hEff->GetNbinsX(); ++iBin) {
    double num = hRec->GetBinContent(iBin);
    double den = hTot->GetBinContent(iBin);
    if (std::abs(den) < 1.e-24) continue;
    double eff = num/den;
    hEff->SetBinContent(iBin, eff);
    hEff->SetBinError(iBin, std::sqrt(eff * (1. - eff) / den));
  }
  hEff->Write();
  return hEff;
}

void G3G4syst(bool bRMS = false) {

  gStyle->SetOptStat(0);
  /// Taking all the histograms from the MC file
  unique_ptr<TFile> input_file[2];
  input_file[0] = utils::make_unique<TFile>(kMCgeant3.data()); //reference
  if(!input_file[0]){
    std::cout << "Impossible to find GEANT3 file" << std::endl;
    exit(1);
  }
  input_file[1] = utils::make_unique<TFile>(kMCgeant4.data());
  if(!input_file[1]){
    std::cout << "Impossible to find GEANT4 file" << std::endl;
    exit(1);
  }
  TFile output_file(kG3G4output.data(),"recreate");

  /// Getting all the histograms
  TH2F  *fITS_TPC[2][2],*fITS_TPC_TOF[2][2],*fTotal[2][2];
  TH1F  *tpcMB[2][2],*tofMB[2][2],*totMB[2][2];
  TH1F  *effTpcMB[2][2],*effTofMB[2][2];

  for(int iFile=0; iFile<2; iFile++){
    TTList* list = (TTList*)input_file[iFile]->Get("nuclei_deuteronsMB_");
    std::cout<<kFilterListNames.data()<<std::endl;
    if(!list){
      std::cout << "Check out the name of the list" << std::endl;
      exit(1);
    }
    TLegend* leg = new TLegend(0.552,0.177,0.923,0.403);
    leg->SetNColumns(2);
    leg->AddEntry((TObject*)nullptr,"TPC","");
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
      TH1F *tpcMB_tmp = (TH1F*)fITS_TPC[iS][iFile]->ProjectionY(Form("tpcMB_%c_%s_tmp",kLetter[iS],kLabels[iFile]));
      tpcMB[iS][iFile] = (TH1F*)tpcMB_tmp->Rebin(kNPtBins,Form("tpcMB_%c_%s",kLetter[iS],kLabels[iFile]),kPtBins);
      std::cout << "         TPC" << std::endl;
      std::cout << "tmp                def" << std::endl;
      for(int i=1; i<=kNPtBins; i++){
        printf("%.8f      %.8f\n", tpcMB_tmp->GetBinLowEdge(i+1),tpcMB[iS][iFile]->GetBinLowEdge(i));
      }
      TH1F *tofMB_tmp = (TH1F*)fITS_TPC_TOF[iS][iFile]->ProjectionY(Form("tofMB_%c_%s_tmp",kLetter[iS],kLabels[iFile]));
      tofMB[iS][iFile] = (TH1F*)tofMB_tmp->Rebin(kNPtBins,Form("tofMB_%c_%s",kLetter[iS],kLabels[iFile]),kPtBins);
      std::cout << "         TOF" << std::endl;
      std::cout << "tmp                def" << std::endl;
      for(int i=1; i<=kNPtBins; i++){
        printf("%.8f      %.8f\n", tofMB_tmp->GetBinLowEdge(i+1),tofMB[iS][iFile]->GetBinLowEdge(i));
      }
      TH1F *totMB_tmp = (TH1F*)fTotal[iS][iFile]->ProjectionY(Form("totMB_%c_%s_tmp",kLetter[iS],kLabels[iFile]));
      totMB[iS][iFile] = (TH1F*)totMB_tmp->Rebin(kNPtBins,Form("totMB_%c_%s",kLetter[iS],kLabels[iFile]),kPtBins);
      std::cout << "tmp                def" << std::endl;
      for(int i=1; i<=kNPtBins; i++){
        printf("%.8f      %.8f\n", totMB_tmp->GetBinLowEdge(i+1),totMB[iS][iFile]->GetBinLowEdge(i));
      }
    }
  }

  TCanvas* cEff[2];
  TLegend* aLeg[2];

  for(int iS=0; iS<2; iS++){
    cEff[iS] = new TCanvas(Form("cEff%c",kLetter[iS]),Form("cEff%c",kLetter[iS]));
    cEff[iS]->cd();
    aLeg[iS] = new TLegend(0.706767,0.189565,0.919799,0.408696);
    aLeg[iS]->SetBorderSize(0);
    for(int iFile=0; iFile<2; iFile++){

      effTpcMB[iS][iFile] = (TH1F*)DoEff(tpcMB[iS][iFile],totMB[iS][iFile],Form("TpcMB_%c_%d",kLetter[iS],iFile));
      effTpcMB[iS][iFile]->SetDirectory(0);
      effTpcMB[iS][iFile]->SetMarkerColor(kColors[iFile]);
      effTpcMB[iS][iFile]->SetMarkerSize(1.2);
      effTpcMB[iS][iFile]->SetLineColor(kColors[iFile]);
      effTpcMB[iS][iFile]->SetMarkerStyle(20);
      effTpcMB[iS][iFile]->Draw("same");

      aLeg[iS]->AddEntry(effTpcMB[iS][iFile],Form("TPC %s",kLabels[iFile]),"PL");

      effTofMB[iS][iFile] = (TH1F*)DoEff(tofMB[iS][iFile],totMB[iS][iFile],Form("TofMB_%c_%d",kLetter[iS],iFile));
      effTofMB[iS][iFile]->SetDirectory(0);
      effTofMB[iS][iFile]->SetMarkerColor(kColors[iFile]);
      effTofMB[iS][iFile]->SetMarkerSize(1.2);
      effTofMB[iS][iFile]->SetLineColor(kColors[iFile]);
      effTofMB[iS][iFile]->SetMarkerStyle(21);
      effTofMB[iS][iFile]->Draw("same");

      aLeg[iS]->AddEntry(effTofMB[iS][iFile],Form("TOF %s",kLabels[iFile]),"PL");
    }
    cEff[iS]->cd();
    aLeg[iS]->Draw();
    output_file.cd();
    cEff[iS]->Write();
  }


  TH1D *systG3G4tpc[2],*systG3G4tof[2];
  TH1D *corrG3G4tpc[2], *corrG3G4tof[2];
  for(int iS=0; iS<2;iS++){
    systG3G4tpc[iS] = new TH1D(Form("systg3g4tpc%c",kLetter[iS]),";#it{p}_{T} (GeV/#it{c});relative error",kNPtBins,kPtBins);
    systG3G4tpc[iS]->GetYaxis()->SetRangeUser(0.,1.);
    corrG3G4tpc[iS] = new TH1D(Form("corrg3g4tpc%c",kLetter[iS]),";#it{p}_{T} (GeV/#it{c});GEANT3/GEANT4",kNPtBins,kPtBins);
    corrG3G4tpc[iS]->GetYaxis()->SetRangeUser(0.5,1.5);
    //
    systG3G4tof[iS] = new TH1D(Form("systg3g4tof%c",kLetter[iS]),";#it{p}_{T} (GeV/#it{c});relative error",kNPtBins,kPtBins);
    systG3G4tof[iS]->GetYaxis()->SetRangeUser(0.,1.);
    corrG3G4tof[iS] = new TH1D(Form("corrg3g4tof%c",kLetter[iS]),";#it{p}_{T} (GeV/#it{c});GEANT3/GEANT4",kNPtBins,kPtBins);
    corrG3G4tof[iS]->GetYaxis()->SetRangeUser(0.7,1.7);
  }
  
  float tpc_tmp;
  for(int iB=1; iB<kNPtBins; iB++){
    for(int iS=0; iS<2; iS++){
      float bin_center = systG3G4tpc[iS]->GetBinCenter(iB);
      //printf("iB: %d iS: %d bin_center: %f\n",iB,iS,bin_center);
      float g3val = effTpcMB[iS][0]->GetBinContent(effTpcMB[iS][0]->FindBin(bin_center));
      float g3err = effTpcMB[iS][0]->GetBinError(effTpcMB[iS][0]->FindBin(bin_center));
      float g4val = effTpcMB[iS][1]->GetBinContent(effTpcMB[iS][1]->FindBin(bin_center));
      float g4err = effTpcMB[iS][1]->GetBinError(effTpcMB[iS][1]->FindBin(bin_center));
      float g3g4corr = g3val/g4val;
      float g3g4corrErr = g3g4corr*TMath::Sqrt((g3err*g3err/(g3val*g3val))+(g4err*g4err/(g4val*g4val)));
      if(bRMS){
        tpc_tmp = TMath::Abs(g3val-g4val)/2/g3val;
      }
      else{
        tpc_tmp = TMath::Abs(g3val-g4val)/TMath::Sqrt(12)/g3val;
      }
      systG3G4tpc[iS]->SetBinContent(iB, tpc_tmp);
      corrG3G4tpc[iS]->SetBinContent(iB,g3g4corr);
      corrG3G4tpc[iS]->SetBinError(iB,g3g4corrErr);
    }
  }
  TF1* fit_function_tpc[2]{
    new TF1("functionM_tpc","pol0",0.6,5.),
    new TF1("functionA_tpc","pol0",0.6,5.) //"[0]+[1]*exp(-[2]*x)"
  };
  corrG3G4tpc[0]->Fit("functionM_tpc","R"); 
  fit_function_tpc[1]->SetParLimits(0,0.1,2.);
  fit_function_tpc[1]->SetParLimits(1,0.,1.);
  fit_function_tpc[1]->SetParLimits(1,0.,5.);
  corrG3G4tpc[1]->Fit("functionA_tpc","R");

  TCanvas* cCorrTPC = new TCanvas("cCorrTPC","cCorrTPC");
  TLegend* legTPC = new TLegend(0.72,0.17,0.90,0.33);
  legTPC->SetBorderSize(0.);

  for(int iS=0; iS<2; iS++){
    systG3G4tpc[iS]->Write();
    corrG3G4tpc[iS]->Write();
    fit_function_tpc[iS]->Write();
    int color = (iS==0) ? kBlack : kBlue;
    plotting::SetHistStyle(corrG3G4tpc[iS],color);
    legTPC->AddEntry(corrG3G4tpc[iS],kNames[iS].data(),"PE");
  }
  cCorrTPC->cd();
  corrG3G4tpc[0]->Draw("PE");
  corrG3G4tpc[1]->Draw("Same");
  legTPC->Draw();
  cCorrTPC->Write();

  float tof_tmp;
  for(int iB=1; iB<kNPtBins; iB++){
    for(int iS=0; iS<2; iS++){
      float bin_center = systG3G4tof[iS]->GetBinCenter(iB);
      //printf("iB: %d iS: %d bin_center: %f\n",iB,iS,bin_center);
      if(bin_center < kTOFminPt) continue;
      float g3val = effTofMB[iS][0]->GetBinContent(effTofMB[iS][0]->FindBin(bin_center));
      float g3err = effTofMB[iS][0]->GetBinError(effTofMB[iS][0]->FindBin(bin_center));
      float g4val = effTofMB[iS][1]->GetBinContent(effTofMB[iS][1]->FindBin(bin_center));
      float g4err = effTofMB[iS][1]->GetBinError(effTofMB[iS][1]->FindBin(bin_center));
      float g3g4corr = g3val/g4val;
      float g3g4corrErr = g3g4corr*TMath::Sqrt((g3err*g3err/(g3val*g3val))+(g4err*g4err/(g4val*g4val)));
      if(bRMS){
        tof_tmp = TMath::Abs(g3val-g4val)/2/g3val;
      }
      else{
        tof_tmp = TMath::Abs(g3val-g4val)/TMath::Sqrt(12)/g3val;
      }
      systG3G4tof[iS]->SetBinContent(iB, tof_tmp);
      corrG3G4tof[iS]->SetBinContent(iB,g3g4corr);
      corrG3G4tof[iS]->SetBinError(iB,g3g4corrErr);
    }
  }
  TF1* fit_function_tof[2]{
    new TF1("functionM_tof","pol0",1.,5.),
    new TF1("functionA_tof","pol0",1.,5.) //"[0]+[1]*exp(-[2]*x)"
  };
  corrG3G4tof[0]->Fit("functionM_tof","R"); 
  fit_function_tof[1]->SetParLimits(0,0.1,2.);
  fit_function_tof[1]->SetParLimits(1,0.,1.);
  fit_function_tof[1]->SetParLimits(1,0.,5.);
  corrG3G4tof[1]->Fit("functionA_tof","R");
  
  TCanvas* cCorrTOF = new TCanvas("cCorrTOF","cCorrTOF");
  TLegend* legTOF = new TLegend(0.72,0.17,0.90,0.33);
  legTOF->SetBorderSize(0.);

  for(int iS=0; iS<2; iS++){
    systG3G4tof[iS]->Write();
    corrG3G4tof[iS]->Write();
    fit_function_tof[iS]->Write();
    int color = (iS==0) ? kBlack : kBlue;
    plotting::SetHistStyle(corrG3G4tof[iS],color);
    legTOF->AddEntry(corrG3G4tof[iS],kNames[iS].data(),"PE");
  }
  cCorrTOF->cd();
  corrG3G4tof[0]->Draw("PE");
  corrG3G4tof[1]->Draw("Same");
  legTOF->Draw();
  cCorrTOF->Write();

  for(int iFile=0; iFile<2; iFile++){
    input_file[iFile]->Close();
  }
  output_file.Close();
}