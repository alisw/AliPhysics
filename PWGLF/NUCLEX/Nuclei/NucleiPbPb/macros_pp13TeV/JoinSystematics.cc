#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLine.h"

#include <vector>

void JoinSystematics(){
  TFile TOF_systematics_file(kSystematicsOutput.data());
  TFile TPC_systematics_file(kSystematicsOutputTPC.data());
  TFile joinsystematics_file(kJoinSystematicsOutput.data(),"recreate");

  const char* syst_names[11] = {"dcaxysyst","dcazsyst","tpcsyst","pidsyst","countsyst","shiftsyst","matsyst","abssyst","siglossyst", "itstpcsyst", "totsyst"};

  const char* R_names[6] = {"dcaxy_R","dcaz_R","tpc_R","pid_R","count_R","shift_R"};

  const char* syst_labels[11] = {"DCA_{xy}","DCA_{z}","TPC clusters","TPC sigmas","Range broadening","Range shifting","Material budget","TPC_TOF Matching","Signal Loss", "ITS-TPC mathing","Total"};

  const char* syst_labels_compact[7] = {"Track selection","Signal extraction","Material budget","Hadronic interaction", "ITS-TPC mathing", "Signal Loss", "Total"};

  const char* syst_names_compact[7] = {"tracksyst","signalsyst","matsyst","abssyst","itstpcsyst","siglossyst", "totsyst"};

  std::vector<int> components[7] = {{0,1,2,3},{4,5},{6},{7},{9},{8},{10}};

  const int syst_colors[11] = {plotting::kHighContrastColors[8],plotting::kHighContrastColors[0],plotting::kHighContrastColors[7],plotting::kHighContrastColors[6],plotting::kHighContrastColors[1],plotting::kHighContrastColors[5],plotting::kHighContrastColors[2],plotting::kHighContrastColors[3],plotting::kHighContrastColors[9],plotting::kHighContrastColors[4]};

  const int compact_colors[7] = {kRed, kOrange+1, kGreen+3, kBlue, kViolet, kAzure+10, kBlack};

  const int order_leg[11] = {2,3,4,5,0,1,6,7,9,8,10};

  TH1F* vec_syst_tof[2][kCentLength][11];
  TH1F* vec_syst_tpc[2][kCentLength][11];
  TH1F* vec_syst_all[2][kCentLength][11];
  TH1F* vec_R_tof[2][kCentLength][6];
  TH1F* vec_R_tpc[2][kCentLength][6];
  TH1F* vec_R_all[2][kCentLength][6];
  TH1F* vec_syst_compact[2][kCentLength][7];
  TH1F* vec_syst_pt_uncorr_tpc[2][kCentLength];
  TH1F* vec_syst_pt_uncorr_tof[2][kCentLength];
  TH1F* vec_syst_pt_uncorr_all[2][kCentLength];
  TH1F* vec_syst_pt_corr_tpc[2][kCentLength];
  TH1F* vec_syst_pt_corr_tof[2][kCentLength];
  TH1F* vec_syst_pt_corr_all[2][kCentLength];
  TH1F* vec_syst_mult_corr_tpc[2][kCentLength];
  TH1F* vec_syst_mult_corr_tof[2][kCentLength];
  TH1F* vec_syst_mult_corr_all[2][kCentLength];
  TH1F* vec_syst_mult_uncorr_tpc[2][kCentLength];
  TH1F* vec_syst_mult_uncorr_tof[2][kCentLength];
  TH1F* vec_syst_mult_uncorr_all[2][kCentLength];

  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = joinsystematics_file.mkdir(kNames[iS].data());
    for (int iC = kCentLength-1; iC>=0; iC--) {
      TDirectory *c_dir = s_dir->mkdir(std::to_string(iC).data());
      TDirectory *tof_dir = c_dir->mkdir("tof");
      TDirectory *tpc_dir = c_dir->mkdir("tpc");
      TDirectory *joined_dir = c_dir->mkdir("joined");
      TDirectory *r_dir = c_dir->mkdir("R");
      for(int iSyst=0; iSyst<11; iSyst++){
        string tof_path = (kNames[iS] + "/" + std::to_string(iC) + "/" + syst_names[iSyst] + "_" + std::to_string(iC)).data();
        vec_syst_tof[iS][iC][iSyst] = (TH1F*) TOF_systematics_file.Get(tof_path.data());
        Requires(vec_syst_tof[iS][iC][iSyst],Form("Missing TOF %s", tof_path.data()));
        tof_dir->cd();
        vec_syst_tof[iS][iC][iSyst]->Write();
        //
        string tpc_path = (kNames[iS] + "/" + std::to_string(iC) + "/" + syst_names[iSyst] + "_" + std::to_string(iC)).data();
        vec_syst_tpc[iS][iC][iSyst] = (TH1F*) TPC_systematics_file.Get(tpc_path.data());
        Requires(vec_syst_tpc[iS][iC][iSyst],Form("Missing TPC %s", tpc_path.data()));
        tpc_dir->cd();
        vec_syst_tpc[iS][iC][iSyst]->Write();
        //
        vec_syst_all[iS][iC][iSyst] = new TH1F(Form("%s_%c%d",syst_names[iSyst],kLetter[iS],iC), Form("%4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); Systematics uncertainties",kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
        for(int iB=1; iB<=kNPtBins; iB++){
          float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
          if (vec_syst_all[iS][iC][iSyst]->GetXaxis()->GetBinCenter(iB) >kCentPtLimits[iC]) continue;
          if(bin_center < 1.){
            vec_syst_all[iS][iC][iSyst]->SetBinContent(iB,vec_syst_tpc[iS][iC][iSyst]->GetBinContent(iB));
          }
          else{
            vec_syst_all[iS][iC][iSyst]->SetBinContent(iB,vec_syst_tof[iS][iC][iSyst]->GetBinContent(iB));
          }
          vec_syst_all[iS][iC][iSyst]->SetBinError(iB,0.);
        }
        joined_dir->cd();
        vec_syst_all[iS][iC][iSyst]->Write();

        if(iSyst>5) continue;
        //
        string R_tof_path = (kNames[iS] + "/" + std::to_string(iC) + "/" + R_names[iSyst] + "_" + std::to_string(iC)).data();
        vec_R_tof[iS][iC][iSyst] = (TH1F*) TOF_systematics_file.Get(R_tof_path.data());
        Requires(vec_R_tof[iS][iC][iSyst],Form("Missing TOF %s", R_tof_path.data()));
        //
        string R_tpc_path = (kNames[iS] + "/" + std::to_string(iC) + "/" + R_names[iSyst] + "_" + std::to_string(iC)).data();
        vec_R_tpc[iS][iC][iSyst] = (TH1F*) TPC_systematics_file.Get(R_tpc_path.data());
        Requires(vec_R_tpc[iS][iC][iSyst],Form("Missing TPC %s", R_tpc_path.data()));
        //
        vec_R_all[iS][iC][iSyst] = new TH1F(Form("%s_%c%d",R_names[iSyst],kLetter[iS],iC), Form("%4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); R",kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
        for(int iB=1; iB<=kNPtBins; iB++){
          float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
          if (vec_R_all[iS][iC][iSyst]->GetXaxis()->GetBinCenter(iB) >kCentPtLimits[iC]) continue;
          if(bin_center < 1.){
            vec_R_all[iS][iC][iSyst]->SetBinContent(iB,vec_R_tpc[iS][iC][iSyst]->GetBinContent(iB));
          }
          else{
            vec_R_all[iS][iC][iSyst]->SetBinContent(iB,vec_R_tof[iS][iC][iSyst]->GetBinContent(iB));
          }
          vec_R_all[iS][iC][iSyst]->SetBinError(iB,0.);
        }
        r_dir->cd();
        vec_R_all[iS][iC][iSyst]->Write();

      }
      for(int iSystCompact = 0; iSystCompact < 7; iSystCompact++){
        vec_syst_compact[iS][iC][iSystCompact] = new TH1F(Form("compact_%s_%c%d",syst_names_compact[iSystCompact],kLetter[iS],iC), Form("%s, %4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); Systematics uncertainties",kNames[iS].data(),kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
        for(int iB=1; iB<=kNPtBins; iB++){
          vec_syst_compact[iS][iC][iSystCompact]->SetBinContent(iB,0.);
          for(int iV = 0; iV < (int)components[iSystCompact].size(); iV++){
            float current_val = vec_syst_compact[iS][iC][iSystCompact]->GetBinContent(iB);
            int place = components[iSystCompact][iV];
            float addendum = vec_syst_all[iS][iC][place]->GetBinContent(iB);
            vec_syst_compact[iS][iC][iSystCompact]->SetBinContent(iB,current_val+addendum*addendum);
          }
          float current_val = vec_syst_compact[iS][iC][iSystCompact]->GetBinContent(iB);
          vec_syst_compact[iS][iC][iSystCompact]->SetBinContent(iB,TMath::Sqrt(current_val));
        }
      }

      // pt uncorrelated uncertainties

      string tof_pt_uncorr_path = (kNames[iS] + "/" + std::to_string(iC) + "/pt_uncorr_syst_" + std::to_string(iC)).data();
      vec_syst_pt_uncorr_tof[iS][iC] = (TH1F*) TOF_systematics_file.Get(tof_pt_uncorr_path.data());
      Requires(vec_syst_pt_uncorr_tof[iS][iC],Form("Missing TOF %s", tof_pt_uncorr_path.data()));
      tof_dir->cd();
      vec_syst_pt_uncorr_tof[iS][iC]->Write();

      string tpc_pt_uncorr_path = (kNames[iS] + "/" + std::to_string(iC) + "/pt_uncorr_syst_" + std::to_string(iC)).data();
      vec_syst_pt_uncorr_tpc[iS][iC] = (TH1F*) TPC_systematics_file.Get(tpc_pt_uncorr_path.data());
      Requires(vec_syst_pt_uncorr_tpc[iS][iC],Form("Missing TPC %s", tpc_pt_uncorr_path.data()));
      tpc_dir->cd();
      vec_syst_pt_uncorr_tpc[iS][iC]->Write();

      vec_syst_pt_uncorr_all[iS][iC] = new TH1F(Form("pt_uncorr_%c%d",kLetter[iS],iC), Form("%s, %4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); Systematics uncertainties",kNames[iS].data(),kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      for(int iB=1; iB<=kNPtBins; iB++){
        float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
        if (vec_syst_pt_uncorr_all[iS][iC]->GetXaxis()->GetBinCenter(iB) >kCentPtLimits[iC]) continue;
        if(bin_center < 1.){
          vec_syst_pt_uncorr_all[iS][iC]->SetBinContent(iB,vec_syst_pt_uncorr_tpc[iS][iC]->GetBinContent(iB));
        }
        else{
          vec_syst_pt_uncorr_all[iS][iC]->SetBinContent(iB,vec_syst_pt_uncorr_tof[iS][iC]->GetBinContent(iB));
        }
        vec_syst_pt_uncorr_all[iS][iC]->SetBinError(iB,0.);
      }
      joined_dir->cd();
      vec_syst_pt_uncorr_all[iS][iC]->Write();

      // pt correlated uncertainties

      string tof_pt_corr_path = (kNames[iS] + "/" + std::to_string(iC) + "/pt_corr_syst_" + std::to_string(iC)).data();
      vec_syst_pt_corr_tof[iS][iC] = (TH1F*) TOF_systematics_file.Get(tof_pt_corr_path.data());
      Requires(vec_syst_pt_corr_tof[iS][iC],Form("Missing TOF %s", tof_pt_corr_path.data()));
      tof_dir->cd();
      vec_syst_pt_corr_tof[iS][iC]->Write();

      string tpc_pt_corr_path = (kNames[iS] + "/" + std::to_string(iC) + "/pt_corr_syst_" + std::to_string(iC)).data();
      vec_syst_pt_corr_tpc[iS][iC] = (TH1F*) TPC_systematics_file.Get(tpc_pt_corr_path.data());
      Requires(vec_syst_pt_corr_tpc[iS][iC],Form("Missing TPC %s", tpc_pt_corr_path.data()));
      tpc_dir->cd();
      vec_syst_pt_corr_tpc[iS][iC]->Write();

      vec_syst_pt_corr_all[iS][iC] = new TH1F(Form("pt_corr_%c%d",kLetter[iS],iC), Form("%s, %4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); Systematics uncertainties",kNames[iS].data(),kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      for(int iB=1; iB<=kNPtBins; iB++){
        float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
        if (vec_syst_pt_corr_all[iS][iC]->GetXaxis()->GetBinCenter(iB) >kCentPtLimits[iC]) continue;
        if(bin_center < 1.){
          vec_syst_pt_corr_all[iS][iC]->SetBinContent(iB,vec_syst_pt_corr_tpc[iS][iC]->GetBinContent(iB));
        }
        else{
          vec_syst_pt_corr_all[iS][iC]->SetBinContent(iB,vec_syst_pt_corr_tof[iS][iC]->GetBinContent(iB));
        }
        vec_syst_pt_corr_all[iS][iC]->SetBinError(iB,0.);
      }
      joined_dir->cd();
      vec_syst_pt_corr_all[iS][iC]->Write();

      // multiplicity correlated and uncorrelated uncertainties
      string tof_corr_path = (kNames[iS] + "/" + std::to_string(iC) + "/mult_corr_syst_" + std::to_string(iC)).data();
      vec_syst_mult_corr_tof[iS][iC] = (TH1F*) TOF_systematics_file.Get(tof_corr_path.data());
      Requires(vec_syst_mult_corr_tof[iS][iC],Form("Missing TOF %s", tof_corr_path.data()));
      tof_dir->cd();
      vec_syst_mult_corr_tof[iS][iC]->Write();

      string tpc_corr_path = (kNames[iS] + "/" + std::to_string(iC) + "/mult_corr_syst_" + std::to_string(iC)).data();
      vec_syst_mult_corr_tpc[iS][iC] = (TH1F*) TPC_systematics_file.Get(tpc_corr_path.data());
      Requires(vec_syst_mult_corr_tpc[iS][iC],Form("Missing TPC %s", tpc_corr_path.data()));
      tpc_dir->cd();
      vec_syst_mult_corr_tpc[iS][iC]->Write();

      vec_syst_mult_corr_all[iS][iC] = new TH1F(Form("mult_corr_%c%d",kLetter[iS],iC), Form("%s, %4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); Systematics uncertainties",kNames[iS].data(),kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      for(int iB=1; iB<=kNPtBins; iB++){
        float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
        if (vec_syst_mult_corr_all[iS][iC]->GetXaxis()->GetBinCenter(iB) >kCentPtLimits[iC]) continue;
        if(bin_center < 1.){
          vec_syst_mult_corr_all[iS][iC]->SetBinContent(iB,vec_syst_mult_corr_tpc[iS][iC]->GetBinContent(iB));
        }
        else{
          vec_syst_mult_corr_all[iS][iC]->SetBinContent(iB,vec_syst_mult_corr_tof[iS][iC]->GetBinContent(iB));;
        }
        vec_syst_mult_corr_all[iS][iC]->SetBinError(iB,0.);
      }
      joined_dir->cd();
      vec_syst_mult_corr_all[iS][iC]->Write();

      string tof_uncorr_path = (kNames[iS] + "/" + std::to_string(iC) + "/mult_uncorr_syst_" + std::to_string(iC)).data();
      vec_syst_mult_uncorr_tof[iS][iC] = (TH1F*) TOF_systematics_file.Get(tof_uncorr_path.data());
      Requires(vec_syst_mult_uncorr_tof[iS][iC],Form("Missing TOF %s", tof_uncorr_path.data()));
      tof_dir->cd();
      vec_syst_mult_uncorr_tof[iS][iC]->Write();

      string tpc_uncorr_path = (kNames[iS] + "/" + std::to_string(iC) + "/mult_uncorr_syst_" + std::to_string(iC)).data();
      vec_syst_mult_uncorr_tpc[iS][iC] = (TH1F*) TPC_systematics_file.Get(tpc_uncorr_path.data());
      Requires(vec_syst_mult_uncorr_tpc[iS][iC],Form("Missing TPC %s", tpc_uncorr_path.data()));
      tpc_dir->cd();
      vec_syst_mult_uncorr_tpc[iS][iC]->Write();

      vec_syst_mult_uncorr_all[iS][iC] = new TH1F(Form("mult_uncorr_%c%d",kLetter[iS],iC), Form("%s, %4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c}); Systematics uncertainties",kNames[iS].data(),kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      for(int iB=1; iB<=kNPtBins; iB++){
        float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
        if (vec_syst_mult_uncorr_all[iS][iC]->GetXaxis()->GetBinCenter(iB) >kCentPtLimits[iC]) continue;
        if(bin_center < 1.){
          vec_syst_mult_uncorr_all[iS][iC]->SetBinContent(iB,vec_syst_mult_uncorr_tpc[iS][iC]->GetBinContent(iB));
        }
        else{
          vec_syst_mult_uncorr_all[iS][iC]->SetBinContent(iB,vec_syst_mult_uncorr_tof[iS][iC]->GetBinContent(iB));
        }
        vec_syst_mult_uncorr_all[iS][iC]->SetBinError(iB,0.);
      }
      joined_dir->cd();
      vec_syst_mult_uncorr_all[iS][iC]->Write();
    }
  }

  for (int iS = 0; iS < 2; ++iS) {
    joinsystematics_file.cd(kNames[iS].data());
    for (int iC = 0; iC < kCentLength; ++iC) {
      TCanvas cSyst("Systematics","Systematics",3200,1800);
      cSyst.DrawFrame(
          0.4,
          0.,
          1.1*kCentPtLimits[iC],
          0.3,
          Form("%s, %4.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});Systematic uncertainties",kNames[iS].data(),kCentLabels[iC][0],kCentLabels[iC][1])
          );
      TLegend leg(0.6,0.56,0.89,0.84);
      for(int iSyst=0; iSyst<11; iSyst++){
        vec_syst_all[iS][iC][iSyst]->SetLineColor(syst_colors[iSyst]);
        iSyst==10 ? vec_syst_all[iS][iC][iSyst]->SetLineWidth(3) : vec_syst_all[iS][iC][iSyst]->SetLineWidth(2);
        vec_syst_all[iS][iC][iSyst]->Draw("same");
        if (iSyst==10) vec_syst_all[iS][iC][iSyst]->SetLineColor(plotting::kHighContrastColors[4]);
      }
      for(int iSyst=0; iSyst<11; iSyst++){
        leg.AddEntry(vec_syst_all[iS][iC][order_leg[iSyst]],syst_labels[order_leg[iSyst]],"l");
      }
      leg.Draw();
      joinsystematics_file.cd(Form("%s/%d",kNames[iS].data(),iC));      
      cSyst.Write();
      if(iS==0 && iC==0){
        cSyst.Print(Form("%splots/systematics.pdf[",kBaseOutputDir.data()),"pdf");
      }
      cSyst.Print(Form("%splots/systematics.pdf",kBaseOutputDir.data()),"pdf");
      if(iS==1 && iC==kCentLength-1){
        cSyst.Print(Form("%splots/systematics.pdf]",kBaseOutputDir.data()),"pdf");
      }

      TCanvas cSystCompact(Form("cCompSyst_%d_%c",iC,kLetter[iS]),Form("cCompSyst_%d_%c",iC,kLetter[iS]),800,600);
      cSystCompact.DrawFrame(0.4,0.,1.1*kCentPtLimits[iC],0.15,";#it{p}_{T} (GeV/#it{c});Realtive systematic uncertainties");
      TPaveText paveCompSys(0.70,0.73,0.88,0.85,"blNDC");
	    paveCompSys.SetBorderSize(0);
	    paveCompSys.SetFillColor(0);
	    paveCompSys.SetTextSize(0.03);
	    paveCompSys.AddText("This work");
	    paveCompSys.AddText("#bf{pp, #sqrt{#it{s}} = 13 TeV}");
      const char* title_tmp = (iC==kCentLength-1) ? "MB" : Form("Mult. class %s",kRomanLabels[iC]);
      paveCompSys.AddText(Form("#bf{%s, %s}",kNames[iS].data(),title_tmp));
	    paveCompSys.Draw();
      if(iS==1 && iC==8){
        paveCompSys.SetY1NDC(0.75);
        paveCompSys.SetY2NDC(0.87);
        gPad->Update();
        gPad->Modified();
      }
      TLegend leg_comp(0.16,0.66,0.63,0.85);
      leg_comp.SetNColumns(2);
      leg_comp.SetTextSize(0.03);
      for(int iSyst=0; iSyst<7; iSyst++){
        vec_syst_compact[iS][iC][iSyst]->SetLineColor(compact_colors[iSyst]);
        iSyst==6 ? vec_syst_compact[iS][iC][iSyst]->SetLineWidth(3) :  vec_syst_compact[iS][iC][iSyst]->SetLineWidth(2);
        vec_syst_compact[iS][iC][iSyst]->Draw("same");
      }
      for(int iSyst=0; iSyst<7; iSyst++){
        leg_comp.AddEntry(vec_syst_compact[iS][iC][iSyst],syst_labels_compact[iSyst],"l");
      }
      leg_comp.Draw();
      if(iS==1 && iC==8){
        leg_comp.SetY1NDC(0.68);
        leg_comp.SetY2NDC(0.87);
        gPad->Update();
        gPad->Modified();
      }
      joinsystematics_file.cd(Form("%s/%d",kNames[iS].data(),iC));      
      cSystCompact.Write();
      if(iS==0 && iC==0){
        cSystCompact.Print(Form("%splots/systematics_compact.pdf[",kBaseOutputDir.data()),"pdf");
      }
      cSystCompact.Print(Form("%splots/systematics_compact.pdf",kBaseOutputDir.data()),"pdf");
      if(iS==1 && iC==kCentLength-1){
        cSystCompact.Print(Form("%splots/systematics_compact.pdf]",kBaseOutputDir.data()),"pdf");
      }

      TCanvas cR("R","R",800,600);
      cR.DrawFrame(
          0.4,
          0.97,
          1.1*kCentPtLimits[iC],
          1.07,
          ";#it{p}_{T} (GeV/#it{c});R");
      TLine l;
      l.SetLineStyle(kDashed);
      l.SetLineColor(kBlack);
      l.DrawLine(0.4,1.,1.1*kCentPtLimits[iC],1.);
      TLegend legR(0.16,0.60,0.60,0.74,"","brNDC");
      legR.SetNColumns(2);
      legR.SetTextSize(0.03);
      for(int iSyst=0; iSyst<6; iSyst++){
        vec_R_all[iS][iC][iSyst]->SetMarkerColor(syst_colors[iSyst]);
        vec_R_all[iS][iC][iSyst]->SetMarkerStyle(20);
        vec_R_all[iS][iC][iSyst]->SetMarkerSize(1);
        vec_R_all[iS][iC][iSyst]->Draw("p same");
        legR.AddEntry(vec_R_all[iS][iC][iSyst],syst_labels[iSyst],"p");
      }
      legR.Draw();
      TPaveText paveR(0.21,0.75,0.39,0.87,"blNDC");
	    paveR.SetBorderSize(0);
	    paveR.SetFillColor(0);
	    paveR.SetTextSize(0.03);
	    paveR.AddText("This work");
	    paveR.AddText("#bf{pp, #sqrt{#it{s}} = 13 TeV}");
      paveR.AddText(Form("#bf{%s, %s}",kNames[iS].data(),title_tmp));
      paveR.Draw();
      joinsystematics_file.cd(Form("%s/%d",kNames[iS].data(),iC));      
      cR.Write();
      if(iS==0 && iC==0){
        cR.Print(Form("%splots/R.pdf[",kBaseOutputDir.data()),"pdf");
      }
      cR.Print(Form("%splots/R.pdf",kBaseOutputDir.data()),"pdf");
      if(iS==1 && iC==kCentLength-1){
        cR.Print(Form("%splots/R.pdf]",kBaseOutputDir.data()),"pdf");
      }
    }
  }
}