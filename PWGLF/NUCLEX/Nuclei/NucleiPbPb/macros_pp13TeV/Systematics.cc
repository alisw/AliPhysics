#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include "src/Plotting.h"

#include <map>
#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <limits.h>
using std::array;
using std::vector;
using std::string;
using std::to_string;

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>

enum eSyst{
  kShift = 0,
  kDCAxy,
  kDCAz,
  kPID,
  kTPC,
  kWidth
};

const float kPtCorr[6] = {0.07, 0.34, 0.58, 0.33, 0.18, 0.57};

void Systematics(bool isTPC = false, bool bDebug = false) {
  TFile input_file(kSpectraOutput.data());
  TFile countsyst_file(kSignalOutput.data());
  TFile shiftsyst_file(kSignalOutput.data());
  TFile matsyst_file(kMaterialOutput.data());
  TFile siglos_file(kSignalLossOutput.data());
  //TFile secsyst_file(kSecondariesOutput.data());
  std::string output_string = (isTPC) ? kSystematicsOutputTPC.data() : kSystematicsOutput.data();
  TFile output_file(output_string.data(),"recreate");

  std::ofstream fDebug;

  if(bDebug){
    fDebug.open(Form("%s/systematic_debug.txt",kBaseOutputDir.data()));
  }

  vector<TH1D*> references(kCentLength,nullptr);
  vector<TH1D*> cutsyst(kCentLength,nullptr);
  vector<TH1D*> dcaxysyst(kCentLength,nullptr);
  vector<TH1D*> dcazsyst(kCentLength,nullptr);
  vector<TH1D*> tpcsyst(kCentLength,nullptr);
  vector<TH1D*> pidsyst(kCentLength,nullptr);
  vector<TH1D*> countsyst(kCentLength,nullptr);
  vector<TH1D*> shiftsyst(kCentLength,nullptr);
  vector<TH1D*> matsyst(kCentLength,nullptr);
  vector<TH1D*> abssyst(kCentLength,nullptr);
  vector<TH1D*> itstpcsyst(kCentLength,nullptr);
  vector<TH1D*> siglossyst(kCentLength,nullptr);
  vector<TH1D*> totsyst(kCentLength,nullptr);
  vector<TH1D*> pt_uncorr_syst(kCentLength,nullptr);
  vector<TH1D*> pt_corr_syst(kCentLength,nullptr);
  vector<TH1D*> mult_corr_syst(kCentLength,nullptr);
  vector<TH1D*> mult_uncorr_syst(kCentLength,nullptr);
  vector<TH1D*> totsystcheck(kCentLength,nullptr);
  vector<TH1D*> totsystptcheck(kCentLength,nullptr);
  vector<TH1D*> totsystmultcheck(kCentLength,nullptr);

  /// R-mult systematics
  vector<TH1D*> dcaxyR(kCentLength,nullptr);
  vector<TH1D*> dcazR(kCentLength,nullptr);
  vector<TH1D*> tpcR(kCentLength,nullptr);
  vector<TH1D*> pidR(kCentLength,nullptr);
  vector<TH1D*> countR(kCentLength,nullptr);
  vector<TH1D*> shiftR(kCentLength,nullptr);

  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* species_dir = output_file.mkdir(kNames[iS].data());//cent_dir->mkdir(kNames[iS].data());
    for (int iC = kCentLength -1; iC >= 0; iC--) {
      TDirectory* cent_dir = species_dir->mkdir(std::to_string(iC).data());

      std::string basepath;
      if(isTPC){
        basepath = kFilterListNames + "/" + kNames[iS] + "/" + std::to_string(iC) + "/TPC/TPCspectra" + kLetter[iS] + std::to_string(iC);
      } else {
        basepath = kFilterListNames + "/" + kNames[iS] + "/" + std::to_string(iC) + "/TOF/TOFspectra" + kLetter[iS] + std::to_string(iC);
      }
      cent_dir->cd();
      TH1D* references_tmp = (TH1D*)input_file.Get(basepath.data());
      Requires(references_tmp,Form("Missing %s",basepath.data()));
      string referenceName = (isTPC) ? Form("TPCspectra_%d",iC) : Form("TOFspectra_%d",iC);
      references[iC] = (TH1D*)references_tmp->Rebin(kNPtBins,referenceName.data(),kPtBins);
      references[iC]->Write("reference");
      auto ptAxis = references[iC]->GetXaxis();
      std::string count_sys_path;
      if(isTPC){
        count_sys_path = kFilterListNames + "/" + kNames[iS] + Form("/Systematic/C_%i/hWidenRangeSystTPC",iC) + kLetter[iS] + std::to_string(iC);
      } else {
        count_sys_path = kFilterListNames + "/" + kNames[iS] + Form("/Systematic/C_%i/hWidenRangeSyst",iC) + kLetter[iS] + std::to_string(iC);
      }
      TH1D* countsyst_tmp = (TH1D*)countsyst_file.Get(count_sys_path.data());
      Requires(countsyst_tmp,"Missing range widening systematic");
      countsyst[iC] = (TH1D*)countsyst_tmp->Rebin(kNPtBins,Form("countsyst_%d",iC),kPtBins);
      for(int iBin=1; iBin <= countsyst[iC]->GetNbinsX(); iBin++){
        countsyst[iC]->SetBinError(iBin,0);
      }
      //species_dir->cd();
      //countsyst_tmp->Write("countsyst_tmp");
      //countsyst[iC]->Write();

      string shift_sys_path;
      if(isTPC){
        shift_sys_path = kFilterListNames + "/" + kNames[iS] + Form("/Systematic/C_%i/hShiftRangeSystTPC",iC) + kLetter[iS] + std::to_string(iC);
      } else {
        shift_sys_path = kFilterListNames + "/" + kNames[iS] + Form("/Systematic/C_%i/hShiftRangeSyst",iC) + kLetter[iS] + std::to_string(iC);
      }
      TH1D* shiftsyst_tmp = (TH1D*)shiftsyst_file.Get(shift_sys_path.data());
      Requires(shiftsyst_tmp,"Missing range shifting systematic");
      shiftsyst[iC] = (TH1D*)shiftsyst_tmp->Rebin(kNPtBins,Form("shiftsyst_%d",iC),kPtBins);
      for(int iBin=1; iBin <= shiftsyst[iC]->GetNbinsX(); iBin++){
        shiftsyst[iC]->SetBinError(iBin,0);
      }
      // species_dir->cd();
      // shiftsyst_tmp->Write("shiftsyst_tmp");
      // shiftsyst[iC]->Write();

      string mat_sys_path;
      mat_sys_path = (isTPC) ? Form("deuterons%ctpc",kLetter[iS]) : Form("deuterons%ctof",kLetter[iS]);
      TH1D* matsyst_tmp = (TH1D*)matsyst_file.Get(mat_sys_path.data());
      Requires(matsyst_tmp,mat_sys_path.data());
      matsyst[iC] = (TH1D*)matsyst_tmp->Rebin(kNPtBins,Form("matsyst_%d",iC),kPtBins);
      // species_dir->cd();
      // matsyst_tmp->Write("matsyst_tmp");
      // matsyst[iC]->Write();

      string sigloss_sys_path;
      sigloss_sys_path = Form("%d/hTotMeanSignalLoss_%d",iC,iC);
      TH1D* sigloss_tmp = (TH1D*)siglos_file.Get(sigloss_sys_path.data());
      Requires(sigloss_tmp,Form("Missing %d/hTotMeanSignalLoss_%d",iC,iC));
      siglossyst[iC] = (TH1D*)sigloss_tmp->Rebin(kNPtBins,Form("siglossyst_%d",iC),kPtBins);
      siglossyst[iC]->GetYaxis()->SetRangeUser(0.,1.);
      siglossyst[iC]->SetOption("histo");
      for(int iB=1; iB<=kNPtBins; iB++){
        if(siglossyst[iC]->GetBinCenter(iB)>kCentPtLimits[iC]){ 
          siglossyst[iC]->SetBinContent(iB,0);
        }
        siglossyst[iC]->SetBinContent(iB,(siglossyst[iC]->GetBinContent(iB)-1)/2);

      }

      cutsyst[iC] = (TH1D*)countsyst[iC]->Clone(("cutsyst_" + std::to_string(iC)).data());
      cutsyst[iC]->Reset();

      dcaxysyst[iC] = (TH1D*)countsyst[iC]->Clone(("dcaxysyst_" + std::to_string(iC)).data());
      dcaxysyst[iC]->Reset();

      dcazsyst[iC] = (TH1D*)countsyst[iC]->Clone(("dcazsyst_" + std::to_string(iC)).data());
      dcazsyst[iC]->Reset();

      tpcsyst[iC] = (TH1D*)countsyst[iC]->Clone(("tpcsyst_" + std::to_string(iC)).data());
      tpcsyst[iC]->Reset();

      pidsyst[iC] = (TH1D*)countsyst[iC]->Clone(("pidsyst_" + std::to_string(iC)).data());
      pidsyst[iC]->Reset();

      abssyst[iC] = (TH1D*)countsyst[iC]->Clone(("abssyst_" + std::to_string(iC)).data());
      abssyst[iC]->Reset();

      itstpcsyst[iC] = (TH1D*)countsyst[iC]->Clone(("itstpcsyst_" + std::to_string(iC)).data());
      itstpcsyst[iC]->Reset();

      totsyst[iC] = (TH1D*)countsyst[iC]->Clone(("totsyst_" + std::to_string(iC)).data());
      totsyst[iC]->Reset();
      
      totsystcheck[iC] = (TH1D*)countsyst[iC]->Clone(("totsystcheck_" + std::to_string(iC)).data());
      totsystcheck[iC]->Reset();

      totsystptcheck[iC] = (TH1D*)countsyst[iC]->Clone(("totsystptcheck_" + std::to_string(iC)).data());
      totsystptcheck[iC]->Reset();

      totsystmultcheck[iC] = (TH1D*)countsyst[iC]->Clone(("totsystmultcheck_" + std::to_string(iC)).data());
      totsystmultcheck[iC]->Reset();

      pt_corr_syst[iC] = (TH1D*)countsyst[iC]->Clone(("pt_corr_syst_" + std::to_string(iC)).data());
      pt_corr_syst[iC]->Reset();

      pt_uncorr_syst[iC] = (TH1D*)countsyst[iC]->Clone(("pt_uncorr_syst_" + std::to_string(iC)).data());
      pt_uncorr_syst[iC]->Reset();

      mult_corr_syst[iC] = (TH1D*)countsyst[iC]->Clone(("mult_corr_syst_" + std::to_string(iC)).data());
      mult_corr_syst[iC]->Reset();

      mult_uncorr_syst[iC] = (TH1D*)countsyst[iC]->Clone(("mult_uncorr_syst_" + std::to_string(iC)).data());
      mult_uncorr_syst[iC]->Reset();

      dcaxyR[iC] = new TH1D(Form("dcaxy_R_%d",iC),";#it{p}_{T} (GeV/#it{c});R",kNPtBins,kPtBins);
      dcaxyR[iC]->GetYaxis()->SetRangeUser(0.95,1.05);
      plotting::SetHistStyle(dcaxyR[iC], plotting::kSpectraColors[0],20,"P");
      dcazR[iC] = new TH1D(Form("dcaz_R_%d",iC),";#it{p}_{T} (GeV/#it{c});R",kNPtBins,kPtBins);
      dcazR[iC]->GetYaxis()->SetRangeUser(0.95,1.05);
      plotting::SetHistStyle(dcazR[iC], plotting::kSpectraColors[1],20,"P");
      tpcR[iC] = new TH1D(Form("tpc_R_%d",iC),";#it{p}_{T} (GeV/#it{c});R",kNPtBins,kPtBins);
      tpcR[iC]->GetYaxis()->SetRangeUser(0.95,1.05);
      plotting::SetHistStyle(tpcR[iC], plotting::kSpectraColors[2],20,"P");
      pidR[iC] = new TH1D(Form("pid_R_%d",iC),";#it{p}_{T} (GeV/#it{c});R",kNPtBins,kPtBins);
      pidR[iC]->GetYaxis()->SetRangeUser(0.95,1.05);
      plotting::SetHistStyle(pidR[iC], plotting::kSpectraColors[3],20,"P");
      countR[iC] = new TH1D(Form("count_R_%d",iC),";#it{p}_{T} (GeV/#it{c});R",kNPtBins,kPtBins);
      countR[iC]->GetYaxis()->SetRangeUser(0.95,1.05);
      plotting::SetHistStyle(countR[iC], plotting::kSpectraColors[4],20,"P");
      shiftR[iC] = new TH1D(Form("shift_R_%d",iC),";#it{p}_{T} (GeV/#it{c});R",kNPtBins,kPtBins);
      shiftR[iC]->GetYaxis()->SetRangeUser(0.95,1.05);
      plotting::SetHistStyle(shiftR[iC], plotting::kSpectraColors[5],20,"P");

      for(int iB=1; iB<kNPtBins; iB++){
        double bin_center = itstpcsyst[iC]->GetBinCenter(iB);
        if (bin_center > kCentPtLimits[iC]) continue;
        if(bin_center<1){
          itstpcsyst[iC]->SetBinContent(iB,0.01);
        } else if (bin_center<2){
          itstpcsyst[iC]->SetBinContent(iB,0.015);
        } else {
          itstpcsyst[iC]->SetBinContent(iB,0.025);
        }
      }

      for(int iB=1; iB<kNPtBins; iB++){
        double bin_center = abssyst[iC]->GetBinCenter(iB);
        if (bin_center > kCentPtLimits[iC]) continue;
        abssyst[iC]->SetBinContent(iB,kAbsSyst[iS]);
      }

      for (auto& syst : kCutNames) {
        if(isTPC){
          basepath = kFilterListNames + syst.first.data() + "%i/" + kNames[iS] + "/" + std::to_string(iC) + "/TPC/TPCspectra" + kLetter[iS] + std::to_string(iC);
        } else {
          basepath = kFilterListNames + syst.first.data() + "%i/" + kNames[iS] + "/" + std::to_string(iC) + "/TOF/TOFspectra" + kLetter[iS] +std::to_string(iC);
        }
        TDirectory* cut_dir = cent_dir->mkdir(syst.first.data());
        vector<TH1D*> variations(syst.second.size(),nullptr);
        vector<TH1D*> sigmas(syst.second.size(),nullptr);
        for (size_t iV = 0; iV < syst.second.size(); ++iV) {
          TH1D* variations_tmp = (TH1D*)input_file.Get(Form(basepath.data(),iV));
          if (!variations_tmp) {
            std::cout << basepath.data() << " is missing." << std::endl;
            return;
          }
          variations[iV] = (TH1D*)variations_tmp->Rebin(kNPtBins,Form("variation_%zu",iV),kPtBins);
          cut_dir->cd();
          variations_tmp->Write(Form("variations_tmp_%zu",iV));
          variations[iV]->SetName(("cut_" + syst.first + "_" + std::to_string(iV)).data());
          plotting::SetHistStyle(variations[iV],plotting::kSpectraColors[iV]);
          variations[iV]->Write();
          sigmas[iV] = (TH1D*)variations[iV]->Clone(("sigma_" + syst.first + "_" + std::to_string(iV)).data());
          sigmas[iV]->Reset();
          sigmas[iV]->SetDrawOption("e");
        }

        vector<double> rms(kNPtBins,0.f);

        TH1D* h_rms = (TH1D*)references[iC]->Clone(Form("rms_%s",syst.first.data()));
        h_rms->GetYaxis()->SetTitle("RMS");
        h_rms->Reset();

        for (int iB = 1; iB <= kNPtBins; ++iB) {
          if (ptAxis->GetBinCenter(iB) < kPtRange[0]||
              ptAxis->GetBinCenter(iB) > kPtRange[1])
            continue;
          if(isTPC){
            if (ptAxis->GetBinCenter(iB) > kTPCmaxPt) continue;
          } else {
            if (ptAxis->GetBinCenter(iB)<kTOFminPt) continue;
            if (ptAxis->GetBinCenter(iB)>kCentPtLimits[iC]) continue;
          }
          const double m0 = references[iC]->GetBinContent(iB);
          const double s0 = references[iC]->GetBinError(iB);
          if(bDebug){
            fDebug << std::endl << std::endl << "pT bin : [" << ptAxis->GetBinLowEdge(iB) << ", " << ptAxis->GetBinLowEdge(iB+1) <<"] GeV/c" << std::endl;
          }
          vector<double> values{m0};
          vector<double> weigths{s0};
          for (size_t iV = 0; iV < syst.second.size(); ++iV) {
            const double m1 = variations[iV]->GetBinContent(iB);
            const double s1 = variations[iV]->GetBinError(iB);
            const double z = utils::zTest(m0,s0,m1,s1);
            if(bDebug){
              fDebug << "State: " << kNames[iS] << std::endl;
              fDebug << "Multiplicity class: " << iC << std::endl;
              fDebug << "Cut name: " << syst.first.data() << std::endl << std::endl;
              fDebug << "cut iV : " << iV << std::endl;
              fDebug << "   m0 : " << m0 << std::endl;
              fDebug << "   s0 : " << s0 << std::endl;
              fDebug << "   m1 : " << m1 << std::endl;
              fDebug << "   s1 : " << s1 << std::endl;
              fDebug << "   m1 / m0 : " << m1/m0 << std::endl;
              fDebug << "   m1 - m0 : " << m1-m0 << std::endl;
              fDebug << "   s1 - s0 : " << std::abs(s1-s0) << std::endl;
              fDebug << "   z : " << z << std::endl;
            }
            sigmas[iV]->SetBinContent(iB,z);
            if (std::abs(z) < 2. && kUseBarlow) {
              variations[iV]->SetBinContent(iB,0.f);
              variations[iV]->SetBinError(iB,0.f);
            } else {
              values.push_back(m1);
              weigths.push_back(s1);
            }
          }
          rms[iB - 1] = TMath::RMS(values.begin(),values.end());
          rms[iB - 1] /= m0;
          if(bDebug){
            fDebug << "values : ";
            for(auto &p : values){
              fDebug << p << ",";
            }
            fDebug << std::endl << "rms : " << rms[iB-1] << std::endl;
            fDebug << std::endl << "rms/m0 : " << rms[iB-1] << std::endl;
          }
          cutsyst[iC]->SetBinContent(iB, cutsyst[iC]->GetBinContent(iB) + rms[iB-1] * rms[iB-1]);
          if (syst.first == "dcaxy") dcaxysyst[iC]->SetBinContent(iB, rms[iB-1] * rms[iB-1]);
          else if (syst.first == "dcaz") dcazsyst[iC]->SetBinContent(iB, rms[iB-1] * rms[iB-1]);
          else if (syst.first == "tpc") tpcsyst[iC]->SetBinContent(iB, rms[iB-1] * rms[iB-1]);
          else if (syst.first == "pid") pidsyst[iC]->SetBinContent(iB, rms[iB-1] * rms[iB-1]);
          h_rms->SetBinContent(iB,rms[iB-1]);
        }

        cut_dir->cd();
        TCanvas cv_variations("cv_variations","cv_variations");
        cv_variations.cd();
        references[iC]->Draw();
        for (auto& var : variations)
          var->Draw("same");
        cv_variations.Write();

        h_rms->Write();

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
      }
      for (int iB = 1; iB <= kNPtBins; ++iB) {
        if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
            ptAxis->GetBinCenter(iB) > kPtRange[1])
          continue;
        if (isTPC) {
          if (ptAxis->GetBinCenter(iB) > kTPCmaxPt) continue;
        } else {
          if (ptAxis->GetBinCenter(iB)<kTOFminPt) continue;
          if (ptAxis->GetBinCenter(iB)>kCentPtLimits[iC]) continue;
        }
        cutsyst[iC]->SetBinContent(iB,sqrt(cutsyst[iC]->GetBinContent(iB)));
        dcaxysyst[iC]->SetBinContent(iB,sqrt(dcaxysyst[iC]->GetBinContent(iB)));
        dcazsyst[iC]->SetBinContent(iB,sqrt(dcazsyst[iC]->GetBinContent(iB)));
        tpcsyst[iC]->SetBinContent(iB,sqrt(tpcsyst[iC]->GetBinContent(iB)));
        pidsyst[iC]->SetBinContent(iB,sqrt(pidsyst[iC]->GetBinContent(iB)));
      }

      if (kSmoothSystematics) {
        if(isTPC){
          cutsyst[iC]->GetXaxis()->SetRange(cutsyst[iC]->FindBin(kPtRange[0]+0.01),cutsyst[iC]->FindBin(kTPCmaxPt-0.01));
          dcaxysyst[iC]->GetXaxis()->SetRange(dcaxysyst[iC]->FindBin(kPtRange[0]+0.01),dcaxysyst[iC]->FindBin(kTPCmaxPt-0.01));
          dcazsyst[iC]->GetXaxis()->SetRange(dcazsyst[iC]->FindBin(kPtRange[0]+0.01),dcazsyst[iC]->FindBin(kTPCmaxPt-0.01));
          tpcsyst[iC]->GetXaxis()->SetRange(tpcsyst[iC]->FindBin(kPtRange[0]+0.01),tpcsyst[iC]->FindBin(kTPCmaxPt-0.01));
          pidsyst[iC]->GetXaxis()->SetRange(pidsyst[iC]->FindBin(kPtRange[0]+0.01),pidsyst[iC]->FindBin(kTPCmaxPt-0.01));
          countsyst[iC]->GetXaxis()->SetRange(countsyst[iC]->FindBin(kPtRange[0]+0.01),countsyst[iC]->FindBin(kTPCmaxPt-0.01));
          shiftsyst[iC]->GetXaxis()->SetRange(shiftsyst[iC]->FindBin(kPtRange[0]+0.01),shiftsyst[iC]->FindBin(kTPCmaxPt-0.01));
          matsyst[iC]->GetXaxis()->SetRange(matsyst[iC]->FindBin(kPtRange[0]+0.01),matsyst[iC]->FindBin(kTPCmaxPt-0.01));
        } else {
          cutsyst[iC]->GetXaxis()->SetRange(cutsyst[iC]->FindBin(kTOFminPt+0.01),cutsyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
          dcaxysyst[iC]->GetXaxis()->SetRange(dcaxysyst[iC]->FindBin(kTOFminPt+0.01),dcaxysyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
          dcazsyst[iC]->GetXaxis()->SetRange(dcazsyst[iC]->FindBin(kTOFminPt+0.01),dcazsyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
          tpcsyst[iC]->GetXaxis()->SetRange(tpcsyst[iC]->FindBin(kTOFminPt+0.01),tpcsyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
          pidsyst[iC]->GetXaxis()->SetRange(pidsyst[iC]->FindBin(kTOFminPt+0.01),pidsyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
          countsyst[iC]->GetXaxis()->SetRange(countsyst[iC]->FindBin(kTOFminPt+0.01),countsyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
          shiftsyst[iC]->GetXaxis()->SetRange(shiftsyst[iC]->FindBin(kTOFminPt+0.01),shiftsyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
          matsyst[iC]->GetXaxis()->SetRange(matsyst[iC]->FindBin(kTOFminPt+0.01),matsyst[iC]->FindBin(kCentPtLimits[iC]-0.01));
        }
        cutsyst[iC]->Smooth(1,"R");
        dcaxysyst[iC]->Smooth(1,"R");
        dcazsyst[iC]->Smooth(1,"R");
        tpcsyst[iC]->Smooth(1,"R");
        pidsyst[iC]->Smooth(1,"R");
        countsyst[iC]->Smooth(1,"R");
        shiftsyst[iC]->Smooth(1,"R");
        matsyst[iC]->Smooth(1,"R");
      }

      for (int iB = 1; iB <= kNPtBins; ++iB) {
        if (ptAxis->GetBinCenter(iB) < kPtRange[0] ||
            ptAxis->GetBinCenter(iB) > kPtRange[1])
          continue;
        if (ptAxis->GetBinCenter(iB)>kCentPtLimits[iC]){
          matsyst[iC]->SetBinContent(iB,0.);
        }
        if( (!isTPC && (ptAxis->GetBinCenter(iB) < kTOFminPt || ptAxis->GetBinCenter(iB) > kCentPtLimits[iC])) || (isTPC && ptAxis->GetBinCenter(iB) > kTPCmaxPt) ){
          cutsyst[iC]->SetBinContent(iB,0.);
          dcaxysyst[iC]->SetBinContent(iB,0.);
          dcazsyst[iC]->SetBinContent(iB,0.);
          tpcsyst[iC]->SetBinContent(iB,0.);
          pidsyst[iC]->SetBinContent(iB,0.);
          matsyst[iC]->SetBinContent(iB,0.);
          abssyst[iC]->SetBinContent(iB,0.);
          countsyst[iC]->SetBinContent(iB,0.);
          shiftsyst[iC]->SetBinContent(iB,0.);
          siglossyst[iC]->SetBinContent(iB,0);
          itstpcsyst[iC]->SetBinContent(iB,0.);
          totsyst[iC]->SetBinContent(iB,0.);
          pt_corr_syst[iC]->SetBinContent(iB,0.);
          pt_uncorr_syst[iC]->SetBinContent(iB,0.);
          mult_corr_syst[iC]->SetBinContent(iB,0.);
          mult_uncorr_syst[iC]->SetBinContent(iB,0.);
          totsystcheck[iC]->SetBinContent(iB,0.);
          totsystptcheck[iC]->SetBinContent(iB,0.);
          totsystmultcheck[iC]->SetBinContent(iB,0.);
        }
        else{

          /// Computing R
          double dcaxyR_val = (1 + dcaxysyst[iC]->GetBinContent(iB))/(1 + dcaxysyst[kCentLength-1]->GetBinContent(iB));
          dcaxyR[iC]->SetBinContent(iB,dcaxyR_val);

          double dcazR_val = (1 + dcazsyst[iC]->GetBinContent(iB))/(1 + dcazsyst[kCentLength-1]->GetBinContent(iB));
          dcazR[iC]->SetBinContent(iB,dcazR_val);

          double tpcR_val = (1 + tpcsyst[iC]->GetBinContent(iB))/(1 + tpcsyst[kCentLength-1]->GetBinContent(iB));
          tpcR[iC]->SetBinContent(iB,tpcR_val);

          double pidR_val = (1 + pidsyst[iC]->GetBinContent(iB))/(1 + pidsyst[kCentLength-1]->GetBinContent(iB));
          pidR[iC]->SetBinContent(iB,pidR_val);

          double countR_val = (1 + countsyst[iC]->GetBinContent(iB))/(1 + countsyst[kCentLength-1]->GetBinContent(iB));
          countR[iC]->SetBinContent(iB,countR_val);

          double shiftR_val = (1 + shiftsyst[iC]->GetBinContent(iB))/(1 + shiftsyst[kCentLength-1]->GetBinContent(iB));
          shiftR[iC]->SetBinContent(iB,shiftR_val);

          double tot =
              utils::Sq(cutsyst[iC]->GetBinContent(iB)) +
              utils::Sq(matsyst[iC]->GetBinContent(iB)) +
              utils::Sq(abssyst[iC]->GetBinContent(iB)) +
              utils::Sq(countsyst[iC]->GetBinContent(iB)) +
              utils::Sq(shiftsyst[iC]->GetBinContent(iB)) +
              utils::Sq(itstpcsyst[iC]->GetBinContent(iB)) +
              utils::Sq(siglossyst[iC]->GetBinContent(iB));
          tot = sqrt(tot);
          totsyst[iC]->SetBinContent(iB,tot);

          double totcheck_val =  
              utils::Sq(dcaxysyst[iC]->GetBinContent(iB)) +
              utils::Sq(dcazsyst[iC]->GetBinContent(iB)) +
              utils::Sq(tpcsyst[iC]->GetBinContent(iB)) +
              utils::Sq(pidsyst[iC]->GetBinContent(iB)) +
              utils::Sq(matsyst[iC]->GetBinContent(iB)) +
              utils::Sq(abssyst[iC]->GetBinContent(iB)) +
              utils::Sq(countsyst[iC]->GetBinContent(iB)) +
              utils::Sq(shiftsyst[iC]->GetBinContent(iB)) +
              utils::Sq(itstpcsyst[iC]->GetBinContent(iB)) +
              utils::Sq(siglossyst[iC]->GetBinContent(iB));
          
          totcheck_val = sqrt(totcheck_val);
          totsystcheck[iC]->SetBinContent(iB,totcheck_val);

          double pt_corr =
              utils::Sq(matsyst[iC]->GetBinContent(iB)) +
              utils::Sq(abssyst[iC]->GetBinContent(iB)) +
              utils::Sq(itstpcsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kDCAxy] * utils::Sq(dcaxysyst[iC]->GetBinContent(iB)) +
              kPtCorr[kDCAz] * utils::Sq(dcazsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kTPC] * utils::Sq(tpcsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kPID] * utils::Sq(pidsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kWidth] * utils::Sq(countsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kShift] * utils::Sq(shiftsyst[iC]->GetBinContent(iB));

          pt_corr = sqrt(pt_corr);
          pt_corr_syst[iC]->SetBinContent(iB,pt_corr);

          double pt_uncorr =
              utils::Sq(siglossyst[iC]->GetBinContent(iB)) +
              TMath::Abs(1 - kPtCorr[kDCAxy]) * utils::Sq(dcaxysyst[iC]->GetBinContent(iB)) +
              TMath::Abs(1 - kPtCorr[kDCAz]) * utils::Sq(dcazsyst[iC]->GetBinContent(iB)) +
              TMath::Abs(1 - kPtCorr[kTPC]) * utils::Sq(tpcsyst[iC]->GetBinContent(iB)) +
              TMath::Abs(1 - kPtCorr[kPID]) * utils::Sq(pidsyst[iC]->GetBinContent(iB)) +
              TMath::Abs(1 - kPtCorr[kWidth]) * utils::Sq(countsyst[iC]->GetBinContent(iB)) +
              TMath::Abs(1 - kPtCorr[kShift]) * utils::Sq(shiftsyst[iC]->GetBinContent(iB));
          
          pt_uncorr = sqrt(pt_uncorr);
          pt_uncorr_syst[iC]->SetBinContent(iB,pt_uncorr);

          double totsystptcheck_val = Sq(pt_corr) + Sq(pt_uncorr);
          totsystptcheck_val = sqrt(totsystptcheck_val);
          totsystptcheck[iC]->SetBinContent(iB,totsystptcheck_val);

          double mult_corr = 
              utils::Sq(matsyst[iC]->GetBinContent(iB)) +
              utils::Sq(abssyst[iC]->GetBinContent(iB)) +
              utils::Sq(itstpcsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kDCAxy] * (1 - TMath::Abs(1 - dcaxyR_val) ) * utils::Sq(dcaxysyst[iC]->GetBinContent(iB)) +
              kPtCorr[kDCAz] * (1 - TMath::Abs(1 - dcazR_val) ) *utils::Sq(dcazsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kTPC] * (1 - TMath::Abs(1 - tpcR_val) ) * utils::Sq(tpcsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kPID] * (1 - TMath::Abs(1 - pidR_val) ) * utils::Sq(pidsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kWidth] * (1 - TMath::Abs(1 - countR_val) ) * utils::Sq(countsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kShift] * (1 - TMath::Abs(1 - shiftR_val) ) * utils::Sq(shiftsyst[iC]->GetBinContent(iB));
          
          mult_corr = sqrt(mult_corr);
          mult_corr_syst[iC]->SetBinContent(iB,mult_corr);

          double mult_uncorr = 
              utils::Sq(siglossyst[iC]->GetBinContent(iB)) +
              kPtCorr[kDCAxy] * TMath::Abs(1 - dcaxyR_val) * utils::Sq(dcaxysyst[iC]->GetBinContent(iB)) +
              kPtCorr[kDCAz] * TMath::Abs(1 - dcazR_val) * utils::Sq(dcazsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kTPC] * TMath::Abs(1 - tpcR_val) * utils::Sq(tpcsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kPID] * TMath::Abs(1 - pidR_val) * utils::Sq(pidsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kWidth] * TMath::Abs(1 - countR_val) * utils::Sq(countsyst[iC]->GetBinContent(iB)) +
              kPtCorr[kShift] * TMath::Abs(1 - shiftR_val) * utils::Sq(shiftsyst[iC]->GetBinContent(iB));
          mult_uncorr = sqrt(mult_uncorr);
          mult_uncorr_syst[iC]->SetBinContent(iB,mult_uncorr);

          double totsystmultcheck_val = sqrt(Sq(mult_corr)+Sq(mult_uncorr));
          totsystmultcheck[iC]->SetBinContent(iB,totsystmultcheck_val);
        }
      }



      TCanvas summary("summary","Summary");
      summary.DrawFrame(0.3,0.,4.1,0.5,";#it{p}_{T} (GeV/#it{c}); Systematics uncertainties");
      TLegend leg (0.6,0.56,0.89,0.84);
      leg.SetBorderSize(0);
      dcaxysyst[iC]->SetLineColor(plotting::kSpectraColors[0]);
      dcaxysyst[iC]->SetLineWidth(2);
      dcaxysyst[iC]->Draw("same");
      leg.AddEntry(dcaxysyst[iC],"dcaxy cut","l");
      dcazsyst[iC]->SetLineColor(plotting::kSpectraColors[1]);
      dcazsyst[iC]->SetLineWidth(2);
      dcazsyst[iC]->Draw("same");
      leg.AddEntry(dcazsyst[iC],"dcaz cut","l");
      tpcsyst[iC]->SetLineColor(plotting::kSpectraColors[2]);
      tpcsyst[iC]->SetLineWidth(2);
      tpcsyst[iC]->Draw("same");
      leg.AddEntry(tpcsyst[iC],"tpc cut","l");
      pidsyst[iC]->SetLineColor(plotting::kSpectraColors[3]);
      pidsyst[iC]->SetLineWidth(2);
      pidsyst[iC]->Draw("same");
      leg.AddEntry(pidsyst[iC],"PID","l");
      countsyst[iC]->SetLineColor(plotting::kSpectraColors[4]);
      countsyst[iC]->SetLineWidth(2);
      countsyst[iC]->Draw("histo same");
      leg.AddEntry(countsyst[iC],"Range broadening","l");
      shiftsyst[iC]->SetLineColor(plotting::kSpectraColors[5]);
      shiftsyst[iC]->SetLineWidth(2);
      shiftsyst[iC]->Draw("histo same");
      leg.AddEntry(shiftsyst[iC],"Range shifting","l");
      matsyst[iC]->SetLineColor(plotting::kSpectraColors[6]);
      matsyst[iC]->SetLineWidth(2);
      matsyst[iC]->Draw("histo same");
      leg.AddEntry(matsyst[iC],"Material budget","l");
      abssyst[iC]->SetLineColor(plotting::kSpectraColors[7]);
      abssyst[iC]->SetLineWidth(2);
      abssyst[iC]->Draw("histo same");
      leg.AddEntry(abssyst[iC],"Hadronic interaction","l");
      itstpcsyst[iC]->SetLineColor(plotting::kSpectraColors[8]);
      itstpcsyst[iC]->SetLineWidth(2);
      itstpcsyst[iC]->Draw("histo same");
      leg.AddEntry(itstpcsyst[iC],"ITS-TPC matching","l");
      siglossyst[iC]->SetLineColor(plotting::kSpectraColors[9]);
      siglossyst[iC]->SetLineWidth(2);
      siglossyst[iC]->Draw("histo same");
      leg.AddEntry(siglossyst[iC],"Signal loss","l");
      totsyst[iC]->SetLineColor(plotting::kSpectraColors[10]);
      totsyst[iC]->SetLineWidth(2);
      totsyst[iC]->Draw("histo same");
      leg.AddEntry(totsyst[iC],"Total","l");
      leg.Draw();

      totsyst[iC]->SetLineColor(kMagenta);
      totsyst[iC]->SetLineWidth(2);

      cent_dir->cd();
      cutsyst[iC]->Write();
      dcaxysyst[iC]->Write();
      dcazsyst[iC]->Write();
      tpcsyst[iC]->Write();
      pidsyst[iC]->Write();
      countsyst[iC]->Write();
      shiftsyst[iC]->Write();
      dcaxyR[iC]->Write();
      dcazR[iC]->Write();
      tpcR[iC]->Write();
      pidR[iC]->Write();
      countR[iC]->Write();
      shiftR[iC]->Write();
      abssyst[iC]->Write();
      matsyst[iC]->Write();
      itstpcsyst[iC]->Write();
      siglossyst[iC]->Write();
      pt_uncorr_syst[iC]->Write();
      pt_corr_syst[iC]->Write();
      mult_corr_syst[iC]->Write();
      mult_uncorr_syst[iC]->Write();
      totsyst[iC]->Write();
      totsystcheck[iC]->Write();
      totsystptcheck[iC]->Write();
      totsystmultcheck[iC]->Write();
      summary.Write();

      // if (kPrintFigures) {
      //   summary.SaveAs((kFiguresFolder + "syst" + kLetter[iS] + std::to_string(iC) + ".eps").data());
      //   summary.SaveAs((kMacrosFolder + "syst" + kLetter[iS] + std::to_string(iC) + ".C").data());
      // }
    }
  }
  output_file.Close();

}
