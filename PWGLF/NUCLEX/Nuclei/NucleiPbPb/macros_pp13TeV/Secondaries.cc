#include "src/Common.h"
#include "src/Utils.h"
#include "src/Plotting.h"
using namespace utils;

#include <TDirectory.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFractionFitter.h>
#include <TH1D.h>
#include <TH3F.h>
#include <TList.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TError.h>
#include <TLine.h>
#include <TMath.h>

void Secondaries(bool isTPC = false, bool use_antideuterons = true, bool bRebin = false, bool use_sec_loose = false) {
  /// Taking all the histograms from the MC and data files
  TFile mc_file(kMCfilename.data());
  const char* mc_name_MB = (kUseIntegratedForMB) ? kMCfilename.data() : kMCfilenameMB.data();
  TFile mc_file_MB(mc_name_MB);
  TFile data_file(kDataFilename.data());
  const char* data_name_MB = (kUseIntegratedForMB) ? kDataFilename.data() : kDataFilenameMB.data();
  TFile data_file_MB(data_name_MB);
  std::string output_string = (isTPC) ? kSecondariesTPCoutput.data() : kSecondariesOutput.data();
  TFile output_file(output_string.data(),"recreate");

  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);
  gErrorIgnoreLevel=kError;

  TH1D* sc_arr[kNPtBins] = {nullptr};
  if(use_sec_loose){
    TTList* sec_list = (TTList*)mc_file_MB.Get("nuclei_deuteronsMB_dcaz3");
    TH3F *secondaries_gen = (isTPC) ? (TH3F*)sec_list->FindObject("fMDCASecondaryTPC") : (TH3F*)sec_list->FindObject("fMDCASecondaryTOF");
    Requires(secondaries_gen,"Missing secondaries");
    for (int iB = 1; iB <= kNPtBins; ++iB) {
      float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
      if (((bin_center < kPtRangeMatCorrectionTPC[0] || bin_center > kPtRangeMatCorrectionTPC[1]) && isTPC) || 
          ((bin_center < kPtRangeMatCorrectionTOF[0] || bin_center > kPtRangeMatCorrectionTOF[1]) && !isTPC)) continue;
      int iBin1 = secondaries_gen->GetYaxis()->FindBin(kPtBins[iB-1]+0.005);
      int iBin2 = secondaries_gen->GetYaxis()->FindBin(kPtBins[iB]-0.005);
      int first_bin = (kUseIntegratedForMB) ? 2 : 1;
      TH1D* sc_tmp = (TH1D*)secondaries_gen->ProjectionZ(Form("sc_tmp_%i",iB),first_bin,kCentBinsArray[kCentLength-1][1],iBin1,iBin2);
      sc_arr[iB] = (TH1D*)sc_tmp->Rebin(kNDCAbins,Form("sc_%i",iB),kDCAbins);
      sc_arr[iB]->SetTitle(Form("%4.1f < p_{T} #leq %4.1f (GeV/#it{c})",kPtBins[iB-1],kPtBins[iB]));
    }
  }

  TObjArray obj(2);
  for (auto list_key : *data_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find("_MV") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;
    if (string(list_key->GetName()).find("_p_selection") != string::npos) continue;
    if (string(list_key->GetName()).find("_sectemplate") != string::npos) continue;

    string list_name = list_key->GetName();
    TTList* mcList = (TTList*)mc_file.Get(list_name.data());
    TTList* dtList = (TTList*)data_file.Get(list_name.data());
    string list_name_MB = list_name; 
    if(!kUseIntegratedForMB) list_name_MB.insert(kFilterListNames.size()-1,"MB");
    TTList* mcListMB = (TTList*)mc_file_MB.Get(list_name_MB.data());
    Requires(mcListMB, Form("MC: %s",list_name_MB.data()));
    TTList* dtListMB = (TTList*)data_file_MB.Get(list_name_MB.data());
    Requires(dtListMB, Form("data: %s",list_name_MB.data()));

    TH3F* primaries = nullptr;
    TH3F* secondaries = nullptr;
    TH3F* data = nullptr;

    TH3F* primaries_MB = nullptr;
    TH3F* secondaries_MB = nullptr;
    TH3F* data_MB = nullptr;

    if(isTPC){
      if(use_antideuterons){
        primaries = (TH3F*)dtList->FindObject("fADCAxyTPC");
        Requires(primaries,"Missing primaries (real anti-deuterons) TPC");
        primaries_MB = (TH3F*)dtListMB->FindObject("fADCAxyTPC");
        Requires(primaries_MB,"Missing primaries (real anti-deuterons) TPC (MB)");
      } else {
        primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTPC");
        Requires(primaries,"Missing primaries TPC");
        primaries_MB = (TH3F*)mcListMB->FindObject("fMDCAPrimaryTPC");
        Requires(primaries_MB,"Missing primaries TPC (MB)");
      }
      if(!use_sec_loose){
        secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTPC");
        Requires(secondaries,"Missing secondaries for TPC");
        secondaries_MB = (TH3F*)mcListMB->FindObject("fMDCASecondaryTPC");
        Requires(secondaries_MB,"Missing secondaries for TPC (MB)");
      }
      data = (TH3F*)dtList->FindObject("fMDCAxyTPC");
      Requires(data,"Missing data for TPC");
      data_MB = (TH3F*)dtListMB->FindObject("fMDCAxyTPC");
      Requires(data_MB,"Missing data for TPC (MB)");
    } else {
      if(use_antideuterons){
        primaries = (TH3F*)dtList->FindObject("fADCAxyTOF");
        Requires(primaries,"Missing primaries (real anti-deuterons) TOF");
        primaries_MB = (TH3F*)dtListMB->FindObject("fADCAxyTOF");
        Requires(primaries_MB,"Missing primaries (real anti-deuterons) TOF (MB)");
      } else {
        primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTOF");
        Requires(primaries,"Missing primaries TOF");
        primaries_MB = (TH3F*)mcListMB->FindObject("fMDCAPrimaryTOF");
        Requires(primaries_MB,"Missing primaries TOF (MB)");
      }
      if(!use_sec_loose){
        secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTOF");
        Requires(secondaries,"Missing secondaries");
        secondaries_MB = (TH3F*)mcListMB->FindObject("fMDCASecondaryTOF");
        Requires(secondaries_MB,"Missing secondaries for TOF (MB)");
      }
      data = (TH3F*)dtList->FindObject("fMDCAxyTOF");
      Requires(data,"Missing data");
      data_MB = (TH3F*)dtListMB->FindObject("fMDCAxyTOF");
      Requires(data_MB,"Missing data for TOF (MB)");
    }

    TDirectory *base_dir = output_file.mkdir(list_key->GetName());
    TDirectory *datadir = base_dir->mkdir("Data");
    TDirectory *primdir = base_dir->mkdir("Primaries");
    TDirectory *secodir = base_dir->mkdir("Secondaries");
    TDirectory *fitdir = base_dir->mkdir("Fit");
    TDirectory *tffdir = base_dir->mkdir("TFractionFitter");
    TDirectory *resdir = base_dir->mkdir("Results");
    TDirectory *ratio2MBdir = base_dir->mkdir("RatioToMB");

    TH1D* hResTFF[kCentLength] = {nullptr};
    TH1D* hRatioToMB[kCentLength] = {nullptr};
    TH1F* ratio_fit_data[kCentLength] = {nullptr};
    for(int iC=0; iC<kCentLength; iC++){
      hResTFF[iC]= new TH1D(Form("hResTFF_%i",iC),";p_{T} GeV/c;Fraction",kNPtBins,kPtBins);
      hResTFF[iC]->GetXaxis()->SetRangeUser(0.6,4.);
      ratio_fit_data[iC] = new TH1F(Form("ratio_fit_data_%d",iC),Form("%2.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});data/fit",kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      hRatioToMB[iC] = new TH1D(Form("hRatioToMB_%d",iC),Form("%2.0f - %2.0f %%;#it{p}_{T} (GeV/#it{c});Ratio to MB",kCentLabels[iC][0],kCentLabels[iC][1]),kNPtBins,kPtBins);
      primdir->mkdir(Form("%d",iC));
      datadir->mkdir(Form("%d",iC));
      tffdir->mkdir(Form("%d",iC));
      resdir->mkdir(Form("%d",iC));
      fitdir->mkdir(Form("%d",iC));
      ratio2MBdir->mkdir(Form("%d",iC));
    }

    for (int iB = 1; iB <= kNPtBins; ++iB) {
      float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
      if (((bin_center < kPtRangeMatCorrectionTPC[0] || bin_center > kPtRangeMatCorrectionTPC[1]) && isTPC) || 
          ((bin_center < kPtRangeMatCorrectionTOF[0] || bin_center > kPtRangeMatCorrectionTOF[1]) && !isTPC)) continue;
      int iBin = data->GetYaxis()->FindBin(bin_center);
      int iBin_MB = data_MB->GetYaxis()->FindBin(bin_center);
      int iBin1_prim = primaries->GetYaxis()->FindBin(kPtBins[iB-1]+0.005);
      int iBin2_prim = primaries->GetYaxis()->FindBin(kPtBins[iB]-0.005);
      int iBin1_prim_MB = primaries_MB->GetYaxis()->FindBin(kPtBins[iB-1]+0.005);
      int iBin2_prim_MB = primaries_MB->GetYaxis()->FindBin(kPtBins[iB]-0.005);
      int iBin1_sec = secondaries->GetYaxis()->FindBin(kPtBins[iB-1]+0.005);
      int iBin2_sec = secondaries->GetYaxis()->FindBin(kPtBins[iB]-0.005);
      int iBin1_sec_MB = secondaries_MB->GetYaxis()->FindBin(kPtBins[iB-1]+0.005);
      int iBin2_sec_MB = secondaries_MB->GetYaxis()->FindBin(kPtBins[iB]-0.005);
      TH1D* pr[kCentLength] = {nullptr};
      TH1D* pr_tmp = nullptr;
      for(int iC=kCentLength; iC--;){
        if(use_antideuterons){
          pr_tmp = (iC==9 && !kUseIntegratedForMB) ? (TH1D*) primaries_MB->ProjectionZ(Form("pr_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin_MB,iBin_MB) : (TH1D*)primaries->ProjectionZ(Form("pr_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin,iBin);
        } else {
          pr_tmp = (iC==9 && !kUseIntegratedForMB) ? (TH1D*)primaries_MB->ProjectionZ(Form("pr_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin1_prim_MB,iBin2_prim_MB) : (TH1D*)primaries->ProjectionZ(Form("pr_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin1_prim,iBin2_prim);
        }
        if(bRebin){
          pr[iC] = (TH1D*)pr_tmp->Rebin(kNDCAbins,Form("pr_%i_%i",iC,iB),kDCAbins);
        } else {
          pr[iC] = (TH1D*)pr_tmp->Clone(Form("pr_%i_%i",iC,iB));
        }
        pr[iC]->SetTitle(Form("Mult: %1.0f - %1.0f %% , %4.1f < p_{T} #leq %4.1f (GeV/#it{c})",kCentLabels[iC][0],kCentLabels[iC][1],kPtBins[iB-1],kPtBins[iB]));
        primdir->cd(Form("%d",iC));
        pr[iC]->Write();
      }
      TH1D *sc = nullptr;
      if(!use_sec_loose){
        int first_bin = (kUseIntegratedForMB) ? 2 : 1;
        TH1D* sc_tmp = (TH1D*)secondaries_MB->ProjectionZ(Form("sc_tmp_%i",iB),first_bin,kCentBinsArray[kCentLength-1][1],iBin1_sec_MB,iBin2_sec_MB);
        if(bRebin){
          sc = (TH1D*)sc_tmp->Rebin(kNDCAbins,Form("sc_%i",iB),kDCAbins);
        } else {
          sc = (TH1D*)sc_tmp->Clone(Form("sc_%i",iB));
        }
        sc->SetTitle(Form("%4.1f < p_{T} #leq %4.1f (GeV/#it{c});counts; DCA_{xy} (cm)",kPtBins[iB-1],kPtBins[iB]));
        secodir->cd();
        sc->Write();
      }
      datadir->cd();
      TH1D* dt[kCentLength] = {nullptr};
      TH1D* dt_tmp = nullptr;
      for(int iC=kCentLength; iC--;){
        dt_tmp = (iC==9 && !kUseIntegratedForMB)? (TH1D*)data->ProjectionZ(Form("dt_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin_MB,iBin_MB,"e") : (TH1D*)data->ProjectionZ(Form("dt_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin,iBin,"e");
        if(bRebin){
          dt[iC] = (TH1D*)dt_tmp->Rebin(kNDCAbins,Form("dt_%i_%i",iC,iB),kDCAbins);
        } else {
          dt[iC] = (TH1D*)dt_tmp->Clone(Form("dt_%i_%i",iC,iB));
        }
        dt[iC]->SetTitle(Form("Mult: %1.0f - %1.0f %% , %4.1f < p_{T} #leq %4.1f (GeV/#it{c})",kCentLabels[iC][0],kCentLabels[iC][1],kPtBins[iB-1],kPtBins[iB]));
        datadir->cd(Form("%i",iC));
        dt[iC]->Write();
      }
      for(int iC=kCentLength; iC--;){
        obj.Add(pr[iC]);
        if(use_sec_loose){
          obj.Add(sc_arr[iB]);
        } else {
          obj.Add(sc);
        }
        TFractionFitter fitter(dt[iC],&obj,"Q");
        fitter.Constrain(0, 0., 1.);
        fitter.Constrain(1, 0., 1.);
        Double_t yieldPri = 0., yieldSec = 1., errorPri = 0., errorSec = 0.;
        char input = 'n';
        TVirtualFitter::SetMaxIterations(10000);
        Int_t result = fitter.Fit();
        if (result == 0) {
          TH1F* hPrimary = (TH1F*)fitter.GetMCPrediction(0);
          TH1F* hSecondary = (TH1F*)fitter.GetMCPrediction(1);
          fitter.GetResult(0, yieldPri, errorPri);
          fitter.GetResult(1, yieldSec, errorSec);
          TH1F* hFit = (TH1F*)fitter.GetPlot();
          Float_t dataIntegral = dt[iC]->Integral();
          hFit->SetLineColor(kGreen + 3);
          hFit->SetLineWidth(3);
          hSecondary->Scale(dataIntegral * yieldSec / hSecondary->Integral());
          hPrimary->Scale(dataIntegral * yieldPri / hPrimary->Integral());
          dt[iC]->SetMarkerStyle(20);
          dt[iC]->SetMarkerSize(0.5);
          dt[iC]->SetMarkerColor(kBlack);

          TCanvas cv(Form("tff_%i_%i",iB,iC),Form("TFractionFitter_%i_%i",iB,iC));
          cv.cd();
          dt[iC]->Scale(1., "width");
          dt[iC]->DrawCopy("PE");

          float sec_cut = 0.12;  
          if (string(list_key->GetName()).find("_dcaxy0") != string::npos) sec_cut = 0.10;
          if (string(list_key->GetName()).find("_dcaxy1") != string::npos) sec_cut = 0.14;
      
          float tot_integral = hFit->Integral(hFit->GetXaxis()->FindBin(-1.*sec_cut+0.005),hFit->GetXaxis()->FindBin(sec_cut-0.005));
          float prim_integral = hPrimary->Integral(hPrimary->GetXaxis()->FindBin(-1.*sec_cut+0.005),hPrimary->GetXaxis()->FindBin(sec_cut-0.005));

          float ratio = prim_integral/tot_integral;

          hFit->Scale(1., "width");
          hFit->DrawCopy("hist same");
          fitdir->cd(Form("%i",iC));
          hFit->Write(Form("fit_%i_%i",iC,iB));
          hSecondary->SetLineColor(kRed);
          hPrimary->SetLineColor(kBlue);
          hPrimary->SetLineWidth(2);
          hSecondary->SetLineWidth(2);
          hSecondary->Scale(1., "width");
          hSecondary->DrawCopy("hist same");
          hPrimary->Scale(1., "width");
          hPrimary->DrawCopy("hist same");
          fitdir->cd(Form("%d",iC));
          hPrimary->Write(Form("pr_%i_%i",iC,iB));
          fitdir->cd(Form("%d",iC));
          hSecondary->Write(Form("sc_%i_%i",iC,iB));
          tffdir->cd(Form("%d",iC));
          TLegend leg (0.6,0.56,0.89,0.84);
          leg.SetBorderSize(0.);
          leg.AddEntry(dt[iC],"Data","pe");
          leg.AddEntry(hPrimary,"Primaries","l");
          leg.AddEntry(hSecondary,"Secondaries","l");
          leg.AddEntry(hFit,"TFF","l");
          leg.AddEntry((TObject*)nullptr,Form("#chi^{2}/NDF = %.2f/%d",fitter.GetChisquare(),fitter.GetNDF()),"");
          leg.Draw();
          cv.Write();
          hResTFF[iC]->SetBinContent(iB, ratio);
          hResTFF[iC]->SetBinError(iB, std::sqrt(ratio * (1. - ratio) / tot_integral));
        } else {
          std::cout << list_key->GetName() << Form(", %4.1f < p_{T} #leq %4.1f (GeV/#it{c}), iC: %d ==> ",kPtBins[iB-1],kPtBins[iB],iC);
          std::cout << " the TFF does not converge." << std::endl;
        }
        delete dt[iC];
        obj.Remove(pr[iC]);
        if(use_sec_loose){
          obj.Remove(sc_arr[iB]);
        } else {
          obj.Remove(sc);
        }
        delete pr[iC];
      }
      secodir->cd();
      if(use_sec_loose){
        delete sc_arr[iB];
      }

    }
    for(int iC=0; iC<kCentLength; iC++){
      resdir->cd(Form("%d",iC));
      /// Building the function used to fit the primary fraction distribution
      //float left_fit_range = (isTPC) ? kPtRangeMatCorrectionTPC[0] : kPtRangeMatCorrectionTOF[0];
      //float right_fit_range = (isTPC) ? kPtRangeMatCorrectionTPC[1] : kPtRangeMatCorrectionTOF[1];
      TF1 fitModel("fitFrac","1/(1+[0]*exp([1]*x))",0.6,4);
      fitModel.SetParLimits(0, 0., 30.);
      fitModel.SetParLimits(1, -5., -2.);
      hResTFF[iC]->Fit(&fitModel,"RQ");
      hResTFF[iC]->Write();
      for(int iB=1; iB<=kNPtBins; iB++){
        float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
        if (((bin_center < kPtRangeMatCorrectionTPC[0] || bin_center > kPtRangeMatCorrectionTPC[1]) && isTPC) || 
          ((bin_center < kPtRangeMatCorrectionTOF[0] || bin_center > kPtRangeMatCorrectionTOF[1]) && !isTPC)) continue;
        float data_val = hResTFF[iC]->GetBinContent(iB);
        float fit_val = fitModel.Eval(hResTFF[iC]->GetBinCenter(iB));
        float data_err = hResTFF[iC]->GetBinError(iB);
        ratio_fit_data[iC]->SetBinContent(iB,data_val/fit_val);
        ratio_fit_data[iC]->SetBinError(iB,data_err/fit_val);
      }
      ratio_fit_data[iC]->GetXaxis()->SetRangeUser(0.6,4.);
      ratio_fit_data[iC]->GetYaxis()->SetRangeUser(0.8,1.2);
      ratio_fit_data[iC]->Write();
      TCanvas* cV = new TCanvas(Form("RatioFitData_%d",iC),Form("cRatioFitData_%d",iC));
      cV->Divide(2);
      cV->cd(1);
      hResTFF[iC]->GetYaxis()->SetRangeUser(0.4,1.2);
      hResTFF[iC]->Draw();
      cV->cd(2);
      ratio_fit_data[iC]->GetYaxis()->SetTitleOffset(1.2);
      ratio_fit_data[iC]->Draw();
      TLine *line_one = new TLine(0.5,1.,2.,1.);
      line_one->SetLineColor(kBlack);
      line_one->SetLineStyle(kDashed);
      line_one->Draw();
      cV->Write();
      /// Ratio between the primary fraction obtained for each multiplicity class and the integrated multiplicity
      for(int iB=1; iB<=kNPtBins; iB++){
        float num = hResTFF[iC]->GetBinContent(iB);
        float num_err = hResTFF[iC]->GetBinError(iB);
        float den = hResTFF[9]->GetBinContent(iB);
        float den_err = hResTFF[9]->GetBinError(iB);
        if(TMath::Abs(num) < 1e-4 || TMath::Abs(den) < 1e-4) continue;
        float ratio = num/den;
        float ratio_err = ratio * TMath::Sqrt((num_err/num)*(num_err/num)+(den_err/den)*(den_err/den));
        hRatioToMB[iC]->SetBinContent(iB,ratio);
        hRatioToMB[iC]->SetBinError(iB,ratio_err);   
      }
      plotting::SetHistStyle(hRatioToMB[iC], plotting::kSpectraColors[iC]);
      hRatioToMB[iC]->GetXaxis()->SetRangeUser(0.6,2);
      hRatioToMB[iC]->GetYaxis()->SetRangeUser(0.8,1.2);
      ratio2MBdir->cd(Form("%d",iC));
      hRatioToMB[iC]->Write();
      TCanvas cRatioToMB(Form("cRatioToMB_%d",iC),Form("cRatioToMB_%d",iC));
      hRatioToMB[iC]->Draw("PE");
      TLine *line_toMB = new TLine(0.6,1.,2.,1.);
      line_toMB->SetLineColor(kBlack);
      line_toMB->SetLineStyle(kDashed);
      line_toMB->Draw();
      cRatioToMB.Write();
    }
  }
}
