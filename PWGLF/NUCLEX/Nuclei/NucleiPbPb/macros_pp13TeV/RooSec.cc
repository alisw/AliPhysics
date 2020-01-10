#include "src/Common.h"
#include "src/Plotting.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TH3F.h>
#include <TAxis.h>
#include <TString.h>
#include <TError.h>
#include <TCanvas.h>
#include <TLine.h>
#include <RooRealVar.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooMsgService.h>
#include <RooAbsReal.h>
#include <RooPolynomial.h>

void RooSec(bool isTPC = false,bool use_antideuterons = false,bool bRebin = false){

  /// Suppressing the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel=kError; // Suppressing warning outputs

  TFile mc_file(kMCfilename.data());
  TFile data_file(kDataFilename.data());
  std::string output_string = (isTPC) ? kSecondariesTPCoutput.data() : kSecondariesOutput.data();
  TFile output_file(output_string.data(),"recreate");

  TF1 fitModel("fitFrac","1/(1+[0]*exp([1]*x))",0.4,6.);
  fitModel.SetParLimits(0, 0, 100000);
  fitModel.SetParLimits(1, -30, 30);

  gStyle->SetOptFit(1111);

  for(auto list_key : *mc_file.GetListOfKeys()){
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find("_MV") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;
    if (string(list_key->GetName()).find("_p_selection") != string::npos) continue;
    if (string(list_key->GetName()).find("_sectemplate") != string::npos) continue;
    if (string(list_key->GetName()).find("_dcaxy0") != string::npos) continue;
    if (string(list_key->GetName()).find("_dcaxy1") != string::npos) continue;

    TTList* mcList = (TTList*)mc_file.Get(list_key->GetName());
    TTList* dtList = (TTList*)data_file.Get(list_key->GetName());

    TH3F *primaries = nullptr;
    TH3F *secondaries = nullptr;
    TH3F *data = nullptr;

    if(isTPC){
      if(use_antideuterons){
        primaries = (TH3F*)dtList->FindObject("fADCAxyTPC");
        Requires(primaries,"Missing primaries (real anti-deuterons) TPC");
      } else {
        primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTPC");
        Requires(primaries,"Missing primaries TPC");
      }
      secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTPC");
      Requires(secondaries,"Missing secondaries TPC");
      data = (TH3F*)dtList->FindObject("fMDCAxyTPC");
      Requires(data,"Missing data TPC");
    } else {
      if(use_antideuterons){
        primaries = (TH3F*)dtList->FindObject("fADCAxyTOF");
        Requires(primaries,"Missing primaries (real anti-deuterons)");
      } else {
        primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTOF");
        Requires(primaries,"Missing primaries");
      }
      secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTOF");
      Requires(secondaries,"Missing secondaries TOF");
      data = (TH3F*)dtList->FindObject("fMDCAxyTOF");
      Requires(data,"Missing data TOF");
    }

    TDirectory *base_dir = output_file.mkdir(list_key->GetName());
    TDirectory *datadir = base_dir->mkdir("Data");
    TDirectory *primdir = base_dir->mkdir("Primaries");
    TDirectory *secodir = base_dir->mkdir("Secondaries");
    TDirectory *roofitter = base_dir->mkdir("RooFitter");
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
      roofitter->mkdir(Form("%d",iC));
      resdir->mkdir(Form("%d",iC));
      ratio2MBdir->mkdir(Form("%d",iC));
    }

    for (int iB = 1; iB <= kNPtBins; ++iB){
      float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
      if (((bin_center < kPtRangeMatCorrectionTPC[0] || bin_center > kPtRangeMatCorrectionTPC[1]) && isTPC) || 
          ((bin_center < kPtRangeMatCorrectionTOF[0] || bin_center > kPtRangeMatCorrectionTOF[1]) && !isTPC)) continue;
      int iBin = data->GetYaxis()->FindBin(bin_center);
      int iBin1_prim = primaries->GetYaxis()->FindBin(kPtBins[iB-1]+0.005);
      int iBin2_prim = primaries->GetYaxis()->FindBin(kPtBins[iB]-0.005);
      int iBin1_sec = secondaries->GetYaxis()->FindBin(kPtBins[iB-1]+0.005);
      int iBin2_sec = secondaries->GetYaxis()->FindBin(kPtBins[iB]-0.005);
      //
      TH1D* pr[kCentLength] = {nullptr};
      TH1D* pr_tmp = nullptr;
      for(int iC=kCentLength; iC--;){
        if(use_antideuterons){
          pr_tmp = (TH1D*)primaries->ProjectionZ(Form("pr_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin,iBin);
        } else {
          pr_tmp = (TH1D*)primaries->ProjectionZ(Form("pr_tmp_%i_%i",iC,iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin1_prim,iBin2_prim);
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
      //
      TH1D *sc = nullptr;
      TH1D* sc_tmp= (TH1D*)secondaries->ProjectionZ(Form("sc_tmp_%i",iB),kCentBinsArray[kCentLength-1][0],kCentBinsArray[kCentLength-1][1],iBin1_sec,iBin2_sec);
      if(bRebin){
        sc = (TH1D*)sc_tmp->Rebin(kNDCAbins,Form("sc_%i",iB),kDCAbins);
      } else {
        sc = (TH1D*)sc_tmp->Clone(Form("sc_%i",iB));
      }
      sc->SetTitle(Form("%4.1f < p_{T} #leq %4.1f (GeV/#it{c});counts; DCA_{xy} (cm)",kPtBins[iB-1],kPtBins[iB]));
      secodir->cd();
      sc->Write();
      //
      TH1D* dt[kCentLength] = {nullptr};
      TH1D* dt_tmp = nullptr;
      for(int iC=kCentLength; iC--;){
        dt_tmp = (TH1D*)data->ProjectionZ(Form("dt_tmp_%i",iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin,iBin,"e");
        if(bRebin){
          dt[iC] = (TH1D*)dt_tmp->Rebin(kNDCAbins,Form("dt_%i_%i",iC,iB),kDCAbins);
        } else {
          dt[iC] = (TH1D*)dt_tmp->Clone(Form("dt_%i_%i",iC,iB));
        }
        dt[iC]->SetTitle(Form("Mult: %1.0f - %1.0f %% , %4.1f < p_{T} #leq %4.1f (GeV/#it{c})",kCentLabels[iC][0],kCentLabels[iC][1],kPtBins[iB-1],kPtBins[iB]));
        datadir->cd(Form("%d",iC));
        dt[iC]->Write();
      }
      //
      float low_limit = pr[0]->GetBinLowEdge(1);
      float up_limit = pr[0]->GetBinLowEdge(pr[0]->GetNbinsX()+1);
      std::unique_ptr<RooRealVar> DCAxy(new RooRealVar("x","DCAxy",low_limit,up_limit));
      DCAxy->setRange("interest",-0.12,0.12);     
      std::unique_ptr<RooDataHist> roo_sc(new RooDataHist(Form("roo_%s",sc->GetName()),Form("%s",sc->GetTitle()),*DCAxy,sc));
      std::unique_ptr<RooHistPdf> roo_sc_pdf(new RooHistPdf("sc_pdf",Form("%s_pdf",roo_sc->GetTitle()),*DCAxy,*roo_sc));
      //std::unique_ptr<RooPolynomial> roo_sc_pdf(new RooPolynomial("sc_pdf",Form("%s_pdf",roo_sc->GetTitle()),*DCAxy,RooArgList()));
      std::unique_ptr<RooDataHist> roo_dt[kCentLength];
      std::unique_ptr<RooDataHist> roo_pr[kCentLength];
      for(int iC=0; iC<kCentLength; iC++){
        
        roo_pr[iC]=utils::make_unique<RooDataHist>(Form("roo_%s",pr[iC]->GetName()),Form("%s",pr[iC]->GetTitle()),*DCAxy,pr[iC]);
        std::unique_ptr<RooHistPdf> roo_pr_pdf(new RooHistPdf("pr_pdf",Form("%s_pdf",roo_pr[iC]->GetTitle()),*DCAxy,*roo_pr[iC]));        
        std::unique_ptr<RooRealVar> f(new RooRealVar("f","Fraction of primaries",0.5,0.0000001,1.));
        std::unique_ptr<RooAddPdf> total(new RooAddPdf("total", "Total template",RooArgList(*roo_pr_pdf,*roo_sc_pdf),RooArgList(*f)));
        //
        roo_dt[iC] = utils::make_unique<RooDataHist>(Form("roo_%s",dt[iC]->GetTitle()),Form("roo_%s",dt[iC]->GetTitle()),*DCAxy,dt[iC]);
        RooPlot* frame = DCAxy->frame();
        roo_dt[iC]->plotOn(frame);
        total->fitTo(*roo_dt[iC],RooFit::Verbose(kFALSE),RooFit::PrintEvalErrors(-1),RooFit::PrintLevel(-1));
        roofitter->cd(Form("%d",iC));
        total->plotOn(frame);
        total->plotOn(frame,RooFit::Components("pr_pdf"),RooFit::LineColor(kGreen+3));
        total->plotOn(frame,RooFit::Components("sc_pdf"),RooFit::LineColor(kRed));
        TString iTitle = Form(" %1.0f - %1.0f %% %1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", kCentLabels[iC][0], kCentLabels[iC][1], kPtBins[iB-1], kPtBins[iB]);
        TString iName = Form("d%i_%i",iC,iB);
        frame->SetTitle(iTitle.Data());
        frame->SetName(iName.Data());
        frame->GetYaxis()->SetTitle(Form("Counts / (%.2f cm)",frame->GetXaxis()->GetBinWidth(1)));
        frame->Write();
        /// Determining the primary fraction
        RooAbsReal* prim_int = roo_pr_pdf->createIntegral(*DCAxy,RooFit::NormSet(*DCAxy),RooFit::Range("interest"));
        RooAbsReal* tot_int = total->createIntegral(*DCAxy,RooFit::NormSet(*DCAxy),RooFit::Range("interest"));
        /// Writing on the final histogram
        float ratio = f->getVal()*prim_int->getVal()/tot_int->getVal();
        hResTFF[iC]->SetBinContent(iB, ratio);
        hResTFF[iC]->SetBinError(iB, std::sqrt(ratio * (1. - ratio) / (tot_int->getVal()*dt[iC]->GetEntries())));
      }
    }
    for(int iC=0; iC<kCentLength; iC++){
      resdir->cd(Form("%d",iC));
      hResTFF[iC]->Fit(&fitModel,"Q");
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
