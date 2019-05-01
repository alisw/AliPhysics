#include "src/Common.h"
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
#include <RooRealVar.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooPlot.h>
#include <RooMsgService.h>
#include <RooAbsReal.h>

void RooSec(){

  /// Suppressing the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  RooMsgService::instance().setSilentMode(true);
  gErrorIgnoreLevel=kError; // Suppressing warning outputs

  TFile mc_file(kMCfilename.data());
  TFile data_file(kDataFilename.data());
  string out_path = kBaseOutputDir + "RooSec.root";
  TFile output_file(out_path.data(),"RECREATE");

  TF1 fitModel("fitFrac","1/(1-[0]*exp([1]*x))",0.4,6.);
  fitModel.SetParLimits(0, -100000, 0);
  fitModel.SetParLimits(1, -30, 30);

  gStyle->SetOptFit(1111);

  const int nDCAbins = 34;
  const double dcabins[35] = {
    -1.30,-1.20,-1.10,-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,
    -0.35,-0.30,-0.25,-0.20,-0.15,-0.10,-0.05, 0.00, 0.05, 0.10,
     0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80,
     0.90, 1.00, 1.10, 1.20, 1.30
  };

  for(auto&& list_key : *mc_file.GetListOfKeys()){
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    TTList* mcList = (TTList*)mc_file.Get(list_key->GetName());
    TTList* dtList = (TTList*)data_file.Get(list_key->GetName());

    TH3F *primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTOF");
    TH3F *secondaries = (TH3F*)mcList->FindObject("fMDCASecondaryTOF");
    TH3F *data = (TH3F*)dtList->FindObject("fMDCAxyTOF");

    TDirectory *root = output_file.mkdir(list_key->GetName());
    TDirectory *datadir = root->mkdir("Data");
    TDirectory *primdir = root->mkdir("Primaries");
    TDirectory *secodir = root->mkdir("Secondaries");
    TDirectory *roofitter = root ->mkdir("RooFitter");
    TDirectory *resdir = root->mkdir("Results");

    TAxis *pt = data->GetYaxis();
    auto pt_labels = *(pt->GetXbins());
    TAxis *cen = data->GetXaxis();
    auto cent_labels = *(cen->GetXbins());
    int n_cent_bins = secondaries->GetNbinsX();
    TH1D* hResTFF[kCentLength] = {nullptr};
    roofitter->cd();
    for(int iC=0; iC<kCentLength; iC++){
      hResTFF[iC]= new TH1D(Form("hResTFF_%i",iC),";p_{T} GeV/c;Fraction",pt->GetNbins(),pt->GetXbins()->GetArray());
      roofitter->mkdir(Form("%d",iC));
    }

    for (int iB = pt->FindBin(kPtRangeMatCorrection[0]); iB <= pt->FindBin(kPtRangeMatCorrection[1]); ++iB){

      TH1D *pr_tmp = primaries->ProjectionZ(Form("pr_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* pr = (TH1D*)pr_tmp->Rebin(nDCAbins,Form("pr_%i",iB),dcabins);
      TH1D *sc_tmp = secondaries->ProjectionZ(Form("sc_tmp%i",iB),1,n_cent_bins,iB,iB);
      TH1D* sc = (TH1D*)sc_tmp->Rebin(nDCAbins,Form("sc_%i",iB),dcabins);
      pr->SetTitle(Form("%4.1f < p_{T} #leq %4.1f",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
      sc->SetTitle(Form("%4.1f < p_{T} #leq %4.1f",pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));

      primdir->cd();
      pr->Write();
      secodir->cd();
      sc->Write();

      float low_limit = pr->GetBinLowEdge(1);
      float up_limit = pr->GetBinLowEdge(pr->GetNbinsX()+1);

      std::unique_ptr<RooRealVar> DCAxy(new RooRealVar("x","DCAxy",low_limit,up_limit));
      DCAxy->setRange("interest",-0.12,0.12);
      std::unique_ptr<RooDataHist> roo_pr(new RooDataHist(Form("roo_%s",pr->GetName()),Form("%s",pr->GetTitle()),*DCAxy,pr));
      std::unique_ptr<RooHistPdf> roo_pr_pdf(new RooHistPdf("pr_pdf",Form("%s_pdf",roo_pr->GetTitle()),*DCAxy,*roo_pr));
      std::unique_ptr<RooDataHist> roo_sc(new RooDataHist(Form("roo_%s",sc->GetName()),Form("%s",sc->GetTitle()),*DCAxy,sc));
      std::unique_ptr<RooHistPdf> roo_sc_pdf(new RooHistPdf("sc_pdf",Form("%s_pdf",roo_sc->GetTitle()),*DCAxy,*roo_sc));
      std::unique_ptr<RooRealVar> f(new RooRealVar("f","Fraction of primaries",0.5,0.0000001,1.));
      std::unique_ptr<RooAddPdf> total(new RooAddPdf("total", "Total template",RooArgList(*roo_pr_pdf,*roo_sc_pdf),RooArgList(*f)));

      TH1D* dt[kCentLength] = {nullptr};
      std::unique_ptr<RooDataHist> roo_dt[kCentLength];
      for(int iC=0; iC<kCentLength; iC++){
        datadir->cd();
        TH1D *dt_tmp = data->ProjectionZ(Form("dt_tmp%i",iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iB,iB);
        dt[iC] = (TH1D*)dt_tmp->Rebin(nDCAbins,Form("dt_%i_%i",iB,iC),dcabins);
        dt[iC]->SetTitle(Form("%1.0f - %1.0f %%  %4.1f < p_{T} #leq %4.1f",cen->GetBinLowEdge(kCentBinsArray[iC][0]),cen->GetBinUpEdge(kCentBinsArray[iC][1]),pt->GetBinLowEdge(iB),pt->GetBinLowEdge(iB+1)));
        dt[iC]->Write();
        roo_dt[iC] = make_unique<RooDataHist>(Form("roo_%s",dt[iC]->GetTitle()),Form("roo_%s",dt[iC]->GetTitle()),*DCAxy,dt[iC]);
        RooPlot* frame = DCAxy->frame();
        roo_dt[iC]->plotOn(frame);
        total->fitTo(*roo_dt[iC],RooFit::Verbose(kFALSE),RooFit::PrintEvalErrors(-1),RooFit::PrintLevel(-1));
        roofitter->cd(Form("%d",iC));
        total->plotOn(frame);
        total->plotOn(frame,RooFit::Components("pr_pdf"),RooFit::LineColor(kGreen+3));
        total->plotOn(frame,RooFit::Components("sc_pdf"),RooFit::LineColor(kRed));
        TString iTitle = Form(" %1.0f - %1.0f %% %1.1f #leq #it{p}_{T} < %1.1f GeV/#it{c}", cent_labels[kCentBinsArray[iC][0]-1], cent_labels[kCentBinsArray[iC][1]], pt_labels[iB-1], pt_labels[iB]);
        TString iName = Form("d%i_%i",iC,iB);
        frame->SetTitle(iTitle.Data());
        frame->SetName(iName.Data());
        frame->GetYaxis()->SetTitle(Form("Counts / (%.2f Gev/#it{c}^{2})",frame->GetXaxis()->GetBinWidth(1)));
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
    resdir->cd();
    for(int iC=0; iC<kCentLength; iC++){
      hResTFF[iC]->Fit(&fitModel,"Q");
      hResTFF[iC]->Write();
    }
    root->Write();
  }
}
