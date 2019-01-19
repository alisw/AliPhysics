#if !defined (__CINT__) || defined (__CLING__)

#include <Riostream.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TString.h>
#include <TList.h>
#include <TMath.h>
#include <TGaxis.h>

#include "AliAnalysisTaskSEDmesonPIDSysProp.h"

#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief macro for estimation of PID Systematic uncertainty propagated from the single track to the D mesons   //
// \main function: EstimateDmesonPIDsyst                                                                        //
// \author: A. M. Barbano, anastasia.maria.barbano@cern.ch                                                      //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                                                    //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

//_____________________________________________________
//METHOD PROTOTYPES
int EstimateDmesonPIDsyst(int ch = AliAnalysisTaskSEDmesonPIDSysProp::kDstoKKpi,
                          int pid = AliAnalysisTaskSEDmesonPIDSysProp::kStrongPID,
                          TString infilename="LHC18a4a2_prop.root",
                          TString suffix="ppMB_kINT7",
                          double ptmin=2.,
                          double ptmax=24.);
void SetStyle();

//_____________________________________________________
//METHOD IMPLEMENTATIONS
int EstimateDmesonPIDsyst(int ch, int pid, TString infilename, TString suffix, double ptmin, double ptmax) {

  SetStyle();

  TString mesonname = "";
  if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kD0toKpi) mesonname = "D0";
  else if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kDplustoKpipi) mesonname = "Dplus";
  else if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kDstartoKpipi) mesonname = "Dstar";
  else if(ch == AliAnalysisTaskSEDmesonPIDSysProp::kDstoKKpi) mesonname = "Ds";

  TString PIDname = "";
  double ymax = 0.05;
  if(pid == AliAnalysisTaskSEDmesonPIDSysProp::kConservativePID) PIDname = "conservativePID";
  else if(pid == AliAnalysisTaskSEDmesonPIDSysProp::kStrongPID) {
    PIDname = "strongPID";
    ymax = 0.1;
  }
  else if(pid == AliAnalysisTaskSEDmesonPIDSysProp::knSigmaPID) PIDname = "nSigmaPID";

  TFile* infile = TFile::Open(infilename.Data());
  if(!infile || !infile->IsOpen()) return 1;
  TDirectoryFile* indir = (TDirectoryFile*)infile->Get(Form("PWGHF_D2H_PIDeffsyst_%s_%s_%s",mesonname.Data(),PIDname.Data(),suffix.Data()));
  cout << Form("PWGHF_D2H_PIDeffsyst_%s_%s_%s",mesonname.Data(),PIDname.Data(),suffix.Data()) << endl;
  if(!indir) {
    cerr << "TDirectoryFile not found, please check the name. Exit" << endl;
    return 2;
  }
  TList* inlist = (TList*)indir->Get(Form("systUnc_%s_%s_%s",mesonname.Data(),PIDname.Data(),suffix.Data()));
  if(!inlist) {
    cerr << "TList not found, please check the name. Exit" << endl;
    return 3;
  }
  TH2F* hSystVPt = (TH2F*)inlist->FindObject("fHistSystPIDEffD");
  TH2F* hPtCorr = (TH2F*)inlist->FindObject("fHistPtDauVsD");
  hSystVPt->SetDirectory(0);
  hPtCorr->SetDirectory(0);
  infile->Close();

  TH1F* hSyst = (TH1F*)hSystVPt->ProjectionX("hSyst");
  hSyst->SetLineWidth(2);
  hSyst->SetLineColor(kRed);
  hSyst->SetMarkerColor(kRed);
  hSyst->SetMarkerStyle(kFullCircle);
  for(int iPt=0; iPt<hSyst->GetNbinsX(); iPt++) {
    TH1F* hTmp = (TH1F*)hSystVPt->ProjectionY("hTmp",iPt+1,iPt+1);
    hSyst->SetBinContent(iPt+1,hTmp->GetMean());
    hSyst->SetBinError(iPt+1,hTmp->GetMeanError());
  }
  
  hSystVPt->GetYaxis()->SetTitleOffset(1.4);
  hSystVPt->GetXaxis()->SetTitleOffset(1.2);
  hSystVPt->GetYaxis()->SetTitleSize(0.045);
  hSystVPt->GetXaxis()->SetTitleSize(0.045);
  hSystVPt->GetYaxis()->SetLabelSize(0.04);
  hSystVPt->GetXaxis()->SetLabelSize(0.04);
  hPtCorr->GetYaxis()->SetTitleOffset(1.4);
  hPtCorr->GetXaxis()->SetTitleOffset(1.2);
  hPtCorr->GetYaxis()->SetTitleSize(0.045);
  hPtCorr->GetXaxis()->SetTitleSize(0.045);
  hPtCorr->GetYaxis()->SetLabelSize(0.04);
  hPtCorr->GetXaxis()->SetLabelSize(0.04);

  TCanvas* cSyst = new TCanvas("cSyst","",1000,500);
  cSyst->Divide(2,1);
  cSyst->cd(1)->DrawFrame(ptmin,0.,ptmax,ptmax,";#it{p}_{T}^{D} (GeV/#it{c});#it{p}_{T}^{daugh} (GeV/#it{c})");
  cSyst->cd(1)->SetLogz();
  hPtCorr->Draw("same colz");
  cSyst->cd(2)->DrawFrame(ptmin,0.,ptmax,ymax,";#it{p}_{T}^{D} (GeV/#it{c});relative systematic uncertainty");
  cSyst->cd(2)->SetLogz();
  hSystVPt->Draw("colz same");
  hSyst->Draw("same");

  TFile outfile(Form("%s_%s_PIDsyst.root",mesonname.Data(),PIDname.Data()),"recreate");
  hPtCorr->Write();
  hSystVPt->Write();
  hSyst->Write();
  outfile.Close();
  
  cSyst->SaveAs(Form("%s_%s_PIDsyst.pdf",mesonname.Data(),PIDname.Data()));
  return 0;
}

//_____________________________________________________
void SetStyle() {
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleSize(0.045,"xyz");
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetOptTitle(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  
  TGaxis::SetMaxDigits(3);
}
