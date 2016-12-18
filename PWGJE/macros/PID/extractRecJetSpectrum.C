#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TStyle.h"

#include "THnSparseDefinitions.h"


//______________________________________________________
TCanvas* createCanvas(TString name, TString title)
{
  TCanvas *canv = new TCanvas(name.Data(), title.Data(), 1824, 168, 700, 500);
  canv->Range(-23.68421,0.311344,205.2632,7.009505);
  canv->SetFillColor(10);
  canv->SetBorderMode(0);
  canv->SetBorderSize(2);
  canv->SetLogy();
  canv->SetGridx(0);
  canv->SetGridy(1);
  canv->SetTickx(1);
  canv->SetTicky(1);
  canv->SetLeftMargin(0.115);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.017);
  canv->SetBottomMargin(0.13);
  canv->SetFrameFillColor(0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderMode(0);
  
  return canv;
}


//______________________________________________________
void setupHisto(TH1D* h)
{
  if (!h)
    return;
  
  h->SetStats(0);
  h->SetFillColor(15);
  h->SetLineWidth(2.);
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle("#it{p}_{T}^{jet, ch} (GeV/#it{c})");
  h->GetXaxis()->SetRange(1,11);
  h->GetXaxis()->SetLabelFont(62);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetTitleFont(62);
  h->GetYaxis()->SetLabelFont(62);
  h->GetYaxis()->SetTitle("#it{N}_{Jets, rec}");
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleFont(62);
  h->GetZaxis()->SetLabelSize(0.05);
  h->GetZaxis()->SetTitleSize(0.06);
}


//______________________________________________________
Int_t extractRecJetSpectrum(TString path = "finalCuts/pp/7TeV/10d_10e_merged.pass2/finalisedSplines/finalMapsAndTail/Jets/nclCut/0differentTrackFilterBits/528",
                            TString fileName = "bhess_PID_Jets_results_LLFit__Pt_2_reg1_regFac1.00_noMuons_idSpectra_jetPt5.0_10.0.root")
{
  TString filePathName = Form("%s/%s", path.Data(), fileName.Data());
  TFile* f = TFile::Open(filePathName.Data(), "READ");
  if (!f) {
    printf("Failed to open file \"%s\"!\n", filePathName.Data());
    return -1;
  }
  
  TH2D* h2 = (TH2D*)f->Get("fh2FFJetPtRec");
  if (!h2) {
    printf("Failed to load histo from file \"%s\"!\n", filePathName.Data());
    f->Close();
    return -1;
  }
  
  TH1D* h = h2->ProjectionY("hRecJetSpectrum");
  setupHisto(h);
  h->SetDirectory(0x0);
  
  f->Close();
  
  TCanvas* canv = createCanvas("canv_recJetSpectrum", "Reconstructed Jet Spectrum");
  canv->cd();
  h->Draw("E1BAR");
  
  
  ClearTitleFromHistoInCanvas(canv);
  
  TString savePathName = Form("%s/RecJetSpectrum_10d_e_pass2.root", path.Data());
  canv->SaveAs(savePathName.Data());
  
  savePathName = savePathName.ReplaceAll(".root", ".pdf");
  canv->SaveAs(savePathName.Data());
    
  return 0;
}
