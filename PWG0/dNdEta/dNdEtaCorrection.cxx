/* $Id$ */

#include "dNdEtaCorrection.h"

#include <TCanvas.h>

//____________________________________________________________________
ClassImp(dNdEtaCorrection)

//____________________________________________________________________
dNdEtaCorrection::dNdEtaCorrection(Char_t* name) {

  fName = TString(name);
  
  hEtaVsVtx_meas  = new TH2F("etaVsVtx_meas", "etaVsVtx_meas" ,80,-20,20,120,-6,6);
  hEtaVsVtx_gene  = new TH2F("etaVsVtx_gene", "etaVsVtx_gene" ,80,-20,20,120,-6,6);
  hEtaVsVtx_corr  = new TH2F("etaVsVtx_corr", "etaVsVtx_corr", 80,-20,20,120,-6,6);

  hEtaVsVtx_ratio = new TH2F("etaVsVtx_ratio","etaVsVtx_ratio",80,-20,20,120,-6,6);

  hEtaVsVtx_meas ->SetXTitle("vtx z [cm]");  hEtaVsVtx_meas ->SetYTitle("#eta"); 
  hEtaVsVtx_gene ->SetXTitle("vtx z [cm]");  hEtaVsVtx_gene ->SetYTitle("#eta"); 
  hEtaVsVtx_corr ->SetXTitle("vtx z [cm]");  hEtaVsVtx_corr ->SetYTitle("#eta"); 
  hEtaVsVtx_ratio->SetXTitle("vtx z [cm]");  hEtaVsVtx_ratio->SetYTitle("#eta"); 

}

//____________________________________________________________________
void
dNdEtaCorrection::Finish() {  

  hEtaVsVtx_ratio->Divide(hEtaVsVtx_meas, hEtaVsVtx_gene, 1,1,"B");
  hEtaVsVtx_corr->Divide(hEtaVsVtx_gene, hEtaVsVtx_meas, 1,1,"B");

  Int_t nBinsVtx = hEtaVsVtx_corr->GetNbinsX();
  Int_t nBinsEta = hEtaVsVtx_corr->GetNbinsY();

  TH2F* tmp = (TH2F*)hEtaVsVtx_corr->Clone("tmp");

  // cut at 0.2
  for (Int_t bx=0; bx<=nBinsVtx; bx++) {
    for (Int_t by=0; by<=nBinsEta; by++) {
      if (tmp->GetBinContent(bx,by)<0.2) {
	hEtaVsVtx_corr->SetBinContent(bx,by,0);
	hEtaVsVtx_corr->SetBinError(bx,by,0);
	
	tmp->SetBinContent(bx,by,0);
      }
      else 
	tmp->SetBinContent(bx,by,1);
    }
  }
  
}

//____________________________________________________________________
void
dNdEtaCorrection::RemoveEdges(Float_t cut, Int_t nBinsVtx, Int_t nBinsEta) {

  // remove edges of correction histogram by removing
  // - bins with content bigger than cut
  // - bins next to bins with zero bin content

  Int_t nBinsX = hEtaVsVtx_corr->GetNbinsX();
  Int_t nBinsY = hEtaVsVtx_corr->GetNbinsY();

  // set bin content to zero for bins with content bigger than cut
  for (Int_t bx=0; bx<=nBinsX; bx++) {
    for (Int_t by=0; by<=nBinsY; by++) {
      if (hEtaVsVtx_corr->GetBinContent(bx,by)>cut) {
	hEtaVsVtx_corr->SetBinContent(bx,by,0);
	hEtaVsVtx_corr->SetBinError(bx,by,0);
      }
    }
  }

  // set bin content to zero for bins next to bins with zero
  TH2F* tmp = (TH2F*)hEtaVsVtx_corr->Clone("tmp");
  tmp->Reset();

  Bool_t done = kFALSE;
  Int_t nBinsVtxCount = 0;
  Int_t nBinsEtaCount = 0;
  while (!done) {
    if (nBinsVtxCount<nBinsVtx)
      for (Int_t bx=0; bx<=nBinsX; bx++) {
	for (Int_t by=0; by<=nBinsY; by++) {
	  if ((hEtaVsVtx_corr->GetBinContent(bx+1,by)==0)||
	      (hEtaVsVtx_corr->GetBinContent(bx-1,by)==0))
	    tmp->SetBinContent(bx,by,1);

	}
      }
    if (nBinsEtaCount<nBinsEta)
      for (Int_t bx=0; bx<=nBinsX; bx++) {
	for (Int_t by=0; by<=nBinsY; by++) {
	  if ((hEtaVsVtx_corr->GetBinContent(bx,by+1)==0)||
	      (hEtaVsVtx_corr->GetBinContent(bx,by-1)==0))
	    tmp->SetBinContent(bx,by,1);
	}
      }
    for (Int_t bx=0; bx<=nBinsX; bx++) {
      for (Int_t by=0; by<=nBinsY; by++) {
	if (tmp->GetBinContent(bx,by)==1) {
	  hEtaVsVtx_corr->SetBinContent(bx,by,0);
	  hEtaVsVtx_corr->SetBinError(bx,by,0);
	}
      }
    }
    nBinsVtxCount++;
    nBinsEtaCount++;
    if ((nBinsVtxCount>=nBinsVtx)&&(nBinsEtaCount>=nBinsEta)) done=kTRUE;
  }
  tmp->Delete();

}

//____________________________________________________________________
Bool_t
dNdEtaCorrection::LoadHistograms(Char_t* fileName, Char_t* dir) {

  TFile* fin = TFile::Open(fileName);  
  
  // add test of file
  // return kFALSE

  hEtaVsVtx_meas  = (TH2F*)fin->Get(Form("%s/etaVsVtx_meas",dir));
  hEtaVsVtx_gene  = (TH2F*)fin->Get(Form("%s/etaVsVtx_gene",dir));
  hEtaVsVtx_corr  = (TH2F*)fin->Get(Form("%s/etaVsVtx_corr",dir));

  hEtaVsVtx_ratio = (TH2F*)fin->Get(Form("%s/etaVsVtx_ratio",dir));

  return kTRUE;
}


//____________________________________________________________________
void
dNdEtaCorrection::SaveHistograms() {

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());
  
  hEtaVsVtx_meas ->Write();
  hEtaVsVtx_gene ->Write();

  if (hEtaVsVtx_corr)
    hEtaVsVtx_corr->Write();

  if (hEtaVsVtx_ratio)
    hEtaVsVtx_ratio->Write();


  gDirectory->cd("../");
}

//____________________________________________________________________
void dNdEtaCorrection::DrawHistograms()
{
  TCanvas* canvas = new TCanvas("dNdEtaCorrection", "dNdEtaCorrection", 800, 800);
  canvas->Divide(2, 2);
  
  canvas->cd(1);
  if (hEtaVsVtx_meas)
    hEtaVsVtx_meas->Draw("COLZ");

  canvas->cd(2);
  if (hEtaVsVtx_gene)
    hEtaVsVtx_gene->Draw("COLZ");

  canvas->cd(3);
  if (hEtaVsVtx_ratio)
    hEtaVsVtx_ratio->Draw("COLZ");

  canvas->cd(4);
  if (hEtaVsVtx_corr)
    hEtaVsVtx_corr->Draw("COLZ");
}
