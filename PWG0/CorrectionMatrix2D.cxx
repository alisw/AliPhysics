#include "CorrectionMatrix2D.h"

//____________________________________________________________________
ClassImp(CorrectionMatrix2D);

//____________________________________________________________________
CorrectionMatrix2D::CorrectionMatrix2D(Char_t* name) {

  fName = TString(name);

}
//________________________________________________________________________
void CorrectionMatrix2D::SetHist(Char_t* title ,Int_t nBinX, Float_t Xmin, Float_t Xmax,
				  Int_t nBinY, Float_t Ymin, Float_t Ymax) 
{
//Char_t* name , Chat_t* title ,   

  hmeas  = new TH2F("meas",  Form("meas_%s",title), nBinX,Xmin,Xmax,nBinY,Ymin,Ymax);
  hgene  = new TH2F("gene",  Form("gene_%s",title), nBinX,Xmin,Xmax,nBinY,Ymin,Ymax);
  hcorr  = new TH2F("corr",  Form("corr_%s",title), nBinX,Xmin,Xmax,nBinY,Ymin,Ymax);
  hratio = new TH2F("ratio", Form("ratio_%s",title), nBinX,Xmin,Xmax,nBinY,Ymin,Ymax);
}


void CorrectionMatrix2D::SetHist(Char_t* title ,Int_t nBinX, Float_t *X,
				                Int_t nBinY, Float_t *Y) 
{

  hmeas  = new TH2F("meas",  Form("meas_%s",title), nBinX,X,nBinY,Y);
  hgene  = new TH2F("gene",  Form("gene_%s",title), nBinX,X,nBinY,Y);
  hcorr  = new TH2F("corr",  Form("corr_%s",title), nBinX,X,nBinY,Y);
  hratio = new TH2F("ratio", Form("ratio_%s",title),nBinX,X,nBinY,Y);
  
}



//________________________________________________________________________
void CorrectionMatrix2D::SetHistTitle(Char_t* titleX, Char_t* titleY) 
{
 
  hmeas ->SetXTitle(titleX);  hmeas ->SetYTitle(titleY); 
  hgene ->SetXTitle(titleX);  hgene ->SetYTitle(titleY); 
  hcorr ->SetXTitle(titleX);  hcorr ->SetYTitle(titleY); 
  hratio->SetXTitle(titleX);  hratio->SetYTitle(titleY); 

}


//____________________________________________________________________
void CorrectionMatrix2D::Finish() {  

if (!hmeas || !hgene)  return;
 
TCanvas *c1 = new TCanvas();
c1->Divide(2,1);
c1->cd(1);
hgene->Draw();
c1->cd(2);
hmeas->Draw();

  hratio->Divide(hmeas, hgene, 1,1,"B");
  hcorr->Divide(hgene, hmeas, 1,1,"B");

  Int_t nBinsVtx = hcorr->GetNbinsX();
  Int_t nBinsEta = hcorr->GetNbinsY();

  TH2F* tmp = (TH2F*)hcorr->Clone("tmp");

  // cut at 0.2
  for (Int_t bx=0; bx<=nBinsVtx; bx++) {
    for (Int_t by=0; by<=nBinsEta; by++) {
      if (tmp->GetBinContent(bx,by)<0.2) {
	hcorr->SetBinContent(bx,by,0);
	hcorr->SetBinError(bx,by,0);
	
	tmp->SetBinContent(bx,by,0);
      }
      else 
	tmp->SetBinContent(bx,by,1);
    }
  }
  
}

//____________________________________________________________________
void
CorrectionMatrix2D::RemoveEdges(Float_t cut, Int_t nBinsXedge, Int_t nBinsYedge) {

  // remove edges of correction histogram by removing 
  // - bins with content less than cut
  // - bins next to bins with zero bin content
  
  Int_t nBinsX = hcorr->GetNbinsX();
  Int_t nBinsY = hcorr->GetNbinsY();

  // set bin content to zero for bins with content smaller cut
  for (Int_t bx=0; bx<=nBinsX; bx++) {
    for (Int_t by=0; by<=nBinsY; by++) {
      if (hcorr->GetBinContent(bx,by)>cut) {
	  hcorr->SetBinContent(bx,by,0);
	  hcorr->SetBinError(bx,by,0);
      }
    }
  }

  // set bin content to zero for bins next to bins with zero
  TH2F* tmp = (TH2F*)hcorr->Clone("tmp");
  tmp->Reset();
  
  Bool_t done = kFALSE;
  Int_t nBinsXCount = 0;
  Int_t nBinsYCount = 0;
  while (!done) {    
    if (nBinsXCount<nBinsXedge) 
      for (Int_t bx=0; bx<=nBinsX; bx++) {
	for (Int_t by=0; by<=nBinsY; by++) {
	  if ((hcorr->GetBinContent(bx+1,by)==0)|| 
	      (hcorr->GetBinContent(bx-1,by)==0))
	    tmp->SetBinContent(bx,by,1);	
	  
	}
      }
    if (nBinsYCount<nBinsYedge) 
      for (Int_t bx=0; bx<=nBinsX; bx++) {
	for (Int_t by=0; by<=nBinsY; by++) {
	  if ((hcorr->GetBinContent(bx,by+1)==0)|| 
	      (hcorr->GetBinContent(bx,by-1)==0))
	    tmp->SetBinContent(bx,by,1);	
	}
      }    
    for (Int_t bx=0; bx<=nBinsX; bx++) {
      for (Int_t by=0; by<=nBinsY; by++) {
	if (tmp->GetBinContent(bx,by)==1) {
	  hcorr->SetBinContent(bx,by,0);
	  hcorr->SetBinError(bx,by,0);
	}
      }
    }
    nBinsXCount++;
    nBinsYCount++;
    if ((nBinsXCount>=nBinsXedge)&&(nBinsYCount>=nBinsYedge)) done=kTRUE;
  }
  tmp->Delete();  

}

//____________________________________________________________________
Bool_t CorrectionMatrix2D::LoadHistograms(Char_t* fileName, Char_t* dir) {

  TFile* fin = TFile::Open(fileName);  
  
  if(!fin) 
  {
    //Info("LoadHistograms",Form(" %s file does not exist",fileName));
  return kFALSE;
  }
  
  if(hgene)  {delete hgene;  hgene=0;}
  if(hcorr)  {delete hcorr;  hcorr=0;}
  if(hratio) {delete hratio; hratio=0;}
  if(hmeas)  {delete hmeas;  hmeas=0;}
  
  hmeas  = (TH2F*)fin->Get(Form("%s/meas",dir));
      if(!hmeas)  Info("LoadHistograms","No meas  hist available");
  hgene  = (TH2F*)fin->Get(Form("%s/gene",dir));
      if(!hgene)  Info("LoadHistograms","No gene  hist available");
  hratio = (TH2F*)fin->Get(Form("%s/ratio",dir));
      if(!hratio) Info("LoadHistograms","No ratio hist available");
  hcorr  = (TH2F*)fin->Get(Form("%s/corr",dir));
      if(!hcorr) {Info("LoadHistograms","No corr  hist available");
      return kFALSE;}
      
  return kTRUE;
}


//____________________________________________________________________
void
CorrectionMatrix2D::SaveHistograms() {

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());
  
  hmeas ->Write();
  hgene ->Write();

  if (hcorr)
    hcorr->Write();

  if (hratio)
    hratio->Write();

  gDirectory->cd("../");
}

//____________________________________________________________________
void CorrectionMatrix2D::DrawHistograms()
{
  TCanvas* canvas = new TCanvas("Correction", "Correction", 800, 800);
  canvas->Divide(2, 2);

  canvas->cd(1);
  if (hmeas)
    hmeas->Draw("COLZ");

  canvas->cd(2);
  if (hgene)
    hgene->Draw("COLZ");

  canvas->cd(3);
  if (hratio)
    hratio->Draw("COLZ");

  canvas->cd(4);
  if (hcorr)
    hcorr->Draw("COLZ");
}


