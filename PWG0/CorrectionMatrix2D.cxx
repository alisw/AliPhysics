#include "CorrectionMatrix2D.h"

//____________________________________________________________________
ClassImp(CorrectionMatrix2D);

//____________________________________________________________________
CorrectionMatrix2D::CorrectionMatrix2D(Char_t* name, Char_t* title,
				       Int_t nBinX, Float_t Xmin, Float_t Xmax,
				       Int_t nBinY, Float_t Ymin, Float_t Ymax) 
  : TNamed(name, title)
{

  fhMeas  = new TH2F(Form("meas_%s",name), Form("meas_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhGene  = new TH2F(Form("gene_%s",name), Form("gene_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhCorr  = new TH2F(Form("corr_%s",name), Form("corr_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhRatio = new TH2F(Form("ratio_%s",name),Form("ratio_%s",title), nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
}

//____________________________________________________________________
CorrectionMatrix2D::CorrectionMatrix2D(Char_t* name,Char_t* title, 
				       Int_t nBinX, Float_t *X, Int_t nBinY, Float_t *Y) 
  : TNamed(name, title)
{
  fhMeas  = new TH2F(Form("meas_%s",name), Form("meas_%s",title),  nBinX, X, nBinY, Y);
  fhGene  = new TH2F(Form("gene_%s",name), Form("gene_%s",title),  nBinX, X, nBinY, Y);
  fhCorr  = new TH2F(Form("corr_%s",name), Form("corr_%s",title),  nBinX, X, nBinY, Y);
  fhRatio = new TH2F(Form("ratio_%s",name),Form("ratio_%s",title), nBinX, X, nBinY, Y);
}


//____________________________________________________________________
CorrectionMatrix2D::~CorrectionMatrix2D() {
  // Destructor
  //
  if (fhMeas)  delete fhMeas;
  if (fhGene)  delete fhCorr;
  if (fhRatio) delete fhRatio;
  if (fhCorr)  delete fhCorr;
}


//________________________________________________________________________
void CorrectionMatrix2D::SetAxisTitles(Char_t* titleX, Char_t* titleY) 
{ 
  fhMeas ->SetXTitle(titleX);  fhMeas ->SetYTitle(titleY);
  fhGene ->SetXTitle(titleX);  fhGene ->SetYTitle(titleY);
  fhCorr ->SetXTitle(titleX);  fhCorr ->SetYTitle(titleY);
  fhRatio->SetXTitle(titleX);  fhRatio->SetYTitle(titleY);
}


//____________________________________________________________________
void CorrectionMatrix2D::Finish() {  

if (!fhMeas || !fhGene)  return;
 

  fhRatio->Divide(fhMeas, fhGene, 1,1,"B");
  fhCorr->Divide(fhGene, fhMeas, 1,1,"B");

  Int_t nBinsVtx = fhCorr->GetNbinsX();
  Int_t nBinsEta = fhCorr->GetNbinsY();

  TH2F* tmp = (TH2F*)fhCorr->Clone("tmp");

  // cut at 0.2
  for (Int_t bx=0; bx<=nBinsVtx; bx++) {
    for (Int_t by=0; by<=nBinsEta; by++) {
      if (tmp->GetBinContent(bx,by)<0.2) {
	fhCorr->SetBinContent(bx,by,0);
	fhCorr->SetBinError(bx,by,0);
	
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
  
  Int_t nBinsX = fhCorr->GetNbinsX();
  Int_t nBinsY = fhCorr->GetNbinsY();

  // set bin content to zero for bins with content smaller cut
  for (Int_t bx=0; bx<=nBinsX; bx++) {
    for (Int_t by=0; by<=nBinsY; by++) {
      if (fhCorr->GetBinContent(bx,by)>cut) {
	  fhCorr->SetBinContent(bx,by,0);
	  fhCorr->SetBinError(bx,by,0);
      }
    }
  }

  // set bin content to zero for bins next to bins with zero
  TH2F* tmp = (TH2F*)fhCorr->Clone("tmp");
  tmp->Reset();
  
  Bool_t done = kFALSE;
  Int_t nBinsXCount = 0;
  Int_t nBinsYCount = 0;
  while (!done) {    
    if (nBinsXCount<nBinsXedge) 
      for (Int_t bx=0; bx<=nBinsX; bx++) {
	for (Int_t by=0; by<=nBinsY; by++) {
	  if ((fhCorr->GetBinContent(bx+1,by)==0)|| 
	      (fhCorr->GetBinContent(bx-1,by)==0))
	    tmp->SetBinContent(bx,by,1);	
	  
	}
      }
    if (nBinsYCount<nBinsYedge) 
      for (Int_t bx=0; bx<=nBinsX; bx++) {
	for (Int_t by=0; by<=nBinsY; by++) {
	  if ((fhCorr->GetBinContent(bx,by+1)==0)|| 
	      (fhCorr->GetBinContent(bx,by-1)==0))
	    tmp->SetBinContent(bx,by,1);	
	}
      }    
    for (Int_t bx=0; bx<=nBinsX; bx++) {
      for (Int_t by=0; by<=nBinsY; by++) {
	if (tmp->GetBinContent(bx,by)==1) {
	  fhCorr->SetBinContent(bx,by,0);
	  fhCorr->SetBinError(bx,by,0);
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
  
  if(fhGene)  {delete fhGene;  fhGene=0;}
  if(fhCorr)  {delete fhCorr;  fhCorr=0;}
  if(fhRatio) {delete fhRatio; fhRatio=0;}
  if(fhMeas)  {delete fhMeas;  fhMeas=0;}
  
  fhMeas  = (TH2F*)fin->Get(Form("%s/meas",dir));
      if(!fhMeas)  Info("LoadHistograms","No meas  hist available");
  fhGene  = (TH2F*)fin->Get(Form("%s/gene",dir));
      if(!fhGene)  Info("LoadHistograms","No gene  hist available");
  fhRatio = (TH2F*)fin->Get(Form("%s/ratio",dir));
      if(!fhRatio) Info("LoadHistograms","No ratio hist available");
  fhCorr  = (TH2F*)fin->Get(Form("%s/corr",dir));
      if(!fhCorr) {Info("LoadHistograms","No corr  hist available");
      return kFALSE;}
      
  return kTRUE;
}


//____________________________________________________________________
void
CorrectionMatrix2D::SaveHistograms() {

  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());
  
  fhMeas ->Write();
  fhGene ->Write();

  if (fhCorr)
    fhCorr->Write();

  if (fhRatio)
    fhRatio->Write();

  gDirectory->cd("../");
}

//____________________________________________________________________
void CorrectionMatrix2D::DrawHistograms()
{
  TCanvas* canvas = new TCanvas("Correction", "Correction", 800, 800);
  canvas->Divide(2, 2);

  canvas->cd(1);
  if (fhMeas)
    fhMeas->Draw("COLZ");

  canvas->cd(2);
  if (fhGene)
    fhGene->Draw("COLZ");

  canvas->cd(3);
  if (fhRatio)
    fhRatio->Draw("COLZ");

  canvas->cd(4);
  if (fhCorr)
    fhCorr->Draw("COLZ");
}


