// ------------------------------------------------------
//
// Class to handle 2d-corrections. 
//
// ------------------------------------------------------
//

/* $Id$ */

#include <TFile.h>
#include <TCanvas.h>

#include <AliLog.h>

#include "CorrectionMatrix2D.h"

//____________________________________________________________________
ClassImp(CorrectionMatrix2D)

//____________________________________________________________________
CorrectionMatrix2D::CorrectionMatrix2D(const CorrectionMatrix2D& c) 
  : TNamed(c)
{
  // copy constructor
  ((CorrectionMatrix2D &)c).Copy(*this);
}

//____________________________________________________________________
CorrectionMatrix2D::CorrectionMatrix2D(Char_t* name, Char_t* title,
				       Int_t nBinX, Float_t Xmin, Float_t Xmax,
				       Int_t nBinY, Float_t Ymin, Float_t Ymax) 
  : TNamed(name, title)
{
  //
  // constructor
  //


  fhMeas  = new TH2F(Form("meas_%s",name), Form("meas_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhGene  = new TH2F(Form("gene_%s",name), Form("gene_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhCorr  = new TH2F(Form("corr_%s",name), Form("corr_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
}

//____________________________________________________________________
CorrectionMatrix2D::CorrectionMatrix2D(Char_t* name,Char_t* title, 
				       Int_t nBinX, Float_t *X, Int_t nBinY, Float_t *Y) 
  : TNamed(name, title) 
{
  //
  // constructor
  //

  fhMeas  = new TH2F(Form("meas_%s",name), Form("meas_%s",title),  nBinX, X, nBinY, Y);
  fhGene  = new TH2F(Form("gene_%s",name), Form("gene_%s",title),  nBinX, X, nBinY, Y);
  fhCorr  = new TH2F(Form("corr_%s",name), Form("corr_%s",title),  nBinX, X, nBinY, Y);
}


//____________________________________________________________________
CorrectionMatrix2D::~CorrectionMatrix2D() {
  //
  // destructor
  //
  if (fhMeas)  delete fhMeas;
  if (fhGene)  delete fhGene;
  if (fhCorr)  delete fhCorr;
}

//____________________________________________________________________
CorrectionMatrix2D &CorrectionMatrix2D::operator=(const CorrectionMatrix2D &c)
{
  // assigment operator

  if (this != &c) 
    ((CorrectionMatrix2D &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
TH1F* CorrectionMatrix2D::Get1DCorrection(Char_t* opt) {
  //
  // integrate the correction over one variable 
  // 

  TH1D* meas1D; 
  TH1D* gene1D; 

  if (strcmp(opt,"x")==0) {
    meas1D = fhMeas->ProjectionX();
    gene1D = fhGene->ProjectionX();      
  }
  if (strcmp(opt,"y")==0) {
    meas1D = fhMeas->ProjectionY();
    gene1D = fhGene->ProjectionY();      
  }
  gene1D->Sumw2();

  gene1D->SetName(Form("corr_1D_%s",fName.Data()));
  gene1D->SetTitle(Form("corr_1D_%s",fName.Data()));

  gene1D->Divide(gene1D, meas1D, 1, 1, "B");
  
  return (TH1F*)gene1D;   
}


//____________________________________________________________________
void
CorrectionMatrix2D::Copy(TObject& c) const 
{
  // copy function

  CorrectionMatrix2D& target = (CorrectionMatrix2D &) c;

  target.fhMeas  = fhMeas;
  target.fhGene  = fhGene;
  target.fhCorr  = fhCorr;
}


//________________________________________________________________________
void CorrectionMatrix2D::SetAxisTitles(Char_t* titleX, Char_t* titleY) 
{ 
  //
  // method for setting the axis titles of the histograms
  //

  fhMeas ->SetXTitle(titleX);  fhMeas ->SetYTitle(titleY);
  fhGene ->SetXTitle(titleX);  fhGene ->SetYTitle(titleY);
  fhCorr ->SetXTitle(titleX);  fhCorr ->SetYTitle(titleY);
}

//____________________________________________________________________
Long64_t CorrectionMatrix2D::Merge(TCollection* list) {
  // Merge a list of CorrectionMatrix2D objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of measured and generated histograms
  TList* collectionMeas = new TList;
  TList* collectionGene = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    CorrectionMatrix2D* entry = dynamic_cast<CorrectionMatrix2D*> (obj);
    if (entry == 0) 
      continue;

    collectionMeas->Add(entry->GetMeasuredHistogram());
    collectionGene->Add(entry->GetGeneratedHistogram());

    count++;
  }
  fhMeas->Merge(collectionMeas);
  fhGene->Merge(collectionGene);

  // is this really faster than just adding the histograms in the list???
  delete collectionMeas;
  delete collectionGene;


  return count+1;
}


//____________________________________________________________________
void CorrectionMatrix2D::Divide() {  
  //
  // divide the histograms to get the correction
  // 

  if (!fhMeas || !fhGene)  return; 

  fhCorr->Divide(fhGene, fhMeas, 1,1,"B");
  
}

//____________________________________________________________________
void
CorrectionMatrix2D::RemoveEdges(Float_t cut, Int_t nBinsXedge, Int_t nBinsYedge) 
{
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
  //
  // loads the histograms from a file
  //
  
  TFile* fin = TFile::Open(fileName);  
  
  if(!fin) {
    //Info("LoadHistograms",Form(" %s file does not exist",fileName));
    return kFALSE;
  }
  
  if(fhGene)  {delete fhGene;  fhGene=0;}
  if(fhCorr)  {delete fhCorr;  fhCorr=0;}
  if(fhMeas)  {delete fhMeas;  fhMeas=0;}
  
  fhMeas  = (TH2F*)fin->Get(Form("%s/meas_%s",dir,fName.Data()));
      if(!fhMeas)  Info("LoadHistograms","No meas  hist available");
  fhGene  = (TH2F*)fin->Get(Form("%s/gene_%s",dir,fName.Data()));
      if(!fhGene)  Info("LoadHistograms","No gene  hist available");
  fhCorr  = (TH2F*)fin->Get(Form("%s/corr_%s",dir,fName.Data()));
      if(!fhCorr) {Info("LoadHistograms","No corr  hist available");
      return kFALSE;}
      
  return kTRUE;
}


//____________________________________________________________________
void
CorrectionMatrix2D::SaveHistograms() {
  //
  // saves the histograms 
  //
  
  gDirectory->mkdir(fName.Data());
  gDirectory->cd(fName.Data());
  
  fhMeas ->Write();
  fhGene ->Write();

  if (fhCorr)
    fhCorr->Write();

  gDirectory->cd("../");
}

//____________________________________________________________________
void CorrectionMatrix2D::DrawHistograms()
{
  //
  // draws all the four histograms on one TCanvas
  //

  TCanvas* canvas = new TCanvas(Form("correction_%s",fName.Data()), 
				Form("correction_%s",fName.Data()), 800, 800);
  canvas->Divide(2, 2);
  
  canvas->cd(1);
  if (fhMeas)
    fhMeas->Draw("COLZ");
  
  canvas->cd(2);
  if (fhGene)
    fhGene->Draw("COLZ");

  canvas->cd(3);
  if (fhCorr)
    fhCorr->Draw("COLZ");

  canvas->cd(4);

  // add: draw here the stat. errors of the correction histogram
  
}


