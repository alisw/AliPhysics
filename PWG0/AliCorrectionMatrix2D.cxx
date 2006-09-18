/* $Id$ */

// ------------------------------------------------------
//
// Class to handle 2d-corrections.
//
// ------------------------------------------------------
//

#include <TH2F.h>

#include <AliLog.h>

#include "AliCorrectionMatrix2D.h"

//____________________________________________________________________
ClassImp(AliCorrectionMatrix2D)

//____________________________________________________________________
AliCorrectionMatrix2D::AliCorrectionMatrix2D() :
  AliCorrectionMatrix()
{
  // default constructor
}

//____________________________________________________________________
AliCorrectionMatrix2D::AliCorrectionMatrix2D(const AliCorrectionMatrix2D& c)
  : AliCorrectionMatrix(c)
{
  // copy constructor
  ((AliCorrectionMatrix2D &)c).Copy(*this);
}

//____________________________________________________________________
AliCorrectionMatrix2D &AliCorrectionMatrix2D::operator=(const AliCorrectionMatrix2D &c)
{
  // assigment operator

  if (this != &c)
    ((AliCorrectionMatrix2D &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
AliCorrectionMatrix2D::AliCorrectionMatrix2D(const Char_t* name, const Char_t* title,
				       Int_t nBinX, Float_t Xmin, Float_t Xmax,
				       Int_t nBinY, Float_t Ymin, Float_t Ymax) 
  : AliCorrectionMatrix(name, title)
{
  //
  // constructor
  //

  fhMeas  = new TH2F(Form("meas_%s",name), Form("meas_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhGene  = new TH2F(Form("gene_%s",name), Form("gene_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhCorr  = new TH2F(Form("corr_%s",name), Form("corr_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);

  fhMeas->Sumw2();
  fhGene->Sumw2();
  fhCorr->Sumw2();
}

//____________________________________________________________________
AliCorrectionMatrix2D::AliCorrectionMatrix2D(const Char_t* name, const Char_t* title,
				       Int_t nBinX, Float_t *X, Int_t nBinY, Float_t *Y) 
  : AliCorrectionMatrix(name, title)
{
  //
  // constructor
  //

  fhMeas  = new TH2F(Form("meas_%s",name), Form("meas_%s",title),  nBinX, X, nBinY, Y);
  fhGene  = new TH2F(Form("gene_%s",name), Form("gene_%s",title),  nBinX, X, nBinY, Y);
  fhCorr  = new TH2F(Form("corr_%s",name), Form("corr_%s",title),  nBinX, X, nBinY, Y);

  fhMeas->Sumw2();
  fhGene->Sumw2();
  fhCorr->Sumw2();
}

//____________________________________________________________________
AliCorrectionMatrix2D::~AliCorrectionMatrix2D()
{
  //
  // destructor
  //

  // histograms already deleted in base class
}

TH2F* AliCorrectionMatrix2D::GetGeneratedHistogram() const
{
  // return generated histogram casted to correct type
  return dynamic_cast<TH2F*> (fhGene);
}

TH2F* AliCorrectionMatrix2D::GetMeasuredHistogram() const
{
  // return measured histogram casted to correct type
  return dynamic_cast<TH2F*> (fhMeas);
}

//____________________________________________________________________
TH1F* AliCorrectionMatrix2D::Get1DCorrection(Char_t* opt, Float_t min, Float_t max)
{
  //
  // integrate the correction over one variable 
  // 

  TH1D* meas1D = 0;
  TH1D* gene1D = 0; 

  if (strcmp(opt,"x")==0) {
    Int_t binMin = GetMeasuredHistogram()->GetYaxis()->FindBin(min);
    Int_t binMax = GetMeasuredHistogram()->GetYaxis()->FindBin(max);

    if (min==0 && max==0) {
      meas1D = GetMeasuredHistogram()->ProjectionX();
      gene1D = GetGeneratedHistogram()->ProjectionX();
    }
    else {
      AliDebug(AliLog::kDebug+1, Form("Getting 1D map. Including y-bins %d to %d \n", binMin, binMax));

      meas1D = GetMeasuredHistogram()->ProjectionX("pm",binMin,binMax);
      gene1D = GetGeneratedHistogram()->ProjectionX("pg",binMin,binMax);
    }
  }
  if (strcmp(opt,"y")==0) {
    Int_t binMin = GetMeasuredHistogram()->GetXaxis()->FindBin(min);
    Int_t binMax = GetMeasuredHistogram()->GetXaxis()->FindBin(max);

    if (min==0 && max==0) {
      meas1D = GetMeasuredHistogram()->ProjectionY();
      gene1D = GetGeneratedHistogram()->ProjectionY();
    }
    else {
      AliDebug(AliLog::kDebug+1, Form("Getting 1D map. Including x-bins %d to %d \n", binMin, binMax));

      meas1D = GetMeasuredHistogram()->ProjectionY("pm", binMin, binMax);
      gene1D = GetGeneratedHistogram()->ProjectionY("pg", binMin, binMax);
    }
  }
  gene1D->Sumw2();

  gene1D->SetName(Form("corr_1D_%s",fName.Data()));
  gene1D->SetTitle(Form("corr_1D_%s",fName.Data()));

  gene1D->Divide(gene1D, meas1D, 1, 1, "B");
  
  return (TH1F*)gene1D;   
}

//____________________________________________________________________
void AliCorrectionMatrix2D::FillMeas(Float_t ax, Float_t ay)
{
  // add value to measured histogram
  GetMeasuredHistogram()->Fill(ax, ay);
}

//____________________________________________________________________
void AliCorrectionMatrix2D::FillGene(Float_t ax, Float_t ay)
{
  // add value to generated histogram
  GetGeneratedHistogram()->Fill(ax, ay);
}

//____________________________________________________________________
Float_t AliCorrectionMatrix2D::GetCorrection(Float_t ax, Float_t ay) const
{
  // returns a value of the correction map
  return fhCorr->GetBinContent(fhCorr->FindBin(ax,ay));
}

//____________________________________________________________________
void AliCorrectionMatrix2D::RemoveEdges(Float_t cut, Int_t nBinsXedge, Int_t nBinsYedge)
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

