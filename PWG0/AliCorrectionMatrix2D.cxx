/* $Id$ */

// ------------------------------------------------------
//
// Class to handle 2d-corrections.
//
// ------------------------------------------------------
//

#include <TH2F.h>
#include <TMath.h>

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

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhMeas  = new TH2F("measured",   Form("%s measured", GetTitle()),   nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhGene  = new TH2F("generated",  Form("%s generated", GetTitle()),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);
  fhCorr  = new TH2F("correction", Form("%s correction", GetTitle()), nBinX, Xmin, Xmax, nBinY, Ymin, Ymax);

  TH1::AddDirectory(oldStatus);

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

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

	fhMeas  = new TH2F("measured",   Form("%s measured",title),   nBinX, X, nBinY, Y);
  fhGene  = new TH2F("generated",  Form("%s generated",title),  nBinX, X, nBinY, Y);
  fhCorr  = new TH2F("correction", Form("%s correction",title), nBinX, X, nBinY, Y);

  TH1::AddDirectory(oldStatus);

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

TH2* AliCorrectionMatrix2D::GetGeneratedHistogram() const
{
  // return generated histogram casted to correct type
  return dynamic_cast<TH2*> (fhGene);
}

TH2* AliCorrectionMatrix2D::GetMeasuredHistogram() const
{
  // return measured histogram casted to correct type
  return dynamic_cast<TH2*> (fhMeas);
}

//____________________________________________________________________
TH1* AliCorrectionMatrix2D::Get1DCorrectionHistogram(const Char_t* opt, Float_t min, Float_t max, Bool_t binomialErrors)
{
  //
  // integrate the correction over one variable 
  // 

  TH1D* meas1D = 0;
  TH1D* gene1D = 0; 

  if (strcmp(opt,"x")==0) {
    Int_t binMin = fhMeas->GetYaxis()->FindBin(min);
    Int_t binMax = fhGene->GetYaxis()->FindBin(max);

    if (min>=max) {
      meas1D = ((TH2F*)fhMeas)->ProjectionX();
      gene1D = ((TH2F*)fhGene)->ProjectionX();
    }
    else {
      Printf("Getting 1D map. Including y-bins %d to %d", binMin, binMax);

      meas1D = ((TH2F*)fhMeas)->ProjectionX(Form("%s_x_pm", GetName()),binMin,binMax);
      gene1D = ((TH2F*)fhGene)->ProjectionX(Form("%s_x_pg", GetName()),binMin,binMax);
    }
  }
  else if (strcmp(opt,"y")==0) {
    Int_t binMin = fhMeas->GetXaxis()->FindBin(min);
    Int_t binMax = fhMeas->GetXaxis()->FindBin(max);

    if (min>=max) {
      meas1D = ((TH2F*)fhMeas)->ProjectionY();
      gene1D = ((TH2F*)fhGene)->ProjectionY();
    }
    else {
      Printf("Getting 1D map. Including x-bins %d to %d \n", binMin, binMax);

      meas1D = ((TH2F*)fhMeas)->ProjectionY(Form("%s_y_pm", GetName()), binMin, binMax);
      gene1D = ((TH2F*)fhGene)->ProjectionY(Form("%s_y_pg", GetName()), binMin, binMax);
    }
  }
  else {
    Printf("ERROR: Invalid option");
    return 0;
  }

  if (!binomialErrors)
  {
    // set the errors on gene manually, and clear the ones on meas.
    gene1D->Sumw2();
    for (Int_t bin=0; bin <= gene1D->GetNbinsX()+1; bin++)
    {
      gene1D->SetBinError(bin, TMath::Sqrt(gene1D->GetBinContent(bin)));
      meas1D->SetBinError(bin, 0);
    }
  }
  
  gene1D->SetName(Form("corr_1D_%s",fName.Data()));
  gene1D->SetTitle(Form("corr_1D_%s",fName.Data()));
 
  TH1* divided = (TH1*) gene1D->Clone(Form("corr_1D_%s",fName.Data()));
  divided->Reset();
  
  divided->Divide(gene1D, meas1D, 1, 1, (binomialErrors) ? "B" : "");

  Printf("%p %p", gene1D, meas1D);
  
  return (TH1F*)divided;   
}

//____________________________________________________________________
void AliCorrectionMatrix2D::FillMeas(Float_t ax, Float_t ay)
{
  // add value to measured histogram
  ((TH2F*)fhMeas)->Fill(ax, ay);
}

//____________________________________________________________________
void AliCorrectionMatrix2D::FillGene(Float_t ax, Float_t ay)
{
  // add value to generated histogram
  ((TH2F*)fhGene)->Fill(ax, ay);
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

//____________________________________________________________________
void AliCorrectionMatrix2D::Rebin(Int_t x, Int_t y)
{
	// rebins the histograms, recalculates the correction

        GetGeneratedHistogram()->Rebin2D(x, y);
        GetMeasuredHistogram()->Rebin2D(x, y);
        GetCorrectionHistogram()->Rebin2D(x, y);
        Divide();
}
