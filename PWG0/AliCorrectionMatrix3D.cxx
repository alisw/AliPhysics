/* $Id$ */

// ------------------------------------------------------
//
// Class to handle 3d-corrections.
//
// ------------------------------------------------------
//

#include <TDirectory.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>

#include <AliLog.h>

#include "AliCorrectionMatrix3D.h"
#include "AliCorrectionMatrix2D.h"
#include "AliPWG0Helper.h"

//____________________________________________________________________
ClassImp(AliCorrectionMatrix3D)

//____________________________________________________________________
AliCorrectionMatrix3D::AliCorrectionMatrix3D() :
  AliCorrectionMatrix()
{
  // default constructor
}

//____________________________________________________________________
AliCorrectionMatrix3D::AliCorrectionMatrix3D(const AliCorrectionMatrix3D& c)
  : AliCorrectionMatrix(c)
{
  // copy constructor
  ((AliCorrectionMatrix3D &)c).Copy(*this);
}

//____________________________________________________________________
AliCorrectionMatrix3D &AliCorrectionMatrix3D::operator=(const AliCorrectionMatrix3D &c)
{
  // assigment operator

  if (this != &c)
    ((AliCorrectionMatrix3D &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
AliCorrectionMatrix3D::AliCorrectionMatrix3D(const Char_t* name, const Char_t* title,
              Int_t nBinX, Float_t Xmin, Float_t Xmax,
              Int_t nBinY, Float_t Ymin, Float_t Ymax,
              Int_t nBinZ, Float_t Zmin, Float_t Zmax)
  : AliCorrectionMatrix(name, title)
{
  //
  // constructor
  //

  Float_t* binLimitsX = new Float_t[nBinX+1];
  for (Int_t i=0; i<=nBinX; ++i)
    binLimitsX[i] = Xmin + (Xmax - Xmin) / nBinX * i;

  Float_t* binLimitsY = new Float_t[nBinY+1];
  for (Int_t i=0; i<=nBinY; ++i)
    binLimitsY[i] = Ymin + (Ymax - Ymin) / nBinY * i;

  Float_t* binLimitsZ = new Float_t[nBinZ+1];
  for (Int_t i=0; i<=nBinZ; ++i)
    binLimitsZ[i] = Zmin + (Zmax - Zmin) / nBinZ * i;

  CreateHists(nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, binLimitsZ);

  delete[] binLimitsX;
  delete[] binLimitsY;
  delete[] binLimitsZ;
}

AliCorrectionMatrix3D::AliCorrectionMatrix3D(const Char_t* name, const Char_t* title,
      Int_t nBinX, Float_t Xmin, Float_t Xmax,
      Int_t nBinY, Float_t Ymin, Float_t Ymax,
      Int_t nBinZ, const Float_t* zbins)
  : AliCorrectionMatrix(name, title)
{
  // constructor with variable bin sizes

  Float_t* binLimitsX = new Float_t[nBinX+1];
  for (Int_t i=0; i<=nBinX; ++i)
    binLimitsX[i] = Xmin + (Xmax - Xmin) / nBinX * i;

  Float_t* binLimitsY = new Float_t[nBinY+1];
  for (Int_t i=0; i<=nBinY; ++i)
    binLimitsY[i] = Ymin + (Ymax - Ymin) / nBinY * i;

  CreateHists(nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, zbins);

  delete[] binLimitsX;
  delete[] binLimitsY;
}

AliCorrectionMatrix3D::AliCorrectionMatrix3D(const Char_t* name, const Char_t* title, TH3F* hBinning)
  : AliCorrectionMatrix(name, title)
{
  // constructor with variable bin sizes (uses binning of hBinning)

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhMeas = (TH3F*)hBinning->Clone("measured");
  fhGene = (TH3F*)hBinning->Clone("generated");
  fhCorr = (TH3F*)hBinning->Clone("correction");

  fhMeas->SetTitle(Form("%s measured",title)  );
  fhGene->SetTitle(Form("%s generated",title ));
  fhCorr->SetTitle(Form("%s correction",title));

  fhMeas->Reset();
  fhGene->Reset();
  fhCorr->Reset();

  TH1::AddDirectory(oldStatus);

  fhMeas->Sumw2();
  fhGene->Sumw2();
  fhCorr->Sumw2();
}

//____________________________________________________________________
void AliCorrectionMatrix3D::CreateHists(Int_t nBinX, const Float_t* binLimitsX,
      Int_t nBinY, const Float_t* binLimitsY,
      Int_t nBinZ, const Float_t* binLimitsZ)
{
  // create the histograms

  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhMeas  = new TH3F("measured",   Form("%s measured",GetTitle()),   nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, binLimitsZ);
  fhGene  = new TH3F("generated",  Form("%s generated",GetTitle()),  nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, binLimitsZ);
  fhCorr  = new TH3F("correction", Form("%s correction",GetTitle()), nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, binLimitsZ);

  fhMeas->Sumw2();
  fhGene->Sumw2();
  fhCorr->Sumw2();

  TH1::AddDirectory(oldStatus);
}


//____________________________________________________________________
AliCorrectionMatrix3D::~AliCorrectionMatrix3D()
{
  //
  // destructor
  //

  // histograms already deleted in base class
}

//____________________________________________________________________
TH3F* AliCorrectionMatrix3D::GetGeneratedHistogram()
{
  // return generated histogram casted to correct type
  return dynamic_cast<TH3F*> (fhGene);
}

//____________________________________________________________________
TH3F* AliCorrectionMatrix3D::GetMeasuredHistogram()
{
  // return measured histogram casted to correct type
  return dynamic_cast<TH3F*> (fhMeas);
}

//____________________________________________________________________
TH3F* AliCorrectionMatrix3D::GetCorrectionHistogram()
{
  // return correction histogram casted to correct type
  return dynamic_cast<TH3F*> (fhCorr);
}

//____________________________________________________________________
AliCorrectionMatrix2D* AliCorrectionMatrix3D::Get2DCorrection(Char_t* opt, Float_t aMin, Float_t aMax)
{
  // returns a 2D projection of this correction

  TString option = opt;

  if (aMin<aMax) {
    if (option.Contains("xy") || option.Contains("yx")) {
      Int_t bMin = fhMeas->GetZaxis()->FindBin(aMin);
      Int_t bMax = fhMeas->GetZaxis()->FindBin(aMax);
      fhMeas->GetZaxis()->SetRange(bMin, bMax);      
    }
    else if (option.Contains("xz") || option.Contains("zx")) {
      Int_t bMin = fhMeas->GetYaxis()->FindBin(aMin);
      Int_t bMax = fhMeas->GetYaxis()->FindBin(aMax);
      fhMeas->GetYaxis()->SetRange(bMin, bMax);      
    }
    else if (option.Contains("yz") || option.Contains("zy")) {
      Int_t bMin = fhMeas->GetXaxis()->FindBin(aMin);
      Int_t bMax = fhMeas->GetXaxis()->FindBin(aMax);
      fhMeas->GetXaxis()->SetRange(bMin, bMax);      
    }
    else {
      AliDebug(AliLog::kWarning, Form("WARNING: unknown projection option %s", opt));
      return 0;
    }
  }
  AliCorrectionMatrix2D* corr2D = new AliCorrectionMatrix2D(Form("%s_%s",GetName(),opt),Form("%s projection %s",GetName(),opt),100,0,100,100,0,100);

  TH2F* meas = (TH2F*) ((TH3F*)fhMeas)->Project3D(opt);
  TH2F* gene = (TH2F*) ((TH3F*)fhGene)->Project3D(opt);

  TH2F* corr = (TH2F*)gene->Clone("corr");
  corr->Reset();

  corr2D->SetGeneratedHistogram(gene);
  corr2D->SetMeasuredHistogram(meas);
  corr2D->SetCorrectionHistogram(corr);

  corr2D->Divide();

  // unzoom
  fhMeas->GetXaxis()->UnZoom();  
  fhMeas->GetYaxis()->UnZoom();  
  fhMeas->GetZaxis()->UnZoom();  

  fhGene->GetXaxis()->UnZoom();
  fhGene->GetYaxis()->UnZoom();
  fhGene->GetZaxis()->UnZoom();

  return corr2D;
}

//____________________________________________________________________
TH1F* AliCorrectionMatrix3D::Get1DCorrectionHistogram(Char_t* opt, Float_t aMin1, Float_t aMax1, Float_t aMin2, Float_t aMax2)
{
  // returns a 1D projection of this correction
  AliDebug(AliLog::kWarning, Form("WARNING: test"));
  
  AliCorrectionMatrix2D* corr2D;
  if (strcmp(opt,"x")==0) {  
    corr2D = Get2DCorrection("xy",aMin1,aMax1);
    return corr2D->Get1DCorrectionHistogram("x",aMin2,aMax2);
  }
  if (strcmp(opt,"y")==0) {
    corr2D = Get2DCorrection("xy",aMin1,aMax1);
    return corr2D->Get1DCorrectionHistogram("y",aMin2,aMax2);
  }  
  if (strcmp(opt,"z")==0) {
    corr2D = Get2DCorrection("yz",aMin1,aMax1);    
    return corr2D->Get1DCorrectionHistogram("x",aMin2,aMax2);
  }  
  AliDebug(AliLog::kWarning, Form("WARNING: unknown projection option %s (should be x,y or z)", opt));  
  
  return 0;
}



//____________________________________________________________________
void AliCorrectionMatrix3D::FillMeas(Float_t ax, Float_t ay, Float_t az)
{
  // add value to measured histogram
  GetMeasuredHistogram()->Fill(ax, ay, az);
}

//____________________________________________________________________
void AliCorrectionMatrix3D::FillGene(Float_t ax, Float_t ay, Float_t az)
{
  // add value to generated histogram
  GetGeneratedHistogram()->Fill(ax, ay, az);
}

//____________________________________________________________________
Float_t AliCorrectionMatrix3D::GetCorrection(Float_t ax, Float_t ay, Float_t az) const
{
  // returns a value of the correction map
  return fhCorr->GetBinContent(fhCorr->FindBin(ax, ay, az));
}

//____________________________________________________________________
//void AliCorrectionMatrix3D::RemoveEdges(Float_t cut, Int_t nBinsXedge, Int_t nBinsYedge, Int_t nBinsZedge)
void AliCorrectionMatrix3D::RemoveEdges(Float_t, Int_t, Int_t, Int_t)
{
  // so what do we do here...
}

//____________________________________________________________________
void AliCorrectionMatrix3D::SaveHistograms()
{
  //
  // saves the histograms
  //

  AliCorrectionMatrix::SaveHistograms();

  if (GetGeneratedHistogram() && GetMeasuredHistogram())
  {
    gDirectory->cd(GetName());
    
    AliPWG0Helper::CreateDividedProjections(GetGeneratedHistogram(), GetMeasuredHistogram(), 0, kFALSE, kTRUE);
    
    gDirectory->cd("..");
  }
}

//____________________________________________________________________
Int_t AliCorrectionMatrix3D::CheckEmptyBins(Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax, Bool_t quiet)
{
  //
  // counts the number of empty Bins in a given region
  //

  TH3F* hist = GetGeneratedHistogram();
  if (!hist)
    return -1;

  Int_t emptyBins = 0;
  for (Int_t x=hist->GetXaxis()->FindBin(xmin); x<=hist->GetXaxis()->FindBin(xmax); ++x)
    for (Int_t y=hist->GetYaxis()->FindBin(ymin); y<=hist->GetYaxis()->FindBin(ymax); ++y)
      for (Int_t z=hist->GetZaxis()->FindBin(zmin); z<=hist->GetZaxis()->FindBin(zmax); ++z)
        if (hist->GetBinContent(x, y, z) == 0)
        {
          if (!quiet)
            printf("Empty bin in %s at vtx = %f, eta = %f, pt = %f\n", GetName(), hist->GetXaxis()->GetBinCenter(x), hist->GetYaxis()->GetBinCenter(y), hist->GetZaxis()->GetBinCenter(z));
          ++emptyBins;
        }

  return emptyBins;
}

//____________________________________________________________________
TH1F* AliCorrectionMatrix3D::PlotBinErrors(Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax)
{
  //
  // makes a 1d plots of the relative bin errors
  //

  TH3F* hist = GetCorrectionHistogram();
  if (!hist)
    return 0;
    
  TH1F* target = new TH1F("relerrors", "relerrors", 100, 0, 10);

  for (Int_t x=hist->GetXaxis()->FindBin(xmin); x<=hist->GetXaxis()->FindBin(xmax); ++x)
    for (Int_t y=hist->GetYaxis()->FindBin(ymin); y<=hist->GetYaxis()->FindBin(ymax); ++y)
      for (Int_t z=hist->GetZaxis()->FindBin(zmin); z<=hist->GetZaxis()->FindBin(zmax); ++z)
        if (hist->GetBinContent(x, y, z) != 0)
          target->Fill(100 * hist->GetBinError(x, y, z) / hist->GetBinContent(x, y, z));

  return target;
}

