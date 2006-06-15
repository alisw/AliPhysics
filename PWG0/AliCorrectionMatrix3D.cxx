/* $Id$ */

// ------------------------------------------------------
//
// Class to handle 3d-corrections.
//
// ------------------------------------------------------
//

#include <TH3F.h>

#include <AliLog.h>

#include "AliCorrectionMatrix3D.h"
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
AliCorrectionMatrix3D::AliCorrectionMatrix3D(const Char_t* name, const Char_t* title,
              Int_t nBinX, Float_t Xmin, Float_t Xmax,
              Int_t nBinY, Float_t Ymin, Float_t Ymax,
              Int_t nBinZ, Float_t Zmin, Float_t Zmax)
  : AliCorrectionMatrix(name, title)
{
  //
  // constructor
  //

  fhMeas  = new TH3F(Form("meas_%s",name), Form("meas_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax, nBinZ, Zmin, Zmax);
  fhGene  = new TH3F(Form("gene_%s",name), Form("gene_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax, nBinZ, Zmin, Zmax);
  fhCorr  = new TH3F(Form("corr_%s",name), Form("corr_%s",title),  nBinX, Xmin, Xmax, nBinY, Ymin, Ymax, nBinZ, Zmin, Zmax);

  fhMeas->Sumw2();
  fhGene->Sumw2();
  fhCorr->Sumw2();
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

  fhMeas  = new TH3F(Form("meas_%s",name), Form("meas_%s",title), nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, zbins);
  fhGene  = new TH3F(Form("gene_%s",name), Form("gene_%s",title), nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, zbins);
  fhCorr  = new TH3F(Form("corr_%s",name), Form("corr_%s",title), nBinX, binLimitsX, nBinY, binLimitsY, nBinZ, zbins);

  delete[] binLimitsX;
  delete[] binLimitsY;

  fhMeas->Sumw2();
  fhGene->Sumw2();
  fhCorr->Sumw2();
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
    AliPWG0Helper::CreateDividedProjections(GetGeneratedHistogram(), GetMeasuredHistogram());
}
