#ifndef ALICORRECTIONMATRIX2D_H
#define ALICORRECTIONMATRIX2D_H

/* $Id$ */

// ------------------------------------------------------
//
// Class to handle 2d-corrections.
//
// ------------------------------------------------------

#include <AliCorrectionMatrix.h>

class TH2F;
class TH1F;

class AliCorrectionMatrix2D : public AliCorrectionMatrix
{
public:
  AliCorrectionMatrix2D();
  AliCorrectionMatrix2D(const AliCorrectionMatrix2D& c);
  AliCorrectionMatrix2D(const Char_t* name, const Char_t* title,
         Int_t nBinX, Float_t Xmin, Float_t Xmax,
         Int_t nBinY, Float_t Ymin, Float_t Ymax);

  AliCorrectionMatrix2D(const Char_t* name, const Char_t* title,
         Int_t nBinX, Float_t *X,
         Int_t nBinY, Float_t *Y);

  virtual ~AliCorrectionMatrix2D();

  AliCorrectionMatrix2D& operator= (const AliCorrectionMatrix2D& c);

  TH2* GetGeneratedHistogram() const;
  TH2* GetMeasuredHistogram() const;

  TH2* GetCorrectionHistogram() {return (TH2*)fhCorr;}

  TH1* Get1DCorrection(const Char_t* opt="x", Float_t min=0, Float_t max=0) {return Get1DCorrectionHistogram(opt,min,max);}
  TH1* Get1DCorrectionHistogram(const Char_t* opt="x", Float_t min=0, Float_t max=0, Bool_t binomialErrors = kFALSE);

  void Rebin(Int_t x = 1, Int_t y = 1);

  void FillMeas(Float_t ax, Float_t ay);
  void FillGene(Float_t ax, Float_t ay);
  Float_t GetCorrection(Float_t ax, Float_t ay) const;

  void RemoveEdges(Float_t cut=2, Int_t nBinsX=0, Int_t nBinsY=0);

protected:
  ClassDef(AliCorrectionMatrix2D,1)

private:
};

#endif

