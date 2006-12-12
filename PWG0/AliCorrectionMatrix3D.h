#ifndef ALICORRECTIONMATRIX3D_H
#define ALICORRECTIONMATRIX3D_H

/* $Id$ */

// ------------------------------------------------------
//
// Class to handle 3d-corrections.
//
// ------------------------------------------------------

#include <AliCorrectionMatrix.h>
#include <AliCorrectionMatrix2D.h>

class TH3F;
class TH2F;
class TH1F;


class AliCorrectionMatrix3D : public AliCorrectionMatrix
{
public:
  AliCorrectionMatrix3D();
  AliCorrectionMatrix3D(const AliCorrectionMatrix3D& c);
  AliCorrectionMatrix3D(const Char_t* name, const Char_t* title,
		     Int_t nBinX, Float_t Xmin, Float_t Xmax,
		     Int_t nBinY, Float_t Ymin, Float_t Ymax,
		     Int_t nBinZ, Float_t Zmin, Float_t Zmax);

  AliCorrectionMatrix3D(const Char_t* name, const Char_t* title,
         Int_t nBinX, Float_t Xmin, Float_t Xmax,
         Int_t nBinY, Float_t Ymin, Float_t Ymax,
         Int_t nBinZ, const Float_t* zbins);

  AliCorrectionMatrix3D(const Char_t* name, const Char_t* title, TH3F* hBinning);
 
  virtual ~AliCorrectionMatrix3D();

  AliCorrectionMatrix3D& operator= (const AliCorrectionMatrix3D& c);

  void CreateHists(Int_t nBinX, const Float_t* binLimitsX,
		   Int_t nBinY, const Float_t* binLimitsY,
		   Int_t nBinZ, const Float_t* binLimitsZ);

  TH3F* GetGeneratedHistogram();
  TH3F* GetMeasuredHistogram();
  TH3F* GetCorrectionHistogram();

  AliCorrectionMatrix2D* Get2DCorrection(Char_t* opt, Float_t aMin, Float_t aMax);
  TH2F* Get2DCorrectionHistogram(Char_t* opt, Float_t aMin, Float_t aMax) {return Get2DCorrection(opt,aMin,aMax)->GetCorrectionHistogram();}
  TH1F* Get1DCorrectionHistogram(Char_t* opt, Float_t aMins1, Float_t aMax1, Float_t aMins2, Float_t aMax2);

  void FillMeas(Float_t ax, Float_t ay, Float_t az);
  void FillGene(Float_t ax, Float_t ay, Float_t az);

  Float_t GetCorrection(Float_t ax, Float_t ay, Float_t az) const;

  void RemoveEdges(Float_t cut=2, Int_t nBinsXedge = 0, Int_t nBinsYedge = 0, Int_t nBinsZedge = 0);

  virtual void SaveHistograms();

  Int_t CheckEmptyBins(Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax, Bool_t quiet = kFALSE);
  TH1F* PlotBinErrors(Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax);


protected:
  ClassDef(AliCorrectionMatrix3D,1)
};

#endif

