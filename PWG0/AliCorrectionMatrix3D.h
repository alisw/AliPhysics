#ifndef ALICORRECTIONMATRIX3D_H
#define ALICORRECTIONMATRIX3D_H

/* $Id$ */

// ------------------------------------------------------
//
// Class to handle 3d-corrections.
//
// ------------------------------------------------------

#include <AliCorrectionMatrix.h>

class TH3F;
class TH1F;

class AliCorrectionMatrix3D : public AliCorrectionMatrix
{
public:
  AliCorrectionMatrix3D();
  AliCorrectionMatrix3D(const AliCorrectionMatrix3D& c);
  AliCorrectionMatrix3D(const Char_t* name, const Char_t* title,
		     Int_t nBinX=10, Float_t Xmin=0., Float_t Xmax=10.,
		     Int_t nBinY=10, Float_t Ymin=0., Float_t Ymax=10.,
		     Int_t nBinZ=10, Float_t Zmin=0., Float_t Zmax=10.);

  AliCorrectionMatrix3D(const Char_t* name, const Char_t* title,
         Int_t nBinX, Float_t Xmin, Float_t Xmax,
         Int_t nBinY, Float_t Ymin, Float_t Ymax,
         Int_t nBinZ, const Float_t* zbins);

  virtual ~AliCorrectionMatrix3D();

  AliCorrectionMatrix3D& operator= (const AliCorrectionMatrix3D& c);

  TH3F* GetGeneratedHistogram();
  TH3F* GetMeasuredHistogram();
  TH3F* GetCorrectionHistogram();

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

