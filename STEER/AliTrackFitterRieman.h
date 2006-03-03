#ifndef ALITRACKFITTERRIEMAN_H
#define ALITRACKFITTERRIEMAN_H

#include "TMatrixDSym.h"

#include "AliTrackFitter.h"

class AliTrackFitterRieman : public AliTrackFitter{
 public:
  AliTrackFitterRieman();
  AliTrackFitterRieman(AliTrackPointArray *array, Bool_t owner = kTRUE);
  AliTrackFitterRieman(const AliTrackFitterRieman &rieman);
  AliTrackFitterRieman &operator =(const AliTrackFitterRieman& rieman);
  virtual ~AliTrackFitterRieman();

  Bool_t Fit(const TArrayI *volIds,const TArrayI *volIdsFit = 0x0,
	     AliAlignObj::ELayerID layerRangeMin = AliAlignObj::kFirstLayer,
	     AliAlignObj::ELayerID layerRangeMax = AliAlignObj::kLastLayer);
  Bool_t GetPCA(const AliTrackPoint &p, AliTrackPoint &p2) const;

  void Reset();
  void AddPoint(Float_t x, Float_t y, Float_t z, Float_t sy, Float_t sz);
  void Update();

  Double_t GetC(); 
  Double_t GetYat(Double_t x);
  Double_t GetZat(Double_t x) const;
  Double_t GetDYat(Double_t x) const;
  Double_t GetDZat(Double_t x) const;
  Bool_t   GetXYZat(Double_t r, Float_t *xyz);

 protected:

  Double_t      fAlpha;     //angle to transform to the fitting coordinate system
  Double_t      fSumXY[9];  //sums for XY part
  Double_t      fSumYY;     //sum for YY part
  Double_t      fSumXZ[9];  //sums for XZ part
  Double_t      fSumZZ;     //sum for ZZ part
  Int_t         fNUsed;     //actual number of space-points used in the fit
  Bool_t        fConv;      // indicates convergation
  Float_t      *fX;
  Float_t      *fY;
  Float_t      *fZ;
  Float_t      *fSy;
  Float_t      *fSz;

 private:

  ClassDef(AliTrackFitterRieman,1)  // Fast fit of helices on ITS RecPoints

};

#endif
