#ifndef ALITRACKFITTERRIEMAN_H
#define ALITRACKFITTERRIEMAN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//
// Class to the track points on the Riemann sphere. Inputs are
// the set of id's (volids) of the volumes in which residuals are
// calculated to construct a chi2 function to be minimized during 
// the alignment procedures
//
//////////////////////////////////////////////////////////////////////////////

#include "AliTrackFitter.h"
#include "AliRieman.h"  
class TTreeSRedirector;
class AliRieman;

class AliTrackFitterRieman : public AliTrackFitter{
 public:
  AliTrackFitterRieman();
  AliTrackFitterRieman(AliTrackPointArray *array, Bool_t owner = kTRUE);
  AliTrackFitterRieman(const AliTrackFitterRieman &rieman);
  AliTrackFitterRieman &operator =(const AliTrackFitterRieman& rieman);
  virtual ~AliTrackFitterRieman();

  Bool_t Fit(const TArrayI *volIds,const TArrayI *volIdsFit = 0x0,
	     AliGeomManager::ELayerID layerRangeMin = AliGeomManager::kFirstLayer,
	     AliGeomManager::ELayerID layerRangeMax = AliGeomManager::kLastLayer);
  Bool_t GetPCA(const AliTrackPoint &p, AliTrackPoint &p2) const;
  void SetMaxDelta(Float_t maxDelta) { fMaxDelta = maxDelta;}
  Float_t GetMaxDelta() const { return fMaxDelta;}
  void  SetCorrection(Bool_t correction){ fBCorrection=correction;}
  Bool_t  GetCorrection() const {return fBCorrection ;}
  void Reset();
  void AddPoint(Float_t x, Float_t y, Float_t z, Float_t sy, Float_t sz);
  Bool_t Update();

  Double_t GetC() const              {return fRieman->GetC();}
  Double_t GetYat(Double_t x) const;
  Double_t GetZat(Double_t x) const;
  Double_t GetDYat(Double_t x) const {return fRieman->GetDYat(x);}
  Double_t GetDZat(Double_t x) const {return fRieman->GetDZat(x);}  
  Double_t GetErrY2at(Double_t x) const;
  Double_t GetErrZ2at(Double_t x) const;

  Bool_t   GetXYZat(Double_t r, Float_t *xyz) const {return fRieman->GetXYZat(r, fAlpha,xyz);}
  AliRieman *GetRieman() const {return fRieman;}

  // Methods needed in order to use the class
  // inside the millipede alignment algorithm
  Double_t   GetAlpha() const { return fAlpha; }
  void SetParam(Int_t i, Double_t par);

  void SetMinPointRadius(Double_t radius) { fMinPointRadius = radius; }
  void SetMaxPointRadius(Double_t radius) { fMaxPointRadius = radius; }
  Double_t GetMinPointRadius() const { return fMinPointRadius; }
  Double_t GetMaxPointRadius() const { return fMaxPointRadius; }

 protected:
  Bool_t        fBCorrection; //add  correction for non-helicity
  Double_t      fAlpha;     //angle to transform to the fitting coordinate system
  Int_t         fNUsed;     //actual number of space-points used in the fit
  Bool_t        fConv;      //indicates convergation
  Float_t       fMaxDelta;  // maximal allowed delta in PCA exported for PCA minimization
  AliRieman    *fRieman;    // rieman fitter
  Double_t      fCorrY[4];  // correction polynom coef
  Double_t      fCorrZ[4];  // correction polynom coef
  Double_t      fMinPointRadius; // minimal allowed radius for the point to be included in the fit
  Double_t      fMaxPointRadius; // maximal allowed radius for the point to be included in the fit
 private:
  Bool_t Begin(Int_t, Int_t) {Reset(); return kTRUE;}
  Bool_t AddPoint(const AliTrackPoint *) {return kTRUE;}
  TTreeSRedirector *fDebugStream;   //!debug streamer
  ClassDef(AliTrackFitterRieman,3)  // Fast fit of helices on ITS RecPoints

};

#endif
