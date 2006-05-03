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
	     AliAlignObj::ELayerID layerRangeMin = AliAlignObj::kFirstLayer,
	     AliAlignObj::ELayerID layerRangeMax = AliAlignObj::kLastLayer);
  Bool_t GetPCA(const AliTrackPoint &p, AliTrackPoint &p2) const;

  void Reset();
  void AddPoint(Float_t x, Float_t y, Float_t z, Float_t sy, Float_t sz);
  void Update();

  Double_t GetC() const; 
  Double_t GetYat(Double_t x) const  {return fRieman->GetYat(x);}
  Double_t GetZat(Double_t x) const  {return fRieman->GetZat(x);}
  Double_t GetDYat(Double_t x) const {return fRieman->GetDYat(x);}
  Double_t GetDZat(Double_t x) const {return fRieman->GetDZat(x);}
  Bool_t   GetXYZat(Double_t r, Float_t *xyz) const {return fRieman->GetXYZat(r, fAlpha,xyz);}
  AliRieman *GetRieman() const {return fRieman;}
 protected:

  Double_t      fAlpha;     //angle to transform to the fitting coordinate system
  Int_t         fNUsed;     //actual number of space-points used in the fit
  Bool_t        fConv;      //indicates convergation
  AliRieman    *fRieman;    // rieman fitter
 private:
  TTreeSRedirector *fDebugStream;   //!debug streamer
  ClassDef(AliTrackFitterRieman,1)  // Fast fit of helices on ITS RecPoints

};

#endif
