#ifndef ALITPCTRACKERPARAM_H
#define ALITPCTRACKERPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice                               */
/* $Id$ */
//-----------------------------------------------------------------------------
//                    TPC Tracking Parameterization Class
//
//   Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it
//-----------------------------------------------------------------------------
#include "alles.h"
#include "AliMagF.h"
#include "AliTPCtrack.h"

class AliTPCtrackerParam {
 public:
  AliTPCtrackerParam(const Int_t coll=0,const Double_t Bz=0.4);
  virtual ~AliTPCtrackerParam();


  Int_t BuildTPCtracks(const TFile *inp, TFile *out,Int_t n=1);

 private:
  Int_t    fColl; // collision code (0: PbPb6000)
  Double_t fBz;   // value of the z component of L3 field (Tesla)
 
  
  AliTPCtrack* BuildTrack(Double_t alpha,Double_t x,Double_t y,Double_t z,
			  Double_t Px,Double_t Py,Double_t Pz,Double_t Pt,
			  Int_t ch,Int_t lab);
  
  Bool_t SelectedTrack(Int_t pdg, Double_t pt, Double_t eta);
  
  Int_t GetBin(Double_t pt,Double_t eta);
  
  TMatrixD GetSmearingMatrix(Double_t* cc, Double_t pt,Double_t eta);
  
  void SmearTrack(Double_t* xx,Double_t* xxsm,TMatrixD cov);

  Double_t LinearInterpolation(Int_t ptBins,Double_t *value,Double_t pt,Double_t eta);
  
  void CookTracks(TObjArray& tarray,TObjArray& newtarray);
  
  
  ClassDef(AliTPCtrackerParam,1)
};

#endif


