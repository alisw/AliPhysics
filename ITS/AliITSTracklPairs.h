#ifndef ALIITSTRACKLPAIRS_H 
#define ALIITSTRACKLPAIRS_H

#include<TObject.h>
/* Copyright(c) 2009-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////
// Helper class for 3D primary vertexing                      //
// Used by AliITSSortTrkl                                     //
// Origin M.Masera (masera@to.infn.it)                        //
////////////////////////////////////////////////////////////////

class AliITSTracklPairs : public TObject {

 public:

  AliITSTracklPairs();
  AliITSTracklPairs(Int_t t1, Int_t t2, Double_t dca, Double_t *coo);
  virtual ~AliITSTracklPairs();
  Int_t GetTrack1() const {return fTrack1;}
  Int_t GetTrack2() const {return fTrack2;}
  Double_t GetDCA() const {return fDCA;}
  void GetCrossCoord(Double_t *cr) const {for(int i=0;i<3;i++)cr[i]=fCross[i];}
  Double_t GetDistance(const AliITSTracklPairs& pair) const;
  Bool_t HasTrack(Int_t tr) const {return ((tr == fTrack1) || (tr == fTrack2));}

 protected:
  Int_t fTrack1;      // first tracklet index
  Int_t fTrack2;      // second tracklet index
  Double_t fDCA;      // DCA
  Double_t fCross[3]; // intersection coordinates

  ClassDef(AliITSTracklPairs,1);
};

#endif
