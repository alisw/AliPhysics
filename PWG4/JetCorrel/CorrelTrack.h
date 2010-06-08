#ifndef CORRELTRACK_H
#define CORRELTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_____________________________________________________
// Container class for global tracks
// Use CorrelKFTrack_t to reconstruct parent with AliKF
//-- Author: Paul Constantin
 
#include "CorrelParticle.h"

class CorrelTrack_t : public CorrelParticle_t {
 public:
  
  CorrelTrack_t();
  CorrelTrack_t(Float_t pt, Float_t p, Float_t e, Float_t m, cPartType_t i, Float_t x, Float_t y, Float_t z);
  CorrelTrack_t(const CorrelTrack_t &p);
  virtual ~CorrelTrack_t() {;}
  virtual CorrelTrack_t* Copy();
  
  Float_t X() const {return fTPCx;}
  Float_t Y() const {return fTPCy;}
  Float_t Z() const {return fTPCz;}
  Float_t Dist(CorrelTrack_t * const trk) const;
  
  void SetTPCEntry(Float_t x, Float_t y, Float_t z) {fTPCx=x; fTPCy=y; fTPCz=z;}
  virtual void Show() const;
  
 private:
  Float_t fTPCx;   // x-coord TPC entrance
  Float_t fTPCy;   // y-coord TPC entrance
  Float_t fTPCz;   // z-coord TPC entrance
};

#endif
