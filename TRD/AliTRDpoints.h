#ifndef ALITRDPOINTS_H
#define ALITRDPOINTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliPoints.h"

const Int_t kNTRpoints = 75;

class AliTRDpoints : public AliPoints {

 public:

  AliTRDpoints();
  AliTRDpoints(const AliTRDpoints &p);
  AliTRDpoints(Int_t nhitsE, Int_t nhitsT);
  virtual ~AliTRDpoints();
  AliTRDpoints &operator=(const AliTRDpoints &p);

  virtual void           Copy(TObject &p);   
  virtual void           Draw(Option_t *option);

  virtual void           SetTRpoints(Int_t n, Float_t *coor);

 protected:

  Float_t          fTRpoints[kNTRpoints];       //  The hits from TR photons
  Int_t            fNTRpoints;                  //  The number of TR photon hits
  TPolyMarker3D   *fTRpolyMarker;               //! Polymarker to draw the photon hits
  
  ClassDef(AliTRDpoints,1)                      // Class to draw TRD hits 

};
#endif
