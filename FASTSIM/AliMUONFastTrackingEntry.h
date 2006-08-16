#ifndef ALIMUONFASTTRACKINGENTRY
#define ALIMUONFASTTRACKINGENTRY


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TClassTable.h>

static const Int_t kSplitP = 5; 
static const Int_t kSplitTheta = 3; 

class AliMUONFastTrackingEntry {
 public:
  AliMUONFastTrackingEntry();
  virtual ~AliMUONFastTrackingEntry(){;}
  Float_t fP;
  Float_t fTheta;
  Float_t fPhi;
  Float_t fMeanp;
  Float_t fMeantheta;
  Float_t fMeanphi;
  Float_t fSigmap;
  Float_t fSigmatheta;
  Float_t fSigmaphi;
  Float_t fSigma1p;
  Float_t fChi2p;
  Float_t fChi2theta;
  Float_t fChi2phi;
  Float_t fAcc[5][3];
  Float_t fEff[5][3];
  Float_t fNormG2;
  Float_t fMeanG2;
  Float_t fSigmaG2;
  ClassDef(AliMUONFastTrackingEntry,1)       
};


#endif
