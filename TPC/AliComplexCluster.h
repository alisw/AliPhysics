#ifndef ALICOMPLEXCLUSTER_H
#define ALICOMPLEXCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"

class AliComplexCluster : public TObject {
public:
  Int_t     fTracks[3];//labels of overlapped tracks
  Float_t   fX ;       //Y of cluster
  Float_t   fY ;       //Z of cluster
  Float_t   fQ ;       //Q of cluster (in ADC counts)
  Float_t   fSigmaX2;  //Sigma Y square of cluster
  Float_t   fSigmaY2;  //Sigma Z square of cluster
  Float_t   fSigmaXY;  //      XY moment 
  Float_t   fArea;     //area of cluster
  Float_t   fMax;     //amplitude at maximum 
public:
  AliComplexCluster() {
    fTracks[0]=fTracks[1]=fTracks[2]=0; 
    fX=fY=fQ=fSigmaX2=fSigmaY2=0.;
  }
  virtual ~AliComplexCluster() {;}
  Bool_t    IsSortable() const;
  Int_t Compare(const TObject *o) const;
  ClassDef(AliComplexCluster,1)  // Cluster manager
};

class AliDigitCluster : public AliComplexCluster {
public:
  Int_t fNx; //number of accepted x bins
  Int_t fNy; //number of accepted y bins
  Float_t fMaxX; //maximum x bin
  Float_t fMaxY; //maximum y bin
public:  
  ClassDef(AliDigitCluster,1)  // Tclusters
};

class AliDifCluster : public AliDigitCluster {
public:
  Float_t fDx; //delta x 
  Float_t fDy; //delta y
  Float_t fAngleX;//x angle
  Float_t fAngleY;//y angle
  Float_t fOrigQ; //original charge
  Int_t fGenerTrack;  //track number of nearest track
  ClassDef(AliDifCluster,1)  // // Cluster manager
};

#endif //ALICOMPLEXCLUSTER_H
