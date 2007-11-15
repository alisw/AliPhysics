#ifndef ALITPCTRACKLET_H
#define ALITPCTRACKLET_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
// A class that contains a tracklet (a track that lives only in a single TPC
// sector).
////


#include "TObject.h"

class TObjArray;
class AliTPCseed;
class AliExternalTrackParam;

class AliTPCTracklet:public TObject {
public:
  AliTPCTracklet();
  AliTPCTracklet(const AliTPCseed *s,Int_t sector);
  AliTPCTracklet(const AliTPCTracklet &t);
  AliTPCTracklet& operator=(const AliTPCTracklet &t);
  virtual ~AliTPCTracklet();

  static TObjArray CreateTracklets(const AliTPCseed *s,
				   Int_t minClusters=0,
				   Int_t maxTracklets=72);

  // Returns the tracklet parametrisation at its outer most cluster.
  AliExternalTrackParam* GetOuter() const {return fOuter;};
  // Returns the tracklet parametrisation at its inner most cluster.
  AliExternalTrackParam* GetInner() const {return fInner;};
  // Returns the tracklet parametrisation at X=0, i.e. the "primary vertex".
  AliExternalTrackParam* GetPrimary() const {return fPrimary;};
  // Returns the sector in which the tracklet lives.
  Int_t GetSector() const {return fSector;}
  // Returns the number of clusters assined to the tracklet.
  Int_t GetNClusters() const {return fNClusters;}
private:
  Int_t fNClusters; // The number of clusters assined to the tracklet.
  Int_t fSector; // The sector this tracklet lives in.
  AliExternalTrackParam *fOuter; // The tracklet parametrisation at its outer most cluster.
  AliExternalTrackParam *fInner; // The tracklet parametrisation at its inner most cluster.
  AliExternalTrackParam *fPrimary; // The tracklet parametrisation at X=0, i.e. the "primary vertex".

  ClassDef(AliTPCTracklet,1)
};

#endif
