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
class AliTPCclusterMI;

#include "TEllipse.h"

class AliTPCTracklet:public TObject {
public: 
  enum TrackType {kKalman,kRiemann,kLinear,kQuadratic};

  AliTPCTracklet();
  AliTPCTracklet(const AliTPCseed *s,Int_t sector,TrackType type=kKalman,
		 Bool_t storeClusters=kFALSE);
  AliTPCTracklet(const TObjArray &clusters,Int_t sector,TrackType type=kKalman,
		 Bool_t storeClusters=kFALSE);
  AliTPCTracklet(const AliTPCTracklet &t);
  AliTPCTracklet& operator=(const AliTPCTracklet &t);
  virtual ~AliTPCTracklet();

  static TObjArray CreateTracklets(const TObjArray &clusters,
				   TrackType type=kKalman,
				   Bool_t storeClusters=kFALSE,
				   Int_t minClusters=0,
				   Int_t maxTracklets=72);

  static TObjArray CreateTracklets(const AliTPCseed *s,
				   TrackType type=kKalman,
				   Bool_t storeClusters=kFALSE,
				   Int_t minClusters=0,
				   Int_t maxTracklets=72);

  static Bool_t PropagateToMeanX(const AliTPCTracklet &t1,
				 const AliTPCTracklet &t2,
				 AliExternalTrackParam *&t1m,
				 AliExternalTrackParam *&t2m);

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
  // Returns the clusters of this tracklet. In case they weren't stored it
  // returns 0.
  AliTPCclusterMI* GetClusters() const {return fClusters;};
  // Test the functionality of the class. Generates some random tracks and
  // refits them into tracklets. 
  static void Test(const char *filename);
  static void RandomND(Int_t ndim,const Double_t *p,const Double_t *c,
		       Double_t *x);
  static TEllipse ErrorEllipse(Double_t x,Double_t y,
			       Double_t sx,Double_t sy,Double_t sxy);
  static inline void SetEdgeCut(Float_t edgeX, Float_t edgeY);
private:
  static Bool_t RejectCluster(AliTPCclusterMI* cl,AliExternalTrackParam * param=0);
  static const Double_t kB2C; //! ugly to have the track parametrised in a way, that constand is allways needed
  static double GetBz(Double_t *xyz);
  static Float_t        fgEdgeCutY; //cut on the edge effect in local Y 
  static Float_t        fgEdgeCutX; //cut on the edge effect in local X 
  void FitLinear(const AliTPCseed *track,Int_t sector,TrackType type);
  void FitKalman(const AliTPCseed *track,Int_t sector);
  void FitRiemann(const AliTPCseed *track,Int_t sector);
  void Quadratic2Helix(Double_t *a,Double_t *ca,
		       Double_t *b,Double_t *cb,
		       Double_t x0,
		       Double_t *p,Double_t *c);
  Bool_t Riemann2Helix(Double_t *a,Double_t *ca,
		       Double_t *b,Double_t *cb,
		       Double_t x0,
		       Double_t *p,Double_t *c);   
  Int_t fNClusters; // The number of clusters assined to the tracklet.
  Int_t fNStoredClusters; // The number of stored clusters.
  AliTPCclusterMI *fClusters; //[fNStoredClusters] The clusters of the track, if stored (otherwise 0)
  Int_t fSector; // The sector this tracklet lives in.
  AliExternalTrackParam *fOuter; // The tracklet parametrisation at its outer most cluster.
  AliExternalTrackParam *fInner; // The tracklet parametrisation at its inner most cluster.
  AliExternalTrackParam *fPrimary; // The tracklet parametrisation at X=0, i.e. the "primary vertex".

  ClassDef(AliTPCTracklet,1)
};


void AliTPCTracklet::SetEdgeCut(Float_t edgeX, Float_t edgeY){
  //
  //
  fgEdgeCutY=edgeY;
  fgEdgeCutX=edgeX;
}

#endif
