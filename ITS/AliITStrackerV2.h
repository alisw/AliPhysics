#ifndef ALIITSTRACKER_H
#define ALIITSTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          ITS tracker
//
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
#include "AliTracker.h"
#include "AliITSrecoV2.h"
#include "AliITStrackV2.h"

class AliITSclusterV2;
class AliITSgeom;
class TFile;


//-------------------------------------------------------------------------
class AliITStrackerV2 : public AliTracker {
public:
  AliITStrackerV2():AliTracker(){}
  AliITStrackerV2(const AliITSgeom *geom) throw (const Char_t *);

  AliCluster *GetCluster(Int_t index) const;
  Int_t Clusters2Tracks(const TFile *in, TFile *out);
  Int_t PropagateBack(const TFile *in, TFile *out) {return 0;}

private:

  Double_t GetEffectiveThickness(Double_t phi, Double_t z) const;

  void  FollowProlongation();
  Int_t TakeNextProlongation();

  void ResetBestTrack() {
     fBestTrack.~AliITStrackV2();
     new(&fBestTrack) AliITStrackV2(fTrackToFollow);
  }

  void ResetTrackToFollow(const AliITStrackV2 &t) {
     fTrackToFollow.~AliITStrackV2();
     new(&fTrackToFollow) AliITStrackV2(t);
  }

class AliITSdetector {
private:
  Double_t fR;    // polar coordinates 
  Double_t fPhi;  // of this detector

public:
  AliITSdetector(){}
  AliITSdetector(Double_t r,Double_t phi) {fR=r; fPhi=phi;}

  void *operator new(size_t s,AliITSdetector *p) {return p;}

  Double_t GetR()   const {return fR;}
  Double_t GetPhi() const {return fPhi;}
};

class AliITSlayer {
  Double_t fR;                // mean radius of this layer
  Double_t fPhiOffset;        // offset of the first detector in Phi
  Int_t fNladders;            // number of ladders
  Double_t fZOffset;          // offset of the first detector in Z
  Int_t fNdetectors;          // detectors/ladder
  AliITSdetector *fDetectors; // array of detectors

  Int_t fN;                   // number of clusters
  AliITSclusterV2 *fClusters[kMaxClusterPerLayer];    // pointers to clusters

  Double_t fZmax;      //       edges
  Double_t fYmin;      //      of  the
  Double_t fYmax;      //      "window"
  Int_t fI;            // index of the current cluster within the "window"

  Int_t FindClusterIndex(Double_t z) const;

public:
  AliITSlayer();
  AliITSlayer(Double_t r, Double_t p, Double_t z, Int_t nl, Int_t nd);
 ~AliITSlayer();
  Int_t InsertCluster(AliITSclusterV2 *c);
  void SelectClusters(Double_t zmin,Double_t zmax,Double_t ymin,Double_t ymax);
  const AliITSclusterV2 *GetNextCluster(Int_t &ci);

  void *operator new(size_t s, AliITSlayer *p) {return p;}

  Double_t GetR() const {return fR;}
  AliITSclusterV2 *GetCluster(Int_t i) const {return fClusters[i];} 
  AliITSdetector &GetDetector(Int_t n) const { return fDetectors[n]; }
  Int_t FindDetectorIndex(Double_t phi, Double_t z) const;
  Double_t GetThickness(Double_t phi, Double_t z) const;

  Int_t InRoad() const ;
};

  Int_t fI;                       // index of the current layer
  static AliITSlayer fLayers[kMaxLayer]; // ITS layers
  AliITStrackV2 fTracks[kMaxLayer]; // track estimations at the ITS layers
  AliITStrackV2 fBestTrack;         // "best" track 
  AliITStrackV2 fTrackToFollow;     // followed track

  Double_t fYV;                   // Y-coordinate of the primary vertex
  Double_t fZV;                   // Z-coordinate of the primary vertex
};



#endif
