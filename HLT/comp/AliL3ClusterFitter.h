// @(#) $Id$

#ifndef AliL3_ClusterFitter
#define AliL3_ClusterFitter

#include "AliL3Modeller.h"
#include "AliL3Transform.h"

class AliL3ModelTrack;
class AliL3TrackArray;
class AliL3SpacePointData;

class AliL3ClusterFitter : public AliL3Modeller {
  
 public:
  AliL3ClusterFitter();
  AliL3ClusterFitter(Char_t *path);
  virtual ~AliL3ClusterFitter();
  
  void Init(Int_t slice,Int_t patch,Int_t *rowrange,AliL3TrackArray *tracks);
  void Init(Int_t slice,Int_t patch);
  void LoadSeeds(Int_t *rowrange,Bool_t offline=kTRUE,Int_t eventnr=0,Float_t zvertex=0.0);
  void LoadLocalSegments();
  void FindClusters();
  void AddClusters();
  void WriteClusters(Bool_t global=kTRUE);
  void WriteTracks(Int_t minHits);
  void SetNmaxOverlaps(Int_t i) {fNmaxOverlaps=i;}
  //void SetChiSqMax(Float_t f,Bool_t overlapping) {fChiSqMax[(Int_t)overlapping] = f;}
  void SetChiSqMax(Float_t f,Int_t lpatch) {fChiSqMax[lpatch] = f;}
  void SetInnerWidthFactor(Float_t y,Float_t z) {fYInnerWidthFactor=y; fZInnerWidthFactor=z;}
  void SetOuterWidthFactor(Float_t y,Float_t z) {fYOuterWidthFactor=y; fZOuterWidthFactor=z;}
  
  Float_t GetYWidthFactor() {return fCurrentPadRow < AliL3Transform::GetLastRow(1) ? fYInnerWidthFactor : fYOuterWidthFactor;}
  Float_t GetZWidthFactor() {return fCurrentPadRow < AliL3Transform::GetLastRow(1) ? fZInnerWidthFactor : fZOuterWidthFactor;}
  AliL3TrackArray *GetSeeds() {return fSeeds;}
  
 private:
  Int_t fNmaxOverlaps;
  Int_t fRowMin;
  Int_t fRowMax;
  Float_t fChiSqMax[3];
  Float_t fYInnerWidthFactor;
  Float_t fZInnerWidthFactor;
  Float_t fYOuterWidthFactor;
  Float_t fZOuterWidthFactor;
  Int_t fFitted;
  Int_t fFailed;
  static Int_t fgBadFitError;
  static Int_t fgFitError;
  static Int_t fgResultError;
  static Int_t fgFitRangeError;
  Bool_t fSeeding;
  Int_t fNMaxClusters;
  Int_t fNClusters;
  Int_t fEvent;
  AliL3TrackArray *fSeeds; //!
  AliL3TrackArray *fProcessTracks; //!
  AliL3SpacePointData *fClusters; //!
  
  void FitClusters(AliL3ModelTrack *track,Int_t *padrange,Int_t *timerange);
  Bool_t CheckCluster(Int_t trackindex);
  Bool_t IsMaximum(Int_t pad,Int_t time);
  Bool_t SetFitRange(AliL3ModelTrack *track,Int_t *padrange,Int_t *timerange);
  void SetClusterfitFalse(AliL3ModelTrack *track);
  void CalculateWeightedMean(AliL3ModelTrack *track,Int_t *padrange,Int_t *timerange);
  
  ClassDef(AliL3ClusterFitter,1) 

};

#endif
