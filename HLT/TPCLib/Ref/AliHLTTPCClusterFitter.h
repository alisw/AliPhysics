// @(#) $Id$

#ifndef AliHLTTPC_ClusterFitter
#define AliHLTTPC_ClusterFitter

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCModeller.h"
#include "AliHLTTPCTransform.h"

class AliHLTTPCModelTrack;
class AliHLTTPCTrackArray;
class AliHLTTPCSpacePointData;

class AliHLTTPCClusterFitter : public AliHLTTPCModeller {
  
 private:
  Int_t fNmaxOverlaps;
  Int_t fRowMin;
  Int_t fRowMax;
  Float_t fChiSqMax[2];
  Float_t fYInnerWidthFactor;
  Float_t fZInnerWidthFactor;
  Float_t fYOuterWidthFactor;
  Float_t fZOuterWidthFactor;
  Int_t fFitted;
  Int_t fFailed;
  static Int_t fBadFitError;
  static Int_t fFitError;
  static Int_t fResultError;
  static Int_t fFitRangeError;
  Bool_t fSeeding;
  Int_t fNMaxClusters;
  Int_t fNClusters;
  Int_t fEvent;
  AliHLTTPCTrackArray *fSeeds; //!
  AliHLTTPCTrackArray *fProcessTracks; //!
  AliHLTTPCSpacePointData *fClusters; //!
  
  void FitClusters(AliHLTTPCModelTrack *track,Int_t *padrange,Int_t *timerange);
  Bool_t CheckCluster(Int_t trackindex);
  Bool_t IsMaximum(Int_t pad,Int_t time);
  Bool_t SetFitRange(AliHLTTPCModelTrack *track,Int_t *padrange,Int_t *timerange);
  void SetClusterfitFalse(AliHLTTPCModelTrack *track);
  
 public:
  AliHLTTPCClusterFitter();
  AliHLTTPCClusterFitter(Char_t *path);
  virtual ~AliHLTTPCClusterFitter();
  
  void Init(Int_t slice,Int_t patch,Int_t *rowrange,AliHLTTPCTrackArray *tracks);
  void Init(Int_t slice,Int_t patch);
  void LoadSeeds(Int_t *rowrange,Bool_t offline=kTRUE,Int_t eventnr=0);
  void LoadLocalSegments();
  void FindClusters();
  void AddClusters();
  void WriteClusters(Bool_t global=kTRUE);
  void WriteTracks(Int_t min_hits);
  void SetNmaxOverlaps(Int_t i) {fNmaxOverlaps=i;}
  void SetChiSqMax(Float_t f,Bool_t overlapping) {fChiSqMax[(Int_t)overlapping] = f;}
  void SetInnerWidthFactor(Float_t y,Float_t z) {fYInnerWidthFactor=y; fZInnerWidthFactor=z;}
  void SetOuterWidthFactor(Float_t y,Float_t z) {fYOuterWidthFactor=y; fZOuterWidthFactor=z;}
  
  Float_t GetYWidthFactor() {return fCurrentPadRow < AliHLTTPCTransform::GetLastRow(1) ? fYInnerWidthFactor : fYOuterWidthFactor;}
  Float_t GetZWidthFactor() {return fCurrentPadRow < AliHLTTPCTransform::GetLastRow(1) ? fZInnerWidthFactor : fZOuterWidthFactor;}
  AliHLTTPCTrackArray *GetSeeds() {return fSeeds;}
  
  
  ClassDef(AliHLTTPCClusterFitter,1) 

};

#endif
