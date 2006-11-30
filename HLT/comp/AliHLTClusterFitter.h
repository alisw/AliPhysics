// @(#) $Id$

#ifndef AliHLT_ClusterFitter
#define AliHLT_ClusterFitter

#include "AliHLTModeller.h"
#include "AliHLTTransform.h"

class AliHLTModelTrack;
class AliHLTTrackArray;
class AliHLTSpacePointData;

class AliHLTClusterFitter : public AliHLTModeller {
  
 public:
  AliHLTClusterFitter();
  AliHLTClusterFitter(Char_t *path);
  virtual ~AliHLTClusterFitter();
  
  void Init(Int_t slice,Int_t patch,Int_t *rowrange,AliHLTTrackArray *tracks);
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
  
  Float_t GetYWidthFactor() const
    {return fCurrentPadRow < AliHLTTransform::GetLastRow(1) ? fYInnerWidthFactor : fYOuterWidthFactor;}
  Float_t GetZWidthFactor() const 
    {return fCurrentPadRow < AliHLTTransform::GetLastRow(1) ? fZInnerWidthFactor : fZOuterWidthFactor;}
  AliHLTTrackArray *GetSeeds() {return fSeeds;}
  
 private:
  Int_t fNmaxOverlaps; // Max number of overlaps
  Int_t fRowMin; // Minimal row number (?)
  Int_t fRowMax; // Maximal row number (?)
  Float_t fChiSqMax[3]; // Maximal chi2 (?)
  Float_t fYInnerWidthFactor; // Inner width factor in Y
  Float_t fZInnerWidthFactor; // Inner width factor in Z
  Float_t fYOuterWidthFactor; // Outer width factor in Y
  Float_t fZOuterWidthFactor; // Outer width factor in Z
  Int_t fFitted; // Code for fitted (?)
  Int_t fFailed; // Code for failed
  static Int_t fgBadFitError; // Bad fit error
  static Int_t fgFitError; // Fit Error
  static Int_t fgResultError; // Result error
  static Int_t fgFitRangeError; // Fit range error
  Bool_t fSeeding; // Seeding (?)
  Int_t fNMaxClusters; // Max number of clusters
  Int_t fNClusters; // umver of clusters
  Int_t fEvent; // Current event
  AliHLTTrackArray *fSeeds; //! Array of seed
  AliHLTTrackArray *fProcessTracks; //! Array of processed tracks
  AliHLTSpacePointData *fClusters; //! Array of clusters
  
  void FitClusters(AliHLTModelTrack *track,Int_t *padrange,Int_t *timerange);
  Bool_t CheckCluster(Int_t trackindex);
  Bool_t IsMaximum(Int_t pad,Int_t time);
  Bool_t SetFitRange(AliHLTModelTrack *track,Int_t *padrange,Int_t *timerange);
  void SetClusterfitFalse(AliHLTModelTrack *track);
  void CalculateWeightedMean(AliHLTModelTrack *track,Int_t *padrange,Int_t *timerange);
  
  ClassDef(AliHLTClusterFitter,1) 

};

typedef AliHLTClusterFitter AliL3ClusterFitter; // for backward compatibility

#endif
