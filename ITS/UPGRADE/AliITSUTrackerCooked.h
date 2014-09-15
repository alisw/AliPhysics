#ifndef ALIITSUTRACKERCooked_H
#define ALIITSUTRACKERCooked_H

//-------------------------------------------------------------------------
//                   The stand-alone ITSU tracker
//    The pattern recongintion based on the "cooked covariance" approach
//-------------------------------------------------------------------------

#include "AliITSUTrackerGlo.h"

class TTree;
class TClonesArray;
class TObjArray;

class AliESDEvent;
class AliCluster;
class AliITSUClusterPix;
class AliITSUTrackCooked;
class AliITSUReconstructor;

//-------------------------------------------------------------------------
class AliITSUTrackerCooked : public AliITSUTrackerGlo {
public:
  enum {
     kNLayers=7, kMaxClusterPerLayer=15000, kMaxSelected=kMaxClusterPerLayer/10
  };
  AliITSUTrackerCooked(AliITSUReconstructor *rec);
  virtual ~AliITSUTrackerCooked();

  // These functions must be implemented 
  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);
  Int_t LoadClusters(TTree *ct);
  void UnloadClusters();
  Bool_t 
  RefitAt(Double_t x, AliITSUTrackCooked *seed, const AliITSUTrackCooked *t);

  AliCluster *GetCluster(Int_t index) const;

  void SetSAonly(Bool_t sa=kTRUE) {fSAonly=sa;}
  Bool_t GetSAonly() const {return fSAonly;}

  // internal helper classes
  class AliITSUlayer;

protected:
  AliITSUTrackerCooked(const AliITSUTrackerCooked&);
  // Other protected functions
  Int_t MakeSeeds();
  Bool_t AddCookedSeed(const Float_t r1[3], Int_t l1, Int_t i1,
                       const Float_t r2[3], Int_t l2, Int_t i2,
                       const AliCluster *c3,Int_t l3, Int_t i3);
  void  FollowProlongation();
  Int_t TakeNextProlongation();
  void ResetTrackToFollow(const AliITSUTrackCooked &t);
  void ResetBestTrack();

private:
  AliITSUTrackerCooked &operator=(const AliITSUTrackerCooked &tr);

  // Data members
  // Internal tracker arrays, layers, modules, etc
  static AliITSUlayer fgLayers[kNLayers];// Layers
    
  TObjArray *fSeeds; // Track seeds

  Int_t fI;                              // index of the current layer
  AliITSUTrackCooked *fBestTrack;        // "best" track 
  AliITSUTrackCooked *fTrackToFollow;    // followed track
  
  Bool_t fSAonly; // kTRUE if the standalone tracking only

  ClassDef(AliITSUTrackerCooked,2)   //ITSU stand-alone tracker
};



// The helper classes
class AliITSUTrackerCooked::AliITSUlayer {
  public:
    AliITSUlayer();
    ~AliITSUlayer();

    void InsertClusters(TClonesArray *clusters, Bool_t seedingLayer, Bool_t sa);
    void SetR(Double_t r) {fR=r;}
    void DeleteClusters();
    void ResetSelectedClusters() {fI=0;}
    void SelectClusters(Float_t phi, Float_t dy, Float_t z, Float_t dz);
    const AliCluster *GetNextCluster(Int_t &i); 
    void ResetTrack(const AliITSUTrackCooked &t);
    Int_t FindClusterIndex(Double_t z) const;
    Float_t GetR() const {return fR;}
    AliCluster *GetCluster(Int_t i) const { return fClusters[i]; } 
    Float_t GetXRef(Int_t i) const { return fXRef[i]; } 
    Float_t GetAlphaRef(Int_t i) const { return fAlphaRef[i]; } 
    Float_t GetClusterPhi(Int_t i) const { return fPhi[i]; } 
    Int_t GetNumberOfClusters() const {return fN;}
    const AliITSUTrackCooked *GetTrack() const {return fTrack;}

  protected:
    AliITSUlayer(const AliITSUlayer&);
    AliITSUlayer &operator=(const AliITSUlayer &tr);  
    Int_t InsertCluster(AliCluster *c);

    Float_t fR;                // mean radius of this layer

    AliCluster *fClusters[kMaxClusterPerLayer]; // All clusters
    Float_t fXRef[kMaxClusterPerLayer];     // x of the reference plane
    Float_t fAlphaRef[kMaxClusterPerLayer]; // alpha of the reference plane
    Float_t fPhi[kMaxClusterPerLayer]; // cluster phi 
    Int_t fN; // Total number of clusters 

    Int_t fIndex[kMaxSelected]; 
    Int_t fNsel;      // Number of selected clusters
    Int_t fI;         // Running index for the selected clusters  

    AliITSUTrackCooked *fTrack; // track estimation at this layer
  };


#endif
