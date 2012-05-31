#ifndef ALIITSTRACKERV2_H
#define ALIITSTRACKERV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                          ITS tracker
//     reads AliITSRecPoint clusters and creates AliITStrackV2 tracks
//           Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
#include "AliTracker.h"
#include "AliITSRecoParam.h"
#include "AliITStrackV2.h"
#include "AliITSgeomTGeo.h"


class AliITSRecPoint;
class AliESDEvent;
class TTree;


//-------------------------------------------------------------------------
class AliITStrackerV2 : public AliTracker {
public:
  AliITStrackerV2();
  AliITStrackerV2(const Char_t *geom);
  ~AliITStrackerV2(){}
  AliCluster *GetCluster(Int_t index) const;
  AliITSRecPoint *GetCluster(Int_t l, Int_t c) const {
    return fgLayers[l].GetCluster(c);
  }
  Int_t GetNumberOfClustersLayer(Int_t n) const {
     return fgLayers[n].GetNumberOfClusters();
  }   
  Int_t LoadClusters(TTree *cf);
  void UnloadClusters();
  Int_t Clusters2Tracks(AliESDEvent *event);
  Int_t PropagateBack(AliESDEvent *event);
  Int_t RefitInward(AliESDEvent *event);
  Bool_t RefitAt(Double_t x, AliITStrackV2 *seed, 
                 const AliITStrackV2 *t, Bool_t extra=kFALSE);
  void SetupFirstPass(Int_t *flags, Double_t *cuts=0);
  void SetupSecondPass(Int_t *flags, Double_t *cuts=0);

  void SetLastLayerToTrackTo(Int_t l=0) {fLastLayerToTrackTo=l;} 
  void SetLayersNotToSkip(Int_t *l);

  void UseClusters(const AliKalmanTrack *t, Int_t from=0) const;

  class AliITSdetector {
  public:
    AliITSdetector():fR(0.),fPhi(0.){}
    AliITSdetector(Double_t r,Double_t phi):fR(r),fPhi(phi){}
    Double_t GetR()   const {return fR;}
    Double_t GetPhi() const {return fPhi;}
  private:
    Double_t fR;    // polar coordinates 
    Double_t fPhi;  // of this detector
  };

  class AliITSlayer {
  public:
    enum {kNsector=5, kMaxClusterPerSector=AliITSRecoParam::kMaxClusterPerLayer/kNsector};
    AliITSlayer();
    AliITSlayer(Double_t r, Double_t p, Double_t z, Int_t nl, Int_t nd);
   ~AliITSlayer();
    Int_t InsertCluster(AliITSRecPoint *c);
    void ResetClusters();
    Int_t SelectClusters(Float_t zmi, Float_t zma, Float_t ymi, Float_t yma);
    const AliITSRecPoint *GetNextCluster(Int_t &ci);
    void ResetRoad();
    Double_t GetRoad() const {return fRoad;}
    Double_t GetR() const {return fR;}
    AliITSRecPoint *GetCluster(Int_t i) const { return fClusters[i]; } 
    AliITSdetector &GetDetector(Int_t n) const { return fDetectors[n]; }
    Int_t FindDetectorIndex(Double_t phi, Double_t z) const;
    Double_t GetThickness(Double_t y, Double_t z, Double_t &x0) const;
    Int_t GetNladders() const {return fNladders;}
    Int_t GetNdetectors() const {return fNdetectors;}
    Int_t GetNumberOfClusters() const;
  protected:
    AliITSlayer(const AliITSlayer&);
    AliITSlayer &operator=(const AliITSlayer &tr);  
    Double_t fR;                // mean radius of this layer
    Double_t fPhiOffset;        // offset of the first detector in Phi
    Int_t fNladders;            // number of ladders
    Double_t fZOffset;          // offset of the first detector in Z
    Int_t fNdetectors;          // detectors/ladder
    AliITSdetector *fDetectors; // array of detectors

    AliITSRecPoint *fClusters[AliITSRecoParam::kMaxClusterPerLayer]; // pointers to clusters
    Int_t fN[kNsector];         // numbers of clusters sector by sector 

    Int_t fIndex[AliITSRecoParam::kMaxClusterPerLayer]; // indexes of selected clusters 
    Int_t fNsel;                       // number of selected clusters

    Double_t fRoad;     // road defined by the cluster density
    Int_t FindClusterIndex(Float_t z, Int_t s) const;
  };

protected:
  AliITStrackerV2(const AliITStrackerV2&);
  void CookLabel(AliKalmanTrack *t,Float_t wrong) const;
  Double_t GetEffectiveThickness(Double_t y, Double_t z) const;
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
  Int_t fI;                              // index of the current layer
  static AliITSlayer fgLayers[AliITSgeomTGeo::kNLayers];// ITS layers
  AliITStrackV2 fTracks[AliITSgeomTGeo::kNLayers];      // track estimations at the ITS layers
  AliITStrackV2 fBestTrack;              // "best" track 
  AliITStrackV2 fTrackToFollow;          // followed track
  Int_t fPass;                           // current pass through the data 
  Int_t fConstraint[2];                  // constraint flags

  Int_t fLayersNotToSkip[AliITSgeomTGeo::kNLayers];     // layer masks
  Int_t fLastLayerToTrackTo;             // the innermost layer to track to

private:
  AliITStrackerV2 &operator=(const AliITStrackerV2 &tr);  
  ClassDef(AliITStrackerV2,1)   //ITS tracker V2
};


inline void AliITStrackerV2::SetupFirstPass(Int_t *flags, Double_t *cuts) {
  // This function sets up flags and cuts for the first tracking pass   
  //
  //   flags[0] - vertex constaint flag                                
  //              negative means "skip the pass"                        
  //              0        means "no constraint"                        
  //              positive means "normal constraint"                    

   fConstraint[0]=flags[0];
   if (cuts==0) return;
}

inline void AliITStrackerV2::SetupSecondPass(Int_t *flags, Double_t *cuts) {
  // This function sets up flags and cuts for the second tracking pass   
  //
  //   flags[0] - vertex constaint flag                                
  //              negative means "skip the pass"                        
  //              0        means "no constraint"                        
  //              positive means "normal constraint"                    

   fConstraint[1]=flags[0];
   if (cuts==0) return;
}

inline void AliITStrackerV2::CookLabel(AliKalmanTrack *t,Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
   Int_t tpcLabel=t->GetLabel();
   if (tpcLabel<0) return;
   AliTracker::CookLabel(t,wrong);
   if (tpcLabel != t->GetLabel()) t->SetLabel(-tpcLabel); 
}

#endif
