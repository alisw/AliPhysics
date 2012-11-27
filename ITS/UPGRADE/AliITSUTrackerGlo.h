#ifndef ALIITSTRACKERU_H
#define ALIITSTRACKERU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                ITS upgrade tracker base class
//-------------------------------------------------------------------------

#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliITSUSeed.h"
#include "AliITSUTrackCond.h"

class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSUClusterPix;
class AliESDtrack;

class TTree;


//-------------------------------------------------------------------------
class AliITSUTrackerGlo : public AliTracker {

  public:
  enum { // info from track extrapolation to layer for cluster check
    kTrXIn ,kTrYIn ,kTrZIn ,kTrPhiIn , // entrance (outer) point on the layer from above 
    kTrXOut,kTrYOut,kTrZOut,kTrPhiOut, // exit (inner) point on the layer
    kTrPhi0, kTrDPhi, kTrZ0, kTrDZ,     // mean phi,dPhi, mean z, dZ (don't change this order)
    kNTrImpData};
  //
  enum {kMissingCluster=0  // no cluster found on this layer
	,kTransportFailed=1  // seed did not reach target layer
	,kRWCheckFailed =2  // failed to rotate the seed to frame of the layer impact point
  };
  enum {kStopSearchOnSensor,kClusterNotMatching,kClusterMatching}; // flags for track-to-cluster checks

  AliITSUTrackerGlo(AliITSUReconstructor* rec);
  virtual ~AliITSUTrackerGlo();

  virtual Int_t          Clusters2Tracks(AliESDEvent *event);
  virtual Int_t          PropagateBack(AliESDEvent *event);
  virtual Int_t          RefitInward(AliESDEvent *event);
  virtual Int_t          LoadClusters(TTree * treeRP=0);
  virtual void           UnloadClusters();
  virtual AliCluster*    GetCluster(Int_t index) const;
  //------------------------------------
  AliITSURecoDet*        GetITSInterface()       const {return fITS;}
  //
  //------------------------------------
  Bool_t                 NeedToProlong(AliESDtrack* estTr);
  void                   Init(AliITSUReconstructor* rec);
  void                   FindTrack(AliESDtrack* esdTr);
  Bool_t                 InitSeed(AliESDtrack *esdTr);
  Int_t                  GetNSeeds(Int_t lr)              const {return fSeedsLr[lr].GetEntriesFast();} //RS TOCHECK
  AliITSUSeed*           GetSeed(Int_t lr, Int_t sID)     const {return (AliITSUSeed*)fSeedsLr[lr].UncheckedAt(sID);} //RS TOCHECK
  Bool_t                 TransportToLayer(AliITSUSeed* seed, Int_t lFrom, Int_t lTo);
  Bool_t                 NeedToKill(AliITSUSeed* seed, Int_t flag);
  void                   KillSeed(Int_t ilr, Int_t id) {} // todo
  Bool_t                 GetRoadWidth(AliITSUSeed* seed, int ilrA);
  Int_t                  CheckCluster(AliITSUSeed* seed, Int_t lr, Int_t clID);
  void                   AddProlongationHypothesis(AliITSUSeed* seed, Int_t lr);
  //
  AliITSUSeed*           NewSeedFromPool(const AliITSUSeed* src=0);
  void                   DeleteLastSeedFromPool()               {fSeedsPool.RemoveLast();}
  void                   ResetSeedTree();  // RS TOCHECK
  //
 private:
  
  AliITSUTrackerGlo(const AliITSUTrackerGlo&);
  AliITSUTrackerGlo &operator=(const AliITSUTrackerGlo &tr);
  //
 protected:
  AliITSUReconstructor*           fReconstructor;  // ITS global reconstructor 
  AliITSURecoDet*                 fITS;            // interface to ITS
  AliESDtrack*                    fCurrESDtrack;   // current esd track in processing
  Double_t                        fCurrMass;       // current track mass
  Double_t                        fTrImpData[kNTrImpData];  // data on track impact on the layer
  //
  // the seeds management to be optimized
  TObjArray*                      fSeedsLr;        // seeds at each layer
  TClonesArray                    fSeedsPool;      // pool for seeds
  //
  AliITSUTrackCond                fTrCond;         // tmp, to be moved to recoparam
  //
  ClassDef(AliITSUTrackerGlo,1)   //ITS upgrade tracker
    
};

//_________________________________________________________________________
inline void AliITSUTrackerGlo::AddProlongationHypothesis(AliITSUSeed* seed, Int_t lr)
{
  // add new seed prolongation hypothesis 
  fSeedsLr[lr].AddLast(seed);
  printf("*** Adding: "); seed->Print();
}


#endif

