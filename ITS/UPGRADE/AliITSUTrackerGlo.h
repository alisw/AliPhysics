#ifndef ALIITSTRACKERU_H
#define ALIITSTRACKERU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#define _FILL_CONTROL_HISTOS_

//-------------------------------------------------------------------------
//                ITS upgrade tracker base class
//-------------------------------------------------------------------------
#include "AliTracker.h"
#include "AliESDEvent.h"
#include "AliITSUSeed.h"
#include "AliITSUTrackCond.h"
#include "AliITSUTrackHyp.h"
#include <TArrayI.h>

class AliITSUReconstructor;
class AliITSURecoDet;
class AliITSUClusterPix;
class AliESDtrack;
class AliITSURecoLayer;
class TTree;


//-------------------------------------------------------------------------
class AliITSUTrackerGlo : public AliTracker {

  public:
  enum {kClus2Tracks,kPropBack,kRefitInw,kNTrackingPhases};   // tracking phases
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
  //
  enum {kDummyLabel=-3141593};
  //
  AliITSUTrackerGlo(AliITSUReconstructor* rec);
  virtual ~AliITSUTrackerGlo();

  virtual Int_t          Clusters2Tracks(AliESDEvent *event);
  virtual Int_t          PropagateBack(AliESDEvent *event);
  virtual Int_t          RefitInward(AliESDEvent *event);
  virtual Int_t          LoadClusters(TTree * treeRP=0);
  virtual void           UnloadClusters();
  virtual AliCluster*    GetCluster(Int_t index) const;
  void                   FlagSplitClusters();
  Bool_t                 ContainsSplitCluster(const AliITSUSeed* seed, Int_t maxSize=99999);
  void                   PrintSeedClusters(const AliITSUSeed* seed, Option_t* option="");
  //
  Int_t                  GetCountPronlongationTrials() const {return fCountProlongationTrials;}
  Int_t                  GetCountITSin()               const {return fCountITSin;}
  Int_t                  GetCountITSout()              const {return fCountITSout;}
  Int_t                  GetCountITSrefit()            const {return fCountITSrefit;}

  //------------------------------------
  AliITSURecoDet*        GetITSInterface()       const {return fITS;}
  //
  //------------------------------------
  Bool_t                 NeedToProlong(AliESDtrack* estTr);
  void                   Init(AliITSUReconstructor* rec);
  void                   FindTrack(AliESDtrack* esdTr, Int_t esdID);
  void                   CreateDefaultTrackCond();
  AliITSUTrackHyp*       InitHypothesis(AliESDtrack *esdTr, Int_t esdID);
  Bool_t                 TransportToLayer(AliITSUSeed* seed, Int_t lFrom, Int_t lTo);
  Bool_t                 TransportToLayer(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo);
  Bool_t                 TransportToLayerX(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo, Double_t xStop);  
  Bool_t                 GoToExitFromLayer(AliITSUSeed* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kTRUE);
  Bool_t                 GoToExitFromLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kTRUE);
  Bool_t                 GoToEntranceToLayer(AliITSUSeed* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kFALSE);
  Bool_t                 GoToEntranceToLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr, Int_t dir, Bool_t check=kFALSE);
  Bool_t                 PropagateSeed(AliITSUSeed *seed, Double_t xToGo, Double_t mass, Double_t maxStep=1.0, Bool_t matCorr=kTRUE);
  Bool_t                 PropagateSeed(AliExternalTrackParam *seed, Double_t xToGo, Double_t mass, Double_t maxStep=1.0, Bool_t matCorr=kTRUE);
  Double_t               RefitTrack(AliExternalTrackParam* trc, Double_t r, Int_t stopCond=0);
  Int_t                  GetTrackingPhase()                 const {return fTrackPhase;}

  //
  void                   KillSeed(AliITSUSeed* seed, Bool_t branch=kFALSE);
  Bool_t                 NeedToKill(AliITSUSeed* seed, Int_t flag);
  Bool_t                 GetRoadWidth(AliITSUSeed* seed, int ilrA);
  Bool_t                 CheckBackwardMatching(AliITSUSeed* seed);
  Int_t                  CheckCluster(AliITSUSeed* seed, Int_t lr, Int_t clID);
  void                   AddProlongationHypothesis(AliITSUSeed* seed, Int_t lr);
  Bool_t                 AddSeedBranch(AliITSUSeed* seed);
  void                   ValidateAllowedBranches(Int_t accMax);
  void                   ValidateAllowedCandidates(Int_t ilr, Int_t accMax);
  void                   FlagSeedClusters(const AliITSUSeed* seed, Bool_t flg);
  //
  AliITSUSeed*           NewSeedFromPool(const AliITSUSeed* src=0);
  void                   ResetSeedsPool();
  void                   MarkSeedFree(AliITSUSeed* seed );

  AliITSUTrackHyp*       GetTrackHyp(Int_t id)               const  {return (AliITSUTrackHyp*)fHypStore.UncheckedAt(id);}
  void                   SetTrackHyp(AliITSUTrackHyp* hyp,Int_t id) {fHypStore.AddAtAndExpand(hyp,id);}
  void                   DeleteLastSeedFromPool()                   {fSeedsPool.RemoveLast();}
  void                   SaveReducedHypothesesTree(AliITSUTrackHyp* dest);
  void                   FinalizeHypotheses();
  Bool_t                 FinalizeHypothesis(AliITSUTrackHyp* hyp);
  void                   UpdateESDTrack(AliITSUTrackHyp* hyp,Int_t flag);
  void                   CookMCLabel(AliITSUTrackHyp* hyp);
  void                   SetTrackingPhase(Int_t p)        {fTrackPhase = p;}
 //
 protected:
  TObject*&              NextFreeSeed();
  //
 private:
  //
  AliITSUTrackerGlo(const AliITSUTrackerGlo&);
  AliITSUTrackerGlo &operator=(const AliITSUTrackerGlo &tr);
  //
 protected:
  AliITSUReconstructor*           fReconstructor;  // ITS global reconstructor 
  AliITSURecoDet*                 fITS;            // interface to ITS, borrowed from reconstructor
  AliESDtrack*                    fCurrESDtrack;   // current esd track in processing
  Int_t                           fCurrESDtrMClb;  // its eventual mc label
  Double_t                        fCurrMass;       // current track mass
  Double_t                        fTrImpData[kNTrImpData];  // data on track impact on the layer
  //
  Int_t                           fCountProlongationTrials;   // number of TPC seeds
  Int_t                           fCountITSin;     // number of successful ITSin 
  Int_t                           fCountITSout;    // number of successful ITSout
  Int_t                           fCountITSrefit;  // number of successful ITSrefit 
  //
  // the seeds management to be optimized
  TObjArray                       fHypStore;       // storage for tracks hypotheses
  Int_t                           fLayerMaxCandidates; //! size of tmp candidates array 
  AliITSUSeed**                   fLayerCandidates;//! array for branches of current track prolongation
  Int_t                           fNBranchesAdded; // number of branches created for current seed in prolongation
  Int_t                           fNCandidatesAdded; // number of candidates added for current seed in prolongation
  AliITSUTrackHyp*                fCurrHyp;        //! hypotheses container for current track
  AliITSUTrackHyp*                fWorkHyp;        //! temporary hypothesis for track finding
  TClonesArray                    fSeedsPool;      //! pool for seeds
  TArrayI                         fFreeSeedsID;    //! array of ID's of freed seeds
  TArrayI                         fESDIndex;       //! array of ID's of freed seeds
  Int_t                           fNFreeSeeds;     //! number of seeds freed in the pool
  Int_t                           fLastSeedID;     //! id of the pool seed on which is returned by the NextFreeSeed method
  Int_t                           fNLrActive;      //! number of active layers
  //
  TObjArray                       fDefTrackConds;  //! default tracking conditions
  AliITSUTrackCond*               fCurrTrackCond;  //! current tracking condition
  Int_t                           fCurrActLrID;    //! current active layer ID being processed (set only when needed, not guaranteed)
  AliITSURecoLayer*               fCurrLayer;      //! current layer being processed  (set only when needed, not guaranteed)
  Int_t                           fTrackPhase;     //! tracking phase
  //
  static const Double_t           fgkToler;        // tracking tolerance
  //
#ifdef  _FILL_CONTROL_HISTOS_
  // this code is only for special histos needed to extract some control parameters
  void BookControlHistos();
  TObjArray* fCHistoArr;
  enum {kHResY=0,kHResYP=10,kHResZ=20,kHResZP=30,kHChi2Cl=40,kHChi2Nrm=50,kHBestInBranch=60,kHBestInCand=70};
  enum {kHClShare=0,kHChiMatchCorr,kHChiMatchFake,kHChiMatchCorrMiss,kHChiMatchFakeMiss,
	kHChiITSSACorr,kHChiITSSAFake}; // custom histos 
  enum {kHistosPhase=100};
  //
#endif
  //
  ClassDef(AliITSUTrackerGlo,1)   //ITS upgrade tracker
    
};

//________________________________________
inline TObject *&AliITSUTrackerGlo::NextFreeSeed()
{
  // return next free slot where the seed can be created
  fLastSeedID = fNFreeSeeds ? fFreeSeedsID.GetArray()[--fNFreeSeeds] : fSeedsPool.GetEntriesFast();
  //  AliInfo(Form("%d",fLastSeedID));
  return fSeedsPool[ fLastSeedID ];
  //
}

//_________________________________________________________________________
inline AliITSUSeed* AliITSUTrackerGlo::NewSeedFromPool(const AliITSUSeed* src)
{
  // create new seed, optionally copying from the source
  AliITSUSeed* sd =  src ? new( NextFreeSeed() ) AliITSUSeed(*src) : new( NextFreeSeed() ) AliITSUSeed();
  sd->SetPoolID(fLastSeedID);
  return sd;
}

//_________________________________________________________________________
inline void AliITSUTrackerGlo::KillSeed(AliITSUSeed* seed, Bool_t branch)
{
  // flag seed as killed, if requested, kill recursively its parents whose sole child is the seed being killed
  seed->Kill();
  seed = (AliITSUSeed*)seed->GetParent();
  if (seed && !seed->DecChildren() && branch) KillSeed(seed,branch);
}

//_________________________________________________________________________
inline void AliITSUTrackerGlo::AddProlongationHypothesis(AliITSUSeed* seed, Int_t lr)
{
  // add new seed prolongation hypothesis 
  fCurrHyp->AddSeed(seed,lr);
}

#endif

