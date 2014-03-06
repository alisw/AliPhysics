#ifndef ALITOFTRACKERV2_H
#define ALITOFTRACKERV2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:  $ */

//------------------------------------------------------------------//
//                                                                  //
// AliTOFtrackerV2 Class                                            //
// Task: Perform association of the ESD tracks to TOF Clusters      //
// and Update ESD track with associated TOF Cluster parameters      //
//                                                                  //
// -- Authors : A. De Caro (Centro Studi e Ricerche E.Fermi)        //
// -- Contacts: Annalisa.De.Caro@cern.ch                            //
//                                                                  //
//------------------------------------------------------------------//

#include "AliTracker.h"
#include "AliTOFcluster.h"
#include "AliESDTOFCluster.h"

class TClonesArray;
class TObjArray;

class AliESDEvent;
class AliESDpid;
class AliTOFRecoParam;

class AliTOFtrackerV2 : public AliTracker {

 public:

 AliTOFtrackerV2(); 

 virtual ~AliTOFtrackerV2();
 virtual void GetPidSettings(AliESDpid *esdPID);
 virtual Int_t Clusters2Tracks(AliESDEvent* /*event*/) {return -1;};
 virtual Int_t PropagateBack(AliESDEvent * const event);
 virtual Int_t RefitInward(AliESDEvent* /*event*/) {return -1;};
 virtual Int_t LoadClusters(TTree * cTree); // Load Clusters
 virtual void  UnloadClusters();// UnLoad Clusters
 Bool_t GetTrackPoint(Int_t index, AliTrackPoint& p) const;
 Int_t GetNumberOfMatchedTOFtracks() const {return fnmatch;}
 virtual AliCluster *GetCluster(Int_t index) const
 {//Int_t index = fWrittenInPos[indexOr]; // indexOr should refer to the clusters written in ESD
   if (index==-1 || index >= fN) return NULL;
   return (AliCluster *) &(fClusters[index]);};

private:
 void Clusterize();
 void MergeClusters(Int_t i,Int_t j);

 enum {kMaxCluster=77777}; //maximal number of the TOF clusters

 AliTOFtrackerV2(const AliTOFtrackerV2 &t); //Copy Ctor 
 AliTOFtrackerV2& operator=(const AliTOFtrackerV2 &source); // ass. op.

 Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
 void  MatchTracks(); // Matching Algorithm 
 void  CollectESD(); // Select starting Set for Matching 
 Float_t CorrectTimeWalk(Float_t dist,Float_t tof) const; // Time Walk correction

 const AliTOFRecoParam* fkRecoParam;    // Pointer to TOF Recon. Pars
 
 Int_t fN;              // Number of Clusters
 Int_t fNseeds;         // Number of track seeds  
 Int_t fNseedsTOF;      // TPC BP tracks
 Int_t fnunmatch;       // Unmatched tracks
 Int_t fnmatch;         // Total matched tracks

 TObjArray* fSeeds;  //! pointer to the TObjArray with ESD tracks
 AliTOFcluster    *fClusters[kMaxCluster];     //! pointers to the TOF cluster
 TClonesArray     *fClustersESD;  //! base line for ESD clusters
 TClonesArray     *fHitsESD;       //! filter list of TOF hits for ESD
 Int_t            fWrittenInPos[kMaxCluster]; //! the position where the cluster is already written
 
 AliESDEvent      *fEvent;    //! pointer to the event

 Int_t fNsteps;                         //! number of propagation steps
 Double_t *fTimesAr[AliPID::kSPECIESC]; //! array to compute expected times for each propagation step
 Float_t  *fTrackPos[4];                //! array to compute pos for each propation step

 ClassDef(AliTOFtrackerV2, 2) // TOF tracker 
};

#endif
