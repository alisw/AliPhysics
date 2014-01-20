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
#include "AliESDTOFcluster.h"

class TClonesArray;
class TObjArray;

class AliESDEvent;
class AliESDpid;

class AliTOFRecoParam;
class AliTOFGeometry;

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
 virtual AliCluster *GetCluster(Int_t index) const
   {if (index==-1 || index >= fN) return NULL;
     return (AliCluster *) (&fClusters[index]);};
 Bool_t GetTrackPoint(Int_t index, AliTrackPoint& p) const;
 Int_t GetNumberOfMatchedTOFtracks() const {return fnmatch;}
 void FillClusterArray(TObjArray* arr) const;
 void Clusterize();

private:

 enum {kMaxCluster=77777}; //maximal number of the TOF clusters

 AliTOFtrackerV2(const AliTOFtrackerV2 &t); //Copy Ctor 
 AliTOFtrackerV2& operator=(const AliTOFtrackerV2 &source); // ass. op.

 Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
 void  MatchTracks(); // Matching Algorithm 
 void  CollectESD(); // Select starting Set for Matching 
 Float_t CorrectTimeWalk(Float_t dist,Float_t tof) const; // Time Walk correction

 const AliTOFRecoParam* fkRecoParam;    // Pointer to TOF Recon. Pars
 AliTOFGeometry*  fGeom;                // Pointer to TOF geometry
 
 Int_t fN;              // Number of Clusters
 Int_t fNseeds;         // Number of track seeds  
 Int_t fNseedsTOF;      // TPC BP tracks
 Int_t fnunmatch;       // Unmatched tracks
 Int_t fnmatch;         // Total matched tracks

 TClonesArray* fTracks; //! pointer to the TClonesArray with TOF tracks
 TObjArray* fSeeds;  //! pointer to the TObjArray with ESD tracks
 AliESDTOFcluster *fClusters; // pointers to the TOF clusters

 ClassDef(AliTOFtrackerV2, 1) // TOF tracker 
};

#endif
