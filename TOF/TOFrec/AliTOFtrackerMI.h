#ifndef ALITOFTRACKERMI_H
#define ALITOFTRACKERMI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------//
//                                                                 //
//   AliTOFtrackerMI Class                                         //
//   Task: Perform association of the ESD tracks to TOF Clusters   //
//   and Update ESD track with associated TOF Cluster parameters   //
//                                                                 //
//-----------------------------------------------------------------//

#include "AliTracker.h"

class TTreeSRedirector;
class TClonesArray;
class TObjArray;

class AliESDEvent;
class AliESDpid;

class AliTOFcluster;
class AliTOFRecoParam;
class AliTOFGeometry;
class AliTOFtrack;

class AliTOFtrackerMI : public AliTracker {

public:

 AliTOFtrackerMI(); 

 //  virtual ~AliTOFtrackerMI() {delete fTOFpid;}
 virtual ~AliTOFtrackerMI();
 virtual void GetPidSettings(AliESDpid *esdPID);
 virtual Int_t Clusters2Tracks(AliESDEvent* /*event*/) {return -1;};
 virtual Int_t PropagateBack(AliESDEvent * const event);
 virtual Int_t RefitInward(AliESDEvent* /*event*/) {return -1;};
 virtual Int_t LoadClusters(TTree *dTree); // Loading Clusters from Digits
 virtual void  UnloadClusters();// UnLoad Clusters
 virtual AliCluster *GetCluster(Int_t index) const
   {if (index==-1 || index >= fN) return NULL;
   return (AliCluster *) fClusters[index];};
 void    GetLikelihood(Float_t dy, Float_t dz, const Double_t *cov, AliTOFtrack * track, Float_t & py, Float_t &pz) const;
 void FillClusterArray(TObjArray* arr) const;

private:

 enum {kMaxCluster=77777}; //maximal number of the TOF clusters

 AliTOFtrackerMI(const AliTOFtrackerMI &t); //Copy Ctor 
 AliTOFtrackerMI& operator=(const AliTOFtrackerMI &source); // ass. op.

 Int_t InsertCluster(AliTOFcluster *c); // Fills TofClusters Array
 Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
 void  MatchTracks(Bool_t mLastStep) const; // Matching Algorithm 
 void  MatchTracksMI(Bool_t mLastStep); // Matching Algorithm 
 void  CollectESD(); // Select starting Set for Matching 
 //void  Init();
 Float_t GetLinearDistances(AliTOFtrack * track, AliTOFcluster *cluster, Float_t distances[5]);

 const AliTOFRecoParam* fkRecoParam;    // Pointer to TOF Recon. Pars
 AliTOFGeometry*  fGeom;                // Pointer to TOF geometry
 AliTOFcluster *fClusters[kMaxCluster]; // pointers to the TOF clusters

 Int_t fN;              // Number of Clusters
 Int_t fNseeds;         // Number of track seeds  
 Int_t fNseedsTOF;      // TPC BP tracks
 Int_t fngoodmatch;     // Correctly matched  tracks
 Int_t fnbadmatch;      // Wrongly matched tracks
 Int_t fnunmatch;       // Unmatched tracks
 Int_t fnmatch;         // Total matched tracks
 
 Float_t fR;            // Intermediate radius in TOF, used in matching
 Float_t fTOFHeigth;    // Inner TOF radius for propagation
 Float_t fdCut;         // Cut on minimum distance track-pad in matching 
 Float_t fDx;           // Pad Size in X   
 Float_t fDy;           // Pad Size in Y (== X  TOF convention)
 Float_t fDz;           // Pad Size in Z 
 TClonesArray* fTracks; //! pointer to the TClonesArray with TOF tracks
 TObjArray* fSeeds;  //! pointer to the TObjArray with ESD tracks
 TTreeSRedirector *fDebugStreamer;     //!debug streamer
 ClassDef(AliTOFtrackerMI, 3) // TOF trackerMI 
};

#endif
