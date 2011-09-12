#ifndef ALITOFTRACKER_H
#define ALITOFTRACKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//----------------------------------------------------------------------//
//                                                                      //
// AliTOFtracker Class                                                  //
// Task: Perform association of the ESD tracks to TOF Clusters          //
// and Update ESD track with associated TOF Cluster parameters          //
//                                                                      //
// -- Authors : S. Arcelli, C. Zampolli (Bologna University and INFN)   //
// -- Contacts: Annalisa.De.Caro@cern.ch                                //
// --         : Chiara.Zampolli@bo.infn.it                              //
// --         : Silvia.Arcelli@bo.infn.it                               //
//                                                                      //

#include "AliTracker.h"

#include "TObject.h"


class TClonesArray;
class TObjArray;

class TH1F;
class TH2F;

class AliESDEvent;
class AliESDpid;

class AliTOFcluster;
class AliTOFRecoParam;
class AliTOFGeometry;

class AliTOFtrackPoint : public TObject {

 public:

  AliTOFtrackPoint() :
    fIndex(0),fDistance(0),fDistanceZ(0),
    fDistanceY(0),fDistanceX(0),fPropRadius(0),fLength(0) { };
  AliTOFtrackPoint(Int_t index,Float_t dist,Float_t distZ,
		   Float_t distY,Float_t distX,Float_t radius,Float_t length) :
    TObject(),
    fIndex(index),fDistance(dist),fDistanceZ(distZ),
    fDistanceY(distY),fDistanceX(distX),fPropRadius(radius),fLength(length) { };
  AliTOFtrackPoint(const AliTOFtrackPoint & source) :
    TObject(source),
    fIndex(source.fIndex),
    fDistance(source.fDistance),
    fDistanceZ(source.fDistanceZ),
    fDistanceY(source.fDistanceY),
    fDistanceX(source.fDistanceX),
    fPropRadius(source.fPropRadius),
    fLength(source.fLength) { };
  AliTOFtrackPoint & operator=(const AliTOFtrackPoint & source)
    { if (this == &source) return *this;
      TObject::operator=(source);
      fDistance=source.fDistance;fDistanceZ=source.fDistanceZ;fDistanceY=source.fDistanceY;fDistanceX=source.fDistanceX;
      fPropRadius=source.fPropRadius;fLength=source.fLength;
      return *this; };

  Int_t Index()       const {return fIndex;} // cluster index
  Float_t Distance()  const {return fDistance;} // distance
  Float_t DistanceZ() const {return fDistanceZ;} // distance, Z component
  Float_t DistanceY() const {return fDistanceY;} // distance, Y  component
  Float_t DistanceX() const {return fDistanceX;} // distance, X  component
  Float_t PropRadius() const {return fPropRadius;} // propagation radius at TOF
  Float_t Length() const {return fLength;} // reconstructed track length at TOF

 private:

  Int_t fIndex; // cluster index 
  Float_t fDistance; // track-cluster distance
  Float_t fDistanceZ; //  Z component of track-cluster distance
  Float_t fDistanceY; //  Y component of track-cluster distance
  Float_t fDistanceX; //  X component of track-cluster distance
  Float_t fPropRadius; // track propagation radius
  Float_t fLength; // receonstructed track length

  //ClassDef(AliTOFtrackPoint, 2) // TOF matchable cluster

}; 

class AliTOFtracker : public AliTracker {

 public:

 AliTOFtracker(); 

 virtual ~AliTOFtracker();
 virtual void GetPidSettings(AliESDpid *esdPID);
 virtual Int_t Clusters2Tracks(AliESDEvent* /*event*/) {return -1;};
 virtual Int_t PropagateBack(AliESDEvent * const event);
 virtual Int_t RefitInward(AliESDEvent* /*event*/) {return -1;};
 virtual Int_t LoadClusters(TTree * cTree); // Load Clusters
 virtual void  UnloadClusters();// UnLoad Clusters
 virtual AliCluster *GetCluster(Int_t index) const
   {if (index==-1 || index >= fN) return NULL;
   return (AliCluster *) fClusters[index];};
 Bool_t GetTrackPoint(Int_t index, AliTrackPoint& p) const;
 void InitCheckHists();
 void SaveCheckHists();
 void FillClusterArray(TObjArray* arr) const;

private:

 enum {kMaxCluster=77777}; //maximal number of the TOF clusters

 AliTOFtracker(const AliTOFtracker &t); //Copy Ctor 
 AliTOFtracker& operator=(const AliTOFtracker &source); // ass. op.

 Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
 void  MatchTracks(Bool_t mLastStep); // Matching Algorithm 
 void  CollectESD(); // Select starting Set for Matching 
 Float_t CorrectTimeWalk(Float_t dist,Float_t tof) const; // Time Walk correction

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
 
 TClonesArray* fTracks; //! pointer to the TClonesArray with TOF tracks
 TObjArray* fSeeds;  //! pointer to the TObjArray with ESD tracks
                     //Digits/Reco QA histos
 TObjArray* fTOFtrackPoints; //! pointer to TObjArray of matchable TOF
			     //track points

 TH2F * fHDigClusMap; //Digits QA, Cluster Map 
 TH1F * fHDigNClus;   //Digits QA, # of clusters on TOF/event
 TH1F * fHDigClusTime;//Digits QA, Cluster Time (ns)
 TH1F * fHDigClusToT; //Digits QA, Cluster ToT (ns)
 TH1F * fHRecNClus; //Reco QA, cluster occupancy in search window
 TH1F * fHRecDist;//Reco QA, track-TOF cluster closest distance (cm)
 TH2F * fHRecSigYVsP;//Reco QA, track error in Y at TOF inner surface (cm)
 TH2F * fHRecSigZVsP; //Reco QA, track error in Z at TOF inner surface (cm)
 TH2F * fHRecSigYVsPWin;//Reco QA, search window size in Y (cm)
 TH2F * fHRecSigZVsPWin;//Reco QA, search window size in X (cm)
 TTree * fCalTree; // Tree for on-the-fly offline Calibration
 // internal variables in tree for on-the-fly TOF Calibration

 Int_t   fIch; //TOF channel number
 Float_t fToT; // Time over Threshold, ns
 Float_t fTime; //TOF time, ps
 Float_t fExpTimePi; // exp time, Pions
 Float_t fExpTimeKa; // exp time, Kaons
 Float_t fExpTimePr; // exp time, Protons

 ClassDef(AliTOFtracker, 6) // TOF tracker 
};

#endif
