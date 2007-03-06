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
//----------------------------------------------------------------------//

#include "AliTracker.h"

#include "AliTOFpidESD.h"

class TClonesArray;

class AliESD;

class AliTOFcluster;
class AliTOFGeometry;
class TH1F;
class TH2F;

class AliTOFtracker : public AliTracker {

enum {kMaxCluster=77777}; //maximal number of the TOF clusters

public:

 AliTOFtracker(AliTOFGeometry* geom, Double_t parPID[2]); 
 AliTOFtracker(const AliTOFtracker &t); //Copy Ctor 
 AliTOFtracker& operator=(const AliTOFtracker &source); // ass. op.

 virtual ~AliTOFtracker();
  //  virtual ~AliTOFtracker() {delete fTOFpid; SaveCheckHists();}
  virtual Int_t Clusters2Tracks(AliESD* /*event*/) {return -1;};
  virtual Int_t PropagateBack(AliESD* event);
  virtual Int_t RefitInward(AliESD* /*event*/) {return -1;};
  virtual Int_t LoadClusters(TTree * cTree); // Load Clusters
  virtual void  UnloadClusters();// UnLoad Clusters
  virtual AliCluster *GetCluster(Int_t /*index*/) const {return NULL;};
  Bool_t GetTrackPoint(Int_t index, AliTrackPoint& p) const;
  void InitCheckHists();
  void SaveCheckHists();

private:

  Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
  void  MatchTracks(Bool_t mLastStep); // Matching Algorithm 
  void  CollectESD(); // Select starting Set for Matching 

  AliTOFGeometry*  fGeom;                 // Pointer to TOF geometry
  AliTOFpidESD*    fTOFpid;               // Pointer to TOF PID
  AliTOFcluster *fClusters[kMaxCluster];  // pointers to the TOF clusters

  Bool_t fHoles;         // flag for Geometry Version(w/wo Holes) 
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
  Float_t fDzMax;        // Upper limit in z for the size of the search window 
  Float_t fDyMax;        // Upper limit in y for the size of the search window 
  TClonesArray* fTracks; //! pointer to the TClonesArray with TOF tracks
  TClonesArray* fSeeds;  //! pointer to the TClonesArray with ESD tracks
  //Digits/Reco QA histos

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

  Int_t fIch; //TOF channel number
  Float_t fToT; // Time over Threshold, ns
  Float_t fTime; //TOF time, ps
  Float_t fExpTimePi; // exp time, Pions
  Float_t fExpTimeKa; // exp time, Kaons
  Float_t fExpTimePr; // exp time, Protons

  ClassDef(AliTOFtracker, 1) // TOF tracker 
};

#endif
