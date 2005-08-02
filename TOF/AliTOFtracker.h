#ifndef ALITOFTRACKER_H
#define ALITOFTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// AliTOFtracker Class
// Task: Perform association of the ESD tracks to TOF Clusters
// and Update ESD track with associated TOF Cluster parameters 
//
// -- Authors : S. Arcelli, C. Zampolli (Bologna University and INFN) 
// -- Contacts: Annalisa.De.Caro@cern.ch
// --         : Chiara.Zampolli@bo.infn.it
// --         : Silvia.Arcelli@bo.infn.it
//--------------------------------------------------------------------

/* $Id$ */

#include "TClonesArray.h"

#include "AliESD.h"
#include "AliTracker.h"

#include "AliTOFpidESD.h"

class AliTOFGeometry;
class AliTOFcluster;

class AliTOFtracker : public AliTracker {

enum {kMaxCluster=77777}; //maximal number of the TOF clusters

public:

 AliTOFtracker(AliTOFGeometry* geom, Double_t parPID[2]); 
 AliTOFtracker(const AliTOFtracker &t); //Copy Ctor 
  virtual ~AliTOFtracker() {delete fTOFpid;}
  virtual Int_t Clusters2Tracks(AliESD* /*event*/) {return -1;};
  virtual Int_t PropagateBack(AliESD* event);
  virtual Int_t RefitInward(AliESD* /*event*/) {return -1;};
  virtual Int_t LoadClusters(TTree *dTree); // Loading Clusters from Digits
  virtual void  UnloadClusters();// UnLoad Clusters
  virtual AliCluster *GetCluster(Int_t /*index*/) const {return NULL;};

private:

  Int_t InsertCluster(AliTOFcluster *c); // Fills TofClusters Array
  Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
  void  MatchTracks(Bool_t mLastStep); // Matching Algorithm 
  void  CollectESD(); // Select starting Set for Matching 
  void  Init();

  AliTOFGeometry*  fGeom;                 // Pointer to TOF geometry
  AliTOFpidESD*    fTOFpid;               // Pointer to TOF PID
  AliTOFcluster *fClusters[kMaxCluster];  // pointers to the TOF clusters

  Bool_t fHoles;         // flag for Geometry Version(w/wo Holes) temporary!
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
  TClonesArray* fSeeds;  //! pointer to the TClonesArray with ESD tracks

  ClassDef(AliTOFtracker, 1) // TOF tracker 
};

#endif


