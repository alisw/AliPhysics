#ifndef ALIRICH_H
#define ALIRICH_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////
//  Manager and hits classes for set:RICH     //
////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliRICHConst.h"
#include "AliRICHChamber.h"
static const int kNCH=7;


class AliRICHHit;
class AliRICHPadHit;
class AliRICHRawCluster;
class AliRICHRecHit;
class AliRICHClusterFinder;
class AliRICHDetect;
class AliRICHChamber;
class AliRICHCerenkov;
class AliRICHSegmentation;
class AliRICHResponse;
class AliRICHEllipse;
class AliRICHGeometry;

class AliRICH : public  AliDetector {
 public:
    AliRICH();
    AliRICH(const char *name, const char *title);
    AliRICH(const AliRICH& RICH);
    virtual       ~AliRICH();
    virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits);
    virtual void   AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
    virtual void   AddPadHit(Int_t *clhits);
    virtual void   AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual void   AddRawCluster(Int_t id, const AliRICHRawCluster& cluster);
    virtual void   AddRecHit(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsx, Int_t* padsy);


    virtual void   BuildGeometry();
    virtual void   CreateGeometry() {}
    virtual void   CreateMaterials() {}
    virtual void   StepManager();
    Int_t          DistancetoPrimitive(Int_t px, Int_t py);
    virtual Int_t  IsVersion() const =0;
//
    TClonesArray  *PadHits() {return fPadHits;}
    TClonesArray  *Cerenkovs() {return fCerenkovs;}
    virtual void   MakeBranch(Option_t *opt=" ");
    void           SetTreeAddress();
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetRawClusters();
    virtual void   ResetRecHits();
    virtual void   FindClusters(Int_t nev,Int_t lastEntry);
    virtual void   Digitise(Int_t nev,Int_t flag,Option_t *opt=" ",Text_t *name=" ");
// 
// Configuration Methods (per station id)
//
// Set Chamber Segmentation Parameters
// id refers to the station and isec to the cathode plane   
// Set Segmentation and Response Model
    virtual void   SetGeometryModel(Int_t id, AliRICHGeometry *geometry);
    virtual void   SetSegmentationModel(Int_t id, AliRICHSegmentation *segmentation);
    virtual void   SetResponseModel(Int_t id, AliRICHResponse *response);
    virtual void   SetNsec(Int_t id, Int_t nsec);
// Set Reconstruction Model
    virtual void   SetReconstructionModel(Int_t id, AliRICHClusterFinder *reconstruction);
// Set source debugging level
    void SetDebugLevel(Int_t level) {fDebugLevel=level;}
// Get source debugging level
    Int_t GetDebugLevel() {return fDebugLevel;}
// Response Simulation
    virtual Int_t   MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id, ResponseType res);
// Return reference to Chamber #id
    virtual AliRICHChamber& Chamber(Int_t id) {return *((AliRICHChamber *) (*fChambers)[id]);}
// Retrieve pad hits for a given Hit
    virtual AliRICHPadHit* FirstPad(AliRICHHit *hit, TClonesArray *clusters);
    virtual AliRICHPadHit* NextPad(TClonesArray *clusters);
// Return pointers to digits 
    TObjArray            *Dchambers() {return fDchambers;}
    Int_t                *Ndch() {return fNdch;}
    virtual TClonesArray *DigitsAddress(Int_t id) {return ((TClonesArray *) (*fDchambers)[id]);}
// Return pointers to rec. hits
    TObjArray            *RecHits() {return fRecHits;}
    Int_t                *Nrechits() {return fNrechits;}
    virtual TClonesArray *RecHitsAddress(Int_t id) {return ((TClonesArray *) (*fRecHits)[id]);}
// Return pointers to reconstructed clusters
    virtual TClonesArray *RawClustAddress(Int_t id) {return ((TClonesArray *) (*fRawClusters)[id]);}    
// Assignment operator
    AliRICH& operator=(const AliRICH& rhs);
    
    
 protected:
    TObjArray            *fChambers;           // List of Tracking Chambers
    Int_t                 fNPadHits;           // Number of clusters
    Int_t                 fNcerenkovs;         // Number of cerenkovs
    TClonesArray         *fPadHits;            // List of clusters
    TObjArray            *fDchambers;          // List of digits
    TClonesArray         *fCerenkovs;          // List of cerenkovs
    Int_t                *fNdch;               // Number of digits
    Text_t               *fFileName;           // Filename for event mixing
    TObjArray            *fRawClusters;        // List of raw clusters
    TObjArray            *fRecHits;            // List of rec. hits
    Int_t                *fNrawch;             // Number of raw clusters
    Int_t                *fNrechits;           // Number of rec hits 
    Int_t                 fDebugLevel;          // Source debugging level

 protected:
    
    ClassDef(AliRICH,1)  //Hits manager for set:RICH
};
#endif















