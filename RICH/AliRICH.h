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

class TFile;

class AliRICHHit;
class AliRICHSDigit;
class AliRICHRawCluster;
class AliRICHRecHit1D;
class AliRICHRecHit3D;
class AliRICHClusterFinder;
class AliRICHDetect;
class AliRICHChamber;
class AliRICHCerenkov;
class AliSegmentation;
class AliRICHResponse;
class AliRICHEllipse;
class AliRICHGeometry;
class AliRICHMerger;

class AliRICH : public  AliDetector {
 public:
    AliRICH();
    AliRICH(const char *name, const char *title);
    AliRICH(const AliRICH& RICH);
    virtual       ~AliRICH();
    virtual void   AddHit(Int_t track, Int_t *vol, Float_t *hits);
    virtual void   AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
    virtual void   AddSDigit(Int_t *clhits);
    virtual void   AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits);
    virtual void   AddRawCluster(Int_t id, const AliRICHRawCluster& cluster);
    virtual void   AddRecHit1D(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsx, Int_t* padsy);
    virtual void   AddRecHit3D(Int_t id, Float_t* rechit);


    virtual void   BuildGeometry();
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual Float_t AbsoCH4(Float_t x);
    virtual Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
    virtual void   StepManager();
    Int_t          DistancetoPrimitive(Int_t px, Int_t py);
    virtual Int_t  IsVersion() const =0;
//
    TClonesArray  *SDigits() {return fSDigits;}
    TClonesArray  *Cerenkovs() {return fCerenkovs;}
    virtual void   MakeBranch(Option_t *opt=" ", char *file=0);
    void           SetTreeAddress();
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetRawClusters();
    virtual void   ResetRecHits1D();
    virtual void   ResetRecHits3D();
    virtual void   FindClusters(Int_t nev,Int_t lastEntry);
    virtual void   Digitise(Int_t nev,Int_t flag,Option_t *opt=" ",Text_t *name=" ");
    virtual void   SDigits2Digits();
// 
// Configuration Methods (per station id)
//
// Set Chamber Segmentation Parameters
// id refers to the station and isec to the cathode plane   
// Set Segmentation and Response Model
    virtual void   SetGeometryModel(Int_t id, AliRICHGeometry *geometry);
    virtual void   SetSegmentationModel(Int_t id, AliSegmentation *segmentation);
    virtual void   SetResponseModel(Int_t id, AliRICHResponse *response);
// Set Reconstruction Model
    virtual void   SetReconstructionModel(Int_t id, AliRICHClusterFinder *reconstruction);
// Set source debugging level
    void SetDebugLevel(Int_t level) {fDebugLevel=level;}
// Set Merger
    virtual void   SetMerger(AliRICHMerger* thisMerger) {fMerger=thisMerger;}  
// Get source debugging level
    Int_t GetDebugLevel() {return fDebugLevel;}
// Response Simulation
    virtual Int_t   Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id, ResponseType res);
// Return reference to Chamber #id
    virtual AliRICHChamber& Chamber(Int_t id) {return *((AliRICHChamber *) (*fChambers)[id]);}
// Retrieve pad hits for a given Hit
    virtual AliRICHSDigit* FirstPad(AliRICHHit *hit, TClonesArray *clusters);
    virtual AliRICHSDigit* NextPad(TClonesArray *clusters);
// Return pointers to digits 
    TObjArray            *Dchambers() {return fDchambers;}
    Int_t                *Ndch() {return fNdch;}
    virtual TClonesArray *DigitsAddress(Int_t id) {return ((TClonesArray *) (*fDchambers)[id]);}
// Return pointers to rec. hits
    TObjArray            *RecHits1D() {return fRecHits1D;}
    Int_t                *Nrechits1D() {return fNrechits1D;}
    virtual TClonesArray *RecHitsAddress1D(Int_t id) {return ((TClonesArray *) (*fRecHits1D)[id]);}
    TObjArray            *RecHits3D() {return fRecHits3D;}
    Int_t                *Nrechits3D() {return fNrechits3D;}
    virtual TClonesArray *RecHitsAddress3D(Int_t id) {return ((TClonesArray *) (*fRecHits3D)[id]);}
    
// Return pointers to reconstructed clusters
    virtual TClonesArray *RawClustAddress(Int_t id) {return ((TClonesArray *) (*fRawClusters)[id]);}    
// Assignment operator
    AliRICH& operator=(const AliRICH& rhs);
    
    
 protected:
    TObjArray            *fChambers;           // List of Tracking Chambers
    Int_t                 fNSDigits;           // Number of clusters
    Int_t                 fNcerenkovs;         // Number of cerenkovs
    TClonesArray         *fSDigits;            // List of clusters
    TObjArray            *fDchambers;          // List of digits
    TClonesArray         *fCerenkovs;          // List of cerenkovs
    Int_t                 fNdch[kNCH];         // Number of digits
    TObjArray            *fRawClusters;        // List of raw clusters
    TObjArray            *fRecHits1D;          // List of rec. hits
    TObjArray            *fRecHits3D;          // List of rec. hits
    Int_t                 fNrawch[kNCH];       // Number of raw clusters
    Int_t                 fNrechits1D[kNCH];   // Number of rec hits 
    Int_t                 fNrechits3D[kNCH];   // Number of rec hits 
    Int_t                 fDebugLevel;         // Source debugging level

    Int_t fCkovNumber;                   // Number of Cerenkov photons
    Int_t fCkovQuarz;                    // Cerenkovs crossing quartz
    Int_t fCkovGap;                      // Cerenkovs crossing gap
    Int_t fCkovCsi;                      // Cerenkovs crossing csi
    Int_t fLostRfreo;                    // Cerenkovs reflected in freon
    Int_t fLostRquar;                    // Cerenkovs reflected in quartz
    Int_t fLostAfreo;                    // Cerenkovs absorbed in freon 
    Int_t fLostAquarz;                   // Cerenkovs absorbed in quartz
    Int_t fLostAmeta;                    // Cerenkovs absorbed in methane
    Int_t fLostCsi;                      // Cerenkovs below csi quantum efficiency 
    Int_t fLostWires;                    // Cerenkovs lost in wires
    Int_t fFreonProd;                    // Cerenkovs produced in freon
    Float_t fMipx;                       // x coord. of MIP
    Float_t fMipy;                       // y coord. of MIP
    Int_t fFeedbacks;                    // Number of feedback photons
    Int_t fLostFresnel;                  // Cerenkovs lost by Fresnel reflection


// Background eent for event mixing
    Text_t *fFileName;           // ! File with background hits
    AliRICHMerger *fMerger;   // ! pointer to merger
    
    ClassDef(AliRICH,1)  //Hits manager for set:RICH
};
#endif















