#ifndef ALIMUON_H
#define ALIMUON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliMUONTriggerCircuit.h" // cp

class AliMUONChamber;
class AliMUONLocalTrigger;
class AliMUONGlobalTrigger;
class AliMUONTriggerCircuit;
class AliMUONTriggerDecision;
class AliSegmentation;
class AliMUONResponse;
class AliMUONHit;
class AliMUONPadHit;
class AliMUONRawCluster;
class AliMUONClusterFinderVS;
class AliMUONReconstHit;
class TVector;
class TObjArray;
class TFile;
class TTree;


class AliMUON : public  AliDetector {
 public:
    AliMUON();
    AliMUON(const char *name, const char *title);
    AliMUON(const AliMUON& rMUON);
    virtual       ~AliMUON();
    virtual void   AddHit(Int_t track , Int_t *vol, Float_t *hits);
    virtual void   AddPadHit(Int_t* clhits);
    virtual void   AddDigits(Int_t id, Int_t* tracks, Int_t* charges,
			     Int_t* digits);
    virtual void   AddRawCluster(Int_t id, const AliMUONRawCluster& clust);
    virtual void   BuildGeometry();
    void           AddGlobalTrigger(Int_t *singlePlus, Int_t *singleMinus,
				    Int_t *singleUndef, Int_t *pairUnlike, 
				    Int_t *pairLike);
    void           AddLocalTrigger(Int_t* ltrigger);
    Int_t          DistancetoPrimitive(Int_t px, Int_t py);
    virtual Int_t  IsVersion() const {return 0;}
    TClonesArray  *PadHits() {return fPadHits;}
    TClonesArray  *LocalTrigger() {return fLocalTrigger;}
    TClonesArray  *GlobalTrigger() {return fGlobalTrigger;}
    virtual void   MakeBranch(Option_t *opt=" ");
    void           SetTreeAddress();
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetTrigger();
    virtual void   ResetRawClusters();
    // Cluster Finding
    virtual void   FindClusters(Int_t event ,Int_t lastEntry);
    // Digitisation 
    virtual void   Digitise(Int_t nev,Int_t bgrEvent, Option_t *opt1=" ",
    			    Option_t *opt2=" ",Text_t *name=" ");
    virtual void   SortTracks(Int_t *tracks,Int_t *charges, Int_t ntr);
// Configuration Methods (per station id)
//
// Set Chamber Segmentation Parameters
// id refers to the station and isec to the cathode plane
// Set Z values for all chambers
    virtual void SetChambersZ(const Float_t *Z);
    virtual void SetChambersZToDefault(void);
    virtual void SetPadSize(Int_t id, Int_t isec, Float_t p1, Float_t p2);
// Set Signal Generation Parameters
    virtual void   SetSigmaIntegration(Int_t id, Float_t p1);
    virtual void   SetChargeSlope(Int_t id, Float_t p1);
    virtual void   SetChargeSpread(Int_t id, Float_t p1, Float_t p2);
    virtual void   SetMaxAdc(Int_t id, Int_t p1);
// Set Segmentation and Response Model
    virtual void   SetSegmentationModel(Int_t id, Int_t isec,
					AliSegmentation *segmentation);
    virtual void   SetResponseModel(Int_t id, AliMUONResponse *response);
    virtual void   SetNsec(Int_t id, Int_t nsec);
// Set Reconstruction Model
    virtual void   SetReconstructionModel(Int_t id, AliMUONClusterFinderVS *reconstruction);
// Set Stepping Parameters
    virtual void   SetMaxStepGas(Float_t p1);
    virtual void   SetMaxStepAlu(Float_t p1);
    virtual void   SetMaxDestepGas(Float_t p1);
    virtual void   SetMaxDestepAlu(Float_t p1);
    virtual void   SetMuonAcc(Bool_t acc=0, Float_t angmin=2, Float_t angmax=9);
// Response Simulation
    virtual void   MakePadHits(Float_t xhit,Float_t yhit, Float_t zhit,
			       Float_t eloss, Float_t tof, Int_t id);
// get Trigger answer
    void   Trigger(Int_t nev);
// Return reference to Chamber #id
    virtual AliMUONChamber& Chamber(Int_t id)
	{return *((AliMUONChamber *) (*fChambers)[id]);}
// Return reference to Circuit #id
    virtual AliMUONTriggerCircuit& TriggerCircuit(Int_t id)
      {return *((AliMUONTriggerCircuit *) (*fTriggerCircuits)[id]);}
// Retrieve pad hits for a given Hit
    virtual AliMUONPadHit* FirstPad(AliMUONHit *hit, TClonesArray *padHits);
    virtual AliMUONPadHit* NextPad(TClonesArray *padHits);
// Return pointers to digits
    TObjArray            *Dchambers() {return fDchambers;}
    Int_t                *Ndch() {return fNdch;}
    virtual TClonesArray *DigitsAddress(Int_t id)
	{return ((TClonesArray *) (*fDchambers)[id]);}
// Return pointers to reconstructed clusters
    TObjArray            *RawClusters() {return fRawClusters;}
    Int_t                *Nrawch() {return fNrawch;} 
    virtual TClonesArray *RawClustAddress(Int_t id) 
	{return ((TClonesArray *) (*fRawClusters)[id]);}

    AliMUONRawCluster    *RawCluster(Int_t ichamber, Int_t icathod,
				     Int_t icluster);
// Copy Operator
    AliMUON& operator = (const AliMUON& rhs);
    
	    
 protected:
    TObjArray            *fChambers;           // List of Tracking Chambers
    TObjArray            *fTriggerCircuits;   // List of Trigger Circuits
    Int_t                 fNPadHits;           // Number of pad hits
    TClonesArray         *fPadHits;            // List of pad hits
    TObjArray            *fDchambers;          // List of digits
    Int_t                *fNdch;               // Number of digits
    TObjArray            *fRawClusters;        // List of raw clusters
    Int_t                *fNrawch;             // Number of raw clusters
    Int_t                 fNLocalTrigger;      // Number of Local Trigger 
    TClonesArray         *fLocalTrigger;       // List of Local Trigger      
    Int_t                 fNGlobalTrigger;     // Number of Global Trigger
    TClonesArray         *fGlobalTrigger;      // List of Global Trigger  

//
    Bool_t   fAccCut;          //Transport acceptance cut
    Float_t  fAccMin;          //Minimum acceptance cut used during transport
    Float_t  fAccMax;          //Minimum acceptance cut used during transport
//  

//  Stepping Parameters
   Float_t fMaxStepGas;      // Maximum step size inside the chamber gas
   Float_t fMaxStepAlu;      // Maximum step size inside the chamber aluminum
   Float_t fMaxDestepGas;    // Maximum relative energy loss in gas
   Float_t fMaxDestepAlu;    // Maximum relative energy loss in aluminum



// Pad Iterator
   Int_t fMaxIterPad;        // Maximum pad index
   Int_t fCurIterPad;        // Current pad index
// Background eent for event mixing
   Text_t *fFileName;           // ! File with background hits
   TTree *fTrH1;                // Hits Tree for background event
   TClonesArray *fHits2;        // List of hits for one track only
   TClonesArray *fPadHits2;     // List of clusters for one track only

   ClassDef(AliMUON,1)  // MUON Detector base class
};
#endif















