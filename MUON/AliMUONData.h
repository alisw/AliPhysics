#ifndef ALIMUONDATA_H
#define ALIMUONDATA_H
//
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// AliMUONData
// Class containing MUON data: hits, digits, rawclusters, globaltrigger, localtrigger, etc ...
// Gines Martinez, Subatech,  September 2003
//

#include "AliLoader.h"

class TClonesArray;
class TNamed;
class TObjArray;
class TTree;


class AliMUONConstants;
class AliMUONRawCluster;
class AliMUONTrack;

//__________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliMUONData                                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliMUONData : public TNamed {
 public:
    AliMUONData();
    AliMUONData(AliLoader * loader, const char* name, const char* title);
    AliMUONData(const AliMUONData& rMUONData);
    virtual ~AliMUONData();  
    virtual void   AddDigit(Int_t id, Int_t* tracks, Int_t* charges,
			     Int_t* digits); 
    virtual void   AddHit(Int_t fIshunt, Int_t track, Int_t iChamber, 
			  Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			  Float_t tof, Float_t momentum, Float_t theta, 
			  Float_t phi, Float_t length, Float_t destep);
    virtual void   AddHit(Int_t fIshunt, Int_t track, Int_t iChamber, 
			  Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			  Float_t tof, Float_t momentum, Float_t theta, 
			  Float_t phi, Float_t length, Float_t destep, 
			  Float_t Xref,Float_t Yref,Float_t Zref);
    virtual void   AddGlobalTrigger(Int_t *singlePlus, Int_t *singleMinus,
				    Int_t *singleUndef, Int_t *pairUnlike, 
				    Int_t *pairLike);
    virtual void   AddLocalTrigger(Int_t* ltrigger);
    virtual void   AddRawCluster(Int_t id, const AliMUONRawCluster& clust);
    virtual void   AddRecTrack(const AliMUONTrack& track);

    TClonesArray*  Hits() {return fHits;}
    TClonesArray*  Digits(Int_t DetectionPlane) 
      {return ( (TClonesArray*) fDigits->At(DetectionPlane) );}
    TClonesArray*  LocalTrigger() {return fLocalTrigger;}
    TClonesArray*  GlobalTrigger() {return fGlobalTrigger;}
    TClonesArray*  RawClusters(Int_t DetectionPlane)
      {return ( (TClonesArray*) fRawClusters->At(DetectionPlane) );}
    TClonesArray*  RecTracks() {return fRecTracks;}

    void           GetTrack(Int_t it) {fLoader->TreeH()->GetEvent(it);}
    Int_t          GetNtracks()       {return (Int_t) fLoader->TreeH()->GetEntries();}
    void           GetCathode(Int_t ic) {fLoader->TreeD()->GetEvent(ic);}
    void           GetRawClusters() {fLoader->TreeR()->GetEvent(0);}
    void           GetTrigger() {fLoader->TreeR()->GetEvent(0);}

    Bool_t        IsRawClusterBranchesInTree();
    Bool_t        IsTriggerBranchesInTree();

    virtual AliLoader* GetLoader() {return fLoader;}
    virtual void       SetLoader(AliLoader * loader) {fLoader=loader;}    
    
    virtual void   Fill(Option_t* opt=" ");
    virtual void   MakeBranch(Option_t *opt=" ");
    virtual void   SetTreeAddress(Option_t *opt=" ");
    
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetTrigger();
    virtual void   ResetRawClusters();
    virtual void   ResetRecTracks();
  
    TTree*         TreeH() {return fLoader->TreeH(); }
    TTree*         TreeD() {return fLoader->TreeD(); }
    TTree*         TreeR() {return fLoader->TreeR(); }
    TTree*         TreeT() {return fLoader->TreeT(); }
    TTree*         TreeP() {return fLoader->TreeP(); }

 private:  
    //descendant classes should
    //use protected interface methods to access these folders

  
 protected: 
    AliLoader*  fLoader;

    // One event in treeH per primary track
    TClonesArray*   fHits;    
    // One event in treeD and one branch per detection plane
    TObjArray*      fDigits; 
    //One event in TreeR/Rawcluster and one branch per tracking detection plane
    TObjArray*      fRawClusters; 
    //! List of Global Trigger One event in TreeR/GlobalTriggerBranch
    TClonesArray*   fGlobalTrigger; 
    //! List of Local Trigger, One event in TreeR/LocalTriggerBranch
    TClonesArray*   fLocalTrigger;  
    // pointer to array of reconstructed tracks
    TClonesArray*   fRecTracks; 

    Int_t           fNhits; //!
    Int_t*          fNdigits;//!
    Int_t*          fNrawclusters;//!
    Int_t           fNglobaltrigger;//!
    Int_t           fNlocaltrigger;//!
    Int_t           fNrectracks; //!

    ClassDef(AliMUONData,1)
 };
#endif
