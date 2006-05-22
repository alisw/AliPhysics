#ifndef ALIMUONDATA_H
#define ALIMUONDATA_H
//
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONData
/// \brief MUON data
///
/// Class containing MUON data: hits, digits, rawclusters, globaltrigger, 
/// localtrigger, etc ...
///
/// Author: Gines Martinez, Subatech,  September 2003

#include <TNamed.h>

class TArrayI;

#include "AliLoader.h"

class TClonesArray;
class TNamed;
class TObjArray;
class TTree;
class TIterator;

class AliMUONConstants;
class AliMUONRawCluster;
class AliMUONTrack;
class AliMUONTriggerTrack;
class AliMUONDigit;
class AliMUONHit;
class AliMUONLocalTrigger;
class AliMUONGlobalTrigger;

//__________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliMUONData                                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliMUONData : public TNamed 
{
  public:
  
  //  enum EChamberIteration { kAllChambers, kTrackingChambers, kTriggerChambers };
  
    AliMUONData();
    AliMUONData(AliLoader * loader, const char* name, const char* title);
    virtual ~AliMUONData();  
    virtual void   AddDigit(Int_t id, Int_t* tracks, Int_t* charges,
			     Int_t* digits); 
    virtual void   AddSDigit(Int_t id, Int_t* tracks, Int_t* charges,
			     Int_t* digits); 
    virtual void   AddDigit(Int_t id, const AliMUONDigit& digit); // use copy constructor
    virtual void   AddSDigit(Int_t id, const AliMUONDigit& digit); // use copy constructor
    virtual void   AddHit(Int_t fIshunt, Int_t track, Int_t iChamber, 
			  Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			  Float_t tof, Float_t momentum, Float_t theta, 
			  Float_t phi, Float_t length, Float_t destep, 
			  Float_t Xref,Float_t Yref,Float_t Zref);
			  // TBR
    virtual void   AddHit2(Int_t fIshunt, Int_t track, Int_t detElemId, 
			  Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			  Float_t tof, Float_t momentum, Float_t theta, 
			  Float_t phi, Float_t length, Float_t destep, 
			  Float_t Xref,Float_t Yref,Float_t Zref);
    
    virtual void   AddGlobalTrigger(const AliMUONGlobalTrigger& trigger); // use copy constructor

    virtual void   AddLocalTrigger(const AliMUONLocalTrigger& trigger); // use copy constructor

    virtual void   AddRawCluster(Int_t id, const AliMUONRawCluster& clust);
    virtual void   AddRecTrack(const AliMUONTrack& track);
    virtual void   AddRecTriggerTrack(const AliMUONTriggerTrack& triggertrack);

    TClonesArray*  Hits() {return fHits;}
    TClonesArray*  Digits(Int_t DetectionPlane) const;
    TClonesArray*  SDigits(Int_t DetectionPlane) const;
    TClonesArray*  LocalTrigger() const;
    TClonesArray*  GlobalTrigger() const;    
    TClonesArray*  RawClusters(Int_t DetectionPlane);
    TClonesArray*  RecTracks() {return fRecTracks;}
    TClonesArray*  RecTriggerTracks() {return fRecTriggerTracks;}

    void           GetTrack(Int_t it) const  {fLoader->TreeH()->GetEvent(it);}
    Int_t          GetNtracks() const      {return (Int_t) fLoader->TreeH()->GetEntries();}
    void           GetDigits() const;
    void           GetSDigits() const {fLoader->TreeS()->GetEvent(0);}
    void           GetRawClusters() const {fLoader->TreeR()->GetEvent(0);}
    void           GetTrigger() const {fLoader->TreeR()->GetEvent(0);}
    void           GetTriggerD() const {fLoader->TreeD()->GetEvent(0);}
    Int_t          GetSplitLevel() const {return fSplitLevel;}
    void           GetRecTracks() const {fLoader->TreeT()->GetEvent(0);}
    void           GetRecTriggerTracks() const {fLoader->TreeT()->GetEvent(0);}

    Bool_t        IsRawClusterBranchesInTree();
    Bool_t        IsDigitsBranchesInTree();
    Bool_t        IsTriggerBranchesInTree();
    Bool_t        IsTriggerBranchesInTreeD();
    Bool_t        IsTrackBranchesInTree();
    Bool_t        IsTriggerTrackBranchesInTree();

    virtual AliLoader* GetLoader() const { return fLoader; }
    virtual void       SetLoader(AliLoader * loader) {fLoader=loader;}    
    
    virtual void   Fill(Option_t* opt=" ");
    virtual void   MakeBranch(Option_t *opt=" ");
    virtual void   SetTreeAddress(Option_t *opt=" ");
    
    void           SetSplitLevel(Int_t SplitLevel) {fSplitLevel=SplitLevel;}
    
    virtual void Print(Option_t* opt="") const;
    
    virtual void   ResetHits();
    virtual void   ResetDigits();
    virtual void   ResetSDigits();
    virtual void   ResetTrigger();
    virtual void   ResetRawClusters();
    virtual void   ResetRecTracks();
    virtual void   ResetRecTriggerTracks();
  
    TTree*         TreeH() {return fLoader->TreeH(); }
    TTree*         TreeD() {return fLoader->TreeD(); }
    TTree*         TreeS() {return fLoader->TreeS(); }
    TTree*         TreeR() {return fLoader->TreeR(); }
    TTree*         TreeT() {return fLoader->TreeT(); }
    TTree*         TreeP() {return fLoader->TreeP(); }

    //    TIterator* CreateDigitIterator(AliMUONData::EChamberIteration type);
    
  protected: 
    AliMUONData(const AliMUONData& rhs);
    AliMUONData& operator=(const AliMUONData& rhs);

    AliLoader*      fLoader;  //!< Detector Loader pointer
    TClonesArray*   fHits;    ///< One event in treeH per primary track
    TObjArray*      fDigits;  ///< One event in treeD and one branch per detection plane
    TObjArray*      fSDigits; ///< One event in treeS and one branch per detection plane
    TObjArray*      fRawClusters; ///< One event in TreeR/Rawcluster and one branch per tracking detection plane
    TClonesArray*   fGlobalTrigger; ///< List of Global Trigger One event in TreeR/GlobalTriggerBranch
    TClonesArray*   fLocalTrigger;  ///< List of Local Trigger, One event in TreeR/LocalTriggerBranch
    TClonesArray*   fRecTracks; ///< pointer to array of reconstructed tracks
    TClonesArray*   fRecTriggerTracks; ///< pointer to array of reconstructed trigger tracks

    Int_t           fNhits;   //!< Number of Hits
    Int_t*          fNdigits; //!< Number of Digits
    Int_t*          fNSdigits;//!< Number of Digits
    Int_t*          fNrawclusters;  //!< Number of Raw Clusters
    Int_t           fNglobaltrigger;//!< Number of Global trigger
    Int_t           fNlocaltrigger; //!< Number of Local trigger
    Int_t           fNrectracks;    //!< Number of reconstructed tracks
    Int_t           fNrectriggertracks; //!< Number of reconstructed tracks
    Int_t           fSplitLevel;   ///< Splitting of branches 0 no spitting (root files are smaller) 1 splitting (larger output files)

    mutable Int_t fCurrentEvent; ///< Current event we're dealing with
    
private:  

  ClassDef(AliMUONData,3) // Data accessor for MUON module
      
};


#endif

