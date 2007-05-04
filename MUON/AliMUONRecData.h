#ifndef ALIMUONRECDATA_H
#define ALIMUONRECDATA_H
//
/*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup rec
/// \class AliMUONRecData
/// \brief Class containing MUON data: hits, digits, rawclusters, globaltrigger, 
/// localtrigger, etc ...
///
//  Author: Gines Martinez, Subatech,  September 2003

#include "AliLoader.h"

#include "AliMUONData.h"

class AliMUONRawCluster;
class AliMUONTrack;
class AliMUONTriggerTrack;

class AliRunLoader;

class TClonesArray;
class TObjArray;
class TTree;


//__________________________________________________________________
/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliMUONData                                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliMUONRecData : public AliMUONData
{
  public:
  
  //  enum EChamberIteration { kAllChambers, kTrackingChambers, kTriggerChambers };
  
    AliMUONRecData();
    AliMUONRecData(AliLoader * loader, const char* name, const char* title);
    AliMUONRecData(const char* galiceFile);
    virtual ~AliMUONRecData();  
    
    virtual void   AddRawCluster(Int_t id, const AliMUONRawCluster& clust);
    virtual void   AddRecTrack(const AliMUONTrack& track);
    virtual void   AddRecTriggerTrack(const AliMUONTriggerTrack& triggertrack);

    TClonesArray*  RawClusters(Int_t DetectionPlane);
    
                    /// Return reconstructed tracks
    TClonesArray*  RecTracks() {return fRecTracks;} 
                   /// Return reconstructed trigger tracks
    TClonesArray*  RecTriggerTracks() {return fRecTriggerTracks;}

    void           GetRawClusters() const;
    void           GetTrigger() const;
    void           GetRecTracks() const;
    void           GetRecTriggerTracks() const;

    Bool_t        IsRawClusterBranchesInTree();
    Bool_t        IsTrackBranchesInTree();
    Bool_t        IsTriggerBranchesInTree();
    Bool_t        IsTriggerTrackBranchesInTree();
    
    virtual void   Fill(Option_t* opt=" ");
    virtual void   MakeBranch(Option_t *opt=" ");
    virtual void   SetDataContainer(Option_t *opt=" ");
    virtual void   SetTreeAddress(Option_t *opt=" ");
    
    virtual void   ResetRawClusters();
    virtual void   ResetRecTracks();
    virtual void   ResetRecTriggerTracks();
  
                   /// Return tree with raw clusters
    TTree*         TreeR() {return fLoader->TreeR(); }
                   /// Return tree with tracks
    TTree*         TreeT() {return fLoader->TreeT(); }

                   // Methods to dump data
    void DumpRecPoints(Int_t event2Check=0, Option_t* opt="full");
    void DumpTracks(Int_t event2Check=0, Option_t* opt="full");
    void DumpRecTrigger(Int_t event2Check=0, Int_t write = 0, Bool_t readFromRP = kTRUE);
    
  protected: 
    /// Not implemented
    AliMUONRecData(const AliMUONRecData& rhs);
    /// Not implemented
    AliMUONRecData& operator=(const AliMUONRecData& rhs);

    TObjArray*      fRawClusters; ///< One event in TreeR/Rawcluster and one branch per tracking detection plane
    TClonesArray*   fRecTracks; ///< pointer to array of reconstructed tracks
    TClonesArray*   fRecTriggerTracks; ///< pointer to array of reconstructed trigger tracks

    Int_t*          fNrawclusters;  //!< Number of Raw Clusters
    Int_t           fNrectracks;    //!< Number of reconstructed tracks
    Int_t           fNrectriggertracks; //!< Number of reconstructed tracks
    Int_t           fSplitLevel;   ///< Splitting of branches 0 no spitting (root files are smaller) 1 splitting (larger output files)

    mutable Int_t fCurrentEvent; ///< Current event we're dealing with
    
private:  
    void   FillOwn(Option_t* opt=" ");
    void   MakeOwnBranch(Option_t *opt=" ");
    void   SetOwnDataContainer(Option_t *opt=" ");
    void   SetOwnTreeAddress(Option_t *opt=" ");
    

  ClassDef(AliMUONRecData,3) // Data accessor for MUON module
      
};
// inline functions


/// Load raw clusters tree
inline void AliMUONRecData::GetRawClusters() const {
  if (fLoader && fLoader->TreeR())
    fLoader->TreeR()->GetEvent(0);
}

/// Load trigger tree
inline void AliMUONRecData::GetTrigger() const {
  if (fLoader && fLoader->TreeR())
    fLoader->TreeR()->GetEvent(0);
}

/// Load reconstructed tracks
inline void AliMUONRecData::GetRecTracks() const {
  if (fLoader && fLoader->TreeT())
    fLoader->TreeT()->GetBranch("MUONTrack")->GetEvent(0);
}

/// Load reconstructed trigger tracks
inline void AliMUONRecData::GetRecTriggerTracks() const {
  if (fLoader && fLoader->TreeT())
    fLoader->TreeT()->GetBranch("MUONTriggerTrack")->GetEvent(0);
}



#endif

