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
/// \brief Class containing MUON data: hits, digits, rawclusters, globaltrigger, 
/// localtrigger, etc ...
///
//  Author: Gines Martinez, Subatech,  September 2003

#include <TNamed.h>

#include "AliLoader.h"

class AliMUONDigit;
class AliMUONLocalTrigger;
class AliMUONRegionalTrigger;
class AliMUONGlobalTrigger;

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

class AliMUONData : public TNamed 
{
  public:
  
  //  enum EChamberIteration { kAllChambers, kTrackingChambers, kTriggerChambers };
  
    AliMUONData();
    AliMUONData(AliLoader * loader, const char* name, const char* title);
    AliMUONData(const char* galiceFile, const char* folderName);
    virtual ~AliMUONData();  
    virtual void   AddSDigit(Int_t id, const AliMUONDigit& digit); // use copy constructor
    virtual void   AddDigit(Int_t id, const AliMUONDigit& digit); // use copy constructor

    
    virtual void   AddGlobalTrigger(const AliMUONGlobalTrigger& trigger); // use copy constructor

    virtual void   AddLocalTrigger(const AliMUONLocalTrigger& trigger); // use copy constructor
 
    virtual void   AddRegionalTrigger(const AliMUONRegionalTrigger& trigger); // use copy constructor


    TClonesArray*  SDigits(Int_t DetectionPlane) const;
    TClonesArray*  Digits(Int_t DetectionPlane) const;
    TClonesArray*  LocalTrigger() const;
    TClonesArray*  RegionalTrigger() const;
    TClonesArray*  GlobalTrigger() const;    

    void           GetSDigits() const;
    void           GetDigits() const;
    void           GetTriggerD() const;

                   /// Return split level
    Int_t          GetSplitLevel() const {return fSplitLevel;}

                   /// Set split level
    void           SetSplitLevel(Int_t SplitLevel) {fSplitLevel=SplitLevel;}

    Bool_t        IsDigitsBranchesInTree();
    Bool_t        IsTriggerBranchesInTreeD();

                       /// Get loader
    virtual AliLoader* GetLoader() const { return fLoader; }
                       /// Set loader
    virtual void       SetLoader(AliLoader * loader) {fLoader=loader;}    
    
    virtual void   Fill(Option_t* opt=" ");
    virtual void   MakeBranch(Option_t *opt=" ");
    virtual void   SetDataContainer(Option_t *opt=" ");
    virtual void   SetTreeAddress(Option_t *opt=" ");
    
    virtual void Print(Option_t* opt="") const;
    
    virtual void   ResetSDigits();
    virtual void   ResetDigits();
    virtual void   ResetTrigger();
  
                   /// Return tree with summable digits
    TTree*         TreeS() {return fLoader->TreeS(); }
                   /// Return tree with digits
    TTree*         TreeD() {return fLoader->TreeD(); }

                   // Methods to dump data
    void DumpSDigits(Int_t event2Check=0, Option_t* opt="tracks");
    void DumpDigits(Int_t event2Check=0, Option_t* opt="tracks");
    
  protected: 
    /// Not implemented
    AliMUONData(const AliMUONData& rhs);
    /// Not implemented
    AliMUONData& operator=(const AliMUONData& rhs);

    AliRunLoader*   fRunLoader; //!< Run loader pointer
    AliLoader*      fLoader;  //!< Detector Loader pointer

    TObjArray*      fSDigits; ///< One event in treeS and one branch per detection plane
    TObjArray*      fDigits;  ///< One event in treeD and one branch per detection plane
    TClonesArray*   fGlobalTrigger; ///< List of Global Trigger One event in TreeR/GlobalTriggerBranch
    TClonesArray*   fLocalTrigger;  ///< List of Local Trigger, One event in TreeR/LocalTriggerBranch
    TClonesArray*   fRegionalTrigger;  ///< List of Regional Trigger, One event in TreeR/LocalTriggerBranch

    Int_t*          fNSdigits;//!< Number of Digits
    Int_t*          fNdigits; //!< Number of Digits
    Int_t           fNglobaltrigger;//!< Number of Global trigger
    Int_t           fNlocaltrigger; //!< Number of Local trigger
    Int_t           fNregionaltrigger; //!< Number of regional trigger
    Int_t           fSplitLevel;   ///< Splitting of branches 0 no spitting (root files are smaller) 1 splitting (larger output files)

    mutable Int_t fCurrentEvent; ///< Current event we're dealing with
    
private:  

  ClassDef(AliMUONData,4) // Data accessor for MUON module
      
};
// inline functions


/// Load sdigits tree
inline void AliMUONData::GetSDigits() const {
  if (fLoader && fLoader->TreeS())
    fLoader->TreeS()->GetEvent(0);
}

/// Load trigger D tree
inline void AliMUONData::GetTriggerD() const {
  if (fLoader && fLoader->TreeD())
    fLoader->TreeD()->GetEvent(0);
}


#endif

