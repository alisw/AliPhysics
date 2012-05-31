/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: $ 

/// \ingroup mptrigger
/// \class AliMpRegionalTrigger
/// \brief The class defines the properties of regional trigger crate
///
/// \author Ch. Finck, Subatech Nantes; I. Hrivnacova, IPN Orsay

#ifndef ALI_MP_REGIONAL_TRIGGER_H
#define ALI_MP_REGIONAL_TRIGGER_H

#include <TObject.h>

#include "AliMpExMap.h"

#include <TObjArray.h>

class AliMpTriggerCrate;
class AliMpLocalBoard;
class AliMpDataStreams;
class TIterator;

class AliMpRegionalTrigger : public  TObject{

  public:
    AliMpRegionalTrigger();
    AliMpRegionalTrigger(const AliMpRegionalTrigger& rhs);
    AliMpRegionalTrigger(TRootIOCtor* ioCtor);
    virtual ~AliMpRegionalTrigger();
    
    // operators
    AliMpRegionalTrigger& operator=(const AliMpRegionalTrigger& rhs);

    // methods
    Bool_t ReadData(const TString& fileName);
    Bool_t ReadData(const AliMpDataStreams& dataStreams);
    
    AliMpTriggerCrate* FindTriggerCrate(TString crateName, Bool_t warn = true) const;
    AliMpLocalBoard*   FindLocalBoard(Int_t localBoardId, Bool_t warn = true) const;

    // method for looping

    TIterator* CreateCrateIterator() const;
    
    TIterator* CreateLocalBoardIterator() const;
    
    Int_t LocalBoardId(Int_t index) const;
    
    Int_t GetNofTriggerCrates() const;

    Int_t GetNofLocalBoards() const;
    
    // ownership
    void SetTriggerCratesOwner(Bool_t owner);

  private:
    Bool_t ReadData(istream& in);

    // data members  
    AliMpExMap  fTriggerCrates; ///< map for trigger crates
    AliMpExMap  fLocalBoardMap; ///< map of local boards (owner of boards)
    TObjArray   fLocalBoardArray; ///< array of local boards (not owner of boards, the map is the owner)
 
  ClassDef(AliMpRegionalTrigger,2) // Regional trigger crate
};


#endif //ALI_MP_REGIONAL__CRATE_H














