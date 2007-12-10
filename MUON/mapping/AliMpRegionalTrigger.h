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

class AliMpRegionalTrigger : public  TObject{

  public:
    AliMpRegionalTrigger();
    AliMpRegionalTrigger(const AliMpRegionalTrigger& rhs);
    AliMpRegionalTrigger(TRootIOCtor* /*ioCtor*/);
    virtual ~AliMpRegionalTrigger();
    
    // operators
    AliMpRegionalTrigger& operator=(const AliMpRegionalTrigger& rhs);

    // methods
    Bool_t ReadData(const TString& fileName = "");
    
    AliMpTriggerCrate* FindTriggerCrate(TString crateName, Bool_t warn = true) const;
    AliMpLocalBoard*   FindLocalBoard(Int_t localBoardId, Bool_t warn = true) const;

    // method for looping
    
    Int_t GetNofTriggerCrates() const;
    AliMpTriggerCrate* GetTriggerCrate(Int_t index) const;
    AliMpTriggerCrate* GetTriggerCrateFast(Int_t index) const;
    TExMapIter GetTriggerCrateItr() const;

    Int_t GetNofLocalBoards() const;
    AliMpLocalBoard* GetLocalBoard(Int_t index) const;
    AliMpLocalBoard* GetLocalBoardFast(Int_t index) const;
    TExMapIter GetLocalBoardItr() const;

  private:
    // data members  
    AliMpExMap  fTriggerCrates; ///< map for trigger crates
    AliMpExMap  fLocalBoards;   ///< map for local boards
 
  ClassDef(AliMpRegionalTrigger,1) // Regional trigger crate
};

/// Return trigger crates iterator
inline TExMapIter AliMpRegionalTrigger::GetTriggerCrateItr() const { 
  return fTriggerCrates.GetIterator(); 
}

/// Return trigger local board iterator
inline TExMapIter AliMpRegionalTrigger::GetLocalBoardItr() const { 
  return fLocalBoards.GetIterator(); 
}
 

#endif //ALI_MP_REGIONAL__CRATE_H














