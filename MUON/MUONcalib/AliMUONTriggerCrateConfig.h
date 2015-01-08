/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 

/// \ingroup calib
/// \class AliMUONTriggerCrateConfig
/// \brief The class defines the configuration of trigger crate
///
/// \author Ch. Finck, Subatech Nantes

#ifndef ALI_MUON_TRIGGER_CRATE_CONFIG_H
#define ALI_MUON_TRIGGER_CRATE_CONFIG_H

#include "AliMpArrayI.h"
#include "AliMpTriggerCrate.h"

#include <TObject.h>
#include <TString.h>
#include "AliMpArrayI.h"

class AliMUONTriggerCrateConfig : public  TObject {

  public:
    AliMUONTriggerCrateConfig(AliMpTriggerCrate* mpTriggerCrate);
    AliMUONTriggerCrateConfig(TRootIOCtor* ioCtor);
    virtual ~AliMUONTriggerCrateConfig();
    
      // set methods
    void SetMask(UShort_t mask);
    void SetMode(UShort_t mode);
    void SetCoinc(UShort_t coinc);

      // get methods
    const Char_t* GetName() const;
    UShort_t GetId()  const;
    UShort_t GetMask() const;
    UShort_t GetMode() const;
    UShort_t GetCoinc() const;
    Int_t  GetNofLocalBoards() const;
    Int_t  GetLocalBoardId(Int_t index) const;
    Bool_t HasLocalBoard(Int_t localBoardId) const;
    Bool_t AddLocalBoard(Int_t localBoardId);
    
    // Only for checking data memebres for backward compatibility
    // These methods should not be called from other code !!!
    Int_t  GetNofLocalBoardsOld() const;
    Int_t  GetLocalBoardIdOld(Int_t index) const;

  private:
    /// Not implemented
    AliMUONTriggerCrateConfig();
    /// Not implemented
    AliMUONTriggerCrateConfig(const AliMUONTriggerCrateConfig& rhs);
    /// Not implemented
    AliMUONTriggerCrateConfig& operator=(const AliMUONTriggerCrateConfig& rhs);

    // data members
    AliMpTriggerCrate* fMpCrate; ///< mapping crate
    UShort_t           fMask;    ///< regional mask
    UShort_t           fMode;    ///< mode operating for crate
    UShort_t           fCoinc;   ///< coincidence mode for crate
    
    // not used data members kept for backward compatibility
    UShort_t     fId;            ///< crate number 
    AliMpArrayI  fLocalBoard;    ///< local board connected to this crate
 
  ClassDef(AliMUONTriggerCrateConfig,2)  // The class collectiong electronics properties of DDL
};

// inline functions

/// Set regional mask
inline void AliMUONTriggerCrateConfig::SetMask(UShort_t mask)   
{ fMask = mask; }

/// Set mode operating for crate
inline void AliMUONTriggerCrateConfig::SetMode(UShort_t mode)   
{ fMode = mode; }

/// Set coincidence mode for crate
inline void AliMUONTriggerCrateConfig::SetCoinc(UShort_t coinc) 
{ fCoinc = coinc; }

/// Return  name
inline const Char_t* AliMUONTriggerCrateConfig::GetName() const
{  return fMpCrate->GetName(); }

/// Return  Id
inline UShort_t AliMUONTriggerCrateConfig::GetId() const
{  return fMpCrate->GetId(); }

/// Return mask
inline UShort_t AliMUONTriggerCrateConfig::GetMask() const
{  return fMask; }

/// Return Mode
inline UShort_t AliMUONTriggerCrateConfig::GetMode() const
{  return fMode; }

/// Return coinc
inline UShort_t AliMUONTriggerCrateConfig::GetCoinc() const
{  return fCoinc; }

#endif














