/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 

/// \ingroup calib
/// \class AliMUONTriggerCrateConfig
/// \brief The class defines the configuration of trigger crate
///
/// \author Ch. Finck, Subatech Nantes

#ifndef ALIMUON_TRIGGER_CRATE_CONFIG_H
#define ALIMUON_TRIGGER_CRATE_CONFIG_H

#include "AliMpArrayI.h"

#include <TNamed.h>
#include <TString.h>
#include "AliMpArrayI.h"

class AliMUONTriggerCrateConfig : public  TNamed {

  public:
  
    AliMUONTriggerCrateConfig();
    AliMUONTriggerCrateConfig(const Char_t* name, UShort_t Id, UShort_t mask, 
                              UShort_t mode, UShort_t coinc);
    virtual ~AliMUONTriggerCrateConfig();
    
      /// get methods
    UShort_t GetId()  const;
    UShort_t GetMask() const;
    UShort_t GetMode() const;
    UShort_t GetCoinc() const;
    Int_t  GetNofLocalBoards() const;
    Int_t  GetLocalBoardId(Int_t index) const;
    Bool_t HasLocalBoard(Int_t localBoardId) const;
    Bool_t AddLocalBoard(Int_t localBoardId);
    
  private:

    /// Not implemented
    AliMUONTriggerCrateConfig(const AliMUONTriggerCrateConfig& rhs);
    /// Not implemented
    AliMUONTriggerCrateConfig& operator=(const AliMUONTriggerCrateConfig& rhs);

    // data members
    UShort_t     fId;         ///< crate number
    UShort_t     fMask;       ///< regional mask
    UShort_t     fMode;       ///< mode operating for crate
    UShort_t     fCoinc;      ///< coincidence mode for crate
    AliMpArrayI  fLocalBoard; ///< local board connected to this crate

  ClassDef(AliMUONTriggerCrateConfig,1)  // The class collectiong electronics properties of DDL
};

// inline functions

/// Return  Id
inline UShort_t AliMUONTriggerCrateConfig::GetId() const
{  return fId; }

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














