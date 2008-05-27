/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$ 

/// \ingroup calib
/// \class AliMUONRegionalTriggerConfig
/// \brief The class defines the properties of regional trigger crate
///
/// \author Ch. Finck, Subatech Nantes; I. Hrivnacova, IPN Orsay

#ifndef ALIMUON_REGIONAL_TRIGGER_CONFIG_H
#define ALIMUON_REGIONAL_TRIGGER_CONFIG_H

#include <TObject.h>

#include "AliMpExMap.h"

#include <TObjArray.h>

class AliMUONTriggerCrateConfig;
class AliMUONLocalBoardConfig;

class AliMUONRegionalTriggerConfig : public  TObject{

  public:
    AliMUONRegionalTriggerConfig();
    AliMUONRegionalTriggerConfig(const AliMUONRegionalTriggerConfig& rhs);
    virtual ~AliMUONRegionalTriggerConfig();
    
    // operators
    AliMUONRegionalTriggerConfig& operator=(const AliMUONRegionalTriggerConfig& rhs);

    // methods
    Int_t ReadData(const TString& fileName = "");
    
    AliMUONTriggerCrateConfig* FindTriggerCrate(TString crateName, Bool_t warn = true) const;

    // method for looping
    
    Int_t GetNofTriggerCrates() const;
    AliMUONTriggerCrateConfig* GetTriggerCrate(Int_t index) const;
    AliMUONTriggerCrateConfig* GetTriggerCrateFast(Int_t index) const;
    TExMapIter GetTriggerCrateItr() const;

  
  private:
    // data members  
    AliMpExMap  fTriggerCrates; ///< map for trigger crates
 
  ClassDef(AliMUONRegionalTriggerConfig,1) // Regional trigger crate config
};

/// Return trigger crates iterator
inline TExMapIter AliMUONRegionalTriggerConfig::GetTriggerCrateItr() const { 
  return fTriggerCrates.GetIterator(); 
}
#endif 














