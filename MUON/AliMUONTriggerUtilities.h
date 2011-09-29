/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup core
/// \class AliMUONTriggerUtilities
/// \brief Utilities for trigger (check if pad is masked)
///
//  Author Diego Stocco

#ifndef ALIMUONTRIGGERUTILITIES_H
#define ALIMUONTRIGGERUTILITIES_H

#include "TObject.h"
#include "TArrayI.h"

class AliMUONCalibrationData;
class AliMUONVDigit;
class AliMpPad;

class AliMUONTriggerUtilities : public TObject
{
public:
  AliMUONTriggerUtilities(AliMUONCalibrationData* calibData);
  ~AliMUONTriggerUtilities();
  
  Bool_t IsMasked(const AliMUONVDigit& digit) const;
  Bool_t IsMasked(const AliMpPad& pad, Int_t detElemId, Int_t cathode) const;

private:
  /// Not implemented
  AliMUONTriggerUtilities(const AliMUONTriggerUtilities& other);
  /// Not implemented
  AliMUONTriggerUtilities& operator=(const AliMUONTriggerUtilities& other);
  
  Bool_t Init();
  Int_t GetArrayIndex(Int_t cathode, Int_t trigCh, Int_t localCircuit) const;
  
  AliMUONCalibrationData* fCalibrationData; //!< pointer to access calib parameters
  TArrayI fTriggerStatusMap; //!< Trigger masks
  
  ClassDef(AliMUONTriggerUtilities,0) // MUON Trigger utilities
};

#endif
