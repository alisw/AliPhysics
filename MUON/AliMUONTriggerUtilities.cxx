/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliLog.h"

#include "AliMUONCalibrationData.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONTriggerCrate.h"
#include "AliMUONTriggerCrateConfig.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONVDigit.h"
#include "AliMUONConstants.h"

#include "AliMpDDLStore.h"
#include "AliMpPad.h"
#include "AliMpLocalBoard.h"

#include "AliMUONTriggerUtilities.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerUtilities)
/// \endcond


//_____________________________________________________________________________
AliMUONTriggerUtilities::AliMUONTriggerUtilities(AliMUONCalibrationData* calibData):
TObject(),
fCalibrationData(calibData),
fTriggerStatusMap(2*AliMUONConstants::NTriggerCh()*AliMUONConstants::NTriggerCircuit())
{
  /// Ctor.
  Init();
}

//_____________________________________________________________________________
AliMUONTriggerUtilities::~AliMUONTriggerUtilities()
{
  /// Destructor. Note we're the owner of some pointers.
  
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerUtilities::Init()
{
  /// Build trigger status map from masks
  AliMUONTriggerCrateStore crates;
  crates.ReadFromFile(fCalibrationData);
  
  AliMUONRegionalTriggerConfig* regionalConfig = fCalibrationData->RegionalTriggerConfig();
  if ( ! regionalConfig ) AliFatal("no valid regional trigger configuration in CDB\n");
  
  // Loop on crates
  AliMUONTriggerCrate* cr = 0x0;
  TIter next ( crates.CreateCrateIterator() );
  while ( ( cr = static_cast<AliMUONTriggerCrate*>(next()) ) ) {
    TObjArray *boards = cr->Boards();
    
    AliMUONTriggerCrateConfig* crateConfig = regionalConfig->FindTriggerCrate(cr->GetName());
    
    if ( ! crateConfig ) AliFatal(Form("Crate %s not present in configuration !!!\n", cr->GetName()));
    
    UShort_t regionalMask = crateConfig->GetMask();
    
    // Loop on boards
    for (Int_t iboard = 1; iboard < boards->GetEntries(); iboard++ ) {      
      Bool_t activeBoard = ( ( regionalMask >> ( iboard - 1) ) & 1 );
      AliMUONLocalTriggerBoard* board = (AliMUONLocalTriggerBoard*)boards->At(iboard);
      Int_t cardNumber = board->GetNumber();
      if ( cardNumber <= 0 ) continue; // interface board are not interested
      AliMUONVCalibParam* localBoardMask = fCalibrationData->LocalTriggerBoardMasks(cardNumber);
      for ( Int_t icath = 0; icath < 2; ++icath ) {
        for ( Int_t ich = 0; ich < 4; ++ich ) {
          Int_t planeIndex = icath * 4 + ich;
          Int_t localMask = ( activeBoard ) ? localBoardMask->ValueAsInt(planeIndex) : 0;
          Int_t arrayIndex = GetArrayIndex(icath, ich, cardNumber);
          fTriggerStatusMap[arrayIndex] = localMask;
        } // loop on chambers
      } // loop on planes
    } // loop on boards
  } // loop on crates
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMUONTriggerUtilities::IsMasked(const AliMUONVDigit& digit) const
{
  /// Check if pad is masked
  Int_t detElemId = digit.DetElemId();
  Int_t localCircuit = digit.ManuId();
  Int_t strip = digit.ManuChannel();
  Int_t cathode = digit.Cathode();
  Int_t trigCh = detElemId/100 - 11;
  
  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(localCircuit);
  Int_t ibitxy = strip;
  if (cathode && localBoard->GetSwitch(AliMpLocalBoard::kZeroAllYLSB)) ibitxy += 8;
  Int_t arrayIndex = GetArrayIndex(cathode, trigCh, localCircuit);
  AliDebug(1,Form("ch %i  cath %i  board %i  strip %i  mask %i\n", trigCh, cathode, localCircuit, strip, (fTriggerStatusMap[arrayIndex] >> ibitxy ) & 0x1));
  return ((( fTriggerStatusMap[arrayIndex] >> ibitxy ) & 0x1 ) == 0 );
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerUtilities::IsMasked(const AliMpPad& pad, Int_t detElemId, Int_t cathode) const
{
  /// Check if pad is masked
  Int_t localCircuit = pad.GetLocalBoardId(0);
  Int_t strip = pad.GetLocalBoardChannel(0);
  Int_t trigCh = detElemId/100 - 11;
  
  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(localCircuit);
  Int_t ibitxy = strip;
  if (cathode && localBoard->GetSwitch(AliMpLocalBoard::kZeroAllYLSB)) ibitxy += 8;
  Int_t arrayIndex = GetArrayIndex(cathode, trigCh, localCircuit);
  AliDebug(1,Form("ch %i  cath %i  board %i  strip %i  mask %i\n", trigCh, cathode, localCircuit, strip, (fTriggerStatusMap[arrayIndex] >> ibitxy ) & 0x1));
  return ((( fTriggerStatusMap[arrayIndex] >> ibitxy ) & 0x1 ) == 0 );
}


//_____________________________________________________________________________
Int_t AliMUONTriggerUtilities::GetArrayIndex(Int_t cathode, Int_t trigCh, Int_t localCircuit) const
{
  /// Get index of array with trigger status map or efficiency
  return
  AliMUONConstants::NTriggerCircuit() * AliMUONConstants::NTriggerCh() * cathode +
  AliMUONConstants::NTriggerCircuit() * trigCh + localCircuit-1;
}
