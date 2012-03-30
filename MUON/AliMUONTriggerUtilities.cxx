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


#include "AliMUONTriggerUtilities.h"

#include "TArrayS.h"

#include "AliLog.h"

#include "AliMUONCalibrationData.h"
//#include "AliMUONTriggerCrateStore.h"
//#include "AliMUONTriggerCrate.h"
//#include "AliMUONTriggerCrateConfig.h"
//#include "AliMUONVCalibParam.h"
//#include "AliMUONRegionalTriggerConfig.h"
//#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONVDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONTriggerElectronics.h"
#include "AliMUONDigitStoreV2R.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONTriggerStoreV1.h"

#include "AliMpDDLStore.h"
#include "AliMpPad.h"
//#include "AliMpLocalBoard.h"
#include "AliMpConstants.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerUtilities)
/// \endcond


//_____________________________________________________________________________
AliMUONTriggerUtilities::AliMUONTriggerUtilities(AliMUONCalibrationData* calibData):
TObject(),
fCalibrationData(calibData),
fTriggerStatusMap(2*AliMUONConstants::NTriggerCh()*AliMUONConstants::NTriggerCircuit()),
fMaskedDigitsStore(new AliMUONDigitStoreV2R())
{
  /// Ctor.
  Init();
}

//_____________________________________________________________________________
AliMUONTriggerUtilities::~AliMUONTriggerUtilities()
{
  /// Destructor. Note we're the owner of some pointers.
  delete fMaskedDigitsStore;
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerUtilities::Init()
{
  /// Build trigger status map from masks
  AliMUONTriggerElectronics trigElectronics(fCalibrationData);
  AliMUONDigitMaker digitMaker(kFALSE);
  AliMUONDigitStoreV2R digitStore;
  AliMUONTriggerStoreV1 triggerStore;
  
  TArrayS xyPatternAll[2]; 	 
  for(Int_t icath=0; icath<AliMpConstants::NofCathodes(); icath++){ 	 
    xyPatternAll[icath].Set(AliMpConstants::NofTriggerChambers()); 	 
    xyPatternAll[icath].Reset(0xFFFF);
  }
  
  // Create a store with all digits in trigger
  for ( Int_t iboard=1; iboard<=AliMpConstants::NofLocalBoards(); iboard++ ) {
    digitMaker.TriggerDigits(iboard, xyPatternAll, digitStore, kFALSE);
  }
  
  // Create trigger with electronics (it applies masks)
  trigElectronics.Digits2Trigger(digitStore, triggerStore);
  
  // Re-compute digits from triggerStore
  // Since the masks were applied in the response,
  // the new store do not contain masked channels
  AliMUONDigitStoreV2R digitStoreMasked;
  digitMaker.TriggerToDigitsStore(triggerStore, digitStoreMasked);
  
  // Loop on non-masked digit store
  // Search for digits in the masked one:
  // if digit is not found, it means it was masked
  TIter next(digitStore.CreateIterator());
  AliMUONVDigit* dig = 0x0;
  while ( ( dig = static_cast<AliMUONVDigit*>(next()) ) ) {
    Int_t cath = dig->Cathode();
    Int_t detElemId = dig->DetElemId();
    Int_t board = dig->ManuId();
    Int_t strip = dig->ManuChannel();
    AliMUONVDigit* currDigit = digitStoreMasked.FindObject(detElemId, board, strip, cath);
    Bool_t isMasked = ( currDigit ) ? kFALSE : kTRUE;
    if ( isMasked ) fMaskedDigitsStore->Add(*((AliMUONVDigit*)dig->Clone()), AliMUONVDigitStore::kDeny);
    else {
      Int_t ich = detElemId/100-11;
      const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, AliMp::GetCathodType(cath));
      AliMpPad pad = seg->PadByIndices(dig->PadX(), dig->PadY(), kTRUE);
      for (Int_t iloc=0; iloc<pad.GetNofLocations(); iloc++) {
        Int_t currBoard = pad.GetLocalBoardId(iloc);
        Int_t arrayIndex = GetArrayIndex(cath, ich, currBoard);
        fTriggerStatusMap[arrayIndex] |= ( 0x1 << strip );
      } // loop on locations (in bending plane we have to fill all copy boards)
    }
  }

//  AliMUONTriggerCrateStore crates;
//  crates.ReadFromFile(fCalibrationData);
//  
//  AliMUONRegionalTriggerConfig* regionalConfig = fCalibrationData->RegionalTriggerConfig();
//  if ( ! regionalConfig ) AliFatal("no valid regional trigger configuration in CDB\n");
//  
//  // Loop on crates
//  AliMUONTriggerCrate* cr = 0x0;
//  TIter next ( crates.CreateCrateIterator() );
//  while ( ( cr = static_cast<AliMUONTriggerCrate*>(next()) ) ) {
//    TObjArray *boards = cr->Boards();
//    
//    AliMUONTriggerCrateConfig* crateConfig = regionalConfig->FindTriggerCrate(cr->GetName());
//    
//    if ( ! crateConfig ) AliFatal(Form("Crate %s not present in configuration !!!\n", cr->GetName()));
//    
//    UShort_t regionalMask = crateConfig->GetMask();
//    
//    // Loop on boards
//    for (Int_t iboard = 1; iboard < boards->GetEntries(); iboard++ ) {      
//      Bool_t activeBoard = ( ( regionalMask >> ( iboard - 1) ) & 1 );
//      AliMUONLocalTriggerBoard* board = (AliMUONLocalTriggerBoard*)boards->At(iboard);
//      Int_t cardNumber = board->GetNumber();
//      if ( cardNumber <= 0 ) continue; // interface board are not interested
//      AliMUONVCalibParam* localBoardMask = fCalibrationData->LocalTriggerBoardMasks(cardNumber);
//      for ( Int_t icath = 0; icath < 2; ++icath ) {
//        for ( Int_t ich = 0; ich < 4; ++ich ) {
//          Int_t planeIndex = icath * 4 + ich;
//          Int_t localMask = ( activeBoard ) ? localBoardMask->ValueAsInt(planeIndex) : 0;
//          Int_t arrayIndex = GetArrayIndex(icath, ich, cardNumber);
//          fTriggerStatusMap[arrayIndex] = localMask;
//        } // loop on chambers
//      } // loop on planes
//    } // loop on boards
//  } // loop on crates
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliMUONTriggerUtilities::IsMasked(const AliMUONVDigit& digit) const
{
  /// Check if pad is masked  
  return IsMasked(digit.DetElemId(), digit.Cathode(), digit.ManuId(), digit.ManuChannel());
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerUtilities::IsMasked(const AliMpPad& pad, Int_t detElemId, Int_t cathode) const
{
  /// Check if pad is masked  
  return IsMasked(detElemId, cathode, pad.GetLocalBoardId(0), pad.GetLocalBoardChannel(0));
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerUtilities::IsMasked(Int_t detElemId, Int_t cathode, Int_t localCircuit, Int_t strip) const
{
  /// Check if pad is masked
  Int_t trigCh = detElemId/100 - 11;
  
//  Int_t ibitxy = strip;
//  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(localCircuit);
//  if ( cathode && localBoard->GetSwitch(AliMpLocalBoard::kZeroAllYLSB) ) ibitxy += 8;
  Int_t arrayIndex = GetArrayIndex(cathode, trigCh, localCircuit);
  Bool_t isMasked = ( ( ( fTriggerStatusMap[arrayIndex] >> strip ) & 0x1 ) == 0 );
  AliDebug(1,Form("detElemId %i  cath %i  board %i  strip %i  is active %i\n", detElemId, cathode, localCircuit, strip, ! isMasked));
  return isMasked;
}


//_____________________________________________________________________________
Int_t AliMUONTriggerUtilities::GetArrayIndex(Int_t cathode, Int_t trigCh, Int_t localCircuit) const
{
  /// Get index of array with trigger status map or efficiency
  return
  AliMUONConstants::NTriggerCircuit() * AliMUONConstants::NTriggerCh() * cathode +
  AliMUONConstants::NTriggerCircuit() * trigCh + localCircuit-1;
}
