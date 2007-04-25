#ifndef ALIMUONDIGITMAKER_H
#define ALIMUONDIGITMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup base
/// \class AliMUONDigitMaker
/// \brief Reading Raw data class for trigger and tracker chambers
///
//  Author: Ch, Finck

#include <TObject.h>
#include "TStopwatch.h"

class TList;
class TArrayS;

class AliRawReader;
class AliMUONData;
class AliMUONDigit;
class AliMUONGlobalTrigger;
class AliMUONLocalTrigger;
class AliMUONTriggerCrateStore;
class AliMUONLocalStruct;

class AliMUONRawStreamTracker;
class AliMUONRawStreamTrigger;

class AliMUONDigitMaker : public TObject 
{
 public:
  AliMUONDigitMaker(Bool_t digit = kTRUE); // Constructor
  virtual ~AliMUONDigitMaker(void); // Destructor
    
  // write raw data
  Int_t  Raw2Digits(AliRawReader* rawReader);

  Int_t  ReadTrackerDDL(AliRawReader* rawReader);
  Int_t  ReadTriggerDDL(AliRawReader* rawReader);
 
         /// Return MUON data
  AliMUONData*   GetMUONData() const {return fMUONData;}
        /// Set MUON data
  void SetMUONData(AliMUONData* data) {fMUONData = data;}
 
  Int_t GetMapping(Int_t buspatchId, UShort_t manuId, 
			  UChar_t channelId, AliMUONDigit* digit );

  Int_t TriggerDigits(Int_t nBoard, TArrayS* xyPattern, TList& digitList );

        /// Set flag to generates scaler event
  void  SetScalerEvent() {fScalerEvent = kTRUE;}

      /// Disable trigger rawdata reading
  void  DisableTrigger() {fTriggerFlag = kFALSE;}

        /// Set Crate array
  void  SetCrateManager(AliMUONTriggerCrateStore* crateManager) {fCrateManager =  crateManager;}

       /// enable only list of digits for the display
  void SetDisplayFlag() { fDisplayFlag = kTRUE;  fDigitFlag = kFALSE;}


 private:
  /// Not implemented
  AliMUONDigitMaker (const AliMUONDigitMaker& rhs); // copy constructor
  /// Not implemented
  AliMUONDigitMaker& operator=(const AliMUONDigitMaker& rhs); // assignment operator

  void GetCrateName(Char_t* name, Int_t iDDL, Int_t iReg) const;

  AliMUONData*     fMUONData;          //!< Data container for MUON subsystem 
  
  Bool_t           fScalerEvent;       //!< flag to generates scaler event

  Bool_t           fDigitFlag;         //!< true for Digit, false for SDigit

  Bool_t           fTriggerFlag;       //!< true for reading also trigger rawdata

  Bool_t           fDisplayFlag;       //!< true for returning digits list to the display

  AliMUONRawStreamTracker* fRawStreamTracker;  //!< pointer of raw stream for tracker
  AliMUONRawStreamTrigger* fRawStreamTrigger;  //!< pointer of raw stream for trigger

  AliMUONDigit*        fDigit;                 //!< pointer to digits

  AliMUONLocalTrigger*  fLocalTrigger;         //!< pointer to local trigger
  AliMUONGlobalTrigger* fGlobalTrigger;        //!< pointer to local trigger

  AliMUONTriggerCrateStore* fCrateManager;     //!< Crate array

  TStopwatch fTrackerTimer;                    //!< time watcher for tracker part
  TStopwatch fTriggerTimer;                    //!< time watcher for trigger part
  TStopwatch fMappingTimer;                    //!< time watcher for mapping-tracker part

  ClassDef(AliMUONDigitMaker,1) // MUON digit maker from rawdata
};
	
#endif
