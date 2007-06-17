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

class TArrayS;

class AliRawReader;
class AliMUONTriggerCrateStore;
class AliMUONLocalStruct;

class AliMUONRawStreamTracker;
class AliMUONRawStreamTrigger;

class AliMUONVDigitStore;
class AliMUONVTriggerStore;

class AliMUONDigitMaker : public TObject 
{
 public:
  AliMUONDigitMaker(); // Constructor
  virtual ~AliMUONDigitMaker(void); // Destructor
    
  // write raw data
  Int_t  Raw2Digits(AliRawReader* rawReader, 
                    AliMUONVDigitStore* digitContainer=0,
                    AliMUONVTriggerStore* triggerStore=0);

  Int_t  ReadTrackerDDL(AliRawReader* rawReader);
  Int_t  ReadTriggerDDL(AliRawReader* rawReader);
  
  Int_t TriggerDigits(Int_t nBoard, TArrayS* xyPattern, 
                      AliMUONVDigitStore& digitStore) const;

        /// Set flag to generates scaler event
  void  SetScalerEvent() { fScalerEvent = kTRUE; }

        /// Set Crate array
  void  SetCrateManager(AliMUONTriggerCrateStore* crateManager) { fCrateManager =  crateManager; }

        /// Set flag whether or not we should generate digits for the trigger
  void  SetMakeTriggerDigits(Bool_t flag = kFALSE) { fMakeTriggerDigits = flag; }

private:
    
  /// Not implemented
  AliMUONDigitMaker (const AliMUONDigitMaker& rhs); // copy constructor
  /// Not implemented
  AliMUONDigitMaker& operator=(const AliMUONDigitMaker& rhs); // assignment operator

  void GetCrateName(Char_t* name, Int_t iDDL, Int_t iReg) const;

private:

  Bool_t           fScalerEvent;       //!< flag to generates scaler event
  Bool_t fMakeTriggerDigits; //!< whether or not we should generate digits for the trigger
  
  AliMUONRawStreamTracker* fRawStreamTracker;  //!< pointer of raw stream for tracker
  AliMUONRawStreamTrigger* fRawStreamTrigger;  //!< pointer of raw stream for trigger

  AliMUONTriggerCrateStore* fCrateManager;     //!< Crate array

  TStopwatch fTrackerTimer;                    //!< time watcher for tracker part
  TStopwatch fTriggerTimer;                    //!< time watcher for trigger part
  TStopwatch fMappingTimer;                    //!< time watcher for mapping-tracker part

  AliMUONVDigitStore* fDigitStore; //!< not owner
  AliMUONVTriggerStore* fTriggerStore; //!< not owner
  
  ClassDef(AliMUONDigitMaker,4) // MUON digit maker from rawdata
};
	
#endif
