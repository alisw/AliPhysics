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

class TArrayS;

class AliRawReader;
class AliMUONLocalStruct;

class AliMUONVRawStreamTracker;
class AliMUONVRawStreamTrigger;

class AliMUONVDigitStore;
class AliMUONVTriggerStore;

class AliMUONLogger;

class AliMUONDigitMaker : public TObject 
{
 public:
  AliMUONDigitMaker(
         Bool_t enableErrorLogger = kTRUE,
         Bool_t useFastTrackerDecoder = kTRUE, Bool_t useFastTriggerDecoder = kTRUE
     ); // Constructor
  virtual ~AliMUONDigitMaker(void); // Destructor
    
  /// Code to indicate readout errors
  enum ErrorCode
  {
    kOK=0,             ///< everything is OK 
    kTrackerBAD=1<<1,  ///< tracker part had readout errors 
    kTriggerBAD=1<<2   ///< trigger part had readout errors 
  };
  
  // write raw data
  Int_t  Raw2Digits(AliRawReader* rawReader, 
                    AliMUONVDigitStore* digitContainer=0,
                    AliMUONVTriggerStore* triggerStore=0);

  Int_t  ReadTrackerDDL(AliRawReader* rawReader);
  Int_t  ReadTriggerDDL(AliRawReader* rawReader);
  Int_t  ReadTriggerDDLFast(AliRawReader* rawReader);
  
  Int_t TriggerDigits(Int_t nBoard, TArrayS* xyPattern, 
                      AliMUONVDigitStore& digitStore) const;

  Bool_t TriggerToDigitsStore(const AliMUONVTriggerStore& triggerStore, 
			      AliMUONVDigitStore& digitStore) const;

        /// Set flag to generates scaler event
  void  SetScalerEvent() { fScalerEvent = kTRUE; }

        /// Set flag whether or not we should generate digits for the trigger
  void  SetMakeTriggerDigits(Bool_t flag = kFALSE) { fMakeTriggerDigits = flag; }

  Bool_t UsingFastTrackerDecoder() const;
  Bool_t UsingFastTriggerDecoder() const;
  void  SetFastTrackerDecoder(Bool_t useFastDecoder);
  void  SetFastTriggerDecoder(Bool_t useFastDecoder);

  /// Return the raw stream object which decodes DDL raw data from tracking stations.
  AliMUONVRawStreamTracker* GetRawStreamTracker() const { return fRawStreamTracker; }

  /// Return the raw stream object which decodes DDL raw data from the trigger system.
  AliMUONVRawStreamTrigger* GetRawStreamTrigger() const { return fRawStreamTrigger; }

  void Print(Option_t* opt="") const;

private:
    
  /// Not implemented
  AliMUONDigitMaker (const AliMUONDigitMaker& rhs); // copy constructor
  /// Not implemented
  AliMUONDigitMaker& operator=(const AliMUONDigitMaker& rhs); // assignment operator

private:
  void CreateRawStreamTracker(Bool_t useFastDecoder);
  void CreateRawStreamTrigger(Bool_t useFastDecoder);

  Bool_t fScalerEvent;       //!< flag to generates scaler event
  Bool_t fMakeTriggerDigits; //!< whether or not we should generate digits for the trigger
  
  AliMUONVRawStreamTracker* fRawStreamTracker; //!< pointer of raw stream for tracker
  AliMUONVRawStreamTrigger* fRawStreamTrigger;  //!< pointer of raw stream for trigger

  AliMUONVDigitStore* fDigitStore; //!< not owner
  AliMUONVTriggerStore* fTriggerStore; //!< not owner

  AliMUONLogger* fLogger; //!< to log messages
  
  ClassDef(AliMUONDigitMaker,6) // MUON digit maker from rawdata
};
	
#endif
