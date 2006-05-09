#ifndef ALIMUONDIGITMAKER_H
#define ALIMUONDIGITMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONDigitMaker
/// \brief Raw data class for trigger and tracker chambers
///
/// Readding Raw data class for trigger and tracker chambers

#include <TObject.h>
#include "TStopwatch.h"


class AliMpBusPatch;
class AliMUONData;
class AliMUONDigit;
class AliMUONGlobalTrigger;
class AliMUONLocalTrigger;

class AliMpSegFactory;
class AliMUONRawStreamTracker;
class AliMUONRawStreamTrigger;

class AliMUONDigitMaker : public TObject 
{
 public:
  AliMUONDigitMaker(AliMUONData* data); // Constructor
  virtual ~AliMUONDigitMaker(void); // Destructor
    
  // write raw data
  Int_t  Raw2Digits(AliRawReader* rawReader);

  Int_t  ReadTrackerDDL(AliRawReader* rawReader);
  Int_t  ReadTriggerDDL(AliRawReader* rawReader);

  AliMUONData*   GetMUONData() const {return fMUONData;}
  void SetMUONData(AliMUONData* data) {fMUONData = data;}
 
  Int_t GetMapping(Int_t buspatchId, UShort_t manuId, 
			  UChar_t channelId, AliMUONDigit* digit );

  void  SetScalerEvent() {fScalerEvent = kTRUE;}

 protected:
  AliMUONDigitMaker();                  // Default constructor
  AliMUONDigitMaker (const AliMUONDigitMaker& rhs); // copy constructor
  AliMUONDigitMaker& operator=(const AliMUONDigitMaker& rhs); // assignment operator

 private:

  AliMUONData*     fMUONData;          //! Data container for MUON subsystem 
  
  AliMpSegFactory* fSegFactory;        //! Mapping segmentation factory

  AliMpBusPatch*   fBusPatchManager;   //! buspatch versus DE's & DDL

  Bool_t           fScalerEvent;       //! flag to generates scaler event

  AliMUONRawStreamTracker* fRawStreamTracker;  //!pointer of raw stream for tracker
  AliMUONRawStreamTrigger* fRawStreamTrigger;  //!pointer of raw stream for trigger

  AliMUONDigit*        fDigit;         //! pointer to digits

  AliMUONLocalTrigger*  fLocalTrigger;  //! pointer to local trigger
  AliMUONGlobalTrigger* fGlobalTrigger;  //! pointer to local trigger

  TStopwatch fTrackerTimer; //!
  TStopwatch fTriggerTimer; //!
  TStopwatch fMappingTimer; //!

  ClassDef(AliMUONDigitMaker,1) // MUON digit maker from rawdata
};
	
#endif
