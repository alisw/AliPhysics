#ifndef ALIMUONRAWREADER_H
#define ALIMUONRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONRawReader
/// \brief Raw data class for trigger and tracker chambers
///
/// Readding Raw data class for trigger and tracker chambers

#include <TObject.h>
#include "TStopwatch.h"

class AliMpBusPatch;
class AliMUONData;
class AliMUONDigit;
class AliMUONDDLTracker;
class AliMUONDDLTrigger;
class AliMUONGlobalTrigger;
class AliRawReader;
class AliMpSegFactory;

class AliMUONRawReader : public TObject 
{
 public:
  AliMUONRawReader(AliMUONData* data); // Constructor
  virtual ~AliMUONRawReader(void); // Destructor
    
  // write raw data
  Int_t   Raw2Digits(AliRawReader* rawReader);

  Int_t ReadTrackerDDL(AliRawReader* rawReader);
  Int_t ReadTriggerDDL(AliRawReader* rawReader);

  AliMUONData*   GetMUONData() {return fMUONData;}

  Int_t GetMapping(Int_t buspatchId, UShort_t manuId, 
			  UChar_t channelId, AliMUONDigit* digit );

  AliMUONGlobalTrigger* GetGlobalTriggerPattern(Int_t gloTrg) const;

  void  SetScalerEvent() {fScalerEvent = kTRUE;}

 protected:
  AliMUONRawReader();                  // Default constructor
  AliMUONRawReader (const AliMUONRawReader& rhs); // copy constructor
  AliMUONRawReader& operator=(const AliMUONRawReader& rhs); // assignment operator

 private:

  AliMUONData*  fMUONData;           //! Data container for MUON subsystem 
 
  AliMpSegFactory* fSegFactory;      //! Mapping segmentation factory
   
  AliMUONDDLTracker* fDDLTracker;    //! DDL tracker class pointers
  AliMUONDDLTrigger* fDDLTrigger;    //! DDL trigger class pointers

  AliMpBusPatch* fBusPatchManager;   //! buspatch versus DE's & DDL

  Bool_t fScalerEvent;               // flag to generates scaler event

  TStopwatch fTrackerTimer; //!
  TStopwatch fTriggerTimer; //!
  TStopwatch fMappingTimer; //!
  
  ClassDef(AliMUONRawReader,0) // MUON cluster reconstructor in ALICE
};
	
#endif
