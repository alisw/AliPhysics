#ifndef ALIMUONRAWWRITER_H
#define ALIMUONRAWWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONRawWriter
/// \brief Raw data class for trigger and tracker chambers
///
/// Writring Raw data class for trigger and tracker chambers

#include <TObject.h>
#include <TClonesArray.h>
#include "AliMpBusPatch.h"
#include "AliMUONSubEventTracker.h"

class TArrayI;
class AliLoader;
class AliMUONData;
class AliMUONDigit;
class AliMUONDDLTracker;
class AliMUONDDLTrigger;
class AliMUONGlobalTrigger;
class AliMUONSubEventTrigger;
class AliRawReader;
class AliMUONGlobalTrigger;
class AliMpSegFactory;

class AliMUONRawWriter : public TObject 
{
 public:
  AliMUONRawWriter(AliLoader* loader, AliMUONData* data); // Constructor
  virtual ~AliMUONRawWriter(void); // Destructor
    
  // write raw data
  Int_t   Digits2Raw();

  AliMUONData*   GetMUONData() {return fMUONData;}

  void AddData(const AliMUONSubEventTracker* event) {
    TClonesArray &temp = *fSubEventArray;
    new(temp[temp.GetEntriesFast()])AliMUONSubEventTracker(*event); 
  }

  // could be private function (public for debugging)
  Int_t GetInvMapping(const AliMUONDigit* digit, Int_t &busPatchId,
		       UShort_t &manuId, UChar_t &channelId);

  Int_t GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg) const;

  void  SetScalerEvent() {fScalerEvent = kTRUE;}

 protected:
  AliMUONRawWriter();                  // Default constructor
  AliMUONRawWriter (const AliMUONRawWriter& rhs); // copy constructor
  AliMUONRawWriter& operator=(const AliMUONRawWriter& rhs); // assignment operator

 private:

  AliMUONData*  fMUONData;           //! Data container for MUON subsystem 
 
  AliLoader*    fLoader;             //! alice loader
 
  AliMpSegFactory* fSegFactory;      //! Mapping segmentation factory

  FILE*         fFile[2];            //! DDL binary file pointer one per 1/2 chamber

  TClonesArray* fSubEventArray;      //! array to sub event tracker
   
  AliMUONDDLTracker* fDDLTracker;    //! DDL tracker class pointers
  AliMUONDDLTrigger* fDDLTrigger;    //! DDL trigger class pointers

  AliMpBusPatch* fBusPatchManager;   //! buspatch versus DE's & DDL

  Bool_t fScalerEvent;               // flag to generates scaler event

  static Int_t fgManuPerBusSwp1B[12];   //! array containing the first manuId for each buspatch st1, Bending
  static Int_t fgManuPerBusSwp1NB[12];  //! array containing the first manuId for each buspatch st1, NBending

  static Int_t fgManuPerBusSwp2B[12];   //! array containing the first manuId for each buspatch st2, Bending
  static Int_t fgManuPerBusSwp2NB[12];  //! array containing the first manuId for each buspatch st2, NBending

  // writing raw data
  Int_t WriteTrackerDDL(Int_t iCh);
  Int_t WriteTriggerDDL();

  ClassDef(AliMUONRawWriter,1) // MUON cluster reconstructor in ALICE
};
	
#endif
