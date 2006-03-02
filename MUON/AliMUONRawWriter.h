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
#include "AliMUONSubEventTracker.h"
#include "TStopwatch.h"

class AliMUONDDLTracker;
class AliMUONDDLTrigger;
class AliMUONData;
class AliMUONDigit;
class AliMUONGlobalTrigger;
class AliMUONSubEventTrigger;
class AliMpBusPatch;
class AliMpSegFactory;

class AliMUONRawWriter : public TObject 
{
 public:
  AliMUONRawWriter(AliMUONData* data); // Constructor
  virtual ~AliMUONRawWriter(void); // Destructor
    
  // write raw data
  Int_t Digits2Raw();

  void  SetScalerEvent() {fScalerEvent = kTRUE;}
  
private:

  void AddData(const AliMUONSubEventTracker& event)
  {
    TClonesArray &temp = *fSubEventArray;
    new(temp[temp.GetEntriesFast()]) AliMUONSubEventTracker(event); 
  }

  Int_t GetBusPatch(const AliMUONDigit& digit);
  
  Int_t GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg) const;

protected:
  AliMUONRawWriter();                  // Default constructor
  AliMUONRawWriter (const AliMUONRawWriter& rhs); // copy constructor
  AliMUONRawWriter& operator=(const AliMUONRawWriter& rhs); // assignment operator

  // writing raw data
  Int_t WriteTrackerDDL(Int_t iCh);
  Int_t WriteTriggerDDL();
  
 private:

  AliMUONData*  fMUONData;           //! Data container for MUON subsystem 
 
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
  
  TStopwatch fTrackerTimer; //!
  TStopwatch fTriggerTimer; //!
  TStopwatch fMappingTimer; //!
  
  AliMpSegFactory* fSegFactory; //!
  
  ClassDef(AliMUONRawWriter,0) // MUON cluster reconstructor in ALICE
};
	
#endif
