#ifndef ALIMUONRAWDATA_H
#define ALIMUONRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup rec
/// \class AliMUONRawData
/// \brief Raw data class for trigger and tracker chambers
///
/// Raw data class for trigger and tracker chambers

#include <TObject.h>
#include "AliMpBusPatch.h"
#include "AliMUONSubEventTracker.h"

class TClonesArray;
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

class AliMUONRawData : public TObject 
{
 public:
  AliMUONRawData(AliLoader* loader); // Constructor
  virtual ~AliMUONRawData(void); // Destructor
    
  // write raw data
  Int_t   Digits2Raw();
  Int_t   Raw2Digits(AliRawReader* rawReader);

  Int_t ReadTrackerDDL(AliRawReader* rawReader);
  Int_t ReadTriggerDDL(AliRawReader* rawReader);

  AliMUONData*   GetMUONData() {return fMUONData;}

  void AddData(const AliMUONSubEventTracker* event) {
    TClonesArray &temp = *fSubEventArray;
    new(temp[temp.GetEntriesFast()])AliMUONSubEventTracker(*event); 
  }

  // could be private function (public for debugging)
  Int_t GetInvMapping(const AliMUONDigit* digit, Int_t &busPatchId,
		       UShort_t &manuId, UChar_t &channelId);

  Int_t GetMapping(Int_t buspatchId, UShort_t manuId, 
			  UChar_t channelId, AliMUONDigit* digit );

  Int_t GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg) const;
  AliMUONGlobalTrigger* GetGlobalTriggerPattern(Int_t gloTrg) const;


 protected:
  AliMUONRawData();                  // Default constructor
  AliMUONRawData (const AliMUONRawData& rhs); // copy constructor
  AliMUONRawData& operator=(const AliMUONRawData& rhs); // assignment operator

 private:

  AliMUONData*  fMUONData;           //! Data container for MUON subsystem 
 
  AliLoader*    fLoader;             //! alice loader
 
  FILE*         fFile[2];            //! DDL binary file pointer one per 1/2 chamber

  TClonesArray* fSubEventArray;      //! array to sub event tracker
   
  AliMUONDDLTracker* fDDLTracker;    //! DDL tracker class pointers
  AliMUONDDLTrigger* fDDLTrigger;    //! DDL trigger class pointers

  AliMpBusPatch* fBusPatchManager;    //! buspatch versus DE's & DDL

  // writing raw data
  Int_t WriteTrackerDDL(Int_t iCh);
  Int_t WriteTriggerDDL();

  ClassDef(AliMUONRawData,1) // MUON cluster reconstructor in ALICE
};
	
#endif
