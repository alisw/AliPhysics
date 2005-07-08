#ifndef ALIMUONRAWDATA_H
#define ALIMUONRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* Raw data class for trigger and tracker chambers */

#include <TObject.h>
#include "AliMUONSubEventTracker.h"

class TClonesArray;
class AliLoader;
class AliMUONData;
class AliMUONDigit;
class AliMUONDDLTracker;
class AliMUONDDLTrigger;
class AliMUONGlobalTrigger;
class AliMUONSubEventTrigger;
class AliRawReader;

class AliMUONRawData : public TObject 
{
 public:
  AliMUONRawData(AliLoader* loader); // Constructor
  virtual ~AliMUONRawData(void); // Destructor
    
  // write raw data
  Int_t   Digits2Raw();
  Int_t   Raw2Digits(AliRawReader* rawReader);

  AliMUONData*   GetMUONData() {return fMUONData;}

  void AddData1(const AliMUONSubEventTracker* event) {
    TClonesArray &temp = *fSubEventArray[0];
    new(temp[temp.GetEntriesFast()])AliMUONSubEventTracker(*event); 
  }

  void AddData2(const AliMUONSubEventTracker* event) {
    TClonesArray &temp = *fSubEventArray[1];
    new(temp[temp.GetEntriesFast()])AliMUONSubEventTracker(*event); 
  }

  void GetDummyMapping(Int_t iCh, Int_t iCath, const AliMUONDigit* digit, Int_t &busPatchId,
		       UShort_t &manuId, UChar_t &channelId);

  void GetInvDummyMapping(Int_t iCh, Int_t buspatchId, UShort_t manuId, 
			  UChar_t channelId, AliMUONDigit* digit );


  Int_t GetGlobalTriggerPattern(const AliMUONGlobalTrigger* gloTrg);

 protected:
  AliMUONRawData();                  // Default constructor
  AliMUONRawData (const AliMUONRawData& rhs); // copy constructor
  AliMUONRawData& operator=(const AliMUONRawData& rhs); // assignment operator

 private:
  static const Int_t fgkDefaultPrintLevel;     // Default print level

  AliMUONData*  fMUONData;           //! Data container for MUON subsystem 
 
  AliLoader*    fLoader;             //! alice loader
 
  FILE*         fFile1;              //! DDL binary file pointer one per 1/2 chamber
  FILE*         fFile2;              //! DDL binary file pointer one per 1/2 chamber

  TClonesArray* fSubEventArray[2];   //! array to sub event tracker
  
  TClonesArray* fSubEventTrigArray[2]; //! array to sub event trigger
 
  AliMUONDDLTracker* fDDLTracker;      //! DDL tracker class pointers
  AliMUONDDLTrigger* fDDLTrigger;      //! DDL trigger class pointers

  // writing raw data
  Int_t WriteTrackerDDL(Int_t iCh);
  Int_t WriteTriggerDDL();

  Int_t ReadTrackerDDL(AliRawReader* rawReader);
  Int_t ReadTriggerDDL(AliRawReader* rawReader);

  ClassDef(AliMUONRawData,1) // MUON cluster reconstructor in ALICE
};
	
#endif
