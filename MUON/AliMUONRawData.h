#ifndef ALIMUONRAWDATA_H
#define ALIMUONRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include "AliMUONSubEventTracker.h"
#include "AliMUONSubEventTrigger.h"

class TClonesArray;
class AliLoader;
class AliMUONData;
class AliMUONDigit;
class AliMUONDDLTracker;
class AliMUONDDLTrigger;

//class AliMUONTriggerDecision;

class AliMUONRawData : public TObject 
{
 public:
  AliMUONRawData(AliLoader* loader); // Constructor
  virtual ~AliMUONRawData(void); // Destructor
    
  // write raw data
  Int_t   WriteRawData();

  AliMUONData*   GetMUONData() {return fMUONData;}

  Int_t GetPrintLevel(void) const {return fPrintLevel;}
  void SetPrintLevel(Int_t printLevel) {fPrintLevel = printLevel;}

  void AddData1(AliMUONSubEventTracker* event) {
    TClonesArray &temp = *fSubEventArray[0];
    new(temp[temp.GetEntriesFast()])AliMUONSubEventTracker(*event); 
  }

  void AddData2(AliMUONSubEventTracker* event) {
    TClonesArray &temp = *fSubEventArray[1];
    new(temp[temp.GetEntriesFast()])AliMUONSubEventTracker(*event); 
  }

  void GetDummyMapping(Int_t iCh, Int_t iCath, const AliMUONDigit* digit, Int_t &busPatchId,
		       UShort_t &manuId, UChar_t &channelId, UShort_t &charge);


 protected:
  AliMUONRawData();                  // Default constructor
  AliMUONRawData (const AliMUONRawData& rhs); // copy constructor
  AliMUONRawData& operator=(const AliMUONRawData& rhs); // assignment operator

 private:
  static const Int_t fgkDefaultPrintLevel;     // Default print level


  AliMUONData*            fMUONData;           //! Data container for MUON subsystem 
/*   AliMUONTriggerDecision* fTrigDec;            //! calculated trigger from digits tmp solution */
/*   AliMUONData*            fTrigData; */

 // print level
  Int_t fPrintLevel;

  // debug
  Int_t fDebug;
  
  // alice loader
  AliLoader* fLoader;

  // DDL binary file pointer one per 1/2 chamber
  FILE* fFile1;
  FILE* fFile2;

  // array to sub event tracker
  TClonesArray* fSubEventArray[2];

  // array to sub event trigger
  TClonesArray* fSubEventTrigArray[2];

  // DDL class pointers
  AliMUONDDLTracker* fDDLTracker;
  AliMUONDDLTrigger* fDDLTrigger;

  // writing raw data
  Int_t WriteTrackerDDL(Int_t iCh);
  Int_t WriteTriggerDDL();

  ClassDef(AliMUONRawData,0) // MUON cluster reconstructor in ALICE
};
	
#endif
