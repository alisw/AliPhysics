#ifndef ALIMUONRAWDATA_H
#define ALIMUONRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
class TClonesArray;
class AliLoader;
class AliMUONData;
class AliMUONDigit;
class AliMUONDDLTracker;

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

 protected:
  AliMUONRawData();                  // Default constructor
  AliMUONRawData (const AliMUONRawData& rhs); // copy constructor
  AliMUONRawData& operator=(const AliMUONRawData& rhs); // assignment operator

 private:
  static const Int_t fgkDefaultPrintLevel;     // Default print level

  Int_t                   fNCh;                // Number of chambers   
  Int_t                   fNTrackingCh;        // Number of tracking chambers*
  Int_t                   fNTriggerCh;         // Number of trigger chambers*

  AliMUONData*            fMUONData;           //! Data container for MUON subsystem 

 // print level
  Int_t fPrintLevel;

  // debug
  Int_t fDebug;
  
  // alice loader
  AliLoader* fLoader;

  // DDL binary file pointer one per 1/2 chamber
  FILE* fFile1;
  FILE* fFile2;

  Int_t WriteDDL(Int_t iCh);

  ClassDef(AliMUONRawData,0) // MUON cluster reconstructor in ALICE
};
	
#endif
