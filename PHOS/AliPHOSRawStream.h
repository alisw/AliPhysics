#ifndef ALIPHOSRAWSTREAM_H
#define ALIPHOSRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to PHOS digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
class TClonesArray ;

// --- AliRoot header files ---
#include "AliAltroRawStream.h"
class AliRawReader;
class AliPHOSConTableDB ;


class AliPHOSRawStream: public TObject {

public :
  
  AliPHOSRawStream(AliRawReader* rawReader);
  

 Bool_t ReadDigits(TClonesArray * digits) ;
 
//PHOS does not need this method
 virtual Bool_t    Next(){return kFALSE ;} ; 
 
 
 Int_t            GetColumn() const {return 0;}
 Int_t            GetModule() const {return 0;}
 Int_t            GetPrevColumn() const {return 0;}
 Int_t            GetPrevModule() const {return 0;}
 Int_t            GetPrevRow() const {return 0;}
 Int_t            GetRow() const {return 0;}
 Int_t            GetSignal() const {return 0;}
 Int_t            GetTime() const {return 0;}
 Bool_t           IsNewColumn() const {return kFALSE; }
 Bool_t           IsNewModule() const {return kFALSE;}
 Bool_t           IsNewRow() const {return kFALSE ;}
 

 void SetConTableDB(AliPHOSConTableDB * ctdb){fctdb = ctdb ;}

 Int_t GetTrigger(void) const {return fTrig ;}
 Bool_t IsLEDevent(void)const {return fTrig == kLED;}
 Bool_t IsPULevent(void)const {return fTrig == kPUL;}
 Bool_t IsPEDevent(void)const {return fTrig == kPED;}
 Bool_t IsWELevent(void)const {return fTrig == kWEL;}
 Bool_t IsNELevent(void)const {return fTrig == kNEL;}


 public:
 enum {	kLED = 129,     // Physics pattern unit for LED events
	kPUL =  33,     // Physics pattern unit for PULSER events
	kPED = 257,     // Physics pattern unit for PEDESTAL events
	kNEL =1029,     // Physics pattern unit for NARROW ELECTRON events
	kWEL =1027,     // Physics pattern unit for WIDE ELECTRON events
	kSOB =2048,     // Pattern unit mask for Start Of Burst trigger or 0 to disable SOB trigger
	kEOB =4096,     // Pattern unit mask for End   Of Burst trigger or 0 to disable EOB trigger

	kPattUnitMarker = 27,  // Equipment type marker for Pattern Unit
	kPattUnitEquipId= 64,  // Equipment ID for Pattern Unit

	kPhosAdcMarker = 22,   // Equipment type marker for PhosAdc (Kurchatov ADC)
	kPhosAdcEquipId= 16,   // Equipment ID for PhosAdc

	kTdcMarker = 26,       // Equipment type marker for Tdc
	kTdcEquipId=128,       // Equipment ID for Tdc

	kChargeAdcMarker = 24, // Equipment type marker for ChargeAdc
	kChargeAdcEquipId= 32, // Equipment ID for ChargeAdc

	kScalerMarker = 25,    // Equipment type marker for Scaler
	kScalerEquipId=256};   // Equipment ID for Scaler


protected :
  AliRawReader*       fRawReader; //! object for reading the raw data
  UChar_t*            fData;      //! raw data
  AliPHOSConTableDB * fctdb ;     //! connection between RAW index and AbsId of crystal 

  Int_t fTrig ; //current trigger

  ClassDef(AliPHOSRawStream, 0)   // class for reading PHOS raw digits
};

#endif
