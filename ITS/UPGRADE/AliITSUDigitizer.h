#ifndef ALIITSUDIGITIZER_H
#define ALIITSUDIGITIZER_H
/* Copyright (c) 1998-2001, ALICE Experiment at CERN, All rights reserved *
 * See cxx source for full Copyright notice                               */

/*
  $Id: AliITSUDigitizer.h 52261 2011-10-23 15:46:57Z hristov $
 */
//////////////////////////////////////////////////////////////////
// Digitizer class for ITS                                      //
//////////////////////////////////////////////////////////////////
class TObjArray;
class TTree;

class AliDigitizationInput;

#include "AliDigitizer.h" // Base class from which this one is derived
#include "AliITSU.h"   // ITS class functions used in inline functions.

class AliITSUDigitizer : public AliDigitizer {
 public:
  AliITSUDigitizer();
  AliITSUDigitizer(AliDigitizationInput* digInput);
  
  virtual ~AliITSUDigitizer();
  virtual Bool_t Init();
  virtual void Digitize(Option_t* opt=0);
  virtual void SetChipActive(Int_t i){if(fModActive) fModActive[i] = kTRUE;}
  virtual void SetChipInActive(Int_t i){if(fModActive) fModActive[i] = kFALSE;}
  virtual void SetByRegionOfInterestFlag(Int_t i=0){fRoif = i;}
  virtual void SetByRegionOfFileNumber(Int_t i=-1){fRoiifile = i;}
  virtual void ClearByRegionOfInterestFlag(){fRoif = 0;}
  //
 private:
  AliITSUDigitizer(const AliITSUDigitizer& dig);
  AliITSUDigitizer& operator=(const AliITSUDigitizer &source);
  AliDigitizationInput* GetDigInput(){return fDigInput;}
  virtual void SetByRegionOfInterest(TTree *ts);
  //
 protected:
  AliITSU   *fITS;      //! local pointer to ITS
  Bool_t    *fModActive;//! flag to indicate which chip to digitize.
  Bool_t     fInit;     //! flag to indecate Initilization when well.
  Int_t      fRoif;     //! Region of interest flag.
  Int_t      fRoiifile; //! The file number with which to determing the region of interest from.
  Bool_t     fFlagFirstEv; //! Flag to control calibration access
  
  ClassDef(AliITSUDigitizer,1) // Task to Digitize ITS from summable hits.
};
#endif
