// -*- mode: C++ -*- 
#ifndef ALIESDHEADER_H
#define ALIESDHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Class AliESDHeader
//   Header data
//   for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

#include "AliVHeader.h"

class AliESDHeader: public AliVHeader {
public:
  AliESDHeader();
  virtual ~AliESDHeader();
  AliESDHeader(const AliESDHeader& header);
  AliESDHeader& operator=(const AliESDHeader& header);

  void      SetTriggerMask(ULong64_t n) {fTriggerMask=n;}
  void      SetOrbitNumber(UInt_t n) {fOrbitNumber=n;}
  void      SetTimeStamp(UInt_t timeStamp){fTimeStamp = timeStamp;}
  void      SetEventType(UInt_t eventType){fEventType = eventType;}
  void      SetEventNumberInFile(Int_t n) {fEventNumberInFile=n;}
  void      SetBunchCrossNumber(UShort_t n) {fBunchCrossNumber=n;}
  void      SetPeriodNumber(UInt_t n) {fPeriodNumber=n;}
  void      SetTriggerCluster(UChar_t n) {fTriggerCluster = n;}

  ULong64_t GetTriggerMask() const {return fTriggerMask;}
  UInt_t    GetOrbitNumber() const {return fOrbitNumber;}
  UInt_t    GetTimeStamp()  const { return fTimeStamp;}
  UInt_t    GetEventType()  const { return fEventType;}
  Int_t     GetEventNumberInFile() const {return fEventNumberInFile;}
  UShort_t  GetBunchCrossNumber() const {return fBunchCrossNumber;}
  UInt_t    GetPeriodNumber() const {return fPeriodNumber;}
  UChar_t   GetTriggerCluster() const {return fTriggerCluster;}

  void      Reset();
  void      Print(const Option_t *opt=0) const;
private:

  // Event Identification
  ULong64_t    fTriggerMask;       // Trigger Type (mask)
  UInt_t       fOrbitNumber;       // Orbit Number
  UInt_t       fTimeStamp;         // Time stamp
  UInt_t       fEventType;         // Type of Event
  UInt_t       fPeriodNumber;      // Period Number
  Int_t        fEventNumberInFile; // running Event count in the file
  UShort_t     fBunchCrossNumber;  // Bunch Crossing Number
  UChar_t      fTriggerCluster;    // Trigger cluster (mask)
  
  ClassDef(AliESDHeader,3)
};

#endif
