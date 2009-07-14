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

#include <TClonesArray.h>
#include "AliVHeader.h"
#include "AliTriggerScalersESD.h"
#include "AliTriggerScalersRecordESD.h"

class AliTriggerScalersESD;
class AliTriggerScalersRecordESD;
class AliTriggerIR;

class AliESDHeader: public AliVHeader {
public:
  AliESDHeader();
  virtual ~AliESDHeader();
  AliESDHeader(const AliESDHeader& header);
  AliESDHeader& operator=(const AliESDHeader& header);
  virtual void Copy(TObject &obj) const;

  void      SetTriggerMask(ULong64_t n) {fTriggerMask=n;}
  void      SetOrbitNumber(UInt_t n) {fOrbitNumber=n;}
  void      SetTimeStamp(UInt_t timeStamp){fTimeStamp = timeStamp;}
  void      SetEventType(UInt_t eventType){fEventType = eventType;}
  void      SetEventNumberInFile(Int_t n) {fEventNumberInFile=n;}
  void      SetBunchCrossNumber(UShort_t n) {fBunchCrossNumber=n;}
  void      SetPeriodNumber(UInt_t n) {fPeriodNumber=n;}
  void      SetTriggerCluster(UChar_t n) {fTriggerCluster = n;}
  Bool_t    AddTriggerIR(const AliTriggerIR* ir);

//************Setters/Getters for Trigger Inputs and TriggerScalersRecordESD
  void SetL0TriggerInputs(UInt_t n) {fL0TriggerInputs=n;}
  void SetL1TriggerInputs(UInt_t n) {fL1TriggerInputs=n;}
  void SetL2TriggerInputs(UShort_t n) {fL2TriggerInputs=n;}
  UInt_t      GetL0TriggerInputs() const {return fL0TriggerInputs;}  
  UInt_t      GetL1TriggerInputs() const {return fL1TriggerInputs;} 
  UShort_t    GetL2TriggerInputs() const {return fL2TriggerInputs;} 
  void SetTriggerScalersRecord(AliTriggerScalersESD *scalerRun) {fTriggerScalers.AddTriggerScalers(scalerRun); }
  const AliTriggerScalersRecordESD* GetTriggerScalersRecord() const {return &fTriggerScalers; }
  const AliTriggerIR* GetTriggerIR(Int_t i) const { return fIRArray[i]; }
//**************************************************************************

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
  Int_t        fEventNumberInFile; // Running Event count in the file
  UShort_t     fBunchCrossNumber;  // Bunch Crossing Number
  UChar_t      fTriggerCluster;    // Trigger cluster (mask)
  UInt_t       fL0TriggerInputs;   //L0 Trigger Inputs 
  UInt_t       fL1TriggerInputs;   //L1 Trigger Inputs
  UShort_t     fL2TriggerInputs;   //L2 Trigger Inputs
  AliTriggerScalersRecordESD fTriggerScalers;  //Trigger counters of triggered classes in event
  enum {kNMaxIR = 3};              // Max number of interaction records (IR)
  AliTriggerIR*  fIRArray[kNMaxIR];// Array with trigger IRs 

  ClassDef(AliESDHeader,6)
};

#endif
