#ifndef ALIITSEVENTHEADER_H
#define ALIITSEVENTHEADER_H

////////////////////////////////////////////////////
//  Base class to define                          //
//  ITS Event Header                              //
//  
//  Origin: E. Crescio crescio@to.infn.it         //
//                                                //
////////////////////////////////////////////////////

#include "AliDetectorEventHeader.h"

typedef enum { kSimulated, kReal, kCalibration1, kCalibration2 }  EventType_t;

class AliITSEventHeader : public AliDetectorEventHeader {
 
 public:


  AliITSEventHeader(const char* name);
  AliITSEventHeader();

  virtual ~AliITSEventHeader() {}
  
  EventType_t GetEventTypeSDD()   const {return fEventTypeSDD;};
  UChar_t     GetL1TriggerType(Int_t i) const {return fL1TriggerType[i];}
  UInt_t      GetOrbitNumber(Int_t i) const  {return fOrbitNumber[i];}
  UShort_t    GetBunchCross(Int_t i)  const  {return fBunchCross[i];}
  UChar_t     GetBlockAttributes(Int_t i) const {return fBlockAttr[i];}
  ULong64_t   GetTriggerClass(Int_t i) const {return fTriggerClass[i];}
  UInt_t      GetStatusBits(Int_t i)  const {return fStatusBits[i];}
  UInt_t      GetMiniEvId(Int_t i)   const {return fMiniEvId[i];}
  UInt_t      GetSubDet(Int_t i)     const {return fSubDet[i];}
  UInt_t      GetVersion(Int_t i)    const {return fVersion[i];}
  Int_t       GetJitterSDD()         const {return fJitterSDD;}

  void SetEventTypeSDD(EventType_t type=kSimulated){fEventTypeSDD=type;}
  void SetL1TriggerType(Int_t i,UChar_t l1trig) {fL1TriggerType[i]=l1trig;}
  void SetOrbitNumber(Int_t i,UInt_t orbitnum) {fOrbitNumber[i]=orbitnum;}
  void SetBunchCross(Int_t i,UShort_t bunchcross) {fBunchCross[i]=bunchcross;}
  void SetBlockAttributes(Int_t i,UChar_t attr) {fBlockAttr[i]=attr;}
  void SetTriggerClass(Int_t i,ULong64_t trigclass){fTriggerClass[i]=trigclass;}
  void SetStatusBits(Int_t i,UInt_t bits) {fStatusBits[i]=bits;}
  void SetMiniEvId(Int_t i,UInt_t minievid) {fMiniEvId[i]=minievid;}
  void SetSubDet(Int_t i,UInt_t subdet) {fSubDet[i]=subdet;}
  void SetVersion(Int_t i,UInt_t version) {fVersion[i]=version;}
  void SetJitterSDD(Int_t jitter) {fJitterSDD=jitter;}
  
 protected:
  
  EventType_t fEventTypeSDD;       //event type for SDD 
  UChar_t   fL1TriggerType[3];     //Level 1 trigger type (0 spd,1 sdd,2 ssd) 
  UInt_t    fOrbitNumber[3];       //Orbit Number (0 spd,1 sdd,2 ssd)
  UShort_t  fBunchCross[3];        //Bunch Crossing Number (0 spd,1 sdd,2 ssd)
  UChar_t   fBlockAttr[3];         //Block Attributes (0 spd,1 sdd,2 ssd)
  ULong64_t fTriggerClass[3];      //Trigger classes (0 spd,1 sdd,2 ssd)
  UInt_t    fStatusBits[3];        //Status Bits (0 spd,1 sdd,2 ssd)
  UInt_t    fMiniEvId[3];          //MiniEvent ID (0 spd,1 sdd,2 ssd)
  UInt_t    fSubDet[3];            //partic. sub-detectors (0 spd,1 sdd,2 ssd)
  UInt_t    fVersion[3];           //Header Version (0 spd,1 sdd,2 ssd)
  Int_t     fJitterSDD;            // SDD jitter between L0 and pascal stop

  ClassDef(AliITSEventHeader,1)  // An Alice ITS event header 

 };


#endif

    
