#ifndef ALIPHOSRAWREADERDATE_H
#define ALIPHOSRAWREADERDATE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading raw data from a date file or event.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliRawReader.h"

#ifdef ALI_DATE
#include "event.h"
#include "equipment.h"
#else
#include "AliPHOSevent.h"
#include "AliPHOSequipment.h"
#endif


class AliPHOSRawReaderDate: public AliRawReader {
  public :
    AliPHOSRawReaderDate(void* event);
    AliPHOSRawReaderDate(const char* fileName, Int_t eventNumber = -1);
    virtual ~AliPHOSRawReaderDate();

    void             RequireHeader(Bool_t required = kTRUE)
      {fRequireHeader = required;};

    virtual UInt_t   GetType() const       { return fEvent ? fEvent->eventType : 0 ; }
    virtual UInt_t   GetRunNumber() const  { return fEvent ? fEvent->eventRunNb: 0 ;}
    virtual const UInt_t* GetEventId() const{return 0; } // ??? fEvent ? fEvent->eventId : 0 ; } //?? 
    virtual const UInt_t* GetTriggerPattern()  const{return 0; } // ??? (!fEvent) ? 0 : fEvent->eventTriggerPattern; } //??
    virtual const UInt_t* GetDetectorPattern() const{return fEvent  ? fEvent->detectorId : 0 ;}
    virtual const UInt_t* GetAttributes() const{return fEvent ? fEvent->eventTypeAttribute : 0 ; } 
    virtual UInt_t   GetLDCId() const {return  0 ; } 
    virtual UInt_t   GetGDCId() const {return  0 ; }

    virtual Int_t    GetEquipmentSize() const {return fEquipment ? fEquipment->rawDataLen : 0 ;  }
    virtual Int_t    GetEquipmentType() const {return fEquipment ? fEquipment->type : 0 ;  }
    virtual Int_t    GetEquipmentId() const   {return fEquipment ? fEquipment->equipmentId : 0 ; }
    virtual const UInt_t* GetEquipmentAttributes() const {return  0 ;}
    virtual Int_t    GetEquipmentElementSize() const {return 0 ;}

    virtual Bool_t   ReadHeader();
    virtual Bool_t   ReadNextData(UChar_t*& data);

    virtual Bool_t   Reset();

    virtual Bool_t   NextEvent();
    virtual Bool_t   RewindEvents();

    virtual Int_t    CheckData() const;
    
    virtual const UInt_t* GetSubEventAttributes() const { return 0;}
  protected :
 
  inline void ChangeOrder(Int_t & dword) {
           dword = dword << 24 | (dword >> 24) & 0xff | (dword & 0xff00) << 8 | (dword & 0xff00) >>8; }
  inline void ChangeOrder(Long_t & dword){
           dword = dword << 24 | (dword >> 24) & 0xff | (dword & 0xff00) << 8 | (dword & 0xff00) >>8; }
  inline void ChangeOrder(ULong_t & dword){
	   dword = dword << 24 | (dword >> 24) & 0xff | (dword & 0xff00) << 8 | (dword & 0xff00) >>8; }
  inline void ChangeOrder(UInt_t & dword){
           dword = dword << 24 | (dword >> 24) & 0xff |	(dword & 0xff00) << 8 | (dword & 0xff00) >>8; }
  inline void ChangeOrder(Short_t & word){
	word = word << 8 | (word >> 8) & 0xff; }
  inline void ChangeOrder(UShort_t & word){
        word = word << 8 | (word >> 8) & 0xff; }
  void SwappEvent(eventHeaderStruct * event) ;

  virtual Bool_t   ReadNext(UChar_t* data, Int_t size);
  
  Bool_t           fRequireHeader; // if false, data without header is accepted
  
  FILE*            fFile;         // DATE file
  eventHeaderStruct* fEvent;      // raw data super event
  eventHeaderStruct* fSubEvent;   // raw data sub event
  equipmentHeaderStruct* fEquipment; // raw data equipment header
  
  UChar_t*         fPosition;     // current position in the raw data
  UChar_t*         fEnd;          // end position of the current data block
  
 private:
    AliPHOSRawReaderDate(const AliPHOSRawReaderDate& rawReader);
    AliPHOSRawReaderDate& operator = (const AliPHOSRawReaderDate& rawReader);
    
    ClassDef(AliPHOSRawReaderDate, 0) // class for reading raw digits from a root file
};
      
#endif
