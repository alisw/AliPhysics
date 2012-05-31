#ifndef ALITPCMONITORDATEFORMAT_H
#define ALITPCMONITORDATEFORMAT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorDateFormat class
////
//// Class for decoding raw data headers in DATE format
//// 
//// Authors: Roland Bramm, 
////          Stefan Kniege, IKF, Frankfurt
////       
/////////////////////////////////////////////////////////////////////////

struct eventHeaderStruct;
struct equipmentHeaderStruct;

#define long32 int
#include "TNamed.h"
using namespace std;

class AliTPCMonitorDateFormat : public TNamed {
 public:
    AliTPCMonitorDateFormat(Char_t* data);
    AliTPCMonitorDateFormat(const  AliTPCMonitorDateFormat &dateformat);
    AliTPCMonitorDateFormat& operator= (const AliTPCMonitorDateFormat& dateformat); 
    ~AliTPCMonitorDateFormat(); 

    //Super Event Header
    Int_t   GetEventSize() const;
    Int_t   GetEventHeaderSize() const;
    Int_t   GetEventHeaderBaseSize() const; 
    Int_t   GetEventID() const;
    Int_t   GetEventLDC() const;
    Int_t   GetEventGDC() const;

    Int_t   GetEventRunID() const;
    Int_t   GetEventVersion() const;
    Int_t   GetEventVersionMajor() const;
    Int_t   GetEventVersionMinor() const;
    Bool_t  IsEventSuperEvent() const;
    Bool_t  IsEventStartOfRun() const;
    Bool_t  IsEventEndOfRun() const;
    Bool_t  IsEventPhysicsEvent() const;
    Bool_t  IsEventSwapped() const;
    Bool_t  IsEventWrongEndian() const;
 
    //Sub Event Header
    void    GotoSubEventHeader();
    void    GotoNextSubEventHeader();
    Bool_t  IsLastSubEventHeader() const;

    Int_t   GetSubEventSize() const;
    Int_t   GetSubEventHeaderSize() const;
    Int_t   GetSubEventLDC() const;
    Int_t   GetSubEventGDC() const;

    Bool_t  IsSubEventSuperEvent();
    Bool_t  IsSubEventStartOfRun() const;
    Bool_t  IsSubEventEndOfRun() const;
    Bool_t  IsSubEventPhysicsEvent() const;

    //Eqipments
    void    GotoFirstEquipment();
    void    GotoNextEquipment();
    Bool_t  IsLastEquipment() const;

    Int_t   GetEquipmentSize() const;
    Int_t   GetEquipmentType() const;
    Int_t   GetEquipmentID() const; 
    Int_t*  GetEquipmentTypeAttribute();
    Int_t   GetEquipmentBasicSize()  const;
    Int_t   GetEquipmentHeaderSize() const; 
    Int_t   GetPayloadSize() const;

    //DATA
    Char_t* GetFirstDataPointer();
    Int_t   GetPosition() const;
    Int_t   GetPositionSubEvent() const;
    
 private:
    Char_t*                fdataPtr;       // pointer to data array (start, will not be changed in event) 
    Char_t*                fsubEventPtr;   // pointer to SubEvent
    Char_t*                fcurrentPtr;    // pointer to current data position (header or data)
    eventHeaderStruct*     fevent;         // event and
    eventHeaderStruct*     fsubEvent;      // subevent structure
    equipmentHeaderStruct* fequipment;     // equipmemnt structure
  
    ClassDef(AliTPCMonitorDateFormat,1);
};

#endif
 
