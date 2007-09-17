#ifndef ALITPCMONITORDATEFORMAT_H
#define ALITPCMONITORDATEFORMAT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitorDateFormat class
//
// Class for decoding and reading raw data headers in DATE format
// 
// Authors: Roland Bramm, 
//          Stefan Kniege, IKF, Frankfurt
//       
/////////////////////////////////////////////////////////////////////////

struct eventHeaderStruct;
struct equipmentHeaderStruct;

#define long32 int
#include "TNamed.h"
using namespace std;

class AliTPCMonitorDateFormat : public TNamed {
 public:
    AliTPCMonitorDateFormat(Char_t* data);
    ~AliTPCMonitorDateFormat();

    //Super Event Header
    Int_t   GetEventSize();
    Int_t   GetEventHeaderSize();
    Int_t   GetEventHeaderBaseSize();
    Int_t   GetEventID();
    Int_t   GetEventLDC();
    Int_t   GetEventGDC();

    Int_t   GetEventRunID();
    Int_t   GetEventVersion();
    Int_t   GetEventVersionMajor();
    Int_t   GetEventVersionMinor();
    Bool_t  IsEventSuperEvent();
    Bool_t  IsEventStartOfRun();
    Bool_t  IsEventEndOfRun();
    Bool_t  IsEventPhysicsEvent();
    Bool_t  IsEventSwapped();
    Bool_t  IsEventWrongEndian();
 
    //Sub Event Header
    void    GotoSubEventHeader();
    void    GotoNextSubEventHeader();
    Bool_t  IsLastSubEventHeader();

    Int_t   GetSubEventSize();
    Int_t   GetSubEventHeaderSize();
    Int_t   GetSubEventLDC();
    Int_t   GetSubEventGDC();

    Bool_t  IsSubEventSuperEvent();
    Bool_t  IsSubEventStartOfRun();
    Bool_t  IsSubEventEndOfRun();
    Bool_t  IsSubEventPhysicsEvent();

    //Eqipments
    void    GotoFirstEquipment();
    void    GotoNextEquipment();
    Bool_t  IsLastEquipment();

    Int_t   GetEquipmentSize();
    Int_t   GetEquipmentType();
    Int_t   GetEquipmentID();
    Int_t*  GetEquipmentTypeAttribute();
    Int_t   GetEquipmentBasicSize();
    Int_t   GetEquipmentHeaderSize();
    Int_t   GetPayloadSize();

    //DATA
    Char_t* GetFirstDataPointer();
    Int_t   GetPosition();
    Int_t   GetPositionSubEvent();
    
 private:
    Char_t*                       fdataPtr;       // pointer to data array (start, will not be changed in event) 
    Char_t*                       fsubEventPtr;   // pointer to SubEvent
    Char_t*                       fcurrentPtr;    // pointer to current data position (header or data)
    eventHeaderStruct*     event;          // event and
    eventHeaderStruct*     subEvent;       // subevent structure
    equipmentHeaderStruct* equipment;      // equipmemnt structure
  
    ClassDef(AliTPCMonitorDateFormat,1);
};

#endif
 
