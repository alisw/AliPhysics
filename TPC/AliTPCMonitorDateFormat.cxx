/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.2  2007/09/17 16:34:54  cvetan
The package was overwriting the rootcint flags. This was fixed by applying the necessary changes in the DATE-dependent parts of the code

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 


////////////////////////////////////////////////////////////////////////
//
// AliTPCMonitorDateFormat class
//
// Class for decoding raw data headers in DATE format
// Reads event and subevent header informations form DATE files
// 
// Authors: Roland Bramm, 
//          Stefan Kniege, IKF, Frankfurt
//       
/////////////////////////////////////////////////////////////////////////

#include "AliTPCMonitorDateFormat.h"
#include "event.h"
#include <iostream>
ClassImp(AliTPCMonitorDateFormat)
//____________________________________________________________________________
AliTPCMonitorDateFormat::AliTPCMonitorDateFormat(Char_t* data): 
  fdataPtr(data),
  fsubEventPtr(data),
  fcurrentPtr(data),
  event((struct eventHeaderStruct*) fdataPtr),
  subEvent(0),
  equipment(0)
{
  // Constructor
}


//____________________________________________________________________________
AliTPCMonitorDateFormat::AliTPCMonitorDateFormat(const AliTPCMonitorDateFormat &dateformat) :
  TNamed(dateformat.GetName(),dateformat.GetTitle()),
  fdataPtr(dateformat.fdataPtr),
  fsubEventPtr(dateformat.fsubEventPtr),
  fcurrentPtr(dateformat.fcurrentPtr),
  event((struct eventHeaderStruct*)dateformat.fdataPtr),
  subEvent(dateformat.subEvent),
  equipment(dateformat.equipment)
{
  // copy constructor
}

//____________________________________________________________________________
AliTPCMonitorDateFormat &AliTPCMonitorDateFormat:: operator= (const AliTPCMonitorDateFormat& dateformat)
{

  // assignment operator 
  if(this!=&dateformat)
    {
      fdataPtr=dateformat.fdataPtr;
      fsubEventPtr=dateformat.fsubEventPtr;
      fcurrentPtr=dateformat.fcurrentPtr;
      event=dateformat.event;
      subEvent=dateformat.subEvent;
      equipment=dateformat.equipment;
    }
  return *this;
}

//____________________________________________________________________________
AliTPCMonitorDateFormat::~AliTPCMonitorDateFormat()
{
  // Destructor
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventSize()
{
  // Return event size
  return (Int_t)event->eventSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventHeaderSize()
{
  // Return event header size
  return (Int_t)event->eventHeadSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventHeaderBaseSize()
{
  // Return event header base size
	return (Int_t)EVENT_HEAD_BASE_SIZE;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventID()
{
  // Return event ID
  return (Int_t)EVENT_ID_GET_NB_IN_RUN(event->eventId);
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventLDC()
{
  // Return event LDC
  return (Int_t)event->eventLdcId;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventGDC()
{
  // return event GDC
  return (Int_t)event->eventGdcId;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventRunID()
{
  // Return run ID
  return (Int_t)event->eventRunNb;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventVersion()
{
  // Return event version
  return (Int_t)event->eventVersion;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventVersionMajor()
{
  // retrun event version (16-32 bit)
  return (Int_t)(GetEventVersion()>>16);
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventVersionMinor()
{
  // return event version (0-15 bit)
  return (Int_t)(GetEventVersion()&0x0000ffff);
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventSuperEvent()
{
  // Check if event ist super event
  Bool_t retval;
  if(TEST_SYSTEM_ATTRIBUTE( event->eventTypeAttribute, ATTR_SUPER_EVENT ) == 1)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventStartOfRun()
{
  // Check if event ist Start of Run
  Bool_t retval;
  if(event->eventType == START_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventEndOfRun()
{
  // Check if event is End of Run
  Bool_t retval;
  if(event->eventType == END_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventPhysicsEvent()
{
  // Check if event is Physics event
  Bool_t retval;
  if(event->eventType == PHYSICS_EVENT)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventSwapped()
{
  // Check if event is swapped
  Bool_t retval;
  if(TEST_SYSTEM_ATTRIBUTE( event->eventTypeAttribute, ATTR_EVENT_SWAPPED ) == 1)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventWrongEndian()
{
  // Check endian  
  Bool_t retval;
  if(EVENT_MAGIC_NUMBER == event->eventMagic)
    retval = false;
  else
    retval = true;
  return retval;
}


//____________________________________________________________________________
void AliTPCMonitorDateFormat::GotoSubEventHeader()
{
  // Set subevent Pointer to sub event
  if(IsEventSuperEvent() ==true){
    fsubEventPtr = fdataPtr+GetEventHeaderSize();
    subEvent = (struct eventHeaderStruct*) (fsubEventPtr);
  }else{
    fsubEventPtr = fdataPtr;
		subEvent = (struct eventHeaderStruct*) (fsubEventPtr);
  }
}

//____________________________________________________________________________
void AliTPCMonitorDateFormat::GotoNextSubEventHeader()
{
  // set subevent pointer to next  sub event 
  fsubEventPtr += GetSubEventSize();
  subEvent = (struct eventHeaderStruct*) (fsubEventPtr);
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsLastSubEventHeader()
{
  // Check if sub event is last sub event
  Bool_t retval;
  Int_t position;
  if(IsEventSuperEvent() ==true){
    position = fsubEventPtr - fdataPtr;
    if( (position+GetSubEventSize()) < GetEventSize() )
			retval = false;
    else
      retval = true;
  }else{
    position = fsubEventPtr - fdataPtr;
    if( (position+GetSubEventSize()) < GetEventSize() )
			retval = false;
    else
      retval = true;
  }
  return retval;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetSubEventSize()
{
  // Return sub event size
  return (Int_t)subEvent->eventSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetSubEventHeaderSize()
{
  // Return sub event header size
  return (Int_t)subEvent->eventHeadSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetSubEventLDC()
{
  // Return sub event LDC
  return (Int_t)subEvent->eventLdcId;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetSubEventGDC()
{
  // return sub event GDC
	return (Int_t)subEvent->eventGdcId;
}


//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventSuperEvent()
{
  // Check if sub event is super event
	Bool_t retval;
	if(TEST_SYSTEM_ATTRIBUTE( subEvent->eventTypeAttribute, ATTR_SUPER_EVENT ) == 1)
		retval = true;
	else
		retval = false;
	return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventStartOfRun()
{
  // Check if sub event is start of run
  Bool_t retval;
  if(subEvent->eventType == START_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventEndOfRun()
{
  // Check if sub event is end of run
  Bool_t retval;
  if(subEvent->eventType == END_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventPhysicsEvent()
{
  // Check if sub event is physics event
  Bool_t retval;
  if(subEvent->eventType == PHYSICS_EVENT)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
void AliTPCMonitorDateFormat::GotoFirstEquipment()
{
  // Set current pointer to first equipment
  fcurrentPtr = fsubEventPtr + GetSubEventHeaderSize();
  equipment = (struct equipmentHeaderStruct*) (fcurrentPtr);
}

//____________________________________________________________________________
void AliTPCMonitorDateFormat::GotoNextEquipment()
{
  // Set current pointer to next equipment
  fcurrentPtr += GetEquipmentSize();
  equipment = (struct equipmentHeaderStruct*) (fcurrentPtr);
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsLastEquipment()
{
  // Check if equipment is last equipment
  Bool_t retval;
  Int_t position;
  position = fcurrentPtr - fsubEventPtr;
  if( (position+GetEquipmentSize()) < GetSubEventSize() )
		  retval = false;
  else
    retval = true;
  return retval; 
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentSize()
{
  // Return equipment size
  return (Int_t)equipment->equipmentSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetPayloadSize()
{
  // Return payload slze
  Int_t retval = 0;
  if(GetEventVersion() < 196610){
    retval = (Int_t)equipment->equipmentSize;
  }else{
    retval = (Int_t)equipment->equipmentSize - GetEquipmentHeaderSize();
  }
  return retval;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentType()
{
  // Return equipment type
  return (Int_t)equipment->equipmentType;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentID()
{
  // Return equipment ID
  return (Int_t)equipment->equipmentId;
}

//____________________________________________________________________________
Int_t* AliTPCMonitorDateFormat::GetEquipmentTypeAttribute()
{
  // Return equipment type attribute
  return (Int_t*)equipment->equipmentTypeAttribute;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentBasicSize()
{
  // Return equipment basic size
  return (Int_t)equipment->equipmentBasicElementSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentHeaderSize()
{
  // Return equipment header size
  return sizeof(struct equipmentHeaderStruct);
}


//____________________________________________________________________________
Char_t *AliTPCMonitorDateFormat::GetFirstDataPointer()
{
  // Return data pointer (after equipment header)
  Char_t *datapointer;
  if(GetEventVersion() < 196610){
    fcurrentPtr += GetEquipmentHeaderSize();
    datapointer = fcurrentPtr;
  }else{
    datapointer = fcurrentPtr + GetEquipmentHeaderSize();
  }
  return datapointer;
 
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetPosition()
{
  // Return current position relative to start of event
  Int_t retval = (Int_t) (fcurrentPtr - fdataPtr);
  return retval;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetPositionSubEvent()
{
  // Return subevent position  relative to start of event
  Int_t retval = (Int_t) (fsubEventPtr - fdataPtr);
  return retval;
}

