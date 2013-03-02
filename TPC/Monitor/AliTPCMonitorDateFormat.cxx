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
Revision 1.3  2007/10/12 13:36:27  cvetan
Coding convention fixes from Stefan

Revision 1.2  2007/09/17 16:34:54  cvetan
The package was overwriting the rootcint flags. This was fixed by applying the necessary changes in the DATE-dependent parts of the code

Revision 1.1  2007/09/17 10:23:31  cvetan
New TPC monitoring package from Stefan Kniege. The monitoring package can be started by running TPCMonitor.C macro located in macros folder.

*/ 


////////////////////////////////////////////////////////////////////////
////
//// AliTPCMonitorDateFormat class
////
//// Class for decoding raw data headers in DATE format
////
//// Reads event and subevent header informations form DATE files
//// 
//// Authors: Roland Bramm, 
////          Stefan Kniege, IKF, Frankfurt
////       
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
  fevent((struct eventHeaderStruct*) fdataPtr),
  fsubEvent(0),
  fequipment(0)
{
  // Constructor
}


//____________________________________________________________________________
AliTPCMonitorDateFormat::AliTPCMonitorDateFormat(const AliTPCMonitorDateFormat &dateformat) :
  TNamed(dateformat.GetName(),dateformat.GetTitle()),
  fdataPtr(dateformat.fdataPtr),
  fsubEventPtr(dateformat.fsubEventPtr),
  fcurrentPtr(dateformat.fcurrentPtr),
  fevent((struct eventHeaderStruct*)dateformat.fdataPtr),
  fsubEvent(dateformat.fsubEvent),
  fequipment(dateformat.fequipment)
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
      fevent=dateformat.fevent;
      fsubEvent=dateformat.fsubEvent;
      fequipment=dateformat.fequipment;
    }
  return *this;
}

//____________________________________________________________________________
AliTPCMonitorDateFormat::~AliTPCMonitorDateFormat()
{
  // Destructor
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventSize()const
{
  // Return event size
  return (Int_t)fevent->eventSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventHeaderSize()const
{
  // Return event header size
  return (Int_t)fevent->eventHeadSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventHeaderBaseSize()const
{
  // Return event header base size
	return (Int_t)EVENT_HEAD_BASE_SIZE;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventID()const
{
  // Return event ID
  return (Int_t)EVENT_ID_GET_NB_IN_RUN(fevent->eventId);
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventLDC()const
{
  // Return event LDC
  return (Int_t)fevent->eventLdcId;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventGDC()const
{
  // return event GDC
  return (Int_t)fevent->eventGdcId;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventRunID()const
{
  // Return run ID
  return (Int_t)fevent->eventRunNb;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventVersion()const
{
  // Return event version
  return (Int_t)fevent->eventVersion; 
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventVersionMajor()const
{
  // retrun event version (16-32 bit)
  return (Int_t)(GetEventVersion()>>16);
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEventVersionMinor()const
{
  // return event version (0-15 bit)
  return (Int_t)(GetEventVersion()&0x0000ffff);
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventSuperEvent()const
{
  // Check if event ist super event
  Bool_t retval;
  if(TEST_SYSTEM_ATTRIBUTE( fevent->eventTypeAttribute, ATTR_SUPER_EVENT ) == 1)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventStartOfRun()const
{
  // Check if event ist Start of Run
  Bool_t retval;
  if(fevent->eventType == START_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventEndOfRun()const
{
  // Check if event is End of Run
  Bool_t retval;
  if(fevent->eventType == END_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventPhysicsEvent()const
{
  // Check if event is Physics event
  Bool_t retval;
  if(fevent->eventType == PHYSICS_EVENT)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventSwapped()const
{
  // Check if event is swapped
  Bool_t retval;
  if(TEST_SYSTEM_ATTRIBUTE( fevent->eventTypeAttribute, ATTR_EVENT_SWAPPED ) == 1)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsEventWrongEndian()const
{
  // Check endian  
  Bool_t retval;
  if(EVENT_MAGIC_NUMBER == fevent->eventMagic)
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
    fsubEvent = (struct eventHeaderStruct*) (fsubEventPtr);
  }else{
    fsubEventPtr = fdataPtr;
		fsubEvent = (struct eventHeaderStruct*) (fsubEventPtr);
  }
}

//____________________________________________________________________________
void AliTPCMonitorDateFormat::GotoNextSubEventHeader()
{
  // set subevent pointer to next  sub event 
  fsubEventPtr += GetSubEventSize();
  fsubEvent = (struct eventHeaderStruct*) (fsubEventPtr);
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsLastSubEventHeader()const
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
Int_t AliTPCMonitorDateFormat::GetSubEventSize()const
{
  // Return sub event size
  return (Int_t)fsubEvent->eventSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetSubEventHeaderSize()const
{
  // Return sub event header size
  return (Int_t)fsubEvent->eventHeadSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetSubEventLDC()const
{
  // Return sub event LDC
  return (Int_t)fsubEvent->eventLdcId;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetSubEventGDC()const
{
  // return sub event GDC
	return (Int_t)fsubEvent->eventGdcId;
}


//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventSuperEvent()
{
  // Check if sub event is super event
	Bool_t retval;
	if(TEST_SYSTEM_ATTRIBUTE( fsubEvent->eventTypeAttribute, ATTR_SUPER_EVENT ) == 1)
		retval = true;
	else
		retval = false;
	return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventStartOfRun()const
{ 
  // Check if sub event is start of run
  Bool_t retval;
  if(fsubEvent->eventType == START_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventEndOfRun()const
{
  // Check if sub event is end of run
  Bool_t retval;
  if(fsubEvent->eventType == END_OF_RUN)
    retval = true;
  else
    retval = false;
  return retval;
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsSubEventPhysicsEvent()const
{
  // Check if sub event is physics event
  Bool_t retval;
  if(fsubEvent->eventType == PHYSICS_EVENT)
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
  fequipment = (struct equipmentHeaderStruct*) (fcurrentPtr);
}

//____________________________________________________________________________
void AliTPCMonitorDateFormat::GotoNextEquipment()
{
  // Set current pointer to next equipment
  fcurrentPtr += GetEquipmentSize();
  fequipment = (struct equipmentHeaderStruct*) (fcurrentPtr);
}

//____________________________________________________________________________
Bool_t AliTPCMonitorDateFormat::IsLastEquipment() const
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
Int_t AliTPCMonitorDateFormat::GetEquipmentSize() const
{
  // Return equipment size
  return (Int_t)fequipment->equipmentSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetPayloadSize() const
{
  // Return payload slze
  Int_t retval = 0;
  if(GetEventVersion() < 196610){
    retval = (Int_t)fequipment->equipmentSize;
  }else{
    retval = (Int_t)fequipment->equipmentSize - GetEquipmentHeaderSize();
  }
  return retval;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentType() const
{
  // Return equipment type
  return (Int_t)fequipment->equipmentType;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentID() const
{
  // Return equipment ID
  return (Int_t)fequipment->equipmentId;
}

//____________________________________________________________________________
Int_t* AliTPCMonitorDateFormat::GetEquipmentTypeAttribute()
{
  // Return equipment type attribute
  return (Int_t*)fequipment->equipmentTypeAttribute;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentBasicSize() const
{
  // Return equipment basic size
  return (Int_t)fequipment->equipmentBasicElementSize;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetEquipmentHeaderSize() const
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
Int_t AliTPCMonitorDateFormat::GetPosition() const
{
  // Return current position relative to start of event
  Int_t retval = (Int_t) (fcurrentPtr - fdataPtr);
  return retval;
}

//____________________________________________________________________________
Int_t AliTPCMonitorDateFormat::GetPositionSubEvent() const
{
  // Return subevent position  relative to start of event
  Int_t retval = (Int_t) (fsubEventPtr - fdataPtr);
  return retval;
}

