/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: $ */

///
/// @file   AliHLTMUONEvent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   29 Sep 2007
/// @brief  Implementation of the AliHLTMUONEvent class.
///
/// The class is used to store all ROOTified data objects from the dHLT chain
/// for a single event together.

#include "AliHLTMUONEvent.h"
#include "AliHLTMUONDecision.h"
#include "AliLog.h"
#include <iostream>

ClassImp(AliHLTMUONEvent);


AliHLTMUONEvent::AliHLTMUONEvent(const AliHLTMUONEvent& event) :
	TObject(event),
	fEventId(event.fEventId),
	fArray()
{
	/// Copy constructor performs a deep copy of the object.
	
	fEventId = event.fEventId;
	fArray.SetOwner(kTRUE);
	TIter next(&event.fArray);
	TObject* obj = NULL;
	while ( (obj = next()) != NULL )
	{
		fArray.Add(obj->Clone());
	}
}


AliHLTMUONEvent& AliHLTMUONEvent::operator = (const AliHLTMUONEvent& event)
{
	/// The assignment operator performs a deep copy of the object.
	
	if (this == &event) return *this;
	TObject::operator = (event);
	fEventId = event.fEventId;
	fArray.Clear();
	TIter next(&event.fArray);
	TObject* obj = NULL;
	while ( (obj = next()) != NULL )
	{
		fArray.Add(obj->Clone());
	}
	return *this;
}


const AliHLTMUONDecision* AliHLTMUONEvent::FindDecision() const
{
	/// Finds the decision object in the event from the list of dHLT objects.
	/// There should only be one such object in the event. If not, then only
	/// the first object found is returned. You will need to manually search
	/// for the other objects.
	/// \returns  The AliHLTMUONDecision object in the event or NULL if none exists.
	
	for (Int_t i = 0; i < fArray.GetEntriesFast(); i++)
	{
		if (fArray[i]->IsA() == AliHLTMUONDecision::Class())
		{
			return static_cast<const AliHLTMUONDecision*>(fArray[i]);
		}
	}
	
	return NULL;
}


void AliHLTMUONEvent::Print(Option_t* option) const
{
	///
	/// Inherited from TObject. Prints the contents of the event objects in fArray.
	/// \param option  This is an option string that is just passed on to individual
	///     objects in the event's fArray list of objects.
	///
	
	std::cout << "################## EVENT: " << fEventId << " ##################" << std::endl;
	for (Int_t i = 0; i < fArray.GetEntriesFast(); i++)
		if (fArray[i] != NULL) fArray[i]->Print(option);
}


void AliHLTMUONEvent::Clear(Option_t* option)
{
	/// Clears the internal array of event objects.
	
	fEventId = AliHLTEventID_t(-1);
	fArray.Clear(option);
}


void AliHLTMUONEvent::Copy(TObject& object) const
{
	/// Deep copy this object to the target object.
	/// \param object  The target object to copy to.
	
	if (object.IsA() != this->IsA())
	{
		AliError(Form("Cannot copy an object of type %s to a type of %s.",
			this->ClassName(), object.ClassName()
		));
		return;
	}
	
	TObject::Copy(object);
	AliHLTMUONEvent& event = static_cast<AliHLTMUONEvent&>(object);
	event.fEventId = fEventId;
	event.fArray.Clear();
	TIter next(&fArray);
	TObject* obj = NULL;
	while ( (obj = next()) != NULL )
	{
		event.fArray.Add(obj->Clone());
	}
}

