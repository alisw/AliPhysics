#ifndef ALIHLTMUONEVENT_H
#define ALIHLTMUONEVENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///
/// @file   AliHLTMUONRootifierComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   9 May 2008
/// @brief  Declaration of class for storing ROOTified data objects for one event.
///

#include "TObjArray.h"
#include "AliHLTMUONDataTypes.h"

class AliHLTMUONDecision;


class AliHLTMUONEvent : public TObject
{
public:

	/// AliHLTMUONEvent class contructor.
	/// \param eventId  The event identifier number for this event. Set to -1 by default.
	AliHLTMUONEvent(AliHLTEventID_t eventId = AliHLTEventID_t(-1))
	  : TObject(), fEventId(eventId), fArray()
	{
		fArray.SetOwner(kTRUE);
	}
	
	/// Copy constructor performs a deep copy of the object.
	AliHLTMUONEvent(const AliHLTMUONEvent& event);
	
	/// The assignment operator performs a deep copy of the object.
	AliHLTMUONEvent& operator = (const AliHLTMUONEvent& event);
	
	/// Default destructor.
	virtual ~AliHLTMUONEvent() {}
	
	/// Returns the corresponding event ID for the data objects.
	AliHLTEventID_t EventID() const { return fEventId; }
	
	/// Returns the array of data objects for this event.
	const TObjArray& Array() const { return fArray; }
	
	/// Returns the array of data objects for this event.
	const TObjArray& DataObjects() const { return fArray; }
	
	/// Finds the decision object in the array of dHLT objects.
	const AliHLTMUONDecision* FindDecision() const;
	
	/**
	 * Overloaded to find the object given by name in the array of event objects.
	 * \param name  The name of the object to find. This will be the class
	 *      name if the object in the array does not overload the GetName
	 *      method.
	 * \returns  The pointer to the found object or NULL if none was found.
	 */
	virtual TObject* FindObject(const char* name) const
	{
		return fArray.FindObject(name);
	}
	
	/**
	 * Overloaded to find the object for which obj->IsEqual() is true in
	 * the array of event objects.
	 * \param obj  The object to compare to.
	 * \returns  The pointer to the found object or NULL if none was found.
	 */
	virtual TObject* FindObject(const TObject* obj) const
	{
		return fArray.FindObject(obj);
	}
	
	/// Adds an object to the event.
	/// \note This method takes ownership of the object.
	void Add(TObject* obj) { fArray.Add(obj); }
	
	/// Inherited method for printing information about all objects in the event.
	virtual void Print(Option_t* option = NULL) const;
	
	/// Inherited method for clearing the event.
	virtual void Clear(Option_t* option = "");
	
	/// Inherited method for deep copying the event.
	virtual void Copy(TObject& object) const;

private:

	AliHLTEventID_t fEventId;  ///< The event ID.
	TObjArray fArray;          ///< Array of event objects.
	
	ClassDef(AliHLTMUONEvent, 3); // Container class for dHLT event results.
};

#endif // ALIHLTMUONEVENT_H
