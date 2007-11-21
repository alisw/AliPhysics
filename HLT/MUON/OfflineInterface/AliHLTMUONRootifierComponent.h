#ifndef ALIHLTMUONROOTIFIERCOMPONENT_H
#define ALIHLTMUONROOTIFIERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// @file   AliHLTMUONRootifierComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   
/// @brief  Component for converting dHLT raw data into ROOT objects.
///

#include "AliHLTProcessor.h"

// Temporary solution for grouping together objects for the same event.
#include "TObjArray.h"

class AliHLTMUONEvent : public TObject
{
public:

	AliHLTMUONEvent(AliHLTEventID_t eventId = AliHLTEventID_t(-1))
		: fEventId(eventId)
	{
		fArray.SetOwner(kTRUE);
	}
	
	virtual ~AliHLTMUONEvent() {}
	
	AliHLTEventID_t EventID() const { return fEventId; }
	const TObjArray& Array() const { return fArray; }
	
	// Takes ownership of the object.
	void Add(TObject* obj) { fArray.Add(obj); }
	
	virtual void Print(Option_t* option = NULL) const;

private:

	AliHLTEventID_t fEventId;  // The event ID.
	TObjArray fArray;          // Array of event objects.
	
	ClassDef(AliHLTMUONEvent, 1); // Container class for dHLT event results.
};


/**
 * Converts dHLT raw data blocks into ROOT objects.
 */
class AliHLTMUONRootifierComponent : public AliHLTProcessor
{
public:

	AliHLTMUONRootifierComponent();
	virtual ~AliHLTMUONRootifierComponent();
	
	virtual const char* GetComponentID();

	virtual void GetInputDataTypes(vector<AliHLTComponentDataType>& list);
	virtual AliHLTComponentDataType GetOutputDataType();

	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

	virtual AliHLTComponent* Spawn();

protected:

	virtual int DoInit(int argc, const char** argv);
	virtual int DoDeinit();
	virtual int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	using AliHLTProcessor::DoEvent;
	
private:

	// Prevent copying of these objects.
	AliHLTMUONRootifierComponent(const AliHLTMUONRootifierComponent& /*object*/);
	AliHLTMUONRootifierComponent& operator = (const AliHLTMUONRootifierComponent& /*object*/);

	ClassDef(AliHLTMUONRootifierComponent, 0); // Converter component of dHLT raw data.
};

#endif // ALIHLTMUONROOTIFIERCOMPONENT_H
