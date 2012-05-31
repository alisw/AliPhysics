// $Id: $
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

///  @file   AliHLTTriggerCounters.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   4 Sep 2010
///  @brief  Implementation of the HLT global trigger counters.
///
/// This implements the trigger counters for the global HLT.
/// These contain the total number of fired triggers by class and their current
/// estimated rate.

#include "AliHLTTriggerCounters.h"
#include "TString.h"
#include "AliLog.h"
#include "TIterator.h"
#include "Riostream.h"
#include <cassert>

ClassImp(AliHLTTriggerCounters);
ClassImp(AliHLTTriggerCounters::AliCounter);


AliHLTTriggerCounters::AliHLTTriggerCounters() :
	AliHLTScalars(AliHLTTriggerCounters::AliCounter::Class(), 128),
	fTimeStamp()
{
	// Default constructor.
}


AliHLTTriggerCounters::AliHLTTriggerCounters(const AliHLTTriggerCounters& obj) :
	AliHLTScalars(obj),
	fTimeStamp(obj.fTimeStamp)
{
	// Copy constructor performs a deep copy.
}


AliHLTTriggerCounters::~AliHLTTriggerCounters()
{
	// Default destructor.
}

AliHLTScalars::AliScalar* AliHLTTriggerCounters::NewScalar(UInt_t i, const char* name, const char* description, Double_t value)
{
	// Creates a new scalar object.
	
	return new (ScalarForConstructor(i)) AliCounter(name, description, 0, value);
}


const AliHLTScalars::AliScalar& AliHLTTriggerCounters::Sentinel() const
{
	// Returns an empty sentinel object.
	
	static AliHLTTriggerCounters::AliCounter sentinel;
	return sentinel;
}


bool AliHLTTriggerCounters::Add(const char* name, const char* description, Double_t rate, ULong64_t value)
{
	// Adds a new counter.

	AliScalar* scalar = NULL;
	bool exists = AliHLTScalars::Add(scalar, name, description, rate);
	assert(scalar != NULL);
	static_cast<AliCounter*>(scalar)->Counter(value);
	return exists;
}

void AliHLTTriggerCounters::Reset()
{
	// Sets all the counter values and rates to zero.
	
	for (UInt_t i = 0; i < NumberOfCounters(); ++i)
	{
		AliCounter* counter = static_cast<AliCounter*>( ScalarUncheckedAt(i) );
		counter->Counter(0);
		counter->Rate(0);
	}
}


void AliHLTTriggerCounters::Copy(TObject& object) const
{
	// Performs a deep copy.
	
	if (object.IsA() != AliHLTTriggerCounters::Class())
	{
		AliError(Form("Cannot copy to an object of type '%s'.", object.ClassName()));
		return;
	}
	AliHLTTriggerCounters* obj = static_cast<AliHLTTriggerCounters*>(&object);
	obj->operator = (*this);
}


void AliHLTTriggerCounters::Print(Option_t* option) const
{
	// Prints the HLT trigger counters.

	TString opt = option;
	if (opt == "compact")
	{
		if (NumberOfCounters() > 0)
		{
			AliCounter* counter = static_cast<AliCounter*>( ScalarUncheckedAt(0) );
			cout << counter->Counter();
		}
		for (UInt_t i = 1; i < NumberOfCounters(); ++i)
		{
			AliCounter* counter = static_cast<AliCounter*>( ScalarUncheckedAt(i) );
			cout << ", " << counter->Counter();
		}
		cout << endl;
		return;
	}
	
	// Calculate the maximum field width required to keep things aligned.
	int fieldwidth = 0;
	for (UInt_t i = 0; i < NumberOfCounters(); ++i)
	{
		AliCounter* counter = static_cast<AliCounter*>( ScalarUncheckedAt(i) );
		int length = strlen(counter->Name()) + strlen(counter->Description()) + 3;
		if (length > fieldwidth) fieldwidth = length;
	}
	if (fieldwidth > 80) fieldwidth = 80;
	
	cout << "HLT trigger counters (" << fTimeStamp.AsString() << "):" << endl;
	for (UInt_t i = 0; i < NumberOfCounters(); ++i)
	{
		AliCounter* counter = static_cast<AliCounter*>( ScalarUncheckedAt(i) );
		TString str = counter->Description();
		str += " (";
		str += counter->Name();
		str += ")";
		cout << setw(fieldwidth) << str.Data() << setw(0)
		     << " = " << counter->Counter() << setw(0)
		     << " (" << counter->Rate() << " Hz)"
		     << endl;
	}
	if (NumberOfCounters() == 0) cout << "(none)" << endl;
}


AliHLTTriggerCounters& AliHLTTriggerCounters::operator = (const AliHLTTriggerCounters& obj)
{
	// Performs a deep copy.
	
	if (this == &obj) return *this;
	Clear();  // Remove existing counters.
	TObject::operator = (obj);
	fTimeStamp = obj.fTimeStamp;
	for (UInt_t i = 0; i < obj.NumberOfCounters(); ++i)
	{
		AliCounter* counter = static_cast<AliCounter*>( obj.ScalarUncheckedAt(i) );
		Add(counter->Name(), counter->Description(), counter->Rate(), counter->Counter());
	}
	return *this;
}


bool AliHLTTriggerCounters::operator == (const AliHLTTriggerCounters& obj) const
{
	// Compares two counter lists to see that they have the same counter values.

	if (NumberOfCounters() != obj.NumberOfCounters()) return false;
	if (fTimeStamp != obj.fTimeStamp) return false;
	
	for (UInt_t i = 0; i < NumberOfCounters(); ++i)
	{
		const AliCounter* a = static_cast<const AliCounter*>( ScalarUncheckedAt(i) );
		const AliCounter* b = static_cast<const AliCounter*>( obj.FindObject(a->Name()) );
		if (b == NULL) return false;
		if (a->Value() != b->Value() or a->Rate() != b->Rate()) return false;
	}
	return true;
}


void AliHLTTriggerCounters::AliCounter::Copy(TObject& object) const
{
	// Performs a deep copy.
	
	if (object.IsA() != AliHLTTriggerCounters::AliCounter::Class())
	{
		AliError(Form("Cannot copy to an object of type '%s'.", object.ClassName()));
		return;
	}
	AliHLTTriggerCounters::AliCounter* obj = static_cast<AliHLTTriggerCounters::AliCounter*>(&object);
	*obj = *this;
}
