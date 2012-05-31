// $Id$
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

///  @file   AliHLTScalars.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   28 Sep 2010
///  @brief  Implementation of the HLT scalars class.
///
/// This implements the scalars class for the HLT, which is a collection of
/// named scalar values. Searching for a named scalar is optimised by using
/// using a hash map.

#include "AliHLTScalars.h"
#include "TString.h"
#include "AliHLTLogging.h" // HLT logging
#include "TClass.h" // for Class_Name macro for logging
#include "TIterator.h"
#include "Riostream.h"
#include <cassert>

ClassImp(AliHLTScalars);
ClassImp(AliHLTScalars::AliScalar);


AliHLTScalars::AliHLTScalars() :
	TObject(),
	fScalars(AliHLTScalars::AliScalar::Class(), 128),
	fMap(TCollection::kInitHashTableCapacity, 2)
{
	// Default constructor.
	
	fMap.SetOwner(kFALSE);
}


AliHLTScalars::AliHLTScalars(const AliHLTScalars& obj) :
	TObject(obj),
	fScalars(obj.fScalars),
	fMap(obj.fMap.GetSize(), obj.fMap.GetRehashLevel())
{
	// Copy constructor performs a deep copy.

	fMap.SetOwner(kFALSE);
	fMap.AddAll(&fScalars);
}


AliHLTScalars::AliHLTScalars(const TClass* cl, Int_t initSize) :
	TObject(),
	fScalars(cl, initSize),
	fMap(TCollection::kInitHashTableCapacity, 2)
{
	// Constructor to be able to specify a custom class for the fScalars TClonesArray object.

	fMap.SetOwner(kFALSE);
}


AliHLTScalars::~AliHLTScalars()
{
	// Default destructor.
	
	Clear();
}


AliHLTScalars::AliScalar* AliHLTScalars::NewScalar(UInt_t i, const char* name, const char* description, Double_t value)
{
	// Creates a new scalar object.
	
	return new (fScalars[i]) AliScalar(name, description, value);
}


const AliHLTScalars::AliScalar& AliHLTScalars::Sentinel() const
{
	// Returns an empty sentinel object.
	
	static AliHLTScalars::AliScalar sentinel;
	return sentinel;
}


bool AliHLTScalars::Add(AliScalar*& scalar, const char* name, const char* description, Double_t value)
{
	// Adds a new scalar (internal version).
	
	scalar = static_cast<AliScalar*>( fMap.FindObject(name) );
	bool exists = scalar != NULL;
	if (not exists)
	{
		scalar = NewScalar(fScalars.GetEntriesFast(), name, description, value);
		fMap.Add(scalar);
	}
	else
	{
		scalar->Value(value);
	}
	return exists;
}


bool AliHLTScalars::Add(const char* name, const char* description, Double_t value)
{
	// Adds a new scalar.

	AliScalar* scalar = NULL;
	return Add(scalar, name, description, value);
}


bool AliHLTScalars::Remove(const char* name)
{
	// Removes a scalar from the list.

	TNamed x(name, "");
	TObject* scalar = fMap.Remove(&x);
	bool existed = scalar != NULL;
	if (existed)
	{
		fScalars.Remove(scalar);
		fScalars.Compress();
	}
	return existed;
}


const AliHLTScalars::AliScalar& AliHLTScalars::GetScalar(const char* name) const
{
	// Fetch the named scalar object.

	const TObject* scalar = fMap.FindObject(name);
	if (scalar != NULL)
	{
		return *static_cast<const AliScalar*>(scalar);
	}
	else
	{
		return Sentinel();
	}
}


AliHLTScalars::AliScalar& AliHLTScalars::GetScalar(const char* name)
{
	// Fetch the named scalar object for editing.

	TObject* scalar = fMap.FindObject(name);
	if (scalar == NULL)
	{
		scalar = NewScalar(fScalars.GetEntriesFast(), name, "", 0);
		fMap.Add(scalar);
	}
	return *static_cast<AliScalar*>(scalar);
}


const AliHLTScalars::AliScalar& AliHLTScalars::GetScalarN(UInt_t n) const
{
	// Fetch the n'th scalar object.

	if (n < NumberOfScalars())
	{
		const TObject* scalar = fScalars.UncheckedAt(Int_t(n));
		return *static_cast<const AliScalar*>(scalar);
	}
	else
	{
		return Sentinel();
	}
}


AliHLTScalars::AliScalar& AliHLTScalars::GetScalarN(UInt_t n)
{
	// Fetch the n'th scalar object for editing.

	TObject* scalar = NULL;
	if (n < NumberOfScalars())
	{
		scalar = fScalars.UncheckedAt(Int_t(n));
	}
	else
	{
		// We have to create all missing scalars since there cannot
		// be gaps in the TClonesArray. This can cause segfaults during
		// ROOT I/O of the class otherwise.
		for (UInt_t i = NumberOfScalars(); i <= n; ++i)
		{
			// Make sure the the name of the scalar is not taken.
			// If it is then find an unused name.
			TString nameToUse = Form("Scalar%d", i);
			if (FindObject(nameToUse.Data()) != NULL)
			{
				UInt_t m = 0;
				do
				{
					nameToUse = Form("Scalar%d_%d", i, m++);
				}
				while (FindObject(nameToUse.Data()) != NULL);
			}
			scalar = NewScalar(i, nameToUse.Data(), "", 0);
			fMap.Add(scalar);
		}
	}
	return *static_cast<AliScalar*>(scalar);
}


void AliHLTScalars::Reset()
{
	// Sets all the scalar values to zero.
	
	for (Int_t i = 0; i < fScalars.GetEntriesFast(); ++i)
	{
		AliScalar* scalar = static_cast<AliScalar*>( fScalars.UncheckedAt(i) );
		scalar->Value(0);
	}
}


void AliHLTScalars::Clear(Option_t* option)
{
	// Clears the array of scalars.

	fMap.Clear();
	fScalars.Delete(option);
}


void AliHLTScalars::Copy(TObject& object) const
{
	// Performs a deep copy.
	
	if (object.IsA() != AliHLTScalars::Class())
	{
		AliHLTLogging log;
		log.LoggingVarargs(kHLTLogError, Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , "Cannot copy to an object of type '%s'.", object.ClassName());
		return;
	}
	AliHLTScalars* obj = static_cast<AliHLTScalars*>(&object);
	obj->operator = (*this);
}


void AliHLTScalars::Print(Option_t* option) const
{
	// Prints the HLT trigger scalars.

	TIter next(&fScalars);
	const AliScalar* scalar = NULL;

	TString opt = option;
	if (opt == "compact")
	{
		scalar = static_cast<const AliScalar*>(next());
		if (scalar != NULL) cout << scalar->Value();
		while ((scalar = static_cast<const AliScalar*>(next())) != NULL)
		{
			cout << ", " << scalar->Value();
		}
		cout << endl;
		return;
	}
	
	// Calculate the maximum field width required to keep things aligned.
	int fieldwidth = 0;
	while ((scalar = static_cast<const AliScalar*>(next())) != NULL)
	{
		int length = strlen(scalar->Name()) + strlen(scalar->Description()) + 3;
		if (length > fieldwidth) fieldwidth = length;
	}
	if (fieldwidth > 80) fieldwidth = 80;
	
	cout << "HLT scalars:" << endl;
	next.Reset();
	while ((scalar = static_cast<const AliScalar*>(next())) != NULL)
	{
		TString str = scalar->Description();
		str += " (";
		str += scalar->Name();
		str += ")";
		cout << setw(fieldwidth) << str.Data() << setw(0)
		     << " = " << scalar->Value() << setw(0)
		     << endl;
	}
	if (fScalars.GetEntriesFast() == 0) cout << "(none)" << endl;
}


AliHLTScalars& AliHLTScalars::operator = (const AliHLTScalars& obj)
{
	// Performs a deep copy.
	
	if (this == &obj) return *this;
	Clear();  // Remove existing scalars.
	TObject::operator = (obj);
	for (Int_t i = 0; i < obj.fScalars.GetEntriesFast(); ++i)
	{
		const AliScalar* x = static_cast<const AliScalar*>(obj.fScalars.UncheckedAt(i));
		new (fScalars[i]) AliScalar(*x);
	}
	fMap.AddAll(&fScalars);
	return *this;
}


Bool_t AliHLTScalars::IsEqual(const TObject *obj) const
{
	// Checks if two sets of scalar lists have the same scalars.

	assert(obj != NULL);
	if (obj->IsA()->GetBaseClass(AliHLTScalars::Class()) == NULL)
	{
		AliHLTLogging log;
		log.LoggingVarargs(kHLTLogError, Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , "Cannot compare object of type '%s'' with an object of type '%s'.",
			      this->ClassName(), obj->ClassName()
			);
		return kFALSE;
	}
	const AliHLTScalars* rhs = static_cast<const AliHLTScalars*>(obj);
	if (fScalars.GetEntriesFast() != rhs->fScalars.GetEntriesFast()) return kFALSE;
	TIter next(&fScalars);
	const AliScalar* scalar = NULL;
	while ((scalar = static_cast<const AliScalar*>(next())) != NULL)
	{
		if (rhs->fScalars.FindObject(scalar->Name()) == NULL) return kFALSE;
	}
	return kTRUE;
}


bool AliHLTScalars::operator == (const AliHLTScalars& obj) const
{
	// Compares two scalar lists to see that they have the same scalar values.

	if (fScalars.GetEntriesFast() != obj.fScalars.GetEntriesFast()) return false;
	TIter next(&fScalars);
	const AliScalar* a = NULL;
	while ((a = static_cast<const AliScalar*>(next())) != NULL)
	{
		const AliScalar* b = static_cast<const AliScalar*>(
				obj.fScalars.FindObject(a->Name())
			);
		if (b == NULL) return false;
		if (a->Value() != b->Value()) return false;
	}
	return true;
}


void AliHLTScalars::AliScalar::Copy(TObject& object) const
{
	// Performs a deep copy.
	
	if (object.IsA() != AliHLTScalars::AliScalar::Class())
	{
		AliHLTLogging log;
		log.LoggingVarargs(kHLTLogError, Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , "Cannot copy to an object of type '%s'.", object.ClassName());
		return;
	}
	AliHLTScalars::AliScalar* obj = static_cast<AliHLTScalars::AliScalar*>(&object);
	*obj = *this;
}
