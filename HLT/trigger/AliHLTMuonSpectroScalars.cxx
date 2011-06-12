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

///  @file   AliHLTMuonSpectroScalars.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   9 Nov 2009
///  @brief  Implementation of the muon spectrometer trigger scalars.
///
/// This implements the trigger scalars for the muon spectrometer that can be used
/// in the HLT global trigger and/or added to the AliESDEvent.

#include "AliHLTMuonSpectroScalars.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "Riostream.h"

ClassImp(AliHLTMuonSpectroScalars);
ClassImp(AliHLTMuonSpectroScalars::AliScalar);


AliHLTMuonSpectroScalars::AliHLTMuonSpectroScalars() :
	TObject(),
	fScalars(AliHLTMuonSpectroScalars::AliScalar::Class(), 128),
	fIndex(128),
	fIndexValid(false)
{
	// Default constructor.
	
	fIndex.SetOwner(kFALSE);
}


AliHLTMuonSpectroScalars::AliHLTMuonSpectroScalars(const AliHLTMuonSpectroScalars& obj) :
	TObject(obj),
	fScalars(obj.fScalars),
	fIndex(obj.fIndex),
	fIndexValid(obj.fIndexValid)
{
	// Copy constructor performs a deep copy.
	
	fIndex.SetOwner(kFALSE);
}


AliHLTMuonSpectroScalars::~AliHLTMuonSpectroScalars()
{
	// Default destructor.
	
	Clear();
}


void AliHLTMuonSpectroScalars::Add(const char* name, const char* description, Double_t value)
{
	// Adds a new scalar.
	
	new (fScalars[fScalars.GetEntriesFast()]) AliScalar(value, name, description);
	fIndexValid = false; // Invalidate the index.
}


bool AliHLTMuonSpectroScalars::Exists(const char* name) const
{
	// Checks if the scalar exists or not.
	
	if (not fIndexValid) MakeIndex();
	AliScalar tmpobj(0, name, "");
	Int_t n = fIndex.BinarySearch(&tmpobj);
	if (n == -1) return false;
	return true;
}


AliHLTMuonSpectroScalars::AliScalar* AliHLTMuonSpectroScalars::GetScalarN(UInt_t n)
{
	// Fetch the n'th scalar object.
	
	if (n >= NumberOfScalars())
	{
		AliError(Form("Value 'n' is out of bounds. Should be in the range [0..%d].", NumberOfScalars()-1));
		return NULL;
	}
	return static_cast<AliScalar*>( fScalars.UncheckedAt(Int_t(n)) );
}


const AliHLTMuonSpectroScalars::AliScalar* AliHLTMuonSpectroScalars::GetScalarN(UInt_t n) const
{
	// Fetch the n'th scalar object.
	
	if (n >= NumberOfScalars())
	{
		AliError(Form("Value 'n' is out of bounds. Should be in the range [0..%d].", NumberOfScalars()-1));
		return NULL;
	}
	return static_cast<const AliScalar*>( fScalars.UncheckedAt(Int_t(n)) );
}


AliHLTMuonSpectroScalars::AliScalar* AliHLTMuonSpectroScalars::GetScalar(const char* name)
{
	// Fetch the named scalar object.
	
	if (not fIndexValid) MakeIndex();
	AliScalar tmpobj(0, name, "");
	Int_t n = fIndex.BinarySearch(&tmpobj);
	if (n == -1)
	{
		AliError(Form("Scalar '%s' could not be found.", name));
		return NULL;
	}
	return static_cast<AliScalar*>( fIndex.UncheckedAt(n) );
}


const AliHLTMuonSpectroScalars::AliScalar* AliHLTMuonSpectroScalars::GetScalar(const char* name) const
{
	// Fetch the named scalar object.
	
	if (not fIndexValid) MakeIndex();
	AliScalar tmpobj(0, name, "");
	Int_t n = fIndex.BinarySearch(&tmpobj);
	if (n == -1)
	{
		AliError(Form("Scalar '%s' could not be found.", name));
		return NULL;
	}
	return static_cast<const AliScalar*>( fIndex.UncheckedAt(n) );
}


Double_t AliHLTMuonSpectroScalars::GetN(UInt_t n) const
{
	// Fetches the n'th scalar value.
	
	const AliScalar* scalar = GetScalarN(n);
	if (scalar == NULL) return 0;
	return scalar->Value();
}


Double_t AliHLTMuonSpectroScalars::Get(const char* name) const
{
	// Fetches the n'th scalar value.
	
	const AliScalar* scalar = GetScalar(name);
	if (scalar == NULL) return 0;
	return scalar->Value();
}


bool AliHLTMuonSpectroScalars::SetN(UInt_t n, Double_t value)
{
	// Sets the n'th scalar value.
	
	AliScalar* scalar = GetScalarN(n);
	if (scalar == NULL) return false;
	scalar->Value(value);
	return true;
}


bool AliHLTMuonSpectroScalars::Set(const char* name, Double_t value)
{
	// Sets the named scalar value.
	
	AliScalar* scalar = GetScalar(name);
	if (scalar == NULL) return false;
	scalar->Value(value);
	return true;
}


bool AliHLTMuonSpectroScalars::IncrementN(UInt_t n, UInt_t count)
{
	// Increments the n'th scalar by a value of 'count'.
	
	AliScalar* scalar = GetScalarN(n);
	if (scalar == NULL) return false;
	scalar->Increment(count);
	return true;
}


bool AliHLTMuonSpectroScalars::Increment(const char* name, UInt_t count)
{
	// Increments the named scalar by a value of 'count'.
	
	AliScalar* scalar = GetScalar(name);
	if (scalar == NULL) return false;
	scalar->Increment(count);
	return true;
}


const char* AliHLTMuonSpectroScalars::Name(UInt_t n) const
{
	// Fetches the n'th scalar's name.
	
	const AliScalar* scalar = GetScalarN(n);
	if (scalar == NULL) return NULL;
	return scalar->Name();
}


const char* AliHLTMuonSpectroScalars::Description(UInt_t n) const
{
	// Fetches the n'th scalar's description.
	
	const AliScalar* scalar = GetScalarN(n);
	if (scalar == NULL) return NULL;
	return scalar->Description();
}


void AliHLTMuonSpectroScalars::Reset()
{
	// Sets all the scalar values to zero.
	
	for (Int_t i = 0; i < fScalars.GetEntriesFast(); ++i)
	{
		AliScalar* scalar = static_cast<AliScalar*>( fScalars.UncheckedAt(i) );
		scalar->Value(0);
	}
}


void AliHLTMuonSpectroScalars::Clear(Option_t* option)
{
	// Clears the array of scalars.
	
	fScalars.Delete(option);
	fIndex.Clear();
}


void AliHLTMuonSpectroScalars::Copy(TObject& object) const
{
	// Performs a deep copy.
	
	if (object.IsA() != AliHLTMuonSpectroScalars::Class())
	{
		AliError(Form("Cannot copy to an object of type '%s'.", object.ClassName()));
		return;
	}
	AliHLTMuonSpectroScalars* obj = static_cast<AliHLTMuonSpectroScalars*>(&object);
	*obj = *this;
}


TObject* AliHLTMuonSpectroScalars::FindObject(const char* name) const
{
	// Finds the scalar object by name.
	
	if (fIndexValid)
	{
		AliScalar tmpobj(0, name, "");
		Int_t n = fIndex.BinarySearch(&tmpobj);
		if (n != -1) return fIndex.UncheckedAt(n);
	}
	else
	{
		return fScalars.FindObject(name);
	}
	return NULL;
}


TObject* AliHLTMuonSpectroScalars::FindObject(const TObject* obj) const
{
	// Finds the scalar object with the same name as obj->GetName().
	
	return FindObject(obj->GetName());
}


void AliHLTMuonSpectroScalars::Print(Option_t* option) const
{
	// Prints the muon spectrometer's HLT trigger scalars.
	
	TString opt = option;
	if (opt == "compact")
	{
		if (NumberOfScalars() > 0)
		cout << GetN(0);
		for (UInt_t i = 1; i < NumberOfScalars(); ++i) cout << ", " << GetN(i);
		cout << endl;
		return;
	}
	
	// Calculate the maximum field width required to keep things aligned.
	int fieldwidth = 0;
	for (Int_t i = 0; i < fScalars.GetEntriesFast(); ++i)
	{
		AliScalar* scalar = static_cast<AliScalar*>( fScalars.UncheckedAt(i) );
		int length = strlen(scalar->Description());
		if (length > fieldwidth) fieldwidth = length;
	}
	if (fieldwidth > 80) fieldwidth = 80;
	
	cout << "HLT muon spectrometer trigger scalars:" << endl;
	for (Int_t i = 0; i < fScalars.GetEntriesFast(); ++i)
	{
		AliScalar* scalar = static_cast<AliScalar*>( fScalars.UncheckedAt(i) );
		cout << setw(fieldwidth) << scalar->Description() << setw(0)
		     << " = " << scalar->Value() << endl;
	}
	if (fScalars.GetEntriesFast() == 0) cout << "(none)" << endl;
}


AliHLTMuonSpectroScalars& AliHLTMuonSpectroScalars::operator = (const AliHLTMuonSpectroScalars& obj)
{
	// Performs a deep copy.
	
	if (this == &obj) return *this;
	fScalars.Delete();
	for (Int_t i = 0; i < obj.fScalars.GetEntriesFast(); ++i)
	{
		AliScalar* scalar = static_cast<AliScalar*>( obj.fScalars.UncheckedAt(i) );
		Add(scalar->Name(), scalar->Description(), scalar->Value());
	}
	MakeIndex();
	return *this;
}


bool AliHLTMuonSpectroScalars::operator == (const AliHLTMuonSpectroScalars& obj) const
{
	// Compares two scalar objects.
	
	if (not fIndexValid) MakeIndex();
	if (fScalars.GetEntriesFast() != obj.fScalars.GetEntriesFast()) return false;
	for (Int_t i = 0; i < obj.fScalars.GetEntriesFast(); ++i)
	{
		AliScalar* scalar1 = static_cast<AliScalar*>( obj.fScalars.UncheckedAt(i) );
		Int_t n = fIndex.BinarySearch(scalar1);
		if (n == -1) return false;
		AliScalar* scalar2 = static_cast<AliScalar*>( fIndex.UncheckedAt(n) );
		if (scalar1->Value() != scalar2->Value()) return false;
	}
	return true;
}


void AliHLTMuonSpectroScalars::MakeIndex() const
{
	// Makes the index fIndex required for faster searching in fScalars.
	
	fIndex.Clear();
	for (Int_t i = 0; i < fScalars.GetEntriesFast(); ++i)
	{
		fIndex.Add(fScalars.UncheckedAt(i));
	}
	fIndex.Sort();
	fIndexValid = true;
}


void AliHLTMuonSpectroScalars::AliScalar::Copy(TObject& object) const
{
	// Performs a deep copy.
	
	if (object.IsA() != AliHLTMuonSpectroScalars::AliScalar::Class())
	{
		AliError(Form("Cannot copy to an object of type '%s'.", object.ClassName()));
		return;
	}
	AliHLTMuonSpectroScalars::AliScalar* obj = static_cast<AliHLTMuonSpectroScalars::AliScalar*>(&object);
	*obj = *this;
}
