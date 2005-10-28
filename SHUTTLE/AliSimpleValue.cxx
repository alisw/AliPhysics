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
Revision 1.1.1.1  2005/09/12 22:11:40  byordano
SHUTTLE package

Revision 1.2  2005/08/30 10:53:23  byordano
some more descriptions added

*/

//
// This class is a simple wrapper of
// all primitive types used in PVSS SCADA system.
//


#include "AliSimpleValue.h"

#include "AliLog.h"
#include <TClass.h>

TObject* AliSimpleValue::BoolHolder::Clone(const char* /*name*/) const {
	return new BoolHolder(fValue);
}

Bool_t AliSimpleValue::BoolHolder::IsEqual(const TObject* obj) const {
	
	if (this == obj) {
		return kTRUE;
	}

	if (BoolHolder::Class() != obj->IsA()) {
		return kFALSE;
	}

	return fValue == ((const BoolHolder*) obj)->fValue;
}

TObject* AliSimpleValue::ByteHolder::Clone(const char* /*name*/) const {
        return new ByteHolder(fValue);
}

Bool_t AliSimpleValue::ByteHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (ByteHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        return fValue == ((const ByteHolder*) obj)->fValue;
}

TObject* AliSimpleValue::IntHolder::Clone(const char* /*name*/) const {
	return new IntHolder(fValue);
}

Bool_t AliSimpleValue::IntHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (IntHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        return fValue == ((const IntHolder*) obj)->fValue;
}

TObject* AliSimpleValue::UIntHolder::Clone(const char* /*name*/) const {
        return new UIntHolder(fValue);
}

Bool_t AliSimpleValue::UIntHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (UIntHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        return fValue == ((const UIntHolder*) obj)->fValue;
}

TObject* AliSimpleValue::FloatHolder::Clone(const char* /*name*/) const {
        return new FloatHolder(fValue);
}

Bool_t AliSimpleValue::FloatHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (FloatHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        return fValue == ((const FloatHolder*) obj)->fValue;
}

TObject* AliSimpleValue::DynBoolHolder::Clone(const char* /*name*/) const {
        return new DynBoolHolder(fSize, fValues);
}

Bool_t AliSimpleValue::DynBoolHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (DynBoolHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

	const DynBoolHolder* other = ((const DynBoolHolder*) obj);
	
	if (fSize != other->fSize) {
		return kFALSE;
	}

        return !memcmp(fValues, other->fValues, fSize * sizeof(Bool_t));
}

TObject* AliSimpleValue::DynByteHolder::Clone(const char* /*name*/) const {
        return new DynByteHolder(fSize, fValues);
}

Bool_t AliSimpleValue::DynByteHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (DynByteHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        const DynByteHolder* other = ((const DynByteHolder*) obj);

        if (fSize != other->fSize) {
                return kFALSE;
        }

        return !memcmp(fValues, other->fValues, fSize * sizeof(Char_t));
}

TObject* AliSimpleValue::DynIntHolder::Clone(const char* /*name*/) const {
        return new DynIntHolder(fSize, fValues);
}

Bool_t AliSimpleValue::DynIntHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (DynIntHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        const DynIntHolder* other = ((const DynIntHolder*) obj);

        if (fSize != other->fSize) {
                return kFALSE;
        }

        return !memcmp(fValues, other->fValues, fSize * sizeof(Int_t));
}

TObject* AliSimpleValue::DynUIntHolder::Clone(const char* /*name*/) const {
        return new DynUIntHolder(fSize, fValues);
}

Bool_t AliSimpleValue::DynUIntHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (DynUIntHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        const DynUIntHolder* other = ((const DynUIntHolder*) obj);

        if (fSize != other->fSize) {
                return kFALSE;
        }

        return !memcmp(fValues, other->fValues, fSize * sizeof(UInt_t));
}

TObject* AliSimpleValue::DynFloatHolder::Clone(const char* /*name*/) const {
        return new DynFloatHolder(fSize, fValues);
}

Bool_t AliSimpleValue::DynFloatHolder::IsEqual(const TObject* obj) const {

        if (this == obj) {
                return kTRUE;
        }

        if (DynFloatHolder::Class() != obj->IsA()) {
                return kFALSE;
        }

        const DynFloatHolder* other = ((const DynFloatHolder*) obj);

        if (fSize != other->fSize) {
                return kFALSE;
        }

        return !memcmp(fValues, other->fValues, fSize * sizeof(Float_t));
}

ClassImp(AliSimpleValue)

AliSimpleValue::AliSimpleValue():
	fHolder(NULL), fType(kInvalid)
{

}

AliSimpleValue::AliSimpleValue(const AliSimpleValue& other):
	TObject(other), fHolder(NULL), fType(other.fType)
{
	if (other.fHolder) {
		fHolder = other.fHolder->Clone();
	}
}

AliSimpleValue::AliSimpleValue(AliSimpleValue::Type type, Int_t size):
	fHolder(NULL), fType(type)
{
	
	switch (type) {
		case kBool:
			fHolder = new BoolHolder();
			break;
		case kByte:
			fHolder = new ByteHolder();
			break;
		case kInt:
			fHolder = new IntHolder();
			break;
		case kUInt:
			fHolder = new UIntHolder();
			break;
		case kFloat:
			fHolder = new FloatHolder();
			break;
		case kDynBool:
			fHolder = new DynBoolHolder(size);
			break;
		case kDynByte:
			fHolder = new DynByteHolder(size);
			break;
		case kDynInt:
			fHolder = new DynIntHolder(size);
			break;
		case kDynUInt:
			fHolder = new DynUIntHolder(size);
			break;
		case kDynFloat:
			fHolder = new DynFloatHolder(size);
			break;
		default:	
			break;
	}
}

AliSimpleValue::AliSimpleValue(Bool_t val) {

	fType = kBool;
	fHolder = new BoolHolder(val);
}

AliSimpleValue::AliSimpleValue(Char_t val) {
	
	fType = kByte;
	fHolder = new ByteHolder(val);
}

AliSimpleValue::AliSimpleValue(Int_t val) {
	
	fType = kInt;
	fHolder = new IntHolder(val);
}

AliSimpleValue::AliSimpleValue(UInt_t val) {
	
	fType = kUInt;
	fHolder = new UIntHolder(val);
}

AliSimpleValue::AliSimpleValue(Float_t val) {

	fType = kFloat;
	fHolder = new FloatHolder(val);
}

AliSimpleValue::AliSimpleValue(Int_t size, const Bool_t* buf) {
	
	fType = kDynBool;
	fHolder = new DynBoolHolder(size, buf);
}

AliSimpleValue::AliSimpleValue(Int_t size, const Char_t* buf) {

        fType = kDynByte;
        fHolder = new DynByteHolder(size, buf);
}

AliSimpleValue::AliSimpleValue(Int_t size, const Int_t* buf) {

        fType = kDynInt;
        fHolder = new DynIntHolder(size, buf);
}

AliSimpleValue::AliSimpleValue(Int_t size, const UInt_t* buf) {

        fType = kDynUInt;
        fHolder = new DynUIntHolder(size, buf);
}

AliSimpleValue::AliSimpleValue(Int_t size, const Float_t* buf) {

        fType = kDynFloat;
        fHolder = new DynFloatHolder(size, buf);
}

AliSimpleValue::~AliSimpleValue() {
	
	if (fHolder) {
		delete fHolder;
	}
}

AliSimpleValue& AliSimpleValue::operator=(const AliSimpleValue& other) {
	
	if (fHolder) {
		delete fHolder;
	}

	fType = other.fType;
	
	if (other.fHolder) {
		fHolder = other.fHolder->Clone();
	} else {
		fHolder = NULL;
	}

	return *this;
}

Bool_t AliSimpleValue::operator==(const AliSimpleValue& other) const {

	if (fType != other.fType) {
		return kFALSE;
	}

	if (!(fHolder && other.fHolder)) {
		return kFALSE;
	}

	return fHolder->IsEqual(other.fHolder);
}

void AliSimpleValue::SetBool(Bool_t val) {

	if (!TypeOk(kBool)) {
		return;
	}

	((BoolHolder*) fHolder)->fValue = val;
}

void AliSimpleValue::SetByte(Char_t val) {

	if (!TypeOk(kByte)) {
		return;
	}

	((ByteHolder*) fHolder)->fValue = val;
}

void AliSimpleValue::SetInt(Int_t val) {

	if (!TypeOk(kInt)) {
		return;
	}

	((IntHolder*) fHolder)->fValue = val;
}

void AliSimpleValue::SetUInt(UInt_t val) {

	if (!TypeOk(kUInt)) {
		return;
	}
	
	((UIntHolder*) fHolder)->fValue = val;
}

void AliSimpleValue::SetFloat(Float_t val) {

	if (!TypeOk(kFloat)) {
		return;
	}

	((FloatHolder*) fHolder)->fValue = val;
}

Bool_t AliSimpleValue::GetBool() const {

	if (!TypeOk(kBool)) {
		return kFALSE;
	}	

	return ((BoolHolder*) fHolder)->fValue;
}

Char_t AliSimpleValue::GetByte() const {

	if (!TypeOk(kByte)) {
		return 0;
	}
	
        return ((ByteHolder*) fHolder)->fValue;
}

Int_t AliSimpleValue::GetInt() const {

	if (!TypeOk(kInt)) {
		return 0;
	}
        return ((IntHolder*) fHolder)->fValue;
}

UInt_t AliSimpleValue::GetUInt() const {

	if (!TypeOk(kUInt)) {
		return 0;
        }

        return ((UIntHolder*) fHolder)->fValue;
}

Float_t AliSimpleValue::GetFloat() const {

	if (!TypeOk(kFloat)) {
		return 0;
	}

        return ((FloatHolder*) fHolder)->fValue;
}

void AliSimpleValue::SetDynBool(Int_t n, Bool_t val) {
	
	if (!TypeOk(kDynBool)) {
		return;
	}

	if (!BoundsOk(n)) {
		return;
	}

	((DynBoolHolder*) fHolder)->fValues[n] = val;
}

void AliSimpleValue::SetDynByte(Int_t n, Char_t val) {

	if (!TypeOk(kDynByte)) {
                return;
        }

        if (!BoundsOk(n)) {
                return;
        }

        ((DynByteHolder*) fHolder)->fValues[n] = val;
}

void AliSimpleValue::SetDynInt(Int_t n, Int_t val) {

        if (!TypeOk(kDynInt)) {
                return;
        }

        if (!BoundsOk(n)) {
                return;
        }

        ((DynIntHolder*) fHolder)->fValues[n] = val;
}

void AliSimpleValue::SetDynUInt(Int_t n, UInt_t val) {

        if (!TypeOk(kDynUInt)) {
                return;
        }

        if (!BoundsOk(n)) {
                return;
        }

        ((DynUIntHolder*) fHolder)->fValues[n] = val;
}

void AliSimpleValue::SetDynFloat(Int_t n, Float_t val) {

        if (!TypeOk(kDynFloat)) {
                return;
        }

        if (!BoundsOk(n)) {
                return;
        }

        ((DynFloatHolder*) fHolder)->fValues[n] = val;
}

Bool_t AliSimpleValue::GetDynBool(Int_t n) const {

	if (!TypeOk(kDynBool)) {
                return kFALSE;
        }

        if (!BoundsOk(n)) {
                return kFALSE;
        }

	return ((DynBoolHolder*) fHolder)->fValues[n];
}

Char_t AliSimpleValue::GetDynByte(Int_t n) const {

        if (!TypeOk(kDynByte)) {
                return 0;
        }

        if (!BoundsOk(n)) {
                return 0;
        }

        return ((DynByteHolder*) fHolder)->fValues[n];
}

Int_t AliSimpleValue::GetDynInt(Int_t n) const {

        if (!TypeOk(kDynInt)) {
                return 0;
        }

        if (!BoundsOk(n)) {
                return 0;
        }

        return ((DynIntHolder*) fHolder)->fValues[n];
}

UInt_t AliSimpleValue::GetDynUInt(Int_t n) const {

        if (!TypeOk(kDynUInt)) {
                return 0;
        }

        if (!BoundsOk(n)) {
                return 0;
        }

        return ((DynUIntHolder*) fHolder)->fValues[n];
}

Float_t AliSimpleValue::GetDynFloat(Int_t n) const {

        if (!TypeOk(kDynFloat)) {
                return 0;
        }

        if (!BoundsOk(n)) {
                return 0;
        }

        return ((DynFloatHolder*) fHolder)->fValues[n];
}

Bool_t AliSimpleValue::TypeOk(AliSimpleValue::Type type) const {

	if (fType != type) {
		AliError(Form("SimpleValue type is not %s!", 
			GetTypeString(type)));
		return kFALSE;
	}

	return kTRUE;
}

Bool_t AliSimpleValue::BoundsOk(Int_t n) const {

	switch (fType) {
		case kDynBool:
		case kDynByte:
		case kDynInt:
		case kDynUInt:
		case kDynFloat: {
			Int_t size = ((DynHolder*) fHolder)->fSize;
			if (n < 0 || n >= size) {
				AliError(Form("Index %d out of bounds!", n));
				return kFALSE;
			}
			return kTRUE;
		}
		case kBool:
		case kByte:
		case kInt:
		case kUInt:
		case kFloat:
			AliError(Form("SimpleValue type %s is not dynamic!",
				GetTypeString(fType)));
			return kFALSE;
		default:
			AliError("Invalid or unknown type!");
			return kFALSE;
	}
}

Int_t AliSimpleValue::GetDynamicSize() const {
	//
	// returns the size of dynamic type or 0 in case of 
	// none dynamic type.
	//

	if (!fHolder) {
		return 0;
	}

	if (!fHolder->IsA()->InheritsFrom(DynHolder::Class())) {
		return 0;
	}

	return ((DynHolder*) fHolder)->fSize;
}

TString AliSimpleValue::ToString() const {
	
	TString result;
	
	result += "Type: ";
	result += GetTypeString(fType);
	
	result += ", Value: ";
	switch (fType) {
		case kBool:
			result += GetBool();
			break;
		case kByte:
			result += (Int_t) GetByte();
			break;
		case kInt:
			result += GetInt();
			break;
		case kUInt:
			result += GetUInt();
			break;
		case kFloat:
			result += GetFloat();
			break;
		case kDynBool: {
				result += "[";
				Int_t size = GetDynamicSize();
				for (Int_t k = 0; k < size; k ++) {
					result += GetDynBool(k);
					if (k + 1 < size) {
						result += ", ";
					}
				}
				result += "]";
			}
			break;
		case kDynByte: {
                                result += "[";
                                Int_t size = GetDynamicSize();
                                for (Int_t k = 0; k < size; k ++) {
                                        result += GetDynByte(k);
                                        if (k + 1 < size) {
                                                result += ", ";
                                        }
                                }
                                result += "]";
                        }
                        break;
		case kDynInt: {
                                result += "[";
                                Int_t size = GetDynamicSize();
                                for (Int_t k = 0; k < size; k ++) {
                                        result += GetDynInt(k);
                                        if (k + 1 < size) {
                                                result += ", ";
                                        }
                                }
                                result += "]";
                        }
                        break;
		case kDynUInt: {
                                result += "[";
                                Int_t size = GetDynamicSize();
                                for (Int_t k = 0; k < size; k ++) {
                                        result += GetDynUInt(k);
                                        if (k + 1 < size) {
                                                result += ", ";
                                        }
                                }
                                result += "]";
                        }
                        break;
		case kDynFloat: {
                                result += "[";
                                Int_t size = GetDynamicSize();
                                for (Int_t k = 0; k < size; k ++) {
                                        result += GetDynFloat(k);
                                        if (k + 1 < size) {
                                                result += ", ";
                                        }
                                }
                                result += "]";
                        }
                        break;
		default:
			result += "Unknown";		
	}	

	return result;
}

Bool_t AliSimpleValue::IsDynamic(AliSimpleValue::Type type) {

	 switch (type) {
                case kDynBool:
                case kDynByte:
                case kDynInt:
                case kDynUInt:
                case kDynFloat:
			return kTRUE;
		default:
			return kFALSE; 
	}
}

Int_t AliSimpleValue::GetSize() const {
	//
	// return the number of bytes used by this value.
	// In case of dynamic type it returns dynamic size multiplied
	// by the size of corresponding primitive type.
	//
	
	return IsDynamic(fType) ? 
		GetDynamicSize() * AliSimpleValue::GetPrimitiveSize(fType): 
		AliSimpleValue::GetPrimitiveSize(fType);
} 

Int_t AliSimpleValue::GetPrimitiveSize(AliSimpleValue::Type type) {
	//
	// returns the number of bytes used by particular primitive type
	// or by the corresponding primitive type in case of dynamic type.
	//


	switch (type) {

		case kBool: 
		case kDynBool: return sizeof(Bool_t);
		case kByte: 
		case kDynByte: return sizeof(Char_t);
		case kInt: 
		case kDynInt: return sizeof(Int_t);
		case kUInt: 
		case kDynUInt: return sizeof(UInt_t);
		case kFloat: 
		case kDynFloat: return sizeof(Float_t);
		default:
			return 0;
	}
} 

const char* AliSimpleValue::GetTypeString(AliSimpleValue::Type type) {

	switch (type) {
		case kBool: return "Bool";
		case kByte: return "Byte";
		case kInt: return "Int";
		case kUInt: return "UInt";
		case kFloat: return "Float";
		case kDynBool: return "DynBool";
		case kDynByte: return "DynByte";
		case kDynInt: return "DynInt";
		case kDynUInt: return "DynUInt";
		case kDynFloat: return "DynFloat";
		default:
			return "Unknown";
	}
}
