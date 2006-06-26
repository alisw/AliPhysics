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

Revision 1.2  2006/06/13 11:19:48  hristov

Coding conventions (Alberto)



Revision 1.1  2006/06/02 14:14:36  hristov

Separate library for CDB (Jan)



Revision 1.2  2006/03/07 07:52:34  hristov

New version (B.Yordanov)



Revision 1.2  2005/11/17 14:43:23  byordano

import to local CVS



Revision 1.1.1.1  2005/10/28 07:33:58  hristov

Initial import as subdirectory in AliRoot



Revision 1.1.1.1  2005/09/12 22:11:40  byordano

SHUTTLE package



Revision 1.2  2005/08/30 10:53:23  byordano

some more descriptions added



*/



//

// This class is a simple wrapper of

// all primitive types used in PVSS SCADA system.

// Usage examples: AliSimpleValue sim(20) -> Int type

// simfl.SetFloat(3.5) -> Float type

// Char_t chararr[4] = {'c','i','a','o'};

// AliSimpleValue simchar(4,chararr); -> DynChar type

//





#include "AliSimpleValue.h"



#include "AliLog.h"

#include <TClass.h>



//______________________________________________________________________

TObject* AliSimpleValue::AliBoolHolder::Clone(const char* /*name*/) const {

// Clone a value



	return new AliBoolHolder(fValue);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliBoolHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj

	

	if (this == obj) {

		return kTRUE;

	}



	if (AliBoolHolder::Class() != obj->IsA()) {

		return kFALSE;

	}



	return fValue == ((const AliBoolHolder*) obj)->GetValue();

}



//______________________________________________________________________

TObject* AliSimpleValue::AliByteHolder::Clone(const char* /*name*/) const {

// Clone a value



        return new AliByteHolder(fValue);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliByteHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliByteHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        return fValue == ((const AliByteHolder*) obj)->GetValue();

}



//______________________________________________________________________

TObject* AliSimpleValue::AliIntHolder::Clone(const char* /*name*/) const {

// Clone a value

	return new AliIntHolder(fValue);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliIntHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliIntHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        return fValue == ((const AliIntHolder*) obj)->GetValue();

}



//______________________________________________________________________

TObject* AliSimpleValue::AliUIntHolder::Clone(const char* /*name*/) const {

// Clone a value

        return new AliUIntHolder(fValue);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliUIntHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliUIntHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        return fValue == ((const AliUIntHolder*) obj)->GetValue();

}



//______________________________________________________________________

TObject* AliSimpleValue::AliFloatHolder::Clone(const char* /*name*/) const {

 // Clone a value

       return new AliFloatHolder(fValue);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliFloatHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliFloatHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        return fValue == ((const AliFloatHolder*) obj)->GetValue();

}



//______________________________________________________________________

TObject* AliSimpleValue::AliDynBoolHolder::Clone(const char* /*name*/) const {

 // Clone a value

       return new AliDynBoolHolder(fSize, fValues);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliDynBoolHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliDynBoolHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



	const AliDynBoolHolder* other = ((const AliDynBoolHolder*) obj);

	

	if (fSize != other->GetSize()) {

		return kFALSE;

	}



        return !memcmp(fValues, other->GetValues(), fSize * sizeof(Bool_t));

}



//______________________________________________________________________

TObject* AliSimpleValue::AliDynByteHolder::Clone(const char* /*name*/) const {

// Clone a value

        return new AliDynByteHolder(fSize, fValues);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliDynByteHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliDynByteHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        const AliDynByteHolder* other = ((const AliDynByteHolder*) obj);



        if (fSize != other->GetSize()) {

                return kFALSE;

        }



        return !memcmp(fValues, other->GetValues(), fSize * sizeof(Char_t));

}



//______________________________________________________________________

TObject* AliSimpleValue::AliDynIntHolder::Clone(const char* /*name*/) const {

// Clone a value

        return new AliDynIntHolder(fSize, fValues);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliDynIntHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



       if (this == obj) {

                return kTRUE;

        }



        if (AliDynIntHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        const AliDynIntHolder* other = ((const AliDynIntHolder*) obj);



        if (fSize != other->GetSize()) {

                return kFALSE;

        }



        return !memcmp(fValues, other->GetValues(), fSize * sizeof(Int_t));

}



//______________________________________________________________________

TObject* AliSimpleValue::AliDynUIntHolder::Clone(const char* /*name*/) const {

// Clone a value

        return new AliDynUIntHolder(fSize, fValues);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliDynUIntHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliDynUIntHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        const AliDynUIntHolder* other = ((const AliDynUIntHolder*) obj);



        if (fSize != other->GetSize()) {

                return kFALSE;

        }



        return !memcmp(fValues, other->GetValues(), fSize * sizeof(UInt_t));

}



//______________________________________________________________________

TObject* AliSimpleValue::AliDynFloatHolder::Clone(const char* /*name*/) const {

        return new AliDynFloatHolder(fSize, fValues);

}



//______________________________________________________________________

Bool_t AliSimpleValue::AliDynFloatHolder::IsEqual(const TObject* obj) const {

// check whether this is equal to obj



        if (this == obj) {

                return kTRUE;

        }



        if (AliDynFloatHolder::Class() != obj->IsA()) {

                return kFALSE;

        }



        const AliDynFloatHolder* other = ((const AliDynFloatHolder*) obj);



        if (fSize != other->GetSize()) {

                return kFALSE;

        }



        return !memcmp(fValues, other->GetValues(), fSize * sizeof(Float_t));

}



//______________________________________________________________________

//______________________________________________________________________



ClassImp(AliSimpleValue)



//______________________________________________________________________

AliSimpleValue::AliSimpleValue():

	fHolder(NULL), fType(kInvalid)

{

// empty constructor

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(const AliSimpleValue& other):

	TObject(other), fHolder(NULL), fType(other.fType)

{

// copy contructor

	if (other.fHolder) {

		fHolder = other.fHolder->Clone();

	}

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(AliSimpleValue::Type type, Int_t size):

	fHolder(NULL), fType(type)

{

// constructor

	

	switch (type) {

		case kBool:

			fHolder = new AliBoolHolder();

			break;

		case kByte:

			fHolder = new AliByteHolder();

			break;

		case kInt:

			fHolder = new AliIntHolder();

			break;

		case kUInt:

			fHolder = new AliUIntHolder();

			break;

		case kFloat:

			fHolder = new AliFloatHolder();

			break;

		case kDynBool:

			fHolder = new AliDynBoolHolder(size);

			break;

		case kDynByte:

			fHolder = new AliDynByteHolder(size);

			break;

		case kDynInt:

			fHolder = new AliDynIntHolder(size);

			break;

		case kDynUInt:

			fHolder = new AliDynUIntHolder(size);

			break;

		case kDynFloat:

			fHolder = new AliDynFloatHolder(size);

			break;

		default:	

			break;

	}

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Bool_t val) {

// contructor



	fType = kBool;

	fHolder = new AliBoolHolder(val);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Char_t val) {

// contructor

	

	fType = kByte;

	fHolder = new AliByteHolder(val);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Int_t val) {

// contructor

	

	fType = kInt;

	fHolder = new AliIntHolder(val);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(UInt_t val) {

// contructor

	

	fType = kUInt;

	fHolder = new AliUIntHolder(val);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Float_t val) {

// contructor



	fType = kFloat;

	fHolder = new AliFloatHolder(val);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Int_t size, const Bool_t* buf) {

// contructor

	

	fType = kDynBool;

	fHolder = new AliDynBoolHolder(size, buf);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Int_t size, const Char_t* buf) {

// contructor



        fType = kDynByte;

        fHolder = new AliDynByteHolder(size, buf);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Int_t size, const Int_t* buf) {

// contructor



        fType = kDynInt;

        fHolder = new AliDynIntHolder(size, buf);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Int_t size, const UInt_t* buf) {

// contructor



        fType = kDynUInt;

        fHolder = new AliDynUIntHolder(size, buf);

}



//______________________________________________________________________

AliSimpleValue::AliSimpleValue(Int_t size, const Float_t* buf) {

// contructor



        fType = kDynFloat;

        fHolder = new AliDynFloatHolder(size, buf);

}



//______________________________________________________________________

AliSimpleValue::~AliSimpleValue() {

// destructor

	

	if (fHolder) {

		delete fHolder;

	}

}



//______________________________________________________________________

AliSimpleValue& AliSimpleValue::operator=(const AliSimpleValue& other) {

// assignment op



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



//______________________________________________________________________

Bool_t AliSimpleValue::operator==(const AliSimpleValue& other) const {

// equality op



	if (fType != other.fType) {

		return kFALSE;

	}



	if (!(fHolder && other.fHolder)) {

		return kFALSE;

	}



	return fHolder->IsEqual(other.fHolder);

}



//______________________________________________________________________

void AliSimpleValue::SetBool(Bool_t val) {

// set value



	if (!TypeOk(kBool)) {

		return;

	}



	((AliBoolHolder*) fHolder)->SetValue(val);

}



//______________________________________________________________________

void AliSimpleValue::SetByte(Char_t val) {

// set value



	if (!TypeOk(kByte)) {

		return;

	}



	((AliByteHolder*) fHolder)->SetValue(val);

}



//______________________________________________________________________

void AliSimpleValue::SetInt(Int_t val) {

// set value



	if (!TypeOk(kInt)) {

		return;

	}



	((AliIntHolder*) fHolder)->SetValue(val);

}



//______________________________________________________________________

void AliSimpleValue::SetUInt(UInt_t val) {

// set value



	if (!TypeOk(kUInt)) {

		return;

	}

	

	((AliUIntHolder*) fHolder)->SetValue(val);

}



//______________________________________________________________________

void AliSimpleValue::SetFloat(Float_t val) {

// set value



	if (!TypeOk(kFloat)) {

		return;

	}



	((AliFloatHolder*) fHolder)->SetValue(val);

}



//______________________________________________________________________

Bool_t AliSimpleValue::GetBool() const {

// get value



	if (!TypeOk(kBool)) {

		return kFALSE;

	}	



	return ((AliBoolHolder*) fHolder)->GetValue();

}



//______________________________________________________________________

Char_t AliSimpleValue::GetByte() const {

// get value



	if (!TypeOk(kByte)) {

		return 0;

	}

	

        return ((AliByteHolder*) fHolder)->GetValue();

}



//______________________________________________________________________

Int_t AliSimpleValue::GetInt() const {

// get value



	if (!TypeOk(kInt)) {

		return 0;

	}

        return ((AliIntHolder*) fHolder)->GetValue();

}



//______________________________________________________________________

UInt_t AliSimpleValue::GetUInt() const {

// get value



	if (!TypeOk(kUInt)) {

		return 0;

        }



        return ((AliUIntHolder*) fHolder)->GetValue();

}



//______________________________________________________________________

Float_t AliSimpleValue::GetFloat() const {

// get value



	if (!TypeOk(kFloat)) {

		return 0;

	}



        return ((AliFloatHolder*) fHolder)->GetValue();

}



//______________________________________________________________________

void AliSimpleValue::SetDynBool(Int_t n, Bool_t val) {

// set dyn value

	

	if (!TypeOk(kDynBool)) {

		return;

	}



	if (!BoundsOk(n)) {

		return;

	}



	((AliDynBoolHolder*) fHolder)->SetValue(n,val);

}



//______________________________________________________________________

void AliSimpleValue::SetDynByte(Int_t n, Char_t val) {

// set dyn value



	if (!TypeOk(kDynByte)) {

                return;

        }



        if (!BoundsOk(n)) {

                return;

        }



        ((AliDynByteHolder*) fHolder)->SetValue(n,val);

}



//______________________________________________________________________

void AliSimpleValue::SetDynInt(Int_t n, Int_t val) {

// set dyn value



        if (!TypeOk(kDynInt)) {

                return;

        }



        if (!BoundsOk(n)) {

                return;

        }



        ((AliDynIntHolder*) fHolder)->SetValue(n,val);

}



//______________________________________________________________________

void AliSimpleValue::SetDynUInt(Int_t n, UInt_t val) {

// set dyn value



        if (!TypeOk(kDynUInt)) {

                return;

        }



        if (!BoundsOk(n)) {

                return;

        }



        ((AliDynUIntHolder*) fHolder)->SetValue(n,val);

}



//______________________________________________________________________

void AliSimpleValue::SetDynFloat(Int_t n, Float_t val) {

// set dyn value



        if (!TypeOk(kDynFloat)) {

                return;

        }



        if (!BoundsOk(n)) {

                return;

        }



        ((AliDynFloatHolder*) fHolder)->SetValue(n,val);

}



//______________________________________________________________________

Bool_t AliSimpleValue::GetDynBool(Int_t n) const {

// get dyn value



	if (!TypeOk(kDynBool)) {

                return kFALSE;

        }



        if (!BoundsOk(n)) {

                return kFALSE;

        }



	return ((AliDynBoolHolder*) fHolder)->GetValue(n);

}



//______________________________________________________________________

Char_t AliSimpleValue::GetDynByte(Int_t n) const {

// get dyn value



        if (!TypeOk(kDynByte)) {

                return 0;

        }



        if (!BoundsOk(n)) {

                return 0;

        }



        return ((AliDynByteHolder*) fHolder)->GetValue(n);

}



//______________________________________________________________________

Int_t AliSimpleValue::GetDynInt(Int_t n) const {

// get dyn value



        if (!TypeOk(kDynInt)) {

                return 0;

        }



        if (!BoundsOk(n)) {

                return 0;

        }



        return ((AliDynIntHolder*) fHolder)->GetValue(n);

}



//______________________________________________________________________

UInt_t AliSimpleValue::GetDynUInt(Int_t n) const {

// get dyn value



        if (!TypeOk(kDynUInt)) {

                return 0;

        }



        if (!BoundsOk(n)) {

                return 0;

        }



        return ((AliDynUIntHolder*) fHolder)->GetValue(n);

}



//______________________________________________________________________

Float_t AliSimpleValue::GetDynFloat(Int_t n) const {

// get dyn value



        if (!TypeOk(kDynFloat)) {

                return 0;

        }



        if (!BoundsOk(n)) {

                return 0;

        }



        return ((AliDynFloatHolder*) fHolder)->GetValue(n);

}



//______________________________________________________________________

Bool_t AliSimpleValue::TypeOk(AliSimpleValue::Type type) const {

// check that AliSimpleValue is of type type



	if (fType != type) {

		AliError(Form("SimpleValue type is not %s!", 

			GetTypeString(type)));

		return kFALSE;

	}



	return kTRUE;

}



//______________________________________________________________________

Bool_t AliSimpleValue::BoundsOk(Int_t n) const {

// Check that n is within bounds of dyn value



	switch (fType) {

		case kDynBool:

		case kDynByte:

		case kDynInt:

		case kDynUInt:

		case kDynFloat: {

			Int_t size = ((AliDynHolder*) fHolder)->GetSize();

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



//______________________________________________________________________

Int_t AliSimpleValue::GetDynamicSize() const {

	//

	// returns the size of dynamic type or 0 in case of 

	// none dynamic type.

	//



	if (!fHolder) {

		return 0;

	}



	if (!fHolder->IsA()->InheritsFrom(AliDynHolder::Class())) {

		return 0;

	}



	return ((AliDynHolder*) fHolder)->GetSize();

}



//______________________________________________________________________

TString AliSimpleValue::ToString() const {

// Print value

	

	TString result;

	

	result += "Type: ";

	result += GetTypeString(fType);

	

	result += ", Value: ";

	switch (fType) {

		case kBool:

			result += GetBool();

			break;

		case kByte:

			result += GetByte();

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



//______________________________________________________________________

Bool_t AliSimpleValue::IsDynamic(AliSimpleValue::Type type) {

// check that type is dynamic



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



//______________________________________________________________________

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



//______________________________________________________________________

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



//______________________________________________________________________

const char* AliSimpleValue::GetTypeString(AliSimpleValue::Type type) {

// return type name correspondyng to type



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

