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
Revision 1.6  2006/08/15 10:50:00  jgrosseo
effc++ corrections (alberto)

Revision 1.5  2006/07/20 09:54:40  jgrosseo
introducing status management: The processing per subdetector is divided into several steps,
after each step the status is stored on disk. If the system crashes in any of the steps the Shuttle
can keep track of the number of failures and skips further processing after a certain threshold is
exceeded. These thresholds can be configured in LDAP.

Revision 1.4  2006/07/04 14:59:57  jgrosseo
revision of AliDCSValue: Removed wrapper classes, reduced storage size per value by factor 2

Revision 1.3  2006/06/12 09:11:16  jgrosseo
coding conventions (Alberto)

Revision 1.2  2006/03/07 07:52:34  hristov
New version (B.Yordanov)

Revision 1.3  2005/11/17 17:47:34  byordano
TList changed to TObjArray

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
// This class is a wrapper of AliDCSMessage.
// These are the messages which form AliDCSProtocol.
// Every message has header and body. The body size is written in the header.
// There are five message types:
//	1) Request - used by the client to form a single request to DCS server
//	2) Count - returned by the server to inidicate the total number of 
//		values which would be sent to the client.
//	3) ResultSet - returned by the server and contains part of values set
//		which forms the server resposen. 
//	4) Error - returned by the server in case of error
//	5) MultiRequest - used by the client to form multi request.
//		This is a request which serves many aliases/dp at the same time
//		For all aliases/dp the same time interval is used.	
// Short description of the schema:
//	The client sends a request (Request or MultiRequest) and the server 
// 	returns:
//		1) Count - the total number of values that the client should 
//			expect.
//		2) ResultSet* - every ResultSet message contains a part
//			of valueSet (some values) which the client should expect
//			The client can wait for ResultMessage until it gets
//			all values (total number) which was returned by the
//			Count message at the beginning of the ResutlSet sereie.
//	In case of error:
//		1) Error - contains the error code and error description
//

#include "AliDCSMessage.h"

#include "AliLog.h"

#include <Bytes.h>
#include <TObjString.h>

#include <ctype.h>
#include <assert.h>

ClassImp(AliDCSMessage)

//______________________________________________________________________
AliDCSMessage::AliDCSMessage():
	fMessage(NULL), fMessageSize(0), fType(kInvalid),
	fStartTime(0), fEndTime(0),
	fRequestString(""), fCount(0),
	fValueType(AliDCSValue::kInvalid),
	fErrorCode(kNoneError), fErrorString(""),
	fRequestStrings()
{
// default constructor
	fValues = new TObjArray();
	fValues->SetOwner(0);

}

//______________________________________________________________________
AliDCSMessage::AliDCSMessage(const char* message, UInt_t size):
        fMessageSize(size), fType(kInvalid),
	fStartTime(0), fEndTime(0),
	fRequestString(""), fCount(0),
	fValueType(AliDCSValue::kInvalid),
	fErrorCode(kNoneError), fErrorString(""),
	fRequestStrings()
{
// default constructor

        fMessage = new char[size];

        memcpy(fMessage, message, size);
	fValues = new TObjArray();
	fValues->SetOwner(0);
}

//______________________________________________________________________
AliDCSMessage::AliDCSMessage(const AliDCSMessage& /*other*/):
	TObject(), fMessage(NULL), fMessageSize(0), fType(kInvalid),
	fStartTime(0), fEndTime(0),
	fRequestString(""), fCount(0),
	fValueType(AliDCSValue::kInvalid),
	fErrorCode(kNoneError), fErrorString(""),
	fRequestStrings()
{
// copy constructor (not implemented)

}

//______________________________________________________________________
AliDCSMessage &AliDCSMessage::operator=(const AliDCSMessage& /*other*/)
{
// assignment operator (not implemented)

return *this;
}

//______________________________________________________________________
AliDCSMessage::~AliDCSMessage() 
{
// destructor

	DestroyMessage();
        DestroyBuffer();
	if(fValues) delete fValues; fValues=0;
}

//______________________________________________________________________
void AliDCSMessage::CreateRequestMessage(RequestType type,
	UInt_t startTime, UInt_t endTime, const char* request)
{
// Create request message

	DestroyMessage();

	fType = AliDCSMessage::kRequest;
	fRequestType = type;
	fStartTime = startTime;
	fEndTime = endTime;
	fRequestString = request;
}

//______________________________________________________________________
void AliDCSMessage::CreateMultiRequestMessage(RequestType type, 
	UInt_t startTime, UInt_t endTime)
{
// Create multi request message

	DestroyMessage();

        fType = AliDCSMessage::kMultiRequest;
	fRequestType = type;
        fStartTime = startTime;
	fEndTime = endTime;
}
	
//______________________________________________________________________
void AliDCSMessage::CreateCountMessage(UInt_t count) 
{
// Create count request message

	DestroyMessage();

	fType = AliDCSMessage::kCount;
	fCount = count;
}

//______________________________________________________________________
void AliDCSMessage::CreateResultSetMessage(AliDCSValue::Type type)
{
  // Create result set message

	DestroyMessage();

	fType = AliDCSMessage::kResultSet;
	fValueType = type;
}

//______________________________________________________________________
void AliDCSMessage::CreateErrorMessage(ErrorCode errorCode, 
	const char* errorString)
{
// Create error message

	DestroyMessage();

	fType = AliDCSMessage::kError;
	fErrorCode = errorCode;
	fErrorString = errorString;
}

/*
void AliDCSMessage::CreateNextMessage() {
	DestroyMessage();

	fType = AliDCSMessage::kNext;
} */

//______________________________________________________________________
void AliDCSMessage::DestroyMessage() 
{
// Destroy message

	fType = kInvalid;
        ClearValues();
        ClearRequestStrings();
}

//______________________________________________________________________
void AliDCSMessage::SetBool(char* buf, Bool_t val) 
{
// Set bool value to buf 

	tobuf(buf, val);
}

//______________________________________________________________________
void AliDCSMessage::SetByte(char* buf, Char_t val) 
{
// Set byte value to buf 

	tobuf(buf, val);
}

//______________________________________________________________________
void AliDCSMessage::SetUByte(char* buf, UChar_t val)
{
// Set ubyte value to buf

	tobuf(buf, val);
}

//______________________________________________________________________
void AliDCSMessage::SetInt(char* buf, Int_t val) 
{
// Set int value to buf 

	tobuf(buf, val);
}

//______________________________________________________________________
void AliDCSMessage::SetUInt(char* buf, UInt_t val) 
{
// Set uint value to buf 

	tobuf(buf, val);
}

//______________________________________________________________________
void AliDCSMessage::SetFloat(char* buf, Float_t val) 
{
// Set float value to buf 

	tobuf(buf, val);
}

//______________________________________________________________________
Bool_t AliDCSMessage::GetBool(const char* buf) 
{
// get bool value from buf 

	Bool_t val;
	char* aBuffer = (char*) buf;

	frombuf(aBuffer, &val);

	return val;
}

//______________________________________________________________________
Char_t AliDCSMessage::GetByte(const char* buf) 
{
// get byte value from buf 

	Char_t val;
	char* aBuffer = (char*) buf;

	frombuf(aBuffer, &val);

	return val;
}

//______________________________________________________________________
UChar_t AliDCSMessage::GetUByte(const char* buf) 
{
// get ubyte value from buf 

	UChar_t val;
	char* aBuffer = (char*) buf;

        frombuf(aBuffer, &val);

	return val;
}

//______________________________________________________________________
Int_t AliDCSMessage::GetInt(const char* buf) 
{
// get int value from buf 

	Int_t val;
	char* aBuffer = (char*) buf;

        frombuf(aBuffer, &val);

	return val;
}

//______________________________________________________________________
UInt_t AliDCSMessage::GetUInt(const char* buf) 
{
// get uint value from buf 

	UInt_t val;
	char* aBuffer = (char*) buf;

        frombuf(aBuffer, &val);

	return val;
}

//______________________________________________________________________
Float_t AliDCSMessage::GetFloat(const char* buf) 
{
// get float value from buf 

	Float_t val;
	char* aBuffer = (char*) buf;

        frombuf(aBuffer, &val);

	return val;
}

//______________________________________________________________________
TString AliDCSMessage::GetString(const char* buf, Int_t maxLen) 
{
// get string from buf 

	for (Int_t k = 0; k < maxLen; k ++) {
		if (buf[k] == 0) {
			return TString(buf);
		}
	}

	return TString(buf, maxLen);
}

//______________________________________________________________________
void AliDCSMessage::StoreHeader() 
{
// store header message
	
	SetUByte(fMessage + ID_OFFSET, 'A');
	SetUByte(fMessage + ID_OFFSET + 1, 'D');

	SetUByte(fMessage + VERSION_OFFSET, 1);

	SetUByte(fMessage + TYPE_OFFSET, fType);

	SetUInt(fMessage + BODY_SIZE_OFFSET, fMessageSize - HEADER_SIZE);
}

//______________________________________________________________________
void AliDCSMessage::StoreRequestMessage() 
{
// store request message
	
	fMessageSize = REQUEST_STRING_OFFSET +
		fRequestString.Length() + 1;
	
	fMessage = new char[fMessageSize];

	StoreHeader();

	SetUByte(fMessage + REQUEST_TYPE_OFFSET, fRequestType);
	SetUInt(fMessage + START_TIME_OFFSET, fStartTime);
	SetUInt(fMessage + END_TIME_OFFSET, fEndTime);
	strcpy(fMessage + REQUEST_STRING_OFFSET, fRequestString.Data());
}

//______________________________________________________________________
void AliDCSMessage::StoreCountMessage() 
{
// store count message

	fMessageSize = COUNT_OFFSET + sizeof(UInt_t);

	fMessage = new char[fMessageSize];

	StoreHeader();

	SetUInt(fMessage + COUNT_OFFSET, fCount);
}

//______________________________________________________________________
void AliDCSMessage::StoreResultSetMessage()
{
// store result set message

  TIter iter(fValues);
  AliDCSValue* aValue;

  UInt_t valueDataSize = 0;
  while ((aValue = (AliDCSValue*) iter.Next())) {
    valueDataSize += aValue->GetSize();
  }

  fMessageSize = VALUES_OFFSET + valueDataSize;

  fMessage = new char[fMessageSize];

  StoreHeader();

  SetUByte(fMessage + SVT_OFFSET, fValueType);
  SetUInt(fMessage + VALUE_COUNT_OFFSET, GetValueCount());

  UInt_t cursor = VALUES_OFFSET;

  iter.Reset();

  if (fValueType == AliDCSValue::kBool) {
    while ((aValue = (AliDCSValue*) iter.Next())) {
      SetBool(fMessage + cursor, aValue->GetBool());
      cursor += 1;
      SetUInt(fMessage + cursor, aValue->GetTimeStamp());
            cursor += sizeof(UInt_t);
    }
  } else if (fValueType == AliDCSValue::kChar) {
    while ((aValue = (AliDCSValue*) iter.Next())) {
      SetByte(fMessage + cursor, aValue->GetChar());
      cursor += sizeof(Char_t);
      SetUInt(fMessage + cursor, aValue->GetTimeStamp());
      cursor += sizeof(UInt_t);
    }
  } else if (fValueType == AliDCSValue::kInt) {
    while ((aValue = (AliDCSValue*) iter.Next())) {
      SetInt(fMessage + cursor, aValue->GetInt());
      cursor += sizeof(Int_t);
      SetUInt(fMessage + cursor, aValue->GetTimeStamp());
      cursor += sizeof(UInt_t);
    }
  } else if (fValueType == AliDCSValue::kUInt) {
    while ((aValue = (AliDCSValue*) iter.Next())) {
      SetUInt(fMessage + cursor, aValue->GetUInt());
      cursor += sizeof(UInt_t);
      SetUInt(fMessage + cursor, aValue->GetTimeStamp());
      cursor += sizeof(UInt_t);
    }
  } else if (fValueType == AliDCSValue::kFloat) {
    while ((aValue = (AliDCSValue*) iter.Next())) {
      SetFloat(fMessage + cursor, aValue->GetFloat());
      cursor += sizeof(Float_t);
      SetUInt(fMessage + cursor, aValue->GetTimeStamp());
      cursor += sizeof(UInt_t);
    }
  } else {
    AliError("Invalid or unknown ValueType!");
    return;
  }	

}

//______________________________________________________________________
void AliDCSMessage::StoreErrorMessage() 
{
// store error message

	fMessageSize = ERROR_STRING_OFFSET + fErrorString.Length() + 1;

	fMessage = new char[fMessageSize];

	StoreHeader();

	SetUByte(fMessage + ERROR_CODE_OFFSET, fErrorCode);
	strcpy(fMessage + ERROR_STRING_OFFSET, fErrorString.Data());
}

//______________________________________________________________________
void AliDCSMessage::StoreMultiRequestMessage() 
{
// store multi request message
	
	UInt_t requestDataSize = 0;

	TIter iter(&fRequestStrings);
	TObjString* anObjString;

	while ((anObjString = (TObjString*) iter.Next())) {
		assert(anObjString->String().Length() <= 255);
		requestDataSize += anObjString->String().Length() + 1;
	}

	fMessageSize = REQUEST_STRINGS_OFFSET + requestDataSize;

	fMessage = new char[fMessageSize];
	
	StoreHeader();

	SetUByte(fMessage + REQUEST_TYPE_OFFSET, fRequestType);
        SetUInt(fMessage + START_TIME_OFFSET, fStartTime);
        SetUInt(fMessage + END_TIME_OFFSET, fEndTime);
	
	iter.Reset();

	UInt_t cursor = REQUEST_STRINGS_OFFSET;

	while ((anObjString = (TObjString*) iter.Next())) {
		UChar_t strLength = anObjString->String().Length();
		SetUByte(fMessage + cursor, strLength);
		cursor += 1;
		strncpy(fMessage + cursor, anObjString->String().Data(), 
			strLength);
		cursor += strLength;
	}
}

//______________________________________________________________________
/*
void AliDCSMessage::StoreNextMessage() {

        fMessageSize = HEADER_SIZE;

        fMessage = new char[fMessageSize];

        StoreHeader(); 
} */

//______________________________________________________________________
Bool_t AliDCSMessage::ValidateHeader(const char* buf) 
{
// validate message header

	if (!(buf[ID_OFFSET] == 'A' && buf[ID_OFFSET + 1] == 'D')) {
		AliError("Bad message ID!");
		return kFALSE;
	}

	if (buf[VERSION_OFFSET] != 1) {
		AliError("Bad message version!");
		return kFALSE;
	}

	Type type = (Type) GetUByte(buf + TYPE_OFFSET);	
	switch (type) {
		case kRequest:
		case kCount:
		case kResultSet:
		case kError:
		case kMultiRequest:
			break;
		default:
			AliError("Unknown message type!");
			return kFALSE;
	}

	UInt_t bodySize = GetInt(buf + BODY_SIZE_OFFSET);
	if (bodySize > MAX_BODY_SIZE) {
		AliError("Too big message body size!");
		return kFALSE;
	} 

	return kTRUE;
}

//______________________________________________________________________
void AliDCSMessage::LoadRequestMessage() 
{
// load request message

	if (fMessageSize < REQUEST_STRING_OFFSET) {
		AliError("Body size is too small for request message!");
		return;
	}

	fRequestType = (RequestType) GetUByte(fMessage + REQUEST_TYPE_OFFSET);

	fStartTime = GetUInt(fMessage + START_TIME_OFFSET);
	fEndTime = GetUInt(fMessage + END_TIME_OFFSET);
	fRequestString = GetString(fMessage + REQUEST_STRING_OFFSET,
		fMessageSize - REQUEST_STRING_OFFSET);

	switch (fRequestType) {
		case kAlias:
		case kDPName:
			fType = kRequest;
			break;
		default:
			AliError("Invalid request type!");
	}
}

//______________________________________________________________________
void AliDCSMessage::LoadCountMessage() 
{
// load count message

	if (fMessageSize < HEADER_SIZE + sizeof(UInt_t)) {
		AliError("Body size is too small for count message!");
		return;
	}

	fCount = GetUInt(fMessage + COUNT_OFFSET);

	fType = kCount;
}

//______________________________________________________________________
void AliDCSMessage::LoadResultSetMessage()
{
  // load result message

  if (fMessageSize < VALUES_OFFSET) {
    AliError("Body size is too small for result set message!");
    return;
  }

  fValueType = (AliDCSValue::Type) GetUByte(fMessage + SVT_OFFSET);
  UInt_t count = GetUInt(fMessage + VALUE_COUNT_OFFSET);

  UInt_t cursor = VALUES_OFFSET;

  if (fValueType == AliDCSValue::kBool) {
    if (VALUES_OFFSET + count + count * sizeof(UInt_t) >
      fMessageSize) {
      AliError("Too many bool values for this buffer size!");
      return;
    }

    for (UInt_t k = 0; k < count; k ++) {
      Bool_t aBool = GetBool(fMessage + cursor);
      cursor += 1;
      UInt_t timeStamp = GetUInt(fMessage + cursor);
      cursor += sizeof(UInt_t);
      fValues->Add(new AliDCSValue(aBool, timeStamp));
    }
  } else if (fValueType == AliDCSValue::kChar) {
    if (VALUES_OFFSET + count + count * sizeof(UInt_t) >
      fMessageSize) {
      AliError("Too many byte values for this buffer size!");
      return;
    }

    for (UInt_t k = 0; k < count; k ++) {
      Char_t aByte = GetByte(fMessage + cursor);
      cursor += sizeof(Char_t);
      UInt_t timeStamp = GetUInt(fMessage + cursor);
      cursor += sizeof(UInt_t);
      fValues->Add(new AliDCSValue(aByte, timeStamp));
    }
  } else if (fValueType == AliDCSValue::kInt) {
    if (VALUES_OFFSET + count * sizeof(Int_t) +
      count * sizeof(UInt_t) > fMessageSize) {
            AliError("Too many int values for this buffer size!");
            return;
    }

    for (UInt_t k = 0; k < count; k ++) {
            Int_t aInt = GetInt(fMessage + cursor);
            cursor += sizeof(Int_t);
            UInt_t timeStamp = GetUInt(fMessage + cursor);
            cursor += sizeof(UInt_t);
            fValues->Add(new AliDCSValue(aInt, timeStamp));
    }

  } else if (fValueType == AliDCSValue::kUInt) {
    if (VALUES_OFFSET + count * sizeof(UInt_t) +
      count * sizeof(UInt_t) > fMessageSize) {
      AliError("Too many uint values for this buffer size!");
      return;
    }

    for (UInt_t k = 0; k < count; k ++) {
      UInt_t aUInt = GetUInt(fMessage + cursor);
      cursor += sizeof(UInt_t);
      UInt_t timeStamp = GetUInt(fMessage + cursor);
      cursor += sizeof(UInt_t);
      fValues->Add(new AliDCSValue(aUInt, timeStamp));
    }
  } else if (fValueType == AliDCSValue::kFloat) {
    if (VALUES_OFFSET + count * sizeof(Float_t) +
      count * sizeof(UInt_t) > fMessageSize) {
      AliError("Too many float values for this buffer size!");
      return;
    }

    for (UInt_t k = 0; k < count; k ++) {
      Float_t aFloat = GetFloat(fMessage + cursor);
      cursor += sizeof(Float_t);
      UInt_t timeStamp = GetUInt(fMessage + cursor);
      cursor += sizeof(UInt_t);
      fValues->Add(new AliDCSValue(aFloat, timeStamp));
    }

  } else {
    AliError("Unknown or invalid value type!");
  }

  fType = kResultSet;
}

//______________________________________________________________________
void AliDCSMessage::LoadErrorMessage()
{
// load error message
	
	if (fMessageSize < ERROR_STRING_OFFSET) {
		AliError("Body size is too small for error message!");
		return;
	}

	fErrorCode = (ErrorCode) GetUByte(fMessage + ERROR_CODE_OFFSET);
	fErrorString = GetString(fMessage + ERROR_STRING_OFFSET,
		fMessageSize - ERROR_STRING_OFFSET);

	switch (fErrorCode) {
		case kUnknownAliasDPName:
                case kInvalidTimeRange:
                case kInvalidBufferSize:
                case kInvalidRequest:
                case kUnsupportedType:
                case kUnknownError:
			fType = kError;
			break;
		default:
			AliError("Invalid error code!");
	}
}

//______________________________________________________________________
void AliDCSMessage::LoadMultiRequestMessage() 
{
// load multi request message
	
	if (fMessageSize - HEADER_SIZE < REQUEST_STRINGS_OFFSET) {
		AliError("Body size is too small for multi request message!");
		return;
	}

	fRequestType = (RequestType) GetUByte(fMessage + REQUEST_TYPE_OFFSET);

        fStartTime = GetUInt(fMessage + START_TIME_OFFSET);
        fEndTime = GetUInt(fMessage + END_TIME_OFFSET);

        switch (fRequestType) { 
                case kAlias:
                case kDPName:
                        fType = kRequest;
                        break; 
                default:
                        AliError("Invalid request type!");
			return;
        }

	UInt_t cursor = REQUEST_STRINGS_OFFSET;
	
	while ((cursor < fMessageSize)) {
		UChar_t strSize = GetUByte(fMessage + cursor);
		cursor += 1;

		if (cursor + strSize > fMessageSize) {
			AliError("Invalid multi request message!");
			return;
		}		

		TObjString* anObjString = new TObjString(
			GetString(fMessage + cursor, strSize));
		fRequestStrings.AddLast(anObjString);

		cursor += strSize;
	}	

	fType = kMultiRequest;
}

//______________________________________________________________________
/*
void AliDCSMessage::LoadNextMessage() {
	
	fType = kNext;
} */

//______________________________________________________________________
void AliDCSMessage::StoreToBuffer() 
{
	// Creates an underlying message buffer which can be sent to the socket.

	DestroyBuffer();
	
	switch (fType) {
		case kRequest: 
			StoreRequestMessage();
			break;
		case kCount:
			StoreCountMessage();
			break;
		case kResultSet:
			StoreResultSetMessage();
			break;
		case kError:
			StoreErrorMessage();
			break;
		case kMultiRequest:
			StoreMultiRequestMessage();
			break;
/*		case kNext:
			StoreNextMessage();
			break; */
		default:
			AliError("Can't store to buffer invalid message!");
	}
}

//______________________________________________________________________
void AliDCSMessage::LoadFromBuffer()
{
	// Reads the underlying message buffer and if it's valid message
	// creates the corresponding message.  
	// If not set the message type kInvalid.
	// This buffer is read from the socket.
	
	DestroyMessage();

	if (!fMessage) {
		AliError("Message buffer is empty! Can't load it.");
		return;
	}
	
	if (fMessageSize < HEADER_SIZE) {
		AliError("Invalid message buffer. Too small for the header!");
		return;
	}

	if (!ValidateHeader(fMessage)) {
		AliError("Invalid message header!");
		return;
	}

	UInt_t bodySize = GetUInt(fMessage + BODY_SIZE_OFFSET);
	if (bodySize > fMessageSize - HEADER_SIZE) {
		AliError("Message size is to small for the message body!");
		return;
	}

	fMessageSize = HEADER_SIZE + bodySize;

	Type aType = (Type) GetUByte(fMessage + TYPE_OFFSET);

	switch (aType) {
		case kRequest:
			LoadRequestMessage();
			break;
		case kCount:
			LoadCountMessage();
			break;
		case kResultSet:
			LoadResultSetMessage();
			break;
		case kError:
			LoadErrorMessage();
			break;
		case kMultiRequest:
			LoadMultiRequestMessage();
			break;
/*		case kNext:
			LoadNextMessage();
			break; */
		default:
			AliError("Invalid message type!");
	}	
}

//______________________________________________________________________
AliDCSMessage::RequestType AliDCSMessage::GetRequestType() const 
{
	// Request and MultiRequest.
	// Returns the request type: alias or dp (Data Point)

	if (!(fType == kRequest || fType == kMultiRequest)) {
		AliError("Invalid AliDCSMessage type!");
		return kNoneType;
	}

	return fRequestType;
}

//______________________________________________________________________
UInt_t AliDCSMessage::GetStartTime() const 
{
	// Request and MultiRequest.
	// Returns the request start time. (begining of the time interval).

        if (!(fType == kRequest || fType == kMultiRequest)) {
                AliError("Invalid AliDCSMessage type!");
		return 0;
        }

	return fStartTime;
}

//______________________________________________________________________
UInt_t AliDCSMessage::GetEndTime() const 
{
        // Request and MultiRequest.
        // Returns the request start time. (end of the time interval).

	
        if (!(fType == kRequest || fType == kMultiRequest)) {
                AliError("Invalid AliDCSMessage type!");
                return 0;
        }

        return  fEndTime;
}

//______________________________________________________________________
TString AliDCSMessage::GetRequestString() const 
{
        // Request.
        // Returns the request string. (alias or dp)

        if (fType != kRequest) {
                AliError("Invalid AliDCSMessage type!");
		return TString("");
        }

	return fRequestString;
}

//______________________________________________________________________
Bool_t AliDCSMessage::AddRequestString(const char* request)
{
        // MultRequest.
        // Add a request to the request set.
	// Returns kFALSE in case of invalid request (too long request string).
	// Otherwise returns kTRUE.
 
	
	if (fType != kMultiRequest) {
		AliError("Invalid AliDCSMessage type!");
		return kFALSE;
	}

	if (strlen(request) > 255) {
		AliError("Alias/dpName is too long! Max size 255.");
		return kFALSE;	
	}

	fRequestStrings.AddLast(new TObjString(request));
	return kTRUE;
}

//______________________________________________________________________
void AliDCSMessage::ClearRequestStrings()
{
        // MultRequest.
	// Clears the request set.
 
	fRequestStrings.Delete();
}

//______________________________________________________________________
void AliDCSMessage::GetRequestStrings(TObjArray& result) const
{
        // MultRequest.
        // Returns all request strings in this message.
	// result: container where the requests are returned. Collection of
	// TObjString.


	if (fType != kMultiRequest) {
		AliError("Invalid AliDCSMessage type!");
		return;
	}

	TIter iter(&fRequestStrings);
	TObjString* anObjString;
	
	while ((anObjString = (TObjString*) iter.Next())) {
		result.AddLast(new TObjString(*anObjString));
	}
}

//______________________________________________________________________
UInt_t AliDCSMessage::GetCount() const 
{
        // Count.
        // Returns the total number of values.


        if (fType != kCount) {
                AliError("Invalid AliDCSMessage type!");
                return 0;
        }

	return fCount;
}

//______________________________________________________________________
AliDCSValue::Type AliDCSMessage::GetValueType() const
{
  // ResultSet.
  // Returns simple value type (see AliDCSValue) for the values
  // in this ResultSet.

  if (fType != kResultSet) {
          AliError("Invalid AliDCSMessage type!");
          return AliDCSValue::kInvalid;
  }

  return fValueType;
}

//______________________________________________________________________
UInt_t AliDCSMessage::GetValueCount() const
{
  // ResultSet.
  // Returns the count of values in this ResultSet.


  if (fType != kResultSet) {
          AliError("Invalid AliDCSMessage type!");
          return 0;
  }

  return fValues->GetEntriesFast();
}

//______________________________________________________________________
UInt_t AliDCSMessage::GetValues(TObjArray* result) const
{
  // ResultSet.
  // Returns the number of values got from the message.
  // result: used to return the values. Collection of AliDCSValue.
  // result must be owner of the AliDCSValues because fVaule is not!
  // creator of the result array and used GetValues to fill it must delete object by himself!

  // TODO do not copy -> corrected?

	if (fType != kResultSet) {
                AliError("Invalid AliDCSMessage type!");
                return 0;
        }

	TIter iter(fValues);
	AliDCSValue* aValue;
	
	while ((aValue = (AliDCSValue*) iter.Next())) {
		result->AddLast(aValue);
	}

	return fValues->GetEntriesFast();
}


//______________________________________________________________________
Bool_t AliDCSMessage::AddValue(AliDCSValue& value)
{
  // Adds value to the ResultSet value list.
  // Returns kFALSE in case of error.
  // Otherwise returns kTRUE;

  if (fType != kResultSet) {
    AliError("Invalid AliDCSMessage type!");
    return kFALSE;
  }

  if (value.GetType() != fValueType) {
    AliError(Form("Can't add value with type %d to this message!", value.GetType()));
    return kFALSE;
  }

  fValues->Add(&value);

  return kTRUE;
}


//______________________________________________________________________
void AliDCSMessage::ClearValues()
{
// clear values array

	if(fValues) fValues->Clear();
}

//______________________________________________________________________
AliDCSMessage::ErrorCode AliDCSMessage::GetErrorCode() const 
{
	//
	// Error.
	// Returns the error code which has this error message.
	//

        if (fType != kError) {
                AliError("Invalid AliDCSMessage type!");
                return kNoneError;
        }

	return fErrorCode;
}

//______________________________________________________________________
TString AliDCSMessage::GetErrorString() const 
{
	//
	// Error.
	// Returns the error string (error description) which has this 
	// error message.
	//

        if (GetType() != kError) {
                AliError("Invalid AliDCSMessage type!");
                return TString("");
        }

	return fErrorString;
}


//______________________________________________________________________
void AliDCSMessage::Print(Option_t* /*option*/) const 
{
// print message

	if (AliLog::GetGlobalDebugLevel() < 2) {
		return;
	}

	TString printString;
	printString += "\n <<AliDCSMessage>>\n";

	printString += " Size: ";
	printString += fMessageSize;
	printString += '\n';

	printString += " Type: ";
	switch (GetType()) {
		case kRequest: {
			printString += "Request\n";

			printString += " RequestType: ";
			if (GetRequestType() == kDPName) {
				printString += "DPName";
			} else {
				printString += "Alias";
			}
			printString += '\n';

			printString += " RequestString: ";
			printString += GetRequestString();
			printString += '\n';
			printString += " StartTime: ";
			printString += GetStartTime();
			printString += '\n';
			printString += " EndTime: ";
			printString += GetEndTime();
			printString += '\n';
			break;
		}

		case kCount: {
			printString += "Count\n";
			printString += " Count: ";
			printString += GetCount();
			printString += '\n';	
			break;
		} 

		case kResultSet: {
			printString += "ResultSet\n";
			printString += " SimpleValueType: ";
			printString += fValueType;
			printString += '\n';
			printString += " ValueCount: ";
			printString += GetValueCount();
			printString += '\n';
			break;
		}

		case kError: {
			printString += "Error\n";
			printString += " ErrorCode: ";
			switch (GetErrorCode()) {
				case AliDCSMessage::kNoneError:
					printString += "NoneError";
					break;
				case AliDCSMessage::kUnknownAliasDPName:
					printString += "UnknownAliasDPName";
					break;
				case AliDCSMessage::kInvalidTimeRange:
					printString += "InvalidTimeRange";
					break;
				case AliDCSMessage::kInvalidBufferSize:
					printString += "InvalidBufferSize";
					break;
				case AliDCSMessage::kInvalidRequest:
					printString += "InvalidRequest";
					break;
				case AliDCSMessage::kUnsupportedType:
					printString += "UnsupportedType";
					break;
				case AliDCSMessage::kUnknownError:
					printString += "UnknownError";
					break;
				default:
					printString += "Invalid";
			}

			printString += '\n';
			printString += " ErrorString: ";
			printString += GetErrorString();
			printString += '\n';
			break;
		}

		case kMultiRequest: {
			printString += "MultiRequest\n";

                        printString += " RequestType: ";
                        if (GetRequestType() == kDPName) {
                                printString += "DPName";
                        } else {
                                printString += "Alias";
                        }
                        printString += '\n';

                        printString += " RequestStrings: ";
			TIter iter(&fRequestStrings);
			TObjString* anObjString;
			while ((anObjString = (TObjString*) iter.Next())) {
				printString += anObjString->String();
				printString += ' ';
			}
                        printString += '\n';

                        printString += " StartTime: ";
                        printString += GetStartTime();
                        printString += '\n';
                        printString += " EndTime: ";
                        printString += GetEndTime();
                        printString += '\n';
			break;
		} 	

/*		case kNext: {
			printString += "Next\n";
			break;
		} */

		default:
			printString += "Invalid\n";
	}

	if (AliLog::GetGlobalDebugLevel() >= 3 && fMessage) {
		PrintBuffer(fMessage, fMessageSize, printString);
	}

	AliDebug(2, printString);
} 

//______________________________________________________________________
Bool_t AliDCSMessage::SetRawHeader(const char* header) 
{
	//
	// Checks if the header buffer represents a valid header message.
	// If so it creates a message buffer with the appropriate body size
	// and returns true.
	// If not returns false.
	// header: header buffer
	//

	if (!ValidateHeader(header)) {
		AliError("Invalid message header!");
		return kFALSE;
	}
	
	DestroyBuffer();

	UInt_t bodySize = GetUInt(header + BODY_SIZE_OFFSET);
	fMessageSize = HEADER_SIZE + bodySize;
	
	fMessage = new char[fMessageSize];
	
	memcpy(fMessage, header, HEADER_SIZE);

	return kTRUE;
}


//______________________________________________________________________
void AliDCSMessage::DestroyBuffer() 
{
	//
	// Destroy the underlying message buffer.
	//
	
	if (fMessage) {
		delete[] fMessage;
		fMessage = NULL;
	}

	fMessageSize = 0;
}

//______________________________________________________________________
void AliDCSMessage::PrintBuffer(const char* buffer, UInt_t size, 
		TString& output)
{
// print buffer

	UInt_t index = 0;

	while (index < size) {
		if (!(index % 16)) {
			output += Form("\n %.4x:", index); 
		} 	

		if (!(index % 8)) {
			output += ' ';
		}

		output += Form(" %.2x", (UChar_t) buffer[index]);

		if (!((index + 1) % 16) || index + 1 == size) {
			if (index + 1 == size) {
				output.Append(' ',3 * (15 - index % 16));
				if (index % 16 < 8) {
					output.Append(' ');
				}
			}

			output.Append(' ', 2);
			for (Int_t k = index % 16; k >= 0; k --) {
				Char_t aChar = buffer[index - k];
				output += isgraph(aChar) ? aChar: '.';
			}
		} 

		index ++;
	}

	output += '\n';
}
