#ifndef ALI_DCS_MESSAGE_H
#define ALI_DCS_MESSAGE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is a wrapper of AliDCSMessage.
// These are the messages which form AliDCSProtocol.
// Used by AliDCSClient to communicate with the DCS Amanda server
//

#include <TString.h>
#include <TObjArray.h>
#include "AliDCSValue.h"

#define HEADER_SIZE 8
#define ID_OFFSET 0
#define VERSION_OFFSET 2
#define TYPE_OFFSET 3
#define BODY_SIZE_OFFSET 4

#define MAX_BODY_SIZE 40000

#define REQUEST_TYPE_OFFSET HEADER_SIZE
#define START_TIME_OFFSET (HEADER_SIZE + 1)
#define END_TIME_OFFSET (HEADER_SIZE + 5)
#define REQUEST_STRING_OFFSET (HEADER_SIZE + 9)

#define REQUEST_STRINGS_OFFSET REQUEST_STRING_OFFSET

#define COUNT_OFFSET HEADER_SIZE

#define SVT_OFFSET HEADER_SIZE
#define VALUE_COUNT_OFFSET (HEADER_SIZE + 1)
#define VALUES_OFFSET (HEADER_SIZE + 5)

#define ERROR_CODE_OFFSET HEADER_SIZE
#define ERROR_STRING_OFFSET (HEADER_SIZE + 1)

class AliDCSMessage: public TObject {
public:
	enum Type {
		kInvalid = 0,
		kRequest = 1,
		kCount = 2,
		kResultSet = 3,
		kError = 4,
		kMultiRequest = 5
//		kNext = 6
	};

	enum RequestType {
		kNoneType = 0,
		kAlias = 1,
		kDPName = 2	
	};

	enum ErrorCode {
		kNoneError = 0,
		kUnknownAliasDPName = 1,
		kInvalidTimeRange = 2,
		kInvalidBufferSize = 3,
		kInvalidRequest = 4,
		kUnsupportedType = 5,
		kUnknownError = 255
	};


	AliDCSMessage();

        AliDCSMessage(const char* buffer, UInt_t size);

        virtual ~AliDCSMessage();


        void CreateRequestMessage(RequestType type, 
		UInt_t startTime, UInt_t endTime, const char* request);

        void CreateMultiRequestMessage(RequestType type, 
		UInt_t startTime, UInt_t endTime);

        void CreateCountMessage(UInt_t count);

        void CreateResultSetMessage(AliDCSValue::Type type);

        void CreateErrorMessage(ErrorCode code, const char* errorString);
	
	//void CreateNextMessage();

	void DestroyMessage();


        Bool_t SetRawHeader(const char* header);

        void StoreToBuffer();

        void LoadFromBuffer();

        void DestroyBuffer();


        Bool_t IsValid() const {return fType != kInvalid;};

        UInt_t GetMessageSize() const {return fMessageSize;};

        char* GetMessage() const {return fMessage;};

        UInt_t GetBodySize() const {return fMessageSize - HEADER_SIZE;};

        char* GetBody() const {return fMessage + HEADER_SIZE;};

        Type GetType() const {return fType;};

	 // RequestType and MultiReuqestType Message getters
        RequestType GetRequestType() const;

        UInt_t GetStartTime() const;

        UInt_t GetEndTime() const;

        TString GetRequestString() const;

        // MultiRequestType Message getters and setters
        void GetRequestStrings(TObjArray& result) const;

        Bool_t AddRequestString(const char* request);

        void ClearRequestStrings();

        // CountType Message getters            
        UInt_t GetCount() const;

        // ResultSetType Message getters ans setters       
        AliDCSValue::Type GetValueType() const;

        UInt_t GetValueCount() const;

        UInt_t GetValues(TObjArray* result) const;

        Bool_t AddValue(AliDCSValue& value); 

        void ClearValues();

        // ErrorType Message getters    
        ErrorCode GetErrorCode() const;

        TString GetErrorString() const;


        virtual void Print(Option_t* option = NULL) const;

        static void PrintBuffer(const char* buf, UInt_t size, TString& output);

private:

	AliDCSMessage(const AliDCSMessage& other); 	
	AliDCSMessage& operator= (const AliDCSMessage& other); 	


	char* fMessage; 	// Array of bytes building the message

	UInt_t fMessageSize; 	// Size of the message array


	Type fType; 		// Message type (request, count...)
	
	//Request message fields
	RequestType fRequestType; 	// Type of request message
	
	UInt_t fStartTime; 		// Start time of query

	UInt_t fEndTime; 		// End time of query

	TString fRequestString; 	// Request string
	
	//Count message fields
	UInt_t fCount; 			// count counter

	//ResultSet message fields
	AliDCSValue::Type fValueType; // Simple value type

	TObjArray* fValues; 		// array of received values

	//Error message fields
	ErrorCode fErrorCode; 		// error code
	
	TString fErrorString; 		// error string

	//MultiRequest message fields
	TObjArray fRequestStrings; 	// multi request string array

	
	// Message setter helpers
	void StoreHeader();

	void StoreRequestMessage();

	void StoreCountMessage();

	void StoreResultSetMessage();

	void StoreErrorMessage();

	void StoreMultiRequestMessage();

	//void StoreNextMessage();


	Bool_t ValidateHeader(const char* buf);

	void LoadRequestMessage();

	void LoadCountMessage();

	void LoadResultSetMessage();

	void LoadErrorMessage();

	void LoadMultiRequestMessage();

	//void LoadNextMessage();

	// Buffer helpers
	static void SetBool(char* buf, Bool_t val);
	
	static void SetByte(char* buf, Char_t val);

	static void SetUByte(char* buf, UChar_t val);

	static void SetInt(char* buf, Int_t val);
		
	static void SetUInt(char* buf, UInt_t val);

	static void SetFloat(char* buf, Float_t val);

	static Bool_t GetBool(const char* buf);

	static Char_t GetByte(const char* buf);

	static UChar_t GetUByte(const char* buf);

	static Int_t GetInt(const char* buf);

	static UInt_t GetUInt(const char* buf);

	static Float_t GetFloat(const char* buf);

	static TString GetString(const char* buf, Int_t maxLen);

	
	ClassDef(AliDCSMessage, 0);
};

#endif
