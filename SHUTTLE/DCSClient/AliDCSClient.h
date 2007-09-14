#ifndef ALI_DCS_CLIENT_H
#define ALI_DCS_CLIENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class represents the AliDCSClient.
// The client used for data retrieval from DCS server.
// For more info see AliDCSClient.cxx
//

#include "AliDCSMessage.h"

class TObjArray;
class TSocket;
class TMap;
class TCollection;

class AliDCSClient: public TObject {
public:

	enum {
		fgkBadState=-1, 	     // Bad state
		fgkInvalidParameter = -2,    // Invalid parameter
		fgkTimeout = -3,	     // Timeout
		fgkBadMessage = -4,	     // Bad message
		fgkCommError = -5,	     // Communication error
		fgkServerError = -6	     // Server error
	};
	
////	friend class AliShuttle;

	AliDCSClient(const char* host, Int_t port, UInt_t timeout = 5000,
                        Int_t retries = 5, Int_t multiSplit = 100);
        virtual ~AliDCSClient();


        Int_t GetDPValues(const char* dpName, UInt_t startTime, UInt_t endTime,
                                TObjArray* result);

        Int_t GetAliasValues(const char* alias, UInt_t startTime,
                                UInt_t endTime, TObjArray* result);

        TMap* GetDPValues(const TSeqCollection* dpList, UInt_t startTime, UInt_t endTime, Int_t startIndex = 0, Int_t endIndex = -1);

        TMap* GetAliasValues(const TSeqCollection* aliasList, UInt_t startTime, UInt_t endTime, Int_t startIndex = 0, Int_t endIndex = -1);

        AliDCSMessage::ErrorCode GetServerErrorCode() const
                { return fServerErrorCode;}

        const TString& GetServerError() const {return fServerError;}


        Bool_t IsConnected();

        void Close();


        static const char* GetErrorString(Int_t code);

private:
	static const char* fgkBadStateString;		// Bad state string
	static const char* fgkInvalidParameterString;	// Invalid parameter string
	static const char* fgkTimeoutString;    	// Timeout string
	static const char* fgkBadMessageString; 	// Bad message string
	static const char* fgkCommErrorString;  	// Communication error string
	static const char* fgkServerErrorString;	// Server error string

	TSocket* fSocket;	// Pointer to the TCP socket client
	TString fHost;  	// server host
	Int_t   fPort;		// server port
	UInt_t fTimeout;	// timeout parameter
	Int_t fRetries;		// number of retries
  	Int_t   fMultiSplit; // number of datapoints that are queried at a time in a multi dp request, if set to 1 forces single requests 
	AliDCSMessage::ErrorCode fServerErrorCode;	// error code
	TString fServerError;	// server error string

	Bool_t Connect();

	Int_t SendBuffer(const char* buffer, Int_t size);

	Int_t ReceiveBuffer(char* buffer, Int_t size);

	Int_t SendMessage(AliDCSMessage& message);

	Int_t ReceiveMessage(AliDCSMessage& message);

	Int_t GetValues(AliDCSMessage::RequestType requestType,
		const char* requestString, UInt_t startTime, UInt_t endTime,
		TObjArray* result);

	TMap* GetValues(AliDCSMessage::RequestType requestType,
		const TSeqCollection* list, UInt_t startTime, UInt_t endTime,
		Int_t startIndex, Int_t endIndex);

	Int_t ReceiveValueSet(TObjArray* result);

	AliDCSClient(const AliDCSClient& other);		// not implemented
	AliDCSClient& operator= (const AliDCSClient& other);	// not implemented

	ClassDef(AliDCSClient, 0);
};

#endif
