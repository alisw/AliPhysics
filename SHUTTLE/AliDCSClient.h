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

class AliDCSClient: public TObject {
public:

	friend class AliShuttle;

	AliDCSClient(const char* host, Int_t port, UInt_t timeout = 5000,
                        Int_t retries = 5);
        virtual ~AliDCSClient();


        Int_t GetDPValues(const char* dpName, UInt_t startTime, UInt_t endTime,
                                TObjArray* result);

        Int_t GetAliasValues(const char* alias, UInt_t startTime,
                                UInt_t endTime, TObjArray* result);

        Int_t GetDPValues(UInt_t startTime, UInt_t endTime, TMap& result);

        Int_t GetAliasValues(UInt_t startTime, UInt_t endTime, TMap& result);


        AliDCSMessage::ErrorCode GetServerErrorCode() const
                { return fServerErrorCode;};

        const TString& GetServerError() const {return fServerError;};


        Bool_t IsConnected();

        void Close();


        static const char* GetErrorString(Int_t code);

private:

	static const Int_t fgkBadState = -1;		// Bad state
	static const Int_t fgkInvalidParameter = -2;	// Invalid parameter
	static const Int_t fgkTimeout = -3;		// Timeout
	static const Int_t fgkBadMessage = -4;		// Bad message
	static const Int_t fgkCommError = -5;		// Communication error
	static const Int_t fgkServerError = -6;		// Server error

	static const char* fgkBadStateString;		// Bad state string
	static const char* fgkInvalidParameterString;	// Invalid parameter string
	static const char* fgkTimeoutString;    	// Timeout string
	static const char* fgkBadMessageString; 	// Bad message string
	static const char* fgkCommErrorString;  	// Communication error string
	static const char* fgkServerErrorString;	// Server error string

	AliDCSClient(const AliDCSClient& other);
	AliDCSClient& operator= (const AliDCSClient& other);

	TSocket* fSocket;	// Pointer to the TCP socket client

	UInt_t fTimeout;	// timeout parameter

	Int_t fRetries;		// number of retries

	AliDCSMessage::ErrorCode fServerErrorCode;	// error code

	TString fServerError;	// server error string


	Int_t SendBuffer(const char* buffer, Int_t size);

	Int_t ReceiveBuffer(char* buffer, Int_t size);

	Int_t SendMessage(AliDCSMessage& message);

	Int_t ReceiveMessage(AliDCSMessage& message);

	Int_t GetValues(AliDCSMessage::RequestType requestType,
		const char* requestString, UInt_t startTime, UInt_t endTime,
		TObjArray* result);

	Int_t GetValues(AliDCSMessage::RequestType requestType,
		UInt_t startTime, UInt_t endTime, TMap& result);

	Int_t ReceiveValueSet(TObjArray* result);


	ClassDef(AliDCSClient, 0);
};

#endif
