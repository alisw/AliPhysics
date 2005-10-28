#ifndef ALI_DCS_CLIENT_H
#define ALI_DCS_CLIENT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class represents the AliDCSClient.
// The client used for data retrieval from DCS server.
//

#include "AliDCSMessage.h"

#include <TSocket.h>


class TList;
class TMap;

class AliDCSClient: public TObject {
public:
	
	static const Int_t fgkBadState = -1;

	static const Int_t fgkInvalidParameter = -2;

	static const Int_t fgkTimeout = -3;

	static const Int_t fgkBadMessage = -4;

	static const Int_t fgkCommError = -5;

	static const Int_t fgkServerError = -6;

	static const char* fgkBadStateString;

	static const char* fgkInvalidParameterString;

	static const char* fgkTimeoutString;

	static const char* fgkBadMessageString; 

	static const char* fgkCommErrorString;

	static const char* fgkServerErrorString;


	AliDCSClient(const char* host, Int_t port, UInt_t timeout = 5000,
                        Int_t retries = 5);
        virtual ~AliDCSClient();


        Int_t GetDPValues(const char* dpName, UInt_t startTime, UInt_t endTime,
                                TList& result);

        Int_t GetAliasValues(const char* alias, UInt_t startTime,
                                UInt_t endTime, TList& result);

        Int_t GetDPValues(UInt_t startTime, UInt_t endTime, TMap& result);

        Int_t GetAliasValues(UInt_t startTime, UInt_t endTime, TMap& result);


        AliDCSMessage::ErrorCode GetServerErrorCode()
                { return fServerErrorCode;};

        const TString& GetServerError() {return fServerError;};


        Bool_t IsConnected();

        void Close();


        static const char* GetErrorString(Int_t code);

private:

	TSocket* fSocket;
	
	UInt_t fTimeout;

	Int_t fRetries;
	
	AliDCSMessage::ErrorCode fServerErrorCode;
	
	TString fServerError;


	Int_t SendBuffer(const char* buffer, Int_t size);

	Int_t ReceiveBuffer(char* buffer, Int_t size);

	Int_t SendMessage(AliDCSMessage& message);

	Int_t ReceiveMessage(AliDCSMessage& message);	

	Int_t GetValues(AliDCSMessage::RequestType requestType,
		const char* requestString, UInt_t startTime, UInt_t endTime,
		TList& result);
	
	Int_t GetValues(AliDCSMessage::RequestType requestType,
		UInt_t startTime, UInt_t endTime, TMap& result);

	Int_t ReceiveValueSet(TList& result);


	ClassDef(AliDCSClient, 0);
};

#endif
