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
Revision 1.6  2006/10/02 16:38:39  jgrosseo
update (alberto):
fixed memory leaks
storing of objects that failed to be stored to the grid before
interfacing of shuttle status table in daq system

Revision 1.5  2006/08/15 10:50:00  jgrosseo
effc++ corrections (alberto)

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

Revision 1.3  2005/08/30 10:53:23  byordano
some more descriptions added

*/

//
// This class represents the AliDCSClient.
// The client used for data retrieval from DCS server.
// There are two way for retrieving data from the server.
//	1) asking for DP (DataPoint) - usually changed frequently.
//	2) asking for Alias (Alias) - alias should be the same through whole
//		experimnet.
//		
// There are two type of read operations:
//	Asking for single alias/dp or asking for set of aliases/dp
//
// In case of ServerError the coresponding error code and 
// error string (description) could be got by GetServerErrorCode() and
// GetServerErrorString()
//

#include "AliDCSClient.h"
#include "AliDCSValue.h"
#include "AliLog.h"

#include <TSocket.h>
#include <TObjArray.h>
#include <TMap.h>
#include <TObjString.h>
#include <TSystem.h>

ClassImp(AliDCSClient)

const char* AliDCSClient::fgkBadStateString = "BadState";
const char* AliDCSClient::fgkInvalidParameterString = "InvalidParameter";
const char* AliDCSClient::fgkTimeoutString = "Timeout";
const char* AliDCSClient::fgkBadMessageString = "BadMessage";
const char* AliDCSClient::fgkCommErrorString = "CommunicationError";
const char* AliDCSClient::fgkServerErrorString = "ServerError";

//______________________________________________________________________
AliDCSClient::AliDCSClient(const char* host, Int_t port, UInt_t timeout,
	Int_t retries):
	fSocket(NULL), fTimeout(timeout), fRetries(retries),
	fServerErrorCode(AliDCSMessage::kNoneError), fServerError("")
{
	// 
	// host: DCS server host
	// port: DCS server port
	// timeout: in case of communication error or socket read/write this 
	// 		timeout will be used before the next try is made. 
	// retries: the number of retries after which the connection is 
	//		is considered as invalid and error is returned.
	//

	Int_t tries = 0;	
	
	while (tries < fRetries) {
		fSocket = new TSocket(host, port);
		if (fSocket->IsValid()) {
			AliDebug(1, Form("Connected to %s:%d", host, port));
			fSocket->SetOption(kNoBlock, 1);
			break;
		}

		AliDebug(1, Form("Connection timeout! tries <%d> ...", tries));

		delete fSocket;
		fSocket = NULL;

		gSystem->Sleep(fTimeout);
		tries ++;
	}
}

//______________________________________________________________________
AliDCSClient::AliDCSClient(const AliDCSClient& /*other*/):
	TObject(), fSocket(NULL), fTimeout(0), fRetries(0),
	fServerErrorCode(AliDCSMessage::kNoneError), fServerError("")

{
// copy constructor (not implemented)

}

//______________________________________________________________________
AliDCSClient &AliDCSClient::operator=(const AliDCSClient& /*other*/)
{
// assignment operator (not implemented)

return *this;
}

//______________________________________________________________________
AliDCSClient::~AliDCSClient() 
{
// destructor

	if (fSocket) {
		Close();
		delete fSocket;
	}
}

//______________________________________________________________________
Int_t AliDCSClient::SendBuffer(const char* buffer, Int_t size) 
{
// send buffer containing the message to the DCS server

	Int_t sentSize = 0;
        Int_t tries = 0;

        while (sentSize < size && tries < fRetries) {

                Int_t sResult = fSocket->Select(TSocket::kWrite, fTimeout);

                if (sResult == 0) {
			AliDebug(1, Form("Timeout! tries <%d> ...", tries));
                        tries ++;
                        continue;

                } else if (sResult < 0) {
			AliDebug(1, Form("Communication error <%d>!", 
					fSocket->GetErrorCode()));
                        return AliDCSClient::fgkCommError;
                }

                sResult = fSocket->SendRaw(buffer + sentSize, size - sentSize,
				 	kDontBlock);

                if (sResult > 0) {
                        sentSize += sResult;
                } else {
			AliDebug(1, Form("Communication error <%d>!",
                                        fSocket->GetErrorCode()));
                        return AliDCSClient::fgkCommError;
                }
        }

	if (tries == fRetries) {
		return AliDCSClient::fgkTimeout;
	}

        return sentSize;
}

//______________________________________________________________________
Int_t AliDCSClient::ReceiveBuffer(char* buffer, Int_t size) 
{
// Receive message from the DCS server and fill buffer

	Int_t receivedSize = 0;
        Int_t tries = 0;

        while (receivedSize < size && tries < fRetries) {

                Int_t sResult = fSocket->Select(TSocket::kRead, fTimeout);

                if (sResult == 0) {
                        AliDebug(1, Form("Timeout! tries <%d> ...", tries));
                        tries ++;
                        continue;

                } else if (sResult < 0) {
                        AliDebug(1, Form("Communication error <%d>",
                                        fSocket->GetErrorCode()));
                        return AliDCSClient::fgkCommError;
                }

                sResult = fSocket->RecvRaw(buffer + receivedSize, 
			size - receivedSize, kDontBlock);

                if (sResult > 0) {
                        receivedSize += sResult;
                } else {
                        AliDebug(1, Form("Communication error <%d>",
                                        fSocket->GetErrorCode()));
                        return AliDCSClient::fgkCommError;
                }
        }

        if (tries == fRetries) {
                return AliDCSClient::fgkTimeout;
        }

        return receivedSize;
}

//______________________________________________________________________
Int_t AliDCSClient::SendMessage(AliDCSMessage& message) 
{
// send message to the DCS server

	message.StoreToBuffer();

	AliDebug(2, "Sending message.\n"); 
	message.Print();

	return SendBuffer(message.GetMessage(), message.GetMessageSize());
}

//______________________________________________________________________
Int_t AliDCSClient::ReceiveMessage(AliDCSMessage& message) 
{
// receive message from the DCS server
	
	char header[HEADER_SIZE];

	Int_t sResult;

	if ((sResult = ReceiveBuffer(header, HEADER_SIZE)) < 0) {
		AliDebug(1, Form("Can't receive message header! Reason: %s", 
			GetErrorString(sResult)));
		return sResult;
	}

	if (!message.SetRawHeader(header)) {
		return AliDCSClient::fgkBadMessage;
	}

	if ((sResult = ReceiveBuffer(message.GetBody(), 
			message.GetBodySize())) < 0) {
		
		AliDebug(1, Form("Can't receive message body! Reason: %s", 
			GetErrorString(sResult)));
		return sResult;
	}

	message.LoadFromBuffer();

	AliDebug(2, "Message received.");
	message.Print();

	return HEADER_SIZE + sResult;
}

//______________________________________________________________________
Int_t AliDCSClient::GetValues(AliDCSMessage::RequestType reqType,
	const char* reqString, UInt_t startTime, UInt_t endTime, TObjArray* result)
{
// get array of DCS values from the DCS server
// reqString: alias name
// startTime, endTime: start time and end time of the query
// result: contains the array of retrieved AliDCSValue's

	if (!IsConnected()) {
		AliError("Not connected!");
		return AliDCSClient::fgkBadState;
	}

	Int_t sResult;
	AliDCSMessage requestMessage;
	requestMessage.CreateRequestMessage(reqType, startTime, endTime,
			reqString);

	if ((sResult = SendMessage(requestMessage)) < 0) {
		AliError(Form("Can't send request message! Reason: %s",
			GetErrorString(sResult)));
		Close();
		return sResult;
	}

	sResult = ReceiveValueSet(result);

	Close();

	return sResult;
}

//______________________________________________________________________
Int_t AliDCSClient::GetValues(AliDCSMessage::RequestType reqType,
	UInt_t startTime, UInt_t endTime, TMap& result)
{
// get array of DCS values from the DCS server
// startTime, endTime: start time and end time of the query
// result: map containing the array of alias names. It will be filled with
// the values retrieved for each alias

	if (!IsConnected()) {
		AliError("Not connected!");
		return AliDCSClient::fgkBadState;
	}

	AliDCSMessage multiRequestMessage;
	multiRequestMessage.CreateMultiRequestMessage(reqType,
			startTime, endTime);

	TObjArray requests;

	TIter iter(&result);
	TObjString* aRequest;

	// copy request strings to temporary TObjArray because
	// TMap doesn't guarantee the order of elements!!!
	while ((aRequest = (TObjString*) iter.Next())) {
		requests.AddLast(aRequest);
		if (!multiRequestMessage.AddRequestString(aRequest->String()))
		{
			return AliDCSClient::fgkInvalidParameter;
		}
	}

	Int_t sResult;
	if ((sResult = SendMessage(multiRequestMessage)) < 0) {
                AliError(Form("Can't send request message! Reason: %s",
                        GetErrorString(sResult)));
                Close();
                return sResult;
        }

	result.SetOwner(0);
	result.Clear();

	TIter reqIter(&requests);
	while ((aRequest = (TObjString*) reqIter.Next())) {
		TObjArray* resultSet = new TObjArray();
		resultSet->SetOwner(1);

		if ((sResult = ReceiveValueSet(resultSet)) < 0) {
			AliError(Form("Can't get values for %s!" ,
				aRequest->String().Data()));

			delete resultSet;
			break;
		}

		result.Add(aRequest, resultSet);
	}

	if (sResult < 0) {
		result.DeleteValues();
		result.Clear();

		requests.Delete();
	} else {
		result.SetOwner(1);
	}

	Close();

	return sResult;
}
	
//______________________________________________________________________
Int_t AliDCSClient::ReceiveValueSet(TObjArray* result)
{
// receive set of values 

	Int_t sResult;

	AliDCSMessage responseMessage;
	if ((sResult = ReceiveMessage(responseMessage)) < 0) {
		AliError(Form("Can't receive response message! Reason: %s",
			GetErrorString(sResult)));
		return sResult;
	}	

	UInt_t valueCount;

	if (responseMessage.GetType() == AliDCSMessage::kCount) {
		valueCount = responseMessage.GetCount();

	} else if (responseMessage.GetType() == AliDCSMessage::kError) {
		fServerErrorCode = responseMessage.GetErrorCode();
		fServerError = responseMessage.GetErrorString();

		return AliDCSClient::fgkServerError;

	} else {
		AliError("Bad message type received!");	
		return AliDCSClient::fgkBadMessage;
	}

	UInt_t receivedValues = 0;

	AliDCSValue::Type valueType = AliDCSValue::kInvalid;

	while (receivedValues < valueCount) {

                AliDCSMessage message;

                if ((sResult = ReceiveMessage(message)) < 0) {
                        AliError(Form("Can't receive message! Reason: %s",
				GetErrorString(sResult)));
                        return sResult;
                }

                if (message.GetType() == AliDCSMessage::kResultSet) {

			if (valueType == AliDCSValue::kInvalid) {
				valueType = message.GetValueType();
			} else {
				if (valueType != message.GetValueType()) {
					AliError("Unexpected value type!");
					return AliDCSClient::fgkBadMessage;
				}
			}

			receivedValues += message.GetValues(result);

                	if (receivedValues > valueCount) {
                        	AliError("Message contains more values than expected!");
	                        return AliDCSClient::fgkBadMessage;
        	        }

		} else if (message.GetType() == AliDCSMessage::kError) {
				fServerErrorCode = 
					responseMessage.GetErrorCode();
				fServerError = responseMessage.GetErrorString();

                		return AliDCSClient::fgkServerError;
		} else {
                        AliError("Bad message type received!");
                        return AliDCSClient::fgkBadMessage;
                }       
        }       
	
	return receivedValues;
}
		
//______________________________________________________________________
Int_t AliDCSClient::GetDPValues(const char* dpName, UInt_t startTime,
				UInt_t endTime, TObjArray* result)
{
	//
	// Reads a values from the server which correspond to this
	// DataPoint (dpName) in time interval (startTime - endTime).
	// result: Collection of AliDCSValue which contains the read values.
	//
	// Returns:
	//	If >= 0 , the number of values read.
	//	if < 0, the error code which has occured during the read.
	//

	return GetValues(AliDCSMessage::kDPName, 
			dpName, startTime, endTime, result);
}

//______________________________________________________________________
Int_t AliDCSClient::GetAliasValues(const char* alias, UInt_t startTime,
				UInt_t endTime, TObjArray* result)
{
	//
        // Reads a values from the server which correspond to this
        // alias (alias) in time interval (startTime - endTime).
        // result: Collection of AliDCSValue which contains the read values.
	//
        // Returns:
        //      If >= 0 , the number of values read.
        //      if < 0, the error code which has occured during the read.
        //

	return GetValues(AliDCSMessage::kAlias, 
			alias, startTime, endTime, result);
}

//______________________________________________________________________
Int_t AliDCSClient::GetDPValues(UInt_t startTime, UInt_t endTime, 
				TMap& result) 
{
	//
        // For every key of 'result' (which must be TObjString) 
	// reads a valueSet. The key represents particular DataPoint to be read.
        // For all DataPoints time interval (startTime - endTime) is used.
        // After the read, the correspoding value for every key is a 
	// TObjArray - collection of AliDCSValue, or result is an empty map in
	// case of error.
	// 
        // Returns:
        //      If >= 0 , the number of values read.
        //      if < 0, the error code which has occured during the read.
        //

	return GetValues(AliDCSMessage::kDPName, startTime, endTime, result);
}

//______________________________________________________________________
Int_t AliDCSClient::GetAliasValues(UInt_t startTime, UInt_t endTime,
				TMap& result)
{
	//
        // For every key of 'result' (which must be TObjString) 
        // reads a valueSet. The key represents particular Alias to be read.
        // For all aliases time interval (startTime - endTime) is used.
        // After the read, the correspoding value for every key is a 
        // TObjArray - collection of AliDCSValue, or result is an empty map in
        // case of error.
        // 
        // Returns:
        //      If >= 0 , the number of values read.
        //      if < 0, the error code which has occured during the read.
        //

	return GetValues(AliDCSMessage::kAlias, startTime, endTime, result);
}

//______________________________________________________________________
Bool_t AliDCSClient::IsConnected() 
{
	//
	// Returns kTRUE if there is a valid connection to the server.
	//

	if (fSocket) {
		return fSocket->IsValid();
	}

	return kFALSE;
}

//______________________________________________________________________
void AliDCSClient::Close() 
{
	//
	// Close the connection.
	//

	if (fSocket) {
		fSocket->Close();
	}
}

//______________________________________________________________________
const char* AliDCSClient::GetErrorString(Int_t code) 
{
	//
	// Returns a short string describing the error code.
	// code: the error code.
	//
	

	switch (code) {
		case AliDCSClient::fgkBadState:
			return AliDCSClient::fgkBadStateString;

		case AliDCSClient::fgkInvalidParameter:
			return AliDCSClient::fgkInvalidParameterString;

		case AliDCSClient::fgkTimeout:
			return AliDCSClient::fgkTimeoutString;

		case AliDCSClient::fgkBadMessage:
			return AliDCSClient::fgkBadMessageString;
	
		case AliDCSClient::fgkCommError:
			return AliDCSClient::fgkCommErrorString;

		case AliDCSClient::fgkServerError:
			return AliDCSClient::fgkServerErrorString;

		default:
			AliErrorGeneral("AliDCSClient::GetErrorString",
				"Unknown error code!");
			return "UnknownCode";
	}
}

