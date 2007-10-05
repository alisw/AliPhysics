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
Revision 1.4  2007/09/14 16:46:14  jgrosseo
1) Connect and Close are called before and after each query, so one can
keep the same AliDCSClient object.
2) The splitting of a query is moved to GetDPValues/GetAliasValues.
3) Splitting interval can be specified in constructor

Revision 1.3  2007/09/11 16:42:02  jgrosseo
starting modifying AliDCSClient to transparently switch between single and multi query
first step: same alidcsclient can be used for several queries

Revision 1.2  2007/06/09 13:01:09  jgrosseo
Switching to retrieval of several DCS DPs at a time (multiDPrequest)

Revision 1.1  2006/11/06 14:22:47  jgrosseo
major update (Alberto)
o) reading of run parameters from the logbook
o) online offline naming conversion
o) standalone DCSclient package

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
	Int_t retries, Int_t multiSplit):
	fSocket(NULL), fHost(host), fPort(port), fTimeout(timeout), fRetries(retries), fMultiSplit(multiSplit),
	fServerErrorCode(AliDCSMessage::kNoneError), fServerError(""), fResultErrorCode(0)
{
	// 
	// host: DCS server host
	// port: DCS server port
	// timeout: in case of communication error or socket read/write this 
	// 		timeout will be used before the next try is made. 
	// retries: the number of retries after which the connection is 
	//		is considered as invalid and error is returned.
	//
}

//______________________________________________________________________
AliDCSClient::~AliDCSClient() 
{
	// destructor

	Close();
}

//______________________________________________________________________
Bool_t AliDCSClient::Connect()
{
	//  connects to the AMANDA server
	
	Close();
	
	Int_t tries = 0;	
	while (tries < fRetries) 
	{
		fSocket = new TSocket(fHost, fPort);
		if (fSocket->IsValid()) 
		{
			AliDebug(1, Form("Connected to %s:%d", fHost.Data(), fPort));
			fSocket->SetOption(kNoBlock, 1);
			return kTRUE;
		}

		AliDebug(1, Form("Connection timeout! tries <%d> ...", tries));

		delete fSocket;
		fSocket = NULL;

		gSystem->Sleep(fTimeout);
		tries ++;
	}
	
	return kFALSE;
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

	Connect();
	
	if (!IsConnected()) {
		AliError("Not connected!");
		return AliDCSClient::fgkBadState;
	}

	Int_t sResult = -1;
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
TMap* AliDCSClient::GetValues(AliDCSMessage::RequestType reqType,
	const TSeqCollection* list, UInt_t startTime, UInt_t endTime,
	Int_t startIndex, Int_t endIndex)
{
	// get array of DCS values from the DCS server
	// list, startTime, endTime: list of dp/alias names, start time and end time of the query
	//
	// Result: map containing the array of alias names. It will be filled with
	// the values retrieved for each alias

	TMap* result = new TMap;
	result->SetOwner(1);

	if (endIndex < 0 || endIndex > list->GetEntries())
		endIndex = list->GetEntries();

	for (Int_t subsetBegin = startIndex; subsetBegin < endIndex; subsetBegin += fMultiSplit)
	{
    		Connect();
    	
	    	if (!IsConnected()) 
		{
    			AliError("Not connected!");
            		delete result;
			fResultErrorCode = fgkBadState;
	    		return 0;
	    	}

		Int_t subsetEnd = subsetBegin + fMultiSplit;
	        if (subsetEnd > endIndex)
        	    subsetEnd = endIndex;
	
		AliDCSMessage requestMessage;
	        if (fMultiSplit == 1)
        	{
	            // single dp request
            
	    		TObjString* aRequest = (TObjString*) list->At(subsetBegin);
        	  	requestMessage.CreateRequestMessage(reqType, startTime, endTime, aRequest->String());
	        }
        	else
	        {
        	    // multi dp request
            
	            requestMessage.CreateMultiRequestMessage(reqType,
			      	startTime, endTime);

	        	for (Int_t i=subsetBegin; i<subsetEnd; i++)
        		{
        			TObjString* aRequest = (TObjString*) list->At(i);
        			if (!requestMessage.AddRequestString(aRequest->String()))
				{
					delete result;
					fResultErrorCode = fgkBadMessage;
        				return 0;
				}
	        	}
        	}
	
    		if ((fResultErrorCode = SendMessage(requestMessage)) < 0)
    		{
                    AliError(Form("Can't send request message! Reason: %s",
                            GetErrorString(fResultErrorCode)));
                    Close();
                    delete result;
                    return 0;
	        }
    
    		for (Int_t i=subsetBegin; i<subsetEnd; i++)
	    	{
    			TObjString* aRequest = (TObjString*) list->At(i);
    		
    			TObjArray* resultSet = new TObjArray();
	    		resultSet->SetOwner(1);
    
    			if ((fResultErrorCode = ReceiveValueSet(resultSet)) < 0) 
	        	{
    				AliError(Form("Can't get values for %s!" ,
    				aRequest->String().Data()));
    
    				delete resultSet;
	    			result->DeleteValues();
    				delete result;
    				return 0;
		    	}

			result->Add(aRequest, resultSet);
		   }

	       Close();

 		   AliInfo(Form("Retrieved entries %d..%d (total %d..%d); E.g. %s has %d values collected",
					subsetBegin, subsetEnd, startIndex, endIndex, list->At(subsetBegin)->GetName(), ((TObjArray*)
					result->GetValue(list->At(subsetBegin)->GetName()))->GetEntriesFast()));
		
	    }

	return result;
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
TMap* AliDCSClient::GetDPValues(const TSeqCollection* dpList, UInt_t startTime, UInt_t endTime,
	Int_t startIndex, Int_t endIndex) 
{
	//
        // For every entry (fron startIndex to endIndex) in dpList (which must be TObjString) 
	// reads a valueSet. The key represents particular DataPoint to be read.
        // For all DataPoints time interval (startTime - endTime) is used.
        // After the read, the correspoding value for every key is a 
	// TObjArray - collection of AliDCSValue, or result is an empty map in
	// case of error.
	// 
        // Returns:
        //      TMap of results, 0 in case of failure
        //

	return GetValues(AliDCSMessage::kDPName, dpList, startTime, endTime, startIndex, endIndex);
}

//______________________________________________________________________
TMap* AliDCSClient::GetAliasValues(const TSeqCollection* aliasList, UInt_t startTime, UInt_t endTime,
	Int_t startIndex, Int_t endIndex)
{
	//
        // For every entry (fron startIndex to endIndex) in dpList (which must be TObjString) 
        // reads a valueSet. The key represents particular Alias to be read.
        // For all aliases time interval (startTime - endTime) is used.
        // After the read, the correspoding value for every key is a 
        // TObjArray - collection of AliDCSValue, or result is an empty map in
        // case of error.
        // 
        // Returns:
        //      TMap of results, 0 in case of failure
        //

	return GetValues(AliDCSMessage::kAlias, aliasList, startTime, endTime, startIndex, endIndex);
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
		delete fSocket;
		fSocket = 0;
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

