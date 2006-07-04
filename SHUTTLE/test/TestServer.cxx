#include "TestServer.h"

#include "AliDCSMessage.h"
#include "AliLog.h"

#include <TTimeStamp.h>

ClassImp(TestServer)

const Int_t TestServer::kBadState;

const Int_t TestServer::kTimeout;

const Int_t TestServer::kBadMessage;

const Int_t TestServer::kCommError;

const Int_t TestServer::kServerError;

TestServer::TestServer(Int_t port, Long_t timeout, Int_t retries):
	fServerSocket(port), fTimeout(timeout), fRetries(retries)
{

}

Int_t TestServer::SendBuffer(TSocket* socket, const char* buffer, Int_t size) {

        Int_t sentSize = 0;
        Int_t tries = 0;

        while (sentSize < size && tries < fRetries) {

                Int_t sResult = socket->Select(TSocket::kWrite, fTimeout);

                if (sResult == 0) {
                        AliDebug(1, Form("Timeout! tries <%d> ...", tries));
                        tries ++;
                        continue;

                } else if (sResult < 0) {
                        AliDebug(1, Form("Communication error <%d>!",
                                        socket->GetErrorCode()));
                        return TestServer::kCommError;
                }

                sResult = socket->SendRaw(buffer + sentSize, size - sentSize,
                                        kDontBlock);

                if (sResult > 0) {
                        sentSize += sResult;
                } else {
                        AliDebug(1, Form("Communication error <%d>!",
                                        socket->GetErrorCode()));
                        return TestServer::kCommError;
                } 
        }

        if (tries == fRetries) {
                return TestServer::kTimeout;
        }

        return sentSize;
}

Int_t TestServer::ReceiveBuffer(TSocket* socket, char* buffer, Int_t size) {

        Int_t readSize = 0;
        Int_t tries = 0;

        while (readSize < size && tries < fRetries) {

                Int_t sResult = socket->Select(TSocket::kRead, fTimeout);

                if (sResult == 0) {
                        AliDebug(1, Form("Timeout! tries <%d> ...", tries));
                        tries ++;
                        continue;

                } else if (sResult < 0) {
                        AliDebug(1, Form("Communication error <%d>!",
                                        socket->GetErrorCode()));
                        return TestServer::kCommError;
                }

                sResult = socket->RecvRaw(buffer + readSize, size - readSize,
                                         kDontBlock);

                if (sResult > 0) {
                        readSize += sResult;
                } else {
                        AliDebug(1, Form("Communication error <%d>!",
                                        socket->GetErrorCode()));
                        return TestServer::kCommError;
                }
        }

        if (tries == fRetries) {
                return TestServer::kTimeout;
        }

        return readSize;
}

Int_t TestServer::SendMessage(TSocket* socket, AliDCSMessage& message) {
	
	message.StoreToBuffer();
	
        return SendBuffer(socket, message.GetMessage(), 
				message.GetMessageSize());
}

Int_t TestServer::ReceiveMessage(TSocket* socket, AliDCSMessage& message) {

	char header[HEADER_SIZE];

        Int_t sResult;

        if ((sResult = ReceiveBuffer(socket, header, HEADER_SIZE)) < 0) {
                return sResult;
        }

        if (!message.SetRawHeader(header)) {
                return TestServer::kBadMessage;
        }

        if ((sResult = ReceiveBuffer(socket, message.GetBody(),
                        message.GetBodySize())) < 0) {

                return sResult;
        }

	message.LoadFromBuffer();	

        return HEADER_SIZE + sResult;
}


void TestServer::Run(Int_t count, Int_t rsSize) {

        if (!fServerSocket.IsValid()) {
                AliError("Invalid server socket!");
                return;
        }


        TSocket* connSocket;
	while (1) {
		
        	AliInfo("Waiting for connection");

	        connSocket = fServerSocket.Accept();

		AliInfo("Client accepted");

	        if (!connSocket) {
        	        AliError("Can't accept connection!");
                	continue;
	        }

        	connSocket->SetOption(kNoBlock, 1);

        	if (!connSocket->IsValid()) {
                	AliError("Invalid connection socket!");
			delete connSocket;
	                continue;
        	}

		Int_t result;

        	AliDCSMessage message;
	        if ((result = ReceiveMessage(connSocket, message)) < 0) {
        	        AliError(Form("Communication failure: ", result));
			delete connSocket;
			continue;
	        }
		
        	//message.Print();

		if (message.GetType() != AliDCSMessage::kRequest) {
			AliError("Bad message type!");
			delete connSocket;
			continue;
		}

		AliDCSMessage countMessage;
		countMessage.CreateCountMessage(count);
		if ((result = SendMessage(connSocket, countMessage)) < 0) {
			AliError(Form("Communication failure: %d", result));
			continue;
		}

			
		Int_t sentValues = 0;
		
		AliInfo(Form("Count: %d", count));

		while (sentValues < count) {
			
			Int_t pSize = rsSize < count - sentValues ? 
					rsSize : count - sentValues;


			TTimeStamp currentTime;
			AliDCSMessage rsMessage;
			rsMessage.CreateResultSetMessage(AliDCSValue::kInt);
				
			for (Int_t k = 0; k < pSize; k ++) {
				AliDCSValue aValue(k, currentTime.GetSec() + k);
				rsMessage.AddValue(aValue);
			}

			if ((result = SendMessage(connSocket, rsMessage)) < 0) {
				AliError(Form("Communication failure: %d", 
						result));
				break;
			}

			sentValues += pSize;
		}
		
		AliInfo(Form("Sent values: %d", sentValues));		
		delete connSocket; 

		//break;
	}
}
