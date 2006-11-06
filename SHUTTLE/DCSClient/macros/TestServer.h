#ifndef TEST_SERVER_H
#define TEST_SERVER_H

#include <TServerSocket.h>

class AliDCSMessage;

class TestServer: public TObject {
public:

        static const Int_t kBadState = -1;

        static const Int_t kTimeout = -2;

        static const Int_t kBadMessage = -3;

        static const Int_t kCommError = -4;

        static const Int_t kServerError = -5;

private:

        TServerSocket fServerSocket;

        Long_t fTimeout;

        Int_t fRetries;


        Int_t SendBuffer(TSocket* socket, const char* buffer, Int_t size);

        Int_t ReceiveBuffer(TSocket* socket, char* buffer, Int_t size);

        Int_t SendMessage(TSocket* socket, AliDCSMessage& message);

        Int_t ReceiveMessage(TSocket* socket, AliDCSMessage& message);

public:

	TestServer(Int_t port, Long_t timeout = 5000, Int_t retries = 5);

	void Run(Int_t count, Int_t rsSize);

	ClassDef(TestServer, 0);
};

#endif
