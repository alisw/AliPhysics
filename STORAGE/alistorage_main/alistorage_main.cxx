#include "AliStorageClientThread.h"
#include "AliStorageServerThread.h"

#include <iostream>

using namespace std;

void *ClientThreadHandle(void*)
{
	cout<<"ALICE Storage Manager -- Starting client thread"<<endl;
	AliStorageClientThread *client = new AliStorageClientThread();

	if(client)
	{
		client->CollectData();
	}
	else
	{
		cout<<"ALICE Storage Manager -- ERROR - failed to start client thread!!"<<endl;
	}
	
	if(client){delete client;}
	return 0;
}

void *ServerThreadHandle(void*)
{
	cout<<"\nALICE Storage Manager -- Starting server thread"<<endl;
	AliStorageServerThread *server = new AliStorageServerThread();

	if(!server)
	{
		cout<<"ALICE Storage Manager -- ERROR - failed to start server thread!!"<<endl;
	}
	if(server){delete server;}
	return 0;
}

int main()
{
	TThread *clientThread = new TThread("clientThread", ClientThreadHandle,NULL);
	
	TThread *serverThread = new TThread("serverThread", ServerThreadHandle,NULL);
	clientThread->Run();
	serverThread->Run();
	
	clientThread->Join();
	serverThread->Kill();//if client thread if finished, server thread is killed	
	return 0;
}
