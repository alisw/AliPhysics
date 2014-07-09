#include "AliStorageEventManager.h"
#include "AliStorageTypes.h"
#include "AliESDEvent.h"

#include "zmq.hpp"
#include <iostream>

using namespace std;
using namespace zmq;

/* parameters:

0 - connect directly to reconstruction socket
1 - connect to Storage Manager and ask for last event

*/

int main(int argc, char **argv)
{
	if(argc<2)
	{
		cout<<"Usage: alifakedisplay <mode>"<<endl;
		cout<<"mode:"<<endl;
		cout<<"0 - connect directly to reconstruction socket"<<endl;
		cout<<"1 - connect to Storage Manager and ask for last event"<<endl;
		return 0;
	}
	
	context_t *context = new context_t();
	socket_t *socket;
	AliStorageEventManager *manager = new AliStorageEventManager();
	AliESDEvent *event;
	
	if(atoi(argv[1])==0)
	{
		socket = new socket_t(*context,ZMQ_SUB);
		socket->setsockopt(ZMQ_SUBSCRIBE,"",0);
		socket->connect(Form("tcp://137.138.93.150:%d",gEventsSubscriberPort));
		while(1)
		{
			event = manager->GetEvent(socket);
			if(event)
			{
				cout<<"Received event. Run:"<<event->GetRunNumber()<<"\t event:"<<event->GetEventNumberInFile()<<endl;
		
				delete event;
			}
			else
			{
				cout<<"NO EVENT"<<endl;
			}
		}
	}
	else if(atoi(argv[1])==1)
	{
		socket = new socket_t(*context,ZMQ_REQ);
		socket->connect(Form("tcp://137.138.93.150:%d",gServerCommunicationPort));
		while(1)
		{
			struct serverRequestStruct *requestMessage = new struct serverRequestStruct;
			requestMessage->messageType = REQUEST_GET_LAST_EVENT;

			manager->Send(requestMessage,socket);			
			event = manager->GetEvent(socket);
			if(event)
			{
				cout<<"Last event - Run:"<<event->GetRunNumber()<<"\t event:"<<event->GetEventNumberInFile()<<endl;
		
				delete event;
			}
			else
			{
				cout<<"NO EVENT"<<endl;
			}
			sleep(1);
		}
	}

	return 0;
}
