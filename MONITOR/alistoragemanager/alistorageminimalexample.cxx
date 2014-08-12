#ifdef ZMQ
#include "AliStorageEventManager.h"
#endif
#include <iostream>
#include <TFile.h>


int main()
{
#ifdef ZMQ
	AliStorageEventManager *manager = AliStorageEventManager::GetEventManagerInstance();
	manager->CreateSocket(EVENTS_SERVER_SUB);
	AliESDEvent *event  = manager->GetEvent(EVENTS_SERVER_SUB);
		
	std::cout<<"Received event:"<<event->GetEventNumberInFile()<<std::endl;
#endif
	return 0;
}
