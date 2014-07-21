#ifndef AliStorageServerThread_H
#define AliStorageServerThread_H

#include "AliStorageEventManager.h"
#include "AliStorageDatabase.h"

class AliStorageServerThread
{
public:
	AliStorageServerThread();
	~AliStorageServerThread();
private:
	//communication via socket
	void StartCommunication();
	AliStorageEventManager *fEventManager;
	
	//handling different requests
	AliESDEvent* GetEvent(struct eventStruct eventStruct);
	AliESDEvent* GetNextEvent(struct eventStruct eventStruct);
	AliESDEvent* GetLastEvent();
	bool MarkEvent(struct eventStruct event);

	//connection to database and storage
	AliStorageDatabase *fDatabase;
	std::string fStoragePath;

	AliStorageServerThread(const AliStorageServerThread&);
	AliStorageServerThread& operator=(const AliStorageServerThread&);
};

#endif
