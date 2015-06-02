#ifndef AliStorageServerThread_H
#define AliStorageServerThread_H

#include "AliZMQManager.h"
#include "AliStorageDatabase.h"

#include <vector>
#include <string>
#include <set>

class AliStorageServerThread
{
public:
	AliStorageServerThread();
	~AliStorageServerThread();
private:
	void StartCommunication();
	bool MarkEvent(struct eventStruct event);
    std::vector<string100> GetTriggerClasses();
    
    
	//connection to database and storage
	AliStorageDatabase *fDatabase;
	std::string fStoragePath;

	AliStorageServerThread(const AliStorageServerThread&);
	AliStorageServerThread& operator=(const AliStorageServerThread&);
};

#endif
