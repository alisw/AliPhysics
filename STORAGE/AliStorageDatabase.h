#ifndef AliStorageDatabase_H
#define AliStorageDatabase_H

#include "AliStorageTypes.h"
#include "AliESDEvent.h"

#include <string>
#include <vector>

#include <TSQLServer.h>

class AliStorageDatabase
{
public:
	AliStorageDatabase();
	~AliStorageDatabase();
	
	void InsertEvent(int runNumber,
			 int eventNumber,
			 char *system,
			 int multiplicity,
			 char *filePath); //more parameters of the event can be added to this method

	bool MarkEvent(struct eventStruct event);
	void RemoveEvent(struct eventStruct event);
	std::string GetFilePath(struct eventStruct event);
	struct eventStruct GetOldestEvent();
	std::vector<serverListStruct> GetList(struct listRequestStruct listStruct);

	//not tested:
	AliESDEvent* GetEvent(struct eventStruct event);
	AliESDEvent* GetNextEvent(struct eventStruct event);
	AliESDEvent* GetLastEvent();
	//
private:
	std::string fHost;
	std::string fPort;
	std::string fDatabase;
	std::string fUID;
	std::string fPassword;
	std::string fTable;

	TSQLServer *fServer;
	std::string fStoragePath;

	 AliStorageDatabase(const AliStorageDatabase&);
	 AliStorageDatabase& operator=(const AliStorageDatabase&);
};

#endif
