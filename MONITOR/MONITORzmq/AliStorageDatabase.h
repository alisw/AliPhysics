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
			 char *filePath,
             ULong64_t triggerMask,
             ULong64_t triggerMaskNext50); //more parameters of the event can be added to this method

	bool MarkEvent(struct eventStruct event);
	bool UpdateEventPath(struct eventStruct event,const char *newPath);
	void RemoveEvent(struct eventStruct event);
    void RemoveEventsWithPath(std::string path);
	std::string GetFilePath(struct eventStruct event);
	struct eventStruct GetOldestEvent();
	std::vector<serverListStruct> GetList(struct listRequestStruct listStruct);
    std::vector<int> GetListOfRuns();   // returns list of run numbers stored in the database
    
	AliESDEvent* GetEvent(struct eventStruct event);
	AliESDEvent* GetNextEvent(struct eventStruct event);
	AliESDEvent* GetLastEvent();
	AliESDEvent* GetPrevEvent(struct eventStruct event);
	AliESDEvent* GetFirstEvent();
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
