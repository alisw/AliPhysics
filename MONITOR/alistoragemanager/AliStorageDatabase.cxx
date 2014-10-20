#include "AliStorageDatabase.h"

#include <iostream>
#include <fstream>

#include <TSQLRow.h>
#include <TSQLResult.h>
#include <TThread.h>
#include <TSystem.h>
#include <TFile.h>

using namespace std;

AliStorageDatabase::AliStorageDatabase() :
	fHost(""),
	fPort(""),
	fDatabase(""),
	fUID(""),
	fPassword(""),
	fTable(""),
	fServer(0),
	fStoragePath("")
{
	TThread::Lock();	
	ifstream configFile (GetConfigFilePath());

	
	if (configFile.is_open())
	{
		string line;
		int from,to;
		while(configFile.good())
		{
			getline(configFile,line);
			from = line.find("\"")+1;
			to = line.find_last_of("\"");
			if(line.find("HOST=")==0)
			{
				fHost=line.substr(from,to-from);
			}
			else if(line.find("PORT=")==0)
			{
				fPort=line.substr(from,to-from);
			}
			else if(line.find("DATABASE=")==0)
			{
				fDatabase=line.substr(from,to-from);
			}
			else if(line.find("USER=")==0)
			{
				fUID=line.substr(from,to-from);
			}
			else if(line.find("PASS=")==0)
			{
				fPassword=line.substr(from,to-from);
			}
			else if(line.find("TABLE=")==0)
			{
				fTable=line.substr(from,to-from);
			}
			else if(line.find("STORAGE_PATH=")==0)
			{
				fStoragePath=line.substr(from,to-from);
			}

		}
		if(configFile.eof())
		{
			configFile.clear();
		}
		configFile.close();
	}
	else
	{
		cout << "DATABASE -- Unable to open file" <<endl;
	}
	TThread::UnLock();

	
	fServer = TSQLServer::Connect(Form("mysql://%s:%s/%s",fHost.c_str(),fPort.c_str(),fDatabase.c_str()),fUID.c_str(),fPassword.c_str());
}

AliStorageDatabase::~AliStorageDatabase(){
  if (fServer) {delete fServer;}
}

void AliStorageDatabase::InsertEvent(int runNumber,
				     int eventNumber,
				     char *system,
				     int multiplicity,
				     char *filePath)
{
  TSQLResult* res;
  res = fServer->Query(Form("replace into %s (run_number,event_number,system,multiplicity,permanent,file_path) values (%d,%d,'%s',%d,0,'%s');",fTable.c_str(),runNumber,eventNumber,system,multiplicity,filePath));
  delete res;

}

bool AliStorageDatabase::MarkEvent(struct eventStruct event)
{  
  TSQLResult* res;
  res = fServer->Query(Form("UPDATE %s SET permanent = 1 WHERE run_number = %d AND event_number = %d;",fTable.c_str(),event.runNumber,event.eventNumber));
  if(res) {
    delete res;
    return 1;
  }
  else {
    delete res;
    return 0;
  }
}

vector<serverListStruct> AliStorageDatabase::GetList(struct listRequestStruct list)
{
	TSQLResult *result = fServer->Query(Form("SELECT * FROM %s WHERE run_number >= %d AND run_number <= %d AND event_number >= %d AND event_number <= %d AND multiplicity >= %d AND multiplicity <= %d AND (permanent = %d OR permanent = %d) AND (system = '%s' OR system = '%s') ORDER BY run_number,event_number;",
						 fTable.c_str(),
						 list.runNumber[0],
						 list.runNumber[1],
						 list.eventNumber[0],
						 list.eventNumber[1],
						 list.multiplicity[0],
						 list.multiplicity[1],
						 list.marked[0],
						 list.marked[1],
						 list.system[0],
						 list.system[1]));

	
	TSQLRow *row;
	vector<serverListStruct> eventsVector;
	
	while((row = result->Next()))
	{
		serverListStruct resultList;

		resultList.runNumber = atoi(row->GetField(0));
		resultList.eventNumber = atoi(row->GetField(1));
		strcpy(resultList.system, row->GetField(2));
		resultList.multiplicity = atoi(row->GetField(3));
		resultList.marked = atoi(row->GetField(4));

		eventsVector.push_back(resultList);
		delete row;
	}
	delete result;
	return eventsVector;
}

AliESDEvent* AliStorageDatabase::GetEvent(struct eventStruct event)
{
	cout<<"database - get event"<<endl;
	string pathToFile = GetFilePath(event);

	if(!strcmp(pathToFile.c_str(),""))
	{
		cout<<"DATABASE -- no such file in database"<<endl;
		return NULL;
	}
	
	TFile *tmpFile = new TFile(pathToFile.c_str(),"read");
	if(!tmpFile)
	{
		cout<<"DATABASE -- couldn't open temp file"<<endl;
		return NULL;
	}
	AliESDEvent *data;
	tmpFile->GetObject(Form("event%d;1",event.eventNumber),data);

	return data;
}

void AliStorageDatabase::RemoveEvent(struct eventStruct event)
{
  TSQLResult* res;
  res = fServer->Query(Form("DELETE FROM %s WHERE run_number = %d AND event_number = %d",fTable.c_str(),event.runNumber,event.eventNumber));
  delete res;
}

void AliStorageDatabase::RemoveEventsWithPath(string path)
{
    TSQLResult *res = fServer->Query(Form("DELETE FROM %s WHERE file_path = \"%s\";",fTable.c_str(),path.c_str()));
    delete res;
}

string AliStorageDatabase::GetFilePath(struct eventStruct event)
{
	TSQLResult *result = fServer->Query(Form("SELECT * FROM %s WHERE run_number = %d AND event_number = %d;",fTable.c_str(),event.runNumber,event.eventNumber));
	TSQLRow *row;
	row = result->Next();
	if(row)
	{
		string path(row->GetField(5));
		delete row;
		return path;
	}
	else
	{
		return "";
	}
}

AliESDEvent* AliStorageDatabase::GetNextEvent(struct eventStruct event)
{
    cout<<"Database:"<<event.runNumber<<"\t"<<event.eventNumber<<endl;
    
	TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number;",fTable.c_str()));

	TSQLRow *row;
	bool isCurrentEvent=false;
	struct eventStruct nextEvent;
	
	while((row = result->Next()))
	{
		if(isCurrentEvent)
		{
			nextEvent.runNumber = atoi(row->GetField(0));
			nextEvent.eventNumber = atoi(row->GetField(1));
			return GetEvent(nextEvent);
		}

		//if current event found
		if(atoi(row->GetField(0))==event.runNumber && atoi(row->GetField(1))==event.eventNumber)
		{
			isCurrentEvent=true;
		}
		
		delete row;
	}

	return NULL;
}

AliESDEvent* AliStorageDatabase::GetPrevEvent(struct eventStruct event)
{
	TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number DESC;",fTable.c_str()));

	TSQLRow *row;
	bool isCurrentEvent=false;
	struct eventStruct nextEvent;
	
	while((row = result->Next()))
	{
		if(isCurrentEvent)
		{
			nextEvent.runNumber = atoi(row->GetField(0));
			nextEvent.eventNumber = atoi(row->GetField(1));
			return GetEvent(nextEvent);
		}

		//if current event found
		if(atoi(row->GetField(0))==event.runNumber && atoi(row->GetField(1))==event.eventNumber)
		{
			isCurrentEvent=true;
		}
		
		delete row;
	}
	delete result;
	return NULL;
}

struct eventStruct AliStorageDatabase::GetOldestEvent()
{
    TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number;",fTable.c_str()));
    
    TSQLRow *row;
    struct eventStruct oldestEvent = {0,0};
    
    if((row = result->Next()))
    {
        oldestEvent.runNumber = atoi(row->GetField(0));
        oldestEvent.eventNumber = atoi(row->GetField(1));
        delete row;
    }
    else
    {
        cout<<"DATABASE -- NO OLDEST EVENT FOUND. Storage may be corrupted."<<endl;
    }
    return oldestEvent;
}

AliESDEvent* AliStorageDatabase::GetLastEvent()
{
	TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number;",fTable.c_str()));

	TSQLRow *row;
	struct eventStruct lastEvent = {0,0};

	while((row = result->Next()))
	{
		lastEvent.runNumber = atoi(row->GetField(0));
		lastEvent.eventNumber = atoi(row->GetField(1));
		delete row;
	}
	cout<<"Last event is:"<<lastEvent.eventNumber<<endl;
	return GetEvent(lastEvent);

}

AliESDEvent* AliStorageDatabase::GetFirstEvent()
{
    cout<<"Database - first"<<endl;
	TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number DESC;",fTable.c_str()));

	TSQLRow *row;
	struct eventStruct firstEvent = {0,0};

	while((row = result->Next()))
	{
		firstEvent.runNumber = atoi(row->GetField(0));
		firstEvent.eventNumber = atoi(row->GetField(1));
		delete row;
	}
	cout<<"First event is:"<<firstEvent.eventNumber<<endl;
	return GetEvent(firstEvent);

}
