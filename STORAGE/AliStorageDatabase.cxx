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
	ifstream configFile (Form("%s/STORAGE/setupStorageDatabase.sh",
				  gSystem->Getenv("ALICE_ROOT")));

	
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

AliStorageDatabase::~AliStorageDatabase(){}

void AliStorageDatabase::InsertEvent(int runNumber,
				     int eventNumber,
				     char *system,
				     int multiplicity,
				     char *filePath)
{
	fServer->Query(Form("replace into %s (run_number,event_number,system,multiplicity,permanent,file_path) values (%d,%d,'%s',%d,0,'%s');",fTable.c_str(),runNumber,eventNumber,system,multiplicity,filePath));

}

bool AliStorageDatabase::MarkEvent(struct eventStruct event)
{
	if(fServer->Query(Form("UPDATE %s SET permanent = 1 WHERE run_number = %d AND event_number = %d;",fTable.c_str(),event.runNumber,event.eventNumber)))
	{
		return 1;
	}
	else
	{
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
	/*
	//deserialize file??
	TTree* tree= new TTree("esdTree", "esdTree");
	data->WriteToTree(tree);
	tree-> Fill();
	AliESDEvent* requestedEvent= new AliESDEvent();	
	requestedEvent->ReadFromTree(tree);
	tree->GetEntry(0);
	delete data;
	delete tree;
	
	
	if(requestedEvent->GetRunNumber()<0)//if file is not in temp storage, check if it is in perm storage
	{
		cout<<"DATABASE -- could not find file in temp storage -- searching in perm"<<endl;
		TFile *permFile = new TFile(Form("%s/permEvents.root",fStoragePath.c_str()),"update");//open/create perm file
	
		if(!permFile)
		{
			cout<<"DATABASE -- Couldn't open perm file"<<endl;
			tmpFile->Close();
			if(tmpFile){delete tmpFile;}
			if(requestedEvent){delete requestedEvent;}
			return NULL;
		}
		TDirectory *runDirectory = permFile->GetDirectory(Form("run%d",event.runNumber));

		if(!runDirectory)
		{
			cout<<"DATABASE -- Couldn't open run directory"<<endl;
			
			if(tmpFile)
			{
				tmpFile->Close();
				delete tmpFile;
			}
			if(permFile)
			{
				permFile->Close();
				delete permFile;
			}
			if(requestedEvent){delete requestedEvent;}
			return NULL;
		}
		
		//if(requestedEvent){delete requestedEvent;}
		AliESDEvent *requestedEventPerm = (AliESDEvent*)runDirectory->Get(Form("event%d",event.eventNumber));

		if(requestedEventPerm->GetRunNumber()<0)
		{
			cout<<"DATABASE -- could not find event in perm storage"<<endl;
			tmpFile->Close();
			permFile->Close();
			if(tmpFile){delete tmpFile;}
			if(permFile){delete permFile;}
			return NULL;
		}
		else
		{
			return requestedEventPerm;
		}
	}
	else
	{
		cout<<"DATABASE -- sending event:"<<requestedEvent->GetRunNumber()<<endl;
		tmpFile->Close();
		if(tmpFile){delete tmpFile;}
		return requestedEvent;
		}*/
}

struct eventStruct AliStorageDatabase::GetOldestEvent()
{
	struct eventStruct oldestEvent = {0,0};

	TSQLResult *result = fServer->Query(Form("SELECT * FROM %s ORDER BY run_number,event_number;",fTable.c_str()));

	TSQLRow *row;	
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

void AliStorageDatabase::RemoveEvent(struct eventStruct event)
{
	fServer->Query(Form("DELETE FROM %s WHERE run_number = %d AND event_number = %d",fTable.c_str(),event.runNumber,event.eventNumber));
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
