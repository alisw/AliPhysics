#ifndef AliStorageTypes_H
#define AliStorageTypes_H

#include <TSystem.h>
#include <TFile.h>

inline const char* GetConfigFilePath(){
    return Form("%s/setupStorageDatabase.sh",gSystem->Getenv("HOME"));
}

enum storageSockets{
	SERVER_COMMUNICATION_REQ=0,
	SERVER_COMMUNICATION_REP,
	CLIENT_COMMUNICATION_REQ,
	CLIENT_COMMUNICATION_REP,
	EVENTS_SERVER_PUB,
	EVENTS_SERVER_SUB,
	XML_PUB,
	ITS_POINTS_PUB,
	ITS_POINTS_SUB,
	NUMBER_OF_SOCKETS
};
	
enum statusType{
	STATUS_WAITING=1,
	STATUS_OK,
	STATUS_ERROR,
	STATUS_DOWN
};

enum requestType{
	REQUEST_CONNECTION=1,
	REQUEST_RECEIVING,
	REQUEST_SAVING,
	REQUEST_CURRENT_SIZE,
	REQUEST_LIST_EVENTS,
	REQUEST_GET_EVENT,
	REQUEST_GET_NEXT_EVENT,
	REQUEST_GET_PREV_EVENT,
	REQUEST_GET_LAST_EVENT,
	REQUEST_GET_FIRST_EVENT,
	REQUEST_MARK_EVENT,
	REQUEST_SET_PARAMS,
	REQUEST_GET_PARAMS,
    REQUEST_GET_TRIGGER_LIST
};

struct clientRequestStruct
{
    clientRequestStruct() :
    messageType(-1),
    maxStorageSize(-1),
    maxOccupation(-1),
    removeEvents(-1),
    eventsInChunk(-1)
    {};
    clientRequestStruct(const clientRequestStruct& crs){
        messageType = crs.messageType;
        maxStorageSize =crs.maxStorageSize;
        maxOccupation = crs.maxOccupation;
        removeEvents = crs.removeEvents;
        eventsInChunk = crs.eventsInChunk;
    }
    
	int messageType;
	int maxStorageSize;
	int maxOccupation;
	int removeEvents;
	int eventsInChunk;
};

struct eventStruct{
	int runNumber;
	int eventNumber;
};

struct listRequestStruct{
	int runNumber[2];
	int eventNumber[2];
	int marked[2];
	int multiplicity[2];
	char system[2][20];
    char triggerClass[100];
};

struct serverRequestStruct
{
    serverRequestStruct():
    messageType(-1),
    eventsRunNumber(-1),
    eventsEventNumber(-1),
    runNumber(),
    eventNumber(),
    marked(),
    multiplicity(),
    system(),
    triggerClass()
    {};
    serverRequestStruct(const serverRequestStruct& src){
        messageType = src.messageType;
        eventsRunNumber = src.eventsRunNumber;
        eventsEventNumber = src.eventsEventNumber;
        runNumber[0] = src.runNumber[0];
        runNumber[1] = src.runNumber[1];
        eventNumber[0] = src.eventNumber[0];
        eventNumber[1] = src.eventNumber[1];
        marked[0] = src.marked[0];
        marked[1] = src.marked[1];
        multiplicity[0] = src.multiplicity[0];
        multiplicity[1] = src.multiplicity[1];
        strcpy(system[0],src.system[0]);
        strcpy(system[1],src.system[1]);
        strcpy(triggerClass,src.triggerClass);
    }
	int messageType;
    
    int eventsRunNumber;
    int eventsEventNumber;

    int runNumber[2];
    int eventNumber[2];
    int marked[2];
    int multiplicity[2];
    char system[2][20];
    char triggerClass[100];
    
//	struct eventStruct event;
//	struct listRequestStruct list;
};

struct recPointsStruct{
  TFile *files[10];
};

typedef struct serverListStruct{
	int runNumber;
	int eventNumber;
	char system[20];
	int multiplicity;
	int marked;
    char triggerClass[100];
}serverListStruct;

struct string100 {
    char data[100];
    string100(const char* init){strcpy(data,init);}
};

#endif
