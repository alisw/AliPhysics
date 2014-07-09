#ifndef AliStorageTypes_H
#define AliStorageTypes_H

//sockets
const int gServerCommunicationPort = 5066;
const int gEventsSubscriberPort = 5024;
const int gClientCommunicationPort = 5088;
const int gClientPublisherPort = 5077;

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
	REQUEST_GET_LAST_EVENT,
	REQUEST_MARK_EVENT,
	REQUEST_SET_PARAMS,
	REQUEST_GET_PARAMS
};

struct clientRequestStruct{
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
};

struct serverRequestStruct{
	int messageType;
	struct eventStruct event;
	struct listRequestStruct list;
};

typedef struct serverListStruct{
	int runNumber;
	int eventNumber;
	char system[20];
	int multiplicity;
	int marked;
}serverListStruct;

#endif
