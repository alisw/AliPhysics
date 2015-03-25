#include "AliStorageEventManager.h"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fstream>



#include <TList.h>
#include <TStreamerInfo.h>
#include <TThread.h>

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTrackPointArray.h"
#include "AliESDfriendTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliTracker.h"

using namespace std;

AliStorageEventManager *AliStorageEventManager::fManagerInstance = 0;

AliStorageEventManager::AliStorageEventManager()
{
    fZmqError = false;
    
    //read config file
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
            if(line.find("STORAGE_SERVER=")==0)
            {
                fStorageServer=line.substr(from,to-from);
            }
            else if(line.find("EVENT_SERVER=")==0)
            {
                fEventServer=line.substr(from,to-from);
            }
            else if(line.find("STORAGE_SERVER_PORT=")==0)
            {
                fStorageServerPort=atoi(line.substr(from,to-from).c_str());
            }
            else if(line.find("EVENT_SERVER_PORT=")==0)
            {
                fEventServerPort=atoi(line.substr(from,to-from).c_str());
            }
            else if(line.find("STORAGE_CLIENT_PORT=")==0)
            {
                fStorageClientPort=atoi(line.substr(from,to-from).c_str());
            }
            else if(line.find("XML_SERVER_PORT=")==0)
            {
                fXmlServerPort=atoi(line.substr(from,to-from).c_str());
            }
            else if(line.find("ITS_POINTS_SERVER_PORT=")==0)
            {
                fItsPointsServerPort=atoi(line.substr(from,to-from).c_str());
                cout<<"ITS port is:"<<fItsPointsServerPort<<endl;
            }
        }
        if(configFile.eof()){configFile.clear();}
        configFile.close();
    }
    else{cout<<"EVENT MANAGER -- Unable to open config file"<<endl;}
    TThread::UnLock();
    
    for(int i=0;i<NUMBER_OF_SOCKETS;i++){fContexts[i] = zmq_ctx_new();}
}
AliStorageEventManager::~AliStorageEventManager()
{
    if(fManagerInstance){delete fManagerInstance;fManagerInstance=0;}
}

AliStorageEventManager* AliStorageEventManager::GetEventManagerInstance()
{
    TThread::Lock();
    if(fManagerInstance==0)
    {
        fManagerInstance = new AliStorageEventManager();
    }
    TThread::UnLock();
    return fManagerInstance;
}

bool AliStorageEventManager::CreateSocket(storageSockets socket)
{
    cout<<"Creating socket:"<<socket<<endl;
    int connectionStatus=-1,optionStatus=-1;
    
    switch(socket)
    {
        case SERVER_COMMUNICATION_REQ:
        {
            fSockets[SERVER_COMMUNICATION_REQ] = zmq_socket(fContexts[SERVER_COMMUNICATION_REQ],ZMQ_REQ);
            connectionStatus = zmq_connect(fSockets[SERVER_COMMUNICATION_REQ],
                                           Form("tcp://%s:%d",fStorageServer.c_str(),fStorageServerPort));
            break;
        }
        case SERVER_COMMUNICATION_REP:
        {
            fSockets[SERVER_COMMUNICATION_REP] = zmq_socket(fContexts[SERVER_COMMUNICATION_REP],ZMQ_REP);
            connectionStatus = zmq_bind(fSockets[SERVER_COMMUNICATION_REP],Form("tcp://*:%d",fStorageServerPort));
            break;
        }
        case CLIENT_COMMUNICATION_REQ:
        {
            fSockets[CLIENT_COMMUNICATION_REQ] = zmq_socket(fContexts[CLIENT_COMMUNICATION_REQ],ZMQ_REQ);
            connectionStatus = zmq_connect(fSockets[CLIENT_COMMUNICATION_REQ],Form("tcp://%s:%d",fStorageServer.c_str(), fStorageClientPort));
            break;
        }
        case CLIENT_COMMUNICATION_REP:
        {
            fSockets[CLIENT_COMMUNICATION_REP] = zmq_socket(fContexts[CLIENT_COMMUNICATION_REP],ZMQ_REP);
            connectionStatus = zmq_bind(fSockets[CLIENT_COMMUNICATION_REP],Form("tcp://*:%d",fStorageClientPort));
            break;
        }
        case EVENTS_SERVER_PUB:
        {
            fSockets[EVENTS_SERVER_PUB] = zmq_socket(fContexts[EVENTS_SERVER_PUB],ZMQ_PUB);
            connectionStatus = zmq_bind(fSockets[EVENTS_SERVER_PUB],Form("tcp://*:%d",fEventServerPort));
            break;
        }
        case EVENTS_SERVER_SUB:
        {
            fSockets[EVENTS_SERVER_SUB] = zmq_socket(fContexts[EVENTS_SERVER_SUB],ZMQ_SUB);
            if(zmq_setsockopt(fSockets[EVENTS_SERVER_SUB],ZMQ_SUBSCRIBE,"",0)!=0)
            {
                cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
                return 0;
            }
            connectionStatus = zmq_connect(fSockets[EVENTS_SERVER_SUB],Form("tcp://%s:%d",fEventServer.c_str(),fEventServerPort));
            break;
        }
        case XML_PUB:
        {
            fSockets[XML_PUB] = zmq_socket(fContexts[XML_PUB],ZMQ_PUB);
            connectionStatus = zmq_bind(fSockets[XML_PUB],Form("tcp://*:%d",fXmlServerPort));
            break;
        }
        case ITS_POINTS_PUB:
        {
            fSockets[ITS_POINTS_PUB] = zmq_socket(fContexts[ITS_POINTS_PUB],ZMQ_PUB);
            connectionStatus = zmq_bind(fSockets[ITS_POINTS_PUB],Form("tcp://*:%d",fItsPointsServerPort));
            break;
        }
        case ITS_POINTS_SUB:
        {
            fSockets[ITS_POINTS_SUB] = zmq_socket(fContexts[ITS_POINTS_SUB],ZMQ_SUB);
            if(zmq_setsockopt(fSockets[ITS_POINTS_SUB],ZMQ_SUBSCRIBE,"",0)!=0)
            {
                cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
                return 0;
            }
            connectionStatus = zmq_connect(fSockets[ITS_POINTS_SUB],Form("tcp://%s:%d",fEventServer.c_str(),fItsPointsServerPort));
            break;
        }
        default:break;
    }
    
    if (connectionStatus!=0) {
        cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
        return 0;
    }
    return 1;
}

void AliStorageEventManager::Send(vector<serverListStruct> list,storageSockets socket)
{
    //send size of the struct first
    int numberOfRecords = list.size();
    zmq_msg_t buffer;
    
    zmqInit(&buffer,sizeof(int));
    memcpy(zmq_msg_data(&buffer),&numberOfRecords,sizeof(int));
    
    zmqSend(&buffer,fSockets[socket],0);
    zmqRecv(&buffer,fSockets[socket],0);
    zmq_msg_close(&buffer);
    
    zmq_msg_t message;
    zmqInit(&message,sizeof(serverListStruct)*numberOfRecords);
    memcpy(zmq_msg_data(&message),reinterpret_cast<void*> (&list[0]), sizeof(serverListStruct)*numberOfRecords);
    
    zmqSend(&message,fSockets[socket],0);
    zmq_msg_close(&message);
}

bool AliStorageEventManager::Send(struct serverRequestStruct *request,storageSockets socket,int timeout)
{
    //check timeout
    fZmqError=false; // clean old errors
    zmqPoll(fSockets[socket],timeout);
    if(fZmqError){return 0;}
    
    
    int sizeOfRequest = sizeof(struct serverRequestStruct)+sizeof(struct listRequestStruct)+sizeof(struct eventStruct);
    
    zmq_msg_t buffer;
    zmqInit(&buffer,sizeOfRequest);
    memcpy(zmq_msg_data(&buffer),request,sizeOfRequest);
    zmqSend(&buffer,fSockets[socket],0);
    zmq_msg_close(&buffer);
    
    if(fZmqError){return 0;} // if there were errors while sending message
    return 1;
}

bool AliStorageEventManager::Send(struct clientRequestStruct *request,storageSockets socket,int timeout)
{
    //check timeout
    fZmqError=false;
    zmqPoll(fSockets[socket],timeout);
    if(fZmqError){return 0;}
    
    //put clientRequestStruct in buffer
    zmq_msg_t buffer;
    zmqInit(&buffer,sizeof(struct clientRequestStruct));
    memcpy(zmq_msg_data(&buffer),request,sizeof(struct clientRequestStruct));
    
    //send buffer
    zmqSend(&buffer,fSockets[socket],0);
    zmq_msg_close(&buffer);
    
    if(fZmqError) // if there were errors while sending message
    {
        if(fSockets[socket]){zmq_close(fSockets[socket]);fSockets[socket]=0;}
        CreateSocket(socket);
        return 0;
    }
    return 1;
}
void AliStorageEventManager::Send(long message,storageSockets socket)
{
    zmq_msg_t buffer;
    zmqInit(&buffer,sizeof(long));
    memcpy(zmq_msg_data(&buffer),&message,sizeof(long));
    zmqSend(&buffer,fSockets[socket],0);
    zmq_msg_close(&buffer);
}
void AliStorageEventManager::Send(bool message,storageSockets socket)
{
    zmq_msg_t buffer;
    zmqInit(&buffer,sizeof(bool));
    memcpy(zmq_msg_data(&buffer),&message,sizeof(bool));
    zmqSend(&buffer,fSockets[socket],0);
    zmq_msg_close(&buffer);
}
void AliStorageEventManager::Send(AliESDEvent *event, storageSockets socket)
{
    if(!event){return;}
    
    TMessage tmess(kMESS_OBJECT);
    tmess.Reset();
    tmess.WriteObject(event);
    TMessage::EnableSchemaEvolutionForAll(kTRUE);
    
    int bufsize = tmess.BufferSize();
    zmq_msg_t buffer;
    zmqInit(&buffer,bufsize);
    memcpy(zmq_msg_data(&buffer),tmess.Buffer(),bufsize);
    zmqSend(&buffer,fSockets[socket],0);
    zmq_msg_close(&buffer);
}

void AliStorageEventManager::SendAsXml(AliESDEvent *event,storageSockets socket)
{
    // to be tested from online reconstruction !!
    
    cout<<"SENDING AS XML"<<endl;
    stringstream bufferStream;
    bufferStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>"<<endl;
    bufferStream << "<ESD xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\">" << endl;
    
    for(int i=0;i<event->GetNumberOfTracks();i++)
    {
        AliESDtrack *track = event->GetTrack(i);
        AliKalmanTrack *ITStrack = track->GetITStrack();
        const AliTrackPointArray *array = track->GetTrackPointArray();
        
        bufferStream << "\t<track mass=\""<<track->GetMass()<<"\"";
        //	bufferStream << "\t<track esdpid=\""<<track->GetESDpid();
        bufferStream << "\t pid=\""<<track->PID()<<"\"";
        bufferStream << "\t energy=\""<<track->E()<<"\"";
        bufferStream << "\t volumeID=\""<<array->GetVolumeID()<<"\">" <<endl;
 
        if(array)
        {
            const float *x = array->GetX();
            const float *y = array->GetY();
            const float *z = array->GetZ();
            int n = array->GetNPoints();
            
            for(int j=0;j<n;j++)
            {
                bufferStream <<"\t\t<point volumeID=\""<<array->GetVolumeID()[j]<<"\">"<<endl;
                bufferStream <<"\t\t\t<x>"<< x[j] <<"</x>"<<endl;
                bufferStream <<"\t\t\t<y>"<< y[j] <<"</y>"<<endl;
                bufferStream <<"\t\t\t<z>"<< z[j] <<"</z>"<<endl;
                bufferStream <<"\t\t</point>"<<endl;
            }
        }
        else cout<<"no array"<<endl;
        
        bufferStream << "\t</track>"<<endl;
    }
    
    bufferStream << "</ESD>"<<endl;
    
    string bufferString = bufferStream.str();
    
    zmq_msg_t buffer;
    zmqInit(&buffer,bufferString.size());
    memcpy (zmq_msg_data(&buffer), bufferString.data(), bufferString.size());
    
    zmqSend(&buffer,fSockets[socket],0);
    zmq_msg_close(&buffer);
    cout<<"xml sent"<<endl;
}

vector<serverListStruct> AliStorageEventManager::GetServerListVector(storageSockets socket, int timeout)
{
    //check timeout
    fZmqError=false;
    zmqPoll(fSockets[socket],timeout);
    if(fZmqError){vector<serverListStruct> emptyVector;return emptyVector;}
    
    //get size of the incomming message
    zmq_msg_t buffer;
    zmqInit(&buffer);
    zmqRecv(&buffer,fSockets[socket],0);
    int numberOfRecords;
    memcpy(&numberOfRecords,&buffer,sizeof(int));
    if(numberOfRecords==0){cout<<"MANAGER -- list is empty"<<endl;}
    
    //send empty message just to keep req-rep order:
    zmqSend(&buffer,fSockets[socket],0);
    
    //get list of events
    zmqRecv(&buffer,fSockets[socket],0);
    void *tmp = zmq_msg_data(&buffer);
    vector<serverListStruct> receivedList(static_cast<serverListStruct*>(tmp),
                                          static_cast<serverListStruct*>(tmp)+numberOfRecords);
    zmq_msg_close(&buffer);
    return receivedList;
}

AliESDEvent* AliStorageEventManager::GetEvent(storageSockets socket,int timeout)
{
    //check timeout
    fZmqError=false;
    zmqPoll(fSockets[socket],timeout);
    if(fZmqError){return NULL;}
    
    //reveive buffer
    zmq_msg_t buffer;
    zmqInit(&buffer);
    zmqRecv(&buffer,fSockets[socket],0);
    if(fZmqError){return NULL;}
    
    //read buffer to TMessage
    TBufferFile *mess = new TBufferFile(TBuffer::kRead,
                                        zmq_msg_size(&buffer)+sizeof(UInt_t),
                                        zmq_msg_data(&buffer));
    mess->InitMap();
    mess->ReadClass();// get first the class stored in message
    mess->SetBufferOffset(sizeof(UInt_t) + sizeof(kMESS_OBJECT));
    mess->ResetMap();
    
    //read ESDEvent from TMessage
    AliESDEvent* data = (AliESDEvent*)(mess->ReadObjectAny(AliESDEvent::Class()));
    if (data)
    {
        data->GetStdContent();
        zmq_msg_close(&buffer);
        return data;
    }
    else
    {
        zmq_msg_close(&buffer);
        return NULL;
    }
}

struct serverRequestStruct* AliStorageEventManager::GetServerStruct(storageSockets socket)
{
    fZmqError=false;
    struct serverRequestStruct *request;
    
    zmq_msg_t buffer;
    zmqInit(&buffer);
    zmqRecv(&buffer,fSockets[socket],0);
    
    if(fZmqError){ //if there were errors when receiving message
        request = new struct serverRequestStruct;
        request->messageType = -1;
    }
    else{request = new struct serverRequestStruct(*(static_cast<struct serverRequestStruct*>(zmq_msg_data(&buffer))));}
    zmq_msg_close(&buffer);
    return request;
}

struct clientRequestStruct* AliStorageEventManager::GetClientStruct(storageSockets socket,int timeout)
{
    //  pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}} ;
    //  if(timeout>=0){if(poll (&items[0], 1, timeout)==0){return NULL;}}
    
    zmq_msg_t buffer;
    zmqInit(&buffer);
    zmqRecv(&buffer,fSockets[socket],0);
    
    struct clientRequestStruct *request = new struct clientRequestStruct(*(static_cast<struct clientRequestStruct*>(zmq_msg_data(&buffer))));
    
    zmq_msg_close(&buffer);
    return request;
}

bool AliStorageEventManager::GetBool(storageSockets socket)
{
    zmq_msg_t buffer;
    zmqInit(&buffer);
    zmqRecv(&buffer,fSockets[socket],0);
    bool result;
    memcpy(&result,&buffer,sizeof(bool));
    zmq_msg_close(&buffer);
    return result;
}

long AliStorageEventManager::GetLong(storageSockets socket)
{
    zmq_msg_t buffer;
    zmqInit(&buffer);
    zmqRecv(&buffer,fSockets[socket],0);
    long result;
    memcpy(&result,&buffer,sizeof(long));
    zmq_msg_close(&buffer);
    return result;
}

// ZMQ methods wrappers:
void AliStorageEventManager::zmqInit(zmq_msg_t *msg,size_t size)
{
    if(size==0){
        if(zmq_msg_init(msg) != 0){
            if(zmq_errno() != EAGAIN){
                cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
                fZmqError=true;
            }
        }
    }
    else{
        if(zmq_msg_init_size(msg,size) != 0){
            if(zmq_errno() != EAGAIN){
                cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
                fZmqError=true;
            }
        }
    }
}

void AliStorageEventManager::zmqSend(zmq_msg_t *msg,void *socket,int flags)
{
    if(zmq_msg_send(msg,socket,flags) != 0){
        if(zmq_errno() != EAGAIN)
        {
            cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
            fZmqError=true;
        }
    }
}

void AliStorageEventManager::zmqRecv(zmq_msg_t *msg,void *socket,int flags)
{
    if(zmq_msg_recv(msg,socket,flags) != 0){
        if(zmq_errno() != EAGAIN)
        {
            cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
            fZmqError=true;
        }
    }
}

void AliStorageEventManager::zmqPoll(void *socket,int timeout)
{
    if(timeout>=0)
    {
        zmq_pollitem_t items[1] =  {{socket,0,ZMQ_POLLIN,0}};
        if(zmq_poll(&items[0],1, timeout)!=0)
        {
            if(zmq_errno() != EAGAIN)
            {
                cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
                fZmqError=true;
            }
        }
    }
}




