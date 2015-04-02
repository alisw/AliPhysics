#include "AliZMQManager.h"

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

AliZMQManager *AliZMQManager::fManagerInstance = 0;

AliZMQManager::AliZMQManager()
{
    //read config file
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
    
    for(int i=0;i<NUMBER_OF_SOCKETS;i++)
    {
        fContexts[i] = zmq_ctx_new();
//        CreateSocket((storageSockets)i);
    }
}
AliZMQManager::~AliZMQManager()
{
    // close sockets and destroy contexts
    for (int i=0; i<NUMBER_OF_SOCKETS; i++)
    {
        zmq_close(fSockets[i]);
        zmq_ctx_destroy(fContexts[i]);
    }
}

AliZMQManager* AliZMQManager::GetInstance()
{
    TThread::Lock();
    if(fManagerInstance==0)
    {
        cout<<"\n\nMANAGER -- creating new instance of ZMQ Manager\n\n"<<endl;
        fManagerInstance = new AliZMQManager();
    }
    TThread::UnLock();
    return fManagerInstance;
}

void AliZMQManager::CreateSocket(storageSockets socket)
{
    int timeout = 5000;
    int linger = -1;
    
    switch (socket)
    {
        case SERVER_COMMUNICATION_REQ:
        {
            // server communication - REQ
            fSockets[SERVER_COMMUNICATION_REQ] = zmq_socket(fContexts[SERVER_COMMUNICATION_REQ],ZMQ_REQ);
            
            if(0 != zmq_setsockopt(fSockets[SERVER_COMMUNICATION_REQ],ZMQ_RCVTIMEO,&timeout,sizeof(int)))
            {cout<<"MANAGER -- create socket 1 -- "<<zmq_strerror(zmq_errno())<<endl;}
            if(0 != zmq_setsockopt(fSockets[SERVER_COMMUNICATION_REQ],ZMQ_LINGER,&linger,sizeof(int)))
            {cout<<"MANAGER -- create socket 1a -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            if(0 != zmq_connect(fSockets[SERVER_COMMUNICATION_REQ],Form("tcp://%s:%d",fStorageServer.c_str(),fStorageServerPort)))
            {cout<<"MANAGER -- create socket 1b -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        case SERVER_COMMUNICATION_REP:
        {
            // server communication - REP
            fSockets[SERVER_COMMUNICATION_REP] = zmq_socket(fContexts[SERVER_COMMUNICATION_REP],ZMQ_REP);
            
            if(0 != zmq_setsockopt(fSockets[SERVER_COMMUNICATION_REP],ZMQ_RCVTIMEO,&timeout,sizeof(int)))
            {cout<<"MANAGER -- create socket 2 -- "<<zmq_strerror(zmq_errno())<<endl;}
            if(0 != zmq_setsockopt(fSockets[SERVER_COMMUNICATION_REP],ZMQ_LINGER,&linger,sizeof(int)))
            {cout<<"MANAGER -- create socket 2a -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            if(0 != zmq_bind(fSockets[SERVER_COMMUNICATION_REP],Form("tcp://*:%d",fStorageServerPort)))
            {cout<<"MANAGER -- create socket 2b -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        case CLIENT_COMMUNICATION_REQ:
        {
            // client communication - REQ
            fSockets[CLIENT_COMMUNICATION_REQ] = zmq_socket(fContexts[CLIENT_COMMUNICATION_REQ],ZMQ_REQ);

            if(0 != zmq_setsockopt(fSockets[CLIENT_COMMUNICATION_REQ],ZMQ_RCVTIMEO,&timeout,sizeof(int)))
            {cout<<"MANAGER -- create socket 3 -- "<<zmq_strerror(zmq_errno())<<endl;}
            if(0 != zmq_setsockopt(fSockets[CLIENT_COMMUNICATION_REQ],ZMQ_LINGER,&linger,sizeof(int)))
            {cout<<"MANAGER -- create socket 3a -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            if(0 != zmq_connect(fSockets[CLIENT_COMMUNICATION_REQ],Form("tcp://%s:%d",fStorageServer.c_str(), fStorageClientPort)))
            {cout<<"MANAGER -- create socket 3b -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        case CLIENT_COMMUNICATION_REP:
        {
            // client communication - REP
            fSockets[CLIENT_COMMUNICATION_REP] = zmq_socket(fContexts[CLIENT_COMMUNICATION_REP],ZMQ_REP);
            
            if(0 != zmq_setsockopt(fSockets[CLIENT_COMMUNICATION_REP],ZMQ_RCVTIMEO,&timeout,sizeof(int)))
            {cout<<"MANAGER -- create socket 4 -- "<<zmq_strerror(zmq_errno())<<endl;}
            if(0 != zmq_setsockopt(fSockets[CLIENT_COMMUNICATION_REP],ZMQ_LINGER,&linger,sizeof(int)))
            {cout<<"MANAGER -- create socket 4a -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            if(0 != zmq_bind(fSockets[CLIENT_COMMUNICATION_REP],Form("tcp://*:%d",fStorageClientPort)))
            {cout<<"MANAGER -- create socket 4b -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        case EVENTS_SERVER_PUB:
        {
            // events publisher
            fSockets[EVENTS_SERVER_PUB] = zmq_socket(fContexts[EVENTS_SERVER_PUB],ZMQ_PUB);
            if(0 != zmq_bind(fSockets[EVENTS_SERVER_PUB],Form("tcp://*:%d",fEventServerPort)))
            {cout<<"MANAGER -- create socket 5 -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            break;
        }
        case EVENTS_SERVER_SUB:
        {
            // events subscriber
            fSockets[EVENTS_SERVER_SUB] = zmq_socket(fContexts[EVENTS_SERVER_SUB],ZMQ_SUB);
            if(0 != zmq_setsockopt(fSockets[EVENTS_SERVER_SUB],ZMQ_SUBSCRIBE,"",0))
            {cout<<"MANAGER -- create socket 6b -- "<<zmq_strerror(zmq_errno())<<endl;}
            if(0 != zmq_connect(fSockets[EVENTS_SERVER_SUB],Form("tcp://%s:%d",fEventServer.c_str(),fEventServerPort)))
            {cout<<"MANAGER -- create socket 6c -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            break;
        }
        case XML_PUB:
        {
            // xml publisher
            fSockets[XML_PUB] = zmq_socket(fContexts[XML_PUB],ZMQ_PUB);
            if(0 != zmq_bind(fSockets[XML_PUB],Form("tcp://*:%d",fXmlServerPort)))
            {cout<<"MANAGER -- create socket 7 -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        case ITS_POINTS_PUB:
        {
            // ITS recpoints publisher
            fSockets[ITS_POINTS_PUB] = zmq_socket(fContexts[ITS_POINTS_PUB],ZMQ_PUB);
            if(0 != zmq_bind(fSockets[ITS_POINTS_PUB],Form("tcp://*:%d",fItsPointsServerPort)))
            {cout<<"MANAGER -- create socket 8 -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        case ITS_POINTS_SUB:
        {// ITS recpoints subscriber
            fSockets[ITS_POINTS_SUB] = zmq_socket(fContexts[ITS_POINTS_SUB],ZMQ_SUB);
            if(zmq_setsockopt(fSockets[ITS_POINTS_SUB],ZMQ_SUBSCRIBE,"",0)!=0)
            {cout<<"MANAGER -- create socket 9 -- "<<zmq_strerror(zmq_errno())<<endl;}
            if(0 != zmq_connect(fSockets[ITS_POINTS_SUB],Form("tcp://%s:%d",fEventServer.c_str(),fItsPointsServerPort)))
            {cout<<"MANAGER -- create socket 9a -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        default:break;
    }
}

bool AliZMQManager::Send(vector<serverListStruct> list,storageSockets socket)
{
    //send size of the struct first
    int numberOfRecords = list.size();
    cout<<"MANAGER -- sending vector with "<<numberOfRecords<<" records"<<endl;
    zmq_msg_t buffer;
    
    if(!zmqInit(&buffer,sizeof(int))){return false;}
    memcpy(zmq_msg_data(&buffer),&numberOfRecords,sizeof(int));
    
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        cout<<"MANAGER -- couldn't send list's size"<<endl;
        zmq_msg_close(&buffer);
        return false;
    }
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        cout<<"MANAGER -- couldn't receive message inside Send vector call. This may cause serious problems!"<<endl;
        zmq_msg_close(&buffer);
        return false;
    }
    if(numberOfRecords==0)
    {
        cout<<"MANAGER -- list size = 0"<<endl;
//        return false;
    }
    zmq_msg_close(&buffer);
    
    zmq_msg_t message;
    if(!zmqInit(&message,sizeof(serverListStruct)*numberOfRecords)){return false;}
    memcpy(zmq_msg_data(&message),reinterpret_cast<void*> (&list[0]), sizeof(serverListStruct)*numberOfRecords);
    
    if(!zmqSend(&message,fSockets[socket],0))
    {
        zmq_msg_close(&message);
        return false;
    }
    zmq_msg_close(&message);
    return true;
}

bool AliZMQManager::Send(struct serverRequestStruct *request,storageSockets socket)
{
    int sizeOfRequest = sizeof(struct serverRequestStruct)+sizeof(struct listRequestStruct)+sizeof(struct eventStruct);
    
    cout<<"MANAGER -- sending serverRequestStruct:"<<request->messageType<<"\t"<<request->list.runNumber[0]<<endl;
    
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeOfRequest)){return false;}
    memcpy(zmq_msg_data(&buffer),request,sizeOfRequest);
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    zmq_msg_close(&buffer);
    return true;
}

bool AliZMQManager::Send(struct clientRequestStruct *request,storageSockets socket)
{
    cout<<"MANAGER -- sending clientRequestStruct:"<<request->messageType<<"\t"<<endl;
    
    //put clientRequestStruct in buffer
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeof(struct clientRequestStruct))){return false;}
    memcpy(zmq_msg_data(&buffer),request,sizeof(struct clientRequestStruct));
    
    //send buffer
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    zmq_msg_close(&buffer);
    return true;
}
bool AliZMQManager::Send(long message,storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeof(long))){return false;}
    memcpy(zmq_msg_data(&buffer),&message,sizeof(long));
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    
    zmq_msg_close(&buffer);
    return true;
}
bool AliZMQManager::Send(bool message,storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeof(bool))){return false;}
    memcpy(zmq_msg_data(&buffer),&message,sizeof(bool));
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    zmq_msg_close(&buffer);
    return true;
}
bool AliZMQManager::Send(AliESDEvent *event, storageSockets socket)
{
    if(!event){cout<<"MANAGER -- no event"<<endl;return false;}
    
    TMessage tmess(kMESS_OBJECT);
    tmess.Reset();
    tmess.WriteObject(event);
    TMessage::EnableSchemaEvolutionForAll(kTRUE);
    
    int bufsize = tmess.BufferSize();
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,bufsize)){return false;}
    memcpy(zmq_msg_data(&buffer),tmess.Buffer(),bufsize);
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    zmq_msg_close(&buffer);
    return true;
}

bool AliZMQManager::SendAsXml(AliESDEvent *event,storageSockets socket)
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
    if(!zmqInit(&buffer,bufferString.size())){return false;}
    memcpy (zmq_msg_data(&buffer), bufferString.data(), bufferString.size());
    
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    zmq_msg_close(&buffer);
    return true;
}

bool AliZMQManager::Get(std::vector<serverListStruct>* &result,storageSockets socket)
{
    //get size of the incomming message
    zmq_msg_t buffer;
    zmqInit(&buffer);
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    int numberOfRecords;
    memcpy(&numberOfRecords,&buffer,sizeof(int));
    cout<<"MANAGER -- number of records:"<<numberOfRecords<<endl;
    //send empty message just to keep req-rep order:
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        cout<<"MANAGER -- couldn't send message inside GetServerListVector. This may cause seriouse problems!"<<endl;
        zmq_msg_close(&buffer);
        return false;
    }

    //get list of events
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }

    if(numberOfRecords==0){
        cout<<"MANAGER -- list is empty"<<endl;
        result = nullptr;
        zmq_msg_close(&buffer);
        return false;
    }
    else
    {
        // read data from buffer:
        void *tmp = zmq_msg_data(&buffer);
        
        // vector's range constructor rebuilding vector from void*
        vector<serverListStruct> tmpVector(static_cast<serverListStruct*>(tmp),
                                            static_cast<serverListStruct*>(tmp)+numberOfRecords);
        
        cout<<"MANAGER -- size of vector:"<<tmpVector.size()<<endl;
        
        // create pointer to this vector:
        result = new vector<serverListStruct>(tmpVector);
        
        cout<<"MANAGER -- size of vector (from pointer:)"<<result->size()<<endl;
        
        zmq_msg_close(&buffer);
        return true;
    }
}

bool AliZMQManager::Get(AliESDEvent* &result, storageSockets socket)
{
    //reveive buffer
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return false;}
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    
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
        cout<<"MANAGER -- received valid event"<<endl;
        data->GetStdContent();
        cout<<"MANAGER -- reading std content:"<<data->GetEventNumberInFile()<<endl;
        zmq_msg_close(&buffer);
        result = data;
        return true;
    }
    else
    {
        zmq_msg_close(&buffer);
        return false;
    }
}

bool AliZMQManager::Get(struct serverRequestStruct* &result, storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return false;}
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    else
    {
        result = new struct serverRequestStruct(*(static_cast<struct serverRequestStruct*>(zmq_msg_data(&buffer))));
        
        cout<<"MANAGER -- received server request:"<<result->messageType<<"\t"<<result->list.runNumber[0]<<endl;
        zmq_msg_close(&buffer);
        return true;
    }
}

bool AliZMQManager::Get(struct clientRequestStruct* &result, storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return false;}
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    else
    {
        result = new struct clientRequestStruct(*(static_cast<struct clientRequestStruct*>(zmq_msg_data(&buffer))));
        cout<<"MANAGER -- received client request:"<<result->messageType<<endl;
        zmq_msg_close(&buffer);
        return true;
    }
}

bool AliZMQManager::Get(long *result, storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return false;}
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    else
    {
        memcpy(result,&buffer,sizeof(bool));
        zmq_msg_close(&buffer);
        return true;
    }
}

bool AliZMQManager::Get(bool *result, storageSockets socket)
{
    zmq_msg_t buffer;

    if(!zmqInit(&buffer)){return false;}
    if(!zmqRecv(&buffer,fSockets[socket],0))
    {
        zmq_msg_close(&buffer);
        return false;
    }
    else
    {
        memcpy(result,&buffer,sizeof(long));
        zmq_msg_close(&buffer);
        return true;
    }
}

void AliZMQManager::RecreateSocket(storageSockets socket)
{
    cout<<"MANAGER -- recreating socket:"<<socket<<endl;
    cout<<"zmq_close:"<<zmq_close(fSockets[socket])<<endl;
//    cout<<"zmq_term:"<<zmq_ctx_destroy(fContexts[socket])<<endl;
//    fContexts[socket] = zmq_ctx_new();
    CreateSocket(socket);
}

// ZMQ methods wrappers:
bool AliZMQManager::zmqInit(zmq_msg_t *msg,size_t size)
{
    if(size==0){
        if(zmq_msg_init(msg) != 0){
            if(zmq_errno() != EAGAIN){
                cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
                return false;
            }
        }
    }
    else{
        if(zmq_msg_init_size(msg,size) != 0){
            if(zmq_errno() != EAGAIN){
                cout<<"MANAGER -- "<<zmq_strerror(zmq_errno())<<endl;
                return false;
            }
        }
    }
    return true;
}

bool AliZMQManager::zmqSend(zmq_msg_t *msg,void *socket,int flags)
{
    if(zmq_msg_send(msg,socket,flags) == -1){
        if(zmq_errno() != EAGAIN) // ignore timeout problems
        {
            if(zmq_errno() == EFSM)// cannot be accomplished in current state
            {}
            cout<<"MANAGER -- zmqSend -- "<<zmq_strerror(zmq_errno())<<endl;
            return false;
        }
    }
    return true;
}

bool AliZMQManager::zmqRecv(zmq_msg_t *msg,void *socket,int flags)
{
    if(zmq_msg_recv(msg,socket,flags) == -1){
        cout<<"MANAGER -- zmqRecv -- "<<zmq_strerror(zmq_errno())<<endl;
        return false;
    }
    return true;
}



