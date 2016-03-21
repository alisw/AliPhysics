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

#include "zmq.h"

using namespace std;

AliZMQManager *AliZMQManager::fManagerInstance = 0;

/// Defaul counstructor of AliZMQManager
///
/// In constructor, config files with hostnames and port number are being read.
/// Contexts for all sockets are being created.
AliZMQManager::AliZMQManager()
	: fStorageServer()
	, fEventServer()
	, fStorageServerPort(0)
	, fStorageClientPort(0)
	, fEventServerPort(0)
	, fXmlServerPort(0)
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
        }
        if(configFile.eof()){configFile.clear();}
        configFile.close();
    }
    else{cout<<"EVENT MANAGER -- Unable to open config file"<<endl;}
    fContext = zmq_ctx_new();
}

/// Default destructor of AliZMQManager.
///
/// All contexts and sockets will be closed/destroyed.
AliZMQManager::~AliZMQManager()
{
}

void AliZMQManager::Close()
{
    // destroy context
    cout<<"\n\nDestroying context\n\n"<<endl;
    zmq_ctx_destroy(fContext);
}

/// Method to get instance of AliZMQManager.
///
/// Only one instance of AliZMQManager will be created (thread-safe).
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

/// Method to create socket with given name.
///
/// Names of sockets are defined in AliStorageTypes.h
/// This method should be extended in case of need for additional sockets.
///
/// \param socket Name of socket (of type storageSockets) to be created.
void AliZMQManager::CreateSocket(storageSockets socket)
{
    int timeout = 5000;
    int linger = -1;
    
    switch (socket)
    {
        case SERVER_COMMUNICATION_REQ:
        {
            // server communication - REQ
            fSockets[SERVER_COMMUNICATION_REQ] = zmq_socket(fContext,ZMQ_REQ);
            
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
            fSockets[SERVER_COMMUNICATION_REP] = zmq_socket(fContext,ZMQ_REP);
            
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
            fSockets[CLIENT_COMMUNICATION_REQ] = zmq_socket(fContext,ZMQ_REQ);

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
            fSockets[CLIENT_COMMUNICATION_REP] = zmq_socket(fContext,ZMQ_REP);
            
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
            fSockets[EVENTS_SERVER_PUB] = zmq_socket(fContext,ZMQ_PUB);
            if(0 != zmq_bind(fSockets[EVENTS_SERVER_PUB],Form("tcp://*:%d",fEventServerPort)))
            {cout<<"MANAGER -- create socket 5 -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            break;
        }
        case EVENTS_SERVER_SUB:
        {
            // events subscriber
            fSockets[EVENTS_SERVER_SUB] = zmq_socket(fContext,ZMQ_SUB);
            if(0 != zmq_setsockopt(fSockets[EVENTS_SERVER_SUB],ZMQ_SUBSCRIBE,"",0))
            {cout<<"MANAGER -- create socket 6b -- "<<zmq_strerror(zmq_errno())<<endl;}
            if(0 != zmq_connect(fSockets[EVENTS_SERVER_SUB],Form("tcp://%s:%d",fEventServer.c_str(),fEventServerPort)))
            {cout<<"MANAGER -- create socket 6c -- "<<zmq_strerror(zmq_errno())<<endl;}
            
            break;
        }
        case XML_PUB:
        {
            // xml publisher
            fSockets[XML_PUB] = zmq_socket(fContext,ZMQ_PUB);
            if(0 != zmq_bind(fSockets[XML_PUB],Form("tcp://*:%d",fXmlServerPort)))
            {cout<<"MANAGER -- create socket 7 -- "<<zmq_strerror(zmq_errno())<<endl;}
            break;
        }
        default:break;
    }
}

/// Method sends vector of structs of type serverListStruct to given socket.
///
/// \param list Vector to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Send(vector<serverListStruct> list,storageSockets socket)
{
    //send size of the struct first
    int numberOfRecords = list.size();
    cout<<"MANAGER -- sending vector with "<<numberOfRecords<<" records"<<endl;
    zmq_msg_t buffer;
    
    if(!zmqInit(&buffer,sizeof(int))){return false;}
    memcpy(zmq_msg_data(&buffer),&numberOfRecords,sizeof(int));
    
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        cout<<"MANAGER -- couldn't send list's size"<<endl;
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus != 1)
    {
        cout<<"MANAGER -- couldn't receive message inside Send vector call. This may cause serious problems!"<<endl;
        zmq_msg_close(&buffer);
        return recvStatus;
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
    sendStatus = zmqSend(&message,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&message);
        return sendStatus;
    }
    zmq_msg_close(&message);
    return 1;
}

/// Method sends vector of structs of type string100 to given socket.
/// string100 is a fixed length string defined in AliStorageTypes.h
/// Use member called "data" to access string stored in it.
///
/// \param list Vector to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Send(vector<string100> list,storageSockets socket)
{
    cout<<"MANAGER -- sending vector of string100:"<<endl;
    vector<string100>::iterator it;
    for (it = list.begin(); it != list.end(); ++it)
    {
        cout<<(*it).data<<endl;
    }
    
    //send size of the set first
    int numberOfRecords = list.size();
    cout<<"MANAGER -- sending set with "<<numberOfRecords<<" records"<<endl;
    zmq_msg_t buffer;
    
    if(!zmqInit(&buffer,sizeof(int))){return 0;}
    memcpy(zmq_msg_data(&buffer),&numberOfRecords,sizeof(int));
    
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        cout<<"MANAGER -- couldn't send list's size"<<endl;
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus != 1)
    {
        cout<<"MANAGER -- couldn't receive message inside Send set call. This may cause serious problems!"<<endl;
        zmq_msg_close(&buffer);
        return recvStatus;
    }
    if(numberOfRecords==0)
    {
        cout<<"MANAGER -- list size = 0"<<endl;
        //        return false;
    }
    zmq_msg_close(&buffer);
    
    zmq_msg_t message;
    if(!zmqInit(&message,sizeof(string100)*numberOfRecords)){return false;}

    memcpy(zmq_msg_data(&message),reinterpret_cast<void*> (&list[0]), sizeof(string100)*numberOfRecords);

    sendStatus = zmqSend(&message,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&message);
        return sendStatus;
    }
    zmq_msg_close(&message);
    return 1;
}

/// Method sends structs of type serverRequestStruct to given socket.
///
/// Define your request as:
///
///     struct serverRequestStruct *request = new struct serverRequestStruct;
///
/// and call this method as:
///
///     manager->Send(request,socket);
///
/// \param request Struct to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Send(struct serverRequestStruct *request,storageSockets socket)
{
    size_t sizeOfRequest = sizeof(struct serverRequestStruct);/*+sizeof(struct listRequestStruct)+sizeof(struct eventStruct);*/
    
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeOfRequest)){return false;}
    memcpy(zmq_msg_data(&buffer),request,sizeOfRequest);
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    zmq_msg_close(&buffer);
    return 1;
}

/// Method sends structs of type clientRequestStruct to given socket.
///
/// Define your request as:
///
///     struct clientRequestStruct *request = new struct clientRequestStruct;
///
/// and call this method as:
///
///     manager->Send(request,socket);
///
/// \param request Struct to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Send(struct clientRequestStruct *request,storageSockets socket)
{
    cout<<"MANAGER -- sending clientRequestStruct:"<<request->messageType<<"\t"<<endl;
    
    //put clientRequestStruct in buffer
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeof(struct clientRequestStruct))){return false;}
    memcpy(zmq_msg_data(&buffer),request,sizeof(struct clientRequestStruct));
    
    //send buffer
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    zmq_msg_close(&buffer);
    return 1;
}

/// Method sends message of type long to given socket.
///
/// \param message Long message to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Send(long message,storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeof(long))){return false;}
    memcpy(zmq_msg_data(&buffer),&message,sizeof(long));
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    
    zmq_msg_close(&buffer);
    return 1;
}

/// Method sends message of type bool to given socket.
///
/// \param message Bool message to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Send(bool message,storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer,sizeof(bool))){return false;}
    memcpy(zmq_msg_data(&buffer),&message,sizeof(bool));
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    zmq_msg_close(&buffer);
    return 1;
}

/// Method sends message of type AliESDEvent to given socket.
///
/// \param event Pointer to AliESDEvent object which is to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Send(AliESDEvent *event, storageSockets socket)
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
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    zmq_msg_close(&buffer);
    return 1;
}

/// Method sends AliESDEvent as an xml to given socket.
///
/// This method extracts some information from AliESDEvent, puts it in a string
/// formatted as an xml and sends it to a given socket.
///
///
/// \param event Pointer to AliESDEvent object which is to be sent
/// \param socket Name of socket to which list shold be sent
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::SendAsXml(AliESDEvent *event,storageSockets socket)
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
    
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        zmq_msg_close(&buffer);
        return sendStatus;
    }
    zmq_msg_close(&buffer);
    return 1;
}

/// Method to get vector of serverListStructs from given socket.
///
/// On of the possible ways to use this method is:
///
/// vector<serverListStruct> *tmpVector;
/// manager->Get(tmpVector,socket);
/// vector<serverListStruct> &receivedList = *tmpVector;
///
/// \param result Adres of pointer to which resulting vector should be written
/// \param socket Name of socket from which message should be received
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Get(vector<serverListStruct>* &result,storageSockets socket)
{
    //get size of the incomming message
    zmq_msg_t buffer;
    zmqInit(&buffer);
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus != 1)
    {
        cout<<"MANAGER -- couldn't receive number of records inside GetServerListVector."<<endl;
        zmq_msg_close(&buffer);
        return recvStatus;
    }
    int numberOfRecords;
    memcpy(&numberOfRecords,&buffer,sizeof(int));
    cout<<"MANAGER -- number of records:"<<numberOfRecords<<endl;
    //send empty message just to keep req-rep order:
    int sendStatus = zmqSend(&buffer,fSockets[socket],0);
    if(sendStatus != 1)
    {
        cout<<"MANAGER -- couldn't send message inside GetServerListVector. This may cause seriouse problems!"<<endl;
        zmq_msg_close(&buffer);
        return sendStatus;
    }

    //get list of events
    recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus != 1)
    {
        zmq_msg_close(&buffer);
        return recvStatus;
    }

    if(numberOfRecords==0){
        cout<<"MANAGER -- list is empty"<<endl;
        result = nullptr;
        zmq_msg_close(&buffer);
        return 0;
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
        
        //delete tmpVector; ?
        zmq_msg_close(&buffer);
        return 1;
    }
}

/// Method to get vector of string100 structs from given socket.
///
/// On of the possible ways to use this method is:
///
/// vector<string100> *tmpVector;
/// manager->Get(tmpVector,socket);
/// vector<string100> &receivedList = *tmpVector;
///
/// \param result Adres of pointer to which resulting vector should be written
/// \param socket Name of socket from which message should be received
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Get(vector<string100>* &result,storageSockets socket)
{
    //get size of the incomming message
    zmq_msg_t buffer;
    zmqInit(&buffer);
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus != 1)
    {
        cout<<"MANAGER -- couldn't receive number of records inside Get set<string>."<<endl;
        zmq_msg_close(&buffer);
        return recvStatus;
    }
    int numberOfRecords;
    memcpy(&numberOfRecords,&buffer,sizeof(int));
    cout<<"MANAGER -- number of records:"<<numberOfRecords<<endl;
    //send empty message just to keep req-rep order:
    if(!zmqSend(&buffer,fSockets[socket],0))
    {
        cout<<"MANAGER -- couldn't send message inside Get set<string>. This may cause seriouse problems!"<<endl;
        zmq_msg_close(&buffer);
        return false;
    }
    
    //get list of events
    recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus != 1)
    {
        zmq_msg_close(&buffer);
        return recvStatus;
    }
    
    if(numberOfRecords==0){
        cout<<"MANAGER -- list is empty"<<endl;
        result = nullptr;
        zmq_msg_close(&buffer);
        return 0;
    }
    else
    {
        // read data from buffer:
        void* tmp = zmq_msg_data(&buffer);
        
        // vector's range constructor rebuilding vector from void*
        vector<string100> newVector(static_cast<string100*>(tmp),
                                    static_cast<string100*>(tmp)+numberOfRecords);
        
        cout<<"MANAGER -- size of vector:"<<newVector.size()<<endl;
        
        cout<<"MANAGER -- received vector:"<<endl;
        vector<string100>::iterator it;
        for (it = newVector.begin(); it != newVector.end(); ++it)
        {
            cout<<(*it).data<<endl;
        }
        
        // create pointer to this vector:
        result = new vector<string100>(newVector);
        cout<<"MANAGER -- size of vector (from pointer:)"<<result->size()<<endl;
        
        //delete tmpSet; ?
        zmq_msg_close(&buffer);
        return 1;
    }
}

/// Method to get AliESDEvent from given socket.
///
/// Usage:
///
/// AliESDEvent *event;
/// manager->Get(event,socket);
///
/// \param result Adres of pointer to which resulting event should be written
/// \param socket Name of socket from which message should be received
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Get(AliESDEvent* &result, storageSockets socket)
{
    //reveive buffer
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return 0;}
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus != 1)
    {
        zmq_msg_close(&buffer);
        return recvStatus;
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
        return 1;
    }
    else
    {
        zmq_msg_close(&buffer);
        result=0;
        return 0;
    }
}

/// Method to get serverRequestStruct from given socket.
///
/// Usage:
///
/// struct serverRequestStruct *request;
/// manager->Get(request,socket);
///
/// \param result Adres of pointer to which resulting struct should be written
/// \param socket Name of socket from which message should be received
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Get(struct serverRequestStruct* &result, storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return 0;}
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus == 1)
    {
        result = new struct serverRequestStruct(*(static_cast<struct serverRequestStruct*>(zmq_msg_data(&buffer))));
        
        cout<<"MANAGER -- received server request:"<<result->messageType<<"\t"<<result->runNumber[0]<<endl;
    }
    zmq_msg_close(&buffer);
    return recvStatus;
}

/// Method to get clientRequestStruct from given socket.
///
/// Usage:
///
/// struct clientRequestStruct *request;
/// manager->Get(request,socket);
///
/// \param result Adres of pointer to which resulting struct should be written
/// \param socket Name of socket from which message should be received
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Get(struct clientRequestStruct* &result, storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return 0;}
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus == 1)
    {
        result = new struct clientRequestStruct(*(static_cast<struct clientRequestStruct*>(zmq_msg_data(&buffer))));
        cout<<"MANAGER -- received client request:"<<result->messageType<<endl;
    }
    zmq_msg_close(&buffer);
    return recvStatus;
}

/// Method to get message of type long from given socket.
///
/// \param result Pointer to which resulting message should be written
/// \param socket Name of socket from which message should be received
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Get(long *result, storageSockets socket)
{
    zmq_msg_t buffer;
    if(!zmqInit(&buffer)){return 0;}
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus == 1)
    {
        memcpy(result,&buffer,sizeof(bool));

    }
    zmq_msg_close(&buffer);
    return recvStatus;
}

/// Method to get message of type bool from given socket.
///
/// \param result Pointer to which resulting message should be written
/// \param socket Name of socket from which message should be received
///
/// \return Returns true in case of success, false in case of failure
int AliZMQManager::Get(bool *result, storageSockets socket)
{
    zmq_msg_t buffer;

    if(!zmqInit(&buffer)){return 0;}
    int recvStatus = zmqRecv(&buffer,fSockets[socket],0);
    if(recvStatus == 1)
    {
        memcpy(result,&buffer,sizeof(long));
    }
    zmq_msg_close(&buffer);
    return recvStatus;
}

/// Method to close and create from scratch given socket.
///
/// This must be used if socket is expecting to receive a message, but you need to send one
/// and vice versa. Use only when you are sure what you are doing, you should design communication
/// pattern yourself depending on what you need.
///
/// \param socket Name of socket to be recreated
void AliZMQManager::RecreateSocket(storageSockets socket)
{
    cout<<"MANAGER -- recreating socket:"<<socket<<endl;
    zmq_close(fSockets[socket]);
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

int AliZMQManager::zmqSend(zmq_msg_t *msg,void *socket,int flags)
{
    if(zmq_msg_send(msg,socket,flags) == -1){
        if(zmq_errno() != EAGAIN) // ignore timeout problems
        {
            cout<<"MANAGER -- zmqSend -- "<<zmq_strerror(zmq_errno())<<endl;
            zmq_close(socket);
            return -1;
        }
        else
        {
            cout<<"MANAGER -- zmqSend -- timeout"<<endl;
            return 0;
        }
    }
    return 1;
}

int AliZMQManager::zmqRecv(zmq_msg_t *msg,void *socket,int flags)
{
    if(zmq_msg_recv(msg,socket,flags) == -1)
    {
        if((zmq_errno() != EAGAIN) && (zmq_errno() != EINTR)) // ignore timeout problems and interrupted system calls
        {
            cout<<"MANAGER -- zmqRecv -- "<<zmq_strerror(zmq_errno())<<endl;
            zmq_close(socket);
            return -1;
        }
        else
        {
            cout<<"MANAGER -- zmqRecv -- timeout"<<endl;
            return 0;
        }
    }
    return 1;
}


