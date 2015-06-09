#include "AliStorageEventManager.h"

#include <iostream>
#include <sstream>
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

#include "zmq.hpp"

using namespace zmq;
using namespace std;

AliStorageEventManager *AliStorageEventManager::fManagerInstance = 0;

AliStorageEventManager::AliStorageEventManager()
{
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
    
    for(int i=0;i<NUMBER_OF_SOCKETS;i++){fContexts[i] = new context_t();}
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


void freeBuff (void *data, void *hint)
{
    //  free(data);
}

bool AliStorageEventManager::CreateSocket(storageSockets socket)
{
    cout<<"Creating socket:"<<socket<<endl;
    
    switch(socket)
    {
        case SERVER_COMMUNICATION_REQ:
        {
            fSockets[SERVER_COMMUNICATION_REQ] =
            new socket_t(*fContexts[SERVER_COMMUNICATION_REQ],ZMQ_REQ);
            try
            {
                fSockets[SERVER_COMMUNICATION_REQ]->connect(Form("tcp://%s:%d",fStorageServer.c_str(),fStorageServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case SERVER_COMMUNICATION_REP:
        {
            fSockets[SERVER_COMMUNICATION_REP] =
            new socket_t(*fContexts[SERVER_COMMUNICATION_REP],ZMQ_REP);
            try
            {
                fSockets[SERVER_COMMUNICATION_REP]->bind(Form("tcp://*:%d",fStorageServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case CLIENT_COMMUNICATION_REQ:
        {
            fSockets[CLIENT_COMMUNICATION_REQ] =
            new socket_t(*fContexts[CLIENT_COMMUNICATION_REQ],ZMQ_REQ);
            try
            {
                fSockets[CLIENT_COMMUNICATION_REQ]->connect(Form("tcp://%s:%d",fStorageServer.c_str(), fStorageClientPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case CLIENT_COMMUNICATION_REP:
        {
            fSockets[CLIENT_COMMUNICATION_REP] =
            new socket_t(*fContexts[CLIENT_COMMUNICATION_REP],ZMQ_REP);
            try
            {
                fSockets[CLIENT_COMMUNICATION_REP]->bind(Form("tcp://*:%d",fStorageClientPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case EVENTS_SERVER_PUB:
        {
            fSockets[EVENTS_SERVER_PUB] =
            new socket_t(*fContexts[EVENTS_SERVER_PUB],ZMQ_PUB);
            try
            {
                fSockets[EVENTS_SERVER_PUB]->bind(Form("tcp://*:%d",fEventServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
            break;
        case EVENTS_SERVER_SUB:
        {
            fSockets[EVENTS_SERVER_SUB] =
            new socket_t(*fContexts[EVENTS_SERVER_SUB],ZMQ_SUB);
            fSockets[EVENTS_SERVER_SUB]->setsockopt(ZMQ_SUBSCRIBE,"",0);
            try
            {
                fSockets[EVENTS_SERVER_SUB]->connect(Form("tcp://%s:%d",fEventServer.c_str(),fEventServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
                
            }
        }
            break;
        case XML_PUB:
        {
            fSockets[XML_PUB] =
            new socket_t(*fContexts[XML_PUB],ZMQ_PUB);
            try
            {
                fSockets[XML_PUB]->bind(Form("tcp://*:%d",fXmlServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
	break;
    case ITS_POINTS_PUB:
        {
            fSockets[ITS_POINTS_PUB] =
            new socket_t(*fContexts[ITS_POINTS_PUB],ZMQ_PUB);
            try
            {
                fSockets[ITS_POINTS_PUB]->bind(Form("tcp://*:%d",fItsPointsServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
            }
        }
	break;
    case ITS_POINTS_SUB:
        {
            fSockets[ITS_POINTS_SUB] =
            new socket_t(*fContexts[ITS_POINTS_SUB],ZMQ_SUB);
            fSockets[ITS_POINTS_SUB]->setsockopt(ZMQ_SUBSCRIBE,"",0);
            try
            {
	      fSockets[ITS_POINTS_SUB]->connect(Form("tcp://%s:%d",fEventServer.c_str(),fItsPointsServerPort));
            }
            catch (const zmq::error_t& e)
            {
                cout<<"MANAGER -- "<<e.what()<<endl;
                return 0;
                
            }
        }
            break;
        default:break;
    }
    return 1;
}

void AliStorageEventManager::Send(vector<serverListStruct> list,storageSockets socket)
{
    //send size of the struct first
    int numberOfRecords = list.size();
    message_t message(20);
    snprintf ((char *)message.data(), 20 ,"%d",numberOfRecords);
    try{
        fSockets[socket]->send(message);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send vector -- "<<e.what()<<endl;
    }
    //if(numberOfRecords==0)return;
    message_t *tmpMessage = new message_t();
    
    try{
        fSockets[socket]->recv(tmpMessage);//empty message just to keep req-rep order
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send vector -- "<<e.what()<<endl;
    }
    // //prepare message with event's list
    // char *buffer = reinterpret_cast<char*> (&list[0]);
    // message_t *reply = new message_t((void*)buffer,
    // 		      sizeof(serverListStruct)*numberOfRecords,0);
    // fSockets[socket]->send(*reply);
    // if(reply){delete reply;}
    
    message_t reply(sizeof(serverListStruct)*numberOfRecords);
    memcpy(reply.data(), reinterpret_cast<const char*> (&list[0]), sizeof(serverListStruct)*numberOfRecords);
    
    try{
        fSockets[socket]->send(reply);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send vector -- "<<e.what()<<endl;
    }
    if(tmpMessage){delete tmpMessage;}
}

bool AliStorageEventManager::Send(struct serverRequestStruct *request,storageSockets socket,int timeout)
{
  
  if(timeout>=0)
    {
      pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}};
      if(poll (&items[0], 1, timeout)==0){
	cout<<"Event Manager -- couldn't send serverRequestStruct -- timeout"<<endl;
	return 0;
      }
    }

    char *buffer = (char*)(request);
    message_t *requestMessage = new message_t((void*)buffer,
                                              sizeof(struct serverRequestStruct)
                                              +sizeof(struct listRequestStruct)
                                              +sizeof(struct eventStruct),freeBuff);
    try{
        fSockets[socket]->send(*requestMessage);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send serverRequestStruct -- "<<e.what()<<endl;
	return 0;
    }
    return 1;
}

bool AliStorageEventManager::Send(struct clientRequestStruct *request,storageSockets socket,int timeout)
{
    if(timeout>=0){
      pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}};
      if(poll (&items[0], 1, timeout)==0){return 0;}
    }
    
    char *buffer = (char*)(request);
    message_t *requestMessage = new message_t((void*)buffer,sizeof(struct clientRequestStruct),freeBuff);
    
    try{
        fSockets[socket]->send(*requestMessage);
    }
    catch (const zmq::error_t& e)
    {
        cout<<"MANAGER -- "<<e.what()<<endl;
        cout<<e.num()<<endl;
        if(fSockets[socket]){delete fSockets[socket];fSockets[socket]=0;}
        
        CreateSocket(socket);
        delete requestMessage;
        return 0;
        
    }
    
    delete requestMessage;
    return 1;
}

void AliStorageEventManager::Send(long message,storageSockets socket)
{
    stringstream streamBuffer;
    streamBuffer << message;
    string stringBuffer = streamBuffer.str();
    char *buffer = (char*)stringBuffer.c_str();
    message_t *replyMessage = new message_t((void*)buffer,sizeof(stringBuffer),freeBuff);
    
    try{
        fSockets[socket]->send(*replyMessage);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send long -- "<<e.what()<<endl;
    }
    delete replyMessage;
    streamBuffer.str(string());
    streamBuffer.clear();
}

void AliStorageEventManager::Send(bool message,storageSockets socket)
{
    char *buffer;
    if(message==true)
    {
        buffer = (char*)("true");
    }
    else
    {
        buffer = (char*)("false");
    }
    message_t *replyMessage = new message_t((void*)buffer,sizeof(buffer),freeBuff);
    try{
        fSockets[socket]->send(*replyMessage);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send bool -- "<<e.what()<<endl;
    }
    delete replyMessage;
}

void AliStorageEventManager::Send(AliESDEvent *event, storageSockets socket)
{
  if(!event){return;}

    TMessage tmess(kMESS_OBJECT);
    tmess.Reset();
    tmess.WriteObject(event);
    //TMessage::EnableSchemaEvolutionForAll(kTRUE);
    
    int bufsize = tmess.Length();
    char* buf = (char*) malloc(bufsize * sizeof(char));
    memcpy(buf, tmess.Buffer(), bufsize);
    
    message_t message((void*)buf, bufsize, freeBuff);
    try{
        fSockets[socket]->send(message);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send AliESDEvent -- "<<e.what()<<endl;
    }
}

void AliStorageEventManager::Send(TFile *file, storageSockets socket)
{
  cout<<"sending tfile to socket:"<<endl;
  TTree *tree;

  if(file)
    {
      file->ls();
      cout<<"1"<<endl;
      //std::string *dfTitle = new std::string(file->GetListOfKeys()->Last()->GetTitle());
      string *name = new string(file->GetListOfKeys()->Last()->GetTitle());
      cout<<"treeTitle:"<<name->data()<<endl;
      
      tree = (TTree*)((TDirectoryFile*)file->Get(name->data()))->Get("TreeR");


	// TDirectoryFile *df = (TDirectoryFile*)file->Get(dfTitle->data());
	//if(df)
	//{
	//	  cout<<"directory file extracted"<<endl;
        //tree = (TTree*)df->Get("TreeR");
	cout<<"2"<<endl;
      tree->Branch("event",&name);
      cout<<"branch created:"<<endl;
      tree->Print();
      tree->Fill();
      cout<<"ttree filled with name"<<endl;
	  // delete df;
	  //}
      tree->Print();
    }
  else
    {
      tree = NULL;
      cout<<"file is empty"<<endl;
    }

  //TMessage::EnableSchemaEvolutionForAll(kTRUE);
  TMessage tmess(kMESS_OBJECT);
  tmess.Reset();
  tmess.WriteObject(tree);
  cout<<"file written to tmessage"<<endl;
 
  int bufsize = tmess.Length();
  char* buf = (char*) malloc(bufsize * sizeof(char));
  memcpy(buf, tmess.Buffer(), bufsize);
  cout<<"messaged copied to buffer"<<endl;
      
  message_t message((void*)buf, bufsize, freeBuff);
  cout<<"message_t created"<<endl;
  try{
    fSockets[socket]->send(message);
  }
  catch(const zmq::error_t &e)
    {
      cout<<"MANAGER -- send TFile -- "<<e.what()<<endl;
    }
  //if(tree){delete tree;}
}

void AliStorageEventManager::Send(struct recPointsStruct *files, storageSockets socket)
{
  for(int i=0;i<10;i++){Send(files->files[i],socket);}
}
void AliStorageEventManager::SendAsXml(AliESDEvent *event,storageSockets socket)
{
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
    message_t message(bufferString.size());
    memcpy (message.data(), bufferString.data(), bufferString.size());
    
    try{
        fSockets[socket]->send(message);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- send send xml -- "<<e.what()<<endl;
    }
    cout<<"xml sent"<<endl;
}

vector<serverListStruct> AliStorageEventManager::GetServerListVector(storageSockets socket, int timeout)
{
  pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}} ;
  if(timeout>=0){if(poll (&items[0], 1, timeout)==0){vector<serverListStruct> emptyVector;return emptyVector;}}

    //get size of the incomming message
    message_t sizeMessage;
    
    try{
        fSockets[socket]->recv(&sizeMessage);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- get vector -- "<<e.what()<<endl;
    }
    int numberOfRecords;
    istringstream iss(static_cast<char*>(sizeMessage.data()));
    iss >> numberOfRecords;
    
    if(numberOfRecords==0){cout<<"MANAGER -- list is empty"<<endl;}
    
    try{
        fSockets[socket]->send(*(new message_t()));//receive empty message just to keep req-rep order
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- get vector -- "<<e.what()<<endl;
    }
    //get list of events
    message_t *response = new message_t(sizeof(serverListStruct)*numberOfRecords);
    try{
        fSockets[socket]->recv(response);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- get vector -- "<<e.what()<<endl;
    }
    
    vector<serverListStruct> receivedList(static_cast<serverListStruct*>(response->data()), static_cast<serverListStruct*>(response->data()) + numberOfRecords);
    
    if (response) {delete response;}
    return receivedList;
}

zmq::message_t* AliStorageEventManager::GetMessage(storageSockets socket,int timeout)
{
  pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}} ;
  if(timeout>=0){
    try{(poll (&items[0], 1, timeout)==0);}
    catch(const zmq::error_t &e){
      cout<<"EVENT MANAGER -- GetEvent():"<<e.what()<<endl;
      return NULL;
    }
  }
  message_t* message = new message_t();

  try
  {
    fSockets[socket]->recv(message);
  }
  catch (const zmq::error_t& e)
  {
    cout<<"MANAGER -- "<<e.what()<<endl;
    return NULL;
  }
  return message;
}

AliESDEvent* AliStorageEventManager::GetEvent(storageSockets socket,int timeout)
{
    pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}} ;
    if(timeout>=0){
	try{(poll (&items[0], 1, timeout)==0);}
	catch(const zmq::error_t &e){
	  cout<<"EVENT MANAGER -- GetEvent():"<<e.what()<<endl;
	  return NULL;
	  }
    }
    message_t* message = new message_t();
    
    try
    {
        fSockets[socket]->recv(message);
    }
    catch (const zmq::error_t& e)
    {
        cout<<"MANAGER -- "<<e.what()<<endl;
        return NULL;
    }
    TBufferFile *mess = new TBufferFile(TBuffer::kRead,
                                        message->size()+sizeof(UInt_t),
                                        message->data());
    mess->InitMap();
    mess->ReadClass();// get first the class stored in message
    mess->SetBufferOffset(sizeof(UInt_t) + sizeof(kMESS_OBJECT));
    mess->ResetMap();
    
    AliESDEvent* data = (AliESDEvent*)(mess->ReadObjectAny(AliESDEvent::Class()));
    
    if (data)
    {
        data->GetStdContent();
        if(message){delete message;}
        return data;
    }
    else
    {
        if(message){delete message;}
        return NULL;
    }
}

TFile* AliStorageEventManager::GetFile(storageSockets socket,int timeout)
{
  cout<<"get file"<<endl;
    pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}} ;
    if(timeout>=0){
	try{(poll (&items[0], 1, timeout)==0);}
	catch(const zmq::error_t &e){
	  cout<<"EVENT MANAGER -- GetFile():"<<e.what()<<endl;
	    return NULL;
	  }
    }
    message_t* message = new message_t();
    cout<<"polling passed"<<endl;

    try
    {
      cout<<"waiting for file on socket:"<<socket<<endl;
        fSockets[socket]->recv(message);
    }
    catch (const zmq::error_t& e)
    {
        cout<<"MANAGER -- "<<e.what()<<endl;
        return NULL;
    }
    cout<<"createing buffer file"<<endl;
    TBufferFile *mess = new TBufferFile(TBuffer::kRead,
					message->size()+sizeof(UInt_t),
                                        message->data());

    //TMessage *mess = new TMessage();
    //mess->SetReadMode();
   
    mess->InitMap();
    mess->ReadClass();
    mess->SetBufferOffset(sizeof(UInt_t) + sizeof(kMESS_OBJECT));
    mess->ResetMap();

    //mess->ReadBuf(message->data(),message->size()+sizeof(UInt_t));
    
    cout<<"reading file from buffer"<<endl;
    TTree* tree = (TTree*)(mess->ReadObjectAny(TTree::Class()));
    //TFile* data = (TFile*)mess->ReadObject(TFile::Class());

    if(tree)
      {
	cout<<"received a tree:"<<endl;
	tree->Print();
	std::string *dfTitle = new std::string();
	tree->SetBranchAddress("event",&dfTitle);
	tree->GetEntry(0);

	cout<<"setting df's name to:"<<dfTitle->data()<<endl;
	TDirectoryFile *df = new TDirectoryFile(dfTitle->data(),dfTitle->data());
	df->Add(tree);
	cout<<"added tree to directory file"<<endl;
	
	TFile *file = new TFile();
	df->Write();

	cout<<"created file:"<<endl;
	file->ls();
	if(message){delete message;}
	if(df){delete df;}
	if(tree){delete tree;}
	return file;
      }
    else
      {
	cout<<"no tree found"<<endl;
	if(message){delete message;}
        return NULL;
      }
}

struct recPointsStruct* AliStorageEventManager::GetFiles(storageSockets socket,int timeout)
{
  struct recPointsStruct *files = new struct recPointsStruct;

  for(int i=0;i<10;i++)
    {
      files->files[i] = GetFile(socket,timeout);
    }
  return files;
}

struct serverRequestStruct* AliStorageEventManager::GetServerStruct(storageSockets socket)
{
    struct serverRequestStruct *request = new struct serverRequestStruct;
    message_t *requestMessage = new message_t();
    try{
        fSockets[socket]->recv(requestMessage);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- get serverRequestStruct -- "<<e.what()<<endl;
	request->messageType = -1;
	return request;
    }
    request = static_cast<struct serverRequestStruct*>(requestMessage->data());
    return request;
}

struct clientRequestStruct* AliStorageEventManager::GetClientStruct(storageSockets socket,int timeout)
{
    pollitem_t items[1] =  {{*fSockets[socket],0,ZMQ_POLLIN,0}} ;
    if(timeout>=0){if(poll (&items[0], 1, timeout)==0){return NULL;}}
    
    struct clientRequestStruct *request = new struct clientRequestStruct;
    message_t *requestMessage = new message_t();
    try{
        fSockets[socket]->recv(requestMessage);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- get clientRequestStruct -- "<<e.what()<<endl;
    }
    request = static_cast<struct clientRequestStruct*>(requestMessage->data());
    return request;
}

bool AliStorageEventManager::GetBool(storageSockets socket)
{
    message_t *response = new message_t();
    try{
        fSockets[socket]->recv(response);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- get bool -- "<<e.what()<<endl;
    }
    char *result = (char*)response->data();
    
    if(!strcmp("true",result)){return true;}
    else{return false;}
}

long AliStorageEventManager::GetLong(storageSockets socket)
{
    message_t *responseMessage = new message_t();
    try{
        fSockets[socket]->recv(responseMessage);
    }
    catch(const zmq::error_t &e)
    {
        cout<<"MANAGER -- get long -- "<<e.what()<<endl;
    }
    
    long result = 0;
    
    if(responseMessage)
    {
        result = (long)atoi(static_cast<char*>(responseMessage->data()));
        delete responseMessage;
    }
    return result;
}

