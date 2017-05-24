// a helper library for using ZMQ with ROOT, focussed on multipart messaging
// blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
// Copyright (C) 2015 Goethe University Frankfurt
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "AliZMQhelpers.h"

#include "zmq.h"
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <vector>

#include "TObjArray.h"
#include "TCollection.h"
#include "TStreamerInfo.h"
#include "TClass.h"
#include "TSystem.h"
#include "TList.h"

#include <netinet/in.h>

using namespace AliZMQhelpers;

using namespace AliZMQhelpers;

//init the shared context to null
void* AliZMQhelpers::gZMQcontext = NULL;

const UInt_t AliZMQhelpers::BaseDataTopic::fgkMagicNumber = CharArr2uint32("O2O2");
const ULong64_t AliZMQhelpers::DataTopic::fgkDataTopicDescription = CharArr2uint64("DataHDR");
const UInt_t AliZMQhelpers::DataTopic::fgkTopicSerialization = CharArr2uint64("NONE");
const ULong64_t AliZMQhelpers::kSerializationROOT = CharArr2uint64("ROOT   ");
const ULong64_t AliZMQhelpers::kSerializationNONE = CharArr2uint64("NONE   ");

const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeStreamerInfos("ROOTSTRI","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeInfo("INFO____","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeConfig("CONFIG__","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeTObject("ROOTTOBJ","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeTH1("ROOTHIST","***\n",0);

const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeStreamerInfos("ROOTSTRI","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeInfo("INFO____","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeConfig("CONFIG__","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeTObject("ROOTTOBJ","***\n",0);
const AliZMQhelpers::DataTopic AliZMQhelpers::kDataTypeTH1("ROOTHIST","***\n",0);

//_______________________________________________________________________________________
void* AliZMQhelpers::alizmq_context()
{
  if (!AliZMQhelpers::gZMQcontext) AliZMQhelpers::gZMQcontext=zmq_ctx_new();
  return AliZMQhelpers::gZMQcontext;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_detach (void *self, const char *endpoints, bool serverish)
{
    assert (self);
    if (!endpoints)
        return 0;

    if (strlen(endpoints)<2)
        return 0;

    //  We hold each individual endpoint here
    char endpoint [256];
    while (*endpoints) {
        const char *delimiter = strchr (endpoints, ',');
        if (!delimiter)
            delimiter = endpoints + strlen (endpoints);
        if (delimiter - endpoints > 255)
            return -1;
        memcpy (endpoint, endpoints, delimiter - endpoints);
        endpoint [delimiter - endpoints] = 0;

        int rc;
        if (endpoint [0] == '@')
            rc = zmq_bind (self, endpoint + 1);
        else
        if (endpoint [0] == '>' || endpoint [0] == '-' || endpoint [0] == '+' )
            rc = zmq_connect (self, endpoint + 1);
        else
        if (serverish)
            rc = zmq_unbind (self, endpoint);
        else
            rc = zmq_disconnect (self, endpoint);

        if (rc == -1)
            return -1;          //  Bad endpoint syntax

        if (*delimiter == 0)
            break;
        endpoints = delimiter + 1;
    }
    return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_attach (void *self, const char *endpoints, bool serverish)
{
    assert (self);
    if (!endpoints)
        return 0;
    if (strlen(endpoints)<2)
        return 0;

    //  We hold each individual endpoint here
    char endpoint [256];
    while (*endpoints) {
        const char *delimiter = strchr (endpoints, ',');
        if (!delimiter)
            delimiter = endpoints + strlen (endpoints);
        if (delimiter - endpoints > 255)
            return -1;
        memcpy (endpoint, endpoints, delimiter - endpoints);
        endpoint [delimiter - endpoints] = 0;

        int rc;
        if (endpoint [0] == '@')
            rc = zmq_bind (self, endpoint + 1);
        else
        if (endpoint [0] == '>' || endpoint [0] == '-' || endpoint [0] == '+' )
            rc = zmq_connect (self, endpoint + 1);
        else
        if (serverish)
            rc = zmq_bind (self, endpoint);
        else
            rc = zmq_connect (self, endpoint);

        if (rc == -1)
            return -1;          //  Bad endpoint syntax

        if (*delimiter == 0)
            break;
        endpoints = delimiter + 1;
    }
    return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_socket_state(void* socket)
{
  int events=0;
  size_t len = sizeof(events);
  zmq_getsockopt(socket, ZMQ_EVENTS, &events, &len);
  return events;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_socket_type(void* socket)
{
  //get the type of the socket
  int type=-1;
  size_t typeLen=sizeof(type);
  int rc=0;
  rc = zmq_getsockopt(socket, ZMQ_TYPE, &type, &typeLen);
  if (rc<0) return rc;
  return type;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_socket_type(std::string config)
{
  if (config.compare(0,3,"PUB")==0) return ZMQ_PUB;
  else if (config.compare(0,3,"SUB")==0) return ZMQ_SUB;
  else if (config.compare(0,3,"REP")==0) return ZMQ_REP;
  else if (config.compare(0,3,"REQ")==0) return ZMQ_REQ;
  else if (config.compare(0,4,"PUSH")==0) return ZMQ_PUSH;
  else if (config.compare(0,4,"PULL")==0) return ZMQ_PULL;
  else if (config.compare(0,6,"DEALER")==0) return ZMQ_DEALER;
  else if (config.compare(0,6,"ROUTER")==0) return ZMQ_ROUTER;
  else if (config.compare(0,6,"STREAM")==0) return ZMQ_STREAM;
  else if (config.compare(0,4,"PAIR")==0) return ZMQ_PAIR;
  else if (config.compare(0,4,"XSUB")==0) return ZMQ_XSUB;
  else if (config.compare(0,4,"XPUB")==0) return ZMQ_XPUB;

  //printf("Invalid socket type %s\n", config.c_str());
  return -1;
}

//_______________________________________________________________________________________
const char* AliZMQhelpers::alizmq_socket_name(int socketType)
{
  switch (socketType)
  {
    case ZMQ_PUB: return "PUB";
    case ZMQ_SUB: return "SUB";
    case ZMQ_REP: return "REP";
    case ZMQ_REQ: return "REQ";
    case ZMQ_PUSH: return "PUSH";
    case ZMQ_PULL: return "PULL";
    case ZMQ_DEALER: return "DEALER";
    case ZMQ_ROUTER: return "ROUTER";
    case ZMQ_STREAM: return "STREAM";
    case ZMQ_PAIR: return "PAIR";
    case ZMQ_XPUB: return "XPUB";
    case ZMQ_XSUB: return "XSUB";
    default: return "INVALID";
  }
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_socket_close(void*& socket, int linger)
{
  zmq_setsockopt(socket, ZMQ_LINGER, &linger, sizeof(linger));
  int rc = zmq_close(socket);
  if (rc>=0) socket = NULL;
  return rc;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_socket_init(void*& socket, void* context, std::string config, int timeout, int highWaterMark)
{
  int rc = 0;
  int zmqSocketMode = 0;
  std::string zmqEndpoints = "";

  size_t configStartPos = config.find_first_not_of(" \t\n");
  size_t configEndPos = config.find_last_not_of(" \t\n");
  if (configStartPos!=std::string::npos && configEndPos!=std::string::npos)
  { config = config.substr(configStartPos,configEndPos-configStartPos+1); }

  if (config.empty()) {
    alizmq_socket_close(socket);
    socket = NULL;
    return 999999;
  }

  std::size_t found = config.find_first_of("@>-+");
  if (found == 0)
  {
    //printf("misformed socket config string %s\n", config.c_str());
    return -1;
  }

  zmqSocketMode = alizmq_socket_type(config);

  if (found!=std::string::npos)
  { zmqEndpoints=config.substr(found,std::string::npos); }

  bool newSocket=true;
  //init the ZMQ socket
  if (socket)
  {
    newSocket=false;
    int lingerValue = 10;
    rc = zmq_setsockopt(socket, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
    if (rc!=0)
    {
      //printf("cannot set linger 0 on socket before closing\n");
      return -2;
    }
    rc = zmq_close(socket);
    if (rc!=0)
    {
      //printf("zmq_close() says: %s\n",zmq_strerror(errno));
      return -3;
    }
  }

  socket  = zmq_socket(context, zmqSocketMode);

  //set socket options
  rc += zmq_setsockopt(socket, ZMQ_RCVHWM, &highWaterMark, sizeof(highWaterMark));
  rc += zmq_setsockopt(socket, ZMQ_SNDHWM, &highWaterMark, sizeof(highWaterMark));
  rc += zmq_setsockopt(socket, ZMQ_RCVTIMEO, &timeout, sizeof(timeout));
  rc += zmq_setsockopt(socket, ZMQ_SNDTIMEO, &timeout, sizeof(timeout));
  if (rc!=0)
  {
    //printf("cannot set socket options\n");
    return -4;
  }

  //by default subscribe to everything if we happen to be SUB
  if (zmqSocketMode == ZMQ_SUB) {
    rc = zmq_setsockopt(socket, ZMQ_SUBSCRIBE, "", 0);
  }

  //connect the socket to the endpoints
  //when reinitializing sometimes it is not possible to bind the same port again fast,
  //we need to retry a few times if we are indeed reconnecting, otherwise we just exit
  int i=100;
  while (i-->0)
  {
    rc = alizmq_attach(socket,  zmqEndpoints.c_str() );
    if ( rc==0 || newSocket ) break;
    usleep(100000);
  }
  if (rc!=0)
  {
    //printf("cannot attach to %s\n",zmqEndpoints.c_str());
    return -5;
  }

  int lingerValue = 0;
  rc += zmq_setsockopt(socket, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
  //printf("socket mode: %s, endpoints: %s\n",alizmq_socket_name(zmqSocketMode), zmqEndpoints.c_str());

  //reset the object containers
  return zmqSocketMode;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_add(aliZMQmsg* message, const std::string& topicString, const std::string& data)
{
  //add a frame to the mesage
  int rc = 0;

  DataTopic topic;
  topic.SetID(topicString.c_str());

  //prepare topic msg
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, sizeof(topic));
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg),&topic,sizeof(topic));

  //prepare data msg
  zmq_msg_t* dataMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( dataMsg, data.size());
  if (rc<0) {
    zmq_msg_close(topicMsg);
    zmq_msg_close(dataMsg);
    delete topicMsg;
    delete dataMsg;
    return -1;
  }
  memcpy(zmq_msg_data(dataMsg),data.data(),data.size());

  //add the frame to the message
  message->push_back(std::make_pair(topicMsg,dataMsg));
  return message->size();
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_add(aliZMQmsg* message, const DataTopic* topic, void* data, int size)
{
  //add a frame to the mesage
  int rc = 0;

  //prepare topic msg
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, sizeof(*topic));
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg),topic,sizeof(*topic));

  //prepare data msg
  zmq_msg_t* dataMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( dataMsg, size);
  if (rc<0) {
    zmq_msg_close(topicMsg);
    zmq_msg_close(dataMsg);
    delete topicMsg;
    delete dataMsg;
    return -1;
  }
  memcpy(zmq_msg_data(dataMsg),data,size);

  //add the frame to the message
  message->push_back(std::make_pair(topicMsg,dataMsg));
  return message->size();
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_add(aliZMQmsg* message, const DataTopic* topic, const std::string& data)
{
  //add a frame to the mesage
  int rc = 0;

  //prepare topic msg
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, sizeof(*topic));
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg),topic,sizeof(*topic));

  //prepare data msg
  zmq_msg_t* dataMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( dataMsg, data.size());
  if (rc<0) {
    zmq_msg_close(topicMsg);
    zmq_msg_close(dataMsg);
    delete topicMsg;
    delete dataMsg;
    return -1;
  }
  memcpy(zmq_msg_data(dataMsg),data.data(),data.size());

  //add the frame to the message
  message->push_back(std::make_pair(topicMsg,dataMsg));
  return message->size();
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_add(aliZMQmsg* message, DataTopic* topic, TObject* object,
                   int compression, aliZMQrootStreamerInfo* streamers)
{
  //add a frame to the mesage
  int rc = 0;

  //we are definitely serializing ROOT objects here...
  topic->SetSerialization(kSerializationROOT);

  //prepare topic msg
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, sizeof(*topic));
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg), topic, sizeof(*topic));

  //prepare data msg
  ZMQTMessage* tmessage = ZMQTMessage::Stream(object, compression, 0, streamers);
  if (!tmessage) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }

  if (streamers) {
    alizmq_update_streamerlist(streamers, tmessage->GetStreamerInfos());
  }

  zmq_msg_t* dataMsg = new zmq_msg_t;
  rc = zmq_msg_init_data( dataMsg, tmessage->Buffer(), tmessage->Length(),
       alizmq_deleteTObject, tmessage);
  if (rc<0) {
    zmq_msg_close(topicMsg);
    zmq_msg_close(dataMsg);
    delete topicMsg;
    delete dataMsg;
    return -1;
  }

  //add the frame to the message
  message->push_back(std::make_pair(topicMsg,dataMsg));

  return message->size();
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_send(aliZMQmsg* message, void* socket, int flagsUser)
{
  int nBytes = 0;
  int rc = 0;
  int flags = flagsUser | ZMQ_SNDMORE;

  for (aliZMQmsg::iterator i=message->begin(); i!=message->end(); ++i)
  {
    zmq_msg_t* topic = i->first;
    zmq_msg_t* data = i->second;

    rc = zmq_msg_send(topic, socket, flags);
    if (rc<0) break;
    nBytes+=rc;

    if (&*i == &*message->rbegin()) flags=flagsUser; //last frame
    rc = zmq_msg_send(data, socket, flags);
    if (rc<0) break;
    nBytes+=rc;
  }
  if (rc<0) nBytes=rc;
  return nBytes;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_send(std::string topicString, std::string data, void* socket, int flags)
{
  int nBytes = 0;
  int rc = 0;

  DataTopic topic;
  topic.SetID(topicString.c_str());

  zmq_msg_t topicMsg;
  zmq_msg_init_size(&topicMsg, sizeof(topic));
  memcpy(zmq_msg_data(&topicMsg), &topic, zmq_msg_size(&topicMsg));
  rc = zmq_msg_send(&topicMsg, socket, ZMQ_SNDMORE);
  if (rc<0)
  {
    //printf("unable to send topicString: %s\n", topicString.c_str());
    zmq_msg_close(&topicMsg);
    return rc;
  }
  nBytes+=rc;

  zmq_msg_t dataMsg;
  zmq_msg_init_size(&dataMsg, data.size());
  memcpy(zmq_msg_data(&dataMsg), data.data(), zmq_msg_size(&dataMsg));
  rc = zmq_msg_send(&dataMsg, socket, flags);
  if (rc<0)
  {
    //printf("unable to send data: %s\n", data.c_str());
    zmq_msg_close(&dataMsg);
    return rc;
  }
  nBytes+=rc;
  return nBytes;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_prepend_streamer_infos(aliZMQmsg* message, aliZMQrootStreamerInfo* streamers)
{
  //prepend the streamer info to the message as first block.
  int rc = 0;

  DataTopic topic = kDataTypeStreamerInfos;
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, sizeof(topic));
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg), &topic, sizeof(topic));

  //prepare data msg
  TObjArray listOfInfos;
  for (aliZMQrootStreamerInfo::const_iterator i=streamers->begin(); i!=streamers->end(); ++i) {
    listOfInfos.Add(*i);
  }
  ZMQTMessage* tmessage = ZMQTMessage::Stream(&listOfInfos, 1); //compress
  zmq_msg_t* dataMsg = new zmq_msg_t;
  rc = zmq_msg_init_data( dataMsg, tmessage->Buffer(), tmessage->Length(),
       alizmq_deleteTObject, tmessage);
  if (rc<0) {
    zmq_msg_close(topicMsg);
    zmq_msg_close(dataMsg);
    delete topicMsg;
    delete dataMsg;
    return -1;
  }

  message->insert(message->begin(),std::make_pair(topicMsg,dataMsg));

  return 0;
}

//_______________________________________________________________________________________
void AliZMQhelpers::alizmq_update_streamerlist_from_object(aliZMQrootStreamerInfo* streamers, TObject* object)
{
  //update the list of streamer infos with the streamers needed for object
  ZMQTMessage* tmessage = ZMQTMessage::Stream(object, 0, 0, streamers);
  if (streamers) {
    alizmq_update_streamerlist(streamers, tmessage->GetStreamerInfos());
  }

  delete tmessage;
}

//_______________________________________________________________________________________
void AliZMQhelpers::alizmq_update_streamerlist(aliZMQrootStreamerInfo* streamers, const TCollection* newStreamers)
{
  //update the list of streamers used
  if (!streamers) return;
  if (!newStreamers) return;

  TIter nextStreamer(newStreamers);
  size_t i = 0;
  while (TVirtualStreamerInfo* info = static_cast<TVirtualStreamerInfo*>(nextStreamer())) {
    const char* name = info->GetName();
    int version = info->GetClassVersion();
    bool found=false;
    for (aliZMQrootStreamerInfo::iterator i=streamers->begin(); i!=streamers->end(); ++i)
    {
      const char* existingName = (*i)->GetName();
      int existingVersion = (*i)->GetClassVersion();
      if (name == existingName && version==existingVersion) {
        found=true;
        break;
      }
    }
    if (!found) {
      streamers->push_back(info);
    }
    ++i;
  }
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_init_streamer_infos(aliZMQmsg::iterator it)
{
  int rc = 0;
  TObject* obj = NULL;
  rc = alizmq_msg_iter_data(it,obj);
  TObjArray* pSchemas = dynamic_cast<TObjArray*>(obj);
  if (!pSchemas) {
    return -1;
  }

  pSchemas->SetOwner(kTRUE);

  for (int i=0; i<pSchemas->GetEntriesFast(); i++) {
    if (pSchemas->At(i)) {
      TStreamerInfo* pSchema=dynamic_cast<TStreamerInfo*>(pSchemas->At(i));
      if (pSchema) {
        int version=pSchema->GetClassVersion();
        TClass* pClass=TClass::GetClass(pSchema->GetName());
        if (pClass) {
          if (pClass->GetClassVersion()==version) {
            //AliDebug(0,Form("skipping schema definition %d version %d to class %s as this is the native version", i, version, pSchema->GetName()));
            continue;
          }
          TObjArray* pInfos=const_cast<TObjArray*>(pClass->GetStreamerInfos());
          if (pInfos /*&& version<pInfos->GetEntriesFast()*/) {
            TVirtualStreamerInfo* pInfo = dynamic_cast<TVirtualStreamerInfo*>(pInfos->At(version));
            if (pInfo==NULL) {
              pSchema->SetClass(pClass);
              pSchema->BuildOld();
              pInfos->AddAtAndExpand(pSchema, version);
              pSchemas->Remove(pSchema);
              printf("adding %s %i\n",pSchema->GetName(),version);
              //AliDebug(0,Form("adding schema definition %d version %d to class %s", i, version, pSchema->GetName()));
            } else {
              if (pInfo && pInfo->GetClassVersion()==version) {
                //AliDebug(0,Form("schema definition %d version %d already available in class %s, skipping ...", i, version, pSchema->GetName()));
              } else {
                //AliError(Form("can not verify version for already existing schema definition %d (%s) version %d: version of existing definition is %d", i, pSchema->GetName(), version, pInfo?pInfo->GetClassVersion():-1));
              }
            }
          } else {
            //AliError(Form("skipping schema definition %d (%s), unable to set version %d in info array of size %d", i, pSchema->GetName(), version, pInfos?pInfos->GetEntriesFast():-1));
          }
        } else {
          //AliError(Form("skipping schema definition %d (%s), unable to find class", i, pSchema->GetName()));
        }
      } else {
        //AliError(Form("skipping schema definition %d, not of TStreamerInfo", i));
      }
    }
  }

  delete pSchemas; //this destroys remaining schemas as pSchemas is set owner
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_send(DataTopic& topic, TObject* object, void* socket, int flags,
                    int compression, aliZMQrootStreamerInfo* streamers)
{
  int nBytes = 0;
  int rc = 0;

  ZMQTMessage* tmessage = ZMQTMessage::Stream(object, compression, 0, streamers);
  zmq_msg_t dataMsg;
  rc = zmq_msg_init_data( &dataMsg, tmessage->Buffer(), tmessage->Length(),
      alizmq_deleteTObject, tmessage);

  if (streamers) {
    alizmq_update_streamerlist(streamers, tmessage->GetStreamerInfos());
  }

  //we are definitely serializing ROOT objects here...
  topic.SetSerialization(kSerializationROOT);

  //then send the object topic
  rc = zmq_send( socket, &topic, sizeof(topic), ZMQ_SNDMORE );
  if (rc<0)
  {
    zmq_msg_close(&dataMsg);
    //printf("unable to send topic: %s %s\n", topic.Description().c_str(), zmq_strerror(errno));
    return rc;
  }
  nBytes += rc;

  //send the object itself
  rc = zmq_msg_send(&dataMsg, socket, flags);
  if (rc<0)
  {
    //printf("unable to send data: %s %s\n", tmessage->GetName(), zmq_strerror(errno));
    zmq_msg_close(&dataMsg);
    return rc;
  }
  nBytes += rc;

  return nBytes;
}

//______________________________________________________________________________
int AliZMQhelpers::alizmq_msg_send(const DataTopic& topic, const std::string& data, void* socket, int flags)
{
  int nBytes = 0;
  int rc = 0;
  rc = zmq_send( socket, &topic, sizeof(topic), ZMQ_SNDMORE );
  if (rc<0)
  {
    //printf("unable to send topic: %s %s\n", topic.Description().c_str(), zmq_strerror(errno));
    return rc;
  }
  nBytes += rc;

  zmq_msg_t dataMsg;
  zmq_msg_init_size(&dataMsg, data.size());
  memcpy(zmq_msg_data(&dataMsg), data.data(), zmq_msg_size(&dataMsg));
  rc = zmq_msg_send(&dataMsg, socket, flags);
  if (rc<0)
  {
    //printf("unable to send data: %s\n", data.c_str());
    zmq_msg_close(&dataMsg);
    return rc;
  }
  nBytes = rc;
  return nBytes;
}

//______________________________________________________________________________
void AliZMQhelpers::alizmq_deleteTObject(void*, void* object)
{
  //delete the TBuffer, for use in zmq_msg_init_data(...) only.
  //printf("deleteObject called! ZMQ just sent and destroyed the message!\n");
  TObject* tobject = static_cast<TObject*>(object);
  delete tobject;
}

//______________________________________________________________________________
void AliZMQhelpers::alizmq_deleteTopic(void*, void* object)
{
  //delete the TBuffer, for use in zmq_msg_init_data(...) only.
  //printf("deleteObject called! ZMQ just sent and destroyed the message!\n");
  DataTopic* topic = static_cast<DataTopic*>(object);
  delete topic;
}


//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_close(aliZMQmsg* message)
{
  int rc = 0;
  for (aliZMQmsg::iterator i=message->begin(); i!=message->end(); ++i)
  {
    int rc1 = zmq_msg_close(i->first);
    delete i->first; i->first=NULL;
    int rc2 = zmq_msg_close(i->second);
    delete (i->second); i->second=NULL;
  }
  message->clear();
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_check(aliZMQmsg::iterator it, const DataTopic& topic)
{
  DataTopic actualTopic;
  alizmq_msg_iter_topic(it, actualTopic);
  if (actualTopic == topic) return 0;
  return 1;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_check_id(aliZMQmsg::iterator it, const DataTopic& topic)
{
  DataTopic actualTopic;
  alizmq_msg_iter_topic(it, actualTopic);
  if (actualTopic.GetID() == topic.GetID()) return 0;
  return 1;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_check_id(aliZMQmsg::iterator it, const std::string& topic)
{
  DataTopic actualTopic;
  alizmq_msg_iter_topic(it, actualTopic);
  std::string topicID = actualTopic.GetID();
  return topicID.compare(0,topic.size(),topic);
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_topic(aliZMQmsg::iterator it, std::string& topic)
{
  zmq_msg_t* message = it->first;
  void* buf = zmq_msg_data(message);
  DataTopic* intopic = DataTopic::Get(buf);
  if (!intopic) return 1;
  char* arr = reinterpret_cast<char*>(&intopic->fDataDescription[1]);
  size_t nbytes = sizeof(&intopic->fDataDescription[1]);
  topic.assign(arr,nbytes);
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_data(aliZMQmsg::iterator it, std::string& data)
{
  zmq_msg_t* message = it->second;
  data.assign((char*)zmq_msg_data(message),zmq_msg_size(message));
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_data(aliZMQmsg::iterator it, void*& buffer, size_t& size)
{
  zmq_msg_t* message = it->second;
  buffer = zmq_msg_data(message);
  size = zmq_msg_size(message);
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_topic(aliZMQmsg::iterator it, DataTopic& topic)
{
  zmq_msg_t* message = it->first;
  void* buf = zmq_msg_data(message);
  DataTopic* intopic = DataTopic::Get(buf);
  if (!intopic) return 1;
  memcpy(&topic, zmq_msg_data(message),std::min(zmq_msg_size(message),sizeof(topic)));
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_data(aliZMQmsg::iterator it, TObject*& object)
{
  zmq_msg_t* message = it->second;
  size_t size = zmq_msg_size(message);
  void* data = zmq_msg_data(message);

  object = ZMQTMessage::Extract(data, size);
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_copy(aliZMQmsg* dst, aliZMQmsg* src)
{
  //copy (append) src to dst
  int numberOfMessages=0;
  for (aliZMQmsg::iterator i=src->begin(); i!=src->end(); ++i)
  {
    int rc=0;
    zmq_msg_t* topicMsg = new zmq_msg_t;
    rc  = zmq_msg_init(topicMsg);
    rc += zmq_msg_copy(topicMsg, i->first);
    if (rc<0) numberOfMessages=-1;

    zmq_msg_t* dataMsg = new zmq_msg_t;
    rc  = zmq_msg_init(dataMsg);
    rc += zmq_msg_copy(dataMsg, i->second);
    if (rc<0) numberOfMessages=-1;

    if (numberOfMessages<0)
    {
      zmq_msg_close(topicMsg);
      zmq_msg_close(dataMsg);
      delete topicMsg;
      delete dataMsg;
      return -1;
    }

    dst->push_back(std::make_pair(topicMsg, dataMsg));
    numberOfMessages++;
  }
  return numberOfMessages;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_recv(aliZMQmsg* message, void* socket, int flags)
{
  int rc = -1;
  int nBytes=0;
  while (true)
  {
    zmq_msg_t* topicMsg = new zmq_msg_t;
    rc = zmq_msg_init(topicMsg);
    rc = zmq_msg_recv(topicMsg, socket, flags);
    if (!zmq_msg_more(topicMsg) || rc<0)
    {
      zmq_msg_close(topicMsg);
      delete topicMsg;
      nBytes=-1;
      break;
    }
    nBytes+=rc;

    zmq_msg_t* dataMsg = new zmq_msg_t;
    rc = zmq_msg_init(dataMsg);
    rc = zmq_msg_recv(dataMsg, socket, flags);
    if (rc<0)
    {
      zmq_msg_close(topicMsg);
      zmq_msg_close(dataMsg);
      delete topicMsg;
      delete dataMsg;
      nBytes=-1;
      break;
    }
    nBytes+=rc;

    message->push_back(std::make_pair(topicMsg,dataMsg));

    int more=0;
    size_t moreLength = sizeof(more);
    rc = zmq_getsockopt(socket, ZMQ_RCVMORE, &more, &moreLength);
    if (!more) break;
  }
  return nBytes;
}

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_recv(aliZMQmsgStr* message, void* socket, int flags)
{
  int rc = -1;
  int nBytes=0;
  while (true)
  {
    zmq_msg_t topicMsg;
    rc = zmq_msg_init(&topicMsg);
    rc = zmq_msg_recv(&topicMsg, socket, flags);
    if (!zmq_msg_more(&topicMsg) || nBytes<0)
    {
      zmq_msg_close(&topicMsg);
      nBytes=-1;
      break;
    }
    nBytes+=rc;

    zmq_msg_t dataMsg;
    rc = zmq_msg_init(&dataMsg);
    rc = zmq_msg_recv(&dataMsg, socket, flags);
    if (nBytes<0)
    {
      zmq_msg_close(&topicMsg);
      zmq_msg_close(&dataMsg);
      nBytes=-1;
      break;
    }
    nBytes+=rc;

    std::string data;
    std::string topic;
    topic.assign(static_cast<char*>(zmq_msg_data(&topicMsg)), zmq_msg_size(&topicMsg));
    data.assign(static_cast<char*>(zmq_msg_data(&dataMsg)), zmq_msg_size(&dataMsg));

    message->push_back(std::make_pair(topic,data));

    rc = zmq_msg_close(&topicMsg);
    rc = zmq_msg_close(&dataMsg);

    int more=0;
    size_t moreLength = sizeof(more);
    rc = zmq_getsockopt(socket, ZMQ_RCVMORE, &more, &moreLength);
    if (!more) break;
  }
  return nBytes;
}

//______________________________________________________________________________
//tokenize a std::string
using namespace std;
vector<string> AliZMQhelpers::TokenizeString(const string input, const string delimiters)
{
  vector<string> output;
  output.reserve(10);
  size_t start = 0;
  size_t end = input.find_first_of(delimiters);
  if (end == string::npos)
  {
    output.push_back(input.substr(start, input.length()));
  }
  else do
  {
    output.push_back(input.substr(start, end-start));
    start = ++end;
    end = input.find_first_of(delimiters, start);
    if (end == string::npos)
    {
      output.push_back(input.substr(start, input.length()));
    }
  } while (end != string::npos);
  return output;
}

//______________________________________________________________________________
//tokenize a string delimited by semi-colons and return a map
//of key value pairs, "KEY=VALUE;key2=value2"
stringMap AliZMQhelpers::ParseParamString(const string paramString)
{
  vector<string> tokens = TokenizeString( paramString, ";");
  stringMap output;
  for (vector<string>::iterator i=tokens.begin(); i!=tokens.end(); ++i)
  {
    if (i->empty()) continue;
    size_t pos = i->find_first_of("=");
    if (pos == string::npos)
    {
      output[*i] = "";
    }
    else
    {
      output[i->substr(0,pos)] = i->substr(pos+1,string::npos);
    }
  }
  return output;
}

//______________________________________________________________________________
//a much faster version of the param string parser - just gives you one value
std::string AliZMQhelpers::GetParamString(const std::string param, const std::string paramstring)
{
  size_t start = paramstring.find(param+"=");
  if (start==std::string::npos) return "";
  start = paramstring.find_first_of("=",start);
  if (start==std::string::npos) return "";
  size_t end = paramstring.find_first_of(";", start);
  start++;
  return paramstring.substr(start,end-start);
}


//______________________________________________________________________________
int AliZMQhelpers::LoadROOTlibs(string libString, bool verbose)
{
  //load specified libraries, return number of loaded libs or negative in case of error
  int nLibs=0;
  vector<string> libs = TokenizeString(libString, ",");
  for (vector<string>::iterator i=libs.begin(); i!=libs.end(); ++i)
  {
    if (verbose) { printf("...loading %s", i->c_str()); }
    if (gSystem->Load(i->c_str())<0) {
      nLibs = -1;
      if (verbose) { printf(" ERROR!\n"); }
      break;
    }
    if (verbose) { printf(" OK\n"); }
    nLibs++;
  }
  return nLibs;
}

//______________________________________________________________________________
bool AliZMQhelpers::Topicncmp(const char* topic,
                              const char* reference,
                              int topicSize,
                              int referenceSize)
{
  for (int i=0; i<((topicSize<referenceSize)?topicSize:referenceSize); i++)
  {
    if (!(topic[i]=='*' || reference[i]=='*' ||
          topic[i]=='\0' || reference[i]=='\0' ||
          topic[i]==reference[i])) {
      return false;
    }
  }
  return true;
}

//______________________________________________________________________________
AliZMQhelpers::BaseDataTopic::BaseDataTopic()
  : fMagicNumber(fgkMagicNumber)
  , fHeaderSize(sizeof(BaseDataTopic))
  , fFlags(0)
  , fBaseHeaderVersion(1)
  , fHeaderDescription(0)
  , fHeaderSerialization(0)
{
}

//______________________________________________________________________________
AliZMQhelpers::BaseDataTopic::BaseDataTopic(UInt_t size, ULong64_t desc, ULong64_t seri)
  : fMagicNumber(fgkMagicNumber)
  , fHeaderSize(size)
  , fFlags(0)
  , fBaseHeaderVersion(1)
  , fHeaderDescription(desc)
  , fHeaderSerialization(seri)
{
}

//__________________________________________________________________________________________________
//helper function to print a hex/ASCII dump of some memory
void AliZMQhelpers::hexDump (const char* desc, const void* addr, int len)
{
  int i;
  unsigned char buff[17];       // stores the ASCII data
  const unsigned char *pc = static_cast<const unsigned char*>(addr);     // cast to make the code cleaner.

  // Output description if given.
  if (desc != NULL)
    printf ("%s:\n", desc);

  // In case of null pointer addr
  if (addr==nullptr) {printf("  nullptr\n"); return;}

  // Process every byte in the data.
  for (i = 0; i < len; i++) {
    // Multiple of 16 means new line (with line offset).
    if ((i % 16) == 0) {
      // Just don't print ASCII for the zeroth line.
      if (i != 0)
        printf ("  %s\n", buff);

      // Output the offset.
      //printf ("  %04x ", i);
      printf ("  %p ", &pc[i]);
    }

    // Now the hex code for the specific character.
    printf (" %02x", pc[i]);

    // And store a printable ASCII character for later.
    if ((pc[i] < 0x20) || (pc[i] > 0x7e))
      buff[i % 16] = '.';
    else
      buff[i % 16] = pc[i];
    buff[(i % 16) + 1] = '\0';
  }

  // Pad out last line if not exactly 16 characters.
  while ((i % 16) != 0) {
    printf ("   ");
    i++;
  }

  // And print the final ASCII bit.
  printf ("  %s\n", buff);
}

//__________________________________________________________________________________________________
void AliZMQhelpers::hexDump (aliZMQmsg* message, size_t maxsize)
{
  for (aliZMQmsg::iterator i=message->begin(); i!=message->end(); ++i)
  {
    zmq_msg_t* topic = i->first;
    zmq_msg_t* data = i->second;

    printf("______________________________________________________________________________\n");
    hexDump("header",zmq_msg_data(topic),zmq_msg_size(topic));
    hexDump("payload",zmq_msg_data(data),maxsize);
  }
}

//__________________________________________________________________________________________________
TObject* AliZMQhelpers::ZMQTMessage::Extract(const void* pBuffer, unsigned bufferSize, unsigned verbosity)
{
   /// Helper function to extract an object from a buffer.
   /// The returned object must be cleaned by the caller
  if (!pBuffer || bufferSize<sizeof(Int_t)) {
    return NULL;
  }

  Int_t firstWord=*((Int_t*)pBuffer);
  if (
      //TODO: these checks are rather braindead as we need to support both
      //little and big endian sizes on receive - think of something better!
      ( firstWord==bufferSize-sizeof(Int_t) || ntohl(firstWord)==bufferSize-sizeof(Int_t) ) &&
      ( firstWord>=34 || ntohl(firstWord)>=34 ) /*thats the minimum size of a streamed TObject*/
      )
  {
    ZMQTMessage msg(const_cast<void*>(pBuffer), bufferSize);
    TClass* objclass=msg.GetClass();
    TObject* pObject=msg.ReadObject(objclass);
    if (pObject && objclass) {
      return pObject;
    } else {
    }
  } else {
  }
  return NULL;
}

//__________________________________________________________________________________________________
AliZMQhelpers::ZMQTMessage* AliZMQhelpers::ZMQTMessage::Stream(TObject* pSrc, Int_t compression, unsigned verbosity, bool enableSchema)
{
  /// Helper function to stream an object into an ZMQTMessage
  /// The returned instance must be cleaned by the caller
  ///
  if (!pSrc) return NULL;

  ZMQTMessage* pMsg=new ZMQTMessage(kMESS_OBJECT);
  if (!pMsg) {
    if (verbosity>0) printf("ZMQTMessage: problem allocating a message\n");
    return NULL;
  }

  pMsg->EnableSchemaEvolution(enableSchema);
  pMsg->SetCompressionLevel(compression);
  pMsg->WriteObject(pSrc);
  if (pMsg->Length()>=0) {
    pMsg->SetLength(); // sets the length to the first (reserved) word

    // does nothing if the level is 0
    pMsg->Compress();

    if (pMsg->CompBuffer()) {
      pMsg->SetLength(); // set once more to have the byte order
    } else {
      if (verbosity>0) printf("ZMQTMessage: no CompBuffer\n");
    }
  } else {
    if (verbosity>0) printf("ZMQTMessage: problem with message length\n");
  }
  return pMsg;
}

