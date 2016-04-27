// blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//see header file for details
//

#include "AliZMQhelpers.h"

#include "zmq.h"
#include <cstring>
#include <cassert>
#include <unistd.h>
#include <vector>

#include "AliHLTDataTypes.h"
#include "TString.h"
#include "TObjString.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "AliHLTMessage.h"
#include "TStreamerInfo.h"
#include "TClass.h"

//init the shared context to null
void* AliZMQhelpers::gZMQcontext = NULL;

//_______________________________________________________________________________________
void* alizmq_context()
{
  if (!AliZMQhelpers::gZMQcontext) AliZMQhelpers::gZMQcontext=zmq_ctx_new();
  return AliZMQhelpers::gZMQcontext;
}

//_______________________________________________________________________________________
int alizmq_detach (void *self, const char *endpoints, bool serverish)
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
int alizmq_attach (void *self, const char *endpoints, bool serverish)
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
int alizmq_socket_state(void* socket)
{
  int events=0;
  size_t len = sizeof(events);
  zmq_getsockopt(socket, ZMQ_EVENTS, &events, &len);
  return events;
}

//_______________________________________________________________________________________
int alizmq_socket_type(void* socket)
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
int alizmq_socket_type(std::string config)
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
const char* alizmq_socket_name(int socketType)
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
int alizmq_socket_close(void*& socket, int linger)
{
  zmq_setsockopt(socket, ZMQ_LINGER, &linger, sizeof(linger));
  int rc = zmq_close(socket);
  if (rc>=0) socket = NULL;
  return rc;
}

//_______________________________________________________________________________________
int alizmq_socket_init(void*& socket, void* context, std::string config, int timeout, int highWaterMark)
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
int alizmq_msg_add(aliZMQmsg* message, const std::string& topic, const std::string& data)
{
  //add a frame to the mesage
  int rc = 0;
  
  //prepare topic msg
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, topic.size());
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg),topic.data(),topic.size());

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
int alizmq_msg_add(aliZMQmsg* message, const AliHLTDataTopic* topic, void* data, int size)
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
int alizmq_msg_add(aliZMQmsg* message, const AliHLTDataTopic* topic, const std::string& data)
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
int alizmq_msg_add(aliZMQmsg* message, const AliHLTDataTopic* topic, TObject* object,
                   int compression, aliZMQrootStreamerInfo* streamers)
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
  memcpy(zmq_msg_data(topicMsg), topic, sizeof(*topic));

  //prepare data msg
  AliHLTMessage* tmessage = AliHLTMessage::Stream(object, compression, 0, streamers);
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
int alizmq_msg_send(aliZMQmsg* message, void* socket, int flagsUser)
{
  int nBytes=0;
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
int alizmq_msg_send(std::string topic, std::string data, void* socket, int flags)
{
  int rc = 0;
  zmq_msg_t topicMsg;
  zmq_msg_init_size(&topicMsg, topic.size());
  memcpy(zmq_msg_data(&topicMsg), topic.data(), zmq_msg_size(&topicMsg));
  rc = zmq_msg_send(&topicMsg, socket, ZMQ_SNDMORE);
  if (rc<0) 
  {
    //printf("unable to send topic: %s\n", topic.c_str());
    zmq_msg_close(&topicMsg);
    return rc;
  }

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
  return rc;
}

//_______________________________________________________________________________________
int alizmq_msg_prepend_streamer_infos(aliZMQmsg* message, aliZMQrootStreamerInfo* streamers)
{
  //prepend the streamer info to the message as first block.
  int rc = 0;

  AliHLTDataTopic topic = kAliHLTDataTypeStreamerInfo;
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
  AliHLTMessage* tmessage = AliHLTMessage::Stream(&listOfInfos, 1); //compress
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
void alizmq_update_streamerlist(aliZMQrootStreamerInfo* streamers, const TObjArray* newStreamers)
{
  //update the list of streamers used
  if (!streamers) return;
  if (!newStreamers) return;

  for (int i=0; i<newStreamers->GetEntriesFast(); i++) {
    TVirtualStreamerInfo* info = const_cast<TVirtualStreamerInfo*> (
      static_cast<const TVirtualStreamerInfo*>((*newStreamers)[i]) );
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
  }
}

//_______________________________________________________________________________________
int alizmq_msg_iter_init_streamer_infos(aliZMQmsg::iterator it)
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
int alizmq_msg_send(const AliHLTDataTopic& topic, TObject* object, void* socket, int flags, 
                    int compression, aliZMQrootStreamerInfo* streamers)
{
  int rc = 0;

  AliHLTMessage* tmessage = AliHLTMessage::Stream(object, compression);
  zmq_msg_t dataMsg;
  rc = zmq_msg_init_data( &dataMsg, tmessage->Buffer(), tmessage->Length(),
      alizmq_deleteTObject, tmessage);
  
  if (streamers) {
    alizmq_update_streamerlist(streamers, tmessage->GetStreamerInfos());
  }

  //then send the object topic
  rc = zmq_send( socket, &topic, sizeof(topic), ZMQ_SNDMORE );
  if (rc<0) 
  {
    zmq_msg_close(&dataMsg);
    //printf("unable to send topic: %s %s\n", topic.Description().c_str(), zmq_strerror(errno));
    return rc;
  }

  //send the object itself
  rc = zmq_msg_send(&dataMsg, socket, flags);
  if (rc<0) 
  {
    //printf("unable to send data: %s %s\n", tmessage->GetName(), zmq_strerror(errno));
    zmq_msg_close(&dataMsg);
    return rc;
  }
  return rc;
}

//______________________________________________________________________________
int alizmq_msg_send(const AliHLTDataTopic& topic, const std::string& data, void* socket, int flags)
{
  int rc = 0;
  rc = zmq_send( socket, &topic, sizeof(topic), ZMQ_SNDMORE );
  if (rc<0) 
  {
    //printf("unable to send topic: %s %s\n", topic.Description().c_str(), zmq_strerror(errno));
    return rc;
  }

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
  return rc;
}

//______________________________________________________________________________
void alizmq_deleteTObject(void*, void* object)
{
  //delete the TBuffer, for use in zmq_msg_init_data(...) only.
  //printf("deleteObject called! ZMQ just sent and destroyed the message!\n");
  TObject* tobject = static_cast<TObject*>(object);
  delete tobject;
}

//______________________________________________________________________________
void alizmq_deleteTopic(void*, void* object)
{
  //delete the TBuffer, for use in zmq_msg_init_data(...) only.
  //printf("deleteObject called! ZMQ just sent and destroyed the message!\n");
  AliHLTDataTopic* topic = static_cast<AliHLTDataTopic*>(object);
  delete topic;
}


//_______________________________________________________________________________________
int alizmq_msg_close(aliZMQmsg* message)
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
int alizmq_msg_iter_check(aliZMQmsg::iterator it, const AliHLTDataTopic& topic)
{
  AliHLTDataTopic actualTopic;
  alizmq_msg_iter_topic(it, actualTopic);
  if (actualTopic == topic) return 0;
  return 1;
}

//_______________________________________________________________________________________
int alizmq_msg_iter_check(aliZMQmsg::iterator it, const std::string& topic)
{
  std::string actualTopic;
  alizmq_msg_iter_topic(it, actualTopic);
  return actualTopic.compare(0,topic.size(),topic);
}

//_______________________________________________________________________________________
int alizmq_msg_iter_topic(aliZMQmsg::iterator it, std::string& topic)
{
  zmq_msg_t* message = it->first;
  topic.assign((char*)zmq_msg_data(message),zmq_msg_size(message));
  return 0;
}

//_______________________________________________________________________________________
int alizmq_msg_iter_data(aliZMQmsg::iterator it, std::string& data)
{
  zmq_msg_t* message = it->second;
  data.assign((char*)zmq_msg_data(message),zmq_msg_size(message));
  return 0;
}

//_______________________________________________________________________________________
int alizmq_msg_iter_topic(aliZMQmsg::iterator it, AliHLTDataTopic& topic)
{
  zmq_msg_t* message = it->first;
  memcpy(&topic, zmq_msg_data(message),std::min(zmq_msg_size(message),sizeof(topic)));
  return 0;
}

//_______________________________________________________________________________________
int alizmq_msg_iter_data(aliZMQmsg::iterator it, TObject*& object)
{
  zmq_msg_t* message = it->second;
  size_t size = zmq_msg_size(message);
  void* data = zmq_msg_data(message);

  object = AliHLTMessage::Extract(data, size);
  return 0;  
}

//_______________________________________________________________________________________
int alizmq_msg_copy(aliZMQmsg* dst, aliZMQmsg* src)
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
int alizmq_msg_recv(aliZMQmsg* message, void* socket, int flags)
{
  int rc = -1;
  int receiveStatus=0;
  while (true)
  {
    zmq_msg_t* topicMsg = new zmq_msg_t;
    rc = zmq_msg_init(topicMsg);
    rc = zmq_msg_recv(topicMsg, socket, flags);
    if (!zmq_msg_more(topicMsg) || rc<0)
    {
      zmq_msg_close(topicMsg);
      delete topicMsg;
      receiveStatus=-1;
      break;
    }
    receiveStatus+=rc;

    zmq_msg_t* dataMsg = new zmq_msg_t;
    rc = zmq_msg_init(dataMsg);
    rc = zmq_msg_recv(dataMsg, socket, flags);
    if (rc<0)
    {
      zmq_msg_close(topicMsg);
      zmq_msg_close(dataMsg);
      delete topicMsg;
      delete dataMsg;
      receiveStatus=-1;
      break;
    }
    receiveStatus+=rc;

    message->push_back(std::make_pair(topicMsg,dataMsg));

    int more=0;
    size_t moreLength = sizeof(more);
    rc = zmq_getsockopt(socket, ZMQ_RCVMORE, &more, &moreLength);
    if (!more) break;
  }
  return receiveStatus;
}

//_______________________________________________________________________________________
int alizmq_msg_recv(aliZMQmsgStr* message, void* socket, int flags)
{
  int rc = -1;
  int receiveStatus=0;
  while (true)
  {
    zmq_msg_t topicMsg;
    rc = zmq_msg_init(&topicMsg);
    rc = zmq_msg_recv(&topicMsg, socket, flags);
    if (!zmq_msg_more(&topicMsg) || receiveStatus<0)
    {
      zmq_msg_close(&topicMsg);
      receiveStatus=-1;
      break;
    }
    receiveStatus+=rc;

    zmq_msg_t dataMsg;
    rc = zmq_msg_init(&dataMsg);
    rc = zmq_msg_recv(&dataMsg, socket, flags);
    if (receiveStatus<0)
    {
      zmq_msg_close(&topicMsg);
      zmq_msg_close(&dataMsg);
      receiveStatus=-1;
      break;
    }
    receiveStatus+=rc;

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
  return receiveStatus;
}

//_______________________________________________________________________________________
TString AliOptionParser::GetFullArgString(int argc, char** argv)
{
  TString argString;
  TString argument="";
  if (argc>0) {
    for (int i=1; i<argc; i++) {
      argument=argv[i];
      if (argument.IsNull()) continue;
      if (!argString.IsNull()) argString+=" ";
      argString+=argument;
    }
  }
  return argString;
}

//______________________________________________________________________________
int AliOptionParser::ProcessOptionString(TString arguments)
{
  //process passed options, return number of processed valid options
  aliStringVec* options = TokenizeOptionString(arguments);
  int nOptions=0;
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //printf("  %s : %s\n", i->first.data(), i->second.data());
    if (ProcessOption(i->first,i->second)<0)
    {
      nOptions=-1;
      break;
    }
    nOptions++;
  }
  delete options; //tidy up

  return nOptions;
}

//______________________________________________________________________________
aliStringVec* AliOptionParser::TokenizeOptionString(const TString strIn)
{
  //options have the form:
  // -option value
  // -option=value
  // -option
  // --option value
  // --option=value
  // --option
  // option=value
  // option value
  // (value can also be a string like 'some string')
  //
  // options can be separated by ' ' arbitrarily combined, e.g:
  //"-option option1=value1 --option2 value2, -option4=\'some string\'"

  //optionRE by construction contains a pure option name as 3rd submatch (without --,-, =)
  //valueRE does NOT match options
  TPRegexp optionRE("(?:(-{1,2})|((?='?[^=]+=?)))"
                    "((?(2)(?:(?(?=')'(?:[^'\\\\]++|\\.)*+'|[^ =]+))(?==?))"
                    "(?(1)[^ =]+(?=[= $])))");
  TPRegexp valueRE("(?(?!(-{1,2}|[^ =]+=))"
                   "(?(?=')'(?:[^'\\\\]++|\\.)*+'"
                   "|[^ =]+))");

  aliStringVec* options = new aliStringVec;

  //first split in lines (by newline) and ignore comments
  TObjArray* lines = strIn.Tokenize("\n\r");
  TIter nextLine(lines);
  while (TObjString* objString = (TObjString*)nextLine())
  {
  TString line = objString->String();
  if (line.BeginsWith("#")) continue;
  if (line.BeginsWith("//")) continue;
  TArrayI pos;
  const TString mods="";
  Int_t start = 0;
  while (1) {
    Int_t prevStart=start;
    TString optionStr="";
    TString valueStr="";

    //check if we have a new option in this field
    Int_t nOption=optionRE.Match(line,mods,start,10,&pos);
    if (nOption>0)
    {
      optionStr = line(pos[6],pos[7]-pos[6]);
      optionStr=optionStr.Strip(TString::kBoth,'\n');
      optionStr=optionStr.Strip(TString::kBoth,'\'');
      optionStr=optionStr.Strip(TString::kLeading,'-');
      start=pos[1]; //update the current character to the end of match
    }

    //check if the next field is a value
    Int_t nValue=valueRE.Match(line,mods,start,10,&pos);
    if (nValue>0)
    {
      valueStr = line(pos[0],pos[1]-pos[0]);
      valueStr=valueStr.Strip(TString::kBoth,'\n');
      valueStr=valueStr.Strip(TString::kBoth,'\'');
      start=pos[1]; //update the current character to the end of match
    }

    //skip empty entries
    if (nOption>0 || nValue>0)
    {
      options->push_back(std::make_pair(optionStr.Data(),valueStr.Data()));
    }

    if (start>=line.Length()-1 || start==prevStart ) break;
  }

  }//while(nextLine())
  lines->Delete();
  delete lines;

  return options;
}

//tokenize a std::string
using namespace std;
vector<string> TokenizeString(const string input, const string delimiters)
{
  vector<string> output;
  output.reserve(10);
  size_t start = 0;
  size_t end = input.find_first_of(delimiters);
  do
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

//tokenize a string delimited by semi-colons and return a map
//of key value pairs, "KEY=VALUE;key2=value2"
stringMap ParseParamString(const string paramString)
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

