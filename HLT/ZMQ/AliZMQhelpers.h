#ifndef __AliZMQhelpers__
#define __AliZMQhelpers__

// a helper library for using ZMQ with ROOT, focussed on multipart messaging
// blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
// some of it might be inspired by czmq.h (http://czmq.zeromq.org)

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

#include <string>
#include <map>
#include "TMessage.h"
struct zmq_msg_t;
class TVirtualStreamerInfo;

namespace AliZMQhelpers
{
struct DataTopic;

extern void* gZMQcontext; //a global ZMQ context

//convenience typedefs:
//define a map of strings
typedef std::map<std::string,std::string> stringMap;
typedef std::pair<zmq_msg_t*, zmq_msg_t*> aliZMQframe;
typedef std::vector<aliZMQframe> aliZMQmsg;
typedef std::vector<std::pair<std::string, std::string> > aliZMQmsgStr;
typedef std::vector<std::pair<std::string, std::string> > aliStringVec;
typedef std::vector<TVirtualStreamerInfo*> aliZMQrootStreamerInfo;

//  Init and bind/connect a ZMQ socket using a string:
//  PUB@tcp://*:123123
//  SUB>tcp://localhost:123123,@tcp://*:454545
//  timeout is in ms, -1 is wait forever
int alizmq_socket_init(void*& socket, void* context, std::string config, int timeout=-1, int highWaterMark=10);
int alizmq_socket_close(void*& socket, int linger=0);
int alizmq_socket_state(void* socket);

//get the global context
void* alizmq_context();

// extract the socket mode from a config string
int alizmq_socket_type(std::string config);
int alizmq_socket_type(void* socket);
const char* alizmq_socket_name(int socketType);

//  --------------------------------------------------------------------------
//  Attach a socket to zero or more endpoints. If endpoints is not null,
//  parses as list of ZeroMQ endpoints, separated by commas, and prefixed by
//  '@' (to bind the socket) or '>' (to attach the socket - alternative: '-'). 
//  Returns 0 if all endpoints were valid, or -1 if there was a syntax error. 
//  If the endpoint does not start with '@' or '>'('-'), the serverish
//  argument defines whether it is used to bind (serverish = true)
//  or connect (serverish = false).
int alizmq_attach (void *self, const char *endpoints, bool serverish=false);
int alizmq_detach (void *self, const char *endpoints, bool serverish=false);

//general multipart messages (aliZMQmsg)
//to access, just iterate over it.
int alizmq_msg_recv(aliZMQmsg* message, void* socket, int flags);
int alizmq_msg_add(aliZMQmsg* message, DataTopic* topic, TObject* object, int compression=0, aliZMQrootStreamerInfo* streamers=NULL);
int alizmq_msg_add(aliZMQmsg* message, const DataTopic* topic, const std::string& data);
int alizmq_msg_add(aliZMQmsg* message, const DataTopic* topic, void* buffer, int size);
int alizmq_msg_add(aliZMQmsg* message, const std::string& topic, const std::string& data);
int alizmq_msg_copy(aliZMQmsg* dst, aliZMQmsg* src);
int alizmq_msg_send(aliZMQmsg* message, void* socket, int flags);
int alizmq_msg_close(aliZMQmsg* message);

//ROOT streamers utilities
//add streamers to a message
int alizmq_msg_prepend_streamer_infos(aliZMQmsg* message, aliZMQrootStreamerInfo* streamers);
//initialize ROOT internals with the incoming streamers so the objects can be decoded
int alizmq_msg_iter_init_streamer_infos(aliZMQmsg::iterator it);
//add new and unique streamers to the list based on the output of the ROOT serializer (newStreamers)
void alizmq_update_streamerlist(aliZMQrootStreamerInfo* streamers, const TCollection* newStreamers);
//this one is slow, use only for init kind of stuff
void alizmq_update_streamerlist_from_object(aliZMQrootStreamerInfo* streamers, TObject* object);

//checking identity of the frame via iterator
int alizmq_msg_iter_check(aliZMQmsg::iterator it, const DataTopic& topic);
int alizmq_msg_iter_check_id(aliZMQmsg::iterator it, const DataTopic& topic);
int alizmq_msg_iter_check_id(aliZMQmsg::iterator it, const std::string& topic);
//helpers for accessing data via iterators
int alizmq_msg_iter_topic(aliZMQmsg::iterator it, std::string& topic);
int alizmq_msg_iter_data(aliZMQmsg::iterator it, std::string& data);
int alizmq_msg_iter_topic(aliZMQmsg::iterator it, DataTopic& topic);
int alizmq_msg_iter_data(aliZMQmsg::iterator it, TObject*& object);
int alizmq_msg_iter_data(aliZMQmsg::iterator it, void*& buffer, size_t& size);

//string messages, no need to close, strings are copied
int alizmq_msg_send(std::string topic, std::string data, void* socket, int flags);
int alizmq_msg_recv(aliZMQmsgStr* message, void* socket, int flags);

//send a single block (one header + payload), ZMQ_SNDMORE should not be used
int alizmq_msg_send(DataTopic& topic, TObject* object, void* socket, int flags, int compression=0, aliZMQrootStreamerInfo* streamers=NULL);
int alizmq_msg_send(const DataTopic& topic, const std::string& data, void* socket, int flags);

//deallocate an object - callback for ZMQ
void alizmq_deleteTObject(void*, void* object);
void alizmq_deleteTopic(void*, void* object);

const int kDataTypefIDsize = 8;
const int kDataTypefOriginSize = 4;
const int kDataTypeTopicSize = kDataTypefIDsize+kDataTypefOriginSize;

//Helper functions
bool Topicncmp(const char* topic, const char* reference, int topicSize=kDataTypeTopicSize, int referenceSize=kDataTypeTopicSize);
UInt_t CharArr2uint32(const char* str);
ULong64_t CharArr2uint64(const char* str);

//helper function to print a hex/ASCII dump of some memory
void hexDump (const char* desc, const void* addr, int len);
void hexDump (aliZMQmsg* message, size_t maxsize=16);

struct BaseDataTopic
{
  static const UInt_t fgkMagicNumber;
  UInt_t fMagicNumber;  // 4 bytes
  UInt_t fHeaderSize;   // 4 bytes
  UInt_t fFlags;        // 4 bytes
  UInt_t fBaseHeaderVersion;  // 4 bytes
  ULong64_t fHeaderDescription; // 8 bytes
  ULong64_t fHeaderSerialization; // 8 bytes
  BaseDataTopic();
  BaseDataTopic(UInt_t size, ULong64_t desc, ULong64_t seri);
  static BaseDataTopic* Get(void* buf) {
    return (*reinterpret_cast<UInt_t*>(buf)==fgkMagicNumber)?
           reinterpret_cast<BaseDataTopic*>(buf):
           NULL;
  }
};

//the data header, describes the data frame
struct DataTopic : public BaseDataTopic
{
  static const ULong64_t fgkDataTopicDescription;
  static const UInt_t fgkTopicSerialization;
  ULong64_t fDataDescription[2]; // 2*8 bytes
  UInt_t fDataOrigin; // 4 bytes
  UInt_t fReserved;   // 4 bytes
  ULong64_t fDataSerialization; // 8 bytes
  ULong64_t fSpecification;     // 8 bytes data specification of the data block
  ULong64_t fPayloadSize;       // 8 bytes

  //ctor
  DataTopic()
    : BaseDataTopic(sizeof(DataTopic), fgkDataTopicDescription, fgkTopicSerialization)
    , fDataDescription()
    , fDataOrigin(0)
    , fReserved(0)
    , fDataSerialization(0)
    , fSpecification(0)
    , fPayloadSize(0)
  {
    fDataDescription[0]=0;
    fDataDescription[1]=0;
  }

  //ctor
  DataTopic(const char* id, const char* origin, int spec )
    : BaseDataTopic(sizeof(DataTopic), fgkDataTopicDescription, fgkTopicSerialization)
    , fDataDescription()
    , fDataOrigin(0)
    , fReserved(0)
    , fDataSerialization(0)
    , fSpecification(spec)
    , fPayloadSize(0)
  {
    fDataDescription[0] = 0;
    fDataDescription[1] = CharArr2uint64(id);
    fDataOrigin = CharArr2uint32(origin);
  }

  bool operator==( const DataTopic& dt )
  {
    bool topicMatch = Topicncmp(dt.GetIDstr(),GetIDstr());
    return topicMatch;
  }

  std::string Description() const
  {
    std::string description(GetIDstr(),
                            sizeof(fDataDescription[1])+sizeof(fDataOrigin));
    description+=" spec:";
    char numstr[21];
    snprintf(numstr, 21, "%llx", fSpecification);
    description+=numstr;
    return description;
  }

  inline std::string GetOrigin() const {
    std::string origin(GetOriginStr(), sizeof(fDataOrigin));
    return origin;
  }
  inline std::string GetID() const {
    std::string id(GetIDstr(), sizeof(fDataDescription[1]));
    return id;
  }
  UInt_t GetSpecification() const {return fSpecification;}
  inline const ULong64_t* GetIDptr() const {return &fDataDescription[1];}
  inline const char* GetIDstr() const {return reinterpret_cast<const char*>(GetIDptr());}
  inline const UInt_t* GetOriginPtr() const {return &fDataOrigin;}
  inline const char* GetOriginStr() const {return reinterpret_cast<const char*>(GetOriginPtr());}
  inline void SetID(ULong64_t id) {fDataDescription[1]=id;}
  inline void SetOrigin(UInt_t origin) {fDataOrigin = origin;}
  inline void SetID(const char* s) {fDataDescription[1]=*reinterpret_cast<const ULong64_t*>(s);}
  inline void SetOrigin(const char* s) {fDataOrigin = *reinterpret_cast<const UInt_t*>(s);}
  inline void SetSpecification(UInt_t spec) {fSpecification=spec;}
  inline void SetSerialization(ULong64_t s) {fDataSerialization=s;}
  static DataTopic* Get(void* buf) {
    BaseDataTopic* bdt = BaseDataTopic::Get(buf);
    return (bdt && bdt->fHeaderDescription==fgkDataTopicDescription)?
            reinterpret_cast<DataTopic*>(buf):NULL;
  }
};

//common data type definitions, compatible with AliHLTDataTypes v25
extern const DataTopic kDataTypeStreamerInfos;
extern const DataTopic kDataTypeInfo;
extern const DataTopic kDataTypeConfig;
extern const DataTopic kDataTypeTObject;
extern const DataTopic kDataTypeTH1;

extern const ULong64_t kSerializationROOT;
extern const ULong64_t kSerializationNONE;

//a general utility to tokenize strings
std::vector<std::string> TokenizeString(const std::string input, const std::string delimiters);
//parse
stringMap ParseParamString(const std::string paramString);
std::string GetParamString(const std::string param, const std::string paramstring);

//load ROOT libraries specified in comma separated string
int LoadROOTlibs(std::string libstring, bool verbose=false);

//______________________________________________________________________________
inline ULong64_t CharArr2uint64(const char* str)
{
	return((ULong64_t) str[0] |
         (str[0] ? ((ULong64_t) str[1] << 8 |
         (str[1] ? ((ULong64_t) str[2] << 16 |
         (str[2] ? ((ULong64_t) str[3] << 24 |
         (str[3] ? ((ULong64_t) str[4] << 32 |
         (str[4] ? ((ULong64_t) str[5] << 40 |
         (str[5] ? ((ULong64_t) str[6] << 48 |
         (str[6] ? ((ULong64_t) str[7] << 56 )
          : 0)) : 0)) : 0)) : 0)) : 0)) : 0)) : 0));
}

//______________________________________________________________________________
inline UInt_t CharArr2uint32(const char* str)
{
	return((UInt_t) str[0] |
         (str[0] ? ((UInt_t) str[1] << 8 |
         (str[1] ? ((UInt_t) str[2] << 16 |
         (str[2] ? ((UInt_t) str[3] << 24)
          : 0)) : 0)) : 0));
}

//______________________________________________________________________________
class ZMQTMessage : public TMessage {
  public:
    ZMQTMessage(UInt_t what = kMESS_ANY, Int_t bufsiz = TBuffer::kInitialSize) : TMessage(what,bufsiz) {}
    ZMQTMessage(void* buf, Int_t len) : TMessage(buf, len) {ResetBit(kIsOwner);}
    static TObject* Extract(const void* pBuffer, unsigned bufferSize, unsigned verbosity=0);
    static ZMQTMessage* Stream(TObject* pSrc, Int_t compression, unsigned verbosity=0, bool enableSchema=kFALSE);
#ifdef AliZMQhelpers_AliHLTMessageFormat
    void SetLength() const
    {
      // Set the message length at the beginning of the message buffer a la AliHLTMessage
      // using native byte order (little endian on x86)
      if (IsWriting()) {
        char *buf = Buffer();
        *((UInt_t*)buf) = (UInt_t)(Length() - sizeof(UInt_t));

        if (CompBuffer()) {
          buf = CompBuffer();
          *((UInt_t*)buf) = (UInt_t)(CompLength() - sizeof(UInt_t));
        }
      }
    }
#endif
};

}  //end namespace AliZMQhelpers

#endif

