
// blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
// some of it might be inspired by czmq.h

#include <string>
#include <map>
#include "TString.h"
struct zmq_msg_t;
struct AliHLTDataTopic;

//convenience typedefs:
//define a map of strings
typedef std::map<std::string,std::string> stringMap;
typedef std::vector<std::pair<zmq_msg_t*, zmq_msg_t*> > aliZMQmsg;
typedef std::vector<std::pair<std::string, std::string> > aliZMQmsgStr;

//  Init and bind/connect a ZMQ socket using a string:
//  PUB@tcp://*:123123
//  SUB>tcp://localhost:123123,@tcp://*:454545
//  timeout is in ms, -1 is wait forever
int alizmq_socket_init(void*& socket, void* context, std::string config, int timeout=-1);

// extract the socket mode from a config string
int alizmq_socket_type(std::string config);
int alizmq_socket_type(void* socket);

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
int alizmq_msg_close(aliZMQmsg* message);
//helpers for accessing data via iterators
int alizmq_msg_iter_topic(aliZMQmsg::iterator it, AliHLTDataTopic& topic);
int alizmq_msg_iter_data(aliZMQmsg::iterator it, TObject*& object);

//string messages, no need to close, strings are copied
int alizmq_msg_send(std::string topic, std::string data, void* socket, int flags);
int alizmq_msg_recv(aliZMQmsgStr* message, void* socket, int flags);

//send a single frame
int alizmq_msg_send(const AliHLTDataTopic& topic, TObject* object, void* socket, int flags, int compression=0);

//deallocate an object - callback for ZMQ
void alizmq_deleteTObject(void*, void* object);

//simple zmq multi part message class
//behaves like a map.
//this is to simplify receiving/sending multipart msgs
//and to handle message destruction automatically
//keep it simple!
//class AliZMQmsg {
//public:
//  AliZMQmsg() {}
//  ~AliZMQmsg() {}
//  int Receive(void* socket) {return 0;}
//  int Send(void* socket) {return 0;}
//  void Add(zmq_msg_t* topic, zmq_msg_t* data) {}
//
//  //define (delegate) iterators
//  typedef aliZMQmsg::iterator iterator;
//  typedef aliZMQmsg::const_iterator const_iterator;
//  iterator begin() { return fMessage.begin(); }
//  iterator end() { return fMessage.end(); }
//private:
//  aliZMQmsg fMessage;
//};

//simple option parser class
class AliOptionParser {
public:
  //call this to parse the args
  int ProcessOptionString(TString arguments);
  int ProcessOptionString(int argc, char** argv);
  //implement this to process one option at a time
  int ProcessOption(TString /*option*/, TString /*value*/) {return 0;}
  //convert argc/argv into a TString of options
  TString GetFullArgString(int argc, char** argv);
private:
  stringMap* TokenizeOptionString(const TString str);
};

