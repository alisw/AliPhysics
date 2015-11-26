#include "zmq.h"
#include "AliHLTDirectHOMERManager.h"
#include "AliHLTHOMERBlockDesc.h"
#include <iostream>
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"
#include <time.h>

int sendWithTopic(char bufferTopic[kAliHLTComponentDataTypeTopicSize], void* bufferPtr, size_t bufferSize, void* mainSocket, void* monitorSocket, int flags=0, bool fOnlyPublish=false);

int main(int argc, char** argv)
{
  int mainReturnCode=0;
  //configuration vars
  Bool_t fSendMultipart = kTRUE;
  Bool_t fOnlyPublish = kFALSE;
  Bool_t fVerbose = kTRUE;
  const char* fZMQendpointREP = "tcp://*:60200";
  const char* fZMQendpointPUB = "tcp://*:60201";
  
  int rc = 0;
  //init the ZMQ sockets
  //the req-rep socket - triggers readout
  //control PUB socket, publishes the same message
  void* fZMQcontext = NULL;
  void* fZMQrep = NULL;
  void* fZMQpub = NULL;
  fZMQcontext = zmq_ctx_new();
  fZMQrep = zmq_socket(fZMQcontext, ZMQ_REP);
  fZMQpub = zmq_socket(fZMQcontext, ZMQ_PUB);
  //set socket options
  int lingerValue = 10;
  rc += zmq_setsockopt(fZMQrep, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
  rc += zmq_setsockopt(fZMQpub, ZMQ_LINGER, &lingerValue, sizeof(lingerValue));
  int highWaterMarkSend = 20;
  rc += zmq_setsockopt(fZMQrep, ZMQ_SNDHWM, &highWaterMarkSend, sizeof(highWaterMarkSend));
  rc += zmq_setsockopt(fZMQpub, ZMQ_SNDHWM, &highWaterMarkSend, sizeof(highWaterMarkSend));
  int highWaterMarkRecv = 20;
  rc += zmq_setsockopt(fZMQrep, ZMQ_RCVHWM, &highWaterMarkRecv, sizeof(highWaterMarkRecv));
  rc += zmq_setsockopt(fZMQpub, ZMQ_RCVHWM, &highWaterMarkRecv, sizeof(highWaterMarkRecv));
  if (rc!=0) {printf("cannot set socket options\n"); return 1;}
  //bind the sockets
  rc = 0;
  rc += zmq_bind(fZMQrep,fZMQendpointREP);
  rc += zmq_bind(fZMQpub,fZMQendpointPUB);
  if (rc!=0) {printf("cannot bind to %s %s\n",fZMQendpointREP,fZMQendpointPUB); return 1;}

  //the HOMER connection
  AliHLTDirectHOMERManager man(argc, argv);
  
  //default void topic
  char voidTopic[kAliHLTComponentDataTypeTopicSize];
  AliHLTComponent::DataType2Topic(kAliHLTVoidDataType, voidTopic);
  
  char printable[100];

  //main loop
  while(1)
  {
    errno=0;
    //create a default request: selection of any data:
    int requestTopicSize = -1; 
    char requestTopic[kAliHLTComponentDataTypeTopicSize];
    memset(requestTopic, '*', kAliHLTComponentDataTypeTopicSize);
    int requestSize = -1;
    char request[kAliHLTComponentDataTypeTopicSize];
    memset(request, '*', kAliHLTComponentDataTypeTopicSize);
    
    //wait for the request
    //if the request is just the topic - send an empty reply, and publish data on
    //pub socket; send data back if request is 2 parts - topic+body
    int64_t more = 0;
    size_t more_size = sizeof(more);
    do
    {
      requestTopicSize = zmq_recv(fZMQrep, requestTopic, kAliHLTComponentDataTypeTopicSize, 0);
      
      //logging
      char printable[kAliHLTComponentDataTypeTopicSize+1]; memset(printable,0,kAliHLTComponentDataTypeTopicSize+1);
      memcpy(printable, requestTopic, kAliHLTComponentDataTypeTopicSize);
      char timebuf[100];
      time_t now = time (0);
      strftime (timebuf, 100, "%Y-%m-%d %H:%M:%S", localtime (&now));
      //printf("\n%s\nrequest: %s size: %i, errno: %s\n", timebuf, printable, requestTopicSize, (requestTopicSize<0)?zmq_strerror(errno):"");
      printf("\n%s\nrequest: %s\n", timebuf, printable);
      
      rc = zmq_getsockopt(fZMQrep, ZMQ_RCVMORE, &more, &more_size);
      if (more==1)  {
        requestSize = zmq_recv(fZMQrep, request, kAliHLTComponentDataTypeTopicSize, 0);
        //printf("  size: %i, rc: %i, errno: %s\n", (int)requestSize, rc, (requestSize<0)?zmq_strerror(errno):"");
        rc = zmq_getsockopt(fZMQrep, ZMQ_RCVMORE, &more, &more_size);
      }
    } while (more==1);
    fOnlyPublish=false;
    if (requestSize <= 0) fOnlyPublish=true;

    //process HOMER blocks, and send back(out) matching ones
    TList* blockList = man.GetNextBlocks();
    TIter next(blockList);
    AliHLTHOMERBlockDesc* blockDesc = NULL;
    while ((blockDesc = (AliHLTHOMERBlockDesc*)next())) 
    {
      char blockTopic[kAliHLTComponentDataTypeTopicSize];
      blockDesc->Topic(blockTopic);
      
      bool isSelected = Topicncmp(requestTopic,blockTopic,requestTopicSize);

      if (fVerbose)
      {
        memset(printable,0,100);
        strncat(printable, blockTopic, kAliHLTComponentDataTypeTopicSize);
        printf("blockTopic: %s %s\n",printable,(isSelected)?"   <---- selected":"");
      }
      
      if (!isSelected) continue;

      //send the lot
      void* bufferPtr = NULL;
      size_t bufferSize = 0;
      if (blockDesc->GetData()) 
      {
        bufferPtr = blockDesc->GetData();
        bufferSize = blockDesc->GetSize();
      }
      sendWithTopic(blockTopic, bufferPtr, bufferSize, fZMQrep, fZMQpub, (fSendMultipart)?ZMQ_SNDMORE:0, fOnlyPublish);
    }
    //after sending all messages send an empty message to terminate
    //send reply
    //sendWithTopic(voidTopic, NULL, 0, fZMQrep, fZMQpub);
    sendWithTopic(voidTopic, NULL, 0, fZMQrep, fZMQpub, 0, fOnlyPublish);
  }

  //destroy ZMQ sockets
  zmq_close(fZMQrep);
  zmq_close(fZMQpub);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}

int sendWithTopic(char bufferTopic[kAliHLTComponentDataTypeTopicSize], void* bufferPtr, size_t bufferSize, void* mainSocket, void* monitorSocket, int flags, bool fOnlyPublish)
{
  //prepare messages to send
  //first frame:  topic
  //second frame: block
  int rc = 0;

  if (!bufferPtr) bufferSize=0;

  //prepare the messages
  //for the control socket, this does not copy the buffer
  zmq_msg_t topicMonitor;
  rc = zmq_msg_init_size(&topicMonitor, kAliHLTComponentDataTypeTopicSize);
  memcpy(zmq_msg_data(&topicMonitor), bufferTopic, kAliHLTComponentDataTypeTopicSize);
  zmq_msg_t messageMonitor;
  rc = zmq_msg_init_size(&messageMonitor, bufferSize);
  memcpy(zmq_msg_data(&messageMonitor), bufferPtr, bufferSize);

  //will be sending on two sockets, so prepare a copy
  zmq_msg_t topic;
  zmq_msg_init(&topic);
  rc = zmq_msg_copy(&topic, &topicMonitor);
  zmq_msg_t message;
  zmq_msg_init(&message);
  if (!fOnlyPublish) rc = zmq_msg_copy(&message, &messageMonitor);  
  
  //send reply
  rc=-1;
  if (mainSocket) rc = zmq_msg_send(&topic, mainSocket, ZMQ_SNDMORE);
  if (rc!=0) zmq_msg_close(&topic);
  rc=-1;
  if (mainSocket) rc = zmq_msg_send(&message, mainSocket, flags);
  if (rc!=0) zmq_msg_close(&message);
  
  //publish also on monitoring socket
  rc=-1;
  if (monitorSocket) rc = zmq_msg_send(&topicMonitor, monitorSocket, ZMQ_SNDMORE);
  if (rc!=0) zmq_msg_close(&topicMonitor);
  rc=-1;
  if (monitorSocket) rc = zmq_msg_send(&messageMonitor, monitorSocket, flags);
  if (rc!=0) zmq_msg_close(&messageMonitor);
 
  return 0;
}

