// a helper library for using ZMQ with ALIROOT, focussed on multipart messaging
// this lib implements the HLT specific interface, for general use cases
// see AliZMQhelpers.h
// blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
// some of it might be inspired by czmq.h (http://czmq.zeromq.org)
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

#include "AliHLTZMQhelpers.h"
#include "AliHLTMessage.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TKey.h"
#include "zmq.h"
#include "AliRawData.h"
#include "TSystem.h"
#include <exception>

using namespace AliZMQhelpers;

//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_iter_check_id(aliZMQmsg::iterator it, const AliHLTDataTopic& topic)
{
  AliHLTDataTopic actualTopic;
  alizmq_msg_iter_topic(it, actualTopic);
  if (actualTopic.GetID() == topic.GetID()) return 0;
  return 1;
}


//_______________________________________________________________________________________
int AliZMQhelpers::alizmq_msg_send(const AliHLTDataTopic& topic, TObject* object, void* socket, int flags, 
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
int AliZMQhelpers::alizmq_msg_send(const AliHLTDataTopic& topic, const std::string& data, void* socket, int flags)
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
int AliZMQhelpers::alizmq_msg_add(aliZMQmsg* message, DataTopic* topic, AliRawData* blob)
{
  //add a frame to the mesage
  int rc = 0;

  topic->SetSerialization(kSerializationNONE);

  //prepare topic msg
  zmq_msg_t* topicMsg = new zmq_msg_t;
  rc = zmq_msg_init_size( topicMsg, sizeof(*topic));
  if (rc<0) {
    zmq_msg_close(topicMsg);
    delete topicMsg;
    return -1;
  }
  memcpy(zmq_msg_data(topicMsg), topic, sizeof(*topic));

  zmq_msg_t* dataMsg = new zmq_msg_t;
  rc = zmq_msg_init_data( dataMsg, blob->GetBuffer(), blob->GetSize(),
       alizmq_deleteTObject, blob);
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

//______________________________________________________________________________
int AliZMQhelpers::alizmq_file_write(AtomicFile& afile, AliHLTDataTopic topic, TObject* object)
{
  //To avoid problems with resource leaks in the HLT the use of non owning TCollection
  //is forbidden.
  //Use AliHLTObjArray/AliHLTList instead.
  //cannot use AliFatal here, so we throw.
  TCollection* collection = dynamic_cast<TCollection*>(object);
  if (collection) {
    if (!collection->IsOwner() &&
        !(strcmp(collection->ClassName(),"AliHLTObjArray")==0 || strcmp(collection->ClassName(),"AliHLTList")==0)) {
      afile.Close();
      throw std::runtime_error("Cannot use a non-owning TCollection in HLT unless it is a AliHLTList or AliHLTObjArray");
    }
  }

  topic.SetSerialization(kSerializationROOT);

  TFile* file = afile.GetFile();
  if (!file) { printf("could not get file\n"); return -1; }
  //TODO: do some sanity checks
  AliRawData topicBlob;
  //const_cast is here because AliRawData interface is broken
  topicBlob.SetBuffer(const_cast<AliHLTDataTopic*>(&topic), sizeof(topic));
  file->WriteObject(&topicBlob, "header");
  file->WriteObject(object, "payload");
  return 0;
}

//______________________________________________________________________________
int AliZMQhelpers::alizmq_file_write(AtomicFile& afile, const AliHLTDataTopic& topic, const void* buf, Int_t len)
{
  TFile* file = afile.GetFile();
  if (!file) { printf("could not get file\n"); return -1; }
  //TODO: do some sanity checks
  AliRawData headerBlob;
  //const_cast is here because AliRawData interface is broken
  headerBlob.SetBuffer(const_cast<AliHLTDataTopic*>(&topic), sizeof(topic));
  file->WriteObject(&headerBlob, "header");
  AliRawData payloadBlob;
  //const_cast is here because AliRawData interface is broken
  payloadBlob.SetBuffer(const_cast<void*>(buf), len);
  file->WriteObject(&payloadBlob, "payload");
  return 0;
}

//______________________________________________________________________________
int AliZMQhelpers::alizmq_file_write(AtomicFile& afile, aliZMQmsg* message, bool deserializeROOTobjects)
{
  for (aliZMQmsg::iterator i=message->begin(); i!=message->end(); ++i)
  {
    AliHLTDataTopic topic;
    alizmq_msg_iter_topic(i,topic);
    if (topic.fDataSerialization==kSerializationROOT && deserializeROOTobjects) {
      TObject* object = NULL;
      alizmq_msg_iter_data(i, object);
      if (!object) {
        printf("alizmq_file_writ error: root object could not be deserialized from message\n");
        return -1;
      }
      alizmq_file_write(afile, topic, object);
      delete object;
    } else {
    void* payloadBuf = NULL;
    size_t payloadSize = 0;
    alizmq_msg_iter_data(i,payloadBuf, payloadSize);
    alizmq_file_write(afile, topic, payloadBuf, payloadSize);
    }
  }

  return 0;
}

//______________________________________________________________________________
int AliZMQhelpers::alizmq_file_read(TFile& file, aliZMQmsg* message)
{
  if (!message) {printf("no message!\n"); return -1;}
  if (file.GetNkeys()%2) {printf("number of keys in the file not even\n");return 1;}
  TList* listOfKeys = file.GetListOfKeys();
  TIter iter(listOfKeys, kIterBackward);
  while (TKey* key = static_cast<TKey*>(iter()))
  {
    std::string objectName = key->GetName();
    if (objectName!="header") continue;
    Short_t cycle = key->GetCycle();
    AliRawData* headerBlob = static_cast<AliRawData*>(key->ReadObj());
    if (!headerBlob) {printf("could not read the header blob for cycle %i\n",cycle); continue;}

    AliHLTDataTopic topic;
    memcpy(&topic,
           headerBlob->GetBuffer(),
           (size_t(headerBlob->GetSize())>sizeof(topic))?sizeof(topic):size_t(headerBlob->GetSize()));
    delete headerBlob;

    TKey* payloadKey = file.GetKey("payload",cycle);
    if (!payloadKey) {printf("no payload found for header %i\n",cycle); continue;}
    string payloadClassName = payloadKey->GetClassName();
    if (payloadClassName.find("AliRawData")!=string::npos) {
      AliRawData* payloadBlob = static_cast<AliRawData*>(payloadKey->ReadObj());
      if (!payloadBlob) {printf("could not read the payload blob for cycle %i\n",cycle); continue;}
      if (alizmq_msg_add(message, &topic, payloadBlob)<0) {delete payloadBlob;}
    } else {
      TObject* object = payloadKey->ReadObj();
      if (!object) {printf("could not extract the root object from key %s with cycle %i \n", key->GetName(), cycle); continue;}
      if (alizmq_msg_add(message, &topic, object)<0) {delete object;}
      delete object;
    }
  }
  return 0;
}

//______________________________________________________________________________
AliZMQhelpers::AtomicFile::AtomicFile(const char* name)
  : targetFileName(name)
  , tempFile(NULL) {
  TString tmpname(name);
  tmpname+=".";
  if (FILE* fd = gSystem->TempFileName(tmpname,gSystem->DirName(name))) {
    tempFile = new TFile(tmpname.Data(),"recreate");
    fclose(fd);
  }
}

//______________________________________________________________________________
AliZMQhelpers::AtomicFile::~AtomicFile()
{
  Close();
}

//______________________________________________________________________________
void AliZMQhelpers::AtomicFile::Close()
{
  if (!tempFile) return;
  TString tempname = tempFile->GetName();
  tempFile->Close(); delete tempFile; tempFile=NULL;
  gSystem->Rename(tempname.Data(),targetFileName.Data());
  targetFileName.Clear();
}
