#ifndef __AliHLTZMQhelpers__
#define __AliHLTZMQhelpers__

// a helper library for using ZMQ with ALIROOT, focussed on multipart messaging
// this lib implements the HLT specific interface, for general use cases
// see AliZMQhelpers.h
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

#include "AliZMQhelpers.h"
#include "AliHLTDataTypes.h"
#include <cstring>
#include "TNamed.h"
#include "TFile.h"
class AliRawData;

namespace AliZMQhelpers
{

//this is just to implement the methods which depend on aliroot HLT definitions
//for ZMQ communication
struct AliHLTDataTopic : public DataTopic
{
  //ctor
  AliHLTDataTopic()
    : DataTopic()
  {
  }

  //copy ctor
  AliHLTDataTopic(const AliHLTComponentDataType& dataType)
    : DataTopic()
  {
    SetID(dataType.fID);
    SetOrigin(dataType.fOrigin);
  }

  //copy ctor
  AliHLTDataTopic(const AliHLTComponentBlockData& blockData)
    : DataTopic()
  {
    SetSpecification(blockData.fSpecification);
    SetID(blockData.fDataType.fID);
    SetOrigin(blockData.fDataType.fOrigin);
    if (strncmp(blockData.fDataType.fID,"ROOT",4)==0) {
      SetSerialization(kSerializationROOT);
    }
  }

  //partial (no fSpecification) copy from AliHLTComponentDataType
  AliHLTDataTopic& operator=( const AliHLTComponentDataType& dataType )
  {
    SetID(dataType.fID);
    SetOrigin(dataType.fOrigin);
    return *this;
  }

  //assignment from a AliHLTComponentBlockData
  AliHLTDataTopic& operator=( const AliHLTComponentBlockData& blockData )
  {
    SetSpecification(blockData.fSpecification);
    SetID(blockData.fDataType.fID);
    SetOrigin(blockData.fDataType.fOrigin);
    if (strncmp(blockData.fDataType.fID,"ROOT",4)==0) {
      SetSerialization(kSerializationROOT);
    }
    return *this;
  }

  bool operator==( const AliHLTDataTopic& dt )
  {
    bool topicMatch =  Topicncmp(dt.GetIDstr(),GetIDstr());
    return topicMatch;
  }

  bool operator==( const AliHLTComponentDataType& dataType)
  {
    AliHLTComponentDataType dt;
    Fill(dt);
    return dt==dataType;
  }

  void Fill(AliHLTComponentDataType& dt)
  {
    memcpy( dt.fID, &fDataDescription[1], kAliHLTComponentDataTypefIDsize );
    memcpy( dt.fOrigin, &fDataOrigin, kAliHLTComponentDataTypefOriginSize );
  }

};

class AtomicFile {
  TString targetFileName;
  TFile* tempFile;
  public:
  AtomicFile(const char* name);
  ~AtomicFile();
  TFile* GetFile() {return tempFile;}
  void Close();
};

int alizmq_msg_iter_check_id(aliZMQmsg::iterator it, const AliHLTDataTopic& topic);
int alizmq_msg_send(const AliHLTDataTopic& topic, TObject* object, void* socket, int flags,
                    int compression=0, aliZMQrootStreamerInfo* streamers=NULL);
int alizmq_msg_send(const AliHLTDataTopic& topic, const std::string& data, void* socket, int flags);
int alizmq_msg_add(aliZMQmsg* message, DataTopic* topic, AliRawData* object);

//file operations
int alizmq_file_write(AtomicFile& file, aliZMQmsg* message, bool deserializeROOTobjects=kTRUE);
int alizmq_file_write(AtomicFile& file, AliHLTDataTopic topic, TObject* object);
int alizmq_file_write(AtomicFile& file, const AliHLTDataTopic& topic, const void* buf, Int_t len);
int alizmq_file_read(TFile& file, aliZMQmsg* message);

}  //end namespace AliZMQhelpers

#endif

