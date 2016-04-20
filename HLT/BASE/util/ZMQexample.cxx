#include "zmq.h"
#include "AliZMQhelpers.h"
#include "AliHLTDataTypes.h"
#include "TObject.h"

void* fChannelIN;
TString fConfigIN;
const char* fUSAGE = "example, needs one argument, e.g. in=PULL@ipc:///tmp/examplesource\n";

int ProcessOptionString(TString arguments);

//_______________________________________________________________________________________
int main(int argc, char** argv)
{
  //process args
  int noptions = ProcessOptionString(AliOptionParser::GetFullArgString(argc,argv));
  if (noptions<=0) 
  {
    printf("%s",fUSAGE);
    return 1;
  }

  alizmq_socket_init(fChannelIN, alizmq_context(), fConfigIN.Data());
  aliZMQmsg message;
  int size = alizmq_msg_recv(&message,fChannelIN,0);
  if (size>0) printf("__________________________\n");
  for (aliZMQmsg::iterator block=message.begin(); block!=message.end(); ++block)
  {
    AliHLTDataTopic header;
    alizmq_msg_iter_topic(block, header);
    TObject* object;
    alizmq_msg_iter_data(block,object);
    
    const char* name = "";
    if (object) {
      name = object->GetName();
      delete object;
    }

    std::string desc = header.Description();
    printf("received: %s name: %s\n", desc.c_str(), name);
  }
  alizmq_msg_close(&message);

  return 0;
}

//______________________________________________________________________________
int ProcessOptionString(TString arguments)
{
  //process passed options
  aliStringVec* options = AliOptionParser::TokenizeOptionString(arguments);
  int nOptions = 0;
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //Printf("  %s : %s", i->first.data(), i->second.data());
    TString option = i->first;
    TString value = i->second;
    if ( option.EqualTo("in") )
    {
      fConfigIN = value;
    }
    else
    {
      nOptions=-1;
      break;
    }
    nOptions++;
  }
  delete options; //tidy up

  return nOptions; 
}

