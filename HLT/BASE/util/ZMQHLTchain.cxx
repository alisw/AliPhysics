/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 */

#include "AliHLTOUT.h"
#include "AliHLTMisc.h"
#include "AliHLTSystem.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTPluginBase.h"
#include "AliHLTConfiguration.h"
#include "AliHLTHOMERLibManager.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTDataSource.h"
#include "AliHLTProcessor.h"
#include "AliHLTDataSink.h"
#include "AliHLTMessage.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTCTPData.h"
#include "AliCDBManager.h"
#include "AliDAQ.h"
#include "TROOT.h"
#include "TServerSocket.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "zmq.h"
#include "AliZMQhelpers.h"
#include <vector>
#include <iostream>
#include <sys/ioctl.h>

using namespace std;

//configuration vars
Bool_t  fVerbose = kFALSE;
TString fZMQconfigIN  = "PULL";
TString fZMQsubscriptionIN = "";
TString fZMQconfigForward  = "PUSH@inproc://source";
Int_t   fZMQmaxQueueSize = 100;
Int_t   fZMQtimeout = -1;

//ZMQ stuff
void* fZMQcontext = NULL;    //ze zmq context
void* fZMQforward = NULL;        //the monitoring socket, here we publish a copy of the data
void* fZMQin  = NULL;        //the in socket - entry point for the data to be merged.

double fSleep = 0.1;  // seconds
const char* cdbpath = "local://$ALICE_ROOT/OCDB";

string fECSstring;
string fINFOstring;
string fConfigMacro;

AliHLTUInt32_t fParticipatingDetectors = 0;
int fRunNumber = 0;

const char* fUSAGE = 
    "ZMQHLTchain: run an HLT chain.\n"
    " -in : data in, zmq config string, e.g. PUSH>tcp://localhost:123123\n"
    " -sleep : how long to sleep between requests (events)\n"
    " -Verbose : be verbose\n"
    " -ZMQtimeout : when to timeout a socket when no reply comes\n"
    " -ZMQmaxQueueSize : max size of the input queue\n"
    " -config : ROOT macro defining the HLT chain.\n"
    "           a \"source\" component is provided, use as first parent\n"
    "           a \"sink\" component needs to be defined (last in chain)\n"
    "\n"
    "           example MINIMAL config is:\n"
    "\n"
    "           void config() {\n"
    "             AliHLTSystem* system = AliHLTPluginBase::GetInstance();\n"
    "             system->LoadComponentLibraries(\"libAliHLTGlobal.so\");\n"
    "             AliHLTConfiguration sink(\"sink\",\"ZMQsink\", \"source\", \"out=PUSH@ipc:///tmp/chainOut\");\n"
    "           }\n"
    ;

Int_t ProcessOptionString(TString arguments);
Int_t InitZMQ();
Int_t Run();
void FetchRunInformation(int& fRunNumber, TString& ecsParamString, TString& detectorList);
void DoRequest(void* socket);

//_______________________________________________________________________________________
Int_t InitZMQ()
{
  //init or reinit stuff
  Int_t rc = 0;
  rc = alizmq_socket_init(fZMQin,  fZMQcontext, fZMQconfigIN.Data(), fZMQtimeout, fZMQmaxQueueSize);
  printf("in:  (%s) %s\n", alizmq_socket_name(rc), fZMQconfigIN.Data());
  rc = alizmq_socket_init(fZMQforward, fZMQcontext, fZMQconfigForward.Data(), fZMQtimeout, fZMQmaxQueueSize);
  printf("out: (%s) %s\n", alizmq_socket_name(rc), fZMQconfigForward.Data());
  return 0;
}

//_______________________________________________________________________________________
AliHLTTriggerMask_t GetTriggerClasses(TObject* obj)
{
  AliHLTTriggerMask_t trigMask = 0;
  AliHLTLogging log;
  long long unsigned int eventID = 0;

  AliHLTGlobalTriggerDecision* decision;
  if(obj && (decision=dynamic_cast<AliHLTGlobalTriggerDecision*>(obj) ) ) 
  {
    const TObject* ctpobj = decision->InputObjects().FindObject("AliHLTCTPData");
    const AliHLTCTPData* ctpdata = (ctpobj != NULL ? dynamic_cast<const AliHLTCTPData*>(ctpobj) : NULL);
    if (ctpdata != NULL)
    {
      trigMask = ctpdata->Triggers();
    }
    else
    {
      printf("no CTPData!\n");
    }
  }

  return trigMask;
}

//_______________________________________________________________________________________
void DoRequest(void* socket)
{
  if (fECSstring.empty())
  {
    if (fVerbose) printf("requesting ECS param\n");
    alizmq_msg_send("ECSPARAM","",socket,0);
  }
  else
  {
    if (fVerbose) printf("requesting some data\n");
    alizmq_msg_send("","",socket,0);
  }
}

int Run()
{
  //make the source configuration
	AliHLTSystem* system = AliHLTPluginBase::GetInstance();  
	system->LoadComponentLibraries("libAliHLTGlobal.so");

  //make a source component
  AliHLTConfiguration reader("source","ZMQsource","",
                             "in=PULL>inproc://source ZMQneverBlock=0 ZMQrequestTimeout=1000000");
  
  //load the chain definition
  gROOT->Macro(fConfigMacro.c_str());

  //build the sink configuration
  system->BuildTaskList("sink");

  //now loop over events
  bool interrupted = false;
	do
	{
    double startTime =  TTimeStamp().AsDouble();

    Int_t nSockets=1;
    zmq_pollitem_t sockets[] = { 
      { fZMQin, 0, ZMQ_POLLIN, 0 },
    };

    Int_t rc = 0;
    errno=0;

    Int_t inType=alizmq_socket_type(fZMQin);
    
    //request first
    if (inType==ZMQ_REQ) DoRequest(fZMQin);

    //wait for the data
    rc = zmq_poll(sockets, nSockets, fZMQtimeout); //poll sockets
    if (rc==-1 && errno==ETERM)
    {
      //this can only happen if the context was terminated, one of the sockets are
      //not valid or operation was interrupted
      Printf("zmq_poll was interrupted! rc = %i, %s", rc, zmq_strerror(errno));
      break;
    }

    //if we time out (waiting for a response) reinit the REQ socket(s)
    if (rc==0)
    {
      if (inType==ZMQ_REQ) {
        if (fVerbose) printf("no reply from %s in %i ms, server died?\n", fZMQconfigIN.Data(), fZMQtimeout);
        rc = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data(), fZMQtimeout, fZMQmaxQueueSize);
      }
    }

    AliHLTTriggerMask_t trgMask = 0;
    AliHLTUInt32_t eventType = 0;

    //data present socket 0 - in
    if (sockets[0].revents & ZMQ_POLLIN)
    {
      aliZMQmsg message;
      //get the incoming data
      alizmq_msg_recv(&message, fZMQin, 0);

      //extract infromation from message
      for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
      {
        //get the trigger mask
        if (alizmq_msg_iter_check(i, kAliHLTDataTypeGlobalTrigger)==0)
        {
          TObject* obj = NULL;
          alizmq_msg_iter_data(i, obj);
          if (obj) { trgMask = GetTriggerClasses(obj); }
          if (fVerbose) printf("block kAliHLTDataTypeGlobalTrigger found\n");
        }
        //get the ECS param string
        else if (alizmq_msg_iter_check(i, kAliHLTDataTypeECSParam)==0)
        {
          alizmq_msg_iter_data(i, fECSstring);
          stringMap ecsParamMap = ParseParamString(fECSstring);
          TString detectorList = ecsParamMap["DETECTOR_LIST"];
          detectorList.ReplaceAll(",", " ");
          fParticipatingDetectors = AliDAQ::DetectorPatternOffline(detectorList)
                                   xor AliDAQ::DetectorPatternOffline("DAQ_TEST");
          if (fVerbose) {
            printf("block kAliHLTDataTypeECSParam found\n");
            printf("ECS string: %s\n", fECSstring.c_str());
          }

        }
        //get the INFO block
        else if (alizmq_msg_iter_check(i, kAliHLTDataTypeInfo))
        {
          if (fVerbose) printf("block kAliHLTDataTypeInfo found\n");
          alizmq_msg_iter_data(i, fINFOstring);
          stringMap infoParams = ParseParamString(fINFOstring);
          int runNumber = atoi(infoParams["run"].c_str());
          if (runNumber!=fRunNumber && fRunNumber>0)
          {
            //if run changes, invalidate all params, new ones will be requested
            fRunNumber = runNumber;
            fECSstring.erase();
            if (fVerbose) printf("run changed! old: %i, new: %i\n", fRunNumber, runNumber);
            interrupted=true;
          }
        }
        //get the event info
        else if (alizmq_msg_iter_check(i, kAliHLTDataTypeEvent))
        {
          AliHLTDataTopic topic;
          alizmq_msg_iter_topic(i, topic);
          eventType = topic.fSpecification;
          if (fVerbose) printf("block kAliHLTDataTypeEvent found, event type = %i\n",eventType);
        }
      }

      //forward the data to the chain, force non blocking mode
      if (fVerbose) printf("forwarding message to %s\n", fZMQconfigForward.Data());
      int rc = alizmq_msg_send(&message, fZMQforward, ZMQ_DONTWAIT);
      if (fVerbose) printf("...forwarded, rc=%i\n",rc);
      //close the message
      if (rc<0) alizmq_msg_close(&message);

      AliHLTSystem* system = AliHLTPluginBase::GetInstance();

      TString opt = "ECS=";
      opt += fECSstring;
      system->ScanOptions(opt.Data());

      AliCDBManager::Instance()->SetRun(fRunNumber);
      AliHLTMisc::Instance().InitMagneticField();		

      int nEvents = 1;
      int stop = 0;

      if (fVerbose) {
        printf("running chain\n");
      }
      AliHLTUInt32_t timestamp = AliHLTUInt32_t( TTimeStamp().AsDouble() );
      system->Run(nEvents, stop, trgMask, timestamp, eventType, fParticipatingDetectors);
    }
	
		double endTime = TTimeStamp().AsDouble();
		double timeDiff = endTime - startTime;
		if (timeDiff < 1.0 and timeDiff < fSleep) {
      gSystem->Sleep(int((fSleep - timeDiff)*1000));
    }
	}
	while (!interrupted);
  return 0;
}

//______________________________________________________________________________
Int_t ProcessOptionString(TString arguments)
{
  //process passed options
  Int_t nOptions=0;
  aliStringVec* options = AliOptionParser::TokenizeOptionString(arguments);
  for (aliStringVec::iterator i=options->begin(); i!=options->end(); ++i)
  {
    //Printf("  %s : %s", i->first.data(), i->second.data());
    const TString& option = i->first; 
    const TString& value = i->second;
    if (option.EqualTo("in"))
    {
      fZMQconfigIN = value;
    }
    else if (option.EqualTo("Verbose"))
    {
      fVerbose=kTRUE;
    }
    else if (option.EqualTo("ZMQmaxQueueSize"))
    {
      fZMQmaxQueueSize=value.Atoi();
    }
    else if (option.EqualTo("ZMQtimeout"))
    {
      fZMQtimeout=value.Atoi();
    }
    else if (option.EqualTo("sleep"))
    {
      fSleep=value.Atoi();
    }
    else if (option.EqualTo("config"))
    {
      fConfigMacro = value;
    }
    else
    {
      Printf("unrecognized option |%s|",option.Data());
      nOptions=-1;
      break;
    }
    nOptions++;
  }
  delete options; //tidy up

  return nOptions; 
}

//_______________________________________________________________________________________
int main(Int_t argc, char** argv)
{
  Int_t mainReturnCode=0;

  //process args
  TString argString = AliOptionParser::GetFullArgString(argc,argv);
  if (ProcessOptionString(argString)<=0)
  {
    printf("%s",fUSAGE);
    return 1;
  }

	AliCDBManager::Instance()->SetDefaultStorage(cdbpath);

  //globally enable schema evolution for serializing ROOT objects
  TMessage::EnableSchemaEvolutionForAll(kTRUE);
  //the context
  fZMQcontext = alizmq_context();

  //init stuff
  if (InitZMQ()<0) {
    Printf("failed init");
    return 1;
  }

  Run();

  //destroy ZMQ sockets
  zmq_close(fZMQin);
  zmq_close(fZMQforward);
  zmq_ctx_destroy(fZMQcontext);
  return mainReturnCode;
}
