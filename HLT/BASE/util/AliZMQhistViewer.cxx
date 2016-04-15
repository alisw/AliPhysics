/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliZMQhistViewer.cxx
@author  Mikolaj Krzewicki (mkrzewic@cern.ch)
*/

#include "zmq.h"
#include <algorithm>
#include <iostream>
#include "AliHLTDataTypes.h"
#include "AliHLTComponent.h"
#include "AliHLTMessage.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TMap.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TDirectory.h"
#include "TList.h"
#include "AliZMQhelpers.h"
#include "TMessage.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH1F.h"
#include <time.h>
#include <string>
#include <map>
#include "TSystem.h"
#include "TStyle.h"
#include "signal.h"
#include "AliZMQhistViewer.h"
#include "TThread.h"

using namespace std;

ClassImpQ(AliZMQhistViewer)

//_______________________________________________________________________________________
AliZMQhistViewer::AliZMQhistViewer():
  fCanvas(NULL),
  fVerbose ( kFALSE ),
  fZMQconfigIN  ( "PULL>tcp://localhost:60211" ),
  fZMQsocketModeIN(1 ),
  fZMQconfigCONFIG("PULL>inproc://viewerConfig"),
  fFilter ( "" ),
  fPollInterval ( 0),
  fPollTimeout ( 100000),
  fSort ( kTRUE),
  fZMQcontext ( NULL),
  fZMQin  ( NULL),
  fZMQconfig( NULL),
  fZMQsleeper(NULL),
  fTrashQueue(),
  fContent(),
  fStatus ( ""),
  fRunNumber ( 0),
  fSelectionRegexp ( NULL),
  fUnSelectionRegexp ( NULL),
  fDrawOptions(),
  fScaleLogX ( kFALSE),
  fScaleLogY ( kFALSE),
  fScaleLogZ ( kFALSE),
  fResetOnRequest ( kFALSE),
  fHistStats (0),
  fAllowResetAtSOR ( kTRUE),
  iterations(0),
  fIncoming(NULL),
  fInfo(),
  fClearCanvas(),
  fUpdateCanvas(true),
  fTerminated(false)
{
  //init stuff
  //ZMQ init
  fZMQcontext = alizmq_context();
  alizmq_socket_init(fZMQconfig, fZMQcontext, fZMQconfigCONFIG.Data());
}

//_______________________________________________________________________________________
int AliZMQhistViewer::Run(void* arg)
{
  if (!fCanvas) fCanvas = new TCanvas();
  fCanvas->SetBorderSize(10);
  if (!fIncoming) fIncoming = new std::vector<ZMQviewerObject>;
  fIncoming->reserve(1000);

  //main loop
  while(!GetTerminated())
  {
    int sleepInterrupt=0;
    errno=0;
    //send a request if we are using REQ
    if (fZMQsocketModeIN==ZMQ_REQ)
    {
      string request;

      if (fSelectionRegexp || fUnSelectionRegexp) 
      {
        if (fSelectionRegexp) request += " select="+fSelectionRegexp->GetPattern();
        if (fUnSelectionRegexp) request += " unselect="+fUnSelectionRegexp->GetPattern();
        if (fResetOnRequest) request += " ResetOnRequest";
        alizmq_msg_send("CONFIG", request, fZMQin, ZMQ_SNDMORE);
      }

      if (fVerbose) Printf("sending request CONFIG %s", request.c_str());
      alizmq_msg_send("", "", fZMQin, 0);
    }

    //wait for the data
    zmq_pollitem_t sockets[] = { 
      { fZMQin, 0, ZMQ_POLLIN, 0 }, 
      { fZMQconfig, 0, ZMQ_POLLIN, 0 }, 
    };
    int rc = zmq_poll(sockets, 2, (fZMQsocketModeIN==ZMQ_REQ)?fPollTimeout:-1);

    if (fVerbose) printf("poll exits with rc %i, \"%s\"\n", rc, (rc<0)?zmq_strerror(errno):0);

    if (rc==-1 && errno==ETERM)
    {
      if (fVerbose) Printf("ZMQ context terminated, jumping out");
      break;
    }

    if (!(sockets[0].revents & ZMQ_POLLIN))
    {
      //server died
      Printf("connection timed out, server %s died?", fZMQconfigIN.Data());
      fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, fZMQconfigIN.Data());
      if (fZMQsocketModeIN < 0) return -1;
      continue;
    }
    else if (sockets[0].revents & ZMQ_POLLIN) //ZMQin
    {
      GetData(fZMQin);
      if (fUpdateCanvas) {
        UpdateCanvas(fCanvas);
        fCanvas->Update();
        gSystem->ProcessEvents();
      }
    }
    else if (sockets[1].revents & ZMQ_POLLIN) //ZMQconfig
    {
      aliZMQmsg message; //only string messages expected here
      alizmq_msg_recv(&message, fZMQconfig, 0);
      for (aliZMQmsg::iterator part=message.begin(); part!=message.end(); ++part)
      {
        //we just take the payloads, in this channel we dont expect anything else than
        //config strings
        string configString;
        alizmq_msg_iter_data(part, configString);
        ProcessOptionString(configString.c_str());
      }
    }
    sleepInterrupt=0;
    //sleep
    rc = zmq_recv(fZMQsleeper, &sleepInterrupt, sizeof(sleepInterrupt), 0);
    if (fVerbose) { printf("    sleeper rc: %i (%s)\n", rc, zmq_strerror(errno)); }
  }//main loop

  if (fVerbose) printf("trying to close the ZMQ socket %s\n", fZMQconfigIN.Data());
  zmq_close(fZMQin); fZMQin=NULL;
  zmq_close(fZMQsleeper); fZMQsleeper=NULL;
  zmq_close(fZMQconfig); fZMQconfig=NULL;
  if (fVerbose) printf("sockets closed\n");

  GetIncoming(NULL); //destroy the buffer
  return 0;
}

//_______________________________________________________________________________________
int AliZMQhistViewer::GetData(void* socket)
{
  //get all data (topic+body), possibly many of them
  aliZMQmsg message;
  alizmq_msg_recv(&message, socket, 0);

  vector<ZMQviewerObject>* incomingObjects = new vector<ZMQviewerObject>;
  incomingObjects->reserve(message.size());

  //process message, deserialize objects, puth them in the container 
  for (aliZMQmsg::iterator i=message.begin(); i!=message.end(); ++i)
  {
    if (alizmq_msg_iter_check(i, "INFO")==0)
    {
      //check if we have a runnumber in the string
      string info;
      alizmq_msg_iter_data(i,info);
      if (fVerbose) Printf("processing INFO %s", info.c_str());
      GetInfo(&info);

      size_t runTagPos = info.find("run");
      if (runTagPos != std::string::npos)
      {
        size_t runStartPos = info.find("=",runTagPos);
        size_t runEndPos = info.find(" ");
        string runString = info.substr(runStartPos+1,runEndPos-runStartPos-1);

        if (fVerbose) printf("received run=%s\n",runString.c_str());

        int runnumber = atoi(runString.c_str());

        if (runnumber!=fRunNumber && fAllowResetAtSOR) 
        {
          if (fVerbose) printf("Run changed, resetting!\n");
          GetIncoming(NULL); //this clears it
          bool clear = true;
          GetClearCanvas(&clear);
        }
        fRunNumber = runnumber; 
      }
      continue;
    }
    if (alizmq_msg_iter_check(i,kAliHLTDataTypeStreamerInfo)==0)
    {
      if (fVerbose) printf("extracting schema\n");
      //alizmq_msg_iter_init_streamer_infos(i);
      continue;
    }

    //incoming ROOT objects
    {
      TObject* tmp = NULL;
      alizmq_msg_iter_data(i, tmp);
      if (!tmp) continue;
      ZMQviewerObject object(tmp);

      if (fVerbose) Printf("--in: %s (%s), %p", object.name.c_str(), tmp->ClassName(), tmp);

      incomingObjects->push_back(object);

    }

  } //for iterator i
  alizmq_msg_close(&message);

  GetIncoming(incomingObjects);

  DataReady(); //emit a signal

  return 0;
}

//______________________________________________________________________________
int AliZMQhistViewer::UpdateCanvas(TCanvas* canvas,
                                   TPRegexp* selectionRegexp,
                                   TPRegexp* unSelectionRegexp)
{
  //use the thread safe access
  string info = GetInfo();
  std::vector<ZMQviewerObject>* incomingObjects = GetIncoming();
  bool forceClearCanvas = GetClearCanvas();
  if (forceClearCanvas && fVerbose) printf("##forceClearCanvas\n"); 
  if (!incomingObjects) return 1;

  if (forceClearCanvas && fContent) {
    canvas->Clear();
    //canvas->Update();
    fTrashQueue.push_back(fContent);
    fContent = NULL;
  }

  if (!fContent) {
    fContent = new vector<ZMQviewerObject>;
    fContent->reserve(incomingObjects->size());
  }

  gStyle->SetOptStat(fHistStats);

  const char* drawOptions = fDrawOptions;
  bool verbose = fVerbose;

  int nNewPlots = 0;
  for (vector<ZMQviewerObject>::iterator incoming=incomingObjects->begin();
      incoming!=incomingObjects->end();
      ++incoming)
  {
    Bool_t selected = kTRUE;
    Bool_t unselected = kFALSE;
    if (selectionRegexp) selected = selectionRegexp->Match(incoming->object->GetName());
    if (unSelectionRegexp) unselected = unSelectionRegexp->Match(incoming->object->GetName());
    if (!selected || unselected)
    {
      if (fVerbose) printf("object %s did not pass selection\n", incoming->object->GetName());
      continue;
    }

    //check if it will replace something
    bool found=false;
    incoming->redraw=true;
    for (vector<ZMQviewerObject>::iterator current=fContent->begin();
        current!=fContent->end();
        ++current)
    {
      if (incoming->name.compare(current->name)==0)
      {
        if (fVerbose) printf("  updating %s\n",incoming->name.c_str());
        found = true;
        if (!(incoming->object)) { printf("WARNING: null pointer in input data"); }
        else { 
          current->SwapObject(*incoming);
          if (verbose) printf("  swapping plot %s\n", incoming->name.c_str());
        }
        break;
      }
    }
    if (!found)
    {
      if (fVerbose) printf("  new object %s\n",incoming->name.c_str());
      ZMQviewerObject tmp(*incoming); incoming->object=NULL;
      fContent->push_back(tmp); //completely new one
      nNewPlots++;
    }
  }

  //re-sort the new list if we have new plots
  if (nNewPlots>0 && fSort) {
    if (fVerbose) printf("  re-sorting and setting new pad numbers\n");
    std::sort(fContent->begin(), fContent->end(), ZMQviewerObjectTitleComparator());
    for (int i=0; i<fContent->size(); i++) {
      (*fContent)[i].pad=i+1;
      (*fContent)[i].redraw=true;
    }
  }

  if (nNewPlots>0)
  {
    forceClearCanvas = true;
  }
  
  canvas->cd();

  //after we clear the canvas, the pads are gone, clear the pad cache as well
  int nDrawables = fContent->size();
  int nPads = CountPads(canvas);

  if (verbose) printf("  nDrawables: %i, nPads: %i\n", nDrawables, nPads);
  if (nDrawables > nPads || forceClearCanvas)
  {
    canvas->SetEditable(kTRUE);
    canvas->Clear();
    canvas->DivideSquare(nDrawables);
    canvas->SetEditable(kFALSE);
    nPads = CountPads(canvas);
    if (verbose) printf("  reorganizing canvas, now %i pads\n",nPads);
  }

  //if we have new plots we need to replot the whole thing as the old plots may have been
  //moved
  //if we have no new plots, we just draw the incoming ones at correct locations
  {
    if (verbose) printf("  drawing (%lu) plots\n", fContent->size());
    for (vector<ZMQviewerObject>::iterator current=fContent->begin(); current!=fContent->end(); ++current)
    {
      if (!(current->redraw)) {
        if (verbose) printf("  not replotting %s\n", current->name.c_str());
        continue;
      }
      if (verbose) printf("  plotting %s at pad %i\n", current->name.c_str(),current->pad);
      TVirtualPad* pad = canvas->cd(current->pad);
      if (!pad) pad = canvas->cd();
      if (fScaleLogX) gPad->SetLogx();
      if (fScaleLogY) gPad->SetLogy();
      if (fScaleLogZ) gPad->SetLogz();
      pad->SetEditable(kFALSE);
      if (current->object)
      {
        if (verbose) printf("  removing previous %s at %p\n", current->object->GetName(), current->previous);
        pad->SetName(current->object->GetName());
        pad->RecursiveRemove(current->previous);
        pad->GetListOfPrimitives()->Add(current->object,drawOptions);
        pad->Modified();
      }
      else { if (verbose) printf("  missing object\n"); }
    }
  }
  canvas->SetTitle(info.c_str());
  
  fTrashQueue.push_back(incomingObjects);
  //garbage collect
  if (fTrashQueue.size()>0)
  {
    if (verbose) printf("    garbage collecting\n");
    vector<ZMQviewerObject>* pieceOfTrash = *(fTrashQueue.begin());
    for (vector<ZMQviewerObject>::iterator i=pieceOfTrash->begin();
                                           i!=pieceOfTrash->end(); 
                                           ++i) 
    {
      if (verbose) printf("    destroying %s, %p\n",i->name.c_str(),i->object);
      delete i->object;
      //if (verbose) printf("destroying previous %s, %p\n",i->name.c_str(),i->previous);
      //delete i->previous;
    }
    fTrashQueue.pop_front();
  }

  return 0;
}

//______________________________________________________________________________
int AliZMQhistViewer::ProcessOption(TString option, TString value)
{
  int ret = 0;
  fZMQcontext = alizmq_context();
  TH1::AddDirectory(kFALSE);
  TDirectory::AddDirectory(kFALSE);

  //process passed options
  if (option.EqualTo("PollInterval") || option.EqualTo("sleep"))
  {
    fPollInterval = round(value.Atof()*1e3);
    if (fPollInterval<0) {fPollInterval = -1;}
  }
  else if (option.EqualTo("PollTimeout") || option.EqualTo("timeout"))
  {
    fPollTimeout = round(value.Atof()*1e3);
  }
  else if (option.EqualTo("ZMQconfigIN") || option.EqualTo("in") )
  {
    string valuestr = value.Data();
    GetZMQconfig(&valuestr);
    fZMQsocketModeIN = alizmq_socket_init(fZMQin, fZMQcontext, value.Data(), -1, 2);
    if (fZMQsocketModeIN < 0) return -1;
  }
  else if (option.EqualTo("Verbose"))
  {
    fVerbose=kTRUE;
  }
  else if (option.EqualTo("select"))
  {
    std::string tmp = value.Data();
    GetSelection(&tmp);
  }
  else if (option.EqualTo("unselect"))
  {
    std::string tmp = value.Data();
    GetUnSelection(&tmp);
  }
  else if (option.EqualTo("ResetOnRequest"))
  {
    fResetOnRequest = kTRUE;
  }
  else if (option.EqualTo("drawoptions"))
  {
    fDrawOptions = value;
  }
  else if (option.EqualTo("logx"))
  {
    fScaleLogX=kTRUE;
  }
  else if (option.EqualTo("logy"))
  {
    fScaleLogY=kTRUE;
  }
  else if (option.EqualTo("logz"))
  {
    fScaleLogZ=kTRUE;
  }
  else if (option.EqualTo("sort"))
  {
    fSort=value.Contains(0)?kFALSE:kTRUE;
  }
  else if (option.EqualTo("histstats"))
  {
    fHistStats = value.Atoi();
    gStyle->SetOptStat(fHistStats);
  }
  else if (option.EqualTo("AllowResetAtSOR"))
  {
    fAllowResetAtSOR = (option.Contains("0")||option.Contains("no"))?kFALSE:kTRUE;
  }
  else
  {
    ret = -1;
  }
  int rc = alizmq_socket_init(fZMQsleeper, fZMQcontext, "PULL@inproc://sleep", fPollInterval);
  return ret; 
}

//______________________________________________________________________________
Int_t AliZMQhistViewer::CountPads(TVirtualPad *pad) {
  //count the number of pads in pad
  if (!pad) return 0;
  Int_t npads = 0;
  TObject *obj;
  TIter next(pad->GetListOfPrimitives());
  while ((obj = next())) {
    if (obj->InheritsFrom(TVirtualPad::Class())) npads++;
  }
  return npads;
}

//______________________________________________________________________________
std::vector<ZMQviewerObject>*
AliZMQhistViewer::GetIncoming(std::vector<ZMQviewerObject>* in)
{
  std::vector<ZMQviewerObject>* tmp = NULL;
  TThread::Lock();
  if (!in) { //POP
    tmp = fIncoming;
    fIncoming = NULL;
  } else { //PUSH, this will destroy anything that was not consumed
    if (fIncoming) {
      for (std::vector<ZMQviewerObject>::iterator i=fIncoming->begin(); i!=fIncoming->end(); ++i)
      {
        if (i->object) { if (fVerbose) printf("destroying non consumed object %s %p\n",i->name.c_str(), i->object); }
        delete i->object; i->object = NULL;
      }
      delete fIncoming;
    }
    fIncoming = in;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
std::string AliZMQhistViewer::GetInfo(std::string* in ) {
  string tmp;
  TThread::Lock();
  if (!in) { 
    tmp = fInfo;
  } else {
    fInfo.assign(*in);
    tmp = fInfo;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
bool AliZMQhistViewer::GetClearCanvas(bool* in ) {
  bool tmp;
  TThread::Lock();
  if (!in) {
    tmp = fClearCanvas;
    fClearCanvas = false;
  } else {
    fClearCanvas = *in;;
    tmp = fClearCanvas;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
bool AliZMQhistViewer::GetTerminated(bool* in) {
  bool tmp;
  TThread::Lock();
  if (!in) {
    tmp = fTerminated;
    fTerminated = false;
  } else {
    fTerminated = *in;
    tmp = fTerminated;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
bool AliZMQhistViewer::GetUpdateCanvas(bool* in) {
  bool tmp;
  TThread::Lock();
  if (!in) {
    tmp = fUpdateCanvas;
    fUpdateCanvas = false;
  } else {
    fUpdateCanvas = *in;
    tmp = fUpdateCanvas;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
int AliZMQhistViewer::GetPollInterval(int* in) {
  int tmp;
  TThread::Lock();
  if (!in) {
    tmp = fPollInterval;
  } else {
    fPollInterval = *in;
    tmp = fPollInterval;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
void AliZMQhistViewer::DataReady()
{
  TThread::Lock();
  Emit("DataReady()");
  TThread::UnLock();
}

//______________________________________________________________________________
string AliZMQhistViewer::GetZMQconfig(string* in)
{
  string tmp;
  TThread::Lock();
  if (!in) {
    tmp = fZMQconfigIN.Data();
  } else {
    fZMQconfigIN = *in;
    tmp = *in;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
TPRegexp* AliZMQhistViewer::GetSelection(string* in)
{
  TPRegexp* tmp = NULL;
  TThread::Lock();
  if (!in) {
    if (fSelectionRegexp) {
      tmp = new TPRegexp(*fSelectionRegexp);
    }
  } else {
    delete fSelectionRegexp;
    fSelectionRegexp = new TPRegexp(in->c_str());
    tmp = fSelectionRegexp;
  }
  TThread::UnLock();
  return tmp;
}

//______________________________________________________________________________
TPRegexp* AliZMQhistViewer::GetUnSelection(string* in)
{
  TPRegexp* tmp = NULL;
  TThread::Lock();
  if (!in) {
    if (fUnSelectionRegexp) {
      tmp = new TPRegexp(*fUnSelectionRegexp);
    }
  } else {
    delete fUnSelectionRegexp;
    fUnSelectionRegexp = new TPRegexp(in->c_str());
    tmp = fUnSelectionRegexp;
  }
  TThread::UnLock();
  return tmp;
}


