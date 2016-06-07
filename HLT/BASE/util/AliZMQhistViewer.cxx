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
#include "TNamed.h"
#include "signal.h"
#include "AliZMQhistViewer.h"
#include "TThread.h"
#include "AliAnalysisDataContainer.h"

#define safename(i) (i->object)?i->object->GetName():"NULL"
#define safetitle(i) (i->object)?i->object->GetTitle():"NULL"

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
        UpdateCanvas(fCanvas,fSelectionRegexp, fUnSelectionRegexp);
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
    if (alizmq_msg_iter_check_id(i, "INFO")==0)
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
    if (alizmq_msg_iter_check_id(i,kAliHLTDataTypeStreamerInfo)==0)
    {
      if (fVerbose) printf("extracting schema\n");
      alizmq_msg_iter_init_streamer_infos(i);
      continue;
    }

    //incoming ROOT objects
    {
      TObject* tmp = NULL;
      alizmq_msg_iter_data(i, tmp);
      if (!tmp) continue;

      std::vector<TObject*> listOfObjects; listOfObjects.reserve(100);
      AliAnalysisDataContainer* analKont = dynamic_cast<AliAnalysisDataContainer*>(tmp);
      if (analKont) {
        if (fVerbose) printf("--in AliAnalysisDataContainer %p\n",analKont);
        GetObjects(analKont, &listOfObjects);
        if (fVerbose) printf("  destroying anal container %p\n",tmp);
        delete tmp;
      } else {
        TCollection* collection = dynamic_cast<TCollection*>(tmp);
        if (collection) {
          if (fVerbose) printf("--in TCollection %p\n",collection);
          GetObjects(collection, &listOfObjects);
        if (fVerbose) printf("  destroying collection %p\n",tmp);
          delete tmp;
        } else {
          if (fVerbose) printf("--in TObject\n");
          listOfObjects.push_back(tmp);
        }
      }

      //add all extracted objects to the list of veiwer objects
      for (std::vector<TObject*>::iterator i=listOfObjects.begin(); i!=listOfObjects.end(); ++i) {
        ZMQviewerObject object(*i);
        if (fVerbose) Printf("  adding: %s (%s), %p", (*i)->GetName(), (*i)->ClassName(), *i);
        incomingObjects->push_back(object);
      }
    }
  } //for iterator i
  alizmq_msg_close(&message);

  GetIncoming(incomingObjects);

  DataReady(); //emit a signal

  return 0;
}

//______________________________________________________________________________
int AliZMQhistViewer::GetObjects(AliAnalysisDataContainer* kont, std::vector<TObject*>* list, const char* prefix)
{
  const char* analName = kont->GetName();
  TObject* analData = kont->GetData();
  std::string name = analName;
  std::string namePrefix = name + "/";
  TCollection* collection = dynamic_cast<TCollection*>(analData);
  if (collection) {
    if (fVerbose) Printf("  have a collection %p",collection);
    const char* collName = collection->GetName();
    GetObjects(collection, list, namePrefix.c_str());
    if (fVerbose) printf("  destroying collection %p\n",collection);
    delete collection;
    kont->SetDataOwned(kFALSE);
  } else { //if (collection)
    TNamed* named = dynamic_cast<TNamed*>(analData);
    name = namePrefix + analData->GetName();
    std::string title = namePrefix + analData->GetTitle();
    if (named) {
      named->SetName(name.c_str());
      named->SetTitle(title.c_str());
    }
    if (fVerbose) Printf("--in (from analysis container): %s (%s), %p",
                         named->GetName(),
                         named->ClassName(),
                         named );
    kont->SetDataOwned(kFALSE);
    list->push_back(analData);
  }
  return 0;
}

//______________________________________________________________________________
int AliZMQhistViewer::GetObjects(TCollection* collection, std::vector<TObject*>* list, const char* prefix)
{
  TIter next(collection);
  while (TObject* tmp = next()) {
    collection->Remove(tmp);
    std::string name = tmp->GetName();
    name = prefix + name;
    if (fVerbose) Printf("--in (from a TCollection): %s (%s), %p",
                         tmp->GetName(), tmp->ClassName(), tmp);
    AliAnalysisDataContainer* analKont = dynamic_cast<AliAnalysisDataContainer*>(tmp);
    if (analKont) {
      if (fVerbose) Printf("  have an analysis container %p",analKont);
      GetObjects(analKont,list,name.c_str());
      if (fVerbose) printf("  destroying anal container %p\n",analKont);
      delete analKont;
    } else {
      TNamed* named = dynamic_cast<TNamed*>(tmp);
      if (named) {
        name = named->GetName();
        name = prefix + name;
        std::string title = named->GetTitle();
        title = prefix + title;
        named->SetName(name.c_str());
        named->SetTitle(title.c_str());
      }
      list->push_back(tmp);
    }
    collection->SetOwner(kTRUE);
  } //while
  return 0;
}

//______________________________________________________________________________
int AliZMQhistViewer::UpdateCanvas(TCanvas* canvas,
                                   TPRegexp* selectionRegexp,
                                   TPRegexp* unSelectionRegexp)
{
  //use the thread safe access
  std::vector<ZMQviewerObject>* incomingObjects = GetIncoming();
  if (!incomingObjects) return 1;
  
  string info = GetInfo();
  bool forceClearCanvas = GetClearCanvas();
  if (forceClearCanvas && fVerbose) printf("##forceClearCanvas\n"); 

  if (forceClearCanvas && fContent) {
    canvas->Clear();
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
  int nPlotsToDraw = 0;
  for (vector<ZMQviewerObject>::iterator incoming=incomingObjects->begin();
      incoming!=incomingObjects->end();
      ++incoming)
  {
    Bool_t selected = kTRUE;
    Bool_t unselected = kFALSE;
    if (selectionRegexp) selected = selectionRegexp->Match(safename(incoming));
    if (unSelectionRegexp) unselected = unSelectionRegexp->Match(safename(incoming));
    if (!selected || unselected)
    {
      if (fVerbose) printf("object %s did not pass selection (%s) && !(%s)\n",
                           safename(incoming), (selectionRegexp)?selectionRegexp->GetPattern().Data():"",
                           (unSelectionRegexp)?unSelectionRegexp->GetPattern().Data():"");
      incoming->redraw = false;
      continue;
    }
    else
    {
      incoming->redraw=true;
    }

    //check if it will replace something
    bool found=false;
    for (vector<ZMQviewerObject>::iterator current=fContent->begin();
        current!=fContent->end();
        ++current)
    {
      if (strcmp(safename(incoming), safename(current))==0)
      {
        if (fVerbose) printf("  updating %s %s\n",safename(incoming), safetitle(incoming));
        found = true;
        if (!(incoming->object)) { printf("WARNING: null pointer in input data"); }
        else {
          if (current->SwapObject(*incoming)) {
            if (verbose) printf("  swapped plot %s %s\n", safename(incoming),safetitle(incoming));
          } else {
            printf("duplicate histogram in incoming message! skipping %s\n", safename(incoming));
          }
        }
        break;
      }
    }
    if (!found)
    {
      if (fVerbose) printf("  new object %s\n",safename(incoming));
      ZMQviewerObject tmp(*incoming); incoming->object=NULL;
      tmp.options = drawOptions;
      fContent->push_back(tmp); //completely new one
      nNewPlots++;
    }
  }

  //MUST reset isnew after updating
  for (vector<ZMQviewerObject>::iterator current=fContent->begin();
      current!=fContent->end();
      ++current)
  {
    current->isnew = false;
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
        if (verbose) printf("  not replotting %s\n", safename(current));
        continue;
      }
      if (verbose) printf("  plotting %s at pad %i\n", safename(current),current->pad);
      TVirtualPad* pad = canvas->cd(current->pad);
      if (!pad) pad = canvas->cd();
      if (fScaleLogX) gPad->SetLogx();
      if (fScaleLogY) gPad->SetLogy();
      if (fScaleLogZ) gPad->SetLogz();
      pad->SetEditable(kFALSE);
      if (current->object)
      {
        if (verbose) printf("  removing previous %s at %p\n", safename(current), current->previous);
        pad->SetName(safename(current));
        pad->RecursiveRemove(current->previous);
        pad->GetListOfPrimitives()->Add(current->object,current->options.c_str());
        pad->Modified(kTRUE);
      }
      else { if (verbose) printf("  missing object\n"); }
    }
  }
  canvas->SetTitle(info.c_str());
  canvas->Update();
  
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
      if (verbose) printf("    destroying %s, %p\n",safename(i),i->object);
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
  if (fZMQsocketModeIN == ZMQ_REQ && fPollInterval<=0) fPollInterval=60000;
  if (fZMQsocketModeIN != ZMQ_REQ ) fPollInterval=0;
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
        if (i->object) { if (fVerbose) printf("destroying non consumed object %s %p\n",safename(i), i->object); }
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
    delete fSelectionRegexp; fSelectionRegexp = NULL;
    if (in->size()>0) fSelectionRegexp = new TPRegexp(in->c_str());
    tmp = fSelectionRegexp;
    if (fVerbose) printf("GetSelection: %s\n",(in)?in->c_str():"NULL");
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
    delete fUnSelectionRegexp; fUnSelectionRegexp = NULL;
    if (in->size()>0) fUnSelectionRegexp = new TPRegexp(in->c_str());
    tmp = fUnSelectionRegexp;
    if (fVerbose) printf("GetUnSelection: %s\n",(in)?in->c_str():"NULL");
  }
  TThread::UnLock();
  return tmp;
}


