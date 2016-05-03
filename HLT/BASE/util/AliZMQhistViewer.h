/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliZMQhistViewer.cxx
@author  Mikolaj Krzewicki (mkrzewic@cern.ch)
*/

#ifndef ALIZMQHISTVIEWER_H
#define ALIZMQHISTVIEWER_H

#include "AliZMQhelpers.h"
#include <vector>
#include <string>
#include <memory>
#include "TThread.h"
#include <deque>
#include "TQObject.h"

class TVirtualPad;
class TString;
class TObject;
class TPad;
class TCanvas;
class TPRegexp;

struct ZMQviewerObject {
  TObject* object;
  TObject* previous;
  int pad;
  bool redraw;
  bool isnew;
  
  ZMQviewerObject() : object(NULL), previous(NULL), pad(-1), redraw(true), isnew(true) {}
  ZMQviewerObject(TObject* o) : object(o), previous(NULL), 
                                pad(-1), redraw(true), isnew(true) {}
  ZMQviewerObject(const ZMQviewerObject& o) : object(o.object), previous(o.previous),
                                              pad(o.pad), 
                                              redraw(o.redraw),
                                              isnew(o.isnew) {}
  
  ZMQviewerObject& operator=(const ZMQviewerObject& o) {
    object = o.object;
    previous = o.previous;
    pad = o.pad;
    redraw = o.redraw;
    isnew = o.isnew;
    return *this;
  }

  bool SwapObject(ZMQviewerObject& from) {
    if (isnew) { return false; }
    if (from.object) { 
      previous = object;
      std::swap(object, from.object);
      isnew = true;
      return true;
    }
    return false;
  }

  ~ZMQviewerObject() {}
};

struct AliZMQhistViewer : public AliOptionParser, public TQObject {
  public:
  //methods
  AliZMQhistViewer();
  virtual ~AliZMQhistViewer() {}
  int Run(void* arg);
  int ProcessOption(TString option, TString value);
  void SetCanvas(TCanvas* canvas) {fCanvas = canvas;}

  //signals
  void DataReady(); //*SIGNAL*

  //thread safe stuff
  int UpdateCanvas(TCanvas* canvas, TPRegexp* sel=NULL, TPRegexp* unsel=NULL);
  int GetData(void* socket);
  std::vector<ZMQviewerObject>* GetIncoming(std::vector<ZMQviewerObject>* in = NULL);
  std::string GetInfo(std::string* in = NULL);
  bool GetClearCanvas(bool* in = NULL);
  bool GetUpdateCanvas(bool* in = NULL);
  bool GetTerminated(bool* in = NULL);
  int GetPollInterval(int* in = NULL);
  std::string GetZMQconfig(std::string* in=NULL);
  TPRegexp* GetSelection(std::string* in=NULL);
  TPRegexp* GetUnSelection(std::string* in=NULL);

  private:
  AliZMQhistViewer(const AliZMQhistViewer&);
  AliZMQhistViewer& operator=(const AliZMQhistViewer&);

  //internal stuff
  static Int_t CountPads(TVirtualPad* pad);

  //configuration vars
  TCanvas* fCanvas;
  Bool_t fVerbose ;
  TString fZMQconfigIN ; 
  int fZMQsocketModeIN;
  TString fZMQconfigCONFIG;

  TString fFilter ;
  int fPollInterval ;
  int fPollTimeout ;
  Bool_t fSort ;

  //internal state
  void* fZMQcontext ;
  void* fZMQin  ;
  void* fZMQconfig;
  void* fZMQsleeper;
  void* fZMQsleeper2;
  std::deque<std::vector<ZMQviewerObject>*> fTrashQueue;
  std::vector<ZMQviewerObject>* fContent;

  TString fStatus ;
  Int_t fRunNumber ;
  TPRegexp* fSelectionRegexp ;
  TPRegexp* fUnSelectionRegexp ;
  TString fDrawOptions;
  Bool_t fScaleLogX ;
  Bool_t fScaleLogY ;
  Bool_t fScaleLogZ ;
  Bool_t fResetOnRequest ;
  Int_t fHistStats ;

  Bool_t fAllowResetAtSOR ;

  ULong64_t iterations;

  //inter thread exchange data
  std::vector<ZMQviewerObject>* fIncoming;
  std::string fInfo;
  bool fClearCanvas;
  bool fUpdateCanvas;
  bool fTerminated;

  ClassDef(AliZMQhistViewer, 0)
};

struct ZMQviewerObjectTitleComparator {
  bool operator()(const ZMQviewerObject& left, const std::string& right) {
    if (!(left.object)) {return true;}
    return right.compare(left.object->GetTitle())>0;
  }

  bool operator()(const ZMQviewerObject& left, const ZMQviewerObject& right) {
    if (!(left.object)) {return true;}
    if (!(right.object)) {return false;}
    return strcmp(right.object->GetTitle(),left.object->GetTitle())>0;
  }
};

#endif
