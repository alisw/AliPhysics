#ifndef ALIMONITORPROCESS_H
#define ALIMONITORPROCESS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObject.h>
#include <TString.h>
#include <TFolder.h>
#include <TFile.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TServerSocket.h>
#include <TTimer.h>
#include <TGrid.h>
#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliTPCParam.h"
#include "AliITSgeom.h"


class AliMonitorProcess : public TObject {
public:
  AliMonitorProcess(const char* alienDir,
		    const char* fileNameGalice = "galice.root");
  virtual ~AliMonitorProcess();

  static const char* GetRevision();

  void             Run();
  void             Stop();
  void             Reset();

  void             ProcessFile(const char* fileName);

  UInt_t           GetRunNumber() {return fRunNumber;};
  UInt_t           GetEventPeriodNumber();
  UInt_t           GetEventOrbitNumber();
  UInt_t           GetEventBunchNumber();

  enum EStatus     {kStopped, kWaiting, kReading, kRecTPC, kRecITS, kRecV0s,
		    kFilling, kUpdating, kWriting, kResetting, 
		    kConnecting, kBroadcasting};
  EStatus          GetStatus() {return fStatus;};
  Bool_t           WillStop() {return fStopping;};
  Bool_t           IsStopped() {return (fStatus == kStopped);};

  Int_t            GetNumberOfEvents() {return fNEvents;};
  Int_t            GetNumberOfClients() {return fSockets.GetEntriesFast();};
  TObjArray*       GetListOfClients() {return &fSockets;};
  Int_t            GetNEventsMin() {return fNEventsMin;};
  void             SetNEventsMin(Int_t nEventsMin) {fNEventsMin = nEventsMin;};
  void             SetWriteHistoList(Bool_t writeHistoList = kTRUE) 
                                         {fWriteHistoList = writeHistoList;};
  
  static const Int_t kgPort;

private:
  Bool_t           CheckForNewFile();
  Bool_t           ProcessFile();
  Int_t            GetNumberOfEvents(const char* fileName);
  Bool_t           ReconstructTPC(AliRawReader* rawReader);
  Bool_t           ReconstructITS(AliRawReader* rawReader);
  Bool_t           ReconstructV0s();

  Bool_t           WriteHistos();
  void             StartNewRun();

  void             CheckForConnections();
  void             BroadcastHistos();

  TGrid*           fGrid;
  AliRunLoader*    fRunLoader;
  AliTPCParam*     fTPCParam;
  AliITSgeom*      fITSgeom;
  TString          fLogicalFileName;
  TString          fFileName;

  UInt_t           fRunNumber;
  UInt_t           fSubRunNumber;
  UInt_t           fEventNumber[2];
  Int_t            fNEvents;
  Int_t            fNEventsMin;
  Bool_t           fWriteHistoList;

  TFolder*         fTopFolder;
  TObjArray        fMonitors;
  TFile*           fFile;
  TTree*           fTree;

  TServerSocket*   fServerSocket;
  TObjArray        fSockets;
  TSocket*         fDisplaySocket;

  EStatus          fStatus;
  Bool_t           fStopping;

  ClassDef(AliMonitorProcess, 0)   // class for performing the monitoring
};
 

#endif









