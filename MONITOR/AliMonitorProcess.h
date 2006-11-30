#ifndef ALIMONITORPROCESS_H
#define ALIMONITORPROCESS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class TFile;
class TFolder;
class TGrid;
class TServerSocket;
class TSocket;
class TTimer;
class TTree;

#include <TObjArray.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>

class AliESD;
class AliITSgeom;
class AliRawReader;
class AliRunLoader;
class AliTPCParam;
class AliLevel3;
class AliHLTHough;


class AliMonitorProcess : public TObject {
public:
  AliMonitorProcess(const char* alienHost,
		    const char* alienDir,
		    const char* selection = "ALL",
		    const char* fileNameGalice = "galice.root");
  virtual ~AliMonitorProcess();

  static const char* GetRevision();

  void             Run();
  void             Stop();
  void             Reset();

  void             ProcessFile(const char* fileName);

  UInt_t           GetRunNumber() const {return fRunNumber;};
  UInt_t           GetEventPeriodNumber() const;
  UInt_t           GetEventOrbitNumber() const;
  UInt_t           GetEventBunchNumber() const;

  enum EStatus     {kStopped, kWaiting, kReading, kRecTPC, kRecITS, kRecV0s,
		    kRecHLT, kFilling, kUpdating, kWriting, kResetting, 
		    kConnecting, kBroadcasting};
  EStatus          GetStatus() const 
    {gSystem->ProcessEvents(); return fStatus;};
  Bool_t           WillStop() const {return fStopping;};
  Bool_t           IsStopped() const {return (fStatus == kStopped);};

  Int_t            GetNumberOfEvents() const {return fNEvents;};
  Int_t            GetNumberOfClients() const 
    {return fSockets.GetEntriesFast();};
  TObjArray*       GetListOfClients() {return &fSockets;};
  Int_t            GetNEventsMin() const {return fNEventsMin;};
  void             SetNEventsMin(Int_t nEventsMin) {fNEventsMin = nEventsMin;};
  void             SetWriteHistoList(Bool_t writeHistoList = kTRUE) 
                                         {fWriteHistoList = writeHistoList;};

  static Int_t     GetPort() {return fgkPort;};
  
private:
  AliMonitorProcess(const AliMonitorProcess& process);
  AliMonitorProcess& operator = (const AliMonitorProcess& process);

  Bool_t           CheckForNewFile();
  Bool_t           ProcessFile();
  Int_t            GetNumberOfEvents(const char* fileName) const;
  Bool_t           ReconstructTPC(AliRawReader* rawReader, AliESD* esd);
  Bool_t           ReconstructITS(AliRawReader* rawReader, AliESD* esd);
  Bool_t           ReconstructV0s(AliESD* esd);
  void             CreateHLT(const char* fileName);
  void             CreateHLTHough(const char* fileName);
  Bool_t           ReconstructHLT(Int_t iEvent);
  Bool_t           ReconstructHLTHough(Int_t iEvent);

  Bool_t           WriteHistos();
  void             StartNewRun();

  void             CheckForConnections();
  void             BroadcastHistos(TSocket* toSocket = NULL);
  void             SetStatus(EStatus status);
  Bool_t           IsSelected(const char* name) const
    {return (fSelection.Contains(name) || fSelection.Contains("ALL"));}

  static const Int_t fgkPort;          // port number for client connections

  TString          fSelection;          // selection of monitor histograms
  TGrid*           fGrid;               // pointer to AliEn
  TString          fAlienDir;           // name of alien directory
  AliRunLoader*    fRunLoader;          // the current run loader
  AliTPCParam*     fTPCParam;           // TPC parameters
  AliITSgeom*      fITSgeom;            // ITS parameters
  TString          fLogicalFileName;    // logical AliEn file name
  TString          fFileName;           // physical file name
  AliLevel3*       fHLT;                // the HLT tracker
  AliHLTHough*      fHLTHough;           // the HLT hough transformer

  UInt_t           fRunNumber;          // current run number
  UInt_t           fSubRunNumber;       // current part (=resets per run)
  UInt_t           fEventNumber[2];     // current event number
  Int_t            fNEvents;            // total number of monitored events
  Int_t            fNEventsMin;         // threshold for writing
  Bool_t           fWriteHistoList;     // write event histos or not

  TFolder*         fTopFolder;          // folder with histos
  TObjArray        fMonitors;           // array of monitor objects
  TFile*           fFile;               // file with tree
  TTree*           fTree;               // monitor tree

  TServerSocket*   fServerSocket;       // socket for client connections
  TObjArray        fSockets;            // array of client sockets
  TSocket*         fDisplaySocket;      // socket for an event display

  EStatus          fStatus;             // current status
  Bool_t           fStopping;           // stop of process requested or not

  class AliMonitorInterruptHandler : public TSignalHandler {
  public:
    AliMonitorInterruptHandler(AliMonitorProcess* process);
    virtual ~AliMonitorInterruptHandler() {};
    virtual Bool_t Notify();
  private:
    AliMonitorProcess* fProcess;       // process to notify

    AliMonitorInterruptHandler(const AliMonitorInterruptHandler& handler); // Not implemented
    AliMonitorInterruptHandler& operator = 
      (const AliMonitorInterruptHandler& handler); // Not implemented

  };

  AliMonitorInterruptHandler* fInterruptHandler;  // interrupt handler

  ClassDef(AliMonitorProcess, 0)   // class for performing the monitoring
};
 

#endif









