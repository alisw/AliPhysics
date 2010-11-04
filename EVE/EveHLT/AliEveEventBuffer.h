//-*- Mode: C++ -*-

// $Id$

#ifndef ALIEVEEVENTBUFFER_H
#define ALIEVEEVENTBUFFER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice     
 */

/** @file   AliEveEventBuffer.h
    @author Svein Lindal
    @date
    @brief  Manager for HOMER in aliroot
*/


class TObjArray;
class TObject;
class TTimer;
class TThread;

#include "TTimer.h"

class AliEveEventBuffer : public TObject{

public:
  
  /** default constructor */
  AliEveEventBuffer();

  /** destructor */
  virtual ~AliEveEventBuffer();

  void SetBufferSize(Int_t bs) { fBufferSize = bs;}
  void SetBusy(Bool_t busy) { fBusy = busy;}
  Bool_t GetBusy() { return fBusy;}
  
  //Navigate the event buffer
  // TObject *  NavigateFwd();
  // TObject *  NavigateBack();
  TObject * NextEvent();
  TObject * Back();
  TObject * Fwd();
  
  void StartBufferMonitor();
  void StopBufferMonitor();
  //Needed for Homer buffer
  virtual void ConnectToSource() = 0;
  //Check if more events are needed in buffer
  void MonitorBuffer();

  static void * BufferThread(void * buffer);
  virtual void WriteToFile() {//Do nothing
    ;
  }

  void CreateBufferThread();

  ULong64_t GetEventId() const { return fEventId[fBIndex[kCurrent]]; }

protected:
  
  enum fBufferIndex {
    kCurrent,
    kLast, 
    kTop,
    kSize
  };

  
  Int_t fBufferSize;
  Int_t fPreBuffer;
  Bool_t fBusy;

  //TClonesArray containing the stored events
  TObjArray * fEventBuffer;

  //Pointer to current event
  TObject * fCurrentEvent;

  //Event buffer indexes
  Int_t fBIndex[kSize];
  
  
  //Add event to buffer
  virtual void AddToBuffer(TObject * event);
  //  virtual void AddToBuffer(TObject * event, ULong64_t eventId);
  virtual TObject * GetEventFromSource() = 0;
  virtual ULong64_t GetEventIdFromSource() { return 0;}


  void FetchEvent();
  

  //Calculate buffer index stuff
  Int_t CalculateDifference(Int_t top, Int_t low);
  Int_t CalculatePrevious(Int_t current);
  Int_t CalculateNext(Int_t current);

  void SetEventId(ULong64_t eventId) { fEventId[fBIndex[kCurrent]] = eventId;}
  
  void SetBufferMonStarted(Bool_t started) {fBufferMonStarted = started;}
  Bool_t GetBufferMonStarted () const { return fBufferMonStarted;}

private:

  /** copy constructor prohibited */
  AliEveEventBuffer(const AliEveEventBuffer&);

  /** assignment operator prohibited */
  AliEveEventBuffer& operator=(const AliEveEventBuffer&);

  

  TObject *  GetNextUnSeen();

  void PrintIndeces();
  void PrintBuffer();

  TTimer * fTimer;

  //Current event id
  ULong64_t * fEventId;
  
  Bool_t fBufferMonStarted;



  ClassDef(AliEveEventBuffer, 0); // Manage connections to HLT data-sources.
};

#endif
