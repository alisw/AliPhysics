// $Id$

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice     
 */

/** @file   AliEveEventBuffer.h
    @author Svein Lindal
    @date
    @brief  Event buffer for HOMER
*/


#ifndef ALIEVEEVENTBUFFER_H
#define ALIEVEEVENTBUFFER_H



#include "TObject.h"
#include "TMutex.h"

class TObjArray;
class TTimer;
class TThread;
class TTimer;


class AliEveEventBuffer : public TObject{

public:
  
  /** default constructor */
  AliEveEventBuffer();

  /** destructor */
  virtual ~AliEveEventBuffer();

  void SetBufferSize(Int_t bs) { fBufferSize = bs;}
    
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
  virtual void WriteToFile(Int_t runnumber) = 0; 

  void CreateBufferThread();

  ULong64_t GetEventId() const { return fEventId[fBIndex[kCurrent]]; }
  void SetEventId(ULong64_t eventId) { fEventId[fBIndex[kCurrent]] = eventId;}

  Int_t LockMutex() { return fMutex->TryLock();}
  Int_t UnLockMutex() { return fMutex->UnLock();}


protected:
  
  enum fBufferIndex {
    kCurrent,
    kLast, 
    kTop,
    kSize
  };

  
  Int_t fBufferSize;//Size of event buffer
  Int_t fPreBuffer;//How many events should be prefetched
  TObjArray * fEventBuffer;   //TClonesArray containing the stored events
  TObject * fCurrentEvent;   //Pointer to current event
  Int_t fBIndex[kSize];   //Event buffer indexes

  
  
  //Add event to buffer
  virtual void AddToBuffer(TObject * event);
  virtual TObject * GetEventFromSource() = 0;
  virtual ULong64_t GetEventIdFromSource() { return 0;}


  void FetchEvent();
  

  //Calculate buffer index stuff
  Int_t CalculateDifference(Int_t top, Int_t low) const;
  Int_t CalculatePrevious(Int_t current) const;
  Int_t CalculateNext(Int_t current) const;

  
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

  TTimer * fTimer;//Timer to loop over buffer monitor

  //Current event id
  ULong64_t * fEventId;//Event id
  
  Bool_t fBufferMonStarted;//Has buffer monitor loop started?

  TThread * fThread; //Thread pointer
  TMutex * fMutex;//Mutex

  ClassDef(AliEveEventBuffer, 0); // Manage connections to HLT data-sources.
};

#endif
