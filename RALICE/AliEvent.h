#ifndef ALIEVENT_H
#define ALIEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <iomanip.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
#include "TDatime.h"
 
#include "AliVertex.h"
 
class AliEvent : public AliVertex
{
 public:
  AliEvent();                             // Default constructor
  AliEvent(Int_t n);                      // Create an event to hold initially n tracks
  ~AliEvent();                            // Default destructor
  void SetDayTime(TDatime& stamp);        // Set the date and time stamp
  void SetRunNumber(Int_t run);           // Set the run number
  void SetEventNumber(Int_t evt);         // Set the event number
  void Reset();                           // Reset all values
  TDatime GetDayTime();                   // Provide the date and time stamp
  Int_t GetRunNumber();                   // Provide the run number
  Int_t GetEventNumber();                 // Provide the event number
  void HeaderInfo();                      // Print the event header information
  void Info(TString f="car");             // Print the event info within coordinate frame f

 protected:
  TDatime fDaytime;      // The date and time stamp
  Int_t fRun;            // The run number
  Int_t fEvent;          // The event number

 ClassDef(AliEvent,1) // Creation and investigation of an Alice physics event.
};
#endif
