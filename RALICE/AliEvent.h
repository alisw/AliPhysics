#ifndef ALIEVENT_H
#define ALIEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliEvent.h,v 1.1 2001/06/06 13:22:48 nick Exp $

#include <iomanip.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
#include "TDatime.h"
 
#include "AliVertex.h"
#include "AliCalorimeter.h"
 
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
  void SetCalCopy(Int_t j);               // (De)activate creation of private copies in fCalorimeters
  Int_t GetCalCopy();                     // Provide CalCopy flag value      
  void AddCalorimeter(AliCalorimeter& c); // Add a calorimeter system to the event
  void AddCalorimeter(AliCalorimeter* c) { AddCalorimeter(*c); }
  Int_t GetNcalorimeters();               // Provide the number of calorimeter systems
  AliCalorimeter* GetCalorimeter(Int_t i);// Provide i-th calorimeter system of the event
  AliCalorimeter* GetCalorimeter(TString name); // Provide calorimeter with name "name"

 protected:
  TDatime fDaytime;         // The date and time stamp
  Int_t fRun;               // The run number
  Int_t fEvent;             // The event number
  Int_t fNcals;             // The number of calorimeter systems 
  TObjArray* fCalorimeters; // Array to hold the pointers to the calorimeter systems
  Int_t fCalCopy;           // Flag to denote creation of private copies in fCalorimeters

 ClassDef(AliEvent,1) // Creation and investigation of an Alice physics event.
};
#endif
