#ifndef ALIEVENT_H
#define ALIEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliEvent.h,v 1.11 2003/10/26 14:53:44 nick Exp $

#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
#include "TDatime.h"
#include "TTimeStamp.h"
 
#include "AliVertex.h"
 
class AliEvent : public AliVertex
{
 public:
  AliEvent();                             // Default constructor
  AliEvent(Int_t n);                      // Create an event to hold initially n tracks
  virtual ~AliEvent();                    // Default destructor
  AliEvent(AliEvent& evt);                // Copy constructor
  virtual void SetOwner(Bool_t own=kTRUE);// Set ownership of all added objects
  void SetDayTime(TTimeStamp& stamp);     // Set the date and time stamp exactly as specified (1 ns accuracy)
  void SetDayTime(TDatime& stamp);        // Set date and time stamp interpreted as local time (1 s accuracy)
  void SetRunNumber(Int_t run);           // Set the run number
  void SetEventNumber(Int_t evt);         // Set the event number
  void SetProjectile(Int_t a,Int_t z,Double_t pnuc,Int_t id=0); // Set projectile A, Z, p per nucleon and id
  Int_t GetProjectileA();                 // Provide A value of the projectile
  Int_t GetProjectileZ();                 // Provide Z value of the projectile
  Double_t GetProjectilePnuc();           // Provide the projectile momentum value per nucleon
  Int_t GetProjectileId();                // Provide the user defined particle ID of the projectile
  void SetTarget(Int_t a,Int_t z,Double_t pnuc,Int_t id=0); // Set target A, Z, p per nucleon and id
  Int_t GetTargetA();                     // Provide A value of the target
  Int_t GetTargetZ();                     // Provide Z value of the target
  Double_t GetTargetPnuc();               // Provide the target momentum value per nucleon
  Int_t GetTargetId();                    // Provide the user defined particle ID of the target
  void Reset();                           // Reset all values
  TTimeStamp GetDayTime();                // Provide the date and time stamp
  Int_t GetRunNumber();                   // Provide the run number
  Int_t GetEventNumber();                 // Provide the event number
  virtual void HeaderData();              // Print the event header information
  virtual void Data(TString f="car");     // Print the event info within coordinate frame f
  void SetDevCopy(Int_t j);               // (De)activate creation of private copies of the devices
  Int_t GetDevCopy();                     // Provide DevCopy flag value      
  void AddDevice(TObject& d);             // Add a device to the event
  void AddDevice(TObject* d) { AddDevice(*d); }
  Int_t GetNdevices();                    // Provide the number of devices
  void ShowDevices();                     // Provide on overview of the available devices
  TObject* GetDevice(Int_t i);            // Provide i-th device of the event
  TObject* GetDevice(TString name);       // Provide device with name "name"

 protected:
  TTimeStamp fDaytime;      // The date and time stamp
  Int_t fRun;               // The run number
  Int_t fEvent;             // The event number
  Int_t fAproj;             // The projectile A value
  Int_t fZproj;             // The projectile Z value
  Double_t fPnucProj;       // The projectile momentum per nucleon
  Int_t fIdProj;            // User defined projectile particle ID
  Int_t fAtarg;             // The target A value
  Int_t fZtarg;             // The target Z value
  Double_t fPnucTarg;       // The target momentum per nucleon
  Int_t fIdTarg;            // User defined target particle ID
  TObjArray* fDevices;      // Array to hold the pointers to the various devices
  Int_t fDevCopy;           // Flag to denote creation of private copies of the devices

 ClassDef(AliEvent,10) // Creation and investigation of an Alice physics event.
};
#endif
