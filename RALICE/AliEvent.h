#ifndef ALIEVENT_H
#define ALIEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id: AliEvent.h,v 1.8 2003/02/25 12:36:28 nick Exp $

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
  virtual ~AliEvent();                    // Default destructor
  AliEvent(AliEvent& evt);                // Copy constructor
  virtual void SetOwner(Bool_t own=kTRUE);// Set ownership of all added objects
  void SetDayTime(TDatime& stamp);        // Set the date and time stamp
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
  TDatime GetDayTime();                   // Provide the date and time stamp
  Int_t GetRunNumber();                   // Provide the run number
  Int_t GetEventNumber();                 // Provide the event number
  void HeaderData();                      // Print the event header information
  void Data(TString f="car");             // Print the event info within coordinate frame f
  void SetCalCopy(Int_t j);               // (De)activate creation of private copies in fCalorimeters
  Int_t GetCalCopy();                     // Provide CalCopy flag value      
  void AddCalorimeter(AliCalorimeter& c); // Add a calorimeter system to the event
  void AddCalorimeter(AliCalorimeter* c) { AddCalorimeter(*c); }
  Int_t GetNcalorimeters();               // Provide the number of calorimeter systems
  void ShowCalorimeters();                // Provide on overview of the available calorimeter systems
  AliCalorimeter* GetCalorimeter(Int_t i);// Provide i-th calorimeter system of the event
  AliCalorimeter* GetCalorimeter(TString name); // Provide calorimeter with name "name"

 protected:
  TDatime fDaytime;         // The date and time stamp
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
  Int_t fNcals;             // The number of calorimeter systems 
  TObjArray* fCalorimeters; // Array to hold the pointers to the calorimeter systems
  Int_t fCalCopy;           // Flag to denote creation of private copies in fCalorimeters

 ClassDef(AliEvent,7) // Creation and investigation of an Alice physics event.
};
#endif
