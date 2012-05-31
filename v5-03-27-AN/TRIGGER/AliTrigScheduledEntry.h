#ifndef ALITRIGSCHEDULEDENTRY_H
#define ALITRIGSCHEDULEDENTRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 04/01/2010

//==============================================================================
//
//   AliTrigScheduledEntry - ABC for scheduled responses of a device that is
//                           able to fire-up single response functions or the 
//                           full device scheduled sequence. The start time is
//                           in arbitrary units and in case it is 0 will not be
//                           considered when ordering by time by schedulers.
// 
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class AliTrigDevice;

class AliTrigScheduledEntry : public TNamed {

protected:
  Int_t             fStartTime;           // Time to fire-up
  AliTrigDevice    *fDevice;              // Device to fire-up
  
private:
  AliTrigScheduledEntry(const AliTrigScheduledEntry &other);
  AliTrigScheduledEntry &operator=(const AliTrigScheduledEntry &other);

public:
  AliTrigScheduledEntry() : TNamed(), fStartTime(0), fDevice(NULL) {}
  AliTrigScheduledEntry(const char *name, AliTrigDevice *device, Int_t start=0);
  virtual ~AliTrigScheduledEntry() {}

  Int_t             GetStartTime() const     {return fStartTime;}
  AliTrigDevice    *GetDevice()    const     {return fDevice;}
  virtual void      FireUp(Int_t time)                                      = 0;
  void              SetStartTime(Int_t time) {fStartTime = time;}

  ClassDef(AliTrigScheduledEntry, 1) // ABC for scheduled responses
};   

//==============================================================================
//
//   AliTrigScheduledResponse - Scheduled device response function. Fires-up a
//                              single response function at a time.
//
//==============================================================================

class AliTrigScheduledResponse : public AliTrigScheduledEntry {

private:
  Int_t             fOutputID;            // Device output to be fired

private:
  AliTrigScheduledResponse(const AliTrigScheduledResponse &other);
  AliTrigScheduledResponse &operator=(const AliTrigScheduledResponse &other);

public:
  AliTrigScheduledResponse() : AliTrigScheduledEntry(), fOutputID(-1) {}
  AliTrigScheduledResponse(const char *name, AliTrigDevice *device, Int_t output, Int_t start=0);
  virtual ~AliTrigScheduledResponse() {}
  
  Int_t             GetOutputID() const      {return fOutputID;}
  virtual void      FireUp(Int_t time);
  
  ClassDef(AliTrigScheduledResponse, 1) // Scheduled response function for a device
};

//==============================================================================
//
//   AliTrigScheduledDevice - Scheduled entry for a full device sequence. Invokes
//                            the device scheduler when firing-up.
//
//==============================================================================

class AliTrigScheduledDevice : public AliTrigScheduledEntry {

private:
  AliTrigScheduledDevice(const AliTrigScheduledDevice &other);
  AliTrigScheduledDevice &operator=(const AliTrigScheduledDevice &other);

public:
  AliTrigScheduledDevice() : AliTrigScheduledEntry() {}
  AliTrigScheduledDevice(const char *name, AliTrigDevice *device, Int_t start=0);
  virtual ~AliTrigScheduledDevice() {}
  
  virtual void      FireUp(Int_t time);
  
  ClassDef(AliTrigScheduledDevice, 1) // Scheduled device replay
};
#endif
