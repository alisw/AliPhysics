#ifndef ALISHUTTLESTATUS_H
#define ALISHUTTLESTATUS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class stores the status of the Shuttle processing for a given run and a given detector
//

#include <TObject.h>

class AliShuttleStatus : public TObject
{
public:
  enum Status {
    kInvalid = 0,
    kStarted,
    kDCSStarted,
    kDCSError,
    kPPStarted,
    kPPError,
    kDone, // final
    kFailed,  // final
    kStoreFailed 
  };

  AliShuttleStatus();
  AliShuttleStatus(const AliShuttleStatus& c);

  ~AliShuttleStatus();

  AliShuttleStatus& operator=(const AliShuttleStatus& c);
  virtual void Copy(TObject& c) const;

  AliShuttleStatus(Status status);

  UInt_t GetTimeStamp() const { return fTimeStamp; }
  void SetTimeStamp(UInt_t timeStamp) { fTimeStamp = timeStamp; }

  Status GetStatus() const { return fStatus; }
  const char* GetStatusName() const { return GetStatusName(fStatus); }
  void SetStatus(Status status);

  Int_t GetCount() const { return fCount; }
  void SetCount(Int_t count) { fCount = count; }
  void IncreaseCount() { fCount++; }

  static const char* GetStatusName(Status status);

protected:
  UInt_t fTimeStamp;    // timestamp of the last change
  Status fStatus;       // status of the processing
  Int_t fCount;         // number of retries

  ClassDef(AliShuttleStatus, 1);
};

#endif
