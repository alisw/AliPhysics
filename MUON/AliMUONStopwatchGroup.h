#ifndef ALIMUONSTOPWATCHGROUP_H
#define ALIMUONSTOPWATCHGROUP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONStopwatchGroup
/// \brief A class to group timers by name
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TStopwatch;
class TMap;

class AliMUONStopwatchGroup : public TObject
{
public:
  
  AliMUONStopwatchGroup();
  AliMUONStopwatchGroup(const AliMUONStopwatchGroup& rhs);
  AliMUONStopwatchGroup& operator=(const AliMUONStopwatchGroup& rhs);

  virtual ~AliMUONStopwatchGroup();

  void Continue(const char* detector, const char* method);

  Double_t CpuTime(const char* detector, const char* method) const;
  
  void Print(Option_t* opt="") const;
  
  Double_t RealTime(const char* detector, const char* method) const;
  
  void Reset();
  
  void Start(const char* detector, const char* method);

  void Stop(const char* detector, const char* method);
    
public:  
  
  TMap* Map(const char* detector) const;
  
  TStopwatch* Stopwatch(const char* detector, const char* method) const;
  
  void CopyTo(AliMUONStopwatchGroup& timers) const;
  
private:
    
  TMap* fTimers; //< internal timers (map from TObjString to TStopwatch*)
  
  ClassDef(AliMUONStopwatchGroup,1) // A timer holder
};


#endif
