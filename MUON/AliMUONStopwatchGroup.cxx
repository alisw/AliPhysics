/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

/// \class AliMUONStopwatchGroup
///
/// A class to group timers by name
/// Typically used to time out some methods, e.g.
///
/// AliMUONStopwatchGroup timers;
///
/// void Class::Method()
/// {
///    timers.Start("Class","Method");
///    ...
///    timers.Stop();
/// }
///
/// and later on :
///
/// timers.Print();
///
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONStopwatchGroup.h"

#include "AliLog.h"
#include <TMap.h>
#include <TObjString.h>
#include <TStopwatch.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONStopwatchGroup)
/// \endcond

//_____________________________________________________________________________
AliMUONStopwatchGroup::AliMUONStopwatchGroup() : TObject(), fTimers(new TMap)
{
  /// Ctor
  fTimers->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONStopwatchGroup::AliMUONStopwatchGroup(const AliMUONStopwatchGroup& other) : TObject(), fTimers(new TMap)
{
  /// Copy ctor
  other.CopyTo(*this);
}

//_____________________________________________________________________________
AliMUONStopwatchGroup& 
AliMUONStopwatchGroup::operator=(const AliMUONStopwatchGroup& other)
{
  /// Assignment
  Reset();
  other.CopyTo(*this);
  return *this;
}

//_____________________________________________________________________________
AliMUONStopwatchGroup::~AliMUONStopwatchGroup()
{
  /// Dtor
  Reset();
  delete fTimers;
}

//_____________________________________________________________________________
void AliMUONStopwatchGroup::Continue(const char* detector, const char* method)
{
  /// Resume a previously stop timer
  TStopwatch* t = Stopwatch(detector,method);
  if (t)
  {
    t->Continue();
  }
  else
  {
    AliError(Form("No timer for %s/%s",detector,method));
  }
}

//_____________________________________________________________________________
void
AliMUONStopwatchGroup::CopyTo(AliMUONStopwatchGroup& other) const
{
  /// Copy this to other
  TIter next(fTimers);
  TObjString* detector;
  
  while ( ( detector = static_cast<TObjString*>(next()) ) )
  {
    TMap* m = static_cast<TMap*>(fTimers->GetValue(detector->String().Data()));
    TMap* otherm = new TMap;
    otherm->SetOwner(kTRUE);
    other.fTimers->Add(new TObjString(detector->String()),otherm);
    TIter next2(m);
    TObjString* method;
    while ( ( method = static_cast<TObjString*>(next2()) ) )
    {
      TStopwatch* timer = static_cast<TStopwatch*>(m->GetValue(method->String().Data()));
      otherm->Add(new TObjString(method->String()),new TStopwatch(*timer));
    }
  }
}

//_____________________________________________________________________________
Double_t AliMUONStopwatchGroup::CpuTime(const char* detector, const char* method) const
{
  /// Return cpu time for a given timer
  TStopwatch* t = Stopwatch(detector,method);
  if (t)
  {
    return t->CpuTime();
  }
  else
  {
    return 0;
  }
}

//_____________________________________________________________________________
TMap*
AliMUONStopwatchGroup::Map(const char* detector) const
{
  /// Return the map for a given "detector"
  return static_cast<TMap*>(fTimers->GetValue(detector));
}

//_____________________________________________________________________________
void AliMUONStopwatchGroup::Print(Option_t* /*opt*/) const
{
  /// Print all the timers we hold
  TIter next(fTimers);
  TObjString* detector;
  
  while ( ( detector = static_cast<TObjString*>(next()) ) )
  {
    cout << detector->String() << endl;
    TMap* m = static_cast<TMap*>(fTimers->GetValue(detector->String().Data()));
    TIter next2(m);
    TObjString* method;
    while ( ( method = static_cast<TObjString*>(next2()) ) )
    {
      TStopwatch* timer = static_cast<TStopwatch*>(m->GetValue(method->String().Data()));
      cout << Form("    %s R:%.2fs C:%.2fs (%d slices)",
                   method->String().Data(),timer->RealTime(),
                   timer->CpuTime(),timer->Counter()-1) << endl;
    }
  }
}

//_____________________________________________________________________________
Double_t 
AliMUONStopwatchGroup::RealTime(const char* detector, const char* method) const
{
  /// Return real time of a given time
  TStopwatch* t = Stopwatch(detector,method);
  if (t)
  {
    return t->RealTime();
  }
  else
  {
    return 0;
  }
}

//_____________________________________________________________________________
void
AliMUONStopwatchGroup::Reset()
{
  /// Reset
  TIter next(fTimers);
  TObjString* detector;
  
  while ( ( detector = static_cast<TObjString*>(next()) ) ) 
  {
    TMap* m = static_cast<TMap*>(fTimers->GetValue(detector->String().Data()));
    m->DeleteAll();
  }
  
  fTimers->DeleteAll();
}

//_____________________________________________________________________________
void 
AliMUONStopwatchGroup::Start(const char* detector, const char* method)
{
  /// Start a given time
  TStopwatch* t = Stopwatch(detector,method);
  if (!t)
  {
    TMap* m = Map(detector);
    if (!m)
    {
      m = new TMap;
      m->SetOwner(kTRUE);
      fTimers->Add(new TObjString(detector),m);
    }      
    t = new TStopwatch;
    t->Start(kTRUE);
    t->Stop();
    m->Add(new TObjString(method),t);
  }
  t->Start(kFALSE);
}

//_____________________________________________________________________________
void 
AliMUONStopwatchGroup::Stop(const char* detector, const char* method)
{
  /// Stop a given timer
  TStopwatch* t = Stopwatch(detector,method);
  if (!t)
  {
    AliError(Form("No timer for %s/%s",detector,method));
  }
  else
  {
    t->Stop();
  }
}

//_____________________________________________________________________________
TStopwatch* 
AliMUONStopwatchGroup::Stopwatch(const char* detector, const char* method) const
{
  /// Return the internal TStopwatch for a given timer
  TMap* m = Map(detector);
  if (m)
  {
    return static_cast<TStopwatch*>(m->GetValue(method));
  }
  else
  {
    return 0x0;
  }
}

