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
//  $Id$

//_________________________________________________________________________
// Class to get organized with the way we're timing our methods...
//
// Typical usage is based on macros (like for AliLog related ones AliDebug...)
//
// The idea is to instrument the code with a few macro calls, and then,
// at the end of the execution, get a printout of *all* the timers, by using
// AliCodeTimer::Instance()->Print()
// instead of getting scattered outputs all over the place.
//
// To time a given method, use :
//
// void ClassA::MethodA(....)
// {
//    AliCodeTimerAuto("")
// }
//
// To get several timers within a same method, use : 
//
// void ClassA::MethodB(...)
// {
//   AliCodeTimerStart("doing something")
//   ....
//   AliCodeTimerStop("doing something")
//
//   AliCodeTimerStart("doing something else")
//   ....
//   AliCodeTimerStop("doing something else")
// }

#include "AliCodeTimer.h"

#include <TMap.h>
#include <TObjString.h>
#include <TStopwatch.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliCodeTimer)
ClassImp(AliCodeTimer::AliPair)
/// \endcond

AliCodeTimer* AliCodeTimer::fgInstance(0x0);

//_____________________________________________________________________________
void
AliCodeTimer::AliPair::Print(Option_t* opt) const
{
  // Print timer information
  cout << opt << Form("%s R:%.4fs C:%.4fs (%d slices)",
                      Name().Data(),Timer()->RealTime(),
                      Timer()->CpuTime(),Timer()->Counter()-1) << endl;
}


//_____________________________________________________________________________
AliCodeTimer::AliCodeTimer() : TObject(), fTimers(new TMap)
{
  /// Ctor
  fTimers->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliCodeTimer::~AliCodeTimer()
{
  /// Dtor
  Reset();
  delete fTimers;
}

//_____________________________________________________________________________
AliCodeTimer*
AliCodeTimer::Instance()
{
  // single instance of this class
  if (!fgInstance) fgInstance = new AliCodeTimer;
  return fgInstance;
}

//_____________________________________________________________________________
void AliCodeTimer::Continue(const char* classname, const char* methodname, 
                            const char* message)
{
  /// Resume a previously stop timer
  TStopwatch* t = Stopwatch(classname,methodname,message);
  if (t)
  {
    t->Continue();
  }
  else
  {
    AliError(Form("No timer for %s/%s/%s",classname,methodname,message));
  }
}

//_____________________________________________________________________________
Double_t AliCodeTimer::CpuTime(const char* classname, 
                               const char* methodname,
                               const char* message) const
{
  /// Return cpu time for a given timer
  TStopwatch* t = Stopwatch(classname,methodname,message);
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
AliCodeTimer::MethodMap(const char* classname) const
{
  /// Return the map for a given "classname"
  return static_cast<TMap*>(fTimers->GetValue(classname));
}

//_____________________________________________________________________________
TObjArray*
AliCodeTimer::MessageArray(const char* classname, const char* methodname) const
{
  /// Return the array for a given AliPair (classname,methodname)
  TMap* m = MethodMap(classname);
  if ( m ) 
  {
    return static_cast<TObjArray*>(m->GetValue(methodname));
  }
  return 0;
}

//_____________________________________________________________________________
void AliCodeTimer::PrintMethod(const char* classname, const char* methodname) const
{
  /// Print all the timers for a given method
  TObjArray* messages = MessageArray(classname,methodname);
  messages->Sort();
  
  cout << "   " << methodname << " ";
  
  if ( messages->GetLast() == 0 ) 
  {
    AliPair* p = static_cast<AliPair*>(messages->First());
    p->Print();
  }
  else
  {
    cout << endl;
    
    TIter next(messages);
    AliPair* p;
  
    while ( ( p = static_cast<AliPair*>(next()) ) ) 
    {
      p->Print("        ");
    }   
  }
}

//_____________________________________________________________________________
void AliCodeTimer::PrintClass(const char* classname) const
{
  /// Print all the timers for a given class
  TMap* methods = MethodMap(classname);
  TIter next(methods);
  TObjString* methodname;
  TObjArray methodNameArray;
  
  while ( ( methodname = static_cast<TObjString*>(next()) ) ) 
  {
    methodNameArray.Add(methodname);
  }
  
  cout << classname << endl;
  
  methodNameArray.Sort();
  
  TIter mnext(&methodNameArray);
  
  while ( ( methodname = static_cast<TObjString*>(mnext()) ) ) 
  {
    PrintMethod(classname,methodname->String().Data());
  }
}
  
//_____________________________________________________________________________
void AliCodeTimer::Print(Option_t* /*opt*/) const
{
  /// Print all the timers we hold
  TIter next(fTimers);
  TObjString* classname;
  TObjArray classNameArray;
  
  while ( ( classname = static_cast<TObjString*>(next()) ) )
  {
    classNameArray.Add(classname);
  }
  
  classNameArray.Sort();
  
  TIter cnext(&classNameArray);
  while ( ( classname = static_cast<TObjString*>(cnext()) ) )
  {
    PrintClass(classname->String().Data());
  }
}

//_____________________________________________________________________________
Double_t 
AliCodeTimer::RealTime(const char* classname, const char* methodname,
                       const char* message) const
{
  /// Return real time of a given time
  TStopwatch* t = Stopwatch(classname,methodname,message);
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
AliCodeTimer::Reset()
{
  /// Reset
  TIter next(fTimers);
  TObjString* classname;
  
  while ( ( classname = static_cast<TObjString*>(next()) ) ) 
  {
    TMap* m = static_cast<TMap*>(fTimers->GetValue(classname->String().Data()));
    m->DeleteAll();
  }
  
  fTimers->DeleteAll();
}

//_____________________________________________________________________________
void 
AliCodeTimer::Start(const char* classname, const char* methodname,
                    const char* message)
{
  /// Start a given time
  TStopwatch* t = Stopwatch(classname,methodname,message);
  if (!t)
  {
    TMap* m = MethodMap(classname);
    if (!m)
    {
      m = new TMap;
      m->SetOwner(kTRUE);
      fTimers->Add(new TObjString(classname),m);
    }      
    TObjArray* messages = MessageArray(classname,methodname);
    if (!messages)
    {
      messages = new TObjArray;
      messages->SetOwner(kTRUE);
      m->Add(new TObjString(methodname),messages);
    }
    t = new TStopwatch;
    t->Start(kTRUE);
    t->Stop();
    messages->Add(new AliPair(new TObjString(message),t));
  }
  t->Start(kFALSE);
}

//_____________________________________________________________________________
void 
AliCodeTimer::Stop(const char* classname, const char* methodname,
                   const char* message)
{
  /// Stop a given timer
  TStopwatch* t = Stopwatch(classname,methodname,message);
  if (!t)
  {
    AliError(Form("No timer for %s/%s/%s",classname,methodname,message));
  }
  else
  {
    t->Stop();
  }
}

//_____________________________________________________________________________
TStopwatch* 
AliCodeTimer::Stopwatch(const char* classname, const char* methodname,
                        const char* message) const
{
  /// Return the internal TStopwatch for a given timer
  TObjArray* a = MessageArray(classname,methodname);
  if ( a ) 
  {
    if (message)
    {
      TIter next(a);
      AliPair* p;
      while ( ( p = static_cast<AliPair*>(next()) ) ) 
      {
        TString s = p->Name();
        if ( s == TString(message) ) 
        {
          return p->Timer();
        }
      }
    }
    else
    {
      return static_cast<TStopwatch*>(a->First());
    }
  }
  return 0x0;
}
