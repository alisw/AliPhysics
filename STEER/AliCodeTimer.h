#ifndef ALICODETIMER_H
#define ALICODETIMER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

///
/// A class to organize TStopwatch timers used to time our code
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TString
#  include "TString.h"
#endif
#ifndef ROOT_TObjString
#  include "TObjString.h"
#endif
#ifndef ALILOG_H
#  include "AliLog.h"
#endif

class TStopwatch;
class TMap;

class AliCodeTimer : public TObject
{
public:
  
  AliCodeTimer();
  virtual ~AliCodeTimer();

  /// Unique instance of this class, which is a singleton
  static AliCodeTimer* Instance();
  
  /// Continue timer(classname,methodname,message)
  void Continue(const char* classname, const char* methodname, const char* message="");

  /// Return the cpu time spent in timer(classname,methodname,message)
  Double_t CpuTime(const char* classname, const char* methodname, const char* message="") const;
  
  /// Print the list of timers we manage
  void Print(Option_t* opt="") const;
  
  /// Return the real time spent in timer(classname,methodname,message)
  Double_t RealTime(const char* classname, const char* methodname, const char* message="") const;
  
  /// Reset all our timers
  void Reset();
  
  /// Start timer(classname,methodname,message)
  void Start(const char* classname, const char* methodname, const char* message="");

  /// Stop timer(classname,methodname,message)
  void Stop(const char* classname, const char* methodname, const char* message="");
    
public:
  
  class AliPair : public TObject
  {
  public:
    AliPair() : TObject(),fName(0), fTimer(0) {}
    // ctor
    AliPair(TObjString* name, TStopwatch* timer) : TObject(), fName(name), fTimer(timer) {}
    virtual ~AliPair() { delete fName; }
    
    /// get name
    TString Name() const { return fName->String(); }
    /// get timer
    TStopwatch* Timer() const { return fTimer; }
    
    /// we are sortable (by name)
    virtual Bool_t IsSortable() const { return kTRUE; }
    /// compare the names
    virtual Int_t Compare(const TObject* object) const
    { return fName->Compare(((const AliPair*)(object))->fName); }

    virtual void Print(Option_t* opt="") const;

private:
    AliPair(const AliPair&);
    AliPair& operator=(const AliPair&);
    
    TObjString* fName; // name of the timer
    TStopwatch* fTimer; // actual timer
    
    ClassDef(AliPair,1) // internal class to hold (string,TStopwatch*) AliPair
  };
    
  class AliAutoPtr
  {
    public:
      
    /// ctor
      AliAutoPtr(const char* classname, const char* methodname, const char* message="") 
      : fA(classname), fB(methodname), fC(message)
      { AliCodeTimer::Instance()->Start(classname,methodname,message); } 

    /// dtor
      ~AliAutoPtr() { AliCodeTimer::Instance()->Stop(fA.Data(),fB.Data(),fC.Data()); }
    
    private:
      TString fA; // first id
      TString fB; // second id
      TString fC; // third id
  };
  
private:  
  
  TMap* MethodMap(const char* classname) const;
  TObjArray* MessageArray(const char* classname, const char* methodname) const;
  TStopwatch* Stopwatch(const char* classname, const char* methodname, const char* message="") const;
  void PrintClass(const char* classname) const;
  void PrintMethod(const char* classname, const char* methodname) const;
  
private:

  AliCodeTimer(const AliCodeTimer& rhs);
  AliCodeTimer& operator=(const AliCodeTimer& rhs);
  
  static AliCodeTimer* fgInstance; //< unique instance
  
  TMap* fTimers; //< internal timers
  
  ClassDef(AliCodeTimer,1) // A timer holder
};

#ifndef LOG_NO_DEBUG

#define AliCodeTimerStartClass(message) AliCodeTimer::Instance()->Start(Class()->GetName(),FUNCTIONNAME(),message);
#define AliCodeTimerStopClass(message) AliCodeTimer::Instance()->Stop(Class()->GetName(),FUNCTIONNAME(),message);
#define AliCodeTimerAutoClass(message) AliCodeTimer::AliAutoPtr aliCodeTimerAliAutoPtrVariable(Class()->GetName(),FUNCTIONNAME(),message);

#define AliCodeTimerStart(message) AliCodeTimer::Instance()->Start(ClassName(),FUNCTIONNAME(),message);
#define AliCodeTimerStop(message) AliCodeTimer::Instance()->Stop(ClassName(),FUNCTIONNAME(),message);
#define AliCodeTimerAuto(message) AliCodeTimer::AliAutoPtr aliCodeTimerAliAutoPtrVariable(ClassName(),FUNCTIONNAME(),message);

#define AliCodeTimerStartGeneral(message) AliCodeTimer::Instance()->Start("General",FUNCTIONNAME(),message);
#define AliCodeTimerStopGeneral(message) AliCodeTimer::Instance()->Stop("General",FUNCTIONNAME(),message);
#define AliCodeTimerAutoGeneral(message) AliCodeTimer::AliAutoPtr aliCodeTimerAliAutoPtrVariable("General",FUNCTIONNAME(),message);

#else

#define AliCodeTimerStartClass(message)
#define AliCodeTimerStopClass(message) 
#define AliCodeTimerAutoClass(message) 

#define AliCodeTimerStart(message) 
#define AliCodeTimerStop(message) 
#define AliCodeTimerAuto(message) 

#define AliCodeTimerStartGeneral(message) 
#define AliCodeTimerStopGeneral(message) 
#define AliCodeTimerAutoGeneral(message) 

#endif

#endif
