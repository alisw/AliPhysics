// @(#) $Id$

#ifndef ALIL3_Stopwatch
#define ALIL3_Stopwatch

#ifndef no_root
#include <TStopwatch.h>
typedef TStopwatch AliHLTStopwatch; 

#else
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include "AliHLTRootTypes.h"

class AliHLTStopwatch
{
 private:
  static clock_t gTicks;
  enum EState { kUndefined, kStopped, kRunning };

  Double_t     fStartRealTime;   //wall clock start time
  Double_t     fStopRealTime;    //wall clock stop time
  Double_t     fStartCpuTime;    //cpu start time
  Double_t     fStopCpuTime;     //cpu stop time
  Double_t     fTotalCpuTime;    //total cpu time
  Double_t     fTotalRealTime;   //total real time
  EState       fState;           //stopwatch state
  Int_t        fCounter;         //number of times the stopwatch was started

 public:
  
  AliHLTStopwatch();
  ~AliHLTStopwatch();
  void        Start(Bool_t reset = kTRUE);
  void        Stop();
  void        Continue();
  Int_t       Counter() const { return fCounter; }
  Double_t    RealTime();
  void        Reset() { ResetCpuTime(); ResetRealTime(); }
  void        ResetCpuTime (Double_t time = 0) { Stop(); fTotalCpuTime  = time; }
  void        ResetRealTime(Double_t time = 0) { Stop(); fTotalRealTime = time; }
  Double_t    CpuTime();
  void        Print(Char_t *opt="") const;
  static Double_t GetRealTime();
  static Double_t GetCPUTime();

  ClassDef(TStopwatch,1)  //A stopwatch which times real and cpu time
};

#endif

typedef AliHLTStopwatch AliL3Stopwatch; // for backward compatibility

#endif
