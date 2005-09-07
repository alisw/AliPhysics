// @(#) $Id$

// Author: C. Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStopwatch.h"
#include "MLUCLog.hpp"

#ifdef no_root

#include "AliHLTTPCStandardIncludes.h"

#if GCCVERSION == 3
using namespace std;
#endif


/** \class AliHLTTPCStopwatch
<pre>
//----------------------------------------------------
// AliHLTTPCStopwatch
//
// Stopwatch class. This class returns the real and cpu time between   
// the start and stop events (taken from Root)
</pre>
*/


ClassImp(AliHLTTPCStopwatch)

clock_t AliHLTTPCStopwatch::gTicks = 0;

AliHLTTPCStopwatch::AliHLTTPCStopwatch()
{
  // Create a stopwatch and start it.
  if (!gTicks) gTicks = (clock_t)sysconf(_SC_CLK_TCK);
  fState         = kUndefined;
  fTotalCpuTime  = 0;
  fTotalRealTime = 0;
  fCounter       = 0;
  Start();
}

AliHLTTPCStopwatch::~AliHLTTPCStopwatch()
{
}

void AliHLTTPCStopwatch::Start(Bool_t reset)
{
  // Start the stopwatch. If reset is kTRUE reset the stopwatch before
  // starting it (including the stopwatch counter).
  // Use kFALSE to continue timing after a Stop() without
  // resetting the stopwatch.

  if (reset) {
    fTotalCpuTime  = 0;
      fTotalRealTime = 0;
      fCounter       = 0;
  }
  if (fState != kRunning) {
    fStartRealTime = GetRealTime();
    fStartCpuTime  = GetCPUTime();
  }
  fState = kRunning;
  fCounter++;
}

void AliHLTTPCStopwatch::Stop()
{
   // Stop the stopwatch.
   fStopRealTime = GetRealTime();
   fStopCpuTime  = GetCPUTime();
   if (fState == kRunning) {
     fTotalCpuTime  += fStopCpuTime  - fStartCpuTime;
     fTotalRealTime += fStopRealTime - fStartRealTime;
   }
   fState = kStopped;
}

void AliHLTTPCStopwatch::Continue()
{
  // Resume a stopped stopwatch. The stopwatch continues counting from the last
  // Start() onwards (this is like the laptimer function).
  
  if (fState == kUndefined){
    cerr << "Error in AliHLTTPCStopwatch::Continue! Stopwatch not started." << endl;
    return;
  }

  if (fState == kStopped) {
    fTotalCpuTime  -= fStopCpuTime  - fStartCpuTime;
    fTotalRealTime -= fStopRealTime - fStartRealTime;
  }

  fState = kRunning;
}

Double_t AliHLTTPCStopwatch::RealTime()
{
  // Return the realtime passed between the start and stop events. If the
  // stopwatch was still running stop it first.
  
  if (fState == kUndefined){
    cerr << "Error in AliHLTTPCStopwatch::RealTime! Stopwatch not started." << endl;
    return -1.0;
  }

  if (fState == kRunning)
    Stop();

  return fTotalRealTime;
}

Double_t AliHLTTPCStopwatch::CpuTime()
{
  // Return the cputime passed between the start and stop events. If the
  // stopwatch was still running stop it first.
  
  if (fState == kUndefined){
    cerr << "Error in AliHLTTPCStopwatch::CpuTime! Stopwatch not started." << endl;
    return -1.0;
  }
  if (fState == kRunning)
    Stop();

  return fTotalCpuTime;
}

Double_t AliHLTTPCStopwatch::GetRealTime()
{
  struct tms cpt;
  Double_t trt =  (Double_t)times(&cpt);
  return trt / (Double_t)gTicks;
}

Double_t AliHLTTPCStopwatch::GetCPUTime()
{
   struct tms cpt;
   times(&cpt);
   return (Double_t)(cpt.tms_utime+cpt.tms_stime) / gTicks;
}

//______________________________________________________________________________
void AliHLTTPCStopwatch::Print(Char_t *opt) const
{
   // Print the real and cpu time passed between the start and stop events.
   // and the number of times (slices) this TStopwatch was called
   // (if this number > 1)

   Double_t  realt = ((AliHLTTPCStopwatch*)(this))->RealTime();

   Int_t  hours = Int_t(realt / 3600);
   realt -= hours * 3600;
   Int_t  min   = Int_t(realt / 60);
   realt -= min * 60;
   Int_t  sec   = Int_t(realt);
   Int_t counter = Counter();
   if (counter <= 1 )
       {
       LOG( MLUCLog::kBenchmark, "AliHLTTPCStopwatch::Print", "Results" )
	   << "Real time " << MLUCLog::kDec << hours << ":" << min << ":" << sec
	   << ", CPU time " << MLUCLog::kPrec << 3 << ((AliHLTTPCStopwatch*)(this))->CpuTime()
	   << ENDLOG;
       }
   //printf("Real time %d:%d:%d, CP time %.3f", hours, min, sec, ((AliHLTTPCStopwatch*)(this))->CpuTime());
   else
       {
       LOG( MLUCLog::kBenchmark, "AliHLTTPCStopwatch::Print", "Results" )
	   << "Real time " << MLUCLog::kDec << hours << ":" << min << ":" << sec
	   << ", CPU time " << MLUCLog::kPrec << 3 << ((AliHLTTPCStopwatch*)(this))->CpuTime()
	   << ", " << counter << " slices" << ENDLOG;
       //printf("Real time %d:%d:%d, CP time %.3f, %d slices", hours, min, sec, ((AliHLTTPCStopwatch*)(this))->CpuTime(),counter);
       }
}

#endif
