// @(#) $Id$

// Author: C. Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3Stopwatch.h"

#ifdef no_root

#include "AliL3StandardIncludes.h"

#if __GNUC__ == 3
using namespace std;
#endif


/** \class AliL3Stopwatch
<pre>
//----------------------------------------------------
// AliL3Stopwatch
//
// Stopwatch class. This class returns the real and cpu time between   
// the start and stop events (taken from Root)
</pre>
*/


ClassImp(AliL3Stopwatch)

clock_t AliL3Stopwatch::gTicks = 0;

AliL3Stopwatch::AliL3Stopwatch()
{
  // Create a stopwatch and start it.
  if (!gTicks) gTicks = (clock_t)sysconf(_SC_CLK_TCK);
  fState         = kUndefined;
  fTotalCpuTime  = 0;
  fTotalRealTime = 0;
  fCounter       = 0;
  Start();
}

AliL3Stopwatch::~AliL3Stopwatch()
{
}

void AliL3Stopwatch::Start(Bool_t reset)
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

void AliL3Stopwatch::Stop()
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

void AliL3Stopwatch::Continue()
{
  // Resume a stopped stopwatch. The stopwatch continues counting from the last
  // Start() onwards (this is like the laptimer function).
  
  if (fState == kUndefined){
    cerr << "Error in AliL3Stopwatch::Continue! Stopwatch not started." << endl;
    return;
  }

  if (fState == kStopped) {
    fTotalCpuTime  -= fStopCpuTime  - fStartCpuTime;
    fTotalRealTime -= fStopRealTime - fStartRealTime;
  }

  fState = kRunning;
}

Double_t AliL3Stopwatch::RealTime()
{
  // Return the realtime passed between the start and stop events. If the
  // stopwatch was still running stop it first.
  
  if (fState == kUndefined){
    cerr << "Error in AliL3Stopwatch::RealTime! Stopwatch not started." << endl;
    return -1.0;
  }

  if (fState == kRunning)
    Stop();

  return fTotalRealTime;
}

Double_t AliL3Stopwatch::CpuTime()
{
  // Return the cputime passed between the start and stop events. If the
  // stopwatch was still running stop it first.
  
  if (fState == kUndefined){
    cerr << "Error in AliL3Stopwatch::CpuTime! Stopwatch not started." << endl;
    return -1.0;
  }
  if (fState == kRunning)
    Stop();

  return fTotalCpuTime;
}

Double_t AliL3Stopwatch::GetRealTime()
{
  struct tms cpt;
  Double_t trt =  (Double_t)times(&cpt);
  return trt / (Double_t)gTicks;
}

Double_t AliL3Stopwatch::GetCPUTime()
{
   struct tms cpt;
   times(&cpt);
   return (Double_t)(cpt.tms_utime+cpt.tms_stime) / gTicks;
}

//______________________________________________________________________________
void AliL3Stopwatch::Print(Char_t *opt) const
{
   // Print the real and cpu time passed between the start and stop events.
   // and the number of times (slices) this TStopwatch was called
   // (if this number > 1)

   Double_t  realt = ((AliL3Stopwatch*)(this))->RealTime();

   Int_t  hours = Int_t(realt / 3600);
   realt -= hours * 3600;
   Int_t  min   = Int_t(realt / 60);
   realt -= min * 60;
   Int_t  sec   = Int_t(realt);
   Int_t counter = Counter();
   if (counter <= 1 )
      printf("Real time %d:%d:%d, CP time %.3f", hours, min, sec, ((AliL3Stopwatch*)(this))->CpuTime());
   else
      printf("Real time %d:%d:%d, CP time %.3f, %d slices", hours, min, sec, ((AliL3Stopwatch*)(this))->CpuTime(),counter);
}

#endif
