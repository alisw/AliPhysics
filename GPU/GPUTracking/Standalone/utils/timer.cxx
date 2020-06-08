//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file timer.cxx
/// \author David Rohr

#include "timer.h"
#ifdef _WIN32
#include <windows.h>
#include <winbase.h>
#elif defined(__MACH__) || defined(__APPLE__)
#include <mach/clock.h>
#include <mach/mach.h>
#else
#include <ctime>
#endif

inline double HighResTimer::GetTime()
{
#ifdef _WIN32
  __int64 istart;
  QueryPerformanceCounter((LARGE_INTEGER*)&istart);
  return ((double)istart);
#elif defined(__MACH__) || defined(__APPLE__)
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  return ((double)mts.tv_sec * 1.0e9 + (double)mts.tv_nsec);
#else
  timespec tv;
  clock_gettime(CLOCK_REALTIME, &tv);
  return ((double)tv.tv_sec * 1.0e9 + (double)tv.tv_nsec);
#endif
}

inline double HighResTimer::GetFrequency()
{
#ifdef _WIN32
  __int64 ifreq;
  QueryPerformanceFrequency((LARGE_INTEGER*)&ifreq);
  return ((double)ifreq);
#else
  return (1.0e9);
#endif
}

void HighResTimer::Start()
{
  StartTime = GetTime();
  running = 1;
}

void HighResTimer::ResetStart()
{
  ElapsedTime = 0;
  Start();
}

void HighResTimer::Stop()
{
  if (running == 0) {
    return;
  }
  running = 0;
  double EndTime = 0;
  EndTime = GetTime();
  ElapsedTime += EndTime - StartTime;
}

void HighResTimer::StopAndStart(HighResTimer& startTimer)
{
  if (running == 0) {
    return;
  }
  running = 0;
  double EndTime = 0;
  EndTime = GetTime();
  ElapsedTime += EndTime - StartTime;
  startTimer.StartTime = EndTime;
  startTimer.running = 1;
}

void HighResTimer::Reset()
{
  ElapsedTime = 0;
  StartTime = 0;
  running = 0;
}

double HighResTimer::GetElapsedTime() { return ElapsedTime / Frequency; }

double HighResTimer::GetCurrentElapsedTime(bool reset)
{
  if (running == 0) {
    double retVal = GetElapsedTime();
    if (reset) {
      Reset();
    }
    return (retVal);
  }
  double CurrentTime = GetTime();
  double retVal = (CurrentTime - StartTime + ElapsedTime) / Frequency;
  if (reset) {
    ElapsedTime = 0;
    Start();
  }
  return (retVal);
}

void HighResTimer::AddTime(double t) { ElapsedTime += t * Frequency; }

double HighResTimer::Frequency = HighResTimer::GetFrequency();
