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

/// \file timer.h
/// \author David Rohr

#ifndef QONMODULE_TIMER_H
#define QONMODULE_TIMER_H

class HighResTimer
{
 public:
  HighResTimer() = default;
  ~HighResTimer() = default;
  void Start();
  void Stop();
  void Reset();
  void ResetStart();
  double GetElapsedTime();
  double GetCurrentElapsedTime(bool reset = false);
  void StopAndStart(HighResTimer& startTimer);
  int IsRunning() { return running; }
  void AddTime(double t);

 private:
  double ElapsedTime = 0.;
  double StartTime = 0.;
  int running = 0;

  static double GetFrequency();
  static double GetTime();
#ifndef GPUCODE
  static double Frequency;
#endif
};

#endif
