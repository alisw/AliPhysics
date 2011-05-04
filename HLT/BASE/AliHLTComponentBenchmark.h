//-*- Mode: C++ -*-
// $Id$
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************


#ifndef ALIHLTCOMPONENTBENCHMARK_H
#define ALIHLTCOMPONENTBENCHMARK_H

#include "TStopwatch.h"
#include "TString.h"

/**
 * @class AliHLTComponentBenchmark
 *
 * AliHLTComponentBenchmark can be used to benchmark HLT compnoents
 */
class AliHLTComponentBenchmark
{
  public:

  AliHLTComponentBenchmark( const char *Name="" );
  ~AliHLTComponentBenchmark(){}

  void Reset();
  void SetName( const char *Name );
  void SetTimer( Int_t i, const char *Name );
  void StartNewEvent();
  void Start( Int_t i );
  void Stop( Int_t i );
  void AddInput( Double_t x );
  void AddOutput( Double_t x );
  const char *GetStatistics();

  private:

  TString fComponentName;// name of the component
  Int_t fNTimers; // n of timers
  TStopwatch fTimers[10]; // the timers
  TString fNames[10]; // timer names
  ULong_t fNEvents; // N events processed
  Double_t fTotalRealTime[10]; // total real time
  Double_t fTotalCPUTime[10]; // total CPU time
  Double_t fTotalInput; // total input size
  Double_t fTotalOutput; // total output size
  TString fStatistics;// string with statistics
};

#endif 
