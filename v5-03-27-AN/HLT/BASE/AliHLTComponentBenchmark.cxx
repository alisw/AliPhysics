// $Id$
// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> 
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************


#include "AliHLTComponentBenchmark.h"

AliHLTComponentBenchmark::AliHLTComponentBenchmark( const char *Name )
  :fComponentName(Name),fNTimers(0),fNEvents(0), fTotalInput(0),fTotalOutput(0), fStatistics()
{
  // !
  Reset();
}

void AliHLTComponentBenchmark::Reset()
{
  // !
  fNEvents = 0;
  fTotalInput = 0;
  fTotalOutput = 0;
  for( int i=0; i<10; i++ ){
    fTimers[i].Reset();
    fTotalRealTime[i] = 0;
    fTotalCPUTime[i] = 0;
  }
  fStatistics = "";
}

void AliHLTComponentBenchmark::SetName( const char *Name )
{
  // !
  fComponentName = Name;
}

void AliHLTComponentBenchmark::SetTimer( Int_t i, const char *Name )
{
  // !
  if( i>=10 ) return;
  if( i>=fNTimers ){
    for( ; fNTimers<=i; fNTimers++ ){      
      fTimers[fNTimers].Reset();
      fTotalRealTime[fNTimers] = 0;
      fTotalCPUTime[fNTimers] = 0;
      fNames[fNTimers] = Form("timer %d",fNTimers);
    }
    fNames[i] = Name;   
  }
}

void AliHLTComponentBenchmark::StartNewEvent()
{
  // !
  fNEvents++;
  for( int i=0; i<10; i++ ){
    fTimers[i].Reset();
  } 
}

void AliHLTComponentBenchmark::Start( Int_t i )
{
  // !
  if( i>=10 ) return;
  fTimers[i].Start();
}

void AliHLTComponentBenchmark::Stop( Int_t i )
{
  // !
  if( i>=10 ) return;
  fTimers[i].Stop();
  fTotalRealTime[i]+= fTimers[i].RealTime();
  fTotalCPUTime[i] += fTimers[i].CpuTime();
}


void AliHLTComponentBenchmark::AddInput( Double_t x )
{
  // !
  fTotalInput+=x;
}

void AliHLTComponentBenchmark::AddOutput( Double_t x )
{
  // !
  fTotalOutput+=x;
}

const char *AliHLTComponentBenchmark::GetStatistics()
{
  // !
  if( fNEvents<=0 ) return fStatistics.Data();
  float ratio = 1;
  if( fTotalInput >0 ) ratio = fTotalOutput / fTotalInput;

  fStatistics = Form("%s, %ld events: in %.1f Kb, out %.1f Kb, ratio %.1f", 
		     fComponentName.Data(), fNEvents, fTotalInput/fNEvents/1024, fTotalOutput/fNEvents/1024, ratio);
  
  if( fNTimers<=0 ) return fStatistics.Data();
  float hz = ( fTotalRealTime[0] > 0 ) ?fNEvents/fTotalRealTime[0] : 0;
  fStatistics+=Form("; Time %.1fms/%.1fHz (real/cpu = ",fTotalRealTime[0]/fNEvents*1.e3,hz);
  
  for( int i=0; i<fNTimers; i++ ){
    if( i>0 ) fStatistics+=", ";
    fStatistics+= Form("%s %.1f/%.1fms",fNames[i].Data(),fTotalRealTime[i]/fNEvents*1.e3, fTotalCPUTime[i]/fNEvents*1.e3  );    
  }
  fStatistics+=")";
  return fStatistics.Data();
}
