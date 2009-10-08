// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliDimIntNotifier.h"
#include <TError.h>
#include <signal.h>

//______________________________________________________________________________
// Full description of AliDimIntNotifier
//

ClassImp(AliDimIntNotifier)

Long_t AliDimIntNotifier::fgMainThreadId = 0;

void AliDimIntNotifier::SetMainThreadId()
{
  fgMainThreadId = TThread::SelfId();
}

AliDimIntNotifier::AliDimIntNotifier(const TString& service) :
  DimUpdatedInfo(service, -1),
  fNotifyLck(kTRUE),
  fLastMessage(-1)
{
  fReThreader.Connect("Timeout()", "AliDimIntNotifier", this, "DimMessage()");
}

void AliDimIntNotifier::StartTimer()
{
  fReThreader.Reset();
  fReThreader.TurnOn();
  pthread_kill((pthread_t)fgMainThreadId, SIGALRM); 
}

void AliDimIntNotifier::StopTimer()
{
  fReThreader.TurnOff();
}

void AliDimIntNotifier::infoHandler()
{
  // Handle DIM message

  fNotifyLck.Lock();
  fLastMessage = getData() ? getInt() : -1;
  if (TThread::SelfId() != fgMainThreadId)
  {
    StartTimer();
  }
  else
  {
    ::Warning("DIMinfoHandler", "DIM message received from CINT thread.");
    DimMessage();
  }
  fNotifyLck.UnLock();
}

void AliDimIntNotifier::infoHandlerTest(Int_t fake)
{
  // Fake handler for testing.

  fNotifyLck.Lock();
  fLastMessage = fake;
  if (TThread::SelfId() != fgMainThreadId)
  {
    StartTimer();
  }
  else
  {
    Warning("infoHandlerTest", "Was called from CINT thread ...");
    DimMessage();
  }
  fNotifyLck.UnLock();
}

void AliDimIntNotifier::DimMessage(Int_t)
{
  StopTimer();
  if (fLastMessage != -1)
  {
    Emit("DimMessage(Int_t)", fLastMessage);
  }
}
