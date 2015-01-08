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

#include "AliOnlineRecoTrigger.h"

#include <TSystem.h>

#include "AliLog.h"
#ifdef ALI_DIM
#include "SORNotifier.h"
#endif

ClassImp(TerminateSignalHandler)
ClassImp(AliOnlineRecoTrigger)

//______________________________________________________________________________________________
Bool_t TerminateSignalHandler::Notify()
{
// Sentd terminate command to the Shuttle trigger

	AliInfo("Terminate signal received ...");
	fTrigger->Terminate();

	return kTRUE;
}

//______________________________________________________________________________________________
AliOnlineRecoTrigger::AliOnlineRecoTrigger():
	fNotified(kFALSE), fTerminate(kFALSE),
	fMutex(), fCondition(&fMutex),
	fQuitSignalHandler(0),
	fInterruptSignalHandler(0)
{
  // Default constructor
  //

  fQuitSignalHandler = new TerminateSignalHandler(this, kSigQuit);
  fInterruptSignalHandler = new TerminateSignalHandler(this, kSigInterrupt);

  gSystem->AddSignalHandler(fQuitSignalHandler);
  gSystem->AddSignalHandler(fInterruptSignalHandler);
}

//______________________________________________________________________________________________
AliOnlineRecoTrigger::~AliOnlineRecoTrigger() 
{
  // destructor

  gSystem->RemoveSignalHandler(fQuitSignalHandler);
  gSystem->RemoveSignalHandler(fInterruptSignalHandler);

  delete fQuitSignalHandler;
  fQuitSignalHandler = 0;

  delete fInterruptSignalHandler;
  fInterruptSignalHandler = 0;
}

//______________________________________________________________________________________________
Bool_t AliOnlineRecoTrigger::Notify() {
	//
	// The method is called automaticly by SORNotifier on "start of run" 
	// notification event from ECS
	//

	fMutex.Lock();

	fNotified = kTRUE;
	fCondition.Signal();

	fMutex.UnLock();

	return kTRUE;
}

//______________________________________________________________________________________________
void AliOnlineRecoTrigger::Terminate() {
	//
	// Stop triggers listen mode and exist from Run()
	// Usually called automaticly by TerminateSignalHandler.
	//

	fTerminate = kTRUE;
	fCondition.Signal();
}

//______________________________________________________________________________________________
Int_t AliOnlineRecoTrigger::Run() {
	//
	// AliOnlineRecoTrigger main loop for asynchronized (listen) mode.
	// It spawns DIM service listener and waits for DAQ "start of run"
	// notification.
	//

	fTerminate = kFALSE;

#ifdef ALI_DIM
	SORNotifier* notifier = new SORNotifier(this);
#endif
	Int_t received=0;
	
	AliInfo("Listening for ECS SOR trigger");
	
	fMutex.Lock();

	while (!(fNotified || fTerminate)) {
	  received=fCondition.TimedWaitRelative(1000*86400); //wait 24h
	  if (received==1) break; // 1 = timeout
	}

	fNotified = kFALSE;
		
	fMutex.UnLock();

	Int_t run = -1;
#ifdef ALI_DIM
	run = notifier->GetRun();
	delete notifier;
#endif

	if (fTerminate) {
	  AliInfo("Terminated.");
	  return -1;
	}

	if (received == 0)
	  {
	    AliInfo("Trigger from ECS received!");
	    if (run < 0) {
	      AliError("Invalid run number received!");
	      return -1;
	    }
	    return run;
	  } else if (received == 1) {
	    AliWarning("Timeout waiting for trigger!");
	    return -1;
	  } else {
	    AliError("Error receiving trigger from ECS!");
	    return -1;
	  }
}
