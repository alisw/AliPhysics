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

/*
 $Log$
 Revision 1.10  2006/08/15 10:50:00  jgrosseo
 effc++ corrections (alberto)

 Revision 1.9  2006/08/08 14:19:29  jgrosseo
 Update to shuttle classes (Alberto)

 - Possibility to set the full object's path in the Preprocessor's and
 Shuttle's  Store functions
 - Possibility to extend the object's run validity in the same classes
 ("startValidity" and "validityInfinite" parameters)
 - Implementation of the StoreReferenceData function to store reference
 data in a dedicated CDB storage.

 Revision 1.8  2006/07/21 07:37:20  jgrosseo
 last run is stored after each run

 Revision 1.7  2006/07/20 09:54:40  jgrosseo
 introducing status management: The processing per subdetector is divided into several steps,
 after each step the status is stored on disk. If the system crashes in any of the steps the Shuttle
 can keep track of the number of failures and skips further processing after a certain threshold is
 exceeded. These thresholds can be configured in LDAP.

 Revision 1.6  2006/07/19 10:09:55  jgrosseo
 new configuration, accesst to DAQ FES (Alberto)

 Revision 1.5  2006/07/10 13:01:41  jgrosseo
 enhanced storing of last sucessfully processed run (alberto)

 Revision 1.4  2006/07/04 14:59:57  jgrosseo
 revision of AliDCSValue: Removed wrapper classes, reduced storage size per value by factor 2

 Revision 1.3  2006/06/12 09:11:16  jgrosseo
 coding conventions (Alberto)

 Revision 1.2  2006/06/06 14:26:40  jgrosseo
 o) removed files that were moved to STEER
 o) shuttle updated to follow the new interface (Alberto)

 Revision 1.1  2006/03/07 07:52:34  hristov
 New version (B.Yordanov)

 Revision 1.5  2005/11/21 09:03:48  byordano
 one more print added

 Revision 1.4  2005/11/20 10:12:37  byordano
 comments added to AliShuttleTrigger

 */


// 
// This class is to deal with DAQ LogBook and DAQ "end of run" notification.
// It has severeal two modes:
// 	1) syncrhnized - Collect(), CollectNew() and CollectAll methods
// 	2) asynchronized - Run() - starts listening for DAQ "end of run"
// 		notification by DIM service.
//

#include "AliShuttleTrigger.h"

#include <TSystem.h>
#include "AliLog.h"
#include "AliShuttleConfig.h"
#include "AliShuttle.h"
#include "DATENotifier.h"

ClassImp(TerminateSignalHandler)

//______________________________________________________________________
TerminateSignalHandler::TerminateSignalHandler(const TerminateSignalHandler& /*other*/):
TSignalHandler(), fTrigger()
{
// copy constructor (not implemented)

}

//______________________________________________________________________
TerminateSignalHandler &TerminateSignalHandler::operator=(const TerminateSignalHandler& /*other*/)
{
// assignment operator (not implemented)

return *this;
}

//______________________________________________________________________________________________
Bool_t TerminateSignalHandler::Notify() 
{
// Sentd terminate command to the Shuttle trigger

	AliInfo("Terminate signal received ...");
	fTrigger->Terminate();

	return kTRUE;
}

//______________________________________________________________________________________________
//______________________________________________________________________________________________

ClassImp(AliShuttleTrigger)

//______________________________________________________________________________________________
AliShuttleTrigger::AliShuttleTrigger(const AliShuttleConfig* config,
		UInt_t timeout, Int_t retries):
	fConfig(config), fShuttle(NULL),
	fNotified(kFALSE), fTerminate(kFALSE),
	fMutex(), fCondition(&fMutex),
	fQuitSignalHandler(this, kSigQuit),
	fInterruptSignalHandler(this, kSigInterrupt)
{
	//
	// config - pointer to the AliShuttleConfig object which represents
	// the configuration
	// mainStorage - pointer to AliCDBStorage for the undelying CDBStorage
	// localStorage (local) CDB storage to be used if mainStorage is unavailable
	//

	fShuttle = new AliShuttle(config, timeout, retries);

	gSystem->AddSignalHandler(&fQuitSignalHandler);
	gSystem->AddSignalHandler(&fInterruptSignalHandler);

}



//______________________________________________________________________
AliShuttleTrigger::AliShuttleTrigger(const AliShuttleTrigger& /*other*/):
	TObject(), fConfig(), fShuttle(NULL),
	fNotified(kFALSE), fTerminate(kFALSE),
	fMutex(), fCondition(&fMutex),
	fQuitSignalHandler(this, kSigQuit),
	fInterruptSignalHandler(this, kSigInterrupt)

{
// copy constructor (not implemented)

}

//______________________________________________________________________
AliShuttleTrigger &AliShuttleTrigger::operator=(const AliShuttleTrigger& /*other*/)
{
// assignment operator (not implemented)

return *this;
}





//______________________________________________________________________________________________
AliShuttleTrigger::~AliShuttleTrigger() 
{
// destructor

	gSystem->RemoveSignalHandler(&fQuitSignalHandler);
	gSystem->RemoveSignalHandler(&fInterruptSignalHandler);

	delete fShuttle;
}

//______________________________________________________________________________________________
Bool_t AliShuttleTrigger::Notify() {
	//
	// Trigger CollectNew() methods in asynchronized (listen) mode.
	// Usually called automaticly by DATENotifier on "end of run" 
	// notification event.
	//

	fMutex.Lock();

	fNotified = kTRUE;
	fCondition.Signal();

	fMutex.UnLock();

	return kTRUE;
}

//______________________________________________________________________________________________
void AliShuttleTrigger::Terminate() {
	//
	// Stop triggers listen mode and exist from Run()
	// Usually called automaticly by TerminateSignalHandler.
	//

	fTerminate = kTRUE;
	fCondition.Signal();
}

//______________________________________________________________________________________________
void AliShuttleTrigger::Run() {
	//
	// AliShuttleTrigger main loop for asynchronized (listen) mode.
	// It spawns DIM service listener and waits for DAQ "end of run"
	// notification. Calls CollectNew() on notification.
	//

	fTerminate = kFALSE;

	DATENotifier* notifier = new DATENotifier(this, "/DATE/LOGBOOK/UPDATE");

	while (1) {
	
		fMutex.Lock();

		while (!(fNotified || fTerminate)) {
			fCondition.Wait();
		}

		fNotified = kFALSE;
		
		fMutex.UnLock();

		if (fTerminate) {
			AliInfo("Terminated.");
			break;		
		}
	
		CollectNew();
	}

	delete notifier;
}

//______________________________________________________________________________________________
Bool_t AliShuttleTrigger::Collect(Int_t run)
{
	//
	// Collects conditions data for the given run.
	//

	return fShuttle->Collect(run);
}

//______________________________________________________________________________________________
Bool_t AliShuttleTrigger::CollectNew()
{
	//
	// Collects conditions data for all UNPROCESSED run written to DAQ LogBook.
	//

	return fShuttle->CollectNew();
}

//______________________________________________________________________________________________
Bool_t AliShuttleTrigger::CollectAll()
{
	//
	// Collects conditions data for all runs (even if they're already done!) written in Shuttle LogBook.
	//

	return fShuttle->CollectAll();
}

