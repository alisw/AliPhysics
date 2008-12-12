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
 Revision 1.15  2007/12/10 18:29:23  acolla
 Some log added to the listen mode

 Revision 1.14  2007/12/07 19:14:36  acolla
 in AliShuttleTrigger:

 Added automatic collection of new runs on a regular time basis (settable from the configuration)

 in AliShuttleConfig: new members

 - triggerWait: time to wait for DIM trigger (s) before starting automatic collection of new runs
 - mode: run mode (test, prod) -> used to build log folder (logs or logs_PROD)

 in AliShuttle:

 - logs now stored in logs/#RUN/DET_#RUN.log

 Revision 1.13  2006/11/16 16:16:48  jgrosseo
 introducing strict run ordering flag
 removed giving preprocessor name to preprocessor, they have to know their name themselves ;-)

 Revision 1.12  2006/10/20 15:22:59  jgrosseo
 o) Adding time out to the execution of the preprocessors: The Shuttle forks and the parent process monitors the child
 o) Merging Collect, CollectAll, CollectNew function
 o) Removing implementation of empty copy constructors (declaration still there!)

 Revision 1.11  2006/10/02 16:38:39  jgrosseo
 update (alberto):
 fixed memory leaks
 storing of objects that failed to be stored to the grid before
 interfacing of shuttle status table in daq system

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
// 	1) synchronized - Collect()
// 	2) asynchronized - Run() - starts listening for DAQ "end of run"
// 		notification by DIM service.
//

#include "AliShuttleTrigger.h"

#include <TSystem.h>
#include <TObjString.h>

#include "AliLog.h"
#include "AliShuttleConfig.h"
#include "AliShuttle.h"
#include "DATENotifier.h"

#include <fstream>

ClassImp(TerminateSignalHandler)
ClassImp(AliShuttleTrigger)

//______________________________________________________________________________________________
Bool_t TerminateSignalHandler::Notify()
{
// Sentd terminate command to the Shuttle trigger

	AliInfo("Terminate signal received ...");
	fTrigger->Terminate();

	return kTRUE;
}

//______________________________________________________________________________________________
AliShuttleTrigger::AliShuttleTrigger(const AliShuttleConfig* config):
	fConfig(config), fShuttle(NULL),
	fNotified(kFALSE), fTerminate(kFALSE),
	fMutex(), fCondition(&fMutex),
	fQuitSignalHandler(0),
	fInterruptSignalHandler(0),
	fLastMailDiskSpace(0)
{
	//
	// config - pointer to the AliShuttleConfig object which represents
	// the configuration
	// mainStorage - pointer to AliCDBStorage for the undelying CDBStorage
	// localStorage (local) CDB storage to be used if mainStorage is unavailable
	//

	if (!fConfig->IsValid()) AliFatal("********** !!!!! Invalid configuration !!!!! **********");
	UInt_t timeout = fConfig->GetDCSTimeOut();
	Int_t retries = fConfig->GetDCSRetries();
	fShuttle = new AliShuttle(config, timeout, retries);

	fQuitSignalHandler = new TerminateSignalHandler(this, kSigQuit);
	fInterruptSignalHandler = new TerminateSignalHandler(this, kSigInterrupt);

	gSystem->AddSignalHandler(fQuitSignalHandler);
	gSystem->AddSignalHandler(fInterruptSignalHandler);

}

//______________________________________________________________________________________________
AliShuttleTrigger::~AliShuttleTrigger() 
{
  // destructor

	gSystem->RemoveSignalHandler(fQuitSignalHandler);
	gSystem->RemoveSignalHandler(fInterruptSignalHandler);

	delete fShuttle;

  delete fQuitSignalHandler;
  fQuitSignalHandler = 0;

  delete fInterruptSignalHandler;
  fInterruptSignalHandler = 0;
}

//______________________________________________________________________________________________
Bool_t AliShuttleTrigger::Notify() {
	//
	// Trigger Collect() methods in asynchronized (listen) mode.
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
	// notification. Calls Collect() on notification.
	//

	fTerminate = kFALSE;

	DATENotifier* notifier = new DATENotifier(this, "/LOGBOOK/SUBSCRIBE/ECS_EOR");

	Int_t nTry=0; 
	Int_t nMaxTry = fConfig->GetMaxRetries()+1;
	Int_t received=0;
	
	AliInfo("Listening for ECS trigger");
	
	while (1) {
	
		fMutex.Lock();

		while (!(fNotified || fTerminate)) {
			received=fCondition.TimedWaitRelative(1000*fConfig->GetTriggerWait());
			if (received==1) break; // 1 = timeout
		}

		fNotified = kFALSE;
		
		fMutex.UnLock();

		if (fTerminate) {
			AliInfo("Terminated.");
			break;		
		}
		
		if (received == 0)
		{
			AliInfo("Trigger from ECS received!");
		} else if (received == 1) {
			AliInfo(Form("Timeout (%d s) waiting for trigger. "
				"Starting collection of new runs!", 
					fConfig->GetTriggerWait()));
		} else {
			AliInfo("Error receiving trigger from ECS!");
			break;
		}
		
		nTry++;
		AliInfo(Form("Received %d triggers so far", nTry));
		
		if (fConfig->GetRunMode() == AliShuttleConfig::kTest)
		{
			if(nTry>=nMaxTry)
			{
				AliInfo(Form("Collect() ran more than %d times -> Exiting!", 
						nMaxTry));
				break;
			}
		}
	
		Collect();
	}

	delete notifier;
}

//______________________________________________________________________________________________
Bool_t AliShuttleTrigger::Collect(Int_t run)
{
	//
	// this function creates a thread that runs the shuttle
	// then it checks if the shuttle is still running by checking the monitoring functions of the shuttle
	//

	// first checking disk space
	Long_t id = 0;
	Long_t bsize = 0;
	Long_t blocks = 0;
	Long_t bfree = 0;

	gSystem->GetFsInfo(fConfig->GetShuttleFileSystem(), &id, &bsize, &blocks, &bfree);

	AliInfo(Form("n. of free blocks = %d, total n. of blocks = %d",bfree,blocks));
	Int_t spaceFree = (Int_t)(((Float_t)bfree/(Float_t)blocks)*100);

	if (spaceFree < fConfig->GetFreeDiskWarningThreshold()) {
		AliWarning(Form("************** Free space left = %d%%, below the Warning Threshold (%d%%)",spaceFree,fConfig->GetFreeDiskWarningThreshold()));
		if (TMath::Abs(time(0) - fLastMailDiskSpace) >= 86400){   // 86400 = n. of seconds in 1 d
			SendMailDiskSpace(fConfig->GetFreeDiskWarningThreshold());
			fLastMailDiskSpace = time(0);  // resetting fLastMailDiskSpace to time(0) = now
		}
		if (spaceFree < fConfig->GetFreeDiskFatalThreshold()){
			AliError(Form("*************** Free space left = %d%%, below the Fatal Threshold (%d%%), terminating....",spaceFree,fConfig->GetFreeDiskFatalThreshold()));
			SendMailDiskSpace(fConfig->GetFreeDiskFatalThreshold());
			fTerminate = kTRUE; // terminating....
		}
	}	

	if (fTerminate) {
		return kFALSE;
	}

	return fShuttle->Collect(run);
}
//______________________________________________________________________________________________
Bool_t AliShuttleTrigger::SendMailDiskSpace(Short_t percentage)
{
	//
	// sends a mail to the shuttle experts in case of free disk space < theshold
	//
	
		
	AliInfo("******************* Sending the Mail!! *********************");
	if (!fConfig->SendMail()) 
		return kTRUE;

	Int_t runMode = (Int_t)fConfig->GetRunMode();
	TString tmpStr;
	if (runMode == 0) tmpStr = " Nightly Test:";
	else tmpStr = " Data Taking:"; 
	void* dir = gSystem->OpenDirectory(fShuttle->GetShuttleLogDir());
	if (dir == NULL)
	{
		if (gSystem->mkdir(fShuttle->GetShuttleLogDir(), kTRUE))
		{
			AliWarning(Form("SendMail - Can't open directory <%s>", fShuttle->GetShuttleLogDir()));
			return kFALSE;
		}

	} else {
		gSystem->FreeDirectory(dir);
	}

	// SHUTTLE responsibles in to
	TString to="";
	TIter iterAdmins(fConfig->GetAdmins(AliShuttleConfig::kGlobal));
	TObjString *anAdmin=0;
	while ((anAdmin = (TObjString*) iterAdmins.Next()))
	{
		to += Form("%s,", anAdmin->GetName());
	}
	if (to.Length() > 0)
	  to.Remove(to.Length()-1);
	AliDebug(2, Form("to: %s",to.Data()));

	// mail body 
  	TString bodyFileName;
  	bodyFileName.Form("%s/mail.body", fShuttle->GetShuttleLogDir());
  	gSystem->ExpandPathName(bodyFileName);

  	ofstream mailBody;
  	mailBody.open(bodyFileName, ofstream::out);

  	if (!mailBody.is_open())
	{
    		AliWarning(Form("Could not open mail body file %s", bodyFileName.Data()));
    		return kFALSE;
  	}

	TString subject;
	TString body;

	subject = Form("%s CRITICAL Disk Space usage exceeds %d%c!",
		       tmpStr.Data(),percentage,'%');
	AliDebug(2, Form("subject: %s", subject.Data()));
	Int_t percentage_used = 100 - percentage; 	

	body = "Dear SHUTTLE experts, \n\n";
	body += "The usage of the disk space on the shuttle machine has overcome \n"; 
	body += Form("the threshold of %d%%. \n \n",percentage_used);
	body += "Please check! \n \n";
	body += "Please do not answer this message directly, it is automatically generated.\n\n";
	body += "Greetings,\n\n \t\t\tthe SHUTTLE\n";

	AliDebug(2, Form("Body : %s", body.Data()));

	mailBody << body.Data();
  	mailBody.close();

	// send mail!
	TString mailCommand = Form("mail -s \"%s\" %s < %s",
						subject.Data(),
						to.Data(),
						bodyFileName.Data());
	AliDebug(2, Form("mail command: %s", mailCommand.Data()));

	Bool_t result = gSystem->Exec(mailCommand.Data());

	return result == 0;
}
