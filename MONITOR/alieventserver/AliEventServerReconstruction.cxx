// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#include <TEnv.h>
#include <TSystem.h>
#include <TThread.h>
#include <TXMLEngine.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliCDBManager.h>
#include <AliGRPPreprocessor.h>
#include <AliReconstruction.h>
#include <AliTPCRecoParam.h>

#include <iostream>
#include <sstream>

#include "AliEventServerUtil.h"
#include "AliEventServerReconstruction.h"
#include "AliStorageEventManager.h"
#include "AliStorageTypes.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTrackPointArray.h"
#include "AliESDfriendTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliTracker.h"

using namespace std;

AliEventServerReconstruction::AliEventServerReconstruction()
	: TQObject(),
	  fAliReco(0),
	  fCDBmanager(0),
	  fCurrentRunId(0),
	  fIsListenning(kFALSE),
	  fSettings(0),
	  fHost(0),
	  fRecoThread(0)
{}

AliEventServerReconstruction::~AliEventServerReconstruction()
{
	Close();
	if(fSettings){delete fSettings;fSettings=0;}
}

void AliEventServerReconstruction::Close()
{
	if(fIsListenning)
	{
		Info("AliRecoServer::Close", "Closing Server");
		StopReconstruction();
		fIsListenning = kFALSE;
	}
}

Bool_t AliEventServerReconstruction::StartReconstruction(Int_t run, const char* input)
{
	fCurrentRunId = run;

	StopReconstruction();

	// re-read settings
	if(fSettings){delete fSettings;}
	fSettings = new TEnv(AliEventServerUtil::GetPathToServerConf());
	
	TString recoBaseDir = fSettings->GetValue("server.saveRecoDir",
						  DEFAULT_SERVER_SAVE_RECO_DIR);
	
	// Create directories and logfile
	TString logFile = Form("%s/log/run%d.log",
			       recoBaseDir.Data(),
			       run);
	
	Info("DoStart","Reconstruction log will be written to %s",logFile.Data());
	if( gSystem->RedirectOutput(logFile.Data())!=0)
	{
		printf(Form("AliRecoServer::StartReconstruction [] Error while trying to redirect output to [%s]. Exiting...", logFile.Data()) );
		return kFALSE;
	}

	gSystem->cd(recoBaseDir.Data());

	TString gdcs;
	if (RetrieveGRP(run,gdcs) <= 0 || gdcs.IsNull()){return kFALSE;}
	  
	gSystem->mkdir(Form("run%d", run));
	gSystem->cd(Form("run%d", run));

	// Create Reco and Reco Thread
	SetupReco(input);
 	fAliReco->InitRun(input);

	fHost = (const char*)Form("%s:%d", fSettings->GetValue("server.host", DEFAULT_SERVER_HOST), fSettings->GetValue("server.port", DEFAULT_SERVER_PORT));
	
	fRecoThread = new TThread("AliEventServerReconstruction",
                              Dispatch,
                              (void*)this);
	fRecoThread->Run();
	fIsListenning = kTRUE;

	return true;
}

void AliEventServerReconstruction::StopReconstruction()
{
	if(!fRecoThread) return;
	delete fRecoThread;
	fRecoThread=0;
	Emit("Stopped()");
  
 	if(fAliReco){delete fAliReco;fAliReco=0;}
 	if(fCDBmanager){delete fCDBmanager;fCDBmanager=0;}
}

void AliEventServerReconstruction::ReconstructionHandle()
{
	TThread::SetCancelAsynchronous();
	TThread::SetCancelOn();
	
	if(!fAliReco) return;
	
	AliESDEvent* event;
	
	AliStorageEventManager *eventManager = AliStorageEventManager::GetEventManagerInstance();
	eventManager->CreateSocket(EVENTS_SERVER_PUB);
	eventManager->CreateSocket(XML_PUB);

	fAliReco->Begin(NULL);
	if (fAliReco->GetAbort() != TSelector::kContinue) return;
	fAliReco->SlaveBegin(NULL);
	if (fAliReco->GetAbort() != TSelector::kContinue) return;

	//******* The loop over events
	Int_t iEvent = 0;
	while (fAliReco->HasNextEventAfter(iEvent))
	{
		// check if process has enough resources 
		if (!fAliReco->HasEnoughResources(iEvent)) break;
		Bool_t status = fAliReco->ProcessEvent(iEvent);
      
		if (status)
		{
			event = fAliReco->GetESDEvent();
			eventManager->Send(event,EVENTS_SERVER_PUB);
			eventManager->SendAsXml(event,XML_PUB);
		}
		else
		{
			fAliReco->Abort("ProcessEvent",TSelector::kAbortFile);
		}
      		
		fAliReco->CleanProcessedEvent();
		iEvent++;
	}
	fAliReco->SlaveTerminate();
	if (fAliReco->GetAbort() != TSelector::kContinue) return;
	fAliReco->Terminate();
	if (fAliReco->GetAbort() != TSelector::kContinue) return;  
}

Int_t AliEventServerReconstruction::RetrieveGRP(UInt_t run, TString &gdc)
{
	if(!fSettings) return (-1);

	// Retrieve GRP entry for given run from aldaqdb.
	TString dbHost = fSettings->GetValue("logbook.host", DEFAULT_LOGBOOK_HOST);
	Int_t dbPort =  fSettings->GetValue("logbook.port", DEFAULT_LOGBOOK_PORT);
	TString dbName =  fSettings->GetValue("logbook.db", DEFAULT_LOGBOOK_DB);
	TString user =  fSettings->GetValue("logbook.user", DEFAULT_LOGBOOK_USER);
	TString password = fSettings->GetValue("logbook.pass", DEFAULT_LOGBOOK_PASS);
	
	Int_t ret=AliGRPPreprocessor::ReceivePromptRecoParameters(run, dbHost.Data(),
								  dbPort, dbName.Data(),
								  user.Data(), password.Data(),
								  Form("local://%s",gSystem->pwd()),
								  gdc);

	if(ret>0) Info("RetrieveGRP","Last run of the same type is: %d",ret);
	else if(ret==0) Warning("RetrieveGRP","No previous run of the same type found");
	else if(ret<0) Error("Retrieve","Error code while retrieving GRP parameters returned: %d",ret);
	return(ret);
}

void AliEventServerReconstruction::SetupReco(const char* input)
{
	if(!fSettings) return;

	//AliTPCRecoParam::SetUseTimeCalibration(kFALSE); //-- !probably should be set from conf file!

	printf(Form("=========================[local://%s/..]===========\n",gSystem->pwd()));

	/* Settings CDB */
	fCDBmanager = AliCDBManager::Instance();
  
	fCDBmanager->SetDefaultStorage(fSettings->GetValue("cdb.defaultStorage", DEFAULT_CDB_STORAGE));
	fCDBmanager->SetSpecificStorage(fSettings->GetValue( "cdb.specificStoragePath1", DEFAULT_CDB_SPEC_STORAGE_PATH1),  
				    fSettings->GetValue( "cdb.specificStorageValue1", DEFAULT_CDB_SPEC_STORAGE_VALUE1));

	fCDBmanager->SetSpecificStorage(fSettings->GetValue( "cdb.specificStoragePath2", DEFAULT_CDB_SPEC_STORAGE_PATH2),  
				    fSettings->GetValue( "cdb.specificStorageValue2", DEFAULT_CDB_SPEC_STORAGE_VALUE2));
  
	fCDBmanager->SetSpecificStorage(fSettings->GetValue( "cdb.specificStoragePath3", DEFAULT_CDB_SPEC_STORAGE_PATH3),  
				    fSettings->GetValue( "cdb.specificStorageValue3", DEFAULT_CDB_SPEC_STORAGE_VALUE3));
  
	/* Reconstruction settings */
	if(fAliReco) delete fAliReco;
  
	AliReconstruction* rec = new AliReconstruction;
	
	// QA options
	rec->SetRunQA(fSettings->GetValue( "qa.runDetectors", DEFAULT_QA_RUN));
	rec->SetRunGlobalQA(fSettings->GetValue( "qa.runGlobal", DEFAULT_QA_RUN_GLOBAL));
	rec->SetQARefDefaultStorage(fSettings->GetValue( "qa.defaultStorage",DEFAULT_QAREF_STORAGE)) ;
	rec->SetRunPlaneEff(fSettings->GetValue( "reco.runPlaneEff", DEFAULT_RECO_RUN_PLANE_EFF));

	// AliReconstruction settings
	rec->SetWriteESDfriend(fSettings->GetValue( "reco.writeESDfriend", DEFAULT_RECO_WRITE_ESDF));
	rec->SetWriteAlignmentData(fSettings->GetValue( "reco.writeAlignment",DEFAULT_RECO_WRITE_ALIGN));
	rec->SetInput(input); // reconstruct data from this input
	rec->SetRunReconstruction(fSettings->GetValue( "reco.detectors", DEFAULT_RECO_DETECTORS));
	rec->SetUseTrackingErrorsForAlignment("ITS"); //-- !should be set from conf file!

	// switch off cleanESD
	rec->SetCleanESD(fSettings->GetValue( "reco.cleanESD",DEFAULT_RECO_CLEAN_ESD));

	fAliReco = rec;
}
