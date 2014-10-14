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
	  fRecoThread(0),
	  fRecoIsRunning(false),
	  fRecoWasInitialized(false)
{}

AliEventServerReconstruction::~AliEventServerReconstruction()
{
	Close();
	if(fSettings){delete fSettings;fSettings=0;}
	if(fAliReco){delete fAliReco;fAliReco=0;}
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
  cout<<"Start of run:"<<run<<endl;
	fCurrentRunId = run;

	// re-read settings
	if(fSettings){delete fSettings;fSettings=0;}
	fSettings = new TEnv(AliEventServerUtil::GetPathToServerConf());
	
	TString recoBaseDir = fSettings->GetValue("server.saveRecoDir",DEFAULT_SERVER_SAVE_RECO_DIR);
	
	// Create directories and logfile
	TString logFile = Form("%s/log/run%d.log",recoBaseDir.Data(),run);
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
	cout<<"Setup reco will be called"<<endl;
	SetupReco(input);
	
	fHost = (const char*)Form("%s:%d", fSettings->GetValue("server.host", DEFAULT_SERVER_HOST), fSettings->GetValue("server.port", DEFAULT_SERVER_PORT));
	
	cout<<"Creating new thread"<<endl;
	fRecoThread = new TThread("AliEventServerReconstruction",Dispatch, (void*)this);
	fRecoThread->Run();
	fIsListenning = kTRUE;
	fRecoIsRunning=true;
	cout<<"Reco started"<<endl;
	return true;
}

bool AliEventServerReconstruction::StopReconstruction()
{
  cout<<"Reco server -- StopPeconstruction() called"<<endl;
  if(!fRecoIsRunning || !fRecoThread)
    {
      cout<<"Reco is not running. No need to stop it."<<endl;
      return true;
    }
  if(!fRecoWasInitialized)
    {
      cout<<"Reco is under initialization. Wait until it's finished"<<endl;
      
      return false;
    }
  cout<<"killing thread"<<endl;
  fRecoThread->Kill();
  cout<<"thread killed"<<endl;
  delete fRecoThread;
  fRecoThread=0;
  cout<<"Reco server -- thread removed"<<endl;
  // Emit("Stopped()");

  cout<<"Reco server -- terminating reconstruction"<<endl;	
  fAliReco->SlaveTerminate();
  if (fAliReco->GetAbort() != TSelector::kContinue) return false;
  fAliReco->Terminate();
  if (fAliReco->GetAbort() != TSelector::kContinue) return false; 
	
  if(fAliReco){delete fAliReco;fAliReco=0;}
  cout<<"Reco server -- deleting CDBManager"<<endl;
  if(fCDBmanager){fCDBmanager->Destroy();fCDBmanager=0;}
  cout<<"Reco server -- recontruction stopped"<<endl;
  fRecoIsRunning=false;
  fRecoWasInitialized=false;
  return true;
}

void AliEventServerReconstruction::ReconstructionHandle()
{
	TThread::SetCancelAsynchronous();
	TThread::SetCancelOn();
	
	if(!fAliReco) return;

	AliStorageEventManager *eventManager = AliStorageEventManager::GetEventManagerInstance();
	eventManager->CreateSocket(EVENTS_SERVER_PUB);
	eventManager->CreateSocket(XML_PUB);
	cout<<"Sockets created"<<endl;
	fAliReco->Begin(NULL);
	cout<<"Reco began"<<endl;
	if (fAliReco->GetAbort() != TSelector::kContinue) return;
	fAliReco->SlaveBegin(NULL);
	cout<<"Slave began"<<endl;
	if (fAliReco->GetAbort() != TSelector::kContinue) return;

	fRecoWasInitialized=true;
	//******* The loop over events
	Int_t iEvent = 0;
	AliESDEvent* event;
	while (fAliReco->HasNextEventAfter(iEvent))
	{
		// check if process has enough resources 
	  cout<<"Event server -- checking resources"<<endl;
		if (!fAliReco->HasEnoughResources(iEvent)) break;
		cout<<"Event server -- resources checked"<<endl;
		Bool_t status = fAliReco->ProcessEvent(iEvent);
		cout<<"Event server -- event processed"<<endl;
      
		if (status)
		{
			event = fAliReco->GetESDEvent();
			cout<<"Event server -- sending event"<<endl;
			eventManager->Send(event,EVENTS_SERVER_PUB);
			cout<<"Event server -- sending event as xml"<<endl;
			eventManager->SendAsXml(event,XML_PUB);
			cout<<"Event server -- xml sent"<<endl;
		}
		else
		{
		  cout<<"Event server -- aborting"<<endl;
			fAliReco->Abort("ProcessEvent",TSelector::kAbortFile);
		}
      		cout<<"Event server -- cleaning event"<<endl;
		fAliReco->CleanProcessedEvent();
		cout<<"Event server -- event cleaned"<<endl;
		iEvent++;
	}
	StopReconstruction();
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
	if(!fAliReco)fAliReco = new AliReconstruction();

	// QA options
	//rec->SetRunQA(fSettings->GetValue( "qa.runDetectors", DEFAULT_QA_RUN));
	//rec->SetRunGlobalQA(fSettings->GetValue( "qa.runGlobal", DEFAULT_QA_RUN_GLOBAL));
	fAliReco->SetQARefDefaultStorage(fSettings->GetValue( "qa.defaultStorage",DEFAULT_QAREF_STORAGE)) ;
	//rec->SetRunPlaneEff(fSettings->GetValue( "reco.runPlaneEff", DEFAULT_RECO_RUN_PLANE_EFF));

	fAliReco->SetRunQA(":");
	fAliReco->SetRunGlobalQA(false);
	fAliReco->SetRunPlaneEff(false);

	// AliReconstruction settings
	fAliReco->SetWriteESDfriend(fSettings->GetValue( "reco.writeESDfriend", DEFAULT_RECO_WRITE_ESDF));
	fAliReco->SetWriteAlignmentData(fSettings->GetValue( "reco.writeAlignment",DEFAULT_RECO_WRITE_ALIGN));
	fAliReco->SetInput(input); // reconstruct data from this input
	fAliReco->SetRunReconstruction(fSettings->GetValue( "reco.detectors", DEFAULT_RECO_DETECTORS));
	fAliReco->SetUseTrackingErrorsForAlignment("ITS"); //-- !should be set from conf file!

	// switch off cleanESD
	fAliReco->SetCleanESD(fSettings->GetValue( "reco.cleanESD",DEFAULT_RECO_CLEAN_ESD));

	// init reco for given run
 	fAliReco->InitRun(input);
}
