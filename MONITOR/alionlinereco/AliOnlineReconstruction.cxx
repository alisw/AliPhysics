// Author:  Mihai Niculesu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *)
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliOnlineReconstruction.h"
#include "AliOnlineReconstructionUtil.h"
#include "AliStorageEventManager.h"

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TTimeStamp.h>

#include <signal.h>
#include <iostream>

using namespace std;

bool gQuit = false;
void GotSignal(int){gQuit = true;}

AliOnlineReconstruction::AliOnlineReconstruction(int run) :
  fRun(run),
  fDataSource(""),
  fSettings(0),
  fAliReco(new AliReconstruction()),
  fCDBmanager(AliCDBManager::Instance())
{
  // make sure that destructor is called when kill signal comes
  struct sigaction sa;
  memset(&sa,0,sizeof(sa));
  sa.sa_handler = GotSignal;
  sigfillset(&sa.sa_mask);
  sigaction(SIGINT,&sa,NULL);

  printf("CDB Lock is %s\n",AliCDBManager::Instance()->GetLock() ? "ON":"OFF");


  fSettings.ReadFile(AliOnlineReconstructionUtil::GetPathToServerConf(), kEnvUser);
  StartOfRun();
  cout<<"after startofrun"<<endl;
}

AliOnlineReconstruction::~AliOnlineReconstruction()
{
  cout<<"AliOnlineReconstruction -- destructor called...";
  if(fAliReco)
    {
      //	fAliReco->SlaveTerminate();
      //	fAliReco->Terminate(); 
      //	delete fAliReco;fAliReco=0;
    }
  if(fCDBmanager){fCDBmanager->Destroy();fCDBmanager=0;}
  cout<<"OK"<<endl;
}

void AliOnlineReconstruction::StartOfRun()
{
  if(strcmp("local",fSettings.GetValue("data.source", DEFAULT_DATA_SOURCE))==0)
    {
      cout<<"Starting Reco for run "<<fRun<<endl;
      fDataSource = Form("mem://%s/run%d", gSystem->Getenv("ONLINERECO_RAWFILES_DIR"), fRun);
    }
  else if(strcmp(fSettings.GetValue("data.source", DEFAULT_DATA_SOURCE),"run")==0)
    {
      cout<<"Starting Reco for GDCs active in current run:"<<fRun<<endl;
      fDataSource = fSettings.GetValue("data.online.source", DEFAULT_DATA_ONLINE_SOURCE);
    }
  else{cout<<"\n\nWrong data source. Quitting\n\n"<<endl;}

  TString recoBaseDir = fSettings.GetValue("server.saveRecoDir",DEFAULT_SERVER_SAVE_RECO_DIR);
  cout<<"Reco base dir:"<<recoBaseDir<<endl;

  // Create directories and logfile
  TString logFile = Form("%s/log/run%d.log",recoBaseDir.Data(),fRun);
  Info("DoStart","Reconstruction log will be written to %s",logFile.Data());
  if( gSystem->RedirectOutput(logFile.Data())!=0)
    {
      printf(Form("AliRecoServer::StartReconstruction [] Error while trying to redirect output to [%s]. Exiting...", logFile.Data()) );
      return;
    }
  gSystem->cd(recoBaseDir.Data());

  TString gdcs;
  if (RetrieveGRP(gdcs) <= 0 || gdcs.IsNull()){return;}

  gSystem->Exec(Form("rm -fr run%d;mkdir run%d;cd run%d",fRun,fRun,fRun));

  SetupReco();
  ReconstructionLoop();
}

int AliOnlineReconstruction::RetrieveGRP(TString &gdc)
{
	// Retrieve GRP entry for given run from aldaqdb.
	TString dbHost = fSettings.GetValue("logbook.host", DEFAULT_LOGBOOK_HOST);
	Int_t   dbPort =  fSettings.GetValue("logbook.port", DEFAULT_LOGBOOK_PORT);
	TString dbName =  fSettings.GetValue("logbook.db", DEFAULT_LOGBOOK_DB);
	TString user =  fSettings.GetValue("logbook.user", DEFAULT_LOGBOOK_USER);
	TString password = fSettings.GetValue("logbook.pass", DEFAULT_LOGBOOK_PASS);

	Int_t ret=AliGRPPreprocessor::ReceivePromptRecoParameters(fRun, dbHost.Data(),
								  dbPort, dbName.Data(),
								  user.Data(), password.Data(),
								  Form("local://%s",gSystem->pwd()),
								  gdc);

	if(ret>0) Info("RetrieveGRP","Last run of the same type is: %d",ret);
	else if(ret==0) Warning("RetrieveGRP","No previous run of the same type found");
	else if(ret<0) Error("Retrieve","Error code while retrieving GRP parameters returned: %d",ret);
	return(ret);
}

void AliOnlineReconstruction::SetupReco()
{
	printf(Form("=========================[local://%s/..]===========\n",gSystem->pwd()));

	/* Settings CDB */
	fCDBmanager->SetDefaultStorage(fSettings.GetValue("cdb.defaultStorage", DEFAULT_CDB_STORAGE));
	fCDBmanager->SetSpecificStorage(fSettings.GetValue( "cdb.specificStoragePath1", DEFAULT_CDB_SPEC_STORAGE_PATH1),  
				    fSettings.GetValue( "cdb.specificStorageValue1", DEFAULT_CDB_SPEC_STORAGE_VALUE1));
	fCDBmanager->SetSpecificStorage(fSettings.GetValue( "cdb.specificStoragePath2", DEFAULT_CDB_SPEC_STORAGE_PATH2),  
				    fSettings.GetValue( "cdb.specificStorageValue2", DEFAULT_CDB_SPEC_STORAGE_VALUE2));
	fCDBmanager->SetSpecificStorage(fSettings.GetValue( "cdb.specificStoragePath3", DEFAULT_CDB_SPEC_STORAGE_PATH3),  
				    fSettings.GetValue( "cdb.specificStorageValue3", DEFAULT_CDB_SPEC_STORAGE_VALUE3));
	/* Reconstruction settings */  

	// QA options
	fAliReco->SetRunQA(fSettings.GetValue("qa.runDetectors",DEFAULT_QA_RUN));
	fAliReco->SetRunGlobalQA(fSettings.GetValue("qa.runGlobal",DEFAULT_QA_RUN_GLOBAL));
	fAliReco->SetQARefDefaultStorage(fSettings.GetValue("qa.defaultStorage",DEFAULT_QAREF_STORAGE)) ;
	fAliReco->SetRunPlaneEff(fSettings.GetValue("reco.runPlaneEff",DEFAULT_RECO_RUN_PLANE_EFF));

	// AliReconstruction settings
	fAliReco->SetWriteESDfriend(fSettings.GetValue( "reco.writeESDfriend",DEFAULT_RECO_WRITE_ESDF));
	fAliReco->SetWriteAlignmentData(fSettings.GetValue( "reco.writeAlignment",DEFAULT_RECO_WRITE_ALIGN));
	fAliReco->SetInput(fDataSource.Data()); // reconstruct data from this input
	fAliReco->SetRunReconstruction(fSettings.GetValue( "reco.detectors", DEFAULT_RECO_DETECTORS));
	fAliReco->SetUseTrackingErrorsForAlignment("ITS"); //-- !should be set from conf file!
	fAliReco->SetCleanESD(fSettings.GetValue( "reco.cleanESD",DEFAULT_RECO_CLEAN_ESD));

	// init reco for given run
 	fAliReco->InitRun(fDataSource.Data());
}

void AliOnlineReconstruction::ReconstructionLoop()
{
	AliStorageEventManager *eventManager = AliStorageEventManager::GetEventManagerInstance();
	eventManager->CreateSocket(EVENTS_SERVER_PUB);
	eventManager->CreateSocket(XML_PUB);

	fAliReco->Begin(NULL);
	if (fAliReco->GetAbort() != TSelector::kContinue) return;
	fAliReco->SlaveBegin(NULL);
	if (fAliReco->GetAbort() != TSelector::kContinue) return;
	
	//******* The loop over events
	Int_t iEvent = 0;
	AliESDEvent* event;
	while (fAliReco->HasNextEventAfter(iEvent) && !gQuit)
	{
		if (!fAliReco->HasEnoughResources(iEvent)) break;
		Bool_t status = fAliReco->ProcessEvent(iEvent);
      
		if (status){
			event = fAliReco->GetESDEvent();
			eventManager->Send(event,EVENTS_SERVER_PUB);
			eventManager->SendAsXml(event,XML_PUB);
		}
		else{
		  cout<<"Event server -- aborting"<<endl;
		  fAliReco->Abort("ProcessEvent",TSelector::kAbortFile);
		}
		cout<<"clean"<<endl;
		fAliReco->CleanProcessedEvent();
		cout<<"iEvent++"<<endl;
		iEvent++;
	}
	cout<<"after while"<<endl;
	//	fAliReco->SlaveTerminate();
	//if (fAliReco->GetAbort() != TSelector::kContinue) return;
	//fAliReco->Terminate();
	//if (fAliReco->GetAbort() != TSelector::kContinue) return; 
}
