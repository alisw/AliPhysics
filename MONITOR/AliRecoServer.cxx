// Author: Mihai Niculescu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#include <TEnv.h>
#include <TSystem.h>
#include <TThread.h>

#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliCDBManager.h>
#include <AliGRPPreprocessor.h>
#include <AliReconstruction.h>
#include <AliTPCRecoParam.h>

#include "AliEventServerUtil.h"
#include "AliRecoServer.h"
#include "AliRecoServerThread.h"

ClassImp(AliRecoServer)
AliRecoServer::AliRecoServer()
  : TQObject(),
  fContext(0),
  fReco(0),
  fCDBman(0),
  fCurrentRunId(0),
  fIsListenning(kFALSE),
  fSettings(0),
  fRecoTh(0)
{
	fContext = new zmq::context_t(1);
}

AliRecoServer::~AliRecoServer()
{
  Close(); // Full Close Server
  delete fSettings;
}

void AliRecoServer::Close()
{
  if(fIsListenning){
  	Info("AliRecoServer::Close", "Closing Server");
  	
  	StopReconstruction();
  	
    delete fContext; fContext=0;
     fIsListenning = kFALSE;
  }

}

const char* AliRecoServer::GetError() const
{}

Bool_t AliRecoServer::IsListenning() const
{
  return fIsListenning;
}

void AliRecoServer::ThreadFinished(Int_t status)
{
	Emit("ThreadFinished(Int_t) ", status); 
}

Bool_t AliRecoServer::StartReconstruction(Int_t run, const char* input)
{
	fCurrentRunId = run;
	
  // stop current reconstruction
  StopReconstruction();

	// re-read settings
	if(fSettings) delete fSettings;
	fSettings = new TEnv(Form("%s/MONITOR/%s", gSystem->Getenv("ALICE_ROOT"), ALIEVENTSERVER_CONF) );
	
	TString recoBaseDir = fSettings->GetValue("server.saveRecoDir", DEFAULT_SERVER_SAVE_RECO_DIR);
	
   // Create directories and logfile
	TString logFile = Form("%s/log/run%d.log",
			       recoBaseDir.Data(),
			       run);
	Info("DoStart","Reconstruction log will be written to %s",logFile.Data());
	if( gSystem->RedirectOutput(logFile.Data())!=0){
		printf(Form("AliRecoServer::StartReconstruction [] Error while trying to redirect output to [%s]. Exiting...", logFile.Data()) );
		return kFALSE;
	}

	gSystem->cd(recoBaseDir.Data());

	TString gdcs;
	if (RetrieveGRP(run,gdcs) <= 0 || gdcs.IsNull()) 
	  return kFALSE;
	  
	gSystem->mkdir(Form("run%d", run));
	gSystem->cd(Form("run%d", run));

  // Create Reco and Reco Thread
  SetupReco(input);
 	fReco->InitRun(input);
 	
	fRecoTh = new AliRecoServerThread(fContext, fReco);
	fRecoTh->Start( Form("%s:%d", fSettings->GetValue("server.host", DEFAULT_SERVER_HOST), fSettings->GetValue("server.port", DEFAULT_SERVER_PORT)) );
	fIsListenning = kTRUE;
}

void AliRecoServer::StopReconstruction()
{
  if(!fRecoTh) return;
  
  fRecoTh->Stop();
  
 	delete fReco; fReco = 0;
 	delete fCDBman; fCDBman = 0;
 	
 	// Emit signal
 	ThreadFinished(0);
}

Int_t AliRecoServer::RetrieveGRP(UInt_t run, TString &gdc)
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

void AliRecoServer::SetupReco(const char* input)
{
	if(!fSettings) return;

  //AliTPCRecoParam::SetUseTimeCalibration(kFALSE); //-- !probably should be set from conf file!

	printf(Form("=========================[local://%s/..]===========\n",gSystem->pwd()));

  /* Settings CDB */
	fCDBman = AliCDBManager::Instance();
  
  fCDBman->SetDefaultStorage(fSettings->GetValue("cdb.defaultStorage", DEFAULT_CDB_STORAGE));
  fCDBman->SetSpecificStorage(fSettings->GetValue( "cdb.specificStoragePath1", DEFAULT_CDB_SPEC_STORAGE_PATH1),  
  																									fSettings->GetValue( "cdb.specificStorageValue1", DEFAULT_CDB_SPEC_STORAGE_VALUE1));

	fCDBman->SetSpecificStorage(fSettings->GetValue( "cdb.specificStoragePath2", DEFAULT_CDB_SPEC_STORAGE_PATH2),  
  																									fSettings->GetValue( "cdb.specificStorageValue2", DEFAULT_CDB_SPEC_STORAGE_VALUE2));
  
	fCDBman->SetSpecificStorage(fSettings->GetValue( "cdb.specificStoragePath3", DEFAULT_CDB_SPEC_STORAGE_PATH3),  
  																									fSettings->GetValue( "cdb.specificStorageValue3", DEFAULT_CDB_SPEC_STORAGE_VALUE3));
  
  /* Reconstruction settings */
  if(fReco) delete fReco;
  
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

	fReco = rec;
}
