// Author:  Mihai Niculesu 2013

/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, all rights reserved. *)
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <TEnv.h>
#include <TSystem.h>

#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include <TTimeStamp.h>
#include <TTimer.h>
#include <iostream>

#include <AliLog.h>
#include <AliReconstruction.h>
#include <AliDimIntNotifier.h>

#include "AliEventServerUtil.h"
#include "AliEventServer.h"
#include "AliEventServerReconstruction.h"

ClassImp(AliEventServer)

using namespace std;

AliEventServer::AliEventServer() :
	fRecoServer(0)
{
	fRecoServer = new AliEventServerReconstruction();
	for(Int_t i=0; i<5; ++i)
	{
		fDimSORListener[i] = 0;
		fDimEORListener[i] = 0;
	}
	FillRunsFromDatabase();
	InitDIMListeners();
}

AliEventServer::~AliEventServer()
{
	for (Int_t i = 0; i < 5; ++i)
	{
		if(fDimSORListener[i]) delete fDimSORListener[i];
		if(fDimEORListener[i]) delete fDimEORListener[i];
		
		fDimSORListener[i] = 0;
		fDimEORListener[i] = 0;
	}
	if(fRecoServer){delete fRecoServer;fRecoServer=0;}
}

void AliEventServer::InitDIMListeners()
{
	// DIM interface.  
	for (Int_t i = 0; i < 5; ++i)
	{
		if (i == 0)
		{
			fDimSORListener[i] = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS");
			fDimEORListener[i] = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS");
		}
		else
		{
			fDimSORListener[i] = new AliDimIntNotifier(Form("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS_%d", i));
			fDimEORListener[i] = new AliDimIntNotifier(Form("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS_%d", i));
		}
    
		fDimSORListener[i]->Connect("DimMessage(Int_t)", "AliEventServer", this, "StartOfRun(Int_t)");
		fDimEORListener[i]->Connect("DimMessage(Int_t)", "AliEventServer", this, "EndOfRun(Int_t)");
	}

}

void AliEventServer::StartOfRun(Int_t run)
{
  cout<<"SOR signal received for run:"<<run<<endl;
  if(run<=0) return;
  fRecoServer->StopReconstruction();
	
  TEnv settings;
  settings.ReadFile(AliEventServerUtil::GetPathToServerConf(), kEnvUser);
    
  TString dataSource = settings.GetValue("data.source", DEFAULT_DATA_SOURCE);
  TString eventSource;
    
  if(dataSource=="local")
    {
      cout<<"Starting Reco for run "<<run<<endl;
      eventSource = Form("mem://%s/run%d", gSystem->Getenv("ONLINERECO_RAWFILES_DIR"), run);
    }
  else if(dataSource=="run")
    {
      cout<<"Starting Reco for GDCs active in current run:"<<run<<endl;
      eventSource = "mem://@*:";
    }
     
  fRecoServer->StartReconstruction(run, eventSource.Data());
}

void AliEventServer::EndOfRun(Int_t run)
{
  cout<<"EOR signal received for run:"<<run<<endl;
  if(run<=0) return;
  fRecoServer->StopReconstruction();
}

void AliEventServer::FillRunsFromDatabase()
{
	TEnv settings;
	settings.ReadFile(AliEventServerUtil::GetPathToServerConf(), kEnvUser);
	
	TString dbHost = settings.GetValue("logbook.host", DEFAULT_LOGBOOK_HOST);
	TString dbPort =  Form("%d", settings.GetValue("logbook.port", DEFAULT_LOGBOOK_PORT));
	TString dbName =  settings.GetValue("logbook.db", DEFAULT_LOGBOOK_DB);
	TString user =  settings.GetValue("logbook.user", DEFAULT_LOGBOOK_USER);
	TString password = settings.GetValue("logbook.pass", DEFAULT_LOGBOOK_PASS);

	TString connStr = Form("mysql://%s:%s/%s", dbHost.Data(), dbPort.Data(), dbName.Data()) ;

	cout<<"connecting to:"<<connStr<<endl;

	TSQLServer* server = TSQLServer::Connect(connStr.Data(), user.Data(), password.Data());
	if (!server)
	{
	  cout<<"ERROR: Could not connect to DAQ Logbook"<<endl;
		return;
	}
	TString sqlQuery;
	TTimeStamp ts;
	sqlQuery.Form("SELECT run FROM logbook WHERE DAQ_time_start > %u AND DAQ_time_end IS NULL AND `partition` REGEXP 'PHYSICS.*'",
		      (UInt_t)ts.GetSec()-86400);
	TSQLResult* result = server->Query(sqlQuery);
	if (!result)
	{
	  cout<<"ERROR: Can't execute query:"<< sqlQuery<<endl;
		return;
	}
	if (result->GetRowCount() != 0)
	{
		for (Int_t iRow = 0; iRow < result->GetRowCount(); iRow++)
		{
			TSQLRow* row = result->Next();
			TString runStr = row->GetField(0);
			if (runStr.IsDigit())
				StartOfRun(runStr.Atoi());
			delete row;
		}
	}
	delete result;

}


/*
void AliEventServer::FinishedReconstruction(Int_t status)
{
  // Slot called on termination of child process.
  Int_t run = fServer->GetRunId();
	
  Info("FinishedReconstruction", "Reconstruction Thread finished \tRunId:%d \tstatus=%d", run, status);

  mIntInt_i i =fRun2PidMap.find(run);
  if (i != fRun2PidMap.end())
    {
      fRunList->RemoveEntry(run);
    
      // clean (remove) run's reconstructed directory
      //gSystem->Exec(Form("rm -rf %s/reco/run%d_%d",gSystem->Getenv("ONLINERECO_BASE_DIR"),run,pid));
      
      if (status == 0)
	{
	  fRunList->AddEntrySort(TString::Format("%-20d -- PROCESSED", run), run);
	}
      else
	{
	  fRunList->AddEntrySort(TString::Format("%-20d -- PROCESSED [%d]", run, status), run);
	}
      fRunList->Layout();
    
    }
  else
    {
      Warning("FinishedReconstruction", "Run number %d not registered.", run);
    }

    }*/
