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

#include <TGButton.h>
#include <TGListBox.h>
#include <TGTab.h>
#include <TGTextEntry.h>
#include <TGToolBar.h>
#include <TG3DLine.h>

#include <AliLog.h>
#include <AliReconstruction.h>

#include "AliEventServerUtil.h"
#include "AliEventServerWindow.h"
#include "AliEventServerPreferencesWindow.h"
#include "AliDimIntNotifier.h"
#include "AliRecoServer.h"

//______________________________________________________________________________
// Full description of AliEventServerWindow
//

ClassImp(AliEventServerWindow)

AliEventServerWindow::AliEventServerWindow() :
  TGMainFrame(gClient->GetRoot(), 400, 400),
  fRunList(0), 
  fStartServButt(0), 
  fStopServButt(0), 
  fExitButt(0),
  fRunRunning(0), 
  fRecoServer(0)
{
	SetCleanup(kDeepCleanup);

	SetupToolbar();
  
	fRunList = new TGListBox(this);
  
	AddFrame(fRunList,  new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
	
  for(Int_t i=0; i<5; ++i)
  {
  	fDimSORListener[i] = 0;
  	fDimEORListener[i] = 0;
  }
  
  Connect("CloseWindow()", "AliEventServerWindow", this, "onExit()");
  SetWindowName("ALICE Event Server");
 
  MapSubwindows();
  Resize(250,300);
  MapWindow();
  
  FillRunsFromDatabase();
  InitDIMListeners();
}

AliEventServerWindow::~AliEventServerWindow()
{
  // Destructor.
	
	for (Int_t i = 0; i < 5; ++i)
	{
		if(fDimSORListener[i]) delete fDimSORListener[i];
		if(fDimEORListener[i]) delete fDimEORListener[i];
		
		fDimSORListener[i] = 0;
		fDimEORListener[i] = 0;
	}
 
}

void AliEventServerWindow::InitDIMListeners()
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
    
    fDimSORListener[i]->Connect("DimMessage(Int_t)", "AliEventServerWindow", this, "StartOfRun(Int_t)");
    fDimEORListener[i]->Connect("DimMessage(Int_t)", "AliEventServerWindow", this, "EndOfRun(Int_t)");
  }

}

void AliEventServerWindow::FillRunsFromDatabase()
{
	TEnv settings(ALIEVENTSERVER_CONF);

	TString dbHost = settings.GetValue("logbook.host", DEFAULT_LOGBOOK_HOST);
	TString dbPort =  Form("%d", settings.GetValue("logbook.port", DEFAULT_LOGBOOK_PORT));
	TString dbName =  settings.GetValue("logbook.db", DEFAULT_LOGBOOK_DB);
	TString user =  settings.GetValue("logbook.user", DEFAULT_LOGBOOK_USER);
	TString password = settings.GetValue("logbook.pass", DEFAULT_LOGBOOK_PASS);

  TSQLServer* server = TSQLServer::Connect(Form("mysql://%s:%s/%s", dbHost.Data(), dbPort.Data(), dbName.Data()), user.Data(), password.Data());
  if (!server) {
    AliWarning("ERROR: Could not connect to DAQ Logbook");
    return;
  }
  TString sqlQuery;
  TTimeStamp ts;
  sqlQuery.Form("SELECT run FROM logbook WHERE DAQ_time_start > %u AND DAQ_time_end IS NULL AND partition REGEXP 'PHYSICS.*'",
    (UInt_t)ts.GetSec()-86400);
  TSQLResult* result = server->Query(sqlQuery);
  if (!result)
  {
    AliWarning( Form("ERROR: Can't execute query <%s>!", sqlQuery.Data()) );
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

void AliEventServerWindow::SetupToolbar()
{
	TGToolBar* mToolBar = new TGToolBar(this);
	mToolBar->AddButton(this, new TGPictureButton(mToolBar, Form("%s/MONITOR/icons/start.png", gSystem->Getenv("ALICE_ROOT")), TOOLBUTTON_START ) );
	mToolBar->AddButton(this, new TGPictureButton(mToolBar, Form("%s/MONITOR/icons/stop.png", gSystem->Getenv("ALICE_ROOT")), TOOLBUTTON_STOP) );
	mToolBar->AddButton(this, new TGPictureButton(mToolBar, Form("%s/MONITOR/icons/preferences.png", gSystem->Getenv("ALICE_ROOT")), TOOLBUTTON_PREFERENCES) );
	mToolBar->AddButton(this, new TGPictureButton(mToolBar, Form("%s/MONITOR/icons/exit.png", gSystem->Getenv("ALICE_ROOT")), TOOLBUTTON_EXIT) );
	
	mToolBar->Connect("Clicked(Int_t)", "AliEventServerWindow", this, "HandleToolBarAction(Int_t)");
	
	AddFrame(mToolBar, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
	AddFrame(new TGHorizontal3DLine(this), new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
}

void AliEventServerWindow::HandleToolBarAction(Int_t id)
{
	if(id==-1) return;
	
	switch(id){
	case TOOLBUTTON_START:{
		onStartServer();
		break;
	}
	case TOOLBUTTON_STOP:{
		onStopServer();
		break;
	}
	case TOOLBUTTON_PREFERENCES:{
		new AliEventServerPreferencesWindow(this, "Settings");
		break;
	}
	case TOOLBUTTON_EXIT:{
		onExit();
		break;
	}
	
	}
	
}

/*
void AliEventServerWindow::FinishedReconstruction(Int_t status)
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

}
*/
//------------------------------------------------------------------------------
// Private methods
//------------------------------------------------------------------------------

void AliEventServerWindow::StartReco(Int_t run)
{
  AliInfo(Form("Starting Reco for run %d", run));
  
  TString eventSource = Form("mem://%s/run%d", gSystem->Getenv("ONLINERECO_RAWFILES_DIR"), run);
  
  if(!fRecoServer) LaunchRecoServer();
    
  fRecoServer->StartReconstruction(run, eventSource.Data());
  
  if(fRecoServer->IsListenning()){
    fRunList->RemoveEntry(run);
    fRunList->AddEntrySort(TString::Format("%-20d -- RUNNING", run), run);
    fRunList->Layout();
  }
}


//------------------------------------------------------------------------------
// Handlers of DIM signals.
//------------------------------------------------------------------------------

void AliEventServerWindow::StartOfRun(Int_t run)
{
	if(run<=0) return;

  // Slot called from DIM handler on start of run.
	AliInfo(Form("called for Run %d ", run));

	fRunList->AddEntrySort(TString::Format("%d", run), run);
	fRunList->Layout();
	gClient->NeedRedraw(fRunList);
	
	StartReco(run);
}

void AliEventServerWindow::EndOfRun(Int_t run)
{
	if(run<=0) return;
	
   // Slot called from DIM handler on stop of run.
	AliInfo(Form("called for Run %d", run) );
	if(fRecoServer) fRecoServer->StopReconstruction();
	
  fRunList->RemoveEntry(run);
  fRunList->Layout();
  gClient->NeedRedraw(fRunList);
}

///------------------------------------------------------------------------------
// Handlers of button signals.
//------------------------------------------------------------------------------

void AliEventServerWindow::onStartServer()
{
  // Slot called from Start button.
  AliInfo("Starting server...");
	if(fRecoServer!=0) StopRecoServer();
	
	LaunchRecoServer();
}

void AliEventServerWindow::onStopServer()
{
  // Slot called from Stop button.
	AliInfo("Closing server...");
	
	StopRecoServer();
}

void AliEventServerWindow::onExit()
{
	AliInfo("Closing server & Exiting...");
	
	StopRecoServer();
	CloseWindow();
	
	gSystem->ExitLoop();
}

void AliEventServerWindow::LaunchRecoServer()
{
 	fRecoServer = new AliRecoServer;		
}

bool AliEventServerWindow::StopRecoServer()
{
	if(fRecoServer==0) return true;

	AliInfo("Closing server and stoping process...");
	
	delete fRecoServer;
	fRecoServer=0;	

	return true;
}
