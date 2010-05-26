// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *)
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliOnlineReco.h"
#include "AliChildProcTerminator.h"
#include "AliDimIntNotifier.h"
#include "AliCDBManager.h"
#include "AliGRPPreprocessor.h"

#include <TTimer.h>

#include <TGListBox.h>
#include <TGButton.h>

#include <TInterpreter.h>

#include <unistd.h>
#include <signal.h>

//______________________________________________________________________________
// Full description of AliOnlineReco
//

ClassImp(AliOnlineReco)

AliOnlineReco::AliOnlineReco() :
  TGMainFrame(gClient->GetRoot(), 400, 400),

  fRunList(0), fAutoRun(0), fStartButt(0), fStopButt(0), fExitButt(0),
  fAutoRunTimer(0), fAutoRunScheduled(0), fAutoRunRunning(0),
  fRun2PidMap(),
  fTestMode(kFALSE)
{
  // Constructor.

  // GUI components.
  fRunList = new TGListBox(this);
  AddFrame(fRunList, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));

  TGHorizontalFrame *hf = new TGHorizontalFrame(this, 1, 20);

  fAutoRun = new TGCheckButton(hf, "AutoRun");
  hf->AddFrame(fAutoRun, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  fAutoRun->Connect("Clicked()", "AliOnlineReco", this, "DoAutoRun()");

  fStartButt = new TGTextButton(hf, "Start");
  hf->AddFrame(fStartButt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  fStartButt->Connect("Clicked()", "AliOnlineReco", this, "DoStart()");

  fStopButt = new TGTextButton(hf, "Stop");
  hf->AddFrame(fStopButt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  fStopButt->Connect("Clicked()", "AliOnlineReco", this, "DoStop()");

  fExitButt = new TGTextButton(hf, "Exit");
  hf->AddFrame(fExitButt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  fExitButt->Connect("Clicked()", "AliOnlineReco", this, "DoExit()");

  AddFrame(hf, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));

  MapSubwindows();
  Layout();
  SetWindowName("Alice Online Reconstruction");

  // DIM interface.  
  for (Int_t i = 0; i < 5; ++i)
  {
    if (i == 0)
    {
      fSOR[i] = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS");
      fEOR[i] = new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS");
    }
    else
    {
      fSOR[i] = new AliDimIntNotifier(Form("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS_%d", i));
      fEOR[i] = new AliDimIntNotifier(Form("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS_%d", i));
    }

    fSOR[i]->Connect("DimMessage(Int_t)", "AliOnlineReco", this, "StartOfRun(Int_t)");
    fEOR[i]->Connect("DimMessage(Int_t)", "AliOnlineReco", this, "EndOfRun(Int_t)");
  }

  const Int_t autoRunDelay = 10; // should go to config
  fAutoRunTimer = new TTimer(autoRunDelay * 1000l);
  fAutoRunTimer->Connect("Timeout()", "AliOnlineReco", this, "AutoRunTimerTimeout()");

  // Signal handlers
  // ROOT's TSignalHAndler works not SIGCHLD ...
  AliChildProcTerminator::Instance()->Connect("ChildProcTerm(Int_t,Int_t)", "AliOnlineReco", this, "ChildProcTerm(Int_t,Int_t)");
}

AliOnlineReco::~AliOnlineReco()
{
  // Destructor.

  delete fAutoRunTimer;
}

Int_t AliOnlineReco::GetLastRun() const
{
  // Returns the last started run.

  return fRun2PidMap.empty() ? 0 : fRun2PidMap.rbegin()->first;
}

Bool_t AliOnlineReco::GetAutoRunMode() const
{
  // Return state of auto-run flag.

  return fAutoRun->IsOn();
}

void AliOnlineReco::SetAutoRunMode(Bool_t ar)
{
  // Set auto-run flag.

  if (ar == fAutoRun->IsOn())
    return;

  fAutoRun->SetState(ar ? kButtonDown : kButtonUp, kTRUE);
}

//------------------------------------------------------------------------------
// Private methods
//------------------------------------------------------------------------------

AliOnlineReco::mIntInt_i AliOnlineReco::FindMapEntryByPid(Int_t pid)
{
  // Find run-to-pid map iterator by pid.
  // Requires iteration over map.

  for (mIntInt_i i = fRun2PidMap.begin(); i != fRun2PidMap.end(); ++i)
  {
    if (i->second == pid)
      return i;
  }

  return fRun2PidMap.end();
}

void AliOnlineReco::StartAliEve(mIntInt_i& mi)
{
  // Start alieve to process run given my the run-pid entry.

  Int_t run = mi->first;

  if (mi->second == 0)
  {
    pid_t pid = fork();
    if (pid == -1)
    {
      perror("DoStart -- Fork failed");
      return;
    }

    if (pid)
    {
      mi->second = pid;
      fRunList->RemoveEntry(run);
      fRunList->AddEntrySort(TString::Format("%-20d -- RUNNING", run), run);
      fRunList->Layout();
    }
    else
    {
      gCINTMutex = 0;

      struct sigaction sac;
      sac.sa_handler = 0;
      sigemptyset(&sac.sa_mask);
      sac.sa_flags = 0;
      sigaction(SIGCHLD, &sac, 0);

      int s;
      if (fTestMode)
      {
	s = execlp("alitestproc", "alitestproc", TString::Format("%d", run).Data(), (char*) 0);
      }
      else
      {
	Int_t procPID = gSystem->GetPid();
	TString logFile = Form("%s/reco/log/run%d_%d.log",
			       gSystem->Getenv("ONLINERECO_BASE_DIR"),
			       run,
			       (Int_t)procPID);
	Info("DoStart","Reconstruction log will be written to %s",logFile.Data());
	gSystem->RedirectOutput(logFile.Data());

	gSystem->cd(Form("%s/reco",gSystem->Getenv("ONLINERECO_BASE_DIR")));

	TString gdcs;
	if (RetrieveGRP(run,gdcs) <= 0 || gdcs.IsNull()) 
	  gSystem->Exit(1);

	gSystem->Setenv("DATE_RUN_NUMBER", Form("%d", run));
	// Setting CDB
// 	AliCDBManager * man = AliCDBManager::Instance();
// 	man->SetDefaultStorage("local:///local/cdb");
// 	man->SetSpecificStorage("GRP/GRP/Data",
// 			      Form("local://%s",gSystem->pwd()));
// 	man->SetSpecificStorage("GRP/CTP/Config",
// 			      Form("local://%s",gSystem->pwd()));
// 	man->SetSpecificStorage("ACORDE/Align/Data",
// 				"local://$ALICE_ROOT/OCDB");

	gSystem->mkdir(Form("run%d_%d", run, (Int_t)procPID));
	gSystem->cd(Form("run%d_%d", run, (Int_t)procPID));

	TString recMacroPath(gSystem->Getenv("ONLINERECO_MACRO"));
	if (recMacroPath.IsNull()) {
	  recMacroPath = "$ALICE_ROOT/MONITOR/rec.C";
	}

	s = execlp("alieve",
		   "alieve",
		   "-q",
		   Form("%s(\"mem://@*:\")", gSystem->ExpandPathName(recMacroPath.Data())),
		   (char*) 0);

	gSystem->Exec(Form("rm -rf %s/reco/run%d_%d",gSystem->Getenv("ONLINERECO_BASE_DIR"),run,(Int_t)procPID));
      }

      if (s == -1)
      {
	perror("execlp failed - this will not end well");
	gSystem->Exit(1);
      }
    }
  }
  else
  {
    Error("DoStart", "Process already running.");
  }
}

void AliOnlineReco::KillPid(Int_t pid)
{
  // Terminate process given by pid.

  // Send terminate signal to process ...

  if (fTestMode)
  {
    kill(pid, SIGTERM);
  }
  else
  {
    // alieve will auto-destruct on SIGUSR1
    kill(pid, SIGUSR1);
  }
}

void AliOnlineReco::StartAutoRunTimer(Int_t run)
{
  // Start timer for given run.
  // If an auto-started run is already active, this call is ignored.
  // If timer is already active, it is restarted.

  if (fAutoRunRunning)
    return;

  fAutoRunTimer->Reset();
  fAutoRunTimer->TurnOn();
  fAutoRunScheduled = run;

  Info("StartAutoRunTimer", "Scheduling run %d for auto-display.", run);
}

void AliOnlineReco::StopAutoRunTimer()
{
  // Stop auto-run timer.

  fAutoRunTimer->TurnOff();
  fAutoRunScheduled = 0;
}

void AliOnlineReco::AutoRunTimerTimeout()
{
  // Slot called on auto-timer timeout.

  Int_t run = fAutoRunScheduled;

  StopAutoRunTimer();

  mIntInt_i i = fRun2PidMap.find(run);

  if (i == fRun2PidMap.end())
  {
    Warning("AutoRunTimerTimeout", "run no longer active.");
    return;
  }

  Info("AutoRunTimerTimeout", "Starting display for run %d.", run);

  StartAliEve(i);
  fAutoRunRunning = run;
}

//------------------------------------------------------------------------------
// Handlers of DIM signals.
//------------------------------------------------------------------------------

void AliOnlineReco::StartOfRun(Int_t run)
{
  // Slot called from DIM handler on start of run.

  mIntInt_i i = fRun2PidMap.find(run);
  if (i == fRun2PidMap.end())
  {
    fRun2PidMap[run] = 0;
    fRunList->AddEntrySort(TString::Format("%d", run), run);
    fRunList->Layout();

    if (fAutoRun->IsOn())
    {
      StartAutoRunTimer(run);
    }
  }
  else
  {
    Error("StartOfRun", "Run %d already registered.", run);
  }
}

void AliOnlineReco::EndOfRun(Int_t run)
{
  // Slot called from DIM handler on stop of run.

  mIntInt_i i = fRun2PidMap.find(run);
  if (i != fRun2PidMap.end())
  {
    Int_t pid = i->second;
    fRunList->RemoveEntry(run);
    fRunList->Layout();
    fRun2PidMap.erase(i);
    if (pid)
    {
      KillPid(pid);
    }
    gClient->NeedRedraw(fRunList);

    if (fAutoRunRunning == run)
    {
      fAutoRunRunning = 0;
    }
  }
  else
  {
    Error("EndOfRun", "Run %d not registered.", run);
  }
}

//------------------------------------------------------------------------------
// Handlers of OS signals.
//------------------------------------------------------------------------------

void AliOnlineReco::ChildProcTerm(Int_t pid, Int_t status)
{
  // Slot called on termination of child process.

  printf("child process termination pid=%d, status=%d...\n", pid, status);

  mIntInt_i i = FindMapEntryByPid(pid);
  if (i != fRun2PidMap.end())
  {
    Int_t run = i->first;
    fRunList->RemoveEntry(run);
    if (status == 0)
    {
      fRunList->AddEntrySort(TString::Format("%-20d -- PROCESSED", run), run);
    }
    else
    {
      fRunList->AddEntrySort(TString::Format("%-20d -- PROCESSED [%d]", run, status), run);
    }
    fRunList->Layout();
    i->second = 0;

    if (fAutoRunRunning == run && fAutoRun->IsOn())
    {
      fAutoRunRunning = 0;
      StartAutoRunTimer(run);
    }
    else
    {
      fAutoRunRunning = 0;
    }
  }
  else
  {
    Warning("ChildProcTerm", "Process with pid=%d not registered.", pid);
  }
}

//------------------------------------------------------------------------------
// Handlers of button signals.
//------------------------------------------------------------------------------

void AliOnlineReco::DoAutoRun()
{
  // Slot called from auto-run check-box.

  Bool_t autoRun = fAutoRun->IsOn();

  if (autoRun)
    fStartButt->SetEnabled(kFALSE);
  else
    fStartButt->SetEnabled(kTRUE);    
}

void AliOnlineReco::DoStart()
{
  // Slot called from Start button.

  Int_t run = fRunList->GetSelected();
  mIntInt_i i = fRun2PidMap.find(run);

  if (i == fRun2PidMap.end())
  {
    Error("DoStart", "no selection");
    return;
  }

  StartAliEve(i);
}

void AliOnlineReco::DoStop()
{
  // Slot called from Stop button.

  Int_t run = fRunList->GetSelected();
  mIntInt_i i = fRun2PidMap.find(run);

  if (i == fRun2PidMap.end())
  {
    Error("DoStop", "no selection");
    return;
  }

  Int_t pid = i->second;
  if (pid)
  {
    KillPid(pid);
  }
  else
  {
    Error("DoStop", "Process not running.");
  }
}

void AliOnlineReco::DoExit()
{
  // Slot called from Exit button.

  gSystem->ExitLoop();
}

void AliOnlineReco::CloseWindow()
{
  // Virtual method called when window-manager close-window button is pressed.

  gSystem->ExitLoop();
}

Int_t AliOnlineReco::RetrieveGRP(UInt_t run, TString &gdc)
{
  // Retrieve GRP entry for given run from aldaqdb.

  Int_t ret=AliGRPPreprocessor::ReceivePromptRecoParameters(run, "aldaqdb", 0, "LOGBOOK", "logbook", "alice",
							    Form("local://%s",gSystem->pwd()),
							    gdc);
  if(ret>0) Info("RetrieveGRP","Last run of the same type is: %d",ret);
  else if(ret==0) Warning("RetrieveGRP","No previous run of the same type found");
  else if(ret<0) Error("Retrieve","Error code while retrieving GRP parameters returned: %d",ret);
  return(ret);
}
