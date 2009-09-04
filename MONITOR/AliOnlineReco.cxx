// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliOnlineReco.h"
#include "AliChildProcTerminator.h"
#include "AliDimIntNotifier.h"

#include <TGListBox.h>
#include <TGButton.h>

#include <unistd.h>
#include <signal.h>

//______________________________________________________________________________
// Full description of AliOnlineReco
//

ClassImp(AliOnlineReco)

AliOnlineReco::AliOnlineReco() :
  TGMainFrame(gClient->GetRoot(), 400, 400),

  fRunList(0), fStartButt(0), fStopButt(0), fXyzzButt(0),

  fSOR(new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_SOR_PHYSICS")),
  fEOR(new AliDimIntNotifier("/LOGBOOK/SUBSCRIBE/DAQ_EOR_PHYSICS")),

  fTestMode(kFALSE)
{
  // GUI components.
  fRunList = new TGListBox(this);
  AddFrame(fRunList, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));

  TGHorizontalFrame *hf = new TGHorizontalFrame(this, 1, 20);

  fStartButt = new TGTextButton(hf, "Start");
  hf->AddFrame(fStartButt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  fStartButt->Connect("Clicked()", "AliOnlineReco", this, "DoStart()");

  fStopButt = new TGTextButton(hf, "Stop");
  hf->AddFrame(fStopButt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  fStopButt->Connect("Clicked()", "AliOnlineReco", this, "DoStop()");

  fXyzzButt = new TGTextButton(hf, "Exit");
  hf->AddFrame(fXyzzButt, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  fXyzzButt->Connect("Clicked()", "AliOnlineReco", this, "DoXyzz()");

  AddFrame(hf, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));

  MapSubwindows();
  Layout();
  SetWindowName("Alice Online Reconstruction");

  // DIM interface.
  fSOR->Connect("DimMessage(Int_t)", "AliOnlineReco", this, "StartOfRun(Int_t)");
  fEOR->Connect("DimMessage(Int_t)", "AliOnlineReco", this, "EndOfRun(Int_t)");

  // Signal handlers
  // ROOT's TSignalHAndler works not SIGCHLD ...
  AliChildProcTerminator::Instance()->Connect("ChildProcTerm(Int_t,Int_t)", "AliOnlineReco", this, "ChildProcTerm(Int_t,Int_t)");
}

AliOnlineReco::mIntInt_i AliOnlineReco::FindMapEntryByPid(Int_t pid)
{
  for (mIntInt_i i = fRun2PidMap.begin(); i != fRun2PidMap.end(); ++i)
  {
    if (i->second == pid)
      return i;
  }

  return fRun2PidMap.end();
}

//------------------------------------------------------------------------------
// Handlers of DIM signals.
//------------------------------------------------------------------------------

void AliOnlineReco::StartOfRun(Int_t run)
{
  mIntInt_i i = fRun2PidMap.find(run);
  if (i == fRun2PidMap.end())
  {
    fRun2PidMap[run] = 0;
    fRunList->AddEntrySort(TString::Format("%d", run), run);
    fRunList->Layout();
  }
  else
  {
    Error("StartOfRun", "Run %d already registered.", run);
  }
}

void AliOnlineReco::EndOfRun(Int_t run)
{
  mIntInt_i i = fRun2PidMap.find(run);
  if (i != fRun2PidMap.end())
  {
    Int_t pid = i->second;
    fRunList->RemoveEntry(run);
    fRunList->Layout();
    fRun2PidMap.erase(i);
    if (pid)
    {
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
    gClient->NeedRedraw(fRunList);
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
  printf("child process termination pid=%d, status=%d...\n", pid, status);

  mIntInt_i i = FindMapEntryByPid(pid);
  if (i != fRun2PidMap.end())
  {
    Int_t run = i->first;
    fRunList->RemoveEntry(run);
    if (status == 0)
    {
      fRunList->AddEntrySort(TString::Format("%-20d -- FINISHED", run), run);
    }
    else
    {
      fRunList->AddEntrySort(TString::Format("%-20d -- CRASHED [%d]", run, status), run);
    }
    fRunList->Layout();
    i->second = 0;
  }
  else
  {
    Warning("ChildProcTerm", "Process with pid=%d not registered.", pid);
  }
}

//------------------------------------------------------------------------------
// Handlers of button signals.
//------------------------------------------------------------------------------

void AliOnlineReco::DoStart()
{
  Int_t run = fRunList->GetSelected();
  mIntInt_i i = fRun2PidMap.find(run);

  if (i == fRun2PidMap.end())
  {
    Error("DoStart", "no selection");
    return;
  }

  if (i->second == 0)
  {
    pid_t pid = fork();
    if (pid == -1)
    {
      perror("DoStart -- Fork failed");
      return;
    }

    if (pid)
    {
      i->second = pid;
      fRunList->RemoveEntry(run);
      fRunList->AddEntrySort(TString::Format("%-20d -- RUNNING", run), run);
      fRunList->Layout();
    }
    else
    {
      int s;
      if (fTestMode)
      {
	s = execlp("alitestproc", "alitestproc", TString::Format("%d", run).Data(), (char*) 0);
      }
      else
      {
	// !!!! Cvetan, add proper args here ...
	s = execlp("alieve", "alieve", TString::Format("%d", run).Data(), (char*) 0);
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

void AliOnlineReco::DoStop()
{
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
  else
  {
    Error("DoStop", "Process not running.");
  }
}

void AliOnlineReco::DoXyzz()
{
  gSystem->ExitLoop();
}

void AliOnlineReco::CloseWindow()
{
  gSystem->ExitLoop();
}
