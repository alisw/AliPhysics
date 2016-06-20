// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include <TObject.h>
#include <TList.h>
#include <TTimer.h>
#include <TROOT.h>

#include <TThread.h>
#include <TMutex.h>
#include <TCondition.h>

class CommandQueue : public TObject
{
public:
  class Command : public TObject
  {
  public:
    TString     fCommand;
    TCondition *fAwakeCond; // awake thread after command executed

    Command(const TString& cmd, TCondition* cond) : fCommand(cmd), fAwakeCond(cond) {}

    ClassDef(Command, 0); // An entry in the command-queue.
  };

  TList  fQueue;
  TMutex fQueueLock;
  TTimer fQueueTimer;
  Bool_t fQueueTimerSet;

  CommandQueue()
  {
    fQueueTimer.Connect("Timeout()", "CommandQueue", this, "ProcessQueue()");
    fQueueTimerSet = kFALSE;
  }

  // Destructor, cleanup missing.

  //--------------------------------

  void RegisterCommand(const TString& cmd, TCondition* cond=0)
  {
    fQueueLock.Lock();
    fQueue.Add(new Command(cmd, cond));
    if (!fQueueTimerSet) {
      // Force execution in main thread when it is free.
      fQueueTimer.Start(0, kTRUE);
      fQueueTimerSet = kTRUE;
    }
    fQueueLock.UnLock();
  }

  virtual void ProcessQueue()
  {
    Int_t n_proccessed = 0;

    while (1) {
      fQueueLock.Lock();
      Command* c = (Command*) fQueue.First();
      fQueue.RemoveFirst();
      fQueueLock.UnLock();

      { // Put this into virtual void ProcessComand(Command*);
	// Need also exception handling.
	gROOT->ProcessLineFast(c->fCommand);
	++n_proccessed;
	if (c->fAwakeCond != 0) {
	  c->fAwakeCond->GetMutex()->Lock();
	  c->fAwakeCond->Signal();
	  c->fAwakeCond->GetMutex()->UnLock();
	}
      }
      delete c;

      fQueueLock.Lock();
      if (fQueue.IsEmpty()) {
	fQueueTimerSet = kFALSE;
	fQueueLock.UnLock();
	break;
      }
      fQueueLock.UnLock();
    }

    printf("CommandQueue::ProcessQueue() processed %d commands.\n", n_proccessed);
  }

  ClassDef(CommandQueue,0); // Command queue for exection in global CINT context.
};

//=============================================================================
// Global variable
//=============================================================================

CommandQueue* g_cmd_queue = 0;

void command_queue()
{
  g_cmd_queue = new CommandQueue;
  printf("Starting command-queue ...\n");
}

// ============================================================================
// TEveUtil side
// ============================================================================

#include <TEve.h>
#include <TEveManager.h>
#include <TEvePointSet.h>

#include <TRandom.h>

void make_crap(void* arg)
{
  Int_t num = 1024;
  TRandom rnd(0);

  TEvePointSet* ps = new TEvePointSet("Testus", num);
  for (Int_t i=0; i<num; ++i) {
    ps->SetNextPoint(rnd.Uniform(-100, 100),
		     rnd.Uniform(-100, 100),
		     rnd.Uniform(-100, 100));
  }
  ps->SetMainColor(kRed);
  printf("make_crap() -> produced TEvePointSet* %p)\n", ps);

  ((CommandQueue*)arg)->RegisterCommand
    (Form("register_crap((TEveElement*)0x%lx)", ps));
}

void register_crap(TEveElement* el)
{
  printf("register_crap(TEveElement* %p)\n", el);
  gEve->AddElement(el);
  gEve->Redraw3D();
}

void test()
{
  TThread* thr = new TThread(make_crap, g_cmd_queue);
  thr->Run();
}
