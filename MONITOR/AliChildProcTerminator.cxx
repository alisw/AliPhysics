// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliChildProcTerminator.h"

#include <sys/wait.h>
#include <signal.h>

//______________________________________________________________________________
// Full description of AliChildProcTerminator
//

ClassImp(AliChildProcTerminator)

AliChildProcTerminator* AliChildProcTerminator::fgTheOne = 0;

AliChildProcTerminator* AliChildProcTerminator::Instance()
{
  if (fgTheOne == 0)
    fgTheOne = new AliChildProcTerminator;
  return fgTheOne;
}

AliChildProcTerminator::AliChildProcTerminator()
{
  struct sigaction sac;
  sac.sa_handler = sig_handler;
  sigemptyset(&sac.sa_mask);
  sac.sa_flags = 0;
  sigaction(SIGCHLD, &sac, 0);
}

void AliChildProcTerminator::sig_handler(int /*sig*/)
{
  int   status;
  pid_t pid = wait(&status);
  Instance()->ChildProcTerm(pid, status);
}

void AliChildProcTerminator::ChildProcTerm(Int_t pid, Int_t status)
{
   Long_t args[2];
   args[0] = (Long_t) pid;
   args[1] = (Long_t) status;

   Emit("ChildProcTerm(Int_t,Int_t)", args);
}
