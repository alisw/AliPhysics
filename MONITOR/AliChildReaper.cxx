// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliChildReaper.h"

#include <sys/wait.h>

//______________________________________________________________________________
// Full description of AliChildReaper
//

ClassImp(AliChildReaper)

AliChildReaper* AliChildReaper::fgTheOne = 0;

AliChildReaper* AliChildReaper::Instance()
{
  if (fgTheOne == 0)
    fgTheOne = new AliChildReaper;
  return fgTheOne;
}

AliChildReaper::AliChildReaper()
{
  struct sigaction sac;
  sac.sa_handler = sig_handler;
  sigemptyset(&sac.sa_mask);
  sac.sa_flags = 0;
  sigaction(SIGCHLD, &sac, 0);
}

void AliChildReaper::sig_handler(int /*sig*/)
{
  int   status;
  pid_t pid = wait(&status);
  Instance()->ChildDeath(pid, status);
}

void AliChildReaper::ChildDeath(Int_t pid, Int_t status)
{
   Long_t args[2];
   args[0] = (Long_t) pid;
   args[1] = (Long_t) status;

   Emit("ChildDeath(Int_t,Int_t)", args);
}
