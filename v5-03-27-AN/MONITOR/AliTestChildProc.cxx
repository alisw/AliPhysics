// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliTestChildProc.h"

#include <TGButton.h>

#include <TSystem.h>

#include <cstdio>
#include <cstdlib>


//==============================================================================
// AliTestChildProc
//==============================================================================

//______________________________________________________________________________
// Full description of AliTestChildProc
//

ClassImp(AliTestChildProc)

AliTestChildProc::AliTestChildProc(Int_t run) :
  TGMainFrame(gClient->GetRoot(), 400, 200)
{
  TGTextButton *b;

  b = new TGTextButton(this, "Exit");
  AddFrame(b, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  b->Connect("Clicked()", "AliTestChildProc", this, "DoExit()");

  b = new TGTextButton(this, "Crash");
  AddFrame(b, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  b->Connect("Clicked()", "AliTestChildProc", this, "DoCrash()");

  b = new TGTextButton(this, "Xyzz");
  AddFrame(b, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  b->Connect("Clicked()", "AliTestChildProc", this, "DoXyzz()");

  MapSubwindows();
  Layout();
  SetWindowName(TString::Format("Test window %d", run));
}

void AliTestChildProc::DoExit()
{
  printf("doing exit ...\n");
  gSystem->ExitLoop();
}

void AliTestChildProc::DoCrash()
{
  printf("doing crash ...\n");
  TObject *p = 0;
  p->Dump();
}

void AliTestChildProc::DoXyzz()
{
  printf("doing xyzz ...\n");
}
