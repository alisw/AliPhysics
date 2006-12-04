// $Header$

#include "RMacro.h"

#include <TSystem.h>
#include <TROOT.h>
#include <G__ci.h>

using namespace Reve;

//______________________________________________________________________
// RMacro
//
// Sub-class of TMacro, overriding Exec to unload the previous verison
// and cleanup after the execution.

ClassImp(RMacro)

RMacro::RMacro() : TMacro() {}

RMacro::RMacro(const RMacro& m) : TMacro(m) {}

RMacro::RMacro(const char* name) :
  TMacro()
{
  if (!name) return;

  fTitle = name;

  char *dot   = (char*)strrchr(name, '.');
  char *slash = (char*)strrchr(name, '/');
  if (dot) *dot = 0;
  if (slash) fName = slash + 1;
  else       fName = name;

  ReadFile(fTitle);
}

/**************************************************************************/

#include <TTimer.h>

void RMacro::Exec(const char* params)
{
  if(Reve::CheckMacro(fTitle.Data()))
    G__unloadfile(fTitle.Data());

  // Copy from TMacro::Exec. Difference is that the file is really placed
  // into the /tmp.
  TString fname = "/tmp/";
  {
    //the current implementation uses a file in the current directory.
    //should be replaced by a direct execution from memory by CINT
    fname += GetName();
    fname += ".C";
    SaveSource(fname);
    //disable a possible call to gROOT->Reset from the executed script
    gROOT->SetExecutingMacro(kTRUE);
    //execute script in /tmp
    TString exec = ".x " + fname;
    TString p = params;
    if (p == "") p = fParams;
    if (p != "")
      exec += "(" + p + ")";
    Int_t exit;
    gROOT->ProcessLine(exec, &exit);
    //enable gROOT->Reset
    gROOT->SetExecutingMacro(kFALSE);
    //delete the temporary file
    gSystem->Unlink(fname);
  }

  G__unloadfile(fname);

  // In case an exception was thrown (which i do not know how to detect
  // the execution of next macros does not succeed.
  // However strange this might seem, this solves the problem.
  TTimer::SingleShot(100, "Reve::RMacro", this, "ResetRoot()");
}

#include <TApplication.h>

void RMacro::ResetRoot()
{
  // printf ("RMacro::ResetRoot doing 'gROOT->Reset()'.\n");
  gROOT->GetApplication()->ProcessLine("gROOT->Reset()");
}
