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

RMacro::RMacro(const char* name, const char* /*title*/) : TMacro(name, "") {}

/**************************************************************************/

void RMacro::Exec(const char* params)
{
  if(Reve::CheckMacro(fName))
    G__unloadfile(fTitle.Data());

  // Copy from TMacro::Exec. Difference is that the file is really placed
  // into the /tmp.
  TString fname = "/tmp/";
  {
    //the current implementation uses a file in the current directory.
    //should be replaced by a direct execution from memory by CINT
    fname += GetName();
    fname += ".Cexec";
    SaveSource(fname);
    //disable a possible call to gROOT->Reset from the executed script
    gROOT->SetExecutingMacro(kTRUE);
    //execute script in /tmp
    TString exec = ".x " + fname;
    TString p = params;
    if (p == "") p = fParams;
    if (p != "")
      exec += "(" + p + ")";
    gROOT->ProcessLine(exec);
    //enable gROOT->Reset
    gROOT->SetExecutingMacro(kFALSE);
    //delete the temporary file
    gSystem->Unlink(fname);
  }

  G__unloadfile(fname);
}
