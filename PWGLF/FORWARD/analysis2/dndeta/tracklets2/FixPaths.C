#ifndef __CINT__
#include <TROOT.h>
#include <TInterpreter.h>
#include <TSystem.h>
#include <TError.h>
#endif

struct FixPaths
{
  static FixPaths* fgInstance;
  FixPaths()
  {
    Printf("Fixing include path");
    gSystem->AddIncludePath("-I${ALICE_PHYSICS}/include");
    Printf("Include path: %s", gSystem->GetIncludePath());
  }
};

#ifndef __CINT__
FixPaths* FixPaths::fgInstance = new FixPaths;
#endif



