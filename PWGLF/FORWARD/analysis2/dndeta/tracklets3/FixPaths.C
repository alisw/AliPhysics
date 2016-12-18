/**
 * @file   FixPaths.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Tue Sep 20 17:10:39 2016
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
#ifndef __CINT__
#include <TROOT.h>
#include <TInterpreter.h>
#include <TSystem.h>
#include <TError.h>
#include <TString.h>
#else 
class TString;
class TSystem;
#endif

/**
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
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
//
// EOF
//



