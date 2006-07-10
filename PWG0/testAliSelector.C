/* $Id$ */

//
// This is an example how to run a selector
//
// This script runs the AliSelector on a chain of ESD files,
// some debug information is printed
//

#include "CreateESDChain.C"

// parameters are
//   dataDir: the directory containing subdirectories that contain the ESD files
//   nRuns: the number of files that should be processed
//   offset: the directory to start with
void testAliSelector(const Char_t* dataDir, Int_t nRuns = 5, Int_t offset = 0)
{
  // load needed libraries
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libESD");
  gSystem->Load("libPWG0base");

  // create chain, CreateESDChain() is defined in CreateESDChain.C
  TChain* chain = CreateESDChain(dataDir, nRuns, offset);

  // enable debugging
  AliLog::SetClassDebugLevel("AliSelector", AliLog::kInfo);

  // run selector on chain
  Long64_t result = chain->Process("AliSelector.cxx+");

  if (result != 0)
  {
    printf("ERROR: Executing process failed with %d.\n", result);
    return;
  }

  printf("Execution complete.\n");
}
