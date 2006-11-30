#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TError.h>
#include <TInterpreter.h>
#include "AliRun.h"
#include "MONITOR/AliMonitorProcess.h"
#include "MONITOR/AliMonitorControl.h"
#endif

void monitor(Bool_t batchMode = kFALSE,
	     const char* selection = "ALL",
	     const char* alienHost = "alien://aliens7.cern.ch:15000/?direct",
	     const char* alienDir = "/alice_mdc/DC")
{
  // load libraries
  if (!gROOT->GetClass("AliLevel3")) {
    gSystem->Load("libAliHLTSrc.so");
    gSystem->Load("libAliHLTMisc.so");
    gSystem->Load("libAliHLTHough.so");
    gSystem->Load("libAliHLTComp.so");
  }
  if (!gROOT->GetClass("AliMonitorProcess")) {
    gSystem->Load("libMONITOR.so");
  }

  // make sure galice.root is there
  if (!gSystem->Which(".", "galice.root")) {
    gAlice->Init("$ALICE_ROOT/MONITOR/galice.C");
    gAlice->GetRunLoader()->Write();
    delete gAlice->GetRunLoader();
  }

  // start the monitoring
  AliMonitorProcess *process = new AliMonitorProcess(alienHost, alienDir,
						     selection);
  if (batchMode) {
    process->Run();
    delete process;
  } else {
    new AliMonitorControl(process);
  }
}
