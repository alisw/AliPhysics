#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TError.h>
#include <TInterpreter.h>
#include "MONITOR/AliMonitorProcess.h"
#include "MONITOR/AliMonitorControl.h"
#endif

void monitor(const char* alienHost = "alien://", const char* alienDir = ".")
{
  // load libraries
  if (strcmp(gSystem->Getenv("ALIHLT_USEPACKAGE"), "ALIROOT") == 0) {
    if (!gROOT->GetClass("AliLevel3")) {
      gSystem->Load("libAliL3Src.so");
      gSystem->Load("libAliL3Misc.so");
      gSystem->Load("libAliL3Hough.so");
      gSystem->Load("libAliL3Comp.so");
    }
  }
  if (!gROOT->GetClass("AliMonitorProcess")) {
    gSystem->Load("libMONITOR.so");
  }

  // make sure galice.root and compression tables are there
  if (!gSystem->Which(".", "galice.root")) {
    gAlice->Init("$ALICE_ROOT/MONITOR/galice.C");
    gAlice->GetRunLoader()->Write();
    delete gAlice->GetRunLoader();
  }
  if (!gSystem->Which(".", "Table0.dat")) {
    gSystem->Exec("cp $ALICE_ROOT/RAW/Table*.dat .");
  }

  // start the monitoring
  AliMonitorProcess *process = new AliMonitorProcess(alienHost, alienDir);
  //  process->Run();
  new AliMonitorControl(process);
}
