#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MONITOR/AliMonitorClient.h"
#include "TSystem.h"
#endif

void client()
{
  // load libraries
  if (!gROOT->GetClass("AliLevel3")) {
    gSystem->Load("libAliL3Src.so");
    gSystem->Load("libAliL3Misc.so");
    gSystem->Load("libAliL3Hough.so");
    gSystem->Load("libAliL3Comp.so");
  }
  if (!gROOT->GetClass("AliMonitorClient")) {
    gSystem->Load("libMONITOR.so");
  }
  new AliMonitorClient;
}
