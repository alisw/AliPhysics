#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MONITOR/AliMonitorClient.h"
#include "TSystem.h"
#endif

void client()
{
  // load libraries
  if (!gROOT->GetClass("AliLevel3")) {
    gSystem->Load("libAliHLTSrc.so");
    gSystem->Load("libAliHLTMisc.so");
    gSystem->Load("libAliHLTHough.so");
    gSystem->Load("libAliHLTComp.so");
  }
  if (!gROOT->GetClass("AliMonitorClient")) {
    gSystem->Load("libMONITOR.so");
  }
  new AliMonitorClient;
}
