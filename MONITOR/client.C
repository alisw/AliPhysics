#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MONITOR/AliMonitorClient.h"
#include "TSystem.h"
#endif

void client()
{
  // load libraries
  if (!gROOT->GetClass("AliLevel3")) {
    gSystem->Load("libAliHLTSrc");
    gSystem->Load("libAliHLTMisc");
    gSystem->Load("libAliHLTHough");
    gSystem->Load("libAliHLTComp");
  }
  if (!gROOT->GetClass("AliMonitorClient")) {
    gSystem->Load("libMONITOR");
  }
  new AliMonitorClient;
}
