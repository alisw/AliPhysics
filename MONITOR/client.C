#if !defined(__CINT__) || defined(__MAKECINT__)
#include "MONITOR/AliMonitorClient.h"
#endif

void client()
{
  if (!gROOT->GetClass("AliMonitorClient")) {
    gSystem->Load("libMONITOR.so");
  }
  new AliMonitorClient;
}
