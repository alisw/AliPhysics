#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
//#include "./AliFlatESDEvent.h"
//#include "./AliFlatESDTrack.h"
//#include "../BASE/AliHLTExternalTrackParam.h"
#endif   


void LoadLibs() {

  gSystem->Load("libHLTbase.so");
  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTTPC.so");
  gSystem->Load("libAliHLTITS.so");
  gSystem->Load("libAliHLTGlobal.so");

  return;
}    
