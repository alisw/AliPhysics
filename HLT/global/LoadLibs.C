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

  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -g"); 

  gSystem->Load("libHLTbase.so");
  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTTPC.so");
  gSystem->Load("libAliHLTITS.so");
  gSystem->Load("libAliHLTGlobal.so");

  return;
}    
