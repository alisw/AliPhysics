// This macro is used to make performace studies with MC info
// usage: aliroot rec-hlt-trd-digits.cxx("/data/run/")    reconstruct local digits file
//    or copy into folder and aliroot rec-hlt-trd.cxx     reconstruct local digits file in pwd

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "TString.h"
#include "TMath.h"

#include "AliHLTSystem.h"
#include "AliLog.h"
#include "AliHLTDataTypes.h"
#include "AliHLTConfiguration.h"

#include <valgrind/callgrind.h>
#include <sys/time.h>
#include "TSystem.h"
#include "AliHLTOfflineInterface.h"
#include "AliRunLoader.h"
#include "AliReconstructor.h"
#endif

#include "initGRP.h"
#include "readCDBentry.h"

int rec_hlt_trd_digits(const TString input = gSystem->pwd());
int main(int argc, char** argv)
{
  if(argc==2) return rec_hlt_trd_digits(argv[1]);
  else return rec_hlt_trd_digits();
}

int rec_hlt_trd_digits(const TString input){

  // Use custom arguments for components?
  Bool_t customArgs=kFALSE;

  // Disable HLT flag?
  Bool_t disableHLTflag=kFALSE;



  ////////////////////////////////

  gSystem->ChangeDirectory(input.Data());

  InitGRP("local://$ALICE_ROOT/OCDB",Form("local://",gSystem->pwd()));
  
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  AliHLTOfflineInterface::SetParamsToComponents(runLoader, NULL);

  AliHLTSystem gHLT;
  AliLog::SetGlobalDebugLevel(1); 
  //AliLog::SetClassDebugLevel("AliTRDrawStream", 11);
  /* enum ETpye {kFatal=0, kError, kWarning, kInfo, kDebug, kMaxType}; */
  gHLT.SetGlobalLoggingLevel((AliHLTComponentLogSeverity)0x7f); 
  /*  enum AliHLTComponentLogSeverity {       
      kHLTLogNone      = 0,
      kHLTLogBenchmark = 0x1,
      kHLTLogDebug     = 0x2,
      kHLTLogInfo      = 0x4,
      kHLTLogWarning   = 0x8,
      kHLTLogError     = 0x10,
      kHLTLogFatal     = 0x20,      
      few important messages not to be filtered out.
      redirected to kHLTLogInfo in AliRoot
      kHLTLogImportant = 0x40,
      special value to enable all messages 
      kHLTLogAll       = 0x7f,
      the default logging filter 
      kHLTLogDefault   = 0x79
      useful           = 0x45
  */

  gHLT.LoadComponentLibraries("libAliHLTUtil.so libAliHLTTRD.so libAliHLTMUON.so libAliHLTGlobal.so libAliHLTTrigger.so");

  // digits publisher
  AliHLTConfiguration pubDigiConf("TRD-DigiP", "AliLoaderPublisher", NULL , "-loader TRDLoader -tree digits -datatype 'ALITREED' 'TRD '");

  TString arg="";
  if(customArgs || disableHLTflag){
    arg = readCDBentry("HLT/ConfigTRD/ClusterizerComponent"); //output_percentage 100 -lowflux -experiment -tailcancellation -faststreamer -yPosMethod LUT
    if(customArgs)
      arg += " -highLevelOutput yes -emulateHLToutput no";
    if(disableHLTflag)
      arg += " -HLTflag no";
  }
  // clusterizer which processes digits
  AliHLTConfiguration cfDigiConf("TRD-DigiCF", "TRDOfflineClusterizer", "TRD-DigiP", arg.Data());

  arg="";
  if(customArgs || disableHLTflag){
    arg = readCDBentry("HLT/ConfigTRD/TrackerV1Component"); //"output_percentage 100 -lowflux -NTimeBins 24";
    if(customArgs)
      arg += " -highLevelOutput yes -emulateHLToutput no";
    if(disableHLTflag)
      arg+=" -HLTflag no";
  }
  // tracker reading the output from the clusterizer which processes the digits
  AliHLTConfiguration trDigiConf("TRD-DigiTR", "TRDOfflineTrackerV1", "TRD-DigiCF", arg.Data());

  // root file writer (with esd friends and MC)
  AliHLTConfiguration writerDigiConf("TRD-DigiEsdFile", "TRDEsdWriter", "TRD-DigiTR", "-concatenate-events -concatenate-blocks");
  
  gHLT.BuildTaskList("TRD-DigiEsdFile");
  gHLT.Run(runLoader->GetNumberOfEvents());
}

