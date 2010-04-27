// This macro is used to profile the HLT::TRD code
// usage: aliroot aliHLTTRDrun.cxx("/data/sim/")            reconstruct raw ddls, /data/run/ must contain subfolders rawX

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>

#include "TString.h"
#include "TMath.h"

#include "AliHLTSystem.h"
#include "AliLog.h"
#include "AliHLTDataTypes.h"
#include "AliHLTConfiguration.h"

#include <valgrind/callgrind.h>
#include <sys/time.h>
#endif

#include "initGRP.h"
#include "readCDBentry.h"

void aliHLTTRDrun(const TString inDir = gSystem->pwd());
int main(int argc, char** argv)
{
  if(argc==2) aliHLTTRDrun(argv[1]);
  else aliHLTTRDrun();
}

void aliHLTTRDrun(const TString inDir)
{   

  // Is the TRD full?
  Bool_t fullTRD=kFALSE;

  // If not use these SMs:
  Int_t TRDmodules[18] = {0,1,7,8,9,10,17,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  // Use custom arguments for components? i.e.: not reading OCDB arguments
  Bool_t customArgs=kFALSE;

  // Disable HLT flag?
  Bool_t disableHLTflag=kFALSE;



  /////////////////////////////////////

  Int_t usedModules=0;
  if(fullTRD){
    usedModules = 18;
    for(int i=0; i<18; i++)
      TRDmodules[i]=i;
  }else{
#if !defined (__CINT__) || defined (__MAKECINT__)
    std::sort((UInt_t*)TRDmodules, ((UInt_t*)TRDmodules) + 18);
#endif
    for(int i=0; i<18; i++)
      if(TRDmodules[i]>-1)usedModules++;
  }

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

  InitGRP("local://$ALICE_ROOT/OCDB","local:///software/data/ppbench"/*inDir.Data()*/);
  //TString sGeomPath = " -geometry "+inDir+"/geometry.root";
  TString inFolder = inDir+"/raw", inFile = "/TRD_", inExt = ".ddl";
  TString sinput = " -datatype 'DDL_RAW ' 'TRD '";
  TString temp1, temp2;
  Int_t spec, startEvent=0, nEvents=240, ddl;  //KR: start=10, nEvents=30;
  for(Int_t Event=startEvent; Event<(nEvents+startEvent); Event++){
    temp1=inFolder;
    temp1+=Event;
    temp1+=inFile;
    for(int module=0; module<usedModules; module++){
      spec=TMath::Power(2,TRDmodules[module]);
      ddl=TRDmodules[module]+1024;
      temp2=temp1;
      temp2+=ddl;
      temp2+=inExt;
      sinput+=" -dataspec ";
      sinput+=spec;
      sinput+=" -datafile ";
      sinput+=temp2;
    }
    sinput+=" -nextevent";
  }
  printf("%s\n",sinput.Data());
  TString sCFArgs = "";
  TString sTrackerArgs = "";

  if(customArgs || disableHLTflag){
    sCFArgs = readCDBentry("HLT/ConfigTRD/ClusterizerComponent"); //output_percentage 100 -lowflux -experiment -tailcancellation -faststreamer -yPosMethod LUT
    sTrackerArgs = readCDBentry("HLT/ConfigTRD/TrackerV1Component"); //"output_percentage 100 -lowflux -NTimeBins 24";
    sCFArgs += ""; // -processTracklets
    sTrackerArgs += ""; // -highLevelOutput yes -emulateHLToutput no
    if(disableHLTflag){
      sCFArgs +=" -HLTflag no";
      sTrackerArgs +=" -HLTflag no";
    }
  }

  // ======== Configuring chain
  // Chain 1
  AliHLTConfiguration Hsource("Hsource", "FilePublisher", 0, sinput);
  AliHLTConfiguration HClusterizer("HClusterizer", "TRDClusterizer", "Hsource", sCFArgs);
  AliHLTConfiguration HWriterCF("HWriterCF", "FileWriter", "HClusterizer", "-directory output/ -datafile cf.out");
  AliHLTConfiguration HClustMultTrig("HClustMultTrig", "TrdClusterMultiplicityTrigger", "HClusterizer", "-MultiplicityThresh 400");

  AliHLTConfiguration HTracker("HTracker", "TRDTrackerV1", "HClusterizer", sTrackerArgs);
  AliHLTConfiguration HCalib("HCalib", "TRDCalibration", "HTracker", "-pushback-period=10 -TrgStr hi -rejectTrgStr");
  AliHLTConfiguration HWriterCalib("HWriterCalib", "ROOTFileWriter", "HCalib", "-directory output/ -datafile calib.root -concatenate-events -concatenate-blocks -write-all-events");

  AliHLTConfiguration HESDMaker("HESDMaker", "GlobalEsdConverter", "HTracker", "-notree");
  AliHLTConfiguration HTrackMerger("HTrackMerger", "GlobalTrackMerger", "HTracker", "");

  AliHLTConfiguration HClHisto("HClHisto", "TRDClusterHisto", "HClusterizer", "-pushback-period=10");
  AliHLTConfiguration HClWriterHisto("HClWriterHisto", "ROOTFileWriter", "HClHisto", "-directory output/ -datafile clHisto.root -concatenate-events -concatenate-blocks");

  AliHLTConfiguration HTrHisto("HTrHisto", "TRDTrackHisto", "HTracker", "-pushback-period=10");
  AliHLTConfiguration HTrWriterHisto("HTrWriterHisto", "ROOTFileWriter", "HTrHisto", "-directory output/ -datafile trHisto.root -concatenate-events -concatenate-blocks");

  AliHLTConfiguration writerOffConf("esdWriter", "TRDEsdWriter", "HTracker", "-datafile AliHLTTRDESDs.root -concatenate-events -concatenate-blocks");

  gHLT.BuildTaskList("HTracker");

  gHLT.Run(nEvents);
}

