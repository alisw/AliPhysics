// usage: aliroot rec-hlt-trd.cxx("/data/real_runXXX/raw_input.root")
//    or: aliroot rec-hlt-trd.cxx("/data/sim_run/raw.root")
//    or: aliroot rec-hlt-trd.cxx("/data/sim_run/") *
//    or copy into folder and aliroot rec-hlt-trd.cxx
//
// (*) here sim_run has as subfolders rawX (be sure to have the last "/" !!)

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <iostream>
#include <stdlib.h>
#include <fstream>

#include "TString.h"
#include "TMath.h"

#include "AliHLTSystem.h"
#include "AliHLTPluginBase.h"
#include "AliLog.h"
#include "AliReconstruction.h"
#include "AliHLTDataTypes.h"
#include "AliHLTConfiguration.h"

#include <valgrind/callgrind.h>
#include <sys/time.h>
#include "TSystem.h"

#endif

const Bool_t fullTRD=kTRUE;

void rec_hlt_trd(const char* input ="./raw.root");
int main(int argc, char** argv)
{
  if(argc>1) rec_hlt_trd(argv[1]);
  else rec_hlt_trd();
}

void rec_hlt_trd(const char* input)
{
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  
  // Is the TRD full?
  Bool_t fullTRD=kFALSE;

  // If not use these SMs:
  Int_t TRDmodules[18] = {0,1,7,8,9,10,17,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  // What chains should be run?
  TString chains="TRD-OffEsdFile";

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
 
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }

  TString GRPDir=input;
  GRPDir.Resize(GRPDir.Last('/')+1);
  if(GRPDir.Length()==0)GRPDir=gSystem->pwd();
  printf("GRP dir: %s\n",GRPDir.Data());

  Int_t usedModules=0;
  if(fullTRD){
    usedModules = 18;
    for(int i=0; i<18; i++)
      TRDmodules[i]=i;
  }else{
    for(int i=0; i<18; i++)
      if(TRDmodules[i]>0)usedModules++;
  }

  TString option="libAliHLTUtil.so libAliHLTTRD.so libAliHLTGlobal.so libAliHLTTrigger.so loglevel=0x7f chains=";
  option+=chains;
  TString nextInput, nextOffInput;

  for (int module = 0; module < usedModules; module++) 
    {
      TString arg, publisher, cf, tr, trOff;
      // raw data publisher components
      publisher.Form("TRD-RP_%02d", module);
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TRD ' -dataspec %i -verbose", module+1024, (int)TMath::Power(2,module));
      AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      // Clusterizer
      cf.Form("TRD-CF_%02d", module);
      AliHLTConfiguration cfConf(cf.Data(), "TRDClusterizer", publisher.Data(), "output_percentage 700 -lowflux -simulation -yPosMethod LUT");

      // Tracker
      tr.Form("TRD-TR_%02d", module);
      arg.Form("output_percentage 100 -lowflux -PIDmethod NN");
      AliHLTConfiguration trConf(tr.Data(), "TRDTrackerV1", cf.Data(), arg.Data());
      
      if (nextInput.Length()>0) nextInput+=" ";
      nextInput+=tr;

      trOff.Form("TRD-TROFF_%02d", module);
      AliHLTConfiguration trOffConf(trOff.Data(), "TRDOfflineTrackerV1", cf.Data(), arg.Data());

      if (nextOffInput.Length()>0) nextOffInput+=" ";
      nextOffInput+=trOff;
      
    }

  // calibration
  AliHLTConfiguration calibConf("TRD-Calib", "TRDCalibration", nextInput.Data(), "");
  AliHLTConfiguration writerCalibConf( "TRD-CalibFile", "ROOTFileWriter", "TRD-Calib", "-directory hlt-trd-calib/ -datafile calib.root");

  // esd converter
  AliHLTConfiguration esdConf("TRD-Esd", "GlobalEsdConverter", nextInput.Data(), "-notree");
  
  // root file writer
  AliHLTConfiguration writerConf("TRD-EsdFile", "EsdCollector", "TRD-Esd", "-directory hlt-trd-esd/");

  // root file writer (with esd friends and some day perhaps MC)
  AliHLTConfiguration writerOffConf("TRD-OffEsdFile", "TRDEsdWriter", nextOffInput.Data(), "-datafile AliHLTTRDESDs.root -concatenate-events -concatenate-blocks");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstructio is switched off 
  //
  AliReconstruction rec;
  rec.SetInput(input);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking(":");
  rec.SetLoadAlignFromCDB(0);
  rec.SetFillESD("");
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);
  rec.SetFillTriggerESD(kFALSE);
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");   
  rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",GRPDir.Data()));
  //rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));

  rec.SetOption("HLT", option);
  rec.Run();
}
