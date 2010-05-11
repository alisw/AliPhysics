// This macro is used to simulate the HLT::TRD reconstruction
// usage: aliroot rec-hlt-trd.cxx("/data/run/raw.root")    reconstruct local raw root file (you might add "alien://" to reconstruct remotely)
//    or copy into folder and aliroot rec-hlt-trd.cxx      reconstruct raw.root in pwd
//
// (*1) here /data/run/ must contain subfolders rawX (be sure to have the last "/" !!)

#if !defined (__CINT__) || defined (__MAKECINT__)

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <algorithm>

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
#include "TFile.h"
#include "TGrid.h"
#include "TGridResult.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliExternalTrackParam.h"
#endif

#include "readCDBentry.h"

int rec_hlt_trd(const TString input ="raw.root", TString outPath=gSystem->pwd());
Int_t ExtractRunNumber(const TString str);
int main(int argc, char** argv)
{
  if(argc==2) return rec_hlt_trd(argv[1]);
  else if(argc==3) return rec_hlt_trd(argv[1],argv[2]);
  else return rec_hlt_trd();
}

int rec_hlt_trd(const TString filename, TString outPath)
{
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  
  // What chains should be run? (usually would be: TRD-OffEsdFile)
  TString chains="TRD-OffEsdFile";

  // cosmics or not
  Bool_t bCosmics=kFALSE;

  // look only in data containing TRD triggers?
  Bool_t useOnlyTRDtrigger=kFALSE;

  // Is the TRD full?
  Bool_t fullTRD=kFALSE;

  // If not use these SMs:
  Int_t TRDmodules[18] = {0,1,7,8,9,10,17,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

  // Use custom arguments for components? i.e.: not reading OCDB arguments
  Bool_t customArgs=kFALSE;

  // Disable HLT flag?
  Bool_t disableHLTflag=kFALSE;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // set paths
  //
  TString dataPath = outPath;
  dataPath.Resize(dataPath.Last('/')+1);
  if(!dataPath.Length()){
    dataPath=gSystem->pwd();
    dataPath+="/";
    outPath.Prepend(dataPath);
  }
  
  TString extInput = filename;
  dataPath = filename;
  dataPath.Resize(dataPath.Last('/')+1);
  if(!dataPath.Length()){
    dataPath=gSystem->pwd();
    dataPath+="/";
    extInput.Prepend(dataPath);
  }

  printf("File path %s\n",dataPath.Data());
  printf("Processing file %s\n",extInput.Data());
  printf("Output to %s\n",outPath.Data());

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return -1;
  }

  Bool_t bRealData=kFALSE;
  if(filename.Contains(".root") && !filename.Contains("raw.root")){
    bRealData = kTRUE;
    printf("processing real data\n");
  }else{
    bRealData = kFALSE;
    printf("processing simulated data\n");
  }

  if(filename.Contains("alien://") || bRealData){
    TGrid::Connect("alien://");
  }

  if(filename.Contains(".root") && !TFile::Open(filename))return -1;

  gSystem->mkdir(outPath.Data());
  gSystem->ChangeDirectory(outPath.Data());

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
  
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

  TString option="libAliHLTUtil.so libAliHLTTRD.so libAliHLTMUON.so libAliHLTGlobal.so libAliHLTTrigger.so loglevel=0x7f chains=";
  option+=chains;
  TString afterTr, afterTrOff, afterCf, afterCal;

  for (int module = 0; module < usedModules; module++) 
    {
      TString arg, publisher, cf, tr, trOff, cal;
      // raw data publisher components
      publisher.Form("TRD-RP_%02d", TRDmodules[module]);
      arg.Form("-minid %d -datatype 'DDL_RAW ' 'TRD ' -dataspec %i -verbose", TRDmodules[module]+1024, (int)TMath::Power(2,TRDmodules[module]));
      AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
      
      // Clusterizer
      arg = "";
      if(customArgs || disableHLTflag){
	arg = readCDBentry("HLT/ConfigTRD/ClusterizerComponent"); //output_percentage 100 -lowflux -experiment -tailcancellation -faststreamer -yPosMethod LUT
	if(customArgs)
	  arg += "";
	if(disableHLTflag)
	  arg += " -HLTflag no";
      }

      cf.Form("TRD-CF_%02d", TRDmodules[module]);
      AliHLTConfiguration cfConf(cf.Data(), "TRDClusterizer", publisher.Data(), arg.Data());

      if (afterCf.Length()>0) afterCf+=" ";
      afterCf+=cf;

      // Tracker
      arg="";
      if(customArgs || disableHLTflag){
	arg = readCDBentry("HLT/ConfigTRD/TrackerV1Component"); //"output_percentage 100 -lowflux -NTimeBins 24";
	if(customArgs)
	  arg += "";
	if(disableHLTflag)
	  arg += " -HLTflag no";
      }

      tr.Form("TRD-TR_%02d", TRDmodules[module]);
      AliHLTConfiguration trConf(tr.Data(), "TRDTrackerV1", cf.Data(), arg.Data());
      
      if (afterTr.Length()>0) afterTr+=" ";
      afterTr+=tr;

      // Offline Tracker (for debug purposes only)
      arg = readCDBentry("HLT/ConfigTRD/TrackerV1Component"); //"output_percentage 100 -lowflux -NTimeBins 24";
      if(customArgs)
	arg += " -highLevelOutput yes -emulateHLToutput no";
      if(disableHLTflag)
	arg+=" -HLTflag no";

      trOff.Form("TRD-TROFF_%02d", TRDmodules[module]);
      AliHLTConfiguration trOffConf(trOff.Data(), "TRDOfflineTrackerV1", cf.Data(), arg.Data());

      if (afterTrOff.Length()>0) afterTrOff+=" ";
      afterTrOff+=trOff;

      // new SM wise calibration
      arg="-takeAllEvents";

      cal.Form("TRD-CalHist_%02d", TRDmodules[module]);
      AliHLTConfiguration calConf(cal.Data(), "TRDCalibHisto", tr.Data(), arg.Data());
      
      if (afterCal.Length()>0) afterCal+=" ";
      afterCal+=cal;
      
    }

  // cluster histogramm
  AliHLTConfiguration histoConf("TRD-ClHisto", "TRDClusterHisto", afterCf.Data(), "");
  AliHLTConfiguration writerHistoConf( "TRD-ClHistoFile", "ROOTFileWriter", "TRD-ClHisto", "-directory hlt-trd-histo/ -datafile histo.root -concatenate-events -concatenate-blocks");

  // new calibration (SM wise)
  AliHLTConfiguration calibFitConf("TRD-CalibFit", "TRDCalibFit", afterCal.Data(), "");
  AliHLTConfiguration writerCalibFitConf( "TRD-CalibFitFile", "ROOTFileWriter", "TRD-CalibFit", "-directory hlt-trd-calib/ -datafile calibFit.root -concatenate-events -concatenate-blocks -write-all-events");

  // old calibration (you may use tr or trOff here)
  AliHLTConfiguration calibConf("TRD-Calib", "TRDCalibration", afterTr.Data(), "-takeAllEvents");
  AliHLTConfiguration writerCalibConf( "TRD-CalibFile", "ROOTFileWriter", "TRD-Calib", "-directory hlt-trd-calib/ -datafile calib.root -concatenate-events -concatenate-blocks -write-all-events");

  // esd converter
  AliHLTConfiguration esdConf("TRD-Esd", "GlobalEsdConverter", afterTr.Data(), "-notree");
  
  // root file writer
  AliHLTConfiguration writerConf("TRD-EsdFile", "EsdCollector", "TRD-Esd", "-directory hlt-trd-esd/");

  // root file writer (with esd friends) (you may use tr or trOff here)
  AliHLTConfiguration writerOffConf("TRD-OffEsdFile", "TRDEsdWriter", afterTr.Data(), "-concatenate-events -concatenate-blocks");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Init CDBManager and trigger
  //
  AliCDBManager * man = AliCDBManager::Instance();
  Int_t run = 0;
  if(bRealData){
    man->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB");
    //man->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB?cacheFold=/lustre/alice/local/alice/data/2009/OCDB"); 
    //man->SetSpecificStorage("GRP/GRP/Data","alien://folder=/alice/data/2009/OCDB?cacheFold=/lustre/alice/local/alice/data/2009/OCDB");
    man->SetSpecificStorage("HLT/*","local://$ALICE_ROOT/OCDB");
    run = ExtractRunNumber(filename);
  }else{
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB"); 
    man->SetSpecificStorage("GRP/GRP/Data", Form("local://%s",dataPath.Data()));
  }
  man->SetRun(run);

  if(bCosmics){
    // no magnetic field
    AliExternalTrackParam::SetMostProbablePt(8.);
  }

  // Find TRD triggers
  TString filestring = filename;
  if(useOnlyTRDtrigger){
    AliCDBEntry *grp_ctp = man->Get("GRP/CTP/Config");
    AliTriggerConfiguration *trg_conf = (AliTriggerConfiguration *)grp_ctp->GetObject();
    trg_conf->Print();
    TObjArray trg_masks = trg_conf->GetClasses(); // Reference!!!
    std::vector<unsigned char> triggerconfs;
    for(Int_t iobj = 0; iobj < trg_masks.GetEntriesFast(); iobj++){
      
      AliTriggerClass * trg_class = (AliTriggerClass*)trg_masks.UncheckedAt(iobj);
      AliTriggerCluster * trg_clust = (AliTriggerCluster *)trg_class->GetCluster();
      
      printf("ioj[%d]\n", iobj); trg_class->Print(0x0);

      // cosmic run 2009                                                                                                                                                                                       
      // if(TString(trg_class->GetName()).Contains("TRD")){                                                                                                                                                    
      //        triggerconfs.push_back(trg_class->GetMask());                                                                                                                                                  
      // }                                                                                                                                                                                                     

      // pp run 2009                                                                                                                                                                                           
      if(TString(trg_class->GetName()).Contains("CINT1B-ABCE-NOPF-ALL")){
        triggerconfs.push_back(trg_class->GetMask());
      }

    }
    
    Int_t itrg = 0;
    printf("Number of Trigger Clusters including TRD: %d\n", (Int_t)triggerconfs.size());
    for(std::vector<unsigned char>::iterator it = triggerconfs.begin(); it < triggerconfs.end(); it++)
      printf("Trigger Mask %d for TRD: %d\n", itrg++, *it);
    filestring += "?EventType=7";
    char triggerbuf[256];
    Int_t triggerval = 0;
    for(std::vector<unsigned char>::iterator it = triggerconfs.begin(); it < triggerconfs.end(); it++)
      triggerval += *it;
    sprintf(triggerbuf, "?Trigger=%d", triggerval);
    filestring += triggerbuf; // This line does the trigger selection. It has to be uncommented if one wants to apply trigger selection
  }
  printf("Filename: %s\n", filestring.Data());

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init and run the reconstruction
  //
  AliReconstruction rec;
  rec.SetInput(filestring.Data());
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking(":");
  rec.SetFillESD("");
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);
  rec.SetFillTriggerESD(kFALSE);

  rec.SetOption("HLT", option);
  rec.Run();

  return 0;
}

Int_t ExtractRunNumber(const TString str){
  TObjArray *ptoks = (TObjArray *)str.Tokenize("?");
  TString path = ((TObjString *)ptoks->UncheckedAt(0))->String();
  TObjArray *toks = (TObjArray *)path.Tokenize("/");
  TString fname = ((TObjString *)(toks->UncheckedAt(toks->GetEntriesFast() - 1)))->String();
  TString rstr = fname(2,9);
  return rstr.Atoi();
}
