/////////////////////////////////////////
// Class for SDD digits preprocessing  //
//                                     //
//                                     //
////////////////////////////////////////

/* $Id$ */

#include "AliITSPreprocessorSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDCSAnalyzerSDD.h"
#include "AliShuttleInterface.h"
#include "AliCDBMetaData.h"
#include "TObjArray.h"
#include "AliLog.h"
#include <TObjString.h>
#include <TSystem.h>
#include <TList.h>

const Int_t AliITSPreprocessorSDD::fgkNumberOfSDD = 260;
const Int_t AliITSPreprocessorSDD::fgkNumberOfChannels = 512;
const TString AliITSPreprocessorSDD::fgkNameHistoPedestals = "hpedestal";
const TString AliITSPreprocessorSDD::fgkNameHistoNoise = "hnoise";
ClassImp(AliITSPreprocessorSDD)


UInt_t AliITSPreprocessorSDD::Process(TMap* dcsAliasMap){

  //preprocessing. 

  TString runType = GetRunType();
  Bool_t retcode;
  Char_t command[100];
  Char_t inpFileName[100];
  AliCDBMetaData *md1= new AliCDBMetaData();
  md1->SetResponsible("Elisabetta Crescio, Francesco Prino");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 30 November 2007"); //root version
  md1->SetComment("This is a test");

  if (runType == "PULSER_RUN"){
    TObjArray calSDD(fgkNumberOfSDD);
    calSDD.SetOwner(kFALSE);
    Float_t baseline,rawnoise,cmn,corn,gain;
    Int_t isgoodan,i,im,is,isgoodmod,basmin,basoff;
    Int_t numOfBadChannels[fgkNumberOfSDD];
  
    TList* sourceList = GetFileSources(kDAQ, "SDD_Calib");
    if (!sourceList){ 
      Log("Error: no sources found for SDD_Calib");
      return 2;
    }

    Int_t ind = 0;
    while (sourceList->At(ind)!=NULL) {
      TObjString* tarId = (TObjString*) sourceList->At(ind);
      TString tarName = GetFile(kDAQ, "SDD_Calib", tarId->GetString().Data());
      if(tarName.Length()==0){
	Log(Form("Baseline tar file from source %d not found.",ind));
	return 2;
      }
      sprintf(command,"tar -xf %s",tarName.Data());
      gSystem->Exec(command);
      ind++;
    }
    delete sourceList;

    for(Int_t imod=0;imod<fgkNumberOfSDD;imod++){
      AliITSCalibrationSDD *cal = new AliITSCalibrationSDD("simulated");
      numOfBadChannels[imod]=0;
      Int_t badch[fgkNumberOfChannels];
      for(Int_t isid=0;isid<=1;isid++){
       sprintf(inpFileName,"./SDDbase_mod%03d_sid%d.data",imod,isid);
       FILE* basFil = fopen(inpFileName,"read");
       if (basFil == 0) {
	 Log(Form("File %s not found.",inpFileName));
	 return 2;
       }      
       fscanf(basFil,"%d %d %d\n",&im,&is,&isgoodmod);
       if(!isgoodmod) cal->SetDead();
       for(Int_t ian=0;ian<(fgkNumberOfChannels/2);ian++){
	 fscanf(basFil,"%d %d %f %d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn,&gain);
	 Int_t ich=ian;
	 if(isid==1) ich+=256;
	 if(!isgoodan){ 
	   Int_t ibad=numOfBadChannels[imod];
	   numOfBadChannels[imod]++;
	   badch[ibad]=ich;
	 }
	 cal->SetBaseline(ich,baseline);
	 cal->SetNoiseAfterElectronics(ich,rawnoise);
	 Int_t iChip=cal->GetChip(ich);
	 Int_t iChInChip=cal->GetChipChannel(ich);
	 cal->SetGain(gain,isid,iChip,iChInChip);
       }
       cal->SetDeadChannels(numOfBadChannels[imod]);
       for(Int_t ibad=0;ibad<numOfBadChannels[imod];ibad++){
	 cal->SetBadChannel(ibad,badch[ibad]);
       }
       fclose(basFil);
      }
      calSDD.Add(cal);
    }
    md1->SetObjectClassName("AliITSCalibration");
    retcode = Store("Calib","CalibSDD",&calSDD,md1, 0, kTRUE);
  }else if(runType == "PHYSICS"){

    TObjArray vdrift(2*fgkNumberOfSDD);
    vdrift.SetOwner(kFALSE);
    Int_t evNumb,polDeg; 
    UInt_t timeStamp;
    Float_t param[4];

    TList* sourceList = GetFileSources(kDAQ, "SDD_Injec");
    if (!sourceList){ 
      Log("Error: no sources found for SDD_Injec");
      return 2;
    }
    Int_t ind = 0;
    while (sourceList->At(ind)!=NULL) {
      TObjString* tarId = (TObjString*) sourceList->At(ind);
      TString tarName = GetFile(kDAQ, "SDD_Injec", tarId->GetString().Data());
      if(tarName.Length()==0){
	Log(Form("Injector tar file from source %d not found.",ind));
	return 2;
      }
      sprintf(command,"tar -xf %s",tarName.Data());
      gSystem->Exec(command);
      ind++;
    }
    delete sourceList;

    for(Int_t imod=0;imod<fgkNumberOfSDD;imod++){
      for(Int_t isid=0;isid<=1;isid++){
	AliITSDriftSpeedArraySDD *arr=new AliITSDriftSpeedArraySDD();
	sprintf(inpFileName,"./SDDinj_mod%03d_sid%d.data",imod,isid);
	FILE* injFil = fopen(inpFileName,"read");
	if (injFil == 0) {
	  Log(Form("File %s not found.",inpFileName));
	  return 2;
	}
	fscanf(injFil,"%d",&polDeg);
	while (!feof(injFil)){
	  fscanf(injFil,"%d %d",&evNumb,&timeStamp);
	  if(feof(injFil)) break;
	  for(Int_t ic=0;ic<4;ic++) fscanf(injFil,"%f",&param[ic]);
	  AliITSDriftSpeedSDD *dsp=new AliITSDriftSpeedSDD(evNumb,timeStamp,polDeg,param);
	  arr->AddDriftSpeed(dsp);
	}
	vdrift.Add(arr);
      }
    }
    md1->SetObjectClassName("AliITSDriftSpeedArraySDD");
    retcode = Store("Calib","DriftSpeedSDD",&vdrift,md1,0, kTRUE);    
  }else{
    // do nothing for other run types
    retcode=1;
  }
  if(retcode){
    // process DCS data
    AliITSDCSAnalyzerSDD *dcs=new AliITSDCSAnalyzerSDD();
    dcs->AnalyzeData(dcsAliasMap);
    TObjArray refDCS(fgkNumberOfSDD);
    refDCS.SetOwner(kFALSE);
    for(Int_t imod=0;imod<fgkNumberOfSDD;imod++){
      AliITSDCSDataSDD *dcsdata=dcs->GetDCSData(imod);
      refDCS.Add(dcsdata);
    }
    
    AliCDBMetaData *mddcs= new AliCDBMetaData();
    mddcs->SetResponsible("Francesco Prino");
    mddcs->SetBeamPeriod(0);
    mddcs->SetAliRootVersion("head 18 December 2007"); //root version
    mddcs->SetComment("This is a test");
    mddcs->SetObjectClassName("AliITSDCSDataSDD");
    Int_t retcodedcs = StoreReferenceData("DCS","DataSDD",&refDCS,mddcs);
    
    if(retcodedcs) return 0;
  }
  return 1;
}
