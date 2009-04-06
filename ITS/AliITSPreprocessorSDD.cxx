/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class for SDD preprocessing             //
// Origin: E.Crescio, Torino, crescio@to.infn.it                 //
//         F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSPreprocessorSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDCSAnalyzerSDD.h"
#include "AliITSHLTforSDD.h"
#include "AliShuttleInterface.h"
#include "AliCDBEntry.h"
#include "AliCDBMetaData.h"
#include "TObjArray.h"
#include "AliLog.h"
#include <TObjString.h>
#include <TSystem.h>
#include <TList.h>

const TString AliITSPreprocessorSDD::fgkNameHistoPedestals = "hpedestal";
const TString AliITSPreprocessorSDD::fgkNameHistoNoise = "hnoise";
ClassImp(AliITSPreprocessorSDD)

//______________________________________________________________________
AliITSPreprocessorSDD::AliITSPreprocessorSDD( AliShuttleInterface* shuttle):
    AliPreprocessor("SDD", shuttle)
{
  // constructor
  AddRunType("PULSER");
  AddRunType("INJECTOR");
  AddRunType("PHYSICS"); 
}

//______________________________________________________________________
UInt_t AliITSPreprocessorSDD::Process(TMap* dcsAliasMap){


  // Get DDL map from OCDB
  AliCDBEntry* entry = GetFromOCDB("Calib", "DDLMapSDD");
  if(!entry){
    Log("DDL map file not found in OCDB.");  
    return 2;
  }
  AliITSDDLModuleMapSDD* ddlmap = (AliITSDDLModuleMapSDD*)entry->GetObject();
  if(!ddlmap){ 
    Log("AliITSDDLModuleMapSDD object not in file.");
    return 2;
  }
  ddlmap->PrintDDLMap();

  //preprocessing
  TString runType = GetRunType();
  Int_t retcode=0;

  if (runType == "PULSER"){
    Log("Process FXS files from PULSER RUN");
    retcode=ProcessPulser(ddlmap);
  }else if(runType== "INJECTOR"){
    Log("Process FXS files from INJECTOR RUN");
    retcode=ProcessInjector(ddlmap);
  }else if(runType== "PHYSICS"){
    retcode=ProcessPhysics();
  }
  if(retcode!=0) return retcode;

  Log("Process DCS data");
  Bool_t retcodedcs =ProcessDCSDataPoints(dcsAliasMap);
  if(retcodedcs) return 0; 
  else return 1;           

}
//______________________________________________________________________
UInt_t AliITSPreprocessorSDD::ProcessPhysics(){
  // Get the HLT status for the PHYSICS run 
  // needeed to define the raw data format  

  AliITSHLTforSDD *hltSDD=new AliITSHLTforSDD();
  TString hltMode = GetRunParameter("HLTmode");
  TSubString firstChar = hltMode(0,1);
  if (firstChar == "C") hltSDD->SetHLTmodeC(kTRUE);
  else hltSDD->SetHLTmodeC(kFALSE);

  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("HLT mode C flag for PHYSICS run");
  Bool_t retCode = Store("Calib","HLTforSDD",hltSDD,md);
  if(retCode) return 0;
  else return 1;
}
//______________________________________________________________________
UInt_t AliITSPreprocessorSDD::ProcessPulser(AliITSDDLModuleMapSDD* ddlmap){
  // Process FXS files from PULSER run (baseline, noise, gain)
  // returns 0 in case of success, 
  //         1 in case of storage error, 
  //         2 in case of error with input files
  TObjArray calSDD(kNumberOfSDD);
  calSDD.SetOwner(kFALSE);

  Float_t baseline,rawnoise,cmn,corn,gain;
  Int_t isgoodan,i,im,is,isgoodmod,basmin,basoff;
  Int_t th,tl;
  Int_t numOfBadChannels[kNumberOfSDD];
  
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
    TString command;
    command.Form("tar -xf %s",tarName.Data());
    gSystem->Exec(command);
    ind++;
  }
  delete sourceList;
  
  for(Int_t iddl=0;iddl<kNumberOfDDL;iddl++){
    for(Int_t imod=0;imod<kModulesPerDDL;imod++){
      Int_t modID=ddlmap->GetModuleNumber(iddl,imod);
      if(modID==-1) continue;
      modID-=240; // to have SDD modules numbering from 0 to 260
      AliITSCalibrationSDD *cal = new AliITSCalibrationSDD("simulated");
      cal->SetAMAt20MHz();            // for runs > 51275 with clock at 20 MHz
      numOfBadChannels[modID]=0;
      Int_t badch[kNumberOfChannels];
      for(Int_t isid=0;isid<=1;isid++){
	TString inpFileName;
	inpFileName.Form("./SDDbase_ddl%02dc%02d_sid%d.data",iddl,imod,isid);

	FILE* basFil = fopen(inpFileName,"read");
	if (basFil == 0) {
	  Log(Form("File %s not found.",inpFileName.Data()));
	  cal->SetBad();
	  continue;
	}
	fscanf(basFil,"%d %d %d\n",&im,&is,&isgoodmod);
	if(!isgoodmod) cal->SetBad();
	fscanf(basFil,"%d\n",&th);
	fscanf(basFil,"%d\n",&tl);
	cal->SetZSLowThreshold(isid,tl);
	cal->SetZSHighThreshold(isid,th);
	for(Int_t ian=0;ian<(kNumberOfChannels/2);ian++){
	  fscanf(basFil,"%d %d %f %d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn,&gain);
	  Int_t ich=ian;
	  if(isid==1) ich+=256;
	  if(!isgoodan){ 
	    Int_t ibad=numOfBadChannels[modID];
	    badch[ibad]=ich;
	    numOfBadChannels[modID]++;
	  }
	  cal->SetBaseline(ich,baseline-basoff);
	  cal->SetNoiseAfterElectronics(ich,rawnoise);
	  cal->SetGain(ich,gain);
	}
	cal->SetDeadChannels(numOfBadChannels[modID]);
	for(Int_t ibad=0;ibad<numOfBadChannels[modID];ibad++){
	  cal->SetBadChannel(ibad,badch[ibad]);
	}
	fclose(basFil);
      }
      Log(Form("Put calib obj for module %d (DDL %d  Carlos %d)",modID,iddl,imod));
      calSDD.AddAt(cal,modID);
    }
  }
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("AliITSCalibrationSDD from PEDESTAL+PULSER runs");
  Bool_t retCode = Store("Calib","CalibSDD",&calSDD,md, 0, kTRUE);
  if(retCode) return 0;
  else return 1;
}
//______________________________________________________________________
UInt_t AliITSPreprocessorSDD::ProcessInjector(AliITSDDLModuleMapSDD* ddlmap){
  // Process FXS files from injector events (INJECTOR or PHYSICS runs)
  // returns 0 in case of success, 
  //         1 in case of storage error, 
  //         2 in case of error with input files
  TObjArray vdrift(2*kNumberOfSDD);
  vdrift.SetOwner(kFALSE);
  Int_t evNumb,polDeg; 
  UInt_t timeStamp;
  Bool_t modSet[2*kNumberOfSDD]; // flag modules with good inj.
  for(Int_t ihyb=0; ihyb<2*kNumberOfSDD; ihyb++) modSet[ihyb]=0;
  Double_t nPt = 0;

  Double_t param[4];    // parameters of poly fit
  Double_t minValP0=4.; // min value for param[0]
  Double_t maxValP0=9.; // max value for param[0]
  Double_t minValP1=0.; // min value for param[1]
  Double_t aveCoef[4]={0.,0.,0.,0.};  // average param for good mod.
  Double_t defCoef[4]={6.53227,0.00128941,-5.14493e-06,0};  // default values for param
  Float_t auxP;

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
    TString command;
    command.Form("tar -xf %s",tarName.Data());
    gSystem->Exec(command);
    ind++;
  }
  delete sourceList;


  for(Int_t iddl=0;iddl<kNumberOfDDL;iddl++){
    for(Int_t imod=0;imod<kModulesPerDDL;imod++){
      Int_t modID=ddlmap->GetModuleNumber(iddl,imod);
      if(modID==-1) continue;
      modID-=240; // to have SDD modules numbering from 0 to 260
      for(Int_t isid=0;isid<=1;isid++){
	AliITSDriftSpeedArraySDD *arr=new AliITSDriftSpeedArraySDD();
	TString inpFileName;
	inpFileName.Form("./SDDinj_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	FILE* injFil = fopen(inpFileName,"read");
	if (injFil == 0) {
	  Log(Form("File %s not found.",inpFileName.Data()));
	  AliITSDriftSpeedSDD *dsp=new AliITSDriftSpeedSDD();
	  arr->AddDriftSpeed(dsp);
	  vdrift.AddAt(arr,2*modID+isid);
	  continue;
	}
	fscanf(injFil,"%d",&polDeg);
	while (!feof(injFil)){
	  fscanf(injFil,"%d %u ",&evNumb,&timeStamp);
	  if(feof(injFil)) break;
	  for(Int_t ic=0;ic<4;ic++){ 
	    fscanf(injFil,"%f ",&auxP);
	    param[ic]=auxP;
	  }

	  if(param[0]>minValP0 && param[0]<maxValP0 && param[1]>minValP1){
	    for(Int_t ic=0;ic<4;ic++) aveCoef[ic]+=param[ic];
	    nPt++;
	    AliITSDriftSpeedSDD *dsp=new AliITSDriftSpeedSDD(evNumb,timeStamp,polDeg,param);
	    arr->AddDriftSpeed(dsp);
	    modSet[2*modID+isid]=1;
	  }
	}
	Log(Form("Put calib obj for hybrid %d (DDL %d  Carlos %d)",2*modID+isid,iddl,imod));
	if(modSet[2*modID+isid]) vdrift.AddAt(arr,2*modID+isid);
      }
    }
  }

  // set drift speed for modules with bad injectors
  for(Int_t ic=0;ic<4;ic++){ 
    if(nPt>0) aveCoef[ic]/=nPt; // mean parameters
    else aveCoef[ic]=defCoef[ic]; // default parameters
  }
  AliITSDriftSpeedSDD *avdsp=new AliITSDriftSpeedSDD(evNumb,timeStamp,polDeg,aveCoef);

  for(Int_t ihyb=0; ihyb<2*kNumberOfSDD; ihyb++){
    if(modSet[ihyb]==0){ 
      AliWarning(Form("No good injector events for mod. %d --> use average values",ihyb/2));
      AliITSDriftSpeedArraySDD *arr=new AliITSDriftSpeedArraySDD();
      arr->AddDriftSpeed(avdsp);
      vdrift.AddAt(arr,ihyb);
    }
  }


  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("AliITSDriftSpeedSDD from injector events");
  Bool_t retCode = Store("Calib","DriftSpeedSDD",&vdrift,md,0, kTRUE);
  if(retCode) return 0;
  else return 1;
}
//______________________________________________________________________
Bool_t AliITSPreprocessorSDD::ProcessDCSDataPoints(TMap* dcsAliasMap){
  // Process DCS data
  AliITSDCSAnalyzerSDD *dcs=new AliITSDCSAnalyzerSDD();
  dcs->AnalyzeData(dcsAliasMap);
  TObjArray refDCS(kNumberOfSDD);
  refDCS.SetOwner(kFALSE);
  for(Int_t imod=0;imod<kNumberOfSDD;imod++){
    AliITSDCSDataSDD *dcsdata=dcs->GetDCSData(imod);
    refDCS.Add(dcsdata);
  }    
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("AliITSDCSDataSDD objects from DCS DB");
  Bool_t retCode = StoreReferenceData("DCS","DataSDD",&refDCS,md);
  return retCode;
}

