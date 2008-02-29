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

  //preprocessing
  TString runType = GetRunType();
  Bool_t retcode;

  if (runType == "PULSER"){
    retcode=ProcessPulser(ddlmap);
  }else if(runType== "INJECTOR"){
    retcode=ProcessInjector(ddlmap);
  }else{
    // do nothing for other run types
    retcode=1;
  }
  if(retcode){
    Bool_t retcodedcs =ProcessDCSDataPoints(dcsAliasMap);
    if(retcodedcs) return 0;
  }
  return 1;
}
//______________________________________________________________________
Bool_t AliITSPreprocessorSDD::ProcessPulser(AliITSDDLModuleMapSDD* ddlmap){
  // Process FXS files from PULSER run (baseline, noise, gain)
  TObjArray calSDD(kNumberOfSDD);
  calSDD.SetOwner(kFALSE);
  Char_t command[100];
  Char_t inpFileName[100];
  Float_t baseline,rawnoise,cmn,corn,gain;
  Int_t isgoodan,i,im,is,isgoodmod,basmin,basoff;
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
    sprintf(command,"tar -xf %s",tarName.Data());
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
      numOfBadChannels[modID]=0;
      Int_t badch[kNumberOfChannels];
      for(Int_t isid=0;isid<=1;isid++){
	sprintf(inpFileName,"./SDDbase_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	FILE* basFil = fopen(inpFileName,"read");
	if (basFil == 0) {
	  Log(Form("File %s not found.",inpFileName));
	  cal->SetBad();
	  continue;
	}
	fscanf(basFil,"%d %d %d\n",&im,&is,&isgoodmod);
	if(!isgoodmod) cal->SetBad();
	for(Int_t ian=0;ian<(kNumberOfChannels/2);ian++){
	  fscanf(basFil,"%d %d %f %d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn,&gain);
	  Int_t ich=ian;
	  if(isid==1) ich+=256;
	  if(!isgoodan){ 
	    Int_t ibad=numOfBadChannels[modID];
	    badch[ibad]=ich;
	    numOfBadChannels[modID]++;
	  }
	  cal->SetBaseline(ich,baseline);
	  cal->SetNoiseAfterElectronics(ich,rawnoise);
	  cal->SetGain(ich,gain);
	}
	cal->SetDeadChannels(numOfBadChannels[modID]);
	for(Int_t ibad=0;ibad<numOfBadChannels[modID];ibad++){
	  cal->SetBadChannel(ibad,badch[ibad]);
	}
	fclose(basFil);
      }
      calSDD.AddAt(cal,modID);
    }
  }
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("AliITSCalibrationSDD from PEDESTAL+PULSER runs");
  Bool_t retCode = Store("Calib","CalibSDD",&calSDD,md, 0, kTRUE);
  return retCode;
}
//______________________________________________________________________
Bool_t AliITSPreprocessorSDD::ProcessInjector(AliITSDDLModuleMapSDD* ddlmap){
  // Process FXS files from injector events (INJECTOR or PHYSICS runs)
  TObjArray vdrift(2*kNumberOfSDD);
  vdrift.SetOwner(kFALSE);
  Char_t command[100];
  Char_t inpFileName[100];
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
  
  for(Int_t iddl=0;iddl<kNumberOfDDL;iddl++){
    for(Int_t imod=0;imod<kModulesPerDDL;imod++){
      Int_t modID=ddlmap->GetModuleNumber(iddl,imod);
      if(modID==-1) continue;
      modID-=240; // to have SDD modules numbering from 0 to 260
      for(Int_t isid=0;isid<=1;isid++){
	AliITSDriftSpeedArraySDD *arr=new AliITSDriftSpeedArraySDD();
	sprintf(inpFileName,"./SDDinj_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	FILE* injFil = fopen(inpFileName,"read");
	if (injFil == 0) {
	  Log(Form("File %s not found.",inpFileName));
	  AliITSDriftSpeedSDD *dsp=new AliITSDriftSpeedSDD();
	  arr->AddDriftSpeed(dsp);
	  vdrift.AddAt(arr,2*modID+isid);
	  continue;
	}
	fscanf(injFil,"%d",&polDeg);
	while (!feof(injFil)){
	  fscanf(injFil,"%d %d",&evNumb,&timeStamp);
	  if(feof(injFil)) break;
	  for(Int_t ic=0;ic<4;ic++) fscanf(injFil,"%f",&param[ic]);
	  AliITSDriftSpeedSDD *dsp=new AliITSDriftSpeedSDD(evNumb,timeStamp,polDeg,param);
	  arr->AddDriftSpeed(dsp);
	}
	vdrift.AddAt(arr,2*modID+isid);
      }
    }
  }
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("AliITSDriftSpeedSDD from injector events");
  Bool_t retCode = Store("Calib","DriftSpeedSDD",&vdrift,md,0, kTRUE);
  return retCode;
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
