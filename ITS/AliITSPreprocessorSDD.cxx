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
#include <TF1.h>
#include <TH1D.h>

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
  }
  if(retcode!=0) return retcode;

  Log("Process DCS data");
  Bool_t retcodedcs =ProcessDCSDataPoints(dcsAliasMap);
  if(retcodedcs) return 0; 
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
  // Read ADC sampling frequency from fee.conf
  Int_t amSamplFreq=40;
  Int_t retfscf;
  FILE* feefil=fopen("fee.conf","r");
  if(feefil){
    retfscf=fscanf(feefil,"%d \n",&amSamplFreq);
    fclose(feefil);
    Log(Form("AM sampling frequency = %d MHz",amSamplFreq));
  }else{
    Log("File fee.conf not found. AM sampling set at 40 MHz by default");
  }
  
  for(Int_t iddl=0;iddl<kNumberOfDDL;iddl++){
    for(Int_t imod=0;imod<kModulesPerDDL;imod++){
      Int_t modID=ddlmap->GetModuleNumber(iddl,imod);
      if(modID==-1) continue;
      modID-=240; // to have SDD modules numbering from 0 to 260
      AliITSCalibrationSDD *cal = new AliITSCalibrationSDD("simulated");
      if(amSamplFreq!=40) cal->SetAMAt20MHz();
      numOfBadChannels[modID]=0;
      Int_t badch[kNumberOfChannels];
      Bool_t sid0ok=kTRUE;
      Bool_t sid1ok=kTRUE;
      for(Int_t isid=0;isid<=1;isid++){
	TString inpFileName;
	inpFileName.Form("./SDDbase_ddl%02dc%02d_sid%d.data",iddl,imod,isid);

	FILE* basFil = fopen(inpFileName,"read");
	if (basFil == 0) {
	  Log(Form("File %s not found.",inpFileName.Data()));
	  if(isid==0){
	    sid0ok=kFALSE;
	    for(Int_t iChip=0; iChip<4; iChip++) cal->SetChipBad(iChip);
	    cal->SetDeadChannels(cal->GetDeadChannels()+256);
	    for(Int_t iAnode=0; iAnode<256; iAnode++) cal->SetBadChannel(iAnode,iAnode);
	  }else{
	    sid1ok=kFALSE;
	    for(Int_t iChip=4; iChip<8; iChip++) cal->SetChipBad(iChip);
	    cal->SetDeadChannels(cal->GetDeadChannels()+256);
	    for(Int_t iAnode=0; iAnode<256; iAnode++) cal->SetBadChannel(iAnode,iAnode+256);
	  }
	  continue;
	}

	retfscf=fscanf(basFil,"%d %d %d\n",&im,&is,&isgoodmod);
	if(!isgoodmod){
	  if(isid==0){
	    sid0ok=kFALSE;
	    for(Int_t iChip=0; iChip<4; iChip++) cal->SetChipBad(iChip);
	  }else{
	    sid1ok=kFALSE;
	    for(Int_t iChip=4; iChip<8; iChip++) cal->SetChipBad(iChip);
	  }
	}
	retfscf=fscanf(basFil,"%d\n",&th);
	retfscf=fscanf(basFil,"%d\n",&tl);
	cal->SetZSLowThreshold(isid,tl);
	cal->SetZSHighThreshold(isid,th);
	for(Int_t ian=0;ian<(kNumberOfChannels/2);ian++){
	  retfscf=fscanf(basFil,"%d %d %f %d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&basmin,&basoff,&rawnoise,&cmn,&corn,&gain);
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
      if(!sid0ok && !sid1ok) cal->SetBad();
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
  UInt_t timeStamp,statusInj;
  Bool_t modSet[2*kNumberOfSDD]; // flag modules with good inj.
  for(Int_t ihyb=0; ihyb<2*kNumberOfSDD; ihyb++) modSet[ihyb]=0;
  Double_t nPtLay3 = 0;
  Double_t nPtLay4 = 0;

  Double_t param[4];    // parameters of poly fit
  Double_t minValP0=5.; // min value for param[0]
  Double_t maxValP0=8.; // max value for param[0]
  Double_t minValP1=0.; // min value for param[1]
  Double_t aveCoefLay3[4]={0.,0.,0.,0.};  // average param for good mod.
  Double_t aveCoefLay4[4]={0.,0.,0.,0.};  // average param for good mod.
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
  Int_t retfscf;
  

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
	  continue;
	}
	retfscf=fscanf(injFil,"%d",&polDeg);
	while (!feof(injFil)){
	  retfscf=fscanf(injFil,"%d %u ",&evNumb,&timeStamp);
	  if(evNumb==-99){
	    statusInj=timeStamp;
	    arr->SetInjectorStatus(statusInj);
	  }else{
	    if(feof(injFil)) break;
	    for(Int_t ic=0;ic<4;ic++){ 
	      retfscf=fscanf(injFil,"%f ",&auxP);
	      param[ic]=auxP;
	    }

	    if(polDeg>=0 && polDeg<=AliITSDriftSpeedSDD::GetMaxPolDeg() && 
	       param[0]>minValP0 && param[0]<maxValP0 && param[1]>minValP1){
	      if(polDeg==3){
		if(modID<kNumberOfSDDLay3){ 
		  for(Int_t ic=0;ic<4;ic++) aveCoefLay3[ic]+=param[ic];
		  nPtLay3++;
		}else{ 
		  for(Int_t ic=0;ic<4;ic++) aveCoefLay4[ic]+=param[ic];
		  nPtLay4++;
		}
	      }
	      AliITSDriftSpeedSDD *dsp=new AliITSDriftSpeedSDD(evNumb,timeStamp,polDeg,param);
	      arr->AddDriftSpeed(dsp);
	      modSet[2*modID+isid]=1;
	    }else{
	      Log(Form("Module %d side %d not accepted, degree=%d, params=%g %g %g %g",modID+240,isid,polDeg,param[0],param[1],param[2],param[3]));
	    }
	  }
	}
	fclose(injFil);
	Log(Form("Put calib obj for hybrid %d (DDL %d  Carlos %d)",2*modID+isid,iddl,imod));
	if(modSet[2*modID+isid]) vdrift.AddAt(arr,2*modID+isid);
      }
    }
  }

  // set drift speed for modules with bad injectors
  for(Int_t ic=0;ic<4;ic++){ 
    if(nPtLay3>0) aveCoefLay3[ic]/=nPtLay3; // mean parameters
    else aveCoefLay3[ic]=defCoef[ic]; // default parameters
    if(nPtLay4>0) aveCoefLay4[ic]/=nPtLay4; // mean parameters
    else aveCoefLay4[ic]=defCoef[ic]; // default parameters
  }
  AliITSDriftSpeedSDD *avdsp3=new AliITSDriftSpeedSDD(evNumb,timeStamp,3,aveCoefLay3);
  AliITSDriftSpeedSDD *avdsp4=new AliITSDriftSpeedSDD(evNumb,timeStamp,3,aveCoefLay4);

  // Check status of golden modules
  Int_t idGoldenMod=-1, idGoldenSide=-1;
  Int_t idGoldenModList[5]={319,319,321,243,243};
  Int_t idGoldenSideList[5]={0,1,0,0,1};
  AliITSDriftSpeedSDD* refSpeed=0x0;
  for(Int_t iGold=0; iGold<5; iGold++){
    Int_t indexG=2*(idGoldenModList[iGold]-240)+idGoldenSideList[iGold];
    if(modSet[indexG]){
      idGoldenMod=idGoldenModList[iGold];
      idGoldenSide=idGoldenSideList[iGold];
      AliITSDriftSpeedArraySDD* arrRef=(AliITSDriftSpeedArraySDD*)vdrift.At(indexG);
      refSpeed=arrRef->GetDriftSpeedObject(0);
      break;
    }
  }
  TList* correctionList=0x0;
  if(idGoldenMod>=240 && idGoldenSide>=0){
    // Get rescaling corrections from OCDB
    AliCDBEntry* entry = GetFromOCDB("Calib", "RescaleDriftSpeedSDD");
    if(!entry){
      Log("RescaleDriftSpeedSDD file not found in OCDB.");  
    }
    TList* fullList=(TList*)entry->GetObject();
    if(!fullList){ 
      Log("TList object not found in file.");
    }
    TString listName=Form("RefMod%d_Side%d",idGoldenMod,idGoldenSide);
    correctionList=(TList*)fullList->FindObject(listName.Data());
    if(!correctionList){
      Log(Form("TList for requested module %d side %d not found",idGoldenMod,idGoldenSide));
    }else{
      Log(Form("Use module %d side %d as reference module",idGoldenMod,idGoldenSide));
      if(refSpeed){
	Log(Form("Drift speed params for golden module = %g %g %g %g",refSpeed->GetDriftSpeedParameter(0),refSpeed->GetDriftSpeedParameter(1),refSpeed->GetDriftSpeedParameter(2),refSpeed->GetDriftSpeedParameter(3)));
      }else{
	AliError("No drift speed object for golden module");
      }
    }
  }
  
  for(Int_t ihyb=0; ihyb<2*kNumberOfSDDLay3; ihyb++){
    AliITSDriftSpeedArraySDD *arr=new AliITSDriftSpeedArraySDD();
    if(modSet[ihyb]==0){ 
      Int_t iBadMod=ihyb/2+240;
      Int_t iBadSide=ihyb%2;
      Bool_t goldenUsed=kFALSE;
      if(correctionList && refSpeed){
	Double_t *params=RescaleDriftSpeedModule(correctionList,iBadMod,iBadSide,refSpeed);
	if(params){
	  AliWarning(Form("No good injector events for mod. %d side %d --> use rescaled values from golden module",iBadMod,iBadSide));
	  AliITSDriftSpeedSDD* dspres=new AliITSDriftSpeedSDD(0,refSpeed->GetEventTimestamp(),3,params);
	  arr->AddDriftSpeed(dspres);
	  arr->SetInjectorStatus(1);
	  goldenUsed=kTRUE;
	}
      }
      if(!goldenUsed){
	AliWarning(Form("No good injector events for mod. %d side %d --> use average values for layer 3",iBadMod,iBadSide));
	arr->AddDriftSpeed(avdsp3);
	arr->SetInjectorStatus(0);
      }
      vdrift.AddAt(arr,ihyb);
    }
  }

  for(Int_t ihyb=2*kNumberOfSDDLay3; ihyb<2*kNumberOfSDD; ihyb++){
    if(modSet[ihyb]==0){ 
      AliWarning(Form("No good injector events for mod. %d side %d --> use average values for layer 4",ihyb/2,ihyb%2));
      AliITSDriftSpeedArraySDD *arr=new AliITSDriftSpeedArraySDD();
      arr->AddDriftSpeed(avdsp4);
      arr->SetInjectorStatus(0);
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
Double_t* AliITSPreprocessorSDD::RescaleDriftSpeedModule(TList* theList,
							 Int_t iBadMod, 
							 Int_t iBadSide,
							 AliITSDriftSpeedSDD* refSpeed)
  const
{
  // Rescale driftSpeed for a drift region starting from values of golden module

  if(!refSpeed) return 0x0;
  TString hisName=Form("hRatioMod%d_Side%d",iBadMod,iBadSide);
  TH1D* h=(TH1D*)theList->FindObject(hisName.Data());
  if(!h) return 0x0;

  TF1* fpoly=new TF1("fpoly","pol3",0.,255.);
  for(Int_t iAnode=0; iAnode<256; iAnode++){
    Double_t vref=refSpeed->GetDriftSpeedAtAnode((Double_t)iAnode);
    Double_t vcorr=h->GetBinContent(iAnode+1)*vref;
    h->SetBinContent(iAnode+1,vcorr);
  }
  h->Fit(fpoly,"RNQ");
  Double_t *params=new Double_t[4];
  for(Int_t iPar=0; iPar<4; iPar++) params[iPar]=fpoly->GetParameter(iPar);
  delete fpoly;
  return params;
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

