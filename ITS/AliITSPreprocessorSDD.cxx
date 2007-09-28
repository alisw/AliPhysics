/////////////////////////////////////////
// Class for SDD digits preprocessing  //
//                                     //
//                                     //
////////////////////////////////////////

/* $Id$ */

#include "AliITSPreprocessorSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliShuttleInterface.h"
#include "AliCDBMetaData.h"
#include "TObjArray.h"
#include "AliLog.h"
#include <TObjString.h>
#include <TSystem.h>
#include <TList.h>

const Int_t AliITSPreprocessorSDD::fgkNumberOfSDD = 260;
const Int_t AliITSPreprocessorSDD::fgkNumberOfChannels = 512;
const Int_t AliITSPreprocessorSDD::fgkNumberOfChannelsPerChip = 64;
const TString AliITSPreprocessorSDD::fgkNameHistoPedestals = "hpedestal";
const TString AliITSPreprocessorSDD::fgkNameHistoNoise = "hnoise";
ClassImp(AliITSPreprocessorSDD)


UInt_t AliITSPreprocessorSDD::Process(TMap*/* dcsAliasMap*/){

  //preprocessing. 

  TObjArray respSDD(fgkNumberOfSDD);
  respSDD.SetOwner(kFALSE);
  Float_t baseline,rawnoise,cmn,gain;
  Int_t isgoodan,i,im,is,isgoodmod;
  Int_t numOfBadChannels[fgkNumberOfSDD];
  //TString pwd = gSystem->pwd();
//  const Char_t* tempDir=AliShuttleInterface::GetShuttleTempDir();
  
  TList* sourceList = GetFileSources(kDAQ, "SDD_Calib");
  if (!sourceList){ 
    Log("Error: no sources found for SDD_Calib");
    return 2;
  }

  Int_t ind = 0;
  Char_t command[100];
  while (sourceList->At(ind)!=NULL) {
    TObjString* tarId = (TObjString*) sourceList->At(ind);
    TString tarName = GetFile(kDAQ, "SDD_Calib", tarId->GetString().Data());
//    gSystem->cd(tempDir);
    sprintf(command,"tar -xf %s",tarName.Data());
    gSystem->Exec(command);
    ind++;
  }
   
  sourceList = GetFileSources(kDAQ, "SDD_Injec");
  if (!sourceList){ 
    Log("Error: no sources found for SDD_Injec");
    return 2;
  }

  ind = 0;
  while (sourceList->At(ind)!=NULL) {
   TObjString* tarId = (TObjString*) sourceList->At(ind);
     TString tarName = GetFile(kDAQ, "SDD_Injec", tarId->GetString().Data());
//    gSystem->cd(tempDir);
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
      Char_t basFileName[100];
//      sprintf(basFileName,"%s/SDDbase_mod%03d_sid%d.data",tempDir,imod,isid);
      sprintf(basFileName,"./SDDbase_mod%03d_sid%d.data",imod,isid);
      FILE* basFil = fopen(basFileName,"read");
      if (basFil == 0) {
	Log(Form("File %s not found.",basFileName));
	return 2;
      }      
      fscanf(basFil,"%d %d %d\n",&im,&is,&isgoodmod);
      if(!isgoodmod) cal->SetDead();
      for(Int_t ian=0;ian<(fgkNumberOfChannels/2);ian++){
	fscanf(basFil,"%d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&rawnoise,&cmn,&gain);
	Int_t ich=ian;
	if(isid==1) ich+=256;
	if(!isgoodan){ 
	  Int_t ibad=numOfBadChannels[imod];
	  numOfBadChannels[imod]++;
	  badch[ibad]=ich;
	}
	cal->SetBaseline(ich,baseline);
	cal->SetNoiseAfterElectronics(ich,rawnoise);
	Int_t iChip=ian/fgkNumberOfChannelsPerChip;
	Int_t iChInChip=ian%fgkNumberOfChannelsPerChip;
	cal->SetGain(gain,isid,iChip,iChInChip);
      }
      cal->SetDeadChannels(numOfBadChannels[imod]);
      for(Int_t ibad=0;ibad<numOfBadChannels[imod];ibad++){
	cal->SetBadChannel(ibad,badch[ibad]);
      }
      fclose(basFil);

      Char_t injFileName[100];
      Int_t evNumb; 
      UInt_t timeStamp;
      Float_t param[4];
      sprintf(injFileName,"./SDDinj_mod%03d_sid%d.data",imod,isid);
      FILE* injFil = fopen(injFileName,"read");
      if (injFil == 0) {
	Log(Form("File %s not found.",basFileName));
	return 2;
      }      
      while (!feof(injFil)){
	fscanf(injFil,"%d %d %",&evNumb,&timeStamp);
	if(feof(injFil)) break;
	for(Int_t ic=0;ic<4;ic++) fscanf(injFil,"%f",&param[ic]);
	cal->SetDriftSpeedParam(isid,param);
      }
    }
    respSDD.Add(cal);
  }

  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSCalibration");
  md1->SetResponsible("Elisabetta Crescio, Francesco Prino");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 5 April 2007"); //root version
  md1->SetComment("This is a test");

  Bool_t retcode = Store("Calib","CalibSDD",&respSDD,md1, 0, kTRUE);

  if(retcode) return 0;
  else return 1;
}
