/////////////////////////////////////////
// Class for SDD digits preprocessing  //
//                                     //
//                                     //
////////////////////////////////////////

#include "AliITSPreprocessorSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliShuttleInterface.h"
#include "AliCDBMetaData.h"
#include "TObjArray.h"
#include "AliLog.h"

const Int_t AliITSPreprocessorSDD::fgkNumberOfSDD = 260;
const Int_t AliITSPreprocessorSDD::fgkNumberOfChannels = 512;
const TString AliITSPreprocessorSDD::fgkNameHistoPedestals = "hpedestal";
const TString AliITSPreprocessorSDD::fgkNameHistoNoise = "hnoise";
ClassImp(AliITSPreprocessorSDD)


UInt_t AliITSPreprocessorSDD::Process(TMap*/* dcsAliasMap*/){

  //preprocessing. 


  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSCalibration");
  md1->SetResponsible("Elisabetta Crescio, Francesco Prino");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 5 April 2007"); //root version
  md1->SetComment("This is a test");

  TObjArray respSDD(fgkNumberOfSDD);
  respSDD.SetOwner(kFALSE);
  
  Char_t filid[20];
  Float_t baseline,rawnoise,cmn,gain;
  Int_t isgoodan,i,im,is,isgoodmod;
  Int_t numOfBadChannels[fgkNumberOfSDD];
  for(Int_t imod=0;imod<fgkNumberOfSDD;imod++){
    AliITSCalibrationSDD *cal = new AliITSCalibrationSDD("simulated");
    numOfBadChannels[imod]=0;
    Int_t badch[fgkNumberOfChannels];
    for(Int_t isid=0;isid<=1;isid++){
      sprintf(filid,"DAQDAm%03ds%d",imod,isid);
      const char* filenamed= GetFile(kDAQ,filid,"GDC");
      FILE* filed = fopen(filenamed,"read");
      if (filed == 0) {
	AliWarning("File not found"); 
	continue;
      }
      
      fscanf(filed,"%d %d %d\n",&im,&is,&isgoodmod);
      if(!isgoodmod) cal->SetDead();
      for(Int_t ian=0;ian<(fgkNumberOfChannels/2);ian++){
	fscanf(filed,"%d %d %f %f %f %f\n",&i,&isgoodan,&baseline,&rawnoise,&cmn,&gain);
	Int_t ich=ian;
	if(isid==1) ich+=256;
	if(!isgoodan){ 
	  Int_t ibad=numOfBadChannels[imod];
	  numOfBadChannels[imod]++;
	  badch[ibad]=ich;
	}
	cal->SetBaseline(ich,baseline);
	cal->SetNoiseAfterElectronics(ich,rawnoise);	
	//	cal->SetGain(gain,isid,ian/4,ian%4);
      }
      cal->SetDeadChannels(numOfBadChannels[imod]);
      for(Int_t ibad=0;ibad<numOfBadChannels[imod];ibad++){
	cal->SetBadChannel(ibad,badch[ibad]);
      }
      fclose(filed);
    }
    respSDD.Add(cal);
  }

  Bool_t retcode = Store("Calib","Data",&respSDD,md1);


  if(retcode) return 0;
  else return 1;
}
