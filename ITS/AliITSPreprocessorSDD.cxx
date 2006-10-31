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
#include "TH1F.h"
#include "AliLog.h"
#include "TFile.h"

const Int_t AliITSPreprocessorSDD::fgkNumberOfSDD = 260;
const Int_t AliITSPreprocessorSDD::fgkNumberOfChannels = 512;
const TString AliITSPreprocessorSDD::fgkNameHistoPedestals = "hpedestal";
const TString AliITSPreprocessorSDD::fgkNameHistoNoise = "hnoise";
ClassImp(AliITSPreprocessorSDD)


UInt_t AliITSPreprocessorSDD::Process(TMap*/* dcsAliasMap*/){

  //preprocessing. 

  UInt_t result = 0;
  const char* filename = GetFile(kDAQ,"PEDESTALS","GDC");
  const char* filenamen= GetFile(kDAQ,"NOISE","GDC");
  const char* filenamed= GetFile(kDAQ,"DEADCHANNELS","GDC");
  TFile* f1 = TFile::Open(filename,"r");
  TFile* f2 = TFile::Open(filenamen,"r");
  Char_t namehisto[20];
  Char_t namehisto2[20];

  FILE* filed = fopen(filenamed,"read");
  Int_t numOfBadChannels[fgkNumberOfSDD];
  Int_t** badCh = new Int_t*[fgkNumberOfSDD];
  
  Char_t row[50];
  Int_t nSDD=0;
  Char_t str[20];
  char dims[1];
  Int_t dim;
  sprintf(str,"MODULE=%d",0);
  while(!feof(filed)){
    fscanf(filed,"%s\n",row);
    if(strcmp(row,str)==0){
      fscanf(filed,"%s %d\n",dims,&dim);
      badCh[nSDD] = new Int_t[dim];
      numOfBadChannels[nSDD]=dim;
      for(Int_t ibad=0;ibad<dim;ibad++){
	fscanf(filed,"%d\n",&badCh[nSDD][ibad]);
      }      
    }
    nSDD++;
    sprintf(str,"MODULE=%d",nSDD);
  }
  


  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSCalibration");
  md1->SetResponsible("Elisabetta Crescio");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head September 2005"); //root version
  md1->SetComment("This is a test");

  TObjArray respSDD(260);
  respSDD.SetOwner(kFALSE);
  
  for(Int_t imod=0;imod<fgkNumberOfSDD;imod++){
    AliITSCalibrationSDD *cal = new AliITSCalibrationSDD("simulated");
    cal->SetDeadChannels(numOfBadChannels[imod]);
    for(Int_t ich=0;ich<numOfBadChannels[imod];ich++){
      cal->SetBadChannel(ich,badCh[imod][ich]);
    }
    sprintf(namehisto,"%s_%d",fgkNameHistoPedestals.Data(),imod);
    sprintf(namehisto2,"%s_%d",fgkNameHistoNoise.Data(),imod);
    TH1F* hbas = (TH1F*)f1->Get(namehisto);
    TH1F* hnoi = (TH1F*)f2->Get(namehisto2);
    for(Int_t ien=0;ien<fgkNumberOfChannels;ien++){
      cal->SetBaseline(ien,hbas->GetBinContent(ien+1));
      cal->SetNoiseAfterElectronics(ien,hnoi->GetBinContent(ien+1));
    }
    respSDD.Add(cal);
  }

  result = Store("Calib","Data",&respSDD,md1);

  for(Int_t i=0;i<fgkNumberOfSDD;i++){
    delete badCh[i];
  }
  delete [] badCh;
  f1->Close();
  f2->Close();
  fclose(filed);
  return result;

}
