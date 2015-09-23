/*
CPV BCM DA for processing physics runs and producing bad channel map for further use at reconstruction time.

Contact: Sergey Evdokimov <sevdokim@cern.ch>
Link: https://twiki.cern.ch/twiki/bin/view/ALICE/CPVda
Reference run: 214340 (/afs/cern.ch/user/s/sevdokim/public/CPV_run214340_standalone.raw)
Run Type:  PHYSICS
DA Type: MON
Number of events needed: ~100000 events
Input files: thr?_??.dat CpvPeds.root PHOSCPVBCMda.cfg
Output files: CpvBadMap.root
Trigger types used: PHYSICS_EVENT
*/

//daqDA
#include "event.h"
#include "monitor.h"
#include "daqDA.h"
//AMORE monitoring framework
#include <AmoreDA.h>

//system
#include <Riostream.h>
#include <stdlib.h>
#include <fstream>
#include <string>

//AliRoot
#include "AliPHOSCpvRawDigiProducer.h"
#include "AliPHOSCpvParam.h"
#include "AliRawReaderDate.h"
#include "AliBitPacking.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGeometry.h"

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include "TSAXParser.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH2S.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom.h"

Double_t FillNoisyMap(TH2F* DigMap, TH2I* BadMap); //returns mean occupancy when all bad channels are excluded
void FillDeadMap(TH2F* PedMap, TH2F* DigMap, TH2I* BadMap,Double_t meanOccupancy);

int main( int argc, char **argv )
{
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  Int_t status,statusPeds=0,print,statusBadMap=0;
  Int_t sigcut=3;
  Int_t minAmpl = 10;//minimal amplitude for consideration, to be read from DAQ DB
  Int_t minOccupancy = 10;//min occupancy for publishing in OCDB, to be read from DAQ DB
  Int_t minClustSize=4;//min cluster size to consider
  Bool_t turbo = kTRUE;

  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // log start of process
  printf("Cpv bad channel map DA program started\n");

  /* report progress */
  daqDA_progressReport(0);

  /* retrieve configuration file from DAQ DB */
  status=daqDA_DB_getFile("PHOSCPVBCMda.cfg", "PHOSCPVBCMda.cfg");
  status = 0;
  if(!status) {
    char buf[500]; 
    FILE * fConf = fopen("PHOSCPVBCMda.cfg","r");
    while(fgets(buf, 500, fConf)){
      if(buf[0]=='#') continue;//comment indicator
      if(strstr(buf,"minOccupancy")){ 
	sscanf(buf,"%*s %d",&minOccupancy);
	cout<<"I read minOccupancy="<<minOccupancy<<" from config file"<<endl;
      }
      if(strstr(buf,"minAmpl")) {
	sscanf(buf,"%*s %d",&minAmpl);
	cout<<"I read minAmpl="<<minAmpl<<" from config file"<<endl;
      }
    }
    fclose(fConf);
  }
  

  /* retrieve pedestal tables from DAQ DB */
  for(int iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL+=2){
    if(iDDL!=4) continue; // only one module with DDL=4 by now
    for (int iCC = 0; iCC<AliPHOSCpvParam::kNRows; iCC++){
      status=daqDA_DB_getFile(Form("thr%d_%02d.dat", iDDL, iCC),Form("thr%d_%02d.dat", iDDL, iCC));
      if(status!=0) {
	printf("cannot retrieve file %s from DAQ DB. Exit.\n", Form("thr%d_%02d.dat", iDDL, iCC));
	//return 11;//error code 11 (cannot retrive thr.dat from DAQ DB)
      }
    }
  }

  /* retrieve pedestals in root format to find dead channels */
  statusPeds=daqDA_DB_getFile("CpvPeds.root", "CpvPeds.root");
  if(statusPeds) {
    printf("cannot retrieve CpvPeds.root from DAQ DB! No dead channels will be found.");
    //return 12; //error code 12 (cannot retrive CpvPeds.root from DAQ DB)
  }
  
  /* retrieve bad map from DAQ DB to see if we have some statistics saved form previous runs */
  statusBadMap=daqDA_DB_getFile("CpvBadMap.root", "CpvBadMap.root");
  
  //digiProducer
  AliPHOSCpvRawDigiProducer* digiProducer = new AliPHOSCpvRawDigiProducer();
  digiProducer->SetTurbo(turbo);
  digiProducer->LoadPedFiles();
  digiProducer->SetCpvMinAmp(minAmpl);

  TH2I *hBadChMap[2*AliPHOSCpvParam::kNDDL];
  for(int i=0;i<2*AliPHOSCpvParam::kNDDL;i++)
    hBadChMap[i]=0x0;

  /* retrieve permanent bad map from DAQ DB */
  status=daqDA_DB_getFile("CpvPermanentBadMap.root","CpvPermanentBadMap.root");
  if(status!=0) {
    printf("cannot retrieve file %s from DAQ DB. \n", "CpvPermanentBadMap.root");
  }
  else{
    TFile *fPBM = TFile::Open("CpvPermanentBadMap.root","r");
    for(int iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL+=2){
      if(iDDL!=4) continue; // only one module with DDL=4 by now
      TH2I* badMap=(TH2I*)fPBM->Get(Form("fBadMap%d",iDDL));
      if(badMap){
	digiProducer->SetPermanentBadMap(badMap,iDDL);
	hBadChMap[iDDL]=(TH2I*)badMap->Clone();
      }
    }
  }

  /* connecting to raw data */
  status=monitorSetDataSource( argv[1] );
  if (status!=0) {
    printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
    return -1;
  }

  /* declare monitoring program */
  status=monitorDeclareMp( __FILE__ );
  if (status!=0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(status));
    return -1;
  }

  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

    /* report progress */
  daqDA_progressReport(5);


  // init event counter
  Int_t iPhysEvnt=0;
  Int_t iTotEvnt =0;

  // Reader
  AliRawReader * reader;


  //digits
  TClonesArray *digits = new TClonesArray("AliPHOSDigit",1);
  digits->SetName("DIGITS");
  
  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //maps of digits and bad channels
  TH2F* hMapOfDig[2*AliPHOSCpvParam::kNDDL]; 
  
  for (int i = 0;i<2*AliPHOSCpvParam::kNDDL;i++){
    hMapOfDig[i]= 0x0;
  }

  //any previously gained statistics?
  TFile *fPreviousStatistics = 0x0;
  if(!statusBadMap) fPreviousStatistics = TFile::Open("CpvBadMap.root");

  for(Int_t iDDL=0;iDDL<2*AliPHOSCpvParam::kNDDL;iDDL++){
    if(!statusBadMap){//we have some statistics from previous runs
      if(fPreviousStatistics->Get(Form("hMapOfDig%d",iDDL)))
	//hMapOfDig[iDDL] = new TH2F(*(TH2F*)(fPreviousStatistics->Get(Form("hMapOfDig%d",iDDL))));
	hMapOfDig[iDDL] = (TH2F*)(fPreviousStatistics->Get(Form("hMapOfDig%d",iDDL))->Clone());
    }
    if(!hMapOfDig[iDDL])
      hMapOfDig[iDDL] = new TH2F(Form("hMapOfDig%d",iDDL),Form("Map of digits with subtructed pedestals, DDL = %d",iDDL),
				 AliPHOSCpvParam::kPadPcX,0,AliPHOSCpvParam::kPadPcX,
				 AliPHOSCpvParam::kPadPcY,0,AliPHOSCpvParam::kPadPcY);
    if(!hBadChMap[iDDL])
      hBadChMap[iDDL] = new TH2I(Form("hBadMap%d",iDDL),Form("Bad Channels Map, DDL= %d",iDDL),
				 AliPHOSCpvParam::kPadPcX,0,AliPHOSCpvParam::kPadPcX,
				 AliPHOSCpvParam::kPadPcY,0,AliPHOSCpvParam::kPadPcY);
      
  }
  //if(!statusBadMap) fPreviousStatistics->Close();

  /* report progress */
  daqDA_progressReport(10);

  /* main loop (infinite) */
  for(;;) { // infinite loop
    struct eventHeaderStruct *event;
    eventTypeType eventT;

    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}

    // get next event
    status=monitorGetEventDynamic((void **)&event);
    if (status==MON_ERR_EOF) { // end of monitoring file has been reached
      printf("End of monitoring file has been reached! \n");
      break;
    }

    if (status!=0) {
      printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
      break;
    }

    // retry if got no event
    if (event==NULL) continue;

    // use event - here, just write event id to result file
    eventT=event->eventType;
    if (eventT==PHYSICS_EVENT) { //we use PHYSICS_EVENT for pedestal not CALIBRATION_EVENT???
      iTotEvnt++;
      reader = new AliRawReaderDate((void*)event);
      digiProducer->LoadNewEvent(reader);
      digiProducer->MakeDigits(digits);
      if(digits->GetEntriesFast()>0) iPhysEvnt++;
      for(Int_t i=0;i<digits->GetEntriesFast();i++)
	{
	  Int_t relId[4];
	  AliPHOSDigit *digit = static_cast<AliPHOSDigit*>( digits->At(i) ) ; 
	  geom->AbsToRelNumbering(digit->GetId(),relId);
	  //cout<<"relId : "<<relId[0]<<" "<<relId[1]<<" "<<relId[2]<<" "<<relId[3]<<endl;
	  if(
	     (relId[0]>=0&&relId[0]<=AliPHOSCpvParam::kNDDL)
	     &&(relId[1]==-1)//cpv
	     &&(relId[2]-1>=0&&relId[2]-1<AliPHOSCpvParam::kPadPcX)
	     &&(relId[3]-1>=0&&relId[3]-1<AliPHOSCpvParam::kPadPcY)
	     )
	    hMapOfDig[AliPHOSCpvParam::Mod2DDL(relId[0])]->Fill(relId[2]-1,relId[3]-1);
	}
      digits->Clear("C");
      delete reader;
    } // if PHYSICS_EVENT

    free(event);

    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }

  Printf(" Received %d events, %d good events",iTotEvnt,iPhysEvnt);
  /* report progress */
  daqDA_progressReport(90);

  TH2F* hPedMap[2*AliPHOSCpvParam::kNDDL];
  for(int iDDL = 0; iDDL< 2*AliPHOSCpvParam::kNDDL; iDDL++) hPedMap[iDDL] = 0x0;
  //prepare ped maps for dead channels search
  TFile* fPeds = TFile::Open("CpvPeds.root");
  if(fPeds){
    for(int iDDL = 0; iDDL< 2*AliPHOSCpvParam::kNDDL; iDDL+=2)
      if(fPeds->Get(Form("fPedMeanMap%d",iDDL)))
	//hPedMap[iDDL] = new TH2F(*(TH2F*)(fPeds->Get(Form("fPedMeanMap%d",iDDL))));
	hPedMap[iDDL] = (TH2F*)(fPeds->Get(Form("fPedMeanMap%d",iDDL))->Clone());
    //fPeds->Close();
  }
  
  //find noisy channels (i.e. channelOccupancy > 10*meanOccupancy)
  TFile *fSave = TFile::Open("CpvBadMap.root","RECREATE");

  setenv("AMORE_DA_NAME","CPV-DAs",1);
  amore::da::AmoreDA* myAmore = new amore::da::AmoreDA(amore::da::AmoreDA::kSender);
  Bool_t isStatisticsEnough=kFALSE;
  for(int iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL+=2){
    if(hMapOfDig[iDDL]->GetEntries()>0) {
      Double_t Occupancy = FillNoisyMap(hMapOfDig[iDDL],hBadChMap[iDDL]);
      if(Occupancy>minOccupancy) FillDeadMap(hPedMap[iDDL],hMapOfDig[iDDL],hBadChMap[iDDL],Occupancy);
      fSave->WriteObject(hMapOfDig[iDDL],Form("hMapOfDig%d",iDDL));
      //send digit maps to amore
      myAmore->Send(Form("hMapOfDig%d",iDDL),hMapOfDig[iDDL]);
      fSave->WriteObject(hBadChMap[iDDL],Form("hBadChMap%d",iDDL));
      cout<< "meanOccupancy in DDL"<<iDDL<<" = "<<Occupancy<<"; minOccupancy = "<<minOccupancy<<endl;
      if(Occupancy>minOccupancy)isStatisticsEnough=kTRUE;
    }
  }
  fSave->Close();
  
  
  if(isStatisticsEnough){//send file to FES if only statistics is enough
    status = daqDA_FES_storeFile("CpvBadMap.root","CPVBADMAP");
    if(status) printf("Failed to store CpvBadMap.root in DAQ FXS!\n");
    //store dummy file in DAQ DB
    TFile *fDummy = TFile::Open("dummy.root","RECREATE");
    fDummy->Close();
    status = daqDA_DB_storeFile("dummy.root","CpvBadMap.root");
    if(status) printf("Failed to store dummy.root as CpvBadMap.root in DAQ DB!\n");
    //send bad map to amore as well
    for(int iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL+=2)
      if(hMapOfDig[iDDL]->GetEntries()>0)
	myAmore->Send(Form("hBadChMap%d",iDDL),hBadChMap[iDDL]);
  }
  else{//store file with current statistics in DAQ DB for further use.
    status = daqDA_DB_storeFile("CpvBadMap.root","CpvBadMap.root");
    if(status) printf("Failed to store CpvBadMap.root in DAQ DB!\n");

  }

  /* report progress */
  daqDA_progressReport(100);


  return status;
}
//==============================================================================
Double_t FillNoisyMap(TH2F* hDigMap, TH2I* hBadChMap){
  if(!hDigMap)return 0;
  if(!hBadChMap)return 0;
  Double_t meanOccupancy = hDigMap->GetEntries()/7680.;
  Double_t nDigits=0,nChannels=0;
  int x,y,iterationNumber=1;
  int nBadChannelsTotal=0,nBadChannelsCurrentIteration=0;
  //1st iteration
  //cout<<"Iteration Number = "<<iterationNumber<<endl;
  for(int ix = 1;ix<=AliPHOSCpvParam::kPadPcX;ix++)
    for(int iy = 1;iy<=AliPHOSCpvParam::kPadPcY;iy++)
      if(hDigMap->GetBinContent(ix,iy)>meanOccupancy*10.+1.){
	nBadChannelsCurrentIteration++;
	hBadChMap->Fill(ix-1,iy-1);
	printf("Noisy channel found! DDL=4 x=%d y=%d\n",ix-1,iy-1);
      }
  //next iterations
  while(nBadChannelsCurrentIteration!=0){
    iterationNumber++;
    //cout<<"Iteration Number = "<<iterationNumber<<endl;
    nDigits=0;nChannels=0;
    nBadChannelsTotal+=nBadChannelsCurrentIteration;
    nBadChannelsCurrentIteration=0;
    //1st -- calculate new mean occupancy excluding already badly marked channels  
    for(int ix = 1;ix<=AliPHOSCpvParam::kPadPcX;ix++)
      for(int iy = 1;iy<=AliPHOSCpvParam::kPadPcY;iy++)
	if(hBadChMap->GetBinContent(ix,iy)==0){
	  nDigits+=hDigMap->GetBinContent(ix,iy);
	  nChannels+=1;
	}
    //cout<<"nDigits = "<<nDigits<<" ; nChannels = "<<nChannels<<endl;
    meanOccupancy=nDigits/nChannels;
    //cout<<"meanOccupancy = "<<meanOccupancy<<endl;
    //2nd -- spot new bad channels
    for(int ix = 1;ix<=AliPHOSCpvParam::kPadPcX;ix++)
      for(int iy = 1;iy<=AliPHOSCpvParam::kPadPcY;iy++)
	if(hBadChMap->GetBinContent(ix,iy)==0&&hDigMap->GetBinContent(ix,iy)>meanOccupancy*10.+1.){
	  nBadChannelsCurrentIteration++;
	  hBadChMap->Fill(ix-1,iy-1);
	  printf("Noisy channel found! DDL=4 x=%d y=%d\n",ix-1,iy-1);
	}
    
  }
  cout<<"Total number of noisy channels (DDL = 4): "<<nBadChannelsTotal<<endl;
  return meanOccupancy;
}
//==============================================================================
void FillDeadMap(TH2F* hPedMap, TH2F* hDigMap, TH2I* hBadChMap, Double_t meanOccupancy){
  if(!hPedMap) return;
  if(!hBadChMap) return;
  if(!hDigMap)return;
  Int_t nDeadTotal = 0;
  for(int ix = 1;ix<=AliPHOSCpvParam::kPadPcX;ix++)
    for(int iy = 1;iy<=AliPHOSCpvParam::kPadPcY;iy++){
      if((hPedMap->GetBinContent(ix,iy) < 1.) && (hBadChMap->GetBinContent(ix,iy)==0)) {
	  nDeadTotal++;
	  hBadChMap->Fill(ix-1,iy-1);
	  printf("Dead channel found! DDL=4 x=%d y=%d\n",ix-1,iy-1);
	}
      if((hDigMap->GetBinContent(ix,iy) < meanOccupancy/10.0-1.0) && (hBadChMap->GetBinContent(ix,iy)==0)) {
	  nDeadTotal++;
	  hBadChMap->Fill(ix-1,iy-1);
	  printf("Dead channel found! DDL=4 x=%d y=%d\n",ix-1,iy-1);
	}
    }
  cout<<"Total number of dead channels (DDL = 4): "<<nDeadTotal<<endl;
}
