/*
CPV GAIN DA for processing physics runs and producing amplitude histograms in every channel for further calibration coefs calulation.

Contact: Sergey Evdokimov <sevdokim@cern.ch>
Link: https://twiki.cern.ch/twiki/bin/view/ALICE/CPVda
Reference run: 214340 (/afs/cern.ch/user/s/sevdokim/public/CPV_run214340_standalone.raw)
Run Type:  PHYSICS
DA Type: MON
Number of events needed: 1M events
Input files: thr?_??.dat CpvBadMap.root PHOSCPVGAINda.cfg
Output files: CpvCalibrSupply.root
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
#include "AliPHOSCpvGainCalibDA.h"
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
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH2S.h"
#include "TH2F.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom.h"


int main( int argc, char **argv )
{
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  Int_t status,statusBadCh=0,statusCalibrSupply=0,print;
  Int_t sigcut=3;
  Int_t minAmpl = 10;//minimal amplitude for consideration, to be read from DAQ DB
  Int_t minOccupancy = 1000;//min occupancy for publishing in OCDB, to be read from DAQ DB
  Int_t minClustSize=4;//min cluster size to consider
  Bool_t turbo = kTRUE;

  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // log start of process
  printf("Cpv gain calibration DA program started\n");

  /* report progress */
  daqDA_progressReport(0);

  /* retrieve configuration file from DAQ DB */
  status=daqDA_DB_getFile("PHOSCPVGAINda.cfg", "PHOSCPVGAINda.cfg");
  if(status==0) {
    char buf[500]; 
    FILE * fConf = fopen("PHOSCPVGAINda.cfg","r");
    while(fgets(buf, 500, fConf)){
      if(buf[0]=='#') continue;//comment indicator
      if(strstr(buf,"minOccupancy")) {
	sscanf(buf,"%*s %d",&minOccupancy);
	cout<<"I read minOccupancy="<<minOccupancy<<" from config file"<<endl;
      }
      if(strstr(buf,"minAmpl")) {
	sscanf(buf,"%*s %d",&minAmpl);
	cout<<"I read minAmpl="<<minAmpl<<" from config file"<<endl;
      }
      if(strstr(buf,"minClustSize")) {
	sscanf(buf,"%*s %d",&minClustSize);
	cout<<"I read minClustSize="<<minClustSize<<" from config file"<<endl;
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
	//return -1;
      }
    }
  }

  /* retrieve Bad Channel Map from DAQ DB */
  statusBadCh=daqDA_DB_getFile("CpvBadMap.root", "CpvBadMap.root");
  if(statusBadCh!=0) printf("Cannot retrieve file CpvBadMap.root from DAQ DB! Bad channels map will not be used!");

  /* retrieve previously collected histograms from DAQ DB */
  statusCalibrSupply=daqDA_DB_getFile("CpvCalibrSupply.root", "CpvCalibrSupply.root");
  if(statusCalibrSupply!=0) printf("Cannot retrieve file CpvCalibrSupply.root from DAQ DB! No previously collected histograms found!");
  
  //digiProducer
  AliPHOSCpvRawDigiProducer* digiProducer = new AliPHOSCpvRawDigiProducer();
  digiProducer->SetTurbo(turbo);
  digiProducer->LoadPedFiles();
  digiProducer->SetCpvMinAmp(minAmpl);

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
  
  //DA object
  AliPHOSCpvGainCalibDA *fDA = new AliPHOSCpvGainCalibDA();
  TFile *fCalibrSupplyRoot=0x0; 
  if(!statusCalibrSupply) fCalibrSupplyRoot = TFile::Open("CpvCalibrSupply.root");
  fDA->InitCalibration(fCalibrSupplyRoot);
  if(!statusBadCh)fDA->SetDeadChannelMapFromFile("CpvBadMap.root");
  fDA->SetMinClustSize(minClustSize);

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
      fDA->FillAmplA0Histos(digits);
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

  /* save collected histos, send files to DBs */
  fDA->WriteA0HistosToFile("CpvCalibrSupplyNew.root");

  //calculate occupancy
  Double_t Occupancy = 0;
  TFile* fSave = TFile::Open("CpvCalibrSupplyNew.root");
  for(Int_t iDDL = 0;iDDL<2*AliPHOSCpvParam::kNDDL; iDDL+=2){
    if(iDDL!=4)continue;
    if(fSave->Get(Form("hEntriesMap%d",iDDL))){
      TH2* hEntries = (TH2*)(fSave->Get(Form("hEntriesMap%d",iDDL)));
      Occupancy = hEntries->GetEntries()/7680.;
    }
  }
  fSave->Close();
  cout<<"Occupancy = "<<Occupancy<<"; minOccupancy = "<<minOccupancy<<endl;
  if(Occupancy>minOccupancy){//if we have enough statistics to calculate calibration
    status = daqDA_FES_storeFile("CpvCalibrSupplyNew.root","CPVAMPLITUDES");
    if(status) printf("Failed to store CpvCalibrSupplyNew.root in DAQ FXS!\n");
    TFile * fDummy = TFile::Open("dummy.root","RECREATE");
    status = daqDA_DB_storeFile("dummy.root","CpvCalibrSupply.root");
    fDummy->Close();
    if(status) printf("Failed to store dummy.root as CpvCalibrSupply.root in DAQ DB!\n");
  }
  else{//store CpvCalibrSupply.root in DAQ DB for future
    status = daqDA_DB_storeFile("CpvCalibrSupplyNew.root","CpvCalibrSupply.root");
    if(status) printf("Failed to CpvCalibrSupplyNew.root in DAQ DB!\n");

  }
  //send pictures to amore
  setenv("AMORE_DA_NAME","CPV-DAs",1);
  TList* histos = fDA->GetQAHistos();
  amore::da::AmoreDA* myAmore = new amore::da::AmoreDA(amore::da::AmoreDA::kSender);
  Int_t iHist = 0;
  while(histos->At(iHist)){
    myAmore->Send(histos->At(iHist)->GetName(),histos->At(iHist));
    iHist++;
  }

  /* report progress */
  daqDA_progressReport(100);


  return 0;
}
