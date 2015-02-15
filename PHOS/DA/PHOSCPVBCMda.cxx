/*
CPV BCM DA for processing physics runs and producing bad channel map for further use at reconstruction time.

Contact: Sergey Evdokimov <sevdokim@cern.ch>
Link: https://twiki.cern.ch/twiki/bin/view/ALICE/CPVda
Reference run: 
Run Type:  PHYSICS
DA Type: MON
Number of events needed: ? events
Input files: thr?_??.dat
Output files: CpvBadMap.root
Trigger types used: PHYSICS_EVENT
*/

//daqDA
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

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
#include "TObject.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom.h"

void FillBadMap(TH2* DigMap, TH2* BadMap);

int main( int argc, char **argv )
{
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  Int_t status,print;
  Int_t sigcut=3;
  Bool_t turbo = kTRUE;

  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  // log start of process
  printf("Cpv bad channel map DA program started\n");

  /* report progress */
  daqDA_progressReport(0);


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

  //digiProducer
  AliPHOSCpvRawDigiProducer* digiProducer = new AliPHOSCpvRawDigiProducer();
  digiProducer->SetTurbo(turbo);
  digiProducer->LoadPedFiles();
  digiProducer->SetCpvMinAmp(0);

  //digits
  TClonesArray *digits = new TClonesArray("AliPHOSDigit",1);
  digits->SetName("DIGITS");
  
  const AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance() ;

  //maps of digits and bad channels
  TH2F* hMapOfDig[2*AliPHOSCpvParam::kNDDL]; 
  TH2I *hBadChMap[2*AliPHOSCpvParam::kNDDL]; 
  for(Int_t iDDL=0;iDDL<2*AliPHOSCpvParam::kNDDL;iDDL++){
    hMapOfDig[iDDL] = new TH2F(Form("hMapOfDig%d",iDDL),Form("Map of digits with substructed pedestals, DDL = %d",iDDL),
			  AliPHOSCpvParam::kPadPcX,0,AliPHOSCpvParam::kPadPcX,
			  AliPHOSCpvParam::kPadPcY,0,AliPHOSCpvParam::kPadPcY);
    hBadChMap[iDDL] = new TH2I(Form("hBadMap%d",iDDL),Form("Bad Channels Map, DDL= %d",iDDL),
			  AliPHOSCpvParam::kPadPcX,0,AliPHOSCpvParam::kPadPcX,
			  AliPHOSCpvParam::kPadPcY,0,AliPHOSCpvParam::kPadPcY);
  }

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

  TFile *fSave = TFile::Open("CpvBadMap.root","RECREATE");

  for(int iDDL = 0; iDDL<2*AliPHOSCpvParam::kNDDL; iDDL+=2){
    if(hMapOfDig[iDDL]->GetEntries()>0) {
      FillBadMap(hMapOfDig[iDDL],hBadChMap[iDDL]);
      fSave->WriteObject(hMapOfDig[iDDL],Form("hMapOfDig%d",iDDL));
      fSave->WriteObject(hBadChMap[iDDL],Form("hBadChMap%d",iDDL));
    }
  }
  fSave->Close();
  //status = daqDA_DB_storeFile("CpvBadMap.root","CpvBadMap.root");
  //if(status) printf("Failed to store CpvBadMap.root in DAQ DB!\n");
  status = daqDA_FES_storeFile("CpvBadMap.root","CpvBadMap.root");
  if(status) printf("Failed to store CpvBadMap.root in DAQ FXS!\n");


  /* report progress */
  daqDA_progressReport(100);


  return status;
}
//==============================================================================
void FillBadMap(TH2* hDigMap, TH2* hBadChMap){
  Double_t meanOccupancy = hDigMap->GetEntries()/7680.;
  Double_t nDigits=0,nChannels=0;
  int x,y,iterationNumber=1;
  int nBadChannelsTotal=0,nBadChannelsCurrentIteration=0;
  //1st iteration
  //cout<<"Iteration Number = "<<iterationNumber<<endl;
  for(int ix = 1;ix<=128;ix++)
    for(int iy = 1;iy<=60;iy++)
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
    for(int ix = 1;ix<=128;ix++)
      for(int iy = 1;iy<=60;iy++)
	if(hBadChMap->GetBinContent(ix,iy)==0){
	  nDigits+=hDigMap->GetBinContent(ix,iy);
	  nChannels+=1;
	}
    //cout<<"nDigits = "<<nDigits<<" ; nChannels = "<<nChannels<<endl;
    meanOccupancy=nDigits/nChannels;
    //cout<<"meanOccupancy = "<<meanOccupancy<<endl;
    //2nd -- spot new bad channels
    for(int ix = 1;ix<=128;ix++)
      for(int iy = 1;iy<=60;iy++)
	if(hBadChMap->GetBinContent(ix,iy)==0&&hDigMap->GetBinContent(ix,iy)>meanOccupancy*10.+1.){
	  nBadChannelsCurrentIteration++;
	  hBadChMap->Fill(ix-1,iy-1);
	  printf("Noisy channel found! DDL=4 x=%d y=%d\n",ix-1,iy-1);
	}
    
  }
  cout<<"Total number of bad channels: "<<nBadChannelsTotal<<endl;
}

