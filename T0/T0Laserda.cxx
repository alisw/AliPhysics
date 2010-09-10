/*
T0 DA for online calibration
 
Contact: Alla.Maevskaya@cern.ch
Run Type: AMPLITUDE_CALIBRATION
DA Type: MON
Number of events needed: 5000 
Input Files: inLaser.dat, external parameters
Output Files: daLaser.root, to be exported to the DAQ FXS
Trigger types used: CALIBRATION_EVENT

*/

#define FILE_OUT "daLaser.root"
#define FILE_IN "inLaser.dat"
#include <daqDA.h>
#include <event.h>
#include <monitor.h>
 
#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include <AliRawReaderDate.h>
#include <AliRawReader.h>
#include <AliT0RawReader.h>

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include "TFile.h"
#include "TKey.h"
#include "TH2S.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
//AMORE
//
#ifdef ALI_AMORE
#include <AmoreDA.h>
#endif


float lcfd, hcfd;
/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  int status;
  
  /* magic line */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 
  
  Int_t nEntries = daqDA_ECS_getTotalIteration(); // usually = 11 = Nb of calibration runs
 cout<<" nEntries "<<nEntries<<endl;
 Int_t nInit=1;  // = 0 all DAC values ; = 1 DAC=0 excluded (default=1)
  // Reading current iteration
  Int_t nIndex = daqDA_ECS_getCurrentIteration();
  if(nIndex<0 || nIndex>nEntries) {printf("\n Failed: nIndex = %d\n",nIndex); return -1 ;}
  

// decode the input line
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
 

  Float_t mipsin[20];
   if(daqDA_DB_getFile(FILE_IN, FILE_IN)){
     printf("Couldn't get input file >>inLaser.dat<< from DAQ_DB !!!\n");
     return -1;
  }

  ifstream filein(FILE_IN,ios::in); 
  Int_t k=0;
  while (k<nEntries ) { 
    filein >>mipsin[k] ;
    cout<<" "<<mipsin[k]<<endl;
    k++; }
  filein>>lcfd;
  filein>>hcfd;
  filein.close();



  /*
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  */

  // Char_t inputFile[256]="";


  /* define data source : this is argument 1 */  
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
  

  /* log start of process */
  printf("T0 monitoring program started\n");  
  
  // Allocation of histograms - start
  
  TH1F *hQTC[24];  TH1F *hLED[24]; TH1F *hCFD[24];
  printf(" CFD low limit %f high %f",lcfd, hcfd);
  
  for(Int_t ic=0; ic<24; ic++) 
    {
      hQTC[ic] = new TH1F(Form("hQTC%d_%d",ic+1,nIndex),"QTC",1000, 500, 10000);
      hLED[ic] = new TH1F(Form("hLED%d_%d",ic+1,nIndex),"LED",125,200,700);
      hCFD[ic] = new TH1F(Form("hCFD%d_%d",ic+1,nIndex),"CFD", Int_t ((hcfd-lcfd)/2), lcfd, hcfd);
    }
 
  // Allocation of histograms - end
  
  // Reading current iteration
   
  Int_t iev=0;
 
  /* main loop (infinite) */
  for(;;) {
    struct eventHeaderStruct *event;
    eventTypeType eventT;
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic((void **)&event);
    if (status==(int)MON_ERR_EOF) {
      printf ("End of File detected\n");
      break; /* end of monitoring file has been reached */
    }
    
    if (status!=0) {
      printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
      break;
    }
    
    /* retry if got no event */
    if (event==NULL) {
      continue;
    }
    
    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    switch (event->eventType){
      
    case START_OF_RUN:
      break;
      
    case END_OF_RUN:
      break;
      
    case CALIBRATION_EVENT:
      //      case PHYSICS_EVENT:
      iev++; 
      
      if(iev==1) printf("First event  %i\n",iev);
      
      
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      
      // Enable the following two lines in case of real-data
      reader->RequireHeader(kTRUE);
      AliT0RawReader *start = new AliT0RawReader(reader, kTRUE);

      // Read raw data
      Int_t allData[106][5];
      for(Int_t i0=0;i0<106;i0++)
      	for(Int_t j0=0;j0<5;j0++)
		allData[i0][j0] = 0;
      
      if(start->Next()){
        for (Int_t i=0; i<106; i++) {
	  for(Int_t iHit=0;iHit<5;iHit++){
	    allData[i][iHit]= start->GetData(i,iHit);
	  }
        }
      }
      else 
	printf("No data for T0 found!!!\n");
 
      
      // Fill the histograms
      
      for (Int_t ik = 0; ik<24; ik+=2)
	{
	  Int_t cc = ik/2;
	  if(allData[ik+25][0]>0 && allData[ik+26][0]>0)
	    hQTC[cc]->Fill((allData[ik+25][0]-allData[ik+26][0]));
	  if(allData[cc+13][0]>0 && allData[cc+1][0]>0)
	    hLED[cc]->Fill(allData[cc+13][0]-allData[cc+1][0]);
	  if(allData[cc+1][0] > 0) hCFD[cc]->Fill(allData[cc+1][0]);
	  //	  printf("%i %i CFD %i LED_CFD %i\n",
	  //	  ik,cc,allData[cc+1][0],allData[cc+1][0]);

	}	   
      
      
      for (Int_t ik = 24; ik<48; ik+=2)
	{
	  Int_t cc = ik/2;
	   if(allData[ik+57][0]>0 && allData[ik+58][0]>0 && allData[cc+45][0]>0)
	     hQTC[cc]->Fill(allData[ik+57][0]-allData[ik+58][0]);
	   
	   if(allData[cc+57][0]>0 && allData[cc+45][0]>0)
	     hLED[cc]->Fill(allData[cc+57][0]-allData[cc+45][0]);
	   if(allData[cc+45][0] > 0) hCFD[cc]->Fill(allData[cc+45][0]);
	   //  printf("%i %i CFD %i LED_CFD %i\n",
	   //  ik,cc,allData[cc+45][0],allData[cc+45][0]);

	}
      
      delete start;
      start = 0x0;
      delete reader;
      reader= 0x0;
      // End of fill histograms

    }
    
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      printf("Number of events processed - %i\n ",iev);
      break;
    }
  }

  printf("After loop, before writing histos\n");
  TH1F* hAmp = new TH1F("hAmpLaser"," Laser amplitude ", 1000, 0.5, 20.5);  
  for(Int_t i=0; i<nIndex; i++)   hAmp->Fill(mipsin[i]); 
  // write a file with the histograms
  TFile *hist=0;
  if(nIndex == 1 )
    hist = new TFile(FILE_OUT,"RECREATE");
  else
    hist = new TFile(FILE_OUT,"UPDATE");
  hist->cd();
  
  for(Int_t j=0;j<24;j++)
    {
      if(hQTC[j]->GetEntries()>0) hQTC[j]->Write();
      if(hLED[j]->GetEntries()>0) hLED[j]->Write();
      if(hCFD[j]->GetEntries()>0) hCFD[j]->Write();
    }
  hAmp->Write();
  
  hist->Close();
  delete hist;
  
  status=0;
  
  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_OUT, "AMPLITUDE_CALIBRATION")) {
    status=-2;
  }
  
  return status;
}

