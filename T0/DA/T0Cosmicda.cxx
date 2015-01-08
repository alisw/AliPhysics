/*
T0 DA for online calibration
 
Contact: Michal.Oledzki@cern.ch
Link: http://users.jyu.fi/~mioledzk/
Run Type: PHYSICS
DA Type: MON
Number of events needed: 400000 
Input Files: inCosmicLaser.dat, external parameters
Output Files: daCosmicLaser.root, to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT

*/

#define FILE_OUT "daCosmic.root"
//
//author Michal Oledzki
//DA for old cosmic runs
//
//
#define FILE_IN "inCosmic.dat"
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
//#include "TProfile.h"
int kCcbx;
float kCclx,kCcmx;
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
  
  if(daqDA_DB_getFile(FILE_IN, FILE_IN)){
     printf("Couldn't get input file >>inCosmic.dat<< from DAQ_DB !!!\n");
     return -1;
  }
  
 
  FILE *inp;
  char c;
  inp = fopen(FILE_IN, "r");  
  if(!inp){
	printf("Input file >>inCosmic.dat<< not found !!!\n");
	return -1;
  }

  while((c=getc(inp))!=EOF) {
    switch(c) {
      case 'a': {fscanf(inp, "%d", &kCcbx ); break;} //N of X bins hCFD1_CFD
      case 'b': {fscanf(inp, "%f", &kCclx ); break;} //Low x hCFD1_CFD
      case 'c': {fscanf(inp, "%f", &kCcmx ); break;} //High x hCFD1_CFD
    }
  }
  fclose(inp);

  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


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
  TH1F *hCFD1minCFD[24]; 

   for(Int_t ic=0; ic<24; ic++) {
      hCFD1minCFD[ic] = new TH1F(Form("CFD1-CFD%d",ic+1),"CFD-CFD",kCcbx,kCclx,kCcmx);
    }

  // Allocation of histograms - end

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
 
//    iev++; 

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

      if(iev==1){
           printf("First event - %i\n",iev);
      }

      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
          
      // Enable the following two lines in case of real-data
      	  reader->RequireHeader(kTRUE);
      AliT0RawReader *start = new AliT0RawReader(reader, kTRUE);

      // Read raw data
      Int_t allData[105][5];
      for(Int_t i0=0;i0<105;i0++)
      	for(Int_t j0=0;j0<5;j0++)
		allData[i0][j0] = 0;
    
       if(start->Next()){
       for (Int_t i=0; i<105; i++) {
	for(Int_t iHit=0;iHit<5;iHit++){
	  allData[i][iHit]= start->GetData(i,iHit);
        }
       }
      }
	else
	printf("No T0 data found!!!\n");

      // Fill the histograms
	
      for (Int_t ik = 0; ik<24; ik++)
         for (Int_t iHt=0; iHt<5; iHt++){
                if(allData[ik+1][iHt]!=0 ){
		  if(ik<12){
			 hCFD1minCFD[ik]->Fill(allData[ik+1][iHt]-allData[1][iHt]);
		  }
                  if(ik>11){
                         hCFD1minCFD[ik]->Fill(allData[ik+45][iHt]-allData[57][iHt]);
                  }
		}
	}
	
     delete start;
	start = 0x0;
     reader->Reset();
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
  // write a file with the histograms
  TFile *hist = new TFile(FILE_OUT,"RECREATE");

  for(Int_t j=0;j<24;j++){
     hCFD1minCFD[j]->Write();
    }
  hist->Close();
  delete hist;

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_OUT, "COSMIC")) {
    status=-2;
  }

  return status;
}

