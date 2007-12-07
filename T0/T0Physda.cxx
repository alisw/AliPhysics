/*
T0 DA for online calibration

Contact: Michal.Oledzki@cern.ch
Link: http://users.jyu.fi/~mioledzk/
Run Type: PHYSICS
DA Type: MON
Number of events needed: 500000 
Input Files: inPhys.dat, external parameters
Output Files: daPhys.root, to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT

*/

#define FILE_OUT "daPhys.root"
#define FILE_IN "inPhys.dat"
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
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
int cbx, ccbx;
float clx,cmx,cclx,ccmx;

/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
//int main(){
  int status;

  /* magic line */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  
  FILE *inp;
  char c;
  inp = fopen(FILE_IN, "r");
  while((c=getc(inp))!=EOF) {
    switch(c) {
      case 'a': {fscanf(inp, "%d", &ccbx ); break;} //N of X bins hCFD1_CFD
      case 'b': {fscanf(inp, "%f", &cclx ); break;} //Low x hCFD1_CFD
      case 'c': {fscanf(inp, "%f", &ccmx ); break;} //High x hCFD1_CFD
      case 'd': {fscanf(inp, "%d", &cbx ); break;} //N of X bins hCFD
      case 'e': {fscanf(inp, "%f", &clx ); break;} //Low x hCFD
      case 'f': {fscanf(inp, "%f", &cmx ); break;} //High x hCFD
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

  TH1F *hCFD1_CFD[24];  
  TH1F *hCFD[24];
   
   for(Int_t ic=0; ic<24; ic++) {
      if(ic<12) {
        hCFD1_CFD[ic] = new TH1F(Form("CFD1-CFD%d",ic+1),"CFD-CFD",ccbx,cclx,ccmx);
        hCFD[ic] = new TH1F(Form("CFD%i",ic+1),"CFD",cbx,clx,cmx);
	}
      if(ic>11){
        hCFD1_CFD[ic] = new TH1F(Form("CFD13-CFD%i",ic+1),"CFD-CFD",ccbx,cclx,ccmx);
        hCFD[ic] = new TH1F(Form("CFD%i",ic+1),"CFD",cbx,clx,cmx);
	}
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
 
    iev++; 

    /* use event - here, just write event id to result file */
    eventT=event->eventType;
   
    switch (event->eventType){

      case START_OF_RUN:
	break;

      case END_OF_RUN:
	break;

      case PHYSICS_EVENT:
      printf(" event number = %i \n",iev);


      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
          
      // Enable the following two lines in case of real-data
   //       reader->LoadEquipmentIdsMap("T0map.txt");
    //     reader->RequireHeader(kFALSE);
     // 	  reader->RequireHeader(kTRUE);
      AliT0RawReader *start = new AliT0RawReader(reader, kTRUE);

      // Read raw data
      Int_t allData[105][5];
      for(Int_t i0=0;i0<105;i0++)
      	for(Int_t j0=0;j0<5;j0++)
		allData[i0][j0] = 0;
 
     
      if(start->Next())
      for (Int_t i=0; i<105; i++) {
	for(Int_t iHit=0;iHit<5;iHit++){
	  allData[i][iHit]= start->GetData(i,iHit);
        }
      }
	else 
	printf("No T0 data found!!\n");

      // Fill the histograms
	
      for (Int_t ik = 0; ik<24; ik+=2)
         for (Int_t iHt=0; iHt<5; iHt++){
		 Int_t cc = ik/2;
                if(allData[cc+1][iHt]!=0 ){
		 hCFD1_CFD[cc]->Fill(allData[cc+1][iHt]-allData[1][iHt]);
		 hCFD[cc]->Fill(allData[cc+13][iHt]);
		}
	}

      for (Int_t ik = 24; ik<48; ik+=2)
         for (Int_t iHt=0; iHt<5; iHt++){
		 Int_t cc = ik/2; 
                if(allData[cc+45][iHt]!=0 ){
                 hCFD1_CFD[cc]->Fill(allData[cc+45][iHt]-allData[57][iHt]);
 		 hCFD[cc]->Fill(allData[cc+45][iHt]);
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
      break;
    }
  }
  printf("After loop, before writing histos\n");
  // write a file with the histograms
  TFile *hist = new TFile(FILE_OUT,"RECREATE");

  for(Int_t j=0;j<24;j++){
     hCFD1_CFD[j]->Write();
     hCFD[j]->Write();	
    }
  hist->Close();
  delete hist;

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_OUT, FILE_OUT)) {
    status=-2;
  }

  return status;
}


