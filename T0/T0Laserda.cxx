/*
T0 DA for online calibration

Contact: Michal.Oledzki@cern.ch
Link: http://users.jyu.fi/~mioledzk/
Run Type: STANDALONE
DA Type: MON
Number of events needed: 400000 
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
//#include "TProfile.h"
int cqbx,cqby,clbx,clby,cbx;
float cqlx,cqmx,cqly,cqmy,cllx,clmx,clly,clmy,clx,cmx;
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
     printf("Couldn't get input file >>inLaser.dat<< from DAQ_DB !!!\n");
     return -1;
  }
  
 
  FILE *inp;
  char c;
  inp = fopen(FILE_IN, "r");  
  if(!inp){
	printf("Input file >>inLaser.dat<< not found !!!\n");
	return -1;
  }

  while((c=getc(inp))!=EOF) {
    switch(c) {
      case 'a': {fscanf(inp, "%d", &cqbx ); break;} //N of X bins hCFD_QTC
      case 'b': {fscanf(inp, "%f", &cqlx ); break;} //Low x hCFD_QTC
      case 'c': {fscanf(inp, "%f", &cqmx ); break;} //High x hCFD_QTC
      case 'd': {fscanf(inp, "%d", &cqby ); break;} //N of Y bins hCFD_QTC
      case 'e': {fscanf(inp, "%f", &cqly ); break;} //Low y hCFD_QTC
      case 'f': {fscanf(inp, "%f", &cqmy ); break;} //High y hCFD_QTC
      case 'g': {fscanf(inp, "%d", &clbx ); break;} //N of X bins hCFD_LED
      case 'h': {fscanf(inp, "%f", &cllx ); break;} //Low x hCFD_LED
      case 'i': {fscanf(inp, "%f", &clmx ); break;} //High x hCFD_LED
      case 'j': {fscanf(inp, "%d", &clby ); break;} //N of Y bins hCFD_LED
      case 'k': {fscanf(inp, "%f", &clly ); break;} //Low y hCFD_LED
      case 'l': {fscanf(inp, "%f", &clmy ); break;} //High y hCFD_LED
      case 'm': {fscanf(inp, "%d", &cbx ); break;}  //N of Y bins hCFD 
      case 'n': {fscanf(inp, "%f", &clx ); break;}  //Low x hCFD
      case 'o': {fscanf(inp, "%f", &cmx ); break;}  //High x hCFD
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
  TH1F *hCFD[24];
  TH2F *hCFDvsQTC[24];
  TH2F *hCFDvsLED[24]; 

   for(Int_t ic=0; ic<24; ic++) {
      hCFDvsQTC[ic] = new TH2F(Form("CFD_QTC%d",ic+1),"CFD_QTC",cqbx,cqlx,cqmx,cqby,cqly,cqmy);
      hCFDvsLED[ic] = new TH2F(Form("CFD_LED%d",ic+1),"CFD_LED",clbx,cllx,clmx,clby,clly,clmy);
      if(ic<12){
	hCFD[ic] = new TH1F(Form("T0_C_%d_CFD",ic+1),"CFD", cbx,clx,cmx);	
      }
      else{
        hCFD[ic] = new TH1F(Form("T0_A_%d_CFD",ic-11),"CFD", cbx,clx,cmx); 
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
		if((allData[cc+1][iHt]-allData[0][0]+5000)!=0 && allData[cc+1][iHt]>0){
		 hCFD[cc]->Fill(allData[cc+1][iHt]-allData[0][0]+5000);
		}
		if(allData[ik+25][iHt]!=0 && allData[ik+26][iHt]!=0 && allData[cc+1][iHt]!=0){
                 hCFDvsQTC[cc]->Fill((allData[ik+25][iHt]-allData[ik+26][iHt]) , (allData[cc+1][iHt]-allData[0][0]+5000));
		} 
                if(allData[cc+13][iHt]!=0 && allData[cc+1][iHt]!=0){
		 hCFDvsLED[cc]->Fill(allData[cc+13][iHt]-allData[cc+1][iHt],allData[cc+1][iHt]-allData[0][0]+5000);
		}
	}

      for (Int_t ik = 24; ik<48; ik+=2)
         for (Int_t iHt=0; iHt<5; iHt++){
		 Int_t cc = ik/2;
                if((allData[cc+45][iHt]-allData[0][0]+5000)!=0 && allData[cc+45][iHt]>0){
                 hCFD[cc]->Fill(allData[cc+45][iHt]-allData[0][0]+5000);
                }
                if(allData[ik+57][iHt]!=0 && allData[ik+58][iHt]!=0 && allData[cc+45][iHt]!=0){
                 hCFDvsQTC[cc]->Fill(allData[ik+57][iHt]-allData[ik+58][iHt],allData[cc+45][iHt]-allData[0][0]+5000);
		}
                if(allData[cc+57][iHt]!=0 && allData[cc+45][iHt]!=0){
                 hCFDvsLED[cc]->Fill(allData[cc+57][iHt]-allData[cc+45][iHt],allData[cc+45][iHt]-allData[0][0]+5000);
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
     hCFDvsQTC[j]->Write();
     hCFDvsLED[j]->Write();
     hCFD[j]->Write();
    }
  hist->Close();
  delete hist;

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_OUT, "LASER")) {
    status=-2;
  }

  return status;
}

