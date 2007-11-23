//extern "C" 

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
#include <AliCDBManager.h>

//ROOT
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
#include "TProfile.h"

/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  int status;
  
  
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
  TH2F *hCFD_QTC[24]; TH2F *hCFD_LED[24]; 
  Char_t  buf1[10], buf2[10],buf3[10],buf4[10];

   for(Int_t ic=0; ic<24; ic++) {
      sprintf(buf1,"CFD_QTC%i",ic+1);
      sprintf(buf2,"CFD_LED%i",ic+1);

      hCFD_QTC[ic] = new TH2F(buf1,"CFD_QTC",8000,0.5,8000.5,2000,10000.5,20000.5);
      hCFD_LED[ic] = new TH2F(buf2,"CFD_LED",1000,-500.0,500.0,100,14600.0,14700.0);
      if(ic<12){
	sprintf(buf3,"T0_C_%i_CFD",ic+1);
	hCFD[ic] = new TH1F(buf3,"CFD", 30000,0., 30000.);	
      }
      else{
        sprintf(buf4,"T0_A_%i_CFD",ic-11);
        hCFD[ic] = new TH1F(buf4,"CFD", 30000,0., 30000.); 
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

//      case CALIBRATION_EVENT:
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
		if((allData[cc+1][iHt]-allData[0][0])>0){
		 hCFD[cc]->Fill(allData[cc+1][iHt]-allData[0][0]);
		}
		if(allData[ik+25][iHt]!=0 && allData[ik+26][iHt]!=0 && allData[cc+1][iHt]!=0){
                 hCFD_QTC[cc]->Fill((allData[ik+25][iHt]-allData[ik+26][iHt]) , (allData[cc+1][iHt]-allData[0][0]));
		 } 
                if(allData[cc+13][iHt]!=0 && allData[cc+1][iHt]!=0){
		 hCFD_LED[cc]->Fill(allData[cc+13][iHt]-allData[cc+1][iHt],allData[cc+1][iHt]-allData[0][0]);
		 }
	}

      for (Int_t ik = 24; ik<48; ik+=2)
         for (Int_t iHt=0; iHt<5; iHt++){
		 Int_t cc = ik/2;
                if((allData[cc+45][iHt]-allData[0][0])>0){
                 hCFD[cc]->Fill(allData[cc+45][iHt]-allData[0][0]);
                }

                if(allData[ik+57][iHt]!=0 && allData[ik+58][iHt]!=0){
                if(allData[cc+45][iHt]!=0)
                 hCFD_QTC[cc]->Fill(allData[ik+57][iHt]-allData[ik+58][iHt],allData[cc+45][iHt]-allData[0][0]);
		 }
                if(allData[cc+57][iHt]!=0 && allData[cc+45][iHt]!=0){
                 hCFD_LED[cc]->Fill(allData[cc+57][iHt]-allData[cc+45][iHt],allData[cc+45][iHt]-allData[0][0]);
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
  Char_t filehist[21];
  sprintf(filehist,"daLaser.root");
  TFile *hist = new TFile(filehist,"RECREATE");

  for(Int_t j=0;j<24;j++){
     hCFD_QTC[j]->Write();
     hCFD_LED[j]->Write();
     hCFD[j]->Write();
    }
  hist->Close();
  delete hist;

  return status;
}

