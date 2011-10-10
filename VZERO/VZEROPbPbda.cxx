/*********************************************************************************
- Contact:    Brigitte Cheynis     b.cheynis@ipnl.in2p3.fr
- Link:
- Raw data test file :          
- Reference run number : 137366	      
- Run Type:   PHYSICS
- DA Type:    MON
- Number of events needed: >=2000
- Input Files:  argument list
- Output Files: FXS file     V0_EqualizationFactors.dat (Channel equalization factors)
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/

/**********************************************************************************
*                                                                                 *
* VZERO Detector Algorithm for extracting channel equalization factors for Pb-Pb  *
*                                                                                 *
*                                                                                 *
***********************************************************************************/

// DATE
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

//AliRoot
#include <AliVZERORawStream.h>
#include <AliRawReaderDate.h>
#include <AliRawReader.h>
#include <AliDAQ.h>

// standard
#include <stdio.h>
#include <stdlib.h>

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>

Int_t GetOfflineChannel(Int_t channel);

/* Main routine --- Arguments: monitoring data source */
      
int main(int argc, char **argv) {

/* magic line from Cvetan */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                    "*",
                    "TStreamerInfo",
                    "RIO",
                    "TStreamerInfo()");
  int status;
  if (argc!=2) {
     printf("Wrong number of arguments\n");
     return -1;
  }
  
//___________________________________________________
// Get parameters from V00DAEqualFactors.config file

  Int_t    kStartClock = 9;  // First clock in the search for max adc
  Int_t    kEndClock = 11;   // Last clock in the search for max adc
  Int_t    kNPreClocks = 6;  // Number of clock before max used in the charge sum
  Int_t    kNPostClocks = 1; // Number of clock after max used in the charge sum

  UShort_t    kTriggerAcc = 64;    // Trigger mask for accepted events (64 = CTA1 & CTC1)
  UShort_t    kTriggerRej = 256;   // Trigger mask for rejected events (256 = CTA2 & CTC2)

  Int_t    kNBins = 10000;
  Float_t  kRange = 0.1;

  status = daqDA_DB_getFile("V00DAEqualFactors.config","./V00DAEqualFactors.config");
  if (status) {
    printf("Failed to get Config file (V00DAEqualFactors.config) from DAQ DB, status=%d\n", status);
    printf("Take default values of parameters for pedestal calculation \n");
  } else {
    /* open the config file and retrieve cuts */
    FILE *fpConfig = fopen("V00DAEqualFactors.config","r");
    int res = fscanf(fpConfig,"%d %d %d %d %u %u %d %f",
		     &kStartClock,&kEndClock,&kNPreClocks,&kNPostClocks,&kTriggerAcc,&kTriggerRej,&kNBins,&kRange);
    if(res!=8) {
      printf("Failed to get values from Config file (V00DAEqualFactors.config): wrong file format - 7 integers and 1 float are expected - \n");
    }
    fclose(fpConfig);
  }
  
  printf("First LHC Clock = %d; Last LHC Clock = %d; N Pre Clock = %d ; N Post Clock = %d; Trigger mask for accepted events = %u; Trigger mask for rejected events = %u; Number of histogram bins = %d; Histogram range = %.3f\n",
	 kStartClock, kEndClock, kNPreClocks, kNPostClocks, kTriggerAcc, kTriggerRej, kNBins, kRange);

  TH1D *fMedian[64];
  for(Int_t j = 0; j < 64; ++j) fMedian[j] = new TH1D(Form("fMedian_%d",j),"Slopes weighted median, channel par channel",kNBins,0,kRange);

  Bool_t fFirst = kTRUE;
  Float_t fPrevTotCharge = 0;
  Float_t fPrevadc[64];
  for(Int_t j = 0; j < 64; ++j) fPrevadc[j] = 0;

//___________________________________________________ 
   
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
  
  /* init counters on events */
  int neventsPhysics=0;
  int neventsTotal=0;


  /* loop on events (infinite) */
  for(;;) {
      struct eventHeaderStruct *event;
      eventTypeType eventT;

      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}   

      /* get next event (blocking call until timeout) */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) {
          printf ("End of File detected\n");
      break; /* end of monitoring file has been reached */
      }
    
      if (status!=0) {
           printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
      break;
      }

      /* retry if got no event */
      if (event==NULL) continue;
               
      /* decode event */
      eventT=event->eventType;
	
      switch (event->eventType){
      
      case START_OF_RUN:
           break;
      
      case END_OF_RUN:
           printf("End Of Run detected\n");
           break;
      
      case PHYSICS_EVENT:
 	     		 
	   AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  
	   AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
	   if (rawStream->Next()) {	
	     UShort_t triggers = rawStream->GetTriggerInputs();
	     if (((triggers & kTriggerAcc) == kTriggerAcc) &&  // Check if the requested trigger(s) is fired
		 ((triggers & kTriggerRej) == 0)) { // Check if the requested trigger(s) is NOT fired
	       neventsPhysics++;

	       Float_t adc[64];
	       Float_t totCharge = 0.;
	       for(Int_t i = 0; i < 64; ++i) {
		 adc[i] = 0.;
		 Float_t maxadc=0.;
		 Int_t imax=-1;
		 for(Int_t j = kStartClock; j <= kEndClock; ++j) {
		   Float_t charge = (Float_t)(rawStream->GetPedestal(i,j));
		   if(charge > maxadc){
		     maxadc = charge;
		     imax   = j;
		   }
		 }

		 if (imax != -1) {
		   Int_t start = imax - kNPreClocks;
		   if (start < 0) start = 0;
		   Int_t end = imax + kNPostClocks;
		   if (end > 20) end = 20;
		   for(Int_t iClock = start; iClock <= end; iClock++) {
		     adc[i] += (Float_t)(rawStream->GetPedestal(i,iClock));
		   }
		 }
		 totCharge += adc[i];
	       }

	       if (fFirst) {
		 fFirst = kFALSE;
		 fPrevTotCharge = totCharge;
		 for(int i = 0; i < 64; ++i) fPrevadc[i] = adc[i];
	       }
	       else {
		 fFirst = kTRUE;
		 Float_t deltaTotCharge = totCharge - fPrevTotCharge;
		 Float_t weight = deltaTotCharge*deltaTotCharge;
		 if (weight > 1) {
		   for(int i = 0; i < 64; ++i) {
		     fMedian[i]->Fill((adc[i]-fPrevadc[i])/deltaTotCharge,weight);
		   }
		 }
	       }

	     }
	   }   // End : if rawstream
           delete rawStream;
           rawStream = 0x0;      
           delete rawReader;
           rawReader = 0x0;	     	 						         
      } // end of switch on event type 
	
      neventsTotal++;
      /* free resources */
      free(event);
	
      /* exit when last event received, no need to wait for TERM signal */
      if (eventT==END_OF_RUN) {
	printf("End Of Run event detected\n");
	break;
      }

  }  // loop over events
  
  printf("%d physics events processed\n",neventsPhysics);
    
//___________________________________________________________________________
//  Computes regression parameters
// charge_i = p0 + charge_tot * p1

  if(neventsPhysics>2000){
    /* open result file to be exported to FES */
    FILE *fp=NULL;
    fp=fopen("./V0_EqualizationFactors.dat","w");
    if (fp==NULL) {
      printf("Failed to open local result file\n");
      return -1;}

    Double_t beta[64];
    Double_t q = 0.5;
    for(int i = 0; i < 64; ++i) fMedian[i]->GetQuantiles(1,&beta[i],&q);

    for(Int_t i=0; i<64; i++) {
      fprintf(fp," %d %.3f\n",GetOfflineChannel(i), beta[i]*64.);				       
      printf(" %d %.3f\n",GetOfflineChannel(i), beta[i]*64.);				       
    }

    /* close local result file and FXS result file*/
    fclose(fp);
  }

//________________________________________________________________________
   
  /* export result file to FES */
  status=daqDA_FES_storeFile("./V0_EqualizationFactors.dat","V00DAEqualFactors");
  if (status)    {
    printf("Failed to export file : %d\n",status);
    return -1; }

  /* store result file into Online DB */
  status=daqDA_DB_storeFile("./V0_EqualizationFactors.dat","V00DAEqualFactors");
  if (status)    {
    printf("Failed to store file into Online DB: %d\n",status);
    return -1; }

  return status;
}

 Int_t GetOfflineChannel(Int_t channel) {

// Channel mapping Online - Offline:
 
 Int_t fOfflineChannel[64] = {39, 38, 37, 36, 35, 34, 33, 32, 
                              47, 46, 45, 44, 43, 42, 41, 40, 
			      55, 54, 53, 52, 51, 50, 49, 48, 
			      63, 62, 61, 60, 59, 58, 57, 56,
			       7,  6,  5,  4,  3,  2,  1,  0, 
			      15, 14, 13, 12, 11, 10,  9,  8,
			      23, 22, 21, 20, 19, 18, 17, 16, 
			      31, 30, 29, 28, 27, 26, 25, 24};
 return	fOfflineChannel[channel];		      
}			      
