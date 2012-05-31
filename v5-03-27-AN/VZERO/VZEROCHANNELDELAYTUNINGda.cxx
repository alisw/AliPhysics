/*********************************************************************************
- Contact:    Brigitte Cheynis     b.cheynis@ipnl.in2p3.fr
- Link:       http
- Raw data test file : 
- Reference run number : 	      
- Run Type:   CHANNEL_DELAY_TUNING
- DA Type:    LDC
- Number of events needed: 500
- Input Files:  argument list
- Output Files: local file  V0_ChannelDelayTuning.dat
                FXS file    V0_ChannelDelayTuning.dat
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/


/**********************************************************************************
*                                                                                 *
* VZERO Detector Algorithm used for tuning FEE parameters                         *
*                                                                                 *
* This program reads data on the LDC                                              *
* It cumulates fBB and fBG flags, populates local "./V0_ChannelDelayTuning.dat"   *            
* file, exports it to the FES, and stores it to DAQ DB                            *
* We have 128 channels instead of 64 as expected for V0 due to the two sets of    *
* charge integrators which are used by the FEE ...                                *
* The program reports about its processing progress.                              *
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


/* Main routine --- Arguments: list of DATE raw data files */
      
int main(int argc, char **argv) {

/* magic line from Cvetan */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                    "*",
                    "TStreamerInfo",
                    "RIO",
                    "TStreamerInfo()");
  int status;
  status = daqDA_DB_getFile("V00da_results","./V0_Pedestals.dat");
  if (status) {
      printf("Failed to get Pedestal file (V0_Pedestals.dat) from DAQ DB, status=%d\n", status);
      return -1;   
  }
  
  Float_t MeanPed[128], SigPed[128], fdump;
  
  /* open the pedestal file and retrieve pedestal mean and sigma */
  FILE *fpPed = fopen("V0_Pedestals.dat","r");
  for(int i=0;i<128;i++){
      fscanf(fpPed,"%f %f %f %f \n",&MeanPed[i],&SigPed[i],&fdump,&fdump);
//      printf("%.3f %.3f \n",MeanPed[i],SigPed[i]);
  }
  fclose(fpPed);

//______________________________________________________________________________
// Get running parameters from V00_ChannelDelayTuning_DA.config file

  Int_t    kNbEventSkipped;   // = 100;   number of events skipped - to be tuned
  Float_t  kSigmaCut;         // = 3.0;   number of sigmas for threshold cut - to be tuned
  
  status = daqDA_DB_getFile("V00_ChannelDelayTuning_DA.config","./V00_ChannelDelayTuning_DA.config");
  if (status) {
      printf("Failed to get config file (V00_ChannelDelayTuning_DA.config) from DAQ DB, status=%d\n", status);
      return -1;   
  }
  /* open the config file and retrieve running parameters */
  FILE *fpConfig = fopen("V00_ChannelDelayTuning_DA.config","r");
  fscanf(fpConfig,"%f %d",&kSigmaCut,&kNbEventSkipped);
  fclose(fpConfig);
  
  printf("Number of events skipped = %d ; Number of sigmas for threshold cut = %f\n",kNbEventSkipped,kSigmaCut);
//______________________________________________________________________________

  Int_t  BBFlag[64];
  Int_t  BGFlag[64];
  Int_t  ChargeEoI  = 0;
  Bool_t Integrator = 0;
  Int_t NHit[64];
  for(Int_t i=0; i<64; i++) {
      BBFlag[i] = 0;
      BGFlag[i] = 0;
      NHit[i]   = 0;
  } 
      
  /* log start of process */
  printf("VZERO DA program started - Channel Delay Tuning \n");  

  /* check that we got some arguments  */
  if (argc<2)   {
      printf("Wrong number of arguments\n");
      return -1;}

  /* open result file to be exported to FES */
  FILE *fp=NULL;
  fp=fopen("./V0_ChannelDelayTuning.dat","a");
  if (fp==NULL) {
      printf("Failed to open result file\n");
      return -1;}

  /* open log file to inform user */
  FILE *flog=NULL;
  flog=fopen("./V00log.txt","w");
  if (flog==NULL) {
      printf("Failed to open log file\n");
      return -1;  }
    
  /* report progress */
  daqDA_progressReport(10);

  /* init counters on events */
  int nevents_physics=0;
  int nevents_total=0;
  int iteration;
  sscanf(argv[1],"%d",&iteration);
  
  /* read the n data files */
  for (int n=2; n<argc; n++) {
  
  /* read the data  */
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
        printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
        return -1; }

  /* report progress */
    daqDA_progressReport(10+50*n/argc);

  /* read the data file */
    for(;;) {
        struct eventHeaderStruct *event;
        eventTypeType eventT;

        /* get next event */
        status=monitorGetEventDynamic((void **)&event);
        if (status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
        if (status!=0) {
            printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
            return -1; }

        /* retry if got no event */
        if (event==NULL) break;
        
        /* decode event */
        eventT=event->eventType;
	
        switch (event->eventType){
      
        case START_OF_RUN:
	     printf("START Of Run detected\n");
             break;
      
        case END_OF_RUN:
	     printf("END Of Run detected\n");
             break;
      
        case PHYSICS_EVENT:
             nevents_physics++;

	     AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  
	     AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
	     rawStream->Next();	
	     if(nevents_physics > kNbEventSkipped){
		for(Int_t i=0; i<64; i++) {
	            //ChargeEoI= rawStream->GetADC(i); takes EoI as central clock
		    // Look for the maximum in the LHC clock train instead of central clock 10
                    ChargeEoI    = 0;
                    Int_t iClock = 0;
                    for(size_t iEvent=0; iEvent<21; iEvent++){
                         if(rawStream->GetPedestal(i,iEvent)>ChargeEoI) 
	                    {ChargeEoI= rawStream->GetPedestal(i,iEvent);
	                    iClock    = iEvent;}
                    }   		    
		    Integrator = rawStream->GetIntegratorFlag(i,iClock);
		    Float_t Threshold = MeanPed[i + 64 * Integrator] + kSigmaCut * SigPed[i + 64 * Integrator];

		    if((float)ChargeEoI>Threshold) {
			NHit[i]+=1;
	        	if(rawStream->GetBBFlag(i,iClock)) BBFlag[i]+=1;       
			if(rawStream->GetBGFlag(i,iClock)) BGFlag[i]+=1; 
		    }		 
          	}    
	     }
             delete rawStream;
             rawStream = 0x0;      
             delete rawReader;
             rawReader = 0x0;	     	 						         
        } // end of switch on event type 
	
        nevents_total++;
        /* free resources */
        free(event);

    }    // end of loop over events
  }      // end of loop over data files 
   
//________________________________________________________________________
//  Computes mean values, dumps them into the output text file
  Float_t fBB, fBG;	
  for(Int_t i=0; i<64; i++) {      
      if(NHit[i] > 0) {
         fBB = (float)BBFlag[i]/(float)NHit[i];
         fBG = (float)BGFlag[i]/(float)NHit[i];
       }else{
      	 fBB = fBG = 0.;
       }
       fprintf(fp," %d %d %f %f\n",iteration,i,fBB,fBG);
       fprintf(flog," it %d ch %d BB %f BG %f BBFlag %d BGFlag %d Hit %d\n",iteration,i,fBB,fBG,BBFlag[i],BGFlag[i],NHit[i]);
  } 
  
//________________________________________________________________________
   
  /* write report */
  fprintf(flog,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);
  printf("Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);
  
  /* close result and log files */
  fclose(fp);
  fclose(flog); 
  
  /* report progress */
  daqDA_progressReport(90);

  /* export result file to FES */
  status=daqDA_FES_storeFile("./V0_ChannelDelayTuning.dat","V00da_ChannelDelayTuning");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }
      
  /* store result file on Online DB */
//   status=daqDA_DB_storeFile("./V0_ChannelDelayTuning.dat","V00da_ChannelDelayTuning");
//   if (status)    {
//       printf("Failed to store file into Online DB: %d\n",status);
//       return -1; }
      
  /* report progress */
  daqDA_progressReport(100);
  
  return status;
}
