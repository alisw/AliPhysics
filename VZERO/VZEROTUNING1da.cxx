/*********************************************************************************
- Contact:    Brigitte Cheynis     b.cheynis@ipnl.in2p3.fr
- Link:       http
- Raw data test file : 
- Reference run number : 	      
- Run Type:   STANDALONE
- DA Type:    LDC
- Number of events needed: 500
- Input Files:  argument list
- Output Files: local file  V0_Tuning1.dat
                FXS file    V0_Tuning1.dat
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/


/**********************************************************************************
*                                                                                 *
* VZERO Detector Algorithm used for tuning FEE parameters                         *
*                                                                                 *
* This program reads data on the LDC                                              *
* It cumulates fBB and fBG flags, populates local "./V0_Tuning1.dat"              *            
* file and exports it to the FES.                                                 *
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


/* Main routine --- Arguments: data file */
      
int main(int argc, char **argv) {

/* magic line from Cvetan */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                    "*",
                    "TStreamerInfo",
                    "RIO",
                    "TStreamerInfo()");
  int status;

  Float_t BBFlag[128];
  Float_t BGFlag[128];
  for(Int_t i=0; i<128; i++) {
      BBFlag[i] = 0.0;
      BGFlag[i] = 0.0;
  } 
      
  /* log start of process */
  printf("VZERO DA program started - Tuning FEE parameters\n");  

  /* check that we got some arguments = list of files */
  if (argc<2)   {
      printf("Wrong number of arguments\n");
      return -1;}

  /* open result file to be exported to FES */
  FILE *fp=NULL;
  fp=fopen("./V0_Tuning1.dat","a");
  if (fp==NULL) {
      printf("Failed to open result file\n");
      return -1;}

  /* open log file to inform user */
  FILE *flog=NULL;
  flog=fopen("./V00log.txt","a");
  if (flog==NULL) {
      printf("Failed to open log file\n");
      return -1;  }
    
  /* report progress */
  daqDA_progressReport(10);

  /* init counters on events */
  int nevents_physics=0;
  int nevents_total=0;

  /* read the data  */
  
  status=monitorSetDataSource( argv[1] );
  if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1; }

  /* report progress */
  daqDA_progressReport(10+80/argc);

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
             break;
      
        case END_OF_RUN:
	     printf("End Of Run detected\n");
             break;
      
        case PHYSICS_EVENT:
             nevents_physics++;
 	     
//              fprintf(flog,"Run #%lu, event size: %lu, BC:%u, Orbit:%u, Period:%u\n",
//                  (unsigned long)event->eventRunNb,
//                  (unsigned long)event->eventSize,
//                  EVENT_ID_GET_BUNCH_CROSSING(event->eventId),
//                  EVENT_ID_GET_ORBIT(event->eventId),
//                  EVENT_ID_GET_PERIOD(event->eventId) );
		 
	     AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  
	     AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
	     rawStream->Next();	
	     for(Int_t i=0; i<64; i++) {
	        if(!rawStream->GetIntegratorFlag(i,10))
		{
                   if(rawStream->GetBBFlag(i,10)) BBFlag[i]++;       // even integrator - channel 0 to 63
		   if(rawStream->GetBGFlag(i,10)) BGFlag[i]++;       // even integrator - channel 0 to 63
		} 
		else 
		{ 		   
		   if(rawStream->GetBBFlag(i,10)) BBFlag[i+64]++;    // odd integrator  - channel 64 to 123
		   if(rawStream->GetBGFlag(i,10)) BGFlag[i+64]++;    // odd integrator  - channel 64 to 123
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

  }  // loop over events
    
//________________________________________________________________________
//  Computes mean values, dumps them into the output text file
	
  for(Int_t i=0; i<128; i++) {
      BBFlag[i] = BBFlag[i]/nevents_physics;
      BGFlag[i] = BGFlag[i]/nevents_physics;
      fprintf(fp," %d %f %f\n",i,BBFlag[i],BGFlag[i]);
      fprintf(flog," %d %f %f\n",i,BBFlag[i],BGFlag[i]);
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
  status=daqDA_FES_storeFile("./V0_Tuning1.dat","V00da_results");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }

  /* report progress */
  daqDA_progressReport(100);
  
  return status;
}
