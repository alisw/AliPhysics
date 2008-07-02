/*********************************************************************************
- Contact:    Brigitte Cheynis     b.cheynis@ipnl.in2p3.fr
- Link:       http
- Raw data test file : 
- Reference run number : 	      
- Run Type:   STANDALONE
- DA Type:    LDC
- Number of events needed: 500
- Input Files:  argument list
- Output Files: local file  V0_Tuning2.dat
                FXS file    V0_Tuning2.dat
- Trigger types used: PHYSICS_EVENT
**********************************************************************************/


/**********************************************************************************
*                                                                                 *
* VZERO Detector Algorithm used for tuning FEE parameters                         *
*                                                                                 *
* This program reads data on the LDC                                              *
* It cumulates mean ADC responses ,  populates local  "./V0_Tuning2.dat"  file    * 
* and exports it to the FES.                             *                        *                    *
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

  Float_t ADC_Mean[128];
  for(Int_t i=0; i<128; i++) {
      ADC_Mean[i] = 0.0;
  } 
      
  /* log start of process */
  printf("VZERO DA program started - Tuning FEE parameters\n");  

  /* check that we got some arguments = list of files */
  if (argc<2)   {
      printf("Wrong number of arguments\n");
      return -1;}
               
  /* open result file to be exported to FES */
  FILE *fp=NULL;
  fp=fopen("./V0_Tuning2.dat","a");
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
  
  /* read the n data files */
  for (int n=1; n<argc; n++) {
  
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
             break;
      
        case END_OF_RUN:
	     printf("End Of Run detected\n");
             break;
      
        case PHYSICS_EVENT:
             nevents_physics++;
 	     		 
	     AliRawReader *rawReader = new AliRawReaderDate((void*)event);
  
	     AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader); 
	     rawStream->Next();	
	     for(Int_t i=0; i<64; i++) {
	        if(!rawStream->GetIntegratorFlag(i,10))
		{
		    ADC_Mean[i]=ADC_Mean[i]+rawStream->GetADC(i);       // even integrator 
		}
		else 
		{ 		   
		    ADC_Mean[i+64]=ADC_Mean[i+64]+rawStream->GetADC(i); // odd integrator
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
	
  for(Int_t i=0; i<128; i++) {
      ADC_Mean[i]=ADC_Mean[i]/nevents_physics;
      fprintf(fp," %d %f\n",i,ADC_Mean[i]);
      fprintf(flog," %d %f\n",i,ADC_Mean[i]);
  } 
  fprintf(fp," End of current iteration \n");
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
  status=daqDA_FES_storeFile("./V0_Tuning2.dat","V00da_results");
  if (status)    {
      printf("Failed to export file : %d\n",status);
      return -1; }

  /* report progress */
  daqDA_progressReport(100);
  
  return status;
}
