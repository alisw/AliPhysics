/*

DAcase2.c

This program connects to the DAQ data source passed as argument
and populates local "./result.txt" file with the ids of events received
during the run.

The program exits when being asked to shut down (daqDA_checkshutdown)
or End of Run event.

Messages on stdout are exported to DAQ log system.

contact: alice-datesupport@cern.ch

*/
extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"
//#include "daqDA.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliRawReaderDate.h"
#include "AliPMDCalibPedestal.h"
#include "AliPMDCalibGain.h"

//ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TTree.h"

/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  
  AliPMDCalibPedestal calibped;
  AliPMDCalibGain calibgain;

  TTree *ped  = new TTree("ped","PMD Pedestal tree");
  TTree *gain = new TTree("gain","PMD Gain tree");

  TH1F::AddDirectory(0);
  
      
  // decoding the events
  
  int status;

  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
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
  printf("DA example case2 monitoring program started\n");  

  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  struct eventHeaderStruct *event;
  eventTypeType eventT;
  Int_t iev=0;

  /* main loop (infinite) */
  for(;;) {
    
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
    if (event==NULL) {
      continue;
    }

    iev++; 

   /* use event - here, just write event id to result file */
    eventT=event->eventType;
    switch (event->eventType){
      
      /* START OF RUN */
    case START_OF_RUN:
      break;
      /* END START OF RUN */
      
    /* END OF RUN */
    case END_OF_RUN:
      break;
      
    case PHYSICS_EVENT:
      printf(" event number = %i \n",iev);
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      calibped.ProcessEvent(rawReader);

      calibgain.ProcessEvent(rawReader);

      delete rawReader;
      rawReader = 0x0;

    }
       
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      calibped.Analyse(ped);
      calibgain.Analyse(gain);

      break;
    }
  }
  
  //write the Run level file   
  TFile * fileRun = new TFile ("outPMDdaRun.root","RECREATE"); 
  TBenchmark *bench = new TBenchmark();
  bench->Start("PMD");
  bench->Stop("PMD");
  bench->Print("PMD");
  fileRun->Close();

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);


  TFile * pedRun = new TFile ("pmd_ped.root","RECREATE"); 
  ped->Write();
  pedRun->Close();

  TFile * gainRun = new TFile ("pmd_calib.root","RECREATE"); 
  gain->Write();
  gainRun->Close();



  /* close result file */
  fclose(fp);


  return status;
}
