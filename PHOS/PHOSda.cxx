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


#include "event.h"
#include "monitor.h"
extern "C" {
#include "daqDA.h"
}

#include <stdio.h>
#include <stdlib.h>

#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliPHOSCalibHistoProducer.h"
#include "AliPHOSRawDecoder.h"


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

  AliRawReader *rawReader = NULL;

  AliPHOSCalibHistoProducer hp(200,0.,200.);
  hp.SetOldRCUFormat(kTRUE);
  hp.SetUpdatingRate(200000);
  
  /* main loop (infinite) */
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
    if (event==NULL) {
      continue;
    }


    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    if (eventT==PHYSICS_EVENT) {
      fprintf(fp,"Run #%lu, event size: %lu, BC:%u, Orbit:%u, Period:%u\n",
	  (unsigned long)event->eventRunNb,
        (unsigned long)event->eventSize,
        EVENT_ID_GET_BUNCH_CROSSING(event->eventId),
        EVENT_ID_GET_ORBIT(event->eventId),
        EVENT_ID_GET_PERIOD(event->eventId)
      );
      
      rawReader = new AliRawReaderDate((void*)event);
      AliPHOSRawDecoder dc(rawReader);
      dc.SubtractPedestals(kTRUE);
      hp.SetRawDecoder(&dc);
      hp.Run();
      
      nevents_physics++;
    }
    nevents_total++;


    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }


  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);


  return status;
}
