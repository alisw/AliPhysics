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

#include <TSystem.h>

#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliPHOSDA2.h"
#include "AliPHOSRawDecoderv1.h"
#include "AliCaloAltroMapping.h"


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

  /* Open mapping files */
  AliAltroMapping *mapping[4];
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/PHOS/mapping/RCU";
  TString path2;
  for(Int_t i = 0; i < 4; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    mapping[i] = new AliCaloAltroMapping(path2.Data());
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
  printf("DA2 (bad channels finding) started.\n");  


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  AliRawReader *rawReader = NULL;

  AliPHOSDA1 da2(2); // DA2 ("Checking for bad channels") for module2
  
  Float_t q[64][56][2];

  Int_t gain = -1;
  Int_t X = -1;
  Int_t Z = -1;

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
      
      for(Int_t iX=0; iX<64; iX++) {
	for(Int_t iZ=0; iZ<56; iZ++) {
	  for(Int_t iGain=0; iGain<2; iGain++) {
	    q[iX][iZ][iGain] = 0.;
	  }
	}
      }

      rawReader = new AliRawReaderDate((void*)event);
      AliPHOSRawDecoderv1 dc(rawReader,mapping);
      dc.SubtractPedestals(kTRUE);
      dc.SetOldRCUFormat(kTRUE);
      
      while(dc.NextDigit()) {

	X = dc.GetRow() - 1;
	Z = dc.GetColumn() - 1;

	if(dc.IsLowGain()) gain = 0;
	else
	  gain = 1;
	
	q[X][Z][gain] = dc.GetSampleQuality();
	
      }
      
      da2.FillQualityHistograms(q);       
      //da1.UpdateHistoFile();
      
      delete rawReader;     
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

  for(Int_t i = 0; i < 4; i++) delete mapping[i];  

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);


  return status;
}
