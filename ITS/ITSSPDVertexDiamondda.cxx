/*
Contact: cvetan.cheshkov@cern.ch
Link: missing
Run Type: PHYSICS
DA Type: MON
Number of events needed: 10000
Input Files:
Output Files:
Trigger types used: PHYSICS
*/

#define OUTPUT_FILE "SPDVertexDiamondDA.root"
#define CDB_STORAGE "local://$ALICE_ROOT"
#define N_EVENTS_AUTOSAVE 50

extern "C" {
#include "daqDA.h"
}

#include "event.h"
#include "monitor.h"

#ifdef ALI_AMORE
#include <AmoreDA.h>
//int amore::da::Updated(char const*) {}
#endif

#include <TPluginManager.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>

#include "AliRawReaderDate.h"
#include "AliCDBManager.h"
#include "AliITSMeanVertexer.h"

int main(int argc, char **argv) {

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()"); 

  int status;
  if (argc<2) {
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
  printf("Vertex-Diamond SPD DA started\n");  

  /* init some counters */
  int nevents_with_vertex = 0;
  int nevents_physics=0;
  int nevents_total=0;

  struct eventHeaderStruct *event;
  eventTypeType eventT;

  // Get run number
  if (getenv("DATE_RUN_NUMBER")==0) {
    printf("DATE_RUN_NUMBER not properly set.\n");
    return -1;
  }
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));

  // Global initializations
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(CDB_STORAGE);
  man->SetRun(runNr);

  // Init mean vertexer
  AliITSMeanVertexer *mv = new AliITSMeanVertexer();
  if (!mv->Init()) {
    printf("Initialization of mean vertexer object failed ! Check the log for details");
    return -1;
  }

  // Initialization of AMORE sender
#ifdef ALI_AMORE
  amore::da::AmoreDA vtxAmore(amore::da::AmoreDA::kSender);
#endif
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

    nevents_total++;
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
      nevents_physics++;
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);

      // Run mean-vertexer reco
      if (mv->Reconstruct(rawReader)) nevents_with_vertex++;

      // Auto save
      if ((nevents_physics%N_EVENTS_AUTOSAVE) == 0)
	mv->WriteVertices(OUTPUT_FILE);

      delete rawReader;
    }
    
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }

  mv->WriteVertices(OUTPUT_FILE);

#ifdef ALI_AMORE
  // send the histos to AMORE pool
  printf("AMORE send status: %d\n",vtxAmore.Send(mv->GetVertexXY()->GetName(),mv->GetVertexXY()));
  printf("AMORE send status: %d\n",vtxAmore.Send(mv->GetVertexZ()->GetName(),mv->GetVertexZ()));
#endif

  delete mv;

  /* write report */
  printf("Run #%s, received %d events with vertex, out of %d physics and out of %d total events\n",getenv("DATE_RUN_NUMBER"),nevents_with_vertex,nevents_physics,nevents_total);

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(OUTPUT_FILE, "VertexDiamond")) {
    status=-2;
  }
  
  return status;
}
