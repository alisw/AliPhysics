/*

PHOSPHYSda.cxx

This program reads the DAQ data files passed as argument using the monitoring library.
It stores the required event information into an output file.
Messages on stdout are exported to DAQ log system.

contact: Hisayuki.Torii@cern.ch

*/

#include "event.h"
#include "monitor.h"
#include "daqDA.h"

#include <stdio.h>
#include <stdlib.h>

#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliPHOSRawDecoder.h"
#include "AliCaloAltroMapping.h"
#include "AliPHOSDApi0mip.h"

/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

  int status;

  /* log start of process */
  printf("PHOSPHYSda program started\n");  

  /* Time out counters*/
  int timeout = 180; // 3mins
  time_t start_time = time(0);

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* Retrieve mapping files from DAQ DB */ 
  const char* mapFiles[4] = {"RCU0.data","RCU1.data","RCU2.data","RCU3.data"};
  for(Int_t iFile=0; iFile<4; iFile++) {
    int failed = daqDA_DB_getFile(mapFiles[iFile], mapFiles[iFile]);
    if(failed) { 
      printf("Cannot retrieve file %s from DAQ DB. Exit.\n",mapFiles[iFile]);
      return -1;
    }
  }

  /* Open mapping files */
  AliAltroMapping *mapping[4];
  TString path = "./";
  path += "RCU";
  TString path2;
  for(Int_t i = 0; i < 4; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    mapping[i] = new AliCaloAltroMapping(path2.Data());
  }

  /* report progress */
  daqDA_progressReport(10);

  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  /* init PHOS reader*/
  AliRawReader *rawReader = NULL;

  /* Prepare DA algorithm */
  int iMod = 0;
  AliPHOSDApi0mip* dapi0mip = new AliPHOSDApi0mip(iMod);

  /* read the data files */
  for(int n=1;n<argc;n++) {
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* report progress */
    /* in this example, indexed on the number of files */
    daqDA_progressReport(10+80*(n-1)/argc);

    /* read the file */
    for(;;) {
      struct eventHeaderStruct *event;
      eventTypeType eventT;

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
      if (status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if (event==NULL) {
        break;
      }

      /* Status report */
      // This is just for debugging
      //if( nevents_total%1000 == 0 ){
      //daqDA_progressReport(10+80*(n-1)/argc+nevents_total/10000);
      //printf(" DEBUGDEBUG:: next event %d\n",nevents_total);
      //}
      
      /* get time stamp */
      time_t t;
      t = event->eventTimestamp;
      dapi0mip->SetTime(t);

      /* Applying into DA algorithm */
      int ix, iz, igain;
      eventT=event->eventType;
      if (eventT==PHYSICS_EVENT || eventT==CALIBRATION_EVENT) {
	// ---------------------------------------------------------------
	// User Defined Function
	// ---------------------------------------------------------------
	dapi0mip->NewEvent();
	rawReader = new AliRawReaderDate((void*)event);
	//AliPHOSRawDecoderv1 dc(rawReader,mapping);
	AliPHOSRawDecoder dc(rawReader,mapping);
	dc.SubtractPedestals(false);
	while(dc.NextDigit()) {
	  //printf(".");
	  ix = dc.GetRow() - 1;
	  iz = dc.GetColumn() - 1;
	  if(dc.IsLowGain()) igain = 0; else igain = 1;
	  if( igain == 1 ){
	    dapi0mip->FillDigit(dc.GetEnergy(),ix,iz);
	  }
	  //dc.GetTime();
	}
	delete rawReader;
	//dapi0mip->FillHist(); // This is not necesarry in this DA algorithm
	dapi0mip->FillTree();
	//if( nevents_total%1000 == 0 ){ // Dump for debugging
	//dapi0mip->Print();
	//}
	// ---------------------------------------------------------------
        nevents_physics++;
      }
      nevents_total++;

      /* Check the time out */
      if( nevents_total%1000 == 0 ){
	time_t current_time = time(0);
	if( current_time - start_time > timeout ){
	  free(event);
	  printf(" Warning: Exit due to the processing time exceed the limitation of %d sec\n",timeout);
	  n = argc;
	  break;
	}
      }
      
      /* free resources */
      free(event);
    }
    
  }

  /* report progress */
  daqDA_progressReport(90);

  /* Delete DA algorithm for saving output into a file */
  delete dapi0mip;

  /* Store output files to the File Exchange Server */
  char localfile[1024];
  char remotefile[1024];
  sprintf(localfile,"AliPHOSDApi0mip_mod%d.root",iMod);
  sprintf(remotefile,"PHOSDApi0mip",iMod);
  daqDA_FES_storeFile(localfile,remotefile);

  /* report progress */
  daqDA_progressReport(100);
  
  return status;
}

