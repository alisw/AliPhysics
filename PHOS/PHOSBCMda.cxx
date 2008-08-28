/*
contact: Boris.Polishchuk@cern.ch
link: http://aliceinfo.cern.ch/static/phpBB3/viewtopic.php?f=4&t=17
reference run: /castor/cern.ch/alice/phos/2007/10/02/13/07000008232001.10.root
run type: STANDALONE
DA type: MON
number of events needed: 1000
input files: RCU0.data  RCU1.data  RCU2.data  RCU3.data
Output files: PHOS_Module2_BCM.root
Trigger types used: CALIBRATION_EVENT or PHYSICS
*/


#include "event.h"
#include "monitor.h"
extern "C" {
#include "daqDA.h"
}

#include <stdio.h>
#include <stdlib.h>

#include <TSystem.h>
#include <TROOT.h>
#include <TPluginManager.h>

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
  printf("DA2 (bad channels search) started.\n");  


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  AliRawReader *rawReader = NULL;

  AliPHOSDA2 da2(2); // DA2 ("Checking for bad channels") for module2
  
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
    
    if (eventT==PHYSICS_EVENT || eventT==CALIBRATION_EVENT) {
      
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
  
  /* Store output files to the File Exchange Server */
  char localfile[128];

  for(Int_t iMod=0; iMod<5; iMod++) {
    sprintf(localfile,"PHOS_Module%d_BCM.root",iMod);
    daqDA_FES_storeFile(localfile,"BAD_CHANNELS");
  }

  return status;
}
