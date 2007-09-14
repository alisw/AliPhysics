/*
*********************************************************
Author:                                                 *
this file provides the detector algorithm for HMPID.    *
*********************************************************
*/
extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliHMPIDRawStream.h"
#include "AliHMPIDCalib.h"
#include "AliRawReaderDate.h"
#include "AliBitPacking.h"
#include "TMath.h"

//ROOT
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH2S.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TTreePlayer.h"


int main(int argc, char **argv){ 

  int status;

  /* log start of process */
  printf("HMPID DA program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* copy locally a file from daq detector config db
  status=daqDA_DB_getFile("myconfig","./myconfig.txt");
  if (status) {
    printf("Failed to get config file : %d\n",status);
    return -1;
  }
  and possibly use it */

  /* report progress */
  daqDA_progressReport(10);


  /* init the pedestal calculation */
  AliHMPIDCalib *pCal=new AliHMPIDCalib();
  //pCal->Init();                    //Init the pedestal calculation
  
  /* init event counter */
  Int_t iEvtNcal=0;

  int n;
  for (n=1;n<argc;n++) {

    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* report progress */
    /* in this example, indexed on the number of files */
    daqDA_progressReport(10+80*n/argc);

    for(;;) { // infinite loop 
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

      /* use event - here, just write event id to result file */
      eventT=event->eventType;

      if (eventT==PHYSICS_EVENT) {                                                //we use PHYSICS_EVENT for pedestal not CALIBRATION_EVENT
	
	iEvtNcal++;

	AliRawReader *reader = new AliRawReaderDate((void*)event);

	// Temporary there. Shall be removed as soon as the equipment ID is set correctly
	// For the moment ddl.map file contains one line which maps
	// the observed eqID=225 to the first HMPID DDL with ID=1536
	//	reader->LoadEquipmentIdsMap("ddl.map");
	
	AliHMPIDRawStream stream(reader);

	while(stream.Next()) {
          Int_t ddl=stream.GetDDLNumber();
            for(Int_t row = 1; row <=AliHMPIDRawStream::kNRows; row++){
              for(Int_t dil = 1; dil <=AliHMPIDRawStream::kNDILOGICAdd; dil++){
                for(Int_t pad = 0; pad < AliHMPIDRawStream::kNPadAdd; pad++){
                  pCal->FillPedestal(ddl,row,dil,pad,stream.GetCharge(ddl,row,dil,pad));
                }//pad
              }//dil
            }//row
	} //raw data loop
	
	delete reader;

      }// if CALIBRATION_EVENT

      /* exit when last event received, no need to wait for TERM signal */
      if (eventT==END_OF_RUN) {
	printf("EOR event detected\n");
	break;    
    
      } // events loop   

      free(event);
    }

  }

  /* write report */
  printf("Run #%s, received %d calibration events\n",getenv("DATE_RUN_NUMBER"),iEvtNcal);

  if (!iEvtNcal) {
    printf("No calibration events have been read. Exiting\n");
    return -1;
  }

  /* report progress */
  daqDA_progressReport(90);

  for(Int_t ddl=0; ddl < 14; ddl++) {
    
    /* Calculate pedestal for the given ddl, if there is no ddl go t next */
    if(!pCal->CalcPedestal(ddl,Form("./HmpidPedDdl%02i.txt",ddl))) continue;
    
    /* store the result file on FES */
    status=daqDA_FES_storeFile(Form("./HmpidPedDdl%02i.txt",ddl),Form("HMPID_DA_Pedestals_ddl=%02i",ddl));
    if (status) {
      printf("Failed to export file : %d\n",status);
      return -1;
    }
  }//ddl

  /* report progress */
  daqDA_progressReport(100);

  return status;
}
