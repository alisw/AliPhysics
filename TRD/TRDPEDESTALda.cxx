/*

- "Contact: r.bailhache@gsi.de"
- "Link: http://www-linux.gsi.de/~bailhach/PEDESTAL/raw.root.date" 
- "Run Type: dedicated pedestals run"
- "DA Type: LDC"
- "Number of events needed: 100"
- "Input Files: no config files, no previous result files, RAW DATA files from DDL = 1024 to DDL = 1041 included"
- "Output Files: trdCalibration.root, trdCalibration.root, no persitent file over runs"
- "Trigger types used: PEDESTAL_RUN"

*/

#define RESULT_FILE "trdCalibration.root"

extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"
#include <stdio.h>
#include <stdlib.h>

//
// Root includes
//
#include <TFile.h>
#include <TStopwatch.h>

//
// AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTRDRawStreamV2.h"
#include "AliCDBManager.h"
//
// AliRoot TRD calib classes
//
#include "AliTRDCalibPadStatus.h"


/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

  int status;


  /* log start of process */
  printf("TRD DA PEDESTAL started\n");  


  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /* copy locally a file from daq detector config db */
  //status=daqDA_DB_getFile("myconfig","./myconfig.txt");
  //if (status) {
  //  printf("Failed to get config file : %d\n",status);
  //  return -1;
  //}
  /* and possibly use it */
  

  /* init some counters */
  int nevents_total=0;
  int nevents      =0;
 
  //Instance of AliCDBManager: needed by AliTRDRawStream
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  // AliTRDCalibPadStatus object
  AliTRDCalibPadStatus calipad = AliTRDCalibPadStatus();
  Bool_t passpadstatus = kTRUE;


  /*see the time*/
  TStopwatch timer;
  timer.Start();


  /* read the data files */
  int n;
  for (n=1;n<argc;n++) {
   
    /* define data source : this is argument i */
    printf("Processing file %s\n",argv[n]);
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* read the file  until EOF */
    for(;;) {
      struct eventHeaderStruct *event;
      
      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) {
	printf("End of File %d detected\n",n);
	break; /* end of monitoring file has been reached */
      }
      if (status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        break;
      }

      /* retry if got no event */
      if (event==NULL) {
        break;
      }

      if(passpadstatus){

	AliRawReader *rawReader = new AliRawReaderDate((void*)event);
	AliTRDRawStreamV2 *trdRawStream = new AliTRDRawStreamV2((AliRawReader *)rawReader);
	if(!calipad.ProcessEvent(trdRawStream,(Bool_t)nevents_total)) passpadstatus = kFALSE;
	nevents++;
	delete trdRawStream;
	delete rawReader;
     
      }

      nevents_total++;

      /* free resources */
      free(event);
    }
  }


  /* report progress */
  printf("%d events processed and %d used\n",nevents_total,nevents);

  /*see the time*/
  timer.Stop();
  timer.Print();

  /* write file in any case to see what happens in case of problems*/
  /*see the time*/
  TStopwatch timer1;
  timer1.Start();
  TFile *fileTRD = new TFile(RESULT_FILE,"recreate");
  calipad.AnalyseHisto();
  calipad.Write("calibpadstatus");
  fileTRD->Close();   
  printf("Wrote local file %s\n",RESULT_FILE);
  /*see the time*/
  timer1.Stop();
  timer1.Print();
  
  /* store the result file on FES */
  status=daqDA_FES_storeFile(RESULT_FILE,RESULT_FILE);
  if (status) {
    printf("Failed to export file : %d\n",status);
    return -1;
  }
  
   
  /* report progress */
  daqDA_progressReport(100);

  
  return status;
}
