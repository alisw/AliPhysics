/*

Contact: r.bailhache@gsi.de
Link:
Reference run: 12170
Run Type: PEDESTAL
DA Type: LDC
Number of events needed: 100
Input Files: TRD raw files
Output Files: trdCalibration.root
Trigger types used: PHYSICS_EVENT

*/

#define RESULT_FILE  "trdPedestal.root"
#define FILE_ID "PADSTATUS"

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
#include "TROOT.h"
#include "TPluginManager.h"
//
// AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTRDrawStream.h"
#include "AliTRDrawStreamBase.h"
#include "AliCDBManager.h"
//
// AliRoot TRD calib classes
//
#include "AliTRDCalibPadStatus.h"


/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

  /* magic line from Rene */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");
  
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
  //AliCDBManager *man = AliCDBManager::Instance();
  //man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //man->SetRun(0);
  // AliTRDCalibPadStatus object
  AliTRDCalibPadStatus calipad = AliTRDCalibPadStatus();
  Bool_t passpadstatus = kTRUE;

  // setting
  // AliTRDrawStream::SetNoDebug();
  AliTRDrawStream::SetNoErrorWarning();
  //AliTRDrawStream::SetForceCleanDataOnly();
  AliTRDrawStream::AllowCorruptedData();
  AliTRDrawStream::DisableStackNumberChecker();
  AliTRDrawStream::DisableStackLinkNumberChecker();
  //AliTRDrawStream::SetSkipCDH();
  //AliTRDrawStream::SetExtraWordsFix();
  //AliTRDrawStream::EnableDebugStream();
  //AliTRDrawStream::SetDumpHead(320);
  //AliTRDrawStream::SetDumpHead(80);

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
      eventTypeType eventT;
      
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

      /* use event - here, just write event id to result file */
      eventT=event->eventType;

      if((eventT==PHYSICS_EVENT) && (passpadstatus)){

	AliRawReader *rawReader = new AliRawReaderDate((void*)event);
	rawReader->Select("TRD");
	AliTRDrawStream *trdRawStream = new AliTRDrawStream((AliRawReader *)rawReader);
	Int_t result = calipad.ProcessEvent((AliTRDrawStreamBase *)trdRawStream,0);
	// 0 error, 1 no input, 2 output
	if(result == 2) nevents++;
	if(result == 0) passpadstatus = kFALSE;
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

  /* write file in any case to see what happens in case of problems*/
  TFile *fileTRD = new TFile(RESULT_FILE,"recreate");
  if(nevents < 30) calipad.Destroy();
  calipad.AnalyseHisto();
  calipad.Write("calibpadstatus");
  fileTRD->Close();   
     
  /* store the result file on FES */
  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  
  return status;
}
