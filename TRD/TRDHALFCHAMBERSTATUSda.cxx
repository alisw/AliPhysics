/*

Contact: r.bailhache@gsi.de
Link:
Reference run: 104892
Run Type: PHYSICS
DA Type: MON
Number of events needed: many
Input Files: TRD raw files
Output Files: trdCalibrationhcs.root
Trigger types used: PHYSICS_EVENT

*/

#define RESULT_FILE  "trdHalfChamberStatus.root"
#define FILE_ID "HALFCHAMBERSTATUS"

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
#include "TSystem.h"
//
// AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
#include "AliTRDrawFastStream.h"
#include "AliTRDrawStreamBase.h"
#include "AliLog.h"
//#include "AliCDBManager.h"

//
//AMORE
//
#include <AmoreDA.h>

//
// AliRoot TRD calib classes
//
#include "AliTRDCalibChamberStatus.h"

//functios, implementation below
void SendToAmoreDB(TObject *o, unsigned long32 runNb);

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
  printf("TRD DA HALFCHAMBERSTATUS started\n");  

  if (argc!=2) {
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

  /* init some counters */
  int nevents_total=0;
  int nevents      =0;
 
  //Instance of AliCDBManager: needed by AliTRDRawStream
  //AliCDBManager *man = AliCDBManager::Instance();
  //man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //man->SetRun(0);
  // AliTRDCalibPadStatus object
  AliTRDCalibChamberStatus calipad = AliTRDCalibChamberStatus();
  calipad.Init();
  unsigned long32 runNb=0;      //run number

  // setting
  AliTRDrawFastStream::DisableSkipData();
  AliLog::SetGlobalLogLevel(AliLog::kFatal); 

  /* read the file  until EOF */
  for(;;) {
    struct eventHeaderStruct *event;
    eventTypeType eventT;
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event */
    status=monitorGetEventDynamic((void **)&event);
    if (status==MON_ERR_EOF) {
      printf("End of File detected\n");
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
    
    if(eventT==PHYSICS_EVENT){
      
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      rawReader->Select("TRD");
      
      calipad.ProcessEvent((AliRawReader *) rawReader,nevents);
      nevents++;
      delete rawReader;
      
    }
    
    nevents_total++;
    
    // get the run number
    runNb = event->eventRunNb;
    
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }
  

  /* report progress */
  printf("%d events processed and %d used\n",nevents_total,nevents);
  
  /* write file in any case to see what happens in case of problems*/
  TFile *fileTRD = new TFile(RESULT_FILE,"recreate");
  calipad.AnalyseHisto();
  calipad.Write("calibchamberstatus");
  fileTRD->Close();   

  /* Send to amore */
  SendToAmoreDB(&calipad,runNb);
     
  /* store the result file on FES */
  status=daqDA_FES_storeFile(RESULT_FILE,FILE_ID);
  if (status) {
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  
  return status;

}
void SendToAmoreDB(TObject *o, unsigned long32 runNb)
{
  //AMORE
  const char *amoreDANameorig=gSystem->Getenv("AMORE_DA_NAME");
  //cheet a little -- temporary solution (hopefully)
  //
  //currently amoreDA uses the environment variable AMORE_DA_NAME to create the mysql
  //table in which the calib objects are stored. This table is dropped each time AmoreDA
  //is initialised. This of course makes a problem if we would like to store different
  //calibration entries in the AMORE DB. Therefore in each DA which writes to the AMORE DB
  //the AMORE_DA_NAME env variable is overwritten.
  
  gSystem->Setenv("AMORE_DA_NAME","TRD-dataHALFCHAMBERSTATUS");
  //
  // end cheet
  TDatime time;
  TObjString info(Form("Run: %u; Date: %s",runNb,time.AsSQLString()));
  
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  Int_t statusDA=0;
  statusDA+=amoreDA.Send("DataHalfChamberStatus",o);
  statusDA+=amoreDA.Send("Info",&info);
  if ( statusDA!=0 )
    printf("TRDQAda: Waring: Failed to write one of the calib objects to the AMORE database\n");
  // reset env var
  if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
}