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
#include "TSystem.h"
//
// AliRoot includes
//
#include "AliRawReader.h"
#include "AliRawReaderDate.h"
//#include "AliTRDrawFastStream.h"
//#include "AliTRDrawStreamBase.h"
#include "AliTRDgeometry.h"
#include "AliCDBManager.h"
#include "AliLog.h"

//
//AMORE
//
#include <AmoreDA.h>

//
// AliRoot TRD calib classes
//
#include "AliTRDCalibPadStatus.h"

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
  unsigned long32 runNb=0;      //run number

  // setting
  //AliTRDrawFastStream::DisableSkipData();
  AliLog::SetGlobalLogLevel(AliLog::kFatal); 

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
	// for debug
	//rawReader->SelectEquipment(-1,1024,1025);
	
	Int_t result = calipad.ProcessEvent3((AliRawReader *) rawReader);
	// 0 error, 1 no input, 2 output
	if(result == 2) nevents++;
	if(result == 0) passpadstatus = kFALSE;
	delete rawReader;
      
      }

      nevents_total++;

      // get the run number
      runNb = event->eventRunNb;

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
void SendToAmoreDB(TObject *calipad, unsigned long32 runNb)
{
  //cheet a little -- temporary solution (hopefully)
  //
  //currently amoreDA uses the environment variable AMORE_DA_NAME to create the mysql
  //table in which the calib objects are stored. This table is dropped each time AmoreDA
  //is initialised. This of course makes a problem if we would like to store different
  //calibration entries in the AMORE DB. Therefore in each DA which writes to the AMORE DB
  //the AMORE_DA_NAME env variable is overwritten.
  
  
  //
  // The reference data are stored in:
  // PadStatus1 for sm-00-01-02
  // PadStatus2 for sm-03-04-05
  // PadStatus3 for sm-06-07-08
  // PadStatus4 for sm-09-10-11
  // PadStatus5 for sm-12-13-14
  // PadStatus6 for sm-15-16-17
  // PadStatus0 if nothing found..means problems
  //

  ///////////////////
  // Find wich LDC
  ///////////////////
  Int_t ldcnumber = -1;
  Int_t sm = -1;
  for (Int_t idet=0; idet<540; idet++) {
    AliTRDCalROC *rocMean  = ((AliTRDCalibPadStatus *) calipad)->GetCalRocMean(idet, kFALSE);
    if ( rocMean )  {
      sm  = AliTRDgeometry::GetSector(idet);  
      if((sm==0) || (sm==1) || (sm==2)) ldcnumber = 1;
      if((sm==3) || (sm==4) || (sm==5)) ldcnumber = 2;
      if((sm==6) || (sm==7) || (sm==8)) ldcnumber = 3;
      if((sm==9) || (sm==10) || (sm==11)) ldcnumber = 4;
      if((sm==12) || (sm==13) || (sm==14)) ldcnumber = 5;
      if((sm==15) || (sm==16) || (sm==17)) ldcnumber = 6;
    }
  }
  const char *amoreDANameorig=gSystem->Getenv("AMORE_DA_NAME");
  
  gSystem->Setenv("AMORE_DA_NAME",Form("TRD-dataQA-%02d-%s",ldcnumber,FILE_ID));
 
  /////////////////////
  // Send the stuff
  ////////////////////
  if (ldcnumber>-1){
    TDatime time;
    TObjString info(Form("Run: %u; Date: %s",runNb,time.AsSQLString()));
    
    amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
    Int_t statusDA=0;
    statusDA+=amoreDA.Send("Pedestals",calipad);
    statusDA+=amoreDA.Send("Info",&info);
    if ( statusDA )
      printf("Warning: Failed to write one of the calib objects to the AMORE database\n");
  }  else {
    printf("Warning: No data found!\n");
  }
  
  // reset env var
  if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
  
}
  


