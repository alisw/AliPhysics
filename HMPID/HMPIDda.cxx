/*

HMPID DA for online calibration

Contact: Levente.Molnar@ba.infn.it, Giacomo.Volpe@ba.infn.it
Link: http://richpc1.ba.infn.it/~levente/Files4Public/ValidateHmpidDA/
Run Type: PEDESTAL -- but we select on the PHYSICS_EVENTS in th HMPIDda.cxx
DA Type: LDC
Number of events needed: 1000 events
Input Files: Raw pedestal file, EXTERNAL config file: HmpidSigmaCut.txt on both HMPID LDCs
Output Files: 14 txt files including pedestal values
Trigger types used: PEDESTAL RUN (selecting on PHYSICS_EVENT)

*/

extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include "fstream.h"
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliHMPIDRawStream.h"
#include "AliHMPIDCalib.h"
#include "AliRawReaderDate.h"
#include "AliBitPacking.h"

//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TObject.h"
#include "TPluginManager.h"


int main(int argc, char **argv){ 

  int status;
  const Char_t         *hmpConfigFile = "HmpDaqDaConfig.txt"; 
  const Int_t               ddlOffset = 1536;
        TString                 hmpIn;
        TString                 feeIn;

  
  /* log start of process */
  printf("HMPID DA program started\n");  

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo","*","TStreamerInfo","RIO","TStreamerInfo()");
  
  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* copy locally a file from daq detector config db */
  hmpIn=Form("./%s",hmpConfigFile);
  status=daqDA_DB_getFile(hmpConfigFile,hmpIn.Data());
  if (status) { printf("Failed to get HMPID config file status: %d\n",status); return -1; }

  /* report progress */
  daqDA_progressReport(10);

  
  
  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  /* set local storage of ped files for Fe2C */
 
  
  /* init the pedestal calculation */
  AliHMPIDCalib *pCal=new AliHMPIDCalib();
  /* Set the number of sigma cuts inside the file HmpidSigmaCut.txt on BOTH LDCs! */
  /* If the file is NOT present then the default cut 3 will be used!*/
  pCal->SetSigCutFromFile(hmpIn);
  
  /* ONLY set this option to kTRUE if you want to create the ADC dsitributions for all 161280 pads!!!!*/  
  /* kTRUE is not suggested for production mode b/c of the memory consumption! */
  pCal->SetWriteHistoPads(kFALSE);               //use this option for default production useage!!!
  //pCal->SetWriteHistoPads(kTRUE);              //only for expert debug
  //pCal->SetWriteHistoPads(kTRUE,kTRUE,13);     //DO NOT USE THIS OPTION!
  
  
  /* init event counter */
  Int_t firstEvt=0;
  Int_t iEvtNcal=firstEvt;                      //Start from 1 not 0!                                                 
  ULong_t runNum=0;
  ULong_t ldcId=0;

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
      
       /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
      struct eventHeaderStruct *event;
      eventTypeType eventT;

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF)                                              /* end of monitoring file has been reached */
      {
        printf("End of monitoring file has been reached! \n");
        break;
        }
      
       
      if (status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
	return -1;
      }

      /* retry if got no event */
      if (event==NULL) {
        //break;
        continue;
      }

      /* use event - here, just write event id to result file */
      eventT=event->eventType;

      if (eventT==PHYSICS_EVENT || eventT==CALIBRATION_EVENT  ) {                 //updated: 18/02/2008 based on http://alice-ecs.web.cern.ch/alice-ecs/runtypes_3.16.html
        runNum=(unsigned long)event->eventRunNb;                                  //assuming that only one run is processed at a time
        ldcId=(unsigned long)event->eventLdcId;
        if(iEvtNcal==firstEvt && pCal->GetWritePads()==kTRUE) 
        {
          if(pCal->GetLargePads()==kFALSE)  pCal->InitFile((Int_t)ldcId);         //The number for iEvtNcal should be the same as for the first value
          else pCal->InitFile((Int_t)runNum);                                     //The number for iEvtNcal should be the same as for the first value
        }
        
        iEvtNcal++;        
	AliRawReader *reader = new AliRawReaderDate((void*)event);
	AliHMPIDRawStream stream(reader);stream.SetTurbo(kTRUE);                  //raw data decoding without error checks SetTurbo(kTRUE)
        while(stream.Next())
          {
             for(Int_t iPad=0;iPad<stream.GetNPads();iPad++) {
             pCal->FillPedestal(stream.GetPadArray()[iPad],stream.GetChargeArray()[iPad]);
              } //pads
          }//while -- loop on det load in one event
          
         for(Int_t iddl=0;iddl<stream.GetNDDL();iddl++){   
                 pCal->FillDDLCnt(iddl,stream.GetnDDLInStream()[iddl],stream.GetnDDLOutStream()[iddl]);                     
           for(Int_t ierr=0; ierr < stream.GetNErrors(); ierr++) {
               pCal->FillErrors(iddl,ierr,stream.GetErrors(iddl,ierr));
               }
          }//err   
          
        pCal->SetRunParams(runNum,stream.GetTimeStamp(),stream.GetLDCNumber());   //Get the last TimeStam read and the LDC ID
        stream.Delete();            
      }// if CALIBRATION_EVENT

      /* exit when last event received, no need to wait for TERM signal */
      if (eventT==END_OF_RUN) {
	printf("EOR event detected\n");
	break;    
    
      } // events loop   

      free(event);
    }

  }//arg

  /* write report */
  printf("HMPID DA processed RUN #%s on LDC#%d, with %d calibration events\n",getenv("DATE_RUN_NUMBER"),ldcId,iEvtNcal);

  if (!iEvtNcal) {
    printf("No calibration events have been read. Exiting\n");
    return -1;
  }

  /* report progress */
  daqDA_progressReport(90);
  
  for(Int_t nDDL=0; nDDL < AliHMPIDRawStream::kNDDL; nDDL++) {   
    feeIn=Form("thr%d.dat",ddlOffset+nDDL);         
    /* Calculate pedestal for the given ddl, if there is no ddl go t next */
    if(pCal->CalcPedestal(nDDL,Form("HmpidPedDdl%02i.txt",nDDL),Form("%s",feeIn.Data()),iEvtNcal)) {
      status=daqDA_DB_storeFile(feeIn.Data(),feeIn.Data());                                               //store a single threshold file for a DDL in DAQ DB  
      if (status) { printf("Failed to store file %s in DAQ DB, status: %d\n",feeIn.Data(),status); }
      status=daqDA_FES_storeFile(Form("HmpidPedDdl%02i.txt",nDDL),Form("HmpidPedDdl%02i.txt",nDDL));      //store a single pedestal file for a DDL
      if (status) { printf("Failed to export file on FES: %d\n",status); }
    }//pedestals and thresholds
    if(pCal->WriteErrors(nDDL,Form("HmpidErrorsDdl%02i.txt",nDDL),iEvtNcal)) {
      status=daqDA_FES_storeFile(Form("HmpidErrorsDdl%02i.txt",nDDL),Form("HmpidErrorsDdl%02i.txt",nDDL));
      if (status) { printf("Failed to export file : %d\n",status); }
    }//errors
    /* to create pedestal file as Paolo uncomment the line */
    //if(!pCal->CalcPedestalPaolo(nDDL,Form("%sHmpidPedDdl%02i.txt",sDaOut.Data(),nDDL),iEvtNcal)) continue; 
  }//nDDL

  if(pCal->GetWritePads()==kTRUE) pCal->CloseFile();
  delete pCal;  
  if (status) return -1;
  
  /* report progress */
  daqDA_progressReport(100);
  
  
  return status;
}

