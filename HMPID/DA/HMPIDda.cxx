/*

HMPID DA for online calibration

Contact:                  Levente.Molnar@cern.ch, Giacomo.Volpe@ba.infn.it
Link:                     https://twiki.cern.ch/twiki/bin/view/ALICE/DAInstructions
Run Type:                 PEDESTAL -- but we select on the PHYSICS_EVENTS in the HMPIDda.cxx
DA Type:                  LDC
Reference Run:             78734
Number of events needed:  ~ 1000 events
Input Files:              Raw pedestal file, EXTERNAL config files: HmpDaqDaConfig.txt, HmpDeadChannelMap.txt on all HMPID LDCs
Output Files:             2 x 14 txt files including pedestal and error values
Trigger types used:       CALIBRATION RUN (but selecting on PHYSICS_EVENT)

*/

extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <fstream>
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

//AMORE
#include <AmoreDA.h>


int main(int argc, char **argv){ 

  int status=0;
  const Char_t         *hmpConfigFile        = "HmpDaqDaConfig.txt"; 
  const Char_t         *hmpDeadChannelMapFile = "HmpDeadChannelMap.txt"; 
  const Int_t               ddlOffset = 1536;
        TString                 hmpIn,hmpIn2;
        TString                 feeIn;

  
  /* log start of process */
  printf("HMP PedestalDa: HMPID DA program started\n");  

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo","*","TStreamerInfo","RIO","TStreamerInfo()");
  
  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("HMP PedestalDa: Wrong number of arguments\n");
    return -1;
  }

  /* copy locally a file from daq detector config db */
  
 hmpIn=Form("./%s",hmpConfigFile);
 status=daqDA_DB_getFile(hmpConfigFile,hmpIn.Data()); 
 if (status) { printf("Failed to get HMPID config file status: %d\n",status);return -1; }
 hmpIn2=Form("./%s",hmpDeadChannelMapFile);
 status=daqDA_DB_getFile(hmpDeadChannelMapFile,hmpIn2.Data()); 
 if (status) { printf("Failed to get HMPID dead channel file status: %d\n",status);return -1; }
  

  /* report progress */
  daqDA_progressReport(10);
 
  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  /* set local storage of ped files for Fe2C */

  
  /* init the pedestal calculation */
  AliHMPIDCalib *pCal=new AliHMPIDCalib();
  /* Set the number of sigma cuts inside the file HmpidSigmaCut.txt on all LDCs! */
  /* If the file is NOT present then the default cut 3 will be used!*/
 pCal->SetSigCutFromFile(hmpIn);
 pCal->SetDeadChannelMapFromFile(hmpIn2);  
  
  /* init event counter */
  Int_t firstEvt=0;
  Int_t iEvtNcal=firstEvt;                      //Start from 1 not 0!                                                 
  ULong_t runNum=0;
  ULong_t ldcId=0;

  int n;
  for (n=1;n<argc;n++) {

    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("HMP PedestalDa: monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
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
        printf("HMP PedestalDa: End of monitoring file has been reached! \n");
        break;
        }
      if (status!=0) {
        printf("HMP PedestalDa: monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
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
	printf("HMP PedestalDa: EOR event detected\n");
	break;    
    
      } // events loop   

      free(event);
    }

  }//arg

  /* write report */
  printf("HMP PedestalDa: HMPID DA processed RUN #%s on LDC#%d, with %d calibration events\n",getenv("DATE_RUN_NUMBER"),ldcId,iEvtNcal);

  if (!iEvtNcal) {
    printf("HMP PedestalDa: No calibration events have been read. Exiting\n");
    return -1;
  }

  /* report progress */
  daqDA_progressReport(90);
  /* send pedestal and error files to DAQ FXS */
  for(Int_t nDDL=0; nDDL < AliHMPIDRawStream::kNDDL; nDDL++) {   
    feeIn=Form("thr%d.dat",ddlOffset+nDDL);         
    /* Calculate pedestal for the given ddl, if there is no ddl go t next */
    if(pCal->CalcPedestal(nDDL,Form("HmpidPedDdl%02i.txt",nDDL),Form("%s",feeIn.Data()),iEvtNcal)) {
      status=daqDA_DB_storeFile(feeIn.Data(),feeIn.Data());                                               //store a single threshold file for a DDL in DAQ DB  
      if (status) { printf("HMP PedestalDa: Failed to store file %s in DAQ DB, status: %d\n",feeIn.Data(),status); }
      status=daqDA_FES_storeFile(Form("HmpidPedDdl%02i.txt",nDDL),Form("HmpidPedDdl%02i.txt",nDDL));      //store a single pedestal file for a DDL
      if (status) { printf("HMP PedestalDa: Failed to export file on FES: %d\n",status); }
    }//pedestals and thresholds
    if(pCal->WriteErrors(nDDL,Form("HmpidErrorsDdl%02i.txt",nDDL),iEvtNcal)) {
      status=daqDA_FES_storeFile(Form("HmpidErrorsDdl%02i.txt",nDDL),Form("HmpidErrorsDdl%02i.txt",nDDL));
      if (status) { printf("HMP PedestalDa: Failed to export file : %d\n",status); }
    }//errors
  }//nDDL

  /* send files to AMORE DB */
  daqDA_progressReport(95);
  Int_t statusAmoreDA=0; 
  amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
  for(Int_t iCh=0; iCh < AliHMPIDParam::kMaxCh; iCh++) {  
  statusAmoreDA+=amoreDA.Send(Form("fPedMeanMap%d",iCh), pCal->GetPedMeanMap(iCh));  
  statusAmoreDA+=amoreDA.Send(Form("fPedSigMap%d",iCh),  pCal->GetPedSigMap(iCh));  
  }
  for(Int_t iCh=0;iCh<=AliHMPIDParam::kMaxCh;iCh++)
   {
    for(Int_t iFee=0;iFee<6;iFee++)
     {
         statusAmoreDA+=amoreDA.Send(Form("f1DPedMean_Ch%d_FEE_%d",iCh,iFee),pCal->GetPedMean(6*iCh+iFee));
         statusAmoreDA+=amoreDA.Send(Form("f1DPedSigma_Ch%d_FEE_%d",iCh,iFee),pCal->GetPedSigma(6*iCh+iFee));
       }
     }

  printf("HMP PedestalDa: num. masked pads: %d, num of currently dead pads: %d.\n",pCal->GetNumMaskedPads(),pCal->GetNumDeadPads()); 
           
  delete pCal;  
  if(status) return -1;
  if(statusAmoreDA) return -1;
  
  /* report progress */
  daqDA_progressReport(100);
  
  
  
  return status;
}

