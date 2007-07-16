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
//#include "AliHMPIDRawStream.h"
#include "AliRawReaderDate.h"
#include "AliBitPacking.h"
#include "TMath.h"

//ROOT
#include "TFile.h"
#include "TKey.h"
#include "TH2S.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom.h"

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


  /* init some counters */
  Float_t SummQ[14][48][11][25], SummQ2[14][48][11][25];
  Bool_t  isDDLOn[14];
  for(Int_t ddl=0;ddl<=13;ddl++) {
    isDDLOn[ddl] = kFALSE;
    for(Int_t row=1;row<=24;row++)
      for(Int_t dil=1;dil<=10;dil++)
	for(Int_t adr=0;adr<=47;adr++)
	  {
	    SummQ[ddl][adr][dil][row]=0;
	    SummQ2[ddl][adr][dil][row]=0;
	  }
  }

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

      if (eventT==CALIBRATION_EVENT) {
	
	iEvtNcal++;

	AliRawReader *pRR = new AliRawReaderDate((void*)event);
	
	pRR->Select("HMPID");//select only one DDL files 
      
	UInt_t RawDataWord=0;
      
	while(pRR->ReadNextInt(RawDataWord)) //raw records loop (in selected DDL files)
	  {	  	    
	    Int_t ddl = pRR->GetDDLID();       

	    Int_t a = AliBitPacking::UnpackWord(RawDataWord,12,17); assert(0<=a&&a<=47);   // 1098 7654 3210 9876 5432 1098 7654 3210 DILOGIC address (0..47)
	    Int_t d = AliBitPacking::UnpackWord(RawDataWord,18,21); assert(1<=d&&d<=10);   // 3322 2222 2222 1111 1111 1000 0000 0000 DILOGIC number  (1..10)      
	    Int_t r = AliBitPacking::UnpackWord(RawDataWord,22,26); assert(1<=r&&r<=24);   // Row number (1..24)
	  
	    Int_t q = AliBitPacking::UnpackWord(RawDataWord, 0,11); assert(0<=q&&q<=4095); // 0000 0rrr rrdd ddaa aaaa qqqq qqqq qqqq Qdc        (0..4095)

	    SummQ[ddl][a][d][r]+=q;
            SummQ2[ddl][a][d][r]+=(q*q);
	    isDDLOn[ddl] = kTRUE;
	  }//raw records loop
	
	delete pRR;

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
    if (!isDDLOn[ddl]) continue;

    ofstream out;
    out.open(Form("./HmpidPedDdl%02i.txt",ddl));

    for(Int_t row=1; row < 25; row++)
      for(Int_t dil=1; dil < 11; dil++)
	for(Int_t adr=0; adr < 48; adr++) {

	  Float_t mean = SummQ[ddl][adr][dil][row]/iEvtNcal;
          Float_t sigma = TMath::Sqrt(SummQ2[ddl][adr][dil][row]/iEvtNcal - (SummQ[ddl][adr][dil][row]/iEvtNcal)*(SummQ[ddl][adr][dil][row]/iEvtNcal));
	  Int_t inhard=((Int_t(mean))<<9)+Int_t(mean+3*sigma);

	  out << Form("%2i %2i %2i %5.2f %5.2f %x\n",row,dil,adr,mean,sigma,inhard);
	}

    /* store the result file on FES */
    status=daqDA_FES_storeFile(Form("./HmpidPedDdl%02i.txt",ddl),Form("HMPID_DA_Pedestals_ddl=%02i",ddl));
    if (status) {
      printf("Failed to export file : %d\n",status);
      return -1;
    }
  }


  /* report progress */
  daqDA_progressReport(100);

  return status;
}
