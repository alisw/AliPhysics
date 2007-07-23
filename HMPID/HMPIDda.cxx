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
  Double_t SummQ[14][48][11][25], SummQ2[14][48][11][25];
  Int_t   nEntries[14][48][11][25];
  Bool_t  isDDLOn[14];
  for(Int_t ddl=0;ddl<=13;ddl++) {
    isDDLOn[ddl] = kFALSE;
    for(Int_t row=1;row<=24;row++)
      for(Int_t dil=1;dil<=10;dil++)
	for(Int_t adr=0;adr<=47;adr++)
	  {
	    SummQ[ddl][adr][dil][row]=0;
	    SummQ2[ddl][adr][dil][row]=0;
	    nEntries[ddl][adr][dil][row]=0;
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

      if (eventT==PHYSICS_EVENT) {
	
	iEvtNcal++;

	AliRawReader *reader = new AliRawReaderDate((void*)event);

	// Temporary there. Shall be removed as soon as the equipment ID is set correctly
	// For the moment ddl.map file contains one line which maps
	// the observed eqID=225 to the first HMPID DDL with ID=1536
	//	reader->LoadEquipmentIdsMap("ddl.map");
	
	AliHMPIDRawStream stream(reader);

	while(stream.Next()) {

	  Int_t ddl = stream.GetDDLNumber();

	  for(Int_t row=1; row < 25; row++)
	    for(Int_t dil=1; dil < 11; dil++)
	      for(Int_t adr=0; adr < 48; adr++) {
		Short_t q = stream.GetCharge(row,dil,adr);
		if (q<0) continue;

		SummQ[ddl][adr][dil][row] += q;
		SummQ2[ddl][adr][dil][row] += (q*q);
		nEntries[ddl][adr][dil][row]++;
		isDDLOn[ddl] = kTRUE;
	      }
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
    if (!isDDLOn[ddl]) continue;

    ofstream out;
    out.open(Form("./HmpidPedDdl%02i.txt",ddl));

    for(Int_t row=1; row < 25; row++)
      for(Int_t dil=1; dil < 11; dil++)
	for(Int_t adr=0; adr < 48; adr++) {

	  Int_t n = nEntries[ddl][adr][dil][row];
	  if (!n) {
	    printf("No data for channel: %d %d %d %d\n",ddl,adr,dil,row);
	    continue;
	  }
	  Double_t mean = SummQ[ddl][adr][dil][row]/n;
	  if ((SummQ2[ddl][adr][dil][row]/n - (SummQ[ddl][adr][dil][row]/n)*(SummQ[ddl][adr][dil][row]/n)) < 0) {
	    printf("Invalid sums: %d %d %d   %d  %f %f\n",row,dil,adr,n,SummQ[ddl][adr][dil][row],SummQ2[ddl][adr][dil][row]);
	    return -1;
	  }
          Float_t sigma = TMath::Sqrt(SummQ2[ddl][adr][dil][row]/n - (mean*mean));
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
