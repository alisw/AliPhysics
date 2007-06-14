//////////////////////////////////////////////////////////////////////////////
// Detector Algorithm for analysis of SDD baseline runs.                    //
//                                                                          //
// Produces ASCII and ROOT output files with:                               //
// 1. anode quality bit                                                     //
// 1. Baseline values                                                       //
// 2. Raw noise                                                             //
// 3. Common mode coefficients                                              //
// 4. Noise corrected for common mode                                       //
// Files are stored locally.                                                //
// The next DAQ-DA step on Test Pulse run writes files to FXS               //
//                                                                          //
// Author: F. Prino (prino@to.infn.it)                                      //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////



extern "C" {
#include "daqDA.h"
}

#include "event.h"
#include "monitor.h"

#include <stdio.h>
#include <stdlib.h>

// ROOT includes
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

// AliRoot includes
#include "AliRawReaderDate.h"
#include "AliITSOnlineSDDBase.h"
#include "AliITSOnlineSDDCMN.h"
#include "AliITSRawStreamSDD.h"
/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

  int status = 0;


  /* log start of process */
  printf("ITS SDD TP algorithm program started\n");  


  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  // Loop over modules has to be added
  const Int_t nSDDmodules=12;
  AliITSOnlineSDDBase **base=new AliITSOnlineSDDBase*[2*nSDDmodules];
  AliITSOnlineSDDCMN **corr=new AliITSOnlineSDDCMN*[2*nSDDmodules];
  TH2F **histo=new TH2F*[2*nSDDmodules];

  Char_t hisnam[20];
  for(Int_t imod=0; imod<nSDDmodules;imod++){
    for(Int_t isid=0;isid<2;isid++){
      Int_t index=2*imod+isid;
      base[index]=new AliITSOnlineSDDBase(imod,isid);
      corr[index]=new AliITSOnlineSDDCMN(imod,isid);
      sprintf(hisnam,"his%03ds%d",imod,isid);
      histo[index]=new TH2F(hisnam,"",256,-0.5,255.5,256,-0.5,255.5);
    }
  }
  
  /* report progress */
  daqDA_progressReport(10);
  Int_t iev;
  for(Int_t iStep=0;iStep<2;iStep++){
    /* init some counters */
    printf("Start Analysis Step %d\n",iStep);
    iev=0;
    if(iStep==1){
      for(Int_t imod=0; imod<nSDDmodules;imod++){
	for(Int_t isid=0;isid<2;isid++){
	  Int_t index=2*imod+isid;
	  corr[index]->Reset();
	}
      }
    }

    /* read the data files */
    int n;
    for (n=1;n<argc;n++) {
   
      status=monitorSetDataSource( argv[n] );
      if (status!=0) {
	printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
	return -1;
      }

      /* report progress */
      /* in this example, indexed on the number of files */
      // Progress report inside the event loop as well?
      daqDA_progressReport(10+80*n/argc);

      /* read the file */
      for(;;) {
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

	iev++; 
	if(iev>10) break;

	/* use event - here, just write event id to result file */
	eventT=event->eventType;
	switch (event->eventType){
      
	  /* START OF RUN */
	case START_OF_RUN:
	  break;
      
	  /* END OF RUN */
	case END_OF_RUN:
	  break;

	  //      case PHYSICS_EVENT:
	  //	break;

	case CALIBRATION_EVENT:
	  // for test raw data
	case PHYSICS_EVENT:
	  printf(" event number = %i \n",iev);
	  AliRawReader *rawReader = new AliRawReaderDate((void*)event);
	  rawReader->RequireHeader(kFALSE);

	  // temp for test raw data
	  rawReader->SelectEquipment(17,101,101);

	  Int_t evtyp=0;
	  while(rawReader->ReadHeader()){
	    const UInt_t *subev = rawReader->GetSubEventAttributes();
	    if(subev[0]==0 && subev[1]==0 && subev[2]==0) evtyp=1; 
	  }

	  rawReader->Reset();
	  for(Int_t imod=0; imod<nSDDmodules;imod++){
	    for(Int_t isid=0; isid<2;isid++){
	      Int_t index=2*imod+isid;
	      histo[index]->Reset();
	    }
	  }
	  AliITSRawStreamSDD s(rawReader);
	  // temp for test raw data
	  rawReader->SelectEquipment(17,101,101);

	  while(s.Next()){
	    // calculation of module etc.
	    if(s.GetCarlosId()<nSDDmodules){
	      Int_t index=2*s.GetCarlosId()+s.GetChannel();
	      histo[index]->Fill(s.GetCoord2(),s.GetCoord1(),s.GetSignal());
	    }
	  }
	  delete rawReader;
	  for(Int_t imod=0; imod<nSDDmodules;imod++){
	    for(Int_t isid=0; isid<2;isid++){
	      Int_t index=2*imod+isid;
	      if(iStep==0) base[index]->AddEvent(histo[index]);    
	      if(iStep==1) corr[index]->AddEvent(histo[index]);    
	    }
	  }

	  /* free resources */
	  free(event);
	}
      }
    
    }
    
    for(Int_t imod=0; imod<nSDDmodules;imod++){
      for(Int_t isid=0; isid<2;isid++){
	Int_t index=2*imod+isid;
	if(iStep==0){
	  base[index]->ValidateAnodes();
	  base[index]->WriteToASCII();
	}
	if(iStep==1){
	  corr[index]->ValidateAnodes();
	  corr[index]->WriteToASCII();
	}
      }
    }  
  }

  /* write report */
  printf("Run #%s, received %d calibration events\n",getenv("DATE_RUN_NUMBER"),iev);

  /* report progress */
  daqDA_progressReport(90);



  TFile *fh=new TFile("histos.root","RECREATE");
  for(Int_t imod=0; imod<nSDDmodules;imod++){
    for(Int_t isid=0; isid<2;isid++){
      Int_t index=2*imod+isid;
      corr[index]->WriteToROOT(fh);
    }
  }
  fh->Close();



  // Example how to store the output file ot DAQ FXS
  //  status=daqDA_FES_storeFile("./result.txt","DAcase1_results");

  /* report progress */
  daqDA_progressReport(100);



  return status;
}
