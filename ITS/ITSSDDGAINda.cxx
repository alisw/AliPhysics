/*
- Contact: - prino@to.infn.it
- Link: - http://www.to.infn.it/~prino/alice/RawData/run11173.date
- Run Type: - PULSER_RUN
- DA Type: - LDC
- Number of events needed: 100
- Input Files: - SDDbase_step1_ddl*c*_sid*.data
- Output Files: - SDDbase_ddl*c*_sid*.data
- Trigger types used: 
*/

//////////////////////////////////////////////////////////////////////////////
// Detector Algorithm for analysis of SDD test pulse runs.                  //
//                                                                          //
// Produces ASCII and ROOT output files with:                               //
// 1. anode quality bit                                                     //
// 1. Baseline values                                                       //
// 2. Raw noise                                                             //
// 3. Common mode coefficients                                              //
// 4. Noise corrected for common mode                                       //
// 5. Gain                                                                  //
// Files are written to FXS                                                 //
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
#include <TSystem.h>
#include <TROOT.h>
#include <TPluginManager.h>

// AliRoot includes
#include "AliRawReaderDate.h"
#include "AliITSOnlineSDDTP.h"
#include "AliITSRawStreamSDD.h"
/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

  int status = 0;

  // line added to solve IO problems
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");

  /* log start of process */
  printf("ITS SDD TEST-PULSE algorithm program started\n");  


  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  Int_t maxNEvents=15; // maximum number of events to be analyzed
  const Int_t kTotDDL=24;
  const Int_t kModPerDDL=12;
  const Int_t kSides=2;

  AliITSOnlineSDDTP **tpan=new AliITSOnlineSDDTP*[kTotDDL*kModPerDDL*kSides];
  TH2F **histo=new TH2F*[kTotDDL*kModPerDDL*kSides];
  Bool_t isFilled[kTotDDL*kModPerDDL*kSides];
  Char_t hisnam[20];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	tpan[index]=new AliITSOnlineSDDTP(iddl,imod,isid,100.);
	sprintf(hisnam,"h%02dc%02ds%d",iddl,imod,isid);
	histo[index]=new TH2F(hisnam,"",256,-0.5,255.5,256,-0.5,255.5);
	isFilled[index]=0;
      }
    }
  }
  
  /* report progress */
  daqDA_progressReport(10);
  Int_t iev=0,iAnalyzedEv=0;
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

      if(iAnalyzedEv>=maxNEvents) break;
      iev++; 
      
      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      switch (event->eventType){
	
	/* START OF RUN */
      case START_OF_RUN:
	break;
	
	/* END OF RUN */
      case END_OF_RUN:
	break;
	
	//      case PHYSICS_EVENT:  // comment this line for test raw data
	//	break;               // comment this line for test raw data

      case CALIBRATION_EVENT:
	break;  // uncomment this line for test raw data
      case PHYSICS_EVENT: // uncomment this line for test raw data
	printf(" event number = %i \n",iev);
	AliRawReader *rawReader = new AliRawReaderDate((void*)event);

	Int_t evtyp=0;
	while(rawReader->ReadHeader()){
	  const UInt_t *subev = rawReader->GetSubEventAttributes();
	  if(subev[0]==0 && subev[1]==0 && subev[2]==0) evtyp=1; 
	}

	rawReader->Reset();
	for(Int_t iddl=0; iddl<kTotDDL;iddl++){
	  for(Int_t imod=0; imod<kModPerDDL;imod++){
	    for(Int_t isid=0;isid<kSides;isid++){
	      Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	      histo[index]->Reset();
	    }
	  }
	}
	AliITSRawStreamSDD s(rawReader);
	
	while(s.Next()){
	  Int_t iDDL=rawReader->GetDDLID();
	  Int_t iCarlos=s.GetCarlosId();
	  if(iDDL>=0 && iDDL<kTotDDL && s.IsCompletedModule()==kFALSE){ 
	    Int_t index=kSides*(kModPerDDL*iDDL+iCarlos)+s.GetChannel(); 
	    histo[index]->Fill(s.GetCoord2(),s.GetCoord1(),s.GetSignal());
	    isFilled[index]=1;
	  }
	}
	delete rawReader;
	for(Int_t iddl=0; iddl<kTotDDL;iddl++){
	  for(Int_t imod=0; imod<kModPerDDL;imod++){
	    for(Int_t isid=0;isid<kSides;isid++){
	      Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	      if(isFilled[index]) tpan[index]->AddEvent(histo[index]);    
	    }
	  }
	}
	
	/* free resources */
	iAnalyzedEv++;
	free(event);
      }
    }
    
  }
    
  TFile *fh=new TFile("SDDgainHistos.root","RECREATE");
  Char_t filnam[100],command[120];
  for(Int_t iddl=0; iddl<kTotDDL;iddl++){
    for(Int_t imod=0; imod<kModPerDDL;imod++){
      for(Int_t isid=0;isid<kSides;isid++){
	Int_t index=kSides*(kModPerDDL*iddl+imod)+isid;
	if(isFilled[index]){
	  tpan[index]->ValidateAnodes();
	  tpan[index]->WriteToASCII();
	  tpan[index]->WriteToROOT(fh);
	  sprintf(filnam,"SDDbase_ddl%02dc%02d_sid%d.data",iddl,imod,isid);
	  sprintf(command,"tar -rf SDDbase_LDC.tar %s",filnam);
	  gSystem->Exec(command);
	}
      }
    }  
  }
  fh->Close();

  /* write report */
  printf("Run #%s, received %d calibration events\n",getenv("DATE_RUN_NUMBER"),iAnalyzedEv);

  /* report progress */
  daqDA_progressReport(90);



  fh->Close();



  status=daqDA_FES_storeFile("./SDDbase_LDC.tar","SDD_Calib");

  /* report progress */
  daqDA_progressReport(100);



  return status;
}
