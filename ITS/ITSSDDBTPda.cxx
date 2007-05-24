/*

DAcase1.c

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

contact: alice-datesupport@cern.ch

*/

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
#include "AliITSOnlineSDDTP.h"
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
  AliITSOnlineSDDTP *left=new AliITSOnlineSDDTP(12,1,411.);
  TH2F* histo = new TH2F("histo","",256,-0.5,255.5,256,-0.5,255.5);
  
  /* report progress */
  daqDA_progressReport(10);


  /* init some counters */
  Int_t iev=0;


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

      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      switch (event->eventType){
      
	/* START OF RUN */
      case START_OF_RUN:
	break;
      
	/* END OF RUN */
      case END_OF_RUN:
	break;

      case PHYSICS_EVENT:
	break;

      case CALIBRATION_EVENT:
	// for test raw data
	//      case PHYSICS_EVENT:
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
	histo->Reset();
	AliITSRawStreamSDD s(rawReader);
	// temp for test raw data
	rawReader->SelectEquipment(17,101,101);

	while(s.Next()){
	  // calculation of module etc.
	  if(s.GetCarlosId()==1&&s.GetChannel()==0){
	    histo->Fill(s.GetCoord2(),s.GetCoord1(),s.GetSignal());
	  }
	}
	delete rawReader;
	left->AddEvent(histo);    
      }

      /* free resources */
      free(event);
    }
    
  }


  /* write report */
  printf("Run #%s, received %d calibration events\n",getenv("DATE_RUN_NUMBER"),iev);

  /* report progress */
  daqDA_progressReport(90);


  TH1F *hval=new TH1F("hval","",256,-0.5,255.5);
  TH1F *hgain=new TH1F("hgain","",256,-0.5,255.5);

  for(Int_t ian=0;ian<256;ian++){
    hgain->SetBinContent(ian+1,left->GetChannelGain(ian));
    hval->SetBinContent(ian+1,float(left->IsAnodeGood(ian)));
    printf("Anode: %d, valid=%d, gain=%f\n",ian,left->IsAnodeGood(ian),left->GetChannelGain(ian));
  }

  hgain->SetMaximum(10);

  hgain->GetXaxis()->SetTitle("Anode number");
  hgain->GetYaxis()->SetTitle("Gain (ADC counts/DAC units)");

  left->ValidateAnodes();
  printf("Anodes validated");
  left->WriteToFXS();

  // Example how to store the output file ot DAQ FXS
  //  status=daqDA_FES_storeFile("./result.txt","DAcase1_results");

  /* report progress */
  daqDA_progressReport(100);

  // temp for debug purposes
  TFile fh("histos.root","RECREATE");
  histo->Write();
  hgain->Write();
  hval->Write();
  fh.Close();

  return status;
}
