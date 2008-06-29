/*
PMD DA for online calibration

contact: basanta@phy.iitb.ac.in
Link:http://www.veccal.ernet.in/~pmd/
Reference run:
Run Type: PHYSICS
DA Type: MON
Number of events needed: 1 million for PB+PB, 200 milion for p+p
Input Files: PMD_PED.root, Configfile
Output Files: PMDGAINS.root, to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT

*/
extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"
//#include "daqDA.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliRawReaderDate.h"
#include "AliPMDCalibPedestal.h"
#include "AliPMDCalibGain.h"

//ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TTree.h"
#include "TROOT.h"
#include "TPluginManager.h"



/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  
    /* magic line from Rene */
    gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					  "*",
					  "TStreamerInfo",
					  "RIO",
					  "TStreamerInfo()");

    Int_t filestatus = -1, totevt = -1;
    Int_t maxevt = -1;

    // Reads the pedestal file and keep the values in memory for subtraction

    AliPMDCalibGain calibgain;
    Int_t pstatus = calibgain.ExtractPedestal();

    if(pstatus == -3) return -3;

    TTree *ic = NULL;

    FILE *fp1 = NULL;

    fp1 = fopen("Configfile","r");

    if (fp1 == NULL)
      {
	printf("**** Configfile doesn't exist, creating the file ****\n");
	fp1 = fopen("Configfile","w");
	filestatus = 0;
	totevt     = 0;
	maxevt     = 2000;
	fprintf(fp1,"%d %d %d\n",filestatus, totevt,maxevt);
      }
    else
      {
	fscanf(fp1,"%d %d %d\n",&filestatus, &totevt,&maxevt);
	//printf("%d %d %d\n",filestatus, totevt, maxevt);
      }
    fclose(fp1);
    
    if (filestatus == 1)
      {
	calibgain.ReadIntermediateFile();
      }

    // decoding the events
    
    int status;

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
    
    /* log start of process */
    printf("PMD GAIN DA - strted generating the gain of a cell\n");  
    
    /* init some counters */
    int nevents_physics=0;
    int nevents_total=0;
    
    struct eventHeaderStruct *event;
    eventTypeType eventT = 0;

    Int_t iev=0;
    
    /* main loop (infinite) */
    for(;;) {
	
	/* check shutdown condition */
	if (daqDA_checkShutdown()) {break;}
	
	/* get next event (blocking call until timeout) */
	status=monitorGetEventDynamic((void **)&event);
	if (status==MON_ERR_EOF) {
	    printf ("End of File detected\n");
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
	
	iev++; 
	
	/* use event - here, just write event id to result file */
	nevents_total++;
	eventT=event->eventType;
	switch (event->eventType){
      
	    /* START OF RUN */
	    case START_OF_RUN:
		break;
		/* END START OF RUN */
		
		/* END OF RUN */
	    case END_OF_RUN:
		break;
		
	    case PHYSICS_EVENT:
		nevents_physics++;
		//if(nevents_physics%100 == 0)printf("Physis Events = %d\n",nevents_physics);
		AliRawReader *rawReader = new AliRawReaderDate((void*)event);
		TObjArray *pmdddlcont = new TObjArray();
		calibgain.ProcessEvent(rawReader, pmdddlcont);

		delete pmdddlcont;
		pmdddlcont = 0x0;
		delete rawReader;
		rawReader = 0x0;
		
	}
       
	/* free resources */
	free(event);
	
    }

    /* exit when last event received, no need to wait for TERM signal */

    ic = new TTree("ic","PMD Gain tree");

    totevt += nevents_physics++;

    fp1 = fopen("Configfile","w+");

    if (totevt < maxevt)
      {
	printf("Required Number of Events not reached\n");
	printf("Number of Events processed = %d\n",totevt);
	printf("Writing the intermediate ASCII file\n");
	calibgain.WriteIntermediateFile();

	filestatus = 1;
	fprintf(fp1,"%d %d %d\n",filestatus,totevt,maxevt);
      }
    else if (totevt >= maxevt)
      {
	printf("Required Number of Events reached = %d\n",totevt);
	calibgain.Analyse(ic);

	TFile * gainRun = new TFile ("PMDGAINS.root","RECREATE"); 
	ic->Write();
	gainRun->Close();

	filestatus = 0;
	totevt     = 0;
	fprintf(fp1,"%d %d %d\n",filestatus,totevt,maxevt);
      }
    fclose(fp1);
    
    delete ic;
    ic = 0;
    

    /* store the result file on FES */
 
    if (filestatus == 0)
      {
	printf("root file is created and getting exported\n");
	status = daqDA_FES_storeFile("PMDGAINS.root","gaincalib");
      }

    if (status) {
      status = -2;
    }


    return status;
}
