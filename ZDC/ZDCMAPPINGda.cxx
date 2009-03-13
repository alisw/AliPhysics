/*

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA to write mapping for the ZDC ADCs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: PHYSICS, STANDALONE_BC, STANDALONE_COSMIC, STANDALONE_CENTRAL, 
	  STANDALONE_MB, STANDALONE_SEMICENTRAL
DA Type: MON
Number of events needed: 1 (SOD is read) 
Input Files:  none
Output Files: ZDCChMapping.dat
Trigger Types Used: different trigger types are used

*/

#define MAPDATA_FILE  "ZDCChMapping.dat"

#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>

// DATE
#include <event.h>
#include <monitor.h>
#include <daqDA.h>

//ROOT
#include <TFile.h>

//AliRoot
#include <AliRawReaderDate.h>
#include <AliRawEventHeaderBase.h>
#include <AliZDCRawStream.h>

int main(int argc, char **argv) {

  const Char_t* tableSOD[]  = {"ALL", "no", "SOD", "all", NULL, NULL};
  monitorDeclareTable(const_cast<char**>(tableSOD));
  
  int status = 0;
  int const kNChannels = 48;

  /* log start of process */
  printf("\nZDC MAPPING program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }
  
  FILE *mapFile4Shuttle;

  /* read the data files */
  int n;
  for(n=1;n<argc;n++){
   
    status=monitorSetDataSource( argv[n] );
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
    monitorSetNowait();
    monitorSetNoWaitNetworkTimeout(1000);

    struct eventHeaderStruct *event;
    eventTypeType eventT;

    Int_t iev = 0;
    Bool_t sodRead = kFALSE;
    while(!sodRead && iev<1000){
 
      /* check shutdown condition */
      if (daqDA_checkShutdown()) {break;}

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if(status==MON_ERR_EOF){
        printf ("End of File detected\n");
        break; /* end of monitoring file has been reached */
      }
      if(status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if(event==NULL) continue;
      
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      // --- Reading event header
      //UInt_t evtype = reader->GetType();
      //printf("\t ZDCMAPPINGda -> run # %d\n",reader->GetRunNumber());
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
        
	
      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      
      Int_t ich=0;
      Int_t adcMod[kNChannels], adcCh[kNChannels], sigCode[kNChannels];
      Int_t det[kNChannels], sec[kNChannels];
      
      if(eventT==START_OF_DATA){
	sodRead = kTRUE;
	rawStreamZDC->SetSODReading(kTRUE);
	
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
        else{
	  while(rawStreamZDC->Next()){
            if((rawStreamZDC->IsChMapping()) && (ich<kNChannels)){
	      adcMod[ich] = rawStreamZDC->GetADCModFromMap(ich);
	      adcCh[ich] = rawStreamZDC->GetADCChFromMap(ich);
	      sigCode[ich] = rawStreamZDC->GetADCSignFromMap(ich);
	      det[ich] = rawStreamZDC->GetDetectorFromMap(ich);
	      sec[ich] = rawStreamZDC->GetTowerFromMap(ich);
	      //
	      fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",
	        ich,adcMod[ich],adcCh[ich],sigCode[ich],det[ich],sec[ich]);
	      //
	      //printf("ZDCMAPPINGda.cxx -> %d mod %d ch %d code %d det %d sec %d\n",
	      //   ich,adcMod[ich],adcCh[ich],sigCode[ich],det[ich],sec[ich]);
	      //
	      ich++;
	    }
	  }
	}
        fclose(mapFile4Shuttle);
      } // SOD event
      else{ 
        if(sodRead){
	  printf("\t SOR read -> exiting from DA\n");
	  break;
	}
	else continue;
      }
      
      iev++; 
    }    
      
    /* free resources */
    //free(event);
  }
  
  /* store the result files on FES */
  status = daqDA_FES_storeFile(MAPDATA_FILE, MAPDATA_FILE);
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  return status;
}
