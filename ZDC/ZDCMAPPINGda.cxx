/*

This program reads the DAQ data files passed as argument using the monitoring library.

It computes the average event size and populates local "./result.txt" file with the 
result.

The program reports about its processing progress.

Messages on stdout are exported to DAQ log system.

DA for ZDC standalone pedestal runs

Contact: Chiara.Oppedisano@to.infn.it
Link: 
Run Type: PHYSICS, STANDALONE_BC, STANDALONE_COSMIC, STANDALONE_CENTRAL, 
	  STANDALONE_MB, STANDALONE_SEMICENTRAL
DA Type: LDC
Number of events needed: no constraint 
Input Files:  
Output Files: ZDCChMapping.dat
Trigger Types Used: 

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


/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {
  
  int status = 0;

  /* log start of process */
  printf("\nZDC MAPPING program started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* open result file */
  FILE *fp=NULL;
  fp=fopen("./result.txt","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }
  
  FILE *mapFile4Shuttle;
  
  /* report progress */
  daqDA_progressReport(10);


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  /* read the data files */
  int n;
  for(n=1;n<argc;n++){
   
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* report progress */
    /* in this example, indexed on the number of files */
    daqDA_progressReport(10+80*n/argc);

    /* read the file */
    for(;;) {
      struct eventHeaderStruct *event;
      eventTypeType eventT;

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if(status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
      if(status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if(event==NULL) {
        break;
      }
      
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      reader->Select("ZDC");
      // --- Reading event header
      //UInt_t evtype = reader->GetType();
      //printf("\n\t ZDCPEDESTALda -> ev. type %d\n",evtype);
      //printf("\t ZDCPEDESTALda -> run # %d\n",reader->GetRunNumber());
      //
      AliZDCRawStream *rawStreamZDC = new AliZDCRawStream(reader);
        

      /* use event - here, just write event id to result file */
      eventT=event->eventType;
      
      Int_t ich=0, adcMod[48], adcCh[48], sigCode[48], det[48], sec[48];
      if(eventT==START_OF_DATA){
	  	
	if(!rawStreamZDC->Next()) printf(" \t No raw data found!! \n");
        else{
	  while(rawStreamZDC->Next()){
            if(rawStreamZDC->IsChMapping()){
	      adcMod[ich] = rawStreamZDC->GetADCModFromMap(ich);
	      adcCh[ich] = rawStreamZDC->GetADCChFromMap(ich);
	      sigCode[ich] = rawStreamZDC->GetADCSignFromMap(ich);
	      det[ich] = rawStreamZDC->GetDetectorFromMap(ich);
	      sec[ich] = rawStreamZDC->GetTowerFromMap(ich);
	      ich++;
	    }
	  }
	}
	// --------------------------------------------------------
	// --- Writing ascii data file for the Shuttle preprocessor
        mapFile4Shuttle = fopen(MAPDATA_FILE,"w");
        for(Int_t i=0; i<ich; i++){
	   fprintf(mapFile4Shuttle,"\t%d\t%d\t%d\t%d\t%d\t%d\n",i,
	     adcMod[i],adcCh[i],sigCode[i],det[i],sec[i]);
	   //
	   //printf("ZDCPEDESTALDA.cxx ->  ch.%d mod %d, ch %d, code %d det %d, sec %d\n",
	   //	   i,adcMod[i],adcCh[i],sigCode[i],det[i],sec[i]);
        }
        fclose(mapFile4Shuttle);
      }

      nevents_total++;

      /* free resources */
      free(event);
    
    }
  }  
  
  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);
  
  /* report progress */
  daqDA_progressReport(90);

  /* store the result files on FES */
  status = daqDA_FES_storeFile(MAPDATA_FILE, MAPDATA_FILE);
  if(status){
    printf("Failed to export file : %d\n",status);
    return -1;
  }

  /* report progress */
  daqDA_progressReport(100);

  return status;
}
