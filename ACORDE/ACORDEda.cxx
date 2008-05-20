/*

ACORDE DA for online calibration

Contact: Pedro.Gonzalez.Zamora@cern.ch
Link: missing
Reference Run: missing
Run Type: PHYSICS
DA Type: LDC
Number of events needed: depending on the run, being run-level
Input Files: ACORDEDa.root
Output Files: ACORDEHistos.root to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT

*/

#define FILE_TOTAL "ACORDEHistos.root"
// DATE
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

//AliRoot
#include "AliACORDERawStream.h"
#include "AliRawReaderDate.h"

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>


#include <stdio.h>
#include <stdlib.h>



/* Main routine
      Arguments: list of DATE raw data files
*/
int main(int argc, char **argv) {

  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");

 int status;

 // AliACORDERawStream reader;

  bool fHitMap[60];
  bool fMultMap[60];
  int eventHits=0;
  int multHits=0;

  UInt_t word;
  int kModules=60;
  
  TH1D *fH1=new TH1D("fHist1", "Single MUON Hits", kModules+1, 0, kModules);
  TH1D *fH2=new TH1D("fHist2", "Total  MUON Hits ", 10, 0, 9);
  TH1D *fH3=new TH1D("fHist3", "Hit Multiplicity (Mult words)", kModules+1, 0, kModules);
  TH1D *fH4=new TH1D("fHist4", "Total Hit Multiplicity (Mult words)", 10, 0, 9);
  





  /* log start of process */
  printf("ACORDE DA  program started\n");  


  /* check that we got some arguments = list of files */
  if (argc<2) 
  {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* open log */

  FILE *fp=NULL;
  fp=fopen("./da.log","a");
  if (fp==NULL) {
    printf("Failed to open file\n");
    return -1;
  }

  
  
  /* report progress */
  daqDA_progressReport(10);


  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;


  /* read the data files */
  int n;
  for (n=1;n<argc;n++) {
   
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* report progress */
    daqDA_progressReport(10+80*n/argc);

    /* read the file */
    for(;;) {
      struct eventHeaderStruct *event;
      eventTypeType eventT;

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);

      if (status==MON_ERR_EOF) break; /* end of monitoring file has been reached */

      if (status!=0) 
      {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if (event==NULL) 
      {
        break;
      }

      /* use event - here, just write event id to result file */
      eventT=event->eventType;

      switch (event->eventType)
      {

      case START_OF_RUN:
	break;
      case END_OF_RUN:
	break;
      case PHYSICS_EVENT:
              nevents_physics++;
              fprintf(fp,"Run #%lu, event size: %lu, BC:%u, Orbit:%u, Period:%u\n",
              (unsigned long)event->eventRunNb,
              (unsigned long)event->eventSize,
              EVENT_ID_GET_BUNCH_CROSSING(event->eventId),
              EVENT_ID_GET_ORBIT(event->eventId),
              EVENT_ID_GET_PERIOD(event->eventId));

              AliRawReader *rawReader = new AliRawReaderDate((void*)event);

              AliACORDERawStream* rawStream = new AliACORDERawStream(rawReader);
             
              rawStream->Reset();
          
              rawStream->Next();
              
              word = rawStream->GetWord(0);
              for(int i=0;i<30;i++){
                       fHitMap[i]=word & 1;
                        word>>=1;
              }

              word = rawStream->GetWord(1);
              for(int i=30;i<60;i++){
                fHitMap[i]=word & 1;
                word>>=1;
              }
               
              word = rawStream->GetWord(2);
              for(int i=0;i<30;i++){
                fMultMap[i]=word & 1;
                word>>=1;
              }

              word = rawStream->GetWord(3);
              for(int i=30;i<60;i++){
                fMultMap[i]=word & 1;
                word>>=1;
              }

              for(int i=0; i<kModules; ++i){
           	   if(fHitMap[i]){
              		fH1->Fill(i);
              		++eventHits;
              	   }	
		   if(fMultMap[i]){
             		fH3->Fill(i);
             		++multHits;
                   }		
               }
               fH2->Fill(eventHits);
               fH4->Fill(multHits);

              delete rawStream;
              rawStream = 0x0;
              delete rawReader;
              rawReader = 0x0;
         } //End Swicht 	

               

     

      nevents_total++;

      /* free resources */
      free(event);
    }   
    
  }


//write root file


  TObjArray Hlist(0);
  Hlist.Add(fH1);
  Hlist.Add(fH2);
  Hlist.Add(fH3);
  Hlist.Add(fH4);

  TFile *histoFile = new TFile(FILE_TOTAL,"recreate");
  Hlist.Write();
  histoFile->Close();
  delete histoFile;



  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);


  /* close result file */
  fclose(fp);


  /* report progress */
  daqDA_progressReport(90);


  /* store the result file on FES */
  status=daqDA_FES_storeFile(FILE_TOTAL,"CALIB");
  if (status) {
    printf("Failed to export file : %d\n",status);
    return -1;
  }


  /* report progress */
  daqDA_progressReport(100);

  
  return status;
}
