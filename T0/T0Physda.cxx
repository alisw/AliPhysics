/*
T0 DA for online calibration
 
Contact: Michal.Oledzki@cern.ch
Link: http://users.jyu.fi/~mioledzk/
Run Type: PHYSICS
DA Type: MON
Number of events needed: 500000 
Input Files: inPhys.dat, external parameters
Output Files: daPhys.root, to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT
------------------------------Alla
Now trigger type changed to CALICRATION_EVENT
to have data during test.
SOULD BE CHANGED BACK BEFORE BEAM
------------------------------- Alla
*/

#define FILE_OUT "daPhys.root"
#define FILE_IN "inPhys.dat"
#include <daqDA.h>
#include <event.h>
#include <monitor.h>
 
#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include <AliRawReaderDate.h>
#include <AliRawReader.h>
#include <AliT0RawReader.h>

//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TString.h"
#include "TH1.h"

int cbx, ccbx, t0bx, npmtA, npmtC;
float clx,cmx,cclx,ccmx, t0lx, t0hx;

/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
//int main(){
  int status;

  /* magic line */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
                                        "*",
                                        "TStreamerInfo",
                                        "RIO",
                                        "TStreamerInfo()");
  
  if(daqDA_DB_getFile(FILE_IN, FILE_IN)){
    printf("Couldn't get input file >>inPhys.dat<< from DAQ_DB !!!\n");
    return -1;
  }
  
  
  FILE *inp;
  char c;
  inp = fopen(FILE_IN, "r");
  if(!inp){
    printf("Input file >>inPhys.dat<< not found !!!\n");
    return -1;
  }
  
  while((c=getc(inp))!=EOF) {
    switch(c) {
    case 'a': {fscanf(inp, "%d", &ccbx ); break;} //N of X bins hCFD1_CFD
    case 'b': {fscanf(inp, "%f", &cclx ); break;} //Low x hCFD1_CFD
    case 'c': {fscanf(inp, "%f", &ccmx ); break;} //High x hCFD1_CF
    case 'd': {fscanf(inp, "%d", &npmtC ); break;} //number of reference PMTC
    case 'e': {fscanf(inp, "%d", &npmtA ); break;} //number of reference PMTA
    case 'f': {fscanf(inp, "%d", &t0bx ); break;} //N of X bins hT0
    case 'g': {fscanf(inp, "%f", &t0lx ); break;} //Low x hT0
    case 'k': {fscanf(inp, "%f", &t0hx ); break;} //High x hT0
    }
  }
  fclose(inp);

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
  printf("T0 monitoring program started\n");  
  
  // Allocation of histograms - start

  TH1F *hCFD1minCFD[24];  
   
  for(Int_t ic=0; ic<24; ic++) {
    hCFD1minCFD[ic] = new TH1F(Form("CFD1minCFD%d",ic+1),"CFD-CFD",ccbx,cclx,ccmx);
  }
  TH1F *hVertex = new TH1F("hVertex","T0 time",t0bx,t0lx,t0hx);
  
 
  // Allocation of histograms - end

  Int_t iev=0;
  /* main loop (infinite) */
  for(;;) {
    struct eventHeaderStruct *event;
    eventTypeType eventT;
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) {break;}
    
    /* get next event (blocking call until timeout) */
    status=monitorGetEventDynamic((void **)&event);
    if (status==(int)MON_ERR_EOF) {
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
    
    /* use event - here, just write event id to result file */
    eventT=event->eventType;
    
    switch (event->eventType){
      
      case START_OF_RUN:
	break;
	
    case END_OF_RUN:
      break;
      
        case PHYSICS_EVENT:
	  //    case CALIBRATION_EVENT:
      iev++;
      
      if(iev==1){
	printf("First event - %i\n",iev);
      }
      
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      
      // Enable the following two lines in case of real-data
      reader->RequireHeader(kTRUE);
      AliT0RawReader *start = new AliT0RawReader(reader, kTRUE);
      
      // Read raw data
      Int_t allData[105][5];
      for(Int_t i0=0;i0<105;i0++)
      	for(Int_t j0=0;j0<5;j0++)
	  allData[i0][j0] = 0;
      
      if(start->Next()){
	for (Int_t i=0; i<105; i++) {
	  for(Int_t iHit=0;iHit<5;iHit++){
	    allData[i][iHit]= start->GetData(i,iHit);
	  }
	}
      }
      
      // Fill the histograms
      Float_t besttimeA=9999999;
      Float_t besttimeC=9999999;
      Float_t time[24]; 
       Float_t meanShift[24];
       for (Int_t ik = 0; ik<24; ik++)
	 { 
	   if(ik<12 && allData[ik+1][0]>0 
	      && (allData[ik+13][0]-allData[ik+1][0]) < 530 ){
	     hCFD1minCFD[ik]->Fill(allData[ik+1][0]-allData[npmtC][0]);
	   }
	   
	   if(ik>11 && allData[ik+45][0]>0 
	      && (allData[ik+57][0]-allData[ik+45][0]) <530 ){
	     hCFD1minCFD[ik]->Fill(allData[ik+45][0]-allData[56+npmtA][0]);
	   }
	   if(iev == 10000) {	
	     meanShift[ik] =  hCFD1minCFD[ik]->GetMean();  
	   }
	 }
      //fill  mean time _ fast reconstruction
      if (iev > 10000 )
	{
	  for (Int_t in=0; in<12; in++)  
	    {
	      time[in] = allData[in+1][0] - meanShift[in]  ;
	      time[in+12] = allData[in+56+1][0] ;
	    }
	  for (Int_t ipmt=0; ipmt<12; ipmt++){
	    if(time[ipmt] > 1 ) {
	      if(time[ipmt]<besttimeC)
		besttimeC=time[ipmt]; //timeC
	    }
	  }
	   for ( Int_t ipmt=12; ipmt<24; ipmt++){
	     if(time[ipmt] > 1) {
	       if(time[ipmt]<besttimeA) 
		 besttimeA=time[ipmt]; //timeA
	     }
	   }
	   if(besttimeA<9999999 &&besttimeC< 9999999) {
	     Float_t t0 =0.001* 24.4 * Float_t( besttimeA+besttimeC)/2.;
	     cout<<" t0 "<<t0<<endl;
	     hVertex->Fill(t0);
	   }
	}
	   
     delete start;
      start = 0x0;
      reader->Reset();
      // End of fill histograms
      
    }
    
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      printf("Number of events processed - %i\n ",iev); 	
      break;
    }
  }
  printf("After loop, before writing histos\n");
  // write a file with the histograms

  TFile *hist = new TFile(FILE_OUT,"RECREATE");

  for(Int_t j=0;j<24;j++){
     hCFD1minCFD[j]->Write();
    }
  hVertex->Write();
  hist->Close();
  delete hist;

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_OUT, "PHYSICS")) {
    status=-2;
  }

  return status;
}


