/*
T0 DA for online calibration
Contact: AllaMaevskaya@cern.ch
Run Type: PHYSICS
DA Type: MON
Number of events needed: 10000 
Input Files: inPhys.dat, external parameters, T0/Calib/Slewing_Walk
Output Files: daPhys.root, to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT

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
#include <AliT0CalibWalk.h>
#include <AliCDBManager.h>
#include <AliCDBEntry.h>
//ROOT
#include "TROOT.h"
#include "TPluginManager.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TString.h"
#include "TH1.h"
#include "TMath.h"


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
  int  kcbx, kt0bx, knpmtA, knpmtC,kcfdbx;
  float kclx,kcmx, kt0lx, kt0hx, kcfdlx, kcfdhx;
  
  while((c=getc(inp))!=EOF) {
    switch(c) {
    case 'a': {fscanf(inp, "%d", &kcbx ); break;} //N of X bins hCFD1minCFD
    case 'b': {fscanf(inp, "%f", &kclx ); break;} //Low x hCFD1minCFD
    case 'c': {fscanf(inp, "%f", &kcmx ); break;} //High x hCFD1minCFD
    case 'd': {fscanf(inp, "%d", &knpmtC ); break;} //number of reference PMTC
    case 'e': {fscanf(inp, "%d", &knpmtA ); break;} //number of reference PMTA
    case 'f': {fscanf(inp, "%d", &kcfdbx ); break;} //N of X bins hCFD&OR
    case 'g': {fscanf(inp, "%f", &kcfdlx ); break;} //Low x  hCFD&OR
    case 'k': {fscanf(inp, "%f", &kcfdhx ); break;} //High x  hCFD&OR
    case 'm': {fscanf(inp, "%d", &kt0bx ); break;} //N of X bins TVDC
    case 'n': {fscanf(inp, "%f", &kt0lx ); break;} //Low x  TVDC
    case 'q': {fscanf(inp, "%f", &kt0hx ); break;} //High x  TVDC
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
  
  // Get run number
  if (getenv("DATE_RUN_NUMBER")==0) {
    printf("DATE_RUN_NUMBER not properly set.\n");
    return -1;
  }
  int runNr = atoi(getenv("DATE_RUN_NUMBER"));

 // Get the necessary OCDB files from the DAQ detector DB
  if (gSystem->AccessPathName("localOCDB/T0/Calib/Slewing_Walk/",kFileExists)) {
    if (gSystem->mkdir("localOCDB/T0/Calib/Slewing_Walk/",kTRUE) != 0) {
      printf("Failed to create directory: localOCDB/T0/Calib/Slewing_Walk/");
      return -1;
    }
  }
  status = daqDA_DB_getFile("T0/Calib/Slewing_Walk/Run0,999999999_v0_s0.root","localOCDB/T0/Calib/Slewing_Walk/Run0_999999999_v0_s0.root");
  if (status) {
    printf("Failed to get  file T0/Calib/Slewing_Walk() from DAQdetDB, status=%d\n", status);
    return -1;
  }
 
  TGraph *gr[24]; TGraph *gramp[24];
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://localOCDB");
  man->SetRun(runNr);
  cout<<" run number "<<runNr<<endl;
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("T0/Calib/Slewing_Walk");
  if(entry) {
    AliT0CalibWalk *fParam = (AliT0CalibWalk*)entry->GetObject();
    for (Int_t i=0; i<24; i++) {
      gr[i] = fParam->GetWalk(i); 
      gramp[i] = fParam->GetQTC(i); 
    }
  }
  Int_t chargeQT0[24], chargeQT1[24];
  Float_t adc ,walk, amp;
 
  // Allocation of histograms - start

  TH1F *hCFD1minCFD[24];  
  TH1F *hCFD[24], *hQT1[24], *hPed[24];  
   
  for(Int_t ic=0; ic<24; ic++) {
    hCFD1minCFD[ic] = new TH1F(Form("CFD1minCFD%d",ic+1),"CFD-CFD",kcbx,kclx,kcmx);
    hCFD[ic] = new TH1F(Form("CFD%d",ic+1),"CFD",kcfdbx,kcfdlx,kcfdhx);
    hQT1[ic] = new TH1F(Form("QT1%d",ic+1),"QT1",kt0bx,kt0lx,kt0hx);
    hPed[ic] = new TH1F(Form("hPed%d",ic+1),"pedestal",500,500,2000);
  }
  TH1F *hVertex = new TH1F("hVertex","TVDC",kt0bx,kt0lx,kt0hx);
  TH1F *hOrA = new TH1F("hOrA","OrA",kcfdbx,kcfdlx,kcfdhx);
  TH1F *hOrC = new TH1F("hOrC","OrC",kcfdbx,kcfdlx,kcfdhx);
  
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
      //	 	 case CALIBRATION_EVENT: for test
      iev++;
      
      if(iev==1){
	printf("First event - %i\n",iev);
      }
      
      // Initalize raw-data reading and decoding
      AliRawReader *reader = new AliRawReaderDate((void*)event);
      
      // Enable the following two lines in case of real-data
      reader->RequireHeader(kTRUE);
      AliT0RawReader *start = new AliT0RawReader(reader, kTRUE);
      // start->SetPrintout(kFALSE);
      // Read raw data
      Int_t allData[250][5];
      for(Int_t i0=0;i0<250;i0++)
      	for(Int_t j0=0;j0<5;j0++)
	  allData[i0][j0] = 0;
      
      if(start->Next()){
	for (Int_t i=0; i<226; i++) {
	  for(Int_t iHit=0;iHit<5;iHit++){
	    allData[i][iHit]= start->GetData(i,iHit);
	  }
	}
      }
      if (allData[50][0] == 0)     continue;
      hVertex->Fill(allData[50][0]);
      hOrA->Fill(allData[51][0]);
      hOrC->Fill(allData[52][0]);
      // Fill the histograms
      walk = adc = amp = -999;
      for (Int_t in=0; in<12;  in++)
	{
	  chargeQT0[in]=allData[2*in+25][0];
	  chargeQT1[in]=allData[2*in+26][0];
	}	
      for (Int_t in=12; in<24;  in++)
	{
	  chargeQT0[in]=allData[2*in+57][0];
	  chargeQT1[in]=allData[2*in+58][0];
	}
      Float_t time[24]; 
      Float_t meanShift[24];
      for (Int_t ik = 0; ik<24; ik++)
	{ 	 
	  if( ( chargeQT0[ik] - chargeQT1[ik])>0)  {
	    adc = chargeQT0[ik] - chargeQT1[ik];
	    //	cout<<ik <<"  "<<adc<<endl;
	  }
	  if(gr[ik] && adc>0) 
	    walk = Int_t(gr[ik]->Eval(Double_t(adc) ) );
	  //	  cout<<ik<<" walk "<<walk<<" "<<allData[ik+1][0]<<endl; 
	  if(ik<12 && allData[ik+1][0]>0  ){
	    if( walk >-100) hCFD[ik]->Fill(allData[ik+1][0] - walk);
	    hQT1[ik]->Fill(chargeQT1[ik]);
	    if ( allData[knpmtC][0]>0 )
	      hCFD1minCFD[ik]->Fill(allData[ik+1][0]-allData[knpmtC][0]);
	  }
	  
	  if(ik>11 && allData[ik+45][0]>0 )
	    {
	      if( walk >-100) hCFD[ik]->Fill(allData[ik+45][0] - walk);
	      hQT1[ik]->Fill(chargeQT1[ik]);
	       if (allData[56+knpmtA][0]>0)
		 hCFD1minCFD[ik]->Fill(allData[ik+45][0]-allData[56+knpmtA][0]);
	    }
	  if(ik<12 && allData[ik+1][0]==0 && adc>0 ) hPed[ik]->Fill(adc);
	  if(ik>11 && allData[ik+45][0]==0 && adc>0 ) hPed[ik]->Fill(adc);

	}
	   
     delete start;
     start = 0x0;
     delete reader;
     reader= 0x0;
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

  TFile hist(FILE_OUT,"RECREATE");

  for(Int_t j=0;j<24;j++){
    hCFD1minCFD[j]->SetDirectory(&hist);
    hCFD1minCFD[j]->Write();
    hCFD[j]->Write();
    hQT1[j]->Write();
    hPed[j]->Write();
  }
  hVertex->Write();
  hOrA->Write();
  hOrC->Write();
  hist.Close();
  //delete hist;

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_OUT, "PHYSICS")) {
    status=-2;
  }

  return status;
}


