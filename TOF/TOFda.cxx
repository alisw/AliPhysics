/*

DAcase2.c

This program connects to the DAQ data source passed as argument
and populates local "./result.txt" file with the ids of events received
during the run.

The program exits when being asked to shut down (daqDA_checkshutdown)
or End of Run event.

Messages on stdout are exported to DAQ log system.

contact: alice-datesupport@cern.ch

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
#include "AliTOFRawStream.h"
#include "AliRawReaderDate.h"
#include "AliTOFGeometry.h"
#include "AliTOFGeometryV5.h"
#include "AliT0digit.h"
#include "AliT0RawReader.h"

//ROOT
#include "TFile.h"
#include "TKey.h"
#include "TH2S.h"
#include "TObject.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TRandom.h"


/* Main routine
      Arguments: 
      1- monitoring data source
*/
int main(int argc, char **argv) {
  
  AliTOFGeometry * geomV5 = new AliTOFGeometryV5();
  AliTOFGeometry * geom = new AliTOFGeometry();

  static const Int_t size = geomV5->NPadXSector()*geomV5->NSectors();
  static const Int_t nbins = 500;
  static const Int_t binmin = -20;
  const Float_t c = 2.99792458E10; //speed of light
  TH1F::AddDirectory(0);
  TH2S * htofPartial = new TH2S("htof","histo with delays", size,-0.5,size*1.-0.5,nbins,binmin-0.5,nbins*1.+binmin-0.5);
  
  TRandom *rand = new TRandom();  //to be used for testing with cosmic data
  rand->SetSeed(0);               //to be used for testing with cosmic data

  //TTree *tree = 0x0;  // tree for T0 decoder
      
  // decoding the events
  
  int status;
  bool decodeStatus;
  if (argc!=2) {
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
  printf("DA example case2 monitoring program started\n");  

  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  struct equipmentStruct *equipment;
  int *eventEnd;
  int *eventData;
  int *equipmentEnd;
  int *equipmentData;
  int *equipmentID;
  struct eventHeaderStruct *event;
  eventTypeType eventT;
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
      printf(" event number = %i \n",iev);
      //if (iev%10000 ==0) printf(" event number = %i \n",iev);
      //if (iev > 50000) break;
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      //      rawReader->RequireHeader(kFALSE);
      //rawReader->LoadEquipmentIdsMap("TOFmap.txt");  //to be used if no exact mapping implemented 
      Int_t meantime = 0;     

            
      //      TTree *tree = new TTree();  // tree for T0 decoder
      AliT0RawReader *rawReaderT0 = new AliT0RawReader(rawReader);
      //      printf("rawReaderT0 = %p \n", rawReaderT0);
      //printf("rawReader = %p \n", rawReader);
      //printf("event = %p \n", event);
      //      AliT0digit *digit = 0x0;
      //tree->SetBranchAddress("T0",&digit);
      
      if (!rawReaderT0->Next())
       printf(" no raw data found!! %i", rawReaderT0->Next());
      Int_t allData[110][5];
      for (Int_t i=0; i<110; i++) {
	allData[i][0]=rawReaderT0->GetData(i,0);
      }
      meantime = allData[49][0];
      printf("time zero (ns) = %i (%f) \n", meantime, meantime*25*1E-3-200);
      
      delete rawReaderT0;
      rawReaderT0 = 0x0;
      rawReader->Reset();
      
      AliTOFRawStream *rawStreamTOF = new AliTOFRawStream(rawReader);
      
      //loop the event data	  	  
      Int_t detectorIndex[5] = {-1, -1, -1, -1, -1};
      
      printf("before T0 \n");
      
      while(rawStreamTOF->Next()){
	for (Int_t ii=0; ii<5; ii++) detectorIndex[ii] = -1;
	
	detectorIndex[0] = (Int_t)rawStreamTOF->GetSector();
	detectorIndex[1] = (Int_t)rawStreamTOF->GetPlate();
	detectorIndex[2] = (Int_t)rawStreamTOF->GetStrip();
	detectorIndex[3] = (Int_t)rawStreamTOF->GetPadZ();
	detectorIndex[4] = (Int_t)rawStreamTOF->GetPadX();
	
	if (detectorIndex[0]==-1 ||
	    detectorIndex[1]==-1 ||
	    detectorIndex[2]==-1 ||
	    detectorIndex[3]==-1 ||
	    detectorIndex[4]==-1) continue;
	else {
	  Float_t tdc = (Float_t)rawStreamTOF->GetTDC();
	  Int_t tof = (Int_t)rawStreamTOF->GetTofBin();
	  Int_t detID[5]={0,0,0,0,0};
	  detectorIndex[0] = rawStreamTOF->GetSector(); 
	  detectorIndex[1] = rawStreamTOF->GetPlate();
	  detectorIndex[2] = rawStreamTOF->GetStrip(); 
	  detectorIndex[3] = rawStreamTOF->GetPadZ(); 
	  detectorIndex[4] = rawStreamTOF->GetPadX(); 
	  Int_t index = rawStreamTOF->GetIndex(detectorIndex);
	  Float_t pos[3];
 	  geomV5->GetPosPar(detectorIndex,pos);
	  Float_t texp=TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])/c*1E9; //expected time in ns
	  Float_t texpBin=(texp*1E3-32)/AliTOFGeometry::TdcBinWidth(); //expected time in number of TDC bin
	  //to be uded with cosmics
	  //Float_t tsim = (Float_t)rand->Gaus(texp,1E-1); //TDC measured time in ns, simulated to be used for testing with cosmic data  
	  //Float_t delta = tsim-texp;  
	  //Int_t deltabin = (Int_t)((delta*1E3-32)/AliTOFGeometry::TdcBinWidth());    //to be used for testing with cosmics data
	  Int_t deltabin = tof-TMath::Nint(texpBin);   //to be used with real data; rounding expected time to Int_t
	  htofPartial->Fill(index,deltabin); //channel index start from 0, bin index from 1
	  //debugging printings
	  printf("sector %i, plate %i, strip %i, padz %i, padx %i \n",detectorIndex[0],detectorIndex[1],detectorIndex[2],detectorIndex[3],detectorIndex[4]);
	  //printf("pos x = %f, pos y = %f, pos z = %f \n",pos[0],pos[1],pos[2]);
	  //printf ("expected time = %f (ns)\n",texp);
	  //printf ("expected time bin = %f (TDC bin)\n",texpBin);
	  printf ("measured time bin = %i (TDC bin)\n",tof);
	  //	  printf("index = %i, deltabin = %i , filling index = %i, and bin = % i\n",index, deltabin, index, deltabin);
	}
      }
      
      delete rawStreamTOF;
      rawStreamTOF = 0x0;
 
      delete rawReader;
      rawReader = 0x0;

    }
       
    /* free resources */
    free(event);
    
    /* exit when last event received, no need to wait for TERM signal */
    if (eventT==END_OF_RUN) {
      printf("EOR event detected\n");
      break;
    }
  }
  
  delete geomV5;
  geomV5 = 0x0;
  delete geom;
  geom = 0x0;

  //write the Run level file   
  TFile * fileRun = new TFile ("outciccioTOFT0daRun.root","RECREATE"); 
  TBenchmark *bench = new TBenchmark();
  bench->Start("t0");
  htofPartial->Write();
  bench->Stop("t0");
  bench->Print("t0");
  fileRun->Close();

  TFile * filetot = new TFile ("outciccioTOFT0daTotal.root","READ"); 
  printf("dopo aver aperto il file in modalita' lettura \n");
  TH2S *htoftot = 0x0;
  //look for the file
  if (!filetot->IsZombie()){
    printf("il file non e' zombie \n");
    TIter next(filetot->GetListOfKeys());
    TKey *key;
    //look for the histogram
    while ((key=(TKey*)next())){
      const char * namekey = key->GetName();
      if (strcmp(namekey,"htoftot")==0){
	printf(" histo found \n");
	htoftot = (TH2S*) filetot->Get("htoftot");
	htoftot->AddDirectory(0);
	htoftot->Add(htofPartial);
	break;
      }
    }
  filetot->Close();
  }
  else {
    printf(" no file found \n");
    htoftot = new TH2S(*htofPartial);
    htoftot->SetName("htoftot");
    htoftot->AddDirectory(0);
  }
  filetot=0x0;
  filetot = new TFile ("outciccioTOFT0daTotal.root","RECREATE");
  filetot->cd();
  htoftot->Write();
  filetot->Close();
  
  delete fileRun;
  delete filetot;
  delete bench;
  delete htofPartial;
  delete htoftot;

  fileRun = 0x0;
  filetot = 0x0;
  bench = 0x0;
  htofPartial = 0x0;
  htoftot = 0x0;

  /* write report */
  fprintf(fp,"Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  /* close result file */
  fclose(fp);


  return status;
}
