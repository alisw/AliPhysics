/*

TOF DA for online calibration

Contact: Chiara.Zampolli@bo.infn.it
Link: www.bo.infn.it/~zampolli
Run Type: PHYSICS
DA Type: MON
Number of events needed: depending on the run, being run-level
Input Files: TOFdaTotal.root, to be updated if existing
Output Files: TOFdaRun.root, TOFdaTotal.root, both to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT

*/

#define FILE_TOTAL "TOFdaTotal.root"
#define FILE_RUN "TOFdaRun.root"

// DATE
#include <daqDA.h>
#include <event.h>
#include <monitor.h>

//AliRoot
#include <AliTOFRawStream.h>
#include <AliRawReaderDate.h>
#include <AliRawReader.h>
#include <AliTOFGeometry.h>
#include <AliT0RawReader.h>
#include <AliDAQ.h>
#include <AliTOFHitData.h>
#include <AliTOFHitDataBuffer.h>

//ROOT
#include <TFile.h>
#include <TKey.h>
#include <TH2S.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
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

  AliTOFGeometry * geom = new AliTOFGeometry();

  static const Int_t size = AliTOFGeometry::NPadXSector()*AliTOFGeometry::NSectors();
  static const Int_t nbins = 500;
  static const Int_t binmin = -20;
  const Float_t c = 2.99792458E10; //speed of light
  TH1F::AddDirectory(0);
  TH2S * htofPartial = new TH2S("htof","histo with delays", size,-0.5,size*1.-0.5,nbins,binmin-0.5,nbins*1.+binmin-0.5);
  
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
  printf("TOF DA started\n");  

  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

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
      AliRawReader *rawReader = new AliRawReaderDate((void*)event);
      //rawReader->RequireHeader(kFALSE);

      //T0 event
      Int_t meantime = 0;     
      AliT0RawReader *rawReaderT0 = new AliT0RawReader(rawReader,kTRUE);
      if (!rawReaderT0->Next()) {
        printf("T0: no raw data found!\n");
      } 
      else {
	/*
        Int_t allData[105][5];
        for (Int_t i=0; i<105; i++) {
	  allData[i][0]=rawReaderT0->GetData(i,0);
	}
        meantime = allData[49][0];
	*/
        //meantime = rawReaderT0->GetData(49,0); //OLD
        meantime = (rawReaderT0->GetData(51,0)+rawReaderT0->GetData(52,0))/2.; //Alla
	//        printf("time zero (ns) = %i (%f) \n", meantime, (meantime*24.4-200)*1E-3);   // debugging purpose
      }
      
      delete rawReaderT0;
      rawReaderT0 = 0x0;
      rawReader->Reset();
      
      //TOF event
      Int_t dummy = -1;
      Int_t Volume[5];
      AliTOFHitData *HitData;
      AliTOFHitDataBuffer DataBuffer;
      AliTOFHitDataBuffer PackedDataBuffer;
      AliTOFRawStream *rawStreamTOF = new AliTOFRawStream(rawReader);
      //      rawReader->ReadHeader();
      rawStreamTOF->ResetBuffers();
      rawStreamTOF->DecodeDDL(0, AliDAQ::NumberOfDdls("TOF") - 1,0);
      Int_t nPDBEntriesToT = 0;
      Int_t nDBEntriesToT = 0;
      for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls("TOF"); iDDL++){
	
	/* read decoded data */
	DataBuffer = rawStreamTOF->GetDataBuffer(iDDL);
	PackedDataBuffer = rawStreamTOF->GetPackedDataBuffer(iDDL);
	
	/* get buffer entries */
	Int_t nDBEntries = DataBuffer.GetEntries();
	Int_t nPDBEntries = PackedDataBuffer.GetEntries();
	nPDBEntriesToT+=nPDBEntries;
	nDBEntriesToT+=nDBEntries;
	
	//for (Int_t iHit = 0; iHit < nDBEntries; iHit++){
	// HitData = DataBuffer->GetHit(iHit);
	  /* store volume information */
	//  rawStreamTOF->EquipmentId2VolumeId(HitData, Volume);
	//}
	/* reset buffer */
	DataBuffer.Reset();
	
	/* read data buffer hits */
	for (Int_t iHit = 0; iHit < nPDBEntries; iHit++){
	  HitData = PackedDataBuffer.GetHit(iHit);
	  /* add volume information */
	  HitData->SetDDLID(iDDL);
	  rawStreamTOF->EquipmentId2VolumeId(HitData, Volume);
	  if (Volume[0]==-1 ||
	      Volume[1]==-1 ||
	      Volume[2]==-1 ||
	      Volume[3]==-1 ||
	      Volume[4]==-1) continue;
	  else {
	    dummy = Volume[3];
	    Volume[3] = Volume[4];
	    Volume[4] = dummy;
	    Int_t tof = (Int_t)((Double_t)HitData->GetTime()*1E3/AliTOFGeometry::TdcBinWidth());
	    Int_t index = geom->GetIndex(Volume);
	    Float_t pos[3];
	    geom->GetPosPar(Volume,pos);
	    Float_t texp=TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])/c*1E9; //expected time in ns
	    Float_t texpBin=(texp*1E3-32)/AliTOFGeometry::TdcBinWidth(); //expected time in number of TDC bin
	    Int_t deltabin = tof-TMath::Nint(texpBin);   //to be used with real data; rounding expected time to Int_t
	    htofPartial->Fill(index,deltabin); //channel index start from 0, bin index from 1
	    //debugging printings
	    //printf("sector %i, plate %i, strip %i, padz %i, padx %i \n",Volume[0],Volume[1],Volume[2],Volume[3],Volume[4]);
	    //printf("pos x = %f, pos y = %f, pos z = %f \n",pos[0],pos[1],pos[2]);
	    //printf ("expected time = %f (ns)\n",texp);
	    //printf ("expected time bin = %f (TDC bin)\n",texpBin);
	    //printf ("measured time bin = %i (TDC bin) with %f (ns) and ACQ bit = %i \n",tof, HitData->GetTime(), HitData->GetACQ());
	    //	  printf("index = %i, deltabin = %i , filling index = %i, and bin = % i\n",index, deltabin, index, deltabin);
	    
	  }
	  /* reset buffer */
	  PackedDataBuffer.Reset();
	}
      }
      //printf(" Packed Hit Buffer Entries = %i \n",nPDBEntriesToT);
      //printf(" Hit Buffer Entries = %i \n",nDBEntriesToT);
   
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
  
  delete geom;
  geom = 0x0;

  //write the Run level file   
  TFile * fileRun = new TFile (FILE_RUN,"RECREATE"); 
  htofPartial->Write();
  fileRun->Close();

  //write the Total file
  TH2S *htoftot = 0x0;
  TFile * filetot = 0x0;
  Bool_t isThere=kFALSE;
  const char *dirname = "./";
  TString filename = FILE_TOTAL;
  if((gSystem->FindFile(dirname,filename))!=NULL){
    isThere=kTRUE;
    printf("%s found \n",FILE_TOTAL);
  }
  if(isThere){

    TFile * filetot1 = new TFile (FILE_TOTAL,"READ"); 
    //look for the file
    if (!filetot1->IsZombie()){
      printf("updating file %s \n",FILE_TOTAL);
      TIter next(filetot1->GetListOfKeys());
      TKey *key;
      //look for the histogram
      while ((key=(TKey*)next())){
	const char * namekey = key->GetName();
	if (strcmp(namekey,"htoftot")==0){
	  printf(" histo found \n");
	  htoftot = (TH2S*) filetot1->Get("htoftot");
	  htoftot->AddDirectory(0);
	  htoftot->Add(htofPartial);
	  break;
	}
      }
    }
    filetot1->Close();
    delete filetot1;
    filetot1=0x0;
  }
  else {
    printf(" no %s file found \n",FILE_TOTAL);
    htoftot = new TH2S(*htofPartial);
    htoftot->SetName("htoftot");
    htoftot->AddDirectory(0);
  }

  filetot = new TFile (FILE_TOTAL,"RECREATE");
  filetot->cd();
  htoftot->Write();
  filetot->Close();
  
  delete fileRun;
  delete filetot;
  delete htofPartial;
  delete htoftot;

  fileRun = 0x0;
  filetot = 0x0;
  htofPartial = 0x0;
  htoftot = 0x0;

  /* write report */
  printf("Run #%s, received %d physics events out of %d\n",getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  status=0;

  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_RUN, "RUNLevel")) {
    status=-2;
  }
  if (daqDA_FES_storeFile(FILE_TOTAL, "DELAYS")) {
    status=-2;
  }

  
  return status;
}

