/*

TOF DA for pulser data

Contact: Chiara.Zampolli@bo.infn.it
Link: www.bo.infn.it/~zampolli
Run Type: PULSER
DA Type: LDC
Number of events needed: 10000
Input Files: TOF<nrun>.raw, where <nrun> is the run number 
Output Files: TOFoutPulserLDC_<LDCid>.root, where <LDCid> is the id of the LDC which is read (2 characters field, e.g. TOFoutPulserLDC_03.root),to be exported to the DAQ FXS
Trigger types used: PHYSICS_EVENT (for the time being)

*/

// DATE
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include <AliTOFRawStream.h>
#include <AliRawReaderDate.h>
#include <AliRawReader.h>
#include <AliTOFGeometry.h>
#include <AliDAQ.h>
#include <AliTOFHitData.h>
#include <AliTOFHitDataBuffer.h>
#include <AliTOFDecoder.h>
#include <AliTOFNoiseConfigHandler.h>

//ROOT
#include <TFile.h>
#include <TKey.h>
#include <TH1I.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
#include "TROOT.h"
#include "TPluginManager.h"
#include "TSAXParser.h"

/* Main routine
      Arguments: list of DATE raw data files
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
  TH1F::AddDirectory(0);
  TH1I * htofPulser = new TH1I("hTOFpulser","histo with signals on TOF during pulser", size,-0.5,size-0.5);
  for (Int_t ibin =1;ibin<=size;ibin++)
    htofPulser->SetBinContent(ibin,-1);

  UInt_t ldcId=99;
  UInt_t ldcIdOLD=99;

  int status;

  /* log start of process */
  printf("TOF DA started\n");  

  /* check that we got some arguments = list of files */
  if (argc<2) {
    printf("Wrong number of arguments\n");
    return -1;
  }

  /* retrieve config file */
  int getConfigFile = daqDA_DB_getFile("TOFNoiseConfig.xml","TOFNoiseConfig.xml");
  if (getConfigFile != 0){
    printf("Failed to retrieve config file from DB! returning...\n");
    return -1;
  }

  AliTOFNoiseConfigHandler* tofHandler = new AliTOFNoiseConfigHandler();
  TSAXParser *parser = new TSAXParser();
  parser->ConnectToHandler("AliTOFNoiseConfigHandler", tofHandler);  
  if (parser->ParseFile("./TOFNoiseConfig.xml") != 0) {
    printf("Failed parsing config file! retunring... \n");
    return -1;
  }

  Int_t debugFlag = tofHandler->GetDebugFlag();
  printf("the debug flag is %i\n",debugFlag);

  /* init some counters */
  int nevents_physics=0;
  int nevents_total=0;

  Int_t nPDBEntriesToT = 0;
  Int_t nDBEntriesToT = 0;
  AliTOFHitData *HitData = 0;
  Int_t dummy = -1;
  Int_t Volume[5];
  for (Int_t i=0;i<5;i++) Volume[i]=-1;
  AliTOFRawStream *rawStreamTOF = new AliTOFRawStream();
  AliTOFDecoder * decoderTOF = new AliTOFDecoder();
  AliTOFHitDataBuffer *DataBuffer[AliDAQ::NumberOfDdls("TOF")];
  AliTOFHitDataBuffer *PackedDataBuffer[AliDAQ::NumberOfDdls("TOF")];
  for (Int_t i=0;i<AliDAQ::NumberOfDdls("TOF");i++) {
    DataBuffer[i]=new AliTOFHitDataBuffer();
    PackedDataBuffer[i]=new AliTOFHitDataBuffer();
  }
  Int_t currentEquipment;
  Int_t currentDDL;
  UChar_t *data = 0x0;
  Int_t nchDDL = 0;
  Int_t nDBEntries = 0;
  Int_t nPDBEntries = 0;

  struct eventHeaderStruct *event;
  eventTypeType eventT;

  /* read the data files */
  int n;
  for (n=1;n<argc;n++) {
   
    status=monitorSetDataSource( argv[n] );
    if (status!=0) {
      printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(status));
      return -1;
    }

    /* read the file */
    for(;;) {

      rawStreamTOF->Clear();
      currentEquipment = 0;
      currentDDL = 0;
      data = 0x0;
      nPDBEntriesToT = 0;
      nDBEntriesToT = 0;
      dummy = -1;
      nchDDL = 0;
      nDBEntries = 0;
      nPDBEntries = 0;
      HitData = 0;

      /* get next event */
      status=monitorGetEventDynamic((void **)&event);
      if (status==MON_ERR_EOF) break; /* end of monitoring file has been reached */
      if (status!=0) {
        printf("monitorGetEventDynamic() failed : %s\n",monitorDecodeError(status));
        return -1;
      }

      /* retry if got no event */
      if (event==NULL)
        break;

      /* use event - here, just write event id to result file */
      eventT=event->eventType;

      if (eventT==PHYSICS_EVENT) {
	//printf ("event %i \n", nevents_physics);

	AliRawReader *rawReader = new AliRawReaderDate((void*)event);
	rawStreamTOF->SetRawReader(rawReader);

	while (rawReader->ReadHeader()) {
	  ldcId = rawReader->GetLDCId();
	  //debugging printings
	  if (debugFlag && ldcId!=ldcIdOLD) {
	    ldcIdOLD = ldcId;
	    printf ("ldcId = %i \n",ldcId);
	  }

	  //memory leak prevention (actually data should be always 0x0 here)
	  if (data != 0x0) delete [] data;
	  
	  //get equipment infos
	  currentEquipment = rawReader->GetEquipmentId();
	  currentDDL = rawReader->GetDDLID();
	  const AliRawDataHeader *currentCDH = (AliRawDataHeader*)rawReader->GetDataHeader();
	  if (currentDDL%2==0)
	    nchDDL = 2160;
	  else
	    nchDDL = 2208;

	  //Int_t * array = new Int_t[nchDDL];
	  Int_t array[nchDDL];
	  for (Int_t ii=0; ii<nchDDL; ii++) array[ii] = 0;
	  decoderTOF->GetArrayDDL(array, currentDDL);

	  for (Int_t i=0;i<nchDDL;i++)
	    if (htofPulser->GetBinContent(array[i]+1)<0) htofPulser->SetBinContent(array[i]+1,0);

	  //debugging printings
	  //if (debugFlag) printf(" Equipment = %i, and DDL = %i \n", currentEquipment,currentDDL); 
	  const Int_t kDataSize = rawReader->GetDataSize();
	  const Int_t kDataWords = kDataSize / 4;
	  data = new UChar_t[kDataSize];
	  decoderTOF->SetDataBuffer(DataBuffer[currentDDL]);
	  decoderTOF->SetPackedDataBuffer(PackedDataBuffer[currentDDL]);
	  //start decoding
	  if (!rawReader->ReadNext(data, kDataSize)) {
	    rawReader->AddMajorErrorLog(AliTOFRawStream::kDDLdataReading);
	    printf("Error while reading DDL data. Go to next equipment \n");
	    delete [] data;
	    data = 0x0;
	    continue;
	  }
	  if (decoderTOF->Decode((UInt_t *)data, kDataWords, currentCDH) == kTRUE) {
	    rawReader->AddMajorErrorLog(AliTOFRawStream::kDDLDecoder,Form("DDL # = %d",currentDDL));
	    printf("Error while decoding DDL # %d: decoder returned with errors \n", currentDDL);
	  }

	  nDBEntries = DataBuffer[currentDDL]->GetEntries();
	  nPDBEntries = PackedDataBuffer[currentDDL]->GetEntries();
	  nPDBEntriesToT+=nPDBEntries;
	  nDBEntriesToT+=nDBEntries;
	  /* reset buffer */
	  DataBuffer[currentDDL]->Reset();

	  /* read data buffer hits */
	  for (Int_t iHit = 0; iHit < nPDBEntries; iHit++){
	    HitData = PackedDataBuffer[currentDDL]->GetHit(iHit);
	    /* add volume information */
	    HitData->SetDDLID(currentDDL);
	    for (Int_t i=0;i<5;i++) Volume[i]=-1;
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
	      Int_t index = geom->GetIndex(Volume);
	      // to check array indexes
	      /*
	      Bool_t found =kFALSE;
	      for (Int_t j=0;j<nchDDL;j++){
		if (index==array[j]) {
		  found = kTRUE;
		  break;
		}
	      }
	      printf ("index = %6d, found = %6d\n",index, (Int_t)found);
	      */
	      //printf ("index = %6d \n",index);
	      htofPulser->Fill(index); //channel index start from 0, bin index from 1
	      //debugging printings
	      //if (debugFlag) printf("sector %2d, plate %1d, strip %2d, padz %1d, padx %2d \n",Volume[0],Volume[1],Volume[2],Volume[3],Volume[4]); // too verbose
	    }
	  }
	  /* reset buffer */
	  PackedDataBuffer[currentDDL]->Reset();
	  delete [] data;
	  data = 0x0;
	  //delete [] array;
	}

	//debugging printings
	//if (debugFlag) {
	//  printf(" Packed Hit Buffer Entries = %i \n",nPDBEntriesToT); // too verbose
	//  printf(" Hit Buffer Entries = %i \n",nDBEntriesToT); // too verbose
	//}

	delete rawReader;
	rawReader = 0x0;
        nevents_physics++;
      }
      nevents_total++;

      /* free resources */
      free(event);
    }

  }

  delete decoderTOF;
  decoderTOF=0x0;

  for (Int_t i=0;i<AliDAQ::NumberOfDdls("TOF");i++){ 
    delete DataBuffer[i];
    delete PackedDataBuffer[i];
  }

  delete rawStreamTOF;
  rawStreamTOF = 0x0;

  delete geom;
  geom = 0x0;

  //write the Run level file   
  char filename[100];
  sprintf(filename,"TOFoutPulserLDC_%02i.root",ldcId);
  TFile * fileRun = new TFile (filename,"RECREATE"); 
  htofPulser->Write();
  fileRun->Close();

  /* write report */
  printf("Run #%s, received %d physics events out of %d\n",
	 getenv("DATE_RUN_NUMBER"),nevents_physics,nevents_total);

  status = 0;

  /* store the result file on FES */
  status=daqDA_FES_storeFile(filename,"PULSER");
  if (status)
    status = -2;

  return status;
}
