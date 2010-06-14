/*

TOF DA for online calibration

Contact: Chiara.Zampolli@bo.infn.it
         Roberto.Preghenella@bo.infn.it

Run Type: PHYSICS
DA Type: MON
Number of events needed:
Input Files: no input
Output Files: TOFdaCalib.root
Event types used: CALIBRATION_EVENT

*/

#define FILE_CALIB "TOFdaCalib.root"

// DATE
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

//ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TPluginManager.h"

//AliRoot
#include "AliLog.h"
#include "AliTOFRawStream.h"
#include "AliRawReaderDate.h"
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliTOFGeometry.h"
#include "AliTOFDecoderV2.h"
#include "AliTOFDecoderSummaryData.h"
#include "AliTOFDRMSummaryData.h"
#include "AliTOFTRMSummaryData.h"
#include "AliTOFChainSummaryData.h"
#include "AliTOFTDCHitBuffer.h"
#include "AliTOFTDCHit.h"

/* Main routine
      Arguments: 
      1- monitoring data source
*/
int 
main(int argc, char **argv) 
{
  
  /* magic line from Rene */
  gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					"*",
					"TStreamerInfo",
					"RIO",
					"TStreamerInfo()");
  

  /* log start of process */
  printf("TOF DA started\n");  
  
  /* check that we got some arguments = list of files */
  if (argc!=2) {
    printf("Wrong number of arguments\n");
    return -1;
  }


  /*
   * INIT 
   */

  /* constants */
  const Int_t nChannels = 157248;
  /* counters and flags */
  Int_t nCalibEvents;
  /* variables */
  Int_t ddl, slot, trm, chain, tdc, channel, index, det[5], dummy;
  /* TOF raw data handling */
  AliTOFRawStream *rawStream = new AliTOFRawStream();
  AliTOFDecoderV2 *decoder = rawStream->GetDecoderV2();
  AliTOFDecoderSummaryData *decodersd;
  AliTOFDRMSummaryData *drmsd;
  AliTOFTRMSummaryData *trmsd;  
  AliTOFChainSummaryData *chainsd;
  AliTOFTDCHitBuffer *hitBuffer;
  AliTOFTDCHit *hit;
  UChar_t *data = 0x0;
  Int_t dataSize;
  Int_t dataWords;
  Int_t currentDDL;

  /* init counters and flags */
  nCalibEvents = 0;
  
  /* open CALIB output file */
  TFile *fileOutCalib = new TFile(FILE_CALIB, "RECREATE"); 
  /* create calib hit histo */
  TH1F *hCalibHit = new TH1F("hCalibHit", "Calibration events;index;N_{hits}/N_{events}", nChannels, 0., nChannels);

  /*
   * ONLINE MONITOR
   */

  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  struct eventHeaderStruct *event;
  int ret;
  /* define monitoring table */
  char *monTable[5] = {
    "ALL", "no",
    "CAL", "yes",
    NULL
  };
  ret = monitorDeclareTable(monTable);
  if (ret != 0) {
    printf("monitorDeclareTable() failed: %s\n", monitorDecodeError(ret));
    return -1;
  }
  /* define data source : this is argument 1 */  
  ret = monitorSetDataSource(argv[1]);
  if (ret != 0) {
    printf("monitorSetDataSource() failed : %s\n",monitorDecodeError(ret));
    return -1;
  }
  /* declare monitoring program */
  ret = monitorDeclareMp("TOFdaCalib");
  if (ret != 0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(ret));
    return -1;
  }
  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);
  
  /* loop over events */
  while (1) {
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) break;
    
    /*
     * GET EVENT
     */

    /* get next event (blocking call until timeout) */
    ret = monitorGetEventDynamic((void **)&event);
    if (ret == MON_ERR_EOF) {
      printf ("End of File detected\n");
      break; /* end of monitoring file has been reached */
    }
    if (ret != 0) {
      printf("monitorGetEventDynamic() failed (ret=%d errno=%d): %s\n", ret, errno, monitorDecodeError(ret));
      break;
    }
    /* retry if got no event */
    if (event==NULL) continue;
    /* check TOF in partecipating detectors */
    if (!TEST_DETECTOR_IN_PATTERN(event->eventDetectorPattern, EVENT_DETECTOR_TOF)) {
      free(event);
      continue;
    }
    /* check event type */
    if (event->eventType != CALIBRATION_EVENT) {
      printf("not a calibration event: %d\n", event->eventType);
      free(event);
      continue;
    }
    /* increment number of calib events */
    nCalibEvents++;

    /*
     * DECODE EVENT
     */

    /* create and setup raw reader */
    AliRawReader *rawReader = new AliRawReaderDate((void *)event);
    rawReader->Reset();
    rawReader->Select("TOF", 0, AliDAQ::NumberOfDdls("TOF") - 1);
    /* setup raw stream */
    rawStream->SetRawReader(rawReader);

    /* loop over DDLs - rawReader->ReadHeader() */
    while (rawReader->ReadHeader()) {

      /* read equipment data */
      dataSize = rawReader->GetDataSize();
      data = new UChar_t[dataSize];
      if (!rawReader->ReadNext(data, dataSize)){
	delete [] data;
	continue;
      }

      /* decode data */
      dataWords = dataSize / 4;
      decoder->Decode((UInt_t *)data, dataWords);
      delete [] data;

      /* read equipment info */
      currentDDL = rawReader->GetDDLID();
      /* read decoder summary data */
      decodersd = decoder->GetDecoderSummaryData();
      /* check DRM header/trailer */
      drmsd = decodersd->GetDRMSummaryData();
      if (!drmsd->GetHeader() || !drmsd->GetTrailer()) continue;
      /* loop over TRM to get hits */
      for (Int_t itrm = 0; itrm < 10; itrm++) {
	trmsd = drmsd->GetTRMSummaryData(itrm);
	/* check header/trailer */
	if (!trmsd->GetHeader() || !trmsd->GetTrailer()) continue;
	/* loop over chains */
	for (Int_t ichain = 0; ichain < 2; ichain++) {
	  chainsd = trmsd->GetChainSummaryData(ichain);
	  /* check header/trailer */
	  if (!chainsd->GetHeader() || !chainsd->GetTrailer()) continue;
	  hitBuffer = chainsd->GetTDCPackedHitBuffer();

	  /*
	   * HIT MANIPULATION
	   */
	  
	  /* loop over hits in buffer */
	  for (Int_t ihit = 0; ihit < hitBuffer->GetEntries(); ihit++) {
	    
	    /* get hit */
	    hit = hitBuffer->GetHit(ihit);
	    
	    /* get channel info */
	    ddl = currentDDL;
	    slot = trmsd->GetSlotID();
	    trm = slot - 3;
	    chain = chainsd->GetChain();
	    tdc = hit->GetTDCID();
	    channel = hit->GetChan();
	    /* get index */
	    rawStream->EquipmentId2VolumeId(ddl, slot, chain, tdc, channel, det);
	    dummy = det[4];
	    det[4] = det[3];
	    det[3] = dummy;
	    /* check valid index */
	    if (det[0] < 0 || det[0] > 17 ||
		det[1] < 0 || det[1] > 5 ||
		det[2] < 0 || det[2] > 18 ||
		det[3] < 0 || det[3] > 1 ||
		det[4] < 0 || det[4] > 47) continue;
	    index = AliTOFGeometry::GetIndex(det);
	    
	    /* fill calib hit histo */
	    hCalibHit->Fill(index);
	    
	  } /* end of loop over hits in buffer */
	} /* end of loop over chains */
      } /* end of loop over TRMs */
    } /* end of loop over DDLs - rawReader->ReadHeader() */
    
    /* delete raw reader */
    delete rawReader;
    /* free event */
    free(event);
    
  } /* end of loop over events */

  /* scale calib hit histo by number of calib events */
  printf("found %d calibration events\n", nCalibEvents);
  hCalibHit->Sumw2();
  if (nCalibEvents > 0)
    hCalibHit->Scale(1. / nCalibEvents);
  
  /* write calib hit histo on CALIB file */
  fileOutCalib->cd();
  hCalibHit->Write();
  fileOutCalib->Close();
  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_CALIB, "CALIB"))
    return -2;

  return 0;
}
