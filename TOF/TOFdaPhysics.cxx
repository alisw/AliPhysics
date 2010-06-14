/*

TOF DA for online calibration

Contact: Chiara.Zampolli@bo.infn.it
         Roberto.Preghenella@bo.infn.it

Run Type: PHYSICS
DA Type: MON
Number of events needed:
Input Files: no input
Output Files: TOFdaHits.root TOFdaReadout.root
Event types used: PHYSICS_EVENT

*/

#define FILE_HITS "TOFdaHits.root"
#define FILE_READOUT "TOFdaReadout.root"

#define READOUT_INFO_HISTO 1

// DATE
#include "event.h"
#include "monitor.h"
#include "daqDA.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

//AliRoot
#include "TROOT.h"
#include "AliTOFRawStream.h"
#include "AliRawReaderDate.h"
#include "AliRawReader.h"
#include "AliDAQ.h"
#include "AliTOFHitData.h"
#include "AliTOFHitDataBuffer.h"
#include "AliTOFDaConfigHandler.h"
#include "AliTOFHitField.h"
#include "AliLog.h"
#include "AliTOFGeometry.h"
#include "AliTOFDecoderV2.h"
#include "AliTOFDecoderSummaryData.h"
#include "AliTOFDRMSummaryData.h"
#include "AliTOFTRMSummaryData.h"
#include "AliTOFChainSummaryData.h"
#include "AliTOFTDCHitBuffer.h"
#include "AliTOFTDCHit.h"
#include "AliTOFTDCErrorBuffer.h"
#include "AliTOFTDCError.h"

//ROOT
#include "TFile.h"
#include "TKey.h"
#include "TH2S.h"
#include "TObject.h"
#include "TMath.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TSAXParser.h"
#include "TTree.h"

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
   * CONFIG
   */
  
  /* retrieve config file */
  int getConfigFile = daqDA_DB_getFile("TOFPhysicsConfig.xml","TOFPhysicsConfig.xml");
  if (getConfigFile != 0){
    printf("Failed to retrieve config file from DB! returning...\n");
    return -1;
  }
  /* parse config file */
  AliTOFDaConfigHandler* tofHandler = new AliTOFDaConfigHandler();
  TSAXParser *parser = new TSAXParser();
  parser->ConnectToHandler("AliTOFDaConfigHandler", tofHandler);  
  if (parser->ParseFile("./TOFPhysicsConfig.xml") != 0) {
    printf("Failed parsing config file! retunring... \n");
    return -1;
  }
  /* setup config params */
  Int_t meanMultiplicity = tofHandler->GetMeanMultiplicity(); /* average expected TOF multiplicity */
  Int_t maxHits = tofHandler->GetMaxHits(); /* max number of hits to be collected */
  printf("current settings:\n");
  printf(" - meanMultiplicity = %d\n", meanMultiplicity);
  printf(" - maxHits          = %d\n", maxHits);
  /* constants */
  const Int_t nChannels = 157248;
  Int_t noiseCheckTrigger = 10; /* first noise check after 10 events */
  Float_t meanChannelRate = (Float_t)meanMultiplicity / (Float_t)nChannels; /* average expected channel rate (hits/event) */
  Float_t noiseThreshold = 10. * meanChannelRate; /* noise threshold (hits/event) */
  Int_t minNoiseHits = 10; /* min number of channel hits to check noise */
  /* counters and flags */
  Int_t nPhysicsEvents, nCollectedPhysicsEvents, totHits;
  Int_t nChHits[nChannels];
  Bool_t inhibitCollection;
  Bool_t noiseFlag[nChannels];
  /* variables */
  Int_t ddl, slot, trm, chain, tdc, channel, index, timebin, totbin, deltaBC, l0l1latency, det[5], dummy;
  Float_t noiseHitThreshold;

  /*
   * INIT 
   */

  /* init counters and flags */
  nPhysicsEvents = 0;
  nCollectedPhysicsEvents = 0;
  totHits = 0;
  inhibitCollection = kFALSE;
  for (Int_t ich = 0; ich < nChannels; ich++) {
    nChHits[ich] = 0;
    noiseFlag[ich] = kFALSE;
  }

  /* TOF raw data handling */
  AliTOFRawStream *rawStream = new AliTOFRawStream();
  AliTOFDecoderV2 *decoder = rawStream->GetDecoderV2();
  AliTOFDecoderSummaryData *decodersd;
  AliTOFDRMSummaryData *drmsd;
  //  AliTOFLTMSummaryData *ltmsd;
  AliTOFTRMSummaryData *trmsd;  
  AliTOFChainSummaryData *chainsd;
  AliTOFTDCHitBuffer *hitBuffer;
  AliTOFTDCHit *hit;
  AliTOFTDCErrorBuffer *errorBuffer;
  AliTOFTDCError *error;
  UShort_t errorFlags;
  UChar_t *data = 0x0;
  Int_t dataSize;
  Int_t dataWords;
  Int_t currentDDL;
  const AliRawDataHeader *currentCDH;
  Int_t currentMiniEventID;
  Int_t currentEventID1;
  Int_t currentL0BCID ;
  Int_t currentBunchID;
  Bool_t skipTRM, skipChain;
  Int_t trmIndex, chainIndex, tdcIndex;
  Double_t chainEff;
  
  /* open HITS output file */
  TFile *fileOutHits = new TFile(FILE_HITS, "RECREATE"); 
  /* create hit field data structure */
  AliTOFHitField *hitField = new AliTOFHitField();
  /* create temporary tree */
  TTree *tempTree = new TTree("tempTree", "temporary tree");
  tempTree->Branch("hit", "AliTOFHitField", &hitField);
  /* create output tree */
  TTree *outTree = new TTree("hitTree", "hit tree");
  outTree->Branch("hit", "AliTOFHitField", &hitField);

  /* open READOUT output file */
  TFile *fileOutReadout = new TFile(FILE_READOUT, "RECREATE"); 
  /* create chain readout efficiency histo */
  TH1F *hChainEfficiency = new TH1F("hChainEfficiency", "Chain efficiency;chain;efficiency", 1440, 0., 1440.);

#if READOUT_INFO_HISTO
  /* create TRM data histo */
  TH1F *hTRMData = new TH1F("hTRMData", "TRM data;TRM;frequency", 720, 0., 720.);
  /* create TRM empty event frequency histo */
  TH1F *hTRMEmptyEvent = new TH1F("hTRMEmptyEvent", "TRM empty event error;TRM;frequency", 720, 0., 720.);
  /* create TRM bad event counter frequency histo */
  TH1F *hTRMBadEventCounter = new TH1F("hTRMBadEventCounter", "TRM bad event counter;TRM;frequency", 720, 0., 720.);
  /* create TRM bad CRC frequency histo */
  TH1F *hTRMBadCRC = new TH1F("hTRMBadCRC", "TRM bad CRC;TRM;frequency", 720, 0., 720.);

  /* create chain data histo */
  TH1F *hChainData = new TH1F("hChainData", "Chain data;chain;frequency", 1440, 0., 1440.);
  /* create chain bad status frequency histo */
  TH1F *hChainBadStatus = new TH1F("hChainBadStatus", "Chain bad status;chain;frequency", 1440, 0., 1440.);
  /* create chain bad event counter frequency histo */
  TH1F *hChainBadEventCounter = new TH1F("hChainBadEventCounter", "Chain bad event counter;chain;status;frequency", 1440, 0., 1440.);

  /* create TDC error frequency histo */
  TH1F *hTDCError = new TH1F("hTDCError", "TDC error;TDC;frequency", 21600, 0., 21600.);
  /* create TDC error flags frequency histo */
  TH2F *hTDCErrorFlags = new TH2F("hTDCErrorFlags", "TDC error flags;TDC;error flag;frequency", 21600, 0., 21600., 15, 0., 15);
#endif
  

  /*
   * ONLINE MONITOR
   */

  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  struct eventHeaderStruct *event;
  int ret;
  /* define monitoring table */
  char *monTable[5] = {
    "ALL", "no",
    "PHY", "yes",
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
  ret = monitorDeclareMp("TOFdaPhysics");
  if (ret != 0) {
    printf("monitorDeclareMp() failed : %s\n",monitorDecodeError(ret));
    return -1;
  }
  /* define wait event timeout - 1s max */
  monitorSetNowait();
  monitorSetNoWaitNetworkTimeout(1000);

  /* variables */
  
  /* loop over events */
  while (1) {
    
    /* check shutdown condition */
    if (daqDA_checkShutdown()) break;
    
    /*
     * NOISE CHECK
     */

    /* check inhibit collection */
    if (!inhibitCollection) {
      /* check number of events and check noise */
      if (nCollectedPhysicsEvents >= noiseCheckTrigger || totHits >= maxHits) {
	noiseHitThreshold = noiseThreshold * nCollectedPhysicsEvents;
	printf("noise check triggered after %d events: threshold is %f hits\n", nCollectedPhysicsEvents, noiseHitThreshold);
	/* loop over all channels */
	for (Int_t ich = 0; ich < nChannels; ich++) {
	  /* check */
	  if (nChHits[ich] < minNoiseHits || noiseFlag[ich] || nChHits[ich] < noiseHitThreshold) continue;
	  printf("channel %d tagged as noisy (%d hits): disabled\n", ich, nChHits[ich]);
	  noiseFlag[ich] = kTRUE;
	  totHits -= nChHits[ich];
	} /* end of loop over all channels */
	/* set new noise check trigger value */
	noiseCheckTrigger *= 10;
      } /* end of noise check */    
      
      /* inhibit hit collection when maximum number of hits exceeded */
      if (totHits >= maxHits) {
	printf("maximum number of hits exceeded (%d): inhibit hit collection\n", maxHits);
	inhibitCollection = kTRUE;
      }
    }
    
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
    if (event->eventType != PHYSICS_EVENT) {
      printf("not a physics event: %d\n", event->eventType);
      free(event);
      continue;
    }
    /* increment number of physics events */
    nPhysicsEvents++;
    if (!inhibitCollection) nCollectedPhysicsEvents++;

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
      currentCDH = rawReader->GetDataHeader();
      currentMiniEventID = currentCDH->GetMiniEventID();
      currentEventID1 = currentCDH->GetEventID1();

      /* read decoder summary data */
      decodersd = decoder->GetDecoderSummaryData();
      
      /* check DRM header/trailer */
      drmsd = decodersd->GetDRMSummaryData();
      if (!drmsd->GetHeader() || !drmsd->GetTrailer()) continue;

      /* get DRM data */
      currentL0BCID = drmsd->GetL0BCID();
      
      /* loop over TRM to get hits */
      for (Int_t itrm = 0; itrm < 10; itrm++) {
	trmsd = drmsd->GetTRMSummaryData(itrm);
	trmIndex = itrm + 10 * currentDDL;
	skipTRM = kFALSE;
	
	/* check header/trailer */
	if (!trmsd->GetHeader() || !trmsd->GetTrailer()) continue;
	
#if READOUT_INFO_HISTO
	/* fill TRM data */
	hTRMData->Fill(trmIndex);
	/* fill TRM empty event */
	if (trmsd->GetEBit() != 0) {
	  hTRMEmptyEvent->Fill(trmIndex);
	  skipTRM = kTRUE;
	}
	/* fill TRM bad event counter */
	if (trmsd->GetEventCounter() != drmsd->GetLocalEventCounter()) {
	  hTRMBadEventCounter->Fill(trmIndex);
	  skipTRM = kTRUE;
	}
	/* fill TRM bad CRC */
	if (trmsd->GetEventCRC() != trmsd->GetDecoderCRC()) {
	  hTRMBadCRC->Fill(trmIndex);
	  skipTRM = kTRUE;
	}
#else
	/* check bad condition and skip TRM */
	if ( trmsd->GetEBit() != 0 ||
	     trmsd->GetEventCounter() != drmsd->GetLocalEventCounter() ||
	     trmsd->GetEventCRC() != trmsd->GetDecoderCRC() ) continue;
#endif
	
	/* loop over chains */
	for (Int_t ichain = 0; ichain < 2; ichain++) {
	  chainsd = trmsd->GetChainSummaryData(ichain);
	  chainIndex = ichain + 2 * itrm + 20 * currentDDL;
	  skipChain = kFALSE;
	  
	  /* check header/trailer */
	  if (!chainsd->GetHeader() || !chainsd->GetTrailer()) continue;
	  currentBunchID = chainsd->GetBunchID();
	  hitBuffer = chainsd->GetTDCPackedHitBuffer();
	  errorBuffer = chainsd->GetTDCErrorBuffer();

#if READOUT_INFO_HISTO
	  /* fill chain data */
	  hChainData->Fill(chainIndex);
	  /* check chain bad status */
	  if (chainsd->GetStatus() != 0) {
	    hChainBadStatus->Fill(chainIndex);
	    skipChain = kTRUE;
	  }
	  /* check chain bad event counter */
	  if (chainsd->GetEventCounter() != drmsd->GetLocalEventCounter()) {	
	    hChainBadEventCounter->Fill(chainIndex);
	    skipChain = kTRUE;
	  }
	  /* fill TDC error frequency histo */
	  for (Int_t ierr = 0; ierr < errorBuffer->GetEntries(); ierr++) {
	    error = errorBuffer->GetError(ierr);
	    tdc = error->GetTDCID();
	    tdcIndex = tdc + 15 * ichain + 30 * itrm + 300 * currentDDL;
	    hTDCError->Fill(tdcIndex); 
	    errorFlags = error->GetErrorFlags();
	    for (Int_t ierflg = 0; ierflg < 15; ierflg++)
	      if (errorFlags & (1 << ierflg))
		hTDCErrorFlags->Fill(tdcIndex, ierflg); 
	  }
#else
	  /* check bad condition and skip chain */
	  if ( chainsd->GetStatus() != 0 ||
	       chainsd->GetEventCounter() != drmsd->GetLocalEventCounter() ) continue;
#endif
	  
	  /*
	   * CHAIN READOUT EFFICIENCY
	   */
	  
	  /* compute number of available channels removing TDCs in error */
	  chainEff = (120. - 8. * errorBuffer->GetEntries()) / 120.;
	  /* fill chain readout efficiency histo */
	  if (!skipTRM && !skipChain)
	    hChainEfficiency->Fill(chainIndex, chainEff);
	  
	  /*
	   * HIT MANIPULATION
	   */
	  
	  /* check inhibit collection */
	  if (inhibitCollection) continue;

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
	    
	    /* check noise flag */
	    if (noiseFlag[index]) continue;
	    /* increment number of channel hits and total hits */
	    nChHits[index]++;
	    totHits++;
	    /* get signal info */
	    timebin = hit->GetHitTime();
	    totbin = hit->GetTOTWidth();
	    deltaBC = currentBunchID - currentEventID1;
	    l0l1latency = currentMiniEventID - currentL0BCID;
	    /* set hit field data */
	    hitField->SetIndex(index);
	    hitField->SetTimeBin(timebin);
	    hitField->SetTOTBin(totbin);
	    hitField->SetDeltaBC(deltaBC);
	    hitField->SetL0L1Latency(l0l1latency);
	    /* fill temp tree */
	    tempTree->Fill();

	  } /* end of loop over hits in buffer */
	} /* end of loop over chains */
      } /* end of loop over TRMs */
    } /* end of loop over DDLs - rawReader->ReadHeader() */
    
    /* delete raw reader */
    delete rawReader;
    /* free event */
    free(event);
    
  } /* end of loop over events */
  
  /* final noise check */
  noiseHitThreshold = noiseThreshold * nCollectedPhysicsEvents;
  printf("final noise check after collectiong %d events: threshold is %f hits\n", nCollectedPhysicsEvents, noiseHitThreshold);
  /* loop over all channels */
  for (Int_t ich = 0; ich < nChannels; ich++) {
    /* check */
    if (nChHits[ich] < minNoiseHits || noiseFlag[ich] || nChHits[ich] < noiseHitThreshold) continue;
    printf("channel %d tagged as noisy (%d hits): disabled\n", ich, nChHits[ich]);
    noiseFlag[ich] = kTRUE;
    totHits -= nChHits[ich];
  } /* end of loop over all channels */
  
  /* copy hits into output tree from temp tree */
  printf("copy hits from temporary tree into output tree\n");
  printf("temporary tree contains %d hits\n", (Int_t)tempTree->GetEntries());
  for (Int_t ihit = 0; ihit < tempTree->GetEntries(); ihit++) {
    /* get entry */
    tempTree->GetEntry(ihit);
    /* check noise flag */
    if (noiseFlag[hitField->GetIndex()]) continue;
    /* fill output tree */
    outTree->Fill();
  } /* end of copy hits into output tree from temp tree */
  printf("output tree contains %d hits\n", (Int_t)outTree->GetEntries());

  /* write output tree on HITS file */
  fileOutHits->cd();
  outTree->Write();
  fileOutHits->Close();
  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_HITS, "HITS"))
    return -2;

#if READOUT_INFO_HISTO
  hTRMData->Sumw2();
  hTRMEmptyEvent->Sumw2();
  hTRMBadEventCounter->Sumw2();
  hTRMBadCRC->Sumw2();

  hChainData->Sumw2();
  hChainBadStatus->Sumw2();
  hChainBadEventCounter->Sumw2();

  hTDCError->Sumw2();
  hTDCErrorFlags->Sumw2();

  /* divide histos */
  if (nPhysicsEvents > 0) {
    hTRMEmptyEvent->Divide(hTRMData);
    hTRMBadEventCounter->Divide(hTRMData);
    hTRMBadCRC->Divide(hTRMData);
    hTRMData->Scale(1. / nPhysicsEvents);

    hChainBadStatus->Divide(hChainData);
    hChainBadEventCounter->Divide(hChainData);
    hChainData->Scale(1. / nPhysicsEvents);

    hTDCError->Scale(1. / nPhysicsEvents);
    hTDCErrorFlags->Scale(1. / nPhysicsEvents);
  }
#endif
  
  /* scale chain efficiency by number of physics events */
  printf("found %d physics events\n", nPhysicsEvents);
  hChainEfficiency->Sumw2();
  if (nPhysicsEvents > 0) hChainEfficiency->Scale(1. / (nPhysicsEvents));
  
  /* write efficiency histo on READOUT file */
  fileOutReadout->cd();
#if READOUT_INFO_HISTO
  hTRMData->Write();
  hTRMEmptyEvent->Write();
  hTRMBadEventCounter->Write();
  hTRMBadCRC->Write();
  
  hChainData->Write();
  hChainBadStatus->Write();
  hChainBadEventCounter->Write();
  
  hTDCError->Write();
  hTDCErrorFlags->Write();
#endif
  hChainEfficiency->Write();
  fileOutReadout->Close();
  /* export file to FXS */
  if (daqDA_FES_storeFile(FILE_READOUT, "READOUT"))
    return -2;

  return 0;
}
