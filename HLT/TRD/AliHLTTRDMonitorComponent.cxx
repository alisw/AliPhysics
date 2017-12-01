// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Felix Rettig, Stefan Kirsch                           *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTRDMonitorComponent.cxx
/// @author Felix Rettig, Stefan Kirsch
/// @date   2012-08-16
/// @brief  The TRD monitoring component
///

#include <cstdlib>
#include "TH1I.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "AliESDtrack.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTLogging.h"
#include "AliHLTTRDMonitorComponent.h"
#include "AliTRDonlineTrackingDataContainer.h"

#define TRDMODULES 18

ClassImp(AliHLTTRDMonitorComponent)

#define LogError( ... ) { HLTError(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("ERROR", __VA_ARGS__); } }
#define LogInfo( ... ) { HLTInfo(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INFO", __VA_ARGS__); } }
#define LogInspect( ... ) { HLTDebug(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INSPECT", __VA_ARGS__); } }
#define LogDebug( ... ) { if (fDebugLevel >= 1) { HLTInfo(__VA_ARGS__); DbgLog("DEBUG", __VA_ARGS__); } }

AliHLTTRDMonitorComponent::AliHLTTRDMonitorComponent() : AliHLTProcessor(),
  fTrackHighPtThreshold(2.3),
  fHistoMode(1),
  fTrackingDataDebugOutput(kFALSE),
  fDebugLevel(0),
  fWriteHistos(kFALSE),
  fEventId(fgkInvalidEventId),
  fTrackingData(NULL),
  fHistArray(NULL),
  fHistEventTypes(NULL),
  fHistTrackletY(NULL),
  fHistTrackletDy(NULL),
  fHistTrackletYDy(NULL),
  fHistTrackletZ(NULL),
  fHistTrackletPID(NULL),
  fHistTrackletsHCId(NULL),
  fHistTrackPt(NULL),
  fHistTrackPID(NULL),
  fHistTrackLayers(NULL),
  fHistTrackLayersHighPt(NULL),
  fHistTracksStack(NULL),
  fHistTrackletTimingStack(NULL),
  fHistTrackingTiming(NULL),
  fHistTriggerContribs(NULL)
{
}

AliHLTTRDMonitorComponent::~AliHLTTRDMonitorComponent() {
}

const char* AliHLTTRDMonitorComponent::GetComponentID() {
  return "TRDMonitorComponent";
}

void AliHLTTRDMonitorComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  list.push_back(kAliHLTDataTypeTObject | kAliHLTDataOriginTRD);
}

AliHLTComponentDataType AliHLTTRDMonitorComponent::GetOutputDataType() {
  return (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTRD);
}

void AliHLTTRDMonitorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  constBase = 10000000;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDMonitorComponent::Spawn() {
  return new AliHLTTRDMonitorComponent;
}

int AliHLTTRDMonitorComponent::Configure(const char* /*arguments*/) {
  return 0;
}

int AliHLTTRDMonitorComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) {
  return 0;
}

int AliHLTTRDMonitorComponent::DoInit(int argc, const char** argv) {

  int iResult = 0;

  do {

    fTrackingData = new AliTRDonlineTrackingDataContainer();
    if (!fTrackingData) {
      iResult = -ENOMEM;
      break;
    }
    fTrackingData->SetGtuPtMultiplier(-1.); // this component does not know about the B-field direction

    fHistArray = new TObjArray(25);
    if(!fHistArray)
      return -ENOMEM;
    fHistArray->SetOwner(kTRUE);

  } while (0);

  if (iResult < 0){

    if (fHistArray) delete fHistArray;
    fHistArray = NULL;

    if (fTrackingData) delete fTrackingData;
    fTrackingData = NULL;

  }

  fHistEventTypes = new TH1I("trdmon_event_types", "event types analysed;event type;abundance", 10, 0, 10);
  fHistEventTypes->GetXaxis()->SetBinLabel(1, "any");
  fHistEventTypes->GetXaxis()->SetBinLabel(2, "data");
  fHistEventTypes->GetXaxis()->SetBinLabel(3, "track(lets) present");
  fHistArray->AddLast(fHistEventTypes);

  fHistTrackletY = new TH1I("trdmon_tracklet_y", "tracklet y-position;y (cm);abundance", 125, -75, 75);
  fHistArray->AddLast(fHistTrackletY);

  fHistTrackletDy = new TH1I("trdmon_tracklet_dy", "tracklet deflection;#Delta y (140 #mu m);abundance", 140, -70, 70);
  fHistArray->AddLast(fHistTrackletDy);

  fHistTrackletYDy = new TH2I("trdmon_tracklet_y_dy", "tracklet y-deflection vs. y-position;y (160 #mu m);#Delta y (140 #mu m);", 256, -4096, 4096, 140, -70, 70);
  fHistArray->AddLast(fHistTrackletYDy);

  fHistTrackletZ = new TH1I("trdmon_tracklet_z", "tracklet z-position;z (padrow);abundance", 16, 0, 16);
  fHistArray->AddLast(fHistTrackletZ);

  fHistTrackletPID = new TH1I("trdmon_tracklet_pid", "tracklet PID;PID (a.u.);abundance", 256, 0, 256);
  fHistArray->AddLast(fHistTrackletPID);

  fHistTrackletsHCId = new TH2F("trdmon_tracklets_hc", "tracklet number by HC;TRD sector;TRD half-chamber",
				18, 0, 18, 60, 0, 60);
  fHistArray->AddLast(fHistTrackletsHCId);

  fHistTrackletTimingStack = new TH2I("trdmon_tracklet_timing_stack", "tracklet timing;TRD stack; time after L0 (us)",
				       90, 0., 90, 270, 0., 9.);
  fHistArray->AddLast(fHistTrackletTimingStack);

  fHistTrackPt = new TH1I("trdmon_track_pt", "p_{T} of TRD online tracks;p_{T}^{TRD, online} (GeV/c);abundance", 200, -20, 20.);
  fHistArray->AddLast(fHistTrackPt);

  fHistTrackPID = new TH1I("trdmon_track_pid", "PID of TRD online tracks;PID^{TRD, online} (a.u.);abundance", 256, 0, 256);
  fHistArray->AddLast(fHistTrackPID);

  fHistTrackLayers = new TH1I("trdmon_track_layers", "contributing layers to TRD online tracks;contributing layers;abundance", 7, 0, 7);
  fHistArray->AddLast(fHistTrackLayers);

  fHistTrackLayersHighPt = new TH1I("trdmon_track_layers_hpt", "contributing layers to TRD online tracks;contributing layers;abundance", 7, 0, 7);
  fHistArray->AddLast(fHistTrackLayersHighPt);

  fHistTracksStack = new TH2F("trdmon_tracks_stack", "tracks by stack;TRD sector;TRD stack",
				18, 0, 18, 5, 0, 5);
  fHistArray->AddLast(fHistTracksStack);

  const Double_t trackingTimesTimeBin = 0.025; // us
  const Double_t trackingTimesMaxTime = 12.;   // us
  fHistTrackingTiming = new TH2I("trdmon_tracking_timing", "tracking timing;;time after L0 (#mu s);",
	       4, 0, 4, (Int_t)(trackingTimesMaxTime/trackingTimesTimeBin), 0., trackingTimesMaxTime);
  fHistTrackingTiming->GetXaxis()->SetBinLabel(1, "tracklet start");
  fHistTrackingTiming->GetXaxis()->SetBinLabel(2, "tracklet end");
  fHistTrackingTiming->GetXaxis()->SetBinLabel(3, "stack done");
  fHistTrackingTiming->GetXaxis()->SetBinLabel(4, "sector done");
  fHistArray->AddLast(fHistTrackingTiming);

  fHistTriggerContribs = new TH2I("trdmon_trigger_contribs", "TRD internal contributions by sector;TRD sector;trigger contribution;",
				  18, 0, 18, 12, 0, 12);
  fHistTriggerContribs->GetYaxis()->SetBinLabel(1, "trg0");
  fHistTriggerContribs->GetYaxis()->SetBinLabel(2, "trg1");
  fHistTriggerContribs->GetYaxis()->SetBinLabel(3, "trg2");
  fHistTriggerContribs->GetYaxis()->SetBinLabel(4, "trg3");
  fHistTriggerContribs->GetYaxis()->SetBinLabel(5, "trg4");
  fHistTriggerContribs->GetYaxis()->SetBinLabel(6, "trg5");
  fHistTriggerContribs->GetYaxis()->SetBinLabel(7, "trg5");
  fHistTriggerContribs->GetYaxis()->SetBinLabel(8, "T");
  fHistArray->AddLast(fHistTriggerContribs);

 vector<const char*> remainingArgs;
 for (int i = 0; i < argc; ++i)
   remainingArgs.push_back(argv[i]);

 if (argc > 0)
   ConfigureFromArgumentString(remainingArgs.size(), &(remainingArgs[0]));

  return 0;
}

int AliHLTTRDMonitorComponent::DoDeinit() {

  if ((fHistoMode == 1) && (fWriteHistos)){
    TFile out("mon_out/mon_hists.root", "RECREATE");
    if (!out.IsZombie()) {
      out.cd();
      UInt_t numHists = fHistArray->GetEntries();
      for (UInt_t iHist = 0; iHist < numHists; ++iHist)
	if (fHistArray->At(iHist))
	  fHistArray->At(iHist)->Write();
      out.Close();
    }
  }

  if (fHistArray) delete fHistArray;
  fHistArray = NULL;

  if (fTrackingData) delete fTrackingData;
  fTrackingData = NULL;

  return 0;
}

int AliHLTTRDMonitorComponent::ScanConfigurationArgument(int argc, const char** argv)
{

  if (argc <= 0)
    return 0;

  UShort_t iArg = 0;
  TString argument(argv[iArg]);

  if (!argument.CompareTo("-write-histograms")){
    LogInfo("writing of histograms enabled.");
    fWriteHistos = kTRUE; // enable histogram writing, for debugging/tuning only!
    return 1;
  }

  if (!argument.CompareTo("-debug")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fDebugLevel = argument.Atoi();
    LogInfo("debug level set to %d.", fDebugLevel);
    return 2;
  }

  return 0;

}

void AliHLTTRDMonitorComponent::DbgLog(const char* prefix, ...){
#ifdef __TRDHLTDEBUG
  AliHLTEventID_t eventNumber = fEventId;
  Int_t fRunNumber = -1;
  printf("TRDHLTGM %s-X-%s: [MON] %s",
 	 (fRunNumber >= 0) ? Form("%06d", fRunNumber) : "XXXXXX",
 	 (eventNumber != fgkInvalidEventId) ? Form("%05llu", eventNumber) : "XXXXX",
 	 (strlen(prefix) > 0) ? Form("<%s> ", prefix) : "");
#endif
  va_list args;
  va_start(args, prefix);
  char* fmt = va_arg(args, char*);
  vprintf(fmt, args);
  printf("\n");
  va_end(args);
}

int AliHLTTRDMonitorComponent::PrepareTRDData() {

  int result = 1;

  fTrackingData->Clear();
  for (const AliHLTComponentBlockData* datablock = GetFirstInputBlock(AliHLTTRDDefinitions::fgkOnlineDataType);
       datablock != NULL;
       datablock = GetNextInputBlock())
    {
      fTrackingData->Decompress(datablock->fPtr, datablock->fSize, kTRUE);
    }

  fTrackingData->PrintSummary("monitor component");

  return result;

}

void AliHLTTRDMonitorComponent::DumpTrackingData(){
  TString trklStr("");
  TString matchStr("");
  UShort_t layerMask;

  if (fTrackingData->GetNumTracklets() + fTrackingData->GetNumTracks() == 0)
    return;

  for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack){
    for (Int_t iTrk = 0; iTrk < fTrackingData->GetNumTracks(iStack); ++iTrk){

      layerMask = fTrackingData->GetTrackLayerMask(iStack, iTrk);
      trklStr = Form("trkl: ");
      for (Short_t iLayer = 5; iLayer >= 0; --iLayer){
	if ((layerMask >> iLayer) & 1)
	  trklStr += Form("0x%08x (%+8.3f)  ",
			  fTrackingData->GetTrackTrackletWord(iStack, iTrk, iLayer),
			  fTrackingData->GetTrackTrackletLocalY(iStack, iTrk, iLayer));
	else
	  trklStr += "---------------------  ";
      } // loop over layers
      trklStr.Remove(trklStr.Length() - 2, 2);

      if (fTrackingDataDebugOutput){

	printf("###DOTDB EV%04llu  GTU TRACK - S%02d-%d  pt: %+7.2f  pid: %3d  lm: 0x%02x %s\n",
	       fEventId,
	       iStack/5, iStack%5, fTrackingData->GetTrackPt(iStack, iTrk), fTrackingData->GetTrackPID(iStack, iTrk),
	       layerMask, trklStr.Data());
      }

      // paranoia checks
      for (Short_t iLayer = 5; iLayer >= 0; --iLayer){
	if (((layerMask >> iLayer) & 1) && (fTrackingData->GetTrackTrackletWord(iStack, iTrk, iLayer) == 0x10001000))
	  LogError("invalid layer mask / tracklet value combination A", "");

	if ((((layerMask >> iLayer) & 1) == 0) && (fTrackingData->GetTrackTrackletWord(iStack, iTrk, iLayer) != 0x10001000))
	  LogError("invalid layer mask / tracklet value combination B", "");
      }

    } // loop over tracks in stack
  } // loop over stacks

}


int AliHLTTRDMonitorComponent::ProcessTRDData(){

  UInt_t numTracklets;
  UInt_t numTracks;

  if (fHistoMode == 0){
    UInt_t numHists = fHistArray->GetEntries();
    for (UInt_t iHist = 0; iHist < numHists; ++iHist)
      if (fHistArray->At(iHist))
	((TH1*) (fHistArray->At(iHist)))->Reset();
  }

  // tracklets
  for (UInt_t iDet = 0; iDet < fkTRDChambers; ++iDet){
    numTracklets = fTrackingData->GetNumTracklets(iDet);
    for (UInt_t iTrkl = 0; iTrkl < numTracklets; ++iTrkl){
      fHistTrackletY->Fill(fTrackingData->GetTrackletLocalY(iDet, iTrkl));
      fHistTrackletDy->Fill(fTrackingData->GetTrackletBinDy(iDet, iTrkl));
      fHistTrackletYDy->Fill(fTrackingData->GetTrackletBinY(iDet, iTrkl), fTrackingData->GetTrackletBinDy(iDet, iTrkl));
      fHistTrackletZ->Fill(fTrackingData->GetTrackletBinZ(iDet, iTrkl));
      fHistTrackletPID->Fill(fTrackingData->GetTrackletPID(iDet, iTrkl));
      Int_t hc = fTrackingData->GetTrackletHCId(iDet, iTrkl);
      fHistTrackletsHCId->Fill(hc/60, hc%60);
    }
  }

  for (UShort_t iStack = 0; iStack < fkTRDStacks; ++iStack){
    // timing
    fHistTrackletTimingStack->Fill(iStack, fTrackingData->GetTrackletStartTime(iStack/5, iStack%5));
    fHistTrackletTimingStack->Fill(iStack, fTrackingData->GetTrackletEndTime(iStack/5, iStack%5));
    fHistTrackingTiming->Fill(0., fTrackingData->GetTrackletStartTime(iStack/5, iStack%5));
    fHistTrackingTiming->Fill(1., fTrackingData->GetTrackletEndTime(iStack/5, iStack%5));
    fHistTrackingTiming->Fill(2., fTrackingData->GetTMUTrackingDoneTime(iStack/5, iStack%5));

    // GTU tracks
    numTracks = fTrackingData->GetNumTracks(iStack);
    fHistTracksStack->Fill(iStack/5, iStack%5, numTracks);
    for (UInt_t iTrk = 0; iTrk < numTracks; ++iTrk){
      Double_t gpt = fTrackingData->GetTrackPt(iStack, iTrk);
      fHistTrackPt->Fill(gpt);
      fHistTrackPID->Fill(fTrackingData->GetTrackPID(iStack, iTrk));
      fHistTrackLayers->Fill(fTrackingData->GetTrackLayerNum(iStack, iTrk));
      if (gpt >= fTrackHighPtThreshold)
	fHistTrackLayersHighPt->Fill(fTrackingData->GetTrackLayerNum(iStack, iTrk));
    } // loop over tracks in stack
  } // loop over stacks

  for (UShort_t iSector = 0; iSector < fkTRDSectors; ++iSector){
    fHistTrackingTiming->Fill(3., fTrackingData->GetSMUTrackingDoneTime(iSector), fkTRDStacksPerSector);

    UInt_t sectorTrgFlags = fTrackingData->GetSectorTrgContribs(iSector);
    for (UShort_t iTrgCtb = 0; iTrgCtb < 12; ++iTrgCtb)
      if ((sectorTrgFlags >> iTrgCtb) & 1)
	fHistTriggerContribs->Fill(iSector, iTrgCtb);
  }

  return kTRUE;
}

int AliHLTTRDMonitorComponent::DoEvent(const AliHLTComponentEventData& evtData,
				       AliHLTComponentTriggerData& /*trigData*/) {

  int iResult = 0;
  fEventId = evtData.fEventID;
  fHistEventTypes->Fill(0.);

  LogDebug("### START DoEvent [event id: %llu, %d blocks, size: %d]",
	   fEventId, evtData.fBlockCnt, evtData.fStructSize);

  if (!IsDataEvent()) {  // process data events only
    LogDebug("### END   DoEvent [event id: %llu, %d blocks, size: %d] (skipped: no data event)",
	     fEventId, evtData.fBlockCnt, evtData.fStructSize);
    return iResult;
  }

  fHistEventTypes->Fill(1.);

  fTrackingData->SetLogPrefix(Form("TRDHLTGM XXXXXX-%05llu: [MON] {TrkDat} ", fEventId));

  do {

    // access to TRD specific data from AliHLTTRDPreprocessorComponent
    if (!PrepareTRDData()){
      LogError("access to TRD data failed. Skipping event...", "");
      break;
    }

    if (fTrackingData->GetNumTracks() + fTrackingData->GetNumTracklets() == 0) {
      LogDebug("no TRD-relevant information, skipping further event processing");
      break;
    }

    fHistEventTypes->Fill(2.);

    // DumpTrackingData();

    if (!ProcessTRDData()) {
      LogError("processing of TRD data failed, skipping further event processing");
      break;
    }

    break;

  } while (1);

  PushBack(fHistArray, (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTRD), 0x3fffff);

  LogDebug("### END   DoEvent [event id: %llu, %d blocks, size: %d]",
		      fEventId, evtData.fBlockCnt, evtData.fStructSize);

  return iResult;
}
