// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTRDPreprocessorComponent.cxx
/// @author Felix Rettig, Stefan Kirsch
/// @date   2012-08-16
/// @brief  A pre-processing component for TRD tracking/trigger data on FEP-level
/// @ingroup alihlt_trd_components

#include <cstdlib>
#include "TList.h"
#include "AliLog.h"
#include "AliHLTDataTypes.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDPreprocessorComponent.h"
#include "AliTRDdigitsManager.h"
#include "AliRawReaderMemory.h"
#include "AliTRDrawStream.h"
#include "AliTRDtrackletWord.h"
#include "AliESDTrdTracklet.h"
#include "AliESDTrdTrack.h"
#include "AliTRDonlineTrackingDataContainer.h"

ClassImp(AliHLTTRDPreprocessorComponent)

#define LogError( ... ) { HLTError(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("ERROR", __VA_ARGS__); } }
#define LogInfo( ... ) { HLTInfo(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INFO", __VA_ARGS__); } }
#define LogInspect( ... ) { HLTDebug(__VA_ARGS__); if (fDebugLevel >= 1) { DbgLog("INSPECT", __VA_ARGS__); } }
#define LogDebug( ... ) { if (fDebugLevel >= 1) { HLTInfo(__VA_ARGS__); DbgLog("DEBUG", __VA_ARGS__); } }

AliHLTTRDPreprocessorComponent::AliHLTTRDPreprocessorComponent() :
  AliHLTProcessor(),
  fDebugLevel(0),
  fEventId(fgkInvalidEventId),
  fTrackletArray(NULL),
  fGtuTrackArray(NULL),
  fRawReaderMem(NULL),
  fDigitsManagerTrd(NULL),
  fRawReaderTrd(NULL),
  fTrackingData(NULL)
{
  // constructor
}

AliHLTTRDPreprocessorComponent::~AliHLTTRDPreprocessorComponent() {
  // destructor
}

const char* AliHLTTRDPreprocessorComponent::GetComponentID() {
  return "TRDPreprocessorComponent";
}

void AliHLTTRDPreprocessorComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTRD);
}

AliHLTComponentDataType AliHLTTRDPreprocessorComponent::GetOutputDataType() {
  return kAliHLTMultipleDataType;
}

int AliHLTTRDPreprocessorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) {
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTObject | kAliHLTDataOriginTRD);
  return tgtList.size();
}

void AliHLTTRDPreprocessorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier ) {
  constBase = 5000000;
  inputMultiplier = 0;
}

void AliHLTTRDPreprocessorComponent::GetOCDBObjectDescription( TMap* const /*targetMap*/) {
}

AliHLTComponent* AliHLTTRDPreprocessorComponent::Spawn(){
  return new AliHLTTRDPreprocessorComponent;
}

int AliHLTTRDPreprocessorComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) {
  return 0;
}

int AliHLTTRDPreprocessorComponent::ReadPreprocessorValues(const char* /*modules*/){
  return 0;
}

int AliHLTTRDPreprocessorComponent::ScanConfigurationArgument(int argc, const char** argv){

  if (argc <= 0)
    return 0;

  UShort_t iArg = 0;
  TString argument(argv[iArg]);

  if (!argument.CompareTo("-debug")){
    if (++iArg >= argc) return -EPROTO;
    argument = argv[iArg];
    fDebugLevel = argument.Atoi();
    LogInfo("debug level set to %d.", fDebugLevel);
    return 2;
  }

  return 0;
}

int AliHLTTRDPreprocessorComponent::DoInit(int argc, const char** argv){

  int iResult = 0;

  do {

    fRawReaderMem = new AliRawReaderMemory;
    if (!fRawReaderMem) {
      iResult=-ENOMEM;
      break;
    }

    fTrackletArray = new TClonesArray("AliTRDtrackletWord", 500);
    if (!fTrackletArray) {
      iResult=-ENOMEM;
      break;
    }
    fTrackletArray->BypassStreamer(kTRUE);  //## needed for performance improvement?

    fGtuTrackArray = new TClonesArray("AliESDTrdTrack", 50);
    if (!fGtuTrackArray){
      iResult=-ENOMEM;
      break;
    }

    fRawReaderTrd = new AliTRDrawStream(fRawReaderMem);
    if (!fRawReaderTrd) {
      iResult=-ENOMEM;
      break;
    }

    if (kFALSE){ // do not create digits manager for reader -> do not process ADC raw data
      fDigitsManagerTrd = new AliTRDdigitsManager();
      if (!fDigitsManagerTrd) {
	iResult=-ENOMEM;
	break;
      }
      fDigitsManagerTrd->CreateArrays();
      fRawReaderTrd->SetDigitsManager(fDigitsManagerTrd);
    }

    fRawReaderTrd->SetDigitsManager(fDigitsManagerTrd);
    fRawReaderTrd->SetTrackletArray(fTrackletArray);
    fRawReaderTrd->SetTrackArray(fGtuTrackArray);

    // Disable raw reader error messages that could flood HLT logbook
    AliLog::SetClassDebugLevel("AliTRDrawStream", 0);
    fRawReaderTrd->SetErrorDebugLevel(AliTRDrawStream::kLinkMonitor, 1);

    fTrackingData = new AliTRDonlineTrackingDataContainer();
    if (!fTrackingData){
      iResult=-ENOMEM;
      break;
    }

  } while (0);

  if (iResult < 0) {

    if (fTrackingData) delete fTrackingData;
    fTrackingData = NULL;

    if (fRawReaderTrd) delete fRawReaderTrd;
    fRawReaderTrd = NULL;

    if (fRawReaderMem) delete fRawReaderMem;
    fRawReaderMem = NULL;

    if (fDigitsManagerTrd) delete fDigitsManagerTrd;
    fDigitsManagerTrd = NULL;

    if (fGtuTrackArray) delete fGtuTrackArray;
    fGtuTrackArray = NULL;

    if (fTrackletArray) delete fTrackletArray;
    fTrackletArray = NULL;

  }

  vector<const char*> remainingArgs;
  for (int i = 0; i < argc; ++i)
    remainingArgs.push_back(argv[i]);

  if (argc > 0)
    ConfigureFromArgumentString(remainingArgs.size(), &(remainingArgs[0]));

  return iResult;
}

int AliHLTTRDPreprocessorComponent::DoDeinit() {

  if (fTrackingData) delete fTrackingData;
  fTrackingData = NULL;

  if (fRawReaderTrd) delete fRawReaderTrd;
  fRawReaderTrd = NULL;

  if (fRawReaderMem) delete fRawReaderMem;
  fRawReaderMem = NULL;

  if (fGtuTrackArray) delete fGtuTrackArray;
  fGtuTrackArray = NULL;

  if (fTrackletArray) delete fTrackletArray;
  fTrackletArray = NULL;

  return 0;
}

//void AliHLTTRDPreprocessorComponent::DbgLog(const char* prefix, const char* msg){
//  AliHLTEventID_t eventNumber = fEventId;
//  int runNumber = -1;
//  printf("TRDGM %s-%s: [PRE] %s%s\n",
// 	 (runNumber >= 0) ? Form("%06d", runNumber) : "XXXXXX",
// 	 (eventNumber != fgkInvalidEventId) ? Form("%05llu", eventNumber) : "XXXXX",
// 	 (strlen(prefix) > 0) ? Form("<%s> ", prefix) : "", msg);
//}


void AliHLTTRDPreprocessorComponent::DbgLog(const char* prefix, ...){
#ifdef __TRDHLTDEBUG
  AliHLTEventID_t eventNumber = fEventId;
  int runNumber = -1;
  printf("TRDHLTGM %s-X-%s: [PRE] %s",
 	 (runNumber >= 0) ? Form("%06d", runNumber) : "XXXXXX",
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


int AliHLTTRDPreprocessorComponent::DoEvent(const AliHLTComponentEventData& hltEventData,
					    AliHLTComponentTriggerData& /*trigData*/) {

  fEventId = hltEventData.fEventID;
  fTrackingData->SetLogPrefix(Form("TRDHLTGM XXXXXX-%05llu: [PRE] {TrkDat} ", hltEventData.fEventID));

  LogDebug("### START DoEvent [event id: %llu, %d blocks, size: %d]",
	   hltEventData.fEventID, hltEventData.fBlockCnt, hltEventData.fStructSize);

  // event processing function
  int iResult = 0;
  UInt_t sourceSectors = 0;

  fTrackingData->Clear();
  fTrackletArray->Clear();
  fGtuTrackArray->Clear();

  if (!IsDataEvent()) { // process data events only
    LogDebug("### END   DoEvent [event id: %llu, %d blocks, size: %d] (skipped: no data event)",
	     hltEventData.fEventID, hltEventData.fBlockCnt, hltEventData.fStructSize);
    return iResult;
  }

  TString infoStr("");

  // #FIXME: Also take care of SOR, EOR, etc...

  // loop over all incoming TRD raw data blocks
  for (const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTRD);
       pBlock != NULL && iResult >= 0;
       pBlock = GetNextInputBlock()) {

    int trdSector = -1;

    // determine sector from block specification
    for (unsigned pos = 0; pos < 8*sizeof(AliHLTUInt32_t); pos++) {
      if (pBlock->fSpecification & (0x1 << pos)) {
	if (trdSector >= 0) {
	  HLTWarning("Cannot uniquely identify DDL number from specification, skipping data block %s 0x%08x",
	      	     DataType2Text(pBlock->fDataType).c_str(),
	      	     pBlock->fSpecification);
	  trdSector = -1;
	  break;
	}
	trdSector = pos;
      }
    }
    if (trdSector < 0) continue;

    // add data block to rawreader
    infoStr += Form("%02d, ", trdSector);
    sourceSectors |= pBlock->fSpecification;
    if(!fRawReaderMem->AddBuffer((UChar_t*) pBlock->fPtr, pBlock->fSize, trdSector + 1024)){
      LogError("Could not add buffer of data block  %s, 0x%08x to rawreader",
	       DataType2Text(pBlock->fDataType).c_str(),
	       pBlock->fSpecification);
      continue;
    }

  } // loop over all incoming TRD raw data blocks

  if (sourceSectors){
    infoStr.Remove(infoStr.Length() - 2, 2);
    LogDebug("preprocessing data from sectors: %s...", infoStr.Data());
  } else {
    LogDebug("### END   DoEvent [event id: %llu, %d blocks, size: %d] (skipping: no TRD data)",
			hltEventData.fEventID, hltEventData.fBlockCnt, hltEventData.fStructSize);
    return iResult;
  }

  // extract header info, TRD tracklets and tracks from raw data
  fRawReaderTrd->ReadEvent();

  // read and store header info
  for (UShort_t iSector = 0; iSector < 18; ++iSector){
    if ((sourceSectors >> iSector) & 1){
      fTrackingData->SetSectorTrgWord(iSector, fRawReaderTrd->GetTriggerFlags(iSector));
      for (UShort_t iStack = 0; iStack < 5; ++iStack)
	fTrackingData->SetStackTrgWord(iSector, iStack, fRawReaderTrd->GetTrkFlags(iSector, iStack));
    }
  }

//  //## dbg only!!!!
//  for (UShort_t iSector = 0; iSector < 18; ++iSector){
//    fTrackingData->SetSectorTrgWord(iSector, 0x2345352 ^ (iSector + 34));
//    for (UShort_t iStack = 0; iStack < 5; ++iStack){
//      ULong64_t dwl = ((ULong64_t)0xaffe << 16) | (ULong64_t)iSector << 8 | iStack;
//      ULong64_t dw = (((ULong64_t)0xbead << 16) | (ULong64_t)iSector << 8 | iStack) | ((ULong64_t)dwl << 32);
//      fTrackingData->SetStackTrgWord(iSector, iStack, dw);
//    }
//  }

  // read and process TRD tracklets
  Int_t trackletIndex[1080] = { 0 };
  TList trklList;
  Int_t iTrkl = 0;
  trklList.SetOwner(kFALSE);
  AliTRDrawStream::SortTracklets(fTrackletArray, trklList, trackletIndex);
  TIter trackletIter(&trklList);
  while (AliTRDtrackletBase* tracklet = (AliTRDtrackletBase*) trackletIter()) {

    // tracklet data to tracking data for transport to other HLT components
    fTrackingData->AddTracklet(tracklet->GetHCId(), tracklet->GetTrackletWord());

//    // conversion to AliESDTrdTracklet only for testing
//    // label -2, we assume it's always raw tracklets here
//    AliESDTrdTracklet esdTracklet(tracklet->GetTrackletWord(), tracklet->GetHCId(), -2);
//    Int_t det = esdTracklet.GetDetector();
//    if (kFALSE){
//      printf("TRDPreproc: TRD tracklet %3d - S%02d-%d-%d (det %3d): 0x%08x  - y=%+5d  dy=%+3d  pid=%3d\n",
// 	     iTrkl, det/30, (det%30)/6, det%6, det,
// 	     esdTracklet.GetTrackletWord(), esdTracklet.GetBinY(), esdTracklet.GetBinDy(), esdTracklet.GetPID());
//    }

    ++iTrkl;
  }

  // read and process GTU tracks
  UInt_t numGtuTracks = fGtuTrackArray->GetEntriesFast();
  AliESDTrdTrack *trdTrack = NULL;
  UInt_t stack;
  Int_t refIndex[6];

  for (UInt_t iTrack = 0; iTrack < numGtuTracks; iTrack++) {
    trdTrack = (AliESDTrdTrack*) ((*fGtuTrackArray)[iTrack]);
    if (trdTrack){

      AliTRDrawStream::AssignTracklets(trdTrack, trackletIndex, refIndex);
      stack = trdTrack->GetStack();

      for (Int_t iLayer = 0; iLayer < 6; ++iLayer) {
	AliESDTrdTracklet *trkl = (refIndex[iLayer] > -1) ? (AliESDTrdTracklet*) trklList.At(refIndex[iLayer]) : 0x0;
	if (trkl)
	  trdTrack->AddTrackletReference(trkl, iLayer);
      }

      fTrackingData->AddTrack(trdTrack);

//      for (Int_t iLayer = 0; iLayer < 6; ++iLayer) {
// 	Int_t det = trdTrack->GetSector()*30 + stack*6 + iLayer;
//
// 	AliESDTrdTracklet *trkl = refIndex[iLayer] > -1 ? (AliESDTrdTracklet*) trklList.At(refIndex[iLayer]) : 0x0;
// 	if (trkl) {
// 	  AliDebugClass(5, Form("adding tracklet with index %i: 0x%08x",
// 				refIndex[iLayer], trkl->GetTrackletWord()));
// 	  if (trkl->GetDetector() != det)
// 	    AliErrorClass(Form("inconsistent assignment of tracklet 0x%08x in det %i to track in %i",
// 			       trkl->GetTrackletWord(), trkl->GetDetector(), det));
// 	  trdTrack->AddTrackletReference(trkl, iLayer);
// 	}
//      }
//
//      UInt_t pid = 0;
//      UShort_t lyrNum = 0;
//      TString trackletInfo("");
//      for (UShort_t iLayer = 0; iLayer < 6; ++iLayer){
// 	AliESDTrdTracklet* trkl = trdTrack->GetTracklet(5 - iLayer);
// 	if (trkl){
// 	  trackletInfo += Form("0x%08x  ", trkl->GetTrackletWord());
// 	  pid += trkl->GetPID();
// 	  lyrNum++;
// 	} else
// 	  trackletInfo += "----------  ";
//      }
//      trackletInfo.Remove(trackletInfo.Length() - 2, 2);
//      pid /= lyrNum;
//
//      Int_t pidDiff = trdTrack->GetPID() - pid;
//
//      UInt_t lm = trdTrack->GetLayerMask();
//      DbgLog("", Form("GTU track %d - S%02d-%d pt=%+.3f lm=(L5)%d%d%d%d%d%d(L0) [%s] [pid %s: %d|%d]",
// 		      iTrack, trdTrack->GetSector(), trdTrack->GetStack(), trdTrack->Pt(),
// 		      (lm >> 5) & 0x1, (lm >> 4) & 0x1, (lm >> 3) & 0x1,
// 		      (lm >> 2) & 0x1, (lm >> 1) & 0x1, (lm >> 0) & 0x1,
// 		      trackletInfo.Data(),
// 		      (TMath::Abs(pidDiff) <= 1) ? "ok" : "mismatch",
// 		      trdTrack->GetPID(), pid));
//

    } else
      LogError("GTU track %d is NULL!\n", iTrack);
  } // loop over GTU tracks

  fRawReaderMem->ClearBuffers();

  if (sourceSectors){
    LogDebug("pushing data for sectors: 0x%05x", sourceSectors);
    // transport of TRD tracklets and GTU tracks via tracking data container
    void* buffer;
    UInt_t bufferSize;
    fTrackingData->PrintSummary("preproc component");
    fTrackingData->Compress(buffer, bufferSize);
    iResult += PushBack(buffer, bufferSize, AliHLTTRDDefinitions::fgkOnlineDataType, sourceSectors);
    free(buffer);
  }

  LogDebug("### END   DoEvent [event id: %llu, %d blocks, size: %d]",
		      hltEventData.fEventID, hltEventData.fBlockCnt, hltEventData.fStructSize);

  return iResult;
}
