// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors:                                                       *
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

/** @file   AliHLTTRDClusterizerComponent.cxx
    @author Theodor Rascanu
    @date   
    @brief  A TRDClusterizer processing component for the HLT. 
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "AliHLTTRDClusterizerComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDClusterizer.h"
#include "AliHLTTRDUtils.h"

#include "AliGeomManager.h"
#include "AliTRDReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliTRDrecoParam.h"
#include "AliTRDrawStreamBase.h"
#include "AliTRDcluster.h"

#include "AliRawReaderMemory.h"

#ifdef HAVE_VALGRIND_CALLGRIND_H
#include <valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION (void)0
#define CALLGRIND_STOP_INSTRUMENTATION (void)0
#endif

#include <cstdlib>
#include <cerrno>
#include <string>

ClassImp(AliHLTTRDClusterizerComponent)
   
AliHLTTRDClusterizerComponent::AliHLTTRDClusterizerComponent()
: AliHLTProcessor(),
  fOutputPercentage(100),
  fOutputConst(0),
  fClusterizer(NULL),
  fRecoParam(NULL),
  fMemReader(NULL),
  fReconstructor(NULL),
  fRecoParamType(-1),
  fRecoDataType(-1),
  fRawDataVersion(2),
  fyPosMethod(1),
  fgeometryFileName(""),
  fProcessTracklets(kFALSE),
  fHLTstreamer(kTRUE),
  fTC(kFALSE),
  fHLTflag(kTRUE),
  fHighLevelOutput(kFALSE),
  fEmulateHLTClusters(kFALSE)
{
  // Default constructor

}

AliHLTTRDClusterizerComponent::~AliHLTTRDClusterizerComponent()
{
  // Destructor
  // Work is Done in DoDeInit()
}


const char* AliHLTTRDClusterizerComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDClusterizer"; // The ID of this component
}

void AliHLTTRDClusterizerComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTRD);
}

AliHLTComponentDataType AliHLTTRDClusterizerComponent::GetOutputDataType()
{
  // Get the output data type
  return kAliHLTMultipleDataType;
}

int AliHLTTRDClusterizerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data type
  tgtList.clear();
  tgtList.push_back(AliHLTTRDDefinitions::fgkClusterDataType);
  tgtList.push_back(AliHLTTRDDefinitions::fgkMCMtrackletDataType);
  return tgtList.size();
}


void AliHLTTRDClusterizerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = fOutputConst;
  inputMultiplier = ((double)fOutputPercentage)*4/100.0;
}

AliHLTComponent* AliHLTTRDClusterizerComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDClusterizerComponent;
};

int AliHLTTRDClusterizerComponent::DoInit( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  int iResult=0;
  
  fReconstructor = new AliTRDReconstructor();
  HLTDebug("TRDReconstructor at 0x%x", fReconstructor);

  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }

  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

  if(!fClusterizer){
    HLTFatal("Clusterizer was not initialized!");
    return -1;
  }

  if(iResult<0) return iResult;

  fMemReader = new AliRawReaderMemory;
  fClusterizer->SetReconstructor(fReconstructor);
  fClusterizer->SetUseLabels(kFALSE);

  if(fReconstructor->IsProcessingTracklets())
    fOutputConst = fClusterizer->GetTrMemBlockSize();

  return iResult;
}

int AliHLTTRDClusterizerComponent::DoDeinit()
{
  // Deinitialization of the component
  delete fMemReader;
  fMemReader = 0;
  delete fClusterizer;
  fClusterizer = 0;
  
  fReconstructor->SetClusters(0x0);
  delete fReconstructor;
  fReconstructor = 0x0;
  return 0;

  if (fRecoParam)
    {
      HLTDebug("Deleting fRecoParam");
      delete fRecoParam;
      fRecoParam = 0;
    }
}

int AliHLTTRDClusterizerComponent::DoEvent( const AliHLTComponentEventData& evtData, 
					    const AliHLTComponentBlockData* blocks, 
					    AliHLTComponent_TriggerData& /*trigData*/, 
					    AliHLTUInt8_t* outputPtr, 
					    AliHLTUInt32_t& size, 
					    vector<AliHLTComponent_BlockData>& outputBlocks )
{
  // Process an event

  if (evtData.fEventID == 10)
    CALLGRIND_START_INSTRUMENTATION;

  if(!IsDataEvent())return 0;

  HLTDebug( "NofBlocks %i", evtData.fBlockCnt );
  // Process an event
  AliHLTUInt32_t totalSize = 0, offset = 0;

  //implement a usage of the following
  //   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
  //   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
  //   void *triggerData = trigData.fData;
  //HLTDebug( "Trigger data received. Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  // Loop over all input blocks in the event
  AliHLTComponentDataType expectedDataType = (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTRD);
  for ( UInt_t iBlock = 0; iBlock < evtData.fBlockCnt; iBlock++ )
    {      
      const AliHLTComponentBlockData &block = blocks[iBlock];
      // lets not use the internal TRD data types here : AliHLTTRDDefinitions::fgkDDLRawDataType
      // which is depreciated - we use HLT global defs instead
      //      if ( block.fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTRD) )
      AliHLTComponentDataType inputDataType = block.fDataType;
      if ( inputDataType != expectedDataType)
	{
	  HLTDebug( "Block # %i/%i; Event 0x%08LX (%Lu) Wrong received datatype: %s - required datatype: %s; Skipping",
		    iBlock, evtData.fBlockCnt,
		    evtData.fEventID, evtData.fEventID, 
		    DataType2Text(inputDataType).c_str(), 
		    DataType2Text(expectedDataType).c_str());
	  if(block.fDataType == kAliHLTDataTypeEOR)
	    CALLGRIND_STOP_INSTRUMENTATION;
	  continue;
	}
      else 
	{
	  HLTDebug("We get the right data type: Block # %i/%i; Event 0x%08LX (%Lu) Received datatype: %s; Block Size: %i",
		   iBlock, evtData.fBlockCnt,
		   evtData.fEventID, evtData.fEventID, 
		   DataType2Text(inputDataType).c_str(),
		   block.fSize);
	}
      
#ifndef NDEBUG
      unsigned long constBase;
      double inputMultiplier;
      GetOutputDataSize(constBase,inputMultiplier);
      if(size<(constBase+block.fSize*inputMultiplier)){
	HLTWarning("Memory Block given might be too small: %i < %i; Event %Lu", size, constBase+block.fSize*inputMultiplier, evtData.fEventID);
      }
#endif

      // fMemReader->Reset();
      fMemReader->SetMemory((UChar_t*) block.fPtr, block.fSize);

      AliHLTUInt32_t spec = block.fSpecification;
      
      Int_t id = AliHLTTRDUtils::GetSM(spec) + 1024;

      fMemReader->SetEquipmentID(id);
      
      fClusterizer->SetMemBlock(outputPtr+offset);
      Bool_t bclustered = fClusterizer->Raw2ClustersChamber(fMemReader);
      if(bclustered)
	{
	  HLTDebug("Clustered successfully");
	}
      else
	{
	  HLTError("Clustering ERROR");
	  return -1;
	}

      AliHLTUInt32_t addedSize;
      if(fReconstructor->IsProcessingTracklets()){
	addedSize = fClusterizer->GetAddedTrSize();
	totalSize += fClusterizer->GetTrMemBlockSize();  //if IsProcessingTracklets() is enabled we always reserve a data block of size GetTrMemBlockSize() for the tracklets
	if (addedSize > 0){
	  // Using low-level interface 
	  // with interface classes
	  if ( totalSize > size )
	    {
	      HLTError("Too much data; Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
		       totalSize, size );
	      return EMSGSIZE;
	    }

	  // Fill block 
	  AliHLTComponentBlockData bd;
	  FillBlockData( bd );
	  bd.fOffset = offset;
	  bd.fSize = addedSize;
	  bd.fSpecification = block.fSpecification;
	  bd.fDataType = AliHLTTRDDefinitions::fgkMCMtrackletDataType;
	  outputBlocks.push_back( bd );
	  HLTDebug( "BD ptr 0x%x, offset %i, size %i, dataType %s, spec 0x%x ", bd.fPtr, bd.fOffset, bd.fSize, DataType2Text(bd.fDataType).c_str(), spec);
	}
	offset = totalSize;
      }

      addedSize = fClusterizer->GetAddedClSize();
      if (addedSize > 0){
	
	Int_t* nTimeBins = (Int_t*)(outputPtr+offset+fClusterizer->GetAddedClSize());
	*nTimeBins = fClusterizer->GetNTimeBins();
	addedSize += sizeof(*nTimeBins);

	totalSize += addedSize;
	if ( totalSize > size )
	  {
	    HLTError("Too much data; Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
		     totalSize, size );
	    return EMSGSIZE;
	  }

	// Fill block 
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset = offset;
	bd.fSize = addedSize;
	bd.fSpecification = block.fSpecification;
	bd.fDataType = AliHLTTRDDefinitions::fgkClusterDataType;
	outputBlocks.push_back( bd );
	HLTDebug( "BD ptr 0x%x, offset %i, size %i, dataType %s, spec 0x%x ", bd.fPtr, bd.fOffset, bd.fSize, DataType2Text(bd.fDataType).c_str(), spec);
	offset = totalSize;
      }
      else{
	HLTDebug("Array of clusters is empty!");
      }
    }
  fReconstructor->SetClusters(0x0);

  size = totalSize;
  HLTDebug("Event is done. size written to the output is %i", size);
  return 0;
}

void AliHLTTRDClusterizerComponent::PrintObject( TClonesArray* inClustersArray)
{
  AliTRDcluster* cluster=0x0;
  
  for (Int_t i=0; i < inClustersArray->GetEntriesFast(); i++){
    cluster = dynamic_cast<AliTRDcluster*>(inClustersArray->At(i));
    HLTDebug("cluster[%i]",i);
    HLTDebug("  PadCol = %i; PadRow = %i; PadTime = %i", cluster->GetPadCol(), cluster->GetPadRow(), cluster->GetPadTime());
    HLTDebug("  Detector = %i, Amplitude = %f, Center = %f", cluster->GetDetector(), cluster->GetQ(), cluster->GetCenter());
    HLTDebug("  LocalTimeBin =  %i; NPads = %i; maskedPosition: %s, status: %s", cluster->GetLocalTimeBin(), cluster->GetNPads(),cluster->GetPadMaskedPosition(),cluster->GetPadMaskedPosition());
  }
  
}

int AliHLTTRDClusterizerComponent::Configure(const char* arguments){
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
      if (argument.CompareTo("output_percentage")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Setting output percentage to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fOutputPercentage=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 
      else if (argument.CompareTo("-geometry")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Setting geometry to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fgeometryFileName=((TObjString*)pTokens->At(i))->GetString();
	continue;
      } 
      else if (argument.CompareTo("-lowflux")==0) {
	fRecoParamType = 0;
	HLTInfo("Low flux reconstruction selected");
	continue;
      }
      else if (argument.CompareTo("-highflux")==0) {
	fRecoParamType = 1;
	HLTInfo("High flux reconstruction selected");
	continue;
      }
      else if (argument.CompareTo("-cosmics")==0) {
	fRecoParamType = 2;
	HLTInfo("Cosmics reconstruction selected");
	continue;
      }
      else if (argument.CompareTo("-simulation")==0) {
	fRecoDataType = 0;
	HLTInfo("Awaiting simulated data");
	continue;
      }
      else if (argument.CompareTo("-experiment")==0) {
	fRecoDataType = 1;
	HLTInfo("Awaiting real data");
	continue;
      }
      else if (argument.CompareTo("-processTracklets")==0) {
	fProcessTracklets = kTRUE;
	HLTInfo("Writing L1 tracklets to output");
	continue;
      }
      else if (argument.CompareTo("-noZS")==0) {
	fOutputPercentage = 10;
	HLTInfo("Awaiting non zero surpressed data");
	continue;
      }
      else if (argument.CompareTo("-HLTflag")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	TString toCompareTo=((TObjString*)pTokens->At(i))->GetString();
	if (toCompareTo.CompareTo("yes")==0){
	  HLTInfo("Setting HLTflag to: %s", toCompareTo.Data());
	  fHLTflag=kTRUE;
	}
	else if (toCompareTo.CompareTo("no")==0){
	  HLTInfo("Setting HLTflag to: %s", toCompareTo.Data());
	  fHLTflag=kFALSE;
	}
	else {
	  HLTError("unknown argument for HLTflag: %s", toCompareTo.Data());
	  iResult=-EINVAL;
	  break;
	}
	continue;
      }
      else if (argument.CompareTo("-faststreamer")==0) {
	fHLTstreamer = kTRUE;
	HLTInfo("Useing fast raw streamer");
	continue;
      }
      else if (argument.CompareTo("-nofaststreamer")==0) {
	fHLTstreamer = kFALSE;
	HLTInfo("Don't use fast raw streamer");
	continue;
      }
      else if (argument.CompareTo("-tailcancellation")==0) {
	fTC = kTRUE;
	HLTInfo("Useing tailcancellation");
	continue;
      }
      else if (argument.CompareTo("-rawver")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Raw data version is: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fRawDataVersion=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 
      else if (argument.CompareTo("-highLevelOutput")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	TString toCompareTo=((TObjString*)pTokens->At(i))->GetString();
	if (toCompareTo.CompareTo("yes")==0){
	  HLTWarning("Setting highLevelOutput to: %s", toCompareTo.Data());
	  fHighLevelOutput=kTRUE;
	}
	else if (toCompareTo.CompareTo("no")==0){
	  HLTInfo("Setting highLevelOutput to: %s", toCompareTo.Data());
	  fHighLevelOutput=kFALSE;
	}
	else {
	  HLTError("unknown argument for highLevelOutput: %s", toCompareTo.Data());
	  iResult=-EINVAL;
	  break;
	}
	continue;
      }
      else if (argument.CompareTo("-emulateHLToutput")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	TString toCompareTo=((TObjString*)pTokens->At(i))->GetString();
	if (toCompareTo.CompareTo("yes")==0){
	  HLTWarning("Setting emulateHLToutput to: %s", toCompareTo.Data());
	  fEmulateHLTClusters=kTRUE;
	}
	else if (toCompareTo.CompareTo("no")==0){
	  HLTInfo("Setting emulateHLToutput to: %s", toCompareTo.Data());
	  fEmulateHLTClusters=kFALSE;
	}
	else {
	  HLTError("unknown argument for emulateHLToutput: %s", toCompareTo.Data());
	  iResult=-EINVAL;
	  break;
	}
	continue;
      }
      else if (argument.CompareTo("-yPosMethod")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	TString toCompareTo=((TObjString*)pTokens->At(i))->GetString();
	if (toCompareTo.CompareTo("COG")==0){
	  HLTInfo("Setting yPosMethod method to: %s", toCompareTo.Data());
	  fyPosMethod=0;
	}
	else if (toCompareTo.CompareTo("LUT")==0){
	  HLTInfo("Setting yPosMethod method to: %s", toCompareTo.Data());
	  fyPosMethod=1;
	}
	else if (toCompareTo.CompareTo("Gauss")==0){
	  HLTInfo("Setting yPosMethod method to: %s", toCompareTo.Data());
	  fyPosMethod=2;
	}
	else {
	  HLTError("unknown argument for yPosMethod: %s", toCompareTo.Data());
	  iResult=-EINVAL;
	  break;
	}
	continue;
      }
      
      else {
	HLTError("unknown argument: %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if(iResult>=0){
    iResult=SetParams();
  }
  return iResult;
}

int AliHLTTRDClusterizerComponent::SetParams()
{
  Int_t iResult=0;
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    HLTError("DefaultStorage is not set in CDBManager");
    return -EINVAL;
  }
  if(AliCDBManager::Instance()->GetRun()<0){
    HLTError("Run Number is not set in CDBManager");
    return -EINVAL;
  }
  HLTInfo("CDB default storage: %s; RunNo: %i", (AliCDBManager::Instance()->GetDefaultStorage()->GetBaseFolder()).Data(), AliCDBManager::Instance()->GetRun());

  if(!AliGeomManager::GetGeometry()){
    if(fgeometryFileName.CompareTo("")==0 || !TFile::Open(fgeometryFileName.Data())){
      HLTInfo("Loading standard geometry file");
      AliGeomManager::LoadGeometry();
    }else{
      HLTWarning("Loading NON-standard geometry file");
      AliGeomManager::LoadGeometry(fgeometryFileName.Data());
    }
    if(!AliGeomManager::GetGeometry()){
      HLTError("Could not load geometry");
      return -EINVAL;
    }
    HLTInfo("Applying Alignment from CDB object");
    AliGeomManager::ApplyAlignObjsFromCDB("TRD");
  }
  else{
    HLTInfo("Geometry Already Loaded!");
  }

  if(fReconstructor->GetRecoParam()){
    fRecoParam = new AliTRDrecoParam(*fReconstructor->GetRecoParam());
    HLTInfo("RecoParam already set!");
  }else{
    if(fRecoParamType == 0){
      HLTDebug("Low flux params init.");
      fRecoParam = AliTRDrecoParam::GetLowFluxParam();
    }
    if(fRecoParamType == 1){
      HLTDebug("High flux params init.");
      fRecoParam = AliTRDrecoParam::GetHighFluxParam();
    }
    if(fRecoParamType == 2){
      HLTDebug("Cosmic Test params init.");
      fRecoParam = AliTRDrecoParam::GetCosmicTestParam();
    }
  }

  if (!fRecoParam)
    {
      HLTError("No reco params initialized. Sniffing big trouble!");
      return -EINVAL;
    }

  if(fTC){fRecoParam->SetTailCancelation(kTRUE); HLTDebug("Enableing Tail Cancelation"); }
  else{fRecoParam->SetTailCancelation(kFALSE); HLTDebug("Disableing Tail Cancelation"); }

  switch(fyPosMethod){
  case 0: fRecoParam->SetGAUS(kFALSE); fRecoParam->SetLUT(kFALSE); break;
  case 1: fRecoParam->SetGAUS(kFALSE); fRecoParam->SetLUT(kTRUE); break;
  case 2: fRecoParam->SetGAUS(kTRUE); fRecoParam->SetLUT(kFALSE); break;
  }

  fRecoParam->SetStreamLevel(AliTRDrecoParam::kClusterizer, 0);
  fReconstructor->SetRecoParam(fRecoParam);

  TString recoOptions="!cw";
  if(fHLTflag)
    recoOptions += ",hlt";
  if(fProcessTracklets) recoOptions += ",tp";
  else  recoOptions += ",!tp";

  HLTDebug("Reconstructor options are: %s",recoOptions.Data());
  fReconstructor->SetOption(recoOptions.Data());

  if (fRecoDataType < 0 || fRecoDataType > 1)
    {
      HLTWarning("No data type selected. Use -simulation or -experiment flag. Defaulting to simulation.");
      fRecoDataType = 0;
    }

  if (fRecoDataType == 0)
    {
      AliTRDrawStreamBase::SetRawStreamVersion(AliTRDrawStreamBase::kTRDsimStream);
      HLTDebug("Data type expected is SIMULATION!");
    }

  if (fRecoDataType == 1)
    {
      AliTRDrawStreamBase::SetRawStreamVersion(AliTRDrawStreamBase::kTRDrealStream);
      HLTDebug("Data type expected is EXPERIMENT!");
    }

#ifndef HAVE_NOT_ALITRD_RAWSTREAM_r39608
  if(fHLTstreamer){
    AliTRDrawStreamBase::SetRawStreamVersion("default");
    HLTDebug("fast rawstreamer used");
  }else{
    AliTRDrawStreamBase::SetRawStreamVersion("FAST");
    HLTDebug("old rawstreamer used");
  }
#else
  if(fHLTstreamer){
      AliTRDrawStreamBase::SetRawStreamVersion("FAST");
      HLTDebug("fast rawstreamer used");  
    }
#endif

  if(!fClusterizer){
    fClusterizer = new AliHLTTRDClusterizer("TRDCclusterizer", "TRDCclusterizer");  
    HLTDebug("TRDClusterizer at 0x%x", fClusterizer);
  }

  fClusterizer->SetRawVersion(fRawDataVersion);

  return iResult;
}

int AliHLTTRDClusterizerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation

  int iResult=0;
  const char* path="HLT/ConfigTRD/ClusterizerComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
  	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
  	iResult=Configure(pString->GetString().Data());
      } else {
  	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("cannot fetch object \"%s\" from CDB", path);
    }
  }

  return iResult;
}

void AliHLTTRDClusterizerComponent::GetOCDBObjectDescription(TMap* const targetMap){
  // Get a list of OCDB object description needed for the particular component
  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigTRD/ClusterizerComponent"), new TObjString("component arguments"));
  targetMap->Add(new TObjString("TRD/Calib/ChamberGainFactor"), new TObjString("gain factor of chambers"));
  targetMap->Add(new TObjString("TRD/Calib/ChamberT0"), new TObjString("T0 of chambers"));
  targetMap->Add(new TObjString("TRD/Calib/ChamberVdrift"), new TObjString("drift velocity of chambers"));
  targetMap->Add(new TObjString("TRD/Calib/DetNoise"), new TObjString("noise of chambers"));
  targetMap->Add(new TObjString("TRD/Calib/LocalGainFactor"), new TObjString("per pad gain factor"));
  targetMap->Add(new TObjString("TRD/Calib/LocalT0"), new TObjString("per pad T0"));
  targetMap->Add(new TObjString("TRD/Calib/LocalVdrift"), new TObjString("per pad drift velocity"));
  targetMap->Add(new TObjString("TRD/Calib/PadNoise"), new TObjString("per pad noise"));
  targetMap->Add(new TObjString("TRD/Calib/PadStatus"), new TObjString("pad status"));
  targetMap->Add(new TObjString("TRD/Calib/PRFWidth"), new TObjString("pad response function"));
}
