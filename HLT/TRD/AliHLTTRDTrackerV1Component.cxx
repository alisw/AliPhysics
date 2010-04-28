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

/** @file   AliHLTTRDTrackerV1Component.cxx
    @author Theodor Rascanu
    @date   
    @brief  A TRDTrackerV1 processing component for the HLT.
*/

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTRDTrackerV1Component.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDTrack.h"
#include "AliHLTTRDUtils.h"

#include "TFile.h"
#include "TChain.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"

#include "AliTRDcalibDB.h"
#include "AliTRDReconstructor.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDrecoParam.h"

#include <cstdlib>
#include <cerrno>
#include <string>

ClassImp(AliHLTTRDTrackerV1Component)

void AliHLTTRDTrackerV1Component::AliHLTTRDESDEvent::CreateStdContent()
{
  TClonesArray* tracksArray = new TClonesArray("AliESDtrack",0);
  tracksArray->SetName(AliESDEvent::fgkESDListName[AliESDEvent::kTracks]);
  AddObject(tracksArray);
  GetStdContent();
}

void AliHLTTRDTrackerV1Component::AliHLTTRDESDEvent::Streamer(TBuffer &/*R__b*/)
{
  AliFatal("class is for internal us only and not for streaming");
}

AliHLTTRDTrackerV1Component::AliHLTTRDTrackerV1Component():
  AliHLTProcessor(),
  fOutputPercentage(100), // By default we copy to the output exactly what we got as input 
  fTracker(NULL),
  fRecoParam(NULL),
  fReconstructor(NULL),
  fESD(NULL),
  fClusterArray(NULL),
  fRecoParamType(-1),
  fNtimeBins(-1),
  fPIDmethod(1),
  fgeometryFileName(""),
  fHLTflag(kTRUE),
  fOutputV1Tracks(kTRUE),
  fHighLevelOutput(kFALSE),
  fEmulateHLTTracks(kFALSE),
  fImproveTracklets(kFALSE)
{
  // Default constructor

}

AliHLTTRDTrackerV1Component::~AliHLTTRDTrackerV1Component()
{
  // Destructor
}

const char* AliHLTTRDTrackerV1Component::GetComponentID()
{
  // Return the component ID const char *
  return "TRDTrackerV1"; // The ID of this component
}

void AliHLTTRDTrackerV1Component::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data  
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back(AliHLTTRDDefinitions::fgkClusterDataType);
}

AliHLTComponentDataType AliHLTTRDTrackerV1Component::GetOutputDataType()
{
  // Get the output data type
  return kAliHLTMultipleDataType;
}

int AliHLTTRDTrackerV1Component::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data types
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTrack | kAliHLTDataOriginTRD);
  tgtList.push_back(AliHLTTRDDefinitions::fgkTracksDataType);
  return tgtList.size();
}

void AliHLTTRDTrackerV1Component::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = fOutputV1Tracks ? 2*((double)fOutputPercentage)/100.0 : 0.5*((double)fOutputPercentage)/100.0;
}

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTRDTrackerV1Component::Spawn()
{
  return new AliHLTTRDTrackerV1Component;
};


int AliHLTTRDTrackerV1Component::DoInit( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  int iResult=0;

  fReconstructor = new AliTRDReconstructor();
  HLTDebug("TRDReconstructor at 0x%x", fReconstructor);
  fESD = new AliHLTTRDESDEvent();
  fESD->CreateStdContent();

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

  if(iResult<0) return iResult;

  fTracker = new AliTRDtrackerV1();
  HLTDebug("TRDTracker at 0x%x", fTracker);
  fTracker->SetReconstructor(fReconstructor);

  fClusterArray = new TClonesArray("AliTRDcluster"); // would be nice to allocate memory for all clusters here.

  return iResult;
}

int AliHLTTRDTrackerV1Component::DoDeinit()
{
  // Deinitialization of the component

  fTracker->SetClustersOwner(kFALSE);
  delete fTracker;
  fTracker = NULL;

  fClusterArray->Delete();
  delete fClusterArray;
  fClusterArray = NULL;
  
  // We need to set clusters in Reconstructor to null to prevent from 
  // double deleting, since we delete TClonesArray by ourself.
  fReconstructor->SetClusters(0x0);
  delete fReconstructor;
  fReconstructor = NULL;
  delete fESD;
  fESD = NULL;
  
  AliTRDcalibDB::Terminate();

  return 0;
}

int AliHLTTRDTrackerV1Component::DoEvent( const AliHLTComponentEventData& evtData, 
					  const AliHLTComponentBlockData* blocks, 
					  AliHLTComponent_TriggerData& /*trigData*/, 
					  AliHLTUInt8_t* outputPtr, 
					  AliHLTUInt32_t& size, 
					  vector<AliHLTComponent_BlockData>& outputBlocks )
{
  // Process an event

  HLTDebug("NofBlocks %i", evtData.fBlockCnt );
  
  AliHLTUInt32_t totalSize = 0, offset = 0;

  AliHLTComponentDataType expectedDataType = AliHLTTRDDefinitions::fgkClusterDataType;
  for ( unsigned long iBlock = 0; iBlock < evtData.fBlockCnt; iBlock++ ) 
    {
      const AliHLTComponentBlockData &block = blocks[iBlock];
      AliHLTComponentDataType inputDataType = block.fDataType;

      if(inputDataType != expectedDataType)
	{
	  HLTDebug( "Block # %i/%i; Event 0x%08LX (%Lu) Wrong received datatype: %s - Skipping",
		    iBlock, evtData.fBlockCnt-1,
		    evtData.fEventID, evtData.fEventID, 
		    DataType2Text(inputDataType).c_str());
	  continue;
	}
      else {
	HLTDebug("We get the right data type: Block # %i/%i; Event 0x%08LX (%Lu) Received datatype: %s; Block Size: %i",
		 iBlock, evtData.fBlockCnt-1,
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

      fESD->Reset();
      //fESD->SetMagneticField(GetBz());

      AliHLTTRDUtils::ReadClusters(fClusterArray, block.fPtr, block.fSize, &fNtimeBins);
      HLTDebug("Reading number of time bins from input block. Setting number of timebins to %d", fNtimeBins);
      AliTRDtrackerV1::SetNTimeBins(fNtimeBins);

      HLTDebug("TClonesArray of clusters: nbEntries = %i", fClusterArray->GetEntriesFast());
      fTracker->LoadClusters(fClusterArray);

      fTracker->Clusters2Tracks(fESD);

      Int_t nTracks = fESD->GetNumberOfTracks();
      HLTInfo("Number of tracks  == %d ==", nTracks);  

      TClonesArray* trdTracks;
      trdTracks = fTracker->GetListOfTracks();
      
      if(fHighLevelOutput){
	if(fEmulateHLTTracks && trdTracks){
	  // TClonesArray* oldArr = trdTracks;
	  trdTracks = new TClonesArray(*trdTracks);
	  AliHLTTRDUtils::EmulateHLTTracks(trdTracks);
	  // if(oldArr->At(0)){
	  //   HLTInfo("Old Track:");
	  //   ((AliTRDtrackV1*)oldArr->At(0))->Print("a");
	  //   HLTInfo("\nNew Track:");
	  //   ((AliTRDtrackV1*)trdTracks->At(0))->Print("a");
	  // }
	}

	TObjString strg;
	strg.String() += fNtimeBins;
	if(trdTracks)
	  PushBack(trdTracks, AliHLTTRDDefinitions::fgkHiLvlTracksDataType, block.fSpecification);
	else{
	  TClonesArray temp("AliTRDtrackV1");
	  PushBack(&temp, AliHLTTRDDefinitions::fgkHiLvlTracksDataType, block.fSpecification);
	}
	PushBack(&strg, AliHLTTRDDefinitions::fgkHiLvlTracksDataType, block.fSpecification);

	if(fEmulateHLTTracks && trdTracks){
	  trdTracks->Delete();
	  delete trdTracks;
	}
      }
      else if(nTracks>0){
	HLTDebug("We have an output ESDEvent: 0x%x with %i tracks", fESD, nTracks);
	AliHLTUInt32_t addedSize = AliHLTTRDUtils::AddESDToOutput(fESD, outputPtr+offset);
	totalSize += addedSize;
	  
	// Fill block 
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	//bd.fPtr = outputPtr;
	bd.fOffset = offset;
	bd.fSize = addedSize;
	bd.fSpecification = block.fSpecification;
	bd.fDataType = kAliHLTDataTypeTrack | kAliHLTDataOriginTRD;
	outputBlocks.push_back( bd );
	HLTDebug("BD ptr 0x%x, offset %i, size %i, datav1Type %s, spec 0x%x ", bd.fPtr, bd.fOffset, bd.fSize, DataType2Text(bd.fDataType).c_str(), bd.fSpecification);
	offset = totalSize;

	if (fOutputV1Tracks && trdTracks){
	  HLTDebug("We have an output array: pointer to trdTracks = 0x%x, nbEntries = %i", trdTracks, trdTracks->GetEntriesFast());
	  
	  addedSize = AliHLTTRDUtils::AddTracksToOutput(trdTracks, outputPtr+offset, fNtimeBins);
	  totalSize += addedSize;
	  
	  // Fill block 
	  FillBlockData( bd );
	  //bd.fPtr = outputPtr;
	  bd.fOffset = offset;
	  bd.fSize = addedSize;
	  bd.fSpecification = block.fSpecification;
	  bd.fDataType = AliHLTTRDDefinitions::fgkTracksDataType;
	  outputBlocks.push_back( bd );
	  HLTDebug("BD ptr 0x%x, offset %i, size %i, dataType %s, spec 0x%x ", bd.fPtr, bd.fOffset, bd.fSize, DataType2Text(bd.fDataType).c_str(), bd.fSpecification);
	  offset = totalSize;
	}
      }

      HLTDebug("totalSize: %i", totalSize);
      
//       if ( totalSize > allocSize )
// 	{
// 	  HLTError("Too much data; Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
// 	  totalSize, size );
// 	  return EMSGSIZE;
// 	}

      //here we are deleting clusters (but not the TClonesArray itself)
      fTracker->UnloadClusters();
      AliTRDReconstructor::SetClusters(0x0);
      fClusterArray->Delete();
      
    }
      
  size = totalSize;
  HLTDebug("Event is done. size written to the output is %i", size);
  return 0;
}

int AliHLTTRDTrackerV1Component::Configure(const char* arguments){
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
      else if (argument.CompareTo("-solenoidBz")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTWarning("argument -solenoidBz is deprecated, magnetic field set up globally (%f)", GetBz());
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
      else if (argument.CompareTo("-outputV1Tracks")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	TString toCompareTo=((TObjString*)pTokens->At(i))->GetString();
	if (toCompareTo.CompareTo("yes")==0){
	  HLTInfo("Setting OutputV1Tracks to: %s", toCompareTo.Data());
	  fOutputV1Tracks=kTRUE;
	}
	else if (toCompareTo.CompareTo("no")==0){
	  HLTInfo("Setting OutputV1Tracks to: %s", toCompareTo.Data());
	  fOutputV1Tracks=kFALSE;
	}
	else {
	  HLTError("unknown argument for OutputV1Tracks: %s", toCompareTo.Data());
	  iResult=-EINVAL;
	  break;
	}
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
	  fEmulateHLTTracks=kTRUE;
	}
	else if (toCompareTo.CompareTo("no")==0){
	  HLTInfo("Setting emulateHLToutput to: %s", toCompareTo.Data());
	  fEmulateHLTTracks=kFALSE;
	}
	else {
	  HLTError("unknown argument for emulateHLToutput: %s", toCompareTo.Data());
	  iResult=-EINVAL;
	  break;
	}
	continue;
      }
      else if (argument.CompareTo("-PIDmethod")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	TString toCompareTo=((TObjString*)pTokens->At(i))->GetString();
	if (toCompareTo.CompareTo("LH")==0){
	  HLTInfo("Setting PID method to: %s", toCompareTo.Data());
	  fPIDmethod=0;
	}
	else if (toCompareTo.CompareTo("NN")==0){
	  HLTInfo("Setting PID method to: %s", toCompareTo.Data());
	  fPIDmethod=1;
	}
	else if (toCompareTo.CompareTo("TM")==0){
	  HLTInfo("Setting PID method to: %s", toCompareTo.Data());
	  fPIDmethod=2;
	}
	else {
	  HLTError("unknown argument for PID method: %s", toCompareTo.Data());
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

int AliHLTTRDTrackerV1Component::SetParams()
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

  if(!fRecoParam)
    {
      HLTError("No reco params initialized. Sniffing big trouble!");
      return -EINVAL;
    }

  switch(fPIDmethod){
  case 0: fRecoParam->SetPIDNeuralNetwork(kFALSE); break;
  case 1: fRecoParam->SetPIDNeuralNetwork(kTRUE); break;
  case 2: fRecoParam->SetPIDNeuralNetwork(kFALSE); break;
  }

  fRecoParam->SetImproveTracklets(fImproveTracklets);

  fRecoParam->SetStreamLevel(AliTRDrecoParam::kTracker, 0);
  fReconstructor->SetRecoParam(fRecoParam);

  TString recoOptions="sa,!cw";
  
  if(fHLTflag)
    recoOptions += ",hlt";

  HLTDebug("Reconstructor options are: %s",recoOptions.Data());
  fReconstructor->SetOption(recoOptions.Data());

  return iResult;
}

int AliHLTTRDTrackerV1Component::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation

  int iResult=0;
  const char* path="HLT/ConfigTRD/TrackerV1Component";
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

int AliHLTTRDTrackerV1Component::ReadPreprocessorValues(const char* modules)
{
  // see header file for class documentation
  
  int iResult = 0;
  TString str(modules);
  if(str.Contains("HLT") || str.Contains("TRD") || str.Contains("GRP")){
  
  }  
  return iResult;
}
