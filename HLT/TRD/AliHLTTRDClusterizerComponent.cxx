// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTRDClusterizerComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  A TRDClusterizer processing component for the HLT. */

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
#include "AliHLTTRDCluster.h"

#include "AliGeomManager.h"
#include "AliTRDReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliHLTTRDClusterizer.h"
#include "AliTRDrecoParam.h"
#include "AliTRDrawStreamBase.h"
#include "AliTRDcluster.h"

#include "AliRawReaderMemory.h"

#ifdef HAVE_VALGRIND_CALLGRIND_H
#include <valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION do { } while (0)
#define CALLGRIND_STOP_INSTRUMENTATION do { } while (0)
#endif

#include <cstdlib>
#include <cerrno>
#include <string>

ClassImp(AliHLTTRDClusterizerComponent)
   
AliHLTTRDClusterizerComponent::AliHLTTRDClusterizerComponent():
  AliHLTProcessor(),
  fOutputPercentage(500),
  fOutputConst(0),
  fClusterizer(NULL),
  fRecoParam(NULL),
  fMemReader(NULL),
  fReconstructor(NULL)
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
  list.push_back( (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTRD) );
}

AliHLTComponent_DataType AliHLTTRDClusterizerComponent::GetOutputDataType()
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
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

AliHLTComponent* AliHLTTRDClusterizerComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDClusterizerComponent;
};

int AliHLTTRDClusterizerComponent::DoInit( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  Int_t iRawDataVersion = 2;
  int i = 0;
  char* cpErr;

  Int_t iRecoParamType = -1; // default will be the low flux

  // the data type will become obsolete as soon as the formats are established
  Int_t iRecoDataType = -1; // default will be simulation
  Int_t iyPosMethod = 1;     // 0=COG 1=LUT 2=Gauss 
  Bool_t bProcessTracklets = kFALSE;
  string geometryFileName = "";
  
  while ( i < argc )
    {
      HLTDebug("argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_percentage" ) )
	{
	  if ( i+1>=argc )
	    {
	      HLTError("Missing output_percentage parameter");
	      return ENOTSUP;
	    }
	  HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
	  fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      HLTError("Cannot convert output_percentage parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }
	  HLTInfo("Output percentage set to %i %%", fOutputPercentage );
	  i += 2;
	}
      else if ( strcmp( argv[i], "-lowflux" ) == 0)
	{
	  iRecoParamType = 0;	  
	  HLTDebug("Low flux reco selected.");
	  i++;
	}
      else if ( strcmp( argv[i], "-highflux" ) == 0)
	{
	  iRecoParamType = 1;	  
	  HLTDebug("High flux reco selected.");
	  i++;
	}
      else if ( strcmp( argv[i], "-cosmics" ) == 0)
	{
	  iRecoParamType = 2;	  
	  HLTDebug("Cosmic test reco selected.");
	  i++;
	}
      // raw data type - sim or experiment
      else if ( strcmp( argv[i], "-simulation" ) == 0)
	{
	  iRecoDataType = 0;
	  i++;
	}
      else if ( strcmp( argv[i], "-experiment" ) == 0)
	{
	  iRecoDataType = 1;
	  i++;
	}
      else if ( strcmp( argv[i], "-rawver" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError("Missing -rawver argument");
	      return ENOTSUP;	      
	    }
	  iRawDataVersion = atoi( argv[i+1] );
	  HLTInfo("Raw data version is %d", iRawDataVersion );	  
	  i += 2;
	  
	}      

      else if ( strcmp( argv[i], "-geometry" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError("Missing -geometry argument");
	      return ENOTSUP;	      
	    }
	  geometryFileName = argv[i+1];
	  HLTInfo("GeomFile storage is %s", geometryFileName.c_str() );	  
	  i += 2;
	} 
      else if ( strcmp( argv[i], "-processTracklets" ) == 0)
	{
	  bProcessTracklets = kTRUE;
 	  i++; 
 	}
      else if ( strcmp( argv[i], "-yPosMethod" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError("Missing -yPosMethod argument");
	      return ENOTSUP;	      
	    }
	  if( strcmp(argv[i], "COG") )
	    iyPosMethod=0;
	  else if( strcmp(argv[i], "LUT") )
	    iyPosMethod=1;
	  else if( strcmp(argv[i], "Gauss") )
	    iyPosMethod=2;
	  else {
	    HLTError("Unknown -yPosMethod argument");
	    return ENOTSUP;	      
	  }
	  i += 2;
	}
      else if ( strcmp( argv[i], "-noZS" ) == 0) //no zero surpression in the input data
	{
	  fOutputPercentage = 100;
 	  i++; 
 	}
      
      else{
	HLTError("Unknown option '%s'", argv[i] );
	return EINVAL;
      }
      
    }

  // THE "REAL" INIT COMES HERE

if (iRecoParamType < 0 || iRecoParamType > 2)
    {
      HLTWarning("No reco param selected. Use -lowflux or -highflux flag. Defaulting to low flux.");
      iRecoParamType = 0;
    }

  if (iRecoParamType == 0)
    {
      fRecoParam = AliTRDrecoParam::GetLowFluxParam();
      HLTDebug("Low flux params init.");
    }

  if (iRecoParamType == 1)
    {
      fRecoParam = AliTRDrecoParam::GetHighFluxParam();
      HLTDebug("High flux params init.");
    }

  if (iRecoParamType == 2)
    {
      fRecoParam = AliTRDrecoParam::GetCosmicTestParam();
      HLTDebug("Cosmic Test params init.");
    }

  if (fRecoParam == 0)
    {
      HLTError("No reco params initialized. Sniffing big trouble!");
      return -1;
    }

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(fRecoParam);
  fReconstructor->SetStreamLevel(0, AliTRDReconstructor::kClusterizer); // default value
  HLTInfo("Not writing clusters. I.e. output is a TClonesArray of clusters");
  TString recoOptions="hlt,!cw,sl_cf_0";
  switch(iRecoDataType){
  case 0: recoOptions += ",tc"; break;
  case 1: recoOptions += ",!tc"; break;
  }
  switch(iyPosMethod){
  case 0: recoOptions += ",!gs,!lut"; break;
  case 1: recoOptions += ",!gs,lut"; break;
  case 2: recoOptions += ",gs,!lut"; break;
  }
  if(bProcessTracklets) recoOptions += ",tp";
  else  recoOptions += ",!tp";

  HLTInfo("Reconstructor options: %s",recoOptions.Data());
  fReconstructor->SetOption(recoOptions.Data());
  
  // init the raw data type to be used...
  // the switch here will become obsolete as soon as the data structures is fixed 
  // both: in sim and reality
  if (iRecoDataType < 0 || iRecoDataType > 1)
    {
      HLTWarning("No data type selected. Use -simulation or -experiment flag. Defaulting to simulation.");
      iRecoDataType = 0;
    }

  if (iRecoDataType == 0)
    {
      AliTRDrawStreamBase::SetRawStreamVersion(AliTRDrawStreamBase::kTRDsimStream);
      HLTDebug("Data type expected is SIMULATION!");
    }

  if (iRecoDataType == 1)
    {
      AliTRDrawStreamBase::SetRawStreamVersion(AliTRDrawStreamBase::kTRDrealStream);
      HLTDebug("Data type expected is EXPERIMENT!");
    }

  // the DATA BASE STUFF
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    HLTError("DefaultStorage is not Set in CDBManager");
    return -1;
  }
  if(AliCDBManager::Instance()->GetRun()<0){
    AliCDBManager *cdb = AliCDBManager::Instance();
    if (cdb)
      {
  	cdb->SetRun(0);
  	HLTWarning("Setting CDB Runnumber to 0. CDB instance 0x%x", cdb);
      }
    else
      {
  	HLTError("Could not get CDB instance", "cdb 0x%x", cdb);
  	return -1;
      }
  }
  HLTInfo("CDB default storage: %s; RunNo: %i", (AliCDBManager::Instance()->GetDefaultStorage()->GetBaseFolder()).Data(), AliCDBManager::Instance()->GetRun());
  
  if(!AliGeomManager::GetGeometry()){
    if(!TFile::Open(geometryFileName.c_str())){
      HLTInfo("Loading standard geometry file");
      AliGeomManager::LoadGeometry();
    }else{
      HLTWarning("Loading non-standard geometry file");
      AliGeomManager::LoadGeometry(geometryFileName.c_str());
    }
    if(!AliGeomManager::GetGeometry()){
      HLTError("Cannot load geometry");
      return EINVAL;
    }
  }
  else{
    HLTInfo("Geometry Already Loaded");
  }  
  
  fMemReader = new AliRawReaderMemory;

  fClusterizer = new AliHLTTRDClusterizer("TRDCclusterizer", "TRDCclusterizer");
  fClusterizer->SetReconstructor(fReconstructor);
  fClusterizer->SetUseLabels(kFALSE);
  fClusterizer->SetRawVersion(iRawDataVersion);

  if(fReconstructor->IsProcessingTracklets())
    fOutputConst = fClusterizer->GetTrMemBlockSize();
  return 0;
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

  if (evtData.fEventID == 1)
    CALLGRIND_START_INSTRUMENTATION;

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
      
      //      fMemReader->Reset();
      fMemReader->SetMemory((UChar_t*) block.fPtr, block.fSize);

      AliHLTUInt32_t spec = block.fSpecification;
      
      Int_t id = 1024;
      
      for ( Int_t ii = 0; ii < 18 ; ii++ ) {
	if ( spec & 0x00000001 ) {
	  id += ii;
	  break;
	}
	spec = spec >> 1 ;
      }

      fMemReader->SetEquipmentID( id ); 
      
      fClusterizer->SetMemBlock(outputPtr+offset);
      Bool_t iclustered = fClusterizer->Raw2ClustersChamber(fMemReader);
      if (iclustered == kTRUE)
	{
	  HLTDebug( "Clustered successfully");
	}
      else
	{
	  HLTError("Clustering ERROR");
	  return -1;
	}

      // put the tree into output
      //fcTree->Print();
      
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
	  //bd.fSpecification = spec;
	  bd.fSpecification = gkAliEventTypeData;
	  bd.fDataType = AliHLTTRDDefinitions::fgkMCMtrackletDataType;
	  outputBlocks.push_back( bd );
	  HLTDebug( "BD fPtr 0x%x, fOffset %i, size %i, dataType %s, spec 0x%x ", bd.fPtr, bd.fOffset, bd.fSize, DataType2Text(bd.fDataType).c_str(), spec);
	}
	offset = totalSize;
      }

      addedSize = fClusterizer->GetAddedClSize();
      if (addedSize > 0){
	// Using low-level interface 
	// with interface classes
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
	//bd.fSpecification = spec;
	bd.fSpecification = gkAliEventTypeData;
	bd.fDataType = AliHLTTRDDefinitions::fgkClusterDataType;
	outputBlocks.push_back( bd );
	HLTDebug( "BD fPtr 0x%x, fOffset %i, size %i, dataType %s, spec 0x%x ", bd.fPtr, bd.fOffset, bd.fSize, DataType2Text(bd.fDataType).c_str(), spec);
	offset = totalSize;
	      
      }
      else 
	HLTWarning("Array of clusters is empty!");
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
