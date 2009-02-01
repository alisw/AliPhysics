// $Id: AliHLTTRDTrackerV1Component.cxx 23618 2008-01-29 13:07:38Z hristov $

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

/** @file   AliHLTTRDTrackerV1Component.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  A TRDTrackerV1 processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTRDTrackerV1Component.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDCluster.h"
#include "AliHLTTRDTrack.h"

#include "TFile.h"
#include "TChain.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliMagF.h"
#include "AliESDfriend.h"

#include "AliTRDcalibDB.h"
#include "AliTRDReconstructor.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDrecoParam.h"

#include <cstdlib>
#include <cerrno>
#include <string>

#ifdef HAVE_VALGRIND_CALLGRIND_H
#include <valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION() do { } while (0)
#define CALLGRIND_STOP_INSTRUMENTATION() do { } while (0)
#endif

ClassImp(AliHLTTRDTrackerV1Component);
    
AliHLTTRDTrackerV1Component::AliHLTTRDTrackerV1Component():
  AliHLTProcessor(),
  fOutputPercentage(100), // By default we copy to the output exactly what we got as input  
  fStrorageDBpath("local://$ALICE_ROOT"),
  fCDB(NULL),
  fGeometryFileName(""),
  fUseHLTClusters(kFALSE),
  fUseHLTTracks(kFALSE),
  fTracker(NULL),
  fRecoParam(NULL),
  fReconstructor(NULL)
{
  // Default constructor

  fGeometryFileName = getenv("ALICE_ROOT");
  fGeometryFileName += "/HLT/TRD/geometry.root";
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
  list.push_back( AliHLTTRDDefinitions::fgkClusterDataType );
}

AliHLTComponent_DataType AliHLTTRDTrackerV1Component::GetOutputDataType()
{
  // Get the output data type
  return AliHLTTRDDefinitions::fgkClusterDataType;
}

void AliHLTTRDTrackerV1Component::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTRDTrackerV1Component::Spawn()
{
  return new AliHLTTRDTrackerV1Component;
};

/**
 * Convert AliTRDtrackV1 to AliHLTTRDTrack 
 * Add HLTTrack to the output, defined by pointer
 * Fill block desctiptors 
 * Return size of the added to ouput objects
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDTrackerV1Component::AddToOutput(TClonesArray* inTrackArray, AliHLTUInt8_t* output)
{
  cout << "\nWriting tracks to the Memory\n ============= \n";
  AliTRDtrackV1* track = 0;
  AliHLTUInt32_t addedSize = 0;
  AliHLTUInt8_t *iterPtr = output;
  AliHLTTRDTrack * outPtr = (AliHLTTRDTrack*)iterPtr;
  
  if (inTrackArray){
    Int_t nbTracks  = inTrackArray->GetEntries();
    for (Int_t iTrack = 0; iTrack<nbTracks; iTrack++){
      AliHLTUInt32_t trackSize=0;
      
      track = dynamic_cast<AliTRDtrackV1*>(inTrackArray->At(iTrack));
      //track->Print();
      
      AliHLTTRDTrack *hltTrack = new (outPtr) AliHLTTRDTrack(track);
      trackSize = hltTrack->GetSize();
      addedSize += trackSize;
      HLTDebug("addedSize %i, trackSize %i", addedSize, trackSize);
      
      iterPtr += trackSize;
      outPtr = (AliHLTTRDTrack*)iterPtr;
    }
  }
  return addedSize;
  
}

int AliHLTTRDTrackerV1Component::DoInit( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  fOutputPercentage = 100;
  int i = 0;
  char* cpErr;
  

  Int_t iRecoParamType = -1; // default will be the low flux
  Int_t iNtimeBins = -1;     // number of time bins for the tracker to use
  Int_t iMagneticField = -1; // magnetic field: 0==OFF and 1==ON
  Bool_t bHLTMode = kTRUE, bWriteClusters = kFALSE;
  
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
	  HLTInfo("Output percentage set to %lu %%", fOutputPercentage );
	  i += 2;
	  
	}

      else if ( !strcmp( argv[i], "-NTimeBins" ) )
	{
	  if ( i+1>=argc )
	    {
	      HLTError("Missing -NTimeBins parameter");
	      return ENOTSUP;
	    }
	  HLTDebug("Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	  iNtimeBins = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      HLTError("Wrong Argument. Cannot convert -NTimeBins parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }	  
	  i += 2;
	  
	}

      else if ( strcmp( argv[i], "-cdb" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError( "Missing -cdb argument");
	      return ENOTSUP;	      
	    }
	  fStrorageDBpath = argv[i+1];
	  HLTInfo("DB storage is %s", fStrorageDBpath.c_str() );	  
	  i += 2;
	  
	}      

      else if ( strcmp( argv[i], "-geometry" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError("Missing -geometry argument");
	      return ENOTSUP;	      
	    }
	  fGeometryFileName = argv[i+1];
	  HLTInfo("GeomFile storage is %s", 
		  fGeometryFileName.c_str() );	  
	  i += 2;
	  
	}      

      // the flux parametrizations
      else if ( strcmp( argv[i], "-lowflux" ) == 0)
	{
	  iRecoParamType = 0;	  
	  HLTDebug("Low flux reco selected.");
	  i++;
	  
	}      
      
      else if ( strcmp( argv[i], "-highflux" ) == 0)
	{
	  iRecoParamType = 1;	  
	  HLTDebug("Low flux reco selected.");
	  i++;
	}      

      else if ( strcmp( argv[i], "-cosmics" ) == 0)
	{
	  iRecoParamType = 2;	  
	  HLTDebug("Cosmic test reco selected.");
	  i++;
	}      

      else if ( strcmp( argv[i], "-magnetic_field_ON" ) == 0)
	{
	  iMagneticField = 1;
	  i++;
	}
      else if ( strcmp( argv[i], "-magnetic_field_OFF" ) == 0)
	{
	  iMagneticField = 0;
	  i++;
	}
      else if ( strcmp( argv[i], "-writeClusters" ) == 0)
	{
	  bWriteClusters = kTRUE;
	  HLTDebug("input clusters are expected to be in a TTree.");
	  i++;
	}
      else if ( strcmp( argv[i], "-offlineMode" ) == 0)
	{
	  bHLTMode=kFALSE;
	  HLTDebug("Using standard offline tracking.");
	  i++;
	}
      else if ( strcmp( argv[i], "-useHLTClusters" ) == 0)
	{
	  fUseHLTClusters = kTRUE;
	  i++;
	  HLTInfo("expecting AliHLTCluster as input");
	}
      else if ( strcmp( argv[i], "-useHLTTracks" ) == 0)
	{
	  fUseHLTTracks = kTRUE;
	  i++;
	  HLTInfo("Using AliHLTTrack to pass data further in the chain");
	}
      else {
	HLTError("Unknown option '%s'", argv[i] );
	return EINVAL;
      }
      
    }

  // THE "REAL" INIT COMES HERE
  // offline condition data base
  fCDB = AliCDBManager::Instance();
  if (!fCDB)
    {
      HLTError("Could not get CDB instance", "fCDB 0x%x", fCDB);
      return -1;
    }
  else
    {
      fCDB->SetRun(0); // THIS HAS TO BE RETRIEVED !!!
      fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
      HLTDebug("CDB instance", "fCDB 0x%x", fCDB);
    }

  // check if the N of time bins make sense
  if (iNtimeBins <= 0)
    {
      HLTError("Sorry. Tracker needs number of time bins. At the moment you have to provide it with -NTimeBins <value>. The simulation always had 24 and the real data 30. Take your pick. Make sure the information is correct. Ask offline to implement how to propagate this information into clusters/cluster tree.");
      return -1;
    }

  if (iNtimeBins < 24 || iNtimeBins > 30)
    {
      HLTWarning("The number of time bins seems to be strange = %d. But okay. Let's try it...", iNtimeBins);
    }

  HLTDebug("The number of time bins = %d.", iNtimeBins);
  AliTRDtrackerV1::SetNTimeBins(iNtimeBins);

  // !!!! THIS IS IMPORTANT
  // init alifield map - temporarly via parameter - should come from a DB or DCS ?
  // !!!! 
  if (iMagneticField < 0)
    {
      iMagneticField = 0;
      HLTWarning("No magnetic field switch stated. Use -magnetic_field_ON or -magnetic_field_OFF flag. Defaulting to OFF = NO MAGNETIC FIELD");
    }
  
  if (!TGeoGlobalMagField::Instance()->IsLocked()) {
    if (iMagneticField == 0)
      {
	// magnetic field OFF
	AliMagF* field = new AliMagF("Maps","Maps",2,0.,0., 10.,AliMagF::k5kGUniform);
	TGeoGlobalMagField::Instance()->SetField(field);
	HLTDebug("Magnetic field is OFF.");
      }
    
    if (iMagneticField == 1)
      {
	// magnetic field ON
	AliMagF* field = new AliMagF("Maps","Maps",2,1.,1., 10.,AliMagF::k5kG);
	TGeoGlobalMagField::Instance()->SetField(field);
	HLTDebug("Magnetic field is ON.");
      }
  }
  else {
    HLTError("Magnetic field is already set and locked, cannot redefine it." );
  }

  // reconstruction parameters
  if (iRecoParamType < 0 || iRecoParamType > 2)
    {
      HLTWarning("No reco param selected. Use -lowflux -highflux -cosmics flags. Defaulting to low flux.");
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
  //   fRecoParam->SetChi2Y(.1);
  //   fRecoParam->SetChi2Z(5.);
  fReconstructor->SetRecoParam(fRecoParam);
  // write clusters [cw] = true
  // track seeding (stand alone tracking) [sa] = true
  // PID method in reconstruction (NN) [nn] = true
  // write online tracklets [tw] = false
  // drift gas [ar] = false
  // sl_tr_0 = StreamLevel_task_Level
  //  fReconstructor->SetOption("sa,!cw,hlt,sl_tr_0");
  TString recoOptions="sa,sl_tr_0";
  
  if (bWriteClusters)
    {
      recoOptions += ",cw";
    } 
  else
    {
      recoOptions += ",!cw";
    }
  if (bHLTMode)
    recoOptions += ",hlt";
  
  fReconstructor->SetOption(recoOptions.Data());
  HLTDebug("Reconstructor options are: %s",recoOptions.Data());
  
  if((AliGeomManager::GetGeometry()) == NULL){
    
    if ( TFile::Open(fGeometryFileName.c_str())) {
      AliGeomManager::LoadGeometry(fGeometryFileName.c_str());
    }
    else {
      HLTError("Cannot load geometry from file %s",fGeometryFileName.c_str());
      return EINVAL;
    }
  }
  else
    HLTInfo("Geometry Already Loaded");
  
  // create the tracker
  fTracker = new AliTRDtrackerV1();
  fTracker->SetReconstructor(fReconstructor);
  HLTDebug("TRDTracker at 0x%x", fTracker);

  if (fTracker == 0)
    {
      HLTError("Unable to create the tracker!");
      return -1;
    }

  return 0;
}

int AliHLTTRDTrackerV1Component::DoDeinit()
{
  // Deinitialization of the component

  fTracker->SetClustersOwner(kFALSE);
  delete fTracker;
  fTracker = 0x0;
  
  // We need to set clusters in Reconstructor to null to prevent from 
  // double deleting, since we delete TClonesArray by ourself in DoEvent.
  fReconstructor->SetClusters(0x0);
  delete fReconstructor;
  fReconstructor = 0x0;
  
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
  Bool_t bWriteClusters = fReconstructor->IsWritingClusters();

  HLTDebug("NofBlocks %lu", evtData.fBlockCnt );
  
  AliHLTUInt32_t totalSize = 0, offset = 0;
  AliHLTUInt32_t dBlockSpecification = 0;

  vector<AliHLTComponent_DataType> expectedDataTypes;
  GetInputDataTypes(expectedDataTypes);
  if (evtData.fEventID == 1)
    CALLGRIND_START_INSTRUMENTATION();
  for ( unsigned long iBlock = 0; iBlock < evtData.fBlockCnt; iBlock++ ) 
    {
      const AliHLTComponentBlockData &block = blocks[iBlock];
      offset = totalSize;
      AliHLTComponentDataType inputDataType = block.fDataType;
      Bool_t correctDataType = kFALSE;
      
      for(UInt_t i = 0; i < expectedDataTypes.size(); i++)
	if( expectedDataTypes.at(i) == inputDataType)
	  correctDataType = kTRUE;
      if (!correctDataType)
	{
	  HLTDebug( "Block # %i/%i; Event 0x%08LX (%Lu) Wrong received datatype: %s - Skipping",
		    iBlock, evtData.fBlockCnt-1,
		    evtData.fEventID, evtData.fEventID, 
		    DataType2Text(inputDataType).c_str());
	  continue;
	}
      else 
	HLTDebug("We get the right data type: Block # %i/%i; Event 0x%08LX (%Lu) Received datatype: %s",
		    iBlock, evtData.fBlockCnt-1,
		    evtData.fEventID, evtData.fEventID, 
		    DataType2Text(inputDataType).c_str());
      
      
      TTree *clusterTree = 0x0;
      TClonesArray* clusterArray = 0x0;
      ReadAndLoadClusters(clusterTree, clusterArray, &block);
      
      // maybe it is not so smart to create it each event? clear is enough ?
      AliESDEvent *esd = new AliESDEvent();
      esd->CreateStdContent();
      fTracker->Clusters2Tracks(esd);

      //here transport the esd tracks further
      Int_t nTracks = esd->GetNumberOfTracks();
      Int_t nTRDTracks = esd->GetNumberOfTrdTracks();
      HLTInfo("Number of tracks  == %d == Number of TRD tracks %d", nTracks, nTRDTracks);  

      TClonesArray* trdTracks = fTracker->GetListOfTracks();
   
      if (trdTracks)
	totalSize += TransportTracks(trdTracks, outputPtr, outputBlocks, offset, dBlockSpecification);
      else 
	HLTDebug("Bad array trdTracks = 0x%x", trdTracks);
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
      delete esd;
      if (bWriteClusters)
	delete clusterTree;
      else{
	//clusterArray->Delete();
	delete clusterArray;
      }
      
    }
  size = totalSize;
  HLTDebug("Event is done. size written to the output is %i", size);
  //  CALLGRIND_STOP_INSTRUMENTATION();
  return 0;
}


/**
 * ReadClusters from the component input and load them to the tracker
 */
//============================================================================
void AliHLTTRDTrackerV1Component::ReadAndLoadClusters(TTree *inClusterTree, 
						      TClonesArray *inClusterArray, const AliHLTComponentBlockData *inBlock)
{
  Bool_t bWriteClusters = fReconstructor->IsWritingClusters();
  const AliHLTComponentBlockData &block = *inBlock;


  TObject *tobjin = 0x0;
  int ibForce = 0; // almost obsolet

  if (bWriteClusters){
    tobjin = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkClusterDataType, "TTree", ibForce);
    HLTDebug("1stBLOCK; Pointer = 0x%x", tobjin);
    inClusterTree = (TTree*)tobjin;
    if (inClusterTree)
      {
	HLTDebug("CLUSTERS; Pointer to TTree = 0x%x Name = %s", inClusterTree, inClusterTree->GetName());
	HLTDebug("TTree of clusters: nbEntries = %i", inClusterTree->GetEntriesFast());
	fTracker->LoadClusters(inClusterTree);
      }
    else
      {
	HLTError("First Input Block is not a TTree 0x%x", tobjin);
      }
  }
  else if (fUseHLTClusters)
    {
      inClusterArray = new TClonesArray("AliTRDcluster"); // would be nice to allocate memory for all clusters here.
      ReadClusters(inClusterArray, block.fPtr, block.fSize);
      HLTDebug("TClonesArray of clusters: nbEntries = %i", inClusterArray->GetEntriesFast());
      fTracker->LoadClusters(inClusterArray);
    }
  else
    {
      tobjin = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkClusterDataType, "TClonesArray", ibForce);
      HLTDebug("1stBLOCK; Pointer = 0x%x", tobjin);
      inClusterArray = (TClonesArray*)tobjin;
      if (inClusterArray)
	{
	  HLTDebug("CLUSTERS; Pointer to TClonesArray = 0x%x Name = %s", inClusterArray, inClusterArray->GetName());
	  HLTDebug("TClonesArray of clusters: nbEntries = %i", inClusterArray->GetEntriesFast());
	  fTracker->LoadClusters(inClusterArray);
	}
      else
	{
	  HLTError("First Input Block not a TClonesArray 0x%x", tobjin);
	}
    }
}

/**
 * Read cluster to the TClonesArray from the memory 
 */
//============================================================================
Int_t AliHLTTRDTrackerV1Component::ReadClusters(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size)
{
  //HLTDebug("\nReading clusters from the Memory\n ============= \n");
  AliHLTTRDCluster * curCluster;
  UInt_t clusterSize = sizeof(AliHLTTRDCluster), curSize = 0;
  Int_t i=0;
  
  curCluster = (AliHLTTRDCluster*) inputPtr;
  while (curSize + clusterSize <= size)
    {
      //  HLTDebug(" fX = %f; fY = %f; fZ = %f", curCluster->fX, curCluster->fY, curCluster->fZ);

      AliTRDcluster* curTRDCluster = new((*outArray)[i]) AliTRDcluster();
      curCluster->ExportTRDCluster(curTRDCluster);
      //      HLTDebug(" fX = %f; fY = %f; fZ = %f", curTRDCluster->GetX(), curTRDCluster->GetY(), curTRDCluster->GetZ());
      curSize += clusterSize; 
      i++;
      curCluster++;
      //cout << " current readed size is " << curSize << "/" << size << endl;
    }
  
  return i;
}
/**
 * Transport tracks to the next component
 * Return Numbers of bytes written to the output
 */
//============================================================================
AliHLTUInt32_t AliHLTTRDTrackerV1Component::TransportTracks(TClonesArray *inTracksArray, AliHLTUInt8_t* output,
						    vector<AliHLTComponent_BlockData>& outputBlocks, AliHLTUInt32_t inOffset, AliHLTUInt32_t inSpec)
{
  Int_t nbTracks=inTracksArray->GetEntriesFast();
  if (nbTracks>0)
    {
      HLTDebug("We have an output array: pointer to inTracksArray = 0x%x, nbEntries = %i", inTracksArray, nbTracks);
      if (fUseHLTTracks){
	// Using low-level interface 
	// with interface classes
	AliHLTUInt32_t addedSize = AddToOutput(inTracksArray, output);
	
	// Fill block 
	AliHLTComponentBlockData bd;
	FillBlockData( bd );

	bd.fPtr = output;
	bd.fOffset = inOffset;
	bd.fSize = addedSize;
	bd.fSpecification = inSpec;
	bd.fDataType = AliHLTTRDDefinitions::fgkTRDSATracksDataType;
	outputBlocks.push_back( bd );
	HLTDebug("BD fPtr 0x%x, fOffset %i, fSize %i, fSpec 0x%x", bd.fPtr, bd.fOffset, bd.fSize, bd.fSpecification);
	
	return addedSize;
      }
      else{
	inTracksArray->BypassStreamer(kFALSE);
	PushBack(inTracksArray, AliHLTTRDDefinitions::fgkTRDSATracksDataType);
	return 0;
      }
	  
    }
  return 0;
  
}
