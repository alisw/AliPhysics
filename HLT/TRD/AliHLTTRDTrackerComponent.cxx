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

/** @file   AliHLTTRDTrackerComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  A TRDTracker processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TFile.h"
#include "TChain.h"

#include "AliHLTTRDTrackerComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliCDBManager.h"

#include "AliTRDReconstructor.h"
#include "AliESDEvent.h"
//#include "AliTRDtrackerHLT.h"
#include "AliTRDtracker.h"
#include "AliTRDCalibraFillHisto.h"
#include "AliMagFMaps.h"
#include "AliTRDcluster.h"
#include "AliESDfriend.h"
#include <cstdlib>
#include <cerrno>
#include <string>

// this is a global object used for automatic component registration, do not use this
AliHLTTRDTrackerComponent gAliHLTTRDTrackerComponent;

ClassImp(AliHLTTRDTrackerComponent);
    
AliHLTTRDTrackerComponent::AliHLTTRDTrackerComponent()
  : AliHLTProcessor()
  , fOutputPercentage(100) // By default we copy to the output exactly what we got as input  
  , fStrorageDBpath("local://$ALICE_ROOT")
  , fCDB(NULL)
  , fField(NULL)
  , fGeometryFileName("")
  , fGeometryFile(NULL)
  , fGeoManager(NULL)
  , fTracker(NULL)
{
  // Default constructor
  fCDB = AliCDBManager::Instance();
  //fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
  fCDB->SetRun(0);

  fGeometryFileName = getenv("ALICE_ROOT");
  fGeometryFileName += "/HLT/TRD/geometry.root";
}

AliHLTTRDTrackerComponent::~AliHLTTRDTrackerComponent()
{
  // Destructor
}

const char* AliHLTTRDTrackerComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDTracker"; // The ID of this component
}

void AliHLTTRDTrackerComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data  
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back( AliHLTTRDDefinitions::fgkClusterDataType );
}

AliHLTComponent_DataType AliHLTTRDTrackerComponent::GetOutputDataType()
{
  // Get the output data type
  return AliHLTTRDDefinitions::fgkClusterDataType;
}

void AliHLTTRDTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTRDTrackerComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDTrackerComponent;
};

int AliHLTTRDTrackerComponent::DoInit( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  fOutputPercentage = 100;
  int i = 0;
  char* cpErr;
  while ( i < argc )
    {
      Logging( kHLTLogDebug, "HLT::TRDTracker::DoInit", "Arguments", "argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_percentage" ) )
	{
	  if ( i+1>=argc )
	    {
	      Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "Missing Argument", "Missing output_percentage parameter");
	      return ENOTSUP;
	    }
	  Logging( kHLTLogDebug, "HLT::TRDTracker::DoInit", "Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	  fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "Wrong Argument", "Cannot convert output_percentage parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }
	  Logging( kHLTLogInfo, "HLT::TRDTracker::DoInit", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
	  i += 2;
	  continue;
	}

      if ( strcmp( argv[i], "-cdb" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "Missing Argument", "Missing -cdb argument");
	      return ENOTSUP;	      
	    }
	  fStrorageDBpath = argv[i+1];
	  Logging( kHLTLogInfo, "HLT::TRDTracker::DoInit", "DB storage set", "DB storage is %s", fStrorageDBpath.c_str() );	  
	  fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
	  i += 2;
	  continue;
	}      

      if ( strcmp( argv[i], "-geometry" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "Missing Argument", "Missing -geometry argument");
	      return ENOTSUP;	      
	    }
	  fGeometryFileName = argv[i+1];
	  Logging( kHLTLogInfo, "HLT::TRDTracker::DoInit", "GeomFile storage set", "GeomFile storage is %s", 
		   fGeometryFileName.c_str() );	  
	  i += 2;
	  continue;
	}      

      Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
      return EINVAL;
    }

  //init alifield map - temporarly fixed - should come from a DB
  fField = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  if (fField)
    AliTracker::SetFieldMap(fField,1);
  else
    Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "Field", "Unable to init the field");
    
  fGeometryFile = TFile::Open(fGeometryFileName.c_str());
  if (fGeometryFile)
    {
      fGeoManager = (TGeoManager *)fGeometryFile->Get("Geometry");
      //fTracker = new AliTRDtrackerHLT(fGeometryFile);
      fTracker = new AliTRDtracker(fGeometryFile);
      //fTracker = new AliTRDtracker(fGeometryFile);
    }
  else
    {
      Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "fGeometryFile", "Unable to open file. FATAL!");
      return -1;
    }

  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
  if (calibra == 0)
    {
      Logging(kHLTLogError, "HLT::TRDTracker::DoInit", "Calibration Histos", "::Instance failed");
      return -1;      
    }
  else
    {
      calibra->SetMITracking(kTRUE);
      calibra->Init2Dhistos();
    }

  return 0;
}

int AliHLTTRDTrackerComponent::DoDeinit()
{
  // Deinitialization of the component

  delete fField;
  fField = 0;

  delete fTracker;
  fTracker = 0;
  
  if (fGeometryFile)
    {
      fGeometryFile->Close();
      delete fGeometryFile;
      fGeometryFile = 0;
    }

  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
  if (calibra)
    {
      // should not write in here!
      calibra->Write2d();
      Logging( kHLTLogInfo, "HLT::TRDTracker::DoDeinit", "CALIBRA", "before destroy");
      calibra->Destroy();
      Logging( kHLTLogInfo, "HLT::TRDTracker::DoDeinit", "CALIBRA", "after destroy");
    }

  return 0;
}

int AliHLTTRDTrackerComponent::DoEvent( const AliHLTComponentEventData & evtData,
					AliHLTComponentTriggerData & trigData )
{
  // Process an event
  
  Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
  Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "BLOCKS", "NofBlocks %lu", evtData.fBlockCnt );

  AliHLTUInt32_t fDblock_Specification = 0;

  //implement a usage of the following
//   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
//   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
//   void *triggerData = trigData.fData;
  Logging( kHLTLogDebug, "HLT::TRDTracker::DoEvent", "Trigger data received", 
	   "Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  AliHLTComponentBlockData *dblock = (AliHLTComponentBlockData *)GetFirstInputBlock( AliHLTTRDDefinitions::fgkClusterDataType );
  if (dblock != 0)
    {
      fDblock_Specification = dblock->fSpecification;
    }
  else
    {
      Logging( kHLTLogWarning, "HLT::TRDTracker::DoEvent", "DATAIN", "First Input Block not found! 0x%x", dblock);
      return -1;
    }

  int ibForce = 0;
  TObject *tobjin = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkClusterDataType, "TTree", ibForce);
  Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "1stBLOCK", "Pointer = 0x%x", tobjin);

  TTree *clusterTree = (TTree*)tobjin;
  if (!clusterTree)
    {
      Logging( kHLTLogWarning, "HLT::TRDTracker::DoEvent", "DATAIN", "First Input Block not a tree! 0x%x", tobjin);
      return -1;
    }

  Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "1stBLOCK", "Pointer = 0x%x Name = %s", clusterTree, clusterTree->GetName());

  while (tobjin != 0)
    {
      if (clusterTree)
	{
	  Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "CLUSTERS", "Pointer = 0x%x Name = %s", clusterTree, clusterTree->GetName());
	  Int_t iNentries = clusterTree->GetEntries();
	  Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "COUNT", "N of tree entries = %d", iNentries);
	  fTracker->LoadClusters(clusterTree);
	}
      else
	{
	  Logging( kHLTLogError, "HLT::TRDTracker::DoEvent", "CLUSTERS", "Tree Pointer = 0x%x", clusterTree);
	}

      tobjin = (TObject *)GetNextInputObject( ibForce );
      Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "nextBLOCK", "Pointer = 0x%x", tobjin);
      clusterTree = (TTree*)tobjin;
    }

  AliTRDReconstructor::SetSeedingOn(kTRUE);

  fTracker->SetAddTRDseeds();

  AliESDfriend *esdFriend = new AliESDfriend();

  AliESDEvent *esd = new AliESDEvent();
  esd->CreateStdContent();
  fTracker->PropagateBack(esd);
  fTracker->RefitInward(esd);
  fTracker->Clusters2Tracks(esd);

  esd->GetESDfriend(esdFriend);
  //here transport the esd tracks further

  Int_t nTracks = esd->GetNumberOfTracks();
  Int_t nTRDTracks = esd->GetNumberOfTrdTracks();
  Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "DONE", "Number of tracks %d Number of TRD tracks %d", nTracks, nTRDTracks);

//   AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
//   calibra->Init2Dhistostrack();

  for (Int_t it = 0; it < nTracks; it++)
    {
      AliESDtrack* track = esd->GetTrack(it);

//       Int_t nCalibObjects = 0;
//       Int_t idx = 0;
//       while (track->GetCalibObject(idx) != 0)
// 	{
// 	  nCalibObjects++;
// 	  idx++;
// 	}
//       Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "DONE", "Track 0x%x NcalibObjects %d", track, nCalibObjects);

      Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "DONE", "Track %d 0x%x Pt %1.2f", it, track, track->Pt());
      PushBack(track, AliHLTTRDDefinitions::fgkTRDSATracksDataType, ++fDblock_Specification);
//       if (calibra->GetMItracking())
//  	{
//  	  calibra->UpdateHistograms(track);
// 	}
    }

  //PushBack(esd, AliHLTTRDDefinitions::fgkTRDSAEsdDataType, fDblock_Specification);
  //PushBack(esdFriend, AliHLTTRDDefinitions::fgkTRDSAEsdDataType, fDblock_Specification);

  //no receiver defined yet(!)
  //Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "DONE", "now deleting");
  delete esd;
  delete esdFriend;

  //Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "DONE", "after delete esd");
  delete clusterTree;

  //Logging( kHLTLogInfo, "HLT::TRDTracker::DoEvent", "DONE", "after delete clusterTree");
  return 0;
}
