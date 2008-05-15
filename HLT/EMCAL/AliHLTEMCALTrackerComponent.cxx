// $Id: AliHLTEMCALTrackerComponent.cxx 23618 2008-01-29 13:07:38Z hristov $

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

/** @file   AliHLTEMCALTrackerComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  A EMCALTracker processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TFolder.h"
#include "TFile.h"

#include "AliHLTEMCALTrackerComponent.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliCDBManager.h"
#include "AliEMCALTracker.h"
#include "AliEMCALReconstructor.h"
#include "AliESDEvent.h"
#include "AliMagFMaps.h"
#include "AliESDfriend.h"

#include <cstdlib>
#include <cerrno>
#include <string>

// this is a global object used for automatic component registration, do not use this
AliHLTEMCALTrackerComponent gAliHLTEMCALTrackerComponent;

ClassImp(AliHLTEMCALTrackerComponent);
    
AliHLTEMCALTrackerComponent::AliHLTEMCALTrackerComponent()
  : AliHLTProcessor()
  , fOutputPercentage(100) // By default we copy to the output exactly what we got as input  
  , fStrorageDBpath("local://$ALICE_ROOT")
  , fCDB(NULL)
  , fField(NULL)
  , fGeometryFileName("")
  , fGeometryFile(NULL)
  , fGeoManager(NULL)
  , fTracker(NULL)
  , fInputFolder(new TFolder("EMCALtrackerFolder", "EMCALtrackerFolder"))
{
  // Default constructor
  fInputFolder->SetOwner(kTRUE);

  fGeometryFileName = getenv("ALICE_ROOT");
  fGeometryFileName += "/HLT/EMCAL/geometry.root";
}

AliHLTEMCALTrackerComponent::~AliHLTEMCALTrackerComponent()
{
  // Destructor
  delete fInputFolder;
  fInputFolder = NULL;
}

const char* AliHLTEMCALTrackerComponent::GetComponentID()
{
  // Return the component ID const char *
  return "EMCALTracker"; // The ID of this component
}

void AliHLTEMCALTrackerComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data  
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back( AliHLTEMCALDefinitions::fgkClusterDataType );
}

AliHLTComponent_DataType AliHLTEMCALTrackerComponent::GetOutputDataType()
{
  // Get the output data type
  return AliHLTEMCALDefinitions::fgkEMCALESDDataType;
}

void AliHLTEMCALTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTEMCALTrackerComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTEMCALTrackerComponent;
};

int AliHLTEMCALTrackerComponent::DoInit( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  fOutputPercentage = 100;
  int i = 0;
  char* cpErr;
  while ( i < argc )
    {
      Logging( kHLTLogDebug, "HLT::EMCALTracker::DoInit", "Arguments", "argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_percentage" ) )
	{
	  if ( i+1>=argc )
	    {
	      Logging(kHLTLogError, "HLT::EMCALTracker::DoInit", "Missing Argument", "Missing output_percentage parameter");
	      return ENOTSUP;
	    }
	  Logging( kHLTLogDebug, "HLT::EMCALTracker::DoInit", "Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	  fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      Logging(kHLTLogError, "HLT::EMCALTracker::DoInit", "Wrong Argument", "Cannot convert output_percentage parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }
	  Logging( kHLTLogInfo, "HLT::EMCALTracker::DoInit", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
	  i += 2;
	  continue;
	}

      if ( strcmp( argv[i], "-cdb" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      Logging(kHLTLogError, "HLT::EMCALTracker::DoInit", "Missing Argument", "Missing -cdb argument");
	      return ENOTSUP;	      
	    }
	  fStrorageDBpath = argv[i+1];
	  Logging( kHLTLogInfo, "HLT::EMCALTracker::DoInit", "DB storage set", "DB storage is %s", fStrorageDBpath.c_str() );	  
	  i += 2;
	  continue;
	}      

      if ( strcmp( argv[i], "-geometry" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      Logging(kHLTLogError, "HLT::EMCALTracker::DoInit", "Missing Argument", "Missing -geometry argument");
	      return ENOTSUP;	      
	    }
	  fGeometryFileName = argv[i+1];
	  Logging( kHLTLogInfo, "HLT::EMCALTracker::DoInit", "GeomFile storage set", "GeomFile storage is %s", 
		   fGeometryFileName.c_str() );	  
	  i += 2;
	  continue;
	}      

      Logging(kHLTLogError, "HLT::EMCALTracker::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
      return EINVAL;
    }

  //init alifield map - temporarly fixed - should come from a DB
  fField = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  if (fField)
    AliTracker::SetFieldMap(fField,1);
  else
    Logging(kHLTLogError, "HLT::EMCALTracker::DoInit", "Field", "Unable to init the field");

  fCDB = AliCDBManager::Instance();
  if (!fCDB)
    {
      Logging(kHLTLogError, "HLT::EMCALCalibration::DoInit", "Could not get CDB instance", "fCDB 0x%x", fCDB);
    }
  else
    {
      fCDB->SetRun(0); // THIS HAS TO BE RETRIEVED !!!
      fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
      Logging(kHLTLogDebug, "HLT::EMCALCalibration::DoInit", "CDB instance", "fCDB 0x%x", fCDB);
    }
    
  fGeometryFile = TFile::Open(fGeometryFileName.c_str());
  if (fGeometryFile)
    {
      fGeoManager = (TGeoManager *)fGeometryFile->Get("Geometry");
      fTracker = new AliEMCALTracker;
    }
  else
    {
      Logging(kHLTLogError, "HLT::EMCALTracker::DoInit", "fGeometryFile", "Unable to open file. FATAL!");
      return -1;
    }

  return 0;
}

int AliHLTEMCALTrackerComponent::DoDeinit()
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

  fInputFolder->Clear();

  return 0;
}

int AliHLTEMCALTrackerComponent::DoEvent( const AliHLTComponentEventData & evtData,
					AliHLTComponentTriggerData & trigData )
{
  // Process an event
  
  Logging( kHLTLogInfo, "HLT::EMCALTracker::DoEvent", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
  Logging( kHLTLogInfo, "HLT::EMCALTracker::DoEvent", "BLOCKS", "NofBlocks %lu", evtData.fBlockCnt );

  AliHLTUInt32_t dBlockSpecification = 0;

  //implement a usage of the following
//   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
//   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
//   void *triggerData = trigData.fData;
  Logging( kHLTLogDebug, "HLT::EMCALTracker::DoEvent", "Trigger data received", 
	   "Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  AliHLTComponentBlockData *dblock = (AliHLTComponentBlockData *)GetFirstInputBlock( AliHLTEMCALDefinitions::fgkClusterDataType );
  if (dblock != 0)
    {
      dBlockSpecification = dblock->fSpecification;
    }
  else
    {
      Logging( kHLTLogWarning, "HLT::EMCALTracker::DoEvent", "DATAIN", "First Input Block not found! 0x%x", dblock);
      return -1;
    }

  fInputFolder->Clear();
  AliESDEvent *esd = 0; // we assume we receive this one from a global merger component
  TTree *clusterTree = 0;
  
  int ibForce = 0;
  TObject *tobjin = 0;

  // here getfirstinput finds the first object with spec type..
  // get the clusters tree
  tobjin = (TObject *)GetFirstInputObject( AliHLTEMCALDefinitions::fgkClusterDataType, "TTree", ibForce);
  while (tobjin != 0)
    {
      tobjin = (TObject *)GetNextInputObject( ibForce );
      Logging( kHLTLogInfo, "HLT::EMCALTracker::DoEvent", "nextBLOCK", "Pointer = 0x%x", tobjin);
      clusterTree = (TTree*)tobjin;
      if (clusterTree != 0)
	{
	  Int_t iNentries = clusterTree->GetEntries();
	  Logging( kHLTLogInfo, "HLT::EMCALTracker::DoEvent", "COUNT", "N of tree entries = %d", iNentries);
	  fTracker->LoadClusters(clusterTree);
	  fInputFolder->Add(clusterTree);      
	}
    }

  // now get the ESD(s) - should in principle be only one...
  tobjin = (TObject *)GetFirstInputObject( AliHLTEMCALDefinitions::fgkESDDataType, "AliESDevent", ibForce);
  while (tobjin != 0)
    {
      esd = (AliESDEvent *)tobjin;
      if (esd != 0)
	{
	  HLTInfo("Got ESDevent");
	  fInputFolder->Add(esd);	  
	}
      tobjin = (TObject *)GetNextInputObject( ibForce );
      Logging( kHLTLogInfo, "HLT::EMCALTracker::DoEvent", "nextBLOCK", "Pointer = 0x%x", tobjin);
      clusterTree = (TTree*)tobjin;
    }

  esd = (AliESDEvent*)fInputFolder->FindObject("AliESDevent");
  if (esd != 0)
    {
      fTracker->PropagateBack(esd);
      
      //here transport the esd tracks further
      Int_t nTracks = esd->GetNumberOfTracks();
      HLTInfo("Number of tracks %d", nTracks);  
      //esd->Print();
      PushBack(esd, AliHLTEMCALDefinitions::fgkEMCALESDDataType);      
    }
  else
    {
      HLTError("No ESD events received!");
    }

  fTracker->UnloadClusters();  
  fInputFolder->Clear();

  HLTDebug("Event done.");
  return 0;
}
