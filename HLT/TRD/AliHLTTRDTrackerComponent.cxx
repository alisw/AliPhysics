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
#include "AliTRDrecoParam.h"
#include "AliESDEvent.h"
//#include "AliTRDtrackerHLT.h"
#include "AliTRDtracker.h"
#include "AliTRDCalibraFillHisto.h"
#include "AliMagF.h"
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
  , fGeometryFileName("")
  , fGeometryFile(NULL)
  , fGeoManager(NULL)
  , fTracker(NULL)
{
  // Default constructor

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
      HLTDebug("argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_percentage" ) )
	{
	  if ( i+1>=argc )
	    {
	      HLTError("Missing output_percentage parameter");
	      return ENOTSUP;
	    }
	  HLTDebug("Arguments", "argv[%d+1] == %s", i, argv[i+1] );
	  fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      HLTError("Cannot convert output_percentage parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }
	  HLTInfo("Output percentage set to %lu %%", fOutputPercentage );
	  i += 2;
	  continue;
	}

      if ( strcmp( argv[i], "-cdb" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError("Missing -cdb argument");
	      return ENOTSUP;	      
	    }
	  fStrorageDBpath = argv[i+1];
	  HLTInfo("DB storage is %s", fStrorageDBpath.c_str() );	  
	  i += 2;
	  continue;
	}      

      if ( strcmp( argv[i], "-geometry" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError("Missing -geometry argument");
	      return ENOTSUP;	      
	    }
	  fGeometryFileName = argv[i+1];
	  HLTInfo("GeomFile storage is %s", fGeometryFileName.c_str() );	  
	  i += 2;
	  continue;
	}      

      HLTError("Unknown option '%s'", argv[i] );
      return EINVAL;
    }


  fCDB = AliCDBManager::Instance();
  if (!fCDB)
    {
      HLTError("Could not get CDB instance", "fCDB 0x%x", fCDB);
    }
  else
    {
      fCDB->SetRun(0); // THIS HAS TO BE RETRIEVED !!!
      fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
      HLTDebug("fCDB 0x%x", fCDB);
    }
    
  fGeometryFile = TFile::Open(fGeometryFileName.c_str());
  if (fGeometryFile)
    {
      fGeoManager = (TGeoManager *)fGeometryFile->Get("Geometry");
      //fTracker = new AliTRDtrackerHLT(fGeometryFile);
      AliTRDrecoParam *fPars = AliTRDrecoParam::GetLowFluxParam();
      //fPars->SetSeeding(kTRUE);
      //fPars->SetStreamLevel(0);
      AliTRDReconstructor reconstructor; reconstructor.SetRecoParam(fPars);
      // write clusters [cw] = true
      // track seeding (stand alone tracking) [sa] = true
      // PID method in reconstruction (NN) [nn] = true
      // write online tracklets [tw] = false
      // drift gas [ar] = false
      reconstructor.SetOption("cw,sa");
      fTracker = new AliTRDtracker(fGeometryFile);
      //fTracker = new AliTRDtracker(fGeometryFile);
    }
  else
    {
      HLTError("Unable to open file. FATAL!");
      return -1;
    }

  AliTRDCalibraFillHisto *calibra = AliTRDCalibraFillHisto::Instance();
  if (calibra == 0)
    {
      HLTError("Calibration Histos ::Instance failed");
      return -1;      
    }
  else
    {
      calibra->Init2Dhistos();
    }

  return 0;
}

int AliHLTTRDTrackerComponent::DoDeinit()
{
  // Deinitialization of the component

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
      calibra->Destroy();
    }

  return 0;
}

int AliHLTTRDTrackerComponent::DoEvent( const AliHLTComponentEventData & /*evtData*/,
					AliHLTComponentTriggerData & /*trigData*/ )
{
  // Process an event
  
  HLTInfo("Output percentage set to %lu %%", fOutputPercentage );
  HLTInfo("NofBlocks %lu", GetNumberOfInputBlocks() );

  AliHLTUInt32_t dBlockSpecification = 0;

  //implement a usage of the following
  //   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
  //   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
  //   void *triggerData = trigData.fData;
  //HLTDebug("Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  AliHLTComponentBlockData *dblock = (AliHLTComponentBlockData *)GetFirstInputBlock( AliHLTTRDDefinitions::fgkClusterDataType );
  if (dblock != 0)
    {
      dBlockSpecification = dblock->fSpecification;
    }
  else
    {
      HLTWarning("First Input Block not found! 0x%x", dblock);
      return -1;
    }

  int ibForce = 0;
  TObject *tobjin = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkClusterDataType, "TTree", ibForce);
  HLTInfo("Pointer = 0x%x", tobjin);

  TTree *clusterTree = (TTree*)tobjin;
  if (!clusterTree)
    {
      HLTWarning("First Input Block not a tree! 0x%x", tobjin);
      return -1;
    }

  HLTInfo("Pointer = 0x%x Name = %s", clusterTree, clusterTree->GetName());

  while (tobjin != 0)
    {
      if (clusterTree)
	{
	  HLTInfo("Pointer = 0x%x Name = %s", clusterTree, clusterTree->GetName());
	  Int_t iNentries = clusterTree->GetEntries();
	  HLTInfo("N of tree entries = %d", iNentries);
	  fTracker->LoadClusters(clusterTree);
	}
      else
	{
	  HLTError("Tree Pointer = 0x%x", clusterTree);
	}

      tobjin = (TObject *)GetNextInputObject( ibForce );
      HLTInfo("Pointer = 0x%x", tobjin);
      clusterTree = (TTree*)tobjin;
    }

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
  HLTInfo( "Number of tracks %d Number of TRD tracks %d", nTracks, nTRDTracks);


  for (Int_t it = 0; it < nTracks; it++)
    {
      AliESDtrack* track = esd->GetTrack(it);
      HLTInfo("Track %d 0x%x Pt %1.2f", it, track, track->Pt());
      PushBack(track, AliHLTTRDDefinitions::fgkTRDSATracksDataType, ++dBlockSpecification);
    }

  delete esd;
  delete esdFriend;

  delete clusterTree;

  return 0;
}
