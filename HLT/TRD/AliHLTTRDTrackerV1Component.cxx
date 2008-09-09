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

#include "TFile.h"
#include "TChain.h"

#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliMagFMaps.h"
#include "AliESDfriend.h"

#include "AliTRDcalibDB.h"
#include "AliTRDReconstructor.h"
#include "AliTRDtrackerV1.h"
#include "AliTRDrecoParam.h"

#include <cstdlib>
#include <cerrno>
#include <string>

// this is a global object used for automatic component registration, do not use this
AliHLTTRDTrackerV1Component gAliHLTTRDTrackerV1Component;

ClassImp(AliHLTTRDTrackerV1Component);
    
AliHLTTRDTrackerV1Component::AliHLTTRDTrackerV1Component()
  : AliHLTProcessor()
  , fOutputPercentage(100) // By default we copy to the output exactly what we got as input  
  , fStrorageDBpath("local://$ALICE_ROOT")
  , fCDB(NULL)
  , fField(NULL)
  , fGeometryFileName("")
  , fGeometryFile(NULL)
  , fTracker(NULL)
  , fRecoParam(NULL)
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
  
  if (iMagneticField == 0)
    {
      // magnetic field OFF
      fField = new AliMagFMaps("Maps","Maps", 2, 0., 10., 1);
      HLTDebug("Magnetic field is OFF.");
    }

  if (iMagneticField == 1)
    {
      // magnetic field ON
      fField = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
      HLTDebug("Magnetic field is ON.");
    }

  if (fField == 0)
    {
      HLTError("Unable to init the field. Trouble at this point.");
      return -1;
    }

  // kTRUE sets the map uniform
  AliTracker::SetFieldMap(fField,kTRUE);

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
  TString recoOptions="sa,sl_cf_0";
  
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
  
  fGeometryFile = 0;
  fGeometryFile = TFile::Open(fGeometryFileName.c_str());
  if (fGeometryFile)
    {
      AliGeomManager::LoadGeometry(fGeometryFileName.c_str());
    }
  else
    {
      HLTError("Unable to open file. FATAL!");
      return -1;
    }
  
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

  delete fField;
  fField = 0x0;

  // fTracker->SetClustersOwner(kFALSE);
  delete fTracker;
  fTracker = 0x0;
  
  // We need to set clusters in Reconstructor to null to prevent from 
  // double deleting, since we delete TClonesArray by ourself in DoEvent.
  fReconstructor->SetClusters(0x0);
  delete fReconstructor;
  fReconstructor = 0x0;
  
  AliTRDcalibDB::Terminate();

  if (fGeometryFile)
    {
      fGeometryFile->Close();
      delete fGeometryFile;
      fGeometryFile = 0x0;
    }

  return 0;
}

int AliHLTTRDTrackerV1Component::DoEvent( const AliHLTComponentEventData & evtData,
					  AliHLTComponentTriggerData & trigData )
{
  // Process an event
  Bool_t bWriteClusters = fReconstructor->IsWritingClusters();

  HLTDebug("Output percentage set to %lu %%", fOutputPercentage );
  HLTDebug("NofBlocks %lu", evtData.fBlockCnt );
  
  AliHLTUInt32_t dBlockSpecification = 0;


  //implement a usage of the following
  //   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
  //   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
  //   void *triggerData = trigData.fData;
  HLTDebug("Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);
  
  AliHLTComponentBlockData *dblock = (AliHLTComponentBlockData *)GetFirstInputBlock( AliHLTTRDDefinitions::fgkClusterDataType );
  if (dblock != 0)
    {
      AliHLTComponentDataType inputDataType = dblock->fDataType;
      HLTDebug( "Event 0x%08LX (%Lu) received datatype: %s",
		evtData.fEventID, evtData.fEventID, 
		DataType2Text(inputDataType).c_str());
      dBlockSpecification = dblock->fSpecification;
    }
  else
    {
      HLTWarning("First Input Block not found! 0x%x", dblock);
      return 0;
    }

  int ibForce = 0;
  
  TObject *tobjin = 0x0;

  TTree *clusterTree = 0x0;
  TClonesArray *clusterArray = 0x0;
  if (bWriteClusters){
    tobjin = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkClusterDataType, "TTree", ibForce);
    HLTDebug("1stBLOCK; Pointer = 0x%x", tobjin);
    clusterTree = (TTree*)tobjin;
  
    //     while (tobjin != 0)
    //         {
    if (clusterTree)
      {
	HLTDebug("CLUSTERS; Pointer to TTree = 0x%x Name = %s", clusterTree, clusterTree->GetName());
	HLTDebug("TTree of clusters: nbEntries = %i", clusterTree->GetEntriesFast());
	fTracker->LoadClusters(clusterTree);
      }
    else
      {
	HLTError("First Input Block not a TTree 0x%x", tobjin);
      }
  }
  else{
    tobjin = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkClusterDataType, "TClonesArray", ibForce);
    HLTDebug("1stBLOCK; Pointer = 0x%x", tobjin);
    clusterArray = (TClonesArray*)tobjin;
    if (clusterArray)
      {
	HLTDebug("CLUSTERS; Pointer to TClonesArray = 0x%x Name = %s", clusterArray, clusterArray->GetName());
	HLTDebug("TClonesArray of clusters: nbEntries = %i", clusterArray->GetEntriesFast());
	Int_t nb = clusterArray->GetEntriesFast();
	for (Int_t i=0; i<nb; i++){
	  AliTRDcluster * cluster = (AliTRDcluster* ) clusterArray->At(i);
	  //HLTDebug("Cluster[%i]: detector %i", i, cluster->GetDetector());
	}
  
	fTracker->LoadClusters(clusterArray);
      }
    else
      {
	HLTError("First Input Block not a TClonesArray 0x%x", tobjin);
      }
  }
  
  // maybe it is not so smart to create it each event? clear is enough ?
  AliESDEvent *esd = new AliESDEvent();
  esd->CreateStdContent();
  fTracker->Clusters2Tracks(esd);

  //here transport the esd tracks further
  Int_t nTracks = esd->GetNumberOfTracks();
  Int_t nTRDTracks = esd->GetNumberOfTrdTracks();
  HLTInfo("Number of tracks  == %d == Number of TRD tracks %d", nTracks, nTRDTracks);  
  //esd->Print();
  //PushBack(esd, AliHLTTRDDefinitions::fgkTRDSAEsdDataType);

  // extract the friend ?
  //   AliESDfriend *esdFriend = new AliESDfriend();
  //   esd->GetESDfriend(esdFriend);
  //   PushBack(esdFriend, AliHLTTRDDefinitions::fgkTRDSAEsdDataType);
  //   delete esdFriend;

   TClonesArray* trdTracks = fTracker->GetListOfTracks();
   
   if (trdTracks)
    {
      nTracks=trdTracks->GetEntriesFast();
      if (nTracks>0)
	{
	  HLTDebug("We have an output array: pointer to trdTracks = 0x%x, nbEntries = %i", trdTracks, nTracks);
      HLTDebug("We have an output array: pointer to trdTracks = 0x%x, nbEntries = %i", trdTracks, nTracks);
      
      trdTracks->BypassStreamer(kFALSE);
      PushBack(trdTracks, AliHLTTRDDefinitions::fgkTRDSATracksDataType);
      //trdTracks->Delete();
      //delete trdTracks;
	}
    }
  else 
    HLTDebug("Bad array trdTracks = 0x%x", trdTracks);


  
  //here we are deleting clusters (but not the TClonesArray itself)
  fTracker->UnloadClusters();
  delete esd;
  if (bWriteClusters)
    delete clusterTree;
  else
    delete clusterArray;

  HLTDebug("Event done.");
  return 0;
}

///////////////////////////////
/*
  consider transporting TRD tracks only as they might have signigicantly smaller size... on the other hand you will need to prodece ESDs at some point...

  // this is for ESDtrack
  //   for (Int_t it = 0; it < nTracks; it++)
  //     {
  //       AliESDtrack* track = esd->GetTrack(it);
  //       HLTInfo("Track %d 0x%x Pt %1.2f", it, track, track->Pt());
  //       //PushBack(track, AliHLTTRDDefinitions::fgkTRDSATracksDataType, ++dBlockSpecification);
  //       PushBack(track, AliHLTTRDDefinitions::fgkTRDSATracksDataType);
  //     }

  // one can do similar things with the TRDtrack
  esd->GetNumberOfTrdTracks();
  and then
  for (i;;)
  AliESDTrdTrack *trdtrack = esd->GetTrdTrack(i)
*/

