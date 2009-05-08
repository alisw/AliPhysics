// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCOfflineTrackerComponent.cxx
    @author Jacek Otwinowski & Matthias Richter
    @date   
    @brief  Wrapper component to the TPC offline tracker
*/

#include "AliHLTTPCOfflineTrackerComponent.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliVParticle.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGeomManager.h"
#ifndef HAVE_NOT_ALIMAGF30848
#include "AliMagF.h"
#else // keeping the <30489 code for backward compatibility, to be merged at some point
#include "AliMagFMaps.h"
#endif
#include "AliTPCReconstructor.h"
#include "AliTPCParam.h"
#include "AliTPCRecoParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCtrackerMI.h"
#include "AliTPCClustersRow.h"
#include "AliTPCseed.h"
#include "AliESDEvent.h"
#include "AliHLTTPCDefinitions.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCOfflineTrackerComponent)

AliHLTTPCOfflineTrackerComponent::AliHLTTPCOfflineTrackerComponent() : AliHLTProcessor(),
fGeometryFileName(""),
fTPCGeomParam(0),
fTracker(0),
fESD(0)
{
  // Default constructor
  fGeometryFileName = getenv("ALICE_ROOT");
  fGeometryFileName += "/HLT/TPCLib/offline/geometry.root";
}

AliHLTTPCOfflineTrackerComponent::~AliHLTTPCOfflineTrackerComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCOfflineTrackerComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCOfflineTracker";
}

void AliHLTTPCOfflineTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // get input data type
  list.push_back(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineClustersDataType*/);
}

AliHLTComponentDataType AliHLTTPCOfflineTrackerComponent::GetOutputDataType()
{
  // create output data type
  return kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineTrackSegmentsDataType*/;
}

void AliHLTTPCOfflineTrackerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // get output data size
  constBase = 2000000;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCOfflineTrackerComponent::Spawn()
{
  // create instance of the component
  return new AliHLTTPCOfflineTrackerComponent;
}

int AliHLTTPCOfflineTrackerComponent::DoInit( int argc, const char** argv )
{
  // init configuration 
  //
  int iResult=0;
#ifdef HAVE_NOT_TPCOFFLINE_REC
  HLTFatal("AliRoot version > v4-13-Release required");
  return -ENOSYS;
#endif

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;

  // loop over input parameters
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    if (argument.CompareTo("-geometry")==0) {
      if ((bMissingParam=(++i>=argc))) break;

      HLTInfo("got \'-geometry\' argument: %s", argv[i]);
      fGeometryFileName = argv[i];
      HLTInfo("Geometry file is: %s", fGeometryFileName.c_str());

      // the remaining arguments are treated as configuration
    } else {
      if (!configuration.IsNull()) configuration+=" ";
      configuration+=argument;
    }
  } // end loop

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult>=0 && !configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

  //
  // initialisation
  //
 
  // Load geometry
  AliGeomManager::LoadGeometry(fGeometryFileName.c_str());
  if((AliGeomManager::GetGeometry()) == 0) {
    HLTError("Cannot load geometry from file %s",fGeometryFileName.c_str());
    iResult=-EINVAL;
  }

  // TPC reconstruction parameters
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetHLTParam();
  if(tpcRecoParam) {
    tpcRecoParam->SetClusterSharing(kTRUE);

    AliTPCReconstructor tpcReconstructor;
    tpcReconstructor.SetRecoParam(tpcRecoParam);
  }
 
  // TPC geometry parameters
  fTPCGeomParam = new AliTPCParamSR;
  if (fTPCGeomParam) {
    fTPCGeomParam->ReadGeoMatrices();
  }

  // Init tracker
  fTracker = new AliTPCtrackerMI(fTPCGeomParam);

  // AliESDEvent event needed by AliTPCtrackerMI
  // output of the component
  fESD = new AliESDEvent();
  if (fESD) {
    fESD->CreateStdContent();
  }

  if (!fTracker || !fESD || !fTPCGeomParam) {
    HLTError("failed creating internal objects");
    iResult=-ENOMEM;
  }

  if (iResult>=0) {
    // read the default CDB entries
    iResult=Reconfigure(NULL, NULL);
  }

  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::DoDeinit()
{
  // deinit configuration

  if(fTPCGeomParam) delete fTPCGeomParam; fTPCGeomParam = 0; 
  if(fTracker) delete fTracker; fTracker = 0; 
  if(fESD) delete fESD; fESD = 0;

  return 0;
}

int AliHLTTPCOfflineTrackerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // tracker function
  HLTInfo("DoEvent processing data");

  int iResult=0;
  TClonesArray *clusterArray=0;
  int slice, patch;

  const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC); 
  if(!pBlock) {
     HLTError("Cannot get first data block 0x%08x ",pBlock);
     iResult=-ENOMEM; return iResult;
  }
  int minSlice=AliHLTTPCDefinitions::GetMinSliceNr(pBlock->fSpecification);
  int maxSlice=AliHLTTPCDefinitions::GetMaxSliceNr(pBlock->fSpecification);
  int minPatch=AliHLTTPCDefinitions::GetMinPatchNr(pBlock->fSpecification);
  int maxPatch=AliHLTTPCDefinitions::GetMaxPatchNr(pBlock->fSpecification);  

  if (fTracker && fESD) {
      // loop over input data blocks: TClonesArrays of clusters
      for (TObject *pObj = (TObject *)GetFirstInputObject(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineClustersDataType*/,"TClonesArray",0);
	 pObj !=0 && iResult>=0;
	 pObj = (TObject *)GetNextInputObject(0)) {
      clusterArray = dynamic_cast<TClonesArray*>(pObj);
      if (!clusterArray) continue;

      HLTInfo("load %d clusters from block %s 0x%08x", clusterArray->GetEntries(), DataType2Text(GetDataType(pObj)).c_str(), GetSpecification(pObj));
      slice=AliHLTTPCDefinitions::GetMinSliceNr(GetSpecification(pObj));
      patch=AliHLTTPCDefinitions::GetMinPatchNr(GetSpecification(pObj));

      if(slice < minSlice) minSlice=slice;
      if(slice > maxSlice) maxSlice=slice;
      if(patch < minPatch) minPatch=patch;
      if(patch > maxPatch) maxPatch=patch;
#ifndef HAVE_NOT_TPCOFFLINE_REC
      fTracker->LoadClusters(clusterArray);
#endif //HAVE_NOT_TPCOFFLINE_REC
 
    clusterArray->Delete();
    }// end loop over input objects

#ifndef HAVE_NOT_TPCOFFLINE_REC
    // Load outer sectors
      fTracker->LoadOuterSectors();
    // Load inner sectors
      fTracker->LoadInnerSectors();
#endif

    // set magnetic field for the ESD, assumes correct initialization of
    // the field map
    fESD->SetMagneticField(AliTracker::GetBz());
  
    // run tracker
    fTracker->Clusters2Tracks(fESD);

    // unload clusters
    fTracker->UnloadClusters();

    Int_t nTracks = fESD->GetNumberOfTracks();
    HLTInfo("Number TPC tracks %d", nTracks);

    // calculate specification from the specification of input data blocks
    AliHLTUInt32_t iSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( minSlice, maxSlice, minPatch, maxPatch );
    HLTInfo("minSlice %d, maxSlice %d, minPatch %d, maxPatch %d", minSlice, maxSlice, minPatch, maxPatch);

    // send data
    PushBack(fESD, kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC, iSpecification);

    // reset ESDs and ESDs friends
    fESD->Reset();

  } else {
    HLTError("component not initialized");
    iResult=-ENOMEM;
  }

  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::Configure(const char* arguments)
{
  // see header file for class documentation
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

      if (argument.CompareTo("-solenoidBz")==0) {
#ifndef HAVE_NOT_ALIMAGF30848
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	// TODO: Matthias 2009-05-08 this changes below look weird to me
	// maybe this was not correctly done by Federico in r 30849
	float SolenoidBz=((TObjString*)pTokens->At(i))->GetString().Atof();
	if (SolenoidBz<kAlmost0Field) SolenoidBz=kAlmost0Field;
	float factor=SolenoidBz/5;
	//
	AliMagF::BMap_t map = AliMagF::k5kG;
	if (SolenoidBz<3.) {
	  map=AliMagF::k2kG;
	  factor=SolenoidBz/2;
	} /*else if (SolenoidBz>=3. && SolenoidBz<4.5) {
	  map=AliMagFMaps::k4kG;
	  factor=SolenoidBz/4;
	  }
	else {
	  map=AliMagFMaps::k5kG;
	  factor=SolenoidBz/5;
	  } */
	// the magnetic field map is not supposed to change
	// field initialization should be done once in the beginning
	// TODO: does the factor need adjustment?
	const AliMagF* currentMap = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
	if (!currentMap) {
	  AliMagF* field = new AliMagF("MagneticFieldMap", "BMap", 2, 1., 1., 10., map);
	  TGeoGlobalMagField::Instance()->SetField(field);
	  HLTInfo("Solenoid Field set to: %f map %d", SolenoidBz, map);
	} else if (currentMap->GetMapType()!=map) {
	  HLTWarning("omitting request to override field map %d with %d", currentMap->GetMapType(), map);
	}
#else // keeping the <30489 code for backward compatibility, to be merged at some point
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	// TODO: check if there is common functionality in the AliMagF* classes
	float SolenoidBz=((TObjString*)pTokens->At(i))->GetString().Atof();
	if (SolenoidBz<kAlmost0Field) SolenoidBz=kAlmost0Field;
	float factor=1.;
	int map=AliMagFMaps::k2kG;
	if (SolenoidBz<3.) {
	  map=AliMagFMaps::k2kG;
	  factor=SolenoidBz/2;
	} else if (SolenoidBz>=3. && SolenoidBz<4.5) {
	  map=AliMagFMaps::k4kG;
	  factor=SolenoidBz/4;
	} else {
	  map=AliMagFMaps::k5kG;
	  factor=SolenoidBz/5;
	}
	// the magnetic field map is not supposed to change
	// field initialization should be done once in the beginning
	// TODO: does the factor need adjustment?
	const AliMagF* currentMap=AliTracker::GetFieldMap();
	if (!currentMap) {
	  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., map);
	  AliTracker::SetFieldMap(field,kFALSE);
	  HLTInfo("Solenoid Field set to: %f map %d", SolenoidBz, map);
	} else if (currentMap->Map()!=map) {
	  HLTWarning("omitting request to override field map %s with %s", currentMap->Map(), map);
	}
#endif
	continue;
      } else {
	HLTError("unknown argument %s", argument.Data());
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
  return iResult;
}

int AliHLTTPCOfflineTrackerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path=kAliHLTCDBSolenoidBz;
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    if (chainId) {} // just to get rid of warning, can not comment argument due to debug message
    HLTDebug("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTDebug("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  
  return iResult;
}
