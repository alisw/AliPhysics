// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Øvrebekk <st05886@alf.uib.no>                   *
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

/** @file   AliHLTITSClusterFinderSSDComponent.cxx
    @author Gaute Øvrebekk <st05886@alf.uib.no>
    @date   
    @brief  Component to run offline clusterfinder for SSD
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSClusterFinderSSDComponent.h" 

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliHLTDataTypes.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include <AliHLTDAQ.h>
#include "AliGeomManager.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSClusterFinderSSDComponent);

AliHLTITSClusterFinderSSDComponent::AliHLTITSClusterFinderSSDComponent()
  :
  fNModules(AliITSgeomTGeo::GetNDetectors(5)*AliITSgeomTGeo::GetNLadders(5) + AliITSgeomTGeo::GetNDetectors(6)*AliITSgeomTGeo::GetNLadders(6)/*1698*/),
  fClusterFinder(NULL),
  fRawReader(NULL),
  fDettype(NULL),
  fClusters(NULL),
  fgeom(NULL),
  fgeomInit(NULL),
  fSeg(NULL)
{
 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // AliITSDetTypeRec::fgkDefaultNModulesSPD private for the moment
//   if (AliITSDetTypeRec::fgkDefaultNModulesSSD!=1698) {
//     HLTWarning("Number of modules has changed (AliITSDetTypeRec::fgkDefaultNModulesSSD)");
//   }
}

AliHLTITSClusterFinderSSDComponent::~AliHLTITSClusterFinderSSDComponent() {
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSClusterFinderSSDComponent::GetComponentID()
{
  // see header file for class documentation

  return "ITSClusterFinderSSD";
}

void AliHLTITSClusterFinderSSDComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD );

}

AliHLTComponentDataType AliHLTITSClusterFinderSSDComponent::GetOutputDataType() {
  // see header file for class documentation
  return  kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD;
}

void AliHLTITSClusterFinderSSDComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation

  constBase = 0;
  inputMultiplier = 100;
}

AliHLTComponent* AliHLTITSClusterFinderSSDComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTITSClusterFinderSSDComponent();
}
	
Int_t AliHLTITSClusterFinderSSDComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  if ( fClusterFinder )
    return -EINPROGRESS;

  fClusters = new TClonesArray*[fNModules]; 
  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    fClusters[iModule] = NULL;
  }

  AliCDBManager* man = AliCDBManager::Instance();
  AliGeomManager::LoadGeometry();
 
  if (!man->IsDefaultStorageSet()){
    HLTError("Default CDB storage has not been set !");
    return -ENOENT;
  }

  //fgeomInit = new AliITSInitGeometry(kvSSD02,2);
  fgeomInit = new AliITSInitGeometry(kvPPRasymmFMD,2);
  //fgeomInit->InitAliITSgeom(fgeom);
  fgeom = fgeomInit->CreateAliITSgeom();
  
  //set dettype
  fDettype = new AliITSDetTypeRec();
  fDettype->SetITSgeom(fgeom);
  fDettype->SetReconstructionModel(2,fClusterFinder);
  fDettype->SetDefaultClusterFindersV2(kTRUE);
  fDettype->SetDefaults();
  //fSeg = new AliITSsegmentationSSD();
  //fSeg->Init();
  //fDettype->SetSegmentationModel(2,fSeg);
  //fDettype->GetCalibration();
    
  fClusterFinder = new AliITSClusterFinderV2SSD(fDettype); 
  fClusterFinder->InitGeometry();
  
  if ( fRawReader )
    return -EINPROGRESS;

  fRawReader = new AliRawReaderMemory();

  return 0;
}

Int_t AliHLTITSClusterFinderSSDComponent::DoDeinit() {
  // see header file for class documentation
    
  if ( fRawReader )
    delete fRawReader;
  fRawReader = NULL;
  
  if ( fDettype )
    delete fDettype;
  fDettype = NULL;
  
  if ( fClusterFinder )
    delete fClusterFinder;
  fClusterFinder = NULL;
  
  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    delete fClusters[iModule];
    fClusters[iModule] = NULL;
  }
  
  if ( fgeomInit )
    delete fgeomInit;
  fgeomInit = NULL;

  return 0;
}

Int_t AliHLTITSClusterFinderSSDComponent::DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/)
{  // see header file for class documentation

  // -- Iterator over Data Blocks --
  const AliHLTComponentBlockData* iter = NULL;
  
  if (!IsDataEvent()) return 0;

  if ( evtData.fBlockCnt<=0 )
      {
	Logging( kHLTLogWarning, "HLT::ITSClusterFinderSSD::DoEvent", "DoEvent", "no blocks in event" );
	return 0;
      }

  // -- Loop over blocks
  for ( iter = GetFirstInputBlock(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD); iter != NULL; iter = GetNextInputBlock() ) {
  
    // -- Debug output of datatype --
    HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, 
	       DataType2Text(iter->fDataType).c_str(), 
	       DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD).c_str());
    
    // -- Check for the correct data type
    if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSSD) )  
      continue;
    
    // -- Get equipment ID out of specification
    AliHLTUInt32_t spec = iter->fSpecification;
  
    if(spec>0x00040000){
      HLTDebug("The Spec is to high for ITS SSD");
    }

    Int_t id = AliHLTDAQ::DdlIDOffset("ITSSSD");                  
    for ( Int_t ii = 0; ii < AliHLTDAQ::NumberOfDdls("ITSSSD") ; ii++ ) {
      if ( spec & 0x00000001 ) {
	id += ii;
	break;
      }
      spec = spec >> 1 ;
    }

    // -- Set equipment ID to the raw reader
    
    if(!fRawReader->AddBuffer((UChar_t*) iter->fPtr, iter->fSize, id)){
      HLTWarning("Could not add buffer");
    }
    
    fClusterFinder->RawdataToClusters(fRawReader,fClusters);

    UInt_t nClusters=0;
    for(int i=0;i<fNModules;i++){
      if(fClusters[i] != NULL){
        nClusters += fClusters[i]->GetEntries(); 
      }
    }
    
    UInt_t bufferSize = nClusters * sizeof(AliHLTITSSpacePointData) + sizeof(AliHLTITSClusterData);
    AliHLTUInt8_t *buffer = new AliHLTUInt8_t[bufferSize];
    AliHLTITSClusterData *outputClusters = reinterpret_cast<AliHLTITSClusterData*>(buffer);
    outputClusters->fSpacePointCnt=nClusters;
    
    int clustIdx=0;
    for(int i=0;i<fNModules;i++){
      if(fClusters[i] != NULL){
        for(int j=0;j<fClusters[i]->GetEntries();j++){
          AliITSRecPoint *recpoint = (AliITSRecPoint*) fClusters[i]->At(j);
          outputClusters->fSpacePoints[clustIdx].fY=recpoint->GetY();
          outputClusters->fSpacePoints[clustIdx].fZ=recpoint->GetZ();
          outputClusters->fSpacePoints[clustIdx].fSigmaY2=recpoint->GetSigmaY2();
          outputClusters->fSpacePoints[clustIdx].fSigmaZ2=recpoint->GetSigmaZ2();
          outputClusters->fSpacePoints[clustIdx].fSigmaYZ=recpoint->GetSigmaYZ();
          outputClusters->fSpacePoints[clustIdx].fQ=recpoint->GetQ();
          outputClusters->fSpacePoints[clustIdx].fNy=recpoint->GetNy();
          outputClusters->fSpacePoints[clustIdx].fNz=recpoint->GetNz();
          outputClusters->fSpacePoints[clustIdx].fLayer=recpoint->GetLayer();
          outputClusters->fSpacePoints[clustIdx].fIndex=recpoint->GetDetectorIndex();// | recpoint->GetPindex() | recpoint->GetNindex();
          outputClusters->fSpacePoints[clustIdx].fTracks[0]=recpoint->GetLabel(0);
          outputClusters->fSpacePoints[clustIdx].fTracks[1]=recpoint->GetLabel(1);
          outputClusters->fSpacePoints[clustIdx].fTracks[2]=recpoint->GetLabel(2);
	  
          clustIdx++;
        }
      }
    }
    
    PushBack(buffer,bufferSize,kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD,iter->fSpecification);  
        
    for (Int_t iModule = 0; iModule < fNModules; iModule++) {       
      if(fClusters[iModule]){delete fClusters[iModule];}
      fClusters[iModule] = NULL;
    }
    
    delete buffer; 

    fRawReader->ClearBuffers();
    
  } //  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {    
  
  return 0;
}
