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

/** @file   AliHLTITSClusterFinderSPDComponent.cxx
    @author Gaute Øvrebekk <st05886@alf.uib.no>
    @date   
    @brief  Component to run offline clusterfinder for SPD
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSClusterFinderSPDComponent.h" 

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliHLTDataTypes.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliHLTITSSpacePointData.h"
#include "AliHLTITSClusterDataFormat.h"
#include <AliHLTDAQ.h>
#include "AliGeomManager.h"
#include "TTree.h"
#include "TBranch.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSClusterFinderSPDComponent);

AliHLTITSClusterFinderSPDComponent::AliHLTITSClusterFinderSPDComponent()
  :
  fNModules(AliITSgeomTGeo::GetNDetectors(1)*AliITSgeomTGeo::GetNLadders(1) + AliITSgeomTGeo::GetNDetectors(2)*AliITSgeomTGeo::GetNLadders(2)/*240*/),
  fClusterFinder(NULL),
  fRawReader(NULL),
  fDettype(NULL),
  fgeom(NULL),
  fgeomInit(NULL)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // AliITSDetTypeRec::fgkDefaultNModulesSPD private for the moment
//   if (AliITSDetTypeRec::fgkDefaultNModulesSPD!=240) {
//     HLTWarning("Number of modules has changed (AliITSDetTypeRec::fgkDefaultNModulesSPD)");
//   }
}

AliHLTITSClusterFinderSPDComponent::~AliHLTITSClusterFinderSPDComponent() {
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSClusterFinderSPDComponent::GetComponentID()
{
  // see header file for class documentation

  return "ITSClusterFinderSPD";
}

void AliHLTITSClusterFinderSPDComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD );

}

AliHLTComponentDataType AliHLTITSClusterFinderSPDComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD;
}

void AliHLTITSClusterFinderSPDComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation

  constBase = 0;
  inputMultiplier = 100;
}

AliHLTComponent* AliHLTITSClusterFinderSPDComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTITSClusterFinderSPDComponent();
}
	
Int_t AliHLTITSClusterFinderSPDComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  if ( fClusterFinder )
    return -EINPROGRESS;
  
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()){
    HLTError("Default CDB storage has not been set !");
    return -ENOENT;
  }

  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }
  //fgeomInit = new AliITSInitGeometry(kvSPD02,2);
  fgeomInit = new AliITSInitGeometry(kvPPRasymmFMD,2);
  //fgeomInit->InitAliITSgeom(fgeom);
  fgeom = fgeomInit->CreateAliITSgeom();
  
  //set dettype
  fDettype = new AliITSDetTypeRec();
  fDettype->SetITSgeom(fgeom);
  fDettype->SetReconstructionModel(0,fClusterFinder);
  fDettype->SetDefaultClusterFindersV2(kTRUE);
  fDettype->SetDefaults();
  
  fClusterFinder = new AliITSClusterFinderV2SPD(fDettype); 
  fClusterFinder->InitGeometry();

  if ( fRawReader )
    return -EINPROGRESS;

  fRawReader = new AliRawReaderMemory();

  return 0;
}

Int_t AliHLTITSClusterFinderSPDComponent::DoDeinit() {
  // see header file for class documentation

  if ( fRawReader )
    delete fRawReader;
  fRawReader = NULL;

  if ( fClusterFinder )
    delete fClusterFinder;
  fClusterFinder = NULL;

  if ( fDettype )
    delete fDettype;
  fDettype = NULL;

  if ( fgeomInit )
    delete fgeomInit;
  fgeomInit = NULL;

  return 0;
}

Int_t AliHLTITSClusterFinderSPDComponent::DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/)
{  // see header file for class documentation

  // -- Iterator over Data Blocks --
  const AliHLTComponentBlockData* iter = NULL;
  
  if (!IsDataEvent()) return 0;

  if ( evtData.fBlockCnt<=0 )
      {
	Logging( kHLTLogWarning, "HLT::ITSClusterFinderSPD::DoEvent", "DoEvent", "no blocks in event" );
	return 0;
      }

  // -- Loop over blocks
  for ( iter = GetFirstInputBlock(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD); iter != NULL; iter = GetNextInputBlock() ) {
  
    // -- Debug output of datatype --
    HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, 
	       DataType2Text(iter->fDataType).c_str(), 
	       DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD).c_str());
    
    // -- Check for the correct data type
    if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITSSPD) )  
      continue;
    
    // -- Get equipment ID out of specification
    AliHLTUInt32_t spec = iter->fSpecification;
  
    if(spec>0x00040000){
      HLTDebug("The Spec is to high for ITS SPD");
    }
    
    Int_t id = AliHLTDAQ::DdlIDOffset("ITSSPD");
    for ( Int_t ii = 0; ii < AliHLTDAQ::NumberOfDdls("ITSSPD") ; ii++ ) {   //number of ddl's
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
    
    //fClusterFinder->RawdataToClusters(fRawReader,fClusters);
    TTree *tree = new TTree();
    fDettype->DigitsToRecPoints(fRawReader,tree,"SPD");

    UInt_t nClusters=0;
    TClonesArray *array=new TClonesArray("AliITSRecPoint",1000);
    TBranch *branch = tree->GetBranch("ITSRecPoints");
    branch->SetAddress(&array);
    for(int ev=0;ev<branch->GetEntries();ev++){
      branch->GetEntry(ev);
      nClusters += array->GetEntries();
    }
    
    UInt_t bufferSize = nClusters * sizeof(AliHLTITSSpacePointData) + sizeof(AliHLTITSClusterData);
    AliHLTUInt8_t *buffer = new AliHLTUInt8_t[bufferSize];
    AliHLTITSClusterData *outputClusters = reinterpret_cast<AliHLTITSClusterData*>(buffer);
    outputClusters->fSpacePointCnt=nClusters;

    int clustIdx=0;
    for(int i=0;i<branch->GetEntries();i++){
      branch->GetEntry(i);
      for(int j=0;j<array->GetEntries();j++){
	AliITSRecPoint *recpoint = (AliITSRecPoint*) array->At(j);
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
    
    PushBack(buffer,bufferSize,kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD,iter->fSpecification);

    array->Delete();
    delete array;
    delete tree;
    delete buffer; 
    fRawReader->ClearBuffers();
    
  } //  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {    
  
  return 0;
}
 
