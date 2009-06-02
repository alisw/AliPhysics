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

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSClusterFinderSPDComponent);

AliHLTITSClusterFinderSPDComponent::AliHLTITSClusterFinderSPDComponent()
  :
  fNModules(240/*AliITSDetTypeRec::fgkDefaultNModulesSPD*/),
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
  return kAliHLTDataTypeTObjArray;
}

void AliHLTITSClusterFinderSPDComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation

  constBase = 0;
  inputMultiplier = 0.3;
}

AliHLTComponent* AliHLTITSClusterFinderSPDComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTITSClusterFinderSPDComponent();
}
	
Int_t AliHLTITSClusterFinderSPDComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  if ( fClusterFinder )
    return EINPROGRESS;

  fClusters = new TClonesArray*[fNModules]; 
  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    fClusters[iModule] = NULL;
  }

  //fgeomInit = new AliITSInitGeometry(kvSPD02,2);
  fgeomInit = new AliITSInitGeometry(kvPPRasymmFMD,2);
  fgeom = fgeomInit->CreateAliITSgeom();
  
  //set dettype
  fDettype = new AliITSDetTypeRec();
  fDettype->SetITSgeom(fgeom);
  fDettype->SetReconstructionModel(0,fClusterFinder);
  fDettype->SetDefaultClusterFindersV2(kTRUE);
  fDettype->GetCalibration();
  fSeg = new AliITSsegmentationSSD();
  fDettype->SetSegmentationModel(0,fSeg);
  
  fClusterFinder = new AliITSClusterFinderV2SPD(fDettype); 
  fClusterFinder->InitGeometry();

  if ( fRawReader )
    return EINPROGRESS;

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

  for (Int_t iModule = 0; iModule < fNModules; iModule++) {
    if (fClusters[iModule]) delete fClusters[iModule];
    fClusters[iModule] = NULL;
  }
  
  if ( fgeomInit )
    delete fgeomInit;
  fgeomInit = NULL;

  return 0;
}

Int_t AliHLTITSClusterFinderSPDComponent::DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/)
{  // see header file for class documentation

  // -- Iterator over Data Blocks --
  const AliHLTComponentBlockData* iter = NULL;
  
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )){
    return 0;
  }

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

    Int_t id = 0;
    for ( Int_t ii = 0; ii < 20 ; ii++ ) {   //number of ddl's
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
    
    /*
    Float_t xyz[3];
    filebuf fb;
    fb.open ("test.txt",ios::out | ios::app);
    ostream os(&fb);
    for(int i=0;i<fNModules;i++){
      if(fClusters[i] != NULL){
	for(int j=0;j<fClusters[i]->GetEntries();j++){
	  AliITSRecPoint *recpoint = (AliITSRecPoint*) fClusters[i]->At(j);
	  recpoint->GetGlobalXYZ(xyz);
	  os<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<endl;
	}
      }
    }
    */

    PushBack(*fClusters,kAliHLTDataTypeTObjArray|kAliHLTDataOriginITSSPD,iter->fSpecification);

    /*
    TClonesArray *clustersOut = new TClonesArray();
          
    for(int i=0;i<fNModules;i++){
      if(fClusters[i] != NULL){
	for(int j=0;j<fClusters[i]->GetEntriesFast();j++) {
	  clustersOut->Add(&fClusters[i][j]);
	  //PushBack(fClusters[i],kAliHLTDataTypeTObjArray|kAliHLTDataOriginITSSPD,iter->fSpecification);
	}
      }
    }
    
    PushBack(clustersOut,kAliHLTDataTypeTObjArray|kAliHLTDataOriginITSSPD,iter->fSpecification);

    Float_t xyz2[3];
    filebuf fb2;
    fb2.open ("test2.txt",ios::out | ios::app);
    ostream os2(&fb2);
    for(int j=0;j<clustersOut->GetEntriesFast();j++){
      AliITSRecPoint *recpoint = (AliITSRecPoint*) clustersOut->At(j);
      recpoint->GetGlobalXYZ(xyz2);
      os2<<xyz2[0]<<" "<<xyz2[1]<<" "<<xyz2[2]<<endl;
    }

    clustersOut->Clear();
    delete clustersOut;
    */
    for (Int_t iModule = 0; iModule < fNModules; iModule++) {  
      if(fClusters[iModule]){delete fClusters[iModule];}
      fClusters[iModule] = NULL;
    }
    
    fRawReader->ClearBuffers();
    
  } //  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {    
  
    //fClusterFinder->RawdataToClusters(fRawReader,&fClusters);
  
  //PushBack( (TObject**) fClusters,kAliHLTDataTypeTObjArray,0x00000000);

  return 0;
}
