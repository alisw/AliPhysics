// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kenneth Aamodt                                        *
 *                  Kalliopi Kanaki                                       *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCKryptonClusterFinderComponent.cxx
    @author Kenneth Aamodt, Kalliopi Kanaki
    @date   
    @brief  The TPC krypton cluster finder processing component
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCKryptonClusterFinderComponent.h"
#include "AliHLTTPCDigitReaderDecoder.h"
#include "AliHLTTPCKryptonClusterFinder.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusters.h"
#include "AliHLTTPCDefinitions.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TH1F.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

AliHLTTPCKryptonClusterFinderComponent gAliHLTTPCKryptonClusterFinderComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCKryptonClusterFinderComponent)

AliHLTTPCKryptonClusterFinderComponent::AliHLTTPCKryptonClusterFinderComponent()
  :
  fKryptonClusterFinder(NULL),
  fReader(NULL),
  fSpecification(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCKryptonClusterFinderComponent::~AliHLTTPCKryptonClusterFinderComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCKryptonClusterFinderComponent::GetComponentID()
{
  // see header file for class documentation
    return "TPCKryptonClusterFinder";
}

void AliHLTTPCKryptonClusterFinderComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCKryptonClusterFinderComponent::GetOutputDataType()
{
  // see header file for class documentation
 return kAliHLTMultipleDataType;
 //  return kAliHLTDataTypeHistogram;
}

int AliHLTTPCKryptonClusterFinderComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  tgtList.push_back(kAliHLTDataTypeHwAddr16);
  return tgtList.size();
}

void AliHLTTPCKryptonClusterFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.  
  constBase = 0;
  inputMultiplier = (100 * 0.4);
}

AliHLTComponent* AliHLTTPCKryptonClusterFinderComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCKryptonClusterFinderComponent();
}
	
int AliHLTTPCKryptonClusterFinderComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  //  HLTFatal("Initializing the kryptonclusterfindercomponent");

  if ( fKryptonClusterFinder )
    return EINPROGRESS;

  fKryptonClusterFinder = new AliHLTTPCKryptonClusterFinder();

  Int_t i = 0;

  while ( i < argc ) {      

    //No arguments so far

    Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  fReader = new AliHLTTPCDigitReaderDecoder();
  fKryptonClusterFinder->SetReader(fReader);
 
  return 0;
}

int AliHLTTPCKryptonClusterFinderComponent::DoDeinit()
{
  // see header file for class documentation

  if ( fKryptonClusterFinder )
    delete fKryptonClusterFinder;
  fKryptonClusterFinder = NULL;
 
  if ( fReader )
    delete fReader;
  fReader = NULL;
    
  return 0;
}

int AliHLTTPCKryptonClusterFinderComponent::DoEvent( const AliHLTComponentEventData& evtData, 
						     const AliHLTComponentBlockData* /*blocks*/, 
					      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, 
					      vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )){
    size=0;
    return 0;
  }

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;

  Int_t slice, patch;

  unsigned long maxPoints = 0;

  AliHLTTPCClusterData* outPtr;

  outPtr = (AliHLTTPCClusterData*)outputPtr;

  for (iter = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC); iter != NULL; iter = GetNextInputBlock()){

    HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
	    evtData.fEventID, evtData.fEventID,
	    DataType2Text(iter->fDataType).c_str(), 
	    DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());

    if (iter->fDataType == AliHLTTPCDefinitions::fgkDDLPackedRawDataType && GetEventCount()<2){
      HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeDDLRaw)!", 
		 DataType2Text(AliHLTTPCDefinitions::fgkDDLPackedRawDataType).c_str(),
		 DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
    }      
     
    if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC) &&
	 iter->fDataType != AliHLTTPCDefinitions::fgkDDLPackedRawDataType ) continue;
      
    fSpecification = iter->fSpecification;
      
    slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );


    fKryptonClusterFinder->SetPatch(patch);

    fKryptonClusterFinder->InitSlice( slice, patch, maxPoints );

    fKryptonClusterFinder->SetOutputArray((AliHLTTPCSpacePointData*)outPtr->fSpacePoints);
    
    fKryptonClusterFinder->SetMaxOutputSize(size-sizeof(AliHLTUInt32_t));

    fKryptonClusterFinder->ReadDataUnsorted(iter->fPtr, iter->fSize);
      
    fKryptonClusterFinder->FindRowClusters();

    fKryptonClusterFinder->FindKryptonClusters();
    
    //    if(fKryptonClusterFinder->GetNKryptonClusters()>0){

    //    HLTFatal("Number of kryptonClusters: %d",(Int_t)fKryptonClusterFinder->GetNKryptonClusters());

      if (sizeof(AliHLTUInt32_t)+fKryptonClusterFinder->GetNKryptonClusters()*sizeof(AliHLTTPCSpacePointData)>=size) {
	HLTFatal("Buffer too small too add more spacepoints: %d of %d byte(s) already used",sizeof(AliHLTUInt32_t)+fKryptonClusterFinder->GetNKryptonClusters()*sizeof(AliHLTTPCSpacePointData) ,size);
	return -ENOSPC;
      }
      
      outPtr->fSpacePointCnt=fKryptonClusterFinder->GetNKryptonClusters();
      
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = 0;
      bd.fSize = sizeof(AliHLTUInt32_t)+fKryptonClusterFinder->GetNKryptonClusters()*sizeof(AliHLTTPCSpacePointData);
      bd.fSpecification = iter->fSpecification;
      bd.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
      outputBlocks.push_back( bd );

      //adding the list of hardware addresses to the outpot buffer

      AliHLTUInt32_t dataOffsetBeforeHW=sizeof(AliHLTUInt32_t)+fKryptonClusterFinder->GetNKryptonClusters()*sizeof(AliHLTTPCSpacePointData);
      AliHLTUInt32_t sizeOfHWArray=fKryptonClusterFinder->fHWAddressVector.size()*sizeof(AliHLTUInt16_t);
      //      cout<<"size of array: "<<sizeOfHWArray<<"     number of entries"<<fKryptonClusterFinder->fHWAddressVector.size()<<endl;

      if (dataOffsetBeforeHW+sizeOfHWArray>=size){
	HLTFatal("Buffer too small too add the HW address list: %d of %d byte(s) already used",dataOffsetBeforeHW+sizeOfHWArray ,size);
	return -ENOSPC;
      }

      AliHLTUInt16_t *outputHWPtr=(AliHLTUInt16_t*)(outputPtr+dataOffsetBeforeHW);
      for(UInt_t hw=0;hw<fKryptonClusterFinder->fHWAddressVector.size();hw++){
	*outputHWPtr = fKryptonClusterFinder->fHWAddressVector[hw];
	outputHWPtr++;
      }
      
      AliHLTComponentBlockData bdHW;
      FillBlockData( bdHW );
      bdHW.fOffset =  dataOffsetBeforeHW;
      bdHW.fSize = sizeOfHWArray;
      bdHW.fSpecification = iter->fSpecification;
      bdHW.fDataType = kAliHLTDataTypeHwAddr16;
      outputBlocks.push_back( bdHW );

      size=dataOffsetBeforeHW+sizeOfHWArray;
      
      /*    }
	    else{
	    size=0;
      }*/
      fReader->Reset();
  }
  return 0;
}

int AliHLTTPCKryptonClusterFinderComponent::Configure(const char* arguments) { 
// see header file for class documentation
  
  int iResult=0;
  if (!arguments) return iResult;
  HLTInfo("parsing configuration string \'%s\'", arguments);
  /*
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
     
      if (argument.CompareTo("-apply-noisemap")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-apply-noisemap\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-plot-side-c")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-plot-side-c\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-plot-side-a")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-plot-side-a\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    } // end for
  
    delete pTokens;
  
  } // end if pTokens
  
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  */
  return iResult;
}

int AliHLTTPCKryptonClusterFinderComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  const char* path="HLT/ConfigTPC";
  if (cdbEntry) path=cdbEntry;
  if (path) {
    HLTInfo("reconfigure from entry %s, chain id %s", path, (chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object: %s", pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }
  return 0;
}
