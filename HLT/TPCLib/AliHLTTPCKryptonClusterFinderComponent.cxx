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
  fClusterFinder(NULL),
  fReader(NULL)
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
  return kAliHLTDataTypeHistogram;
}

int AliHLTTPCKryptonClusterFinderComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  return tgtList.size();
}

void AliHLTTPCKryptonClusterFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.  
  constBase = 0;
  inputMultiplier = (6 * 0.4);
}

AliHLTComponent* AliHLTTPCKryptonClusterFinderComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCKryptonClusterFinderComponent();
}
	
int AliHLTTPCKryptonClusterFinderComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  if ( fClusterFinder )
    return EINPROGRESS;

  fClusterFinder = new AliHLTTPCKryptonClusterFinder();

  Int_t i = 0;

  while ( i < argc ) {      

    //No arguments so far

    Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  fReader = new AliHLTTPCDigitReaderDecoder();
  fClusterFinder->SetReader(fReader);
 
  return 0;
}

int AliHLTTPCKryptonClusterFinderComponent::DoDeinit()
{
  // see header file for class documentation

  if ( fClusterFinder )
    delete fClusterFinder;
  fClusterFinder = NULL;
 
  if ( fReader )
    delete fReader;
  fReader = NULL;
    
  return 0;
}

int AliHLTTPCKryptonClusterFinderComponent::DoEvent( const AliHLTComponentEventData& evtData, 
						     const AliHLTComponentBlockData* blocks, 
						     AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
						     AliHLTUInt32_t& size, 
						     vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  //  == OUTdatatype pointer
  AliHLTTPCClusterData* outPtr;

  AliHLTUInt8_t* outBPtr;
  UInt_t offset, mysize, nSize, tSize = 0;

  outBPtr = outputPtr;
  outPtr = (AliHLTTPCClusterData*)outBPtr;

  Int_t slice, patch, row[2];
  unsigned long maxPoints, realPoints = 0;

  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;


      HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	       evtData.fEventID, evtData.fEventID, 
	       DataType2Text( iter->fDataType).c_str(), 
	       DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
      
      if (iter->fDataType == AliHLTTPCDefinitions::fgkDDLPackedRawDataType &&
	  GetEventCount()<2) {
	HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeDDLRaw)!",
		   DataType2Text(AliHLTTPCDefinitions::fgkDDLPackedRawDataType).c_str(),
		   DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
      }
      
      if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC) &&
	   iter->fDataType != AliHLTTPCDefinitions::fgkDDLPackedRawDataType ) continue;
      
      
      slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
      row[0] = AliHLTTPCTransform::GetFirstRow( patch );
      row[1] = AliHLTTPCTransform::GetLastRow( patch );


      fClusterFinder->SetPatch(patch);

      outPtr = (AliHLTTPCClusterData*)outBPtr;

      maxPoints = (size-tSize-sizeof(AliHLTTPCClusterData))/sizeof(AliHLTTPCSpacePointData);

      fClusterFinder->InitSlice( slice, patch, row[0], row[1], maxPoints );
	
      fClusterFinder->ReadDataUnsorted(iter->fPtr, iter->fSize);
      
      fClusterFinder->FindNormalClusters();

      fClusterFinder->FindKryptonClusters();

      realPoints = fClusterFinder->GetNumberOfClusters();
	
    
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fSpecification = iter->fSpecification;
      bd.fDataType = AliHLTTPCDefinitions::fgkClustersDataType;
      //AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
      outputBlocks.push_back( bd );
  
      tSize += mysize;
      outBPtr += mysize;
      outPtr = (AliHLTTPCClusterData*)outBPtr;
  
  
      if ( tSize > size )
	{
	  Logging( kHLTLogFatal, "HLT::TPCClusterFinder::DoEvent", "Too much data", 
		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
		   tSize, size );
	  return EMSGSIZE;
	}
  
      size = tSize;
    }
  return 0;
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
