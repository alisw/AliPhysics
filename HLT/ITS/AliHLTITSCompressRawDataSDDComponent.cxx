// $Id: AliHLTTPCClusterFinderComponent.cxx 27006 2008-07-01 09:21:45Z richterm $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTITSCompressRawDataSDDComponent.cxx
    @author Jochen Thaeder <thaeder@kip.uni-heidelberg.de>
    @date   
    @brief  Component to run data compression for SDD
*/

#if __GNUC__>= 3
using namespace std;
#endif


#include "AliHLTITSCompressRawDataSDDComponent.h" 

#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include <sys/time.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSCompressRawDataSDDComponent);

AliHLTITSCompressRawDataSDDComponent::AliHLTITSCompressRawDataSDDComponent()
  :
  fDataCompressor(NULL),
  fRawReader(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTITSCompressRawDataSDDComponent::~AliHLTITSCompressRawDataSDDComponent() {
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSCompressRawDataSDDComponent::GetComponentID()
{
  // see header file for class documentation

  return "ITSDataCompressorSDD";
}

void AliHLTITSCompressRawDataSDDComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITS );

}

AliHLTComponentDataType AliHLTITSCompressRawDataSDDComponent::GetOutputDataType() {
  // see header file for class documentation
  return (kAliHLTDataTypeDDLRaw| kAliHLTDataOriginITS);
}

void AliHLTITSCompressRawDataSDDComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation

  constBase = 0;
  inputMultiplier = 0.05;
}

AliHLTComponent* AliHLTITSCompressRawDataSDDComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTITSCompressRawDataSDDComponent();
}
	
Int_t AliHLTITSCompressRawDataSDDComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  if ( fDataCompressor )
    return EINPROGRESS;

  fDataCompressor = new AliITSCompressRawDataSDD();

  if ( fRawReader )
    return EINPROGRESS;

  fRawReader = new AliRawReaderMemory();

  return 0;
}

Int_t AliHLTITSCompressRawDataSDDComponent::DoDeinit() {
  // see header file for class documentation

  if ( fRawReader )
    delete fRawReader;
  fRawReader = NULL;

  if ( fDataCompressor )
    delete fDataCompressor;
  fDataCompressor = NULL;

  return 0;
}

Int_t AliHLTITSCompressRawDataSDDComponent::DoEvent( const AliHLTComponentEventData& evtData, 
						     const AliHLTComponentBlockData* blocks, 
						     AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
						     AliHLTUInt32_t& size, 
						     vector<AliHLTComponentBlockData>& outputBlocks ) {
  // see header file for class documentation

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )){
    size=0;
    return 0;
  }

  // -- Iterator over Data Blocks --
  const AliHLTComponentBlockData* iter = NULL;

  // -- Ptr to output shm region --
  AliHLTUInt8_t* outShmPtr = outputPtr;

  // -- Initialize out sizes
  UInt_t offset    = 0;   // offset of current outblock
  UInt_t mySize    = 0;   // out size produced from current block
  UInt_t totalSize = 0;   // total out size of this event
  UInt_t availSize = 0;   // still availible out size for this event

  // -- Loop over blocks
  for ( ULong_t ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {

    iter = blocks+ndx;

    mySize = 0;
    offset = totalSize;
    availSize = size - totalSize;

    // -- Debug output of datatype --
    HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
	     evtData.fEventID, evtData.fEventID, 
	     DataType2Text(iter->fDataType).c_str(), 
	     DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITS).c_str());
    
    // -- Check for the correct data type
    if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginITS) ) { 
      HLTError("WRONG FORMAT");
      continue;
    }
    
    // -- Set RawReader
    fRawReader->SetMemory( (UChar_t*) iter->fPtr, iter->fSize );

    // -- Get equipment ID out of specification
    AliHLTUInt32_t spec = iter->fSpecification;

    Int_t id = 256;
    for ( Int_t ii = 0; ii < 24 ; ii++ ) {
      if ( spec & 0x00000001 ) {
	id += ii;
	break;
      }
      spec = spec >> 1 ;
    }
    
    // -- Set equipment ID to the raw reader
    fRawReader->SetEquipmentID( id ); 
    
    // -- Set raw reader
    fDataCompressor->SetRawReader( fRawReader );

    // -- Set ptr to output shm
    fDataCompressor->SetPointerToData( (UChar_t*) outShmPtr );

    // -- Set availible outputspace
    fDataCompressor->SetSize( (UInt_t) availSize );

    // -- Compress event
    mySize = fDataCompressor->CompressEvent( (UChar_t*) iter->fPtr );

    // -- Fill output blocks
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset = offset;
    bd.fSize = mySize;
    bd.fSpecification = iter->fSpecification;
    bd.fDataType = iter->fDataType;
    outputBlocks.push_back( bd );

    // -- Increase size counters
    totalSize += mySize;

    // -- Increase output shm ptr
    outShmPtr += mySize;
    
    // -- Check if data was written over allowed buffer
    if ( totalSize > size ) {
      HLTFatal( "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",  totalSize, size );
      return EMSGSIZE;
    }
    
  } //  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ ) {    

  // -- Set total output size
  size = totalSize;
    
  return 0;
}
