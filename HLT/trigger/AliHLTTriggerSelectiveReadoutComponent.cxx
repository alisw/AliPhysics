//-*- Mode: C++ -*-
// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTTriggerSelectiveReadoutComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  Component for the Selective Readout Trigger
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTriggerSelectiveReadoutComponent.h"
#include "AliHLTTPCDefinitions.h"

#include <cstdlib>
#include <cerrno>


// ** This is a global object used for automatic component registration, do not use this
AliHLTTriggerSelectiveReadoutComponent gAliHLTTriggerSelectivereadoutComponent;

ClassImp(AliHLTTriggerSelectiveReadoutComponent)
    
AliHLTTriggerSelectiveReadoutComponent::AliHLTTriggerSelectiveReadoutComponent() :
  fDetector(),
  fEnableThresholdSize(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  memset( fThreshold, 0, (fkNThreshold * sizeof(AliHLTUInt32_t)) );
}

AliHLTTriggerSelectiveReadoutComponent::~AliHLTTriggerSelectiveReadoutComponent() {
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTriggerSelectiveReadoutComponent::GetComponentID() {
  // see header file for class documentation
  return "TriggerSelectiveReadout"; 
}

void AliHLTTriggerSelectiveReadoutComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) {
  // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw );
}

AliHLTComponentDataType AliHLTTriggerSelectiveReadoutComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTMultipleDataType;
}

Int_t AliHLTTriggerSelectiveReadoutComponent::GetOutputDataTypes( AliHLTComponentDataTypeList& tgtList) {
  // see header file for class documentation

  tgtList.clear(); 
  tgtList.push_back( kAliHLTDataTypeDDLRaw );
  tgtList.push_back( kAliHLTDataTypeDDL );
}

void AliHLTTriggerSelectiveReadoutComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation

  constBase = sizeof( AliHLTEventDDL );

  if ( ! fDetector.CompareTo("TPC") )
    inputMultiplier = 216.0;
  else
    HLTError("Output DataSize for Detector '%s' has not been implemented", fDetector.Data() );
}

// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTriggerSelectiveReadoutComponent::Spawn() {
  // see header file for class documentation

  return new AliHLTTriggerSelectiveReadoutComponent;
}

Int_t AliHLTTriggerSelectiveReadoutComponent::DoInit( int argc, const char** argv ) {
  // see header file for class documentation
    
  Int_t iResult = 0;
  TString argument = "";
  TString parameter = "";
  Int_t bMissingParam=0;
  
  for ( Int_t ii=0; ii<argc && iResult>=0; ii++ ) {
  
    argument = argv[ii];

    if ( argument.IsNull() ) continue;

    // -detector
    if ( argument.CompareTo("-detector") == 0) {

      if ( ( bMissingParam=( ++ii >= argc ) ) ) break;
      
      fDetector = argv[ii];  
      fDetector.Remove(TString::kLeading, ' ');

      fDetector = parameter;
      HLTInfo( "Detector has been set to '%s'.", fDetector.Data() );
    }
    
    // -enableThreshold
    else if ( ! argument.CompareTo("-enableThreshold") ) {

      if ( ( bMissingParam=( ++ii >= argc ) ) ) break;
      
      parameter = argv[ii];  
      parameter.Remove(TString::kLeading, ' ');

      if ( ! parameter.CompareTo("size") ) {
	fEnableThresholdSize = kTRUE;
	HLTInfo( "Threshold on data 'size' has been enabled." );
      }
      else {
	HLTError( "This threshold '%s' has not been implemented.", parameter.Data() );
	iResult = -EPROTO;
      }
    } 

    // -enableThreshold
    else if ( ! argument.CompareTo("-treshold") ) {
       HLTError("Not yet implemented '%s'", argument.Data() );
    } 
    
    // - unknow parameter
    else {
      iResult = -EINVAL;
      HLTError("Unknown argument '%s'", argument.Data() );
    }

  } // for ( Int_t ii=0; ii<argc && iResult>=0; ii++ ) {
    
  if ( bMissingParam ) {
    HLTError( "Missing parameter for argument '%s'.", argument.Data() );
    iResult = -EPROTO;
  }

  return iResult;
}
    

Int_t AliHLTTriggerSelectiveReadoutComponent::DoDeinit() {
  // see header file for class documentation

  return 0;
}

Int_t AliHLTTriggerSelectiveReadoutComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData ) {
  // see header file for class documentation

  // ** Create empty Readoutlist
  AliHLTEventDDL readoutList;
  memset( &readoutList, 0, sizeof(AliHLTEventDDL) );

  const AliHLTComponentBlockData* iter = NULL;
  Int_t ddlId = 0;
  AliHLTUInt8_t patch = 0;

  // ** Loop over all input blocks and specify which data format should be read - only select Raw Data
  iter = GetFirstInputBlock( kAliHLTDataTypeDDLRaw );
  
  while ( iter != NULL ) {
    
    // ** Check if block has the origin as specified in fDetector
    if ( fDetector.CompareTo( iter->fDataType.fOrigin ) == 0 ) {
      iter = GetNextInputBlock();
      continue;
    }

    // ** Check if block has only the CDH Header ( 32 Bytes = gkAliHLTCommonHeaderCount *sizeof(AliHLTUInt32_t) )
    if ( iter->fSize <= ( gkAliHLTCommonHeaderCount * sizeof(AliHLTUInt32_t) ) ) {
      iter = GetNextInputBlock();
      continue;
    }

    //
    // ** Check for detector specifications
    //

    // ** Check for TPC specifc specification
    if ( fDetector.CompareTo("TPC") == 0 ) {
      
      // ** Get DDL ID
      AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
	
      if (patch < 2) ddlId = 768 + (2 * slice) + patch;
      else ddlId = 838 + (4*slice) + patch;
	
      HLTDebug ( "Input Data - TPC - Slice/Patch: %d/%d - EquipmentID : %d.", slice, patch, ddlId );
    } // if ( fDetector.CompareTo("TPC") == 0 ) {
    
    // ** Check for other detector specifc specification
    else {
      HLTError("Detector '%s' has not been implemented yet.", fDetector.Data() );
      iter = GetNextInputBlock();
      continue;
    }
      
    //
    // ** Check if threshold should be considerd
    //
    
    if ( fEnableThresholdSize ) {

      // ** Check for TPC threshold
      if ( ! fDetector.CompareTo("TPC") ) {
	if ( patch < fkNThreshold ){ 

	  // ** if datablock to big, do not put it in HLT readout, enable for DAQ readout
	  if ( ( fThreshold[patch] !=0 ) && ( iter->fSize > fThreshold[patch] ) ) {
	    //----------------------SetEnableDDL( &readoutList, ddlId );
	    iter = GetNextInputBlock();
	    continue;
	  }
	  
	}
      } // if ( ! fDetector.CompareTo("TPC") ) {
      
      // ** Check for other detector threshold
      else
	HLTError("Check for size threshold as not been implemented yet for Detector '%s'.", fDetector.Data() );
      
    } // if ( fEnableThresholdSize ) {

    // ** PushBack data block
    PushBack( iter->fPtr , iter->fSize, iter->fDataType, iter->fSpecification);
      
    // ** Get next input block, with the same specification as defined in GetFirstInputBlock()
    iter = GetNextInputBlock();
    
  } //  while ( iter != NULL ) {
    
  // ** PushBack readout list
  PushBack( &readoutList, sizeof(AliHLTEventDDL), kAliHLTDataTypeDDL, (AliHLTUInt32_t) 0 );
  
  return 0;
}
