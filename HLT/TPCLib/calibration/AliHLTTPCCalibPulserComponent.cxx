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

/** @file   AliHLTTPCCalibPulserComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  A pulser calibration component for the TPC.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTransform.h"

#include "AliHLTTPCCalibPulserComponent.h"

#include "AliRawDataHeader.h"
#include "AliRawReaderMemory.h"
#include "AliTPCRawStream.h"

#ifndef HAVE_NOT_ALITPCCALIBPULSER
#include "AliTPCCalibPulser.h"
#endif // HAVE_NOT_ALITPCCALIBPULSER

#include <stdlib.h>
#include <errno.h>
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCalibPulserComponent)

AliHLTTPCCalibPulserComponent::AliHLTTPCCalibPulserComponent()
  :
  fRawReader(NULL),
  fRawStream(NULL),
  fCalibPulser(NULL),
  // Note: initialization of min and max seems to be in the wrong order but is on
  // purpose in order to work correctly with the conditional in DoEvent line 233
  //  if ( patch < fMinPatch ) fMinPatch =  patch;
  //  if ( patch > fMaxPatch ) fMaxPatch =  patch;
  fMinPatch(5),
  fMaxPatch(0),
  fSpecification(0) ,
  fEnableAnalysis(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCCalibPulserComponent::~AliHLTTPCCalibPulserComponent() {
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCCalibPulserComponent::GetComponentID() {
  // see header file for class documentation

  return "TPCCalibPulser";
}

void AliHLTTPCCalibPulserComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation

  list.clear(); 
  list.push_back( AliHLTTPCDefinitions::fgkDDLPackedRawDataType );
}

AliHLTComponentDataType AliHLTTPCCalibPulserComponent::GetOutputDataType() {
  // see header file for class documentation

  return AliHLTTPCDefinitions::fgkCalibPulserDataType;
}

void AliHLTTPCCalibPulserComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation

  // XXX TODO: Find more realistic values.  
  constBase = 0;
  inputMultiplier = (2.0);
}

AliHLTComponent* AliHLTTPCCalibPulserComponent::Spawn() {
  // see header file for class documentation

  return new AliHLTTPCCalibPulserComponent();
}  


Int_t AliHLTTPCCalibPulserComponent::ScanArgument( Int_t argc, const char** argv ) {
  // see header file for class documentation

  Int_t iResult = 0;
  TString argument = "";
  TString parameter = "";

  if ( !argc ) 
    return -EINVAL;

  argument = argv[iResult];
  
  if ( argument.IsNull() ) 
    return -EINVAL;

  // -rcuformat
  if ( argument.CompareTo("-rcuformat") == 0 ) {

    if ( ++iResult >= argc  ) {
      iResult = -EPROTO;
    }
    else {
      parameter = argv[1];
      if ( parameter.CompareTo("old") == 0 ) {
        HLTInfo( "RCU Format is set to old." );
      }
      else if ( parameter.CompareTo("new") == 0 ) {
        HLTInfo( "RCU Format is set to new." );
      }
      else {
	HLTError( "Cannot convert rcu format specifier '%s'.", argv[1] );
	iResult = -EPROTO;
      }
    } 
  }
  else if ( argument.CompareTo("-enableanalysis") == 0 ) {
    HLTInfo( "Analysis before shipping data to FXS enabled." );
    fEnableAnalysis = kTRUE;
  }
  else {
    iResult = -EINVAL;
  }
      
  return iResult;
}

Int_t AliHLTTPCCalibPulserComponent::InitCalibration() {
  // see header file for class documentation
    
#ifndef HAVE_NOT_ALITPCCALIBPULSER
  // ** Create pulser calibration
  if ( fCalibPulser )
    return EINPROGRESS;
  
  fCalibPulser = new AliTPCCalibPulser();

  // **  Create AliRoot Memory Reader
  if (fRawReader)
    return EINPROGRESS;

  fRawReader = new AliRawReaderMemory();

  return 0;
#else // HAVE_NOT_ALITPCCALIBPULSER
#warning AliTPCCalibPulser not available in this AliRoot version - AliHLTTPCCalibPulserComponent not functional
  HLTFatal("AliTPCCalibPulser  not available - check your build");
  return -ENODEV;
#endif //HAVE_NOT_ALITPCCALIBPULSER
}

Int_t AliHLTTPCCalibPulserComponent::DeinitCalibration() {
  // see header file for class documentation

  if ( fRawReader )
    delete fRawReader;
  fRawReader = NULL;

#ifndef HAVE_NOT_ALITPCCALIBPULSER
  if ( fCalibPulser ) {
    delete fCalibPulser;
  }
#endif
  fCalibPulser = NULL;

  return 0;
}

Int_t AliHLTTPCCalibPulserComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/, 
							 AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation
  
  const AliHLTComponentBlockData* iter = NULL;

  AliHLTUInt8_t slice=0, patch=0;
  Int_t ddlId = 0;
    
  // ** Loop over all input blocks and specify which data format should be read - only select Raw Data
  iter = GetFirstInputBlock( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC);
  
  while ( iter != NULL ) {
    
    // ** Print Debug output which data format was received
    //    char tmp1[14], tmp2[14];
    //    DataType2Text( iter->fDataType, tmp1 );
    //    DataType2Text( AliHLTTPCDefinitions::fgkDDLPackedRawDataType, tmp2 );
    //     HLTDebug ( "Event received - Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
    // 	       evtData.fEventID, evtData.fEventID, tmp1, tmp2 );

    // ** Get DDL ID in order to tell the memory reader which slice/patch to use
    slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

    if (patch < 2) ddlId = 768 + (2 * slice) + patch;
    else ddlId = 838 + (4*slice) + patch;

    HLTDebug ( "Input Raw Data - Slice/Patch: %d/%d - EquipmentID : %d.", slice, patch, ddlId );

    // ** Get min and max patch, used for output specification
    if ( patch < fMinPatch ) fMinPatch =  patch;
    if ( patch > fMaxPatch ) fMaxPatch =  patch;

    // ** Init TPCRawStream
    fRawReader->SetMemory( reinterpret_cast<UChar_t*>( iter->fPtr ), iter->fSize );
    fRawReader->SetEquipmentID(ddlId);

    fRawStream = new AliTPCRawStream( fRawReader );

#ifndef HAVE_NOT_ALITPCCALIBPULSER
    // ** Process actual Pulser Calibration - Fill histograms
    fCalibPulser->ProcessEvent( fRawStream );
#else //!HAVE_NOT_ALITPCCALIBPULSER
    HLTFatal("AliTPCCalibPulser  not available - check your build");
    return -ENODEV;
#endif //HAVE_NOT_ALITPCCALIBPULSER
  
    // ** Delete TPCRawStream
    if ( fRawStream )
      delete fRawStream;
    fRawStream = NULL;    

    // ** Get next input block, with the same specification as defined in GetFirstInputBlock()
    iter = GetNextInputBlock();

  } //  while ( iter != NULL ) {
    
  // ** Get output specification
  fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( slice, slice, fMinPatch, fMaxPatch );

  // ** PushBack data to shared memory ... 
  PushBack( (TObject*) fCalibPulser, AliHLTTPCDefinitions::fgkCalibPulserDataType, fSpecification);
  
  return 0;
} // Int_t AliHLTTPCCalibPulserComponent::ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData ) {


Int_t AliHLTTPCCalibPulserComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, 
						    AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation
    
#ifndef HAVE_NOT_ALITPCCALIBPULSER
  if ( fEnableAnalysis )
    fCalibPulser->Analyse();
#else //!HAVE_NOT_ALITPCCALIBPULSER
  HLTFatal("AliTPCCalibPulser  not available - check your build");
  return -ENODEV;
#endif //HAVE_NOT_ALITPCCALIBPULSER
  
  // ** PushBack data to FXS ...
  PushToFXS( (TObject*) fCalibPulser, "TPC", "Pulser" ) ;
  
  return 0;
} // Int_t AliHLTTPCCalibPulserComponent::ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData ) {
