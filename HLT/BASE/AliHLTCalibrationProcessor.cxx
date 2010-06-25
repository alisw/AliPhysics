/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  Sebastian Bablok                                      *
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

/** 
 * @file   AliHLTCalibrationProcessor.cxx
 * @author Jochen Thaeder, Sebastian Bablok
 * @date 
 * @brief  Base class of HLT calibration components.  */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTCalibrationProcessor.h"
#include "AliHLTMemoryFile.h"
#include "AliHLTReadoutList.h"

#include <cstdlib>
#include <cerrno>
#include <string.h>
#include <TObjString.h>
#include <TFile.h>

ClassImp(AliHLTCalibrationProcessor);


const AliHLTUInt32_t AliHLTCalibrationProcessor::fgkFXSProtocolHeaderSize = 204;
const AliHLTUInt32_t AliHLTCalibrationProcessor::fgkFXSProtocolHeaderVersion = 1;

/*
 * ################ Constructor / Destructor ####################
 */

AliHLTCalibrationProcessor::AliHLTCalibrationProcessor() : 
  fEventModulo(0),
  fUseCorruptEvents(kFALSE),
  fEventCounter(0),
  fDDLNumber(),
  fDummy(0) {

  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCalibrationProcessor::~AliHLTCalibrationProcessor() {
  // see header file for class documentation
}

/*
 * ######################## InitCalibration #####################
 */

Int_t AliHLTCalibrationProcessor::DoInit( int argc, const char** argv ) {
  // see header file for class documentation

  Int_t iResult = 0;
  TString argument = "";
  TString parameter = "";
  Int_t bMissingParam=0;
  
  for ( Int_t ii=0; ii<argc && iResult>=0; ii++ ) {
    argument = argv[ii];

    if ( argument.IsNull() ) continue;

    // -eventmodulo
    if ( argument.CompareTo("-eventmodulo") == 0 ) {
      if ( ( bMissingParam=( ++ii >= argc ) ) ) break;
      parameter = argv[ii];  
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	fEventModulo = parameter.Atoi();
	HLTInfo("EventModulo is set to %d.", fEventModulo);
      } 
      else {
	HLTError("Cannot convert EventModulo specifier '%s'.", argv[ii]);
	iResult = -EINVAL;
	break;
      }
    } 
    // -usecorruptevents
    else if ( argument.CompareTo( "-usecorruptevents" ) == 0 ) {
      HLTInfo( "Passing through of corrupted events enabled." );
      fUseCorruptEvents = kTRUE;
    } 
    // not known by calibration processor
    else {
      if ( ( iResult = ScanArgument( argc-ii, &argv[ii] ) ) == -EINVAL ) {
	HLTError( "unknown argument %s", argument.Data() );
	break;
      } 
      else if ( iResult == -EPROTO ) {
	bMissingParam = 1;
	break;
      } 
      else if ( iResult >= 0 ) {
	ii += iResult;
	iResult = 0;
      }
    }
  }
  
  if ( bMissingParam ) {
    HLTError( "missing parameter for argument %s", argument.Data() );
    iResult = -EPROTO;
  }

  if ( iResult >= 0 ) {
    iResult = InitCalibration();
  }

  // Reset the DDLNumberList
  memset( fDDLNumber, 0, gkAliHLTFXSHeaderfDDLNumberSize);

  return iResult;
}

Int_t AliHLTCalibrationProcessor::InitCalibration() {
  // see header file for class documentation
  
  // fDummy is just touched here to avoid coding convention violation RC11. 
  // The function can not be declared  const since it is just the default implementation, 
  // overloaded virtual function might not be const.
  fDummy = 0;

  return 0;
}


Int_t AliHLTCalibrationProcessor::ScanArgument(int argc, const char** argv) {
  // see header file for class documentation
  
  // there are no other arguments than the standard ones
  if ( argc == 0 && argv == NULL) {
    // this is just to get rid of the warning "unused parameter"
  }
  
  // fDummy is just touched here to avoid coding convention violation RC11. 
  // The function can not be declared  const since it is just the default implementation, 
  // overloaded virtual function might not be const.
  fDummy = 0;

  return -EINVAL;
}

/*
 * ######################## DeinitCalibration #####################
 */

Int_t AliHLTCalibrationProcessor::DoDeinit() {
  // see header file for class documentation
 
  int iResult = DeinitCalibration();
 
  return iResult;
}

Int_t AliHLTCalibrationProcessor::DeinitCalibration() {
  // see header file for class documentation

  // fDummy is just touched here to avoid coding convention violation RC11. 
  // The function can not be declared  const since it is just the default implementation, 
  // overloaded virtual function might not be const.
  fDummy = 0;

  return 0;
}

/*
 * ######################## DoEvent #####################
 */

Int_t AliHLTCalibrationProcessor::DoEvent( const AliHLTComponentEventData& evtData,
	       const AliHLTComponentBlockData* blocks, 
	       AliHLTComponentTriggerData& trigData,
	       AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size,
	       vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation

  Int_t iResult = 0;

  const AliHLTComponentBlockData* blkSOR = NULL;
  blkSOR = GetFirstInputBlock( kAliHLTDataTypeSOR );
  const AliHLTComponentBlockData* blkEOR = NULL;
  blkEOR = GetFirstInputBlock( kAliHLTDataTypeEOR );

  // ** if event Type is not SOR or EOR -> fill DDLNumber list and process data
  if ( ! blkEOR  && !blkSOR ) {
    // ** ProcessData
    iResult = ProcessCalibration( evtData, blocks, trigData, outputPtr, size, outputBlocks );
    fEventCounter++; 
  }  

  // - if event Type is EOR or event modulo is set -> ship data to FXS
  // - use event modulo ship data also during the run 
  if ( ( fEventModulo > 0 && fEventCounter == fEventModulo ) || blkEOR ){
    HLTDebug ( "Ship Data to FXS!" );
    iResult = ShipDataToFXS( evtData, blocks, trigData, outputPtr, size, outputBlocks );
    fEventCounter = 0;
  }

  return iResult;
}

/*
 * ######################## ProcessCalibration #####################
 */

Int_t AliHLTCalibrationProcessor::ProcessCalibration( const AliHLTComponentEventData& evtData,
						      const AliHLTComponentBlockData* /*blocks*/, 
						      AliHLTComponentTriggerData& trigData,
						      AliHLTUInt8_t* /*outputPtr*/, 
						      AliHLTUInt32_t& /*size*/,
						      vector<AliHLTComponentBlockData>& /*outputBlocks*/ ) {
  // see header file for class documentation

  // we just forward to the high level method, all other parameters already
  // have been stored internally
  return ProcessCalibration( evtData, trigData );
}

Int_t AliHLTCalibrationProcessor::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/, 
						      AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  HLTFatal( "no processing method implemented" );
  return -ENOSYS;
}

/*
 * ######################## ShipDataToFXS #####################
 */

Int_t AliHLTCalibrationProcessor::ShipDataToFXS( const AliHLTComponentEventData& evtData,
						 const AliHLTComponentBlockData* /*blocks*/, 
						 AliHLTComponentTriggerData& trigData,
						 AliHLTUInt8_t* /*outputPtr*/, 
						 AliHLTUInt32_t& /*size*/,
						 vector<AliHLTComponentBlockData>& /*outputBlocks*/ ) {
  // see header file for class documentation

  // we just forward to the high level method, all other parameters already
  // have been stored internally
  return ShipDataToFXS( evtData, trigData );
}


Int_t AliHLTCalibrationProcessor::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, 
						 AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation

  HLTFatal( "no processing method implemented" );
  return -ENOSYS;
}

/*
 * ######################## CreateFXSHeader #####################
 */

Int_t AliHLTCalibrationProcessor::CreateFXSHeader( AliHLTFXSHeader &pHeader, const char* pDetector, const char* pFileID, const AliHLTReadoutList* pDDLList ) {
  // see header file for class documentation

  Int_t iResult = 0;

  // ** Fill header version
  pHeader.fHeaderVersion = AliHLTCalibrationProcessor::fgkFXSProtocolHeaderVersion;

  // ** Fill run number
  pHeader.fRunNumber = GetRunNo(); 
  
  // ** Fill origin
  HLTDebug( "FXS Header Detector  size max %i - actual %i .",gkAliHLTFXSHeaderfOriginSize, strlen( pDetector ) );
  strncpy ( pHeader.fOrigin, pDetector, gkAliHLTFXSHeaderfOriginSize ) ; 

  // To take care if fileIDs which are longer than gkAliHLTFXSHeaderfOriginSize, write one 0 is cheaper than an if.
  pHeader.fOrigin[gkAliHLTFXSHeaderfOriginSize-1] = 0;

  // ** Fill file ID
  HLTDebug( "FXS Header FileID size max %i - actual %i .",gkAliHLTFXSHeaderfFileIDSize, strlen( pFileID ) );
  strncpy ( pHeader.fFileID, pFileID, gkAliHLTFXSHeaderfFileIDSize ) ; 

  // To take care if fileIDs which are longer than gkAliHLTFXSHeaderfFileIDSize, write one 0 is cheaper than an if.
  pHeader.fFileID[gkAliHLTFXSHeaderfFileIDSize-1] = 0;

  // ** Fill DDL number

  // -- if component provides list, convert to fDDLNumber
  if ( pDDLList ) {
    // use user list
    
    AliHLTReadoutList::EDetectorId detid = pDDLList->GetFirstUsedDetector();
    Int_t wordNdx = AliHLTReadoutList::GetFirstWord(detid);
    Int_t wordCount = AliHLTReadoutList::GetWordCount(detid);
    
    if (pDDLList->GetFirstUsedDetector(detid) != AliHLTReadoutList::kNoDetector or wordNdx < 0)
    {
      HLTError("DDLIDs for minimum of TWO detectors ( %s, %s ) set, this function works only for ONE detector.",
          AliHLTReadoutList::DetectorIdToString(detid),
          AliHLTReadoutList::DetectorIdToString(pDDLList->GetFirstUsedDetector(detid))
        );
      iResult = -1;
    }
    else
    {
      // check word for word 
      for ( Int_t ndx = 0; ndx < wordCount; ndx++ ) {
	AliHLTUInt32_t word = pDDLList->Buffer()->fList[wordNdx+ndx];
	
	// set only 4 bit into one Char_t
	for ( Int_t charNdx = 0; charNdx < 8; charNdx++) {
	  fDDLNumber[(8*ndx)+charNdx] = (Char_t) word & 0x0000000F;
	  word = word >> 4;
	}
      }
    }
  } //   if ( pDDLList ) {

  // -- fill header with ascii chars
  for (Int_t ndx = 0; ndx < gkAliHLTFXSHeaderfDDLNumberSize; ndx++ ){
    Int_t numberToChar = (Int_t) fDDLNumber[ndx];
    // Get ASCII
    if ( numberToChar > 9 ) numberToChar += 55;
    else numberToChar += 48;
    
    pHeader.fDDLNumber[ndx] = (Char_t) numberToChar;
  }
  
  return iResult;
}  // Int_t AliHLTCalibrationProcessor::CreateFXSHeader( AliHLTXSHeader &pHeader, const char* pDetector, const char* pFileID, const AliHLTReadoutList* pDDLList ) {

/*
 * ######################## PushToFXS #####################
 */

Int_t AliHLTCalibrationProcessor::PushToFXS(TObject* pObject, const char* pDetector, const char* pFileID, const AliHLTReadoutList* pDDLList ) {
  // see header file for class documentation

  Int_t iResult = 0;
  
  AliHLTFXSHeader pHeader;

  CreateFXSHeader( pHeader, pDetector, pFileID, pDDLList );
  
  if ( pObject ) {
    
    AliHLTMemoryFile* pMemFile = CreateMemoryFile( kAliHLTDataTypeFXSCalib, kAliHLTVoidDataSpec );
    if ( pMemFile ) {
      
      iResult = pMemFile->WriteHeaderBuffer( (const char*) &pHeader, AliHLTCalibrationProcessor::fgkFXSProtocolHeaderSize );

      if ( iResult ) {
	HLTError( "Buffer size to small - for header!" ); 
      }
      else {
	iResult = Write( pMemFile, pObject );
	if ( !iResult ) {
	  HLTError( "Buffer size to small - for data!" ); 
	}
	else {
	  iResult = CloseMemoryFile( pMemFile );  
	}
      }
    } 
    else {
      iResult = -ENOMEM;
    }
  } 
  else {
    iResult=-EINVAL;
  }

  return iResult;

} // Int_t AliHLTCalibrationProcessor::PushToFXS(TObject* pObject, const char* detector, const char* pFileID, const AliHLTReadoutList* pDDLList ) {

Int_t AliHLTCalibrationProcessor::PushToFXS( void* pBuffer, int iSize, const char* pDetector, const char* pFileID, const AliHLTReadoutList* pDDLList ) {
  // see header file for class documentation

  Int_t iResult = 0;
  
  AliHLTFXSHeader pHeader;

  CreateFXSHeader( pHeader, pDetector, pFileID, pDDLList );
  
  iResult = PushBack( pBuffer, iSize, kAliHLTDataTypeFXSCalib, kAliHLTVoidDataSpec, (void*) (&pHeader), AliHLTCalibrationProcessor::fgkFXSProtocolHeaderSize );

  return iResult;

} // Int_t AliHLTCalibrationProcessor::PushToFXS(void* pBuffer, int iSize, const char* pDdetector, const char* pFileID, const AliHLTReadoutList* pDDLList ) {

