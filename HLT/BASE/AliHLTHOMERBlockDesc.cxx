/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
 *          for The ALICE HLT Project.                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTHOMERBlockDesc.cxx
    @author Jochen Thaeder
    @date   
    @brief  Container for HOMER Blocks
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTMessage.h"

#include "TMath.h"
#include "TClass.h"


ClassImp(AliHLTHOMERBlockDesc)

//##################################################################################
AliHLTHOMERBlockDesc::AliHLTHOMERBlockDesc() :
  fData(NULL),
  fSize(0),
  fIsTObject(kFALSE),
  fIsRawData(kFALSE),
  fMessage(NULL),
  fClassName(),
  fDetector(),
  fSubDetector(),
  fSubSubDetector(),
  fSpecification(0),
  fDataType(),
  fHasSubDetectorRange(kFALSE),
  fHasSubSubDetectorRange(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

//##################################################################################
AliHLTHOMERBlockDesc::AliHLTHOMERBlockDesc( void * data, ULong_t size, TString origin, 
					    TString dataType, ULong_t specification ) :
  fData(data),
  fSize(size),
  fIsTObject(kFALSE),
  fIsRawData(kFALSE),
  fMessage(NULL),
  fClassName(),
  fDetector(origin),
  fSubDetector(),
  fSubSubDetector(), 
  fSpecification(specification),
  fDataType(dataType),
  fHasSubDetectorRange(kFALSE),
  fHasSubSubDetectorRange(kFALSE) {
  // see header file for class documentation
  
  // -- Set block parameters
  SetBlockParameters();

} 

//##################################################################################
AliHLTHOMERBlockDesc::~AliHLTHOMERBlockDesc() {
  // see header file for class documentation

  if ( fMessage != NULL )
    delete fMessage;
  fMessage = NULL;
}

//##################################################################################
void AliHLTHOMERBlockDesc::SetBlock( void * data, ULong_t size, TString origin, 
				     TString dataType, ULong_t specification ) {
  // see header file for class documentation

  fData = data;
  fSize = size;
  fDetector = origin; 
  fDataType = dataType;
  fSpecification = specification; 

  // -- Set block parameters
  SetBlockParameters();

  return;
}

//##################################################################################
void AliHLTHOMERBlockDesc::SetBlockParameters() {
  // see header file for class documentation

  Int_t iError = 0;

  // ---- SET CLASS NAME, DATA CONTENTS ----

  // -- Check if block contains raw data
  if ( ! CheckIfRawData() ) {

    // -- Check if block contains TObject  
    if ( ! CheckIfTObject() ) {

      // -- Contains arbitrary data type
      
      // **** TPC ****
      if ( ! fDetector.CompareTo("TPC") ) {
	
	if ( ! fDataType.CompareTo("CLUSTERS") )
	  fClassName = "AliHLTTPCSpacePoints";
	else 
	  iError = 1;
      }

      // **** TRD ****
      
      else if ( ! fDetector.CompareTo("TRD") ) {
	iError = 1;
      }

      // **** PHOS ****
      
      else if ( ! fDetector.CompareTo("PHOS") ) {
	iError = 1;
      }

      // **** MUON ****
      
      else if ( ! fDetector.CompareTo("MUON") ) {
	iError = 1;
      }

      // **** OTHER ****

      else {
	iError = 1;
      }
    }
  } // if ( ! CheckIfRawData() ) {
  
  // -- Check if Data Type has been defined
  if ( iError ) {
    //AliError( Form("The data type %s for the detector %s has not been defined yet.", 
    //		   fDataType.Data(), fDetector.Data()) );
  }
  
  // ---- SET SPECIFICATIONS ----
  // -- Convert Specification to "SubDetector/SubSubDetector"

  // **** TPC **** ( has special treatment )

  if ( ! fDetector.CompareTo("TPC") ) {

    Int_t minPatch  = (fSpecification & 0x000000FF);
    Int_t maxPatch  = (fSpecification & 0x0000FF00) >> 8;
    Int_t minSector = (fSpecification & 0x00FF0000) >> 16 ;
    Int_t maxSector = (fSpecification & 0xFF000000) >> 24;
    
    fSubDetector += minSector;
    fSubSubDetector += minPatch;
  
    // -- check for ranges
    if ( minSector != maxSector )
      fHasSubDetectorRange = kTRUE;

    if ( minPatch != maxPatch )
      fHasSubSubDetectorRange = kTRUE;
  }
  
  // **** OTHER DETECTORS ****
  
  else {
    
    if ( fSpecification ) {
      // find the max bin which is set to 1
      fSubDetector += TMath::FloorNint( TMath::Log2(fSpecification) );

      // -- check for ranges
      if ( TMath::Log2(fSpecification) != ( (Double_t) TMath::FloorNint( TMath::Log2(fSpecification) ) ) )
	fHasSubDetectorRange = kTRUE;

    }
  }

  return;
}
  
//##################################################################################
Bool_t AliHLTHOMERBlockDesc::CheckIfTObject() {
  // see header file for class documentation

  // -- Check Length - First 32 bit word of payload contains the length
  UInt_t len = *( (UInt_t*) fData  );
  
  if ( len != ( fSize - sizeof(UInt_t) ) ) 
    return fIsTObject;
    
  // -- set AliHLTMessage
  if ( fMessage ) 
    delete fMessage;
  fMessage = NULL;    

  fMessage = new AliHLTMessage( fData, fSize );
  
  // -- Check if TMessage payload is TObject
  if ( fMessage->What() == kMESS_OBJECT ) {
    fClassName = fMessage->GetClass()->GetName();
    fIsTObject = kTRUE;
  }

  return fIsTObject;
}

//##################################################################################
Bool_t AliHLTHOMERBlockDesc::CheckIfRawData() {
  // see header file for class documentation

  if ( ! fDataType.CompareTo("DDL_RAW") )
    fIsRawData = kTRUE;

  return fIsRawData;
}
 
//##################################################################################
TObject* AliHLTHOMERBlockDesc::GetTObject() {
  // see header file for class documentation
  
 if ( fMessage )
   return (TObject*) fMessage->ReadObject( fMessage->GetClass() );
 else
   return NULL;
}
