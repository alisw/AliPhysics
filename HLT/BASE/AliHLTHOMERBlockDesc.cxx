// $Id$
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

#include "AliHLTHOMERBlockDesc.h"

#include "AliHLTMessage.h"

#include "TMath.h"
#include "TClass.h"

ClassImp(AliHLTHOMERBlockDesc)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor 
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
AliHLTHOMERBlockDesc::AliHLTHOMERBlockDesc() :
  fData(NULL),
  fSize(0),
  fBlockName(),
  fIsTObject(kFALSE),
  fIsRawData(kFALSE),
  fMessage(NULL),
  fTObject(NULL),
  fClassName(),
  fDataType(),
  fDetector(),
  fSpecification(0),
  fSubDetector(0),
  fSubSubDetector(0),
  fHasSubDetectorRange(kFALSE),
  fHasSubSubDetectorRange(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

//##################################################################################
AliHLTHOMERBlockDesc::~AliHLTHOMERBlockDesc() {
  // see header file for class documentation

  if ( fMessage != NULL )
    delete fMessage;
  fMessage = NULL;

  if ( fData )
    delete [] fData;
  fData = NULL;

  if ( fTObject ) 
    delete fTObject;
  fTObject = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Data Handling - Setter - public
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
void AliHLTHOMERBlockDesc::SetBlock( void * data, ULong_t size, TString origin, 
				     TString dataType, ULong_t specification ) {
  // see header file for class documentation

  fData = new Char_t[size];
  memcpy( fData, data, size);
  
  fSize = size;
  fDetector = origin; 
  fDataType = dataType;
  fSpecification = specification; 

  fBlockName.Form("%s_%s_0x%08lX", fDetector.Data(), fDataType.Data(), fSpecification ); 

  // -- Set block parameters
  SetBlockParameters();

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Data Handling - private
 * --------------------------------------------------------------------------------- 
 */

//##################################################################################
void AliHLTHOMERBlockDesc::SetBlockParameters() {
  // see header file for class documentation

  //Int_t iResult = 0;

  // ---- SET SPECIFICATIONS ----
  // ----------------------------

  // **** TPC **** ( has special treatment )
  if ( ! fDetector.CompareTo("TPC") ) {
    
    Int_t minPatch  = (fSpecification & 0x000000FF);
    Int_t maxPatch  = (fSpecification & 0x0000FF00) >> 8;
    Int_t minSector = (fSpecification & 0x00FF0000) >> 16 ;
    Int_t maxSector = (fSpecification & 0xFF000000) >> 24;
    
    fSubDetector = minSector;
    fSubSubDetector = minPatch;
    
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
      fSubDetector = TMath::FloorNint( TMath::Log2(fSpecification) );
      
      // -- check for ranges
      if ( TMath::Log2(fSpecification) != 
	   static_cast<Double_t>(TMath::FloorNint(TMath::Log2(fSpecification))) )
	fHasSubDetectorRange = kTRUE;
    }
  }

  // ---- SET CLASS NAME, DATA CONTENTS ----
  // ---------------------------------------

  // -- Check if block contains raw data
  if ( CheckIfRawData() )
    return;

  // -- Check if block contains TObject  
  if ( CheckIfTObject() )
    return;

  // -- Contains arbitrary data type
      
  // **** TPC ****
  if ( ! fDetector.CompareTo("TPC") ) {
    
    if ( ! fDataType.CompareTo("CLUSTERS") )
      fClassName = "AliHLTTPCSpacePoints";
    //else 
    //  iResult = -1;
  }
  /*
  // **** TRD ****
  else if ( ! fDetector.CompareTo("TRD") ) {
    iResult = -1;
  }
  
  // **** PHOS ****
  else if ( ! fDetector.CompareTo("PHOS") ) {
    iResult = -1;
  }
  
  // **** MUON ****
  else if ( ! fDetector.CompareTo("MUON") ) {
    iResult = -1;
  }
  
  // **** OTHER ****
  else {
    iResult = -1;
  }
  */
  
  // -- Check if classname has been defined
  //  if ( iResult < 0 ) {
    //   AliWarning( Form("The classname for data type %s for the detector %s has not been defined yet.", 
    //	     fDataType.Data(), fDetector.Data()) );
  // }
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
  if ( fMessage->What() == kMESS_OBJECT and fMessage->GetClass() != NULL and fMessage->GetClass() != (TClass*)-1)
  {
    fClassName = fMessage->GetClass()->GetName();
    fIsTObject = kTRUE;
    
    fTObject = fMessage->ReadObject( fMessage->GetClass() );
  }
  
  fMessage->Reset();

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
void AliHLTHOMERBlockDesc::ComponentDataType(AliHLTComponentDataType& dataType)
{
  dataType.fStructSize = sizeof(AliHLTComponentDataType);
  memcpy(dataType.fID, fDataType.Data(), kAliHLTComponentDataTypefIDsize);
  memcpy(dataType.fOrigin, fDetector.Data(), kAliHLTComponentDataTypefOriginSize);
}

