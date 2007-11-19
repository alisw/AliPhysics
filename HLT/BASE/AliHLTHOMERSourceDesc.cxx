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

/** @file   AliHLTHOMERSourceDesc.cxx
    @author Jochen Thaeder
    @date   
    @brief  Container for HOMER Sources
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliHLTHOMERSourceDesc.h"


ClassImp(AliHLTHOMERSourceDesc)

//##################################################################################
AliHLTHOMERSourceDesc::AliHLTHOMERSourceDesc() :
  fSelected( kFALSE ),
  fHostname(),
  fPort(),
  fSourceName(),
  fClassName(),
  fDetector(),
  fSubDetector(),
  fSubSubDetector(),
  fSpecification(),
  fDataType() {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


//##################################################################################
AliHLTHOMERSourceDesc::AliHLTHOMERSourceDesc( TString hostname, Int_t port ) :
  fSelected( kFALSE ),
  fHostname( hostname ),
  fPort( port ),
  fSourceName(),
  fClassName(),
  fDetector(),
  fSubDetector(),
  fSubSubDetector(), 
  fSpecification(),
  fDataType() {
  // see header file for class documentation

}


//##################################################################################
AliHLTHOMERSourceDesc::~AliHLTHOMERSourceDesc() {
  // see header file for class documentation
}

