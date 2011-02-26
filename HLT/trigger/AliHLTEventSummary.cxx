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

/** @file   AliHLTEventSummary.cxx
    @author Jochen Thaeder
    @date   
    @brief  Summary class for run statistics, merges all detectors
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTEventSummary.h"

ClassImp(AliHLTEventSummary)
    
  AliHLTEventSummary::AliHLTEventSummary() : 
    fRejected(kTRUE),
    fRunNumber(0),
    fRunType(0),
    fTriggerClass(0),
    fDetectorArray(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  fDetectorArray = new TObjArray;
  fDetectorArray->SetOwner();
}

AliHLTEventSummary::~AliHLTEventSummary() {
  // see header file for class documentation

  if ( fDetectorArray )
    delete fDetectorArray;
  fDetectorArray = NULL;
}

Bool_t AliHLTEventSummary::AddDetector( TObject * obj ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;
  
  if ( fDetectorArray ) fDetectorArray->Add( obj );
  else bResult = kFALSE;
  
  return bResult;
}
