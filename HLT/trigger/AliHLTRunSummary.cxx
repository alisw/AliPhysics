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

/** @file   AliHLTRunSummary.cxx
    @author Jochen Thaeder
    @date   
    @brief  Summary class for a run, merges all detectors
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTRunSummary.h"

ClassImp(AliHLTRunSummary)
    
  AliHLTRunSummary::AliHLTRunSummary() :
    fNEvents(0),
    fNEventsRejected(0),
    fRunNumber(0),
    fRunType(0),
    fTriggerClasses(),
    fDetectorArray(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  memset( fTriggerClasses, 0, (gkNCTPTriggerClasses * sizeof(ULong_t)) );
  
  fDetectorArray = new TObjArray;
  fDetectorArray->SetOwner();
}

AliHLTRunSummary::~AliHLTRunSummary() {
  // see header file for class documentation

  if ( fDetectorArray )
    delete fDetectorArray;
  fDetectorArray = NULL;

}

Bool_t AliHLTRunSummary::AddTriggerClass( Int_t ndx ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  if ( ndx < gkNCTPTriggerClasses ) fTriggerClasses[ndx]++;
  else  bResult = kFALSE;

  return bResult;
}

Bool_t AliHLTRunSummary::AddDetector( TObject * obj ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;
  
  if ( fDetectorArray ) fDetectorArray->Add( obj );
  else bResult = kFALSE;
  
  return bResult;
}
