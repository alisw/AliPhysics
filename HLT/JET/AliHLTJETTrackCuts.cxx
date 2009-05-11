//-*- Mode: C++ -*-
// $Id: AliHLTJETTrackCuts.cxx  $
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

/** @file   AliHLTJETTrackCuts.h
    @author Jochen Thaeder
    @date   
    @brief  Cuts for jet input tracks
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTJETTrackCuts.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETTrackCuts)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETTrackCuts::AliHLTJETTrackCuts(const Char_t* name, const Char_t* title )
  : 
  AliAnalysisCuts(name, title),
  fPtMin(0.) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETTrackCuts::~AliHLTJETTrackCuts() {
  // see header file for class documentation

}

/*
 * ---------------------------------------------------------------------------------
 *                                   Selection
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Bool_t AliHLTJETTrackCuts::IsSelected( TObject *obj ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  if ( ! strcmp(obj->ClassName(),"TParticle") )
    bResult = IsSelected( dynamic_cast<TParticle*> (obj));
  else if ( ! strcmp(obj->ClassName(),"AliESDtrack") )
    bResult = IsSelected( dynamic_cast<AliESDtrack*> (obj));
  else {
    HLTError("Unknown object type %s", obj->ClassName() );
    bResult = kFALSE;
  }
  
   HLTError("Unknown object dd type %s", obj->ClassName() );

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETTrackCuts::IsSelected( TParticle *particle ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  Int_t   status  = particle->GetStatusCode();
  Int_t   pdg     = TMath::Abs( particle->GetPdgCode() );
  
  // -- Skip non-final state particles (status != 1), neutrinos (12,14,16)
  if ( (status != 1) || (pdg == 12 || pdg == 14 || pdg == 16) )
    bResult = kFALSE;
  else
    HLTInfo("Is selected !");

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETTrackCuts::IsSelected( AliESDtrack */*esdTrack*/ ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  return bResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Setter
 * ---------------------------------------------------------------------------------
 */
