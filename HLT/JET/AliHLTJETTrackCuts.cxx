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
#include <TDatabasePDG.h>

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
  fChargedOnly(kFALSE),
  fPtMin(0.0),
  fEtaMin(-0.9),
  fEtaMax(0.9),
  fPhiMin(0.0),
  fPhiMax(6.3) {
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
  
  if ( obj->IsA() == TParticle::Class() )
    bResult = IsSelected( static_cast<TParticle*> (obj));
  else if ( obj->IsA() == AliESDtrack::Class() )
    bResult = IsSelected( static_cast<AliESDtrack*> (obj));
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

  // ----------------------------------
  //  Applied before in AliHLTMCEvent:
  //  - Is Physical Primary  
  //   stack->IsPhysicalPrimary(iter)
  //  - final state
  //   particle->GetNDaughters() == 0
  // ----------------------------------

  Bool_t bResult = kTRUE;

  Int_t   status  = particle->GetStatusCode();
  Int_t   pdgCode = TMath::Abs( particle->GetPdgCode() );
  
  // -- Skip non-final state particles (status != 1), neutrinos (12,14,16)
  if ( (status != 1) || (pdgCode == 12 || pdgCode == 14 || pdgCode == 16) )
    bResult = kFALSE;

  // -- Charged particles only
  if ( ! particle->GetPDG() )
    bResult = kFALSE;
  else {
    if ( fChargedOnly && !particle->GetPDG()->Charge() )
      bResult = kFALSE;
  }

  // -- cut on min Pt
  if ( particle->Pt() < fPtMin )
    bResult = kFALSE;

  // -- cut on eta acceptance
  if ( ( particle->Eta() < fEtaMin ) || ( particle->Eta() > fEtaMax ) )
    bResult = kFALSE;

  // -- cut on phi acceptance
  if ( ( particle->Phi() < fPhiMin ) || ( particle->Phi() > fPhiMax ) )
    bResult = kFALSE;

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETTrackCuts::IsSelected( AliESDtrack *esdTrack ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  // -- cut on min Pt
  if ( esdTrack->Pt() < fPtMin )
    bResult = kFALSE;

  // -- cut on eta acceptance
  if ( ( esdTrack->Eta() < fEtaMin ) || ( esdTrack->Eta() > fEtaMax ) )
    bResult = kFALSE;

  // -- cut on phi acceptance
  if ( ( esdTrack->Phi() < fPhiMin ) || ( esdTrack->Phi() > fPhiMax ) )
    bResult = kFALSE;

  return bResult;
}
