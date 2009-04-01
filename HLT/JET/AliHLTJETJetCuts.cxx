//-*- Mode: C++ -*-
// $Id: AliHLTJETJetCuts.cxx  $
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

/** @file   AliHLTJETJetCuts.h
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

#include "AliHLTJETJetCuts.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETJetCuts)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETJetCuts::AliHLTJETJetCuts(const Char_t* name, const Char_t* title )
  : 
  AliAnalysisCuts(name, title),
  fEtMin(0.) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETJetCuts::~AliHLTJETJetCuts() {
  // see header file for class documentation

}

/*
 * ---------------------------------------------------------------------------------
 *                                   Selection
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Bool_t AliHLTJETJetCuts::IsSelected( TObject *obj ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  if ( ! strcmp(obj->ClassName(),"AliHLTJETConeJetCandidate") )
    bResult = IsSelected( dynamic_cast<AliHLTJETConeJetCandidate*> (obj));
  else if ( ! strcmp(obj->ClassName(),"AliAODJet") )
    bResult = IsSelected( dynamic_cast<AliAODJet*> (obj));
  else {
    HLTError("Unknown object type %s", obj->ClassName() );
    bResult = kFALSE;
  }
  
   HLTError("Unknown object dd type %s", obj->ClassName() );

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETJetCuts::IsSelected( AliHLTJETConeJetCandidate* jet ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  // -- cut on min Pt
  if ( jet->GetEt() < fEtMin )
    bResult = kFALSE;

  return bResult;
}

// #################################################################################
Bool_t AliHLTJETJetCuts::IsSelected( AliAODJet* jet ) {
  // see header file for class documentation

  Bool_t bResult = kTRUE;

  // -- cut on min Pt
  if ( jet->Pt() < fEtMin )
    bResult = kFALSE;

  return bResult;
}


#if 0
/*
Int_t AliHLTJetFinder::CleanJetCandidates(){
  // see header file for class documentation

  Int_t iResult = 0;
  
  fJets->SortJets();
  
  Int_t nIntitalJetCandidates = fJets->GetNJets();

  for ( Int_t iter = 0; iter < nIntitalJetCandidates; iter++ ) {

    if ( ! fJets->IsJet( iter ) )
	continue;

    AliHLTJetFinderJetCandidate* jet = fJets->GetJet( iter );

    for ( Int_t compareIter = 0; compareIter < iter ; compareIter++ ) {

      if ( ! fJets->IsJet( compareIter ) )
	continue;

      if ( ! fJets->IsJet( iter ) )
	break;

      
      AliHLTJetFinderJetCandidate* compareJet = fJets->GetJet( compareIter );

      Double_t distance2 = AliHLTJetDefinitions::GetDistance2( jet->GetEta(), jet->GetPhi(), 
							       compareJet->GetEta(), compareJet->GetPhi() );

      // -- check if Jet is close to another one
      if ( distance2 > ( fDistanceCutJet*fDistanceCutJet) ) 
	continue;

      // -- one has to go
      if ( jet->GetPt() >= compareJet->GetPt() )
	fJets->RemoveJet( compareIter );
      
      else {
	fJets->RemoveJet( iter );
	break;
      }

    } // for ( Int_t compareIter = 0; compareIter < iter ; compareIter++ ) {
    
  } // for ( Int_t iter = 0; iter < fNJetCandidates; iter++ ) {
  
  fJets->CompressJets();

  return iResult;
}

*/
#endif

