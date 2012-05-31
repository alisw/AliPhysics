//-*- Mode: C++ -*-
// $Id: AliHLTJets.cxx  $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTJets.h
    @author Jochen Thaeder
    @date   
    @brief  Container holding produced Jets
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "TLorentzVector.h"

#include "AliHLTJets.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJets)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTJets::AliHLTJets() :
  fComment(""),
  fCurrentJetIndex(-1),
  fNAODJets(0),
  fAODJets(new TClonesArray( "AliAODJet", 50 )) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

//##################################################################################
AliHLTJets::~AliHLTJets() {
  // see header file for class documentation

  if ( fAODJets ){
    fAODJets->Clear();
    delete fAODJets;
  }
  fAODJets = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                                   Initialize / Reset
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJets::Reset() {
  // see header file for class documentation  

  fAODJets->Clear();
  fNAODJets = 0;

  return;
}

// #################################################################################
void AliHLTJets::ResetEvent() {
  // see header file for class documentation

  fCurrentJetIndex = -1; 
  
  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Getter
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliAODJet* AliHLTJets::GetJet( Int_t iter ) const { 
  // see header file for class documentation
  
  if ( iter > fNAODJets )
    return NULL;
  else
    return reinterpret_cast<AliAODJet*>((*fAODJets)[iter]); 
}

//##################################################################################
AliAODJet* AliHLTJets::NextJet() { 
  // see header file for class documentation
  
  fCurrentJetIndex++;
  return GetJet( fCurrentJetIndex );
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Setter
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJets::AddJet( Float_t eta, Float_t phi, Float_t pt, Float_t et ) {
  // see header file for class documentation
  
  if ( pt == 0. ) {
    HLTError("Jet Pt=0, eta=%f, phi=%f", eta, phi);
    return;
  }

  // -- create TLorentzVector
  TLorentzVector v;
  v.SetPtEtaPhiE( pt, eta, phi, et );

  // -- add AliAODJet
  new ((*fAODJets)[fNAODJets]) AliAODJet(v);
  fNAODJets++;

  return;
}

//##################################################################################
void AliHLTJets::AddJet( AliAODJet* jet ) {
  // see header file for class documentation

  if ( jet->Pt() == 0. ) {
    HLTError("Jet Pt=0, eta=%f, phi=%f", jet->Eta(), jet->Phi());
    return;
  }

  // -- create TLorentzVector
  TLorentzVector v;
  v.SetPtEtaPhiE( jet->Pt(), jet->Eta(), jet->Phi(), jet->E() );

  // -- add AliAODJet
  new ((*fAODJets)[fNAODJets]) AliAODJet(v);
  fNAODJets++;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Helper
 * ---------------------------------------------------------------------------------
 */

/** Sort Jets with decreasing Et */
void AliHLTJets::Sort() {
  // see header file for class documentation

  fAODJets->Sort();
  ResetEvent();

  return;
}
