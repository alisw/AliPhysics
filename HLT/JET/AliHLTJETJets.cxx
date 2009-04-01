//-*- Mode: C++ -*-
// $Id: AliHLTJETJets.cxx  $
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

/** @file   AliHLTJETJets.h
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

#include "AliHLTJETJets.h"
#include "TLorentzVector.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETJets)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliHLTJETJets::AliHLTJETJets() :
  fNAODJets(0),
  fAODJets(new TClonesArray( "AliAODJet", 20 ) ) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

//##################################################################################
AliHLTJETJets::~AliHLTJETJets() {
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
void AliHLTJETJets::Reset() {
  // see header file for class documentation  

  fAODJets->Clear();
  fNAODJets = 0;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Getter
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliAODJet* AliHLTJETJets::GetJet( Int_t iter )  { 
  // see header file for class documentation
  
  if ( iter > fNAODJets )
    return NULL;
  else
    return reinterpret_cast<AliAODJet*>((*fAODJets)[iter]); 
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Setter
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETJets::AddJet( AliHLTJETConeJetCandidate* jet ) {
  // see header file for class documentation
  
  // -- create TLorentzVector
  TLorentzVector v;
  v.SetPtEtaPhiE( jet->GetPt(),jet->GetEta(),jet->GetPhi(), jet->GetEt() );

  // -- add AliAODJet
  new ((*fAODJets)[fNAODJets]) AliAODJet(v);
  fNAODJets++;

  return;
}
