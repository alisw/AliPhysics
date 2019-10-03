//-*- Mode: C++ -*-
// $Id: AliHLTJETFastJetHeader.cxx  $
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

/** @file   AliHLTJETFastJetHeader.cxx
    @author Jochen Thaeder
    @date   
    @brief  Header of the FastJet finder interface
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTJETJetCuts.h"

#include "AliHLTJETReaderHeader.h"
#include "AliHLTJETFastJetHeader.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETFastJetHeader)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETFastJetHeader::AliHLTJETFastJetHeader() 
  :
  AliJetHeader("AliHLTJETFastJetHeader"),
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fGhostArea(0.05),
  fActiveAreaRepeats(1),
  fAreaType(fastjet::active_area), 
  fPtMin(5.),
  fJetDefinition(NULL),
  fGhostedAreaSpec(NULL),
  fAreaDefinition(NULL),
  fRangeDefinition(NULL),
  fReaderHeader(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
  //  SetRapRange(rapmin, rapmax);
}

// #################################################################################
AliHLTJETFastJetHeader::~AliHLTJETFastJetHeader() {
  // see header file for class documentation
  
  if ( fJetDefinition ) 
    delete fJetDefinition;
  fJetDefinition = NULL;

  if ( fGhostedAreaSpec)
    delete fGhostedAreaSpec;
  fGhostedAreaSpec = NULL;

  if ( fAreaDefinition )
    delete fAreaDefinition;
  fAreaDefinition = NULL;
   
  if ( fRangeDefinition )
    delete fRangeDefinition;
  fRangeDefinition = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Initialize
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETFastJetHeader::Initialize() {
  // see header file for class documentation
  
  HLTInfo(" -= AliHLTJETFastJetHeader =- " );

  Int_t iResult = 0;

  AliHLTJETReaderHeader*  readerHeader = dynamic_cast<AliHLTJETReaderHeader*>(fReaderHeader);

  if ( !fJetCuts ) {
    HLTError("No HLT jet cuts set!");
    return -1;
  }
  
  // -- Set min pt
  fPtMin = fJetCuts->GetMinEt();

  // -- Set Algorithm
  AliHLTJETBase::JetAlgorithmType_t GetJetAlgorithm();

  if ( readerHeader->GetJetAlgorithm() == AliHLTJETBase::kKt )
    fAlgorithm = fastjet::kt_algorithm;
  else if ( readerHeader->GetJetAlgorithm() == AliHLTJETBase::kAntiKt )
    fAlgorithm = fastjet::antikt_algorithm;
  else {
    HLTError("No algorithm found!");
    return -1;
  }

  // -- Create an object that represents your choice of jet algorithm,
  //    and the associated parameters
  if ( fJetDefinition ) 
    delete fJetDefinition;
  fJetDefinition = new fastjet::JetDefinition(fAlgorithm, readerHeader->GetConeRadius(),
					      fRecombScheme, fStrategy);

  // -- Create the object that holds info about ghosts
  if ( fGhostedAreaSpec ) 
    delete fGhostedAreaSpec;
  fGhostedAreaSpec = new fastjet::GhostedAreaSpec(readerHeader->GetFiducialEtaMax(),
						  fActiveAreaRepeats, fGhostArea);
  
  // -- Create area definition
  if ( fAreaDefinition ) 
    delete fAreaDefinition;
  // xxx  fAreaDefinition = new fastjet::AreaDefinition(fAreaType, fGhostedAreaSpec);
  fAreaDefinition = new fastjet::AreaDefinition(fAreaType, (*fGhostedAreaSpec) );
  
  // -- Set the rapididty, phi range within which to study the background 
  if ( fRangeDefinition )
    delete fRangeDefinition;
  fRangeDefinition = new fastjet::RangeDefinition( readerHeader->GetFiducialEtaMin()+readerHeader->GetConeRadius(),
						   readerHeader->GetFiducialEtaMax()-readerHeader->GetConeRadius(),
						   readerHeader->GetFiducialPhiMin(),
						   readerHeader->GetFiducialPhiMax() );

  // -- Check initialization
  if ( !fJetDefinition || !fGhostedAreaSpec || 
       !fAreaDefinition || !fRangeDefinition ) {
    HLTError("Initializing FastJet failed!");
    iResult = kFALSE;
  }

  // -- Create comment
  if ( !iResult ) {
    fComment = "Running FastJet algorithm with the following setup.";
    fComment+= "\n-Jet definition:\n  ";
    fComment+= TString(fJetDefinition->description());
    fComment+= ". \n-Area definition:\n  ";
    fComment+= TString(fAreaDefinition->description());
    fComment+= ". \n-Strategy adopted by FastJet:\n  ";
  }
  
  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                                      Helper
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
void AliHLTJETFastJetHeader::PrintParameters() const {
  // see header file for class documentation

  HLTInfo( "FastJet algorithm  parameters:" );
  
  HLTInfo( "-- Jet Definition --- " );
#if 0
  HLTInfo( "R %f ",fRparam );
  HLTInfo( "Jet Algorithm %s", fAlgorithm ); 
  HLTInfo( "Strategy " << fStrategy );  
  HLTInfo( "Recombination Scheme " << fRecombScheme ); 
  
  HLTInfo( "-- Ghosted Area Spec parameters --- " );
  HLTInfo( "Ghost Eta Max " << fGhostEtaMax );
  HLTInfo( "Ghost Area " << fGhostArea );
  HLTInfo( "Active Area Repeats " << fActiveAreaRepeats );
  
  HLTInfo( "-- Area Definition parameters --- " );
  HLTInfo( "Area Type " << fAreaType ); 
  
  HLTInfo( "-- Cluster Sequence Area parameters --- " );
  HLTInfo( "pt min " << fPtMin ); 
  
  HLTInfo( "-- Range Definition parameters --- " );
  HLTInfo( " bkg rapidity range from  " << fRapMin << " to " << fRapMax );
  HLTInfo( " bkg phi range from " << fPhiMin << " to " << fPhiMax );
#endif
}
