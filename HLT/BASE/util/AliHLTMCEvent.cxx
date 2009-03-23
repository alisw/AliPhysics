//-*- Mode: C++ -*-
// $Id: AliHLTMCEvent.cxx  $
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

/** @file   AliHLTMCEvent.cxx
    @author Jochen Thaeder
    @date   
    @brief  Container class for an AliMCEvent
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif



#include "AliHLTMCEvent.h"
#include "AliStack.h"

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMCEvent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTMCEvent::AliHLTMCEvent()
  : 
  fNParticles(0),
  fCurrentParticleIndex(-1),
  fCurrentParticle(NULL),
  fStack( new TClonesArray("TParticle", 1000) ) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTMCEvent::AliHLTMCEvent( Int_t iNumberTracks)
  :
  fNParticles(0),
  fCurrentParticleIndex(-1),
  fCurrentParticle(NULL),
  fStack( new TClonesArray("TParticle", iNumberTracks) ) {
  // see header file for class documentation

}

// #################################################################################
AliHLTMCEvent::AliHLTMCEvent( AliMCEvent *pMCEvent )
  :
  fNParticles(0),
  fCurrentParticleIndex(-1),
  fCurrentParticle(NULL),
  fStack(NULL) {
  // see header file for class documentation

  FillMCEvent( pMCEvent );
}

// #################################################################################
AliHLTMCEvent::~AliHLTMCEvent() {
  // see header file for class documentation

  if ( fStack ) {
    fStack->Delete();
    delete fStack;
  }
  fStack = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Getter
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
TParticle* AliHLTMCEvent::Particle( Int_t iParticle ) {
  // see header file for class documentation
 
  if ( iParticle >= fNParticles || !fStack  )
    return NULL;
 
  return (TParticle*) (*fStack)[iParticle];
}

// #################################################################################
TParticle* AliHLTMCEvent::NextParticle() {
  // see header file for class documentation

  fCurrentParticleIndex++;
  return Particle( fCurrentParticleIndex );
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Setter
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
void AliHLTMCEvent::AddParticle( const TParticle* particle ) {
  // see header file for class documentation
  
  new( (*fStack) [fNParticles] ) TParticle( *particle );
  fNParticles++;

  return;
}

// #################################################################################
void AliHLTMCEvent::FillMCEvent( AliMCEvent *pMCEvent ) {
  // see header file for class documentation
  
  AliStack *stack = pMCEvent->Stack();

  // -- Create local stack
  if ( stack && stack->GetNtrack() > 0 )
    fStack = new TClonesArray("TParticle", stack->GetNtrack() );
  else
    return;
  
  // -- Loop over off-line stack and fill local stack
  for (Int_t iterStack = 0; iterStack < stack->GetNtrack(); iterStack++) {

    TParticle *particle = stack->Particle(iterStack);
    if ( !particle) {
      printf( "Error reading particle %i out of %i \n", iterStack,stack->GetNtrack() );
      HLTError( "Error reading particle %i out of %i \n", iterStack,stack->GetNtrack() );
      continue;
    }

    // ----------------
    // -- Apply cuts         --> Do be done better XXX
    // ----------------

    // -- primary
    if ( !(stack->IsPhysicalPrimary(iterStack)) )
      continue;
      
    // -- final state
    if ( particle->GetNDaughters() != 0 )
      continue;

    // -- particle in DB
    TParticlePDG * particlePDG = particle->GetPDG();
    if ( ! particlePDG ) {
      particlePDG = TDatabasePDG::Instance()->GetParticle( particle->GetPdgCode() );

      if ( ! particlePDG ) {
	HLTError("Particle %i not in PDG database", particle->GetPdgCode() );
	continue;
      }
    }
  
    // -- only charged particles
    //  if ( !(particle->GetPDG()->Charge()) )
    //continue;

    // -- Add particle after cuts
    AddParticle ( particle );
  }

  Compress();

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Helper
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
void AliHLTMCEvent::Compress() {
  // see header file for class documentation

  fStack->Compress();
}

// #################################################################################
void AliHLTMCEvent::Reset() {
  // see header file for class documentation

  fCurrentParticleIndex = -1; 
  fCurrentParticle      = NULL;
}
