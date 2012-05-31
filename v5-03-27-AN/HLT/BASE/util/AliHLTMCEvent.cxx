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

#include "TList.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "AliGenCocktailEventHeader.h"

#include "AliHLTMCEvent.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMCEvent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTMCEvent::AliHLTMCEvent() : 
  fCurrentParticleIndex(-1),
  fNParticles(0),
  fStack( NULL ),
  fCurrentGenJetIndex(-1),
  fNGenJets(0),
  fGenJets(NULL),
  fHasParticleCuts(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTMCEvent::AliHLTMCEvent( Bool_t applyParticleCuts ) : 
  fCurrentParticleIndex(-1),
  fNParticles(0),
  fStack( NULL ),
  fCurrentGenJetIndex(-1),
  fNGenJets(0),
  fGenJets(NULL),
  fHasParticleCuts(applyParticleCuts) {
  // see header file for class documentation
}

// #################################################################################
AliHLTMCEvent::~AliHLTMCEvent() {
  // see header file for class documentation

  if ( fStack ) {
    fStack->Delete();
    delete fStack;
  }
  fStack = NULL;

  if ( fGenJets ) {
    fGenJets->Delete();
    delete fGenJets;
  }
  fGenJets = NULL;
}

/*
 * ---------------------------------------------------------------------------------
 *                               Setter - public
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTMCEvent::FillMCEvent( AliStack *stack, AliHeader *header ) {
  // see header file for class documentation

  Int_t iResult = 0;

  if ( stack ) {
    if ( (iResult = FillMCTracks(stack)) )
      HLTError("Error filling particles" );
  }
  else {
    HLTError("Error reading stack" );
    iResult = -2;
  }

  if ( header ) {
    AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*> (header->GenEventHeader());
    if ( pythiaGenHeader ) {
      if ( (iResult = FillMCJets(pythiaGenHeader)) )
	HLTError("Error filling jets" );
    }
    else {
      HLTError("Error reading PythiaHeader" );
      iResult = -2;
    }
  }
  else {
    HLTError("Error reading header" );
    iResult = -2;
  }
    
  Compress();

  return iResult;;
}

// #################################################################################
Int_t AliHLTMCEvent::FillMCEvent( AliMCEvent *pMCEvent ) {
  // see header file for class documentation
  
  Int_t iResult = 0;
  
  if ( (iResult = FillMCTracks(pMCEvent->Stack())) )
    HLTError("Error filling particles" );

  if ( (iResult = FillMCJets(GetPythiaEventHeader(pMCEvent))) )
    HLTError("Error filling jets" );
  
  Compress();

  return iResult;
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

// #################################################################################
AliAODJet* AliHLTMCEvent::GenJet( Int_t iJet ) const {
  // see header file for class documentation
 
  if ( iJet >= fNGenJets || !fGenJets  )
    return NULL;
 
  return reinterpret_cast<AliAODJet*> ((*fGenJets)[iJet]);
}

// #################################################################################
AliAODJet* AliHLTMCEvent:: NextGenJet() {
  // see header file for class documentation

  fCurrentGenJetIndex++;
  return GenJet( fCurrentGenJetIndex );
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Helper
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
void AliHLTMCEvent::Compress() {
  // see header file for class documentation

  if (fStack)
    fStack->Compress();
  if (fGenJets)
    fGenJets->Compress();
}

// #################################################################################
void AliHLTMCEvent::Reset() {
  // see header file for class documentation

  fCurrentParticleIndex = -1; 
  fCurrentGenJetIndex = -1; 
}

/*
 * ---------------------------------------------------------------------------------
 *                               Setter - private
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTMCEvent::FillMCTracks( AliStack* stack ) {
  // see header file for class documentation
  
  Int_t iResult = 0;

  // -- Create local stack
  if ( stack && stack->GetNtrack() > 0 )
    fStack = new TClonesArray("TParticle", stack->GetNtrack() );
  else {
    HLTError( "Error creating local stack" );
    iResult = -EINPROGRESS;
  }

  // -- Loop over off-line stack and fill local stack
  for (Int_t iterStack = 0; !iResult &&iterStack < stack->GetNtrack(); iterStack++) {

    TParticle *particle = stack->Particle(iterStack);
    if ( !particle) {
      HLTError( "Error reading particle %i out of %i", iterStack, stack->GetNtrack() );
      iResult = -EINPROGRESS;
      continue;
    }

    // ----------------
    // -- Apply cuts
    // ----------------

    if ( fHasParticleCuts ) {
      
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
	  iResult = -EINPROGRESS;
	  continue;
	}
      }
    }
    
    // -- Add particle after cuts
    AddParticle ( particle );
  }
  
  return iResult;
}

// #################################################################################
void AliHLTMCEvent::AddParticle( const TParticle* particle ) {
  // see header file for class documentation
  
  new( (*fStack) [fNParticles] ) TParticle( *particle );
  fNParticles++;

  return;
}

// #################################################################################
Int_t AliHLTMCEvent::FillMCJets( AliGenPythiaEventHeader* header ) {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Check if Header is present
  if ( !header ) {
    HLTError( "Error no PythiaHeader present" );
    iResult = -EINPROGRESS;
    return iResult;
  }

  // -- Create jet array  
  if ( header->NTriggerJets() > 0 )
    fGenJets = new TClonesArray("AliAODJet", header->NTriggerJets());
  else 
    return iResult;

  // -- Loop over jets in header and fill local array
  for (Int_t iterJet = 0; iterJet < header->NTriggerJets() && !iResult; iterJet++) {

    // -- Add jet
    AddGenJet(header, iterJet);

  } // for (Int_t iterJet = 0; iterJet < header->NTriggerJets() && !iResult; iterJet++) {

  HLTDebug("Pythia Jets found: %d", fNGenJets );
  
  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                           Pythia jets - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
AliGenPythiaEventHeader* AliHLTMCEvent::GetPythiaEventHeader(AliMCEvent *mcEvent) {
  // see header file for class documentation

  Int_t iResult = 0;

  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*> (mcEvent->GenEventHeader());

  if ( !pythiaGenHeader ) {

    // -- is cocktail header
    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(mcEvent->GenEventHeader());
    if ( !genCocktailHeader ) {
      HLTError("Error: Unknown header type (not Pythia or Cocktail)");
      iResult = -1;
    }

    if ( !iResult ) {
      // -- Get Header
      TList* headerList = genCocktailHeader->GetHeaders();

      for (Int_t iter = 0; iter < headerList->GetEntries(); iter++ ) {
	pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(iter));
	if ( pythiaGenHeader )
	  break;
      }
    }  
  }
  
  if ( pythiaGenHeader && !iResult ) 
    return pythiaGenHeader;
  else {
    HLTError("PythiaHeader not found!");
    return NULL;
  }
}

//##################################################################################
void AliHLTMCEvent::AddGenJet( AliGenPythiaEventHeader* header, Int_t iterJet ) {
  // see header file for class documentation

  Float_t pJet[] = {0., 0., 0., 0.};          // jet 4-vector ( x y z )

  // -- Get jet
  header->TriggerJet(iterJet, pJet);

  // -- Create TLorentzVector
  TLorentzVector v(pJet);
  
  // -- Add AliAODJet
  new( (*fGenJets) [fNGenJets] ) AliAODJet(v);
  fNGenJets++;

  return;
}
