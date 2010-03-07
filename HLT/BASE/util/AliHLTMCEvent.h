//-*- Mode: C++ -*-

// $Id:  $

#ifndef ALIHLTMCEVENT_H
#define ALIHLTMCEVENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTMCEvent.h
    @author Jochen Thaeder
    @date   
    @brief  Container class for an AliMCEvent
*/

#include "TObject.h"
#include "TParticle.h"
#include "TClonesArray.h"

#include "AliGenPythiaEventHeader.h"
#include "AliMCEvent.h"
#include "AliAODJet.h"
#include "AliStack.h"
#include "AliHeader.h"

#include "AliHLTLogging.h"

/**
 * @class  AliHLTMCEvent
 * Container class for off-line AliMCEvent, as the class AliMCEvent 
 * is coupled to TTrees and can not be easily copied and send via 
 * the HLT chain.
 *
 * This class as the complete functionality as the off-line class,
 * only the ones needed so far.
 *
 * There is no extra "stack" class, the stack is included in this class 
 * a TClonesArray of TParticles
 *
 */

class AliHLTMCEvent : public TObject, public AliHLTLogging  {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTMCEvent();

  /** standard constructor */
  AliHLTMCEvent( Bool_t applyParticleCuts );
  
  /** destructor */
  virtual ~AliHLTMCEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                               Setter - public
   * ---------------------------------------------------------------------------------
   */

  /** Apply particle cuts to the this event */
  void ApplyParticleCuts() { fHasParticleCuts = kTRUE; }

  /** Fill the off-line MC event in to this class
   *  @param stack  ptr to AliStack
   *  @param header ptr to AliHeader
   *  @return 0 on sucess, <0 on error
   */
  Int_t FillMCEvent( AliStack *stack, AliHeader *header );

  /** Fill the off-line MC event in to this class
   *  @param pMCEvent ptr to off-line AliMCEvent
   *  @return 0 on sucess, <0 on error
   */
  Int_t FillMCEvent( AliMCEvent *pMCEvent );

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  // -- Particles
  // --------------

  /** Get number of particles
   *  @return number of particles
   */
  Int_t GetNumberOfTracks() const { return fNParticles; }

  /** Return particle at index iParticle
   *  @param iParticle Particle index in local stack
   *  @return ptr on success, NULL on failure
   */
  TParticle* Particle( Int_t iParticle );
  
  /** Return next particle
   *  @return ptr on success, NULL on failure
   */
  TParticle* NextParticle();

  /** Get Index of current particle 
   *  @return Index of current particle
   */
  Int_t GetIndex() { return fCurrentParticleIndex; }

  /** Check if particle cuts have been applied */
  Bool_t HasParticleCuts() { return fHasParticleCuts; }

  // -- Generated jets
  // -------------------

  /** Get number of generated jets
   *  @return number of generated jets
   */
  Int_t GetNumberOfGenJets() const { return fNGenJets; }

  /** Return generated jets at index iJet
   *  @param iJet Generated jet index in local array
   *  @return ptr on success, NULL on failure
   */
  AliAODJet* GenJet( Int_t iJet ) const;
  
  /** Return next generated jet
   *  @return ptr on success, NULL on failure
   */
  AliAODJet* NextGenJet() ;

  /*
   * ---------------------------------------------------------------------------------
   *                                    Helper
   * ---------------------------------------------------------------------------------
   */

  /** Compress the TClonesArray fStack */
  void Compress();

  /** Reset index for next particle */
  void Reset();

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTMCEvent (const AliHLTMCEvent&);

  /** assignment operator prohibited */
  AliHLTMCEvent& operator= (const AliHLTMCEvent&);

  /*
   * ---------------------------------------------------------------------------------
   *                               Setter - private
   * ---------------------------------------------------------------------------------
   */

  // -- Particles
  // --------------
 
  /** Fill the particles in into the local stack
   *  @param stack ptr to off-line AliStack
   *  @return         0 on sucess, <0 on error
   */
  Int_t FillMCTracks( AliStack *stack );
  
  /** Add a particle to this container
   *  @param particle ptr to TParticle classs
   */
  void AddParticle( const TParticle* particle);
  
  // -- Generated jets
  // -------------------
  
  /** Fill the jets into local array
   *  @param stack ptr to off-line AliGenPythiaEventHeader
   *  @return         0 on sucess, <0 on error
   */
  Int_t FillMCJets( AliGenPythiaEventHeader* header );

  /** Add a jet to his container 
   *  @param stack    ptr to off-line AliGenPythiaEventHeader
   *  @param iterJet  idx of jet in event
   */
  void AddGenJet( AliGenPythiaEventHeader* header, Int_t iterJet );

  /*
   * ---------------------------------------------------------------------------------
   *                           Pythia jets - private
   * ---------------------------------------------------------------------------------
   */

  /** Retrieve pythia event header out of an AliMCEvent 
   *  @param mcEvent   ptr to AliMCEvent
   *  @return          ptr on sucess, NULL on error   
   */
  AliGenPythiaEventHeader* GetPythiaEventHeader(AliMCEvent *mcEvent);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // -- Particles
  // --------------

  /** Current particle */
  Int_t           fCurrentParticleIndex;         // see above

  /** Number of particles */
  Int_t           fNParticles;                   // see above

  /** Stack of particles [TParticle]*/
  TClonesArray   *fStack;                        // see above

  // -- Generated jets
  // -------------------

  /** Current generated jet */
  Int_t           fCurrentGenJetIndex;           // see above

  /** Number of generated jets */
  Int_t           fNGenJets;                     // see above
  
  /** Array of generated jets [AliAODJet]*/
  TClonesArray   *fGenJets;                      // see above

  // -- Status Flags
  // -----------------

  /** Particle cuts have been applied 
   *  - Is primary
   *  - Is final state
   *  - Is known to PDG datebase
   */
  Bool_t fHasParticleCuts;                       // see above

  ClassDef(AliHLTMCEvent, 1)

};
#endif

