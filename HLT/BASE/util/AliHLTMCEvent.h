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

#include "AliHLTLogging.h"
#include "AliMCEvent.h"

/**
 * @class  AliHLTMCEvent
 * Container class for off-line AliMCEvent, as this class is coupled to 
 * TTrees an can not be easily copied and send via the HLT chain.
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

  /** constructor 
   *  @param iNumberTracks number of initial tracks
   */
  AliHLTMCEvent( Int_t iNumberTracks );
 
  /** constructor 
   *  @param pMCEvent ptr to off-line AliMCEvent
   */
  AliHLTMCEvent( AliMCEvent *pMCEvent );

  /** destructor */
  virtual ~AliHLTMCEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

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

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Add a prticle to this container
   *  @param particle ptr to TParticle classs
   */
  void AddParticle( const TParticle* particle);
  
  /** Fill the off-line MC event in to this class
   *  @param pMCEvent ptr to off-line AliMCEvent
   */
  void FillMCEvent( AliMCEvent *pMCEvent );
  
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
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Number of particles */
  Int_t           fNParticles;                   // see above

  /** Current particle */
  Int_t           fCurrentParticleIndex;         // see above

  /** Ptr to current particle */
  TParticle      *fCurrentParticle;              // see above

  /** Stack of particles */
  TClonesArray   *fStack;                        // see above

  ClassDef(AliHLTMCEvent, 1)

};
#endif

