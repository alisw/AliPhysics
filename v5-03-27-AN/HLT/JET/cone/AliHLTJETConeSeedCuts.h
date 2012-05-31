//-*- Mode: C++ -*-

// $Id: AliHLTJETConeSeedCuts.h $

#ifndef ALIHLTJETCONESEEDCUTS_H
#define ALIHLTJETCONESEEDCUTS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETConeSeedCuts.h
    @author Jochen Thaeder
    @date   
    @brief  Cuts for jet input tracks
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TParticle.h"

#include "AliAnalysisCuts.h"
#include "AliESDtrack.h"

#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

/**
 * @class  AliHLTJETConeSeedCuts
 * Cuts for seed MC tracks and ESD tracks 
 *
 * @ingroup alihlt_jet_cone
 */

class AliHLTJETConeSeedCuts : public AliAnalysisCuts, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTJETConeSeedCuts(const Char_t* name = "AliHLTJETConeSeedCuts", 
			const Char_t* title = "");

  /** destructor */
  virtual ~AliHLTJETConeSeedCuts();

  /*
   * ---------------------------------------------------------------------------------
   *                                   Selection
   * ---------------------------------------------------------------------------------
   */
  
  /** Select track
      @param obj esd track or particle
      @return kTRUE if selected, kFALSE otherwise
  */
  Bool_t IsSelected( TObject* obj );

  /** Select track
      @param particle particle
      @return kTRUE if selected, kFALSE otherwise
  */
  Bool_t IsSelected( TParticle* particle );

  /** Select track
      @param esdTrack esd track 
      @return kTRUE if selected, kFALSE otherwise
  */
  Bool_t IsSelected( AliESDtrack* esdTrack );
  
  /** Select track
      Not implemented
      @return kTRUE 
  */
  Bool_t IsSelected( TList* /*list*/ ) { return kTRUE; }

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set cut on min pt */
  void SetMinPt( Float_t f ) { fPtMin = f; }

  /** Set cut on eta acceptance */
  void SetEtaRange( Float_t etaMin, Float_t etaMax ) { fEtaMin = etaMin; fEtaMax = etaMax; }

  /** Set cut on phi acceptance */
  void SetPhiRange( Float_t phiMin, Float_t phiMax ) { fPhiMin = phiMin; fPhiMax = phiMax; }

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETConeSeedCuts (const AliHLTJETConeSeedCuts&);

  /** assignment operator prohibited */
  AliHLTJETConeSeedCuts& operator= (const AliHLTJETConeSeedCuts&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Minimum pt */
  Float_t            fPtMin;                    // see above

  /** Minimum eta */
  Float_t            fEtaMin;                   // see above

  /** Maximum eta */
  Float_t            fEtaMax;                   // see above

  /** Minimum phi */
  Float_t            fPhiMin;                   // see above

  /** Maximum phi */
  Float_t            fPhiMax;                   // see above

  ClassDef(AliHLTJETConeSeedCuts, 1)

};
#endif

