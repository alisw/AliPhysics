//-*- Mode: C++ -*-

// $Id:  $

#ifndef ALIHLTJETTRACKCUTS_H
#define ALIHLTJETTRACKCUTS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETTrackCuts.h
    @author Jochen Thaeder
    @date   
    @brief  Cuts for jet input tracks
*/

#include "AliAnalysisCuts.h"
#include "AliHLTLogging.h"

#include "AliESDtrack.h"
#include "TParticle.h"

/**
 * @class  AliHLTJETTrackCuts
 * Cuts for MC tracks and ESD tracks 
 *
 * @ingroup alihlt_jet
 */

class AliHLTJETTrackCuts : public AliAnalysisCuts, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETTrackCuts(const Char_t* name = "AliHLTJETTrackCuts", const Char_t* title = "");

  /** destructor */
  virtual ~AliHLTJETTrackCuts();

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
  

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETTrackCuts (const AliHLTJETTrackCuts&);

  /** assignment operator prohibited */
  AliHLTJETTrackCuts& operator= (const AliHLTJETTrackCuts&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Minimum pt  */
  Double_t           fPtMin;                   // see above

  ClassDef(AliHLTJETTrackCuts, 1)

};
#endif

