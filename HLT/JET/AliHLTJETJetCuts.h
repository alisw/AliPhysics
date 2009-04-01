//-*- Mode: C++ -*-

// $Id: AliHLTJETJetCuts.h $

#ifndef ALIHLTJETJETCUTS_H
#define ALIHLTJETJETCUTS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETJetCuts.h
    @author Jochen Thaeder
    @date   
    @brief  Cuts for jet input tracks
*/

#include "AliAODJet.h"

#include "AliAnalysisCuts.h"
#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

#include "AliHLTJETConeJetCandidate.h"

/**
 * @class  AliHLTJETJetCuts
 * Cuts for jets
 *
 * @ingroup alihlt_jet
 */

class AliHLTJETJetCuts : public AliAnalysisCuts, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETJetCuts(const Char_t* name = "AliHLTJETJetCuts", const Char_t* title = "");

  /** destructor */
  virtual ~AliHLTJETJetCuts();

  /*
   * ---------------------------------------------------------------------------------
   *                                   Selection
   * ---------------------------------------------------------------------------------
   */
  
  /** Select jet
      @param obj AliHLTJETConeJetCandidate or AliAODJet
      @return kTRUE if selected, kFALSE otherwise
  */
  Bool_t IsSelected( TObject* obj );

  /** Select jet
      @param jet AliHLTJETConeJetCandidate jet
      @return kTRUE if selected, kFALSE otherwise
  */
  Bool_t IsSelected( AliHLTJETConeJetCandidate* jet );

  /** Select jet
      @param jet AliAODJet jet
      @return kTRUE if selected, kFALSE otherwise
  */
  Bool_t IsSelected( AliAODJet* jet );
  
  /** Select jet
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
  void SetMinEt( Float_t f ) { fEtMin = f; }

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETJetCuts (const AliHLTJETJetCuts&);

  /** assignment operator prohibited */
  AliHLTJETJetCuts& operator= (const AliHLTJETJetCuts&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Minimum Et */
  Float_t            fEtMin;                    // see above

#if 0
  /** Distance between to jets */
  Float_t            fDistanceCut;
#endif

  ClassDef(AliHLTJETJetCuts, 1)

};
#endif

