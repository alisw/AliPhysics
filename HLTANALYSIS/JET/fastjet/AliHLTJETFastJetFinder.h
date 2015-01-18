//-*- Mode: C++ -*-

// $Id: AliHLTJETFastJetFinder.h  $

#ifndef ALIHLTJETFASTJETFINDER_H
#define ALIHLTJETFASTJETFINDER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETFastJetFinder.h
    @author Jochen Thaeder
    @date   
    @brief  FastJet finder interface
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "AliJetFinder.h"

#include "AliHLTJets.h"
#include "AliHLTLogging.h"

#include "AliHLTJETBase.h"

/**
 * @class  AliHLTJETFastJetFinder
 * FastJet Interface for the fasjet package ( v.2.4.1 )
 *
 * <b>Usage in off-line</b><br>
 *  * Initialization phase :
 *      <pre>jetFinder->Init();</pre> 
 *  * Set the input event via the reader
 *      <pre>jetReader->SetInputEvent( ... )</pre>
 *  * Process one event (contains reset per event)
 *      <pre>jetFinder->ProcessEvent();</pre>
 * 
 * <b>Usage in on-line</b><br>
 *  * Initialization phase :
 *      <pre>jetFinder->Initialize();</pre> 
 *  * Set the input event via the reader
 *      <pre>jetReader->SetInputEvent( ... )</pre>
 *  * Process one event 
 *      * Fill input vector (contains reset per event)
 *           <pre>jetReader->FillVectorXXX();</pre> 
 *           Where XXX is has to be replaced by MC, ESD or AOD, 
 *           depending, on the input object 
 *      * Process one event (contains reset per event)
 *          <pre>jetFinder->ProcessHLTEvent();</pre>
 *   
 * @ingroup alihlt_jet_fastjet
 */

class AliHLTJETFastJetFinder : public AliJetFinder, public AliHLTLogging {
 
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETFastJetFinder();

  /** destructor */
  virtual ~AliHLTJETFastJetFinder();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Initialize
   * ---------------------------------------------------------------------------------
   */
  
  /** Initialize the jet finder and the search grid
   *  ONLY for use in off-line ... it inherits from a virtual void ?!?!
   */
  void Init() { Initialize(); }

  /** Initialize the jet finder and the search grid
   *  @return 0 on success, < 0 on failure
   */
  Int_t Initialize();

  /** Reset for next event */
  void Reset();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set ptr to output container */
  void SetOutputJets( AliHLTJets* jets ) { fJets = jets; }

  /*
   * ---------------------------------------------------------------------------------
   *                                      Process
   * ---------------------------------------------------------------------------------
   */
  
  /** Process one event
   *  ONLY for use in off-line ... it inherits from a virtual Bool_t
   *  @return kTRUE on success, kFALSE on failure
   */
  Bool_t ProcessEvent();
  
  /** Process one event
   *  @return kTRUE on success, kFALSE on failure
   */
  Bool_t ProcessHLTEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                                      Helper
   * ---------------------------------------------------------------------------------
   */

  /** Print found jets */
  void PrintJets(vector<fastjet::PseudoJet> &jets, fastjet::ClusterSequenceArea &clust_seq);

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETFastJetFinder (const AliHLTJETFastJetFinder&);

  /** assignment operator prohibited */
  AliHLTJETFastJetFinder& operator= (const AliHLTJETFastJetFinder&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Process - private
   * ---------------------------------------------------------------------------------
   */

  /** Find jets, fill jets and apply jet cuts in one event
   *  @return 0 on success, < 0 on failure
   */
  Int_t FindFastJets(); 

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Input vector for fastJet */
  vector<fastjet::PseudoJet>  *fInputVector;     //! transient

  /** Container of AliAODJets */
  AliHLTJets                  *fJets;            //! transient

  ClassDef(AliHLTJETFastJetFinder,1)
};

#endif
