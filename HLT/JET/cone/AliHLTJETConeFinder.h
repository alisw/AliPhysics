//-*- Mode: C++ -*-

// $Id: AliHLTJETConeFinder.h  $

#ifndef ALIHLTJETCONEFINDER_H
#define ALIHLTJETCONEFINDER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETConeFinder.h
    @author Jochen Thaeder
    @date   
    @brief  Jet cone finder
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliJetFinder.h"

#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

#include "AliHLTJETJets.h"
#include "AliHLTJETConeGrid.h"

/**
 * @class  AliHLTJETConeFinder
 * ConeFinder for jet finder
 *
 * <b>Usage in off-line</b><br>
 *  * Initialization phase :
 *      <pre>jetFinder->Init();</pre> 
 *  * Set the input event via the reader
 *      <pre>jetReader->SetInputEvent( ... )</pre>
 *  * Process one event
 *      <pre>jetFinder->ProcessEvent();</pre>
 *
 * <b>Usage in on-line</b><br>
 *  * Initialization phase :
 *      <pre>jetFinder->Initialize();</pre> 
 *  * Reset
        <pre>jetFinder->Reset();</pre>
 *  * Set the input event via the reader
 *      <pre>jetReader->SetInputEvent( ... )</pre>
 *  * Process one event
 *      <pre>jetReader->FillGridXXX();</pre> 
 *       Where XXX is has to be replaced by MC, ESD or AOD, depending
 *       on the input object
 *      <pre>jetFinder->ProcessConeEvent();</pre>
 *   
 * @ingroup alihlt_jet_cone
 */

class AliHLTJETConeFinder : public AliJetFinder, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETConeFinder();

  /** destructor */
  virtual ~AliHLTJETConeFinder();

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
  void SetOutputJets( AliHLTJETJets* jets ) { fJets = jets; }
  
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
  Bool_t ProcessConeEvent();

  /** Find jets in one event
   *  ONLY for use in off-line ... it inherits from a virtual void ?!?!
   */
  void FindJets() { FindConeJets(); }

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETConeFinder (const AliHLTJETConeFinder&);

  /** assignment operator prohibited */
  AliHLTJETConeFinder& operator= (const AliHLTJETConeFinder&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Process - private
   * ---------------------------------------------------------------------------------
   */

  /** Sort seed descending in pt
   *  if header::fUseLeading only is set, keep only seed with highest pt
   *  @return 0 on success, < 0 on failure
   */ 
  Int_t FindConeLeading();

  /** Find jets in one event
   *  @return 0 on success, < 0 on failure
   */
  Int_t FindConeJets();

  /** Fill jets in output container and apply cuts
   *  @return 0 on success, < 0 on failure
   */
  Int_t FillConeJets();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Grid for cone finder */
  AliHLTJETConeGrid           *fGrid;           //! transient

  AliHLTJETJets               *fJets;           //! transient

  ClassDef(AliHLTJETConeFinder, 1)

};
#endif

