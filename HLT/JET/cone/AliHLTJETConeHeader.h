//-*- Mode: C++ -*-

// $Id: AliHLTJETConeHeader.h $

#ifndef ALIHLTJETCONEHEADER_H
#define ALIHLTJETCONEHEADER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETConeHeader.h
    @author Jochen Thaeder
    @date   
    @brief  Header of the cone finder
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliJetHeader.h"

#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

#include "AliHLTJETJetCuts.h"

/**
 * @class  AliHLTJETConeHeader
 * Header of the cone finder
 *
 * @ingroup alihlt_jet_cone
 */

class AliHLTJETConeHeader : public AliJetHeader, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETConeHeader();

  /** destructor */
  virtual ~AliHLTJETConeHeader();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Initialize
   * ---------------------------------------------------------------------------------
   */

  /** Initialize the jet header
   *  @return 0 on success, < 0 on failure
   */
  Int_t Initialize();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set Analysis Cuts
   *  @param cuts ptr to AliHLTJETJetCuts 
   */
  void SetJetCuts( AliHLTJETJetCuts* cuts )  { fJetCuts = cuts; }
 
  /** Set flag to use only leading seed */
  void SetUseLeading( Bool_t b ) { fUseLeading = b; }

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Analysis Cuts
   *  @return ptr to AliHLTJETJetCuts 
   */
  AliHLTJETJetCuts* GetJetCuts()     { return fJetCuts; }

  /** Get flag to use only leading seed 
   *  @return   if kTRUE, only leading seed is used
   */
  Bool_t GetUseLeading() { return fUseLeading; }

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETConeHeader (const AliHLTJETConeHeader&);

  /** assignment operator prohibited */
  AliHLTJETConeHeader& operator= (const AliHLTJETConeHeader&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Cuts on jet selection */
  AliHLTJETJetCuts          *fJetCuts;              //! transient
  
  /** if kTRUE, only leading seed is used */
  Bool_t                     fUseLeading;           // see above

  /*
    Int_t                      fgFinderType;
  */

  ClassDef(AliHLTJETConeHeader, 1)

};
#endif

