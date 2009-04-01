//-*- Mode: C++ -*-

// $Id: AliHLTJETJets.h  $

#ifndef ALIHLTJETJETS_H
#define ALIHLTJETJETS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETJets.h
    @author Jochen Thaeder
    @date   
    @brief  Container holding produced Jets
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TClonesArray.h"

#include "AliAODJet.h"

#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

#include "AliHLTJETConeJetCandidate.h"

/**
 * @class AliHLTJETJets
 * This class contains AliAODJets and comments
 *
 * @ingroup alihlt_jet
 */

class AliHLTJETJets : public TObject, public AliHLTLogging {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Constructor */
  AliHLTJETJets();
  
  /** Destructor */
  ~AliHLTJETJets();

  /*
   * ---------------------------------------------------------------------------------
   *                                   Initialize / Reset
   * ---------------------------------------------------------------------------------
   */

  /** Reset output array */
  void Reset();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */
  
  /* Get number of AODJets */
  Int_t         GetNAODJets() { return fNAODJets; }

  /** Get ptr to AODJets */
  TClonesArray* GetAODJets()  { return fAODJets; }

  /** Get AODHet with idx iter */
  AliAODJet*    GetJet( Int_t iter );


  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */

  /** Add Jet to Container 
   * @param  jetCandidate  ptr to AliHLTJetFinderJetSeed
   */
  void AddJet( AliHLTJETConeJetCandidate* jet );

  ///////////////////////////////////////////////////////////////////////////////////
  
 private:
 
  /** copy constructor prohibited */
  AliHLTJETJets(const AliHLTJETJets&);
  
  /** assignment operator prohibited */
  AliHLTJETJets& operator=(const AliHLTJETJets&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Number of AOD jets */
  Int_t                   fNAODJets;          // see above

  /** Array of AOD jets */
  TClonesArray*           fAODJets;        // see above

  ClassDef(AliHLTJETJets, 1)
};
#endif
