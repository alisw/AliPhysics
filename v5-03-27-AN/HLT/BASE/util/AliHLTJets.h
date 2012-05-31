//-*- Mode: C++ -*-

// $Id: AliHLTJets.h  $

#ifndef ALIHLTJETS_H
#define ALIHLTJETS_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJets.h
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
#include "TString.h"

#include "AliAODJet.h"

#include "AliHLTLogging.h"

/**
 * @class AliHLTJets
 * This class contains AliAODJets and comments
 *
 * @ingroup alihlt_jet
 */

class AliHLTJets : public TObject, public AliHLTLogging {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Constructor */
  AliHLTJets();
  
  /** Destructor */
  ~AliHLTJets();

  /*
   * ---------------------------------------------------------------------------------
   *                                   Initialize / Reset
   * ---------------------------------------------------------------------------------
   */

  /** Reset output array */
  void Reset();

  /** Reset output array */
  void ResetEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Comment */
  TString       GetComment()          const { return fComment; }
  
  /** Get number of AliAODJets */
  Int_t         GetNAODJets()         const { return fNAODJets; }

  /** Get ptr to TClonesArray of AliAODJets */
  TClonesArray* GetAODJets()          const { return fAODJets; }

  /** Get AliAODJet with idx iter */
  AliAODJet*    GetJet( Int_t iter )  const;

  /** Get next AliAODJet 
   *  @return Ptr to Jet, NULL if no next jet is present
   */  
  AliAODJet*    NextJet();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */

  /** Set Comment
   * @param   comment  arguments of jet finder
   */
  void SetComment( const Char_t* comment ) { fComment = TString(comment); }

  /** Set Comment
   * @param   comment  arguments of jet finder
   */
  void SetComment( TString comment ) { fComment = TString(comment); }

  /** Add Jet to Container 
   *  @param  eta     Jet eta
   *  @param  phi     Jet phi
   *  @param  pt      Jet pt
   *  @param  et      Jet et
   */
  void AddJet( Float_t eta, Float_t phi, Float_t pt, Float_t et );

  /** Add Jet to Container 
   *  @param  jet     Ptr to AliAODJet
   */
  void AddJet( AliAODJet* jet );

  /*
   * ---------------------------------------------------------------------------------
   *                                     Helper
   * ---------------------------------------------------------------------------------
   */

  /** Sort Jets with decreasing Et */
  void Sort();

  ///////////////////////////////////////////////////////////////////////////////////
  
 private:
 
  /** copy constructor prohibited */
  AliHLTJets(const AliHLTJets&);
  
  /** assignment operator prohibited */
  AliHLTJets& operator=(const AliHLTJets&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** String containing comment */
  TString          fComment;           // see above 

  /** Current Jet index */
  Int_t            fCurrentJetIndex;   //! transient

  /** Number of AOD jets */
  Int_t            fNAODJets;          // see above

  /** Array of AOD jets */
  TClonesArray    *fAODJets;           // see above
  
  ClassDef(AliHLTJets, 1)
};
#endif
