//-*- Mode: C++ -*-

// $Id: AliHLTJETFastJetHeader.h  $

#ifndef ALIHLTJETFASTJETHEADER_H
#define ALIHLTJETFASTJETHEADER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETFastJetHeader.h
    @author Jochen Thaeder
    @date   
    @brief  Header of the FastJet finder interface  
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "fastjet/PseudoJet.hh"

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"

#include "AliJetHeader.h"
#include "AliJetReaderHeader.h"
#include "AliHLTLogging.h"

#include "AliHLTJETBase.h"

/**
 * @class  AliHLTJETFastJetHeader
 * FastJet Interface for the fasjet package ( v.2.4.1 )
 *
 *   
 * @ingroup alihlt_jet_fastjet
 */

class AliHLTJETFastJetHeader : public AliJetHeader, public AliHLTLogging {
 
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETFastJetHeader();

  /** destructor */
  virtual ~AliHLTJETFastJetHeader();

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

  /** Set Reader Header
   *  @param rh ptr to Analysis Cuts
   */
  void SetReaderHeader(AliJetReaderHeader* rh ) { fReaderHeader = rh; }

  /** Set Analysis Cuts
   *  @param cuts ptr to AliHLTJETJetCuts 
   */
  void SetJetCuts( AliHLTJETJetCuts* cuts )     { fJetCuts = cuts; }

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Min Pt for jets */
  Double_t                     GetPtMin()             const {return fPtMin;}

  /** Get Analysis Cuts
   *  @return ptr to AliHLTJETJetCuts 
   */
  AliHLTJETJetCuts*            GetJetCuts()                 { return fJetCuts; }

  // -- FastJet Getters
  // --------------------

  /** Get jet definition */
  fastjet::JetDefinition*      GetJetDefinition()           {return fJetDefinition;}

  /** Get area definition */
  fastjet::AreaDefinition*     GetAreaDefinition()          {return fAreaDefinition;}

  /** Get range definition */
  fastjet::RangeDefinition*    GetRangeDefinition()         {return fRangeDefinition;}

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Print fastjet parameters */
  void PrintParameters() const;

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETFastJetHeader (const AliHLTJETFastJetHeader&);

  /** assignment operator prohibited */
  AliHLTJETFastJetHeader& operator= (const AliHLTJETFastJetHeader&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // -- fastjet::JetDefinition parameters
  // --------------------------------------
  
  /** Jet Algorithm 
   *  Default : fastjet::kt_algorithm
   */
  fastjet::JetAlgorithm        fAlgorithm;         // see above

  /** Finding strategy 
   *  Default : fastjet::Best
   */
  fastjet::Strategy            fStrategy;          // see above

  /** Recombination scheme 
   *  Default : fastjet::BIpt_scheme;
   */
  fastjet::RecombinationScheme fRecombScheme;      // see above

  // -- fastjet::GhostedAreaSpec parameters
  // ----------------------------------------

  /** Ghost area */
  Double_t                     fGhostArea;         // see above

  /** Ghost area - active repeats */
  Int_t                        fActiveAreaRepeats; // see above

  // -- fastjet::AreaDefinition parameters
  // ---------------------------------------

  /** Area Definition */
  fastjet::AreaType            fAreaType;          // see above

  // -- fastjet::ClusterSequenceArea options parameters
  // ----------------------------------------------------

  /** Jet pt > ptmin */
  Double_t                     fPtMin;             // see above

  // -- fastjet classes
  // --------------------
  
  /** Jet definition */
  fastjet::JetDefinition      *fJetDefinition;     //! transient

  /** Ghost definition */
  fastjet::GhostedAreaSpec    *fGhostedAreaSpec;   //! transient

  /** Area definition */
  fastjet::AreaDefinition     *fAreaDefinition;    //! transient

  /** Range definition */
  fastjet::RangeDefinition    *fRangeDefinition;   //! transient

  // -- Ptr to classes
  // -------------------

  /** Ptr to jet reader header */
  AliJetReaderHeader          *fReaderHeader;      //! transient

  /** Cuts on jet selection */
  AliHLTJETJetCuts            *fJetCuts;           //! transient


  ClassDef(AliHLTJETFastJetHeader,1)
};
 
#endif
