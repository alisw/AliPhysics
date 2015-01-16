//-*- Mode: C++ -*-

// $Id:  $

#ifndef ALIHLTJETREADER_H
#define ALIHLTJETREADER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETReader.h
    @author Jochen Thaeder
    @date   
    @brief  Reader for jet finder
*/

#ifdef HAVE_FASTJET
#include "fastjet/PseudoJet.hh"
#endif

#include "AliJetReader.h"
#include "AliJetReaderHeader.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"

#include "AliHLTLogging.h"
#include "AliHLTMCEvent.h"

#include "AliHLTJETBase.h"
#include "AliHLTJETTrackCuts.h"
#include "AliHLTJETReaderHeader.h"

#include "AliHLTJETConeSeedCuts.h"
#include "AliHLTJETConeGrid.h"

/**
 * @class  AliHLTJETReader
 * This class is the reader class for the JetFinder in the HLT
 * It implements the reading of ESDs and MCs. AOD reading is 
 * not yet implemented
 * <br>
 * Usage :<br>
 * - Initilize() // Initializes the reader dependent of the algorithm
 *
 * @ingroup alihlt_jet
 */

class AliHLTJETReader : public AliJetReader, public AliHLTLogging {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETReader();

  /** destructor */
  virtual ~AliHLTJETReader();

  /*
   * ---------------------------------------------------------------------------------
   *                                Initialize / Reset
   * ---------------------------------------------------------------------------------
   */

  /** Initialize reader 
   *  Calls AliHLTJETReaderHeader::Initialize
   *  and the private Initialize methods
   *  @return 0 on success, otherwise <0
   */
  Int_t Initialize();
 
  /** Reset the event */
  void ResetEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set pointer to input event
   *  Needs "useMC" flag for running in analysis task only
   *  @param esd an AliESDEvent
   *  @param aod an AliAODEvent
   *  @param mc an AliHLTMCEvent
   */
  void SetInputEvent(const TObject* esd, const TObject* aod, const TObject* mc);

  /** Set number of jet candates = seeds */
  void SetNJetCandidates( Int_t i ) { fNJetCandidates = i; }

  /*
   * ---------------------------------------------------------------------------------
   *                            Fastjet Reader functionality
   * ---------------------------------------------------------------------------------
   */

#ifdef HAVE_FASTJET
  /** Fill tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillVector();

  /** Fill MC tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillVectorMC();

  /** Fill HLT MC tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillVectorHLTMC();

  /** Fill ESD tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillVectorESD();

  /** Fill AOD tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillVectorAOD();
#endif

  /*
   * ---------------------------------------------------------------------------------
   *                               Grid functionality
   * ---------------------------------------------------------------------------------
   */
  
  /** Fill tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillGrid();

  /** Fill MC tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillGridMC();

  /** Fill HLT MC tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillGridHLTMC();

  /** Fill ESD tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillGridESD();

  /** Fill AOD tracks in momentum array 
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillGridAOD();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Ptr to AliHLTJETReaderHeader
   *  @return ptr to AliHLTJETReaderHeader
   */
  AliHLTJETReaderHeader*      GetReaderHeader() const { return dynamic_cast<AliHLTJETReaderHeader*>(fReaderHeader);}

#ifdef HAVE_FASTJET
  /** Get Ptr to input vector of Fastjet
   *  @return ptr to input vector of Fastjet
   */
  vector<fastjet::PseudoJet>* GetVector()             { return fMomentumVector; }
#endif

  /** Get Ptr to grid of cone finder
   *  @return ptr to grid of cone finder
   */
  AliHLTJETConeGrid*          GetGrid()               { return fGrid; }

  /** Get number of jet candates = seeds */
  Int_t                       GetNJetCandidates()     { return fNJetCandidates; }

  /** Get ptr to jet candiates = seeds for cone finder */
  TClonesArray*               GetJetCandidates()      { return fJetCandidates; }

  /*
   * ---------------------------------------------------------------------------------
   *                                     Seeds
   * ---------------------------------------------------------------------------------
   */

  /** Add new seed
   *  @param aEtaPhi     eta and phi of the seed
   *  @param aGridIdx    indeces in the grid
   *  @param coneRadius  coneRadius
   */
  void AddSeed( const Float_t* aEtaPhi, const Int_t* aGridIdx, 
		const Float_t coneRadius );

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETReader (const AliHLTJETReader&);

  /** assignment operator prohibited */
  AliHLTJETReader& operator= (const AliHLTJETReader&);

  /*
   * ---------------------------------------------------------------------------------
   *                         Initialize - private
   * ---------------------------------------------------------------------------------
   */
  
  /** Initialize reader for the FFSC cone jet finder
   *  @return 0 on success, otherwise <0
   */
  Int_t InitializeFFSC();

#ifdef HAVE_FASTJET
  /** Initialize reader for the fastjet jet finders
   *  @return 0 on success, otherwise <0
   */
  Int_t InitializeFastjet();
#endif

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // -- Input
  // ----------

  /** ESD event */
  AliESDEvent                 *fESD;            //! transient

  /** off-line MC event */
  AliMCEvent                  *fMC;             //! transient

  /** on-line MC event */
  AliHLTMCEvent               *fHLTMC;          //! transient

  /** AOD event */
  AliAODEvent                 *fAOD;            //! transient

  // -- Particle structures
  // ------------------------

#ifdef HAVE_FASTJET
  /** Vector of fastjet momemtum entries */
  vector<fastjet::PseudoJet>  *fMomentumVector; //! transient
#endif

  /** Grid for cone finder */
  AliHLTJETConeGrid           *fGrid;           //! transient

  // -- Output 
  // -----------

  /** Number of jet candates = seeds */
  Int_t                        fNJetCandidates; // see above

  /** Jet candiates = seeds for cone finder */
  TClonesArray                *fJetCandidates;  //! transient

  // -- Cuts
  // ---------

  /** Ptr to seed cuts */
  AliHLTJETConeSeedCuts       *fSeedCuts;       //! transient

  /** Ptr to track cuts */
  AliHLTJETTrackCuts          *fTrackCuts;      //! transient

  ClassDef(AliHLTJETReader, 1)

};
#endif

