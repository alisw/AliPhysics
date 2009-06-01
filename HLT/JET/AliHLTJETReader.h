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
 * Reader for jet finder
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


  void SetTrackCuts( AliHLTJETTrackCuts * cuts) {fTrackCuts = cuts; }

  /** Initialize reader for cone jet finder
   *  Calls AliHLTJETReaderHeader::Initialize
   *  @return 0 on success, otherwise <0
   */
  Int_t Initialize();
 
  /** Reset the event */
  void ResetEvent();

  /*
   * ---------------------------------------------------------------------------------
   *                            Fastjet Reader functionality
   * ---------------------------------------------------------------------------------
   */

#ifdef HAVE_FASTJET
  /** Fill tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArrayFast();

  /** Fill MC tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArrayFastMC();

  /** Fill ESD tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArrayFastESD();

  /** Fill AOD tracks in fastjet momemtum vector
   *  @return kTRUE on success, otherwise kFALSE
   */
  Bool_t FillMomentumArrayFastAOD();
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
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set pointer to input event
   *  @param esd an AliESDEvent
   *  @param aod an AliAODEvent
   *  @param mc an AliHLTMCEvent
   */
  void SetInputEvent(TObject* esd, TObject* aod, TObject* mc);

  /** Set number of jet candates = seeds */
  void SetNJetCandidates( Int_t i ) { fNJetCandidates = i; }

  /*
   * ---------------------------------------------------------------------------------
   *                                     Getter
   * ---------------------------------------------------------------------------------
   */

  /** Get Ptr to AliHLTJETReaderHeader
   *  @return ptr to AliHLTJETReaderHeader
   */
  AliHLTJETReaderHeader* GetReaderHeader() { return dynamic_cast<AliHLTJETReaderHeader*>(fReaderHeader);}

#ifdef HAVE_FASTJET
  /** Get Ptr to input vector of Fastjet
   *  @return ptr to input vector of Fastjet
   */
  vector<fastjet::PseudoJet>* GetMomentumVectorFast() { return fMomentumVector; }
#endif

  /** Get Ptr to grid of cone finder
   *  @return ptr to grid of cone finder
   */
  AliHLTJETConeGrid* GetGrid() { return fGrid; }
  
  /** Get number of jet candates = seeds */
  Int_t         GetNJetCandidates() { return fNJetCandidates; }

  /** Get ptr to jet candiates = seeds for cone finder */
  TClonesArray* GetJetCandidates()  { return fJetCandidates; }

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
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** ESD event */
  AliESDEvent                 *fESD;            //! transient

  /** MC event */
  AliHLTMCEvent               *fMC;             //! transient

  /** AOD event */
  AliAODEvent                 *fAOD;            //! transient

#ifdef HAVE_FASTJET
  /** Vector of fastjet momemtum entries */
  vector<fastjet::PseudoJet>  *fMomentumVector; //! transient
#endif

  /** Grid for cone finder */
  AliHLTJETConeGrid           *fGrid;           //! transient

  /** Number of jet candates = seeds */
  Int_t                        fNJetCandidates; // see above

  /** Jet candiates = seeds for cone finder */
  TClonesArray                *fJetCandidates;  //! transient

  /** Ptr to seed cuts */
  AliHLTJETConeSeedCuts       *fSeedCuts;       //! transient

  /** Ptr to track cuts */
  AliHLTJETTrackCuts          *fTrackCuts;      //! transient

  ClassDef(AliHLTJETReader, 1)

};
#endif

