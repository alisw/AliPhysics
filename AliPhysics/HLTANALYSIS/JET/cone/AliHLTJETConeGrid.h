//-*- Mode: C++ -*-

// $Id: AliHLTJETConeGrid.h  $

#ifndef ALIHLTJETCONEGRID_H
#define ALIHLTJETCONEGRID_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETConeGrid.h
    @author Jochen Thaeder
    @date   
    @brief  Eta-Phi grid of the cone finder
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "TClonesArray.h"
#include "TParticle.h"

#include "AliESDtrack.h"

#include "AliHLTLogging.h"
#include "AliHLTJETBase.h"

/**
 * @class  AliHLTJETConeGrid
 * Eta-Phi grid of the cone finder
 *
 * @ingroup alihlt_jet_cone
 */

class AliHLTJETConeGrid : public AliHLTLogging, public TObject  {
  
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** standard constructor */
  AliHLTJETConeGrid();

  /** destructor */
  virtual ~AliHLTJETConeGrid();

  /*
   * ---------------------------------------------------------------------------------
   *                                    Process
   * ---------------------------------------------------------------------------------
   */

  /** Fill particle into grid -> into cell
   *  retrieve the grid indeces and get (eta,phi,pt) of the particle
   *  @param esdTrack ptr to track
   *  @param aEtaPhi  array to be filled with (eta,phi,pt) of the particle
   *  @param aGridIdx array to be filled with grid indeces
   *  @return 0 on sucess, < 0 for error
   */
  Int_t FillTrack( TParticle* particle, const Float_t* aEtaPhi, Int_t* aGridIdx );

  /** Fill track into grid -> into cell
   *  retrieve the grid indeces and get (eta,phi,pt) of the track
   *  @param esdTrack ptr to track
   *  @param aEtaPhi  array to be filled with (eta,phi,pt) of the track
   *  @param aGridIdx array to be filled with grid indeces
   *  @return 0 on sucess, < 0 for error
   */
  Int_t FillTrack( AliESDtrack* esdTrack, const Float_t* aEtaPhi, Int_t* aGridIdx );

  /*
   * ---------------------------------------------------------------------------------
   *                                   Initialize / Reset
   * ---------------------------------------------------------------------------------
   */

  /** Initialize grid 
   *  @return 0 on sucess, < 0 for error
   */
  Int_t Initialize();

  /** Reset grid */
  void Reset();

  /*
   * ---------------------------------------------------------------------------------
   *                                     Setter
   * ---------------------------------------------------------------------------------
   */
  
  /** Set eta range parameters */
  void SetEtaRange( Float_t etaMin, Float_t etaMax, Float_t etaRange ) 
  { fEtaMin = etaMin; fEtaMax = etaMax; fEtaRange = etaRange; }

  /** Set eta range parameters */
  void SetPhiRange( Float_t phiMin, Float_t phiMax, Float_t phiRange ) 
  { fPhiMin = phiMin; fPhiMax = phiMax; fPhiRange = phiRange; }

  /** Set binning */
  void SetBinning( Float_t etaBinning, Float_t phiBinning ) 
  { fEtaBinning = etaBinning; fPhiBinning = phiBinning; }
  
  /** Set cone radius */
  void SetConeRadius( Float_t coneRadius) { fConeRadius = coneRadius; }

  /*
   * ---------------------------------------------------------------------------------
   *                             Helper - public
   * ---------------------------------------------------------------------------------
   */

  /** Next cell
   *  Returns idx for the next cell in around seed
   *  @return idx, -1 if false
   */
  Int_t NextCell();
  
  /** Seed cell iterator
   *  Sets cell iterator for the NextCell method for
   *  a given seed
   *  @param etaIdx Eta index of seed
   *  @param phiIdx Phi index of seed
   */
  void SetCellIter( const Int_t etaIdx, const Int_t phiIdx );


  /** Check if there is an object at cellIdx in fGrid
   *  @param   cellIdx    CellIdx where there coulf be an object
   *  @return             ptr to cell, NULL if empty
   */
  TObject* UncheckedAt( Int_t cellIdx ) { return fGrid->UncheckedAt(cellIdx); }
  

  //  reinterpret_cast<AliHLTJETConeEtaPhiCell*>((*fGrid)[cellIdx])
  //if ( ( iResult = jet->AddCell(fGrid->GetCell(cellIdx)) ) ) {



  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTJETConeGrid (const AliHLTJETConeGrid&);

  /** assignment operator prohibited */
  AliHLTJETConeGrid& operator= (const AliHLTJETConeGrid&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Helper - private
   * ---------------------------------------------------------------------------------
   */

  /** Get Cell index out of eta and phi
   *  @param aEtaPhi  input structure containing eta,phi
   *                  [0] = eta,    
   *                  [1] = phi
   *  @param aGridIdx output structure containing grid index
   *                  [0] = 1D  index in primary region
   *                  [1] = 1D  index in outter region
   *                  [2] = eta index in primary region
   *                  [3] = phi index in primary region
   *                  [4] = phi index in outter region
   *  @return  0 on sucess, 1 for outter phi present, <0 for error
   */
  Int_t GetCellIndex( const Float_t* aEtaPhi, Int_t* aGridIdx );

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // -- Cells

  /** Search Grid */
  TClonesArray*  fGrid;                    //! transient

  // -- Grid boundaries in eta and phi - set via setter

  /** Minimum eta */
  Float_t        fEtaMin;                  // see above

  /** Maximum eta */
  Float_t        fEtaMax;                  // see above

  /** Minimum phi */
  Float_t        fPhiMin;                  // see above

  /** Maximum phi */
  Float_t        fPhiMax;                  // see above

  /** Total eta range */
  Float_t        fEtaRange;                // see above

  /** Total phi range 
   *  abs(phimin)+phiMax + 2*coneRadius;
   */
  Float_t        fPhiRange;                // see above

  // -- Grid parameter - set via setter
  
  /** Binning in eta */
  Float_t        fEtaBinning;              // see above

  /** Binning in phi */
  Float_t        fPhiBinning;              // see above

  // -- Number of bins in eta and phi - set via Initialize()

  /** Number of grid bins in eta */
  Int_t          fEtaNGridBins;            // see above

  /** Number of grid bins in phi */
  Int_t          fPhiNGridBins;            // see above

  /** Total number of bins */
  Int_t          fNBins;                   // see above

  // -- Bins around R - set via Initialize()

  /** Number of grid bins in eta in R */
  Int_t          fEtaNRBins;               // see above

  /** Number of grid bins in phi in R */
  Int_t          fPhiNRBins;               // see above

  // -- Cell iterator - set in SetCellIter()

  /** Eta Idx - current */
  Int_t          fEtaIdxCurrent;           // see above

  /** Eta Idx - min */
  Int_t          fEtaIdxMin;               // see above

  /** Eta Idx - max */
  Int_t          fEtaIdxMax;               // see above

  /** Phi Idx - current */
  Int_t          fPhiIdxCurrent;           // see above

  /** Phi Idx - min */
  Int_t          fPhiIdxMin;               // see above

  /** Phi Idx - max */
  Int_t          fPhiIdxMax;               // see above

  // -- Cone radius - set via setter

  /** Cone radius */
  Float_t        fConeRadius;              // see above

  ClassDef(AliHLTJETConeGrid, 1)

};
#endif

