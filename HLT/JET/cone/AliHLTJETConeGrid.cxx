//-*- Mode: C++ -*-
// $Id: AliHLTJETConeGrid.cxx  $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTJETConeGrid.cxx
    @author Jochen Thaeder
    @date   
    @brief  Eta-Phi grid of the cone finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTJETConeGrid.h"
#include "AliHLTJETConeEtaPhiCell.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETConeGrid)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETConeGrid::AliHLTJETConeGrid()
  : 
  fGrid(NULL),
  fEtaMin(-0.9),
  fEtaMax(0.9),
  fPhiMin(0.0),
  fPhiMax(6.3),
  fEtaRange(1.8),
  fPhiRange(6.3),
  fEtaBinning(0.05),
  fPhiBinning(0.05),
  fEtaNGridBins(-1),
  fPhiNGridBins(-1),
  fNBins(-1),
  fEtaNRBins(-1),
  fPhiNRBins(-1),
  fEtaIdxCurrent(0),
  fEtaIdxMin(0),
  fEtaIdxMax(0),
  fPhiIdxCurrent(0),
  fPhiIdxMin(0),
  fPhiIdxMax(0),
  fConeRadius(0.0) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETConeGrid::~AliHLTJETConeGrid() {
  // see header file for class documentation
 
  if ( fGrid ) {
    fGrid->Clear("C");
    delete fGrid;
  }
  fGrid = NULL;
 
}

/*
 * ---------------------------------------------------------------------------------
 *                                   Setup / Reset
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETConeGrid::Initialize() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Set total N bins in eta and phi
  fEtaNGridBins = TMath::CeilNint( fEtaRange / fEtaBinning );
  fPhiNGridBins = TMath::CeilNint( fPhiRange / fPhiBinning );

  // -- Set total number of bins
  fNBins = fEtaNGridBins * fPhiNGridBins;

  // -- Set cells around R in eta and phi
  fEtaNRBins = TMath::CeilNint( fConeRadius / fEtaBinning );
  fPhiNRBins = TMath::CeilNint( fConeRadius / fPhiBinning );

  HLTInfo(" -= Grid =- ");
  HLTInfo(" NGridBins (%d,%d)", fEtaNGridBins, fPhiNGridBins);
  HLTInfo(" NRBins    (%d,%d)", fEtaNRBins   , fPhiNRBins);
  HLTInfo(" NBins      %d", fNBins );

  fGrid = new TClonesArray("AliHLTJETConeEtaPhiCell", fNBins );
  
  if ( ! fGrid ) {
    HLTError( "Error: Setup search grid with size %d .", fNBins );
    iResult = 1;
  }

  return iResult;
}

// #################################################################################
void AliHLTJETConeGrid::Reset() { 
  // see header file for class documentation

  if ( fGrid )
    fGrid->Clear("C");

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                    Process
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETConeGrid::FillTrack( TParticle* particle, const Float_t* aEtaPhi, Int_t* aGridIdx ) {
  // see header file for class documentation

  Int_t iResult = 0;

  // ---------------------------
  // -- Get Cell Indices
  // ---------------------------
  iResult = GetCellIndex( aEtaPhi, aGridIdx );
  if ( iResult < 0 ) {
    HLTError("Error getting cell index");
    return iResult;
  }

  // ---------------------------
  // -- Fill track in primary region
  // ---------------------------
  
  // -- Create new cell and add track to cell
  if (! fGrid->UncheckedAt( aGridIdx[kIdxPrimary] ) ) {
    new( (*fGrid ) [aGridIdx[kIdxPrimary]] ) AliHLTJETConeEtaPhiCell( aGridIdx[kIdxEtaPrimary], 
								      aGridIdx[kIdxPhiPrimary], 
								      particle );
  }
  // -- Add track to existing cell
  else {
    (reinterpret_cast<AliHLTJETConeEtaPhiCell*> ((*fGrid)[aGridIdx[kIdxPrimary]]))->AddTrack(particle);
  }

  // ---------------------------
  // -- Fill track in outter region
  // ---------------------------
  
  // -- if it has to be filled
  if ( iResult == 1 ) {

    // -- Create new cell and add track to cell
    if (! fGrid->UncheckedAt( aGridIdx[kIdxOutter] ) ) {
      new( (*fGrid) [aGridIdx[kIdxOutter]] ) AliHLTJETConeEtaPhiCell( aGridIdx[kIdxEtaPrimary], 
								      aGridIdx[kIdxPhiOutter], 
								      particle );
    }
    // -- Add track to existing cell
    else {
      (reinterpret_cast<AliHLTJETConeEtaPhiCell*> ((*fGrid)[aGridIdx[kIdxOutter]]))->AddTrack(particle);
    }
  }

  return 0;
}

//##################################################################################
Int_t AliHLTJETConeGrid::FillTrack( AliESDtrack *esdTrack, const Float_t* aEtaPhi, Int_t* aGridIdx ) {
  // see header file for class documentation

  Int_t iResult = 0;
 
  // ---------------------------
  // -- Get Cell Indices
  // ---------------------------

  iResult = GetCellIndex( aEtaPhi, aGridIdx );
  if ( iResult < 0 ) {
    return iResult;
  }
    
  // ---------------------------
  // -- Fill track in primary region
  // ---------------------------
  
  // -- Create new cell and add track to cell
  if (! fGrid->UncheckedAt( aGridIdx[kIdxPrimary] ) ) {
    new( (*fGrid ) [aGridIdx[kIdxPrimary]] ) AliHLTJETConeEtaPhiCell( aGridIdx[kIdxEtaPrimary], 
								      aGridIdx[kIdxPhiPrimary], 
								      esdTrack );
  }
  // -- Add track to existing cell
  else {
    (reinterpret_cast<AliHLTJETConeEtaPhiCell*> ((*fGrid)[aGridIdx[kIdxPrimary]]))->AddTrack(esdTrack);
  }
   
  // ---------------------------
  // -- Fill track in outter region
  // ---------------------------
  
  // -- if it has to be filled
  if ( iResult == 1 ) {
    
    // -- Create new cell and add track to cell
    if (! fGrid->UncheckedAt( aGridIdx[kIdxOutter] ) ) {
      new( (*fGrid) [aGridIdx[kIdxOutter]] ) AliHLTJETConeEtaPhiCell( aGridIdx[kIdxEtaPrimary], 
								      aGridIdx[kIdxPhiOutter], 
								      esdTrack );
    }
    // -- Add track to existing cell
    else {
      (reinterpret_cast<AliHLTJETConeEtaPhiCell*> ((*fGrid)[aGridIdx[kIdxOutter]]))->AddTrack(esdTrack);
    }
  }
  
  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Helper - public
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETConeGrid::NextCell() {
  // see header file for class documentation

  Int_t cellIdx = 0;

  ++fPhiIdxCurrent;
  
  if ( fPhiIdxCurrent > fPhiIdxMax ) {
    fPhiIdxCurrent = fPhiIdxMin;

    ++fEtaIdxCurrent;
    
    if ( fEtaIdxCurrent > fEtaIdxMax )
      cellIdx = -1;
  }

  if ( cellIdx != -1 ) 
    cellIdx = fEtaIdxCurrent + ( fPhiIdxCurrent * fEtaNGridBins );

  if ( cellIdx > fNBins ) {
    HLTError("Idx out of bound (%d,%d) -> %d", fEtaIdxCurrent, fPhiIdxCurrent, cellIdx );
    HLTError("MAX %d,%d - %d", fEtaNGridBins, fPhiNGridBins, fNBins );

    cellIdx = -1;
  }

  return cellIdx;
}

// #################################################################################
void AliHLTJETConeGrid::SetCellIter( const Int_t etaIdx, const Int_t phiIdx ) {
  // see header file for class documentation

  fEtaIdxMax = etaIdx + fEtaNRBins;
  fEtaIdxMin = etaIdx - fEtaNRBins;
  fEtaIdxCurrent = fEtaIdxMin;

  fPhiIdxMax = phiIdx + fPhiNRBins;
  fPhiIdxMin = phiIdx - fPhiNRBins;
  fPhiIdxCurrent = fPhiIdxMin - 1;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Helper - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETConeGrid::GetCellIndex( const Float_t* aEtaPhi, Int_t* aGridIdx ) {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Prime is relative to (0,0) in grid which is (-fEtaMax,-fConeRadius)
  Float_t etaPrime = fEtaMax + aEtaPhi[kIdxEta];
  Float_t phiPrime = aEtaPhi[kIdxPhi] + fConeRadius;

  HLTDebug("eta    : %f - phi :%f", aEtaPhi[kIdxEta], aEtaPhi[kIdxPhi] );
  HLTDebug("eta\'   : %f - phi\':%f", etaPrime, phiPrime );

  // -- Index in 2D (idxEta,idxPhi) 
  aGridIdx[kIdxEtaPrimary] = TMath::FloorNint( etaPrime / fEtaBinning );
  aGridIdx[kIdxPhiPrimary] = TMath::FloorNint( phiPrime / fPhiBinning );
 
  // -- Index in 1D 
  aGridIdx[kIdxPrimary] = aGridIdx[kIdxEtaPrimary] + ( aGridIdx[kIdxPhiPrimary] * fEtaNGridBins );
  
  // -- Boundery Check 2D
  if ( (aGridIdx[kIdxEtaPrimary] < 0) || 
       (aGridIdx[kIdxPhiPrimary] < 0) || 
       (aGridIdx[kIdxEtaPrimary] >= fEtaNGridBins) || 
       (aGridIdx[kIdxPhiPrimary] >= fPhiNGridBins) ) {
    HLTError ( "Index out of range: idxEta %d - max %d, idxPhi %d - max %d", 
	       aGridIdx[kIdxEtaPrimary], fEtaNGridBins, 
	       aGridIdx[kIdxPhiPrimary], fPhiNGridBins);
    iResult = -1;
  }

  // -- Boundery Check 1D
  if ( (aGridIdx[kIdxPrimary] < 0) ||
       (aGridIdx[kIdxPrimary] >= fNBins) ) {
    HLTError( "Index out of range: 1D idx %d - max %d", 
	      aGridIdx[kIdxPrimary], fNBins );
    iResult = -2;
  }
  
  HLTDebug( "idxEta %d - max %d, idxPhi %d - max %d", 
	    aGridIdx[kIdxEtaPrimary], fEtaNGridBins, 
	    aGridIdx[kIdxPhiPrimary], fPhiNGridBins);
  HLTDebug( "1D idx %d - max %d", aGridIdx[kIdxPrimary], fNBins );
  
  if ( iResult ) 
    return iResult;

  // -- check if to map from border region to outter region
  Float_t phiOutterPrime = 0.0;
  
  // -- upper border
  if ( aEtaPhi[kIdxPhi] > ( fPhiMax - fConeRadius ) ) {
    phiOutterPrime = phiPrime - fPhiMax;
    iResult = 1;  

    HLTDebug("eta    : %f - phiOutter :%f .", aEtaPhi[kIdxEta], aEtaPhi[kIdxPhi] - fPhiMax );
    HLTDebug("eta\'   : %f - phiOutter\':%f .", etaPrime, phiOutterPrime );
  }
  // -- lower border
  else if ( aEtaPhi[kIdxPhi] < fConeRadius ) {
    phiOutterPrime = phiPrime + fPhiMax;
    iResult = 1;

    HLTDebug("eta    : %f - phiOutter :%f .", aEtaPhi[kIdxEta], aEtaPhi[kIdxPhi] + fPhiMax );
    HLTDebug("eta\'   : %f - phiOutter\':%f .", etaPrime, phiOutterPrime );
  }
  
  // -- if outter phi present
  if ( iResult == 1 ) {

    // -- Index in 2D (idxEta,idxPhiOutter) 
    aGridIdx[kIdxPhiOutter] = TMath::FloorNint( phiOutterPrime / fPhiBinning );

    // -- Index in 1D
    aGridIdx[kIdxOutter] = aGridIdx[kIdxEtaPrimary] + ( aGridIdx[kIdxPhiOutter] * fEtaNGridBins );
  
    // -- Boundery Check 2D
    if ( aGridIdx[kIdxPhiOutter] < 0 || aGridIdx[kIdxPhiOutter] >= fPhiNGridBins ) {
      HLTError( "Index out of range ( Outter Phi ): idxPhiOutter %d - max %d", 
		aGridIdx[kIdxPhiOutter], fPhiNGridBins);
      iResult = -3;
    }
    
    // -- Boundery Check 1D
    if ( aGridIdx[kIdxOutter] >= fNBins ) {
      HLTError( "Index out of range ( Outter Phi ): 1D idx %d - max %d", 
		aGridIdx[kIdxOutter], fNBins );
      iResult = -4;
    }

    HLTDebug( "idxEta %d - max %d, idxPhiOutter %d - max %d", 
             aGridIdx[kIdxEtaPrimary], fEtaNGridBins, aGridIdx[kIdxPhiOutter], fPhiNGridBins);
    HLTDebug( "1D idx %d - max %d", aGridIdx[kIdxOutter], fNBins );
  }

  return iResult;
}
