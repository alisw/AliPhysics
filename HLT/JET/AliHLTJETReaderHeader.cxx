//-*- Mode: C++ -*-
// $Id: AliHLTJETReaderHeader.cxx  $
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

/** @file   AliHLTJETReaderHeader.cxx
    @author Jochen Thaeder
    @date   
    @brief   ReaderHeader for jet finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#include "AliHLTJETReaderHeader.h"
#include "AliHLTJETTrackCuts.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETReaderHeader)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETReaderHeader::AliHLTJETReaderHeader()
  : 
  AliJetReaderHeader("AliHLTJETReaderHeader"),
  fTrackCuts(NULL),
  fSeedCuts(NULL),
  fGridEtaBinning(0.0),
  fGridPhiBinning(0.0),
  fGridEtaRange(0.0),
  fGridPhiRange(0.0),
  fAlgorithm(AliHLTJETBase::kFFSCSquareCell),
  fConeRadius(0.0),
  fUseMC(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETReaderHeader::~AliHLTJETReaderHeader() {
  // see header file for class documentation

}

/*
 * ---------------------------------------------------------------------------------
 *                                   Initialize
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETReaderHeader::Initialize() {
  // see header file for class documentation

  Int_t iResult = 0;

  // -- Set eta and phi range for grid
  fGridPhiRange = fFiducialPhiMin + 
    fFiducialPhiMax + ( 2.0 * fConeRadius );

  fGridEtaRange = TMath::Abs( fFiducialEtaMin ) + fFiducialEtaMax;

  HLTInfo(" -= ReaderHeader =- ");
  HLTInfo(" Cone radius      %f", fConeRadius );
  HLTInfo(" Grid eta binning %f", fGridEtaBinning );
  HLTInfo(" Grid phi binning %f", fGridPhiBinning );
  HLTInfo(" Grid eta range   %f", fGridEtaRange );
  HLTInfo(" Grid phi range   %f", fGridPhiRange );
  HLTInfo(" Algorithm        %s", AliHLTJETBase::fgkJetAlgorithmType[fAlgorithm] );

  if (fUseMC) { HLTInfo(" Use Kinematics   TRUE"); }
  else { HLTInfo( " Use Kinematics   FALSE"); }

  if ( ! fTrackCuts ) {
    HLTError("No track cuts set in reader header");
    iResult = -EINPROGRESS;
  }
  else {
    fTrackCuts->SetEtaRange( fFiducialEtaMin, fFiducialEtaMax );
    fTrackCuts->SetPhiRange( fFiducialPhiMin, fFiducialPhiMax );
    HLTInfo(" -= TrackCuts =- " );
  }

  // Check for seeds only cone based algoritms
  if ( fAlgorithm >= AliHLTJETBase::kFFSCSquareCell ) {
    if ( !fSeedCuts ) {
      HLTError("No seed cuts set in reader header");
      iResult = -EINPROGRESS;
    }
    else {
      fSeedCuts->SetEtaRange( fFiducialEtaMin+fConeRadius, 
			      fFiducialEtaMax-fConeRadius );
      fSeedCuts->SetPhiRange( fFiducialPhiMin, fFiducialPhiMax );
      HLTInfo(" -= SeedCuts =- " );
    }
  }
  
  return iResult;
}

