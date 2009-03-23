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

/** @file   AliHLTJETReader.cxx
    @author Jochen Thaeder
    @date   
    @brief   Reader for jet finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTJETReader.h"

#include "AliHLTJETTrackCuts.h"

#include "TLorentzVector.h"
#include "TParticle.h"
#include "TParticlePDG.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETReader)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
// #################################################################################
AliHLTJETReader::AliHLTJETReader()
  : 
  AliJetReader(),
  fESD(NULL), 
  fMC(NULL),
  fAOD(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETReader::~AliHLTJETReader() {
  // see header file for class documentation

}

/*
 * ---------------------------------------------------------------------------------
 *                               Reader functionality
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArray() {
  // see header file for class documentation

  Bool_t bResult = kFALSE;

  if ( fESD )
    bResult = FillMomentumArrayESD();
  else if ( fMC )
    bResult = FillMomentumArrayMC();
  else if ( fAOD )
    bResult = FillMomentumArrayAOD();
  
  return bResult;
}

// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArrayMC() {
  // see header file for class documentation

  if ( ! fMC ) {
    HLTError( "No MC Event!" );
    return kFALSE;
  }

  Int_t nTracks = 0;

  // -- Clear Array
  ClearArray(); 

  // -- Get cuts
  AliHLTJETTrackCuts *trackCuts = reinterpret_cast<AliHLTJETTrackCuts*>(GetReaderHeader()->GetAnalysisCuts()) ;
						
  //---------------------------------------------------------------------------
  //  XXX fAliHeader = fMCEvent->Header();
  //---------------------------------------------------------------------------

  TLorentzVector p4;

  TParticle* particle = NULL;
  
  // -- Loop over particles
  // ------------------------
  while ( (particle = fMC->NextParticle() ) ) {

    // -- Apply cuts 
    if ( ! trackCuts->IsSelected(particle) )
      continue;

    // -- Fill vector
    p4 = TLorentzVector( particle->Px(), particle->Py(), particle->Pz(), particle->Energy() );
        
    // -- new TClonesArray entry
    new ( (*fMomentumArray)[nTracks] ) TLorentzVector(p4);
  
    nTracks++;    
    
  } // while ( (particle = fMC->NextParticle() ) ) {

  HLTInfo(" Number of selected tracks %d \n", nTracks);

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArrayESD() {
  // see header file for class documentation

  return kTRUE;
}

// #################################################################################
Bool_t AliHLTJETReader::FillMomentumArrayAOD() {
  // see header file for class documentation

  return kFALSE;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Setter
 * ---------------------------------------------------------------------------------
 */


// #################################################################################
void AliHLTJETReader::SetInputEvent(TObject* esd, TObject* aod, TObject* mc) {
  // see header file for class documentation

  if ( esd )
    fESD = dynamic_cast<AliESDEvent*> (esd);
  else if ( aod )
    fAOD = dynamic_cast<AliAODEvent*> (aod);
  else if ( mc )
    fMC = dynamic_cast<AliHLTMCEvent*> (mc);
  
  return;
}
