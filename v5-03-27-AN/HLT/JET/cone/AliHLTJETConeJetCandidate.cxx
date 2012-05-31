//-*- Mode: C++ -*-
// $Id: AliHLTJETConeJetCandidate.cxx  $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTJETConeJetCandidate.cxx
    @author Jochen Thaeder
    @date   
    @brief  Jet candidate of the cone finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliHLTJETConeJetCandidate.h"
#include "AliHLTJETConeEtaPhiCell.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETConeJetCandidate)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
  AliHLTJETConeJetCandidate::AliHLTJETConeJetCandidate() :
    fSeedCellIdx(0),
    fSeedEtaIdx(0),
    fSeedPhiIdx(0),
    fSeedEta(0.),
    fSeedPhi(0.),
    fSeedPt(0.),
    fEta(0.),
    fPhi(0.),
    fPt(0.),
    fNTracks(0),
    fUseWholeCell(kTRUE),
    fConeRadius2(0.) {
    // see header file for class documentation
    // or
    // refer to README to build package
    // or
    // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


//##################################################################################
AliHLTJETConeJetCandidate::AliHLTJETConeJetCandidate( const Float_t* aEtaPhi, 
						      const Int_t* aGridIdx, 
						      Float_t coneRadius,
						      Bool_t useWholeCell ) :
  fSeedCellIdx(aGridIdx[kIdxPrimary]),
  fSeedEtaIdx(aGridIdx[kIdxEtaPrimary]),
  fSeedPhiIdx(aGridIdx[kIdxPhiPrimary]),
  fSeedEta(aEtaPhi[kIdxEta]),
  fSeedPhi(aEtaPhi[kIdxPhi]),
  fSeedPt(aEtaPhi[kIdxPt]),
  fEta(0.0),
  fPhi(0.0),
  fPt(0.0),
  fNTracks(0),
  fUseWholeCell(useWholeCell),
  fConeRadius2(coneRadius*coneRadius) {
  // see header file for class documentation
}

//##################################################################################
AliHLTJETConeJetCandidate::~AliHLTJETConeJetCandidate() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Process 
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETConeJetCandidate::AddCell( AliHLTJETConeEtaPhiCell* cell ) {
  // see header file for class documentation
  
  // -- use whole cell
  // -------------------
  if ( fUseWholeCell ) {
    
    fPt  += cell->GetPt();
    fPhi += cell->GetPhi();
    fEta += cell->GetEta();
    fNTracks += cell->GetNTracks();
    
    HLTDebug("Cell : eta: %f - phi: %f - pt: %f - nTracks: %d .", 
	     cell->GetEta(), cell->GetPhi(), cell->GetPt(), cell->GetNTracks() );

  } // if ( fUseWholeCell ) {

  // -- radius compared to every track
  // -----------------------------------
  else {
    
    HLTDebug("Check NTracks %d.", cell->GetNTracks());

    TObjArray* trackList = cell->GetTrackList();

    // -- Loop over all tracks in cell
    for ( Int_t iter = 0; iter < cell->GetNTracks(); iter++ ) {

      // -- MC particles
      if ( cell->GetTrackType() == kTrackMC ) {

	TParticle* particle = reinterpret_cast<TParticle*> ((*trackList)[iter]);
	
	// -- Check if in cone
	if ( ! InCone( particle->Eta(), particle->Phi() ) )
	  continue;
	
	fPt  += particle->Pt();
	fPhi += particle->Phi();
	fEta += particle->Eta();
	++fNTracks;
      
	HLTDebug("Particle : eta: %f - phi: %f - pt: %f .", 
		 particle->Eta(), particle->Phi(), particle->Pt() );
      
      }
      else if ( cell->GetTrackType() == kTrackESD ) {

	AliESDtrack* esdTrack = reinterpret_cast<AliESDtrack*> ((*trackList)[iter]);
	
	// -- Check if in cone
	if ( ! InCone( esdTrack->Eta(), esdTrack->Phi() ) )
	  continue;
	
	fPt  += esdTrack->Pt();
	fPhi += esdTrack->Phi();
	fEta += esdTrack->Eta();
	++fNTracks;
      
	HLTDebug("ESDTrack : eta: %f - phi: %f - pt: %f .", 
		 esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt() );
      }
      else {
	HLTError("Not implemented yet...");
      }

    } // for ( Int_t iter = 0; iter > cell->GetNTracks(); iter++ ) {
  }
  
  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                                Sort of JetCandidates 
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTJETConeJetCandidate::Compare( const TObject* obj) const {
  // see header file for class documentation
  
  if (this == obj ) 
    return 0;
  
  AliHLTJETConeJetCandidate * cand = dynamic_cast<AliHLTJETConeJetCandidate*>(const_cast<TObject*>(obj));
  if (cand) {
    if ( fSeedPt < cand->GetSeedPt() ) 
      return 1;
    else 
      return -1;
  }
  else 
    return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                             Helper - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Float_t AliHLTJETConeJetCandidate::GetDistance2( const Float_t eta1, const Float_t phi1, 
						 const Float_t eta2, const Float_t phi2) {
  // see header file for class documentation

  return ( (eta1-eta2)*(eta1-eta2) ) + ( (phi1-phi2)*(phi1-phi2) );
}

//##################################################################################
Bool_t AliHLTJETConeJetCandidate::InCone( Float_t eta, Float_t phi ) {
  // see header file for class documentation

  if ( GetDistance2(fSeedEta,fSeedPhi,eta,phi) <= fConeRadius2 )
    return kTRUE;
  else
    return kFALSE;
}
