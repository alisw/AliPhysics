//-*- Mode: C++ -*-
// $Id: AliHLTJETConeEtaPhiCell.cxx  $
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

/** @file   AliHLTJETConeEtaPhiCell.cxx
    @author Jochen Thaeder
    @date   
    @brief  Cell in eta-phi space of the cone finder
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliHLTJETConeEtaPhiCell.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETConeEtaPhiCell)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
  AliHLTJETConeEtaPhiCell::AliHLTJETConeEtaPhiCell() :
    fEtaIdx(0),
    fPhiIdx(0),
    fPt(0.),
    fEta(0.),
    fPhi(0.),
    fNTracks(0),
    fTrackList(NULL),
    fTrackType(0) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

//##################################################################################
AliHLTJETConeEtaPhiCell::AliHLTJETConeEtaPhiCell( Int_t etaIdx, Int_t phiIdx, 
						  AliESDtrack* track ) :
    fEtaIdx(etaIdx),
    fPhiIdx(phiIdx),
    fPt(0.),
    fEta(0.),
    fPhi(0.),
    fNTracks(0), 
    fTrackList( new TObjArray(20)), // XXXXXX 20
    fTrackType( kTrackESD ) {  
  // see header file for class documentation

  HLTDebug("New cell for ESD etaIdx %d - phiIdx %d", etaIdx, phiIdx );

  AddTrack( track );
}

//##################################################################################
AliHLTJETConeEtaPhiCell::AliHLTJETConeEtaPhiCell( Int_t etaIdx, Int_t phiIdx, 
						  TParticle* particle ) :
    fEtaIdx(etaIdx),
    fPhiIdx(phiIdx),
    fPt(0.),
    fEta(0.),
    fPhi(0.),
    fNTracks(0), 
    fTrackList( new TObjArray(20)), // XXXXXX 20
    fTrackType( kTrackMC ) {  
  // see header file for class documentation
  
  HLTDebug("New cell for MC etaIdx %d - phiIdx %d", etaIdx, phiIdx );
  
  AddTrack( particle );
}


//##################################################################################
AliHLTJETConeEtaPhiCell::~AliHLTJETConeEtaPhiCell() {
  // see header file for class documentation

  if ( fTrackList )
    delete fTrackList;
  fTrackList = NULL;
}

//##################################################################################
void AliHLTJETConeEtaPhiCell::Clear( Option_t* /*option*/ ) {
  // see header file for class documentation

  if ( fTrackList )
    delete fTrackList;
  fTrackList = NULL;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                                     Process 
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTJETConeEtaPhiCell::AddTrack( AliESDtrack* track  ){
  // see header file for class documentation

  Float_t weight = 1.0; // **  XXXX
  
  fEta   += (weight * track->Eta());
  fPhi   += (weight * track->Phi());
  fPt    += (weight * TMath::Abs(track->Pt()));

  fTrackList->Add( (TObject*) track );
  ++fNTracks;

  return;
}

//##################################################################################
void AliHLTJETConeEtaPhiCell::AddTrack( TParticle* particle ){
  // see header file for class documentation

  Float_t weight = 1.0; // **  XXXX

  fEta   += (weight * particle->Eta());
  fPhi   += (weight * particle->Phi());
  fPt    += (weight * TMath::Abs( particle->Pt()));
  
  fTrackList->Add( (TObject*) particle );
  ++fNTracks;

  return;
}
