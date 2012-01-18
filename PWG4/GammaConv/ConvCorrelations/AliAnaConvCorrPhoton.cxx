/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Author: Svein Lindal <slindal@fys.uio.no>                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliAnaConvCorrPhoton.cxx
/// @author Svein Lindal
/// @brief  Class used to fill calculate correlation between photons - tracks
 
#include "AliAnaConvCorrPhoton.h"
#include "AliAODTrack.h"
#include "TClonesArray.h"
#include "AliAODConversionPhoton.h"

#include <iostream>
// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;
ClassImp(AliAnaConvCorrPhoton)

//________________________________________________________________________
AliAnaConvCorrPhoton::AliAnaConvCorrPhoton() :
AliAnaConvCorrBase("photon_hadron_corr", "Photon corr"), 
  fSkipDecayParticles(kFALSE),
  fDecayOnly(kFALSE)
{
  //consctructor
}
//________________________________________________________________________
AliAnaConvCorrPhoton::AliAnaConvCorrPhoton(TString name, TString title = "Photon corr") :
AliAnaConvCorrBase(name, title),
fSkipDecayParticles(kFALSE),
fDecayOnly(kFALSE)
{
  //consctructor
}

//________________________________________________________________________________
AliAnaConvCorrPhoton::~AliAnaConvCorrPhoton() {
  //destructor
}

// //________________________________________________________________________________
// void AliAnaConvCorrPhoton::Process(const TClonesArray * photons, const TClonesArray * tracks, Bool_t isolated = kFALSE) {
//   //Process list of photons and correlate w tracks
//   for(Int_t ig = 0; ig < photons->GetEntriesFast(); ig++) {

// 	AliAODConversionParticle * photon = static_cast<AliAODConversionParticle*>(photons->UncheckedAt(ig));

// 	Int_t tIDs[4] = {-1, -1, -1, -1};
// 	tIDs[0] =  photon->GetLabel(0);
// 	tIDs[1] =  photon->GetLabel(1);
// 	CorrelateWithTracks(photon, tracks, tIDs, isolated);
		
//   }
// }
