/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliUtilityDimuonSource.cxx 47782 2011-02-24 18:37:31Z martinez $ */

//-----------------------------------------------------------------------------
/// \class AliUtilityDimuonSource
/// MC utilities to classify di-muons in MC
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#include "AliUtilityDimuonSource.h"

// ROOT includes
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TPDGCode.h"
#include "TMCProcess.h"
#include "TMath.h"

// STEER includes
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"

// PWGmuon includes
#include "AliAnalysisMuonUtility.h"

/// \cond CLASSIMP
ClassImp(AliUtilityDimuonSource) // Class implementation in ROOT context
/// \endcond

//_________________________________________________________
AliUtilityDimuonSource::AliUtilityDimuonSource() :
TObject()
{
  /// Default constructor
}


//_________________________________________________________
AliUtilityDimuonSource::~AliUtilityDimuonSource()
{
  /// Default destructor
}

//________________________________________________________________________
Int_t AliUtilityDimuonSource::GetParticleType ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  //
  /// Get particle type from matched MC track
  //

  if ( fUtilityMuonAncestor.IsUnidentified(track,mcEvent) ) return kUnidentified;
  if ( fUtilityMuonAncestor.IsHadron(track,mcEvent) ) return kRecoHadron;
  if ( fUtilityMuonAncestor.IsSecondaryMu(track,mcEvent) ) return kSecondaryMu;
  if ( fUtilityMuonAncestor.IsBeautyChainMu(track,mcEvent) ) return kBeautyChainMu;
  if ( fUtilityMuonAncestor.IsCharmChainMu(track,mcEvent) ) return kCharmChainMu;

  return kOtherMu;
}

//________________________________________________________________________
Int_t AliUtilityDimuonSource::GetCommonAncestor ( const AliVParticle* track1, const AliVParticle* track2, const AliMCEvent* mcEvent, UInt_t mask ) const
{
  /// Get common ancestor of muon pairs

  const AliVParticle* mcTracks[2] = {track1, track2};
  Int_t labels[2] = {-10,-10};
  for ( Int_t itrack=0; itrack<2; itrack++ ) {
    if ( ! AliAnalysisMuonUtility::IsMCTrack(mcTracks[itrack]) ) {
      Int_t label = mcTracks[itrack]->GetLabel();
      if ( label < 0 ) return -2; // No MC info
      mcTracks[itrack] = mcEvent->GetTrack(label);
    }
    labels[itrack] = mcTracks[itrack]->GetMother();
  }

  Bool_t stopAtQuark = ( mask & (kStopAtQuark|kStopAtLightQuark) );
  Int_t quarkAbsPdg = (mask&kStopAtLightQuark) ? 3 : 9;

  Int_t foundCommon = -1;
  if ( labels[0] == labels[1] ) foundCommon = labels[0];
  else {
    // labels[0] has to be the older
    if ( labels[0] > labels[1] ) {
      Int_t tmpVal = labels[0];
      labels[0] = labels[1];
      labels[1] = tmpVal;
    }
    while ( labels[0] >= 0 && labels[1] >= 0 ) {
      AliVParticle* part = mcEvent->GetTrack(labels[1]);
      Int_t absPdg = TMath::Abs(part->PdgCode());
      if ( stopAtQuark && (absPdg <= quarkAbsPdg || absPdg == 21) ) break;
      // Do not check quarks or gluons,
      // otherwise you can have problems with, e.g.:
      // g -> c cbar
      // beam_p -> c cbar

      Int_t currMother = part->GetMother();
      if ( currMother == labels[0] ) {
        foundCommon = currMother;
        break;
      }
      if ( currMother < labels[0] ) {
        labels[1] = labels[0];
        labels[0] = currMother;
      }
      else labels[1] = currMother;
    }
  }

  if ( foundCommon >= 0 ) {
    // Reject common ancestor if it is a quark or gluon
//     (the protection on absPdg in the while loop is not enough
//     when the two particles both come from one quark)
    AliVParticle* part = mcEvent->GetTrack(foundCommon);
//    if ( foundCommon == 0 && AliAnalysisMuonUtility::GetStatusCode(part) == 21 ) return -1; // Pythia initial particle
    Int_t absPdgCommon = TMath::Abs(part->PdgCode());
//    if  ( absPdgCommon < 10 || absPdgCommon == 21 ) {
//      if ( stopAtQuark ) return -1;
//    }
//    else {
      // Reject the di-quarks in pythia
      Int_t tmpVal = absPdgCommon/1000;
      if ( tmpVal >= 1 && tmpVal <= 5 ) {
        tmpVal = absPdgCommon%100;
        if ( tmpVal == 1 || tmpVal == 3 ) return -1;
      }
//    }
  }

  return foundCommon;
}

//________________________________________________________________________
TString AliUtilityDimuonSource::GetPairType ( const AliVParticle* track1, const AliVParticle* track2, const AliMCEvent* mcEvent )
{
  // Get pair type
  // CAVEAT: this is slow in a double loop since it recomputes the particle type of each track
  Int_t partType1 = GetParticleType(track1,mcEvent);
  Int_t partType2 = GetParticleType(track2,mcEvent);
  Int_t commonAncestor = GetCommonAncestor(track1,track2,mcEvent);
  return GetPairType(partType1,partType2,commonAncestor,mcEvent);
}

//________________________________________________________________________
TString AliUtilityDimuonSource::GetPairType ( Int_t partType1, Int_t partType2, Int_t commonAncestor, const AliMCEvent* mcEvent ) const
{
  /// Get pair type

  if ( partType1 == kUnidentified || partType2 == kUnidentified ) return "AtLeastOneUnknown";
  if ( partType1 == kSecondaryMu || partType2 == kSecondaryMu ) return "AtLeastOneSecondary";

  if ( partType1 == kRecoHadron || partType2 == kRecoHadron ) return "AtLeastOneHadron";

  if ( commonAncestor < 0 ) {
    // No common ancestor found
    if ( partType1 == kCharmChainMu ) {
      if ( partType2 == kCharmChainMu ) return "DDuncorr";
      if ( partType2 == kBeautyChainMu ) return "BDuncorr";
    }
    else if ( partType1 == kBeautyChainMu ) {
      if ( partType2 == kCharmChainMu ) return "BDuncorr";
      if ( partType2 == kBeautyChainMu ) return "BBuncorr";
    }
    return "Uncorrelated";
  }

  // Particles have a common ancestors
  Int_t ancestorPdg = mcEvent->GetTrack(commonAncestor)->PdgCode();
  Int_t absPdg = TMath::Abs(ancestorPdg);
  Bool_t isQuark = ( absPdg < 10 );
  Bool_t isGluon = ( absPdg == 21 );

  if ( partType1 == kCharmChainMu && partType2 == kCharmChainMu ) {
    if ( isQuark || isGluon ) return "DDpair";
    return "DHsame";
  }

  if ( partType1 == kBeautyChainMu && partType2 == kBeautyChainMu ) {
    if ( isQuark || isGluon ) return "BBpair";
    return "BHsame";
  }

  if ( isQuark ) return "LightQuarkFragmentation";
  if ( isGluon ) return "GluonFragmentation";

  TParticlePDG* partPdg = TDatabasePDG::Instance()->GetParticle(ancestorPdg);
  if ( partPdg ) {
    if ( partPdg->AntiParticle() ) return "LightPartDecay";

    TString partName = partPdg->GetName();
    partName.ReplaceAll("/","");
    return partName;
  }

  return Form("PDG_%i",ancestorPdg);
}
