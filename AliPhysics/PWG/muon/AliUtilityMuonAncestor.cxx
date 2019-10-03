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

/* $Id: AliUtilityMuonAncestor.cxx 47782 2011-02-24 18:37:31Z martinez $ */

//-----------------------------------------------------------------------------
/// \class AliUtilityMuonAncestor
/// Static utilities to get the muon ancestor in MC
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#include "AliUtilityMuonAncestor.h"

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
ClassImp(AliUtilityMuonAncestor) // Class implementation in ROOT context
/// \endcond

//_________________________________________________________
AliUtilityMuonAncestor::AliUtilityMuonAncestor() :
TObject(),
fPx(0.),
fPy(0.),
fPz(0.),
fMask(0),
fAncestor(-999)
{
  /// Default constructor
}


//_________________________________________________________
AliUtilityMuonAncestor::~AliUtilityMuonAncestor()
{
  /// Default destructor
}

//_________________________________________________________
AliUtilityMuonAncestor::AliUtilityMuonAncestor(const AliUtilityMuonAncestor& obj) :
TObject(obj),
fPx(obj.fPx),
fPy(obj.fPy),
fPz(obj.fPz),
fMask(obj.fMask),
fAncestor(obj.fAncestor)
{
  /// Copy constructor
}

//_________________________________________________________
AliUtilityMuonAncestor& AliUtilityMuonAncestor::operator=(const AliUtilityMuonAncestor& obj)
{
  /// Copy operator
  if ( this != &obj ) {
    TObject::operator=(obj);
    fPx = obj.fPx;
    fPy = obj.fPy;
    fPz = obj.fPz;
    fMask = obj.fMask;
    fAncestor = obj.fAncestor;
  }
  return *this;
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::BuildAncestor ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Build ancestor
  
  // If track is the same as the one in memory, do not re-build the track ancestor
  if ( track->Px() == fPx && track->Py() == fPy && track->Pz() == fPz ) return kTRUE;
  fPx = track->Px();
  fPy = track->Py();
  fPz = track->Pz();
  fMask = 0;
  fAncestor = -999;
  if ( ! mcEvent ) return kFALSE;
  
  const AliVParticle* mcParticle = 0x0;
  
  if ( AliAnalysisMuonUtility::IsMCTrack(track) ) mcParticle = track;
  else {
    Int_t trackLabel = track->GetLabel();
    if ( trackLabel < 0 ) return kFALSE;
    mcParticle = mcEvent->GetTrack(trackLabel);
  }
  
  // Track is MC (or matches a MC track)
  SETBIT(fMask,kIsID);
  
  Int_t recoPdg = mcParticle->PdgCode();

  if ( TMath::Abs(recoPdg) == 13 ) SETBIT(fMask,kIsMuon);
  else return kTRUE; // Track is not a muon. Do not check ancestors
  
  Int_t imother = mcParticle->GetMother();

  while ( imother >= 0 ) {
    const AliVParticle* part = mcEvent->GetTrack(imother);
    
    Int_t absPdg = TMath::Abs(part->PdgCode());
    
    // This is a quark
    if ( absPdg < 10 ) return kTRUE;
    
    fAncestor = imother;
    Bool_t isPrimary = AliAnalysisMuonUtility::IsPrimary(part, mcEvent);
    
    if ( isPrimary ) {
      if ( absPdg == 15 ) SETBIT(fMask,kHasTauParent);
      else {
        Int_t mpdg = absPdg%100000;
        if ( mpdg >= 100 && mpdg < 10000 ) {
          Int_t flv  = Int_t ( mpdg / TMath::Power(10, Int_t(TMath::Log10(mpdg) )));
          if ( flv < 4 ) SETBIT(fMask,kHasLightParent);
          else if ( flv >= 6 ) continue;
          else {
            TParticlePDG* partPdg = TDatabasePDG::Instance()->GetParticle(part->PdgCode());
            if ( partPdg && ! partPdg->AntiParticle() ) SETBIT(fMask,kHasQuarkoniumParent);
            else if ( flv == 4 ) SETBIT(fMask,kHasCharmParent);
            else SETBIT(fMask,kHasBeautyParent);
          }
        } // absPdg within 100 and 10000
      }
    } // is primary
    else {
      UInt_t mcProcess = AliAnalysisMuonUtility::GetMCProcess(part);
      if ( mcProcess == kPHadronic ||
           ( mcProcess == 0 && part->Zv() < -90. ) ) {
        // The MC process is not well computed in the AODs of old MC productions
        // In this case, declare the particle as "secondary" if it is produced inside the front absorber
        SETBIT(fMask,kIsSecondary);
//        return kTRUE;
      }
    } // is secondary
    
    imother = part->GetMother();
    
  } // loop on mothers
  return kTRUE;
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::CheckAncestor ( const AliVParticle* track, const AliMCEvent* mcEvent, Int_t ancestorPdg, Bool_t matchAbsPdg )
{
  /// Check ancestor
  Int_t pdgCode = GetAncestorPdg(track, mcEvent);
  if ( matchAbsPdg ) {
    pdgCode = TMath::Abs(pdgCode);
    ancestorPdg = TMath::Abs(ancestorPdg);
  }
  return ( pdgCode == ancestorPdg );
}

//_________________________________________________________
Int_t AliUtilityMuonAncestor::GetAncestor ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Return ancestor (compute it if necessary)
  BuildAncestor(track,mcEvent);
  return fAncestor;
}

//_________________________________________________________
Int_t AliUtilityMuonAncestor::GetAncestorPdg ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Return ancestor Pdg
  if ( BuildAncestor(track,mcEvent) ) {
    if ( fAncestor >= 0 ) return mcEvent->GetTrack(fAncestor)->PdgCode();
  }
  return 0;
}

//_________________________________________________________
Long64_t AliUtilityMuonAncestor::GetMask ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Return mask
  BuildAncestor(track,mcEvent);
  return fMask;
}


//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsBeautyMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from beauty decays
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) && TESTBIT(mask,kHasBeautyParent) & ! TESTBIT(mask,kHasLightParent) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsBeautyChainMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from beauty decays
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) && TESTBIT(mask,kHasBeautyParent) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsBJpsiMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon B->J/psi decays
  if ( IsBeautyMu(track,mcEvent) ) {
    Int_t imother = track->GetMother();
    if ( imother >= 0 ) return ( mcEvent->GetTrack(imother)->PdgCode() == 443 );
  }
  return kFALSE;
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsCharmMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from charm decays
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) && TESTBIT(mask,kHasCharmParent) && ! TESTBIT(mask,kHasBeautyParent) && ! TESTBIT(mask,kHasLightParent) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsCharmChainMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from charm decays
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) && TESTBIT(mask,kHasCharmParent) && ! TESTBIT(mask,kHasBeautyParent) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsDecayMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from light hadron decays
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) && TESTBIT(mask,kHasLightParent) && ! TESTBIT(mask,kIsSecondary) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsHadron ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Reconstructed hadron
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && ! TESTBIT(mask,kIsMuon) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsMuon ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Track is muon
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsQuarkoniumMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Mu from quarkonium decay
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) && TESTBIT(mask,kHasQuarkoniumParent) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsSecondaryMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from secondary decays in absorber
  Long64_t mask = GetMask(track,mcEvent);
  return ( TESTBIT(mask,kIsID) && TESTBIT(mask,kIsMuon) && TESTBIT(mask,kIsSecondary) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsUnidentified ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Unidentified muon
  Long64_t mask = GetMask(track,mcEvent);
  return ( ! TESTBIT(mask,kIsID) );
}

//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsWBosonMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from W boson decays
  Long64_t mask = GetMask(track,mcEvent);
  return ( IsMuon(track,mcEvent) && ! TESTBIT(mask,kHasTauParent) && CheckAncestor(track,mcEvent,24) );
}
          
//_________________________________________________________
Bool_t AliUtilityMuonAncestor::IsZBosonMu ( const AliVParticle* track, const AliMCEvent* mcEvent )
{
  /// Muon from Z boson decays
  Long64_t mask = GetMask(track,mcEvent);
  return ( IsMuon(track,mcEvent) && ! TESTBIT(mask,kHasTauParent) && CheckAncestor(track,mcEvent,kZ0) );
}
