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

/* $Id: AliAnalysisMuonUtility.cxx 47782 2011-02-24 18:37:31Z martinez $ */

//-----------------------------------------------------------------------------
/// \class AliAnalysisMuonUtility
/// Static tilities for muon analysis
/// The class allows to treat AODs and ESDs
/// as well as MC AODs and MC in a transparent way
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#include "AliAnalysisMuonUtility.h"

// ROOT includes
#include "TAxis.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"

// STEER includes
#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliVVertex.h"
#include "AliLog.h"

// CORRFW includes
#include "AliCFGridSparse.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisMuonUtility) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsAODEvent ( const AliVEvent* event )
{
  /// Check if event is from ESD or AOD
  return ( event->IsA() == AliAODEvent::Class() );
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsAODTrack ( const AliVParticle* track )
{
  /// Check if track is from ESD or AOD
  return ( track->IsA() == AliAODTrack::Class() );
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsMuonTrack ( const AliVParticle* track )
{
  /// Check if track has muon tracker info
  return ( IsAODTrack(track) ) ? static_cast<const AliAODTrack*>(track)->IsMuonTrack() : static_cast<const AliESDMuonTrack*>(track)->ContainTrackerData();
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsMuonGhost ( const AliVParticle* track )
{
  /// Check if track has trigger info
  return ( IsAODTrack(track) ) ? kFALSE : ( ! static_cast<const AliESDMuonTrack*>(track)->ContainTrackerData() );
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetRabs ( const AliVParticle* track )
{
  /// Get Rabs
  return ( IsAODTrack(track) ) ? static_cast<const AliAODTrack*>(track)->GetRAtAbsorberEnd() : static_cast<const AliESDMuonTrack*>(track)->GetRAtAbsorberEnd();
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetThetaAbsDeg ( const AliVParticle* track )
{
  /// Get Theta at absorber end (in deg)
  return TMath::ATan( GetRabs(track) / 505. ) * TMath::RadToDeg();
}

//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetMatchTrigger ( const AliVParticle* track )
{
  /// Get match trigger
  return IsAODTrack(track) ? static_cast<const AliAODTrack*>(track)->GetMatchTrigger() : static_cast<const AliESDMuonTrack*>(track)->GetMatchTrigger();
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetChi2perNDFtracker ( const AliVParticle* track )
{
  /// Get Theta at absorber end (in deg)
  return IsAODTrack(track) ? static_cast<const AliAODTrack*>(track)->Chi2perNDF() : static_cast<const AliESDMuonTrack*>(track)->GetNormalizedChi2();
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetXatVertex ( const AliVParticle* track )
{
  /// Get X at vertex
  Double_t coor = 0.;
  if ( IsAODTrack(track) ) {
    Double_t vtxPos[3];
    static_cast<const AliAODTrack*>(track)->GetXYZ(vtxPos);
    coor = vtxPos[0];
  }
  else coor = static_cast<const AliESDMuonTrack*>(track)->GetNonBendingCoor();
  
  return coor;
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetYatVertex ( const AliVParticle* track )
{
  /// Get Y at vertex
  Double_t coor = 0.;
  if ( IsAODTrack(track) ) {
    Double_t vtxPos[3];
    static_cast<const AliAODTrack*>(track)->GetXYZ(vtxPos);
    coor = vtxPos[1];
  }
  else coor = static_cast<const AliESDMuonTrack*>(track)->GetBendingCoor();
  
  return coor;
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetZatVertex ( const AliVParticle* track )
{
  /// Get Z at vertex
  Double_t coor = 0.;
  if ( IsAODTrack(track) ) {
    Double_t vtxPos[3];
    static_cast<const AliAODTrack*>(track)->GetXYZ(vtxPos);
    coor = vtxPos[2];
  }
  else coor = static_cast<const AliESDMuonTrack*>(track)->GetZ();
  
  return coor;
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetXatDCA ( const AliVParticle* track )
{
  /// Get X at DCA
  return IsAODTrack(track) ? static_cast<const AliAODTrack*>(track)->XAtDCA() : static_cast<const AliESDMuonTrack*>(track)->GetNonBendingCoorAtDCA();
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetYatDCA ( const AliVParticle* track )
{
  /// Get Y at DCA
  return IsAODTrack(track) ? static_cast<const AliAODTrack*>(track)->YAtDCA() : static_cast<const AliESDMuonTrack*>(track)->GetBendingCoorAtDCA();
}

//________________________________________________________________________
UInt_t AliAnalysisMuonUtility::GetMUONTrigHitsMapTrk ( const AliVParticle* track )
{
  /// Get hit pattern in trigger chambers from tracker track extrapolation
  return ( IsAODTrack(track) ) ? const_cast<AliAODTrack*>(static_cast<const AliAODTrack*>(track))->GetMUONTrigHitsMapTrk() : static_cast<const AliESDMuonTrack*>(track)->GetHitsPatternInTrigChTrk();
}

//________________________________________________________________________
UInt_t AliAnalysisMuonUtility::GetMUONTrigHitsMapTrg ( const AliVParticle* track )
{
  /// Get hit pattern in trigger chambers from tracker track extrapolation
  return ( IsAODTrack(track) ) ? const_cast<AliAODTrack*>(static_cast<const AliAODTrack*>(track))->GetMUONTrigHitsMapTrg() : static_cast<const AliESDMuonTrack*>(track)->GetHitsPatternInTrigCh();
}

//________________________________________________________________________
TString AliAnalysisMuonUtility::GetFiredTriggerClasses ( const AliVEvent* event )
{
  /// Check if track is from ESD or AOD
  return ( IsAODEvent(event) ) ? static_cast<const AliAODEvent*>(event)->GetFiredTriggerClasses() : static_cast<const AliESDEvent*>(event)->GetFiredTriggerClasses();
}

//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetNTracks ( const AliVEvent* event )
{
  //
  /// Return the number of tracks in event
  //
  return ( IsAODEvent(event) ) ? static_cast<const AliAODEvent*>(event)->GetNTracks() : static_cast<const AliESDEvent*>(event)->GetNumberOfMuonTracks();
}


//________________________________________________________________________
AliVParticle* AliAnalysisMuonUtility::GetTrack ( Int_t itrack, const AliVEvent* event )
{
  //
  /// Get the current track
  //
  AliVParticle* track = 0x0;
  if ( IsAODEvent(event) ) track = static_cast<const AliAODEvent*>(event)->GetTrack(itrack);
  else {
    AliESDEvent* esdEvent = const_cast<AliESDEvent*>(static_cast<const AliESDEvent*> (event));
    track = esdEvent->GetMuonTrack(itrack);
  }
  return track;
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::MuonMass2()
{
  /// A usefull constant
  return 1.11636129640000012e-02;
}

//________________________________________________________________________
TLorentzVector AliAnalysisMuonUtility::GetTrackPair ( const AliVParticle* track1, const AliVParticle* track2 )
{
  //
  /// Get track pair
  //
  const AliVParticle* tracks[2] = {track1, track2};
  
  TLorentzVector vec[2];
  for ( Int_t itrack=0; itrack<2; ++itrack ) {
    Double_t trackP = tracks[itrack]->P();
    Double_t energy = TMath::Sqrt(trackP*trackP + MuonMass2());
    vec[itrack].SetPxPyPzE(tracks[itrack]->Px(), tracks[itrack]->Py(), tracks[itrack]->Pz(), energy);
  }
  
  TLorentzVector vecPair = vec[0] + vec[1];
  return vecPair;
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsMCEvent ( const AliVEvent* event, const AliMCEvent* mcEvent )
{
  //
  /// Contains MC info
  //
  return ( mcEvent || ( IsAODEvent(event) && static_cast<const AliAODEvent*>(event)->FindListObject(AliAODMCParticle::StdBranchName()) ) );
}


//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsAODMCTrack( const AliVParticle* mcParticle )
{
  /// Check if track is from AOD MC
  return ( mcParticle->IsA() == AliAODMCParticle::Class() );
}


//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetNMCTracks ( const AliVEvent* event, const AliMCEvent* mcEvent )
{
  //
  /// Return the number of MC tracks in event
  //
  Int_t nMCtracks = 0;
  if ( mcEvent ) nMCtracks = mcEvent->GetNumberOfTracks();
  else if ( IsAODEvent(event) ) {
    TClonesArray* mcArray = (TClonesArray*)static_cast<const AliAODEvent*>(event)->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if ( mcArray ) nMCtracks = mcArray->GetEntries();
  }
  return nMCtracks;
}

//________________________________________________________________________
AliVParticle* AliAnalysisMuonUtility::GetMCTrack ( Int_t trackLabel, const AliVEvent* event, const AliMCEvent* mcEvent )
{
  //
  /// MC information can be provided by the MC input handler
  /// (mostly when analyising ESDs) or can be found inside AODs
  /// This method returns the correct one
  //
  AliVParticle* mcTrack = 0x0;
  if ( mcEvent ) mcTrack = mcEvent->GetTrack(trackLabel);
  else if ( IsAODEvent(event) ) {
    TClonesArray* mcArray = (TClonesArray*)static_cast<const AliAODEvent*>(event)->FindListObject(AliAODMCParticle::StdBranchName());
    if ( mcArray ) mcTrack =  (AliVParticle*)mcArray->At(trackLabel);
  }
  if ( ! mcTrack ) AliWarningClass(Form("No track with label %i!", trackLabel));
  return mcTrack;
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetMCVertexZ ( const AliVEvent* event, const AliMCEvent* mcEvent )
{
  /// Get MC vertex Z
  Double_t vz = 0.;
  if ( mcEvent ) vz = mcEvent->GetPrimaryVertex()->GetZ();
  else if ( IsAODEvent(event) ) {
    AliAODMCHeader* aodMCHeader = static_cast<AliAODMCHeader*> (static_cast<const AliAODEvent*>(event)->FindListObject(AliAODMCHeader::StdBranchName()));
    vz = aodMCHeader->GetVtxZ();
  }
  else AliErrorClass("No MC event found");
  return vz;
}

//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetMotherIndex ( const AliVParticle* mcParticle )
{
  //
  /// Return the mother index
  //
  return ( IsAODMCTrack(mcParticle) ) ? static_cast<const AliAODMCParticle*>(mcParticle)->GetMother() : static_cast<const AliMCParticle*>(mcParticle)->GetMother();
}

//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetDaughterIndex ( const AliVParticle* mcParticle, Int_t idaughter )
{
  //
  /// Return the daughter index
  /// idaughter can be:
  /// 0 -> first daughter
  /// 1 -> last daughter
  //
  if ( idaughter < 0 || idaughter > 1 ) {
    AliErrorClass(Form("Requested daughter %i Daughter index can be either 0 (first) or 1 (last)", idaughter));
    return -1;
  }
  
  if ( IsAODMCTrack(mcParticle) ) return static_cast<const AliAODMCParticle*>(mcParticle)->GetDaughter(idaughter);
  
  if ( idaughter == 0 ) return static_cast<const AliMCParticle*>(mcParticle)->GetFirstDaughter();
  else return static_cast<const AliMCParticle*>(mcParticle)->GetLastDaughter();
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsPrimary ( const AliVParticle* mcParticle, const AliMCEvent* mcEvent )
{
  /// Check if the particle is primary
  
  Bool_t isPrimary = kFALSE;
  if ( mcEvent ) {
    // First get the index of the current particle in the stack.
    // For this: get the mother, and get its daughter.
    // Since the mother can have many daughters, you can come up to a "sister"
    // of the particle. Nevertheless, if it is primary, then also the particle itself should be.
    Int_t imother = static_cast<const AliMCParticle*> (mcParticle)->GetMother();
    if ( imother < 0 ) isPrimary = kTRUE;
    else if ( static_cast<const AliMCParticle*>(mcEvent->GetTrack(imother))->GetFirstDaughter() < const_cast<AliMCEvent*>(mcEvent)->GetNumberOfPrimaries() ) isPrimary = kTRUE;
  }
  else isPrimary = static_cast<const AliAODMCParticle*>(mcParticle)->IsPrimary();
  return isPrimary;
}


//________________________________________________________________________
AliVVertex* AliAnalysisMuonUtility::GetVertexSPD ( const AliVEvent* event )
{
  //
  /// Get vertex SPD
  //
  
  AliVVertex* primaryVertex = ( IsAODEvent(event) ) ? (AliVVertex*)static_cast<const AliAODEvent*>(event)->GetPrimaryVertexSPD() : (AliVVertex*)static_cast<const AliESDEvent*>(event)->GetPrimaryVertexSPD();
  return primaryVertex;
}


//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetPassNumber ( const AliInputEventHandler* eventHandler )
{
  /// Get pass number from event

  // At present there is no straightforward way to get the pass number.
  // The pass number is usually written in the path to the ESDs/AODs
  // but this won't work for:
  // - MC data (no pass info available)
  // - local ESDs/AODs
  
  TString filePath = "";
  const AliVEvent* event = eventHandler->GetEvent();  
  if ( IsAODEvent(event) ) {
    // In AODs, the header contains the path to the input ESD event
    // However, sometimes there is not the FULL path of the ESDs.
    // In that case we cannot extract the pass number.
    // To increase the possibility of getting the pass number,
    // try first to find the info in the AOD header
    // (which is a priori safer because it works even on local copies of AODs)
    // and if it does not work, directly check the path to the AOD
    filePath = static_cast<const AliAODEvent*> (event)->GetHeader()->GetESDFileName();
    Int_t passNumber = GetPassNumber(filePath.Data());
    if ( passNumber < 0 ) AliWarningClass("Check again with the AOD path");
    else return passNumber;
  }
  
  filePath = eventHandler->GetTree()->GetCurrentFile()->GetName();
  return GetPassNumber(filePath.Data());
}


//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetPassNumber ( const char* str )
{
  //
  /// Get pass number from string
  //
  
  TString currName(str);
  currName.ToUpper();
  Int_t idx = currName.Index("PASS");
  if ( idx >= 0 ) {
    // Remove all of the words before PASS (and remove the word PASS itself)
    currName.Remove(0,idx+4);
    // Cut the word from the end untill only the digits are left
    // (In this way we can extract the pass number even if it has more than one digit)
    while ( ! currName.IsDigit() && currName.Length() > 0 ) {
      currName.Remove(currName.Length()-1);
    }
    
    if ( currName.Length() > 0 ) return currName.Atoi();
  }
  
  AliWarningClass(Form("Cannot find pass number in: %s", str));
  return -1;
}


//_______________________________________________________________________
Bool_t AliAnalysisMuonUtility::SetSparseRange(AliCFGridSparse* gridSparse,
                                        Int_t ivar, TString labelName,
                                        Double_t varMin, Double_t varMax,
                                        TString option)
{
  //
  /// Set range in a smart way.
  /// Allows to select a bin from the label.
  /// Check the bin limits.
  //
  
  option.ToUpper();
  Int_t minVarBin = -1, maxVarBin = -1;
  TAxis* axis = gridSparse->GetAxis(ivar);
  
  if ( ! axis ) {
    printf("Warning: Axis %i not found in %s", ivar, gridSparse->GetName());
    return kFALSE;
  }
  
  if ( ! labelName.IsNull() ) {
    minVarBin = axis->FindBin(labelName.Data());
    maxVarBin = minVarBin;
    if ( minVarBin < 1 ) {
      printf("Warning: %s: label %s not found. Nothing done", gridSparse->GetName(), labelName.Data());
      return kFALSE;
    }
  }
  else if ( option.Contains( "USEBIN" ) ) {
    minVarBin = (Int_t)varMin;
    maxVarBin = (Int_t)varMax;
  }
  else {
    minVarBin = axis->FindBin(varMin);
    maxVarBin = axis->FindBin(varMax);
  }
  
  if ( axis->GetFirst() == minVarBin && axis->GetLast() == maxVarBin ) return kFALSE;
  
  gridSparse->SetRangeUser(ivar, axis->GetBinCenter(minVarBin), axis->GetBinCenter(maxVarBin));
  
  return kTRUE;
}
