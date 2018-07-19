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
#include "TRegexp.h"

// STEER includes
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODTZERO.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDTZERO.h"
#include "AliVEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliVVertex.h"
#include "AliStack.h"

// CORRFW includes
#include "AliCFGridSparse.h"
#include "assert.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisMuonUtility) // Class implementation in ROOT context
/// \endcond


Bool_t AliAnalysisMuonUtility::fUseSmearedTracks = kFALSE;
const char* AliAnalysisMuonUtility::fSmearedTrackListName = "smearedMuonTracks";

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
  /// Get the normalized chi2 of the tracker track
  return IsAODTrack(track) ? static_cast<const AliAODTrack*>(track)->Chi2perNDF() : static_cast<const AliESDMuonTrack*>(track)->GetNormalizedChi2();
}

//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetChi2MatchTrigger ( const AliVParticle* track )
{
  /// Get the normalized chi2 of the tracker-trigger track matching
  return IsAODTrack(track) ? static_cast<const AliAODTrack*>(track)->GetChi2MatchTrigger() : static_cast<const AliESDMuonTrack*>(track)->GetChi2MatchTrigger();
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
Bool_t AliAnalysisMuonUtility::IsTrkChamberHit( Int_t chamber, const AliVParticle* track )
{
  /// Check if the given tracking chamber has been fired
  if (chamber < 0 || chamber > 9) return kFALSE;
  return IsAODTrack(track) ? static_cast<const AliAODTrack*>(track)->HitsMuonChamber(chamber) : static_cast<const AliESDMuonTrack*>(track)->IsInMuonClusterMap(chamber);
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
Int_t AliAnalysisMuonUtility::GetMuonTrigDevSign ( const AliVParticle* track )
{
  /// Get trigger deviation sign
  return ( IsAODTrack(track) ) ? static_cast<const AliAODTrack*>(track)->GetMuonTrigDevSign() : static_cast<const AliESDMuonTrack*>(track)->GetMuonTrigDevSign();
}

//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetLoCircuit ( const AliVParticle* track )
{
  /// Get local board
  Int_t loCircuit = 0;
  if ( IsAODTrack (track) ) {
    UInt_t pattern = static_cast<const AliAODTrack*>(track)->GetMUONTrigHitsMapTrg();
    loCircuit = AliESDMuonTrack::GetCrossedBoard(pattern);
  }
  else loCircuit = static_cast<const AliESDMuonTrack*>(track)->LoCircuit();
  
  return loCircuit;
}

//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetNTracks ( const AliVEvent* event )
{
  //
  /// Return the number of tracks in event
  //
  if ( fUseSmearedTracks ) {
    TObjArray* smearedTrackList = static_cast<TObjArray*>(event->FindListObject(fSmearedTrackListName));
    if ( smearedTrackList ) return smearedTrackList->GetEntriesFast();
    else return 0;
  }

  return ( IsAODEvent(event) ) ? static_cast<const AliAODEvent*>(event)->GetNumberOfTracks() : static_cast<const AliESDEvent*>(event)->GetNumberOfMuonTracks();
}

//________________________________________________________________________
AliVParticle* AliAnalysisMuonUtility::GetTrack ( Int_t itrack, const AliVEvent* event )
{
  //
  /// Get the current track
  //
  if ( fUseSmearedTracks ) {
    TObjArray* smearedTrackList = static_cast<TObjArray*>(event->FindListObject(fSmearedTrackListName));
    if ( smearedTrackList ) return static_cast<AliVParticle*>(smearedTrackList->At(itrack));
    else return NULL;
  }

  if ( IsAODEvent(event) ) return static_cast<const AliAODEvent*>(event)->GetTrack(itrack);

  return const_cast<AliESDEvent*>(static_cast<const AliESDEvent*> (event))->GetMuonTrack(itrack);
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
Bool_t AliAnalysisMuonUtility::IsAODMCTrack( const AliVParticle* mcParticle )
{
  /// Check if track is from AOD MC
  return ( mcParticle->IsA() == AliAODMCParticle::Class() );
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsMCTrack( const AliVParticle* mcParticle )
{
  /// Check if track is MC track
  return ( mcParticle->IsA() == AliAODMCParticle::Class() || mcParticle->IsA() == AliMCParticle::Class() );
}


//________________________________________________________________________
Double_t AliAnalysisMuonUtility::GetMCVertexZ ( const AliVEvent* event, const AliMCEvent* mcEvent )
{
  /// Get MC vertex Z
  Double_t vz = 0.;
  if ( IsAODEvent(event) ) {
    AliAODMCHeader* aodMCHeader = static_cast<AliAODMCHeader*> (static_cast<const AliAODEvent*>(event)->FindListObject(AliAODMCHeader::StdBranchName()));
    vz = aodMCHeader->GetVtxZ();
  }
  else if ( mcEvent ) vz = mcEvent->GetPrimaryVertex()->GetZ();
  else AliErrorClass("No MC event found");
  return vz;
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::IsPrimary ( const AliVParticle* mcParticle, const AliMCEvent* mcEvent )
{
  /// Check if the particle is primary
  
  if ( IsAODMCTrack(mcParticle) ) return static_cast<const AliAODMCParticle*>(mcParticle)->IsPrimary();
  else if ( mcEvent ) {
    return ( mcParticle->GetLabel() < const_cast<AliMCEvent*>(mcEvent)->GetNumberOfPrimaries() );
  }
  else AliWarningClass("There is no MC info");
  return kFALSE;
}

//________________________________________________________________________
UInt_t AliAnalysisMuonUtility::GetMCProcess ( const AliVParticle* mcParticle )
{
  /// Get MC process
  /// WARNING: the method may fail when running on AODs for old MC production,
  /// the information was propagated to the AOD level only recently
  UInt_t mcProcess = 0;
  if ( IsAODMCTrack(mcParticle) ) {
    mcProcess = const_cast<AliAODMCParticle*>(static_cast<const AliAODMCParticle*>(mcParticle))->GetMCProcessCode();
  }
  else if ( mcParticle->IsA() == AliMCParticle::Class() ) {
    mcProcess = static_cast<const AliMCParticle*>(mcParticle)->Particle()->GetUniqueID();
  }
  else AliWarningClass("There is no MC info");
  return mcProcess;
}

//________________________________________________________________________
UInt_t AliAnalysisMuonUtility::GetStatusCode ( const AliVParticle* mcParticle )
{
  /// Get particle status code
  UInt_t statusCode = 0;
  if ( IsAODMCTrack(mcParticle) ) {
    statusCode = static_cast<const AliAODMCParticle*>(mcParticle)->GetStatus();
  }
  else if ( mcParticle->IsA() == AliMCParticle::Class() ) {
    statusCode = static_cast<const AliMCParticle*>(mcParticle)->Particle()->GetStatusCode();
  }
  else AliWarningClass("There is no MC info");
  return statusCode;
}


//________________________________________________________________________
TString AliAnalysisMuonUtility::GetPassName ( const AliVEventHandler* eventHandler )
{
  /// Get pass name from event
  
  // At present there is no straightforward way to get the pass name.
  // The pass name is usually written in the path to the ESDs/AODs
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
    AliAODHeader * header = dynamic_cast<AliAODHeader*>(event->GetHeader());
    assert(header && "Not a standard AOD");
    filePath = header->GetESDFileName();
    TString passName = GetPassName(filePath.Data());
    if ( passName.IsNull() ) AliWarningClass("Check again with the AOD path");
    else return passName;
  }
  
  filePath = eventHandler->GetTree()->GetCurrentFile()->GetName();
  return GetPassName(filePath.Data());
}


//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetPassNumber ( const AliVEventHandler* eventHandler )
{
  /// Get pass number from event

  // At present there is no straightforward way to get the pass number.
  // The pass number is usually written in the path to the ESDs/AODs
  // but this won't work for:
  // - MC data (no pass info available)
  // - local ESDs/AODs
  
  TString passName = GetPassName(eventHandler);
  return GetPassNumber(passName);
}


//________________________________________________________________________
TString AliAnalysisMuonUtility::GetPassName ( const char* str )
{
  //
  /// Get pass name from string
  //
  
  TString currName(str);
  TObjArray* array = currName.Tokenize("/");
  
  TString passName = "";
  for ( Int_t ientry=0; ientry<array->GetEntries(); ientry++ ) {
    TString currToken = static_cast<TObjString*>(array->At(ientry))->GetString();
    TString checkToken(currToken);
    checkToken.ToUpper();
    if ( checkToken.Contains("PASS") || checkToken.Contains("MUON_CALO") ) {
      passName = currToken;
      break;
    }
  }
  
  delete array;
  
  if ( passName.IsNull() ) AliWarningClass(Form("Cannot find pass name in: %s", str));
  
  return passName;
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

//________________________________________________________________________
Int_t AliAnalysisMuonUtility::GetRunNumber ( const char* str )
{
  /// Get run number string
  TString runNumber = GetRunNumberAsString(str);
  if ( runNumber.IsNull() ) return -1;
  return runNumber.Atoi();
}

//________________________________________________________________________
TString AliAnalysisMuonUtility::GetRunNumberAsString ( const char* str )
{
  /// Get run number in string (as string)
  TString found = "";
  TString queryString(str);
  for ( Int_t ndigits=9; ndigits>=6; ndigits-- ) {
    TString sre = "";
    for ( Int_t idigit=0;idigit<ndigits; idigit++ ) sre += "[0-9]";
    found = queryString(TRegexp(sre.Data()));
    if ( ! found.IsNull() ) break;
  }
  return found;
}

//________________________________________________________________________
TString AliAnalysisMuonUtility::GetTrackHistory ( const AliVParticle* track, const AliMCEvent* mcEvent, Bool_t verbose )
{
  //
  /// Get string containing particle history
  /// Useful when debugging MC
  //
  TString trackHistory = "";
  if ( ! mcEvent ) return trackHistory;
  Int_t fakeMother = 999999999;
  Int_t imother = IsMCTrack(track) ? fakeMother : track->GetLabel();

  while ( imother >= 0 ) {
    const AliVParticle* part = ( imother == fakeMother ) ? track : mcEvent->GetTrack(imother);
    if ( ! part ) break; // In principle not needed...but for some old MC it breaks sometimes
    TParticlePDG* partPdg = TDatabasePDG::Instance()->GetParticle(part->PdgCode());
    TString pname = ( partPdg ) ? partPdg->GetName() : Form("%i",part->PdgCode());
    if ( ! trackHistory.IsNull() ) trackHistory.Append(" <- ");
    if ( imother != fakeMother ) trackHistory.Append(Form("%i ", imother));
    trackHistory.Append(Form("(%s)", pname.Data()));
    if ( verbose ) trackHistory.Append(Form(" [vz %g  mc %i  status %i]", part->Zv(), GetMCProcess(part), GetStatusCode(part)));
    imother = part->GetMother();
  }
  return trackHistory;
}

//________________________________________________________________________
Bool_t AliAnalysisMuonUtility::EAGetTZEROFlags(const AliVEvent* event, Bool_t& backgroundFlag, Bool_t& pileupFlag, Bool_t& satelliteFlag)
{
  // get the TZERO decisions
  // return false if there's no tzero information in this event
  
  Bool_t rv(kFALSE);
  
  if ( event->IsA() == AliESDEvent::Class() )
  {
    const AliESDTZERO* tzero = static_cast<AliESDEvent*>(const_cast<AliVEvent*>(event))->GetESDTZERO();
    if ( tzero )
    {
      backgroundFlag = tzero->GetBackgroundFlag();
      pileupFlag = tzero->GetPileupFlag();
      satelliteFlag = tzero->GetSatellite();
      rv = kTRUE;
    }
  }
  else if ( event->IsA() == AliAODEvent::Class() )
  {
    AliAODTZERO* tzero = static_cast<const AliAODEvent*>(event)->GetTZEROData();
    if ( tzero )
    {
      backgroundFlag = tzero->GetBackgroundFlag();
      pileupFlag = tzero->GetPileupFlag();
      satelliteFlag = tzero->GetSatellite();
      rv = kTRUE;
    }
  }
  else
  {
    AliErrorClass(Form("Unknown class for the event = %s",event->ClassName()));
  }
  
  return rv;
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
  
  if ( axis->GetFirst() == minVarBin && axis->GetLast() == maxVarBin && axis->TestBit(TAxis::kAxisRange) ) return kFALSE;
  
  gridSparse->SetRangeUser(ivar, axis->GetBinCenter(minVarBin), axis->GetBinCenter(maxVarBin));
  
  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisMuonUtility::SetUseSmearedTracks ( Bool_t useSmearedTracks, Bool_t verbose)
{
  //
  /// Set flag to use semared tracks instead of standard ones
  /// Requires the use of AliTaskMuonTrackSmearing
  //
  fUseSmearedTracks = useSmearedTracks;
  if ( verbose ) AliInfoClass(Form("Use smeared tracks set to %i",fUseSmearedTracks));
}
