#ifndef ALIANALYSISMUONUTILITY_H
#define ALIANALYSISMUONUTILITY_H

/* $Id: AliAnalysisMuonUtility.h 47782 2011-02-24 18:37:31Z martinez $ */ 

//
// Utilities for muon analysis
//
// Author: Diego Stocco
//

#include "TObject.h"
#include "TString.h"

class TLorentzVector;
class AliVEventHandler;
class AliVEvent;
class AliMCEvent;
class AliVParticle;
class AliVVertex;
class AliCFGridSparse;

class AliAnalysisMuonUtility : public TObject {
 public:
  
  // Utility methods for CF container
  static Bool_t SetSparseRange(AliCFGridSparse* gridSparse,
                               Int_t ivar, TString labelName,
                               Double_t varMin, Double_t varMax,
                               TString option = "");
  
  // Transparently handle ESD/AOD
  static Bool_t IsAODEvent ( const AliVEvent* event );
  static Int_t GetNTracks ( const AliVEvent* event );
  static AliVParticle* GetTrack ( Int_t itrack, const AliVEvent* event );
  
  static Bool_t IsAODTrack ( const AliVParticle* track );
  static Bool_t IsMuonTrack ( const AliVParticle* track );
  static Bool_t IsMuonGhost ( const AliVParticle* track );
  static Double_t GetRabs ( const AliVParticle* track );
  static Double_t GetThetaAbsDeg ( const AliVParticle* track );
  static Int_t GetMatchTrigger ( const AliVParticle* track );
  static Bool_t MatchApt ( const AliVParticle* track ) { return GetMatchTrigger(track) >= 1; }
  static Bool_t MatchLpt ( const AliVParticle* track ) { return GetMatchTrigger(track) >= 2; }
  static Bool_t MatchHpt ( const AliVParticle* track ) { return GetMatchTrigger(track) >= 3; }
  static Double_t GetChi2perNDFtracker ( const AliVParticle* track );
  static Double_t GetChi2MatchTrigger ( const AliVParticle* track );
  static Double_t GetXatVertex ( const AliVParticle* track );
  static Double_t GetYatVertex ( const AliVParticle* track );
  static Double_t GetZatVertex ( const AliVParticle* track );
  static Double_t GetXatDCA ( const AliVParticle* track );
  static Double_t GetYatDCA ( const AliVParticle* track );
  static Double_t GetZatDCA ( const AliVParticle* track ) { return GetZatVertex(track); }
  static Bool_t IsTrkChamberHit( Int_t chamber, const AliVParticle* track );
  static UInt_t GetMUONTrigHitsMapTrk ( const AliVParticle* track );
  static UInt_t GetMUONTrigHitsMapTrg ( const AliVParticle* track );
  static Int_t GetMuonTrigDevSign ( const AliVParticle* track );
  static Int_t GetLoCircuit ( const AliVParticle* track );
  static TLorentzVector GetTrackPair ( const AliVParticle* track1, const AliVParticle* track2 );

  
  // Transparently handle MC
  static Double_t GetMCVertexZ ( const AliVEvent* event, const AliMCEvent* mcEvent );
  
  static Bool_t IsAODMCTrack ( const AliVParticle* mcParticle );
  static Bool_t IsMCTrack ( const AliVParticle* mcParticle );
  static Bool_t IsPrimary ( const AliVParticle* mcParticle, const AliMCEvent* mcEvent );
  static UInt_t GetMCProcess ( const AliVParticle* mcParticle );
  static UInt_t GetStatusCode ( const AliVParticle* mcParticle );
  
  static Bool_t EAGetTZEROFlags(const AliVEvent* event, Bool_t& backgroundFlag, Bool_t& pileupFlag, Bool_t& satelliteFlag);

  // A useful constant
  static Double_t MuonMass2();
  
  // Utilities for ESD/AOD
  static Int_t GetPassNumber ( const AliVEventHandler* eventHandler );
  static TString GetPassName ( const AliVEventHandler* eventHandler );

  // Path utilities
  static Int_t GetPassNumber ( const char* str );
  static TString GetPassName ( const char* str );
  static Int_t GetRunNumber ( const char* str );
  static TString GetRunNumberAsString ( const char* str );

  // Utilities for MC
  static TString GetTrackHistory ( const AliVParticle* track, const AliMCEvent* mcEvent, Bool_t verbose = kFALSE );

  // Utilities for smeared tracks
  static void SetUseSmearedTracks ( Bool_t useSmearedTracks, Bool_t verbose = kTRUE );
  /// Get name of smeared track list
  static const char* GetSmearedTrackListName () { return fSmearedTrackListName; }

private:
  static Bool_t fUseSmearedTracks; //!<! Flag to use smeared tracks
  static const char* fSmearedTrackListName; //!<! Name of the smeared track objects in the InputEvent

  ClassDef(AliAnalysisMuonUtility, 0);
};

#endif
