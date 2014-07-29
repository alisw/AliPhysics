#ifndef ALICALOTRACKAODREADER_H
#define ALICALOTRACKAODREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class for reading data (AODs) in order to do prompt gamma or other particle
// identification and correlations.
// This part is commented: Mixing analysis can be done, input AOD with events
// is opened in the AliCaloTrackReader::Init()
//
//
// -- Author: Gustavo Conesa (INFN-LNF)

class AliAODEvent;

#include "AliCaloTrackReader.h" 

class AliCaloTrackAODReader : public AliCaloTrackReader {
	
public: 
	
              AliCaloTrackAODReader() ;   // ctor
  
  virtual    ~AliCaloTrackAODReader() {;} // virtual dtor
  
  Bool_t      CheckForPrimaryVertex() const ;
  
  AliVEvent * GetOriginalInputEvent() const { return fOrgInputEvent; }
  
  Bool_t      SelectTrack(AliVTrack* track, Double_t* pTrack);

  ULong_t     GetTrackFilterMask()               const { return fTrackFilterMask      ; }
  void        SetTrackFilterMask(ULong_t bit)          { fTrackFilterMask = bit       ; }
  
  ULong_t     GetTrackFilterMaskComplementary()  const { return fTrackFilterMaskComplementary      ; }
  void        SetTrackFilterMaskComplementary(ULong_t bit) {   fTrackFilterMaskComplementary = bit ; }
  
  void        SwitchOnAODHybridTrackSelection()        { fSelectHybridTracks  = kTRUE  ; }
  void        SwitchOffAODHybridTrackSelection()       { fSelectHybridTracks  = kFALSE ; }
  
  void        SwitchOnAODPrimaryTrackSelection()       { fSelectPrimaryTracks = kTRUE  ; }
  void        SwitchOffAODPrimaryTrackSelection()      { fSelectPrimaryTracks = kFALSE ; }
  
  void        SwitchOnAODTrackSharedClusterSelection() { fSelectFractionTPCSharedClusters = kTRUE  ; }
  void        SwitchOffAODTrackSharedClusterSelection(){ fSelectFractionTPCSharedClusters = kFALSE ; }
  
  void        SetTPCSharedClusterFraction(Float_t fr)  { fCutTPCSharedClustersFraction = fr   ; }
  Float_t     GetTPCSharedClusterFraction() const      { return fCutTPCSharedClustersFraction ; }
  
  void        SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ;
  
private:
  
  AliVEvent * fOrgInputEvent;                   //! Original input event, not from filtering
  
  Bool_t      fSelectHybridTracks;              // Select CTS tracks of type hybrid
  Bool_t      fSelectPrimaryTracks;             // Select CTS tracks of type primary
  ULong_t     fTrackFilterMask;                 // Track selection bit, for AODs (any difference with track status?)
  ULong_t     fTrackFilterMaskComplementary;    // Complementary Track selection bit, for AODs in case hybrid option selected
  Bool_t      fSelectFractionTPCSharedClusters; // Accept only TPC tracks with over a given fraction of shared clusters
  Float_t     fCutTPCSharedClustersFraction;    // Fraction of TPC shared clusters to be accepted.

  
  AliCaloTrackAODReader(              const AliCaloTrackAODReader & r) ; // cpy ctor
  AliCaloTrackAODReader & operator = (const AliCaloTrackAODReader & r) ; // cpy assignment
  
  ClassDef(AliCaloTrackAODReader,7)
  
} ;

#endif //ALICALOTRACKAODREADER_H



