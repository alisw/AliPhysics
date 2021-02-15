#ifndef ALICALOTRACKAODREADER_H
#define ALICALOTRACKAODREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloTrackAODReader
/// \ingroup CaloTrackCorrelationsBase
/// \brief Class for event, clusters and tracks filtering and preparation for the AOD analysis.
///
/// Class for accessing/filtering AOD data. Most of the job is done in the mother class
/// here only very specific methods of the AOD format are implemented.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________

class AliAODEvent;

#include "AliCaloTrackReader.h" 

class AliCaloTrackAODReader : public AliCaloTrackReader {
	
public: 
	
              AliCaloTrackAODReader() ;   // ctor
  
  virtual    ~AliCaloTrackAODReader() {;} // virtual dtor
  
  Bool_t      CheckForPrimaryVertex()  const ;
    
  TClonesArray   * GetAODMCParticles() const ;
  AliAODMCHeader * GetAODMCHeader   () const ;
  AliGenEventHeader* GetGenEventHeader() const ;

  TList *     GetCreateControlHistograms() ;
    
  TObjString *GetListOfParameters() ;
  
  AliVEvent * GetOriginalInputEvent()            const { return fOrgInputEvent        ; }
  
  Bool_t      SelectTrack(AliVTrack* track, Double_t* pTrack);

  Int_t       GetTrackID(AliVTrack* track) ;
  
  ULong_t     GetTrackFilterMask()               const { return fTrackFilterMask      ; }
  void        SetTrackFilterMask(ULong_t bit)          { fTrackFilterMask = bit       ; }
  
  ULong_t     GetTrackFilterMaskComplementary()  const { return fTrackFilterMaskComplementary       ; }
  void        SetTrackFilterMaskComplementary(ULong_t bit) {    fTrackFilterMaskComplementary = bit ; }
  
  void        SwitchOnAODHybridTrackSelection()        { fSelectHybridTracks  = kTRUE  ; }
  void        SwitchOffAODHybridTrackSelection()       { fSelectHybridTracks  = kFALSE ; }
  
  void        SwitchOnAODPrimaryTrackSelection()       { fSelectPrimaryTracks = kTRUE  ; }
  void        SwitchOffAODPrimaryTrackSelection()      { fSelectPrimaryTracks = kFALSE ; }
  
  void        SwitchOnAODTrackSharedClusterSelection() { fSelectFractionTPCSharedClusters = kTRUE  ; }
  void        SwitchOffAODTrackSharedClusterSelection(){ fSelectFractionTPCSharedClusters = kFALSE ; }
  
  void        SetTPCSharedClusterFraction(Float_t fr)  { fCutTPCSharedClustersFraction = fr   ; }
  Float_t     GetTPCSharedClusterFraction() const      { return fCutTPCSharedClustersFraction ; }
  
  void        SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ;
  
  void        Print(const Option_t * opt) const;

private:
  
  AliVEvent * fOrgInputEvent;                   //!<! Original input event, not from filtering
  
  Bool_t      fSelectHybridTracks;              ///< Select CTS tracks of type hybrid
 
  Bool_t      fSelectPrimaryTracks;             ///< Select CTS tracks of type primary
  
  ULong_t     fTrackFilterMask;                 ///< Track selection bit, for AODs (any difference with track status?)
  
  ULong_t     fTrackFilterMaskComplementary;    ///< Complementary Track selection bit, for AODs in case hybrid option selected
  
  Bool_t      fSelectFractionTPCSharedClusters; ///< Accept only TPC tracks with over a given fraction of shared clusters
  
  Float_t     fCutTPCSharedClustersFraction;    ///< Fraction of TPC shared clusters to be accepted.

  TH1F  *     fhCTSAODTrackCutsPt[8];           //!<! control histogram on the different CTS tracks selection cuts, pT
  TH2F  *     fhCTSAODTrackCutsPtCen[8];        //!<! control histogram on the different CTS tracks selection cuts, pT vs centrality
  TH1F  *     fhCTSAODTrackCutsPtSignal[8];     //!<! control histogram on the different CTS tracks selection cuts, pT. Embedded signal
  TH2F  *     fhCTSAODTrackCutsPtCenSignal[8];  //!<! control histogram on the different CTS tracks selection cuts, pT vs centrality. Embedded signal
  
  /// Copy constructor not implemented.
  AliCaloTrackAODReader(              const AliCaloTrackAODReader & r) ; 
  
  /// Assignment operator not implemented.
  AliCaloTrackAODReader & operator = (const AliCaloTrackAODReader & r) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliCaloTrackAODReader,9) ;
  /// \endcond

} ;

#endif //ALICALOTRACKAODREADER_H



