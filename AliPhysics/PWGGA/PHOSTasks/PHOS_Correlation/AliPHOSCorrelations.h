#ifndef AliPHOSCorrelations_cxx
#define AliPHOSCorrelations_cxx

/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for pi0-hadron correlation whis PHOS detector.
// Author:  Daniil Ponomarenko <Daniil.Ponomarenko@cern.ch>
// 20-Feb-2015

class TClonesArray;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliESDEvent ;
class AliESDtrack ;
class AliAODTrack ;
class AliAnalysisUtils;
class AliAODInputHandler;
class AliESDInputHandler;

#ifndef ROOT_TAliAnalysisTaskSE
#include "AliAnalysisTaskSE.h"
#endif

#ifndef ROOT_TArrayD
#include "TArrayD.h"
#endif

class AliPHOSCorrelations : public AliAnalysisTaskSE 
{
public:
  enum EventSelection       { kTotal, kEvent, kEventHandler, 
                              kTriggerMaskSelection, kHasVertex, kHasCentrality, 
                              kHasPHOSClusters, kHasTPCTracks, kPHOSEvent, 
                              kMBEvent, kTotalSelected, kHasAbsVertex } ;
  enum PID                  { kPidAll, kPidCPV, kPidDisp, kPidBoth} ;
  enum PbPbTriggerSelection { kCentAndSCent, kCent, kSCent } ;


public:
  AliPHOSCorrelations() ;
  AliPHOSCorrelations(const char *name) ;
  virtual ~AliPHOSCorrelations() ;

  virtual void   UserCreateOutputObjects() ;
  virtual void   UserExec(Option_t *option) ;

  void SetPeriodName(const TString str)                                           { fPeriod = str                        ; }
  TString GetPeriod()                                                       const { return fPeriod                       ; }

  void SetCentralityEstimator(const char * centr)                                 { fCentralityEstimator = centr         ; }
  void SetEventMixingRPBinning(const Int_t nBins)                                 { fNEMRPBins = nBins                   ; }
  void SetEventMixingVtxBinning(const Int_t nBins) ;
  void SetMaxAbsVertexZ(const Float_t z)                                          { fMaxAbsVertexZ = z                   ; }
 
  void SwitchOnPionEfficiency()                                                   { fUseEfficiency = kTRUE               ; }
  void SwitchOffPionEfficiency()                                                  { fUseEfficiency = kFALSE              ; }

  void EnableTOFCut(const Bool_t enable, const Double_t TOFCut)                   { fTOFCutEnabled=enable; fTOFCut=TOFCut; }
  
  void SetSigmaWidth(const Double_t sigmaWidth)                                   { fNSigmaWidth = sigmaWidth            ; }
  void SetMassMeanParametrs(const Double_t par[2])  ;
  void SetMassSigmaParametrs(const Double_t par[4]) ;
  void SetMassWindow(const Double_t min, const Double_t max)                      { fMassInvMeanMin = min; fMassInvMeanMax = max    ; }
  void SetPtAssocBins(TArrayD * arr)                                              { fAssocBins.Set(arr->GetSize(), arr->GetArray()) ; } 
  void SetCentralityBinning(const TArrayD& edges, const TArrayI& nMixed) ;
  void SetCentralityBorders(const Double_t& downLimit , const Double_t& upLimit) ;
  
  void SwitchOnMassParametrisation()                                              { fUseMassWindowParametrisation = true ; }
  void SwitchOffMassParametrisation()                                             { fUseMassWindowParametrisation = false; }

  ULong_t  GetTrackStatus()                                                 const { return fTrackStatus          ; }
  void     SetTrackStatus(ULong_t bit)                                            { fTrackStatus = bit           ; }   

  ULong_t  GetTrackFilterMask()                                             const { return fTrackFilterMask      ; }
  void     SetTrackFilterMask(ULong_t bit)                                        { fTrackFilterMask = bit       ; }

  void     SwitchOnTrackHitSPDSelection()                                         { fSelectSPDHitTracks = kTRUE  ; }
  void     SwitchOffTrackHitSPDSelection()                                        { fSelectSPDHitTracks = kFALSE ; }

  void     SwitchOnAODTrackSharedClusterSelection()                               { fSelectFractionTPCSharedClusters = kTRUE  ; }
  void     SwitchOffAODTrackSharedClusterSelection()                              { fSelectFractionTPCSharedClusters = kFALSE ; }

  Float_t  GetTPCSharedClusterFraction()                                    const { return fCutTPCSharedClustersFraction ; }
  void     SetTPCSharedClusterFraction(Float_t fr)                                { fCutTPCSharedClustersFraction = fr   ; }

  void     SwitchOnAODHybridTrackSelection()                                      { fSelectHybridTracks = kTRUE         ; } 
  void     SwitchOffAODHybridTrackSelection()                                     { fSelectHybridTracks = kFALSE        ; } 

  void     ShowTaskInfo()  ;

  void     SetEventPlaneMethod(TString m)                                         { fEventPlaneMethod = m               ; }
  TString  GetEventPlaneMethod()                                            const { return fEventPlaneMethod            ; }

  void SetTriggerSelectionInPbPb(PbPbTriggerSelection s)                          { fTriggerSelectionPbPb = s           ; }
  

protected: 
  AliPHOSCorrelations           ( const AliPHOSCorrelations& ) ;                                 // not implemented
  AliPHOSCorrelations& operator=( const AliPHOSCorrelations& ) ;                                 // not implemented

  // Filling hists.
  void FillHistogram( const char * key,Double_t x ) const ;                                      // Fill 1D histogram witn name key
  void FillHistogram( const char * key,Double_t x, Double_t y ) const ;                          // Fill 2D histogram witn name key
  void FillHistogram( const char * key,Double_t x, Double_t y, Double_t z ) const ;              // Fill 3D histogram witn name key
  void FillHistogram( const char * key,Double_t x, Double_t y, Double_t z, Double_t w ) const ;  // Fill 3D histogram witn name key

  // Setup hists.
  void SetHistPtNumTrigger( const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax ) const ; // Set trigger's number of histograms (1-5).
  void SetHistPtAssoc     ( const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax ) const ; // Set pt associated of histograms (1-5).
  void SetHistMass        ( const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax ) const ; // Set mass histograms.
  void SetHistEtaPhi      ( const Int_t& ptMult, const Double_t& ptMin, const Double_t& ptMax ) const ; // Set hists with track's (pt depend.) and cluster's angle distributions.
  void SetHistPHOSClusterMap() ;                                                                        // XZE distribution in PHOS.

  // Logical and debug.
  void LogProgress    ( int step ) ;
  void LogSelection   ( const int& step , const int& internalRunNumber ) const ;

  // Step 1(done once):
  Int_t ConvertToInternalRunNumber(const Int_t& run) const ;                                     // Convert run number to local number. 
  void  SetESDTrackCuts() ;                                                                      // AliESDtrack cuts ( for esd data )

  // Step 2: Preparation variables for new event
  void ZeroingVariables() ;
  void SetGeometry();                                                                            // Initialize the PHOS geometry

  // Step 3: Event trigger selection
  Bool_t RejectTriggerMaskSelection() ;                                                          // Select event trigger and reject.

  // Step 4: Vertex
  void   SetVertex() ;                                                                           // Find vertex of event.
  Bool_t RejectEventVertex()         const ;

  void   SetVertexBinning() ;                                                                    // Define vertex bins by their edges
  Int_t  GetVertexBin(const TVector3&  vertexVector) ;                                           // Find vertex bin
  UInt_t GetNumberOfVertexBins()     const { return fNVtxZBins ; }                               // Get number of vertex bins.

  // Step 5: Centrality
  void   SetCentrality() ;                                                                       // Find centrality of event.
  Bool_t RejectEventCentrality()     const; 

  Int_t  GetCentralityBin(const Float_t& centralityV0M) ;                                        // Find centrality bin.
  UInt_t GetNumberOfCentralityBins() const { return fCentEdges.GetSize()-1 ; }                   // Get number of centrality bins.

  // Step 6: Reaction Plane
  void   EvalReactionPlane() ;                                                                   // Find RP of event.
  Int_t  GetRPBin() ;                                                                            // Return RP (rad).
  UInt_t GetNumberOfRPBins()         const { return fNEMRPBins ; }                               // Get number of RP bins.

  // Step 7: Event Photons (PHOS Clusters) selection
  void SelectPhotonClusters() ;

  // Step 8: Event Associated particles (TPC Tracks) selection
  void SelectAccosiatedTracks() ;

  // Step 9: Fill TPC's track mask and control bining hists.
  void FillTrackEtaPhi()             const;                                                      // Distribution by track's angles.
  void FillEventBiningProperties()   const ;                                                     // Fill fCentBin, fEMRPBin, fVtxBin.

  // Step 10: Extract one most energetic pi0 candidate in this event. 
  void SelectTriggerPi0ME() ;                                                                    // Select most energetic Pi0 in event.

  void  TestPi0ME(const Int_t& ipid, const TLorentzVector& p12, const Int_t& modCase) ;          // Compare Pi0 particles and remember most energetic in current event.
 
  void  SetMEExists(const Int_t pid)                        { fMEExists[pid] = true     ; }
  void  SetMEPhi(const Int_t pid, const Double_t phi)       { fMEPhi[pid] = phi         ; }
  void  SetMEEta(const Int_t pid, const Double_t eta)       { fMEEta[pid] = eta         ; }
  void  SetMEPt( const Int_t pid, const Double_t pT)        { fMEPt[pid] = pT           ; }
  void  SetMEModCase(const Int_t pid, const Int_t modcase)  { fMEModCase[pid] = modcase ; }

  Bool_t   GetMEExists(const Int_t pid)               const { return fMEExists[pid]     ; }
  Double_t GetMEPhi(const Int_t pid)                  const { return fMEPhi[pid]        ; }
  Double_t GetMEEta(const Int_t pid)                  const { return fMEEta[pid]        ; }
  Double_t GetMEPt(const Int_t pid)                   const { return fMEPt[pid]         ; }
  Int_t    GetMEModCase(const Int_t pid)              const { return fMEModCase[pid]    ; }

  // Step 11: Start correlation analysis.
  void ConsiderPi0s() ;                       // Consider the most energetic Pi0 in this event with all tracks of this event.
  void ConsiderPi0s_MBSelection() ;           // Consider the most energetic Pi0 in this event with all tracks of this event using MB events.
  
  void ConsiderPi0sMix() ;                    // Use MIX for catch mass peck.
  void ConsiderTracksMix() ;                  // Consider the most energetic Pi0 in this event with all tracks from MIXing pull.

  void UpdatePhotonLists() ;                  // Fill photons in MIXing pull.
  void UpdateTrackLists() ;                   // Fill Tracks in MIXing pull.


  Bool_t TestMass(const Double_t& m, const Double_t& pt)  const ;                                // Check if mair in pi0 peak window.

  Double_t MassMeanFunction(const Double_t &pt)           const ;                                // Parametrization mean of mass window.
  Double_t MassSigmaFunction(const Double_t &pt)          const ;                                // Parametrization sigma of mass window.
  Double_t GetAssocBin(const Double_t& pt)                const ;                                // Calculates bin for current associated particle pT.
  Double_t GetEfficiency(const Double_t& pt)              const ;                                // Return Pi0 efficiency for current pT (PID: both2core only).
  Int_t GetModCase(const Int_t &mod1, const Int_t &mod2)  const ;                                // Produce part of module neme for pTetaPhi histogram.

  TList* GetCaloPhotonsPHOSList(const UInt_t vtxBin, const UInt_t centBin, const UInt_t rpBin) ; // Return photons from PHOS list from previous events.
  TList* GetTracksTPCList(const UInt_t vtxBin, const UInt_t centBin, const UInt_t rpBin) ;       // Return tracks from TPC list from previous events.

  Bool_t SelectESDTrack(const AliESDtrack * t) const ;                                           // Estimate if this track can be used for the RP calculation.
  Bool_t SelectAODTrack(const AliAODTrack * t) const ;                                           // Estimate if this track can be used for the RP calculation.

  AliAnalysisUtils* GetAnalysisUtils() ;


private:
  //General Data members
  AliPHOSGeometry *   fPHOSGeo ;                        //! Geometry
  TList *   fOutputContainer ;                          //! Output histograms container 

  AliVEvent   *           fEvent;                       //! Current event
  AliESDEvent *           fEventESD;                    //! Current event, if ESD.
  AliAODEvent *           fEventAOD;                    //! Current event, if AOD.
  AliInputEventHandler *  fEventHandler;                //! Event trigger bit.

  TClonesArray *  fCaloPhotonsPHOS ;                    //! PHOS photons in current event
  TClonesArray *  fTracksTPC ;                          //! TPC tracks in current event
  TObjArray *     fCaloPhotonsPHOSLists;                //! array of TList, Containers for events with PHOS photons
  TObjArray *     fTracksTPCLists;                      //! array of TList, Containers for events with TPC tracks

  TString   fPeriod;                                    //  Event period name
  Int_t     fRunNumber;                                 //! Run number
  Int_t     fInternalRunNumber ;                        //! Current internal run number 

  PbPbTriggerSelection fTriggerSelectionPbPb;           //  Selection trigger mask in central PbPb collisions.
  Bool_t    isREALEvent;                                //! PHOS event trigger.
  Bool_t    fMBEvent;                                   //! MB event trigger.

  // Binning [vtx, centrality, reaction-plane]
  Int_t     fNVtxZBins;                                 // Number of Z vertex bins
  TArrayD   fVtxEdges;                                  // Vertex bin Lower edges
  TString   fCentralityEstimator;                       // Centrality estimator ("V0M", "ZNA")
  TArrayD   fCentEdges;                                 // Centrality Bin Lower edges
  TArrayI   fCentNMixed;                                // Number of mixed events for each centrality bin
  TString   fEventPlaneMethod;                          // Name of event plane method, by default "Q"
  UInt_t    fNEMRPBins;                                 // Binning of Reaction plane
  TArrayD   fAssocBins;                                 // Assoc Pt Bin Lower edges  

  Double_t  fVertex[3];                                 //! Event vertex
  TVector3  fVertexVector;                              //! The same
  Int_t     fVtxBin;                                    //! Vertex bin

  Float_t   fCentrality ;                               //! Centrality of the current event
  Int_t     fCentBin ;                                  //! Current centrality bin

  Bool_t    fHaveTPCRP ;                                //! Is TPC RP defined?
  Float_t   fRP ;                                       //! Reaction plane calculated with full TPC
  Int_t     fEMRPBin;                                   //! Event Mixing Reaction Plane Bin

  // Behavior / cuts
  Float_t   fMaxAbsVertexZ;                             // Maximum distence Z component of vertix in cm
  Double_t  fCentralityLowLimit;                        // Ignore Centrality less % 
  Double_t  fCentralityHightLimit;                      // Ignore Centrality over % 

  Double_t  fMinClusterEnergy;                          // Min energy PHOS's cluster
  Double_t  fMinBCDistance;                             // Min distance to nearest bad channel
  Int_t     fMinNCells;                                 // Min count of Cells in cluster
  Double_t  fMinM02;                                    // Min size of M02 in claster
  Bool_t    fTOFCutEnabled;                             // Use time of flight or not?
  Double_t  fTOFCut;                                    // Max time of flight

  Bool_t   fUseMassWindowParametrisation;               // Use parametrization? (Else use fixed mass borders)
  Double_t fMassInvMeanMin ;                            // Mass Pi0 minimum window value
  Double_t fMassInvMeanMax ;                            // Mass Pi0 maximum window value
  Double_t fNSigmaWidth;                                // Width in sigma (*N). If fNSigmaWidth = 0 code will use window fMassInvMeanMin fMassInvMeanMax

  // Funktion of mass window parametrs: [mass, pt]
  Double_t  fMassMean[2];                               // Mass mean parametrisation
  Double_t  fMassSigma[3];                              // Mass sigma parametrisation

  // ME Pi0 selection veriables ([n] = pid).
  Bool_t    fMEExists[4];                               // Does trigger Pi0 candidate exists?
  Double_t  fMEPhi[4];                                  // Phi of ME Pi0 candidate
  Double_t  fMEEta[4];                                  // Eta of ME Pi0 candidate
  Double_t  fMEPt[4];                                   // pT of ME Pi0 candidate
  Int_t     fMEModCase[4];                              // Pair of modules where photons are observed

  Bool_t    fUseEfficiency ;                            // Use efficiensy correction during analysis

  AliESDtrackCuts *  fESDtrackCuts ;                    // Track cut
  
  Bool_t    fSelectHybridTracks ;                       // Select CTS tracks of type hybrid (only for AODs)

  ULong_t   fTrackStatus        ;                       // the statusflag contains a set of flags  which show which parts of the tracking algorithm were successful.
  ULong_t   fTrackFilterMask    ;                       // The filterbits encode combinations of cuts on the AOD, such as the hybrid track selection.
  Bool_t    fSelectSPDHitTracks ;                       // Ensure that track hits SPD layers
  Bool_t    fSelectFractionTPCSharedClusters ;          // Accept only TPC tracks with over a given fraction of shared clusters
  Float_t   fCutTPCSharedClustersFraction    ;          // Fraction of TPC shared clusters to be accepted.

  ClassDef(AliPHOSCorrelations, 5);                     // PHOS analysis task
};

#endif

