#ifndef AliPHOSCorrelations_cxx
#define AliPHOSCorrelations_cxx

/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for identifion PHOS cluster from Pi0 and extracting pi0-hadron correlation.
// Author: 	Daniil Ponomarenko <Daniil.Ponomarenko@cern.ch>
// 20-Sept-2014

class TClonesArray;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliTriggerAnalysis;
class AliESDEvent ;
class AliPIDResponse;
class AliPHOSCalibData ;
class AliESDCaloCluster ;
class AliESDEvent ;
class AliESDtrack ;
class AliAODTrack ;
class AliVCluster ;
class AliAnalysisUtils;
class AliEPFlattener;
class AliAODInputHandler;
class AliESDInputHandler;


#include "TArrayD.h"
#include "AliAnalysisTaskSE.h"

class AliPHOSCorrelations : public AliAnalysisTaskSE 
{
public:
  enum Period               { kUndefinedPeriod, kLHC10h, kLHC11h, kLHC13 } ;
  enum EventSelection       { kTotal, kEvent, kEventHandler, 
                              kTriggerMaskSelection, kHasVertex, kHasCentrality, 
                              kHasPHOSClusters, kHasTPCTracks, kPHOSEvent, 
                              kMBEvent, kTotalSelected, kHasAbsVertex } ;
  enum HibridCheckVeriable  { kOnlyHibridTracks, kWithOutHibridTracks, kAllTracks } ;
  enum PID                  { kPidAll, kPidCPV, kPidDisp, kPidBoth} ;


public:
  AliPHOSCorrelations() ;
  AliPHOSCorrelations(const char *name) ;
  AliPHOSCorrelations(const char *name, Period period ) ;
  virtual ~AliPHOSCorrelations() ;

  virtual void   UserCreateOutputObjects() ;
  virtual void   UserExec(Option_t *option) ;

  void SetPeriod(Period period)                                                   { fPeriod = period;                      }
  void SetCentralityEstimator(const char * centr)                                 { fCentralityEstimator = centr;          }
  void SetEventMixingRPBinning(UInt_t nBins)                                      { fNEMRPBins = nBins;                    }
  void SetMaxAbsVertexZ(Float_t z)                                                { fMaxAbsVertexZ = z;                    }
  void SetSigmaWidth(Double_t sigmaWidth)                                         { fSigmaWidth = sigmaWidth;              }
  void SetUseEfficiency(Bool_t useEff)                                            { fUseEfficiency = useEff;               }
  void SetHibridGlobalCheking(Int_t hibridCheck)                                  { fCheckHibridGlobal = hibridCheck;      }
  void EnableTOFCut(Bool_t enable, Double_t TOFCut)                               { fTOFCutEnabled=enable; fTOFCut=TOFCut; }
  void SetMassMeanParametrs(Double_t par[2])  ;
  void SetMassSigmaParametrs(Double_t par[4]) ;
  void SetPtAssocBins(TArrayD * arr)                                              { fAssocBins.Set(arr->GetSize(), arr->GetArray());    } 
  void SetMassWindow(Double_t massMean, Double_t massSigma)                       { fMassInvMean = massMean; fMassInvSigma = massSigma; }
  void SetCentralityBinning(const TArrayD& edges, const TArrayI& nMixed) ;
  void SetCentralityBorders (double down, double up) ;
  

protected: 
  AliPHOSCorrelations           ( const AliPHOSCorrelations& ) ;                                 // not implemented
  AliPHOSCorrelations& operator=( const AliPHOSCorrelations& ) ;                                 // not implemented

  // Filling hists.
  void FillHistogram( const char * key,Double_t x ) const ;                                      // Fill 1D histogram witn name key
  void FillHistogram( const char * key,Double_t x, Double_t y ) const ;                          // Fill 2D histogram witn name key
  void FillHistogram( const char * key,Double_t x, Double_t y, Double_t z ) const ;              // Fill 3D histogram witn name key
  void FillHistogram( const char * key,Double_t x, Double_t y, Double_t z, Double_t w ) const ;  // Fill 3D histogram witn name key

  // Setup hists.
  void SetHistPtNumTrigger( Int_t  ptMult, Double_t ptMin, Double_t ptMax ) ;                    // Set massive of histograms (1-5).
  void SetHistPtAssoc     ( Int_t  ptMult, Double_t ptMin, Double_t ptMax ) ;                    // Set massive of histograms (1-5).
  void SetHistMass        ( Int_t  ptMult, Double_t ptMin, Double_t ptMax ) ;                    // Set other histograms.
  void SetHistEtaPhi() ;                                                                         // Set hists, with track's and cluster's angle distributions.
  void SetHistPHOSClusterMap() ;                                                                 // XZE distribution in PHOS.

  // Logical and debug.
  void LogProgress    ( int step ) ;
  void LogSelection   ( int step , int internalRunNumber ) ;


 
  // Step 1(done once):
  Int_t ConvertToInternalRunNumber(Int_t run) ;                                                  // Convert run number to local number. 
  void SetESDTrackCuts() ;                                                                       // AliESDtrack cuts ( for esd data )

  // Step 2: Preparation variables for new event
  void ZeroingVariables() ;
  void SetGeometry();                                                                            // Initialize the PHOS geometry


  // Step 3: Event trigger selection
  Bool_t RejectTriggerMaskSelection() ;                                                          // Select event trigger and reject.

  // Step 4: Vertex
  void   SetVertex() ;                                                                           // Find vertex of event.
  Bool_t RejectEventVertex() ;

  // Step 5: Centrality
  void   SetCentrality() ;                                                                       // Find centrality of event.
  Bool_t RejectEventCentrality() ; 

  Int_t  GetCentralityBin(Float_t centralityV0M) ;                                               // Find centrality bin.
  UInt_t GetNumberOfCentralityBins() const { return fCentEdges.GetSize()-1 ; }                   // Get number of centrality bins.

  // Step 6: Reaction Plane
  void  EvalReactionPlane() ;                                                                    // Find RP of event.
  Int_t GetRPBin() ;                                                                             // Return RP (rad).

  // Step 7: Event Photons (PHOS Clusters) selection
  virtual void SelectPhotonClusters() ;

  // Step 8: Event Associated particles (TPC Tracks) selection
  void SelectAccosiatedTracks() ;

  // Step 9: Fill TPC's track mask
  void FillTrackEtaPhi() ;                                                                       // Distribution by track's angles.

  // Step 10: Extract one most energetic pi0 candidate in this event. 
  void SelectTriggerPi0ME() ;                                                                    // Select most energetic Pi0 in event.

  void  TestPi0ME(Int_t ipid, TLorentzVector p12, Int_t modCase) ;                              // Compare Pi0 particles and remember most energetic in current event.
 
  void  SetMEExists(const Int_t pid)                        { fMEExists[pid] = true     ; }
  void  SetMEPhi(const Int_t pid, const Double_t phi)       { fMEPhi[pid] = phi         ; }
  void  SetMEEta(const Int_t pid, const Double_t eta)       { fMEEta[pid] = eta         ; }
  void  SetMEPt(const Int_t pid, const Double_t pT)         { fMEPt[pid] = pT           ; }
  void  SetMEModCase(const Int_t pid, const Int_t modcase)  { fMEModCase[pid] = modcase ; }

  Bool_t      GetMEExists(const Int_t pid)    const { return fMEExists[pid]   ; }
  Double_t    GetMEPhi(const Int_t pid)       const { return fMEPhi[pid]      ; }
  Double_t    GetMEEta(const Int_t pid)       const { return fMEEta[pid]      ; }
  Double_t    GetMEPt(const Int_t pid)        const { return fMEPt[pid]       ; }
  Int_t       GetMEModCase(const Int_t pid)   const { return fMEModCase[pid]  ; }

  // Step 11: Start correlation analysis.
  void ConsiderPi0s() ;                       // Consider the most energetic Pi0 in this event with all tracks of this event.
  void ConsiderPi0s_MBSelection() ;           // Consider the most energetic Pi0 in this event with all tracks of this event using MB events.
  
  void ConsiderPi0sMix() ;                    // Use MIX for catch mass peck.
  void ConsiderTracksMix() ;                  // Consider the most energetic Pi0 in this event with all tracks from MIXing pull.

  void UpdatePhotonLists() ;                  // Fill photons in MIXing pull.
  void UpdateTrackLists() ;                   // Fill Tracks in MIXing pull.



  Bool_t TestMass(Double_t m, Double_t pt) ;                                                     // Check if mair in pi0 peak window.

  Double_t MassMeanFunktion(Double_t &pt) const ;                                                // Parametrization mean of mass window.
  Double_t MassSigmaFunktion(Double_t &pt) const ;                                               // Parametrization sigma of mass window.

  Double_t GetAssocBin(Double_t pt) const ;                                                      //Calculates bin for current associated particle pT.

  Double_t GetEfficiency(Double_t pt) const ;                                                    // Return Pi0 efficiency for current pT (PID: both2core only).

  Int_t GetModCase(Int_t &mod1, Int_t &mod2) const ;                                             // Produce part of module neme for pTetaPhi histogram.

  TList* GetCaloPhotonsPHOSList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin) ;                   // Return photons from PHOS list from previous events.
  TList* GetTracksTPCList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin) ;                         // Return tracks from TPC list from previous events.

  Bool_t SelectESDTrack(AliESDtrack * t) const ;                                                 // Estimate if this track can be used for the RP calculation.
  Bool_t SelectAODTrack(AliAODTrack * t) const ;                                                 // Estimate if this track can be used for the RP calculation.

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

  Int_t     fRunNumber;                                 //! Run number
  Int_t     fInternalRunNumber ;                        //! Current internal run number 
  Period    fPeriod;                                    //! kUndefinedPeriod, kLHC10h, kLHC11h, kLHC13

  Bool_t    fPHOSEvent;                                 //! PHOS event trigger.
  Bool_t    fMBEvent;                                   //! MB event trigger.

  // Binning [vtx, centrality, reaction-plane]
  Int_t     fNVtxZBins;                                 // Number of Z vertex bins
  TArrayD   fCentEdges;                                 //! Centrality Bin Lower edges
  TArrayI   fCentNMixed;                                // Number of mixed events for each centrality bin
  UInt_t    fNEMRPBins;                                 // Binning of Reaction plane
  TArrayD   fAssocBins;                                 //! Assoc Pt Bin Lower edges  

  Double_t  fVertex[3];                                 //! Event vertex
  TVector3  fVertexVector;                              //! The same
  Int_t     fVtxBin;                                    //! Vertex bin

  TString   fCentralityEstimator;                       //! Centrality estimator ("V0M", "ZNA")
  Float_t   fCentrality ;                               //! Centrality of the current event
  Int_t     fCentBin ;                                  //! Current centrality bin

  Bool_t    fHaveTPCRP ;                                //! Is TPC RP defined?
  Float_t   fRP ;                                       //! Reaction plane calculated with full TPC
  Int_t     fEMRPBin;                                   //! Event Mixing Reaction Plane Bin

  // Behavior / cuts
  Float_t   fMaxAbsVertexZ;                             // Maximum distence Z component of vertix in cm
  Double_t  fCentralityLowLimit;                        // Ignore Centrality less % 
  Double_t  fCentralityHightLimit;                      // Ignore Centrality over % 

  AliESDtrackCuts *   fESDtrackCuts;                    // Track cut
  Int_t     fCheckHibridGlobal ;                        // For checking/dischecking/passingcheck: t->IsHybridGlobalConstrainedGlobal();

  Double_t  fMinClusterEnergy;                          // Min energy PHOS's cluster
  Double_t  fMinBCDistance;                             // Min distance to nearest bad channel
  Int_t     fMinNCells;                                 // Min count of Cells in cluster
  Double_t  fMinM02;                                    // Min size of M02 in claster
  Bool_t    fTOFCutEnabled;                             // Use time of flight or not?
  Double_t  fTOFCut;                                    // Max time of flight

  Double_t fMassInvMean ;                               // Mass Pi0
  Double_t fMassInvSigma ;                              // Mass width Pi0
  Double_t fSigmaWidth;                                 // Width in sigma (*N). If fSigmaWidth = 0 code will use fMassInvMean+/-fMassInvSigma

  // Funktion of mass window parametrs: [mass, pt]
  Double_t  fMassMean[2];                               // Mass mean parametrisation
  Double_t  fMassSigma[4];                              // Mass sigma parametrisation

  // ME Pi0 selection veriables ([n] = pid).
  Bool_t    fMEExists[4];                               // Does trigger Pi0 candidate exists?
  Double_t  fMEPhi[4];                                  // Phi of ME Pi0 candidate
  Double_t  fMEEta[4];                                  // Eta of ME Pi0 candidate
  Double_t  fMEPt[4];                                   // pT of ME Pi0 candidate
  Int_t     fMEModCase[4];                              // Pair of modules where photons are observed

  Bool_t    fUseEfficiency ;                            // Use efficiensy correction during analysis

  ClassDef(AliPHOSCorrelations, 2);                     // PHOS analysis task
};

#endif
