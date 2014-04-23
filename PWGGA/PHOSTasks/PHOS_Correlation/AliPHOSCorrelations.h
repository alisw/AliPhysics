#ifndef AliPHOSCorrelations_cxx
#define AliPHOSCorrelations_cxx

/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for identified PHOS cluster from pi0 and take korrelation betwen hadron-pi0 angel's.
/// Authors:   Daniil Ponomarenko (Daniil.Ponomarenko@cern.ch)
//                    Dmitry Blau
// 07-Feb-2014

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


#include "TArrayD.h"
#include "AliAnalysisTaskSE.h"

class AliPHOSCorrelations : public AliAnalysisTaskSE 
{
public:
  enum Period { kUndefinedPeriod, kLHC10h, kLHC11h, kLHC13 };
  enum EventSelection { kTotal, kInternalTriggerMaskSelection, kHasVertex, kHasAbsVertex, kHasCentrality,  kCentUnderUpperBinUpperEdge, kHasPHOSClusters, kHasTPCTracks, kTotalSelected };
  enum HibridCheckVeriable { kOnlyHibridTracks, kWithOutHibridTracks, kAllTracks };
  enum TriggerSelection { kNoSelection, kCentralInclusive, kCentralExclusive, kSemiCentralInclusive, kSemiCentralExclusive, kMBInclusive, kMBExclusive };


public:
  AliPHOSCorrelations();
  AliPHOSCorrelations(const char *name, Period period );
  virtual ~AliPHOSCorrelations();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
//  virtual void   Terminate(Option_t *);

  void SetHibridGlobalCheking(Int_t hibridCheck = kAllTracks) {fCheckHibridGlobal = hibridCheck; }
  void SetCentralityBinning(const TArrayD& edges, const TArrayI& nMixed);
  void SetInternalTriggerSelection(TriggerSelection selection) { fInternalTriggerSelection = selection; }
  void EnableTOFCut(Bool_t enable = kTRUE, Double_t TOFCut = 100.e-9){fTOFCutEnabled=enable; fTOFCut=TOFCut;}
  void SetMassWindow(Double_t massMean = 0.135, Double_t massSigma = 0.01) { fMassInvMean = massMean; fMassInvSigma = massSigma; }
  void SetSigmaWidth(Double_t sigmaWidth= 0) { fSigmaWidth = sigmaWidth; }
  void SetPeriod(Period period);
  void SetCentralityBorders (double down = 0., double up = 90.) ;
  void SetPtAssocBins(TArrayD * arr){fAssocBins.Set(arr->GetSize(), arr->GetArray()) ;} 

  void SetCentralityEstimator(const char * centr) {fCentralityEstimator = centr;}
  void SetEventMixingRPBinning(UInt_t nBins) { fNEMRPBins = nBins; }
  void SetMaxAbsVertexZ(Float_t z) { fMaxAbsVertexZ = z; }
  
protected: 

  AliPHOSCorrelations(const AliPHOSCorrelations&);        // not implemented
  AliPHOSCorrelations& operator=(const AliPHOSCorrelations&); // not implemented
  
  // Histograms and trees.
    void SetHistPtAssoc();                      // Set massive of histograms (1-5).
    void SetHistCutDistribution();              // Set other histograms.
    void SetHistEtaPhi();                       // Set hists, with track's and cluster's angle distributions.
    void FillTrackEtaPhi();                     // Distribution by track's angles.
    void SetHistPHOSClusterMap();       // XZE distribution in PHOS.
    void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key.
    void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key.

    void SetESDTrackCuts(); // AliESDtrack cuts ( for esd data )
    
    Bool_t TestMass(Double_t m, Double_t pt) ;
    Double_t GetAssocBin(Double_t pt) ;

    Int_t ConvertToInternalRunNumber(Int_t run);

    void FillTriggerProgress();
    Bool_t RejectTriggerMaskSelection(); 

    void    SetVertex();
    Bool_t RejectEventVertex();

    void  SetCentrality();   // Find centrality of event.
    Bool_t RejectEventCentrality(); 
    

    Int_t     GetCentralityBin(Float_t centralityV0M);
    UInt_t  GetNumberOfCentralityBins() { return fCentEdges.GetSize()-1; }

    void EvalReactionPlane();                               // Find RP of event.
    void EvalV0ReactionPlane();                           // Find RP of event.
    Int_t GetRPBin(); // Return RP (rad).

    Double_t ApplyFlattening(Double_t phi, Double_t c) ;       //Apply centrality-dependent flattening.
    Double_t ApplyFlatteningV0A(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening.
    Double_t ApplyFlatteningV0C(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening.

    
    virtual void SelectPhotonClusters();
    void SelectAccosiatedTracks();

    void ConsiderPi0s();
    void ConsiderPi0sMix();           // MIX for catch Mass
    void ConsiderTracksMix();       // MIX for catch Yeild
    
    TList* GetCaloPhotonsPHOSList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin);
    TList* GetTracksTPCList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin);

    void UpdatePhotonLists();
    void UpdateTrackLists();


    Bool_t SelectESDTrack(AliESDtrack * t) const; //estimate if this track can be used for the RP calculation
    Bool_t SelectAODTrack(AliAODTrack * t) const; //estimate if this track can be used for the RP calculation

    // Logical and debug.
    void LogProgress(int step);
    void LogSelection(int step, int internalRunNumber);

  // Set / Get parametrs
    void SetManualV0EPCalc(Bool_t manCalc = true) {fManualV0EPCalc = manCalc;}

    AliAnalysisUtils* GetAnalysisUtils();
 
private:
  // Geometry
    AliPHOSGeometry* fPHOSGeo;
  // Make output histograms / conteiners.
    TList * fOutputContainer;              //final histogram / tree container     
   
  // cluster cut variables:
    Double_t fMinClusterEnergy;
    Double_t fMinBCDistance;  //distance to nearest bad channel
    Int_t    fMinNCells;
    Double_t fMinM02;
    Bool_t fTOFCutEnabled;
    Double_t fTOFCut;

  // Binning, [vtx, centrality, reaction-plane]
    Int_t   fNVtxZBins;
    TArrayD fCentEdges;                 // Centrality Bin Lower edges.
    TArrayI fCentNMixed;                // Number of mixed events for each centrality bin.
    UInt_t  fNEMRPBins;                  // Binning of Reaction plane.
    TArrayD fAssocBins;                 //  Assoc Pt Bin Lower edges.

  // Control variables
    Int_t fCheckHibridGlobal;      // For checking/dischecking/passingcheck: t->IsHybridGlobalConstrainedGlobal();

  // Behavior / cuts
    Period fPeriod;
    TriggerSelection fInternalTriggerSelection;
    Float_t fMaxAbsVertexZ;       // in cm.
    Bool_t fManualV0EPCalc;

    Double_t fCentCutoffDown;   // Ignore Centrality less %. (def = 0%)
    Double_t fCentCutoffUp;     // Ignore Centrality over %. (def = 90%)

    Double_t fMassInvMean ;      //
    Double_t fMassInvSigma ;      // 
    Double_t fSigmaWidth;       // 0 = wide

    AliVEvent* fEvent;          //! Current event
    AliESDEvent* fEventESD;     //! Current event, if ESD.
    AliAODEvent* fEventAOD;     //! Current event, if AOD.
    AliESDtrackCuts *fESDtrackCuts;     // Track cut

    Int_t fRunNumber;           //! run number
    Int_t fInternalRunNumber ;  //!Current internal run number

    TProfile* fMultV0;                  // object containing VZERO calibration information
    Float_t fV0Cpol,fV0Apol;            // loaded by OADB
    Float_t fMeanQ[9][2][2];    // and recentering
    Float_t fWidthQ[9][2][2];   //       
    TString fEPcalibFileName; 

    Double_t fVertex[3];          //!
    TVector3 fVertexVector;       //!
    Int_t fVtxBin;                //!

    TString fCentralityEstimator; //! Centrality estimator ("V0M", "ZNA")
    Float_t fCentrality ;         //! Centrality of the current event
    Int_t   fCentBin ;            //! Current centrality bin

    Bool_t fHaveTPCRP ; //! Is TPC RP defined?
    Float_t fRP ;       //! Reaction plane calculated with full TPC
    Int_t fEMRPBin;     //! Event Mixing Reaction Plane Bin

    TClonesArray * fCaloPhotonsPHOS ;      //! PHOS photons in current event
    TClonesArray * fTracksTPC ;            //! TPC Tracks in current event

    TObjArray * fCaloPhotonsPHOSLists;  //! array of TList, Containers for events with PHOS photons
    TObjArray * fTracksTPCLists;        //! array of TList, Containers for events with PHOS photons

  ClassDef(AliPHOSCorrelations, 2);    // PHOS analysis task
};

#endif
