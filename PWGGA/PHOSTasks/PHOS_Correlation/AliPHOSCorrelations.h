#ifndef AliPHOSCorrelations_cxx
#define AliPHOSCorrelations_cxx

/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for identified PHOS cluster from pi0 and extracting pi0-hadron correlation.
// Authors: 	Daniil Ponomarenko <Daniil.Ponomarenko@cern.ch>
// 		Dmitry Blau <Dmitry.Blau@cern.ch>
// 09-Jul-2014

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
  enum Period { kUndefinedPeriod, kLHC10h, kLHC11h, kLHC13 };
  enum EventSelection { kTotal, kEvent, kEventHandler, kTriggerMaskSelection, kHasVertex, kHasCentrality, kHasPHOSClusters, kHasTPCTracks, kPHOSEvent, kMBEvent, kTotalSelected, kHasAbsVertex };
  enum HibridCheckVeriable { kOnlyHibridTracks, kWithOutHibridTracks, kAllTracks };
  enum PID { kPidAll, kPidCPV, kPidDisp, kPidBoth};


public:
  AliPHOSCorrelations();
  AliPHOSCorrelations(const char *name);
  AliPHOSCorrelations(const char *name, Period period );
  virtual ~AliPHOSCorrelations();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
//  virtual void   Terminate(Option_t *);

  void SetHibridGlobalCheking(Int_t hibridCheck = kAllTracks) {fCheckHibridGlobal = hibridCheck; }
  void SetAnalysisAlgoritmForReal(TString algoritm = "ME") {algoritm.Contains("ME")?fUseMEAlgoritmForReal = true:fUseMEAlgoritmForReal = false;}
  void SetAnalysisAlgoritmForMix(TString algoritm = "ME") {algoritm.Contains("ME")?fUseMEAlgoritmForMix = true:fUseMEAlgoritmForMix = false;}
  void SetCentralityBinning(const TArrayD& edges, const TArrayI& nMixed);
  void EnableTOFCut(Bool_t enable = kTRUE, Double_t TOFCut = 100.e-9){fTOFCutEnabled=enable; fTOFCut=TOFCut;}
  void SetMassWindow(Double_t massMean = 0.135, Double_t massSigma = 0.01) { fMassInvMean = massMean; fMassInvSigma = massSigma; }
  void SetSigmaWidth(Double_t sigmaWidth= 0) { fSigmaWidth = sigmaWidth; }
  void SetMassMeanParametrs(Double_t p0 = -20.9476, Double_t p1 = 0.1300) {fMassMeanP0 = p0; fMassMeanP1 = p1;}   // from mass fit
  void SetMassSigmaParametrs(Double_t p0 = 0.001, Double_t p1 = -0.0000001, Double_t p2 = -0.06, Double_t p3 = -0.01) {fMassSigmaP0 = p0; fMassSigmaP1 = p1; fMassSigmaP2 = p2; fMassSigmaP3 = p3;}    // from mass fit
  void SetPeriod(Period period) { fPeriod = period; }
  void SetCentralityBorders (double down = 0., double up = 90.) ;
  void SetUseMoreCorrFunctions(Bool_t makeForPHOS = false, Bool_t makeForTPC = false) {fMakePHOSModulesCorrFunctions = makeForPHOS; fMakeTPCHalfBarrelCorrFunctions = makeForTPC; }
  void SetUseEfficiency(Bool_t useEff = true) {fUseEfficiency = useEff;}
  void SetPtAssocBins(TArrayD * arr){fAssocBins.Set(arr->GetSize(), arr->GetArray()) ;} 

  void SetCentralityEstimator(const char * centr) {fCentralityEstimator = centr;}
  void SetEventMixingRPBinning(UInt_t nBins) { fNEMRPBins = nBins; }
  void SetMaxAbsVertexZ(Float_t z) { fMaxAbsVertexZ = z; }
  
protected: 

  AliPHOSCorrelations(const AliPHOSCorrelations&);        // not implemented
  AliPHOSCorrelations& operator=(const AliPHOSCorrelations&); // not implemented
  
  // Histograms and trees.
    void SetHistPtNumTrigger(Int_t  ptMult, Double_t ptMin, Double_t ptMax);      // Set massive of histograms (1-5).
    void SetHistPtAssoc(Int_t  ptMult, Double_t ptMin, Double_t ptMax);           // Set massive of histograms (1-5).
    void SetHistMass(Int_t  ptMult, Double_t ptMin, Double_t ptMax);              // Set other histograms.
    void SetHistEtaPhi();                       // Set hists, with track's and cluster's angle distributions.
    void SetHistPHOSClusterMap();               // XZE distribution in PHOS.
    void FillHistogram(const char * key,Double_t x) const ;                                     //Fill 1D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y) const ;                         //Fill 2D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ;             //Fill 3D histogram witn name key
    void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key

    void SetESDTrackCuts(); // AliESDtrack cuts ( for esd data )
    
    Bool_t TestMass(Double_t m, Double_t pt) ;
    Double_t MassMeanFunktion(Double_t &pt) const ;
    Double_t MassSigmaFunktion(Double_t &pt) const ;

    Double_t GetAssocBin(Double_t pt) const ;   //Calculates bin of associated particle pt.
    Double_t GetEfficiency(Double_t pt) const ; // Return Pi0 efficiency for current pT.

    Int_t GetModCase(Int_t &mod1, Int_t &mod2) const; // Produce part of module neme for pTetaPhi histogram in mixed events.

    Int_t ConvertToInternalRunNumber(Int_t run);

    Bool_t RejectTriggerMaskSelection();    // Select event trigger and reject.

    void    SetVertex();
    Bool_t RejectEventVertex();

    void  SetCentrality();   // Find centrality of event.
    Bool_t RejectEventCentrality(); 
    
    Int_t     GetCentralityBin(Float_t centralityV0M);
    UInt_t  GetNumberOfCentralityBins() { return fCentEdges.GetSize()-1; }

    void EvalReactionPlane();   // Find RP of event.
    void EvalV0ReactionPlane(); // Find RP of event.
    Int_t GetRPBin();           // Return RP (rad).

    Double_t ApplyFlattening(Double_t phi, Double_t c) ;    // Apply centrality-dependent flattening.
    Double_t ApplyFlatteningV0A(Double_t phi, Double_t c) ; // Apply centrality-dependent flattening.
    Double_t ApplyFlatteningV0C(Double_t phi, Double_t c) ; // Apply centrality-dependent flattening.

    void ZeroingVariables();
    
    virtual void SelectPhotonClusters();
    void SelectAccosiatedTracks();

     void FillTrackEtaPhi();        // Distribution by track's angles.

    void SelectTriggerPi0ME();      //Select most energetic Pi0 in event.

    void ConsiderPi0s();            // Consider all Pi0 with all tracks in same event.
    void ConsiderTracksMix();       // Consider all Pi0 in this event with tracks from MIXing pull.

    void ConsiderPi0sME();             // Consider the most energetic Pi0 in this event with all tracks of this event.
    void ConsiderPi0sME_MBSelection(); // Consider the most energetic Pi0 in this event with all tracks of this event using MB events.
    void ConsiderTracksMixME();        // Consider the most energetic Pi0 in this event with all tracks from MIXing pull.

    void ConsiderPi0sMix();           // MIX for catch Mass

    void TestPi0ME(Int_t ipid, TLorentzVector p12, Int_t modCase);  // Compare Pi0 particles and save most energetic.
    Int_t CheckTriggerEta(Double_t eta);                            // Return 1 if eta>=0, else 2.
    
    TList* GetCaloPhotonsPHOSList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin);
    TList* GetTracksTPCList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin);

    void UpdatePhotonLists();   // Fill photons in MIXing pull.
    void UpdateTrackLists();    // Fill Tracks in MIXing pull.

    void SetGeometry();

    Bool_t SelectESDTrack(AliESDtrack * t) const; //estimate if this track can be used for the RP calculation
    Bool_t SelectAODTrack(AliAODTrack * t) const; //estimate if this track can be used for the RP calculation

    // Logical and debug.
    void LogProgress(int step);
    void LogSelection(int step, int internalRunNumber);

  // Set / Get parametrs
    void SetManualV0EPCalc(Bool_t manCalc = true) {fManualV0EPCalc = manCalc;}

    void SetMEExists(const Int_t pid) {fMEExists[pid] = true;}
    Bool_t GetMEExists(const Int_t pid) const {return fMEExists[pid];}

    void SetMEPhi(const Int_t pid, const Double_t phi) {fMEPhi[pid] = phi;}
    Double_t GetMEPhi(const Int_t pid) const {return fMEPhi[pid];}

    void SetMEEta(const Int_t pid, const Double_t eta) {fMEEta[pid] = eta;}
    Double_t GetMEEta(const Int_t pid) const {return fMEEta[pid];}

    void SetMEPt(const Int_t pid, const Double_t pT) {fMEPt[pid] = pT;}
    Double_t GetMEPt(const Int_t pid) const {return fMEPt[pid];}

    void SetMEModCase(const Int_t pid, const Int_t modcase) {fMEModCase[pid] = modcase;}
    Int_t GetMEModCase(const Int_t pid) const {return fMEModCase[pid];}


    AliAnalysisUtils* GetAnalysisUtils();
 
private:
  // Geometry
    AliPHOSGeometry* fPHOSGeo;
  // Make output histograms/conteiners.
    TList * fOutputContainer;              //final histogram / tree container     
   
  // cluster cut variables:
    Double_t fMinClusterEnergy; // Min energy PHOS's cluster.
    Double_t fMinBCDistance;    // Min distance to nearest bad channel
    Int_t    fMinNCells;        // Min count of Cells in cluster.
    Double_t fMinM02;           // Min size of M02 in claster.
    Bool_t fTOFCutEnabled;      // Use time of flight or not?
    Double_t fTOFCut;           // Max time of flight.

  // Binning, [vtx, centrality, reaction-plane]
    Int_t   fNVtxZBins;
    TArrayD fCentEdges;                 // Centrality Bin Lower edges.
    TArrayI fCentNMixed;                // Number of mixed events for each centrality bin.
    UInt_t  fNEMRPBins;                 // Binning of Reaction plane.
    TArrayD fAssocBins;                 //  Assoc Pt Bin Lower edges.

  // Control variables
    Bool_t fUseMEAlgoritmForReal;        // Use common or ME algoritm for analysis real events.
    Bool_t fUseMEAlgoritmForMix;         // Use common or ME algoritm for analysis mixed events.
    Bool_t fUseEfficiency;                // Use efficiensy during analysis.
    Bool_t fMakePHOSModulesCorrFunctions; // Turn on filling module Phi/Eta/Pt distribution.
    Bool_t fMakeTPCHalfBarrelCorrFunctions; // Turn on filling half barrel TPC distribution.
    Int_t fCheckHibridGlobal;      // For checking/dischecking/passingcheck: t->IsHybridGlobalConstrainedGlobal();

  // Event selection
    Bool_t fPHOSEvent;              // PHOS event trigger.
    Bool_t fMBEvent;                // MB event trigger.

  // Behavior / cuts
    Period fPeriod;
    Float_t fMaxAbsVertexZ;       // Maximum distence Z component of vertix in cm.
    Bool_t fManualV0EPCalc;       //

    Double_t fCentCutoffDown;   // Ignore Centrality less %. (def = 0%)
    Double_t fCentCutoffUp;     // Ignore Centrality over %. (def = 90%)

    Double_t fMassInvMean ;     // Mass Pi0.
    Double_t fMassInvSigma ;    // Mass width Pi0.
    Double_t fSigmaWidth;       // *N sigma. if 0 will use fMassInvMean+/-fMassInvSigma. Else will calculate using function.

  // Funktion of window mass parametrs: [mass, pt]
    Double_t fMassMeanP0;
    Double_t fMassMeanP1;
    Double_t fMassSigmaP0;
    Double_t fMassSigmaP1;
    Double_t fMassSigmaP2;
    Double_t fMassSigmaP3;

    AliVEvent* fEvent;          //! Current event
    AliESDEvent* fEventESD;     //! Current event, if ESD.
    AliAODEvent* fEventAOD;     //! Current event, if AOD.
    AliInputEventHandler *fEventHandler; //! Event trigger bit.
    AliESDtrackCuts *fESDtrackCuts;     // Track cut

    Int_t fRunNumber;           //! run number
    Int_t fInternalRunNumber ;  //! Current internal run number

    TProfile* fMultV0;          // Object containing VZERO calibration information
    Float_t fV0Cpol,fV0Apol;    // oaded by OADB
    Float_t fMeanQ[9][2][2];    // and recentering
    Float_t fWidthQ[9][2][2];   //       
    TString fEPcalibFileName;   //

    Double_t fVertex[3];          //! Event vertex.
    TVector3 fVertexVector;       //! The same.
    Int_t fVtxBin;                //! Vertex bin.

    TString fCentralityEstimator; //! Centrality estimator ("V0M", "ZNA")
    Float_t fCentrality ;         //! Centrality of the current event
    Int_t   fCentBin ;            //! Current centrality bin

    Bool_t fHaveTPCRP ; //! Is TPC RP defined?
    Float_t fRP ;       //! Reaction plane calculated with full TPC
    Int_t fEMRPBin;     //! Event Mixing Reaction Plane Bin

    // ME Pi0 selection veriables. 1...4 = all...both pid.
    Bool_t fMEExists[4];    // Does trigger Pi0 exists?
    Double_t fMEPhi[4];     // Phi ME pi0.
    Double_t fMEEta[4];     // Eta ME Pi0.
    Double_t fMEPt[4];      // pT ME Pi0.
    Int_t fMEModCase[4];    // Pair of modules where photons are observed.

    TClonesArray * fCaloPhotonsPHOS ;      //! PHOS photons in current event
    TClonesArray * fTracksTPC ;            //! TPC Tracks in current event

    TObjArray * fCaloPhotonsPHOSLists;  //! array of TList, Containers for events with PHOS photons
    TObjArray * fTracksTPCLists;        //! array of TList, Containers for events with PHOS photons

  ClassDef(AliPHOSCorrelations, 2);    // PHOS analysis task
};

#endif
