/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskFlowModes_pPb_H
#define AliAnalysisTaskFlowModes_pPb_H

#include <vector>

#include "AliAnalysisTaskSE.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliESDpid.h"
#include "AliFlowBayesianPID.h"

#include "TComplex.h"
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"



class AliAnalysisTaskFlowModes_pPb : public AliAnalysisTaskSE
{
    public:
    enum    RunMode {kTest, kFillWeights, kFull}; // task running mode (NOT GRID MODE)
    enum    AnalType {kAOD, kESD}; // tag for analysis type
    enum    PartSpecies {kUnknown, kCharged, kPion, kKaon, kProton}; // list of all particle species of interest
    enum    ColSystem {kPP, kPbPb, kPPb};
    struct FlowPart // representation of selected particle (species independent) storing only basic properties for flow calculations
    {
                FlowPart(Double_t dPt = 0, Double_t dPhi = 0, Double_t dEta = 0, Short_t iCharge = 0, PartSpecies sSpecies = kUnknown, Double_t dMass = 0, Double_t dPx = 0, Double_t dPy = 0, Double_t dPz = 0) :
                  pt(dPt), px(dPx), py(dPy), pz(dPz), phi(dPhi), eta(dEta), mass(dMass), charge(iCharge), species(sSpecies) {} // constructor

        void    PrintPart() const { printf("pt %g | px %g | py %g| pz %g | phi %g | eta %g | mass %g | charge %d | species %d \n",pt,px,py,pz,phi,eta,mass,charge,species); } // print struct members

        Double_t pt,px,py,pz,phi,eta,mass;
        Short_t charge;
        PartSpecies species;
    };
                            AliAnalysisTaskFlowModes_pPb(); // constructor
                            AliAnalysisTaskFlowModes_pPb(const char *name); // named (primary) constructor
    virtual                 ~AliAnalysisTaskFlowModes_pPb(); // destructor

    virtual void            UserCreateOutputObjects(); //
    virtual void            UserExec(Option_t* option); // main methond - called for each event
    virtual void            Terminate(Option_t* option); // called after all events are processed
    //pPb setter
    void CheckPileUp(){fCheckPileUp = kTRUE;}
    void UseSPDPileUpCuts(){fUsePileUpSPD = kTRUE;}
    void SetPileUpSPDParams(Int_t minVtxPileUpContrSPD, Float_t minPileUpZdistSPD)
    {
        fModifySPDDefaultParams = kTRUE;
        fMinVtxPileUpContrSPD = minVtxPileUpContrSPD;
        fMinPileUpZdistSPD = minPileUpZdistSPD;
    }
    void SetRejectOutOfBunchPileUp(){fRejectOutOfBunchPileUp = kTRUE;}
    // analysis setters
    void                    SetColisionSystem(ColSystem colSystem = kPPb) {fColSystem = colSystem; }
    void                    SetRunMode(RunMode mode = kFull) { fRunMode = mode; }
    void                    SetNumEventsAnalyse(Short_t num) { fNumEventsAnalyse = num; }
    void		            SetAnalysisType(AnalType type = kAOD) { fAnalType = type; }
    void                    SetSampling(Bool_t sample = kTRUE) { fSampling = sample; }
    void                    SetFillQAhistos(Bool_t fill = kTRUE) { fFillQA = fill; }
    void                    SetProcessCharged(Bool_t filter = kTRUE) { fProcessCharged = filter; }
    void                    SetProcessPID(Bool_t filter = kTRUE, Bool_t PIDbayesian = kFALSE) {
                                fProcessPID = filter;
                                if(PIDbayesian){fPIDbayesian = PIDbayesian; fPIDnsigma = kFALSE;}else{fPIDbayesian = kFALSE; fPIDnsigma = kTRUE;}
                              }
    // flow related setters
    void                    SetFlowRFPsPtMin(Float_t pt) { fCutFlowRFPsPtMin = pt; }
    void                    SetFlowRFPsPtMax(Float_t pt) { fCutFlowRFPsPtMax = pt; }
    void                    SetFlowDoFourCorrelations(Bool_t four = kTRUE) { fCutFlowDoFourCorrelations = four; }
    void                    SetFlowDoOnlyMixedCorrelations(Bool_t b = kFALSE) { fDoOnlyMixedCorrelations = b; }
    void                    SetFlowFillWeights(Bool_t weights = kTRUE) { fFlowFillWeights = weights; }
    void                    SetUseNUAWeigthsFile(const char* file) { fFlowNUAWeightsPath = file; if(fFlowNUAWeightsPath){fFlowUseNUAWeights = kTRUE;} }
    void                    SetUseNUEWeigthsFile(const char* file) { fFlowNUEWeightsPath = file; if(fFlowNUEWeightsPath){fFlowUseNUEWeights = kTRUE;} }

    // events setters
    void                    SetMultEstimator(const char* mult = "CHARGED") { fMultEstimator = mult; }
    void                    SetTrigger(Short_t trigger = 0) { fTrigger = trigger; }
    void		              SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
    void                    SetCentralityRange(Bool_t kCentRange){fFullCentralityRange = kCentRange;}
    // track setters
    void                    SetChargedEtaMax(Double_t eta) { fCutChargedEtaMax = eta; }
    void                    SetChargedPtMax(Double_t pt) { fCutChargedPtMax = pt; }
    void                    SetChargedPtMin(Double_t pt) { fCutChargedPtMin = pt; }
    void                    SetChargedDCAzMax(Double_t dcaz) {  fCutChargedDCAzMax = dcaz; }
    void                    SetChargedDCAxyMax(Double_t dcaxy) {  fCutChargedDCAxyMax = dcaxy; }
    void                    SetChargedNumTPCclsMin(UShort_t tpcCls) { fCutChargedNumTPCclsMin = tpcCls; }
    void		              SetMaxChi2perTPCcls(Double_t MaxChi2 = 4){ fMaxChi2perTPCcls = MaxChi2;}
    void                    SetChargedTrackFilterBit(UInt_t filter) { fCutChargedTrackFilterBit = filter; }
    // PID (pi,K,p) setters
    
    void                    SetPIDUseAntiProtonOnly(Bool_t use = kTRUE) { fCutPIDUseAntiProtonOnly = use; }
    void                    SetPIDNumSigmasPionMax(Double_t numSigmas) { fCutPIDnSigmaPionMax = numSigmas; }
    void                    SetPIDNumSigmasKaonMax(Double_t numSigmas) { fCutPIDnSigmaKaonMax = numSigmas; }
    void                    SetPIDNumSigmasProtonMax(Double_t numSigmas) { fCutPIDnSigmaProtonMax = numSigmas; }
    void                    SetPIDNumSigmasCombinedNoTOFrejection(Bool_t reject = kTRUE) { fCutPIDnSigmaCombinedNoTOFrejection = reject; }
    void                    SetPIDnsigmaCombination(Int_t Comb =2){fPIDnsigmaCombination = Comb;}
    void		    SetExtraPileUpCut(){fExtraPileUp = kTRUE;}
    void                    SetPositivelyChargedRef(Bool_t Pos=kFALSE){fPositivelyChargedRef = Pos;}
    void                    SetNegativelyChargedRef(Bool_t Neg=kFALSE){fNegativelyChargedRef = Neg;}
    void                    SetPositivelyChargedPOI(Bool_t Pos=kFALSE){fPositivelyChargedPOI = Pos;}
    void                    SetNegativelyChargedPOI(Bool_t Neg=kFALSE){fNegativelyChargedPOI = Neg;}
    void                    SetBayesianProbability(Double_t prob=0.9){fParticleProbability = prob;}
    void                    SetPriors(Float_t centr = 0); // set Noferini's favourite priors for Bayesian PID (requested if Bayesian PID is used)
    AliESDpid&              GetESDpid() {return fESDpid;}
    Bool_t                  TPCTOFagree(const AliVTrack *track);
    
   private:
    
    //pPb
    Bool_t                  fCheckPileUp; 
    Bool_t                  fUsePileUpSPD;
    Bool_t                  fModifySPDDefaultParams;
    Int_t                   fMinVtxPileUpContrSPD;
    Float_t                 fMinPileUpZdistSPD;
    Bool_t                  fRejectOutOfBunchPileUp;
     // properties
     AliAODEvent*            fEventAOD; //! AOD event countainer
     AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
     AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
     TFile*                  fFlowNUAWeightsFile; //! source file containing weights
     TFile*                  fFlowNUEWeightsFile; //! source file containing weights
     Bool_t                  fInit; // initialization check
     Short_t                 fIndexSampling; // sampling index (randomly generated)
     Short_t                 fIndexCentrality; // centrality bin index (based on centrality est. or number of selected tracks)
     Short_t                 fEventCounter; // event counter (used for local test runmode purpose)
     Short_t                 fNumEventsAnalyse; // [50] number of events to be analysed / after passing selection (only in test mode)
     Int_t                   fRunNumber; // [-1] run number obtained from AliVHeader
     Bool_t                  fExtraPileUp; // extra pile-up cuts 
     // array lenghts & constants
     const Double_t          fPDGMassPion; // [DPGMass] DPG mass of charged pion
     const Double_t          fPDGMassKaon; // [DPGMass] DPG mass of charged kaon
     const Double_t          fPDGMassProton; // [DPGMass] DPG mass of proton
    
     // selected POIs containers
     std::vector<FlowPart>*  fVectorCharged; //! container for selected charged particles
     std::vector<FlowPart>*  fVectorPion; //! container for selected pion candidates
     std::vector<FlowPart>*  fVectorKaon; //! container for selected kaon candidates
     std::vector<FlowPart>*  fVectorProton; //! container for selected proton candidates
    
    //cuts & selection: analysis
    RunMode                 fRunMode; // running mode (not grid related)
    AnalType                fAnalType; // analysis type: AOD / ESD
    Bool_t                  fSampling;      // Do random sampling ? (estimation of vn stat. uncertanity)
    Bool_t                  fFillQA; //[kTRUE] flag for filling the QA plots
    Bool_t                  fProcessCharged; // flag for processing charged tracks (both RPF and POIs)
    Bool_t                  fProcessPID; // flag for processing PID tracks (pi,K,p)
    AliFlowBayesianPID*     fBayesianResponse; //! Baysian response with all the TOF tuning (using fESDpid)

    // cuts & selection: flow related
    Float_t                 fCutFlowRFPsPtMin; // [0] (GeV/c) min pT threshold for RFPs particle for reference flow
    Float_t                 fCutFlowRFPsPtMax; // [0] (GeV/c) max pT threshold for RFPs particle for reference flow
    const Float_t           fFlowPOIsPtMin; // [0] (GeV/c) min pT treshold for POIs for differential flow
    const Float_t           fFlowPOIsPtMax; // [20] (GeV/c) max pT treshold for POIs for differential flow
    const Int_t             fFlowCentMin; // [0] min range for centrality/multiplicity histos
    const Int_t             fFlowCentMax; // [150] min range for centrality/multiplicity histos
    const Int_t             fFlowCentNumBins; // [150] min range for centrality/multiplicity histos
    Bool_t                  fCutFlowDoFourCorrelations; // [kTRUE] flag for processing <4>
    Bool_t                  fDoOnlyMixedCorrelations; // [kTRUE] flag if I only want to analyse mixed harmonics
    Bool_t                  fFlowFillWeights; //[kFALSE] flag for filling weights
    Bool_t                  fFlowUseNUAWeights; //[kFALSE] flag for using the previously filled NUA weights (NOTE: this is turned on only when path to file is applied via fFlowWeightsPath)
    Bool_t                  fFlowUseNUEWeights; //[kFALSE] flag for using the previously filled NUE weights (NOTE: this is turned on only when path to file is applied via fFlowWeightsPath)
    TString                 fFlowNUAWeightsPath;//[] path to source root file with weigthts (if empty unit weights are applied) "alice/cern.ch/user/n/nmohamma/CorrectionMaps/fb96/NUACorrectionMap.root"
   
    TString                 fFlowNUEWeightsPath; //[] path to source root file with weigthts (if empty unit weights are applied) e.g. "alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_LHC16kl.root"
    Bool_t                  fPositivelyChargedRef; //for same charged reference particle studies
    Bool_t                  fNegativelyChargedRef; //for same charged reference particle studies
    Bool_t                  fPositivelyChargedPOI; //for like sign reference particles and POIs studies
    Bool_t                  fNegativelyChargedPOI; //for unlike sign reference particle and POIs studies
    
    
    //cuts & selection: events
    Float_t                    fPVtxCutZ; // (cm) PV z cut
    ColSystem                  fColSystem; // collision system: pp or PbPb both with Run2 data (2016 data for pp and 2015 for PbPb, 2016 data for pPb)
    TString                 fMultEstimator; // [''] multiplicity estimator (suported: ''/Charged,VOA,V0C,V0M,CL0,CL1,ZNA,ZNC)
    Short_t                 fTrigger; // physics selection trigger
    Bool_t                  fFullCentralityRange; // flag for running over the full centrality range (0-100%). Otherwise runs from 50-100%
    
    //cuts & selection: tracks
    Float_t                 fCutChargedEtaMax; // (-) Maximum pseudorapidity range
    Float_t                 fCutChargedPtMax; // (GeV/c) Maximal track pT
    Float_t                 fCutChargedPtMin; // (GeV/c) Minimal track pT
    Float_t                 fCutChargedDCAzMax; // (cm) Maximal DCA-z cuts for tracks (pile-up rejection suggested for LHC16)
    Float_t                 fCutChargedDCAxyMax; // (cm) Maximal DCA-xy cuts for tracks (pile-up rejection suggested for LHC16)
    UInt_t                  fCutChargedTrackFilterBit; // (-) tracks filter bit
    UShort_t                fCutChargedNumTPCclsMin;  // (-) Minimal number of TPC clusters used for track reconstruction
    Double_t                fMaxChi2perTPCcls; // max chi2 per TPC clusters
    
    // cuts & selection: PID selection
    AliESDpid               fESDpid; //! pid obj
    Bool_t                  fCutPIDUseAntiProtonOnly; // [kFALSE] check proton PID charge to select AntiProtons only
    Bool_t                  fPIDnsigma; // [kTRUE] default pid method
    Int_t                   fPIDnsigmaCombination; // 1,2,3 for PID nsigma combinations 1,2,3
    Bool_t                  fPIDbayesian; // [kFALSE] bayesian pid method
    Double_t                fCutPIDnSigmaPionMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for pion candidates
    Double_t                fCutPIDnSigmaKaonMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for kaon candidates
    Double_t                fCutPIDnSigmaProtonMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for proton candidates
    Double_t                fCutPIDnSigmaTPCRejectElectron; // [3] number of TPC nSigma for electron rejection
    Bool_t                  fCutPIDnSigmaCombinedNoTOFrejection; // [kFALSE] flag for rejection candidates in TPC+TOF pt region if TOF is not available (if true and no TOF, only TPC is used)
    Float_t                 fCurrCentr; // current centrality used for set the priors
    Double_t                fParticleProbability; // Minimum Bayesian probability
    
    static const Int_t      fNumCentralities = 7; // number of centrality ranges for NUE maps (for now...)
    static const Int_t      fgkPIDptBin = 32; // pT bins for priors
    Float_t                 fC[fgkPIDptBin][5],fBinLimitPID[fgkPIDptBin]; // pt bin limit and priors
    static const Short_t    fFlowNumHarmonicsMax = 10; // maximum harmonics length of flow vector array
    static const Short_t    fFlowNumWeightPowersMax = 10; // maximum weight power length of flow vector array
    static const Short_t    fFlowPOIsPtNumBins = 100; // number of pT bins for POIs
    
    static const Short_t    fiNumIndexQA = 2; // QA indexes: 0: before cuts // 1: after cuts
   
    const static Short_t    fNumSamples = 10; // overall number of samples (from random sampling) used
    const static Int_t      fNumHarmonics = 5; // number of harmonics
    const static Int_t      fNumMixedHarmonics = 4; // number of mixed harmonics: 4{psi2}, 6{psi3} and 5{psi2,3}
    static Int_t            fHarmonics[fNumHarmonics]; // values of used harmonics
    static Int_t            fMixedHarmonics[fNumMixedHarmonics]; // values of used harmonics
    const static Int_t      fNumEtaGap = 3; // number of harmonics
    static Double_t         fEtaGap[fNumEtaGap]; // values of used harmonics

    
    
    // output lists
    TList*      fQAEvents; //! events list
    TList*      fQACharged; //! charged tracks list
    TList*      fQAPID; //! pi,K,p list
    TList*      fFlowWeights; //! list for flow weights
    TList*      fFlowRefs; //! list for flow of reference particles
    TList*      fFlowCharged; //! list for flow of charged particles
    TList*      fFlowPID; //! list for flow of PID (pi,K,p) particles
    
    // histograms & profiles
    
    // Flow
    TH3D*           fh3BeforeNUAWeightsRefs; //! distribution of Refs particles for estimating weight purpose (phi,eta,vtx_z)
    TH3D*           fh3BeforeNUAWeightsCharged; //! distribution of Charged POIs particles for estimating weight purpose (phi,eta,vtx_z)
    TH3D*           fh3BeforeNUAWeightsPion; //! distribution of Pion POIs particles for estimating weight purpose (phi,eta,vtx_z)
    TH3D*           fh3BeforeNUAWeightsKaon; //! distribution of Kaon POIs particles for estimating weight purpose (phi,eta,vtx_z)
    TH3D*           fh3BeforeNUAWeightsProton; //! distribution of Proton POIs particles for estimating weight purpose (phi,eta,vtx_z)
    
    TH3D*           fh3AfterNUAWeightsRefs; //! distribution of Refs particles after applying the weights (phi,eta,vtx_z)
    TH3D*           fh3AfterNUAWeightsCharged; //! distribution of Charged POIs particles after applying the weights (phi,eta,vtx_z)
    TH3D*           fh3AfterNUAWeightsPion; //! distribution of Pion POIs particles after applying the weights (phi,eta,vtx_z)
    TH3D*           fh3AfterNUAWeightsKaon; //! distribution of Kaon POIs particles after applying the weights (phi,eta,vtx_z)
    TH3D*           fh3AfterNUAWeightsProton; //! distribution of Proton POIs particles after applying the weights (phi,eta,vtx_z)
    
    TH1D*           fhBeforeNUEWeightsRefs; //! distribution of Refs particles for estimating weight purpose (pt)
    TH1D*           fhBeforeNUEWeightsCharged; //! distribution of Charged POIs particles for estimating weight purpose (pt)
    TH1D*           fhBeforeNUEWeightsPion; //! distribution of Pion POIs particles for estimating weight purpose (pt)
    TH1D*           fhBeforeNUEWeightsKaon; //! distribution of Kaon POIs particles for estimating weight purpose (pt)
    TH1D*           fhBeforeNUEWeightsProton; //! distribution of Proton POIs particles for estimating weight purpose (pt)
    
    TH1D*           fhAfterNUEWeightsRefs; //! distribution of Refs particles after applying the weights (pt)
    TH1D*           fhAfterNUEWeightsCharged; //! distribution of Charged POIs particles after applying the weights (pt)
    TH1D*           fhAfterNUEWeightsPion; //! distribution of Pion POIs particles after applying the weights (pt)
    TH1D*           fhAfterNUEWeightsKaon; //! distribution of Kaon POIs particles after applying the weights (pt)
    TH1D*           fhAfterNUEWeightsProton; //! distribution of Proton POIs particles after applying the weights (pt)
    
    TH3D*           fh3NUAWeightRefsPlus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightChargedPlus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightPionPlus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightKaonPlus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightProtonPlus; //! container for loading weights for given run
    
    TH3D*           fh3NUAWeightRefsMinus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightChargedMinus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightPionMinus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightKaonMinus; //! container for loading weights for given run
    TH3D*           fh3NUAWeightProtonMinus; //! container for loading weights for given run
    
    // Events
    TH2D*           fhEventSampling; //! distribution of sampled events (based on randomly generated numbers)
    TH1D*           fhEventCentrality; //! distribution of event centrality
    TH2D*           fh2EventCentralityNumSelCharged; //! distribution of event centrality vs number of selected charged tracks
    TH1D*           fhEventCounter; //! counter following event selection
    // Charged
    TH2D*           fh2RefsMult; //!multiplicity distribution of selected RFPs
    TH2D*           fh2RefsPt; //! pt distribution of selected RFPs
    TH2D*           fh2RefsEta; //! pt distribution of selected RFPs
    TH2D*           fh2RefsPhi; //! pt distribution of selected RFPs
    TH1D*           fhChargedCounter; //! counter following charged track selection
    // PID
    TH2D*           fh2PIDPionMult; //! multiplicity distribution of selected pions
    TH2D*           fh2PIDPionPt; //! pt distribution of selected pions
    TH2D*           fh2PIDPionPhi; //! phi distribution of selected pions
    TH2D*           fh2PIDPionEta; //! eta distribution of selected pions
    TH1D*           fhPIDPionCharge; //! charge distribution of selected pions
    
    TH2D*           fh2PIDKaonMult; //! multiplicity distribution of selected pions
    TH2D*           fh2PIDKaonPt; //! pt distribution of selected kaons
    TH2D*           fh2PIDKaonPhi; //! phi distribution of selected kaons
    TH2D*           fh2PIDKaonEta; //! eta distribution of selected kaons
    TH1D*           fhPIDKaonCharge; //! charge distribution of selected pions
    
    TH2D*           fh2PIDProtonMult; //! multiplicity distribution of selected pions
    TH2D*           fh2PIDProtonPt; //! pt distribution of selected protons
    TH2D*           fh2PIDProtonPhi; //! phi distribution of selected protons
    TH2D*           fh2PIDProtonEta; //! eta distribution of selected protons
    TH1D*           fhPIDProtonCharge; //! charge distribution of selected pions
    
    TH2D*           fh2PIDPionTPCdEdx; //! TPC dEdx response of selected pions
    TH2D*           fh2PIDPionTOFbeta; //! TOF beta of selected pions
    
    TH2D*           fh2PIDKaonTPCdEdx; //! TPC dEdx response of selected pions
    TH2D*           fh2PIDKaonTOFbeta; //! TOF beta of selected pions
    
    TH2D*           fh2PIDProtonTPCdEdx; //! TPC dEdx response of selected pions
    TH2D*           fh2PIDProtonTOFbeta; //! TOF beta of selected pions
    
    TH2D*           fh2PIDPionTPCnSigmaPion; //! TPC nSigma vs pT for selected pions (pion hypothesis)
    TH2D*           fh2PIDPionTOFnSigmaPion; //! TOF nSigma vs pT for selected pions (pion hypothesis)
    TH2D*           fh2PIDPionTPCnSigmaKaon; //! TPC nSigma vs pT for selected pions (kaon hypothesis)
    TH2D*           fh2PIDPionTOFnSigmaKaon; //! TOF nSigma vs pT for selected pions (kaon hypothesis)
    TH2D*           fh2PIDPionTPCnSigmaProton; //! TPC nSigma vs pT for selected pions (proton hypothesis)
    TH2D*           fh2PIDPionTOFnSigmaProton; //! TOF nSigma vs pT for selected pions (proton hypothesis)
    TH2D*           fh2PIDKaonTPCnSigmaPion; //! TPC nSigma vs pT for selected kaons (pion hypothesis)
    TH2D*           fh2PIDKaonTOFnSigmaPion; //! TOF nSigma vs pT for selected kaons (pion hypothesis)
    TH2D*           fh2PIDKaonTPCnSigmaKaon; //! TPC nSigma vs pT for selected kaons (kaon hypothesis)
    TH2D*           fh2PIDKaonTOFnSigmaKaon; //! TOF nSigma vs pT for selected kaons (kaon hypothesis)
    TH2D*           fh2PIDKaonTPCnSigmaProton; //! TPC nSigma vs pT for selected kaons (proton hypothesis)
    TH2D*           fh2PIDKaonTOFnSigmaProton; //! TOF nSigma vs pT for selected kaons (proton hypothesis)
    TH2D*           fh2PIDProtonTPCnSigmaPion; //! TPC nSigma vs pT for selected protons (pion hypothesis)
    TH2D*           fh2PIDProtonTOFnSigmaPion; //! TOF nSigma vs pT for selected protons (pion hypothesis)
    TH2D*           fh2PIDProtonTPCnSigmaKaon; //! TPC nSigma vs pT for selected protons (kaon hypothesis)
    TH2D*           fh2PIDProtonTOFnSigmaKaon; //! TOF nSigma vs pT for selected protons (kaon hypothesis)
    TH2D*           fh2PIDProtonTPCnSigmaProton; //! TPC nSigma vs pT for selected protons (proton hypothesis)
    TH2D*           fh2PIDProtonTOFnSigmaProton; //! TOF nSigma vs pT for selected protons (proton hypothesis)
    
    TH1F*           fhNUEWeightRefsPlus; //! containers for loading weights for given run
    TH1F*           fhNUEWeightRefsMinus; //! container for loading weights for given run

    TH1F*           fhNUEWeightChargedPlus; //! container for loading weights for given run
    TH1F*           fhNUEWeightChargedMinus; //! container for loading weights for given run

    TH1F*           fhNUEWeightPionPlus; //! container for loading weights for given run
    TH1F*           fhNUEWeightPionMinus; //! container for loading weights for given run

    TH1F*           fhNUEWeightKaonPlus; //! container for loading weights for given run
    TH1F*           fhNUEWeightKaonMinus; //! container for loading weights for given run

    TH1F*           fhNUEWeightProtonPlus; //! container for loading weights for given run
    TH1F*           fhNUEWeightProtonMinus; //! container for loading weights for given run

    TComplex                fFlowVecQpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
    TComplex                fFlowVecQneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
    TComplex                fFlowVecPpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]; // flow vector array for flow calculation
    TComplex                fFlowVecPneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]; // flow vector array for flow calculation
    TComplex                fFlowVecS[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]; // flow vector array for flow calculation
    

    // histograms & profiles

    TProfile*       fpMeanQxRefsPos[fNumEtaGap][fNumHarmonics]; //! average of Qx (vs. centrality) for Refs
    TProfile*       fpMeanQxRefsNeg[fNumEtaGap][fNumHarmonics]; //! average of Qx (vs. centrality) for Refs
    TProfile*       fpMeanQyRefsPos[fNumEtaGap][fNumHarmonics]; //! average of Qy (vs. centrality) for Refs
    TProfile*       fpMeanQyRefsNeg[fNumEtaGap][fNumHarmonics]; //! average of Qy (vs. centrality) for Refs

    TProfile*       fpRefsCor2[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2> correlations for RFPs
    TProfile2D*     fp2ChargedCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for Charged tracks POIs: POIs in Eta>0
    TProfile2D*     fp2ChargedCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for Charged tracks POIs: POIs in Eta<0
    TProfile2D*     fp2PionCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for pion POIs: POIs in Eta>0
    TProfile2D*     fp2PionCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for pion POIs: POIs in Eta>0
    TProfile2D*     fp2KaonCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for kaon POIs: POIs in Eta>0
    TProfile2D*     fp2KaonCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for kaon POIs: POIs in Eta>0
    TProfile2D*     fp2ProtonCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for proton POIs: POIs in Eta>0
    TProfile2D*     fp2ProtonCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for proton POIs: POIs in Eta>0

    //TProfile*       fpRefsCor4[fNumSamples][fNumHarmonics]; //! <4> correlations for RFPs
    //TProfile2D*     fp2ChargedCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for Charged tracks POIs
    //TProfile2D*     fp2PionCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for pion POIs
    //TProfile2D*     fp2KaonCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for kaon POIs
    //TProfile2D*     fp2ProtonCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for proton POIs
    
    //Mixed harmonics:
    TProfile*       fpMixedRefsCor4[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <4> correlations for RFPs
    TProfile*       fpMixedRefsCor6[fNumSamples][fNumEtaGap]; //! <6> correlations for RFPs
    TProfile2D*     fpMixedChargedCor3Pos[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for Charged tracks POIs: POIs in Eta>0
    TProfile2D*     fpMixedChargedCor3Neg[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for Charged tracks POIs: POIs in Eta<0
    TProfile2D*     fpMixedPionCor3Pos[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for pion POIs: POIs in Eta>0
    TProfile2D*     fpMixedPionCor3Neg[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for pion POIs: POIs in Eta<0
    TProfile2D*     fpMixedKaonCor3Pos[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for kaon POIs: POIs in Eta>0
    TProfile2D*     fpMixedKaonCor3Neg[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for kaon POIs: POIs in Eta<0
    TProfile2D*     fpMixedProtonCor3Pos[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for proton POIs: POIs in Eta>0
    TProfile2D*     fpMixedProtonCor3Neg[fNumSamples][fNumEtaGap][fNumMixedHarmonics]; //! <3'> correlations for proton POIs: POIs in Eta<0
    TProfile2D*     fpMixedChargedCor4Pos[fNumSamples][fNumEtaGap]; //! <4'> correlations for Charged tracks POIs: POIs in Eta>0
    TProfile2D*     fpMixedChargedCor4Neg[fNumSamples][fNumEtaGap]; //! <4'> correlations for Charged tracks POIs: POIs in Eta<0
    TProfile2D*     fpMixedPionCor4Pos[fNumSamples][fNumEtaGap]; //! <4'> correlations for pion POIs: POIs in Eta>0
    TProfile2D*     fpMixedPionCor4Neg[fNumSamples][fNumEtaGap]; //! <4'> correlations for pion POIs: POIs in Eta<0
    TProfile2D*     fpMixedKaonCor4Pos[fNumSamples][fNumEtaGap]; //! <4'> correlations for kaon POIs: POIs in Eta>0
    TProfile2D*     fpMixedKaonCor4Neg[fNumSamples][fNumEtaGap]; //! <4'> correlations for kaon POIs: POIs in Eta<0
    TProfile2D*     fpMixedProtonCor4Pos[fNumSamples][fNumEtaGap]; //! <4'> correlations for proton POIs: POIs in Eta>0
    TProfile2D*     fpMixedProtonCor4Neg[fNumSamples][fNumEtaGap]; //! <4'> correlations for proton POIs: POIs in Eta<0
    
    // QA: events
    TH1D*           fhQAEventsPVz[fiNumIndexQA]; //!
    TH1D*           fhQAEventsNumContrPV[fiNumIndexQA]; //!
    TH1D*           fhQAEventsNumSPDContrPV[fiNumIndexQA]; //!
    TH1D*           fhQAEventsDistPVSPD[fiNumIndexQA]; //!
    TH1D*           fhQAEventsSPDresol[fiNumIndexQA]; //!
    TH2D*           fhQAEventsCentralityOutliers[fiNumIndexQA]; //!
    TH2D*           fhQAEventsPileUp[fiNumIndexQA]; //!
    TH2D*           fhEventsMultTOFFilterbit32[fiNumIndexQA]; //!
    // QA: charged tracks
    TH1D*           fhQAChargedMult[fiNumIndexQA];       //! number of AOD charged tracks distribution
    TH1D*           fhQAChargedPt[fiNumIndexQA];         //! pT dist of charged tracks
    TH1D*           fhQAChargedEta[fiNumIndexQA];        //! eta dist of charged tracks
    TH1D*           fhQAChargedPhi[fiNumIndexQA];        //! phi dist of charged tracks
    TH1D*           fhQAChargedCharge[fiNumIndexQA];     //! charge dist of charged tracks
    TH1D*           fhQAChargedFilterBit[fiNumIndexQA];  //! filter bit distribution of charged tracks
    TH1D*           fhQAChargedNumTPCcls[fiNumIndexQA];  //! dist of track number of TPC clusters
    TH1D*           fhQAChargedDCAxy[fiNumIndexQA];      //! dist of Charged DCA in transverse plane
    TH1D*           fhQAChargedDCAz[fiNumIndexQA];       //! dist of charged DCA in z coordinate
    // QA: PID tracks
    TH1D*           fhQAPIDTPCstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
    TH1D*           fhQAPIDTOFstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
    TH2D*           fhQAPIDTPCdEdx[fiNumIndexQA];    //! TPC PID information
    TH2D*           fhQAPIDTOFbeta[fiNumIndexQA];    //! TOF PID information
  
    TH3D*           fh3PIDPionTPCTOFnSigmaPion[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected pions (pion hypothesis)
    TH3D*           fh3PIDPionTPCTOFnSigmaKaon[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected pions (kaon hypothesis)
    TH3D*           fh3PIDPionTPCTOFnSigmaProton[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected pions (proton hypothesis)
    
    TH3D*           fh3PIDKaonTPCTOFnSigmaPion[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected kaons (pion hypothesis)
    TH3D*           fh3PIDKaonTPCTOFnSigmaKaon[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected kaons (kaon hypothesis)
    TH3D*           fh3PIDKaonTPCTOFnSigmaProton[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected kaons (proton hypothesis)

    TH3D*           fh3PIDProtonTPCTOFnSigmaPion[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected protons (pion hypothesis)
    TH3D*           fh3PIDProtonTPCTOFnSigmaKaon[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected protons (kaon hypothesis)
    TH3D*           fh3PIDProtonTPCTOFnSigmaProton[fiNumIndexQA]; //! TPC nSigma vs TOF nSigma vs pT for selected protons (proton hypothesis)

    
    // functions:
    Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
    void                    ListParameters(); // list all task parameters
    
    Bool_t                  EventSelection(); // main method for event selection (specific event selection is applied within)
    Bool_t                  IsEventSelected_PbPb(); // event selection for LHC2015 PbPb data
    Bool_t                  IsEventSelected_pp(); // event selection for LHC2016 MB pp data
    Bool_t                  IsEventSelected_pPb(); // event selection for LHC2016 pPb data
    void                    FillEventsQA(const Short_t iQAindex); // filling QA plots related to event selection
    Short_t                 GetSamplingIndex(); // returns sampling index based on sampling selection (number of samples)
    Short_t                 GetCentralityIndex(); // returns centrality index based centrality estimator or number of selected tracks
    Double_t                GetWDist(const AliAODVertex* v0, const AliAODVertex* v1); // gets the distance between the two vertices
    Bool_t                  ProcessEvent(); // main (envelope) method for processing events passing selection
    
    Bool_t                  Filtering(); // main (envelope) method for filtering all POIs in event
    void                    FilterCharged(); // charged tracks filtering
    Bool_t                  IsChargedSelected(const AliAODTrack* track = 0x0); // charged track selection
    void                    FillQARefs(const Short_t iQAindex, const AliAODTrack* track = 0x0); // filling QA plots for RFPs selection
    void                    FillQACharged(const Short_t iQAindex, const AliAODTrack* track = 0x0); // filling QA plots for charged track selection
    void                    FilterPID(); // pi,K,p filtering
    PartSpecies             IsPIDSelected(const AliAODTrack* track); // PID tracks selections
    void                    FillPIDQA(const Short_t iQAindex, const AliAODTrack* track = 0x0, const PartSpecies species = kUnknown); // filling pi,K,p QA histograms
    // Flow related methods
    void                    DoFlowRefs(const Short_t iEtaGapIndex = 0); // Estimate <2> for reference flow
    void                    DoFlowCharged(const Short_t iEtaGapIndex = 0); // Estimate <2'> for pt diff. flow of charged hadrons
    void                    DoFlowPID(const Short_t iEtaGapIndex = 0, const PartSpecies species = kUnknown); // Estimate <2'> for pt diff. flow of PID (pi,K,p) hadrons
    void                    FillRefsVectors(const Short_t iEtaGapIndex = 0); // fill flow vector Q with RFPs for reference flow
    void                    FillPOIsVectors(const Short_t iEtaGapIndex = 0, const PartSpecies species = kUnknown, const Short_t iMassIndex = 0); // fill flow vectors p,q and s with POIs (for given species) for differential flow calculations
    Short_t                 GetPOIsPtBinIndex(const Double_t pt); // return pT bin index based on momenta value
    void                    ResetRFPsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]); // set values to TComplex(0,0,0) for given array
    void                    ResetPOIsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]); // set values to TComplex(0,0,0) for given array
    void                    ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]); // printf all values of given Flow vector array
    
    TComplex                Q(const Short_t n, const Short_t p);
    TComplex                QGapPos(const Short_t n, const Short_t p);
    TComplex                QGapNeg(const Short_t n, const Short_t p);
    TComplex                P(const Short_t n, const Short_t p, const Short_t pt);
    TComplex                PGapPos(const Short_t n, const Short_t p, const Short_t pt);
    TComplex                PGapNeg(const Short_t n, const Short_t p, const Short_t pt);
    TComplex                S(const Short_t n, const Short_t p, const Short_t pt);
    
    TComplex                Two(const Short_t n1, const Short_t n2); // Two particle reference correlation calculations (no eta gap)
    TComplex                TwoGap(const Short_t n1, const Short_t n2); // Two particle reference correlation calculations (with eta gap)
    TComplex                TwoDiff(const Short_t n1, const Short_t n2, const Short_t pt); // Two particle diff. correlation calculations (no eta gap)
    TComplex                TwoDiffGapPos(const Short_t n1, const Short_t n2, const Short_t pt); // Two particle diff. correlation calculations (with eta gap)
    TComplex                TwoDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t pt); // Two particle diff. correlation calculations (with eta gap)
    
    TComplex                ThreeDiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3,const Short_t pt); // Three particle diff. correlation calculations (with eta gap);
    TComplex                ThreeDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3,const Short_t pt); // Three particle diff. correlation calculations (with eta gap);
    
    TComplex                Four(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (no eta gap)
    TComplex                FourDiff(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t pt); // Four particle reference correlation calculations (no eta gap)
    TComplex                FourGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (with eta gap); n1+n2 = n3+n4;  n1, n2 from P & n3, n4 from M
    TComplex                FourGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (with eta gap); n1+n2 = n3+n4;  n1, n2 from M & n3, n4 from P
    TComplex                Four13DiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t pt); // Four particle diff. correlation calculations (with eta gap); n1 = n2+n3+n4;  n1 from P & n2, n3, n4 from M
    TComplex                Four13DiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t pt); // Four particle diff. correlation calculations (with eta gap); n1 = n2+n3+n4;  n1 from M & n2, n3, n4 from P
    TComplex                SixGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t n5, const Short_t n6); // Six particle reference correlation calculations (with eta gap); n1 + n2 + n3 = n4 + n5 + n6 ;n1, n2, n3 from P and n4, n5 and n6 from M
    TComplex                SixGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t n5, const Short_t n6); // Six particle reference correlation calculations (with eta gap); n1 + n2 + n3 = n4 + n5 + n6; n1, n2, n3 from M and n4, n5 and n6 from P


    AliAnalysisTaskFlowModes_pPb(const AliAnalysisTaskFlowModes_pPb&); // not implemented
    AliAnalysisTaskFlowModes_pPb& operator=(const AliAnalysisTaskFlowModes_pPb&); // not implemented

    ClassDef(AliAnalysisTaskFlowModes_pPb, 4);
};

#endif
