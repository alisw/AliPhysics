#ifndef ALIANALYSISTASKMINIJET_H
#define ALIANALYSISTASKMINIJET_H

// Two-particle correlations using all particles over pt threshold
// Extract mini-jet yield and fragmentation properties via Delta-Phi histograms
// Can use ESD or AOD, reconstructed and Monte Carlo data as input
// Author: eva.sicking@cern.ch

class TList;
class TH1F;
class TH2F;
class TProfile;
class THnSparse;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include <vector>

class AliAnalysisTaskMinijet : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMinijet(const char *name="<default name>");
  virtual ~AliAnalysisTaskMinijet();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t *);
  virtual void SetCuts(AliESDtrackCuts* cuts){fCuts = cuts;}
  
  void         SetUseMC(Bool_t useMC=kTRUE, Bool_t mcOnly=kFALSE)    {fUseMC = useMC; fMcOnly=mcOnly;}
  void         SetAnalyseOnlyPrimaries(Bool_t analysePrimOnly)       {fAnalysePrimOnly = analysePrimOnly;} // not used anymore
  void         SetPtRange(Float_t ptMin, Float_t ptMax)              {fPtMin = ptMin; fPtMax = ptMax; }
  void         SetTriggerPtCut(Float_t triggerPtCut)                 {fTriggerPtCut = triggerPtCut;}  
  void         SetAssociatePtCut(Float_t associatePtCut)             {fAssociatePtCut = associatePtCut;} 
  void         SetModeEsdAod(Int_t mode)                             {fMode = mode;}
  void         SetTriggerMask(Int_t triggerType)                     {fTriggerType = triggerType;}
  void         SetFilterBit(Int_t filterBit)                         {fFilterBit = filterBit;}
  void         SetMaxVertexZ(Float_t vertexZCut)                     {fVertexZCut = vertexZCut;}
  void         SetMaxEta(Float_t etaCut)                             {fEtaCut = etaCut;}
  void         SetMaxEtaSeed(Float_t etaCutSeed)                     {fEtaCutSeed = etaCutSeed;}
  void         SetSelectParticles(Int_t selectParticles)             {fSelectParticles = selectParticles;}
  void         SetSelectParticlesAssoc(Int_t selectParticlesAssoc)   {fSelectParticlesAssoc = selectParticlesAssoc;}
  void         SetCheckSDD(Bool_t checkSDD, Int_t selOption)         {fCheckSDD = checkSDD; fSelOption = selOption;}
  void         SetCorrStrangeness(Bool_t corrStrangeness)            {fCorrStrangeness = corrStrangeness;}
  void         SetThreeParticleCorrelation(Bool_t threeParticleCorr) {fThreeParticleCorr = threeParticleCorr;}
  void         SetRejectCorrupted(Bool_t rejectChunks, Int_t nTPC)   {fRejectChunks = rejectChunks; fNTPC = nTPC;}

    void         SetCentralityMethod(TString centralityMethod)                   {fCentralityMethod = centralityMethod;}
    
    
 private:

  Int_t ReadEventESD         (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Float_t> &strangnessWeight,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventESDRecMcProp(std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Float_t> &strangnessWeight,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventESDMC       (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Float_t> &strangnessWeight,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  
  Int_t ReadEventAOD         (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,  
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Float_t> &strangnessWeight,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventAODRecMcProp(std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Float_t> &strangnessWeight,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventAODMC       (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,  
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Float_t> &strangnessWeight,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  
  void  Analyse         (const std::vector<Float_t> &pt, 
			 const std::vector<Float_t> &eta, 
			 const std::vector<Float_t> &phi, 
			 const std::vector<Short_t> &charge, 
			 const std::vector<Float_t> &strangnessWeight,
			 const Int_t ntacks, const Int_t ntacklets=0, 
			 const Int_t nAll=0, const Int_t step=0);
  
  Bool_t                 SelectParticlePlusCharged(const Short_t charge, const Int_t pdg, const Bool_t prim);
  Bool_t                 SelectParticle(const Short_t charge, const Int_t pdg, const Bool_t prim);
  Bool_t                 CheckEvent(const Bool_t recVertex);
  const Double_t*        CreateLogAxis(const Int_t nbins, const Double_t xmin, const Double_t xmax); 
  Bool_t                 CheckLikeSign(const Short_t chargeEventAxis, const Short_t chargeOthers);


  Bool_t       fUseMC;                      // flag for Monte Carlo usages
  Bool_t       fMcOnly;                     // flag defines, if only MC data is used in analysis or also reconstructed data
  Double_t     fBSign;                      // magnetic field
  Bool_t       fAnalysePrimOnly;            // flag for analysis of primaries only (also in reconstructed data)
  Float_t      fPtMin;                      // set lower limit for pt acceptance for mutliplicity defintion
  Float_t      fPtMax;                      // set upper limit for pt acceptance for mutliplicity defintion
  AliESDtrackCuts* fCuts;                   // List of cuts for ESDs
  Float_t      fTriggerPtCut;               // cut on particle pt used as event axis
  Float_t      fAssociatePtCut;             // cut on particle pt used for correlations
  Int_t        fMode;                       // ESD(=0) of AOD(=1) reading 
  Int_t        fTriggerType;                // sets trigger -> AliVEvent::kMB, AliVEvent::kHighMult
  Int_t        fFilterBit;                  // Filter bit written in ESD filter, select track type
  Float_t      fVertexZCut;                 // vertex cut
  Float_t      fEtaCut;                     // eta acceptance cut
  Float_t      fEtaCutSeed;                 // eta acceptance cut for seed
  Int_t        fSelectParticles;            // only in cas of MC: use also neutral particles or not 
  Int_t        fSelectParticlesAssoc;       // only in cas of MC: use also neutral particles or not 
  Bool_t       fCheckSDD;                   // check if SDD was in read out partition (needed for LHC11a)
  Int_t        fSelOption;                  // 0 = use hit in SDD for event selection, 1 = use trigger for event selection
  Bool_t       fCorrStrangeness;            // for data correction -> Pythia simulations underestimate contamination from strangness
  Bool_t       fThreeParticleCorr;          // perform three particle correlation
  Bool_t       fRejectChunks;               // rejection of chunks in which no ITS tracks are reconstructed
  Int_t        fNTPC;                       // track number limit for rejection decision.

  AliESDEvent *fESDEvent;                   //! esd event
  AliAODEvent *fAODEvent;                   //! aod event
  Int_t        fNMcPrimAccept;              // global variable for mc multiplucity
  Int_t        fNRecAccept;                 // global variable for rec multiplucity
  Float_t      fNRecAcceptStrangeCorr;                 // global variable for rec multiplucity
  Int_t        fNMcPrimAcceptTracklet;      // global variable for mc multiplucity
  Int_t        fNRecAcceptTracklet;         // global variable for rec multiplucity
  Float_t      fVzEvent;                    // global variable for rec vertex position
  Double_t     fMeanPtRec;                  // global variable for rec mean pt
  Double_t     fLeadingPtRec;               // global variable for rec mean pt

  TList	     *fHists;                       // output list
  TH1F       *fStep;                        // how many events have passed which correction step
  TH1F       *fEventStat;                   // how many events are accepted by trigger, vertex selection, 1 track in acceptance (for real data)
  TH1F       *fHistPt;                      // Pt spectrum ESD
  TH1F       *fHistPtMC;                    // Pt spectrum MC
  TH2F       *fNContrNtracklets;            // control histogram for vertex->nContributers and number of tracklets
  TH2F       *fNContrNtracks;               // control histogram for vertex->nContributers and number of tracks
  TH2F       *fCorruptedChunks;             // control histogram: TPC tracks versus ITS-TPC-tracks
  TH2F       *fCorruptedChunksAfter;        // control histogram: TPC tracks versus ITS-TPC-tracks

  TH2F       *fNmcNch;                      // N mc - N ch rec
  TProfile   *fPNmcNch;                     // N mc - N ch rec
  TH2F       *fNmcNchVtx;                   // N mc - N ch rec for events with reconstructed vertex
  TH2F       *fNmcNchVtxStrangeCorr;        // N mc - N ch rec for events with reconstructed vertex + strangeness correction
  TProfile   *fPNmcNchVtx;                  // N mc - N ch rec for events with reconstructed vertex
  TH2F       *fNmcNchTracklet;              // N mc - N ch rec
  TProfile   *fPNmcNchTracklet;             // N mc - N ch rec
  TH2F       *fNmcNchVtxTracklet;           // N mc - N ch rec for events with reconstructed vertex
  TProfile   *fPNmcNchVtxTracklet;          // N mc - N ch rec for events with reconstructed vertex
  TH2F       *fChargedPi0;                  // charged versus charged+Pi0
  TH1F       *fVertexCheck;                 // check which fraction of events has vtx_rec but no good vtx_mc
  TH1F       *fPropagateDca;                // check of AliAODtrack::PropagateToDca

  THnSparse  *fMapSingleTrig[8];            //! multi-dim histo for trigger track properties
  THnSparse  *fMapPair[8];                  //! multi-dim histo for pair properties
  THnSparse  *fMapEvent[8];                 //! multi-dim histo for event properties
  THnSparse  *fMapAll[8];                   //! multi-dim histo for properties of all analysed tracks
  THnSparse  *fMapThree[8];                 //! multi-dim histo for properties of three particle correlations
  
  TH1F       * fVertexZ[8];                 // z of vertex
  TH1F       * fNcharge[8];                 // pt
  TH1F       * fPt[8];                      // pt
  TH1F       * fEta[8];                     // eta
  TH1F       * fPhi[8];                     // phi
  TH1F       * fDcaXY[8];                   // dca xy direction
  TH1F       * fDcaZ[8];                    // dca z direction

  TH1F       * fPtSeed[8];                  // pt of seed (event axis)
  TH1F       * fEtaSeed[8];                 // eta of seed 
  TH1F       * fPhiSeed[8];                 // phi of seed

  TH1F       * fPtOthers[8];                // pt of all other particels used in dEtadPhi
  TH1F       * fEtaOthers[8];               // eta of all other particels used in dEtadPhi
  TH1F       * fPhiOthers[8];               // phi of all other particels used in dEtadPhi
  TH2F       * fPtEtaOthers[8];             // pt-eta of all other particels used in dEtadPhi


  TH2F       * fPhiEta[8];                  // eta - phi
  TH2F       * fDPhiDEtaEventAxis[8];       // correlation dEta-dPhi towards event axis
  TH2F       * fDPhiDEtaEventAxisSeeds[8];  // correlation dEta-dPhi towards event axis of trigger particles
  TH1F       * fTriggerNch[8];              // number of triggers with accepted-track number
  TH2F       * fTriggerNchSeeds[8];         // number of triggers with accepted-track number
  TH1F       * fTriggerTracklet[8];         // number of triggers with accepted-tracklet number
  TH2F       * fNch07Nch[8];                // nCharged with pT>fTriggerPtCut vs nCharged
  TProfile   * fPNch07Nch[8];               // nCharged with pT>fTriggerPtCut vs nCharged
  
  TH2F       * fNch07Tracklet[8];           // nCharged with pT>fTriggerPtCut vs nTracklet
  TH2F       * fNchTracklet[8];             // nCharged vs nTracklet
  TProfile   * fPNch07Tracklet[8];           // nCharged with pT>fTriggerPtCut vs nTracklet

  TH1F       * fDPhiEventAxis[8];           // delta phi of associate tracks to event axis
  TH2F       * fDPhi1DPhi2[8];              // dPhi1 versus dPhi2: three particle correlation test
    
    TString fCentralityMethod;        //centrality pA
 
  AliAnalysisTaskMinijet(const AliAnalysisTaskMinijet&); // not implemented
  AliAnalysisTaskMinijet& operator=(const AliAnalysisTaskMinijet&); // not implemented
  
  ClassDef(AliAnalysisTaskMinijet, 1); // example of analysis
};

#endif
