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


 private:

  Int_t ReadEventESD         (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventESDRecMcProp(std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventESDMC       (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  
  Int_t ReadEventAOD         (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,  
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventAODRecMcProp(std::vector<Float_t> &pt,  std::vector<Float_t> &eta,
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  Int_t ReadEventAODMC       (std::vector<Float_t> &pt,  std::vector<Float_t> &eta,  
			      std::vector<Float_t> &phi,  std::vector<Short_t> &charge,
			      std::vector<Int_t> &nTracksTracklets, const Int_t step);
  
  void  Analyse         (const std::vector<Float_t> &pt, 
			 const std::vector<Float_t> &eta, 
			 const std::vector<Float_t> &phi, 
			 const std::vector<Short_t> &charge, 
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

  AliESDEvent *fESDEvent;                   //! esd event
  AliAODEvent *fAODEvent;                   //! aod event
  Int_t        fNMcPrimAccept;              // global variable for mc multiplucity
  Int_t        fNRecAccept;                 // global variable for rec multiplucity
  Float_t      fVzEvent;                    // global variable for rec vertex position
  Double_t     fMeanPtRec;                  // global variable for rec mean pt
  Double_t     fLeadingPtRec;               // global variable for rec mean pt
  
  TList	     *fHists;                       // output list
  TH1F       *fStep;                        // how many events have passed which correction step
  TH1F       *fHistPt;                      // Pt spectrum ESD
  TH1F       *fHistPtMC;                    // Pt spectrum MC

  TH2F       *fNmcNch;                      // N mc - N ch rec
  TProfile   *fPNmcNch;                     // N mc - N ch rec
  TH2F       *fNmcNchVtx;                   // N mc - N ch rec for events with reconstructed vertex
  TProfile   *fPNmcNchVtx;                  // N mc - N ch rec for events with reconstructed vertex
  TH2F       *fChargedPi0;                  // charged versus charged+Pi0

  THnSparse   *fMapSingleTrig[6];           //! multi-dim histo for trigger track properties
  THnSparse   *fMapPair[6];                 //! multi-dim histo for pair properties
  THnSparse   *fMapEvent[6];                //! multi-dim histo for event properties
  THnSparse   *fMapAll[6];                  //! multi-dim histo for properties of all analysed tracks
  
  TH1F       * fVertexZ[6];                 // z of vertex
  TH1F       * fNcharge[6];                 // pt
  TH1F       * fPt[6];                      // pt
  TH1F       * fEta[6];                     // eta
  TH1F       * fPhi[6];                     // phi
  TH1F       * fDcaXY[6];                   // dca xy direction
  TH1F       * fDcaZ[6];                    // dca z direction

  TH1F       * fPtSeed[6];                  // pt of seed (event axis)
  TH1F       * fEtaSeed[6];                 // eta of seed 
  TH1F       * fPhiSeed[6];                 // phi of seed

  TH1F       * fPtOthers[6];                // pt of all other particels used in dEtadPhi
  TH1F       * fEtaOthers[6];               // eta of all other particels used in dEtadPhi
  TH1F       * fPhiOthers[6];               // phi of all other particels used in dEtadPhi
  TH2F       * fPtEtaOthers[6];             // pt-eta of all other particels used in dEtadPhi


  TH2F       * fPhiEta[6];                  // eta - phi
  TH2F       * fDPhiDEtaEventAxis[6];       // correlation dEta-dPhi towards event axis
  TH2F       * fDPhiDEtaEventAxisSeeds[6];  // correlation dEta-dPhi towards event axis of trigger particles
  TH1F       * fTriggerNch[6];              // number of triggers with accepted-track number
  TH2F       * fTriggerNchSeeds[6];         // number of triggers with accepted-track number
  TH1F       * fTriggerTracklet[6];         // number of triggers with accepted-tracklet number
  TH2F       * fNch07Nch[6];                // nCharged with pT>fTriggerPtCut vs nCharged
  TProfile   * fPNch07Nch[6];               // nCharged with pT>fTriggerPtCut vs nCharged
  
  TH2F       * fNch07Tracklet[6];           // nCharged with pT>fTriggerPtCut vs nTracklet
  TH2F       * fNchTracklet[6];             // nCharged vs nTracklet
  TProfile   * fPNch07Tracklet[6];           // nCharged with pT>fTriggerPtCut vs nTracklet

  TH1F       * fDPhiEventAxis[6];           // delta phi of associate tracks to event axis
  
  AliAnalysisTaskMinijet(const AliAnalysisTaskMinijet&); // not implemented
  AliAnalysisTaskMinijet& operator=(const AliAnalysisTaskMinijet&); // not implemented
  
  ClassDef(AliAnalysisTaskMinijet, 1); // example of analysis
};

#endif
