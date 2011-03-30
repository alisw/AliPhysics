#ifndef ALIANALYSISTASKMINIJET_H
#define ALIANALYSISTASKMINIJET_H

// analysis task performing mini jet analysis (use ESD or AOD as input)
// Author: eva.sicking@cern.ch

class TList;
class TH1F;
class TH2F;
class TProfile;

class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskMinijet : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMinijet(const char *name="<default name>");
  virtual ~AliAnalysisTaskMinijet();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t* option);
  virtual void   Terminate(Option_t *);
  
  Int_t LoopESD  (Float_t **pt, Float_t **eta, Float_t **phi, Short_t **charge, Int_t **nTracksTracklets);
  Int_t LoopESDMC(Float_t **pt, Float_t **eta, Float_t **phi, Short_t **charge,  Int_t **nTracksTracklets);
  Int_t LoopAOD  (Float_t **pt, Float_t **eta, Float_t **phi, Short_t **charge,  Int_t **nTracksTracklets);
  Int_t LoopAODMC(Float_t **pt, Float_t **eta, Float_t **phi, Short_t **charge,  Int_t **nTracksTracklets);
  void  Analyse  (Float_t* pt, Float_t* eta, Float_t* phi,  Short_t *charge, Int_t ntacks, Int_t ntacklets=0,Int_t nAll=0, Int_t mode=0);
  void  CleanArrays(Float_t* pt, Float_t* eta, Float_t* phi,  Short_t *charge, Int_t* nTracksTracklets=0);
  Bool_t SelectParticlePlusCharged(Short_t charge, Int_t pdg, Bool_t prim);
  Bool_t SelectParticle(Short_t charge, Int_t pdg, Bool_t prim);


  void UseMC(Bool_t useMC=kTRUE, Bool_t mcOnly=kFALSE)    {fUseMC = useMC; fMcOnly=mcOnly;}

  virtual void   SetCuts(AliESDtrackCuts* cuts)           {fCuts = cuts;}

  void   SetRadiusCut(Float_t radiusCut)                  {fRadiusCut = radiusCut;}  
  void   SetTriggerPtCut(Float_t triggerPtCut)            {fTriggerPtCut = triggerPtCut;}  
  void   SetAssociatePtCut(Float_t associatePtCut)        {fAssociatePtCut = associatePtCut;} 
  void   SetEventAxis(Int_t leadingOrRandom)              {fLeadingOrRandom = leadingOrRandom;}  
  void   SetMode(Int_t mode)                              {fMode = mode;}
  void   SetMaxVertexZ(Float_t vertexZCut)                {fVertexZCut = vertexZCut;}
  void   SetMaxEta(Float_t etaCut)                        {fEtaCut = etaCut;}
  void   SetMaxEtaSeed(Float_t etaCutSeed)                {fEtaCutSeed = etaCutSeed;}

  void   SelectParticles(Int_t selectParticles)           {fSelectParticles = selectParticles;}
  void   SelectParticlesAssoc(Int_t selectParticlesAssoc) {fSelectParticlesAssoc = selectParticlesAssoc;}


 private:

  Bool_t       fUseMC;
  Bool_t       fMcOnly;
  AliESDtrackCuts* fCuts;                   // List of cuts for ESDs
  Float_t      fRadiusCut;                  // radius cut 
  Float_t      fTriggerPtCut;               // cut on particle pt used as event axis
  Float_t      fAssociatePtCut;             // cut on particle pt used for correlations
  Int_t        fLeadingOrRandom;            // event axis:leading track or random track
  Int_t        fMode;                       // ESD(=0) of AOD(=1) reading 
  Float_t      fVertexZCut;                 // vertex cut
  Float_t      fEtaCut;                     // eta acceptance cut
  Float_t      fEtaCutSeed;                 // eta acceptance cut for seed
  Int_t        fSelectParticles;            // only in cas of MC: use also neutral particles or not 
  Int_t        fSelectParticlesAssoc;       // only in cas of MC: use also neutral particles or not 

  AliESDEvent *fESDEvent;                   //! esd event
  AliAODEvent *fAODEvent;                   //! aod event
  Int_t        fNMcPrimAccept;              // global variable for mc multiplucity
  Float_t      fVzEvent;                    // global variable for rec vertex position
  
  TList	     *fHists;                      // output list
  TH1F       *fHistPt;                     // Pt spectrum ESD
  TH1F       *fHistPtMC;                   // Pt spectrum MC
  TH2F       *fNmcNch;                     // N mc - N ch rec
  TProfile   *pNmcNch;                     // N mc - N ch rec
  TH2F       *fChargedPi0;                 // charged versus charged+Pi0
  TH1F       * fVertexZ[4];                 // z of vertex

  TH1F       * fPt[4];                      // pt
  TH1F       * fEta[4];                     // et
  TH1F       * fPhi[4];                     // phi
  TH1F       * fDcaXY[4];                   // dca xy direction
  TH1F       * fDcaZ[4];                    // dca z direction

  TH1F       * fPtSeed[4];                  // pt of seed (event axis)
  TH1F       * fEtaSeed[4];                 // eta of seed 
  TH1F       * fPhiSeed[4];                 // phi of seed

  TH1F       * fPtOthers[4];                // pt of all other particels used in dEtadPhi
  TH1F       * fEtaOthers[4];               // eta of all other particels used in dEtadPhi
  TH1F       * fPhiOthers[4];               // phi of all other particels used in dEtadPhi
  TH2F       * fPtEtaOthers[4];             // pt-eta of all other particels used in dEtadPhi


  TH2F       * fPhiEta[4];                  // eta - phi
  TH2F       * fDPhiDEtaEventAxis[4];       // correlation dEta-dPhi towards event axis
  TH2F       * fDPhiDEtaEventAxisSeeds[4];  // correlation dEta-dPhi towards event axis of trigger particles
  TH1F       * fTriggerNch[4];              // number of triggers with accepted-track number
  TH2F       * fTriggerNchSeeds[4];         // number of triggers with accepted-track number
  TH1F       * fTriggerTracklet[4];         // number of triggers with accepted-tracklet number
  TH2F       * fNch07Nch[4];                // nCharged with pT>fTriggerPtCut vs nCharged
  TProfile   * pNch07Nch[4];                // nCharged with pT>fTriggerPtCut vs nCharged
  TH2F       * fNch07Tracklet[4];           // nCharged with pT>fTriggerPtCut vs nTracklet
  TH2F       * fNchTracklet[4];             // nCharged vs nTracklet
  TProfile   * pNch07Tracklet[4];           // nCharged with pT>fTriggerPtCut vs nTracklet

  TH1F       * fDPhiEventAxis[4];           // delta phi of associate tracks to event axis
  TH1F       * fDPhiEventAxisNchBin[4][150];// delta phi of associate tracks to event axis per Nch bin
  TH1F       * fDPhiEventAxisNchBinTrig[4][150];// "" for all possoble trigger particles

  TH1F       * fDPhiEventAxisTrackletBin[4][150]; // delta phi of associate tracks to event axis per Nch bin
  TH1F       * fDPhiEventAxisTrackletBinTrig[4][150]; // "" for all possible trigger particles

  AliAnalysisTaskMinijet(const AliAnalysisTaskMinijet&); // not implemented
  AliAnalysisTaskMinijet& operator=(const AliAnalysisTaskMinijet&); // not implemented
  
  ClassDef(AliAnalysisTaskMinijet, 1); // example of analysis
};

#endif
