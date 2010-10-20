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
  
  Int_t LoopESD  (Float_t **pt, Float_t **eta, Float_t **phi, Int_t **numbers);
  Int_t LoopESDMC(Float_t **pt, Float_t **eta, Float_t **phi);
  Int_t LoopAOD  (Float_t **pt, Float_t **eta, Float_t **phi, Int_t **numbers);
  Int_t LoopAODMC(Float_t **pt, Float_t **eta, Float_t **phi);
  void  Analyse  (Float_t* pt, Float_t* eta, Float_t* phi, Int_t ntacks, Int_t ntacklets=0,Int_t mode=0);
  void  CleanArrays(Float_t* pt, Float_t* eta, Float_t* phi, Int_t* numbers=0);

  void UseMC(Bool_t useMC=kTRUE) { fUseMC = useMC;}
  virtual void   SetCuts(AliESDtrackCuts* cuts)    {fCuts = cuts;}

  void   SetRadiusCut(Float_t radiusCut)           {fRadiusCut = radiusCut;}  
  void   SetTriggerPtCut(Float_t triggerPtCut)     {fTriggerPtCut = triggerPtCut;}  
  void   SetAssociatePtCut(Float_t associatePtCut) {fAssociatePtCut = associatePtCut;} 
  void   SetEventAxis(Int_t leadingOrRandom)       {fLeadingOrRandom = leadingOrRandom;}  
  void   SetMode(Int_t mode)                       {fMode=mode;}
  void   SetMaxVertexZ(Float_t vertexZCut)         {fVertexZCut=vertexZCut;}


 private:
  Bool_t       fUseMC;
  AliESDtrackCuts* fCuts;                   // List of cuts for ESDs
  Float_t      fRadiusCut;                  // radius cut 
  Float_t      fTriggerPtCut;               // cut on particle pt used as event axis
  Float_t      fAssociatePtCut;             // cut on particle pt used for correlations
  Int_t        fLeadingOrRandom;            // event axis:leading track or random track
  Int_t        fMode;                       // ESD(=0) of AOD(=1) reading 
  Float_t      fVertexZCut;                 // vertex cut

  AliESDEvent *fESDEvent;                   //! esd event
  AliAODEvent *fAODEvent;                   //! aod event
  
  TList	      *fHists;                      // output list
  TH1F        *fHistPt;                     // Pt spectrum ESD
  TH1F        *fHistPtMC;                   // Pt spectrum MC
  TH1F       * fVertexZ[4];                 // z of vertex
  TH1F       * fPt[4];                      // pt
  TH1F       * fEta[4];                     // eta
  TH1F       * fPhi[4];                     // phi
  TH2F       * fDPhiDEtaEventAxis[4];       // correlation dEta-dPhi towards event axis
  TH1F       * fTriggerNch[4];              // number of triggers with accepted-track number
  TH1F       * fTriggerTracklet[4];         // number of triggers with accepted-tracklet number
  TH2F       * fNch07Nch[4];                // nCharged with pT>fTriggerPtCut vs nCharged
  TProfile   * pNch07Nch[4];                // nCharged with pT>fTriggerPtCut vs nCharged
  TH2F       * fNch07Tracklet[4];           // nCharged with pT>fTriggerPtCut vs nTracklet
  TH2F       * fNchTracklet[4];             // nCharged vs nTracklet
  TProfile   * pNch07Tracklet[4];           // nCharged with pT>fTriggerPtCut vs nTracklet

  TH1F       * fDPhiEventAxisNchBin[4][150];// delta phi of associate tracks to event axis per Nch bin
  TH1F       * fDPhiEventAxisNchBinTrig[4][150];// "" for all possoble trigger particles

  TH1F       * fDPhiEventAxisTrackletBin[4][150]; // delta phi of associate tracks to event axis per Nch bin
  TH1F       * fDPhiEventAxisTrackletBinTrig[4][150]; // "" for all possible trigger particles

  AliAnalysisTaskMinijet(const AliAnalysisTaskMinijet&); // not implemented
  AliAnalysisTaskMinijet& operator=(const AliAnalysisTaskMinijet&); // not implemented
  
  ClassDef(AliAnalysisTaskMinijet, 1); // example of analysis
};

#endif
