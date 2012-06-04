/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliAnalysisTaskPhiFlow:
// origin: Redmer Alexander Bertens (rbertens@nikhef.nl)
// analyis task for phi-meson reconstruction and determination of V2
// handles aod's and esd's transparantly

#ifndef ALIANALYSISTASKPHIFLOW_H
#define ALIANALYSISTASKPHIFLOW_H

class TH1F;
class TH1F;
class TH2F;
class TProfile;
class AliESDEvent;
class AliESDtrackCuts;
class AliFlowTrackCuts;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowBayesianPID;

#include "AliAnalysisTaskSE.h"

enum ESDEventStats_t
{
   kNUnlikePairs = 0,
   kNLikeNPairs,
   kNLikePPairs,
   kNUnlikeKPairs,
   kNLikeNKPairs,
   kNLikePKPairs,
   kNStats = kNLikePKPairs,
};

class AliAnalysisTaskPhiFlow : public AliAnalysisTaskSE
{
public:
   AliAnalysisTaskPhiFlow();
   AliAnalysisTaskPhiFlow(const char *name);
   virtual ~AliAnalysisTaskPhiFlow();

   void                                 SetEnableDebugMode() {fDebug = kTRUE; };
   TH1F*                                BookHistogram(const char * name);
   TH2F*                                BookPIDHistogram(const char * name);
   TH1F*                                InitPtSpectraHistograms(Int_t i);
   TH1F*                                BookPtHistogram(const char* name);
   void                                 AddPhiIdentificationOutputObjects();
   virtual void                         UserCreateOutputObjects();
   template <typename T> Double_t       InvariantMass(const T* track1, const T* track2) const;
   template <typename T> Double_t       DeltaDipAngle(const T* track1, const T* track2) const;
   template <typename T> Bool_t         CheckDeltaDipAngle(const T* track1, const T* track2) const;
   template <typename T> Bool_t         CheckCandidateEtaPtCut(const T* track1, const T* track2) const;
   void                                 SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod);
   void                                 SetVertexZ(Float_t z) { fVertexRange = z; };
   void                                 SetMaxDeltaDipAngleAndPt(Float_t a, Float_t pt) { fDeltaDipAngle = a;
                                                                                          fDeltaDipPt = pt;
                                                                                          fApplyDeltaDipCut = kTRUE; };
   template <typename T> Bool_t         EventCut(T* event);
   template <typename T> void           PlotVZeroMultiplcities(const T* event) const;
   template <typename T> Bool_t         CheckVertex(const T* event) const ;
   template <typename T> Bool_t         CheckCentrality(T* event);
   void                                 InitializeBayesianPID(AliESDEvent* event);
   void                                 InitializeBayesianPID(AliAODEvent* event);
   template <typename T> Bool_t         PassesTPCbayesianCut(T* track) const;
   Bool_t                               PassesStrictKaonCuts(AliESDtrack* track) const;
   Bool_t                               PassesStrictKaonCuts(AliAODTrack* track) const;
   Bool_t                               IsKaon(AliESDtrack* track) const;
   Bool_t                               IsKaon(AliAODTrack* track) const;
   template <typename T> Double_t       PhiPt(const T* track_1, const T* track_2) const;
   template <typename T> void           PtSelector(Int_t _track_type, const T* track_1, const T* track_2) const;
   template <typename T> Bool_t         PhiTrack(T* track) const;
   template <typename T> void           SetNullCuts(T* esd);
   void                                 PrepareFlowEvent(Int_t iMulti);
   virtual void                         UserExec(Option_t *option);
   virtual void                         Terminate(Option_t *);
   void                                 SetPOICuts(AliFlowTrackCuts *cutsPOI) { fPOICuts = cutsPOI; }
   void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
   void                                 SetRequireStrictKaonCuts() { fStrictKaonCuts= kTRUE; }
   void                                 SetRequireTPCStandAloneKaons() { fRequireTPCStandAlone = kTRUE; }
   void                                 SetOlderThanPbPbPass2() { fOldTrackParam = kTRUE; }
   void                                 SetPIDConfiguration(Double_t prob[7]) { for(Int_t i = 0; i < 7; i++) fPIDConfig[i] = prob[i]; }
   void                                 SetPOIDCAXYZ(Double_t a, Double_t b) {fDCAXY = a; fDCAZ = b; fDCA = kTRUE; }
   void                                 SetCandidateEtaAndPt(Double_t mineta, Double_t maxeta, Double_t minpt, Double_t maxpt) { fCandidateMinEta = mineta, 
                                                                                                                               fCandidateMaxEta = maxeta,
                                                                                                                               fCandidateMinPt = minpt, 
                                                                                                                               fCandidateMaxPt = maxpt,
                                                                                                                               fCandidateEtaPtCut = kTRUE;}
   void                                 SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass);

private:

   Bool_t               fDebug; //! enable debug mode
   Bool_t               fAODAnalysis; // set aod analysis
   Int_t                fMassBins; // mass bins
   Double_t             fMinMass; // mass range
   Double_t             fMaxMass; // mass range
   AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles
   AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
   AliPIDResponse       *fPIDResponse; //! pid response object
   AliFlowEvent         *fFlowEvent; //! flow events (one for each inv mass band)
   AliFlowBayesianPID   *fBayesianResponse; //!PID response object
   TObjArray            *fCandidates; // candidate array
   Bool_t               fOldTrackParam; // set to true if data is older than pbpb pass 2 production
   Bool_t               fRequireTPCStandAlone; // set TPC standalone cut for kaon selection
   Bool_t               fStrictKaonCuts; // require strict kaon cuts
   Bool_t               fCandidateEtaPtCut; // set eta and pt cut for candidate tracks and combinatorial background
   Double_t             fCandidateMinEta; // minimum eta for candidates
   Double_t             fCandidateMaxEta; // maximum eta for candidates
   Double_t             fCandidateMinPt; // minimum pt for candidates
   Double_t             fCandidateMaxPt; // maximum pt for candidates
   Double_t             fPIDConfig[7]; // set cutoff for bayesian probability
   Double_t             fCentrality; // event centrality
   AliESDEvent          *fESD;    //! ESD object
   AliAODEvent          *fAOD;    //! AOD oject
   TList                *fOutputList; // ! Output list
   TH1F                 *fEventStats; // ! Histogram for event statistics
   TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
   TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
   TH2F                 *fNOPID;//! QA histogram of TPC response of all charged particles
   TH2F                 *fPIDk;//! QA histogram of TPC response of kaons
   TH1F                 *fInvMNP03; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP03; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN03; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP36; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP36; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN36; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP69; //!  Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP69; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN69; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP912; //!  Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP912; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN912; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP1215; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP1215; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN1215; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP1518; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP1518; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN1518; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP1821; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP1821; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN1821; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP2124; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP2124; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN2124; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP2427; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP2427; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN2427; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP2730; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP2730; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN2730; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP3035; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP3035; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN3035; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP3540; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP3540; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN3540; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP4045; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP4045; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN4045; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP4550; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP4550; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN4550; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP5055; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP5055; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN5055; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP5560; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP5560; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN5560; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP6065; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP6065; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN6065; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fInvMNP6570; //! Invariant mass of unlike sign kaon pairs
   TH1F                 *fInvMPP6570; //! Invariant mass of like sign (++) kaon pairs
   TH1F                 *fInvMNN6570; //! Invariant mass of like sign (--) kaon pairs
   TH1F                 *fPtSpectra03; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra36; //! kaon pair Pt spectrum 
   TH1F                 *fPtSpectra69; //!  kaon pair Pt spectrum
   TH1F                 *fPtSpectra912; //!  kaon pair Pt spectrum
   TH1F                 *fPtSpectra1215; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra1518; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra1821; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra2124; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra2427; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra2730; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra3035; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra3540; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra4045; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra4550; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra5055; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra5560; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra6065; //! kaon pair Pt spectrum
   TH1F                 *fPtSpectra6570; //! kaon pair Pt spectrum
   TH1F                 *fPtP; //! QA histogram of p_t distribution of positive particles
   TH1F                 *fPtN; //! QA histogram of p_t distribution of negative particles
   TH1F                 *fPtKP; //! QA histogram of p_t distribution of positive kaons
   TH1F                 *fPtKN; //! QA histogram of p_t distribution of negative kaons
   Double_t             fCentralityMin; // lower bound of cenrality bin
   Double_t             fCentralityMax; // upper bound of centrality bin
   const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)
   AliFlowTrackCuts     *fPOICuts; // cuts for particles of interest (flow package)
   Float_t              fVertexRange; // absolute value of maximum distance of vertex along the z-axis
   TH1F                 *fPhi; //! QA plot of azimuthal distribution of tracks used for event plane estimation
   TH1F                 *fPt; //! QA plot of p_t sectrum of tracks used for event plane estimation
   TH1F                 *fEta; //! QA plot of eta distribution of tracks used for event plane estimation
   TH1F                 *fVZEROA; //! QA plot vzeroa multiplicity (all tracks in event)
   TH1F                 *fVZEROC; //! QA plot vzeroc multiplicity (all tracks in event)
   TH1F                 *fTPCM; //! QA plot TPC multiplicity (tracks used for event plane estimation)
   Float_t              fDeltaDipAngle; // absolute value of delta dip angle to be excluded
   Float_t              fDeltaDipPt; // upper value of pt range in which delta dip angle must be applied
   Bool_t               fApplyDeltaDipCut; // enforce delta dip cut
   Double_t             fDCAXY; // POI DCA XY
   Double_t             fDCAZ; // POI DCA Z
   Bool_t               fDCA; // force propagation of DCA for POI's AOD ONLY
   TH1F                 *fDCAXYQA; //! qa plot of dca xz
   TH1F                 *fDCAZQA; //!qa plot of dca z

   AliAnalysisTaskPhiFlow(const AliAnalysisTaskPhiFlow&); // Not implemented
   AliAnalysisTaskPhiFlow& operator=(const AliAnalysisTaskPhiFlow&); // Not implemented

   void                 MakeTrack(Double_t, Double_t, Double_t, Double_t, Int_t , Int_t[]) const;
     ClassDef(AliAnalysisTaskPhiFlow, 3);

};

#endif



