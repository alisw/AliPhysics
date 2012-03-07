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

   enum PIDtype
   {
      kTPC = 0,
      kCombined,
   };

   TH1F*                                BookHistogram(const char * name);
   TH1F*                                BookDPhiPsiHistogram(const char * name);
   TH2F*                                BookPIDHistogram(const char * name);
   TH1F*                                InitPtSpectraHistograms(Int_t i);
   TProfile*                            BookV2Profile(const char * name, Bool_t pt, Bool_t cos);
   TH1F*                                BookPtHistogram(const char* name);
   void                                 AddPhiIdentificationOutputObjects();
   virtual void                         UserCreateOutputObjects();
   template <typename T> void           PairLoss(const T* track1, const T* track2) const;
   template <typename T> Double_t       InvariantMass(const T* track1, const T* track2) const;
   template <typename T> Double_t       DeltaDipAngle(const T* track1, const T* track2) const;
   template <typename T> Bool_t         CheckDeltaDipAngle(const T* track1, const T* track2) const;
   template <typename T> Bool_t         CheckCandidateEtaPtCut(const T* track1, const T* track2) const;
   void                                 SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod);
   void                                 SetVertexZ(Float_t z) { fVertexRange = z; };
   void                                 SetMaxDeltaDipAngleAndPt(Float_t a, Float_t pt) { fDeltaDipAngle = a;
                                                                                          fDeltaDipPt = pt;
                                                                                          fApplyDeltaDipCut = kTRUE; };
   void                                 EventPlanePtCut(Bool_t a = kFALSE) { fSetEventPlanePtCut = a; };
   void                                 SetEventPlanePtCut(Float_t c) { fEventPlanePtCut = c; };
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
   void                                 EventPlane(const AliESDEvent* event);
   void                                 EventPlane(const AliAODEvent* event);
   void                                 EventPlaneResolution(Bool_t random,
                                                        Double_t cossumQ1,
                                                        Double_t sinsumQ1,
                                                        Double_t cossumQ2,
                                                        Double_t sinsumQ2,
                                                        Double_t cossum,
                                                        Double_t sinsum,
                                                        Double_t cossumQwest,
                                                        Double_t sinsumQwest,
                                                        Double_t cossumQeast,
                                                        Double_t sinsumQeast) const;
   Double_t                             Chi(Double_t res) const;
   Double_t                             ResolutionAsFunctionOfChi(Double_t chi) const;
   Double_t                             EventPlaneStar(Bool_t random,
                                                        Double_t cossumQ1,
                                                        Double_t sinsumQ1,
                                                        Double_t cossumQ2,
                                                        Double_t sinsumQ2,
                                                        Double_t cossum,
                                                        Double_t sinsum,
                                                        Double_t cossumQwest,
                                                        Double_t sinsumQwest,
                                                        Double_t cossumQeast,
                                                        Double_t sinsumQeast) const;
   template <typename T> void           EllipticFlow(Double_t* const v2, const T* track_1, const T* track_2) const;
   template <typename T> Bool_t         EventPlaneTrack(T* track) const;
   template <typename T> Bool_t         PhiTrack(T* track) const;
   void                                 FlowFinish(Double_t* const v2) const;
   template <typename T> void           EllipticFlowSin(Double_t* const v2Sin, const T* track_1, const T* track_2) const;
   void                                 FlowFinishSin(Double_t* const v2Sin) const;
   template <typename T> void           SetNullCuts(T* esd);
   void                                 PrepareFlowEvent(Int_t iMulti);
   virtual void                         UserExec(Option_t *option);
   virtual void                         Terminate(Option_t *);
   void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fRPCuts = cutsRP; }
   void                                 SetPOICuts(AliFlowTrackCuts *cutsPOI) { fPOICuts = cutsPOI; }
   void                                 SetIdentificationType(PIDtype type) { fPIDtype = type; }
   void                                 SetRequireStrictKaonCuts() { fStrictKaonCuts= kTRUE; }
   void                                 SetRequireTPCStandAloneKaons() { fRequireTPCStandAlone = kTRUE; }
   void                                 SetOlderThanPbPbPass2() { fOldTrackParam = kTRUE; }
   void                                 SetBayesianProbability(Double_t prob) { fParticleProbability = prob; }
   void                                 SetMassRanges(Double_t flowBands[2][30]) { for (int i = 0; i != 2; ++i) for (int j = 0; j != 30; ++j) fFlowBands[i][j] = flowBands[i][j]; }
   void                                 SetEtaRanges(Double_t mina, Double_t maxa, Double_t minb, Double_t maxb) { 
                                                                                        fEtaMinA = mina;
                                                                                        fEtaMaxA = maxa; 
                                                                                        fEtaMinB = minb; 
                                                                                        fEtaMaxB = maxb; }
   void                                 SetFlowRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
   void                                 SetCandidateEtaAndPt(Double_t mineta, Double_t maxeta, Double_t minpt, Double_t maxpt) { fCandidateMinEta = mineta, 
                                                                                                                               fCandidateMaxEta = maxeta,
                                                                                                                               fCandidateMinPt = minpt, 
                                                                                                                               fCandidateMaxPt = maxpt,
                                                                                                                               fCandidateEtaPtCut = kTRUE;}

private:

   Bool_t               fAODAnalysis; // set aod analysis
   Double_t             fFlowBands[2][30]; // Array containing the boundaries (upper , lower) of invariant mass bands, gev / c^2
   Double_t             fEtaMinA, fEtaMaxA, fEtaMinB, fEtaMaxB; // upper and lower bounds of included eta regions
   AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles (event plane method)
   AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
   AliFlowEvent         *fFlowEvent[30]; //! flow events (one for each inv mass band)
   AliFlowBayesianPID   *fBayesianResponse; //!PID response object
   Bool_t               fOldTrackParam; // set to true if data is older than pbpb pass 2 production
   Bool_t               fRequireTPCStandAlone; // set TPC standalone cut for kaon selection
   Bool_t               fStrictKaonCuts; // require strict kaon cuts
   Bool_t               fCandidateEtaPtCut; // set eta and pt cut for candidate tracks and combinatorial background
   Double_t             fCandidateMinEta; // minimum eta for candidates
   Double_t             fCandidateMaxEta; // maximum eta for candidates
   Double_t             fCandidateMinPt; // minimum pt for candidates
   Double_t             fCandidateMaxPt; // maximum pt for candidates
   Double_t             fParticleProbability; // set cutoff for bayesian probability
   Double_t             fCentrality; // event centrality
   AliESDEvent          *fESD;    //! ESD object
   AliAODEvent          *fAOD;    //! AOD oject
   TList                *fOutputList; // ! Output list
   TH1F                 *fEventStats; // ! Histogram for event statistics
   TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
   TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
   TH2F                 *fNOPID;//! QA histogram of TPC response of all charged particles
   TH2F                 *fPIDk;//! QA histogram of TPC response of kaons
   TH2F                 *fPIDDeltaDip; //! QA histogram of TPC signal after delta dip cut
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
   TH1F                 *fDeltaPhiPsiNP03; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP36; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP69; //!  Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP912; //!  Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP1215; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP1518; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP1821; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP2124; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP2427; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP2730; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP3035; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP3540; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP4045; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP4550; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP5055; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP5560; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP6065; //! Delta Phi Psi of unlike sign kaon pairs
   TH1F                 *fDeltaPhiPsiNP6570; //! Delta Phi Psi of unlike sign kaon pairs
   TProfile             *fProfV2; //! charged particle v2 (cos terms, event plane method)
   TProfile             *fProfV2Sin; //! charged particle v2 (sin terms, event plane method)
   TProfile             *fProfV2InvM03; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM36; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM69; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM912; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM1215; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM1518; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM1821; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM2124; //! cos component of v2 (eevent plane method)
   TProfile             *fProfV2InvM2427; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM2730; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM3035; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM3540; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM4045; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM4550; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM5055; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM5560; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM6065; //! cos component of v2 (event plane method)
   TProfile             *fProfV2InvM6570; //! cos component of v2 (event plane method)
   TProfile             *fProfV2SinInvM03; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM36; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM69; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM912; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM1215; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM1518; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM1821; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM2124; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM2427; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM2730; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM3035; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM3540; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM4045; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM4550; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM5055; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM5560; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM6065; //! sin component of v2 (event plane method)
   TProfile             *fProfV2SinInvM6570; //! sin component of v2 (event plane method)
   TH1F                 *fPtSpectra03; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra36; //! Kaon Pt spectrum 
   TH1F                 *fPtSpectra69; //!  Kaon Pt spectrum
   TH1F                 *fPtSpectra912; //!  Kaon Pt spectrum
   TH1F                 *fPtSpectra1215; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra1518; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra1821; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra2124; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra2427; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra2730; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra3035; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra3540; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra4045; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra4550; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra5055; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra5560; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra6065; //! Kaon Pt spectrum
   TH1F                 *fPtSpectra6570; //! Kaon Pt spectrum
   TH1F                 *fEventPlaneSTAR; //! Distribution of the orientation of the event plane, taken from STAR
   TProfile             *fEventPlaneResolutionRandom; //! event plane resolution from randomized subevents (event plane method)
   TProfile             *fEventPlaneResolutionEta; //! event plane resolution from eta separated subevents (event plane method)
   TH1F                 *fPtP; //! QA histogram of p_t distribution of positive particles
   TH1F                 *fPtN; //! QA histogram of p_t distribution of negative particles
   TH1F                 *fPtKP; //! QA histogram of p_t distribution of positive kaons
   TH1F                 *fPtKN; //! QA histogram of p_t distribution of negative kaons
   Double_t             fCentralityMin; // lower bound of cenrality bin
   Double_t             fCentralityMax; // upper bound of centrality bin
   const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)
   AliFlowTrackCuts     *fRPCuts; // cuts for reference paricles (flow package)
   AliFlowTrackCuts     *fPOICuts; // cus for particles of interest (flow package)
   Float_t              fVertexRange; // absolute value of maximum distance of vertex along the z-axis
   Double_t             fQx; // dummy variable for the Q_x vector (event plane method)
   Double_t             fQy; // dummy variable for the Q_y vector (event plane method)
   TH1F                 *fEventPlane; //! Distribution of the orientation of the event palen (event plane method)
   TH1F                 *fPhi; //! QA plot of azimuthal distribution of tracks used for event plane estimation
   TH1F                 *fPt; //! QA plot of p_t sectrum of tracks used for event plane estimation
   TH1F                 *fEta; //! QA plot of eta distribution of tracks used for event plane estimation
   TH1F                 *fVZEROA; //! QA plot vzeroa multiplicity (all tracks in event)
   TH1F                 *fVZEROC; //! QA plot vzeroc multiplicity (all tracks in event)
   TH1F                 *fTPCM; //! QA plot TPC multiplicity (tracks used for event plane estimation)
   PIDtype              fPIDtype; // selector to indicate which PID process will be used
   Float_t              fDeltaDipAngle; // absolute value of delta dip angle to be excluded
   Float_t              fDeltaDipPt; // upper value of pt range in which delta dip angle must be applied
   Bool_t               fApplyDeltaDipCut; // enforce delta dip cut
   TH2F                 *fPairLoss; //! pair loss histo
   Float_t              fEventPlanePtCut; // Pt cut on tracks used for event plane estimation
   Bool_t               fSetEventPlanePtCut; // kTRUE for Pt cut on event plane estimation tracks

   AliAnalysisTaskPhiFlow(const AliAnalysisTaskPhiFlow&); // Not implemented
   AliAnalysisTaskPhiFlow& operator=(const AliAnalysisTaskPhiFlow&); // Not implemented

   AliFlowCandidateTrack*  MakeTrack(Double_t, Double_t, Double_t, Double_t, Int_t , Int_t[]) const;

   ClassDef(AliAnalysisTaskPhiFlow, 2);

};

#endif


