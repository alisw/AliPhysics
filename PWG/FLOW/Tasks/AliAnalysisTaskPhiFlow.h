#ifndef ALIANALYSISTASKPHIFLOW_H
#define ALIANALYSISTASKPHIFLOW_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliAnalysisTaskPhiFlow:
// origin: Redmer Alexander Bertens (rbertens@nikhef.nl)
// analyis task for phi-meson reconstruction and determination of V2

class TH1F;
class TH1F;
class TH2F;
class TProfile;
class AliESDEvent;
class AliESDtrackCuts;
class AliFlowTrackCuts;
class AliFlowEvent;
class AliFlowCandidateTrack;

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


   TH1F*          BookHistogram(const char * name);
   TH2F*          BookPIDHistogram(const char * name);
   TProfile*      BookV2Profile(const char * name, Bool_t pt, Bool_t cos);
   TH1F*          BookPtHistogram(const char* name);
   void           AddPhiIdentificationOutputObjects();
   virtual void   UserCreateOutputObjects();
   void           PairLoss(const AliESDtrack* track1, const AliESDtrack* track2) const;
   Double_t       InvariantMass(const AliESDtrack* track_1, const AliESDtrack* track_2) const;
   Double_t       DeltaDipAngle(const AliESDtrack* track1, const AliESDtrack* track2) const;
   Bool_t         CheckDeltaDipAngle(const AliESDtrack* track1, const AliESDtrack* track2) const;
   void           SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod);
   void           SetVertexZ(Float_t z)  { fVertexRange = z; };
   void           SetMaxDeltaDipAngleAndPt(Float_t a, Float_t pt) {fDeltaDipAngle = a; fDeltaDipPt = pt; fApplyDeltaDipCut = kTRUE; };
   void           EventPlanePtCut(Bool_t a = kFALSE) { fSetEventPlanePtCut = a; };
   void           SetEventPlanePtCut(Float_t c) {fEventPlanePtCut = c; } ;
   Bool_t         EventCut( AliESDEvent* esd) const;
   void           PlotVZeroMultiplcities(const AliESDEvent* esd) const;
   Bool_t         CheckVertex(const AliESDEvent* esd) const ;
   Bool_t         CheckCentrality(AliESDEvent* esd) const;
   Bool_t         IsKaon(const AliESDtrack* track) const;
   Double_t       PhiPt(const AliESDtrack* track_1, const AliESDtrack* track_2) const;
   void           PtSelector(Int_t _track_type, const AliESDtrack* track_1, const AliESDtrack* track_2) const;
   void           EventPlane(const AliESDEvent* esd);
   void           EventPlaneResolution(        Bool_t random,
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
   Double_t       Chi(Double_t res) const;
   Double_t       ResolutionAsFunctionOfChi(Double_t chi) const;
   Double_t       EventPlaneStar(       Bool_t random,
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
   void           EllipticFlow(Double_t* const v2, const AliESDtrack* track_1, const AliESDtrack* track_2) const;
   Bool_t         EventPlaneTrack(AliESDtrack* track) const;
   Bool_t         PhiTrack(AliESDtrack* track) const;
   void           FlowFinish(Double_t* const v2) const;
   void           EllipticFlowSin(Double_t* const v2Sin, const AliESDtrack* track_1, const AliESDtrack* track_2) const;
   void           FlowFinishSin(Double_t* const v2Sin) const;
   void           SetNullCuts(AliESDEvent* esd);
   void           PrepareFlowEvent(Int_t iMulti);
   virtual void   UserExec(Option_t *option);
   virtual void   Terminate(Option_t *);
   void           SetRPCuts(AliFlowTrackCuts *cutsRP) { fRPCuts = cutsRP; }
   void           SetPOICuts(AliFlowTrackCuts *cutsPOI) { fPOICuts = cutsPOI; }
   void           SetKaonCuts(AliFlowTrackCuts *cutsKaon) { fKaonCuts = cutsKaon; }
   void           SetIdentificationType(PIDtype type) {fPIDtype = type;}
   void           SetMassRanges(Double_t flowBands[2][30]) { for (int i = 0; i != 2; ++i) for (int j = 0; j != 30; ++j) fFlowBands[i][j] = flowBands[i][j]; }
   void           SetEtaRanges(Double_t mina, Double_t maxa, Double_t minb, Double_t maxb) { fEtaMinA = mina; fEtaMaxA = maxa; fEtaMinB = minb; fEtaMaxB = maxb; }
   void           SetFlowRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }

private:

   Double_t             fFlowBands[2][30]; // Array containing the boundaries (upper , lower) of invariant mass bands, gev / c^2
   Double_t             fEtaMinA, fEtaMaxA, fEtaMinB, fEtaMaxB; // upper and lower bounds of included eta regions
   AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles (event plane method)
   AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
   AliFlowEvent         *fFlowEvent[30]; //! flow events (one for each inv mass band)
   AliESDEvent          *fESD;    //! ESD object
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
   AliFlowTrackCuts     *fKaonCuts; // dummy flowtrackcuts object for bayesian pid method
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

   AliAnalysisTaskPhiFlow(const AliAnalysisTaskPhiFlow&);
   AliAnalysisTaskPhiFlow& operator=(const AliAnalysisTaskPhiFlow&);

   AliFlowCandidateTrack*  MakeTrack(Double_t, Double_t, Double_t, Double_t, Int_t , Int_t[]) const;

   ClassDef(AliAnalysisTaskPhiFlow, 1);

};

#endif


