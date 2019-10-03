#ifndef ALIANALYSISTASKFLOWITSTPCTOFQCSP_H
#define ALIANALYSISTASKFLOWITSTPCTOFQCSP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron Flow ITS TPC TOF                  //
//                                                                    //
//  Author: Andrea Dubla (Utrecht University)                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidTPC;
class AliHFEpidQAmanager;
class AliCFManager;
class AliFlowTrackCuts;
class AliFlowTrack;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowEventSimple;
class AliCentrality;
class AliPIDResponse;
class AliSelectNonHFE;
#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskFlowITSTPCTOFQCSP : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskFlowITSTPCTOFQCSP();
    AliAnalysisTaskFlowITSTPCTOFQCSP(const char *name);
    virtual ~AliAnalysisTaskFlowITSTPCTOFQCSP();
    
    void                                 SetEnableDebugMode() {fDebug = kTRUE; };
    void                                 SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod); //select centrality
    void                                 CheckCentrality(AliAODEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....
    void                                 SelectPhotonicElectron(Int_t itrack,const AliAODTrack *track,Float_t fTPCnSigma,Double_t evPlAngV0, Bool_t &fFlagPhotonicElec, Double_t weightEPflat, Double_t multev);
    void                                 SelectPhotonicElectronMethod(Bool_t dca){fDCA = dca;}
    void                                 SetInvariantMassCut(Double_t invmass) {fInvmassCut = invmass;};
    void                                 SetPtMinAssoCut(Double_t ptminimumasso) {fptminAsso = ptminimumasso;};
    void                                 SetpTCuttrack(Double_t ptcutmin, Double_t ptcutmax);
    void                                 SetTrigger(Int_t trig) {fTrigger = trig;};
    void                                 SetTPCS(Int_t sig) {fTPCS = sig;};
    void                                 SetVz(Int_t ver) {fVz = ver;};
    //    void                                 SetQAPIDSparse(Bool_t qapidsparse) {fQAPIDSparse = qapidsparse;};
    void                                 SetOpeningAngleflag(Bool_t opang){fOP_angle = opang;};
    void                                 SetOpeningAngleCut(Double_t opanglecut) {fOpeningAngleCut = opanglecut;};
    template <typename T> void           PlotVZeroMultiplcities(const T* event) const;
    template <typename T> void           SetNullCuts(T* aod);
    void                                 PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const;
    virtual void                         UserCreateOutputObjects();
    virtual void                         UserExec(Option_t *option);
    virtual void                         Terminate(Option_t *);
    void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
    void                                 SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void                                 SetIDCuts(Double_t minTOFnSigma, Double_t maxTOFnSigma, Double_t minITSnsigmalowpt, Double_t maxITSnsigmalowpt,Double_t minITSnsigmahighpt, Double_t maxITSnsigmahighpt, Double_t minTPCnsigmalowpt, Double_t maxTPCnsigmalowpt,Double_t minTPCnsigmahighpt, Double_t maxTPCnsigmahighpt );
    
    void                                 SetMultCorrelationCut(Bool_t multcut) {fMultCut = multcut;};
    void                                 SetHistoForCentralityFlattening(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef=0.,Int_t switchTRand=0);
    Bool_t                               IsEventSelectedForCentrFlattening(Float_t centvalue);
    void                                 SetAssoITSRefit(Bool_t itsref) {fAssoITSRefit = itsref;};
    void                                 SetAssoTPCCluster(Int_t tpc_clust) {fAssoTPCCluster = tpc_clust;};
    void                                 SetPhiCut(Bool_t phicut){fPhiCut = phicut;};
    
    void                                 SetHistoForEPFlattWeights(TH1D *h);
    Double_t                             GiveMeWeight(Double_t EP);
    void                                 SetEPWeight(Bool_t epw){EPweights = epw;};
    
    void                                 SetTPCPID(AliHFEpidTPC *pidcorr){ftpcpid = pidcorr;};
    void                                 SetMultCorrectionTheo(Bool_t mulcorr){multCorrection = mulcorr;}
    
    void                                 SetEtaMinPos(Double_t minetapos){fEtaMinimumPositive = minetapos;}
    void                                 SetEtaMinNeg(Double_t minetaneg){fEtaMinimumNegative = minetaneg;}
    
    void                                 SetHistoForCentralityFlattening_Bis(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef);
    Bool_t                               IsEventSelectedForCentrFlattening_Bis(Double_t centvalue);
    
    void                                 SetCentralityMine(Bool_t CentFlatMine){fCentFlatMine = CentFlatMine;}
    void                                 SetEtaRange(Double_t etaminimum, Double_t etamaximum);
    
    
    
    AliHFEpid *GetPID() const { return fPID; };
    
private:
    
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    
    Bool_t               fDebug; //! enable debug mode
    AliAODEvent          *fAOD; //AOD object
    AliVEvent            *fVevent;			//ESD object
    TList                *fOutputList;//output list
    AliHFEcuts           *fCuts; //Cut Collection
    Bool_t               fIdentifiedAsOutInz;//Out Of Range in z
    Bool_t               fPassTheEventCut;//Pass The Event Cut
    AliCFManager         *fCFM;//!Correction Framework Manager
    AliHFEpid            *fPID;//PID
    AliHFEpidTPC         *ftpcpid;// for TPC mult/eta correction, are done in the HFE class
    AliHFEpidQAmanager   *fPIDqa; //! PID QA manager
    AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles
    AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
    AliFlowEvent         *fFlowEvent; //! flow events Inclusive e
    const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)
    Double_t             fCentrality; // event centrality for QA
    Double_t             fCentralityMin; // lower bound of cenrality bin
    Double_t             fCentralityMax; // upper bound of centrality bin
    Double_t             fInvmassCut;//invariant mass cut value
    Double_t             fpTCutmin;//pt cut value min
    Double_t             fpTCutmax; //pt cut value max
    Int_t                fTrigger; //invariant mass cut value
    TH1F                 *fPhi; //! QA plot of azimuthal distribution of tracks used for event plane estimation
    TH1F                 *fEta; //! QA plot of eta distribution of tracks used for event plane estimation
    TH1F                 *fVZEROA; //! QA plot vzeroa multiplicity (all tracks in event)
    TH1F                 *fVZEROC; //! QA plot vzeroc multiplicity (all tracks in event)
    TH1F                 *fTPCM; //! QA plot TPC multiplicity (tracks used for event plane estimation)
    TH1F                 *fNoEvents; //!no of events
    TH1F                 *fInclusiveElecPt; //!Inclusive elec pt
    TH2F                 *fTPCnsigma;//!TPC n sigma vs p
    TH2F                 *fTPCnsigmaAft;//!TPC n sigma vs p after HFE pid
    TH2F                 *fITSnsigmaAft;//!ITS n sigma vs p after HFE pid
    TH2F                 *fTPCnsigmaVSptAft;//!lodviow
    TH2F                 *fTOFns;//!track TOFnSigma
    TH2F                 *fTOFnsAft;//!track TOFnSigma after PID
    TH2F                 *fTOFBetaAft;//!track TOFnSigma after PID2
    TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
    TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
    TH1F                 *fInvmassLS1; //!LS Invmass for all rec par
    TH1F                 *fInvmassULS1;//!ULS Invmass for all rec par
    TH1F                 *fPhotoElecPt;//!photonic elec pt
    TH1F                 *fSemiInclElecPt;//!Semi inclusive ele pt
    TH1F                 *fULSElecPt; //! ULS elec Pt
    TH1F                 *fLSElecPt;//! LS elec pt
    TH2F                 *fMultCorAfterCuts; //! QA profile global and tpc multiplicity after outlier cut
    TH2F                 *fMultvsCentr; //! QA profile of centralty vs multiplicity
    TProfile             *fSubEventDPhiv2;
    TH1D                 *EPVzA;//!ep v0a
    TH1D                 *EPVzC;//!ep v0c
    TH1D                 *EPTPC;//!ep tpc
    THnSparseF           *fV2Phi;//! v2 analysis of EP-V0
    TH1D                 *fvertex;//!vertex
    TH2F                 *fMultCorBeforeCuts; //! QA profile global and tpc multiplicity after outlier cut
    THnSparseF           *fQAPid;             //! v2 analysis of EP-V0
    Double_t             fminITSnsigmaLowpT;  //ID cuts its
    Double_t             fmaxITSnsigmaLowpT;  //ID cuts its
    Double_t             fminITSnsigmaHighpT;  //ID cuts its
    Double_t             fmaxITSnsigmaHighpT;  //ID cuts its
    Double_t             fminTPCnsigmaLowpT;  //ID cuts tpc
    Double_t             fmaxTPCnsigmaLowpT;  //ID cuts tpc
    Double_t             fminTPCnsigmaHighpT;  //ID cuts tpc
    Double_t             fmaxTPCnsigmaHighpT;  //ID cuts tpc
    // Bool_t               fQAPIDSparse;       //QAPIDSPARSE
    Double_t             fminTOFnSigma;  //ID cuts tof
    Double_t             fmaxTOFnSigma;//ID cuts tof
    THnSparseF           *fQAPidSparse;             //! v2 analysis of EP-V0
    Int_t                fTPCS; //tpc signal cluster
    Int_t                fVz; //vertex range
    Double_t             fOpeningAngleCut; //openingAngle cut value
    Bool_t               fOP_angle; //to shitch on and off the op_angle cut
    TH1F                 *fOpeningAngleLS;  //!opening angle for LS pairs
    TH1F                 *fOpeningAngleULS; //!opening angle for ULS pairs
    AliSelectNonHFE      *fNonHFE;//new elienos stuff
    Bool_t               fDCA;//selection PHelectron
    TH2F                 *fITSnsigma;//!TPC n sigma vs p
    TH2F                 *fITSnsigmaElect;//!TPC n sigma vs p
    TH2F                 *fTPCnsigmaAftITSTOF; //!TPC n sigma vs p
    TH2F                 *fTPCnsigmaAftTOF; //!TPC n sigma vs p
    TH2F                 *fITSnsigmaAftTOF;//!jbd
    TH2F                 *fITSvsTOF;//!TPC n sigma vs p
    TH2F                 *fTPCvsITS;//!TPC n sigma vs p
    TH2F                 *fTPCvsITSafterTOF;//!TPC n sigma vs p
    TH2F                 *fTPCvsTOF;//!jbd
    Bool_t               fMultCut;//for mult correlationcut
    TH2F                 *fMultCorAfterCentrBeforeCuts; //! QA profile global and tpc multiplicity after outlier cut
    TH2F                 *fMultCorAfterVZTRKComp;//!after cent comp
    TH1F                 *fCentralityBeforePileup;//!cent chneck
    TH1F                 *fCentralityAfterVZTRK;//!cent chneck2
    TH1F                 *fCentralityAfterCorrCut;//!cent chneck2
    TH2F                 *fMultCorAfterCorrCut;//!after cent comp
    TH1D                 *EPVz;//!v0cep
    TH1D                 *EPTPCp;//!tpcep
    TH1D                 *EPTPCn;//!tpcen
    TProfile             *fSubEventDPhiv2new;//!evrr
    THnSparseF           *fV2Phivzerotot;//! v2 analysis of EP-V0
    TH1F                 *fHistCentrDistr;//-> isto for Centr Flat
    TH1F                 *fCentralityNoPassForFlattening; //!QA histogram of events that do not pass centrality cut for flattening
    Double_t              fptminAsso;//minassopt
    THnSparseF           *fSparsephipsiULS;//! ULSSparse
    THnSparseF           *fSparsephipsiLS;//! LSSparse
    THnSparseF           *fSparseMassULS;//!ssss
    THnSparseF           *fSparseMassLS;//!ssssss
    Int_t                 fAssoTPCCluster;//asso tpc cluster
    Bool_t                fAssoITSRefit;//asso its refit
    Bool_t                fPhiCut;//Phi cut to simulate emcal acc
    
    TH1D                 *fHistEPDistrWeight;// isto for Centr Flat
    Bool_t               EPweights;//for mult correlationcut
    TH1D                 *EPVzAftW;//!v0cep
    
    Bool_t                multCorrection;//Flag to activate mult/etacorrection
    Double_t              fEtaMinimumPositive;//for reso
    Double_t              fEtaMinimumNegative;//for reso
    TH1F                 *fCentralityAll;//!centall
    Bool_t               fCentFlatMine; //for purity evaluation
    Double_t             fmineta;  //eta cuts
    Double_t             fmaxeta;//eta cuts

    
    AliAnalysisTaskFlowITSTPCTOFQCSP(const AliAnalysisTaskFlowITSTPCTOFQCSP&); // not implemented
    AliAnalysisTaskFlowITSTPCTOFQCSP& operator=(const AliAnalysisTaskFlowITSTPCTOFQCSP&); // not implemented
    
    ClassDef(AliAnalysisTaskFlowITSTPCTOFQCSP, 2); //!example of analysis
};

#endif




